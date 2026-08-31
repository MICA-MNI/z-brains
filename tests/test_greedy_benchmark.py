"""Tests for the greedy benchmark foundation (objective + greedy loop).

Covers Step 1 (greedy loop locks the best candidate, no processing) and Step 2
(balanced dataset x disease x feature objective). Heavy processing is monkeypatched
out so the loop logic is tested in isolation.

Run: python tests/test_greedy_benchmark.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
import types
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd

import zbrains_staged as st
from results.greedy_auroc import balanced_macro_auroc


# --- Step 2: balanced objective ---------------------------------------------
def test_balanced_objective_weights_cells_equally():
    # 2 diseases x 2 features with very unequal row counts; a flat mean is
    # dominated by FCD/T1w (many rows), the balanced macro is not.
    rows = []
    rows += [{"dataset": "MICs", "lesion_type": "FCD", "feature": "T1w", "auc": 0.50}] * 30
    rows += [{"dataset": "MICs", "lesion_type": "FCD", "feature": "FLAIR", "auc": 0.75}] * 2
    rows += [{"dataset": "MICs", "lesion_type": "TLE", "feature": "T1w", "auc": 0.75}] * 2
    rows += [{"dataset": "MICs", "lesion_type": "TLE", "feature": "FLAIR", "auc": 0.80}] * 2
    detail = pd.DataFrame(rows)
    flat = float(detail["auc"].mean())                       # ~0.544, dominated by the big FCD/T1w cell
    balanced = balanced_macro_auroc(detail, group_cols=("dataset", "lesion_type", "feature"))
    assert abs(balanced - np.mean([0.50, 0.75, 0.75, 0.80])) < 1e-9  # cell-mean-of-means = 0.70
    assert abs(balanced - flat) > 0.05                       # balanced not dominated by the large cell
    # feature restriction works
    only_t1w = balanced_macro_auroc(detail, features=["T1w"], group_cols=("dataset", "lesion_type", "feature"))
    assert abs(only_t1w - np.mean([0.50, 0.75])) < 1e-9


# --- Step 1: greedy loop locks the best candidate (no processing) -----------
def test_greedy_loop_locks_best(monkeypatch):
    # Stub out all heavy processing; the loop should still lock the arg-max.
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: f"/fake/{config['normalization']}",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    stages = [st.Stage("subject_norm", [
        ("whitestripe", {"subject_norm": "whitestripe"}),
        ("wmmean", {"subject_norm": "wmmean"}),
    ])]
    scores = {"none": 0.50, "whitestripe": 0.70, "wmmean": 0.60}

    def mock_score(config, output_dirs):
        return scores[config["normalization"]]

    best, history = st.run_staged_optimization(
        cohorts=[fake_cohort], env=None, stages=stages,
        score_fn=mock_score, verbose=False,
    )
    assert best["normalization"] == "whitestripe"
    winners = history[history["is_winner"] == True]  # noqa: E712
    assert (winners["candidate"] == "WINNER:whitestripe").any()


def test_greedy_loop_keeps_baseline_when_no_improvement(monkeypatch):
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: f"/fake/{config['normalization']}",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    stages = [st.Stage("subject_norm", [
        ("whitestripe", {"subject_norm": "whitestripe"}),
        ("wmmean", {"subject_norm": "wmmean"}),
    ])]
    scores = {"none": 0.90, "whitestripe": 0.70, "wmmean": 0.60}  # baseline best

    best, history = st.run_staged_optimization(
        cohorts=[fake_cohort], env=None, stages=stages,
        score_fn=lambda c, d: scores[c["normalization"]], verbose=False,
    )
    assert best["normalization"] == "none"  # implicit "keep" wins


def test_greedy_resume_starts_after_last_locked_stage(monkeypatch):
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: f"/fake/{config['method']}",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    subject_stage = st.Stage("subject_norm", [
        ("whitestripe", {"subject_norm": "whitestripe"}),
        ("wmmean", {"subject_norm": "wmmean"}),
    ])
    scoring_stage = st.Stage("scoring", [
        ("gp", {"method": "wscore", "wscore_covariate_model": "gaussian_process"}),
    ])

    with tempfile.TemporaryDirectory() as td:
        history_path = Path(td) / "history.csv"

        st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[subject_stage],
            score_fn=lambda c, d: {"none": 0.50, "whitestripe": 0.70, "wmmean": 0.60}[c["subject_norm"]],
            verbose=False, history_path=history_path,
        )

        calls = []

        def resumed_score(config, output_dirs):
            calls.append(dict(config))
            assert config["subject_norm"] == "whitestripe"
            assert config["method"] == "wscore"
            return 0.80

        best, history = st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[subject_stage, scoring_stage],
            score_fn=resumed_score, verbose=False, history_path=history_path,
        )

    assert len(calls) == 1
    assert best["subject_norm"] == "whitestripe"
    assert best["method"] == "wscore"
    assert (history["candidate"] == "WINNER:gp").any()


def test_greedy_resume_skips_recorded_mid_stage_candidate(monkeypatch):
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: f"/fake/{config['subject_norm']}",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    stage = st.Stage("subject_norm", [
        ("whitestripe", {"subject_norm": "whitestripe"}),
        ("wmmean", {"subject_norm": "wmmean"}),
    ])

    baseline, _ = st.validate(dict(st.SIMPLEST_BASELINE))
    whitestripe, _ = st.validate(dict(st.SIMPLEST_BASELINE, subject_norm="whitestripe"))

    with tempfile.TemporaryDirectory() as td:
        history_path = Path(td) / "history.csv"
        partial_history = [
            st._history_row("baseline", "baseline", baseline, 0.50, 0.5, is_winner=True),
            st._history_row("subject_norm", "whitestripe", whitestripe, 0.70, 0.5),
        ]
        st._write_history(partial_history, history_path)

        calls = []

        def resumed_score(config, output_dirs):
            calls.append(config["subject_norm"])
            assert config["subject_norm"] != "whitestripe"  # proves skip, not rerun
            return 0.60

        best, history = st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[stage],
            score_fn=resumed_score, verbose=False, history_path=history_path,
        )

    assert calls == ["wmmean"]
    assert best["subject_norm"] == "whitestripe"
    assert (history["candidate"] == "WINNER:whitestripe").any()


# --- Step 3: T1w/FLAIR sample-surface + self-norm resolver -------------------
def test_t1w_flair_resolver_maps_to_existing_models():
    from zbrains import analysis as zba
    for surf in ["midthickness", "white", "swm1"]:
        for sn, suffix in [(None, ""), ("swm2", "_swm2"), ("owncortex", "_owncortex")]:
            cfg = dict(st.SIMPLEST_BASELINE,
                       sample_surface=surf, t1w_flair_self_norm=sn)
            merged, _ = st.validate(cfg)
            idm = merged["intensity_depth_model"]
            assert idm == f"sample_{surf}{suffix}", (surf, sn, idm)
            assert idm in zba._INTENSITY_DEPTH_TRANSFORMS  # a real, dispatchable transform
    # baseline uses the white surface by default
    merged, _ = st.validate(dict(st.SIMPLEST_BASELINE))
    assert merged["sample_surface"] == "white"
    assert merged["intensity_depth_model"] == "sample_white"
    assert merged["quant_sample_surface"] == "white"
    raw, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="raw"))
    assert raw["intensity_depth_model"] == "raw"
    assert raw["quant_sample_surface"] is None
    # the driver keys are NOT zbrains axes (stripped before validate_configuration)
    axis = st._axis_config(dict(st.SIMPLEST_BASELINE))
    assert "sample_surface" not in axis and "t1w_flair_self_norm" not in axis
    assert "cortical_smoothing" not in axis and "hippocampal_smoothing" not in axis


def test_t1w_flair_resolver_distinct_cache_keys():
    # sampling surface and self-norm produce distinct config-cache keys (via idm)
    base = dict(st.SIMPLEST_BASELINE)
    a, _ = st.validate(dict(base, sample_surface="white"))
    b, _ = st.validate(dict(base, sample_surface="swm1"))
    c, _ = st.validate(dict(base, sample_surface="white", t1w_flair_self_norm="swm2"))
    keys = {st._config_key(a, 0.5), st._config_key(b, 0.5), st._config_key(c, 0.5)}
    assert len(keys) == 3


def test_smoothing_distinct_cache_and_processing_signature():
    base, _ = st.validate(dict(st.SIMPLEST_BASELINE))
    coarse, _ = st.validate(dict(
        st.SIMPLEST_BASELINE,
        cortical_smoothing=10,
        hippocampal_smoothing=5,
    ))
    fine, _ = st.validate(dict(
        st.SIMPLEST_BASELINE,
        cortical_smoothing=2,
        hippocampal_smoothing=1,
    ))
    keys = {st._config_key(base, 0.5), st._config_key(coarse, 0.5), st._config_key(fine, 0.5)}
    assert len(keys) == 3
    # The smoothing-convention tag (smfwhm) is ALWAYS present so pre-`-fwhm`
    # (sigma-smoothed) bases can never be reused; the smoothing pair is also
    # always appended so the 0/0 baseline cannot reuse old no-token 5/2 bases.
    conv = st._SMOOTHING_CONVENTION
    assert st._processing_signature(base) == f"{conv}_smoothctx0hip0"
    assert st._processing_signature(coarse) == f"{conv}_smoothctx10hip5"
    assert st._processing_signature(fine) == f"{conv}_smoothctx2hip1"


def test_default_scoring_stage_has_kernels_but_no_uncertainty_arms():
    scoring = [stage for stage in st.DEFAULT_STAGES if stage.name == "scoring"][0]
    labels = [label for label, _overrides in scoring.candidates]
    assert "zscore" in labels
    assert "wscore_demographics" in labels
    assert "knn_half_controls" in labels
    assert "gp_rbf_ard" in labels
    assert "gp_matern32_ard" in labels
    assert "gp_matern52_ard" in labels
    assert "gp_rational_quadratic_ard" in labels
    # uncertainty arms moved OUT of scoring into the conditional gp_uncertainty stage
    assert not any("uncertainty" in label for label in labels)

    by_label = dict(scoring.candidates)
    assert by_label["zscore"]["method"] == "zscore"
    assert by_label["gp_matern32_ard"]["wscore_covariate_model"] == "gaussian_process_matern32"
    assert by_label["gp_matern52_ard"]["wscore_covariate_model"] == "gaussian_process_matern52"
    assert by_label["gp_rational_quadratic_ard"]["wscore_covariate_model"] == "gaussian_process_rational_quadratic"


def test_gp_uncertainty_stage_is_conditional_and_kernel_agnostic():
    stage = [s for s in st.DEFAULT_STAGES if s.name == "gp_uncertainty"][0]
    # comes AFTER scoring
    names = [s.name for s in st.DEFAULT_STAGES]
    assert names.index("gp_uncertainty") == names.index("scoring") + 1
    # candidates set ONLY the variance percentile (they inherit the winning kernel)
    for label, overrides in stage.candidates:
        assert set(overrides) == {"prediction_variance_percentile"}
        assert "uncertainty_p" in label
    pcts = {ov["prediction_variance_percentile"] for _l, ov in stage.candidates}
    assert pcts == set(st.GPR_UNCERTAINTY_PERCENTILES)
    # condition: runs for any GP kernel, skips for zscore/wscore-linear/knn
    cond = stage.condition
    for model in st.GPR_KERNEL_OPTIONS.values():
        assert cond({"method": "wscore", "wscore_covariate_model": model}) is True
    assert cond({"method": "zscore", "wscore_covariate_model": "linear"}) is False
    assert cond({"method": "wscore", "wscore_covariate_model": "linear"}) is False
    assert cond({"method": "wscore", "wscore_covariate_model": "knn"}) is False


def test_gp_uncertainty_stage_skipped_when_non_gp_wins(monkeypatch):
    # end-to-end: when the scoring winner is NOT a GP kernel, the conditional
    # gp_uncertainty stage is skipped (its winner row is a keep-skip, no arms run).
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: "/fake",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    ran = []

    def mock_score(config, output_dirs):
        # reward w-score (linear) so a NON-GP model wins scoring
        if config.get("prediction_variance_percentile") is not None:
            ran.append(config["prediction_variance_percentile"])   # would only fire if arms ran
        s = 0.50
        if config["method"] == "wscore" and config["wscore_covariate_model"] == "linear":
            s += 0.05
        return s

    best, history = st.run_staged_optimization(
        cohorts=[fake_cohort], env=None, score_fn=mock_score, verbose=False,
        fail_fast=False,  # bare mock cohort can't build a site harmonizer -> skip that arm
    )
    assert best["method"] == "wscore" and best["wscore_covariate_model"] == "linear"
    assert best.get("prediction_variance_percentile") is None      # no filter applied
    assert ran == []                                               # uncertainty arms never evaluated
    # the stage still yields exactly one winner row (a keep-skip)
    unc = history[(history["stage"] == "gp_uncertainty") & (history["is_winner"] == True)]  # noqa: E712
    assert len(unc) == 1
    assert "skipped" in str(unc.iloc[0]["candidate"])


# --- Step 6b: ADC/FA/qT1 sampling -> quant_sample_surface --------------------
def test_sample_surface_resolves_both_params():
    # one sample_surface drives BOTH T1w/FLAIR intensity_depth_model AND the
    # ADC/FA/qT1 quant_sample_surface (plain sample, no self-norm).
    merged, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="white"))
    assert merged["intensity_depth_model"] == "sample_white"   # T1w/FLAIR
    assert merged["quant_sample_surface"] == "white"           # ADC/FA/qT1
    # self-norm applies to T1w/FLAIR only, NOT to the quant surface
    merged2, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="swm1",
                                  t1w_flair_self_norm="swm2"))
    assert merged2["intensity_depth_model"] == "sample_swm1_swm2"
    assert merged2["quant_sample_surface"] == "swm1"           # no _swm2 for quant
    # raw override: both off
    raw, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="raw"))
    assert raw["intensity_depth_model"] == "raw" and raw["quant_sample_surface"] is None
    # baseline: white surface on
    base, _ = st.validate(dict(st.SIMPLEST_BASELINE))
    assert base["intensity_depth_model"] == "sample_white" and base["quant_sample_surface"] == "white"


def test_analysis_kwargs_forwards_quant_sample_surface():
    merged, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="white"))
    kwargs = st._analysis_kwargs(merged, 0.5)
    assert kwargs.get("quant_sample_surface") == "white"
    # analyze() actually accepts the param (backend plumbing landed)
    import inspect
    from zbrains.dataset import zbdataset
    assert "quant_sample_surface" in inspect.signature(zbdataset.analyze).parameters


def test_analysis_kwargs_skips_quant_sample_surface_without_required_companions():
    merged, _ = st.validate(dict(st.SIMPLEST_BASELINE, sample_surface="white"))
    mics_like = ["T1w", "T1w*blur", "FA", "ADC", "thickness", "FLAIR", "FLAIR*blur", "qT1", "qT1*blur"]
    kwargs = st._analysis_kwargs(merged, 0.5, features=mics_like)
    assert "quant_sample_surface" not in kwargs

    companioned = [
        "T1w", "T1w*blur", "FLAIR", "FLAIR*blur",
        "FA", "FA*blur", "ADC", "ADC*blur", "qT1", "qT1*blur",
    ]
    kwargs2 = st._analysis_kwargs(merged, 0.5, features=companioned)
    assert kwargs2.get("quant_sample_surface") == "white"


# --- Full staged driver runs end-to-end (mock processing) -------------------
def test_full_default_stages_mock_run(monkeypatch):
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(
        st, "_analyze_candidate",
        lambda cohort, config, env, thr, verbose=False, **k: "/fake",
    )
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")

    def mock_score(config, output_dirs):
        s = 0.50
        if config["subject_norm"] == "whitestripe":
            s += 0.02
        if str(config.get("intensity_depth_model", "raw")).startswith("sample_white"):
            s += 0.03
        if config["t1w_flair_self_norm"] == "swm2":
            s += 0.01
        if config["blur_depth_model"] == "boundary_gradient":
            s += 0.01
        if config["method"] == "wscore" and config["wscore_covariate_model"] == "gaussian_process":
            s += 0.02
        return s

    best, history = st.run_staged_optimization(
        cohorts=[fake_cohort], env=None, score_fn=mock_score, verbose=False,
        fail_fast=False,  # bare mock cohort can't build a site harmonizer -> skip that arm
    )
    winners = history[history["is_winner"] == True]  # noqa: E712  (baseline + 9 stages)
    assert len(winners) == 1 + len(st.DEFAULT_STAGES)
    # locked choices reflect the mock rewards (dataset_norm stays none: ravel/nyul
    # would drop the whitestripe reward, so the greedy keeps subject_norm=whitestripe)
    assert best["subject_norm"] == "whitestripe" and best["dataset_norm"] == "none"
    assert best["intensity_depth_model"] == "sample_white_swm2"
    assert best["blur_depth_model"] == "boundary_gradient"
    assert best["method"] == "wscore" and best["wscore_covariate_model"] == "gaussian_process"
    assert best["cortical_smoothing"] == 0 and best["hippocampal_smoothing"] == 0


def test_default_score_fn_objective_modes(monkeypatch):
    import results.greedy_auroc as ga
    # within leg: pooled_macro_auroc -> (_, "WITHIN"); control leg: a real df concat
    monkeypatch.setattr(ga, "pooled_macro_auroc", lambda *a, **k: (None, "WITHIN"))
    monkeypatch.setattr(ga, "per_patient_vs_control_auc",
                        lambda root, name, **k: pd.DataFrame(
                            {"auc": [0.6], "dataset": [name], "lesion_type": ["FCD"], "feature": ["T1w"]}))
    # balanced_macro_auroc: 0.80 for the within (string) detail, 0.60 for the control df
    monkeypatch.setattr(ga, "balanced_macro_auroc",
                        lambda detail, **k: 0.80 if isinstance(detail, str) else 0.60)
    config = {"method": "wscore"}
    dirs = {"MICs": "/x", "NOEL": "/y"}
    cs = {"MICs": [("s", "ses-01")], "NOEL": [("s", "ses-01")]}

    assert st.default_score_fn(config, dirs, control_subjects_by_name=cs, objective_mode="within") == 0.80
    assert st.default_score_fn(config, dirs, control_subjects_by_name=cs, objective_mode="control") == 0.60
    # 'both' is now an ALIAS of 'within' (within-subject, TLE+FCD together)
    assert st.default_score_fn(config, dirs, control_subjects_by_name=cs, objective_mode="both") == 0.80
    # 'both' == within regardless of whether held-out controls exist
    assert st.default_score_fn(config, dirs, control_subjects_by_name=None, objective_mode="both") == 0.80
    # 'control' with no controls -> warns + falls back to within
    assert st.default_score_fn(config, dirs, control_subjects_by_name=None, objective_mode="control") == 0.80
    # invalid mode rejected
    try:
        st.default_score_fn(config, dirs, objective_mode="bogus")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_reference_cv_score_fn_reports_per_fold_mean_and_std(monkeypatch):
    # reference_cv_score_fn evaluates the SINGLE-SHOT default_score_fn once per
    # reference fold (each fold has its OWN analysis root) and locks on the MEAN,
    # stashing the per-fold vector + mean/SD for the history CSV -- for EVERY mode.
    import results.greedy_auroc as ga
    monkeypatch.setattr(ga, "pooled_macro_auroc", lambda *a, **k: (None, "WITHIN"))
    # per_patient_vs_control_auc value depends on that fold's control subset:
    # 1 held-out ctrl -> 0.6, else 0.8.
    def fake_ppc(root, name, *, control_subjects=None, **k):
        val = 0.6 if len(control_subjects) == 1 else 0.8
        return pd.DataFrame({"auc": [val], "dataset": [name],
                             "lesion_type": ["FCD"], "feature": ["T1w"]})
    monkeypatch.setattr(ga, "per_patient_vs_control_auc", fake_ppc)
    monkeypatch.setattr(ga, "balanced_macro_auroc",
                        lambda detail, **k: 0.9 if isinstance(detail, str)
                        else float(detail["auc"].mean()))   # within leg passes a str stub
    config = {"method": "wscore"}
    fold_dirs = {"MICs": {0: "/f0", 1: "/f1"}}               # a distinct root per fold
    cs = {"MICs": {0: [("s1", "ses-01")],                    # fold 0: 1 ctrl -> 0.6
                   1: [("s2", "ses-01"), ("s3", "ses-01")]}}  # fold 1: 2 ctrls -> 0.8

    val = st.reference_cv_score_fn(config, fold_dirs, control_subjects_by_name=cs,
                                   objective_mode="control")
    assert abs(val - 0.7) < 1e-9                              # locks on the MEAN of the folds
    stats = st._LAST_CONTROL_FOLD_STATS
    assert stats is not None and stats["n_folds"] == 2
    assert sorted(stats["aurocs"]) == [0.6, 0.8]
    assert abs(stats["std"] - 0.1414213562) < 1e-6           # sample SD across folds (ddof=1)
    # per-feature decomposition (averaged across folds): the single feature "T1w"
    # here == the overall mean since it's the only feature present.
    assert abs(stats["per_feature"]["T1w"] - 0.7) < 1e-9

    # within mode ALSO folds now (5-values-for-all): each fold root -> 0.9 -> mean 0.9
    val_w = st.reference_cv_score_fn(config, fold_dirs, control_subjects_by_name=cs,
                                     objective_mode="within")
    assert abs(val_w - 0.9) < 1e-9 and st._LAST_CONTROL_FOLD_STATS["n_folds"] == 2

    # single-shot default_score_fn still works with a FLAT control list (one value)
    v_single = st.default_score_fn(config, {"MICs": "/f1"},
                                   control_subjects_by_name={"MICs": [("s2", "ses-01"), ("s3", "ses-01")]},
                                   objective_mode="control")
    assert abs(v_single - 0.8) < 1e-9


def test_default_score_fn_within_disease_filters_and_disease_control(monkeypatch):
    import results.greedy_auroc as ga
    # within detail has an FCD row (0.9) and a TLE row (0.6)
    detail = pd.DataFrame({"dataset": ["MICs", "MICs"], "lesion_type": ["FCD", "TLE"],
                           "feature": ["T1w", "T1w"], "auc": [0.9, 0.6]})
    monkeypatch.setattr(ga, "pooled_macro_auroc", lambda *a, **k: (None, detail))
    monkeypatch.setattr(ga, "per_patient_vs_other_disease_auc",
                        lambda root, name, **k: pd.DataFrame(
                            {"dataset": [name], "lesion_type": ["FCD"], "feature": ["T1w"], "auc": [0.7]}))
    # balanced_macro faked as the mean AUC of whatever detail survives filtering
    monkeypatch.setattr(ga, "balanced_macro_auroc",
                        lambda d, **k: float(d["auc"].mean()) if (d is not None and not d.empty) else float("nan"))
    config = {"method": "wscore"}
    dirs = {"MICs": "/x"}

    assert st.default_score_fn(config, dirs, objective_mode="within") == 0.75       # both diseases
    assert st.default_score_fn(config, dirs, objective_mode="within_fcd") == 0.90   # FCD rows only
    assert st.default_score_fn(config, dirs, objective_mode="within_tle") == 0.60   # TLE rows only
    assert st.default_score_fn(config, dirs, objective_mode="disease_control") == 0.70


def test_cross_evaluate_builds_matrix_and_picks_robust_winner(monkeypatch):
    # analyze once per config (stubbed); score under every mode from a per-config table.
    # reference_cv=False exercises the legacy single-root scoring path (matrix-assembly
    # logic is identical either way; the per-fold plumbing is covered separately).
    monkeypatch.setattr(st, "analyze_config_across_cohorts",
                        lambda *a, **k: ({"MICs": "/x"}, {"MICs": [("s", "ses-01")]}))

    def fake_score(config, dirs, *, control_subjects_by_name=None, objective_mode="both"):
        return config["_scores"][objective_mode]

    monkeypatch.setattr(st, "default_score_fn", fake_score)
    # A: strong on within, weak on control/disease -> mean 0.78
    # B: strong on control/disease, weak on within -> mean 0.72  (less robust)
    A = {"_scores": {"within_tle": 0.9, "within_fcd": 0.9, "within": 0.9,
                     "control": 0.6, "disease_control": 0.6}}
    B = {"_scores": {"within_tle": 0.6, "within_fcd": 0.6, "within": 0.6,
                     "control": 0.9, "disease_control": 0.9}}
    matrix = st.cross_evaluate_configs({"A": A, "B": B}, cohorts=[object()], env=None,
                                       reference_cv=False)

    assert list(matrix.index) == ["A", "B"]
    for m in st.CROSS_EVAL_MODES:
        assert m in matrix.columns
    assert matrix.loc["A", "within_tle"] == 0.9 and matrix.loc["A", "control"] == 0.6
    assert abs(matrix.loc["A", "mean"] - 0.78) < 1e-9
    assert abs(matrix.loc["B", "mean"] - 0.72) < 1e-9
    assert matrix["mean"].idxmax() == "A"          # A is the more robust overall winner
    assert matrix.loc["A", "std"] > 0              # varies across metrics


def test_reference_cv_drives_default_objective_and_records_folds(monkeypatch):
    # cv_scope="full" -> even a within-mode run is reference-CV'd: the driver calls
    # analyze_config_across_cohorts(reference_cv=True) and scores via
    # reference_cv_score_fn, recording per-fold + per-feature columns.
    seen = {}

    def fake_analyze(config, cohorts, env, thr, *, needs_controls, built_bases=None,
                     reprocess_controls=False, verbose=False, reference_cv=False):
        seen["reference_cv"] = reference_cv
        seen["needs_controls"] = needs_controls
        fold_dirs = {"MICs": {0: "/f0", 1: "/f1", 2: "/f2"}}
        pairs = {"MICs": {0: [("c0", "ses-01")], 1: [("c1", "ses-01")], 2: [("c2", "ses-01")]}}
        return fold_dirs, pairs

    def fake_single(config, odk, *, control_subjects_by_name=None, objective_mode="both"):
        overall = {"/f0": 0.60, "/f1": 0.70, "/f2": 0.80}[odk["MICs"]]
        # stash a per-feature split like the real default_score_fn does
        st._LAST_FEATURE_AUROC = {"T1w": overall, "FLAIR": overall - 0.10}
        return overall

    monkeypatch.setattr(st, "analyze_config_across_cohorts", fake_analyze)
    monkeypatch.setattr(st, "default_score_fn", fake_single)   # reference_cv_score_fn calls this per fold
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")

    with tempfile.TemporaryDirectory() as tmp:
        hp = os.path.join(tmp, "hist.csv")
        best, hist = st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[], objective_mode="within",
            score_fn=None, verbose=False, resume=False, history_path=hp, cv_scope="full")

    assert seen["reference_cv"] is True            # cv_scope=full -> within reference-CV'd
    assert seen["needs_controls"] is False         # within: no held-out-control negatives
    base = hist[hist["candidate"] == "baseline"].iloc[0]
    assert abs(float(base["auroc"]) - 0.70) < 1e-9         # locked on the per-fold MEAN
    assert int(base["control_n_folds"]) == 3               # within-mode now folds too
    assert abs(float(base["control_fold_mean"]) - 0.70) < 1e-9
    assert float(base["control_fold_std"]) > 0             # spread across folds recorded
    # per-FEATURE columns (averaged across folds) appear alongside the overall mean
    assert abs(float(base["auroc_T1w"]) - 0.70) < 1e-9     # mean of 0.6/0.7/0.8
    assert abs(float(base["auroc_FLAIR"]) - 0.60) < 1e-9   # mean of 0.5/0.6/0.7
    assert "auroc_T1w_std" in hist.columns and float(base["auroc_T1w_std"]) > 0


def test_cv_scope_gates_which_modes_use_reference_cv():
    # The scope helper: default reference-CVs ONLY 'control' (its negatives are
    # controls, which also train the normative model -> leakage to hold out);
    # "full" does every mode; "off" does none.
    assert st._mode_uses_reference_cv("control", "control_disease") is True
    # disease_control's negatives are other-disease patients (never in the control
    # fit) -> no leakage -> single-shot under the default scope.
    assert st._mode_uses_reference_cv("disease_control", "control_disease") is False
    assert st._mode_uses_reference_cv("both", "control_disease") is False
    assert st._mode_uses_reference_cv("within", "control_disease") is False
    assert st._mode_uses_reference_cv("within_tle", "control_disease") is False
    for m in ("within", "within_tle", "control", "disease_control"):
        assert st._mode_uses_reference_cv(m, "full") is True
        assert st._mode_uses_reference_cv(m, "off") is False
    try:
        st._mode_uses_reference_cv("within", "bogus")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_default_scope_within_is_single_shot_control_is_reference_cv(monkeypatch):
    # cv_scope="control_disease" (DEFAULT): within -> single-shot (one value,
    # reference_cv=False, still records per-feature); control -> reference-CV.
    calls = {}

    def fake_analyze(config, cohorts, env, thr, *, needs_controls, built_bases=None,
                     reprocess_controls=False, verbose=False, reference_cv=False):
        calls["reference_cv"] = reference_cv
        if reference_cv:
            return ({"MICs": {0: "/f0", 1: "/f1"}},
                    {"MICs": {0: [("c0", "ses-01")], 1: [("c1", "ses-01")]}})
        return {"MICs": "/root"}, {}                     # single-root shape

    def fake_single(config, odk, *, control_subjects_by_name=None, objective_mode="both"):
        st._LAST_FEATURE_AUROC = {"T1w": 0.66}
        return 0.66

    monkeypatch.setattr(st, "analyze_config_across_cohorts", fake_analyze)
    monkeypatch.setattr(st, "default_score_fn", fake_single)
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")

    # within -> single-shot (default scope): analyze called with reference_cv=False,
    # n_folds=1, but per-feature still recorded.
    with tempfile.TemporaryDirectory() as tmp:
        _b, hist = st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[], objective_mode="within",
            score_fn=None, verbose=False, resume=False,
            history_path=os.path.join(tmp, "w.csv"))
    assert calls["reference_cv"] is False
    base = hist[hist["candidate"] == "baseline"].iloc[0]
    assert abs(float(base["auroc"]) - 0.66) < 1e-9
    assert int(base["control_n_folds"]) == 1 and float(base["control_fold_std"]) == 0.0
    assert abs(float(base["auroc_T1w"]) - 0.66) < 1e-9      # per-feature recorded single-shot too

    # control -> reference-CV under the same default scope
    with tempfile.TemporaryDirectory() as tmp:
        st.run_staged_optimization(
            cohorts=[fake_cohort], env=None, stages=[], objective_mode="control",
            score_fn=None, verbose=False, resume=False,
            history_path=os.path.join(tmp, "c.csv"))
    assert calls["reference_cv"] is True


def test_feature_breakdown_decomposes_overall_by_feature():
    # _feature_breakdown decomposes the overall balanced AUROC into one value per
    # feature, using the SAME cell-averaging (balanced_macro_auroc restricted).
    from results.greedy_auroc import balanced_macro_auroc
    detail = pd.DataFrame({
        "dataset": ["MICs", "MICs", "NOEL", "NOEL"],
        "lesion_type": ["FCD", "FCD", "TLE", "TLE"],
        "feature": ["T1w", "FLAIR", "T1w", "FLAIR"],
        "auc": [0.9, 0.5, 0.7, 0.3],
    })
    gc = ("dataset", "lesion_type", "feature")
    fb = st._feature_breakdown(detail, None, gc)
    assert set(fb) == {"T1w", "FLAIR"}
    # each == balanced_macro_auroc restricted to that feature ...
    assert abs(fb["T1w"] - balanced_macro_auroc(detail, features={"T1w"}, group_cols=gc)) < 1e-12
    assert abs(fb["FLAIR"] - balanced_macro_auroc(detail, features={"FLAIR"}, group_cols=gc)) < 1e-12
    # ... T1w cells (0.9, 0.7) -> 0.8 ; FLAIR (0.5, 0.3) -> 0.4
    assert abs(fb["T1w"] - 0.8) < 1e-9 and abs(fb["FLAIR"] - 0.4) < 1e-9
    # a whitelist restricts which features appear
    assert set(st._feature_breakdown(detail, {"T1w"}, gc)) == {"T1w"}
    # degenerate details -> empty (no crash)
    assert st._feature_breakdown(pd.DataFrame(), None, gc) == {}
    assert st._feature_breakdown("WITHIN", None, gc) == {}


def test_disease_control_mode_does_not_need_controls():
    # ONLY 'control' triggers the expensive held-out-control K-fold scoring. Its
    # negatives ARE controls (also the normative training data), so they must be
    # held out. within_*/both (within-subject) and disease_control (negatives =
    # other-disease patients, never in the control fit) do NOT need held-out controls.
    assert st._CONTROL_MODES == ("control",)
    assert "both" not in st._CONTROL_MODES
    assert "disease_control" not in st._CONTROL_MODES
    assert "within_tle" not in st._CONTROL_MODES and "within_fcd" not in st._CONTROL_MODES


def test_fail_fast_stops_on_candidate_failure(monkeypatch):
    # A candidate whose evaluation raises (e.g. a RAVEL SIGSEGV) must STOP the run
    # under fail_fast=True (record the error, then raise), and must be SKIPPED
    # (lock the incumbent) under fail_fast=False.
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(st, "_analyze_candidate",
                        lambda cohort, config, env, thr, verbose=False, **k: "/fake")
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")

    def boom_score(config, output_dirs):
        if config["subject_norm"] == "whitestripe":
            raise RuntimeError("simulated RAVEL SIGSEGV")
        return 0.5

    stage = st.Stage("subject_norm", [("whitestripe", {"subject_norm": "whitestripe"})])

    with tempfile.TemporaryDirectory() as td:
        hp = os.path.join(td, "h.csv")
        try:
            st.run_staged_optimization(
                cohorts=[fake_cohort], env=None, stages=[stage], score_fn=boom_score,
                verbose=False, resume=False, history_path=hp, fail_fast=True)
            assert False, "expected RuntimeError under fail_fast=True"
        except RuntimeError as exc:
            assert "FAILED" in str(exc) and "whitestripe" in str(exc)
        df = pd.read_csv(hp)                       # the error row was checkpointed
        row = df[df["candidate"] == "whitestripe"]
        assert len(row) == 1 and row.iloc[0]["error"] == row.iloc[0]["error"]  # non-NaN

    # fail_fast=False -> the failed arm is skipped and the incumbent (none) is locked
    best, _hist = st.run_staged_optimization(
        cohorts=[fake_cohort], env=None, stages=[stage], score_fn=boom_score,
        verbose=False, resume=False, fail_fast=False)
    assert best["subject_norm"] == "none"


def test_errored_candidate_is_retried_on_resume(monkeypatch):
    # An error row must NOT count as "done" -- a resumed run retries it (so a fixed
    # crash re-runs instead of being permanently skipped).
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: None)
    monkeypatch.setattr(st, "_analyze_candidate",
                        lambda cohort, config, env, thr, verbose=False, **k: "/fake")
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")
    stage = st.Stage("subject_norm", [("whitestripe", {"subject_norm": "whitestripe"})])
    attempts = {"n": 0}

    def flaky_score(config, output_dirs):
        if config["subject_norm"] == "whitestripe":
            attempts["n"] += 1
            if attempts["n"] == 1:
                raise RuntimeError("transient crash (1st attempt)")
            return 0.9                              # succeeds on retry
        return 0.5

    with tempfile.TemporaryDirectory() as td:
        hp = os.path.join(td, "h.csv")
        try:
            st.run_staged_optimization(cohorts=[fake_cohort], env=None, stages=[stage],
                score_fn=flaky_score, verbose=False, resume=False, history_path=hp, fail_fast=True)
            assert False
        except RuntimeError:
            pass
        # resume: the errored whitestripe candidate is retried (not skipped) and now wins
        best, _hist = st.run_staged_optimization(cohorts=[fake_cohort], env=None, stages=[stage],
            score_fn=flaky_score, verbose=False, resume=True, history_path=hp, fail_fast=True)
        assert attempts["n"] == 2                   # it was retried, not skipped
        assert best["subject_norm"] == "whitestripe"


def test_default_score_fn_partial_preserves_control_subjects_signature():
    # run_staged_optimization binds objective_mode via functools.partial; the
    # control_subjects_by_name param must remain visible so evaluate() still passes
    # held-out controls to the control/both legs.
    import functools
    import inspect
    bound = functools.partial(st.default_score_fn, objective_mode="both")
    assert "control_subjects_by_name" in inspect.signature(bound).parameters


def test_run_staged_gates_control_scoring_by_objective_mode(monkeypatch):
    # ONLY 'control' triggers held-out-control scoring; 'within' and 'both' (now an
    # alias of within) must NOT (the gate is driven by objective_mode in the driver,
    # independent of score_fn).
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset",
                        lambda cohort, excluded: types.SimpleNamespace(name="controls"))
    monkeypatch.setattr(st, "_analyze_candidate",
                        lambda cohort, config, env, thr, verbose=False, **k: "/fake")
    scored = []
    monkeypatch.setattr(st, "_score_heldout_controls",
                        lambda cohort, *a, **k: scored.append(cohort.score_name))
    monkeypatch.setattr(st, "_control_pairs", lambda ds: [("s", "ses-01")])
    fake_cohort = types.SimpleNamespace(name="MICs", score_name="MICs")

    st.run_staged_optimization(cohorts=[fake_cohort], env=None, stages=[],
                               score_fn=lambda c, d: 0.5, objective_mode="within",
                               verbose=False, resume=False)
    assert scored == []

    scored.clear()
    st.run_staged_optimization(cohorts=[fake_cohort], env=None, stages=[],
                               score_fn=lambda c, d: 0.5, objective_mode="both",
                               verbose=False, resume=False)
    assert scored == []                             # 'both' is now within -> no controls

    scored.clear()
    st.run_staged_optimization(cohorts=[fake_cohort], env=None, stages=[],
                               score_fn=lambda c, d: 0.5, objective_mode="control",
                               verbose=False, resume=False)
    assert scored == ["MICs"]                       # only 'control' scores held-out controls


def _refcv_mocks(monkeypatch, shared_patient, control_ds, folds, subset_factory,
                 dataset_norm="none"):
    monkeypatch.setattr(st, "staged_output_directory", lambda *a, **k: "/out/")
    monkeypatch.setattr(st, "_processing_signature", lambda *a, **k: "sig")
    monkeypatch.setattr(st, "_base_signature", lambda *a, **k: "bsig")
    monkeypatch.setattr(st, "_pipeline_settings", lambda *a, **k: {})
    monkeypatch.setattr(st, "_analysis_kwargs", lambda *a, **k: {})
    monkeypatch.setattr(st, "_control_feature_availability", lambda *a, **k: {})
    monkeypatch.setattr(st, "_control_pairs",
                        lambda ds: [("p1", "ses-01")] if ds is shared_patient
                        else [("c0", "ses-01"), ("c1", "ses-01"), ("c2", "ses-01"), ("c3", "ses-01")])
    monkeypatch.setattr(st, "_control_folds", lambda ds, **k: folds)
    monkeypatch.setattr(st, "_subset_control_dataset", subset_factory)
    monkeypatch.setattr(st.zb, "processed_base_directory_for", lambda *a, **k: "/base")
    monkeypatch.setattr(st.zb, "symlink_processed_outputs", lambda *a, **k: None)


def test_reference_cv_clones_patient_and_returns_per_fold_roots(monkeypatch):
    # _analyze_reference_cv must NOT validate/analyze the SHARED patient dataset
    # (mutates .features/.valid_subjects) -- it clones it per fold -- and it returns
    # per-fold roots + held-out pairs. Train-fold reference keeps a control name
    # (FITS), the held-out fold gets a non-control name (APPLIES frozen).
    shared_calls = []

    class _SharedPatient:
        name = "patients"
        demographics = types.SimpleNamespace(
            data=pd.DataFrame({"participant_id": ["p1"], "session_id": ["ses-01"]}))
        def validate(self, **k): shared_calls.append("validate")
        def analyze(self, **k): shared_calls.append("analyze")

    shared_patient = _SharedPatient()
    clone_ops = []
    names_seen = []

    class _CloneDS:
        def __init__(self, tag): self.tag = tag
        def validate(self, **k): clone_ops.append((self.tag, "validate"))
        def analyze(self, **k): clone_ops.append((self.tag, "analyze"))

    def subset_factory(orig, pairs, name=None):
        names_seen.append(name if name is not None else getattr(orig, "name", "controls"))
        tag = "patient" if orig is shared_patient else ("heldout" if name == "heldout" else "reference")
        return _CloneDS(tag)

    control_ds = types.SimpleNamespace(name="controls")
    folds = {0: [("c0", "ses-01")], 1: [("c1", "ses-01")]}
    _refcv_mocks(monkeypatch, shared_patient, control_ds, folds, subset_factory)

    cohort = types.SimpleNamespace(output_dir_prefix="/pfx", score_name="MICs",
                                   patient_dataset=shared_patient, base_pipeline_settings={})
    cfg = {"normalization": "none", "dataset_norm": "none", "site_harmonization": "none"}
    roots, subj = st._analyze_reference_cv([(cohort, control_ds, "")], cfg, None, 0.5,
                                           needs_controls=True, verbose=False)

    assert shared_calls == []                                  # shared patient never touched
    assert set(roots["MICs"]) == {0, 1}                        # per-fold roots
    assert {tuple(v) for v in subj["MICs"].values()} == {(("c0", "ses-01"),), (("c1", "ses-01"),)}
    # each fold analyzed a patient clone + held-out clone (needs_controls=True)
    assert ("patient", "analyze") in clone_ops and ("heldout", "analyze") in clone_ops
    # held-out fold was created with the non-control name "heldout" (-> applies frozen)
    assert "heldout" in names_seen


def test_reference_cv_builds_harmonizer_per_fold(monkeypatch):
    # site_harmonization != none -> a NEW harmonizer is fit PER FOLD on that fold's
    # reference (no 2nd-order harmonizer leak); each fold's scoring uses its own.
    class _SharedPatient:
        name = "patients"
        demographics = types.SimpleNamespace(
            data=pd.DataFrame({"participant_id": ["p1"], "session_id": ["ses-01"]}))
    shared_patient = _SharedPatient()
    seen_harmonizers = []

    class _CloneDS:
        def __init__(self, tag): self.tag = tag
        def validate(self, **k): pass
        def analyze(self, **k):
            if self.tag == "patient":
                seen_harmonizers.append(k.get("site_harmonizer"))

    def subset_factory(orig, pairs, name=None):
        return _CloneDS("patient" if orig is shared_patient
                        else ("heldout" if name == "heldout" else "reference"))

    control_ds = types.SimpleNamespace(name="controls")
    folds = {0: [("c0", "ses-01")], 1: [("c1", "ses-01")]}
    _refcv_mocks(monkeypatch, shared_patient, control_ds, folds, subset_factory)
    built, exported = [], []

    def fake_build(config, cohorts):
        h = types.SimpleNamespace(add_site_features=lambda *a, **k: None)
        built.append(h)
        return h
    monkeypatch.setattr(st, "_build_site_harmonizer", fake_build)
    monkeypatch.setattr(st, "_export_fold_control_features", lambda *a, **k: exported.append(1))

    cohort = types.SimpleNamespace(output_dir_prefix="/pfx", score_name="MICs",
                                   patient_dataset=shared_patient, base_pipeline_settings={})
    cfg = {"normalization": "none", "dataset_norm": "none", "site_harmonization": "combat"}
    st._analyze_reference_cv([(cohort, control_ds, "")], cfg, None, 0.5,
                             needs_controls=True, verbose=False)

    assert len(built) == 2                       # ONE harmonizer per fold
    assert len(exported) == 2                    # reference_k exported once per fold (1 cohort x 2 folds)
    assert len(seen_harmonizers) == 2 and seen_harmonizers[0] is not seen_harmonizers[1]  # per-fold


def test_ensure_fold_base_fits_reference_applies_heldout_and_patients(monkeypatch):
    # Per-fold dataset-norm base: reference_k (control-named) FITS; held-out + patients
    # are processed to APPLY the frozen fit. All three get processed into one fold base.
    processed = []

    class _DS:
        def __init__(self, tag): self.tag = tag
        def process(self, **k): processed.append(self.tag)

    reference_k = _DS("reference"); heldout_k = _DS("heldout")
    patient = _DS("patient")
    tmp = tempfile.mkdtemp()
    cohort = types.SimpleNamespace(name="MICs", patient_dataset=patient,
                                   output_dir_prefix=tmp, base_pipeline_settings={})
    base = os.path.join(tmp, "zbrains_WBRAVEL_reffold0"); os.makedirs(base, exist_ok=True)
    monkeypatch.setattr(st, "_base_signature", lambda *a, **k: "bsig")
    monkeypatch.setattr(st, "_pipeline_settings", lambda *a, **k: {})
    monkeypatch.setattr(st, "_control_pairs",
                        lambda ds: [(f"{ds.tag}0", "ses-01")])   # 1 pair per ds, tag-distinct
    monkeypatch.setattr(st.zb, "processed_base_directory_for", lambda *a, **k: base)
    monkeypatch.setattr(st.zb, "processed_maps_complete", lambda *a, **k: False)  # force processing
    monkeypatch.setattr(st.zb, "base_is_marked_complete", lambda *a, **k: False)

    out = st._ensure_fold_base(cohort, {"normalization": "whitestripeRavel", "dataset_norm": "ravel"},
                               None, reference_k, heldout_k, "exclX", 0, verbose=False)
    assert out == base
    assert processed == ["reference", "heldout", "patient"]   # fit first, then apply to held-out + patients


def test_score_heldout_controls_parallel_folds(monkeypatch):
    # The K control-CV folds are independent analyze() passes scored CONCURRENTLY
    # in threads: every non-skipped fold runs (output-identical), skipped folds warn
    # + are excluded from the returned map, and env.num_threads>1 uses >1 worker.
    import threading, time
    import warnings as _w

    pairs = [(f"hc{i}", "ses-01") for i in range(6)]
    # 3 partitioning folds + 1 fold whose complement has <2 training controls (skip)
    folds = {0: pairs[0:2], 1: pairs[2:4], 2: pairs[4:6], 3: pairs[0:5]}
    recorded = []   # (heldout_pairs, thread_name) per analyze() call

    class _StubDS:
        def __init__(self, subset): self.pairs = tuple(map(tuple, subset))
        def validate(self, **k): pass
        def analyze(self, **k):
            recorded.append((self.pairs, threading.current_thread().name))
            time.sleep(0.03)   # widen the window so concurrent folds overlap

    monkeypatch.setattr(st, "staged_output_directory", lambda *a, **k: "/out")
    monkeypatch.setattr(st, "_processing_signature", lambda *a, **k: "sig")
    monkeypatch.setattr(st, "_pipeline_settings", lambda *a, **k: {})
    monkeypatch.setattr(st, "_analysis_kwargs", lambda *a, **k: {})
    monkeypatch.setattr(st, "_control_feature_availability", lambda *a, **k: {})
    monkeypatch.setattr(st, "_control_pairs", lambda ds: pairs)
    monkeypatch.setattr(st, "_control_folds", lambda ds, **k: folds)
    monkeypatch.setattr(st, "_subset_control_dataset", lambda ds, subset: _StubDS(subset))
    monkeypatch.setattr(st.zb, "processed_base_directory_for", lambda *a, **k: "/base")

    cohort = types.SimpleNamespace(output_dir_prefix="/pfx", score_name="MICs")
    config = {"normalization": "none"}

    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        scored = st._score_heldout_controls(
            cohort, config, types.SimpleNamespace(num_threads=4), 0.5, object(), "")

    # the 3 valid folds were scored (fold 3 skipped: only 1 training control) ...
    assert set(scored) == {0, 1, 2}
    assert {tuple(v) for v in scored.values()} == {
        tuple(map(tuple, pairs[0:2])), tuple(map(tuple, pairs[2:4])), tuple(map(tuple, pairs[4:6]))}
    assert len(recorded) == 3                                   # analyze ran once per scored fold
    assert len({t for _, t in recorded}) > 1                   # ... concurrently (>1 thread)
    assert any("fold skipped" in str(w.message) for w in caught)   # skipped fold warned

    # env.num_threads=1 -> serial fallback: all folds scored, a single thread
    recorded.clear()
    scored2 = st._score_heldout_controls(
        cohort, config, types.SimpleNamespace(num_threads=1), 0.5, object(), "")
    assert set(scored2) == {0, 1, 2} and len(recorded) == 3
    assert len({t for _, t in recorded}) == 1


if __name__ == "__main__":
    class _MP:
        """Minimal monkeypatch stand-in that RESTORES patched attributes after each
        test (so a test patching a module-level function can't leak into the next)."""
        def __init__(self):
            self._undo = []
        def setattr(self, obj, name, val):
            self._undo.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        def undo(self):
            for obj, name, orig in reversed(self._undo):
                setattr(obj, name, orig)
            self._undo = []
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        mp = _MP()
        try:
            fn(mp) if "monkeypatch" in fn.__code__.co_varnames else fn()
        finally:
            mp.undo()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} GREEDY-BENCHMARK TESTS PASSED")
