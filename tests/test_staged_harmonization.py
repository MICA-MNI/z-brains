"""Driver-level wiring of the site_harmonization stage (two-pass orchestration +
keying), with heavy analysis mocked out. The ComBat math + SiteHarmonizer are
covered in test_harmonization.py; the real analyze_dataset injection is validated
by the heavy 2-cohort smoke.

Run: python tests/test_staged_harmonization.py  (or under pytest).
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd

import zbrains_staged as st
import zbrains_benchmark as zb


def _two_cohorts():
    def _demo(n, tag):
        data = pd.DataFrame({"participant_id": [f"sub-{tag}{i}" for i in range(n)],
                             "session_id": ["ses-01"] * n,
                             "AGE": np.linspace(25, 55, n), "SEX": ([0, 1] * n)[:n]})
        return types.SimpleNamespace(data=data, normative_columns=["AGE", "SEX"])
    def _cohort(name):
        ctrl = types.SimpleNamespace(demographics=_demo(12, name[:2]))
        return types.SimpleNamespace(name=name, score_name=name, control_dataset=ctrl,
                                     patient_dataset=None)
    return [_cohort("MICs"), _cohort("NOEL")]


def _install_mocks(monkeypatch, calls):
    monkeypatch.setattr(st, "_ensure_base_processed", lambda *a, **k: None)
    monkeypatch.setattr(st, "_resolve_excluded_pairs", lambda *a, **k: frozenset())
    monkeypatch.setattr(st, "_restricted_control_dataset", lambda cohort, excluded: cohort.control_dataset)
    monkeypatch.setattr(st, "_score_heldout_controls", lambda *a, **k: None)  # heavy; stub out

    _HK = "T1w|L_32k_white_regional"

    def fake_analyze(cohort, config, env, thr, *, control_dataset=None,
                     exclusion_signature="", verbose=False,
                     site_harmonizer=None, export_control_features=None):
        rng = np.random.default_rng(abs(hash(cohort.score_name)) % 2**32)
        shift = 0.0 if cohort.score_name == "MICs" else 4.0     # a site effect
        if export_control_features is not None:
            # EXPORT pass: stash synthetic control features for the shared map key
            mat = rng.normal(shift, 1.0, (12, 20))
            demo = cohort.control_dataset.demographics.data
            subs = list(zip(demo["participant_id"], demo["session_id"]))
            export_control_features[_HK] = (mat, subs, demo)
            calls.append((cohort.score_name, "export"))
            return f"/fake/{cohort.score_name}"
        calls.append((cohort.score_name, "score",
                      site_harmonizer is not None and site_harmonizer.covers(_HK)))
        return f"/fake/{cohort.score_name}"

    monkeypatch.setattr(st, "_analyze_candidate", fake_analyze)


def test_none_arm_is_single_pass_no_harmonizer(monkeypatch):
    calls = []
    _install_mocks(monkeypatch, calls)
    cohorts = _two_cohorts()
    stage = st.Stage("site_harmonization", [])   # only the implicit keep = none
    best, _ = st.run_staged_optimization(
        cohorts=cohorts, env=None, stages=[stage],
        score_fn=lambda c, d: 0.5, verbose=False, resume=False,
    )
    # exactly one (score) call per cohort; no export pass; harmonizer covers=... never True
    assert all(c[1] == "score" for c in calls)
    assert [c[0] for c in calls] == ["MICs", "NOEL"]
    assert all(c[2] is False for c in calls)     # site_harmonizer is None -> covers path not hit


def test_combat_arm_runs_two_pass_and_passes_pooled_harmonizer(monkeypatch):
    calls = []
    _install_mocks(monkeypatch, calls)
    cohorts = _two_cohorts()
    # score the combat arm directly as the baseline (stages=[] scores only start_config)
    start = dict(st.SIMPLEST_BASELINE, site_harmonization="combat")
    best, _ = st.run_staged_optimization(
        cohorts=cohorts, env=None, stages=[], start_config=start,
        score_fn=lambda c, d: 0.5, verbose=False, resume=False,
    )
    exports = [c for c in calls if c[1] == "export"]
    scores = [c for c in calls if c[1] == "score"]
    # PASS 1: one export per cohort; PASS 2: one score per cohort
    assert sorted(c[0] for c in exports) == ["MICs", "NOEL"]
    assert sorted(c[0] for c in scores) == ["MICs", "NOEL"]
    # in the scoring pass the harmonizer is present AND covers the shared 2-site map
    assert all(c[2] is True for c in scores)


def test_heldout_control_scoring_and_control_subjects_reach_objective(monkeypatch):
    calls = []
    _install_mocks(monkeypatch, calls)
    scored = set()

    # _score_heldout_controls now RETURNS the {fold_id: [held-out pairs]} map, and
    # that fold map (not a flat _control_pairs list) is what reaches the objective.
    def fake_score_heldout(cohort, *a, **k):
        scored.add(cohort.score_name)
        demo = cohort.control_dataset.demographics.data
        pairs = list(zip(demo["participant_id"], demo["session_id"]))   # 12 controls
        return {0: pairs[:6], 1: pairs[6:]}                             # 2 folds, 6 each
    monkeypatch.setattr(st, "_score_heldout_controls", fake_score_heldout)
    cohorts = _two_cohorts()
    captured = {}

    # a control-aware score_fn (like the real default_score_fn) receives the fold map
    def control_aware_score(config, output_dirs, *, control_subjects_by_name=None):
        captured["controls"] = control_subjects_by_name
        return 0.5

    st.run_staged_optimization(
        cohorts=cohorts, env=None, stages=[], start_config=dict(st.SIMPLEST_BASELINE),
        score_fn=control_aware_score, verbose=False, resume=False,
    )
    # held-out controls were scored for BOTH cohorts...
    assert scored == {"MICs", "NOEL"}
    # ...and each cohort's fold map (totaling its 12 controls) reached the objective
    assert set(captured["controls"]) == {"MICs", "NOEL"}
    for name in ("MICs", "NOEL"):
        folds = captured["controls"][name]
        assert set(folds) == {0, 1}
        assert sum(len(v) for v in folds.values()) == 12


def test_control_folds_are_subject_level_and_partition():
    # 6 subjects, one longitudinal (2 sessions) -> 7 (pid,ses) pairs.
    pids = [f"sub-{i}" for i in range(6)] + ["sub-0"]
    sess = ["ses-01"] * 6 + ["ses-02"]
    data = pd.DataFrame({"participant_id": pids, "session_id": sess,
                         "AGE": np.linspace(20, 60, len(pids)),
                         "SEX": [0, 1, 0, 1, 0, 1, 0]})
    ds = types.SimpleNamespace(demographics=types.SimpleNamespace(
        data=data, normative_columns=["AGE", "SEX"]))
    folds = st._control_folds(ds, k=3, seed=0)
    all_pairs = st._control_pairs(ds)
    flat = [p for v in folds.values() for p in v]
    assert sorted(flat) == sorted(all_pairs)          # every pair covered...
    assert len(flat) == len(set(flat))                # ...exactly once
    # both sessions of the longitudinal subject travel together (never straddle)
    for v in folds.values():
        s0 = [p for p in v if p[0] == "sub-0"]
        assert len(s0) in (0, 2), folds
    # complement of any fold (= the train set) never contains a held-out subject
    for v in folds.values():
        held_pids = {p[0] for p in v}
        train = [p for p in all_pairs if tuple(p) not in set(map(tuple, v))]
        assert held_pids.isdisjoint({p[0] for p in train})


def test_control_folds_balance_feature_availability():
    # 10 controls, constant age/sex (to isolate the feature dimension); half have
    # FLAIR, half don't. Folds must spread the FLAIR subjects evenly so no fold's
    # train reference is starved of FLAIR.
    n = 10
    pids = [f"sub-{i}" for i in range(n)]
    data = pd.DataFrame({"participant_id": pids, "session_id": ["ses-01"] * n,
                         "AGE": [30] * n, "SEX": [0] * n})
    ds = types.SimpleNamespace(demographics=types.SimpleNamespace(
        data=data, normative_columns=["AGE", "SEX"]))
    avail = {(p, "ses-01"): frozenset({"T1w", "FLAIR"} if i < 5 else {"T1w"})
             for i, p in enumerate(pids)}
    flair_subj = {p for (p, s), feats in avail.items() if "FLAIR" in feats}

    folds = st._control_folds(ds, k=5, seed=0, feature_availability=avail)
    all_pairs = st._control_pairs(ds)
    # 5 FLAIR across 5 folds -> exactly one FLAIR subject held out per fold
    for fid, test_pairs in folds.items():
        n_flair = sum(1 for (p, s) in test_pairs if p in flair_subj)
        assert n_flair == 1, (fid, test_pairs, n_flair)
    # every fold's train complement still contains FLAIR subjects (never starved)
    for test_pairs in folds.values():
        tset = set(map(tuple, test_pairs))
        train = [p for p in all_pairs if tuple(p) not in tset]
        assert sum(1 for (p, s) in train if p in flair_subj) == 4


def test_keying_baseline_identical_and_arms_distinct():
    base = dict(st.SIMPLEST_BASELINE)
    legacy = {k: v for k, v in base.items()
              if k not in ("site_harmonization", "scoring_site_covariate")}
    assert st._config_key(legacy, 0.5) == st._config_key(base, 0.5)
    m0, _ = st.validate(base)
    d0 = st.staged_output_directory(m0, "PFX", 0.5)
    dirs = {d0}
    for arm in st.HARMONIZATION_ARMS:
        m, w = st.validate(dict(base, site_harmonization=arm))
        assert not w
        d = st.staged_output_directory(m, "PFX", 0.5)
        assert d not in dirs and f"_harm{arm}" in d
        dirs.add(d)
    # harmonization must NOT change the processed base (analysis-level only)
    b0 = zb.processed_base_directory_for(m0["normalization"], "PFX", st._processing_signature(m0, ""))
    mc, _ = st.validate(dict(base, site_harmonization="combat"))
    bc = zb.processed_base_directory_for(mc["normalization"], "PFX", st._processing_signature(mc, ""))
    assert b0 == bc


if __name__ == "__main__":
    class _MP:
        def setattr(self, obj, name, val):
            setattr(obj, name, val)
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        fn(_MP()) if "monkeypatch" in fn.__code__.co_varnames else fn()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} STAGED-HARMONIZATION TESTS PASSED")
