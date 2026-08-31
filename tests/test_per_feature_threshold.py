"""Regression tests for rank-based per-feature correlation exclusions.

Covers the pipeline resolver, the filesystem-safe directory label, the
hashability of the greedy driver's config key, and the benchmark's per-feature
threshold selection (Step 0). These are pure-logic tests -- no image processing.

Run directly (``python tests/test_per_feature_threshold.py``) or under pytest.
"""
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import pandas as pd


def legacy_resolve_threshold_scalar_and_mapping():
    from zbrains.analysis import _resolve_control_correlation_threshold as R
    assert R(0.5, "T1w") == 0.5
    assert R({"T1w": 0.7}, "T1w") == 0.7
    assert R({"T1w": 0.7}, "t1w") == 0.7                 # case-insensitive
    assert R({"T1w": 0.7}, "FLAIR") == float("-inf")     # miss -> keep all
    assert R({"T1w": None}, "T1w") == float("-inf")      # None -> keep all
    assert R({"T1w*blur": 0.6}, "T1w*blur") == 0.6       # blur companion


def legacy_threshold_label_and_output_directory():
    import re
    import zbrains_benchmark as zb
    assert zb.threshold_label(0.7) == "corr0p7"
    lbl = zb.threshold_label({"T1w": 0.7, "FLAIR": None, "T1w*blur": 0.6})
    assert lbl == "corrpf_FLAIRoff_T1w0p7_T1wblur0p6"
    assert re.fullmatch(r"[0-9A-Za-z_]+", lbl), "directory label has unsafe chars"
    # an unwieldy mapping collapses to a stable short hash
    big = {f"feat{i}": 0.7 for i in range(40)}
    assert zb.threshold_label(big).startswith("corrpf_")
    assert len(zb.threshold_label(big)) <= 25
    cfg = dict(normalization="none", use_curvature_covariates=False, predictive_wscore=False,
               wscore_distribution="gaussian", wscore_preprocessing="none",
               wscore_covariate_model="linear", wscore_surface_smoothing_iterations=0,
               blur_depth_model="mean", intensity_depth_model="raw",
               control_correlation_filter=True)
    d = zb.output_directory_for(cfg, "/tmp/pfx", {"T1w": 0.7, "FLAIR": 0.5})
    assert "corrpf_FLAIR0p5_T1w0p7" in d


def legacy_config_key_hashable_with_mapping():
    import zbrains_staged as st
    base = dict(st.SIMPLEST_BASELINE)
    k1 = st._config_key(base, {"T1w": 0.7, "FLAIR": 0.5})
    k2 = st._config_key(base, {"FLAIR": 0.5, "T1w": 0.7})  # order-independent
    assert hash(k1) == hash(k2) and k1 == k2


def _synthetic_arm_details():
    """One frame per FILTER_ARMS label with hand-designed per-feature winners."""
    import results.benchmark_normalization_methods as B
    arms = dict(B.FILTER_ARMS)

    def auc_for(label):
        ov = arms[label]
        t = ov["control_correlation_threshold"] if ov.get("control_correlation_filter") else None
        base = 0.0 if t is None else t
        return {
            "T1w": 0.9 - abs(base - 0.7),                        # peak at 0.7
            "FLAIR": 0.8 if t is None else 0.8 - t,              # best off
            "ADC": 0.7 - abs(base - 0.5),                        # peak at 0.5
            "FA": 0.6 if (t is None or t == 0.3) else 0.5,       # tie off vs 0.3
            "thickness": 0.6 if t in (0.3, 0.5) else 0.5,        # tie 0.3 vs 0.5
        }

    return {lab: pd.DataFrame([{"feature": f, "auc": a} for f, a in auc_for(lab).items()])
            for lab, _ in B.FILTER_ARMS}


def legacy_select_per_feature_thresholds():
    import results.benchmark_normalization_methods as B
    with tempfile.TemporaryDirectory() as td:
        best_filter, pft = B.select_per_feature_thresholds(_synthetic_arm_details(), Path(td))
        assert (Path(td) / "step0_per_feature_threshold_auc.csv").exists()
    assert pft["T1w"] == 0.7
    assert pft["FLAIR"] is None                    # off
    assert pft["ADC"] == 0.5
    assert pft["FA"] is None                       # tie off vs 0.3 -> least aggressive (off)
    assert pft["thickness"] == 0.3                 # tie 0.3 vs 0.5 -> lower threshold
    assert best_filter["control_correlation_filter"] is True
    carried = best_filter["control_correlation_threshold"]
    assert set(carried) == {"T1w", "ADC", "thickness"}  # off features omitted


def legacy_off_sentinel_keeps_all_controls_including_degenerate():
    """-inf (per-feature "off") must keep every control -- including a control
    with an undefined mean correlation -- so a carried-forward "off" feature is
    identical to its measured corr_off arm (control_correlation_filter=False)."""
    import numpy as np
    from zbrains.analysis import _filter_reference_controls_by_correlation as F
    base = np.arange(10, dtype=float)
    # subj0, subj1 correlated; subj2 constant -> zero variance -> NaN mean corr
    data = np.vstack([base, 1.01 * base, np.full(10, 3.0)])
    subs = [("sub-a", "ses-01"), ("sub-b", "ses-01"), ("sub-c", "ses-01")]
    # a finite threshold drops the degenerate (NaN-correlation) control
    _, kept_finite, _, _ = F(data, subs, threshold=0.5, verbose=False)
    assert ("sub-c", "ses-01") not in kept_finite and len(kept_finite) == 2
    # the -inf "off" sentinel keeps ALL controls (byte-for-byte == filter off)
    _, kept_off, _, _ = F(data, subs, threshold=float("-inf"), verbose=False)
    assert len(kept_off) == 3 and set(kept_off) == set(subs)


def legacy_all_off_collapses_to_filter_disabled():
    import results.benchmark_normalization_methods as B
    arms = dict(B.FILTER_ARMS)
    details = {
        lab: pd.DataFrame([{"feature": "T1w",
                            "auc": 0.1 if arms[lab].get("control_correlation_filter") else 0.9}])
        for lab, _ in B.FILTER_ARMS
    }
    with tempfile.TemporaryDirectory() as td:
        best_filter, _ = B.select_per_feature_thresholds(details, Path(td))
    assert best_filter == {"control_correlation_filter": False}


# --- RANK-based (drop-bottom-fraction) control correlation filter ---------------
def test_resolve_quantile_scalar_none_and_mapping():
    from zbrains.analysis import _resolve_control_correlation_quantile as Q
    assert Q(None, "FA") is None                    # off -> keep all
    assert Q(0.1, "FA") == 0.1
    assert Q({"FA": 0.2}, "FA") == 0.2
    assert Q({"FA": 0.2}, "fa") == 0.2              # case-insensitive
    assert Q({"FA": 0.2}, "T1w") is None           # miss -> off for that feature
    assert Q({"FA": None}, "FA") is None           # explicit None -> off


def test_quantile_label_and_output_directory():
    import zbrains_benchmark as zb
    assert zb.quantile_label(0.05) == "corrq0p05"
    assert zb.quantile_label(0.2) == "corrq0p2"
    cfg = dict(normalization="none", use_curvature_covariates=False, predictive_wscore=False,
               wscore_distribution="gaussian", wscore_preprocessing="none",
               wscore_covariate_model="linear", wscore_surface_smoothing_iterations=0,
               blur_depth_model="mean", intensity_depth_model="raw",
               control_correlation_filter=True)
    d05 = zb.output_directory_for(cfg, "/tmp/pfx", 0.05)
    d10 = zb.output_directory_for(cfg, "/tmp/pfx", 0.10)
    assert "corrq0p05" in d05 and "corrq0p1" in d10 and d05 != d10
    assert "corr0p6" not in d05


def test_rank_filter_drops_fixed_fraction_not_whole_feature():
    """RANK mode removes the bottom X% independently for each feature."""
    import numpy as np
    from zbrains.analysis import _filter_reference_controls_by_correlation as F
    rng = np.random.RandomState(0)
    base = rng.randn(80)
    # 40 controls with increasing idiosyncratic noise -> a spread of correlations
    data = np.array([base + (0.02 * i) * rng.randn(80) for i in range(40)])
    subs = [(f"sub-{i:02d}", "ses-01") for i in range(40)]
    counts = {}
    for q in (0.05, 0.10, 0.20):
        _, kept, _, _ = F(data, subs, quantile=q, min_controls=2, verbose=False)
        counts[q] = len(kept)
        assert 2 <= len(kept) < 40                       # never collapses, always prunes some
    assert counts[0.05] >= counts[0.10] >= counts[0.20]  # more aggressive -> fewer kept
    assert abs((40 - counts[0.10]) - 4) <= 2             # ~10% dropped


def test_rank_filter_floor_guarantees_minimum_controls():
    import numpy as np
    from zbrains.analysis import _filter_reference_controls_by_correlation as F
    rng = np.random.RandomState(1)
    base = rng.randn(60)
    data = np.array([base + (0.05 * i) * rng.randn(60) for i in range(12)])
    subs = [(f"s{i}", "ses-01") for i in range(12)]
    # quantile 0.5 would keep ~6, but the floor (min 10, capped at n) protects the reference
    _, kept, _, _ = F(data, subs, quantile=0.5, min_controls=10, verbose=False)
    assert len(kept) >= 10


def test_outlier_stage_uses_rank_quantile_arms():
    import zbrains_staged as zs
    stage = {s.name: s for s in zs.DEFAULT_STAGES}["outlier_detection"]
    arms = dict(stage.candidates)
    assert {"corr_q05", "corr_q10", "corr_q20"}.issubset(arms)
    assert not any(lbl.startswith("corr_0.") for lbl in arms)   # absolute arms removed
    for lbl, frac in [("corr_q05", 0.05), ("corr_q10", 0.10), ("corr_q20", 0.20)]:
        assert arms[lbl]["outlier_method"] == "correlation"
        assert arms[lbl]["control_correlation_filter"] is True
        assert arms[lbl]["control_correlation_quantile"] == frac
        _, warns = zs.validate(dict(zs.SIMPLEST_BASELINE, **arms[lbl]))
        assert not warns


def test_analysis_kwargs_disables_second_analysis_filter_for_staged_arm():
    import zbrains_staged as zs
    thr = zs.zb.DEFAULT_CONTROL_CORRELATION_THRESHOLD
    cfg = dict(zs.SIMPLEST_BASELINE, outlier_method="correlation",
               control_correlation_filter=True,
               control_correlation_quantile=0.1)
    kwargs = zs._analysis_kwargs(cfg, thr)
    assert kwargs["control_correlation_filter"] is False
    assert "control_correlation_quantile" not in kwargs
    # The non-staged API still supports an explicit analysis-time quantile.
    legacy = dict(zs.SIMPLEST_BASELINE, control_correlation_filter=True,
                  control_correlation_quantile=0.1)
    assert zs._analysis_kwargs(legacy, thr)["control_correlation_quantile"] == 0.1
    assert "control_correlation_quantile" not in zs._analysis_kwargs(dict(zs.SIMPLEST_BASELINE), thr)


def test_config_key_separates_quantile_arms():
    import zbrains_staged as zs
    base = dict(zs.SIMPLEST_BASELINE, outlier_method="correlation",
                control_correlation_filter=True)
    keys = {
        zs._config_key(dict(base, control_correlation_quantile=q), 0.6)
        for q in (0.05, 0.10, 0.20)
    }
    assert len(keys) == 3


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"[ok] {name}")
    print("\nALL PER-FEATURE TESTS PASSED")
