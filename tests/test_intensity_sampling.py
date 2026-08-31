"""Unit tests for the plain intensity SAMPLING intensity_depth_models.

Covers: the 4-depth column selection, the swm2 and own-cortex self-norms,
consistency with the legacy white_swm2_zscore / white_surface_robust_z, validation
of the new names, and unique output-directory naming per model.

Run: python tests/test_intensity_sampling.py  (or under pytest)
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np

from zbrains import analysis as zba
import zbrains_benchmark as zb

SAMPLE_NAMES = [
    "sample_midthickness", "sample_white", "sample_swm1",
    "sample_midthickness_swm2", "sample_white_swm2", "sample_swm1_swm2",
    "sample_midthickness_owncortex", "sample_white_owncortex", "sample_swm1_owncortex",
]


def _profile(n=60, seed=0):
    rng = np.random.default_rng(seed)
    # columns: [midthickness, white, SWM1, SWM2], deliberately distinct
    return np.column_stack([
        rng.normal(30, 3, n),   # midthickness (GM)
        rng.normal(40, 4, n),   # white (boundary)
        rng.normal(50, 5, n),   # SWM1
        rng.normal(55, 6, n),   # SWM2
    ])


def test_all_sample_models_registered_and_valid():
    for name in SAMPLE_NAMES:
        assert name in zba.INTENSITY_DEPTH_MODELS
        assert name in zba._INTENSITY_DEPTH_TRANSFORMS
        assert zba._normalize_intensity_depth_model(name) == name
    # dash form normalizes to underscore canonical
    assert zba._normalize_intensity_depth_model("sample-white") == "sample_white"


def test_plain_samplers_select_correct_column():
    v = _profile()
    assert np.allclose(zba._INTENSITY_DEPTH_TRANSFORMS["sample_midthickness"](v), v[:, 0])
    assert np.allclose(zba._INTENSITY_DEPTH_TRANSFORMS["sample_white"](v), v[:, 1])
    assert np.allclose(zba._INTENSITY_DEPTH_TRANSFORMS["sample_swm1"](v), v[:, 2])


def test_swm2_selfnorm_matches_formula_and_legacy():
    v = _profile(seed=1)
    swm2 = v[:, 3]
    center = np.median(swm2)
    scale = 1.482602218505602 * np.median(np.abs(swm2 - center))
    expected = (v[:, 1] - center) / scale
    got = zba._INTENSITY_DEPTH_TRANSFORMS["sample_white_swm2"](v)
    assert np.allclose(got, expected)
    # white+swm2 == the legacy white_swm2_zscore
    assert np.allclose(got, zba._white_swm2_zscore(v))
    # and the midthickness/swm1 variants anchor to the SAME swm2 reference
    assert np.allclose(zba._INTENSITY_DEPTH_TRANSFORMS["sample_swm1_swm2"](v),
                       (v[:, 2] - center) / scale)


def test_owncortex_selfnorm_matches_robust_standardize_and_legacy():
    v = _profile(seed=2)
    got_white = zba._INTENSITY_DEPTH_TRANSFORMS["sample_white_owncortex"](v)
    assert np.allclose(got_white, zba._robust_standardize(v[:, 1]))
    assert np.allclose(got_white, zba._white_surface_robust_z(v))  # legacy equivalence
    got_mid = zba._INTENSITY_DEPTH_TRANSFORMS["sample_midthickness_owncortex"](v)
    assert np.allclose(got_mid, zba._robust_standardize(v[:, 0]))
    # own-cortex is centered near 0 with unit-ish robust spread
    assert abs(float(np.median(got_mid))) < 1e-6


def test_unique_output_directory_per_model():
    base = dict(normalization="none", use_curvature_covariates=False, predictive_wscore=False,
                wscore_distribution="gaussian", wscore_preprocessing="none",
                wscore_covariate_model="linear", wscore_surface_smoothing_iterations=0,
                blur_depth_model="mean", control_correlation_filter=False)
    dirs = {}
    for idm in ["raw"] + SAMPLE_NAMES:
        cfg = dict(base, intensity_depth_model=idm)
        dirs[idm] = zb.output_directory_for(cfg, "/tmp/pfx", 0.5)
    # all distinct
    assert len(set(dirs.values())) == len(dirs)
    # raw has no idm suffix; the samplers all carry an "idm" token
    assert "idm" not in dirs["raw"].rsplit("/", 2)[-2]
    for name in SAMPLE_NAMES:
        assert "idm" in dirs[name]


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} INTENSITY-SAMPLING TESTS PASSED")
