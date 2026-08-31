"""Unit tests for the three literature-grounded single-value blur metrics.

All reach SWM2 (2mm) and together cover all 4 depths [mid, white, SWM1, SWM2]:
  gray_white_contrast    = 100*(WM - GM)/(0.5*(WM+GM)), GM=mid, WM=mean(SWM1,SWM2)  [pctsurfcon; Salat 2009]
  boundary_gradient      = OLS slope of I vs mm-depth [0,d,d+1,d+2] over all 4       [Antel 2003/Hong 2014,2017]
  juxtacortical_gradient = 0.5*(SWM1+SWM2) - white                                   [transmantle; MELD sub-boundary]

All are one value per vertex (routed as ordinary features, not the 2-component
gradient_flattening path). Run: python tests/test_blur_gradients.py (or pytest).
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


def _profile(n=50, seed=0):
    rng = np.random.default_rng(seed)
    return np.column_stack([
        rng.normal(30, 3, n),  # midthickness (GM)
        rng.normal(40, 4, n),  # white (boundary)
        rng.normal(50, 5, n),  # SWM1
        rng.normal(55, 6, n),  # SWM2
    ])


def test_registered_and_validated():
    for name in ("gray_white_contrast", "boundary_gradient", "juxtacortical_gradient"):
        assert name in zba.BLUR_DEPTH_MODELS
        assert zba._normalize_blur_depth_model(name) == name
    # aliases + dash forms
    assert zba._normalize_blur_depth_model("boundary") == "boundary_gradient"
    assert zba._normalize_blur_depth_model("gray-white-gradient") == "boundary_gradient"
    assert zba._normalize_blur_depth_model("juxta") == "juxtacortical_gradient"
    assert zba._normalize_blur_depth_model("swm1_gradient") == "juxtacortical_gradient"
    # gray-white contrast aliases (pctsurfcon / GWC / dash forms)
    assert zba._normalize_blur_depth_model("gwc") == "gray_white_contrast"
    assert zba._normalize_blur_depth_model("pctsurfcon") == "gray_white_contrast"
    assert zba._normalize_blur_depth_model("grey-white-contrast") == "gray_white_contrast"


def test_gray_white_contrast_matches_pctsurfcon_formula():
    v = _profile(seed=3)
    got = zba._load_blur_single_gradient(v, "unused", "sub-x", "ses-01", "L", "gray_white_contrast")
    gm = v[:, 0]                                     # GM = midthickness
    wm = 0.5 * (v[:, 2] + v[:, 3])                   # WM = mean(SWM1, SWM2), 1-2mm deep anchor
    expected = 100.0 * (wm - gm) / (0.5 * (wm + gm))
    assert got.shape == (v.shape[0],)
    assert np.allclose(got, expected)


def test_gray_white_contrast_drops_toward_zero_when_blurred():
    # A crisp cortex (clear GM/WM offset) has large |GWC|; a blurred junction
    # (GM ~ deep WM) collapses |GWC| toward zero.
    crisp = np.array([[30.0, 40.0, 50.0, 56.0]])     # WM=mean(50,56)=53 vs GM(mid)=30
    blurred = np.array([[40.0, 40.5, 41.0, 41.2]])   # GM ~ WM
    g_crisp = zba._load_blur_single_gradient(crisp, "u", "s", "e", "L", "gray_white_contrast")
    g_blur = zba._load_blur_single_gradient(blurred, "u", "s", "e", "L", "gray_white_contrast")
    assert abs(g_crisp[0]) > abs(g_blur[0])
    # WM-bright (T1w-like) profile yields positive contrast
    assert g_crisp[0] > 0


def test_gray_white_contrast_near_zero_denominator_is_safe():
    # subject-normalized images can have a near-zero local mean; must not raise/inf.
    # WM=mean(SWM1,SWM2)=-1, GM=+1 -> 0.5*(WM+GM)=0 exactly.
    v = np.array([[1.0, 0.0, -1.0, -1.0]])
    got = zba._load_blur_single_gradient(v, "u", "s", "e", "L", "gray_white_contrast")
    assert np.isfinite(got).all()


def test_juxtacortical_gradient_is_subboundary_band_minus_white():
    v = _profile(seed=1)
    got = zba._load_blur_single_gradient(v, "unused", "sub-x", "ses-01", "L", "juxtacortical_gradient")
    assert got.shape == (v.shape[0],)
    # 1-2mm sub-boundary band mean minus the white boundary; no distance division
    assert np.allclose(got, 0.5 * (v[:, 2] + v[:, 3]) - v[:, 1])


def test_boundary_gradient_is_ols_profile_slope_over_all_four_depths(monkeypatch):
    v = _profile(seed=2)
    dist = np.linspace(0.1, 2.0, v.shape[0])  # small distances are now well-conditioned (no floor)
    monkeypatch.setattr(zba, "_load_gray_white_surface_distance",
                        lambda *a, **k: dist.copy())
    got = zba._load_blur_single_gradient(v, "mica", "sub-x", "ses-01", "L", "boundary_gradient")
    # Closed form must equal an independent per-vertex OLS slope of I vs mm-depth
    # x = [0, d, d+1, d+2] over y = [mid, white, SWM1, SWM2].
    expected = np.array([
        np.polyfit([0.0, d, d + 1.0, d + 2.0], v[i, :4], 1)[0]
        for i, d in enumerate(dist)
    ])
    assert got.shape == (v.shape[0],)
    assert np.isfinite(got).all()
    assert np.allclose(got, expected, atol=1e-9)


def test_boundary_gradient_flattens_toward_zero_when_blurred():
    # A crisp GM->WM profile has a steep slope; a flat (blurred) profile ~0.
    crisp = np.array([[30.0, 40.0, 50.0, 55.0]])
    flat = np.array([[40.0, 40.2, 40.1, 40.3]])
    dist = np.array([1.0])
    orig = zba._load_gray_white_surface_distance
    try:
        zba._load_gray_white_surface_distance = lambda *a, **k: dist.copy()
        g_crisp = zba._load_blur_single_gradient(crisp, "m", "s", "e", "L", "boundary_gradient")
        g_flat = zba._load_blur_single_gradient(flat, "m", "s", "e", "L", "boundary_gradient")
    finally:
        zba._load_gray_white_surface_distance = orig
    assert abs(g_crisp[0]) > abs(g_flat[0])


def test_unique_output_directory_for_blur_models():
    base = dict(normalization="none", use_curvature_covariates=False, predictive_wscore=False,
                wscore_distribution="gaussian", wscore_preprocessing="none",
                wscore_covariate_model="linear", wscore_surface_smoothing_iterations=0,
                intensity_depth_model="raw", control_correlation_filter=False)
    dirs = {}
    for blur in ["gray_white_contrast", "boundary_gradient", "juxtacortical_gradient",
                 "mean", "mean_slope_rms", "gradient_flattening"]:
        cfg = dict(base, blur_depth_model=blur)
        dirs[blur] = zb.output_directory_for(cfg, "/tmp/pfx", 0.5)  # must not KeyError
    assert len(set(dirs.values())) == len(dirs)  # all distinct


if __name__ == "__main__":
    # tiny monkeypatch shim so the module runs without pytest
    class _MP:
        def setattr(self, obj, name, val):
            setattr(obj, name, val)
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        fn = globals()[name]
        fn(_MP()) if "monkeypatch" in fn.__code__.co_varnames else fn()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} BLUR-GRADIENT TESTS PASSED")
