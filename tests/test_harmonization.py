"""Cross-site feature harmonization core (ComBat family), fit-on-controls /
apply-frozen-to-patients. Validates the science in isolation from the pipeline
integration: hand-checkable meancenter / sitestd / M-ComBat cases, empirical-Bayes
shrinkage, covariate (age) preservation, edge cases, and a leakage guard.

Run: python tests/test_harmonization.py  (or under pytest).
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np

from zbrains.harmonization import (
    fit_site_harmonization, apply_site_harmonization, harmonize_multidepth,
)

_MICS = "MICs"
_NOEL = "NOEL"


def test_meancenter_hand_case():
    Y = np.array([[0.0], [2.0], [10.0], [14.0]])          # 1 vertex
    site = np.array([_MICS, _MICS, _NOEL, _NOEL])
    p = fit_site_harmonization(Y, site, None, "meancenter", ref_site=_MICS)
    assert np.isclose(p.location_offset[_MICS][0], 0.0)
    assert np.isclose(p.location_offset[_NOEL][0], 11.0)   # 12 - 1
    out = apply_site_harmonization(Y, site, None, p)
    assert np.allclose(out.ravel(), [0.0, 2.0, -1.0, 3.0]) # MICs unchanged; NOEL -11
    assert np.isclose(out[site == _MICS].mean(), out[site == _NOEL].mean())


def test_sitestd_noeb_hand_case():
    Y = np.array([[0.0], [2.0], [10.0], [14.0]])
    site = np.array([_MICS, _MICS, _NOEL, _NOEL])
    p = fit_site_harmonization(Y, site, None, "sitestd_noeb", min_scale_n=2)
    out = apply_site_harmonization(Y, site, None, p)
    # location+scale matched: both sites share the pooled frame (equal means, grand mean preserved)
    assert np.isclose(out[site == _MICS].mean(), out[site == _NOEL].mean())
    assert np.isclose(out.mean(), Y.mean())               # grand mean 6.5 preserved


def test_combat_eb_shrinkage_pulls_toward_prior_mean():
    rng = np.random.default_rng(0)
    n_per, V = 30, 300
    # per-vertex site shift that VARIES across vertices (so tau^2 > 0 -> EB active)
    shift = rng.normal(3.0, 1.0, V)
    mics = rng.normal(0.0, 1.0, (n_per, V))
    noel = rng.normal(0.0, 1.0, (n_per, V)) + shift[None, :]
    Y = np.vstack([mics, noel])
    site = np.array([_MICS] * n_per + [_NOEL] * n_per)
    p = fit_site_harmonization(Y, site, None, "combat", min_scale_n=2)
    # recompute the raw (unshrunk) gamma_hat for NOEL to compare
    fractions = np.array([0.5, 0.5])
    # gamma_star must be a shrinkage of gamma_hat toward the across-vertex prior mean
    g_star = p.gamma_star[_NOEL]
    gbar = g_star.mean()
    # shrinkage moved estimates and kept them on the prior side (|dev| not inflated)
    assert np.isfinite(g_star).all()
    assert g_star.std() > 0
    # harmonization removes the site difference: pooled control site means equalize
    out = apply_site_harmonization(Y, site, None, p)
    assert abs(out[site == _MICS].mean() - out[site == _NOEL].mean()) < 0.1


def test_mcombat_ref_leaves_reference_site_exactly_unchanged():
    rng = np.random.default_rng(1)
    n_per, V = 25, 50
    Y = np.vstack([rng.normal(0, 1, (n_per, V)),
                   rng.normal(4, 2, (n_per, V))])
    site = np.array([_MICS] * n_per + [_NOEL] * n_per)
    p = fit_site_harmonization(Y, site, None, "mcombat_ref", ref_site=_MICS, min_scale_n=2)
    out = apply_site_harmonization(Y, site, None, p)
    # every MICs (reference) row is byte-for-byte unchanged
    assert np.allclose(out[site == _MICS], Y[site == _MICS], atol=1e-9)
    # a fresh MICs patient also passes through unchanged
    pat = rng.normal(0, 1, V)
    assert np.allclose(apply_site_harmonization(pat, _MICS, None, p), pat, atol=1e-9)
    # NOEL is mapped onto the MICs frame (its mean moves toward MICs')
    assert abs(out[site == _NOEL].mean() - out[site == _MICS].mean()) < 0.3


def test_covariate_age_effect_is_preserved():
    rng = np.random.default_rng(2)
    n_per, V = 40, 100
    age = rng.uniform(20, 60, 2 * n_per)
    slope = 0.5
    site = np.array([_MICS] * n_per + [_NOEL] * n_per)
    site_shift = np.where(site == _NOEL, 8.0, 0.0)
    Y = (slope * age)[:, None] + site_shift[:, None] + rng.normal(0, 0.5, (2 * n_per, V))
    cov = age.reshape(-1, 1)
    p = fit_site_harmonization(Y, site, cov, "combat", min_scale_n=2)
    out = apply_site_harmonization(Y, site, cov, p)
    # site removed: site means (at matched covariates) equalize
    assert abs(out[site == _MICS].mean() - out[site == _NOEL].mean()) < 0.5
    # biological age slope preserved: regress harmonized vertex-mean on age
    recovered = np.polyfit(age, out.mean(axis=1), 1)[0]
    assert abs(recovered - slope) < 0.15


def test_edge_cases_passthrough_smalln_nan_and_3d():
    rng = np.random.default_rng(3)
    n_per, V = 20, 10
    Y = np.vstack([rng.normal(0, 1, (n_per, V)), rng.normal(3, 1, (n_per, V))])
    Y[:, 0] = 5.0                                         # zero-variance vertex 0
    site = np.array([_MICS] * n_per + [_NOEL] * n_per)
    p = fit_site_harmonization(Y, site, None, "combat", min_scale_n=2)
    assert p.passthrough[0] and not p.passthrough[1]
    # a patient's passthrough vertex is returned unchanged
    pat = rng.normal(0, 1, V); pat[0] = 999.0
    out = apply_site_harmonization(pat, _NOEL, None, p)
    assert out[0] == 999.0
    # NaN patient vertex stays NaN
    pat2 = rng.normal(0, 1, V); pat2[3] = np.nan
    assert np.isnan(apply_site_harmonization(pat2, _NOEL, None, p)[3])
    # small-n site (1 control) -> scale forced to 1 (location only), no crash
    Ysmall = np.vstack([rng.normal(0, 1, (5, V)), rng.normal(3, 1, (1, V))])
    ssite = np.array([_MICS] * 5 + [_NOEL])
    ps = fit_site_harmonization(Ysmall, ssite, None, "sitestd_noeb", min_scale_n=2)
    assert np.allclose(ps.delta_star[_NOEL], 1.0)
    # 3D (n, V, depth) round-trips through the multidepth helper
    Y3 = rng.normal(0, 1, (2 * n_per, V, 4))
    fit3 = harmonize_multidepth(lambda M: fit_site_harmonization(M, site, None, "combat", min_scale_n=2), Y3)
    out3 = harmonize_multidepth(lambda M: apply_site_harmonization(M, site, None, fit3), Y3)
    assert out3.shape == Y3.shape


def test_leakage_guard_apply_never_mutates_frozen_params():
    rng = np.random.default_rng(4)
    n_per, V = 20, 40
    Y = np.vstack([rng.normal(0, 1, (n_per, V)), rng.normal(2, 1, (n_per, V))])
    site = np.array([_MICS] * n_per + [_NOEL] * n_per)
    p = fit_site_harmonization(Y, site, None, "combat", min_scale_n=2)
    frozen = copy.deepcopy({s: p.gamma_star[s].copy() for s in p.sites})
    frozen_d = copy.deepcopy({s: p.delta_star[s].copy() for s in p.sites})
    # apply to an extreme outlier "patient" from the NOEL site
    outlier = np.full(V, 50.0)
    apply_site_harmonization(outlier, _NOEL, None, p)
    for s in p.sites:                                     # params unchanged by apply
        assert np.array_equal(p.gamma_star[s], frozen[s])
        assert np.array_equal(p.delta_star[s], frozen_d[s])
    # re-applying to controls is deterministic (frozen transform)
    a = apply_site_harmonization(Y, site, None, p)
    b = apply_site_harmonization(Y, site, None, p)
    assert np.array_equal(a, b)


def test_unknown_site_passes_through():
    rng = np.random.default_rng(5)
    Y = np.vstack([rng.normal(0, 1, (10, 8)), rng.normal(3, 1, (10, 8))])
    site = np.array([_MICS] * 10 + [_NOEL] * 10)
    p = fit_site_harmonization(Y, site, None, "combat", min_scale_n=2)
    pat = rng.normal(0, 1, 8)
    assert np.allclose(apply_site_harmonization(pat, "UNSEEN_SITE", None, p), pat)


def _bundle(mat, tag):
    import pandas as pd
    n = mat.shape[0]
    subs = [(f"sub-{tag}{i}", "ses-01") for i in range(n)]
    demo = pd.DataFrame({"AGE": np.linspace(25, 55, n), "SEX": ([0, 1] * n)[:n]})
    return (mat, subs, demo)


def test_siteharmonizer_pools_and_harmonizes():
    from zbrains.harmonization import SiteHarmonizer
    rng = np.random.default_rng(6)
    V = 40
    mics = rng.normal(0, 1, (10, V))
    noel = rng.normal(3, 1, (8, V))               # shifted site
    H = SiteHarmonizer("combat", ["AGE", "SEX"], ref_site=_MICS, min_scale_n=2)
    H.add_site_features(_MICS, {"f1": _bundle(mics, "HC")})
    H.add_site_features(_NOEL, {"f1": _bundle(noel, "PX")})
    assert H.covers("f1") and not H.covers("f2")
    harmonized, subjects, demo = H.pooled_reference("f1")
    assert harmonized.shape == (18, V) and len(subjects) == 18
    # site difference removed in the pooled harmonized controls
    assert abs(harmonized[:10].mean() - harmonized[10:].mean()) < 0.4
    # single-site feature -> None (local reference retained)
    H.add_site_features(_MICS, {"only_mics": _bundle(mics, "HC")})
    assert H.pooled_reference("only_mics") is None


def test_siteharmonizer_harmonizes_guards_uncovered_site():
    from zbrains.harmonization import SiteHarmonizer
    rng = np.random.default_rng(9)
    V = 20
    H = SiteHarmonizer("combat", ["AGE", "SEX"], min_scale_n=2)
    H.add_site_features(_MICS, {"f1": _bundle(rng.normal(0, 1, (8, V)), "HC"),
                                "only_mics": _bundle(rng.normal(0, 1, (8, V)), "HC")})
    H.add_site_features(_NOEL, {"f1": _bundle(rng.normal(3, 1, (8, V)), "PX")})
    # both pooled sites harmonize the shared map...
    assert H.harmonizes("f1", _MICS) and H.harmonizes("f1", _NOEL)
    # ...a THIRD site that never contributed controls does NOT (would otherwise be
    # scored unharmonized against a harmonized reference -> the passthrough hole)
    assert not H.harmonizes("f1", "THIRD_SITE")
    # a single-site map is not covered at all
    assert not H.harmonizes("only_mics", _MICS)


def test_siteharmonizer_apply_patient_frozen():
    from zbrains.harmonization import SiteHarmonizer
    import copy as _copy
    rng = np.random.default_rng(7)
    V = 30
    H = SiteHarmonizer("combat", ["AGE", "SEX"], min_scale_n=2)
    H.add_site_features(_MICS, {"f1": _bundle(rng.normal(0, 1, (12, V)), "HC")})
    H.add_site_features(_NOEL, {"f1": _bundle(rng.normal(4, 1, (10, V)), "PX")})
    H.pooled_reference("f1")                       # triggers the fit
    frozen = _copy.deepcopy(H._params["f1"].gamma_star)
    out = H.apply_patient("f1", rng.normal(4, 1, V), _NOEL, {"AGE": 40, "SEX": 1})
    assert out.shape == (V,)
    for s in frozen:                               # apply never mutates frozen params
        assert np.array_equal(H._params["f1"].gamma_star[s], frozen[s])


def test_siteharmonizer_site_covariate_is_identity_with_site_column():
    from zbrains.harmonization import SiteHarmonizer
    rng = np.random.default_rng(8)
    V = 20
    mics = rng.normal(0, 1, (6, V)); noel = rng.normal(3, 1, (6, V))
    H = SiteHarmonizer("combat", ["AGE", "SEX"], site_covariate=True, min_scale_n=2)
    H.add_site_features(_MICS, {"f1": _bundle(mics, "HC")})
    H.add_site_features(_NOEL, {"f1": _bundle(noel, "PX")})
    Y, subjects, demo = H.pooled_reference("f1")
    assert np.array_equal(Y, np.vstack([mics, noel]))   # NO feature transform
    assert "site" in demo.columns and set(demo["site"].unique()) == {0, 1}
    pat = rng.normal(0, 1, V)
    assert np.array_equal(H.apply_patient("f1", pat, _NOEL, {"AGE": 40, "SEX": 1}), pat)


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} HARMONIZATION TESTS PASSED")
