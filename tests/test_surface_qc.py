"""Unit tests for the pre-normalization bilateral surface QC (zbrains.surface_qc).

Covers the acceptance tests from the spec, adapted to the corrected design:
bilateral joint scaling, cross-fitted robust PCA, one combined max-anomaly
statistic ``A`` with a single bootstrap-calibrated ``qc_alpha`` threshold.

Synthetic control cohorts are written as real ``.func.gii`` files so the whole
load -> mask -> joint-scale -> cross-fitted PCA -> anomaly -> calibrate path runs.
Run directly (``python tests/test_surface_qc.py``) or under pytest.
"""
from __future__ import annotations

import hashlib
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd
import nibabel as nib

from zbrains import surface_qc as sq

NVERT = 500
N_CONTROLS = 24


def _basis(nvert: int, n_basis: int = 5) -> np.ndarray:
    t = np.linspace(0.0, np.pi, nvert)
    return np.vstack([np.sin((k + 1) * t) + 0.3 * np.cos((k + 1) * t) for k in range(n_basis)])


def make_normal_map(rng, basis) -> np.ndarray:
    coeff = rng.normal(1.0, 0.25, size=basis.shape[0])
    return coeff @ basis + rng.normal(0.0, 0.05, basis.shape[1]) + 10.0


def write_gifti(path: Path, arr: np.ndarray) -> None:
    img = nib.gifti.GiftiImage()
    img.add_gifti_data_array(nib.gifti.GiftiDataArray(np.asarray(arr, dtype=np.float32)))
    nib.save(img, str(path))


def build_cohort(tmp: Path, feature="T1w", n=N_CONTROLS, nvert=NVERT, seed=0, mutate=None):
    rng = np.random.default_rng(seed)
    basis = _basis(nvert)
    rows = []
    for i in range(n):
        subject = f"sub-HC{i:03d}"
        for hemi in ("L", "R"):
            arr = make_normal_map(rng, basis)
            if mutate is not None:
                arr = mutate(i, hemi, arr)
            path = tmp / f"{subject}_hemi-{hemi}_{feature}.func.gii"
            write_gifti(path, arr)
            rows.append({"subject": subject, "feature": feature, "hemisphere": hemi, "file_path": str(path)})
    return pd.DataFrame(rows)


# --- pure helpers -----------------------------------------------------------
def test_robust_z_and_combined_max():
    v = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 100.0])  # MAD > 0 (primary path)
    z = sq.robust_z(v)
    assert z[-1] == z.max() and z[-1] > 5 and np.all(z[:-1] < 2)
    C = np.array([[1.0, 2.0, 0.5], [0.1, 0.2, 0.3]])
    A = sq.combined_anomaly_statistic(C)
    assert A[0] == 2.0 and abs(A[1] - 0.3) < 1e-9


def test_robust_z_zero_mad_no_blowup():
    # >50% share a value (MAD=0): must fall back to std, not divide by EPSILON.
    v = np.array([0.0] * 17 + [1e-4] * 3)
    z = sq.robust_z(v)
    assert np.all(np.isfinite(z)) and np.max(np.abs(z)) < 10.0
    assert np.allclose(sq.robust_z(np.full(12, 5.0)), 0.0)  # truly constant -> 0


def test_calibration_robust_to_high_contamination():
    # 40% of subjects are outliers in one component: a contaminated MinCovDet
    # covariance must NOT inflate the threshold so far that nothing flags.
    rng = np.random.default_rng(0)
    C = rng.normal(size=(25, 3))
    C[:10, 0] += 10.0  # 40% clustered outliers in component 0
    A, thr, pval, status = sq.calibrate_anomaly_threshold(C, 0.01, n_boot=8000, random_seed=1)
    assert np.isfinite(thr) and thr < 6.0          # threshold not exploded
    assert int(np.sum(A > thr)) >= 10              # the 10 obvious outliers flag


def test_numeric_group_by_does_not_crash():
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=20)
        manifest["site"] = 1  # numeric acquisition-group column (realistic CSV input)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01,
                              group_by=["site"], minimum_controls=15, make_figures=False)
    assert len(r.table) == N_CONTROLS


def test_calibrate_bounds_and_threshold():
    rng = np.random.default_rng(0)
    comps = rng.normal(size=(200, 3))
    A, thr, pval, status = sq.calibrate_anomaly_threshold(comps, 0.01, n_boot=5000, random_seed=1)
    assert np.isfinite(thr) and thr > 0
    assert np.all((pval > 0) & (pval <= 1))
    assert A.shape[0] == 200 and "bootstrap" in status
    # analytic path is bounded too
    _, thr2, pv2, st2 = sq.calibrate_anomaly_threshold(comps, 0.01, calibration="analytic")
    assert np.isfinite(thr2) and np.all((pv2 > 0) & (pv2 <= 1)) and "analytic" in st2


# --- bilateral joint scaling (correction's key property) --------------------
def test_bilateral_joint_scaling_invariance_and_unilateral_sensitivity():
    rng = np.random.default_rng(3)
    basis = _basis(NVERT)
    L = make_normal_map(rng, basis)
    R = make_normal_map(rng, basis)
    zL, zR, _, _ = sq.robust_bilateral_scale(L, R)
    # scaling BOTH hemispheres by the same factor / offset is invariant
    zL5, zR5, _, _ = sq.robust_bilateral_scale(5.0 * L, 5.0 * R)
    zLo, zRo, _, _ = sq.robust_bilateral_scale(L + 100.0, R + 100.0)
    assert np.allclose(zL, zL5, atol=1e-6) and np.allclose(zR, zR5, atol=1e-6)
    assert np.allclose(zL, zLo, atol=1e-6) and np.allclose(zR, zRo, atol=1e-6)
    # attenuating ONE hemisphere is NOT hidden by joint scaling
    zLu, zRu, _, _ = sq.robust_bilateral_scale(0.5 * L, R)
    assert not np.allclose(zLu, zL, atol=1e-3)


def test_bilateral_concatenation_dimensions():
    rng = np.random.default_rng(4)
    basis = _basis(NVERT)
    left = [make_normal_map(rng, basis) for _ in range(N_CONTROLS)]
    right = [make_normal_map(rng, basis) for _ in range(N_CONTROLS)]
    mask_l = sq.build_common_vertex_mask(left)
    mask_r = sq.build_common_vertex_mask(right)
    zl, zr, _, _ = sq.robust_bilateral_scale(left[0][mask_l], right[0][mask_r])
    assert np.concatenate([zl, zr]).size == int(mask_l.sum()) + int(mask_r.sum())


def test_robust_pca_fallback_documented():
    rng = np.random.default_rng(5)
    X = rng.normal(size=(24, 400))
    model = sq.fit_robust_pca(X, variance_explained=0.9, max_components=15, random_seed=42)
    assert "robpca_equiv" in model.status  # never silent about using the fallback
    od, sd = model.distances(X[:3])
    assert od.shape == (3,) and sd.shape == (3,)


# --- full pipeline ----------------------------------------------------------
def test_normal_cohort_low_false_positives_and_reproducible():
    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=10)
        r1 = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01, make_figures=False)
        r2 = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01, make_figures=False)
    assert int(r1.table["flagged"].sum()) <= 1          # clean cohort, alpha=1%
    pd.testing.assert_series_equal(r1.table["anomaly_score"], r2.table["anomaly_score"])
    pd.testing.assert_series_equal(r1.table["anomaly_pvalue"], r2.table["anomaly_pvalue"])


def test_localized_artifact_high_reconstruction_distance_and_flagged():
    def mutate(i, hemi, arr):
        if i == 0:
            arr = arr.copy(); arr[50:70] += 25.0
        return arr

    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=11, mutate=mutate)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.05, make_figures=False)
    t = r.table.set_index("subject")
    assert t["orthogonal_distance_bilateral"].idxmax() == "sub-HC000"
    assert bool(t.loc["sub-HC000", "flagged"])


def test_global_scaling_detected_by_raw_intensity():
    def mutate(i, hemi, arr):
        return arr * 5.0 if i == 0 else arr  # both hemispheres scaled

    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=12, mutate=mutate)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.05, make_figures=False)
    t = r.table.set_index("subject")
    # global scaling is invisible to the bilaterally-scaled spatial map but the
    # raw-intensity component must catch it
    assert t.loc["sub-HC000", "raw_median_left"] > 1.5 * t["raw_median_left"].drop("sub-HC000").median()
    assert t.loc["sub-HC000", "z_raw"] == t["z_raw"].max()
    assert bool(t.loc["sub-HC000", "flagged"])


def test_unilateral_attenuation_detected():
    def mutate(i, hemi, arr):
        if i == 0 and hemi == "L":
            return 0.5 * arr  # one-sided global attenuation
        return arr

    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=18, mutate=mutate)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.05, make_figures=False)
    t = r.table.set_index("subject")
    assert t["left_right_median_difference"].idxmax() == "sub-HC000"
    assert bool(t.loc["sub-HC000", "flagged"])


def test_unilateral_structural_failure_excludes_both_hemispheres():
    def mutate(i, hemi, arr):
        if i == 0 and hemi == "L":
            return np.full_like(arr, 7.0)  # constant one hemisphere
        return arr

    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=13, mutate=mutate)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01,
                              make_figures=False, write_exclusion_list=True)
    row = r.table[r.table["subject"] == "sub-HC000"].iloc[0]
    assert bool(row["flagged"])
    excl = r.exclusions[r.exclusions["subject"] == "sub-HC000"]
    assert len(excl) == 1 and excl.iloc[0]["hemisphere_action"] == "remove_both"
    assert excl.iloc[0]["left_file"] and excl.iloc[0]["right_file"]


def test_constant_map_detected_before_pca():
    def mutate(i, hemi, arr):
        return np.full_like(arr, 3.0) if i == 0 else arr

    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=14, mutate=mutate)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01, make_figures=False)
    row = r.table[r.table["subject"] == "sub-HC000"].iloc[0]
    assert "constant_map" in row["flag_reason"]
    assert not np.isfinite(row["orthogonal_distance_bilateral"])  # excluded from PCA


def test_inconsistent_vertex_count_flagged():
    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=15)
        bad = manifest[(manifest["subject"] == "sub-HC000") & (manifest["hemisphere"] == "L")].iloc[0]
        write_gifti(Path(bad["file_path"]), np.random.default_rng(0).normal(size=NVERT - 100))
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01, make_figures=False)
    row = r.table[r.table["subject"] == "sub-HC000"].iloc[0]
    assert "unexpected_vertex_count" in row["flag_reason"]


def test_missing_site_does_not_crash():
    with tempfile.TemporaryDirectory() as td:
        manifest = build_cohort(Path(td), seed=16)
        r = sq.run_surface_qc(manifest, output_dir=None, qc_alpha=0.01,
                              covariates=["age", "sex", "site"], cross_fit_covariates=True,
                              group_by=["site"], make_figures=False)
    assert len(r.table) == N_CONTROLS


def _hash(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def test_inputs_never_modified_and_outputs_written():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        manifest = build_cohort(tmp, seed=17)
        before = {p: _hash(p) for p in manifest["file_path"]}
        out = tmp / "qc_out"
        sq.run_surface_qc(manifest, output_dir=out, qc_alpha=0.01,
                          make_figures=True, write_exclusion_list=True)
        after = {p: _hash(p) for p in manifest["file_path"]}
        assert before == after
        assert (out / "pre_normalization_surface_qc.tsv").exists()
        assert (out / "pre_normalization_surface_qc_summary.json").exists()
        assert (out / "pre_normalization_surface_qc_exclusions.tsv").exists()


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} SURFACE-QC TESTS PASSED")
