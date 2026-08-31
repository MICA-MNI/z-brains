"""Fix D(a): the GP prediction-uncertainty filter must mask vertices the normative
model predicts POORLY (high leave-one-out residual dispersion), NOT the highest
control-amplitude cortex (which is where FCD gray-white blurring lives). This test
builds controls with two vertex classes of MATCHED amplitude but opposite model
reliability and asserts the filter masks the unreliable ones.

Run: python tests/test_gp_uncertainty.py  (or under pytest).
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd
import nibabel as nib

from zbrains.gaussian_process import calculate_gaussian_process_maps


def test_gp_uncertainty_masks_poorly_predicted_not_high_amplitude():
    rng = np.random.default_rng(0)
    n_ctrl = 20
    age = rng.uniform(20, 60, n_ctrl)
    sex = rng.integers(0, 2, n_ctrl).astype(float)
    demo_ref = pd.DataFrame({"AGE": age, "SEX": sex})
    demo_pat = pd.Series({"AGE": 40.0, "SEX": 1.0})

    n_each = 15
    # PREDICTABLE vertices: strong linear age effect (high amplitude) + small noise
    # -> GP predicts held-out controls well -> LOW LOO dispersion.
    slope = 3.0
    predictable = slope * (age[:, None] - age.mean()) + rng.normal(0, 0.3, (n_ctrl, n_each))
    # UNPREDICTABLE vertices: pure noise, matched (actually larger) amplitude, no
    # demographic structure -> GP cannot predict -> HIGH LOO dispersion.
    unpredictable = rng.normal(0, predictable.std() * 1.2, (n_ctrl, n_each))
    reference_data = np.hstack([predictable, unpredictable])       # (n_ctrl, 30)
    patient_data = np.zeros(2 * n_each)                            # patient = grand mean

    with tempfile.TemporaryDirectory() as td:
        out = str(Path(td) / "gp.func.gii")
        calculate_gaussian_process_maps(
            reference_data=reference_data, patient_data=patient_data,
            demographics_ref=demo_ref, demographics_pat=demo_pat, output_file=out,
            normative_columns=["AGE", "SEX"], verbose=False,
            gaussian_process_kernel="rbf_ard", wscore_covariate_model="gaussian_process",
            min_reference_subjects=None, reference_vertex_covariates=None,
            patient_vertex_covariates=None, wscore_preprocessing="none",
            blur_depth_model="mean", intensity_depth_model="raw",
            wscore_surface_smoothing_iterations=0, wscore_fit_cache=None,
            surface_smoother=lambda s, n: s, prediction_variance_percentile=50.0,
        )
        scores = np.asarray(nib.load(out).darrays[0].data, dtype=float)

    assert scores.shape == (2 * n_each,)
    masked = ~np.isfinite(scores)
    # the unpredictable (high-LOO) half should be preferentially masked; the
    # predictable high-amplitude half should be preferentially KEPT.
    frac_masked_unpredictable = masked[n_each:].mean()
    frac_masked_predictable = masked[:n_each].mean()
    assert frac_masked_unpredictable > frac_masked_predictable, (
        frac_masked_predictable, frac_masked_unpredictable)
    # and the amplitude-vs-mask correlation is NOT what the old control-SD^2 filter
    # would produce (it would mask the high-amplitude predictable vertices).
    assert frac_masked_predictable < 0.5


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} GP-UNCERTAINTY TESTS PASSED")
