"""kNN nonparametric normative scoring (k = half the controls), added alongside
z-score / w-score / Gaussian process.

Covers: registration as a wscore covariate model (routed, NOT a GP kernel), the
greedy scoring stage option, the k=half-controls rule, per-vertex 1-D and joint
multi-depth abnormality (injected outliers dominate), leave-one-out calibration,
the degenerate-control guard, and the end-to-end calculate_knn_maps writing a
GIFTI + model JSON with a compatible return dict.

Run: python tests/test_knn_scoring.py  (or under pytest).
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd

from zbrains import analysis as zba
from zbrains.knn_normative import knn_abnormality, calculate_knn_maps
import zbrains_staged as zs
import zbrains_benchmark as zb


def test_knn_registered_as_covariate_model_not_gp_kernel():
    assert "knn" in zba.WSCORE_COVARIATE_MODELS
    assert zba._normalize_wscore_covariate_model("knn") == "knn"
    assert zba.GAUSSIAN_PROCESS_KERNELS.get("knn") is None  # won't route to GP


def test_scoring_stage_has_knn_and_validates():
    labels = [l for l, _ in {s.name: s for s in zs.DEFAULT_STAGES}["scoring"].candidates]
    assert "knn_half_controls" in labels
    cfg = dict(zs.SIMPLEST_BASELINE, method="wscore", wscore_covariate_model="knn")
    merged, warns = zs.validate(cfg)
    assert not warns
    # distinct analysis dir (knn appears as the covariate-model token)
    assert "knn" in zb.output_directory_for(zs._axis_config(merged), "PFX", 0.5)


def test_k_is_half_the_controls():
    rng = np.random.default_rng(1)
    for n_ref in (4, 20, 51):
        ref = rng.normal(0, 1, (n_ref, 3))
        pat = rng.normal(0, 1, 3)
        # calculate_knn_maps sets k=max(1, n_ref//2); verify via the model json
        with tempfile.TemporaryDirectory() as td:
            out = str(Path(td) / "m.func.gii")
            demo_ref = pd.DataFrame({"AGE": rng.integers(20, 60, n_ref),
                                     "SEX": rng.integers(0, 2, n_ref)})
            demo_pat = pd.Series({"AGE": 40, "SEX": 1})
            res = calculate_knn_maps(
                ref, pat, demo_ref, demo_pat, out, ["AGE", "SEX"], False, "knn",
                None, None, None, "none", "boundary_gradient", "raw", 0, None,
                lambda s, n: s, None)
            assert res["k_neighbours"] == max(1, n_ref // 2), (n_ref, res["k_neighbours"])


def test_outlier_vertex_dominates_1d_and_multidepth():
    rng = np.random.default_rng(2)
    n_ref = 20
    k = n_ref // 2
    ref = rng.normal(0, 1, (n_ref, 6))
    pat = rng.normal(0, 1, 6)
    pat[0] = 8.0                                  # gross 1-D outlier
    s = knn_abnormality(ref, pat, k)
    assert np.all(np.isfinite(s))
    assert s[0] > 3 * np.max(s[1:])
    # multi-depth: joint over 4 depths; abnormal depth-vector flagged
    refd = rng.normal(0, 1, (n_ref, 4, 4))
    patd = rng.normal(0, 1, (4, 4))
    patd[0] += 6.0
    sd = knn_abnormality(refd, patd, k)
    assert sd[0] > 2 * np.max(sd[1:])


def test_typical_patient_scores_are_moderate():
    # a patient drawn from the control distribution should not be flagged strongly
    rng = np.random.default_rng(3)
    n_ref, n_vert = 30, 200
    ref = rng.normal(0, 1, (n_ref, n_vert))
    pat = rng.normal(0, 1, n_vert)
    s = knn_abnormality(ref, pat, n_ref // 2)
    # calibrated against the control LOO null, a typical subject sits near the null
    assert np.nanmedian(np.abs(s)) < 5.0


def test_degenerate_controls_give_zero_not_false_signal():
    flat = np.ones((6, 4))
    pat = np.array([1.0, 1.0, 1.0, 9.0])
    s = knn_abnormality(flat, pat, 3)
    assert np.all(s == 0.0)


def test_calculate_knn_maps_writes_gifti_and_model_json():
    import nibabel as nib
    rng = np.random.default_rng(4)
    n_ref, n_vert = 16, 50
    ref = rng.normal(0, 1, (n_ref, n_vert))
    pat = rng.normal(0, 1, n_vert)
    pat[10] = 7.0
    demo_ref = pd.DataFrame({"AGE": rng.integers(20, 60, n_ref),
                             "SEX": rng.integers(0, 2, n_ref)})
    demo_pat = pd.Series({"AGE": 35, "SEX": 0})
    with tempfile.TemporaryDirectory() as td:
        out = str(Path(td) / "sub-x" / "knn.func.gii")
        res = calculate_knn_maps(
            ref, pat, demo_ref, demo_pat, out, ["AGE", "SEX"], False, "knn",
            None, None, None, "none", "boundary_gradient", "raw", 0, None,
            lambda s, n: s, None)
        assert res["wscore_file"] == out and Path(out).exists()
        vals = np.asarray(nib.load(out).darrays[0].data)
        assert vals.shape == (n_vert,)
        assert np.nanargmax(vals) == 10                      # the outlier vertex
        model = json.loads(Path(out.replace(".func.gii", "_model.json")).read_text())
        assert model["method"] == "knn" and model["k_neighbours"] == n_ref // 2


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} kNN-SCORING TESTS PASSED")
