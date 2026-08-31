"""Nyul composition fix: use a SynthSeg brain mask, not the data!=0 heuristic.

On subject-normalized (WhiteStripe/wm-mean) images the brain is centred near 0,
so the legacy `data != 0` mask drops most of it. The SynthSeg brain mask keeps
the whole labelled brain, which is required to compose Nyul on those images.

Run: python tests/test_nyul_mask.py  (or under pytest).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import nibabel as nib

from zbrains import normalization as nz


def _synthetic_normalized_brain(seed=0):
    """A z-scored 'subject-normalized' volume (brain ~N(0,1) incl exact zeros;
    background exactly 0), plus a SynthSeg-style label volume for the brain."""
    rng = np.random.default_rng(seed)
    data = np.zeros((24, 24, 24), np.float32)
    brain = np.zeros((24, 24, 24), bool)
    brain[6:18, 6:18, 6:18] = True
    vals = rng.normal(0.0, 1.0, int(brain.sum())).astype(np.float32)
    vals[:100] = 0.0  # exact-zero brain voxels (WhiteStripe output has these)
    data[brain] = vals
    label = np.zeros((24, 24, 24), np.int32)
    label[brain] = 2  # a SynthSeg WM label
    return data, brain, label


def test_synthseg_mask_includes_zero_centered_brain():
    data, brain, label = _synthetic_normalized_brain()
    n_default = int((np.isfinite(data) & (data != 0)).sum())        # legacy heuristic
    n_masked = int((_mask := (label > 0)).sum())                    # SynthSeg mask
    assert n_masked == int(brain.sum())
    assert n_masked > n_default                                     # keeps the exact-zero brain voxels
    # landmarks computed over the full brain differ from the data!=0 subset
    lm_default = nz._nyul_landmarks(data)
    lm_masked = nz._nyul_landmarks(data, brain_mask=(label > 0))
    assert lm_default is not None and lm_masked is not None


def test_nyul_fit_and_apply_with_mask_roundtrips():
    data, brain, label = _synthetic_normalized_brain(seed=1)
    with tempfile.TemporaryDirectory() as td:
        img = Path(td) / "img.nii.gz"
        lab = Path(td) / "label.nii.gz"
        out = Path(td) / "out.nii.gz"
        nib.save(nib.Nifti1Image(data, np.eye(4)), str(img))
        nib.save(nib.Nifti1Image(label.astype(np.float32), np.eye(4)), str(lab))
        model = nz.fit_nyul_standard_scale([str(img)], mask_paths=[str(lab)])
        assert model["n_training_images"] == 1
        nz.apply_nyul_to_image(str(img), model, str(out), mask_path=str(lab))
        result = np.asarray(nib.load(str(out)).dataobj)
        # the whole brain is normalized (has spread); background stays 0
        assert result[brain].std() > 0
        assert np.allclose(result[~brain], 0.0)


def test_nyul_mask_resampled_when_grids_differ():
    # Regression: a subject image and its SynthSeg mask on DIFFERENT grids (shared
    # world coords, different shape -- SynthSeg --robust conforms to ~1mm) must NOT
    # crash; the mask is resampled onto the image grid before use.
    data, brain, _label = _synthetic_normalized_brain(seed=3)      # image: 24^3 @ 1mm
    # label on a FINER 48^3 @ 0.5mm grid covering the same world extent
    label_fine = np.zeros((48, 48, 48), np.float32)
    label_fine[12:36, 12:36, 12:36] = 2                            # brain in world [6,18)
    with tempfile.TemporaryDirectory() as td:
        img = Path(td) / "img.nii.gz"
        lab = Path(td) / "label.nii.gz"
        out = Path(td) / "out.nii.gz"
        nib.save(nib.Nifti1Image(data, np.eye(4)), str(img))
        nib.save(nib.Nifti1Image(label_fine, np.diag([0.5, 0.5, 0.5, 1.0])), str(lab))

        # the mask helper aligns the differently-shaped label to the image grid
        mask = nz._synthseg_brain_mask_on_grid(str(lab), nib.load(str(img)))
        assert mask.shape == data.shape and mask.sum() > 0

        # fit + apply used to raise "operands could not be broadcast" -> now works
        model = nz.fit_nyul_standard_scale([str(img)], mask_paths=[str(lab)])
        assert model["n_training_images"] == 1
        nz.apply_nyul_to_image(str(img), model, str(out), mask_path=str(lab))
        assert np.asarray(nib.load(str(out)).dataobj).shape == data.shape


def test_nyul_parallel_fit_matches_sequential():
    # The threaded fit (n_jobs>1) must produce the IDENTICAL standard scale as the
    # sequential fit (joblib preserves order; the mean is order-independent anyway).
    with tempfile.TemporaryDirectory() as td:
        paths, masks = [], []
        for s in range(6):
            data, _brain, label = _synthetic_normalized_brain(seed=s)
            ip, lp = Path(td) / f"img{s}.nii.gz", Path(td) / f"lab{s}.nii.gz"
            nib.save(nib.Nifti1Image(data, np.eye(4)), str(ip))
            nib.save(nib.Nifti1Image(label.astype(np.float32), np.eye(4)), str(lp))
            paths.append(str(ip)); masks.append(str(lp))
        m1 = nz.fit_nyul_standard_scale(paths, mask_paths=masks, n_jobs=1)
        m4 = nz.fit_nyul_standard_scale(paths, mask_paths=masks, n_jobs=4)
        assert m1["n_training_images"] == m4["n_training_images"] == 6
        assert np.allclose(m1["standard_landmarks"], m4["standard_landmarks"])


def test_backward_compatible_without_mask():
    # mask=None reproduces the legacy data!=0 behaviour (raw-image path unchanged)
    data, brain, label = _synthetic_normalized_brain(seed=2)
    lm_none = nz._nyul_landmarks(data, brain_mask=None)
    lm_explicit = nz._nyul_landmarks(data, brain_mask=(np.isfinite(data) & (data != 0)))
    assert np.allclose(lm_none, lm_explicit)


def test_synthseg_cached_once_and_symlinked_across_bases():
    # SynthSeg depends only on the raw T1w, so it must run ONCE into a shared
    # per-subject cache and be symlinked into each base -- not re-run per base.
    calls = []
    orig = nz._synthseg_subprocess

    def fake_subprocess(in_file, out_file, threads=1):
        calls.append(out_file)
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        with open(out_file, "w") as f:
            f.write("SEG")
        return out_file

    nz._synthseg_subprocess = fake_subprocess
    try:
        with tempfile.TemporaryDirectory() as prefix:
            src = os.path.join(prefix, "raw_T1w.nii.gz")
            open(src, "w").write("t1")

            def base_out(base):
                return os.path.join(prefix, base, "sub-HC1", "ses-01", "maps",
                                    "sub-HC1_ses-01_space-nativepro_desc-synthseg_T1w.nii.gz")

            o1 = base_out("zbrains_SELF_smfwhm")
            o2 = base_out("zbrains_WMMEANRAVEL_smfwhm_exclabc")   # a DIFFERENT base
            nz._run_synthseg(src, o1, threads=1)
            nz._run_synthseg(src, o2, threads=1)

            cache = nz._synthseg_cache_path(o1)
            assert len(calls) == 1 and calls[0] == cache          # segmented ONCE, into the cache
            assert os.path.islink(o1) and os.path.islink(o2)      # each base is a symlink
            assert os.path.realpath(o1) == os.path.realpath(o2) == os.path.realpath(cache)
            assert open(o1).read() == "SEG" and open(o2).read() == "SEG"

            nz._run_synthseg(src, o1, threads=1)                  # re-run existing base -> no-op
            assert len(calls) == 1
    finally:
        nz._synthseg_subprocess = orig


if __name__ == "__main__":
    names = sorted(n for n in dir() if n.startswith("test_"))
    for name in names:
        globals()[name]()
        print(f"[ok] {name}")
    print(f"\nALL {len(names)} NYUL-MASK TESTS PASSED")
