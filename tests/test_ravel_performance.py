"""Regression tests for result-preserving RAVEL execution optimizations."""

import os
import shutil
from pathlib import Path

import nibabel as nib
import numpy as np

from zbrains import normalization as nz


def _save_nifti(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(np.asarray(data), np.eye(4)), str(path))
    return str(path)


def _make_subject_inputs(base, micapipe, participant="sub-01", session="ses-01"):
    bids_id = f"{participant}_{session}"
    maps = Path(base) / participant / session / "maps"
    maps.mkdir(parents=True, exist_ok=True)
    image = np.arange(27, dtype=np.float32).reshape(3, 3, 3)
    csf = np.ones((3, 3, 3), dtype=np.float32)
    _save_nifti(
        maps / f"{bids_id}_space-nativepro_desc-whitestripe_T1w.nii.gz",
        image,
    )
    _save_nifti(
        maps / f"{bids_id}_space-nativepro_desc-synthsegCsf_T1w.nii.gz",
        csf,
    )

    anat = Path(micapipe) / participant / session / "anat"
    xfm = Path(micapipe) / participant / session / "xfm"
    anat.mkdir(parents=True, exist_ok=True)
    xfm.mkdir(parents=True, exist_ok=True)
    _save_nifti(anat / f"{bids_id}_space-nativepro_T1w.nii.gz", image)
    prefix = (
        xfm
        / f"{bids_id}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN"
    )
    Path(f"{prefix}_0GenericAffine.mat").write_text("affine", encoding="utf-8")
    Path(f"{prefix}_1Warp.nii.gz").write_bytes(b"warp")
    Path(f"{prefix}_1InverseWarp.nii.gz").write_bytes(b"inverse-warp")


def test_shared_mni_cache_reuses_content_across_fold_bases(tmp_path, monkeypatch):
    reference = _save_nifti(
        tmp_path / "reference.nii.gz",
        np.ones((3, 3, 3), dtype=np.float32),
    )
    micapipe = tmp_path / "micapipe"
    base_a = tmp_path / "zbrains_method_reffold0"
    base_b = tmp_path / "zbrains_method_reffold1"
    _make_subject_inputs(base_a, micapipe)
    _make_subject_inputs(base_b, micapipe)

    calls = []

    def fake_forward(moving, output, fixed, transforms, inversions, **kwargs):
        calls.append((moving, kwargs["interpolator"]))
        shutil.copy2(moving, output)
        return output

    monkeypatch.setattr(nz, "_reference_image_path", lambda: reference)
    monkeypatch.setattr(nz, "_forward_transform_ravel_cli", fake_forward)

    first = nz._ensure_mni_inputs(
        "sub-01",
        "ses-01",
        str(base_a),
        str(micapipe),
        ["T1w"],
        verbose=False,
    )
    second = nz._ensure_mni_inputs(
        "sub-01",
        "ses-01",
        str(base_b),
        str(micapipe),
        ["T1w"],
        verbose=False,
    )

    assert len(calls) == 2  # image + CSF once, not once per fold
    assert (
        first["modalities"]["T1w"]["mni_whitestripe"]
        == second["modalities"]["T1w"]["mni_whitestripe"]
    )
    assert "ravel_mni_cache" in first["modalities"]["T1w"]["mni_whitestripe"]

    changed_input = (
        base_b
        / "sub-01"
        / "ses-01"
        / "maps"
        / "sub-01_ses-01_space-nativepro_desc-whitestripe_T1w.nii.gz"
    )
    _save_nifti(changed_input, np.full((3, 3, 3), 42, dtype=np.float32))
    nz._ensure_mni_inputs(
        "sub-01",
        "ses-01",
        str(base_b),
        str(micapipe),
        ["T1w"],
        verbose=False,
    )
    assert len(calls) == 4  # changed content invalidates both coupled cache files


def test_forward_transform_cli_preserves_transform_order_and_thread_budget(
    tmp_path,
    monkeypatch,
):
    output = tmp_path / "out.nii.gz"
    captured = {}

    def fake_run(command, check, env):
        captured["command"] = command
        captured["check"] = check
        captured["env"] = env
        _save_nifti(
            command[command.index("-o") + 1],
            np.full((2, 2, 2), 1.25, dtype=np.float64),
        )

    monkeypatch.setattr(nz.shutil, "which", lambda name: "/fake/antsApplyTransforms")
    monkeypatch.setattr(nz.subprocess, "run", fake_run)
    nz._forward_transform_ravel_cli(
        "moving.nii.gz",
        str(output),
        "fixed.nii.gz",
        ["warp.nii.gz", "affine.mat"],
        [False, True],
        interpolator="NearestNeighbor",
        threads=3,
        atomic=True,
    )

    transformed = nib.load(output)
    assert transformed.get_data_dtype() == np.dtype(np.float32)
    assert np.all(np.asarray(transformed.dataobj) == np.float32(1.25))
    assert captured["command"][-4:] == [
        "-t",
        "warp.nii.gz",
        "-t",
        "[affine.mat,1]",
    ]
    assert captured["command"][captured["command"].index("-n") + 1] == "NearestNeighbor"
    assert captured["env"]["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] == "3"


def test_preloaded_model_produces_identical_patient_output(tmp_path, monkeypatch):
    output_directory = tmp_path / "output"
    ravel_dir = output_directory / "ravel"
    ravel_dir.mkdir(parents=True)
    brain_mask = np.ones((2, 2, 2), dtype=bool)
    control_mask = np.array([True, True, True, True, False, False, False, False])
    beta_x = np.linspace(0.1, 0.8, 8, dtype=np.float32)[None, :]
    beta_w = np.linspace(0.2, 0.9, 8, dtype=np.float32)[None, :]
    np.savez_compressed(
        nz.ravel_model_path(str(output_directory), "T1w"),
        brain_mask=brain_mask,
        control_mask=control_mask,
        beta_x=beta_x,
        beta_w=beta_w,
        model_version=np.array(["chunked_gram_v1"]),
    )
    reference = _save_nifti(
        tmp_path / "reference.nii.gz",
        np.ones(brain_mask.shape, dtype=np.float32),
    )
    patient_input = _save_nifti(
        tmp_path / "patient_input.nii.gz",
        np.linspace(1, 8, 8, dtype=np.float32).reshape(brain_mask.shape),
    )

    outputs = {}

    def fake_inputs(participant, session, *args, **kwargs):
        mni_output = str(tmp_path / f"{session}_mni.nii.gz")
        native_output = str(tmp_path / f"{session}_native.nii.gz")
        outputs[session] = mni_output
        return {
            "reference": reference,
            "native_reference": reference,
            "mni_to_native": [],
            "mni_to_native_invert": [],
            "modalities": {
                "T1w": {
                    "mni_whitestripe": patient_input,
                    "mni_ravel": mni_output,
                    "native_ravel": native_output,
                }
            },
        }

    monkeypatch.setattr(nz, "_reference_image_path", lambda: reference)
    monkeypatch.setattr(nz, "_ensure_mni_inputs", fake_inputs)
    monkeypatch.setattr(nz, "_inverse_transform_ravel_cli", lambda *args, **kwargs: args[1])

    nz.apply_ravel_model_to_subject(
        "sub-01",
        "loaded-from-disk",
        str(output_directory),
        "unused",
        ["T1w"],
        verbose=False,
    )
    from_disk = np.asarray(nib.load(outputs["loaded-from-disk"]).dataobj)

    preloaded = nz.load_ravel_models(str(output_directory), ["T1w"])
    monkeypatch.setattr(
        nz,
        "load_ravel_models",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("preloaded application must not reload the model")
        ),
    )
    nz.apply_ravel_model_to_subject(
        "sub-01",
        "preloaded",
        str(output_directory),
        "unused",
        ["T1w"],
        verbose=False,
        models=preloaded,
    )
    reused = np.asarray(nib.load(outputs["preloaded"]).dataobj)

    assert np.array_equal(from_disk, reused)
    assert all(
        not value.flags.writeable
        for value in preloaded["T1w"].values()
    )


def test_parallel_control_preparation_keeps_fitted_model_bit_identical(
    tmp_path,
    monkeypatch,
):
    shape = (6, 6, 6)
    reference = _save_nifti(
        tmp_path / "fit_reference.nii.gz",
        np.ones(shape, dtype=np.float32),
    )
    subjects = [(f"sub-{index:02d}", "ses-01") for index in range(3)]
    input_paths = {}
    csf_paths = {}
    base_data = np.arange(np.prod(shape), dtype=np.float32).reshape(shape)
    for index, subject in enumerate(subjects):
        input_paths[subject] = _save_nifti(
            tmp_path / "inputs" / f"{subject[0]}.nii.gz",
            base_data * np.float32(index + 1) + np.float32(index * 0.25),
        )
        csf_paths[subject] = _save_nifti(
            tmp_path / "csf" / f"{subject[0]}.nii.gz",
            np.ones(shape, dtype=np.float32),
        )

    def fake_inputs(participant, session, output_directory, *args, **kwargs):
        subject = (participant, session)
        subject_dir = Path(output_directory) / participant / session
        subject_dir.mkdir(parents=True, exist_ok=True)
        return {
            "native_reference": reference,
            "mni_to_native": [],
            "mni_to_native_invert": [],
            "modalities": {
                "T1w": {
                    "mni_whitestripe": input_paths[subject],
                    "mni_csf": csf_paths[subject],
                    "mni_ravel": str(subject_dir / "corrected_mni.nii.gz"),
                    "native_ravel": str(subject_dir / "corrected_native.nii.gz"),
                }
            },
        }

    monkeypatch.setattr(nz, "_reference_image_path", lambda: reference)
    monkeypatch.setattr(nz, "_ensure_mni_inputs", fake_inputs)
    monkeypatch.setattr(nz, "_inverse_transform_ravel_cli", lambda *args, **kwargs: args[1])

    serial = tmp_path / "serial"
    parallel = tmp_path / "parallel"
    nz.fit_and_apply_ravel_to_controls(
        subjects,
        str(serial),
        "unused",
        ["T1w"],
        n_jobs=1,
        threads=1,
        total_threads=1,
        tmp_dir=str(tmp_path / "serial_scratch"),
        verbose=False,
    )
    nz.fit_and_apply_ravel_to_controls(
        subjects,
        str(parallel),
        "unused",
        ["T1w"],
        n_jobs=3,
        threads=1,
        total_threads=3,
        tmp_dir=str(tmp_path / "parallel_scratch"),
        verbose=False,
    )

    with np.load(nz.ravel_model_path(str(serial), "T1w")) as serial_model, np.load(
        nz.ravel_model_path(str(parallel), "T1w")
    ) as parallel_model:
        assert serial_model.files == parallel_model.files
        for key in serial_model.files:
            assert np.array_equal(serial_model[key], parallel_model[key]), key

    for participant, session in subjects:
        serial_data = np.asarray(
            nib.load(serial / participant / session / "corrected_mni.nii.gz").dataobj
        )
        parallel_data = np.asarray(
            nib.load(parallel / participant / session / "corrected_mni.nii.gz").dataobj
        )
        assert np.array_equal(serial_data, parallel_data)


def test_ravel_scratch_prefers_explicit_then_environment(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    environment = tmp_path / "environment"
    monkeypatch.setenv("ZBRAINS_RAVEL_TMPDIR", str(environment))
    assert nz._resolve_ravel_tmp_parent("/unused", str(explicit)) == str(explicit)
    assert nz._resolve_ravel_tmp_parent("/unused") == str(environment)
    assert explicit.is_dir()
    assert environment.is_dir()
