import glob
import json
import os
import shutil
import subprocess
import tempfile

import nibabel as nib
import numpy as np
from joblib import Parallel, delayed


SYNTHSEG_WM_LABELS = (2, 41)
SYNTHSEG_CSF_LABELS = (4, 5, 14, 15, 24, 43, 44)


def _as_float32_image(data, reference_img):
    header = reference_img.header.copy()
    header.set_data_dtype(np.float32)
    return nib.Nifti1Image(data.astype(np.float32), reference_img.affine, header)


def _write_json(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def _coerce_n_jobs(n_jobs, n_items):
    if n_items <= 0:
        return 1
    if n_jobs is None:
        return 1

    n_jobs = int(n_jobs)
    if n_jobs == -1:
        requested = os.cpu_count() or 1
    else:
        requested = max(1, n_jobs)
    return min(requested, n_items)


def _find_raw_anat_files(raw_dir, participant_id, session_id):
    anat_dir = os.path.join(raw_dir, participant_id, session_id, "anat")
    if not os.path.exists(anat_dir):
        raise ValueError(f"Raw anatomy directory missing: {anat_dir}")

    anat_files = glob.glob(os.path.join(anat_dir, "*.nii*"))
    t1w_files = sorted(
        f for f in anat_files
        if (
            ("t1w" in os.path.basename(f).lower()
             and "acq-" not in os.path.basename(f).lower()
             and "desc-" not in os.path.basename(f).lower())
            or "unit1" in os.path.basename(f).lower()
        )
    )
    flair_files = sorted(
        f for f in anat_files
        if (
            "flair" in os.path.basename(f).lower()
            and "acq-" not in os.path.basename(f).lower()
            and "desc-" not in os.path.basename(f).lower()
        )
    )

    if not t1w_files or not flair_files:
        raise ValueError(f"Missing T1w or FLAIR in {anat_dir}")
    return t1w_files[-1], flair_files[-1]


def _reorient_lpi_and_std(in_file, out_file):
    import nibabel.orientations as nio

    img = nib.load(in_file)
    orig_ornt = nio.io_orientation(img.affine)
    targ_ornt = nio.axcodes2ornt(("L", "P", "I"))
    transform = nio.ornt_transform(orig_ornt, targ_ornt)
    img_orient = img.as_reoriented(transform)
    img_can = nib.as_closest_canonical(img_orient)
    nib.save(img_can, out_file)


def _n4_and_resample(in_file, out_file):
    import ants
    import nibabel.processing

    img = ants.image_read(in_file)
    img_n4 = ants.n4_bias_field_correction(img)

    with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
        ants.image_write(img_n4, tmp.name)
        tmp_name = tmp.name

    try:
        nib_img = nib.load(tmp_name)
        resampled_img = nibabel.processing.resample_to_output(
            nib_img,
            voxel_sizes=(1.0, 1.0, 1.0),
        )
        nib.save(resampled_img, out_file)
    finally:
        if os.path.exists(tmp_name):
            os.remove(tmp_name)


def _run_synthseg(in_file, out_file, threads=1):
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return out_file

    subprocess.run(
        [
            "lamareg",
            "synthseg",
            "--i",
            in_file,
            "--o",
            out_file,
            "--robust",
            "--threads",
            str(threads),
            "--cpu",
        ],
        check=True,
    )
    return out_file


def _find_micapipe_flair_transform(micapipe_dir, participant_id, session_id):
    bids_id = f"{participant_id}_{session_id}"
    xfm_dir = os.path.join(micapipe_dir, participant_id, session_id, "xfm")
    if not os.path.isdir(xfm_dir):
        raise ValueError(f"Missing micapipe transform directory: {xfm_dir}")

    affine_path = os.path.join(
        xfm_dir,
        f"{bids_id}_from-flair_to-nativepro_mode-image_desc-affine_0GenericAffine.mat",
    )
    if os.path.exists(affine_path):
        return affine_path

    candidates = sorted(
        glob.glob(
            os.path.join(
                xfm_dir,
                f"{bids_id}_from-flair_to-nativepro*_0GenericAffine.mat",
            )
        )
    )
    if candidates:
        return candidates[-1]

    raise ValueError(f"Missing FLAIR->nativepro transform files in {xfm_dir}")


def _label_mask(label_path, labels, output_path):
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return output_path

    label_img = nib.load(label_path)
    label_data = np.asarray(label_img.dataobj)
    mask = np.isin(np.rint(label_data).astype(np.int16), labels)
    nib.save(_as_float32_image(mask.astype(np.float32), label_img), output_path)
    return output_path


def _estimate_whitestripe(data, labels, wm_labels=SYNTHSEG_WM_LABELS, min_voxels=100, width=0.05):
    label_int = np.rint(labels).astype(np.int16)
    wm_mask = np.isin(label_int, wm_labels)
    wm_values = data[wm_mask & np.isfinite(data) & (data != 0)]

    if wm_values.size < min_voxels:
        raise ValueError(
            f"WhiteStripe requires at least {min_voxels} WM voxels; found {wm_values.size}"
        )

    lo, hi = np.percentile(wm_values, [1, 99])
    wm_values = wm_values[(wm_values >= lo) & (wm_values <= hi)]
    if wm_values.size < min_voxels:
        raise ValueError("WhiteStripe WM estimate is too sparse after outlier trimming")

    hist, edges = np.histogram(wm_values, bins=512)
    mode_idx = int(np.argmax(hist))
    mode = float((edges[mode_idx] + edges[mode_idx + 1]) / 2.0)

    mode_cdf = float(np.mean(wm_values <= mode))
    lower_pct = max(0.0, mode_cdf - width) * 100.0
    upper_pct = min(1.0, mode_cdf + width) * 100.0
    lower, upper = np.percentile(wm_values, [lower_pct, upper_pct])
    stripe = wm_values[(wm_values >= lower) & (wm_values <= upper)]
    if stripe.size < min_voxels:
        q40, q60 = np.percentile(wm_values, [40, 60])
        stripe = wm_values[(wm_values >= q40) & (wm_values <= q60)]
    if stripe.size < min_voxels:
        stripe = wm_values

    robust_sd = float(np.std(wm_values))
    mean = float(np.mean(stripe))
    sd = float(np.std(stripe, ddof=1)) if stripe.size > 1 else robust_sd
    if not np.isfinite(sd) or sd <= 0:
        sd = robust_sd if robust_sd > 0 else 1.0

    return {
        "wm_voxels": int(wm_values.size),
        "stripe_voxels": int(stripe.size),
        "mode": mode,
        "mode_cdf": mode_cdf,
        "width": width,
        "mean": mean,
        "std": sd,
    }


def whitestripe_normalize_image(
    image_path,
    synthseg_label_path,
    output_path,
    stats_path=None,
    wm_labels=SYNTHSEG_WM_LABELS,
):
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return output_path

    img = nib.load(image_path)
    data = np.asarray(img.dataobj, dtype=np.float32)
    labels = np.asarray(nib.load(synthseg_label_path).dataobj)

    if data.shape != labels.shape:
        raise ValueError(
            f"Image and SynthSeg labels have different shapes: {data.shape} vs {labels.shape}"
        )

    stats = _estimate_whitestripe(data, labels, wm_labels=wm_labels)
    brain_mask = np.isfinite(data) & (data != 0)
    normalized = np.zeros(data.shape, dtype=np.float32)
    normalized[brain_mask] = (data[brain_mask] - stats["mean"]) / stats["std"]

    nib.save(_as_float32_image(normalized, img), output_path)
    if stats_path:
        _write_json(stats_path, stats)
    return output_path


def _prepare_from_nativepro(
    participant_id,
    session_id,
    output_dir,
    micapipe_dir,
    threads=1,
):
    bids_id = f"{participant_id}_{session_id}"
    subject_micapipe_dir = os.path.join(micapipe_dir, participant_id, session_id)

    source_paths = {
        "T1w": os.path.join(subject_micapipe_dir, "anat", f"{bids_id}_space-nativepro_T1w.nii.gz"),
        "FLAIR": os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-flair.nii.gz"),
    }
    outputs = {}

    for modality, source_path in source_paths.items():
        if not os.path.exists(source_path):
            raise ValueError(f"Missing {modality} source image: {source_path}")

        label_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_{modality}.nii.gz")
        ws_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_{modality}.nii.gz")
        stats_path = os.path.join(output_dir, f"{bids_id}_desc-whitestripe_{modality}.json")
        wm_mask_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthsegWm_{modality}.nii.gz")
        csf_mask_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthsegCsf_{modality}.nii.gz")

        _run_synthseg(source_path, label_path, threads=threads)
        whitestripe_normalize_image(source_path, label_path, ws_path, stats_path=stats_path)
        _label_mask(label_path, SYNTHSEG_WM_LABELS, wm_mask_path)
        _label_mask(label_path, SYNTHSEG_CSF_LABELS, csf_mask_path)

        outputs[modality] = {
            "source": source_path,
            "whitestripe": ws_path,
            "synthseg": label_path,
            "wm_mask": wm_mask_path,
            "csf_mask": csf_mask_path,
            "stats": stats_path,
        }

    return outputs


def _prepare_from_raw(
    participant_id,
    session_id,
    raw_dir,
    output_dir,
    micapipe_dir,
    tmp_dir,
    threads=1,
    verbose=True,
):
    import ants

    bids_id = f"{participant_id}_{session_id}"
    t1_raw, flair_raw = _find_raw_anat_files(raw_dir, participant_id, session_id)
    os.makedirs(tmp_dir, exist_ok=True)

    preprocessed = {}
    for modality, raw_path in {"T1w": t1_raw, "FLAIR": flair_raw}.items():
        reo_path = os.path.join(tmp_dir, f"{bids_id}_{modality}_reo.nii.gz")
        n4_path = os.path.join(tmp_dir, f"{bids_id}_{modality}_n4.nii.gz")
        label_path = os.path.join(tmp_dir, f"{bids_id}_{modality}_synthseg.nii.gz")
        ws_path = os.path.join(tmp_dir, f"{bids_id}_{modality}_whitestripe.nii.gz")
        stats_path = os.path.join(output_dir, f"{bids_id}_desc-whitestripe_{modality}.json")

        if not os.path.exists(reo_path):
            if verbose:
                print(f"  Reorienting {modality} for WhiteStripe normalization...")
            _reorient_lpi_and_std(raw_path, reo_path)
        if not os.path.exists(n4_path):
            if verbose:
                print(f"  Running N4 bias correction for {modality}...")
            _n4_and_resample(reo_path, n4_path)
        _run_synthseg(n4_path, label_path, threads=threads)
        whitestripe_normalize_image(n4_path, label_path, ws_path, stats_path=stats_path)

        preprocessed[modality] = {
            "raw": raw_path,
            "n4": n4_path,
            "whitestripe": ws_path,
            "synthseg": label_path,
            "stats": stats_path,
        }

    native_ref = os.path.join(
        micapipe_dir,
        participant_id,
        session_id,
        "anat",
        f"{bids_id}_space-nativepro_T1w.nii.gz",
    )
    if not os.path.exists(native_ref):
        raise ValueError(f"Missing nativepro reference T1w image: {native_ref}")

    outputs = {}
    t1_ws_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_T1w.nii.gz")
    t1_label_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_T1w.nii.gz")

    if not os.path.exists(t1_ws_native) or not os.path.exists(t1_label_native):
        if verbose:
            print("  Registering WhiteStripe T1w and SynthSeg labels to nativepro...")
        reg = ants.registration(
            fixed=ants.image_read(native_ref),
            moving=ants.image_read(preprocessed["T1w"]["whitestripe"]),
            type_of_transform="Affine",
        )
        ants.image_write(reg["warpedmovout"], t1_ws_native)
        label_native = ants.apply_transforms(
            fixed=ants.image_read(native_ref),
            moving=ants.image_read(preprocessed["T1w"]["synthseg"]),
            transformlist=reg["fwdtransforms"],
            interpolator="nearestNeighbor",
        )
        ants.image_write(label_native, t1_label_native)

    flair_ws_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_FLAIR.nii.gz")
    flair_label_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_FLAIR.nii.gz")

    if not os.path.exists(flair_ws_native) or not os.path.exists(flair_label_native):
        if verbose:
            print("  Mapping WhiteStripe FLAIR and SynthSeg labels to nativepro...")
        affine_path = _find_micapipe_flair_transform(micapipe_dir, participant_id, session_id)
        fixed = ants.image_read(native_ref)
        flair_native = ants.apply_transforms(
            fixed=fixed,
            moving=ants.image_read(preprocessed["FLAIR"]["whitestripe"]),
            transformlist=[affine_path],
            interpolator="linear",
        )
        ants.image_write(flair_native, flair_ws_native)
        flair_label = ants.apply_transforms(
            fixed=fixed,
            moving=ants.image_read(preprocessed["FLAIR"]["synthseg"]),
            transformlist=[affine_path],
            interpolator="nearestNeighbor",
        )
        ants.image_write(flair_label, flair_label_native)

    for modality, ws_path, label_path in (
        ("T1w", t1_ws_native, t1_label_native),
        ("FLAIR", flair_ws_native, flair_label_native),
    ):
        wm_mask_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthsegWm_{modality}.nii.gz")
        csf_mask_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthsegCsf_{modality}.nii.gz")
        _label_mask(label_path, SYNTHSEG_WM_LABELS, wm_mask_path)
        _label_mask(label_path, SYNTHSEG_CSF_LABELS, csf_mask_path)
        outputs[modality] = {
            "source": preprocessed[modality]["raw"],
            "whitestripe": ws_path,
            "synthseg": label_path,
            "wm_mask": wm_mask_path,
            "csf_mask": csf_mask_path,
            "stats": preprocessed[modality]["stats"],
        }

    return outputs


def prepare_t1w_flair_whitestripe(
    participant_id,
    session_id,
    output_dir,
    micapipe_dir,
    tmp_dir,
    raw_dir=None,
    threads=1,
    verbose=True,
):
    """
    Prepare per-subject T1w and FLAIR WhiteStripe volumes and SynthSeg masks.

    Raw BIDS scans are preferred when ``raw_dir`` is supplied. Otherwise the
    nativepro micapipe T1w and FLAIR maps are normalized in place.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    if raw_dir:
        return _prepare_from_raw(
            participant_id=participant_id,
            session_id=session_id,
            raw_dir=raw_dir,
            output_dir=output_dir,
            micapipe_dir=micapipe_dir,
            tmp_dir=tmp_dir,
            threads=threads,
            verbose=verbose,
        )

    if verbose:
        print("  raw_data_directory not set; WhiteStripe normalizing micapipe nativepro volumes.")
    return _prepare_from_nativepro(
        participant_id=participant_id,
        session_id=session_id,
        output_dir=output_dir,
        micapipe_dir=micapipe_dir,
        threads=threads,
    )


def requested_ravel_modalities(features):
    bases = {str(feature).lower().replace("*blur", "") for feature in features}
    if {"t1w", "flair"} & bases:
        return ["T1w", "FLAIR"]
    return []


def get_normalized_modality_path(subject_output_dir, subject_micapipe_dir, modality):
    bids_id = "_".join(subject_output_dir.rstrip(os.sep).split(os.sep)[-2:])
    maps_dir = os.path.join(subject_output_dir, "maps")

    if modality == "T1w":
        return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-ravel_T1w.nii.gz")
    elif modality == "FLAIR":
        return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-ravel_FLAIR.nii.gz")

    raise ValueError(f"Unsupported normalized modality: {modality}")


def _reference_image_path():
    return os.path.join(os.path.dirname(__file__), "data", "MNI152_T1_2mm.nii.gz")


def _subject_maps_dir(output_directory, participant_id, session_id):
    return os.path.join(output_directory, participant_id, session_id, "maps")


def _ravel_dir(output_directory):
    return os.path.join(output_directory, "ravel")


def ravel_model_path(output_directory, modality):
    return os.path.join(_ravel_dir(output_directory), f"{modality}_ravel_model.npz")


def ravel_models_exist(output_directory, modalities):
    return all(os.path.exists(ravel_model_path(output_directory, modality)) for modality in modalities)


def _safe_copy_transform(path, destination):
    if path and os.path.exists(path) and os.path.abspath(path) != os.path.abspath(destination):
        shutil.copy2(path, destination)
    return destination


def _ensure_mni_inputs(
    participant_id,
    session_id,
    output_directory,
    micapipe_dir,
    modalities,
    verbose=True,
):
    import ants

    bids_id = f"{participant_id}_{session_id}"
    maps_dir = _subject_maps_dir(output_directory, participant_id, session_id)
    ravel_subject_dir = os.path.join(maps_dir, "ravel")
    os.makedirs(ravel_subject_dir, exist_ok=True)

    reference_path = _reference_image_path()
    reference_img = ants.image_read(reference_path)
    native_ref = os.path.join(
        micapipe_dir,
        participant_id,
        session_id,
        "anat",
        f"{bids_id}_space-nativepro_T1w.nii.gz",
    )
    if not os.path.exists(native_ref):
        raise ValueError(f"Missing nativepro T1w reference: {native_ref}")

    transform_json = os.path.join(ravel_subject_dir, f"{bids_id}_ravel_transforms.json")
    native_to_mni = []
    mni_to_native = []

    if os.path.exists(transform_json):
        with open(transform_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        native_to_mni = payload.get("native_to_mni", [])
        mni_to_native = payload.get("mni_to_native", [])

    t1_mni = os.path.join(ravel_subject_dir, f"{bids_id}_space-MNI152_desc-whitestripe_T1w.nii.gz")
    if not native_to_mni or not all(os.path.exists(p) for p in native_to_mni) or not os.path.exists(t1_mni):
        if verbose:
            print(f"  Registering {participant_id}/{session_id} WhiteStripe T1w to MNI for RAVEL...")
        t1_native = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-whitestripe_T1w.nii.gz")
        if not os.path.exists(t1_native):
            raise ValueError(f"Missing WhiteStripe T1w before RAVEL: {t1_native}")

        outprefix = os.path.join(ravel_subject_dir, f"{bids_id}_from-nativepro_to-MNI_")
        reg = ants.registration(
            fixed=reference_img,
            moving=ants.image_read(t1_native),
            type_of_transform="Affine",
            outprefix=outprefix,
        )
        ants.image_write(reg["warpedmovout"], t1_mni)

        native_to_mni = [
            _safe_copy_transform(path, os.path.join(ravel_subject_dir, os.path.basename(path)))
            for path in reg["fwdtransforms"]
        ]
        mni_to_native = [
            _safe_copy_transform(path, os.path.join(ravel_subject_dir, os.path.basename(path)))
            for path in reg["invtransforms"]
        ]
        _write_json(
            transform_json,
            {
                "native_to_mni": native_to_mni,
                "mni_to_native": mni_to_native,
                "reference": reference_path,
            },
        )

    outputs = {
        "reference": reference_path,
        "native_reference": native_ref,
        "native_to_mni": native_to_mni,
        "mni_to_native": mni_to_native,
        "modalities": {},
    }

    fixed = ants.image_read(reference_path)
    for modality in modalities:
        native_ws = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-whitestripe_{modality}.nii.gz")
        native_csf = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-synthsegCsf_{modality}.nii.gz")
        mni_ws = os.path.join(ravel_subject_dir, f"{bids_id}_space-MNI152_desc-whitestripe_{modality}.nii.gz")
        mni_csf = os.path.join(ravel_subject_dir, f"{bids_id}_space-MNI152_desc-synthsegCsf_{modality}.nii.gz")
        native_ravel = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-ravel_{modality}.nii.gz")
        mni_ravel = os.path.join(ravel_subject_dir, f"{bids_id}_space-MNI152_desc-ravel_{modality}.nii.gz")

        if not os.path.exists(native_ws):
            raise ValueError(f"Missing WhiteStripe {modality} before RAVEL: {native_ws}")
        if not os.path.exists(native_csf):
            raise ValueError(f"Missing SynthSeg CSF mask for RAVEL: {native_csf}")

        if not os.path.exists(mni_ws):
            transformed = ants.apply_transforms(
                fixed=fixed,
                moving=ants.image_read(native_ws),
                transformlist=native_to_mni,
                interpolator="linear",
            )
            ants.image_write(transformed, mni_ws)

        if not os.path.exists(mni_csf):
            transformed = ants.apply_transforms(
                fixed=fixed,
                moving=ants.image_read(native_csf),
                transformlist=native_to_mni,
                interpolator="nearestNeighbor",
            )
            ants.image_write(transformed, mni_csf)

        outputs["modalities"][modality] = {
            "native_whitestripe": native_ws,
            "native_csf": native_csf,
            "mni_whitestripe": mni_ws,
            "mni_csf": mni_csf,
            "mni_ravel": mni_ravel,
            "native_ravel": native_ravel,
        }

    return outputs


def _load_brain_vector(path, brain_mask):
    data = np.asarray(nib.load(path).dataobj, dtype=np.float32)
    if data.shape != brain_mask.shape:
        raise ValueError(f"RAVEL image shape mismatch for {path}: {data.shape} vs {brain_mask.shape}")
    return data[brain_mask].astype(np.float32)


def _write_brain_vector(vector, brain_mask, reference_img, output_path):
    data = np.zeros(brain_mask.shape, dtype=np.float32)
    data[brain_mask] = vector.astype(np.float32)
    nib.save(_as_float32_image(data, reference_img), output_path)
    return output_path


def _control_mask_from_csf(csf_paths, brain_mask, min_voxels=100):
    csf_vectors = []
    for path in csf_paths:
        csf = np.asarray(nib.load(path).dataobj) > 0.5
        if csf.shape != brain_mask.shape:
            raise ValueError(f"RAVEL CSF mask shape mismatch for {path}: {csf.shape} vs {brain_mask.shape}")
        csf_vectors.append(csf[brain_mask])

    csf_matrix = np.vstack(csf_vectors)
    control_mask = np.all(csf_matrix, axis=0)
    if int(control_mask.sum()) >= min_voxels:
        return control_mask, "intersection"

    threshold = max(1, int(np.ceil(0.5 * len(csf_vectors))))
    control_mask = np.sum(csf_matrix, axis=0) >= threshold
    if int(control_mask.sum()) >= min_voxels:
        return control_mask, f"present_in_at_least_{threshold}_subjects"

    control_mask = np.any(csf_matrix, axis=0)
    if int(control_mask.sum()) >= min_voxels:
        return control_mask, "union"

    raise ValueError("RAVEL CSF control mask has too few voxels")


def _inverse_transform_ravel(mni_path, native_path, native_reference, mni_to_native):
    import ants

    corrected = ants.apply_transforms(
        fixed=ants.image_read(native_reference),
        moving=ants.image_read(mni_path),
        transformlist=mni_to_native,
        interpolator="linear",
    )
    ants.image_write(corrected, native_path)
    return native_path


def fit_and_apply_ravel_to_controls(
    subjects,
    output_directory,
    micapipe_dir,
    modalities,
    k=1,
    verbose=True,
    n_jobs=1,
):
    """
    Fit RAVEL on the control cohort and write RAVEL-corrected nativepro volumes.

    This follows the core RAVEL correction used by jfortin1/RAVEL: build a
    voxel-by-subject matrix, estimate unwanted factors from CSF control voxels
    with SVD, regress those factors out, then add back the voxelwise cohort mean.
    """
    if not subjects:
        raise ValueError("Cannot fit RAVEL without control subjects")

    os.makedirs(_ravel_dir(output_directory), exist_ok=True)
    modalities = list(modalities)

    def ensure_subject_mni_inputs(subject):
        return subject, _ensure_mni_inputs(
            subject[0],
            subject[1],
            output_directory,
            micapipe_dir,
            modalities,
            verbose=verbose,
        )

    ravel_n_jobs = _coerce_n_jobs(n_jobs, len(subjects))
    if verbose and ravel_n_jobs > 1:
        print(f"Preparing control RAVEL MNI inputs with {ravel_n_jobs} parallel jobs...")

    if ravel_n_jobs == 1:
        mni_inputs = dict(ensure_subject_mni_inputs(subject) for subject in subjects)
    else:
        mni_inputs = dict(
            Parallel(n_jobs=ravel_n_jobs)(
                delayed(ensure_subject_mni_inputs)(subject)
                for subject in subjects
            )
        )

    reference_path = _reference_image_path()
    reference_img = nib.load(reference_path)
    brain_data = np.asarray(reference_img.dataobj)
    brain_mask = np.isfinite(brain_data) & (brain_data > 0)

    for modality in modalities:
        if verbose:
            print(f"Fitting RAVEL model for {modality} on {len(subjects)} control subjects...")

        image_paths = [
            mni_inputs[subject]["modalities"][modality]["mni_whitestripe"]
            for subject in subjects
        ]
        csf_paths = [
            mni_inputs[subject]["modalities"][modality]["mni_csf"]
            for subject in subjects
        ]
        V = np.column_stack([_load_brain_vector(path, brain_mask) for path in image_paths])
        control_mask, control_strategy = _control_mask_from_csf(csf_paths, brain_mask)
        mean = V.mean(axis=1).astype(np.float32)
        V_centered = (V - mean[:, None]).astype(np.float32)
        Vc = V_centered[control_mask, :]

        U, singular_values, vh = np.linalg.svd(Vc, full_matrices=False)
        k_eff = min(int(k), vh.shape[0], max(1, len(subjects)))
        Z = vh.T[:, :k_eff].astype(np.float32)
        gamma = (np.linalg.pinv(Z) @ V_centered.T).astype(np.float32)
        fitted = (Z @ gamma).T.astype(np.float32)
        V_norm = (V - fitted).astype(np.float32)

        model_path = ravel_model_path(output_directory, modality)
        np.savez_compressed(
            model_path,
            brain_mask=brain_mask,
            control_mask=control_mask,
            U=U[:, :k_eff].astype(np.float32),
            singular_values=singular_values[:k_eff].astype(np.float32),
            Z=Z,
            gamma=gamma,
            mean=mean,
            subjects=np.array([f"{pid}_{sid}" for pid, sid in subjects]),
        )
        _write_json(
            model_path.replace(".npz", ".json"),
            {
                "modality": modality,
                "k": k_eff,
                "reference": reference_path,
                "n_subjects": len(subjects),
                "n_brain_voxels": int(brain_mask.sum()),
                "n_control_voxels": int(control_mask.sum()),
                "control_mask_strategy": control_strategy,
                "subjects": [f"{pid}_{sid}" for pid, sid in subjects],
            },
        )

        for column, subject in enumerate(subjects):
            paths = mni_inputs[subject]["modalities"][modality]
            _write_brain_vector(V_norm[:, column], brain_mask, reference_img, paths["mni_ravel"])
            _inverse_transform_ravel(
                paths["mni_ravel"],
                paths["native_ravel"],
                mni_inputs[subject]["native_reference"],
                mni_inputs[subject]["mni_to_native"],
            )


def apply_ravel_model_to_subject(
    participant_id,
    session_id,
    output_directory,
    micapipe_dir,
    modalities,
    verbose=True,
):
    os.makedirs(_ravel_dir(output_directory), exist_ok=True)
    modalities = list(modalities)
    mni_inputs = _ensure_mni_inputs(
        participant_id,
        session_id,
        output_directory,
        micapipe_dir,
        modalities,
        verbose=verbose,
    )

    reference_img = nib.load(_reference_image_path())

    for modality in modalities:
        model_path = ravel_model_path(output_directory, modality)
        if not os.path.exists(model_path):
            raise ValueError(
                f"Missing RAVEL model for {modality}: {model_path}. "
                "Process the control dataset first."
            )

        model = np.load(model_path, allow_pickle=True)
        brain_mask = model["brain_mask"].astype(bool)
        control_mask = model["control_mask"].astype(bool)
        U = model["U"].astype(np.float32)
        singular_values = model["singular_values"].astype(np.float32)
        gamma = model["gamma"].astype(np.float32)
        mean = model["mean"].astype(np.float32)

        paths = mni_inputs["modalities"][modality]
        V = _load_brain_vector(paths["mni_whitestripe"], brain_mask)
        V_centered = (V - mean).astype(np.float32)
        Vc = V_centered[control_mask]
        denom = np.where(singular_values > np.finfo(np.float32).eps, singular_values, np.inf)
        z = ((Vc @ U) / denom).astype(np.float32)
        corrected = (V - (z @ gamma)).astype(np.float32)

        if verbose:
            print(f"  Applying fitted RAVEL model to {participant_id}/{session_id} {modality}...")

        _write_brain_vector(corrected, brain_mask, reference_img, paths["mni_ravel"])
        _inverse_transform_ravel(
            paths["mni_ravel"],
            paths["native_ravel"],
            mni_inputs["native_reference"],
            mni_inputs["mni_to_native"],
        )
