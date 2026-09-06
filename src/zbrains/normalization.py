import glob
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import threading
from contextlib import contextmanager

import fcntl

import nibabel as nib
import numpy as np
from joblib import Parallel, delayed


SYNTHSEG_WM_LABELS = (2, 41)
SYNTHSEG_CSF_LABELS = (4, 5, 14, 15, 24, 43, 44)
RAVEL_CHUNK_SIZE = 50000
_ANTS_LOCK = threading.Lock()


def _configure_native_threads(threads):
    threads = max(1, int(threads or 1))
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(threads)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = str(threads)
    return threads


def _as_float32_image(data, reference_img):
    header = reference_img.header.copy()
    header.set_data_dtype(np.float32)
    return nib.Nifti1Image(data.astype(np.float32), reference_img.affine, header)


def _write_json(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def _write_json_atomic(path, payload):
    """Write JSON without exposing a partially-written cache sidecar."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = (
        f"{path}.tmp-{os.getpid()}-{threading.get_ident()}"
    )
    try:
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


@contextmanager
def _exclusive_file_lock(lock_path):
    """Serialize writers to a cache entry across threads and Linux processes."""
    os.makedirs(os.path.dirname(lock_path), exist_ok=True)
    with open(lock_path, "a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _sha256(path, block_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_identity(path, include_hash=False):
    stat = os.stat(path)
    identity = {
        "path": os.path.realpath(path),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if include_hash:
        identity["sha256"] = _sha256(path)
    return identity


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


def _find_raw_anat_files(raw_dir, participant_id, session_id, modalities):
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

    paths = {}
    if "T1w" in modalities:
        if not t1w_files:
            raise ValueError(f"Missing T1w in {anat_dir}")
        paths["T1w"] = t1w_files[-1]
    if "FLAIR" in modalities:
        if not flair_files:
            raise ValueError(f"Missing FLAIR in {anat_dir}")
        paths["FLAIR"] = flair_files[-1]
    return paths


def _n4_and_resample(in_file, out_file, threads=1):
    import nibabel.processing

    _configure_native_threads(threads)

    with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
        tmp_name = tmp.name
    with _ANTS_LOCK:
        import ants

        img = ants.image_read(in_file)
        img_n4 = ants.n4_bias_field_correction(img)
        ants.image_write(img_n4, tmp_name)

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


def _synthseg_cache_path(out_file):
    """Base-independent shared cache path for a SynthSeg label. The per-base label
    path contains a ``zbrains_<base>`` component (one per normalization/smoothing/
    exclusion arm); swapping it for a fixed ``synthseg_cache`` yields a path keyed
    only by subject + modality, so the SAME raw T1w/FLAIR is segmented ONCE and
    reused across every base. Returns None if no ``zbrains_*`` component is found
    (then SynthSeg runs directly, preserving legacy behaviour)."""
    parts = os.path.abspath(out_file).split(os.sep)
    for i, part in enumerate(parts):
        if part.startswith("zbrains_") or part.startswith("zbrains-"):
            parts[i] = "synthseg_cache"
            return os.sep.join(parts)
    return None


def _symlink_label(link_path, target_path):
    """Point ``link_path`` at ``target_path`` (absolute symlink), replacing any
    stale link/file. No-op if the link already resolves to the target."""
    link_path = os.path.abspath(link_path)
    target_path = os.path.abspath(target_path)
    os.makedirs(os.path.dirname(link_path), exist_ok=True)
    if os.path.islink(link_path) or os.path.exists(link_path):
        try:
            if os.path.realpath(link_path) == target_path:
                return link_path
            os.remove(link_path)
        except OSError:
            pass
    try:
        os.symlink(target_path, link_path)
    except FileExistsError:
        pass
    return link_path


def _synthseg_subprocess(in_file, out_file, threads=1):
    # Isolate the SynthSeg (lamareg) subprocess from the parent's environment:
    # a relative PYTHONPATH (e.g. "src") or a deleted parent CWD makes the child
    # Python fail at startup ("init_fs_encoding: failed to get the Python codec").
    # SynthSeg has its own package install and needs neither our PYTHONPATH nor
    # our CWD, so strip PYTHONPATH, pin a UTF-8 locale, and run from a stable dir.
    out_dir = os.path.dirname(os.path.abspath(out_file))
    os.makedirs(out_dir, exist_ok=True)
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.setdefault("LC_ALL", "C.UTF-8")
    env.setdefault("LANG", "C.UTF-8")

    subprocess.run(
        [
            "lamareg",
            "synthseg",
            "--i",
            os.path.abspath(in_file),
            "--o",
            os.path.abspath(out_file),
            "--robust",
            "--threads",
            str(threads),
            "--cpu",
        ],
        check=True,
        env=env,
        cwd=out_dir,
    )
    return out_file


def _run_synthseg(in_file, out_file, threads=1):
    """Segment ``in_file`` -> ``out_file``, reusing a base-independent per-subject
    cache. SynthSeg depends only on the raw nativepro T1w/FLAIR (identical across
    all normalization/smoothing/exclusion bases), so it is computed ONCE into a
    shared ``synthseg_cache`` and symlinked into each base -- eliminating the ~10x
    redundant re-segmentation of the same image across greedy arms."""
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return out_file
    cache = _synthseg_cache_path(out_file)
    if cache is None:
        return _synthseg_subprocess(in_file, out_file, threads=threads)
    if not (os.path.exists(cache) and os.path.getsize(cache) > 0):
        _synthseg_subprocess(in_file, cache, threads=threads)
    return _symlink_label(out_file, cache)


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


def _find_micapipe_nativepro_mni_transforms(micapipe_dir, participant_id, session_id):
    bids_id = f"{participant_id}_{session_id}"
    xfm_dir = os.path.join(micapipe_dir, participant_id, session_id, "xfm")
    if not os.path.isdir(xfm_dir):
        raise ValueError(f"Missing micapipe transform directory: {xfm_dir}")

    affine_path = os.path.join(
        xfm_dir,
        f"{bids_id}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
    )
    warp_path = os.path.join(
        xfm_dir,
        f"{bids_id}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
    )
    inverse_warp_path = os.path.join(
        xfm_dir,
        f"{bids_id}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1InverseWarp.nii.gz",
    )

    if not os.path.exists(affine_path):
        candidates = sorted(
            glob.glob(os.path.join(xfm_dir, f"{bids_id}_from-nativepro*_to-MNI152*_0GenericAffine.mat"))
        )
        if candidates:
            affine_path = candidates[-1]
    if not os.path.exists(warp_path):
        candidates = sorted(
            glob.glob(os.path.join(xfm_dir, f"{bids_id}_from-nativepro*_to-MNI152*_1Warp.nii.gz"))
        )
        if candidates:
            warp_path = candidates[-1]
    if not os.path.exists(inverse_warp_path):
        candidates = sorted(
            glob.glob(os.path.join(xfm_dir, f"{bids_id}_from-nativepro*_to-MNI152*_1InverseWarp.nii.gz"))
        )
        if candidates:
            inverse_warp_path = candidates[-1]

    missing = [
        path for path in (affine_path, warp_path, inverse_warp_path)
        if not os.path.exists(path)
    ]
    if missing:
        raise ValueError(
            "Missing micapipe nativepro<->MNI152 SyN transform file(s): "
            + ", ".join(missing)
        )

    return {
        "native_to_mni": [warp_path, affine_path],
        "native_to_mni_invert": [False, False],
        "mni_to_native": [affine_path, inverse_warp_path],
        "mni_to_native_invert": [True, False],
    }


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

    q25, q75 = np.percentile(wm_values, [25, 75])
    iqr_sd = float((q75 - q25) / 1.349) if q75 > q25 else 0.0
    trimmed_sd = float(np.std(wm_values, ddof=1)) if wm_values.size > 1 else 0.0
    mean = float(np.mean(stripe))
    sd = iqr_sd if np.isfinite(iqr_sd) and iqr_sd > 0 else trimmed_sd
    if not np.isfinite(sd) or sd <= 0:
        sd = 1.0

    return {
        "wm_voxels": int(wm_values.size),
        "stripe_voxels": int(stripe.size),
        "mode": mode,
        "mode_cdf": mode_cdf,
        "width": width,
        "mean": mean,
        "std": sd,
        "scale_method": "trimmed_wm_iqr",
        "trimmed_wm_iqr_sd": iqr_sd,
        "trimmed_wm_std": trimmed_sd,
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
    label_img = nib.load(synthseg_label_path)
    labels = np.asarray(label_img.dataobj)

    if data.shape != labels.shape:
        # SynthSeg (esp. --robust) conforms its output to ~1mm, so the label
        # volume can land on a different grid than the source image even though
        # both share world coordinates. Resample the segmentation onto the image
        # grid (nearest-neighbor) so the WM mask aligns with the intensities.
        import nibabel.processing

        resampled = nibabel.processing.resample_from_to(
            label_img, (img.shape, img.affine), order=0
        )
        labels = np.asarray(resampled.dataobj)

    stats = _estimate_whitestripe(data, labels, wm_labels=wm_labels)
    brain_mask = np.isfinite(data) & (data != 0)
    normalized = np.zeros(data.shape, dtype=np.float32)
    normalized[brain_mask] = (data[brain_mask] - stats["mean"]) / stats["std"]

    nib.save(_as_float32_image(normalized, img), output_path)
    if stats_path:
        _write_json(stats_path, stats)
    return output_path


def _labels_on_image_grid(label_img, img):
    """Return SynthSeg labels resampled (nearest-neighbor) onto the image grid."""
    labels = np.asarray(label_img.dataobj)
    if labels.shape != img.shape:
        import nibabel.processing

        labels = np.asarray(
            nibabel.processing.resample_from_to(
                label_img, (img.shape, img.affine), order=0
            ).dataobj
        )
    return labels


def _estimate_wm_mean(data, labels, wm_labels=SYNTHSEG_WM_LABELS, min_voxels=100):
    """Robust mean/SD of the WHOLE white-matter tissue (FCM/tissue-mean spirit).

    Where WhiteStripe centers on the mode of a narrow WM *stripe*, this uses the
    trimmed mean/SD of the entire SynthSeg WM mask -- a full-tissue WM landmark.
    """
    label_int = np.rint(labels).astype(np.int16)
    wm_values = data[np.isin(label_int, wm_labels) & np.isfinite(data) & (data != 0)]
    if wm_values.size < min_voxels:
        raise ValueError(f"WM-mean requires >= {min_voxels} WM voxels; found {wm_values.size}")
    lo, hi = np.percentile(wm_values, [1, 99])
    trimmed = wm_values[(wm_values >= lo) & (wm_values <= hi)]
    if trimmed.size < min_voxels:
        trimmed = wm_values
    mean = float(np.mean(trimmed))
    sd = float(np.std(trimmed, ddof=1)) if trimmed.size > 1 else 1.0
    if not np.isfinite(sd) or sd <= 0:
        sd = 1.0
    return {"wm_voxels": int(wm_values.size), "mean": mean, "std": sd,
            "scale_method": "trimmed_whole_wm_mean_sd"}


def wm_mean_normalize_image(image_path, synthseg_label_path, output_path,
                            stats_path=None, wm_labels=SYNTHSEG_WM_LABELS):
    """WM-mean normalization: (x - WM_mean) / WM_sd over the whole WM tissue."""
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return output_path
    img = nib.load(image_path)
    data = np.asarray(img.dataobj, dtype=np.float32)
    labels = _labels_on_image_grid(nib.load(synthseg_label_path), img)
    stats = _estimate_wm_mean(data, labels, wm_labels=wm_labels)
    brain_mask = np.isfinite(data) & (data != 0)
    normalized = np.zeros(data.shape, dtype=np.float32)
    normalized[brain_mask] = (data[brain_mask] - stats["mean"]) / stats["std"]
    nib.save(_as_float32_image(normalized, img), output_path)
    if stats_path:
        _write_json(stats_path, stats)
    return output_path


# --------------------------------------------------------------------------- #
# Nyul-Udupa piecewise-linear histogram standardization (cohort-learned scale) #
# --------------------------------------------------------------------------- #
NYUL_PERCENTILES = (1.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 99.0)
NYUL_OUTPUT_RANGE = (0.0, 100.0)


def _synthseg_brain_mask(label_path):
    """Boolean brain mask from a SynthSeg label volume (any labelled voxel)."""
    return np.asarray(nib.load(label_path).dataobj) > 0


def _synthseg_brain_mask_on_grid(label_path, img):
    """Boolean brain mask (any SynthSeg label > 0) resampled (nearest-neighbor)
    onto ``img``'s grid. SynthSeg (esp. --robust) conforms its output to ~1mm, so
    the label volume can land on a DIFFERENT grid than the subject image (they
    share world coordinates but not shape); aligning it here prevents the
    ``mask & data`` shape mismatch that otherwise crashes the Nyul fit."""
    return _labels_on_image_grid(nib.load(label_path), img) > 0


def _nyul_landmarks(data, percentiles=NYUL_PERCENTILES, brain_mask=None):
    """Landmark percentiles over the brain. ``brain_mask`` (e.g. from SynthSeg) is
    REQUIRED to compose Nyul on subject-normalized (near-zero-centred) images,
    where the legacy ``data != 0`` heuristic would drop most of the brain."""
    if brain_mask is not None:
        values = data[brain_mask & np.isfinite(data)]
    else:
        values = data[np.isfinite(data) & (data != 0)]
    if values.size == 0:
        return None
    return np.percentile(values, percentiles).astype(float)


def _nyul_mapped_landmarks(path, mask_path, percentiles, output_range):
    """Per-image mapped landmark vector for the Nyul fit (or None if unusable).
    Factored out so the cohort fit can compute images concurrently."""
    img = nib.load(path)
    data = np.asarray(img.dataobj, dtype=np.float32)
    brain_mask = _synthseg_brain_mask_on_grid(mask_path, img) if mask_path is not None else None
    landmarks = _nyul_landmarks(data, percentiles, brain_mask=brain_mask)
    if landmarks is None:
        return None
    lo, hi = landmarks[0], landmarks[-1]
    if not np.isfinite(hi - lo) or hi <= lo:
        return None
    return (landmarks - lo) / (hi - lo) * (output_range[1] - output_range[0]) + output_range[0]


def fit_nyul_standard_scale(image_paths, percentiles=NYUL_PERCENTILES,
                            output_range=NYUL_OUTPUT_RANGE, mask_paths=None, n_jobs=1):
    """Learn the standard histogram scale (mean landmark positions) from a cohort.

    Each image's landmark percentiles are linearly mapped so its first/last
    landmarks hit ``output_range``; the standard scale is the mean of those
    mapped landmark vectors across images (classic Nyul-Udupa training).
    ``mask_paths`` (parallel to ``image_paths``) supplies SynthSeg brain masks so
    the fit composes correctly on subject-normalized images. ``n_jobs`` > 1
    computes the per-image landmarks concurrently -- THREADING (not loky): the
    per-image work (gzip image load, label resample, np.percentile) releases the
    GIL, and threads avoid the fork + numpy re-init crash that loky triggers here.
    """
    n = len(image_paths)
    masks = mask_paths if mask_paths is not None else [None] * n
    jobs = max(1, min(int(n_jobs or 1), n)) if n else 1
    if jobs == 1:
        mapped_all = [
            _nyul_mapped_landmarks(image_paths[i], masks[i], percentiles, output_range)
            for i in range(n)
        ]
    else:
        mapped_all = Parallel(n_jobs=jobs, backend="threading")(
            delayed(_nyul_mapped_landmarks)(image_paths[i], masks[i], percentiles, output_range)
            for i in range(n)
        )
    mapped = [m for m in mapped_all if m is not None]
    if not mapped:
        raise ValueError("Nyul fit failed: no usable control images")
    return {
        "method": "nyul_udupa_piecewise_linear",
        "percentiles": list(percentiles),
        "output_range": list(output_range),
        "standard_landmarks": np.mean(mapped, axis=0).tolist(),
        "n_training_images": len(mapped),
    }


def apply_nyul_to_image(image_path, model, output_path, mask_path=None):
    """Piecewise-linearly map an image onto the learned standard histogram scale.

    ``mask_path`` (a SynthSeg label volume) defines the brain; required when the
    input is a subject-normalized (near-zero-centred) image so the ``data != 0``
    heuristic does not drop most of the brain."""
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        return output_path
    img = nib.load(image_path)
    data = np.asarray(img.dataobj, dtype=np.float32)
    percentiles = model["percentiles"]
    standard = np.asarray(model["standard_landmarks"], dtype=float)
    if mask_path is not None:
        brain_mask = _synthseg_brain_mask_on_grid(mask_path, img) & np.isfinite(data)
    else:
        brain_mask = np.isfinite(data) & (data != 0)
    landmarks = _nyul_landmarks(data, percentiles, brain_mask=brain_mask)
    normalized = np.zeros(data.shape, dtype=np.float32)
    if landmarks is not None:
        # strictly increasing landmarks required for interpolation
        landmarks = np.maximum.accumulate(landmarks)
        landmarks[1:] += np.arange(1, landmarks.size) * np.finfo(np.float32).eps
        normalized[brain_mask] = np.interp(data[brain_mask], landmarks, standard).astype(np.float32)
    nib.save(_as_float32_image(normalized, img), output_path)
    return output_path


def _nyul_dir(output_directory):
    return os.path.join(output_directory, "nyul")


def nyul_model_path(output_directory, modality):
    return os.path.join(_nyul_dir(output_directory), f"{modality}_nyul_model.json")


def nyul_models_exist(output_directory, modalities):
    for modality in modalities:
        if not os.path.exists(nyul_model_path(output_directory, modality)):
            return False
    return True


def _nativepro_source_paths(participant_id, session_id, micapipe_dir):
    bids_id = f"{participant_id}_{session_id}"
    subject_dir = os.path.join(micapipe_dir, participant_id, session_id)
    return {
        "T1w": os.path.join(subject_dir, "anat", f"{bids_id}_space-nativepro_T1w.nii.gz"),
        "FLAIR": os.path.join(subject_dir, "maps", f"{bids_id}_space-nativepro_map-flair.nii.gz"),
    }


def _nyul_output_path(output_directory, participant_id, session_id, modality, output_desc="nyul"):
    bids_id = f"{participant_id}_{session_id}"
    maps_dir = os.path.join(output_directory, participant_id, session_id, "maps")
    os.makedirs(maps_dir, exist_ok=True)
    return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{output_desc}_{modality}.nii.gz")


def _nyul_source_path(participant_id, session_id, output_directory, micapipe_dir, modality,
                      input_desc=None):
    """Input image for Nyul. ``input_desc=None`` reads the raw nativepro volume
    (legacy single-mode). A subject-norm label reads the desc-<input_desc>
    intermediate written by the subject-level normalization phase."""
    if input_desc in (None, "raw"):
        return _nativepro_source_paths(participant_id, session_id, micapipe_dir).get(modality)
    bids_id = f"{participant_id}_{session_id}"
    maps_dir = _subject_maps_dir(output_directory, participant_id, session_id)
    return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{input_desc}_{modality}.nii.gz")


def _synthseg_label_for(participant_id, session_id, output_directory, modality):
    """SynthSeg label written by the subject-level phase, or None (legacy)."""
    bids_id = f"{participant_id}_{session_id}"
    maps_dir = _subject_maps_dir(output_directory, participant_id, session_id)
    path = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-synthseg_{modality}.nii.gz")
    return path if os.path.exists(path) else None


def fit_and_apply_nyul_to_controls(subjects, output_directory, micapipe_dir, modalities,
                                   verbose=True, input_desc=None, output_desc="nyul", n_jobs=1):
    """Fit the Nyul standard scale on controls, then normalize every control.

    Reads the desc-<input_desc> subject-norm intermediate (or raw nativepro when
    ``input_desc`` is None) and, when a SynthSeg label is available, restricts
    landmark estimation to the labelled brain so the fit composes on
    subject-normalized (near-zero-centred) images. Writes desc-<output_desc>.
    ``n_jobs`` > 1 parallelizes BOTH the fit (per-image landmarks) and the apply
    (per-image piecewise mapping) with THREADING -- the heavy per-image work
    (gzip load, label resample, np.interp, gzip save) releases the GIL.
    """
    os.makedirs(_nyul_dir(output_directory), exist_ok=True)
    for modality in modalities:
        source_by_subject = {}
        mask_by_subject = {}
        for pid, sid in subjects:
            src = _nyul_source_path(pid, sid, output_directory, micapipe_dir, modality, input_desc)
            if src and os.path.exists(src):
                source_by_subject[(pid, sid)] = src
                mask_by_subject[(pid, sid)] = _synthseg_label_for(pid, sid, output_directory, modality)
        if not source_by_subject:
            if verbose:
                print(f"  Skipping Nyul fit for {modality}; no control sources found.")
            continue
        if verbose:
            print(f"  Fitting Nyul standard scale for {modality} on "
                  f"{len(source_by_subject)} controls (input desc-{input_desc or 'raw'})...")
        keys = list(source_by_subject.keys())
        model = fit_nyul_standard_scale(
            [source_by_subject[k] for k in keys],
            mask_paths=[mask_by_subject[k] for k in keys],
            n_jobs=n_jobs,
        )
        _write_json(nyul_model_path(output_directory, modality), model)

        def _apply_one(key):
            apply_nyul_to_image(
                source_by_subject[key], model,
                _nyul_output_path(output_directory, key[0], key[1], modality, output_desc),
                mask_path=mask_by_subject[key])

        apply_jobs = max(1, min(int(n_jobs or 1), len(keys)))
        if apply_jobs == 1:
            for key in keys:
                _apply_one(key)
        else:
            Parallel(n_jobs=apply_jobs, backend="threading")(
                delayed(_apply_one)(key) for key in keys
            )


def apply_nyul_model_to_subject(participant_id, session_id, output_directory, micapipe_dir,
                                modalities, verbose=True, input_desc=None, output_desc="nyul"):
    """Apply the fitted Nyul standard scale to one subject (e.g. a patient)."""
    for modality in modalities:
        model_file = nyul_model_path(output_directory, modality)
        if not os.path.exists(model_file):
            raise ValueError(f"Missing Nyul model for {modality}: {model_file}. "
                             "Process the control cohort first.")
        with open(model_file, "r", encoding="utf-8") as stream:
            model = json.load(stream)
        src = _nyul_source_path(participant_id, session_id, output_directory, micapipe_dir,
                                modality, input_desc)
        if not src or not os.path.exists(src):
            raise ValueError(f"Missing {modality} source for Nyul: {src}")
        mask = _synthseg_label_for(participant_id, session_id, output_directory, modality)
        apply_nyul_to_image(src, model,
                            _nyul_output_path(output_directory, participant_id, session_id, modality, output_desc),
                            mask_path=mask)


def _prepare_from_nativepro(
    participant_id,
    session_id,
    output_dir,
    micapipe_dir,
    modalities,
    threads=1,
):
    bids_id = f"{participant_id}_{session_id}"
    subject_micapipe_dir = os.path.join(micapipe_dir, participant_id, session_id)

    all_source_paths = {
        "T1w": os.path.join(subject_micapipe_dir, "anat", f"{bids_id}_space-nativepro_T1w.nii.gz"),
        "FLAIR": os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-flair.nii.gz"),
    }
    source_paths = {
        modality: all_source_paths[modality]
        for modality in modalities
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
    modalities,
    threads=1,
    verbose=True,
):
    _configure_native_threads(threads)
    import ants

    bids_id = f"{participant_id}_{session_id}"
    raw_paths = _find_raw_anat_files(raw_dir, participant_id, session_id, modalities)
    os.makedirs(tmp_dir, exist_ok=True)

    preprocessed = {}
    for modality, raw_path in raw_paths.items():
        native_ws_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_{modality}.nii.gz")
        native_label_path = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_{modality}.nii.gz")
        stats_path = os.path.join(output_dir, f"{bids_id}_desc-whitestripe_{modality}.json")

        if os.path.exists(native_ws_path) and os.path.exists(native_label_path):
            if verbose:
                print(f"  Existing nativepro WhiteStripe/SynthSeg outputs found for {modality}; skipping raw preprocessing.")
            preprocessed[modality] = {
                "raw": raw_path,
                "n4": None,
                "whitestripe": native_ws_path,
                "synthseg": native_label_path,
                "stats": stats_path,
            }
            continue

        n4_path = os.path.join(tmp_dir, f"{bids_id}_{modality}_n4.nii.gz")
        label_path = os.path.join(output_dir, f"{bids_id}_desc-preprocSynthseg_{modality}.nii.gz")
        ws_path = os.path.join(output_dir, f"{bids_id}_desc-preprocWhitestripe_{modality}.nii.gz")

        if not os.path.exists(n4_path):
            if verbose:
                print(f"  Running N4 bias correction for {modality}...")
            _n4_and_resample(raw_path, n4_path, threads=threads)
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
    native_outputs = {}

    if "T1w" in modalities:
        t1_ws_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_T1w.nii.gz")
        t1_label_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_T1w.nii.gz")

        if not os.path.exists(t1_ws_native) or not os.path.exists(t1_label_native):
            if verbose:
                print("  Registering WhiteStripe T1w and SynthSeg labels to nativepro...")
            with _ANTS_LOCK:
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

        native_outputs["T1w"] = (t1_ws_native, t1_label_native)

    if "FLAIR" in modalities:
        flair_ws_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-whitestripe_FLAIR.nii.gz")
        flair_label_native = os.path.join(output_dir, f"{bids_id}_space-nativepro_desc-synthseg_FLAIR.nii.gz")

        if not os.path.exists(flair_ws_native) or not os.path.exists(flair_label_native):
            if verbose:
                print("  Mapping WhiteStripe FLAIR and SynthSeg labels to nativepro...")
            affine_path = _find_micapipe_flair_transform(micapipe_dir, participant_id, session_id)
            with _ANTS_LOCK:
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

        native_outputs["FLAIR"] = (flair_ws_native, flair_label_native)

    for modality, (ws_path, label_path) in native_outputs.items():
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
    modalities=None,
    threads=1,
    verbose=True,
):
    """
    Prepare per-subject T1w/FLAIR WhiteStripe volumes and SynthSeg masks.

    Only modalities requested by ``modalities`` are required. Raw BIDS scans are
    preferred when ``raw_dir`` is supplied. Otherwise the nativepro micapipe
    images are normalized in place.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    modalities = list(modalities or ("T1w", "FLAIR"))

    if raw_dir:
        return _prepare_from_raw(
            participant_id=participant_id,
            session_id=session_id,
            raw_dir=raw_dir,
            output_dir=output_dir,
            micapipe_dir=micapipe_dir,
            tmp_dir=tmp_dir,
            modalities=modalities,
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
        modalities=modalities,
        threads=threads,
    )


def prepare_t1w_flair_wmmean(participant_id, session_id, output_dir, micapipe_dir,
                             modalities=None, threads=1, verbose=True):
    """WM-mean normalize nativepro T1w/FLAIR (SynthSeg WM mask + whole-WM tissue mean).

    Mirrors the WhiteStripe nativepro path but rescales by the trimmed mean/SD of
    the entire white-matter tissue rather than the WM stripe mode; output is
    ``desc-wmmean``.
    """
    os.makedirs(output_dir, exist_ok=True)
    modalities = list(modalities or ("T1w", "FLAIR"))
    bids_id = f"{participant_id}_{session_id}"
    source_paths = _nativepro_source_paths(participant_id, session_id, micapipe_dir)

    outputs = {}
    for modality in modalities:
        source_path = source_paths.get(modality)
        if not source_path or not os.path.exists(source_path):
            raise ValueError(f"Missing {modality} source image: {source_path}")
        label_path = os.path.join(
            output_dir, f"{bids_id}_space-nativepro_desc-synthseg_{modality}.nii.gz")
        out_path = os.path.join(
            output_dir, f"{bids_id}_space-nativepro_desc-wmmean_{modality}.nii.gz")
        stats_path = os.path.join(output_dir, f"{bids_id}_desc-wmmean_{modality}.json")
        _run_synthseg(source_path, label_path, threads=threads)
        wm_mean_normalize_image(source_path, label_path, out_path, stats_path=stats_path)
        outputs[modality] = {"source": source_path, "wmmean": out_path, "synthseg": label_path}
    return outputs


def ensure_synthseg_csf(participant_id, session_id, output_dir, micapipe_dir,
                        modalities=None, threads=1, verbose=True):
    """Guarantee desc-synthseg label + desc-synthsegCsf mask exist for each modality.

    RAVEL needs the SynthSeg CSF control region regardless of the subject-level
    normalization. The WhiteStripe path already writes these, but the identity
    (none) and wmmean paths do not write the CSF mask, so RAVEL composed on those
    subject norms calls this to fill the gap (idempotent; skips existing files).
    """
    os.makedirs(output_dir, exist_ok=True)
    modalities = list(modalities or ("T1w", "FLAIR"))
    bids_id = f"{participant_id}_{session_id}"
    source_paths = _nativepro_source_paths(participant_id, session_id, micapipe_dir)
    for modality in modalities:
        source_path = source_paths.get(modality)
        if not source_path or not os.path.exists(source_path):
            raise ValueError(f"Missing {modality} source image for SynthSeg CSF: {source_path}")
        label_path = os.path.join(
            output_dir, f"{bids_id}_space-nativepro_desc-synthseg_{modality}.nii.gz")
        csf_mask_path = os.path.join(
            output_dir, f"{bids_id}_space-nativepro_desc-synthsegCsf_{modality}.nii.gz")
        if not (os.path.exists(label_path) and os.path.getsize(label_path) > 0):
            if verbose:
                print(f"  Running SynthSeg for RAVEL CSF ({modality}, {bids_id})...")
            _run_synthseg(source_path, label_path, threads=threads)
        if not (os.path.exists(csf_mask_path) and os.path.getsize(csf_mask_path) > 0):
            _label_mask(label_path, SYNTHSEG_CSF_LABELS, csf_mask_path)


def requested_ravel_modalities(features):
    bases = {str(feature).lower().replace("*blur", "") for feature in features}
    modalities = []
    if "t1w" in bases:
        modalities.append("T1w")
    if "flair" in bases:
        modalities.append("FLAIR")
    return modalities


def normalize_normalization_mode(normalization):
    mode = str(normalization or "ravel").lower().replace("-", "").replace("_", "")
    if mode in {"ravel"}:
        return "ravel"
    if mode in {"whitestripe", "ws"}:
        return "whitestripe"
    if mode in {"none", "self", "selfnorm", "selfnormalized", "raw", "identity"}:
        return "none"
    if mode in {"wmmean", "fcm", "fcmwm", "wmmeannorm", "tissuemean"}:
        return "wmmean"
    if mode in {"nyul", "nyuludupa", "hist", "histmatch", "piecewiselinear"}:
        return "nyul"
    raise ValueError(
        "Unsupported T1w/FLAIR normalization mode "
        f"'{normalization}'. Use 'none' (self-normalized), 'whitestripe', "
        "'wmmean', 'nyul', or 'ravel'."
    )


# --- Composable normalization (subject-level THEN dataset-level) -------------
# A composition pairs a subject-level intensity normalization (per-subject:
# none / whitestripe / wmmean) with a dataset-level harmonization fit across the
# control cohort (none / ravel / nyul). The canonical desc- label for a pair is
# the subject label alone when dataset_norm == "none", otherwise the CamelCase
# join "<subject><Dataset>" (e.g. whitestripeRavel, noneNyul, wmmeanRavel). The
# subject-level phase writes desc-<subject_norm>; the dataset-level phase reads
# that intermediate and writes desc-<composite>; extraction reads desc-<composite>.
_SUBJECT_NORMS = ("none", "whitestripe", "wmmean")
_DATASET_NORMS = ("none", "ravel", "nyul")


def compose_normalization_label(subject_norm, dataset_norm):
    """Canonical desc- label for a (subject_norm, dataset_norm) composition."""
    s = normalize_normalization_mode(subject_norm)
    d_raw = str(dataset_norm or "none").lower().replace("-", "").replace("_", "")
    d = "none" if d_raw in ("", "none") else normalize_normalization_mode(d_raw)
    if s not in _SUBJECT_NORMS:
        raise ValueError(f"subject_norm must be one of {_SUBJECT_NORMS}, got {subject_norm!r}")
    if d not in _DATASET_NORMS:
        raise ValueError(f"dataset_norm must be one of {_DATASET_NORMS}, got {dataset_norm!r}")
    if d == "none":
        return s
    # Preserve the legacy single-mode desc- names for the two historical pairs so
    # existing desc-ravel / desc-nyul bases and example scripts keep working.
    if s == "whitestripe" and d == "ravel":
        return "ravel"
    if s == "none" and d == "nyul":
        return "nyul"
    return f"{s}{d.capitalize()}"


def decompose_normalization_label(label):
    """Inverse of :func:`compose_normalization_label`.

    Accepts composite CamelCase labels (whitestripeRavel), the subject-only
    labels (none/whitestripe/wmmean), and the legacy single-mode aliases
    ("ravel" == whitestripe+ravel, "nyul" == none+nyul).
    """
    raw = str(label or "none")
    for d in ("ravel", "nyul"):
        suffix = d.capitalize()
        if raw.endswith(suffix) and len(raw) > len(suffix):
            head = raw[: -len(suffix)]
            try:
                s = normalize_normalization_mode(head)
            except ValueError:
                continue
            if s in _SUBJECT_NORMS:
                return s, d
    mode = normalize_normalization_mode(raw)
    if mode in _SUBJECT_NORMS:
        return mode, "none"
    if mode == "ravel":  # legacy single-mode: WhiteStripe then RAVEL
        return "whitestripe", "ravel"
    if mode == "nyul":  # legacy single-mode: Nyul on raw nativepro
        return "none", "nyul"
    raise ValueError(f"Cannot decompose normalization label {label!r}")


def resolve_normalization_desc(normalization):
    """Return the on-disk desc- string for a (possibly composite) label."""
    try:
        return normalize_normalization_mode(normalization)
    except ValueError:
        # Composite label (e.g. whitestripeRavel): validate then use verbatim.
        decompose_normalization_label(normalization)
        return str(normalization)


def _link_or_copy(source_path, dest_path):
    """Symlink ``source_path`` to ``dest_path`` (copy if symlinks are unavailable)."""
    if os.path.islink(dest_path) or os.path.exists(dest_path):
        try:
            if os.path.samefile(dest_path, source_path):
                return dest_path
        except OSError:
            pass
        os.remove(dest_path)
    try:
        os.symlink(os.path.abspath(source_path), dest_path)
    except OSError:
        shutil.copy2(source_path, dest_path)
    return dest_path


def prepare_t1w_flair_identity(
    participant_id,
    session_id,
    output_dir,
    micapipe_dir,
    modalities=None,
    verbose=True,
):
    """Expose raw nativepro T1w/FLAIR as ``desc-none`` maps (no intensity normalization).

    This is the "self-normalization" processing mode: neither WhiteStripe nor
    RAVEL is applied. The raw micapipe nativepro T1w/FLAIR volumes are linked into
    the subject ``maps`` directory under the ``desc-none`` name so the existing
    volume-to-surface sampling picks them up unchanged. Any intensity
    standardization for these non-quantitative modalities is expected to happen
    downstream at analysis time (spatial/robust z-score, depth self-normalization,
    interhemispheric asymmetry, etc.).
    """
    os.makedirs(output_dir, exist_ok=True)
    modalities = list(modalities or ("T1w", "FLAIR"))
    bids_id = f"{participant_id}_{session_id}"
    subject_micapipe_dir = os.path.join(micapipe_dir, participant_id, session_id)
    source_paths = {
        "T1w": os.path.join(subject_micapipe_dir, "anat", f"{bids_id}_space-nativepro_T1w.nii.gz"),
        "FLAIR": os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-flair.nii.gz"),
    }

    outputs = {}
    for modality in modalities:
        if modality not in source_paths:
            raise ValueError(f"Unsupported identity-normalization modality: {modality}")
        source_path = source_paths[modality]
        if not os.path.exists(source_path):
            raise ValueError(f"Missing {modality} source image: {source_path}")

        dest_path = os.path.join(
            output_dir, f"{bids_id}_space-nativepro_desc-none_{modality}.nii.gz"
        )
        _link_or_copy(source_path, dest_path)
        outputs[modality] = {"source": source_path, "none": dest_path}
        if verbose:
            print(f"  Using raw (self-normalized) nativepro {modality}: {dest_path}")

    return outputs


def get_normalized_modality_path(subject_output_dir, subject_micapipe_dir, modality, normalization="ravel"):
    bids_id = "_".join(subject_output_dir.rstrip(os.sep).split(os.sep)[-2:])
    maps_dir = os.path.join(subject_output_dir, "maps")
    desc = resolve_normalization_desc(normalization)

    if modality == "T1w":
        return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{desc}_T1w.nii.gz")
    elif modality == "FLAIR":
        return os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{desc}_FLAIR.nii.gz")

    raise ValueError(f"Unsupported normalized modality: {modality}")


def _reference_image_path():
    return os.path.join(os.path.dirname(__file__), "data", "MNI152_T1_0.8mm.nii.gz")


def _subject_maps_dir(output_directory, participant_id, session_id):
    return os.path.join(output_directory, participant_id, session_id, "maps")


def _ravel_dir(output_directory):
    return os.path.join(output_directory, "ravel")


def _shared_ravel_mni_subject_dir(
    output_directory,
    participant_id,
    session_id,
    input_desc,
    source_root=None,
):
    """Return a fold/base-independent cache directory for RAVEL MNI inputs.

    Subject-level normalized images do not depend on the control fold, exclusion
    arm, or smoothing setting.  Replacing the ``zbrains_*`` base component with a
    fixed cache component lets those identical inputs be transformed only once.
    Direct callers whose output path has no zbrains component retain the legacy
    per-output-directory layout.
    """
    parts = os.path.abspath(output_directory).split(os.sep)
    for i, part in enumerate(parts):
        if part.startswith("zbrains_") or part.startswith("zbrains-"):
            cache_desc = "".join(
                char if char.isalnum() or char in "._-" else "_"
                for char in str(input_desc)
            )
            parts[i] = "ravel_mni_cache"
            cache_root = os.sep.join(parts[: i + 1])
            source_key = hashlib.sha256(
                os.path.realpath(source_root or output_directory).encode("utf-8")
            ).hexdigest()[:12]
            return os.path.join(
                cache_root,
                source_key,
                cache_desc,
                participant_id,
                session_id,
            )
    return None


def _content_identity(path):
    """Content identity that remains stable across identical copied inputs."""
    return {
        "size": int(os.path.getsize(path)),
        "sha256": _sha256(path),
    }


def _ravel_mni_cache_payload(
    native_image,
    native_csf,
    reference_path,
    native_to_mni,
    native_to_mni_invert,
):
    return {
        "cache_version": 1,
        "native_image": _content_identity(native_image),
        "native_csf": _content_identity(native_csf),
        "reference": _file_identity(reference_path),
        "native_to_mni": [
            _file_identity(path) for path in native_to_mni
        ],
        "native_to_mni_invert": [bool(value) for value in native_to_mni_invert],
    }


def _ravel_mni_cache_is_valid(
    manifest_path,
    expected_payload,
    image_path,
    csf_path,
    reference_shape,
):
    if not (
        _image_matches_reference_grid(image_path, reference_shape)
        and _image_matches_reference_grid(csf_path, reference_shape)
        and os.path.exists(manifest_path)
    ):
        return False
    try:
        with open(manifest_path, "r", encoding="utf-8") as stream:
            return json.load(stream) == expected_payload
    except (OSError, ValueError, TypeError):
        return False


def ravel_model_path(output_directory, modality):
    return os.path.join(_ravel_dir(output_directory), f"{modality}_ravel_model.npz")


def ravel_models_exist(output_directory, modalities):
    for modality in modalities:
        path = ravel_model_path(output_directory, modality)
        if not os.path.exists(path):
            return False
        try:
            with np.load(path, allow_pickle=True) as model:
                if "beta_w" not in model or "model_version" not in model:
                    return False
        except Exception:
            return False
    return True


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
    threads=1,
    verbose=True,
    input_desc="whitestripe",
    output_desc="ravel",
):
    bids_id = f"{participant_id}_{session_id}"
    maps_dir = _subject_maps_dir(output_directory, participant_id, session_id)
    ravel_subject_dir = os.path.join(maps_dir, "ravel")
    os.makedirs(ravel_subject_dir, exist_ok=True)

    reference_path = _reference_image_path()
    reference_shape = nib.load(reference_path).shape
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
    transform_source = "micapipe_nativepro_to_mni152_0.8mm_syn"
    native_to_mni = []
    native_to_mni_invert = []
    mni_to_native = []
    mni_to_native_invert = []

    if os.path.exists(transform_json):
        with open(transform_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        if payload.get("source") == transform_source and payload.get("reference") == reference_path:
            native_to_mni = payload.get("native_to_mni", [])
            native_to_mni_invert = payload.get("native_to_mni_invert", [])
            mni_to_native = payload.get("mni_to_native", [])
            mni_to_native_invert = payload.get("mni_to_native_invert", [])

    transforms_ready = (
        native_to_mni
        and mni_to_native
        and all(os.path.exists(p) for p in native_to_mni)
        and all(os.path.exists(p) for p in mni_to_native)
    )
    if not transforms_ready:
        transforms = _find_micapipe_nativepro_mni_transforms(micapipe_dir, participant_id, session_id)
        native_to_mni = transforms["native_to_mni"]
        native_to_mni_invert = transforms["native_to_mni_invert"]
        mni_to_native = transforms["mni_to_native"]
        mni_to_native_invert = transforms["mni_to_native_invert"]
        _write_json(
            transform_json,
            {
                "source": transform_source,
                "native_to_mni": native_to_mni,
                "native_to_mni_invert": native_to_mni_invert,
                "mni_to_native": mni_to_native,
                "mni_to_native_invert": mni_to_native_invert,
                "reference": reference_path,
            },
        )

    outputs = {
        "reference": reference_path,
        "native_reference": native_ref,
        "native_to_mni": native_to_mni,
        "native_to_mni_invert": native_to_mni_invert,
        "mni_to_native": mni_to_native,
        "mni_to_native_invert": mni_to_native_invert,
        "modalities": {},
    }

    shared_cache_dir = _shared_ravel_mni_subject_dir(
        output_directory,
        participant_id,
        session_id,
        input_desc,
        source_root=micapipe_dir,
    )
    if shared_cache_dir:
        os.makedirs(shared_cache_dir, exist_ok=True)

    for modality in modalities:
        native_ws = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{input_desc}_{modality}.nii.gz")
        native_csf = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-synthsegCsf_{modality}.nii.gz")
        mni_input_dir = shared_cache_dir or ravel_subject_dir
        mni_ws = os.path.join(mni_input_dir, f"{bids_id}_space-MNI152_desc-{input_desc}_{modality}.nii.gz")
        mni_csf = os.path.join(mni_input_dir, f"{bids_id}_space-MNI152_desc-synthsegCsf_{modality}.nii.gz")
        native_ravel = os.path.join(maps_dir, f"{bids_id}_space-nativepro_desc-{output_desc}_{modality}.nii.gz")
        mni_ravel = os.path.join(ravel_subject_dir, f"{bids_id}_space-MNI152_desc-{output_desc}_{modality}.nii.gz")

        if not os.path.exists(native_ws):
            raise ValueError(f"Missing subject-normalized ({input_desc}) {modality} before RAVEL: {native_ws}")
        if not os.path.exists(native_csf):
            raise ValueError(f"Missing SynthSeg CSF mask for RAVEL: {native_csf}")

        if shared_cache_dir:
            manifest_path = os.path.join(
                shared_cache_dir,
                f"{bids_id}_desc-{input_desc}_{modality}_ravel-mni-cache.json",
            )
            expected_payload = _ravel_mni_cache_payload(
                native_ws,
                native_csf,
                reference_path,
                native_to_mni,
                native_to_mni_invert,
            )
            with _exclusive_file_lock(f"{manifest_path}.lock"):
                if not _ravel_mni_cache_is_valid(
                    manifest_path,
                    expected_payload,
                    mni_ws,
                    mni_csf,
                    reference_shape,
                ):
                    _forward_transform_ravel_cli(
                        native_ws,
                        mni_ws,
                        reference_path,
                        native_to_mni,
                        native_to_mni_invert,
                        interpolator="Linear",
                        threads=threads,
                        atomic=True,
                    )
                    _forward_transform_ravel_cli(
                        native_csf,
                        mni_csf,
                        reference_path,
                        native_to_mni,
                        native_to_mni_invert,
                        interpolator="NearestNeighbor",
                        threads=threads,
                        atomic=True,
                    )
                    _write_json_atomic(manifest_path, expected_payload)
                elif verbose:
                    print(
                        f"  Reusing cached RAVEL MNI inputs for "
                        f"{bids_id} {modality}."
                    )
        else:
            if not _image_matches_reference_grid(mni_ws, reference_shape):
                _forward_transform_ravel_cli(
                    native_ws,
                    mni_ws,
                    reference_path,
                    native_to_mni,
                    native_to_mni_invert,
                    interpolator="Linear",
                    threads=threads,
                )

            if not _image_matches_reference_grid(mni_csf, reference_shape):
                _forward_transform_ravel_cli(
                    native_csf,
                    mni_csf,
                    reference_path,
                    native_to_mni,
                    native_to_mni_invert,
                    interpolator="NearestNeighbor",
                    threads=threads,
                )

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


def _iter_index_chunks(indices, chunk_size=RAVEL_CHUNK_SIZE):
    for start in range(0, len(indices), chunk_size):
        yield indices[start:start + chunk_size]


def _progress_interval(total_steps, updates=10):
    return max(1, int(np.ceil(max(1, total_steps) / float(updates))))


def _maybe_print_progress(verbose, label, current, total, interval):
    if verbose and (current == 1 or current == total or current % interval == 0):
        print(f"  {label}: {current}/{total}")


def _matrix_from_vectors(vectors, indices):
    return np.vstack([np.asarray(vector[indices], dtype=np.float32) for vector in vectors])


def _create_brain_vector_memmaps(image_paths, brain_mask, tmp_dir, prefix):
    vectors = []
    paths = []
    for i, path in enumerate(image_paths):
        data = np.asarray(nib.load(path).dataobj, dtype=np.float32)
        if data.shape != brain_mask.shape:
            raise ValueError(f"RAVEL image shape mismatch for {path}: {data.shape} vs {brain_mask.shape}")
        mmap_path = os.path.join(tmp_dir, f"{prefix}_{i}.dat")
        vector = np.memmap(mmap_path, dtype=np.float32, mode="w+", shape=(int(brain_mask.sum()),))
        vector[:] = data[brain_mask].astype(np.float32)
        vector.flush()
        vectors.append(vector)
        paths.append(mmap_path)
    return vectors, paths


def _write_corrected_memmaps_to_images(corrected_vectors, brain_mask, reference_img, output_paths):
    for vector, output_path in zip(corrected_vectors, output_paths):
        _write_brain_vector(np.asarray(vector, dtype=np.float32), brain_mask, reference_img, output_path)


def _close_memmaps(vectors):
    for vector in vectors:
        try:
            vector.flush()
        except Exception:
            pass
        mmap_obj = getattr(vector, "_mmap", None)
        if mmap_obj is not None:
            try:
                mmap_obj.close()
            except Exception:
                pass


def _image_matches_reference_grid(image_path, reference_shape):
    if not os.path.exists(image_path):
        return False
    try:
        return nib.load(image_path).shape == tuple(reference_shape)
    except Exception:
        return False


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


def _inverse_transform_ravel(
    mni_path,
    native_path,
    native_reference,
    mni_to_native,
    mni_to_native_invert=None,
    threads=1,
):
    _configure_native_threads(threads)
    import ants

    with _ANTS_LOCK:
        corrected = ants.apply_transforms(
            fixed=ants.image_read(native_reference),
            moving=ants.image_read(mni_path),
            transformlist=mni_to_native,
            whichtoinvert=mni_to_native_invert,
            interpolator="linear",
        )
        ants.image_write(corrected, native_path)
    return native_path


def _ants_apply_transform_args(transformlist, whichtoinvert=None):
    whichtoinvert = whichtoinvert or [False] * len(transformlist)
    args = []
    for transform, invert in zip(transformlist, whichtoinvert):
        if invert:
            args.extend(["-t", f"[{transform},1]"])
        else:
            args.extend(["-t", transform])
    return args


def _thread_limited_environment(threads):
    threads = max(1, int(threads or 1))
    env = os.environ.copy()
    env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(threads)
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        env[variable] = str(threads)
    return env


def _atomic_nifti_path(path):
    suffix = ".nii.gz" if path.endswith(".nii.gz") else os.path.splitext(path)[1]
    stem = path[: -len(suffix)] if suffix else path
    return f"{stem}.tmp-{os.getpid()}-{threading.get_ident()}{suffix}"


def _cast_nifti_in_place(path, dtype):
    """Cast a transform result while retaining its output grid and metadata."""
    image = nib.load(path)
    if image.get_data_dtype() == np.dtype(dtype):
        return path
    data = np.asarray(image.dataobj, dtype=dtype)
    header = image.header.copy()
    header.set_data_dtype(dtype)
    cast_path = _atomic_nifti_path(path + ".cast.nii.gz")
    try:
        nib.save(nib.Nifti1Image(data, image.affine, header), cast_path)
        os.replace(cast_path, path)
    finally:
        if os.path.exists(cast_path):
            os.remove(cast_path)
    return path


def _apply_ravel_transform_cli(
    moving_path,
    output_path,
    fixed_path,
    transformlist,
    whichtoinvert=None,
    interpolator="Linear",
    threads=1,
    atomic=False,
    output_dtype=None,
):
    """Apply an existing ANTs transform without serializing independent jobs."""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    write_path = _atomic_nifti_path(output_path) if atomic else output_path
    ants_apply = shutil.which("antsApplyTransforms")
    try:
        if ants_apply:
            cmd = [
                ants_apply,
                "-d",
                "3",
                "-i",
                moving_path,
                "-r",
                fixed_path,
                "-o",
                write_path,
                "-n",
                interpolator,
                "-z",
                "1",
                "-v",
                "0",
                "--float",
                "0",
                "-e",
                "0",
                "-f",
                "0",
            ]
            cmd.extend(_ants_apply_transform_args(transformlist, whichtoinvert))
            subprocess.run(
                cmd,
                check=True,
                env=_thread_limited_environment(threads),
            )
        else:
            # ANTsPy is retained as a compatibility fallback.  It is not safe to
            # invoke concurrently in this process, hence the narrow fallback lock.
            import ants

            python_interpolator = {
                "Linear": "linear",
                "NearestNeighbor": "nearestNeighbor",
            }.get(interpolator, interpolator)
            with _ANTS_LOCK:
                _configure_native_threads(threads)
                transformed = ants.apply_transforms(
                    fixed=ants.image_read(fixed_path),
                    moving=ants.image_read(moving_path),
                    transformlist=transformlist,
                    whichtoinvert=whichtoinvert,
                    interpolator=python_interpolator,
                )
                ants.image_write(transformed, write_path)

        if output_dtype is not None:
            _cast_nifti_in_place(write_path, output_dtype)
        if atomic:
            os.replace(write_path, output_path)
        return output_path
    finally:
        if atomic and os.path.exists(write_path):
            os.remove(write_path)


def _forward_transform_ravel_cli(
    native_path,
    mni_path,
    mni_reference,
    native_to_mni,
    native_to_mni_invert=None,
    interpolator="Linear",
    threads=1,
    atomic=False,
):
    return _apply_ravel_transform_cli(
        native_path,
        mni_path,
        mni_reference,
        native_to_mni,
        whichtoinvert=native_to_mni_invert,
        interpolator=interpolator,
        threads=threads,
        atomic=atomic,
        output_dtype=np.float32,
    )


def _inverse_transform_ravel_cli(
    mni_path,
    native_path,
    native_reference,
    mni_to_native,
    mni_to_native_invert=None,
    threads=1,
    output_dtype=None,
):
    return _apply_ravel_transform_cli(
        mni_path,
        native_path,
        native_reference,
        mni_to_native,
        whichtoinvert=mni_to_native_invert,
        interpolator="Linear",
        threads=threads,
        output_dtype=output_dtype,
    )


def _resolve_ravel_tmp_parent(output_directory, tmp_dir=None):
    """Choose scratch for transient RAVEL memmaps without changing their dtype."""
    candidate = (
        tmp_dir
        or os.environ.get("ZBRAINS_RAVEL_TMPDIR")
        or os.environ.get("SLURM_TMPDIR")
    )
    if candidate:
        candidate = os.path.abspath(os.fspath(candidate))
        os.makedirs(candidate, exist_ok=True)
        return candidate
    return _ravel_dir(output_directory)


def fit_and_apply_ravel_to_controls(
    subjects,
    output_directory,
    micapipe_dir,
    modalities,
    k=1,
    verbose=True,
    n_jobs=1,
    threads=1,
    input_desc="whitestripe",
    output_desc="ravel",
    total_threads=None,
    tmp_dir=None,
):
    """
    Fit RAVEL on the control cohort and write RAVEL-corrected nativepro volumes.

    This follows the core RAVEL correction used by jfortin1/RAVEL, but streams
    voxels in chunks. A single subject-level nuisance basis is estimated from
    the CSF control region, then applied consistently to all brain voxels.
    """
    if not subjects:
        raise ValueError("Cannot fit RAVEL without control subjects")

    os.makedirs(_ravel_dir(output_directory), exist_ok=True)
    modalities = list(modalities)
    total_threads = max(1, int(total_threads or threads or 1))
    ravel_tmp_parent = _resolve_ravel_tmp_parent(output_directory, tmp_dir=tmp_dir)

    def ensure_subject_mni_inputs(subject):
        return subject, _ensure_mni_inputs(
            subject[0],
            subject[1],
            output_directory,
            micapipe_dir,
            modalities,
            threads=threads,
            verbose=verbose,
            input_desc=input_desc,
            output_desc=output_desc,
        )

    ravel_n_jobs = _coerce_n_jobs(n_jobs, len(subjects))
    if verbose and ravel_n_jobs > 1:
        print(f"Preparing control RAVEL MNI inputs with {ravel_n_jobs} parallel jobs...")

    if ravel_n_jobs == 1:
        mni_inputs = {}
        input_interval = _progress_interval(len(subjects))
        for i, subject in enumerate(subjects, start=1):
            if verbose:
                pid, sid = subject
                print(f"Preparing RAVEL MNI inputs for {pid}/{sid} ({i}/{len(subjects)})...")
            subject_key, subject_inputs = ensure_subject_mni_inputs(subject)
            mni_inputs[subject_key] = subject_inputs
            _maybe_print_progress(verbose, "Prepared RAVEL MNI inputs", i, len(subjects), input_interval)
    else:
        mni_inputs = dict(
            Parallel(n_jobs=ravel_n_jobs, backend="threading")(
                delayed(ensure_subject_mni_inputs)(subject)
                for subject in subjects
            )
        )

    # Preserve the pre-optimization parent-process thread budget for the model
    # fit itself. Per-transform subprocess limits must not leak into NumPy's
    # deterministic fitting path.
    _configure_native_threads(total_threads)
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
        control_mask, control_strategy = _control_mask_from_csf(csf_paths, brain_mask)
        n_subjects = len(subjects)
        n_brain_voxels = int(brain_mask.sum())
        brain_indices = np.arange(n_brain_voxels)
        control_indices = np.flatnonzero(control_mask)
        X = np.ones((n_subjects, 1), dtype=np.float64)
        x_pinv = np.linalg.pinv(X)
        residual_maker = np.eye(n_subjects, dtype=np.float64) - X @ x_pinv

        with tempfile.TemporaryDirectory(
            prefix=f"ravel_{modality}_",
            dir=ravel_tmp_parent,
            ignore_cleanup_errors=True,
        ) as ravel_tmp:
            image_vectors, _ = _create_brain_vector_memmaps(
                image_paths,
                brain_mask,
                ravel_tmp,
                f"{modality}_input",
            )

            if verbose:
                print(
                    f"  Estimating {modality} RAVEL nuisance factors from "
                    f"{len(control_indices)} CSF voxels in chunks of {RAVEL_CHUNK_SIZE}..."
                )

            gram = np.zeros((n_subjects, n_subjects), dtype=np.float64)
            control_chunks = list(_iter_index_chunks(control_indices))
            control_interval = _progress_interval(len(control_chunks))
            for i, control_chunk in enumerate(control_chunks, start=1):
                Y = _matrix_from_vectors(image_vectors, control_chunk).astype(np.float64)
                R = residual_maker @ Y
                gram += R @ R.T
                _maybe_print_progress(
                    verbose,
                    f"{modality} CSF Gram chunks",
                    i,
                    len(control_chunks),
                    control_interval,
                )

            if verbose:
                print(f"  Solving {modality} subject-level nuisance eigensystem...")
            eigvals, eigvecs = np.linalg.eigh(gram)
            order = np.argsort(eigvals)[::-1]
            eigvals = eigvals[order]
            eigvecs = eigvecs[:, order]
            k_eff = min(int(k), eigvecs.shape[1], max(1, n_subjects - X.shape[1]))
            if k_eff < 1:
                raise ValueError("RAVEL requires at least one estimable unwanted factor")
            W = eigvecs[:, :k_eff].astype(np.float64)
            design = np.column_stack([X, W]).astype(np.float64)
            design_pinv = np.linalg.pinv(design)

            beta_x = np.zeros((X.shape[1], n_brain_voxels), dtype=np.float32)
            beta_w = np.zeros((k_eff, n_brain_voxels), dtype=np.float32)
            corrected_paths = [
                os.path.join(ravel_tmp, f"{modality}_corrected_{i}.dat")
                for i in range(n_subjects)
            ]
            corrected_vectors = [
                np.memmap(path, dtype=np.float32, mode="w+", shape=(n_brain_voxels,))
                for path in corrected_paths
            ]

            if verbose:
                print(f"  Applying {modality} RAVEL correction to brain voxels in chunks...")

            brain_chunks = list(_iter_index_chunks(brain_indices))
            brain_interval = _progress_interval(len(brain_chunks))
            for i, brain_chunk in enumerate(brain_chunks, start=1):
                Y = _matrix_from_vectors(image_vectors, brain_chunk).astype(np.float64)
                coef = design_pinv @ Y
                beta_x[:, brain_chunk] = coef[:X.shape[1], :].astype(np.float32)
                beta_w[:, brain_chunk] = coef[X.shape[1]:, :].astype(np.float32)
                Y_corrected = Y - W @ coef[X.shape[1]:, :]
                for subject_index, corrected_vector in enumerate(corrected_vectors):
                    corrected_vector[brain_chunk] = Y_corrected[subject_index, :].astype(np.float32)
                _maybe_print_progress(
                    verbose,
                    f"{modality} brain correction chunks",
                    i,
                    len(brain_chunks),
                    brain_interval,
                )

            for corrected_vector in corrected_vectors:
                corrected_vector.flush()

            model_path = ravel_model_path(output_directory, modality)
            np.savez_compressed(
                model_path,
                brain_mask=brain_mask,
                control_mask=control_mask,
                X=X.astype(np.float32),
                W=W.astype(np.float32),
                beta_x=beta_x,
                beta_w=beta_w,
                eigvals=eigvals[:k_eff].astype(np.float32),
                subjects=np.array([f"{pid}_{sid}" for pid, sid in subjects]),
                chunk_size=np.array([RAVEL_CHUNK_SIZE], dtype=np.int32),
                model_version=np.array(["chunked_gram_v1"]),
            )
            _write_json(
                model_path.replace(".npz", ".json"),
                {
                    "modality": modality,
                    "model_version": "chunked_gram_v1",
                    "k": k_eff,
                    "reference": reference_path,
                    "n_subjects": len(subjects),
                    "n_brain_voxels": n_brain_voxels,
                    "n_control_voxels": int(control_mask.sum()),
                    "control_mask_strategy": control_strategy,
                    "chunk_size": RAVEL_CHUNK_SIZE,
                    "subjects": [f"{pid}_{sid}" for pid, sid in subjects],
                },
            )

            output_paths = [
                mni_inputs[subject]["modalities"][modality]["mni_ravel"]
                for subject in subjects
            ]
            if verbose:
                print(f"  Writing {modality} corrected MNI images for {len(subjects)} controls...")
            _write_corrected_memmaps_to_images(
                corrected_vectors,
                brain_mask,
                reference_img,
                output_paths,
            )
            _close_memmaps(corrected_vectors)
            _close_memmaps(image_vectors)
            del corrected_vectors
            del image_vectors

            inverse_jobs = max(1, min(total_threads, len(subjects)))
            inverse_threads = max(1, total_threads // inverse_jobs)
            if verbose:
                print(
                    f"  Inverse-transforming {modality} RAVEL images with "
                    f"{inverse_jobs} parallel jobs ({inverse_threads} ITK thread/job)..."
                )

            def inverse_subject_to_native(subject):
                paths = mni_inputs[subject]["modalities"][modality]
                _inverse_transform_ravel_cli(
                    paths["mni_ravel"],
                    paths["native_ravel"],
                    mni_inputs[subject]["native_reference"],
                    mni_inputs[subject]["mni_to_native"],
                    mni_inputs[subject]["mni_to_native_invert"],
                    threads=inverse_threads,
                )
                return subject

            if inverse_jobs == 1:
                inverse_interval = _progress_interval(len(subjects))
                for i, subject in enumerate(subjects, start=1):
                    inverse_subject_to_native(subject)
                    _maybe_print_progress(
                        verbose,
                        f"{modality} inverse transforms to nativepro",
                        i,
                        len(subjects),
                        inverse_interval,
                    )
            else:
                completed_inverse = 0
                inverse_interval = _progress_interval(len(subjects))
                # THREADING, not loky: inverse_subject_to_native only shells out to
                # antsApplyTransforms (subprocess.run, which releases the GIL), so it
                # needs no numpy in the worker. A loky (fork) worker instead re-imports
                # numpy and dies with "CPU dispatcher tracer already initlized" ->
                # SIGSEGV. Threads avoid the fork + numpy double-init entirely (this is
                # the same backend the MNI-input prep above uses).
                for _ in Parallel(n_jobs=inverse_jobs, backend="threading", return_as="generator")(
                    delayed(inverse_subject_to_native)(subject)
                    for subject in subjects
                ):
                    completed_inverse += 1
                    _maybe_print_progress(
                        verbose,
                        f"{modality} inverse transforms to nativepro",
                        completed_inverse,
                        len(subjects),
                        inverse_interval,
                    )


def load_ravel_models(output_directory, modalities):
    """Load each fitted model once into read-only arrays for cohort application."""
    loaded = {}
    for modality in modalities:
        model_path = ravel_model_path(output_directory, modality)
        if not os.path.exists(model_path):
            raise ValueError(
                f"Missing RAVEL model for {modality}: {model_path}. "
                "Process the control dataset first."
            )
        with np.load(model_path, allow_pickle=True) as model:
            if "beta_w" not in model:
                raise ValueError(
                    f"RAVEL model for {modality} uses the old dense format. "
                    "Reprocess the control dataset to fit a chunked RAVEL model."
                )
            values = {
                "brain_mask": model["brain_mask"].astype(bool),
                "control_mask": model["control_mask"].astype(bool),
                "beta_x": model["beta_x"].astype(np.float32),
                "beta_w": model["beta_w"].astype(np.float32),
            }
        for value in values.values():
            value.setflags(write=False)
        loaded[modality] = values
    return loaded


def apply_ravel_model_to_subject(
    participant_id,
    session_id,
    output_directory,
    micapipe_dir,
    modalities,
    threads=1,
    verbose=True,
    input_desc="whitestripe",
    output_desc="ravel",
    models=None,
):
    os.makedirs(_ravel_dir(output_directory), exist_ok=True)
    modalities = list(modalities)
    mni_inputs = _ensure_mni_inputs(
        participant_id,
        session_id,
        output_directory,
        micapipe_dir,
        modalities,
        threads=threads,
        verbose=verbose,
        input_desc=input_desc,
        output_desc=output_desc,
    )

    reference_img = nib.load(_reference_image_path())
    loaded_models = models or load_ravel_models(output_directory, modalities)

    for modality in modalities:
        if modality not in loaded_models:
            raise ValueError(f"No preloaded RAVEL model was supplied for {modality}")
        model = loaded_models[modality]
        brain_mask = model["brain_mask"]
        control_mask = model["control_mask"]
        beta_x = model["beta_x"]
        beta_w = model["beta_w"]

        paths = mni_inputs["modalities"][modality]
        V = _load_brain_vector(paths["mni_whitestripe"], brain_mask)
        control_indices = np.flatnonzero(control_mask)
        X_subject = np.ones((1, beta_x.shape[0]), dtype=np.float32)
        control_residual = (
            V[control_indices][None, :]
            - X_subject @ beta_x[:, control_indices]
        )
        beta_w_control = beta_w[:, control_indices]
        nuisance_gram = beta_w_control @ beta_w_control.T
        W_subject = (control_residual @ beta_w_control.T) @ np.linalg.pinv(nuisance_gram)
        corrected = np.zeros_like(V, dtype=np.float32)

        brain_chunks = list(_iter_index_chunks(np.arange(V.size)))
        brain_interval = _progress_interval(len(brain_chunks))
        if verbose:
            print(
                f"  Applying fitted RAVEL model to {participant_id}/{session_id} "
                f"{modality} in {len(brain_chunks)} brain chunks..."
            )
        for i, brain_chunk in enumerate(brain_chunks, start=1):
            corrected[brain_chunk] = (
                V[brain_chunk][None, :]
                - W_subject @ beta_w[:, brain_chunk]
            ).ravel().astype(np.float32)
            _maybe_print_progress(
                verbose,
                f"{participant_id}/{session_id} {modality} RAVEL chunks",
                i,
                len(brain_chunks),
                brain_interval,
            )

        _write_brain_vector(corrected, brain_mask, reference_img, paths["mni_ravel"])
        if verbose:
            print(f"  Transforming {participant_id}/{session_id} {modality} RAVEL image back to nativepro...")
        _inverse_transform_ravel_cli(
            paths["mni_ravel"],
            paths["native_ravel"],
            mni_inputs["native_reference"],
            mni_inputs["mni_to_native"],
            mni_inputs["mni_to_native_invert"],
            threads=threads,
            output_dtype=np.float32,
        )
