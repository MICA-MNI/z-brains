"""k-nearest-neighbour nonparametric normative scoring.

A distance-based abnormality score that sits alongside the parametric normative
models (z-score, w-score, Gaussian process). For each surface vertex the control
cohort defines an empirical distribution; a patient vertex is scored by how far it
sits from its ``k`` nearest control observations, with ``k`` set to HALF the
available controls (``k = max(1, n_controls // 2)``). Age/sex are removed first by
per-vertex linear residualization (so the score is a nonparametric analogue of the
w-score), and the raw kNN distance is calibrated against the controls' own
leave-one-out kNN distances to yield a comparable, unit-scaled abnormality map.

Design notes
------------
* Nonparametric: no Gaussian assumption on the control distribution (robust to the
  skewed/multimodal control distributions some depth features have).
* Naturally multivariate: multi-depth features are scored on the joint depth vector
  per vertex (Euclidean kNN over depths), single-value features are scored in 1-D.
* Calibrated: the patient's mean distance to its k nearest controls is standardized
  (robust z) against the null of every control's leave-one-out mean kNN distance,
  so the output is on a comparable scale to the other scorers and higher == more
  abnormal (unsigned magnitude; downstream lesion detection uses magnitude).
"""
from __future__ import annotations

import json
import os

import nibabel as nib
import numpy as np

EPS = 1e-8
_VERTEX_CHUNK = 2048  # bounds the O(n_ref^2 * chunk) leave-one-out pairwise memory


def _standardize_covariates(demo_ref, demo_pat):
    """Standardize covariates by the control mean/SD; return design matrices
    ``X_ref`` (n_ref x 1+p) and ``x_pat`` (1+p,) with a leading intercept."""
    demo_ref = np.asarray(demo_ref, dtype=float)
    demo_pat = np.asarray(demo_pat, dtype=float).reshape(-1)
    mean = demo_ref.mean(axis=0)
    sd = demo_ref.std(axis=0)
    sd = np.where(np.isfinite(sd) & (sd > EPS), sd, 1.0)
    z_ref = (demo_ref - mean) / sd
    z_pat = (demo_pat - mean) / sd
    X_ref = np.column_stack([np.ones(demo_ref.shape[0]), z_ref])
    x_pat = np.concatenate([[1.0], z_pat])
    return X_ref, x_pat


def _residualize(reference_data, patient_data, X_ref, x_pat):
    """Remove the linear demographic trend per column. ``reference_data`` is
    (n_ref, M), ``patient_data`` is (M,); returns control + patient residuals with
    non-finite control entries imputed to the per-column control median so distances
    stay well defined."""
    Y = np.asarray(reference_data, dtype=float)
    y_pat = np.asarray(patient_data, dtype=float)
    finite = np.isfinite(Y)
    if not finite.all():
        col_median = np.nanmedian(np.where(finite, Y, np.nan), axis=0)
        col_median = np.where(np.isfinite(col_median), col_median, 0.0)
        Y = np.where(finite, Y, col_median[None, :])
    beta, *_ = np.linalg.lstsq(X_ref, Y, rcond=None)   # (1+p, M)
    resid_ref = Y - X_ref @ beta
    resid_pat = y_pat - x_pat @ beta                   # (M,)
    return resid_ref, resid_pat


def _knn_mean_distance_to_set(query, pool, k):
    """Mean Euclidean distance from each query row to its ``k`` nearest pool rows.

    ``query`` (nq, d), ``pool`` (npool, d). Returns (nq,). ``k`` is clamped to the
    pool size. Vectorized; caller chunks vertices to bound memory."""
    # pairwise squared distances (nq, npool)
    d2 = (np.square(query[:, None, :] - pool[None, :, :])).sum(axis=2)
    k = int(min(max(k, 1), pool.shape[0]))
    # k smallest per row via partial sort
    idx = np.argpartition(d2, kth=k - 1, axis=1)[:, :k]
    knn = np.take_along_axis(d2, idx, axis=1)
    return np.sqrt(np.maximum(knn, 0.0)).mean(axis=1)


def knn_abnormality(reference_residuals, patient_residuals, k):
    """Calibrated per-vertex kNN abnormality (>=~0; NaN where degenerate).

    ``reference_residuals`` is (n_ref, n_vert) or (n_ref, n_vert, d); patient is
    (n_vert,) or (n_vert, d). ``k`` neighbours (already = half the controls)."""
    ref = np.asarray(reference_residuals, dtype=float)
    pat = np.asarray(patient_residuals, dtype=float)
    if ref.ndim == 2:
        ref = ref[:, :, None]
        pat = pat[:, None]
    n_ref, n_vert, d = ref.shape
    if n_ref < 2:
        raise ValueError("kNN scoring requires at least two reference subjects")
    k_pat = int(min(max(1, k), n_ref))          # patient vs all controls
    k_loo = int(min(max(1, k), n_ref - 1))      # control vs the OTHER controls

    scores = np.full(n_vert, np.nan, dtype=float)
    for start in range(0, n_vert, _VERTEX_CHUNK):
        stop = min(start + _VERTEX_CHUNK, n_vert)
        R = ref[:, start:stop, :]               # (n_ref, c, d)
        P = pat[start:stop, :]                   # (c, d)
        c = stop - start
        # patient mean distance to its k nearest controls, per vertex
        d2_pat = np.square(P[None, :, :] - R).sum(axis=2)        # (n_ref, c)
        kp = min(k_pat, n_ref)
        part = np.partition(d2_pat, kth=kp - 1, axis=0)[:kp, :]
        Dp = np.sqrt(np.maximum(part, 0.0)).mean(axis=0)         # (c,)
        # leave-one-out control null: each control's mean distance to its k
        # nearest OTHER controls
        d2_cc = np.square(R[:, None, :, :] - R[None, :, :, :]).sum(axis=3)  # (n_ref,n_ref,c)
        diag = np.arange(n_ref)
        d2_cc[diag, diag, :] = np.inf           # exclude self
        kl = min(k_loo, n_ref - 1)
        part_cc = np.partition(d2_cc, kth=kl - 1, axis=1)[:, :kl, :]
        Dc = np.sqrt(np.maximum(part_cc, 0.0)).mean(axis=1)      # (n_ref, c)
        center = np.median(Dc, axis=0)                            # (c,)
        mad = np.median(np.abs(Dc - center[None, :]), axis=0)
        scale = 1.482602218505602 * mad
        good = np.isfinite(scale) & (scale > EPS)
        chunk_scores = np.full(c, np.nan, dtype=float)
        chunk_scores[good] = (Dp[good] - center[good]) / scale[good]
        # degenerate (all controls identical at this vertex): no signal -> 0
        chunk_scores[~good & np.isfinite(Dp)] = 0.0
        scores[start:stop] = chunk_scores
    return scores


def calculate_knn_maps(
    reference_data,
    patient_data,
    demographics_ref,
    demographics_pat,
    output_file,
    normative_columns,
    verbose,
    wscore_covariate_model,
    min_reference_subjects,
    reference_vertex_covariates,
    patient_vertex_covariates,
    wscore_preprocessing,
    blur_depth_model,
    intensity_depth_model,
    wscore_surface_smoothing_iterations,
    wscore_fit_cache,
    surface_smoother,
    prediction_variance_percentile=None,
):
    """kNN nonparametric normative scoring; save a GIFTI abnormality map.

    Mirrors :func:`gaussian_process.calculate_gaussian_process_maps`'s interface so
    it plugs into ``calculate_wscore_maps`` unchanged. Multi-depth features are
    scored on the joint depth vector per vertex; single-value features in 1-D.
    ``k = max(1, n_controls // 2)`` (half the available controls).
    """
    if reference_vertex_covariates is not None or patient_vertex_covariates is not None:
        raise ValueError("kNN normative model supports demographic covariates only; "
                         "disable use_curvature_covariates")

    reference_data = np.asarray(reference_data, dtype=float)
    patient_data = np.asarray(patient_data, dtype=float)
    if reference_data.ndim < 2:
        raise ValueError("reference_data must be (n_ref, n_vert[, n_depths])")
    if reference_data.shape[0] < 2:
        raise ValueError("kNN scoring requires at least two reference subjects")

    demo_ref = demographics_ref[normative_columns].to_numpy(dtype=float)
    demo_pat = np.asarray(demographics_pat[normative_columns], dtype=float).reshape(-1)
    demographic_mask = np.isfinite(demo_ref).all(axis=1)
    demo_ref = demo_ref[demographic_mask]
    reference_data = reference_data[demographic_mask]
    if not np.isfinite(demo_pat).all():
        raise ValueError("Patient demographics must be complete for kNN scoring")
    n_ref = reference_data.shape[0]
    if n_ref < max((min_reference_subjects or 2), 2):
        raise ValueError("Too few controls with complete demographics for kNN scoring")

    n_vert = reference_data.shape[1]
    n_depths = reference_data.shape[2] if reference_data.ndim > 2 else 1
    input_was_multidepth = reference_data.ndim > 2

    # Optional per-subject spatial normalization (matches the other scorers).
    def _spatial(values, robust):
        one_d = values.ndim == 1
        m = values[None, :] if one_d else values
        if robust:
            center = np.nanmedian(m, axis=1, keepdims=True)
            scale = 1.482602218505602 * np.nanmedian(np.abs(m - center), axis=1, keepdims=True)
        else:
            center = np.nanmean(m, axis=1, keepdims=True)
            scale = np.nanstd(m, axis=1, keepdims=True)
        scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
        out = (m - center) / scale
        return out[0] if one_d else out

    if wscore_preprocessing != "none":
        robust = wscore_preprocessing == "spatial_robust_z"
        flat_ref = reference_data.reshape(n_ref, -1)
        reference_data = _spatial(flat_ref, robust).reshape(reference_data.shape)
        patient_data = _spatial(patient_data.reshape(-1), robust).reshape(patient_data.shape)

    # Residualize age/sex per column, then score.
    X_ref, x_pat = _standardize_covariates(demo_ref, demo_pat)
    resid_ref, resid_pat = _residualize(
        reference_data.reshape(n_ref, -1), patient_data.reshape(-1), X_ref, x_pat
    )
    if input_was_multidepth:
        resid_ref = resid_ref.reshape(n_ref, n_vert, n_depths)
        resid_pat = resid_pat.reshape(n_vert, n_depths)

    k = max(1, n_ref // 2)   # HALF the available controls
    scores = knn_abnormality(resid_ref, resid_pat, k)

    smoothing_applied = bool(wscore_surface_smoothing_iterations > 0 and scores.size == 32492)
    if smoothing_applied:
        scores = surface_smoother(scores, wscore_surface_smoothing_iterations)

    # kNN has no per-vertex predictive-variance analogue; the uncertainty filter is
    # a no-op for this scorer (documented; the GP rung owns that behavior).
    if prediction_variance_percentile is not None and verbose:
        print("    NOTE: prediction_variance_percentile is ignored by the kNN scorer.")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    data_array = nib.gifti.GiftiDataArray(
        data=np.asarray(scores, dtype=np.float32), intent="NIFTI_INTENT_NORMAL")
    nib.save(nib.gifti.GiftiImage(darrays=[data_array]), output_file)

    model_info = {
        "method": "knn",
        "score": "loo_calibrated_knn_distance",
        "k_neighbours": int(k),
        "k_rule": "half_of_controls",
        "wscore_covariate_model": wscore_covariate_model,
        "wscore_preprocessing": wscore_preprocessing,
        "blur_depth_model": blur_depth_model if input_was_multidepth else None,
        "intensity_depth_model": intensity_depth_model,
        "n_depths_scored_jointly": int(n_depths),
        "reference_count": int(n_ref),
        "wscore_surface_smoothing_iterations": int(wscore_surface_smoothing_iterations),
        "wscore_surface_smoothing_applied": smoothing_applied,
    }
    with open(output_file.replace(".func.gii", "_model.json"), "w", encoding="utf-8") as s:
        json.dump(model_info, s, indent=2)
    if verbose:
        print(f"    Saved kNN abnormality map (k={k} = half of {n_ref} controls): {output_file}")

    return {
        "wscore_file": output_file,
        "mean_wscore": float(np.nanmean(scores)),
        "std_wscore": float(np.nanstd(scores)),
        "normative_data": None,
        "reference_count": int(n_ref),
        "reference_count_min": int(n_ref),
        "reference_count_max": int(n_ref),
        "vertex_covariate_count": 0,
        "uses_prediction_uncertainty": False,
        "wscore_distribution": "knn_loo_calibrated",
        "wscore_preprocessing": wscore_preprocessing,
        "wscore_covariate_model": wscore_covariate_model,
        "wscore_surface_smoothing_iterations": int(wscore_surface_smoothing_iterations),
        "k_neighbours": int(k),
    }
