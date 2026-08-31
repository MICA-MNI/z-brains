import os
import json
import warnings
from functools import lru_cache
from collections.abc import Mapping

import numpy as np
import nibabel as nib
import pandas as pd
from scipy import stats
from scipy import sparse
from scipy.spatial import cKDTree

from zbrains.gaussian_process import (
    calculate_gaussian_process_csv,
    calculate_gaussian_process_maps,
)


WSCORE_DISTRIBUTIONS = {
    "gaussian",
    "gaussian_mad",
    "gaussian_winsor",
    "student_t",
    "empirical",
    "shash",
}
WSCORE_PREPROCESSING = {"none", "spatial_zscore", "spatial_robust_z"}
BLUR_DEPTH_MODELS = {
    # Legacy reductions kept for backward-compat / older tests (no longer swept by
    # the greedy): "mean" (plain depth-average; no FCD precedent), "mean_slope_rms",
    # "gradient_flattening".
    "mean", "mean_slope_rms", "gradient_flattening",
    # --- The THREE literature-grounded, prior-used blur metrics swept by the
    # greedy. Each is ONE per-vertex value from the four-depth profile
    # [midthickness, white, SWM1(~1mm), SWM2(~2mm)] (+ mid->white distance for the
    # profile slope). Every metric samples down to SWM2 (2mm sub-boundary) -- the
    # empirically most-useful depth -- and together they cover all 4 depths + dist.
    # Each probes a DISTINCT construct of gray-white junction integrity (normalized
    # contrast RATIO vs per-mm profile RATE vs sub-boundary LEVEL step), so the
    # sweep stays minimal and non-overfitting:
    #
    #   gray_white_contrast   = 100*(WM - GM)/(0.5*(WM + GM)),  GM=I_mid,
    #                           WM=mean(I_SWM1, I_SWM2)   [depths: mid, SWM1, SWM2]
    #       Percent gray-white matter contrast (FreeSurfer pctsurfcon; Salat et al.
    #       2009, NeuroImage 48:21-28). Used for FCD in MELD (Spitzer/Adler 2022,
    #       Brain 145:3859) and Thesen et al. 2011 (PLoS ONE 6:e16430). GM at
    #       midthickness (~50% depth, closest proxy for pctsurfcon's 30% sample);
    #       WM anchor extended in-family from pctsurfcon's single 1mm sample to a
    #       robust 1-2mm mean. Self-normalizing -> bias-field robust. FCD-II
    #       blurring collapses |contrast| toward zero.
    #
    #   boundary_gradient     = OLS slope of I vs mm-depth [0,d,d+1,d+2]
    #                           over [mid,white,SWM1,SWM2]  [all 4 depths + dist]
    #       Perpendicular (surface-normal) intensity-profile slope; generalizes the
    #       2-point GM/WM gradient of the MNI/Bernasconi lineage (Antel et al. 2003,
    #       NeuroImage 19:1748-1759; Hong et al. 2014, Neurology 83:48-55; Hong 2017
    #       multisurface profiling) to a least-squares slope over the known-spaced
    #       GM->2mm-WM column. Measures the RATE of the transition per mm; FCD
    #       flattens the profile -> reduced |slope|.
    #
    #   juxtacortical_gradient = 0.5*(I_SWM1 + I_SWM2) - I_white
    #                                                     [depths: white, SWM1, SWM2]
    #       Signed level/step of the 1-2mm sub-boundary WM band vs the white
    #       boundary. Intensity-domain surface reduction MOTIVATED BY juxtacortical
    #       WM signal abnormality / transmantle sign (Barkovich et al. 1997; ILAE
    #       Blumcke et al. 2011, Epilepsia 52:158-174; extension image Huppertz
    #       et al. 2005, Epilepsy Res 67:35-50), built on MELD's fixed sub-boundary
    #       sampling (Adler et al. 2017, NeuroImage:Clin 14:18-27) with a 1-2mm
    #       band-mean as robust deep-WM aggregation. Signal staying gray-like 1-2mm
    #       into WM (large step) flags gray-into-white extension.
    "gray_white_contrast",
    "boundary_gradient",
    "juxtacortical_gradient",
}
INTENSITY_DEPTH_MODELS = {
    "raw",
    "white_swm1_direction_cosine",
    "multisurface_median_abs_dominant",
    # Self-normalized white-surface INTENSITY models (feature meaning preserved;
    # apples-to-apples with WhiteStripe, which also yields a normalized intensity):
    "white_swm_zscore",          # (white - median(SWM)) / (1.4826*MAD(SWM))  [WhiteStripe analogue, median]
    "white_swm_mode_zscore",     # (white - mode(SWM))   / (1.4826*MAD(SWM))  [WhiteStripe-faithful, mode]
    "white_swm2_zscore",         # (white - median(SWM2))/ (1.4826*MAD(SWM2)) [deep-WM reference]
    "white_swm_ratio",           # white / median(SWM)
    "white_gmwm_contrast",       # (white - WM_ref) / (GM_ref - WM_ref)       [within-subject GM->WM unit]
    "white_percentile_scale",    # (white - p50) / (p95 - p5)
    "white_surface_robust_z",    # (white - median(white)) / (1.4826*MAD(white))
    # Plain intensity SAMPLERS (no blur/contrast, just a level): two staged axes
    # -- sampling SURFACE (midthickness/white/swm1 = 4-depth columns 0/1/2) x
    # within-subject SELF-NORM (none / swm2-anchor / own-cortex distribution).
    "sample_midthickness", "sample_white", "sample_swm1",
    "sample_midthickness_swm2", "sample_white_swm2", "sample_swm1_swm2",
    "sample_midthickness_owncortex", "sample_white_owncortex", "sample_swm1_owncortex",
}
GAUSSIAN_PROCESS_KERNELS = {
    "gaussian_process": "rbf_ard",
    "gaussian_process_rbf_isotropic": "rbf_isotropic",
    "gaussian_process_matern32": "matern32_ard",
    "gaussian_process_matern52": "matern52_ard",
    "gaussian_process_rational_quadratic": "rational_quadratic_ard",
}
WSCORE_COVARIATE_MODELS = {
    "linear",
    "quadratic_age",
    "age_sex_interaction",
    "quadratic_age_sex_interaction",
    "knn",   # nonparametric kNN normative scorer (k = half the controls)
    *GAUSSIAN_PROCESS_KERNELS,
}
_SHASH_QUANTILE_PROBABILITIES = np.array([0.10, 0.25, 0.50, 0.75, 0.90])


def _normalize_wscore_distribution(distribution):
    """Validate and normalize a W-score residual distribution name."""
    normalized = str(distribution).strip().lower().replace("-", "_")
    aliases = {
        "normal": "gaussian",
        "mad": "gaussian_mad",
        "robust": "gaussian_mad",
        "winsor": "gaussian_winsor",
        "winsorized": "gaussian_winsor",
        "t": "student_t",
        "rank": "empirical",
        "rank_normal": "empirical",
        "sinh_arcsinh": "shash",
        "sinh_arcsinh_normal": "shash",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in WSCORE_DISTRIBUTIONS:
        choices = ", ".join(sorted(WSCORE_DISTRIBUTIONS))
        raise ValueError(
            f"Unsupported W-score distribution '{distribution}'. Choose one of: {choices}."
        )
    return normalized


def _normalize_wscore_preprocessing(preprocessing):
    """Validate preprocessing applied independently to every subject map."""
    normalized = str(preprocessing).strip().lower().replace("-", "_")
    aliases = {
        "off": "none",
        "raw": "none",
        "zscore": "spatial_zscore",
        "spatial_standard_z": "spatial_zscore",
        "spatial_z": "spatial_robust_z",
        "robust_zscore": "spatial_robust_z",
        "spatial_robust_zscore": "spatial_robust_z",
        "robust_spatial_z": "spatial_robust_z",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in WSCORE_PREPROCESSING:
        choices = ", ".join(sorted(WSCORE_PREPROCESSING))
        raise ValueError(
            f"Unsupported W-score preprocessing '{preprocessing}'. Choose one of: {choices}."
        )
    return normalized


def _normalize_blur_depth_model(model):
    """Validate how multi-depth blur maps are reduced to one surface map."""
    normalized = str(model).strip().lower().replace("-", "_")
    aliases = {
        "average": "mean",
        "legacy": "mean",
        "rms": "mean_slope_rms",
        "mean_slope": "mean_slope_rms",
        "flattening": "gradient_flattening",
        "gradient_flat": "gradient_flattening",
        "gw_swm1_gradient_flattening": "gradient_flattening",
        "boundary": "boundary_gradient",
        "gw_gradient": "boundary_gradient",
        "gray_white_gradient": "boundary_gradient",
        "juxtacortical": "juxtacortical_gradient",
        "juxta": "juxtacortical_gradient",
        "swm1_gradient": "juxtacortical_gradient",
        "gwc": "gray_white_contrast",
        "grey_white_contrast": "gray_white_contrast",
        "gray_white_pct": "gray_white_contrast",
        "pctsurfcon": "gray_white_contrast",
        "salat_contrast": "gray_white_contrast",
        "contrast": "gray_white_contrast",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in BLUR_DEPTH_MODELS:
        choices = ", ".join(sorted(BLUR_DEPTH_MODELS))
        raise ValueError(
            f"Unsupported blur depth model '{model}'. Choose one of: {choices}."
        )
    return normalized


def _blur_depth_mean_and_slope(values):
    """Return depth mean and equal-spaced OLS slope along the last axis."""
    values = np.asarray(values, dtype=float)
    if values.ndim < 2 or values.shape[-1] < 2:
        raise ValueError(
            "Multi-depth blur scoring requires at least two depth measurements"
        )
    depth = np.arange(values.shape[-1], dtype=float)
    depth -= np.mean(depth)
    slope = np.tensordot(values, depth / np.sum(depth**2), axes=([-1], [0]))
    return np.mean(values, axis=-1), slope


def _blur_gradient_flattening_components(values, gray_white_distance):
    """Return distance-corrected boundary and white-to-SWM1 gradients.

    ``values`` must be ordered [midthickness, white, SWM1, SWM2]. The first
    component is |white - midthickness| divided by the physical surface
    separation. The second is |SWM1 - white| over the fixed 1-mm interval.
    """
    values = np.asarray(values, dtype=float)
    distance = np.asarray(gray_white_distance, dtype=float)
    if values.ndim != 2 or values.shape[1] != 4:
        raise ValueError(
            "gradient_flattening requires four depth arrays ordered as "
            "midthickness, white, SWM1, SWM2"
        )
    if distance.shape != (values.shape[0],):
        raise ValueError(
            "gray-white distance must contain one value per surface vertex"
        )
    finite_positive = distance[np.isfinite(distance) & (distance > 0)]
    fallback = float(np.median(finite_positive)) if finite_positive.size else 1.0
    distance = np.where(
        np.isfinite(distance) & (distance > 0), distance, fallback
    )
    distance_floor = (
        float(np.percentile(finite_positive, 1.0))
        if finite_positive.size
        else np.finfo(float).eps
    )
    safe_distance = np.maximum(distance, max(distance_floor, np.finfo(float).eps))
    boundary_gradient = np.abs(values[:, 1] - values[:, 0]) / safe_distance
    swm1_gradient = np.abs(values[:, 2] - values[:, 1])
    return np.column_stack([boundary_gradient, swm1_gradient])


def _load_gray_white_surface_distance(
    micapipe_directory,
    participant_id,
    session_id,
    hemi,
    resolution="32k",
):
    """Load physical midthickness-to-white distance on an fsLR surface."""
    surface_dir = os.path.join(
        micapipe_directory, participant_id, session_id, "surf"
    )
    prefix = (
        f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_"
        f"surf-fsLR-{resolution}_label-"
    )
    middle_file = os.path.join(surface_dir, f"{prefix}midthickness.surf.gii")
    white_file = os.path.join(surface_dir, f"{prefix}white.surf.gii")
    if not os.path.exists(middle_file) or not os.path.exists(white_file):
        raise FileNotFoundError(
            "Missing midthickness/white surfaces required by gradient_flattening: "
            f"{participant_id}/{session_id}/hemi-{hemi}"
        )
    middle = np.asarray(nib.load(middle_file).darrays[0].data, dtype=float)
    white = np.asarray(nib.load(white_file).darrays[0].data, dtype=float)
    if middle.shape != white.shape or middle.ndim != 2 or middle.shape[1] != 3:
        raise ValueError(
            "Midthickness and white surfaces must have matching (vertices, 3) coordinates"
        )
    return np.sqrt(np.sum((middle - white) ** 2, axis=1))


def _load_blur_gradient_flattening_components(
    values,
    micapipe_directory,
    participant_id,
    session_id,
    hemi,
    resolution="32k",
):
    distance = _load_gray_white_surface_distance(
        micapipe_directory,
        participant_id,
        session_id,
        hemi,
        resolution=resolution,
    )
    return _blur_gradient_flattening_components(values, distance)


def _load_bilateral_blur_gradient_flattening_components(
    left,
    right,
    micapipe_directory,
    participant_id,
    session_id,
    resolution="32k",
):
    """Self-normalize both hemispheres, then calculate blur gradients."""
    left_distance = _load_gray_white_surface_distance(
        micapipe_directory,
        participant_id,
        session_id,
        "L",
        resolution=resolution,
    )
    right_distance = _load_gray_white_surface_distance(
        micapipe_directory,
        participant_id,
        session_id,
        "R",
        resolution=resolution,
    )
    return _bilateral_blur_gradient_flattening_components(
        left, right, left_distance, right_distance
    )


def _bilateral_blur_gradient_flattening_components(
    left,
    right,
    left_gray_white_distance,
    right_gray_white_distance,
):
    """Apply bilateral median-absolute scaling before gradient extraction."""
    scaled_left, scaled_right = _bilateral_multisurface_median_abs_scale(
        left, right
    )
    left_components = _blur_gradient_flattening_components(
        scaled_left,
        left_gray_white_distance,
    )
    right_components = _blur_gradient_flattening_components(
        scaled_right,
        right_gray_white_distance,
    )
    return left_components, right_components


def _load_blur_single_gradient(values, micapipe_directory, participant_id,
                               session_id, hemi, model, resolution="32k"):
    """One per-vertex value for the three literature-grounded blur metrics.

    ``values`` is the (vertices, 4) raw four-depth profile ordered
    [midthickness, white, SWM1(~1mm), SWM2(~2mm)]. Returns a 1-D (vertices,) array
    so it feeds the normative model as an ordinary single-value feature (never the
    multi-component path). All three are SIGNED; the control w-score resolves the
    per-modality sign (WM-bright T1w vs WM-dark FLAIR/qT1). Every metric now samples
    down to SWM2 (2mm sub-boundary), which is empirically the most useful depth for
    FCD; the deeper WM is also less partial-volume-contaminated than a 1mm-only
    sample.

    gray_white_contrast   [depths: mid, SWM1, SWM2]
        Percent gray-white matter contrast (FreeSurfer ``pctsurfcon``;
        Salat et al. 2009, NeuroImage 48:21-28): ``100*(WM - GM)/(0.5*(WM + GM))``
        with GM = midthickness (~50% depth; closest available proxy for FreeSurfer's
        30%-depth GM sample) and WM = mean(SWM1, SWM2) -- an in-family extension of
        pctsurfcon's single 1mm anchor to a robust 1-2mm deep-WM anchor (pctsurfcon
        itself samples 1mm only; MELD (Spitzer/Adler 2022) samples multiple
        sub-boundary WM depths, supporting the deeper anchor). Used for FCD in MELD
        and Thesen 2011; FCD-II junction blurring collapses |contrast| toward zero.
    boundary_gradient     [depths: mid, white, SWM1, SWM2 + dist]
        Perpendicular (surface-normal) intensity-PROFILE SLOPE: the OLS slope of
        intensity vs physical mm-depth over all four samples at depths
        [0, d, d+1, d+2] with d = dist(mid->white) (white->SWM1->SWM2 are fixed ~1mm
        steps). Generalizes the 2-point boundary gradient of the MNI/Bernasconi
        lineage (Antel et al. 2003 gradient-magnitude; Hong et al. 2014 perpendicular
        intensity gradient; Hong et al. 2017 multisurface intra-/sub-cortical
        profiling) to a least-squares slope over the known-spaced GM->2mm-WM column.
        |slope| large = sharp junction; toward zero = flattened profile (blurring).
    juxtacortical_gradient [depths: white, SWM1, SWM2]
        ``0.5*(I_SWM1 + I_SWM2) - I_white`` -- signed level/step of the 1-2mm
        sub-boundary WM band relative to the white boundary. Intensity-domain surface
        reduction motivated by juxtacortical WM signal abnormality / the transmantle
        sign (Barkovich et al. 1997; ILAE classification Blumcke et al. 2011; Huppertz
        2005 extension image), built on MELD's fixed sub-boundary offset sampling
        (Adler 2017) with the 1-2mm band-mean as a robust deep-WM aggregation. Signal
        that stays gray-like 1-2mm into WM (large step) flags gray-into-white
        extension. Fixed ~1mm offsets, so no distance division.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] < 4:
        raise ValueError(
            "blur metric requires the four-depth profile [mid, white, SWM1, SWM2]"
        )
    if model == "gray_white_contrast":
        gm = values[:, 0]                        # midthickness ~ GM
        wm = 0.5 * (values[:, 2] + values[:, 3])  # mean(SWM1, SWM2): 1-2mm deep-WM anchor
        denom = 0.5 * (wm + gm)
        # Floor the denominator away from zero (subject-normalized images can have
        # near-zero local means) without flipping its sign.
        safe_denom = np.where(np.abs(denom) < 1e-6,
                              np.where(denom < 0, -1e-6, 1e-6), denom)
        return 100.0 * (wm - gm) / safe_denom
    if model == "boundary_gradient":
        # OLS slope of intensity vs physical mm-depth x=[0, d, d+1, d+2] over the
        # four samples [mid, white, SWM1, SWM2]. Closed form (verified equal to
        # np.polyfit to ~5e-13): num/den with den = 3d^2+6d+11 >= 11 > 0, so no
        # distance floor is needed (d = per-vertex mid->white separation >= 0).
        distance = _load_gray_white_surface_distance(
            micapipe_directory, participant_id, session_id, hemi, resolution=resolution
        )
        d = np.asarray(distance, dtype=float)
        d = np.where(np.isfinite(d) & (d >= 0.0), d, 0.0)
        num = (4.0 * (d * values[:, 1] + (d + 1.0) * values[:, 2] + (d + 2.0) * values[:, 3])
               - (3.0 * d + 3.0) * (values[:, 0] + values[:, 1] + values[:, 2] + values[:, 3]))
        den = 3.0 * d * d + 6.0 * d + 11.0
        return num / den
    if model == "juxtacortical_gradient":
        return 0.5 * (values[:, 2] + values[:, 3]) - values[:, 1]
    raise ValueError(f"unsupported single-gradient blur model '{model}'")


def _normalize_intensity_depth_model(model):
    """Validate locally self-normalized T1w/FLAIR intensity analysis."""
    normalized = str(model).strip().lower().replace("-", "_")
    aliases = {
        "none": "raw",
        "off": "raw",
        "white_swm1": "white_swm1_direction_cosine",
        "direction_cosine": "white_swm1_direction_cosine",
        "depth_cosine": "white_swm1_direction_cosine",
        "multisurface": "multisurface_median_abs_dominant",
        "multisurface_intensity": "multisurface_median_abs_dominant",
        "median_abs_dominant": "multisurface_median_abs_dominant",
        "swm_zscore": "white_swm_zscore",
        "wm_zscore": "white_swm_zscore",
        "whitestripe_self": "white_swm_zscore",
        "swm_mode_zscore": "white_swm_mode_zscore",
        "wm_mode_zscore": "white_swm_mode_zscore",
        "swm2_zscore": "white_swm2_zscore",
        "deep_wm_zscore": "white_swm2_zscore",
        "swm_ratio": "white_swm_ratio",
        "wm_ratio": "white_swm_ratio",
        "gmwm_contrast": "white_gmwm_contrast",
        "gm_wm_contrast": "white_gmwm_contrast",
        "percentile_scale": "white_percentile_scale",
        "white_robust_z": "white_surface_robust_z",
        "robust_intensity": "white_surface_robust_z",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in INTENSITY_DEPTH_MODELS:
        choices = ", ".join(sorted(INTENSITY_DEPTH_MODELS))
        raise ValueError(
            f"Unsupported intensity depth model '{model}'. Choose one of: {choices}."
        )
    return normalized


def _white_swm1_direction_cosine(values, floor_percentile=10.0):
    """Return stabilized white-minus-SWM1 contrast over four-depth energy.

    The final axis must contain [midthickness, white, SWM1, SWM2]. The local
    energy denominator is floored at its within-map percentile to avoid noisy
    ratios in vertices where every depth is close to zero.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim < 2 or values.shape[-1] != 4:
        raise ValueError(
            "white_swm1_direction_cosine requires exactly four depth arrays "
            "ordered as midthickness, white, SWM1, SWM2"
        )
    denominator = np.sqrt(np.sum(values**2, axis=-1))
    finite_positive = denominator[np.isfinite(denominator) & (denominator > 0)]
    floor = (
        float(np.percentile(finite_positive, floor_percentile))
        if finite_positive.size
        else np.finfo(float).eps
    )
    denominator = np.maximum(denominator, max(floor, np.finfo(float).eps))
    return (values[..., 1] - values[..., 2]) / denominator


def _bilateral_multisurface_median_abs_scale(left, right):
    """Divide both four-depth hemispheres by one robust subject-level scale."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.ndim != 2 or right.ndim != 2 or left.shape[1] != 4 or right.shape[1] != 4:
        raise ValueError(
            "multisurface_median_abs_dominant requires four depth arrays in both hemispheres"
        )
    combined = np.concatenate([left.ravel(), right.ravel()])
    combined = combined[np.isfinite(combined)]
    scale = np.median(np.abs(combined)) if combined.size else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return left / scale, right / scale


def _superficial_wm_reference(values):
    """Robust (center, scale) of a subject's own superficial-WM samples.

    ``values`` is the per-hemisphere four-depth profile ordered as
    [midthickness, white, SWM1, SWM2]. The SWM1/SWM2 columns lie inside white
    matter, so their pooled distribution is the surface-sampled analogue of the
    WhiteStripe "white stripe": a within-subject white-matter intensity landmark
    estimated from the subject's own data, with no atlas, SynthSeg, or histogram
    mode required. Returns a robust center (median) and scale (1.4826*MAD, with a
    std fallback).
    """
    values = np.asarray(values, dtype=float)
    swm = np.concatenate([values[..., 2].ravel(), values[..., 3].ravel()])
    swm = swm[np.isfinite(swm)]
    if swm.size == 0:
        return 0.0, 1.0
    center = float(np.median(swm))
    scale = 1.482602218505602 * float(np.median(np.abs(swm - center)))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = float(np.std(swm)) if swm.size > 1 else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return center, scale


def _require_four_depths(values, name):
    values = np.asarray(values, dtype=float)
    if values.ndim < 2 or values.shape[-1] != 4:
        raise ValueError(
            f"{name} requires four depth arrays ordered as "
            "midthickness, white, SWM1, SWM2"
        )
    return values


def _white_swm_zscore(values):
    """Self-normalized white-surface INTENSITY (surface WhiteStripe analogue).

    Standardizes the white-surface value by the subject's own superficial-WM
    center/scale: ``(white - median(SWM)) / (1.4826 * MAD(SWM))``. This mirrors
    WhiteStripe's within-subject WM standardization but is computed per subject
    from the surface-sampled profile. The output is still an intensity (WM-
    standardized brightness of the white surface), so the feature's meaning is
    unchanged -- it is directly comparable to WhiteStripe, not a contrast.
    """
    values = _require_four_depths(values, "white_swm_zscore")
    center, scale = _superficial_wm_reference(values)
    return (values[..., 1] - center) / scale


def _white_swm_ratio(values):
    """Self-normalized white-surface INTENSITY as a ratio to the subject's WM.

    ``white / median(SWM)`` -- puts each subject's superficial WM at ~1.0,
    cancelling global per-subject/scanner intensity scale while preserving the
    intensity interpretation (relative brightness). A minimal, robust,
    within-subject WM normalization.
    """
    values = _require_four_depths(values, "white_swm_ratio")
    swm = np.concatenate([values[..., 2].ravel(), values[..., 3].ravel()])
    swm = swm[np.isfinite(swm)]
    reference = float(np.median(swm)) if swm.size else 1.0
    if not np.isfinite(reference) or abs(reference) <= np.finfo(float).eps:
        reference = 1.0
    return values[..., 1] / reference


def _white_surface_robust_z(values):
    """Self-normalized white-surface INTENSITY by its own robust distribution.

    ``(white - median(white)) / (1.4826 * MAD(white))`` across the subject's
    cortex. A landmark-free, within-subject standardization of the intensity map
    that removes global brightness/contrast differences; the output remains an
    intensity. Uses the white column of the four-depth profile.
    """
    values = _require_four_depths(values, "white_surface_robust_z")
    return _robust_standardize(values[..., 1])


def _robust_standardize(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values
    center = float(np.median(finite))
    scale = 1.482602218505602 * float(np.median(np.abs(finite - center)))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = float(np.std(finite)) if finite.size > 1 else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return (values - center) / scale


def _histogram_mode(values, bins=256):
    """Robust histogram mode of a 1-D sample (trimmed to the 1-99 percentile)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    lo, hi = np.percentile(values, [1.0, 99.0])
    trimmed = values[(values >= lo) & (values <= hi)]
    if trimmed.size < 2 or not np.isfinite(hi - lo) or hi <= lo:
        return float(np.median(values))
    histogram, edges = np.histogram(trimmed, bins=bins)
    peak = int(np.argmax(histogram))
    return float((edges[peak] + edges[peak + 1]) / 2.0)


def _white_swm_mode_zscore(values):
    """Self-normalized white-surface INTENSITY, WhiteStripe-FAITHFUL variant.

    Uses the histogram MODE of the subject's superficial-WM samples as the
    center (WhiteStripe centers on the WM intensity mode, not the median),
    scaled by the robust WM spread: ``(white - mode(SWM)) / (1.4826*MAD(SWM))``.
    This is the closest surface-sampled analogue of the actual WhiteStripe
    estimator. Output is an intensity.
    """
    values = _require_four_depths(values, "white_swm_mode_zscore")
    swm = np.concatenate([values[..., 2].ravel(), values[..., 3].ravel()])
    swm = swm[np.isfinite(swm)]
    if swm.size == 0:
        return values[..., 1]
    center = _histogram_mode(swm)
    scale = 1.482602218505602 * float(np.median(np.abs(swm - np.median(swm))))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = float(np.std(swm)) if swm.size > 1 else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return (values[..., 1] - center) / scale


def _white_swm2_zscore(values):
    """Self-normalized white-surface INTENSITY vs the DEEPEST WM (SWM2 only).

    Uses only the 2 mm superficial-WM depth as the reference, which is the most
    lesion-free / least partial-volumed internal tissue:
    ``(white - median(SWM2)) / (1.4826*MAD(SWM2))``. Output is an intensity.
    """
    values = _require_four_depths(values, "white_swm2_zscore")
    swm2 = values[..., 3]
    finite = swm2[np.isfinite(swm2)]
    if finite.size == 0:
        return values[..., 1]
    center = float(np.median(finite))
    scale = 1.482602218505602 * float(np.median(np.abs(finite - center)))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = float(np.std(finite)) if finite.size > 1 else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return (values[..., 1] - center) / scale


def _white_gmwm_contrast(values):
    """Self-normalized white-surface INTENSITY in within-subject GM->WM units.

    Places the white-surface value on the subject's own gray-to-white intensity
    axis: ``(white - WM_ref) / (GM_ref - WM_ref)`` with GM_ref = median of the
    midthickness (gray) samples and WM_ref = median of the superficial-WM
    samples. Result is ~0 in WM and ~1 in cortical gray -- a biologically
    anchored, subject-invariant intensity unit (still an intensity, not a
    depth contrast).
    """
    values = _require_four_depths(values, "white_gmwm_contrast")
    gm = values[..., 0]
    swm = np.concatenate([values[..., 2].ravel(), values[..., 3].ravel()])
    gm_ref = float(np.median(gm[np.isfinite(gm)])) if np.isfinite(gm).any() else 0.0
    wm_ref = float(np.median(swm[np.isfinite(swm)])) if np.isfinite(swm).any() else 0.0
    denominator = gm_ref - wm_ref
    if not np.isfinite(denominator) or abs(denominator) <= np.finfo(float).eps:
        denominator = 1.0
    return (values[..., 1] - wm_ref) / denominator


def _white_percentile_scale(values):
    """Self-normalized white-surface INTENSITY by its own robust range.

    ``(white - p50(white)) / (p95(white) - p5(white))`` -- a robust min-max style
    scaling using the subject's own white-surface intensity distribution. Output
    is an intensity.
    """
    values = _require_four_depths(values, "white_percentile_scale")
    white = values[..., 1]
    finite = white[np.isfinite(white)]
    if finite.size == 0:
        return white
    p5, p50, p95 = np.percentile(finite, [5.0, 50.0, 95.0])
    spread = float(p95 - p5)
    if not np.isfinite(spread) or spread <= np.finfo(float).eps:
        spread = 1.0
    return (white - float(p50)) / spread


# ---------------------------------------------------------------------------
# Plain intensity SAMPLERS (no depth-difference / blur). Two orthogonal choices:
#   * sampling surface -> column of the 4-depth profile [midthickness, white,
#     SWM1, SWM2]: midthickness=0 (cortical GM), white=1 (GM/WM boundary),
#     swm1=2 (superficial WM ~1mm);
#   * within-subject self-normalization reference: none (plain level), swm2
#     (anchor to the ~2mm deep-WM column), or owncortex (the map's own robust
#     median/MAD across cortex).
# The output is a plain per-vertex INTENSITY level (hyper/hypo); the downstream
# normative model does the cross-subject z/w-scoring.
# ---------------------------------------------------------------------------
_SAMPLE_DEPTH_COLUMN = {"midthickness": 0, "white": 1, "swm1": 2}


def _swm2_reference(values):
    """Robust center/scale from the deep superficial-WM (SWM2) column."""
    swm2 = np.asarray(values[..., 3], dtype=float)
    finite = swm2[np.isfinite(swm2)]
    if finite.size == 0:
        return 0.0, 1.0
    center = float(np.median(finite))
    scale = 1.482602218505602 * float(np.median(np.abs(finite - center)))
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = float(np.std(finite)) if finite.size > 1 else 1.0
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = 1.0
    return center, scale


def _make_surface_sampler(column, self_norm):
    """Return a 4-depth transform: sample one surface, optionally self-normalized."""
    def _transform(values):
        values = _require_four_depths(values, "surface sampler")
        x = values[..., column]
        if self_norm is None:
            return x
        if self_norm == "swm2":
            center, scale = _swm2_reference(values)
            return (x - center) / scale
        if self_norm == "owncortex":
            return _robust_standardize(x)
        raise ValueError(f"unknown self_norm '{self_norm}'")
    return _transform


# Name -> per-hemisphere 4-depth transform. multisurface_median_abs_dominant is
# handled separately (it is a bilateral transform), so it is not in this map.
_INTENSITY_DEPTH_TRANSFORMS = {
    "white_swm1_direction_cosine": _white_swm1_direction_cosine,
    "white_swm_zscore": _white_swm_zscore,
    "white_swm_mode_zscore": _white_swm_mode_zscore,
    "white_swm2_zscore": _white_swm2_zscore,
    "white_swm_ratio": _white_swm_ratio,
    "white_gmwm_contrast": _white_gmwm_contrast,
    "white_percentile_scale": _white_percentile_scale,
    "white_surface_robust_z": _white_surface_robust_z,
}
for _surf, _col in _SAMPLE_DEPTH_COLUMN.items():
    _INTENSITY_DEPTH_TRANSFORMS[f"sample_{_surf}"] = _make_surface_sampler(_col, None)
    _INTENSITY_DEPTH_TRANSFORMS[f"sample_{_surf}_swm2"] = _make_surface_sampler(_col, "swm2")
    _INTENSITY_DEPTH_TRANSFORMS[f"sample_{_surf}_owncortex"] = _make_surface_sampler(_col, "owncortex")


def _normalize_wscore_covariate_model(model):
    normalized = str(model).strip().lower().replace("-", "_")
    aliases = {
        "quadratic": "quadratic_age",
        "age2": "quadratic_age",
        "interaction": "age_sex_interaction",
        "quadratic_interaction": "quadratic_age_sex_interaction",
        "gp": "gaussian_process",
        "gpr": "gaussian_process",
        "gp_isotropic": "gaussian_process_rbf_isotropic",
        "gp_rbf_isotropic": "gaussian_process_rbf_isotropic",
        "gp_matern32": "gaussian_process_matern32",
        "gp_matern52": "gaussian_process_matern52",
        "gp_rq": "gaussian_process_rational_quadratic",
        "gp_rational_quadratic": "gaussian_process_rational_quadratic",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in WSCORE_COVARIATE_MODELS:
        choices = ", ".join(sorted(WSCORE_COVARIATE_MODELS))
        raise ValueError(
            f"Unsupported W-score covariate model '{model}'. Choose one of: {choices}."
        )
    return normalized


def _augment_wscore_covariates(X_ref, X_pat, normative_columns, model):
    """Add control-scaled age curvature and/or age-by-sex terms."""
    model = _normalize_wscore_covariate_model(model)
    X_ref = np.asarray(X_ref, dtype=float)
    X_pat = np.asarray(X_pat, dtype=float).reshape(-1)
    names = list(normative_columns)
    if model == "linear" or model in GAUSSIAN_PROCESS_KERNELS:
        return X_ref, X_pat, names

    lower_names = [str(name).lower() for name in names]
    if "age" not in lower_names:
        raise ValueError(f"Covariate model '{model}' requires an AGE column")
    age_index = lower_names.index("age")
    age_mean = float(np.mean(X_ref[:, age_index]))
    age_sd = float(np.std(X_ref[:, age_index]))
    if not np.isfinite(age_sd) or age_sd <= np.finfo(float).eps:
        raise ValueError("AGE has zero variance; cannot fit nonlinear age terms")
    reference_age_z = (X_ref[:, age_index] - age_mean) / age_sd
    patient_age_z = (X_pat[age_index] - age_mean) / age_sd

    ref_terms = [X_ref]
    pat_terms = [X_pat]
    if model in {"quadratic_age", "quadratic_age_sex_interaction"}:
        ref_terms.append(reference_age_z[:, None] ** 2)
        pat_terms.append(np.array([patient_age_z**2]))
        names.append("AGE_Z2")
    if model in {"age_sex_interaction", "quadratic_age_sex_interaction"}:
        if "sex" not in lower_names:
            raise ValueError(f"Covariate model '{model}' requires a SEX column")
        sex_index = lower_names.index("sex")
        ref_terms.append((reference_age_z * X_ref[:, sex_index])[:, None])
        pat_terms.append(np.array([patient_age_z * X_pat[sex_index]]))
        names.append("AGE_Z_X_SEX")
    return np.column_stack(ref_terms), np.concatenate(pat_terms), names


def _spatial_robust_z(values):
    """Median/MAD normalize each subject across vertices without data leakage."""
    values = np.asarray(values, dtype=float)
    one_dimensional = values.ndim == 1
    matrix = values[None, :] if one_dimensional else values
    center = np.nanmedian(matrix, axis=1, keepdims=True)
    scale = 1.482602218505602 * np.nanmedian(
        np.abs(matrix - center), axis=1, keepdims=True
    )
    valid_scale = np.isfinite(scale) & (scale > np.finfo(float).eps)
    scale = np.where(valid_scale, scale, 1.0)
    normalized = (matrix - center) / scale
    return normalized[0] if one_dimensional else normalized


def _spatial_zscore(values):
    """Mean/SD normalize every subject independently across map locations."""
    values = np.asarray(values, dtype=float)
    one_dimensional = values.ndim == 1
    matrix = values[None, :] if one_dimensional else values
    center = np.nanmean(matrix, axis=1, keepdims=True)
    scale = np.nanstd(matrix, axis=1, keepdims=True)
    valid_scale = np.isfinite(scale) & (scale > np.finfo(float).eps)
    scale = np.where(valid_scale, scale, 1.0)
    normalized = (matrix - center) / scale
    return normalized[0] if one_dimensional else normalized


@lru_cache(maxsize=1)
def _fslr32k_averaging_operator():
    """Build a one-ring fsLR-32k averaging operator including each vertex."""
    surface_file = os.path.join(os.path.dirname(__file__), "data", "fsLR-32k.L.surf.gii")
    surface = nib.load(surface_file)
    faces = next(
        np.asarray(array.data, dtype=np.int64)
        for array in surface.darrays
        if int(array.intent) == 1009
    )
    n_vertices = 32492
    rows = np.concatenate(
        [faces[:, 0], faces[:, 1], faces[:, 2], faces[:, 1], faces[:, 2], faces[:, 0]]
    )
    cols = np.concatenate(
        [faces[:, 1], faces[:, 2], faces[:, 0], faces[:, 0], faces[:, 1], faces[:, 2]]
    )
    adjacency = sparse.coo_matrix(
        (np.ones(rows.size), (rows, cols)), shape=(n_vertices, n_vertices)
    ).tocsr()
    adjacency.data[:] = 1.0
    adjacency = adjacency + sparse.eye(n_vertices, format="csr")
    return sparse.diags(1.0 / np.asarray(adjacency.sum(axis=1)).ravel()) @ adjacency


def _smooth_fslr32k_scores(values, iterations):
    """Diffuse scores over fsLR neighbours while preserving missing vertices."""
    result = np.asarray(values, dtype=float).copy()
    averaging = _fslr32k_averaging_operator()
    for _ in range(int(iterations)):
        finite = np.isfinite(result)
        numerator = averaging @ np.where(finite, result, 0.0)
        denominator = averaging @ finite.astype(float)
        result = np.divide(
            numerator,
            denominator,
            out=np.full_like(result, np.nan),
            where=denominator > 0,
        )
    return result


@lru_cache(maxsize=1)
def _shash_shape_lookup():
    """Build a lookup from robust quantile shape statistics to SHASH shape."""
    skew_values = np.linspace(-2.5, 2.5, 101)
    tail_values = np.exp(np.linspace(np.log(0.30), np.log(3.0), 101))
    skew_grid, tail_grid = np.meshgrid(skew_values, tail_values, indexing="ij")

    normal_quantiles = stats.norm.ppf(_SHASH_QUANTILE_PROBABILITIES)
    shash_quantiles = np.sinh(
        (np.arcsinh(normal_quantiles)[None, None, :] + skew_grid[..., None])
        / tail_grid[..., None]
    )
    iqr = shash_quantiles[..., 3] - shash_quantiles[..., 1]
    bowley_skew = (
        shash_quantiles[..., 3]
        + shash_quantiles[..., 1]
        - 2.0 * shash_quantiles[..., 2]
    ) / iqr
    outer_ratio = (shash_quantiles[..., 4] - shash_quantiles[..., 0]) / iqr

    # Fixed scaling gives skewness and tail weight similar influence in the
    # nearest-neighbour quantile match.
    lookup_points = np.column_stack(
        [bowley_skew.ravel() / 0.25, np.log(outer_ratio.ravel()) / 0.50]
    )
    return cKDTree(lookup_points), skew_grid.ravel(), tail_grid.ravel()


def _fit_shash_quantile_parameters(residuals, min_observations=8):
    """Fit SHASH location, scale, skew, and tail weight by quantile matching.

    The parameterization is

        Z = sinh(delta * asinh((Y - location) / scale) - epsilon),

    where Z is standard normal.  Consequently, transformed patient residuals
    are directly comparable to conventional Gaussian W-scores.
    """
    residuals = np.asarray(residuals, dtype=float)
    if residuals.ndim == 1:
        residuals = residuals[:, None]
    if residuals.ndim != 2:
        raise ValueError("SHASH residuals must be a subject-by-measure matrix")

    n_measures = residuals.shape[1]
    parameters = np.full((n_measures, 4), np.nan, dtype=float)
    fit_success = np.zeros(n_measures, dtype=bool)
    counts = np.isfinite(residuals).sum(axis=0)

    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", category=RuntimeWarning)
        quantiles = np.nanpercentile(
            residuals,
            100.0 * _SHASH_QUANTILE_PROBABILITIES,
            axis=0,
        )
        fallback_scale = np.nanstd(residuals, axis=0)

    # Gaussian fallback keeps sparse or degenerate measures scoreable whenever
    # their ordinary residual SD is valid.
    parameters[:, 0] = 0.0
    parameters[:, 1] = fallback_scale
    parameters[:, 2] = 0.0
    parameters[:, 3] = 1.0

    empirical_iqr = quantiles[3] - quantiles[1]
    empirical_outer = quantiles[4] - quantiles[0]
    valid = (
        (counts >= int(min_observations))
        & np.isfinite(quantiles).all(axis=0)
        & np.isfinite(empirical_iqr)
        & (empirical_iqr > np.finfo(float).eps)
        & np.isfinite(empirical_outer)
        & (empirical_outer > empirical_iqr)
    )
    if not np.any(valid):
        return parameters, fit_success

    bowley_skew = (
        quantiles[3, valid] + quantiles[1, valid] - 2.0 * quantiles[2, valid]
    ) / empirical_iqr[valid]
    outer_ratio = empirical_outer[valid] / empirical_iqr[valid]
    query_points = np.column_stack(
        [bowley_skew / 0.25, np.log(outer_ratio) / 0.50]
    )
    lookup_tree, skew_grid, tail_grid = _shash_shape_lookup()
    _, lookup_indices = lookup_tree.query(query_points)
    fitted_skew = skew_grid[lookup_indices]
    fitted_tail = tail_grid[lookup_indices]

    normal_quantiles = stats.norm.ppf(_SHASH_QUANTILE_PROBABILITIES)
    fitted_standard_quantiles = np.sinh(
        (np.arcsinh(normal_quantiles)[None, :] + fitted_skew[:, None])
        / fitted_tail[:, None]
    )
    fitted_standard_iqr = (
        fitted_standard_quantiles[:, 3] - fitted_standard_quantiles[:, 1]
    )
    fitted_scale = empirical_iqr[valid] / fitted_standard_iqr
    fitted_location = quantiles[2, valid] - fitted_scale * fitted_standard_quantiles[:, 2]
    valid_fit = (
        np.isfinite(fitted_location)
        & np.isfinite(fitted_scale)
        & (fitted_scale > np.finfo(float).eps)
    )

    valid_indices = np.flatnonzero(valid)
    fitted_indices = valid_indices[valid_fit]
    parameters[fitted_indices, 0] = fitted_location[valid_fit]
    parameters[fitted_indices, 1] = fitted_scale[valid_fit]
    parameters[fitted_indices, 2] = fitted_skew[valid_fit]
    parameters[fitted_indices, 3] = fitted_tail[valid_fit]
    fit_success[fitted_indices] = True
    return parameters, fit_success


def _shash_to_normal_score(value, location, scale, skew, tail_weight):
    """Map a value from a fitted SHASH distribution to latent normal space."""
    if (
        not np.isfinite(value)
        or not np.isfinite(location)
        or not np.isfinite(scale)
        or scale <= 0
        or not np.isfinite(skew)
        or not np.isfinite(tail_weight)
        or tail_weight <= 0
    ):
        return np.nan
    standardized = (value - location) / scale
    with np.errstate(over="ignore", invalid="ignore"):
        latent = np.sinh(tail_weight * np.arcsinh(standardized) - skew)
    # Approximately the finite range of a double-precision normal quantile.
    return float(np.clip(latent, -8.0, 8.0))


def _residual_standard_error(residuals, design):
    """Return OLS residual standard error plus degrees of freedom metadata."""
    residuals = np.asarray(residuals, dtype=float).reshape(-1)
    design = np.asarray(design, dtype=float)
    if design.ndim != 2 or design.shape[0] != residuals.size:
        return np.nan, 0, 0

    rank = int(np.linalg.matrix_rank(design))
    df = int(design.shape[0] - rank)
    if df <= 0:
        return np.nan, df, rank

    rss = float(np.sum(residuals ** 2))
    if not np.isfinite(rss):
        return np.nan, df, rank
    return float(np.sqrt(rss / df)), df, rank


def _prediction_uncertainty_denominator(residual_standard_error, design, patient_design):
    """Return full prediction-interval denominator and patient leverage."""
    if not np.isfinite(residual_standard_error) or residual_standard_error <= 0:
        return np.nan, np.nan

    design = np.asarray(design, dtype=float)
    patient_design = np.asarray(patient_design, dtype=float).reshape(-1)
    if design.ndim != 2 or design.shape[1] != patient_design.size:
        return np.nan, np.nan

    try:
        xtx_inv = np.linalg.pinv(design.T @ design)
        leverage = float(patient_design @ xtx_inv @ patient_design)
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    if not np.isfinite(leverage):
        return np.nan, np.nan
    leverage = max(leverage, 0.0)
    return float(residual_standard_error * np.sqrt(1.0 + leverage)), leverage


def calculate_wscore_maps(
    reference_data,
    patient_data,
    demographics_ref,
    demographics_pat,
    output_file,
    normative_columns=['age', 'sex'],
    verbose=True,
    min_reference_subjects=None,
    reference_vertex_covariates=None,
    patient_vertex_covariates=None,
    use_prediction_uncertainty=False,
    wscore_distribution="gaussian",
    wscore_preprocessing="none",
    wscore_covariate_model="linear",
    wscore_surface_smoothing_iterations=0,
    blur_depth_model="mean_slope_rms",
    intensity_depth_model="raw",
    wscore_fit_cache=None,
    prediction_variance_percentile=None,
):
    """
    Calculate W-scores for patient data against reference data using normative modeling and save as GIFTI file.

    Parameters:
    -----------
    reference_data : np.ndarray
        Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
    patient_data : np.ndarray
        Patient data with shape (n_vertices,) or (n_vertices, n_depths)
    demographics_ref : pd.DataFrame
        Demographics data for reference subjects with normative columns
    demographics_pat : pd.Series or pd.DataFrame
        Demographics data for patient (single row if DataFrame)
    output_file : str
        Path to save the W-score map
    normative_columns : list, default=['age', 'sex']
        List of demographic columns to use for normative modeling
    verbose : bool, default=True
        If True, prints processing information
    use_prediction_uncertainty : bool, default=False
        If True, divide by the full OLS predictive standard deviation for a
        new subject, using residual standard error and patient leverage.
    wscore_distribution : str, default="gaussian"
        Residual scoring model. Supported values are ``gaussian``,
        ``gaussian_mad``, ``gaussian_winsor``, ``student_t``, ``empirical``,
        and ``shash``.
    wscore_preprocessing : {"none", "spatial_zscore", "spatial_robust_z"}, default="none"
        Optional per-subject map normalization before normative regression.
    wscore_covariate_model : str, default="linear"
        Demographic design: linear, quadratic age, age-by-sex interaction, or both.
    wscore_surface_smoothing_iterations : int, default=0
        Number of one-ring diffusion steps applied to fsLR-32k cortical
        W-scores after scoring. Other surface densities are left unchanged.
    blur_depth_model : {"mean_slope_rms", "mean", "gradient_flattening"}, default="mean_slope_rms"
        Reduction used for multi-depth blur data. ``mean_slope_rms`` fits the
        depth mean and equal-spaced depth slope independently, then combines
        their W-scores using a signed root-mean-square magnitude. ``mean``
        preserves the legacy depth-average behavior. ``gradient_flattening``
        first divides all four depths in both hemispheres by one subject-level
        median absolute intensity, independently fits gray-white and white-SWM1
        gradient magnitudes, and returns their negative mean W-score. Positive
        values therefore indicate a flatter-than-normal transition.
    intensity_depth_model : {"raw", "white_swm1_direction_cosine", "multisurface_median_abs_dominant"}, default="raw"
        Provenance label for T1w/FLAIR intensity maps derived from companion
        four-depth data before this normative model is fitted.
    
    Returns:
    --------
    dict
        Dictionary containing W-score statistics, normative data, and file path
    """
    wscore_distribution = _normalize_wscore_distribution(wscore_distribution)
    wscore_preprocessing = _normalize_wscore_preprocessing(wscore_preprocessing)
    wscore_covariate_model = _normalize_wscore_covariate_model(wscore_covariate_model)
    blur_depth_model = _normalize_blur_depth_model(blur_depth_model)
    intensity_depth_model = _normalize_intensity_depth_model(intensity_depth_model)
    wscore_surface_smoothing_iterations = int(wscore_surface_smoothing_iterations)
    if wscore_surface_smoothing_iterations < 0:
        raise ValueError("wscore_surface_smoothing_iterations must be non-negative")
    gaussian_process_kernel = GAUSSIAN_PROCESS_KERNELS.get(
        wscore_covariate_model
    )
    if gaussian_process_kernel is not None:
        if wscore_distribution != "gaussian":
            raise ValueError(
                "Gaussian-process models use their predictive Gaussian "
                "distribution; set wscore_distribution='gaussian'"
            )
        return calculate_gaussian_process_maps(
            reference_data=reference_data,
            patient_data=patient_data,
            demographics_ref=demographics_ref,
            demographics_pat=demographics_pat,
            output_file=output_file,
            normative_columns=normative_columns,
            verbose=verbose,
            gaussian_process_kernel=gaussian_process_kernel,
            wscore_covariate_model=wscore_covariate_model,
            min_reference_subjects=min_reference_subjects,
            reference_vertex_covariates=reference_vertex_covariates,
            patient_vertex_covariates=patient_vertex_covariates,
            wscore_preprocessing=wscore_preprocessing,
            blur_depth_model=blur_depth_model,
            intensity_depth_model=intensity_depth_model,
            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
            wscore_fit_cache=wscore_fit_cache,
            surface_smoother=_smooth_fslr32k_scores,
            prediction_variance_percentile=prediction_variance_percentile,
        )
    if wscore_covariate_model == "knn":
        if wscore_distribution != "gaussian":
            raise ValueError(
                "The kNN normative scorer uses its own leave-one-out calibration; "
                "set wscore_distribution='gaussian'"
            )
        from zbrains.knn_normative import calculate_knn_maps
        return calculate_knn_maps(
            reference_data=reference_data,
            patient_data=patient_data,
            demographics_ref=demographics_ref,
            demographics_pat=demographics_pat,
            output_file=output_file,
            normative_columns=normative_columns,
            verbose=verbose,
            wscore_covariate_model=wscore_covariate_model,
            min_reference_subjects=min_reference_subjects,
            reference_vertex_covariates=reference_vertex_covariates,
            patient_vertex_covariates=patient_vertex_covariates,
            wscore_preprocessing=wscore_preprocessing,
            blur_depth_model=blur_depth_model,
            intensity_depth_model=intensity_depth_model,
            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
            wscore_fit_cache=wscore_fit_cache,
            surface_smoother=_smooth_fslr32k_scores,
            prediction_variance_percentile=prediction_variance_percentile,
        )
    cache_reference_data = reference_data
    cache_reference_covariates = reference_vertex_covariates
    input_was_multidepth = len(reference_data.shape) > 2

    if len(reference_data) == 0:
        raise ValueError("No reference data provided")
    
    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate W-scores")
    
    if demographics_ref.shape[0] != reference_data.shape[0]:
        raise ValueError("Demographics and reference data must have same number of subjects")
    
    # Multi-depth blur maps are represented internally as two adjacent
    # pseudo-surfaces (mean, slope). The usual normative model is then fit to
    # each component independently and their W-scores are combined below.
    blur_component_vertices = None
    blur_component_names = None
    blur_component_combination = None
    intensity_component_vertices = None
    intensity_component_count = 0
    if len(reference_data.shape) > 2:
        if len(patient_data.shape) <= 1:
            raise ValueError(
                "Patient blur data must contain the same depth dimension as controls"
            )
        if reference_data.shape[2] != patient_data.shape[1]:
            raise ValueError(
                "Control and patient blur data must have the same number of depths"
            )
        if intensity_depth_model == "multisurface_median_abs_dominant":
            if reference_data.shape[2] != 4:
                raise ValueError(
                    "multisurface_median_abs_dominant requires exactly four depth arrays"
                )
            intensity_component_vertices = reference_data.shape[1]
            intensity_component_count = reference_data.shape[2]
            reference_data = np.transpose(reference_data, (0, 2, 1)).reshape(
                reference_data.shape[0], -1
            )
            patient_data = np.asarray(patient_data).T.reshape(-1)
            if reference_vertex_covariates is not None:
                reference_vertex_covariates = np.concatenate(
                    [reference_vertex_covariates] * intensity_component_count,
                    axis=1,
                )
                patient_vertex_covariates = np.concatenate(
                    [patient_vertex_covariates] * intensity_component_count,
                    axis=0,
                )
        elif blur_depth_model == "mean":
            reference_data = np.mean(reference_data, axis=2)
            patient_data = np.mean(patient_data, axis=1)
        elif blur_depth_model == "gradient_flattening":
            if reference_data.shape[2] != 2:
                raise ValueError(
                    "gradient_flattening expects two precomputed components: "
                    "gray-white and white-SWM1 gradient magnitudes"
                )
            blur_component_vertices = reference_data.shape[1]
            blur_component_names = [
                "distance_corrected_gray_white_gradient_magnitude",
                "white_to_swm1_gradient_magnitude",
            ]
            blur_component_combination = (
                "-(gray_white_gradient_wscore + white_swm1_gradient_wscore)/2"
            )
            reference_data = np.transpose(reference_data, (0, 2, 1)).reshape(
                reference_data.shape[0], -1
            )
            patient_data = np.asarray(patient_data).T.reshape(-1)
            if reference_vertex_covariates is not None:
                reference_vertex_covariates = np.concatenate(
                    [reference_vertex_covariates, reference_vertex_covariates],
                    axis=1,
                )
                patient_vertex_covariates = np.concatenate(
                    [patient_vertex_covariates, patient_vertex_covariates], axis=0
                )
        else:
            reference_mean, reference_slope = _blur_depth_mean_and_slope(
                reference_data
            )
            patient_mean, patient_slope = _blur_depth_mean_and_slope(patient_data)
            blur_component_vertices = reference_mean.shape[1]
            blur_component_names = ["mean", "equal_spaced_ols_slope"]
            blur_component_combination = (
                "sign(mean_wscore)*sqrt((mean_wscore^2+slope_wscore^2)/2)"
            )
            reference_data = np.concatenate(
                [reference_mean, reference_slope], axis=1
            )
            patient_data = np.concatenate([patient_mean, patient_slope])
            if reference_vertex_covariates is not None:
                reference_vertex_covariates = np.concatenate(
                    [reference_vertex_covariates, reference_vertex_covariates],
                    axis=1,
                )
                patient_vertex_covariates = np.concatenate(
                    [patient_vertex_covariates, patient_vertex_covariates], axis=0
                )
    elif len(patient_data.shape) > 1:
        patient_data = np.mean(patient_data, axis=1)

    preprocess = None
    if wscore_preprocessing == "spatial_zscore":
        preprocess = _spatial_zscore
    elif wscore_preprocessing == "spatial_robust_z":
        preprocess = _spatial_robust_z
    if preprocess is not None:
        if blur_component_vertices is None and intensity_component_vertices is None:
            reference_data = preprocess(reference_data)
            patient_data = preprocess(patient_data)
        else:
            component_vertices = (
                intensity_component_vertices
                if intensity_component_vertices is not None
                else blur_component_vertices
            )
            component_count = (
                intensity_component_count
                if intensity_component_vertices is not None
                else 2
            )
            for component in range(component_count):
                component_slice = slice(
                    component * component_vertices,
                    (component + 1) * component_vertices,
                )
                reference_data[:, component_slice] = preprocess(
                    reference_data[:, component_slice]
                )
                patient_data[component_slice] = preprocess(
                    patient_data[component_slice]
                )
    
    n_vertices = reference_data.shape[1]
    
    if reference_vertex_covariates is not None:
        if patient_vertex_covariates is None:
            raise ValueError("patient_vertex_covariates is required when reference_vertex_covariates is provided")
        reference_vertex_covariates = np.asarray(reference_vertex_covariates, dtype=float)
        if reference_vertex_covariates.ndim == 2:
            reference_vertex_covariates = reference_vertex_covariates[:, :, None]
        if reference_vertex_covariates.shape[:2] != reference_data.shape[:2]:
            raise ValueError(
                "reference_vertex_covariates must have shape "
                f"{reference_data.shape[:2]} or {reference_data.shape[:2]} + (n_covariates,), "
                f"got {reference_vertex_covariates.shape}"
            )
        patient_vertex_covariates = np.asarray(patient_vertex_covariates, dtype=float)
        if patient_vertex_covariates.ndim == 1:
            patient_vertex_covariates = patient_vertex_covariates[:, None]
        if patient_vertex_covariates.shape != reference_vertex_covariates.shape[1:]:
            raise ValueError(
                "patient_vertex_covariates must have shape "
                f"{reference_vertex_covariates.shape[1:]}, got {patient_vertex_covariates.shape}"
            )
    n_vertex_predictors = reference_vertex_covariates.shape[2] if reference_vertex_covariates is not None else 0
    
    # Prepare demographics data
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    
    # Convert to numpy arrays
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    
    # Remove subjects with missing demographic data
    demo_mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[demo_mask]
    ref_data_clean = reference_data[demo_mask]
    if reference_vertex_covariates is not None:
        reference_vertex_covariates = reference_vertex_covariates[demo_mask]
    
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")

    X_ref, X_pat, effective_normative_columns = _augment_wscore_covariates(
        X_ref,
        X_pat,
        normative_columns,
        wscore_covariate_model,
    )
    n_predictors = X_ref.shape[1]

    # Format: [intercept, coef_1, ..., vertex covariates, residual SD].
    normative_data = np.zeros((n_vertices, n_predictors + n_vertex_predictors + 2))
    
    # Add intercept term
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    min_reference_subjects = (
        int(min_reference_subjects)
        if min_reference_subjects is not None
        else max(X_ref_with_intercept.shape[1] + n_vertex_predictors + 1, 2)
    )
    reference_counts = np.zeros(n_vertices, dtype=int)
    residual_standard_errors = np.full(n_vertices, np.nan, dtype=float)
    residual_dfs = np.full(n_vertices, np.nan, dtype=float)
    residual_ranks = np.full(n_vertices, np.nan, dtype=float)
    residual_matrix = np.full(ref_data_clean.shape, np.nan, dtype=float)
    all_zero_vertices = np.zeros(n_vertices, dtype=bool)
    
    # Fit normative model for each vertex
    for i in range(n_vertices):
        vertex_subject_mask = np.isfinite(ref_data_clean[:, i])
        if reference_vertex_covariates is not None:
            vertex_subject_mask &= np.isfinite(reference_vertex_covariates[:, i, :]).all(axis=1)
        reference_counts[i] = int(vertex_subject_mask.sum())
        if reference_counts[i] < min_reference_subjects:
            normative_data[i, :] = np.nan
            continue

        vertex_data = ref_data_clean[vertex_subject_mask, i]
        vertex_design = X_ref_with_intercept[vertex_subject_mask]
        if reference_vertex_covariates is not None:
            vertex_design = np.hstack([vertex_design, reference_vertex_covariates[vertex_subject_mask, i, :]])

        # Skip vertices with all zeros
        if np.all(vertex_data == 0):
            # Initialize with zeros and set std to 1 to avoid division by zero
            normative_data[i, :] = 0
            normative_data[i, -1] = 1
            residual_standard_errors[i] = 1.0
            residual_dfs[i] = max(vertex_design.shape[0] - np.linalg.matrix_rank(vertex_design), 0)
            residual_ranks[i] = np.linalg.matrix_rank(vertex_design)
            residual_matrix[vertex_subject_mask, i] = 0.0
            all_zero_vertices[i] = True
            continue

        # Fit linear regression: y = intercept + coef_1*x_1 + coef_2*x_2 + ... + coef_n*x_n
        coefficients = np.linalg.lstsq(vertex_design, vertex_data, rcond=None)[0]

        # Calculate residuals and their standard deviation
        predicted = vertex_design @ coefficients
        residuals = vertex_data - predicted
        residual_matrix[vertex_subject_mask, i] = residuals
        std_residuals = np.std(residuals)
        residual_standard_error, residual_df, residual_rank = _residual_standard_error(residuals, vertex_design)
        residual_standard_errors[i] = residual_standard_error
        residual_dfs[i] = residual_df
        residual_ranks[i] = residual_rank
        if not np.isfinite(std_residuals) or std_residuals <= 0:
            normative_data[i, :] = np.nan
            continue
        
        # Store all coefficients plus std_residuals
        normative_data[i, :-1] = coefficients
        normative_data[i, -1] = std_residuals

    shash_parameters = None
    shash_fit_success = np.zeros(n_vertices, dtype=bool)
    if wscore_distribution == "shash":
        cache_key = (
            "surface_shash",
            id(cache_reference_data),
            id(demographics_ref),
            id(cache_reference_covariates),
            tuple(normative_columns),
            wscore_preprocessing,
            wscore_covariate_model,
            blur_depth_model,
        )
        cached_fit = wscore_fit_cache.get(cache_key) if wscore_fit_cache is not None else None
        cache_matches = (
            cached_fit is not None
            and cached_fit["reference_data"] is cache_reference_data
            and cached_fit["demographics_ref"] is demographics_ref
            and cached_fit["reference_covariates"] is cache_reference_covariates
        )
        if cache_matches:
            shash_parameters = cached_fit["parameters"]
            shash_fit_success = cached_fit["fit_success"]
        else:
            shash_parameters, shash_fit_success = _fit_shash_quantile_parameters(
                residual_matrix
            )
            shash_parameters[all_zero_vertices] = np.array([0.0, 1.0, 0.0, 1.0])
            shash_fit_success[all_zero_vertices] = True
            if wscore_fit_cache is not None:
                wscore_fit_cache[cache_key] = {
                    "reference_data": cache_reference_data,
                    "demographics_ref": demographics_ref,
                    "reference_covariates": cache_reference_covariates,
                    "parameters": shash_parameters,
                    "fit_success": shash_fit_success,
                }

    residual_location = np.zeros(n_vertices, dtype=float)
    residual_scale = normative_data[:, -1].copy()
    if wscore_distribution in {
        "gaussian_mad",
        "gaussian_winsor",
        "student_t",
    }:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if wscore_distribution == "gaussian_mad":
                residual_location = np.nanmedian(residual_matrix, axis=0)
                residual_scale = 1.482602218505602 * np.nanmedian(
                    np.abs(residual_matrix - residual_location[None, :]), axis=0
                )
            elif wscore_distribution == "gaussian_winsor":
                bounds = np.nanpercentile(residual_matrix, [5.0, 95.0], axis=0)
                winsorized = np.clip(residual_matrix, bounds[0], bounds[1])
                residual_location = np.nanmean(winsorized, axis=0)
                residual_scale = np.nanstd(winsorized, axis=0)
            else:
                quartiles = np.nanpercentile(residual_matrix, [25.0, 50.0, 75.0], axis=0)
                residual_location = quartiles[1]
                residual_scale = (quartiles[2] - quartiles[0]) / (
                    2.0 * stats.t.ppf(0.75, 5)
                )
        invalid_scale = ~np.isfinite(residual_scale) | (
            residual_scale <= np.finfo(float).eps
        )
        residual_location[invalid_scale] = 0.0
        residual_scale[invalid_scale] = normative_data[invalid_scale, -1]

    
    intensity_dominant_null_threshold = None
    intensity_dominant_calibration_factor = 1.0
    if intensity_component_vertices is not None:
        # Selecting the maximum of four W-scores inflates its null magnitude.
        # Estimate one control-derived calibration factor so 1.96 retains its
        # familiar two-sided 5% interpretation without changing map ranking.
        with np.errstate(divide="ignore", invalid="ignore"):
            gaussian_control_scores = residual_matrix / normative_data[:, -1][None, :]
        gaussian_control_scores = gaussian_control_scores.reshape(
            gaussian_control_scores.shape[0],
            intensity_component_count,
            intensity_component_vertices,
        )
        control_dominant = np.nanmax(np.abs(gaussian_control_scores), axis=1)
        finite_null = control_dominant[np.isfinite(control_dominant)]
        if finite_null.size:
            intensity_dominant_null_threshold = float(
                np.percentile(finite_null, 95.0)
            )
            if intensity_dominant_null_threshold > np.finfo(float).eps:
                intensity_dominant_calibration_factor = (
                    intensity_dominant_null_threshold / 1.96
                )

    # Calculate W-scores for patient
    wscores = np.zeros(n_vertices)
    prediction_leverages = np.full(n_vertices, np.nan, dtype=float)
    prediction_denominators = np.full(n_vertices, np.nan, dtype=float)
    
    for i in range(n_vertices):
        # Get coefficients and standard deviation
        coefficients = normative_data[i, :-1]  # All coefficients (intercept + predictors)
        std_residuals = normative_data[i, -1]  # Standard deviation of residuals
        
        if not np.isfinite(std_residuals) or std_residuals <= 0:
            wscores[i] = np.nan
            continue

        # Calculate expected value: intercept + demographic and vertex covariate terms.
        patient_design = np.r_[1.0, X_pat[:n_predictors]]
        expected = float(coefficients[0] + np.dot(coefficients[1:1 + n_predictors], X_pat[:n_predictors]))
        if patient_vertex_covariates is not None:
            patient_covariates = patient_vertex_covariates[i]
            if not np.isfinite(patient_covariates).all():
                wscores[i] = np.nan
                continue
            covariate_start = 1 + n_predictors
            expected += float(np.dot(coefficients[covariate_start:covariate_start + n_vertex_predictors], patient_covariates))
            patient_design = np.r_[patient_design, patient_covariates]

        denominator = std_residuals
        leverage = np.nan
        if use_prediction_uncertainty:
            vertex_subject_mask = np.isfinite(ref_data_clean[:, i])
            if reference_vertex_covariates is not None:
                vertex_subject_mask &= np.isfinite(reference_vertex_covariates[:, i, :]).all(axis=1)
            vertex_design = X_ref_with_intercept[vertex_subject_mask]
            if reference_vertex_covariates is not None:
                vertex_design = np.hstack([vertex_design, reference_vertex_covariates[vertex_subject_mask, i, :]])
            gaussian_denominator, leverage = _prediction_uncertainty_denominator(
                residual_standard_errors[i],
                vertex_design,
                patient_design,
            )
            prediction_leverages[i] = leverage
            if not np.isfinite(gaussian_denominator) or gaussian_denominator <= 0:
                wscores[i] = np.nan
                continue
            if wscore_distribution == "gaussian":
                denominator = gaussian_denominator
            elif wscore_distribution == "shash":
                denominator = shash_parameters[i, 1] * np.sqrt(1.0 + leverage)
            elif wscore_distribution != "empirical":
                denominator = residual_scale[i] * np.sqrt(1.0 + leverage)
            prediction_denominators[i] = denominator
        
        # Calculate W-score
        patient_residual = patient_data[i] - expected
        if wscore_distribution == "shash":
            location, scale, skew, tail_weight = shash_parameters[i]
            if use_prediction_uncertainty:
                scale = denominator
            wscores[i] = _shash_to_normal_score(
                patient_residual,
                location,
                scale,
                skew,
                tail_weight,
            )
        elif wscore_distribution in {"gaussian_mad", "gaussian_winsor"}:
            scale = residual_scale[i]
            if use_prediction_uncertainty:
                scale *= np.sqrt(1.0 + leverage)
            wscores[i] = (patient_residual - residual_location[i]) / scale
        elif wscore_distribution == "student_t":
            scale = residual_scale[i]
            if use_prediction_uncertainty:
                scale *= np.sqrt(1.0 + leverage)
            probability = stats.t.cdf(
                (patient_residual - residual_location[i]) / scale,
                df=5,
            )
            wscores[i] = stats.norm.ppf(np.clip(probability, 1e-12, 1.0 - 1e-12))
        elif wscore_distribution == "empirical":
            empirical_residual = patient_residual
            if use_prediction_uncertainty:
                empirical_residual /= np.sqrt(1.0 + leverage)
            control_residuals = residual_matrix[:, i]
            control_residuals = control_residuals[np.isfinite(control_residuals)]
            rank = np.sum(control_residuals <= empirical_residual)
            probability = (rank + 0.5) / (control_residuals.size + 1.0)
            wscores[i] = stats.norm.ppf(np.clip(probability, 1e-12, 1.0 - 1e-12))
        else:
            wscores[i] = patient_residual / denominator

    dominant_depth_indices = None
    if intensity_component_vertices is not None:
        component_wscores = wscores.reshape(
            intensity_component_count, intensity_component_vertices
        )
        dominant_depth_indices = np.argmax(np.abs(component_wscores), axis=0)
        wscores = np.take_along_axis(
            component_wscores, dominant_depth_indices[None, :], axis=0
        )[0]
        wscores = wscores / intensity_dominant_calibration_factor
    elif blur_component_vertices is not None:
        component_wscores = wscores.reshape(2, blur_component_vertices)
        if blur_depth_model == "gradient_flattening":
            # Positive values specifically encode lower-than-normal gradients.
            wscores = -np.mean(component_wscores, axis=0)
        else:
            magnitude = np.sqrt(np.mean(component_wscores**2, axis=0))
            # Preserve a useful direction for signed clinical maps while retaining
            # the exact RMS magnitude evaluated by the lesion benchmark.
            wscores = np.copysign(magnitude, component_wscores[0])

    output_vertices = (
        intensity_component_vertices
        if intensity_component_vertices is not None
        else blur_component_vertices
        if blur_component_vertices is not None
        else n_vertices
    )
    surface_smoothing_applied = bool(
        wscore_surface_smoothing_iterations > 0 and output_vertices == 32492
    )
    if surface_smoothing_applied:
        wscores = _smooth_fslr32k_scores(
            wscores, wscore_surface_smoothing_iterations
        )

    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save W-score map
    wscore_data = nib.gifti.gifti.GiftiDataArray(
        data=wscores.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    wscore_gii = nib.gifti.GiftiImage(darrays=[wscore_data])
    nib.save(wscore_gii, output_file)

    dominant_depth_file = None
    if dominant_depth_indices is not None:
        dominant_depth_file = (
            output_file.replace(".func.gii", "_dominant-depth.shape.gii")
            if output_file.endswith(".func.gii")
            else f"{output_file}_dominant-depth.shape.gii"
        )
        dominant_depth_data = nib.gifti.gifti.GiftiDataArray(
            data=dominant_depth_indices.astype(np.float32),
            intent="NIFTI_INTENT_NORMAL",
        )
        nib.save(
            nib.gifti.GiftiImage(darrays=[dominant_depth_data]),
            dominant_depth_file,
        )

    if reference_vertex_covariates is not None:
        reference_count_min = int(np.nanmin(reference_counts))
        reference_count_max = int(np.nanmax(reference_counts))
    else:
        reference_count_min = len(X_ref)
        reference_count_max = len(X_ref)

    model_info = {
        "method": "wscore",
        "wscore_distribution": wscore_distribution,
        "wscore_preprocessing": wscore_preprocessing,
        "wscore_covariate_model": wscore_covariate_model,
        "wscore_surface_smoothing_iterations": wscore_surface_smoothing_iterations,
        "wscore_surface_smoothing_applied": surface_smoothing_applied,
        "blur_depth_model": (
            blur_depth_model
            if input_was_multidepth and intensity_component_vertices is None
            else None
        ),
        "blur_depth_components": (
            blur_component_names
            if blur_component_vertices is not None
            else ["mean"]
            if input_was_multidepth and intensity_component_vertices is None
            else []
        ),
        "blur_depth_combination": (
            blur_component_combination
            if blur_component_vertices is not None
            else "mean_across_depths"
            if input_was_multidepth and intensity_component_vertices is None
            else None
        ),
        "blur_depth_self_normalization": (
            "bilateral_median_absolute_intensity_across_all_four_depths"
            if blur_depth_model == "gradient_flattening"
            and blur_component_vertices is not None
            else None
        ),
        "intensity_depth_model": intensity_depth_model,
        "intensity_depth_components": (
            ["midthickness", "white", "SWM1", "SWM2"]
            if intensity_component_vertices is not None
            else []
        ),
        "intensity_depth_combination": (
            "signed_depth_wscore_with_maximum_absolute_value"
            if intensity_component_vertices is not None
            else None
        ),
        "intensity_dominant_calibration_basis": (
            "global_control_gaussian_residual_maximum_p95_to_1p96"
            if intensity_component_vertices is not None
            else None
        ),
        "intensity_dominant_uncalibrated_control_p95": intensity_dominant_null_threshold,
        "intensity_dominant_calibration_factor": (
            intensity_dominant_calibration_factor
            if intensity_component_vertices is not None
            else None
        ),
        "dominant_depth_file": dominant_depth_file,
        "normative_columns": list(normative_columns),
        "effective_normative_columns": effective_normative_columns,
        "uses_prediction_uncertainty": bool(use_prediction_uncertainty),
        "prediction_denominator": (
            "shash_scale_times_sqrt_1_plus_patient_leverage"
            if use_prediction_uncertainty and wscore_distribution == "shash"
            else "residual_standard_error_times_sqrt_1_plus_patient_leverage"
            if use_prediction_uncertainty
            else "robust_scale_times_sqrt_1_plus_patient_leverage"
            if use_prediction_uncertainty and wscore_distribution in {
                "gaussian_mad", "gaussian_winsor", "student_t"
            }
            else "empirical_control_cdf_after_leverage_adjustment"
            if use_prediction_uncertainty and wscore_distribution == "empirical"
            else "shash_scale"
            if wscore_distribution == "shash"
            else "median_absolute_deviation"
            if wscore_distribution == "gaussian_mad"
            else "winsorized_residual_sd"
            if wscore_distribution == "gaussian_winsor"
            else "student_t_iqr_scale"
            if wscore_distribution == "student_t"
            else "empirical_control_cdf"
            if wscore_distribution == "empirical"
            else "residual_sd"
        ),
        "uses_vertex_covariates": bool(n_vertex_predictors > 0),
        "vertex_covariate_count": int(n_vertex_predictors),
        "vertex_covariate_names": [
            f"gyrification_curvature_{idx + 1}" for idx in range(n_vertex_predictors)
        ],
        "uses_reference_vertex_mask": False,
        "reference_count": int(len(X_ref)),
        "reference_count_min": int(reference_count_min),
        "reference_count_max": int(reference_count_max),
        "min_reference_subjects": int(min_reference_subjects),
    }
    if wscore_distribution == "shash":
        model_info.update(
            {
                "shash_parameterization": (
                    "z=sinh(tail_weight*asinh((residual-location)/scale)-skew)"
                ),
                "shash_fit_method": "quantile_matching_10_25_50_75_90",
                "shash_fit_vertices": int(shash_fit_success.sum()),
                "shash_fallback_vertices": int((~shash_fit_success).sum()),
            }
        )
    if use_prediction_uncertainty and np.isfinite(prediction_leverages).any():
        model_info.update(
            {
                "prediction_leverage_min": float(np.nanmin(prediction_leverages)),
                "prediction_leverage_max": float(np.nanmax(prediction_leverages)),
                "prediction_denominator_min": float(np.nanmin(prediction_denominators)),
                "prediction_denominator_max": float(np.nanmax(prediction_denominators)),
                "residual_standard_error_min": float(np.nanmin(residual_standard_errors)),
                "residual_standard_error_max": float(np.nanmax(residual_standard_errors)),
                "residual_df_min": int(np.nanmin(residual_dfs)),
                "residual_df_max": int(np.nanmax(residual_dfs)),
                "design_rank_min": int(np.nanmin(residual_ranks)),
                "design_rank_max": int(np.nanmax(residual_ranks)),
            }
        )
    if output_file.endswith(".func.gii"):
        model_file = output_file.replace(".func.gii", "_model.json")
    else:
        model_file = f"{output_file}_model.json"
    with open(model_file, "w") as f:
        json.dump(model_info, f, indent=2)
    
    if verbose:
        print(f"    Saved W-score map: {output_file}")
        if n_vertex_predictors > 0:
            print(
                "      Model includes "
                f"{n_vertex_predictors} vertex-wise gyrification covariate(s); "
                f"reference subjects per vertex: {reference_count_min}-{reference_count_max}"
            )
        if use_prediction_uncertainty:
            print("      Predictive W-score distribution includes patient leverage")
        if wscore_distribution == "shash":
            print(
                "      SHASH residual model: "
                f"{int(shash_fit_success.sum())}/{n_vertices} vertices fitted"
            )
        print(f"      Saved W-score model metadata: {model_file}")
    
    return {
        'wscore_file': output_file,
        'mean_wscore': np.nanmean(wscores),
        'std_wscore': np.nanstd(wscores),
        'normative_data': normative_data,
        'reference_count': len(X_ref),
        'reference_count_min': reference_count_min,
        'reference_count_max': reference_count_max,
        'vertex_covariate_count': n_vertex_predictors,
        'uses_prediction_uncertainty': bool(use_prediction_uncertainty),
        'wscore_distribution': wscore_distribution,
        'wscore_preprocessing': wscore_preprocessing,
        'wscore_covariate_model': wscore_covariate_model,
        'wscore_surface_smoothing_iterations': wscore_surface_smoothing_iterations,
        'wscore_surface_smoothing_applied': surface_smoothing_applied,
        'blur_depth_model': (
            blur_depth_model
            if input_was_multidepth and intensity_component_vertices is None
            else None
        ),
        'intensity_depth_model': intensity_depth_model,
        'dominant_depth_file': dominant_depth_file,
        'dominant_depth_indices': dominant_depth_indices,
        'shash_parameters': shash_parameters,
        'shash_fit_success': shash_fit_success,
        'residual_standard_errors': residual_standard_errors,
        'residual_dfs': residual_dfs,
        'residual_ranks': residual_ranks,
        'model_file': model_file,
    }

def calculate_wscore_csv(
    reference_data,
    patient_data,
    demographics_ref,
    demographics_pat,
    output_file,
    normative_columns=['age', 'sex'],
    verbose=True,
    use_prediction_uncertainty=False,
    wscore_distribution="gaussian",
    wscore_preprocessing="none",
    wscore_covariate_model="linear",
    wscore_fit_cache=None,
):
    """
    Calculate W-scores for subcortical data using normative modeling and save as CSV file.
    
    Parameters:
    -----------
    reference_data : pd.DataFrame
        Combined reference data from all subjects
    patient_data : pd.DataFrame
        Patient data (single subject)
    demographics_ref : pd.DataFrame
        Demographics data for reference subjects
    demographics_pat : pd.Series or pd.DataFrame
        Demographics data for patient
    output_file : str
        Path to save the W-score CSV
    normative_columns : list, default=['age', 'sex']
        List of demographic columns to use for normative modeling
    verbose : bool, default=True
        If True, prints processing information
    use_prediction_uncertainty : bool, default=False
        If True, divide by the full OLS predictive standard deviation for a
        new subject, using residual standard error and patient leverage.
    wscore_distribution : str, default="gaussian"
        Distribution used for control residuals.
    wscore_preprocessing : {"none", "spatial_zscore", "spatial_robust_z"}, default="none"
        Optional robust normalization across structures within each subject.
    wscore_covariate_model : str, default="linear"
        Optional quadratic-age and/or age-by-sex demographic design.
    
    Returns:
    --------
    dict
        Dictionary containing W-score statistics and file path
    """
    wscore_distribution = _normalize_wscore_distribution(wscore_distribution)
    wscore_preprocessing = _normalize_wscore_preprocessing(wscore_preprocessing)
    wscore_covariate_model = _normalize_wscore_covariate_model(wscore_covariate_model)
    gaussian_process_kernel = GAUSSIAN_PROCESS_KERNELS.get(
        wscore_covariate_model
    )
    if gaussian_process_kernel is not None:
        if wscore_distribution != "gaussian":
            raise ValueError(
                "Gaussian-process models use their predictive Gaussian "
                "distribution; set wscore_distribution='gaussian'"
            )
        return calculate_gaussian_process_csv(
            reference_data=reference_data,
            patient_data=patient_data,
            demographics_ref=demographics_ref,
            demographics_pat=demographics_pat,
            output_file=output_file,
            normative_columns=normative_columns,
            verbose=verbose,
            gaussian_process_kernel=gaussian_process_kernel,
            wscore_covariate_model=wscore_covariate_model,
            wscore_preprocessing=wscore_preprocessing,
            wscore_fit_cache=wscore_fit_cache,
        )
    cache_reference_data = reference_data

    # Get structure columns (exclude metadata columns)
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)

    if wscore_preprocessing in {"spatial_zscore", "spatial_robust_z"}:
        reference_data = reference_data.copy()
        patient_data = patient_data.copy()
        transform = (
            _spatial_zscore
            if wscore_preprocessing == "spatial_zscore"
            else _spatial_robust_z
        )
        reference_data.loc[:, structures] = transform(
            reference_data[structures].to_numpy(dtype=float)
        )
        patient_data.loc[:, structures] = transform(
            patient_data[structures].to_numpy(dtype=float)
        )
    
    # Prepare demographics data
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    
    # Convert to numpy arrays
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    
    # Remove subjects with missing demographic data
    mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[mask]
    ref_data_clean = reference_data.loc[mask]
    
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")

    X_ref, X_pat, effective_normative_columns = _augment_wscore_covariates(
        X_ref,
        X_pat,
        normative_columns,
        wscore_covariate_model,
    )
    
    # Add intercept term
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    
    # Calculate W-scores for each structure
    patient_wscores = {}
    structure_stats = {}
    
    for structure in structures:
        if structure in reference_data.columns and structure in patient_data.columns:
            ref_values = ref_data_clean[structure].values
            pat_value = patient_data[structure].iloc[0]
            structure_mask = np.isfinite(ref_values)
            ref_values = ref_values[structure_mask]
            structure_design = X_ref_with_intercept[structure_mask]
            if ref_values.size == 0:
                patient_wscores[structure] = np.nan
                continue
            
            # Skip if all reference values are zero
            if np.all(ref_values == 0):
                patient_wscores[structure] = 0
                structure_stats[structure] = {
                    'intercept': 0, 'coefficients': [0] * X_ref.shape[1],
                    'std_residuals': 1,
                    'residual_standard_error': 1,
                    'residual_df': max(X_ref_with_intercept.shape[0] - np.linalg.matrix_rank(X_ref_with_intercept), 0),
                    'design_rank': np.linalg.matrix_rank(X_ref_with_intercept),
                    'prediction_denominator': 1,
                    'prediction_leverage': 0,
                    'count': len(ref_values),
                    'wscore_distribution': wscore_distribution,
                }
                continue
            
            # Fit linear regression
            coefficients = np.linalg.lstsq(structure_design, ref_values, rcond=None)[0]
            
            # Calculate residuals and their standard deviation
            predicted = structure_design @ coefficients
            residuals = ref_values - predicted
            std_residuals = np.std(residuals)
            residual_standard_error, residual_df, residual_rank = _residual_standard_error(
                residuals,
                structure_design,
            )
            shash_parameters = None
            shash_fit_success = False
            if wscore_distribution == "shash":
                cache_key = (
                    "csv_shash",
                    id(cache_reference_data),
                    id(demographics_ref),
                    tuple(normative_columns),
                    wscore_preprocessing,
                    wscore_covariate_model,
                    structure,
                )
                cached_fit = wscore_fit_cache.get(cache_key) if wscore_fit_cache is not None else None
                cache_matches = (
                    cached_fit is not None
                    and cached_fit["reference_data"] is cache_reference_data
                    and cached_fit["demographics_ref"] is demographics_ref
                )
                if cache_matches:
                    shash_parameters = cached_fit["parameters"]
                    shash_fit_success = cached_fit["fit_success"]
                else:
                    fitted_parameters, fitted_success = _fit_shash_quantile_parameters(
                        residuals
                    )
                    shash_parameters = fitted_parameters[0]
                    shash_fit_success = bool(fitted_success[0])
                    if wscore_fit_cache is not None:
                        wscore_fit_cache[cache_key] = {
                            "reference_data": cache_reference_data,
                            "demographics_ref": demographics_ref,
                            "parameters": shash_parameters,
                            "fit_success": shash_fit_success,
                        }

            residual_location = 0.0
            residual_scale = std_residuals
            if wscore_distribution == "gaussian_mad":
                residual_location = float(np.nanmedian(residuals))
                residual_scale = float(
                    1.482602218505602
                    * np.nanmedian(np.abs(residuals - residual_location))
                )
            elif wscore_distribution == "gaussian_winsor":
                lower, upper = np.nanpercentile(residuals, [5.0, 95.0])
                winsorized = np.clip(residuals, lower, upper)
                residual_location = float(np.nanmean(winsorized))
                residual_scale = float(np.nanstd(winsorized))
            elif wscore_distribution == "student_t":
                lower, residual_location, upper = np.nanpercentile(
                    residuals, [25.0, 50.0, 75.0]
                )
                residual_scale = float(
                    (upper - lower) / (2.0 * stats.t.ppf(0.75, 5))
                )
            if not np.isfinite(residual_scale) or residual_scale <= np.finfo(float).eps:
                residual_location = 0.0
                residual_scale = std_residuals
            
            if std_residuals == 0:
                std_residuals = 1  # Avoid division by zero
                if not np.isfinite(residual_standard_error) or residual_standard_error <= 0:
                    residual_standard_error = 1
            
            # Predict expected value for patient
            patient_design = np.r_[1.0, X_pat]
            expected = coefficients[0] + np.dot(coefficients[1:], X_pat)
            denominator = (
                residual_scale
                if wscore_distribution in {
                    "gaussian_mad", "gaussian_winsor", "student_t"
                }
                else std_residuals
            )
            leverage = np.nan
            if use_prediction_uncertainty:
                denominator, leverage = _prediction_uncertainty_denominator(
                    residual_standard_error,
                    structure_design,
                    patient_design,
                )
                if not np.isfinite(denominator) or denominator <= 0:
                    patient_wscores[structure] = np.nan
                    structure_stats[structure] = {
                        'intercept': coefficients[0],
                        'coefficients': coefficients[1:].tolist(),
                        'std_residuals': std_residuals,
                        'residual_standard_error': residual_standard_error,
                        'residual_df': residual_df,
                        'design_rank': residual_rank,
                        'prediction_denominator': denominator,
                        'prediction_leverage': leverage,
                        'count': len(ref_values)
                    }
                    continue
                if wscore_distribution in {
                    "gaussian_mad", "gaussian_winsor", "student_t"
                }:
                    denominator = residual_scale * np.sqrt(1.0 + leverage)
            
            # Calculate W-score
            patient_residual = pat_value - expected
            if wscore_distribution == "shash":
                location, scale, skew, tail_weight = shash_parameters
                if use_prediction_uncertainty and np.isfinite(leverage):
                    scale *= np.sqrt(1.0 + leverage)
                    denominator = scale
                wscore = _shash_to_normal_score(
                    patient_residual,
                    location,
                    scale,
                    skew,
                    tail_weight,
                )
            elif wscore_distribution in {"gaussian_mad", "gaussian_winsor"}:
                wscore = (patient_residual - residual_location) / denominator
            elif wscore_distribution == "student_t":
                standardized = (patient_residual - residual_location) / denominator
                probability = stats.t.cdf(standardized, df=5)
                wscore = stats.norm.ppf(np.clip(probability, 1e-12, 1.0 - 1e-12))
            elif wscore_distribution == "empirical":
                empirical_residual = patient_residual
                if use_prediction_uncertainty:
                    empirical_residual /= np.sqrt(1.0 + leverage)
                rank = np.sum(residuals <= empirical_residual)
                probability = (rank + 0.5) / (residuals.size + 1.0)
                wscore = stats.norm.ppf(np.clip(probability, 1e-12, 1.0 - 1e-12))
            else:
                wscore = patient_residual / denominator
            patient_wscores[structure] = wscore
            
            structure_stats[structure] = {
                'intercept': coefficients[0],
                'coefficients': coefficients[1:].tolist(),
                'std_residuals': std_residuals,
                'residual_standard_error': residual_standard_error,
                'residual_df': residual_df,
                'design_rank': residual_rank,
                'prediction_denominator': denominator,
                'prediction_leverage': leverage,
                'count': len(ref_values),
                'wscore_distribution': wscore_distribution,
                'shash_parameters': (
                    {
                        'location': float(shash_parameters[0]),
                        'scale': float(shash_parameters[1]),
                        'skew': float(shash_parameters[2]),
                        'tail_weight': float(shash_parameters[3]),
                    }
                    if shash_parameters is not None
                    else None
                ),
                'shash_fit_success': shash_fit_success,
            }
                
    # Create output directory and save CSV
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    wscore_df = pd.DataFrame([patient_wscores])
    wscore_df.to_csv(output_file, index=False)
    
    if verbose:
        print(f"    Saved subcortical W-score CSV: {output_file}")
    
    return {
        'wscore_file': output_file,
        'wscores': patient_wscores,
        'structure_stats': structure_stats,
        'uses_prediction_uncertainty': bool(use_prediction_uncertainty),
        'wscore_distribution': wscore_distribution,
        'wscore_preprocessing': wscore_preprocessing,
        'wscore_covariate_model': wscore_covariate_model,
        'effective_normative_columns': effective_normative_columns,
    }

def calculate_zscore_maps(
    reference_data,
    patient_data,
    output_file,
    analysis='regional',
    verbose=True,
    min_reference_subjects=2,
    reference_vertex_covariates=None,
    patient_vertex_covariates=None,
    wscore_preprocessing="none",
):
    """
    Calculate z-scores for patient data against reference data and save as GIFTI file.
    
    Parameters:
    -----------
    reference_data : np.ndarray
        Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
    patient_data : np.ndarray
        Patient data with shape (n_vertices,) or (n_vertices, n_depths)
    output_file : str
        Path to save the z-score map
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing z-score statistics and file path
    """
    if len(reference_data) == 0:
        raise ValueError("No reference data provided")

    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate z-scores")

    # Per-subject spatial self-normalization (the stage-1 "self_normalization"
    # axis). This must run in the z-score path too -- it is a within-subject map
    # transform, not a w-score-only step -- otherwise sweeping wscore_preprocessing
    # under the simplest (z-score) baseline would be a silent no-op.
    if wscore_preprocessing in ("spatial_zscore", "spatial_robust_z"):
        _preprocess = (
            _spatial_zscore
            if wscore_preprocessing == "spatial_zscore"
            else _spatial_robust_z
        )

        def _self_normalize(array, is_reference):
            array = np.asarray(array, dtype=float)
            if is_reference:
                # (subjects, vertices[, depths]) -> normalize each subject over all locations
                if array.ndim <= 2:
                    return _preprocess(array)
                shape = array.shape
                return _preprocess(array.reshape(shape[0], -1)).reshape(shape)
            # patient: (vertices[, depths]) -> one subject over all locations
            if array.ndim == 1:
                return _preprocess(array)
            shape = array.shape
            return _preprocess(array.reshape(-1)).reshape(shape)

        reference_data = _self_normalize(reference_data, True)
        patient_data = _self_normalize(patient_data, False)

    if reference_vertex_covariates is not None:
        if patient_vertex_covariates is None:
            raise ValueError("patient_vertex_covariates is required when reference_vertex_covariates is provided")
        if len(reference_data.shape) > 2:
            reference_data = np.mean(reference_data, axis=2)
        if len(patient_data.shape) > 1:
            patient_data = np.mean(patient_data, axis=1)

        reference_vertex_covariates = np.asarray(reference_vertex_covariates, dtype=float)
        if reference_vertex_covariates.ndim == 2:
            reference_vertex_covariates = reference_vertex_covariates[:, :, None]
        if reference_vertex_covariates.shape[:2] != reference_data.shape[:2]:
            raise ValueError(
                "reference_vertex_covariates must have shape "
                f"{reference_data.shape[:2]} or {reference_data.shape[:2]} + (n_covariates,), "
                f"got {reference_vertex_covariates.shape}"
            )
        patient_vertex_covariates = np.asarray(patient_vertex_covariates, dtype=float)
        if patient_vertex_covariates.ndim == 1:
            patient_vertex_covariates = patient_vertex_covariates[:, None]
        if patient_vertex_covariates.shape != reference_vertex_covariates.shape[1:]:
            raise ValueError(
                "patient_vertex_covariates must have shape "
                f"{reference_vertex_covariates.shape[1:]}, got {patient_vertex_covariates.shape}"
            )

        n_vertices = reference_data.shape[1]
        zscore_to_save = np.full(n_vertices, np.nan, dtype=np.float32)
        ref_mean = np.full(n_vertices, np.nan, dtype=np.float32)
        ref_std = np.full(n_vertices, np.nan, dtype=np.float32)
        n_vertex_predictors = reference_vertex_covariates.shape[2]
        min_reference_subjects = max(int(min_reference_subjects), n_vertex_predictors + 2)
        reference_counts = np.zeros(n_vertices, dtype=int)

        for vertex in range(n_vertices):
            vertex_subject_mask = (
                np.isfinite(reference_data[:, vertex])
                & np.isfinite(reference_vertex_covariates[:, vertex, :]).all(axis=1)
            )
            reference_counts[vertex] = int(vertex_subject_mask.sum())
            if reference_counts[vertex] < int(min_reference_subjects):
                continue

            vertex_reference = reference_data[vertex_subject_mask, vertex].astype(float)
            vertex_design = np.hstack([
                np.ones((reference_counts[vertex], 1), dtype=float),
                reference_vertex_covariates[vertex_subject_mask, vertex, :],
            ])
            coefficients = np.linalg.lstsq(vertex_design, vertex_reference, rcond=None)[0]
            residuals = vertex_reference - vertex_design @ coefficients
            residual_std = float(np.std(residuals))
            if not np.isfinite(residual_std) or residual_std <= 0:
                continue

            patient_covariates = patient_vertex_covariates[vertex]
            if not np.isfinite(patient_covariates).all():
                continue
            patient_design = np.r_[1.0, patient_covariates]
            expected = float(patient_design @ coefficients)
            zscore_to_save[vertex] = (patient_data[vertex] - expected) / residual_std
            ref_mean[vertex] = float(np.mean(vertex_reference))
            ref_std[vertex] = residual_std
    else:
        # Calculate reference statistics
        ref_mean = np.mean(reference_data, axis=0)
        ref_std = np.std(reference_data, axis=0)

        # Avoid division by zero
        ref_std[ref_std == 0] = 1.0

        # Calculate z-score
        zscore = (patient_data - ref_mean) / ref_std

        # For multi-depth data (blur features), average across depths
        if len(zscore.shape) > 1:
            mean_zscore = np.mean(zscore, axis=1)  # Average across depths
            zscore_to_save = mean_zscore
        else:
            zscore_to_save = zscore
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save z-score map
    zscore_data = nib.gifti.gifti.GiftiDataArray(
        data=zscore_to_save.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    zscore_gii = nib.gifti.GiftiImage(darrays=[zscore_data])
    nib.save(zscore_gii, output_file)
    
    if verbose:
        print(f"    Saved z-score map: {output_file}")
    
    return {
        'zscore_file': output_file,
        'mean_zscore': np.nanmean(zscore_to_save),
        'std_zscore': np.nanstd(zscore_to_save),
        'reference_mean': ref_mean,
        'reference_std': ref_std,
        'reference_count': len(reference_data),
        'reference_count_min': int(np.nanmin(reference_counts)) if reference_vertex_covariates is not None else len(reference_data),
        'reference_count_max': int(np.nanmax(reference_counts)) if reference_vertex_covariates is not None else len(reference_data),
        'vertex_covariate_count': reference_vertex_covariates.shape[2] if reference_vertex_covariates is not None else 0,
    }

def calculate_zscore_csv(reference_data, patient_data, output_file, verbose=True):
    """
    Calculate z-scores for subcortical data and save as CSV file.
    
    Parameters:
    -----------
    reference_data : pd.DataFrame
        Combined reference data from all subjects
    patient_data : pd.DataFrame
        Patient data (single subject)
    output_file : str
        Path to save the z-score CSV
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing z-score statistics and file path
    """
    # Get structure columns (exclude metadata columns)
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)
    
    # Calculate statistics for each structure
    structure_stats = {}
    for structure in structures:
        if structure in reference_data.columns:
            ref_values = reference_data[structure].values
            structure_stats[structure] = {
                'mean': np.mean(ref_values),
                'std': np.std(ref_values) if np.std(ref_values) > 0 else 1.0,
                'count': len(ref_values)
            }
    
    # Calculate z-scores for patient
    patient_zscores = {}
    for structure in structures:
        if structure in patient_data.columns and structure in structure_stats:
            pat_value = patient_data[structure].iloc[0]
            ref_mean = structure_stats[structure]['mean']
            ref_std = structure_stats[structure]['std']
            
            zscore = (pat_value - ref_mean) / ref_std
            patient_zscores[structure] = zscore
    
    # Create output directory and save CSV
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    zscore_df = pd.DataFrame([patient_zscores])
    zscore_df.to_csv(output_file, index=False)
    
    if verbose:
        print(f"    Saved subcortical z-score CSV: {output_file}")
    
    return {
        'zscore_file': output_file,
        'zscores': patient_zscores,
        'structure_stats': structure_stats
    }


def _subject_average_correlations(reference_data):
    """Return each subject's average correlation to the other controls."""
    data = np.asarray(reference_data, dtype=np.float64)
    if data.ndim == 3:
        data = np.nanmean(data, axis=2)
    if data.ndim != 2:
        raise ValueError(f"Expected reference data with 2 or 3 dimensions, got {data.shape}")

    n_subjects = data.shape[0]
    mean_correlations = np.full(n_subjects, np.nan, dtype=float)
    if n_subjects < 2:
        return mean_correlations

    finite = np.isfinite(data)
    with np.errstate(invalid="ignore"):
        row_means = np.divide(
            np.where(finite, data, 0.0).sum(axis=1),
            finite.sum(axis=1),
            out=np.full(n_subjects, np.nan, dtype=float),
            where=finite.sum(axis=1) > 0,
        )
    centered = np.where(finite, data - row_means[:, None], 0.0)
    norms = np.sqrt(np.sum(centered * centered, axis=1))
    denominator = np.outer(norms, norms)
    correlations = np.full((n_subjects, n_subjects), np.nan, dtype=float)
    np.divide(
        centered @ centered.T,
        denominator,
        out=correlations,
        where=denominator > 0,
    )
    np.fill_diagonal(correlations, np.nan)
    finite_correlations = np.isfinite(correlations)
    mean_correlations = np.divide(
        np.where(finite_correlations, correlations, 0.0).sum(axis=1),
        finite_correlations.sum(axis=1),
        out=np.full(n_subjects, np.nan, dtype=float),
        where=finite_correlations.sum(axis=1) > 0,
    )
    return mean_correlations


def _resolve_control_correlation_quantile(quantile, feature):
    """Per-feature drop-fraction for the RANK-based control correlation filter.

    ``quantile`` is either ``None`` (filter off for this feature -> keep every
    control), a scalar fraction in [0, 1) applied to every feature, or a mapping
    ``feature -> fraction`` (per-feature; absent/``None`` -> off for that feature).
    Returns ``None`` when off, else the float fraction of the lowest-correlation
    controls to drop within this feature. Lookup is case-insensitive to match the
    ``dataset.features`` strings (``"FLAIR"``, ``"T1w*blur"``, ...).
    """
    if quantile is None:
        return None
    if isinstance(quantile, Mapping):
        key = str(feature).lower()
        for name, value in quantile.items():
            if str(name).lower() == key:
                return None if value is None else float(value)
        return None
    return float(quantile)


def _filter_reference_controls_by_correlation(
    reference_data,
    reference_subjects,
    reference_demographics_df=None,
    reference_vertex_covariates=None,
    quantile=None,
    min_controls=10,
    verbose=True,
    context="reference map",
):
    """RANK-based control QC: drop the bottom ``quantile`` fraction of controls by
    their mean map correlation to the other controls, *within this feature*.

    The cut adapts to each feature's own correlation distribution, so a fixed
    fraction -- never all -- is removed. This prevents an inherently low-correlation
    feature (e.g. FA) from losing its entire normative reference (which would drop
    the feature out of the objective). A floor of ``min_controls`` (capped at the
    number of finite-correlation controls) guarantees a viable reference. NaN
    correlations rank as the worst and are dropped first.

    ``quantile is None`` means "off for this feature": keep EVERY control (including
    controls with an undefined mean correlation), while still populating the
    asymmetry QC set -- byte-for-byte equivalent to disabling the filter.
    """
    reference_subjects = list(reference_subjects)
    if len(reference_subjects) < 2:
        return reference_data, reference_subjects, reference_demographics_df, reference_vertex_covariates

    mean_correlations = _subject_average_correlations(reference_data)
    if quantile is None:
        # "off" (per-feature or global): keep every control.
        keep_mask = np.ones(len(reference_subjects), dtype=bool)
        filter_desc = "keep-all"
    else:
        finite = np.isfinite(mean_correlations)
        n_finite = int(finite.sum())
        if n_finite < 2:
            keep_mask = finite.copy()   # nothing rankable -> keep the finite ones
        else:
            cut = float(np.quantile(mean_correlations[finite], float(quantile)))
            keep_mask = finite & (mean_correlations >= cut)
            floor = min(int(min_controls), n_finite)
            if int(keep_mask.sum()) < floor:
                # Quantile cut over-pruned (ties/small N): keep the top-``floor``
                # controls by correlation instead, so the reference stays viable.
                order = np.argsort(np.where(finite, mean_correlations, -np.inf))[::-1]
                keep_mask = np.zeros(len(reference_subjects), dtype=bool)
                keep_mask[order[:floor]] = True
        filter_desc = f"drop-bottom {float(quantile):.0%}"
    removed = [
        (subject, mean_correlations[i])
        for i, subject in enumerate(reference_subjects)
        if not keep_mask[i]
    ]

    if verbose:
        kept_count = int(keep_mask.sum())
        removed_count = len(reference_subjects) - kept_count
        print(
            f"    Control correlation filter for {context}: "
            f"kept {kept_count}/{len(reference_subjects)} controls "
            f"({filter_desc})"
        )
        if removed_count:
            preview = ", ".join(
                f"{pid}/{sid}={corr:.3f}" if np.isfinite(corr) else f"{pid}/{sid}=nan"
                for (pid, sid), corr in removed[:5]
            )
            suffix = f", ... and {removed_count - 5} more" if removed_count > 5 else ""
            print(f"      Removed low-correlation controls: {preview}{suffix}")

    filtered_data = reference_data[keep_mask]
    filtered_subjects = [
        subject for subject, keep in zip(reference_subjects, keep_mask) if keep
    ]

    filtered_demographics = reference_demographics_df
    if reference_demographics_df is not None:
        filtered_demographics = reference_demographics_df.iloc[keep_mask].reset_index(drop=True)

    filtered_covariates = reference_vertex_covariates
    if reference_vertex_covariates is not None:
        filtered_covariates = reference_vertex_covariates[keep_mask]

    return filtered_data, filtered_subjects, filtered_demographics, filtered_covariates


def _filter_reference_controls_to_subjects(
    reference_data,
    reference_subjects,
    allowed_subjects,
    reference_demographics_df=None,
    reference_vertex_covariates=None,
    verbose=True,
    context="asymmetry map",
):
    """Apply a control QC set learned from corresponding regional maps."""
    allowed_subjects = set(allowed_subjects)
    keep_mask = np.asarray(
        [subject in allowed_subjects for subject in reference_subjects],
        dtype=bool,
    )
    filtered_subjects = [
        subject for subject, keep in zip(reference_subjects, keep_mask) if keep
    ]
    if verbose:
        print(
            f"    Regional control correlation mask for {context}: "
            f"kept {int(keep_mask.sum())}/{len(reference_subjects)} controls"
        )

    filtered_demographics = reference_demographics_df
    if reference_demographics_df is not None:
        filtered_demographics = reference_demographics_df.iloc[keep_mask].reset_index(drop=True)

    filtered_covariates = reference_vertex_covariates
    if reference_vertex_covariates is not None:
        filtered_covariates = reference_vertex_covariates[keep_mask]

    return (
        reference_data[keep_mask],
        filtered_subjects,
        filtered_demographics,
        filtered_covariates,
    )

def _load_curvature_values(micapipe_directory, participant_id, session_id, hemi, resolution):
    bids_id = f"{participant_id}_{session_id}"
    curv_file = os.path.join(
        micapipe_directory,
        participant_id,
        session_id,
        "maps",
        f"{bids_id}_hemi-{hemi}_surf-fsLR-{resolution}_label-curv.func.gii",
    )
    if not os.path.exists(curv_file):
        return None
    img = nib.load(curv_file)
    if not img.darrays:
        return None
    return np.asarray(img.darrays[0].data, dtype=np.float32)

def load_subject_gyrification_covariates(
    micapipe_directory,
    participant_id,
    session_id,
    hemi,
    resolution,
    analysis="regional",
):
    """
    Return per-vertex curvature covariates from micapipe curvature maps.

    Regional maps use the curvature value from the requested hemisphere. For
    asymmetry maps, both left and right curvature values are included as two
    vertex-wise covariates so the W-score model can adjust for each subject's
    local folding pattern on both sides.
    """
    if analysis == "asymmetry":
        curv_lh = _load_curvature_values(micapipe_directory, participant_id, session_id, "L", resolution)
        curv_rh = _load_curvature_values(micapipe_directory, participant_id, session_id, "R", resolution)
        if curv_lh is None or curv_rh is None or curv_lh.shape != curv_rh.shape:
            return None
        return np.column_stack([curv_lh, curv_rh]).astype(np.float32)

    curv = _load_curvature_values(micapipe_directory, participant_id, session_id, hemi, resolution)
    if curv is None:
        return None
    return curv[:, None].astype(np.float32)

def load_reference_gyrification_covariates(
    micapipe_directory,
    subjects,
    hemi,
    resolution,
    analysis="regional",
    verbose=True,
):
    covariates = []
    loaded_subjects = []
    for participant_id, session_id in subjects:
        subject_covariates = load_subject_gyrification_covariates(
            micapipe_directory,
            participant_id,
            session_id,
            hemi,
            resolution,
            analysis=analysis,
        )
        if subject_covariates is None:
            if verbose:
                print(f"    Warning: Missing curvature map for gyrification covariates: {participant_id}/{session_id}")
            continue
        covariates.append(subject_covariates)
        loaded_subjects.append((participant_id, session_id))

    return (np.asarray(covariates, dtype=np.float32) if covariates else np.array([]), loaded_subjects)

def load_reference_surface_data(
    reference_subjects,
    output_directory,
    file_suffix,
    analysis='regional',
    verbose=True,
    data_transform=None,
    bilateral_data_transform=None,
    subject_data_transform=None,
    bilateral_subject_data_transform=None,
):
    """
    Load reference surface data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry maps
    verbose : bool, default=True
        If True, prints loading information
    data_transform : callable, optional
        Applied independently to each hemisphere after loading all GIFTI
        arrays and before regional/asymmetry construction.
    bilateral_data_transform : callable, optional
        Applied jointly to left/right arrays. Regional loading reads the
        companion hemisphere so normalization uses both hemispheres.
    subject_data_transform : callable, optional
        Called as ``transform(data, participant_id, session_id, hemi)`` before
        regional or asymmetry construction. This supports transforms requiring
        subject-specific geometry, such as physical gray-white distance.
    bilateral_subject_data_transform : callable, optional
        Called as ``transform(left, right, participant_id, session_id)``. This
        supports transforms requiring both hemispheres and subject geometry.
    
    Returns:
    --------
    tuple: (np.ndarray, list)
        - Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        ref_bids_id = f"{ref_pid}_{ref_sid}"
        
        if analysis == 'asymmetry':
            # For asymmetry, we need both left and right hemisphere data
            ref_file_lh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
            )
            ref_file_rh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
            )
            
            if os.path.exists(ref_file_lh) and os.path.exists(ref_file_rh):
                try:
                    ref_img_lh = nib.load(ref_file_lh)
                    ref_img_rh = nib.load(ref_file_rh)
                    
                    data_lh = np.column_stack(
                        [np.asarray(array.data) for array in ref_img_lh.darrays]
                    )
                    data_rh = np.column_stack(
                        [np.asarray(array.data) for array in ref_img_rh.darrays]
                    )
                    if bilateral_subject_data_transform is not None:
                        data_lh, data_rh = bilateral_subject_data_transform(
                            data_lh, data_rh, ref_pid, ref_sid
                        )
                    elif subject_data_transform is not None:
                        data_lh = subject_data_transform(
                            data_lh, ref_pid, ref_sid, "L"
                        )
                        data_rh = subject_data_transform(
                            data_rh, ref_pid, ref_sid, "R"
                        )
                    elif bilateral_data_transform is not None:
                        data_lh, data_rh = bilateral_data_transform(
                            data_lh, data_rh
                        )
                    elif data_transform is not None:
                        data_lh = data_transform(data_lh)
                        data_rh = data_transform(data_rh)
                    elif data_lh.shape[1] == 1:
                        data_lh = data_lh[:, 0]
                        data_rh = data_rh[:, 0]

                    asymmetry_data = compute_asymmetry(data_lh, data_rh)
                    reference_data.append(asymmetry_data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))
                        
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference files {ref_file_lh} or {ref_file_rh}: {e}")
        else:
            # Regional analysis - load single hemisphere data
            ref_file = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix}"
            )
            
            if os.path.exists(ref_file):
                try:
                    ref_img = nib.load(ref_file)
                    
                    data = np.column_stack(
                        [np.asarray(array.data) for array in ref_img.darrays]
                    )
                    if bilateral_subject_data_transform is not None:
                        other_suffix = (
                            file_suffix.replace("hemi-L", "hemi-R")
                            if "hemi-L" in file_suffix
                            else file_suffix.replace("hemi-R", "hemi-L")
                        )
                        other_file = os.path.join(
                            output_directory,
                            ref_pid,
                            ref_sid,
                            "maps",
                            "cortex",
                            f"{ref_bids_id}_{other_suffix}",
                        )
                        other_img = nib.load(other_file)
                        other_data = np.column_stack(
                            [np.asarray(array.data) for array in other_img.darrays]
                        )
                        if "hemi-L" in file_suffix:
                            data, _ = bilateral_subject_data_transform(
                                data, other_data, ref_pid, ref_sid
                            )
                        else:
                            _, data = bilateral_subject_data_transform(
                                other_data, data, ref_pid, ref_sid
                            )
                    elif subject_data_transform is not None:
                        subject_hemi = "L" if "hemi-L" in file_suffix else "R"
                        data = subject_data_transform(
                            data, ref_pid, ref_sid, subject_hemi
                        )
                    elif bilateral_data_transform is not None:
                        other_suffix = (
                            file_suffix.replace("hemi-L", "hemi-R")
                            if "hemi-L" in file_suffix
                            else file_suffix.replace("hemi-R", "hemi-L")
                        )
                        other_file = os.path.join(
                            output_directory,
                            ref_pid,
                            ref_sid,
                            "maps",
                            "cortex",
                            f"{ref_bids_id}_{other_suffix}",
                        )
                        other_img = nib.load(other_file)
                        other_data = np.column_stack(
                            [np.asarray(array.data) for array in other_img.darrays]
                        )
                        if "hemi-L" in file_suffix:
                            data, _ = bilateral_data_transform(data, other_data)
                        else:
                            _, data = bilateral_data_transform(other_data, data)
                    elif data_transform is not None:
                        data = data_transform(data)
                    elif data.shape[1] == 1:
                        data = data[:, 0]
                    reference_data.append(data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))
                        
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed cortical reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return np.array(reference_data) if reference_data else np.array([]), successfully_loaded_subjects

def load_reference_hippocampal_data(reference_subjects, output_directory, file_suffix, analysis='regional', verbose=True):
    """
    Load reference hippocampal data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry maps
    verbose : bool, default=True
        If True, prints loading information
    
    Returns:
    --------
    tuple: (np.ndarray, list)
        - Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
        For asymmetry analysis, returns left-right asymmetry maps
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        if analysis == 'asymmetry':
            # For asymmetry, we need both left and right hemisphere data
            ref_file_lh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
            )
            ref_file_rh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
            )
            
            if os.path.exists(ref_file_lh) and os.path.exists(ref_file_rh):
                try:
                    ref_img_lh = nib.load(ref_file_lh)
                    ref_img_rh = nib.load(ref_file_rh)
                    
                    data_lh = ref_img_lh.darrays[0].data
                    data_rh = ref_img_rh.darrays[0].data
                    
                    # Compute asymmetry
                    asymmetry_data = compute_asymmetry(data_lh, data_rh)
                    reference_data.append(asymmetry_data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                    
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference files {ref_file_lh} or {ref_file_rh}: {e}")
        else:
            # Regional analysis - load single hemisphere data
            ref_file = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix}"
            )
            
            if os.path.exists(ref_file):
                try:
                    ref_img = nib.load(ref_file)
                    reference_data.append(ref_img.darrays[0].data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed hippocampal reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return np.array(reference_data) if reference_data else np.array([]), successfully_loaded_subjects

def load_reference_subcortical_data(reference_subjects, output_directory, file_suffix, analysis='regional', verbose=True):
    """
    Load reference subcortical data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry for subcortical structures
    verbose : bool, default=True
        If True, prints loading information
    
    Returns:
    --------
    tuple: (pd.DataFrame, list)
        - Combined reference data from all subjects
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
        For asymmetry analysis, returns asymmetry values for paired structures
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        ref_bids_id = f"{ref_pid}_{ref_sid}"
        ref_file = os.path.join(
            output_directory,
            ref_pid, ref_sid, "maps", "subcortical",
            f"{ref_bids_id}_{file_suffix}"
        )
        
        if os.path.exists(ref_file):
            try:
                ref_df = pd.read_csv(ref_file)
                
                if analysis == 'asymmetry':
                    # Compute asymmetry for paired structures
                    asymmetry_data = {}
                    for col in ref_df.columns:
                        if col.startswith('L'):
                            right_col = col.replace('L', 'R')
                            if right_col in ref_df.columns:
                                left_val = ref_df[col].iloc[0]
                                right_val = ref_df[right_col].iloc[0]
                                avg_val = (left_val + right_val) / 2
                                asym_val = (left_val - right_val) / avg_val if avg_val > 0 else 0
                                asymmetry_data[col] = asym_val
                    
                    # Convert to DataFrame
                    ref_df = pd.DataFrame([asymmetry_data])
                
                reference_data.append(ref_df)
                successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                
            except Exception as e:
                if verbose:
                    print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed subcortical reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return (pd.concat(reference_data, ignore_index=True) if reference_data else pd.DataFrame(),
            successfully_loaded_subjects)

def analyze_dataset(
    dataset,
    reference,
    method='zscore',
    output_directory=None,
    verbose=True,
    gyrification_matching=None,
    use_curvature_covariates=True,
    predictive_wscore=False,
    wscore_distribution="gaussian",
    wscore_preprocessing="none",
    wscore_covariate_model="linear",
    wscore_surface_smoothing_iterations=0,
    blur_depth_model="mean_slope_rms",
    intensity_depth_model="raw",
    quant_sample_surface=None,
    control_correlation_filter=False,
    control_correlation_quantile=None,
    prediction_variance_percentile=None,
    export_control_features=None,
    site_harmonizer=None,
    site=None,
    scoring_site_covariate=False,
):
    """
    Analyze a dataset by comparing it to a reference dataset using specified method.
    
    Parameters:
    -----------
    dataset : zbdataset
        Dataset to analyze (patient dataset)
    reference : zbdataset
        Reference dataset to compare against (e.g., control dataset)
    method : str, default='zscore'
        Analysis method to use. Options: 'zscore', 'wscore'
    normative_columns : list, optional
        List of demographic columns to use for normative modeling (only used for wscore)
        Default: ['age', 'sex']
    output_directory : str, optional
        Directory where score maps will be saved. If None, uses validation output directory
    verbose : bool, default=True
        If True, prints detailed information about the analysis process
    use_curvature_covariates : bool, default=True
        For cortical surfaces, use micapipe curvature maps during scoring.
        W-score models include each subject's vertex-wise curvature as an
        additional regression covariate. Z-score maps residualize each vertex
        against the same per-subject curvature covariates before scoring.
    gyrification_matching : bool, optional
        Deprecated alias for use_curvature_covariates.
    predictive_wscore : bool, default=False
        If True, W-score maps divide by the predictive standard deviation for
        a new patient rather than only the residual SD from controls.
    wscore_distribution : str, default="gaussian"
        Distribution fitted to control residuals for W-scoring.
    wscore_preprocessing : {"none", "spatial_zscore", "spatial_robust_z"}, default="none"
        Optional within-subject robust spatial normalization.
    wscore_covariate_model : str, default="linear"
        Optional quadratic-age and/or age-by-sex demographic design.
    wscore_surface_smoothing_iterations : int, default=0
        Post-score one-ring smoothing steps for fsLR-32k cortex maps.
    blur_depth_model : {"mean_slope_rms", "mean", "gradient_flattening"}, default="mean_slope_rms"
        Multi-depth blur reduction used for cortical W-score maps.
    intensity_depth_model : {"raw", "white_swm1_direction_cosine", "multisurface_median_abs_dominant"}, default="raw"
        For fsLR-32k T1w/FLAIR white-surface analysis, optionally replace raw
        intensity with stabilized white-minus-SWM1 depth contrast.
    control_correlation_filter : bool, default=False
        If True, each feature/map drops the bottom ``control_correlation_quantile``
        fraction of controls by their average correlation to the other controls
        for that feature/map (rank-based, adaptive per feature).
    control_correlation_quantile : float or dict or None, default=None
        Drop-fraction for the rank-based control filter. Scalar (same fraction for
        every feature), per-feature mapping (absent/None -> off for that feature),
        or None (off).

    Returns:
    --------
    dict
        Dictionary containing analysis results for each feature and region
    """
    if method not in ['zscore', 'wscore']:
        raise ValueError(f"Method '{method}' not supported. Currently supports 'zscore' and 'wscore'.")

    if gyrification_matching is not None:
        use_curvature_covariates = bool(gyrification_matching)
    else:
        use_curvature_covariates = bool(use_curvature_covariates)
    control_correlation_filter = bool(control_correlation_filter)
    # ``control_correlation_quantile`` is the RANK-based drop-fraction (drop the
    # bottom fraction of controls by correlation, per feature). Scalar (one fraction
    # for every feature), per-feature mapping, or ``None`` (off). Coerce values.
    if isinstance(control_correlation_quantile, Mapping):
        control_correlation_quantile = {
            str(name): (None if value is None else float(value))
            for name, value in control_correlation_quantile.items()
        }
    elif control_correlation_quantile is not None:
        control_correlation_quantile = float(control_correlation_quantile)
    wscore_distribution = _normalize_wscore_distribution(wscore_distribution)
    wscore_preprocessing = _normalize_wscore_preprocessing(wscore_preprocessing)
    wscore_covariate_model = _normalize_wscore_covariate_model(wscore_covariate_model)
    blur_depth_model = _normalize_blur_depth_model(blur_depth_model)
    intensity_depth_model = _normalize_intensity_depth_model(intensity_depth_model)
    wscore_surface_smoothing_iterations = int(wscore_surface_smoothing_iterations)
    if wscore_surface_smoothing_iterations < 0:
        raise ValueError("wscore_surface_smoothing_iterations must be non-negative")
    wscore_fit_cache = {}
    
    # Check that features match
    if set(dataset.features) != set(reference.features):
        raise ValueError("Feature sets don't match between patient and reference datasets.")
    if intensity_depth_model != "raw":
        lower_features = {str(feature).lower() for feature in dataset.features}
        missing_companions = [
            f"{base}*blur"
            for base in ("t1w", "flair")
            if base in lower_features and f"{base}*blur" not in lower_features
        ]
        if missing_companions:
            raise ValueError(
                f"intensity_depth_model='{intensity_depth_model}' requires companion "
                "features: " + ", ".join(missing_companions)
            )
    if quant_sample_surface not in (None, "none", "raw"):
        lower_features = {str(feature).lower() for feature in dataset.features}
        missing_quant = [
            f"{base}*blur"
            for base in ("adc", "fa", "qt1")
            if base in lower_features and f"{base}*blur" not in lower_features
        ]
        if missing_quant:
            raise ValueError(
                f"quant_sample_surface='{quant_sample_surface}' requires companion "
                "features: " + ", ".join(missing_quant)
            )

    # Check if demographics data is available for wscore
    if method == 'wscore':
        if not hasattr(dataset, 'demographics') or not hasattr(reference, 'demographics'):
            raise ValueError("Demographics data required for both datasets when using wscore method.")
        
        if not hasattr(dataset.demographics, 'normative_columns') or not hasattr(reference.demographics, 'normative_columns'):
            raise ValueError("Normative columns must be specified in dataset and reference demographics for wscore method.")
        
        # Check if required columns exist
        for col in dataset.demographics.normative_columns:
            if col not in dataset.demographics.data.columns:
                raise ValueError(f"Column '{col}' not found in dataset demographics.")
            if col not in reference.demographics.data.columns:
                raise ValueError(f"Column '{col}' not found in reference demographics.")
    
    # Use default output directory if not provided
    if output_directory is None:
        raise ValueError("output_directory must be specified for saving score maps.")
    
    if verbose:
        print(f"Analyzing dataset {dataset.name} against reference {reference.name} using {method} method...")
        print(f"Features to analyze: {', '.join(dataset.features)}")
        if dataset.cortex and reference.cortex:
            print(f"Per-vertex curvature covariates: {'on' if use_curvature_covariates else 'off'}")
        if method == 'wscore':
            print(f"Normative columns: {', '.join(dataset.demographics.normative_columns)}")
            print(f"W-score residual distribution: {wscore_distribution}")
            print(f"W-score preprocessing: {wscore_preprocessing}")
            print(f"W-score covariate model: {wscore_covariate_model}")
            print(
                "Cortical W-score smoothing iterations: "
                f"{wscore_surface_smoothing_iterations}"
            )
            print(f"Blur depth model: {blur_depth_model}")
            print(f"T1w/FLAIR intensity depth model: {intensity_depth_model}")
            print(f"Predictive W-score uncertainty: {'on' if predictive_wscore else 'off'}")
        if not control_correlation_filter:
            print("Control correlation filter: off")
        elif isinstance(control_correlation_quantile, Mapping):
            desc = ", ".join(
                f"{name}={'off' if value is None else format(float(value), '.0%')}"
                for name, value in sorted(control_correlation_quantile.items())
            )
            print(f"Control correlation filter: on (per-feature drop-bottom: {desc})")
        elif control_correlation_quantile is None:
            print("Control correlation filter: on (keep-all; no quantile set)")
        else:
            print(
                "Control correlation filter: on "
                f"(drop-bottom {control_correlation_quantile:.0%} per feature)"
            )
    
    # Get valid subjects from both datasets
    if not hasattr(dataset, 'valid_subjects'):
        raise ValueError("Dataset does not have valid_subjects attribute. Please validate the dataset first.")
    if not hasattr(reference, 'valid_subjects'):
        raise ValueError("Reference dataset does not have valid_subjects attribute. Please validate the reference dataset first.")
    
    # Print summary of valid subjects
    if verbose:
        print(f"Patient dataset: {len(dataset.valid_subjects['base'])} base subjects")
        print(f"Reference dataset: {len(reference.valid_subjects['base'])} base subjects")
        
        # Print feature-specific counts
        for feature in dataset.features:
            if feature in dataset.valid_subjects:
                pat_count = len(dataset.valid_subjects[feature]['all']) if 'all' in dataset.valid_subjects[feature] else 0
                ref_count = len(reference.valid_subjects[feature]['all']) if feature in reference.valid_subjects and 'all' in reference.valid_subjects[feature] else 0
                print(f"Feature {feature}: {pat_count} patient subjects, {ref_count} reference subjects")
    
    # Feature mapping for different naming conventions
    feature_mapping = {
        "thickness": {"output": "thickness"},
        "flair": {"output": "FLAIR"},
        "adc": {"output": "ADC"},
        "fa": {"output": "FA"},
        "qt1": {"output": "qT1"},  # Changed from "T1map" to "qT1"
        "fmri": {"output": "fMRI"},
        "t1w": {"output": "T1w"},
        "sa": {"output": "SA"}
    }
    
    # Define analysis parameters
    resolutions = ["32k", "5k"]
    labels = ["midthickness", "white"]
    hemispheres = ["L", "R"]
    
    # Store analysis results
    analysis_results = {
        "cortical": {},
        "hippocampal": {},
        "subcortical": {},
        "blur": {}
    }
    
    # Get demographics data for wscore. Cross-site harmonization (fit/apply) and
    # the export pass ALSO need AGE/SEX per subject even when method='zscore', so
    # build the maps whenever harmonizing/exporting too.
    _need_demographics = (method == 'wscore'
                          or site_harmonizer is not None
                          or export_control_features is not None)
    ref_demo_dict = {}
    pat_demo_dict = {}
    if _need_demographics:
        for _, row in reference.demographics.data.iterrows():
            key = (row['participant_id'], row.get('session_id', 'ses-01'))
            ref_demo_dict[key] = row
        for _, row in dataset.demographics.data.iterrows():
            key = (row['participant_id'], row.get('session_id', 'ses-01'))
            pat_demo_dict[key] = row

    # Get normative columns for W-score calculation
    normative_columns = dataset.demographics.normative_columns if method == 'wscore' else None
    # site_covariate arm: add a binary 'site' column to the normative design.
    if scoring_site_covariate and normative_columns is not None:
        normative_columns = list(normative_columns) + ['site']
    
    # Process each feature
    for feature in dataset.features:
        if verbose:
            print(f"\nAnalyzing feature: {feature}")
        
        feat_lower = feature.lower()
        is_blur = feat_lower.endswith("*blur")
        is_fmri = feat_lower == "fmri"
        # Get base feature name and mapping
        if is_blur:
            base_feat = feat_lower.replace("*blur", "")
            if base_feat in feature_mapping:
                output_feat = feature_mapping[base_feat]["output"] + "*blur"
            else:
                output_feat = base_feat + "*blur"
        else:
            if feat_lower in feature_mapping:
                output_feat = feature_mapping[feat_lower]["output"]
            else:
                output_feat = feat_lower
        
        # 1. CORTICAL ANALYSIS
        if dataset.cortex and reference.cortex and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing cortical data for {feature}...")
            
            # Get valid subjects for this feature and structure
            dataset_cortical_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_cortical_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_cortical_subjects or not reference_cortical_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for cortical {feature}")
                continue
                
            cortical_results = {}
            cortical_regional_qc_subjects = {}
            
            for analysis in ['regional', 'asymmetry']:
                for hemi in hemispheres:
                    for resolution in resolutions:
                        for label in labels:
                            map_key = f"{hemi}_{resolution}_{label}_{analysis}"
                            
                            file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"

                            # Depth-model dispatch: T1w/FLAIR use intensity_depth_model
                            # (sampling + optional self-norm); ADC/FA/qT1 use
                            # quant_sample_surface (PLAIN sample only, no self-norm).
                            # Both consume the 4-depth *blur companion and yield one
                            # per-vertex value at the white surface (32k).
                            depth_model = None
                            if resolution == "32k" and label == "white":
                                if feat_lower in {"t1w", "flair"} and intensity_depth_model != "raw":
                                    depth_model = intensity_depth_model
                                elif (feat_lower in {"adc", "fa", "qt1"}
                                      and quant_sample_surface not in (None, "none", "raw")):
                                    depth_model = f"sample_{quant_sample_surface}"
                            use_intensity_depth_metric = depth_model is not None
                            source_file_suffix = file_suffix
                            surface_data_transform = None
                            bilateral_surface_data_transform = None
                            active_intensity_depth_model = "raw"
                            if use_intensity_depth_metric:
                                source_file_suffix = (
                                    f"hemi-{hemi}_surf-fsLR-32k_label-midthickness_"
                                    f"feature-{output_feat}*blur_smooth-"
                                    f"{dataset.cortical_smoothing}mm.func.gii"
                                )
                                if depth_model == "multisurface_median_abs_dominant":
                                    bilateral_surface_data_transform = (
                                        _bilateral_multisurface_median_abs_scale
                                    )
                                else:
                                    surface_data_transform = _INTENSITY_DEPTH_TRANSFORMS[
                                        depth_model
                                    ]
                                active_intensity_depth_model = depth_model
                                if verbose:
                                    print(
                                        f"    Using {depth_model} for "
                                        f"{feature} {map_key}"
                                    )
                            
                            # Load reference data and track which subjects actually loaded
                            reference_data, successfully_loaded_subjects = load_reference_surface_data(
                                reference_cortical_subjects,
                                output_directory,
                                source_file_suffix,
                                analysis,
                                verbose,
                                data_transform=surface_data_transform,
                                bilateral_data_transform=bilateral_surface_data_transform,
                            )
                            
                            if len(reference_data) == 0:
                                if verbose:
                                    print(f"    Warning: No reference data found for {map_key}")
                                continue

                            reference_gyrification_covariates = None
                            if use_curvature_covariates:
                                reference_gyrification_covariates, gyr_subjects = load_reference_gyrification_covariates(
                                    reference.micapipe_directory,
                                    successfully_loaded_subjects,
                                    hemi,
                                    resolution,
                                    analysis=analysis,
                                    verbose=verbose,
                                )
                                if len(gyr_subjects) == 0:
                                    if verbose:
                                        print(f"    Warning: No curvature data found for gyrification adjustment in {map_key}")
                                    continue
                                gyr_index = {subject: i for i, subject in enumerate(gyr_subjects)}
                                ref_indices = [
                                    i for i, subject in enumerate(successfully_loaded_subjects)
                                    if subject in gyr_index
                                ]
                                reference_data = reference_data[ref_indices]
                                successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]
                                reference_gyrification_covariates = np.asarray(
                                    [reference_gyrification_covariates[gyr_index[subject]] for subject in successfully_loaded_subjects],
                                    dtype=np.float32,
                                )
                            
                            ref_demographics_df = None
                            # Get demographics for reference subjects (wscore, or when
                            # harmonizing/exporting -- the harmonizer needs AGE/SEX).
                            if _need_demographics:
                                ref_demographics = []
                                valid_ref_subjects = []
                                # Now iterate over subjects that ACTUALLY loaded, not all subjects
                                for ref_pid, ref_sid in successfully_loaded_subjects:
                                    key = (ref_pid, ref_sid)
                                    if key in ref_demo_dict:
                                        ref_demographics.append(ref_demo_dict[key])
                                        valid_ref_subjects.append(key)
                                
                                if len(ref_demographics) == 0:
                                    if verbose:
                                        print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                                    continue
                                
                                ref_demographics_df = pd.DataFrame(ref_demographics)
                                # Filter reference data to match demographics
                                # Now indices will match because we're iterating over successfully_loaded_subjects
                                ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                                reference_data = reference_data[ref_indices]
                                successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]
                                if reference_gyrification_covariates is not None:
                                    reference_gyrification_covariates = reference_gyrification_covariates[ref_indices]

                            if control_correlation_filter:
                                if analysis == 'regional':
                                    (
                                        reference_data,
                                        successfully_loaded_subjects,
                                        ref_demographics_df,
                                        reference_gyrification_covariates,
                                    ) = _filter_reference_controls_by_correlation(
                                        reference_data,
                                        successfully_loaded_subjects,
                                        reference_demographics_df=ref_demographics_df,
                                        reference_vertex_covariates=reference_gyrification_covariates,
                                        quantile=_resolve_control_correlation_quantile(
                                            control_correlation_quantile, feature
                                        ),
                                        verbose=verbose,
                                        context=f"cortical {feature} {map_key}",
                                    )
                                    cortical_regional_qc_subjects[(resolution, label, hemi)] = set(
                                        successfully_loaded_subjects
                                    )
                                else:
                                    allowed_subjects = (
                                        cortical_regional_qc_subjects.get((resolution, label, 'L'), set())
                                        & cortical_regional_qc_subjects.get((resolution, label, 'R'), set())
                                    )
                                    (
                                        reference_data,
                                        successfully_loaded_subjects,
                                        ref_demographics_df,
                                        reference_gyrification_covariates,
                                    ) = _filter_reference_controls_to_subjects(
                                        reference_data,
                                        successfully_loaded_subjects,
                                        allowed_subjects,
                                        reference_demographics_df=ref_demographics_df,
                                        reference_vertex_covariates=reference_gyrification_covariates,
                                        verbose=verbose,
                                        context=f"cortical {feature} {map_key}",
                                    )
                                if len(reference_data) < 2:
                                    if verbose:
                                        print(
                                            f"    Warning: fewer than 2 controls remain after "
                                            f"correlation filtering for {map_key}; skipping"
                                        )
                                    continue

                            # --- Cross-site harmonization hook (analysis-level) ---
                            _hk = f"{feature}|{map_key}"
                            if export_control_features is not None:
                                # Pass 1: stash this site's control features; skip scoring.
                                export_control_features[_hk] = (
                                    reference_data, successfully_loaded_subjects, ref_demographics_df)
                                continue
                            if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                # Pass 2: replace the local reference with the POOLED
                                # (both-sites) harmonized control reference.
                                _bundle = site_harmonizer.pooled_reference(_hk)
                                if _bundle is not None:
                                    reference_data, successfully_loaded_subjects, ref_demographics_df = _bundle
                                    reference_gyrification_covariates = None

                            # Process patient data
                            patient_scores = []
                            for pat_pid, pat_sid in dataset_cortical_subjects:
                                pat_bids_id = f"{pat_pid}_{pat_sid}"

                                if analysis == 'asymmetry':
                                    # For asymmetry, we need both left and right hemisphere data
                                    pat_file_lh = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{source_file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
                                    )
                                    pat_file_rh = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{source_file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                                    )
                                    
                                    if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                        try:
                                            pat_img_lh = nib.load(pat_file_lh)
                                            pat_img_rh = nib.load(pat_file_rh)

                                            if bilateral_surface_data_transform is not None:
                                                pat_depths_lh = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img_lh.darrays]
                                                )
                                                pat_depths_rh = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img_rh.darrays]
                                                )
                                                pat_data_lh, pat_data_rh = (
                                                    bilateral_surface_data_transform(
                                                        pat_depths_lh, pat_depths_rh
                                                    )
                                                )
                                                pat_data = compute_asymmetry(
                                                    pat_data_lh, pat_data_rh
                                                )
                                            elif surface_data_transform is not None:
                                                pat_depths_lh = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img_lh.darrays]
                                                )
                                                pat_depths_rh = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img_rh.darrays]
                                                )
                                                pat_data_lh = surface_data_transform(pat_depths_lh)
                                                pat_data_rh = surface_data_transform(pat_depths_rh)
                                                pat_data = compute_asymmetry(
                                                    pat_data_lh, pat_data_rh
                                                )
                                            
                                            # Handle multi-depth data (blur features)
                                            elif len(pat_img_lh.darrays) > 1:
                                                pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                                pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                                
                                                for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                                    pat_data_lh[:, e] = darray_lh.data
                                                    pat_data_rh[:, e] = darray_rh.data
                                                
                                                # Compute asymmetry
                                                pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                            else:
                                                pat_data_lh = pat_img_lh.darrays[0].data
                                                pat_data_rh = pat_img_rh.darrays[0].data
                                                
                                                # Compute asymmetry
                                                pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)

                                            # Cross-site: harmonize this patient into the
                                            # pooled reference space with FROZEN params.
                                            if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                                _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                                if _pd is not None:
                                                    pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                            reference_vertex_covariates = None
                                            patient_vertex_covariates = None
                                            if use_curvature_covariates:
                                                patient_vertex_covariates = load_subject_gyrification_covariates(
                                                    dataset.micapipe_directory,
                                                    pat_pid,
                                                    pat_sid,
                                                    hemi,
                                                    resolution,
                                                    analysis=analysis,
                                                )
                                                if patient_vertex_covariates is None:
                                                    if verbose:
                                                        print(f"    Warning: Missing curvature map for patient {pat_pid}/{pat_sid}; skipping gyrification-adjusted scoring")
                                                    continue
                                                reference_vertex_covariates = reference_gyrification_covariates

                                            # Create score output file with analysis type in filename
                                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                            score_file = os.path.join(
                                                score_dir, 
                                                f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                            )
                                            
                                            # Calculate scores
                                            if method == 'zscore':
                                                score_result = calculate_zscore_maps(
                                                    reference_data, pat_data, score_file, analysis, verbose,
                                                    reference_vertex_covariates=reference_vertex_covariates,
                                                    patient_vertex_covariates=patient_vertex_covariates,
                                                    wscore_preprocessing=wscore_preprocessing,
                                                )
                                            else:  # wscore
                                                # Get patient demographics
                                                pat_key = (pat_pid, pat_sid)
                                                if pat_key not in pat_demo_dict:
                                                    if verbose:
                                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                    continue
                                                
                                                pat_demographics = pat_demo_dict[pat_key]
                                                score_result = calculate_wscore_maps(
                                                    reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                    score_file, normative_columns, verbose,
                                                    reference_vertex_covariates=reference_vertex_covariates,
                                                    patient_vertex_covariates=patient_vertex_covariates,
                                                    use_prediction_uncertainty=predictive_wscore,
                                                    wscore_distribution=wscore_distribution,
                                                    wscore_preprocessing=wscore_preprocessing,
                                                    wscore_covariate_model=wscore_covariate_model,
                                                    wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                                    blur_depth_model=blur_depth_model,
                                                    intensity_depth_model=active_intensity_depth_model,
                                                    wscore_fit_cache=wscore_fit_cache,
                                                    prediction_variance_percentile=prediction_variance_percentile,
                                                )
                                            
                                            score_result['subject'] = (pat_pid, pat_sid)
                                            score_result['analysis'] = analysis
                                            patient_scores.append(score_result)
                                            
                                        except Exception as e:
                                            if verbose:
                                                print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                                    else:
                                        if verbose:
                                            print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                                
                                else:
                                    # Regional analysis - load single hemisphere data
                                    pat_file = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{source_file_suffix}"
                                    )
                                    
                                    if os.path.exists(pat_file):
                                        try:
                                            pat_img = nib.load(pat_file)
                                            if bilateral_surface_data_transform is not None:
                                                other_source_suffix = (
                                                    source_file_suffix.replace("hemi-L", "hemi-R")
                                                    if "hemi-L" in source_file_suffix
                                                    else source_file_suffix.replace("hemi-R", "hemi-L")
                                                )
                                                other_pat_file = os.path.join(
                                                    output_directory,
                                                    pat_pid,
                                                    pat_sid,
                                                    "maps",
                                                    "cortex",
                                                    f"{pat_bids_id}_{other_source_suffix}",
                                                )
                                                other_pat_img = nib.load(other_pat_file)
                                                pat_depths = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img.darrays]
                                                )
                                                other_depths = np.column_stack(
                                                    [np.asarray(array.data) for array in other_pat_img.darrays]
                                                )
                                                if "hemi-L" in source_file_suffix:
                                                    pat_data, _ = bilateral_surface_data_transform(
                                                        pat_depths, other_depths
                                                    )
                                                else:
                                                    _, pat_data = bilateral_surface_data_transform(
                                                        other_depths, pat_depths
                                                    )
                                            elif surface_data_transform is not None:
                                                pat_depths = np.column_stack(
                                                    [np.asarray(array.data) for array in pat_img.darrays]
                                                )
                                                pat_data = surface_data_transform(pat_depths)
                                            else:
                                                pat_data = pat_img.darrays[0].data

                                            # Cross-site: harmonize this patient into the
                                            # pooled reference space with FROZEN params.
                                            if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                                _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                                if _pd is not None:
                                                    pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                            reference_vertex_covariates = None
                                            patient_vertex_covariates = None
                                            if use_curvature_covariates:
                                                patient_vertex_covariates = load_subject_gyrification_covariates(
                                                    dataset.micapipe_directory,
                                                    pat_pid,
                                                    pat_sid,
                                                    hemi,
                                                    resolution,
                                                    analysis=analysis,
                                                )
                                                if patient_vertex_covariates is None:
                                                    if verbose:
                                                        print(f"    Warning: Missing curvature map for patient {pat_pid}/{pat_sid}; skipping gyrification-adjusted scoring")
                                                    continue
                                                reference_vertex_covariates = reference_gyrification_covariates

                                            # Create score output file with analysis type in filename
                                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                            score_file = os.path.join(
                                                score_dir, 
                                                f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                            )
                                            
                                            # Calculate scores
                                            if method == 'zscore':
                                                score_result = calculate_zscore_maps(
                                                    reference_data, pat_data, score_file, analysis, verbose,
                                                    reference_vertex_covariates=reference_vertex_covariates,
                                                    patient_vertex_covariates=patient_vertex_covariates,
                                                    wscore_preprocessing=wscore_preprocessing,
                                                )
                                            else:  # wscore
                                                # Get patient demographics
                                                pat_key = (pat_pid, pat_sid)
                                                if pat_key not in pat_demo_dict:
                                                    if verbose:
                                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                    continue
                                                
                                                pat_demographics = pat_demo_dict[pat_key]
                                                score_result = calculate_wscore_maps(
                                                    reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                    score_file, normative_columns, verbose,
                                                    reference_vertex_covariates=reference_vertex_covariates,
                                                    patient_vertex_covariates=patient_vertex_covariates,
                                                    use_prediction_uncertainty=predictive_wscore,
                                                    wscore_distribution=wscore_distribution,
                                                    wscore_preprocessing=wscore_preprocessing,
                                                    wscore_covariate_model=wscore_covariate_model,
                                                    wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                                    blur_depth_model=blur_depth_model,
                                                    intensity_depth_model=active_intensity_depth_model,
                                                    wscore_fit_cache=wscore_fit_cache,
                                                    prediction_variance_percentile=prediction_variance_percentile,
                                                )
                                            
                                            score_result['subject'] = (pat_pid, pat_sid)
                                            score_result['analysis'] = analysis
                                            patient_scores.append(score_result)
                                            
                                        except Exception as e:
                                            if verbose:
                                                print(f"    Warning: Could not process patient file {pat_file}: {e}")
                            
                            cortical_results[map_key] = {
                                f'patient_{method}s': patient_scores
                            }

            analysis_results["cortical"][feature] = cortical_results
        
        # 2. HIPPOCAMPAL ANALYSIS
        if dataset.hippocampus and reference.hippocampus and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing hippocampal data for {feature}...")
            
            # Determine resolution string from dataset
            hipp_res = "8k" if getattr(dataset, "hippunfold_version", 1) == 2 else "0p5mm"
            
            # Get valid subjects for this feature and structure
            dataset_hippocampal_subjects = dataset.valid_subjects[feature]['structures']['hippocampus'] if feature in dataset.valid_subjects else []
            reference_hippocampal_subjects = reference.valid_subjects[feature]['structures']['hippocampus'] if feature in reference.valid_subjects else []
            
            if not dataset_hippocampal_subjects or not reference_hippocampal_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for hippocampal {feature}")
                continue
                
            hippocampal_results = {}
            hippocampal_regional_qc_subjects = {}
            
            for analysis in ['regional', 'asymmetry']:
                for hemi in hemispheres:
                    map_key = f"{hemi}_hippocampus_{analysis}"  # Include analysis in map_key
                    file_suffix = f"hemi-{hemi}_den-{hipp_res}_label-hipp_midthickness_feature-{output_feat}_smooth-{dataset.hippocampal_smoothing}mm.func.gii"
                    
                    # Load reference data
                    reference_data, successfully_loaded_subjects = load_reference_hippocampal_data(
                        reference_hippocampal_subjects, output_directory, file_suffix, analysis, verbose
                    )
                    
                    if len(reference_data) == 0:
                        if verbose:
                            print(f"    Warning: No reference data found for {map_key}")
                        continue
                    
                    ref_demographics_df = None
                    # Get demographics for reference subjects (wscore, or when
                    # harmonizing/exporting -- the harmonizer needs AGE/SEX).
                    if _need_demographics:
                        ref_demographics = []
                        valid_ref_subjects = []
                        for ref_pid, ref_sid in successfully_loaded_subjects:
                            key = (ref_pid, ref_sid)
                            if key in ref_demo_dict:
                                ref_demographics.append(ref_demo_dict[key])
                                valid_ref_subjects.append(key)

                        if len(ref_demographics) == 0:
                            if verbose:
                                print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                            continue

                        ref_demographics_df = pd.DataFrame(ref_demographics)
                        # Filter reference data to match demographics
                        ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                        reference_data = reference_data[ref_indices]
                        successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]

                    if control_correlation_filter:
                        if analysis == 'regional':
                            (
                                reference_data,
                                successfully_loaded_subjects,
                                ref_demographics_df,
                                _,
                            ) = _filter_reference_controls_by_correlation(
                                reference_data,
                                successfully_loaded_subjects,
                                reference_demographics_df=ref_demographics_df,
                                quantile=_resolve_control_correlation_quantile(
                                    control_correlation_quantile, feature
                                ),
                                verbose=verbose,
                                context=f"hippocampal {feature} {map_key}",
                            )
                            hippocampal_regional_qc_subjects[hemi] = set(
                                successfully_loaded_subjects
                            )
                        else:
                            allowed_subjects = (
                                hippocampal_regional_qc_subjects.get('L', set())
                                & hippocampal_regional_qc_subjects.get('R', set())
                            )
                            (
                                reference_data,
                                successfully_loaded_subjects,
                                ref_demographics_df,
                                _,
                            ) = _filter_reference_controls_to_subjects(
                                reference_data,
                                successfully_loaded_subjects,
                                allowed_subjects,
                                reference_demographics_df=ref_demographics_df,
                                verbose=verbose,
                                context=f"hippocampal {feature} {map_key}",
                            )
                        if len(reference_data) < 2:
                            if verbose:
                                print(
                                    f"    Warning: fewer than 2 controls remain after "
                                    f"correlation filtering for {map_key}; skipping"
                                )
                            continue
                    
                    # --- Cross-site harmonization hook (analysis-level) ---
                    _hk = f"{feature}|{map_key}"
                    if export_control_features is not None:
                        export_control_features[_hk] = (
                            reference_data, successfully_loaded_subjects, ref_demographics_df)
                        continue
                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                        _bundle = site_harmonizer.pooled_reference(_hk)
                        if _bundle is not None:
                            reference_data, successfully_loaded_subjects, ref_demographics_df = _bundle

                    # Process patient data
                    patient_scores = []
                    for pat_pid, pat_sid in dataset_hippocampal_subjects:
                        if analysis == 'asymmetry':
                            # For asymmetry, we need both left and right hemisphere data
                            pat_file_lh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"

                            )
                            pat_file_rh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                            )
                            
                            if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                try:
                                    pat_img_lh = nib.load(pat_file_lh)
                                    pat_img_rh = nib.load(pat_file_rh)
                                    
                                    pat_data_lh = pat_img_lh.darrays[0].data
                                    pat_data_rh = pat_img_rh.darrays[0].data
                                    
                                    # Compute asymmetry
                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)

                                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                        _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                        if _pd is not None:
                                            pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "hippocampus")
                                    score_file = os.path.join(
                                        score_dir, 
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose,
                                            use_prediction_uncertainty=predictive_wscore,
                                            wscore_distribution=wscore_distribution,
                                            wscore_preprocessing=wscore_preprocessing,
                                            wscore_covariate_model=wscore_covariate_model,
                                            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                            blur_depth_model=blur_depth_model,
                                            wscore_fit_cache=wscore_fit_cache,
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                            else:
                                if verbose:
                                    print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                        
                        else:
                            # Regional analysis - load single hemisphere data
                            pat_file = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix}"
                            )
                            
                            if os.path.exists(pat_file):
                                try:
                                    pat_img = nib.load(pat_file)
                                    pat_data = pat_img.darrays[0].data

                                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                        _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                        if _pd is not None:
                                            pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "hippocampus")
                                    score_file = os.path.join(
                                        score_dir, 
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose,
                                            use_prediction_uncertainty=predictive_wscore,
                                            wscore_distribution=wscore_distribution,
                                            wscore_preprocessing=wscore_preprocessing,
                                            wscore_covariate_model=wscore_covariate_model,
                                            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                            blur_depth_model=blur_depth_model,
                                            wscore_fit_cache=wscore_fit_cache,
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient file {pat_file}: {e}")
                    
                    hippocampal_results[map_key] = {
                        f'patient_{method}s': patient_scores
                    }
            
            analysis_results["hippocampal"][feature] = hippocampal_results
        
        # 3. BLUR ANALYSIS
        if is_blur:
            if verbose:
                print(f"  Processing blur depth data for {feature}...")
            
            # For blur features, use cortical subjects from the base feature
            base_feature = feat_lower.replace("*blur", "")
            dataset_blur_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_blur_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_blur_subjects or not reference_blur_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for blur {feature}")
                continue
                
            blur_results = {}
            blur_regional_qc_subjects = {}

            blur_bilateral_subject_data_transform = None
            blur_subject_data_transform = None
            if blur_depth_model == "gradient_flattening":
                def blur_bilateral_subject_data_transform(left, right, pid, sid):
                    return _load_bilateral_blur_gradient_flattening_components(
                        left,
                        right,
                        reference.micapipe_directory,
                        pid,
                        sid,
                        resolution="32k",
                    )
            elif blur_depth_model in ("gray_white_contrast", "boundary_gradient", "juxtacortical_gradient"):
                def blur_subject_data_transform(data, pid, sid, transform_hemi):
                    return _load_blur_single_gradient(
                        data, reference.micapipe_directory, pid, sid,
                        transform_hemi, blur_depth_model, resolution="32k",
                    )

            for analysis in ['regional', 'asymmetry']:
                for hemi in hemispheres:
                    map_key = f"{hemi}_blur_depths_{analysis}"  # Include analysis in map_key
                    label = "midthickness"
                    resolution = "32k"
                    
                    file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"
                    
                    # Load reference data
                    reference_data, successfully_loaded_subjects = load_reference_surface_data(
                        reference_blur_subjects,
                        output_directory,
                        file_suffix,
                        analysis,
                        verbose,
                        subject_data_transform=blur_subject_data_transform,
                        bilateral_subject_data_transform=(
                            blur_bilateral_subject_data_transform
                        ),
                    )
                    
                    if len(reference_data) == 0:
                        if verbose:
                            print(f"    Warning: No reference data found for {map_key}")
                        continue
                    
                    if len(reference_data) == 1:
                        if verbose:
                            print(f"    Warning: only one subject in reference data for {map_key}, skipping")
                        continue

                    reference_gyrification_covariates = None
                    if use_curvature_covariates:
                        reference_gyrification_covariates, gyr_subjects = load_reference_gyrification_covariates(
                            reference.micapipe_directory,
                            successfully_loaded_subjects,
                            hemi,
                            resolution,
                            analysis=analysis,
                            verbose=verbose,
                        )
                        if len(gyr_subjects) == 0:
                            if verbose:
                                print(f"    Warning: No curvature data found for gyrification adjustment in {map_key}")
                            continue
                        gyr_index = {subject: i for i, subject in enumerate(gyr_subjects)}
                        ref_indices = [
                            i for i, subject in enumerate(successfully_loaded_subjects)
                            if subject in gyr_index
                        ]
                        reference_data = reference_data[ref_indices]
                        successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]
                        reference_gyrification_covariates = np.asarray(
                            [reference_gyrification_covariates[gyr_index[subject]] for subject in successfully_loaded_subjects],
                            dtype=np.float32,
                        )
                    
                    ref_demographics_df = None
                    # Get demographics for reference subjects (wscore, or when
                    # harmonizing/exporting -- the harmonizer needs AGE/SEX).
                    if _need_demographics:
                        ref_demographics = []
                        valid_ref_subjects = []
                        for ref_pid, ref_sid in successfully_loaded_subjects:
                            key = (ref_pid, ref_sid)
                            if key in ref_demo_dict:
                                ref_demographics.append(ref_demo_dict[key])
                                valid_ref_subjects.append(key)

                        if len(ref_demographics) == 0:
                            if verbose:
                                print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                            continue

                        ref_demographics_df = pd.DataFrame(ref_demographics)
                        # Filter reference data to match demographics
                        ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                        reference_data = reference_data[ref_indices]
                        successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]
                        if reference_gyrification_covariates is not None:
                            reference_gyrification_covariates = reference_gyrification_covariates[ref_indices]

                    if control_correlation_filter:
                        if analysis == 'regional':
                            (
                                reference_data,
                                successfully_loaded_subjects,
                                ref_demographics_df,
                                reference_gyrification_covariates,
                            ) = _filter_reference_controls_by_correlation(
                                reference_data,
                                successfully_loaded_subjects,
                                reference_demographics_df=ref_demographics_df,
                                reference_vertex_covariates=reference_gyrification_covariates,
                                quantile=_resolve_control_correlation_quantile(
                                    control_correlation_quantile, feature
                                ),
                                verbose=verbose,
                                context=f"blur {feature} {map_key}",
                            )
                            blur_regional_qc_subjects[hemi] = set(
                                successfully_loaded_subjects
                            )
                        else:
                            allowed_subjects = (
                                blur_regional_qc_subjects.get('L', set())
                                & blur_regional_qc_subjects.get('R', set())
                            )
                            (
                                reference_data,
                                successfully_loaded_subjects,
                                ref_demographics_df,
                                reference_gyrification_covariates,
                            ) = _filter_reference_controls_to_subjects(
                                reference_data,
                                successfully_loaded_subjects,
                                allowed_subjects,
                                reference_demographics_df=ref_demographics_df,
                                reference_vertex_covariates=reference_gyrification_covariates,
                                verbose=verbose,
                                context=f"blur {feature} {map_key}",
                            )
                        if len(reference_data) < 2:
                            if verbose:
                                print(
                                    f"    Warning: fewer than 2 controls remain after "
                                    f"correlation filtering for {map_key}; skipping"
                                )
                            continue
                    
                    # --- Cross-site harmonization hook (analysis-level) ---
                    _hk = f"{feature}|{map_key}"
                    if export_control_features is not None:
                        export_control_features[_hk] = (
                            reference_data, successfully_loaded_subjects, ref_demographics_df)
                        continue
                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                        _bundle = site_harmonizer.pooled_reference(_hk)
                        if _bundle is not None:
                            reference_data, successfully_loaded_subjects, ref_demographics_df = _bundle
                            reference_gyrification_covariates = None

                    # Process patient data
                    patient_scores = []
                    for pat_pid, pat_sid in dataset_blur_subjects:
                        if analysis == 'asymmetry':
                            # For asymmetry, we need both left and right hemisphere data
                            pat_file_lh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"

                            )
                            pat_file_rh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                            )
                            
                            if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                try:
                                    pat_img_lh = nib.load(pat_file_lh)
                                    pat_img_rh = nib.load(pat_file_rh)
                                    
                                    # Handle multi-depth data
                                    pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                    pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                    
                                    for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                        pat_data_lh[:, e] = darray_lh.data
                                        pat_data_rh[:, e] = darray_rh.data

                                    if blur_depth_model == "gradient_flattening":
                                        (
                                            pat_data_lh,
                                            pat_data_rh,
                                        ) = _load_bilateral_blur_gradient_flattening_components(
                                            pat_data_lh,
                                            pat_data_rh,
                                            dataset.micapipe_directory,
                                            pat_pid,
                                            pat_sid,
                                            resolution="32k",
                                        )
                                    elif blur_depth_model in ("gray_white_contrast", "boundary_gradient", "juxtacortical_gradient"):
                                        pat_data_lh = _load_blur_single_gradient(
                                            pat_data_lh, dataset.micapipe_directory,
                                            pat_pid, pat_sid, "L", blur_depth_model,
                                            resolution="32k",
                                        )
                                        pat_data_rh = _load_blur_single_gradient(
                                            pat_data_rh, dataset.micapipe_directory,
                                            pat_pid, pat_sid, "R", blur_depth_model,
                                            resolution="32k",
                                        )

                                    # Compute asymmetry
                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)

                                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                        _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                        if _pd is not None:
                                            pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                    reference_vertex_covariates = None
                                    patient_vertex_covariates = None
                                    if use_curvature_covariates:
                                        patient_vertex_covariates = load_subject_gyrification_covariates(
                                            dataset.micapipe_directory,
                                            pat_pid,
                                            pat_sid,
                                            hemi,
                                            resolution,
                                            analysis=analysis,
                                        )
                                        if patient_vertex_covariates is None:
                                            if verbose:
                                                print(f"    Warning: Missing curvature map for patient {pat_pid}/{pat_sid}; skipping gyrification-adjusted scoring")
                                            continue
                                        reference_vertex_covariates = reference_gyrification_covariates

                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                    score_file = os.path.join(
                                        score_dir,
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )

                                    # Calculate scores (will automatically average across depths)
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose,
                                            reference_vertex_covariates=reference_vertex_covariates,
                                            patient_vertex_covariates=patient_vertex_covariates,
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose,
                                            reference_vertex_covariates=reference_vertex_covariates,
                                            patient_vertex_covariates=patient_vertex_covariates,
                                            use_prediction_uncertainty=predictive_wscore,
                                            wscore_distribution=wscore_distribution,
                                            wscore_preprocessing=wscore_preprocessing,
                                            wscore_covariate_model=wscore_covariate_model,
                                            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                            blur_depth_model=blur_depth_model,
                                            wscore_fit_cache=wscore_fit_cache,
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    score_result[f'avg_{method}_file'] = score_file
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                            else:
                                if verbose:
                                    print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                        
                        else:
                            # Regional analysis - load single hemisphere data
                            pat_file = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix}"
                            )
                            
                            if os.path.exists(pat_file):
                                try:
                                    pat_img = nib.load(pat_file)
                                    pat_data = np.zeros(shape=(pat_img.darrays[0].data.shape[0], len(pat_img.darrays)))
                                    for e, darray in enumerate(pat_img.darrays):
                                        pat_data[:, e] = darray.data

                                    if blur_depth_model == "gradient_flattening":
                                        other_suffix = (
                                            file_suffix.replace("hemi-L", "hemi-R")
                                            if hemi == "L"
                                            else file_suffix.replace("hemi-R", "hemi-L")
                                        )
                                        other_pat_file = os.path.join(
                                            output_directory,
                                            pat_pid,
                                            pat_sid,
                                            "maps",
                                            "cortex",
                                            f"{pat_pid}_{pat_sid}_{other_suffix}",
                                        )
                                        other_img = nib.load(other_pat_file)
                                        other_data = np.column_stack(
                                            [
                                                np.asarray(array.data)
                                                for array in other_img.darrays
                                            ]
                                        )
                                        left_data = pat_data if hemi == "L" else other_data
                                        right_data = other_data if hemi == "L" else pat_data
                                        left_components, right_components = (
                                            _load_bilateral_blur_gradient_flattening_components(
                                                left_data,
                                                right_data,
                                                dataset.micapipe_directory,
                                                pat_pid,
                                                pat_sid,
                                                resolution="32k",
                                            )
                                        )
                                        pat_data = (
                                            left_components
                                            if hemi == "L"
                                            else right_components
                                        )
                                    elif blur_depth_model in ("gray_white_contrast", "boundary_gradient", "juxtacortical_gradient"):
                                        pat_data = _load_blur_single_gradient(
                                            pat_data, dataset.micapipe_directory,
                                            pat_pid, pat_sid, hemi, blur_depth_model,
                                            resolution="32k",
                                        )

                                    if site_harmonizer is not None and site_harmonizer.harmonizes(_hk, site):
                                        _pd = pat_demo_dict.get((pat_pid, pat_sid))
                                        if _pd is not None:
                                            pat_data = site_harmonizer.apply_patient(_hk, pat_data, site, _pd)

                                    reference_vertex_covariates = None
                                    patient_vertex_covariates = None
                                    if use_curvature_covariates:
                                        patient_vertex_covariates = load_subject_gyrification_covariates(
                                            dataset.micapipe_directory,
                                            pat_pid,
                                            pat_sid,
                                            hemi,
                                            resolution,
                                            analysis=analysis,
                                        )
                                        if patient_vertex_covariates is None:
                                            if verbose:
                                                print(f"    Warning: Missing curvature map for patient {pat_pid}/{pat_sid}; skipping gyrification-adjusted scoring")
                                            continue
                                        reference_vertex_covariates = reference_gyrification_covariates
                                    
                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                    score_file = os.path.join(
                                        score_dir,
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores (will automatically average across depths)
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose,
                                            reference_vertex_covariates=reference_vertex_covariates,
                                            patient_vertex_covariates=patient_vertex_covariates,
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose,
                                            reference_vertex_covariates=reference_vertex_covariates,
                                            patient_vertex_covariates=patient_vertex_covariates,
                                            use_prediction_uncertainty=predictive_wscore,
                                            wscore_distribution=wscore_distribution,
                                            wscore_preprocessing=wscore_preprocessing,
                                            wscore_covariate_model=wscore_covariate_model,
                                            wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                            blur_depth_model=blur_depth_model,
                                            wscore_fit_cache=wscore_fit_cache,
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    score_result[f'avg_{method}_file'] = score_file
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient file {pat_file}: {e}")
                    
                    blur_results[map_key] = {
                        f'patient_{method}s': patient_scores
                    }
            
            analysis_results["blur"][feature] = blur_results
        
        # 4. SUBCORTICAL ANALYSIS
        if dataset.subcortical and reference.subcortical and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing subcortical data for {feature}...")
            
            # Get valid subjects for this feature and structure
            dataset_subcortical_subjects = dataset.valid_subjects[feature]['structures']['subcortical'] if feature in dataset.valid_subjects else []
            reference_subcortical_subjects = reference.valid_subjects[feature]['structures']['subcortical'] if feature in reference.valid_subjects else []
            
            if not dataset_subcortical_subjects or not reference_subcortical_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for subcortical {feature}")
                continue
            
            subcortical_results = {}
            
            for analysis in ['regional', 'asymmetry']:
                map_key = f"subcortical_{analysis}"  # Include analysis in map_key
                
                if feat_lower in ["thickness", "volume"]:
                    file_suffix = "feature-volume.csv"
                else:
                    file_suffix = f"feature-{output_feat}.csv"
                
                # Load reference data
                reference_data, successfully_loaded_subjects = load_reference_subcortical_data(
                    reference_subcortical_subjects, output_directory, file_suffix, analysis, verbose
                )
                
                if reference_data.empty:
                    if verbose:
                        print(f"    Warning: No reference data found for {map_key}")
                    continue
                
                # Get demographics for reference subjects (if using wscore)
                if method == 'wscore':
                    ref_demographics = []
                    valid_ref_indices = []
                    for i, (ref_pid, ref_sid) in enumerate(successfully_loaded_subjects):
                        key = (ref_pid, ref_sid)
                        if key in ref_demo_dict:
                            ref_demographics.append(ref_demo_dict[key])
                            valid_ref_indices.append(i)
                    
                    if len(ref_demographics) == 0:
                        if verbose:
                            print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                        continue
                    
                    ref_demographics_df = pd.DataFrame(ref_demographics)
                    # Filter reference data to match demographics
                    reference_data = reference_data.iloc[valid_ref_indices]
                
                # Process patient data
                patient_scores = []
                for pat_pid, pat_sid in dataset_subcortical_subjects:
                    pat_bids_id = f"{pat_pid}_{pat_sid}"
                    pat_file = os.path.join(
                        output_directory,
                        pat_pid, pat_sid, "maps", "subcortical",
                        f"{pat_bids_id}_{file_suffix}"
                    )
                    
                    if os.path.exists(pat_file):
                        try:
                            pat_df = pd.read_csv(pat_file)
                            
                            # For asymmetry analysis, compute left-right differences
                            if analysis == 'asymmetry':
                                # Create asymmetry data by computing differences between left and right structures
                                asymmetry_data = {}
                                for col in pat_df.columns:
                                    if col.startswith('L'):
                                        right_col = col.replace('L', 'R')
                                        if right_col in pat_df.columns:
                                            left_val = pat_df[col].iloc[0]
                                            right_val = pat_df[right_col].iloc[0]
                                            # Compute asymmetry: (L - R) / ((L + R) / 2)
                                            avg_val = (left_val + right_val) / 2
                                            asym_val = (left_val - right_val) / avg_val if avg_val > 0 else 0
                                            asymmetry_data[col] = asym_val
                                
                                # Convert to DataFrame
                                pat_df = pd.DataFrame([asymmetry_data])
                            
                            # Create score output file with analysis type in filename
                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "subcortical")
                            score_file = os.path.join(
                                score_dir, 
                                f"{pat_bids_id}_{file_suffix.replace('.csv', f'_analysis-{analysis}.csv')}"

                            )
                            
                            # Calculate scores
                            if method == 'zscore':
                                score_result = calculate_zscore_csv(
                                    reference_data, pat_df, score_file, verbose
                                )
                            else:  # wscore
                                # Get patient demographics
                                pat_key = (pat_pid, pat_sid)
                                if pat_key not in pat_demo_dict:
                                    if verbose:
                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                    continue
                                
                                pat_demographics = pat_demo_dict[pat_key]
                                score_result = calculate_wscore_csv(
                                    reference_data, pat_df, ref_demographics_df, pat_demographics, 
                                    score_file, normative_columns, verbose,
                                    use_prediction_uncertainty=predictive_wscore,
                                    wscore_distribution=wscore_distribution,
                                    wscore_preprocessing=wscore_preprocessing,
                                    wscore_covariate_model=wscore_covariate_model,
                                    wscore_fit_cache=wscore_fit_cache,
                                )
                            
                            score_result['subject'] = (pat_pid, pat_sid)
                            score_result['analysis'] = analysis
                            patient_scores.append(score_result)
                            
                        except Exception as e:
                            if verbose:
                                print(f"    Warning: Could not process patient file {pat_file}: {e}")
                
                subcortical_results[map_key] = {
                    f'patient_{method}s': patient_scores
                }
            
            analysis_results["subcortical"][feature] = subcortical_results
    
        # 5. fMRI ANALYSIS
        if is_fmri:
            if verbose:
                print(f"  Processing fMRI data for {feature}...")
            
            # Get valid subjects for this feature
            dataset_fmri_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_fmri_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_fmri_subjects or not reference_fmri_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for fMRI {feature}")
                continue
            
            fmri_features = ['rmssd', 'timescales', 'alff', 'falff']
            
            for fmri_feat in fmri_features:
                if verbose:
                    print(f"    Processing fMRI feature: {fmri_feat}")
                
                fmri_results = {}
                fmri_regional_qc_subjects = {}
                
                for analysis in ['regional', 'asymmetry']:
                    for hemi in hemispheres:
                        for resolution in resolutions:
                            for label in ['midthickness']:
                                map_key = f"{hemi}_{resolution}_{label}_{analysis}"
                                
                                file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{fmri_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"
                                
                                # Load reference data and track which subjects actually loaded
                                reference_data, successfully_loaded_subjects = load_reference_surface_data(
                                    reference_fmri_subjects, output_directory, file_suffix, analysis, verbose
                                )
                                
                                if len(reference_data) == 0:
                                    if verbose:
                                        print(f"    Warning: No reference data found for {map_key}")
                                    continue
                                
                                ref_demographics_df = None
                                # Get demographics for reference subjects (if using wscore)
                                if method == 'wscore':
                                    ref_demographics = []
                                    valid_ref_subjects = []
                                    # Now iterate over subjects that ACTUALLY loaded, not all subjects
                                    for ref_pid, ref_sid in successfully_loaded_subjects:
                                        key = (ref_pid, ref_sid)
                                        if key in ref_demo_dict:
                                            ref_demographics.append(ref_demo_dict[key])
                                            valid_ref_subjects.append(key)
                                    
                                    if len(ref_demographics) == 0:
                                        if verbose:
                                            print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                                        continue
                                    
                                    ref_demographics_df = pd.DataFrame(ref_demographics)
                                    # Filter reference data to match demographics
                                    # Now indices will match because we're iterating over successfully_loaded_subjects
                                    ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                                    reference_data = reference_data[ref_indices]
                                    successfully_loaded_subjects = [successfully_loaded_subjects[i] for i in ref_indices]

                                if control_correlation_filter:
                                    if analysis == 'regional':
                                        (
                                            reference_data,
                                            successfully_loaded_subjects,
                                            ref_demographics_df,
                                            _,
                                        ) = _filter_reference_controls_by_correlation(
                                            reference_data,
                                            successfully_loaded_subjects,
                                            reference_demographics_df=ref_demographics_df,
                                            quantile=_resolve_control_correlation_quantile(
                                                control_correlation_quantile, fmri_feat
                                            ),
                                            verbose=verbose,
                                            context=f"fMRI {fmri_feat} {map_key}",
                                        )
                                        fmri_regional_qc_subjects[(resolution, label, hemi)] = set(
                                            successfully_loaded_subjects
                                        )
                                    else:
                                        allowed_subjects = (
                                            fmri_regional_qc_subjects.get((resolution, label, 'L'), set())
                                            & fmri_regional_qc_subjects.get((resolution, label, 'R'), set())
                                        )
                                        (
                                            reference_data,
                                            successfully_loaded_subjects,
                                            ref_demographics_df,
                                            _,
                                        ) = _filter_reference_controls_to_subjects(
                                            reference_data,
                                            successfully_loaded_subjects,
                                            allowed_subjects,
                                            reference_demographics_df=ref_demographics_df,
                                            verbose=verbose,
                                            context=f"fMRI {fmri_feat} {map_key}",
                                        )
                                    if len(reference_data) < 2:
                                        if verbose:
                                            print(
                                                f"    Warning: fewer than 2 controls remain after "
                                                f"correlation filtering for {map_key}; skipping"
                                            )
                                        continue
                                
                                # Process patient data
                                patient_scores = []
                                for pat_pid, pat_sid in dataset_fmri_subjects:
                                    pat_bids_id = f"{pat_pid}_{pat_sid}"
                                    
                                    if analysis == 'asymmetry':
                                        # For asymmetry, we need both left and right hemisphere data
                                        pat_file_lh = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
                                        )
                                        pat_file_rh = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                                        )
                                        
                                        if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                            try:
                                                pat_img_lh = nib.load(pat_file_lh)
                                                pat_img_rh = nib.load(pat_file_rh)
                                                
                                                # Handle multi-depth data (blur features)
                                                if len(pat_img_lh.darrays) > 1:
                                                    pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                                    pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                                    
                                                    for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                                        pat_data_lh[:, e] = darray_lh.data
                                                        pat_data_rh[:, e] = darray_rh.data
                                                    
                                                    # Compute asymmetry
                                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                                else:
                                                    pat_data_lh = pat_img_lh.darrays[0].data
                                                    pat_data_rh = pat_img_rh.darrays[0].data
                                                    
                                                    # Compute asymmetry
                                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                                
                                                # Create score output file with analysis type in filename
                                                score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                                score_file = os.path.join(
                                                    score_dir, 
                                                    f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                                )
                                                
                                                # Calculate scores
                                                if method == 'zscore':
                                                    score_result = calculate_zscore_maps(
                                                        reference_data, pat_data, score_file, analysis, verbose
                                                    )
                                                else:  # wscore
                                                    # Get patient demographics
                                                    pat_key = (pat_pid, pat_sid)
                                                    if pat_key not in pat_demo_dict:
                                                        if verbose:
                                                            print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                        continue
                                                    
                                                    pat_demographics = pat_demo_dict[pat_key]
                                                    score_result = calculate_wscore_maps(
                                                        reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                        score_file, normative_columns, verbose,
                                                        use_prediction_uncertainty=predictive_wscore,
                                                        wscore_distribution=wscore_distribution,
                                                        wscore_preprocessing=wscore_preprocessing,
                                                        wscore_covariate_model=wscore_covariate_model,
                                                        wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                                        blur_depth_model=blur_depth_model,
                                                        wscore_fit_cache=wscore_fit_cache,
                                                    )
                                                
                                                score_result['subject'] = (pat_pid, pat_sid)
                                                score_result['analysis'] = analysis
                                                patient_scores.append(score_result)
                                                
                                            except Exception as e:
                                                if verbose:
                                                    print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                                        else:
                                            if verbose:
                                                print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                                    
                                    else:
                                        # Regional analysis - load single hemisphere data
                                        pat_file = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix}"
                                        )
                                        
                                        if os.path.exists(pat_file):
                                            try:
                                                pat_img = nib.load(pat_file)
                                                pat_data = pat_img.darrays[0].data
                                                
                                                # Create score output file with analysis type in filename
                                                score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                                score_file = os.path.join(
                                                    score_dir, 
                                                    f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                                )
                                                
                                                # Calculate scores
                                                if method == 'zscore':
                                                    score_result = calculate_zscore_maps(
                                                        reference_data, pat_data, score_file, analysis, verbose
                                                    )
                                                else:  # wscore
                                                    # Get patient demographics
                                                    pat_key = (pat_pid, pat_sid)
                                                    if pat_key not in pat_demo_dict:
                                                        if verbose:
                                                            print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                        continue
                                                    
                                                    pat_demographics = pat_demo_dict[pat_key]
                                                    score_result = calculate_wscore_maps(
                                                        reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                        score_file, normative_columns, verbose,
                                                        use_prediction_uncertainty=predictive_wscore,
                                                        wscore_distribution=wscore_distribution,
                                                        wscore_preprocessing=wscore_preprocessing,
                                                        wscore_covariate_model=wscore_covariate_model,
                                                        wscore_surface_smoothing_iterations=wscore_surface_smoothing_iterations,
                                                        blur_depth_model=blur_depth_model,
                                                        wscore_fit_cache=wscore_fit_cache,
                                                    )
                                                
                                                score_result['subject'] = (pat_pid, pat_sid)
                                                score_result['analysis'] = analysis
                                                patient_scores.append(score_result)
                                                
                                            except Exception as e:
                                                if verbose:
                                                    print(f"    Warning: Could not process patient file {pat_file}: {e}")
                                
                                fmri_results[map_key] = {
                                    f'patient_{method}s': patient_scores
                                }
                
                analysis_results["cortical"][fmri_feat] = fmri_results

    if verbose:
        print(f"\nAnalysis complete! {method.upper()} maps saved to {method}_maps directories")
        
        # Print summary statistics
        for region_type, region_results in analysis_results.items():
            if region_results:
                print(f"\n{region_type.capitalize()} analysis:")
                for feature, feature_results in region_results.items():
                    if region_type == "subcortical":
                        if f'patient_{method}s' in feature_results:
                            print(f"  {feature}: {len(feature_results[f'patient_{method}s'])} subjects analyzed")
                    else:
                        map_count = sum(len(maps[f'patient_{method}s']) for maps in feature_results.values())
                        print(f"  {feature}: {map_count} {method} maps generated")
    
    return analysis_results

def compute_asymmetry(data_lh, data_rh):
    """Compute asymmetry.

    Parameters
    ----------
    x_lh
        Left hemisphere data. Shape (n_subjects, n_points) or (n_points,).
    x_rh
        Right hemisphere data. Shape (n_subjects, n_points) or (n_points,).

    Returns
    -------
    s
        Output data. Shape (n_subjects, n_points) or (n_points,).
    """

    den = data_lh + data_rh
    den *= 0.5
    return np.divide(data_lh - data_rh, den, out=den, where=den > 0)
