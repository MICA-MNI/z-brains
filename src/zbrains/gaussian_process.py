"""Vectorized demographic Gaussian-process normative modeling."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import linalg, optimize


EPS = np.finfo(float).eps

GP_KERNEL_LABELS = {
    "rbf_ard": "ARD_RBF_plus_white_noise",
    "rbf_isotropic": "isotropic_RBF_plus_white_noise",
    "matern32_ard": "ARD_Matern32_plus_white_noise",
    "matern52_ard": "ARD_Matern52_plus_white_noise",
    "rational_quadratic_ard": "ARD_rational_quadratic_plus_white_noise",
}


def normalize_gaussian_process_kernel(kernel_name: str) -> str:
    """Validate a Gaussian-process covariance family name."""

    normalized = str(kernel_name).strip().lower().replace("-", "_")
    aliases = {
        "rbf": "rbf_ard",
        "ard_rbf": "rbf_ard",
        "isotropic_rbf": "rbf_isotropic",
        "matern32": "matern32_ard",
        "matern_3_2": "matern32_ard",
        "matern52": "matern52_ard",
        "matern_5_2": "matern52_ard",
        "rational_quadratic": "rational_quadratic_ard",
        "rq": "rational_quadratic_ard",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in GP_KERNEL_LABELS:
        choices = ", ".join(sorted(GP_KERNEL_LABELS))
        raise ValueError(
            f"Unsupported Gaussian-process kernel '{kernel_name}'. "
            f"Choose one of: {choices}."
        )
    return normalized


@dataclass
class GaussianProcessGroup:
    """Vertices sharing the same set of finite reference subjects."""

    vertices: np.ndarray
    subjects: np.ndarray
    center: np.ndarray
    scale: np.ndarray
    alpha: np.ndarray
    cholesky: np.ndarray


@dataclass
class GaussianProcessFit:
    """Reusable fit for independent output GPs with one shared ARD kernel."""

    demographic_center: np.ndarray
    demographic_scale: np.ndarray
    standardized_demographics: np.ndarray
    kernel_name: str
    kernel_shape_parameter: float | None
    signal_variance: float
    noise_variance: float
    length_scales: np.ndarray
    groups: list[GaussianProcessGroup]
    constant_vertices: np.ndarray
    constant_values: np.ndarray
    reference_counts: np.ndarray
    optimization_success: bool
    optimization_message: str
    optimization_vertices: int
    loo_score_sd: float
    loo_abs_gt_1p96: float

    def metadata(self, normative_columns: list[str]) -> dict[str, object]:
        return {
            "kernel": GP_KERNEL_LABELS[self.kernel_name],
            "kernel_name": self.kernel_name,
            "kernel_shape_parameter": self.kernel_shape_parameter,
            "normative_columns": list(normative_columns),
            "demographic_center": self.demographic_center.tolist(),
            "demographic_scale": self.demographic_scale.tolist(),
            "signal_variance": float(self.signal_variance),
            "noise_variance": float(self.noise_variance),
            "length_scales": self.length_scales.tolist(),
            "optimization_success": bool(self.optimization_success),
            "optimization_message": self.optimization_message,
            "optimization_vertices": int(self.optimization_vertices),
            "loo_score_sd": float(self.loo_score_sd),
            "loo_abs_gt_1p96": float(self.loo_abs_gt_1p96),
            "predictive_variance": "latent_posterior_variance_plus_white_noise",
        }


def _standardize_demographics(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    center = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
    return (values - center) / scale, center, scale


def _kernel(
    left: np.ndarray,
    right: np.ndarray,
    signal_variance: float,
    length_scales: np.ndarray,
    kernel_name: str,
    shape_parameter: float | None,
) -> np.ndarray:
    difference = (
        left[:, None, :] - right[None, :, :]
    ) / length_scales[None, None, :]
    squared_distance = np.sum(difference * difference, axis=2)
    if kernel_name in {"rbf_ard", "rbf_isotropic"}:
        correlation = np.exp(-0.5 * squared_distance)
    elif kernel_name == "matern32_ard":
        scaled_distance = np.sqrt(3.0 * np.maximum(squared_distance, 0.0))
        correlation = (1.0 + scaled_distance) * np.exp(-scaled_distance)
    elif kernel_name == "matern52_ard":
        scaled_distance = np.sqrt(5.0 * np.maximum(squared_distance, 0.0))
        correlation = (
            1.0 + scaled_distance + 5.0 * squared_distance / 3.0
        ) * np.exp(-scaled_distance)
    elif kernel_name == "rational_quadratic_ard":
        alpha = max(float(shape_parameter or 1.0), 1e-8)
        correlation = (1.0 + squared_distance / (2.0 * alpha)) ** (-alpha)
    else:
        raise ValueError(f"Unsupported Gaussian-process kernel: {kernel_name}")
    return signal_variance * correlation


def _factor_kernel(
    demographics: np.ndarray,
    signal_variance: float,
    noise_variance: float,
    length_scales: np.ndarray,
    kernel_name: str,
    shape_parameter: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    covariance = _kernel(
        demographics,
        demographics,
        signal_variance,
        length_scales,
        kernel_name,
        shape_parameter,
    )
    covariance.flat[:: covariance.shape[0] + 1] += noise_variance + 1e-6
    cholesky = linalg.cholesky(covariance, lower=True, check_finite=False)
    return covariance, cholesky


def _unpack_kernel_parameters(
    log_parameters: np.ndarray,
    n_features: int,
    kernel_name: str,
) -> tuple[float, float, np.ndarray, float | None]:
    signal, noise = np.exp(log_parameters[:2])
    length_count = 1 if kernel_name == "rbf_isotropic" else n_features
    raw_length_scales = np.exp(log_parameters[2 : 2 + length_count])
    length_scales = (
        np.repeat(raw_length_scales, n_features)
        if length_count == 1
        else raw_length_scales
    )
    shape_parameter = (
        float(np.exp(log_parameters[2 + length_count]))
        if kernel_name == "rational_quadratic_ard"
        else None
    )
    return float(signal), float(noise), length_scales, shape_parameter


def _select_optimization_outputs(reference_data: np.ndarray, limit: int) -> np.ndarray:
    finite = np.isfinite(reference_data).all(axis=0)
    scale = np.nanstd(reference_data, axis=0)
    candidates = np.flatnonzero(finite & np.isfinite(scale) & (scale > EPS))
    if candidates.size <= limit:
        return candidates
    positions = np.linspace(0, candidates.size - 1, limit, dtype=int)
    return candidates[positions]


def _optimize_kernel(
    demographics: np.ndarray,
    reference_data: np.ndarray,
    max_outputs: int,
    max_iterations: int,
    kernel_name: str,
) -> tuple[np.ndarray, bool, str, int, float, float]:
    output_indices = _select_optimization_outputs(reference_data, max_outputs)
    n_features = demographics.shape[1]
    length_count = 1 if kernel_name == "rbf_isotropic" else n_features
    initial_values = [0.8, 0.2, *np.ones(length_count)]
    bounds = [
        (np.log(0.05), np.log(5.0)),
        (np.log(1e-3), np.log(2.0)),
        *[(np.log(0.1), np.log(10.0)) for _ in range(length_count)],
    ]
    if kernel_name == "rational_quadratic_ard":
        initial_values.append(1.0)
        bounds.append((np.log(0.05), np.log(20.0)))
    initial = np.log(np.asarray(initial_values, dtype=float))
    if output_indices.size == 0:
        return initial, False, "No complete variable outputs; used defaults", 0, np.nan, np.nan

    values = reference_data[:, output_indices]
    center = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    standardized = (values - center) / scale

    def objective(log_parameters: np.ndarray) -> float:
        signal, noise, length_scales, shape_parameter = _unpack_kernel_parameters(
            log_parameters,
            n_features,
            kernel_name,
        )
        try:
            _, cholesky = _factor_kernel(
                demographics,
                signal,
                noise,
                length_scales,
                kernel_name,
                shape_parameter,
            )
            alpha = linalg.cho_solve(
                (cholesky, True),
                standardized,
                check_finite=False,
            )
        except linalg.LinAlgError:
            return 1e30
        quadratic = 0.5 * np.mean(np.sum(standardized * alpha, axis=0))
        log_determinant = np.sum(np.log(np.diag(cholesky)))
        return float(
            quadratic
            + log_determinant
            + 0.5 * demographics.shape[0] * np.log(2.0 * np.pi)
        )

    result = optimize.minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": int(max_iterations), "ftol": 1e-7},
    )
    parameters = result.x if np.isfinite(result.fun) else initial

    signal, noise, length_scales, shape_parameter = _unpack_kernel_parameters(
        parameters,
        n_features,
        kernel_name,
    )
    try:
        covariance, cholesky = _factor_kernel(
            demographics,
            signal,
            noise,
            length_scales,
            kernel_name,
            shape_parameter,
        )
        alpha = linalg.cho_solve(
            (cholesky, True),
            standardized,
            check_finite=False,
        )
        inverse = linalg.cho_solve(
            (cholesky, True),
            np.eye(covariance.shape[0]),
            check_finite=False,
        )
        inverse_diagonal = np.maximum(np.diag(inverse), EPS)
        loo_scores = alpha / np.sqrt(inverse_diagonal[:, None])
        loo_score_sd = float(np.std(loo_scores))
        loo_abs_gt_1p96 = float(np.mean(np.abs(loo_scores) > 1.96))
    except linalg.LinAlgError:
        loo_score_sd = np.nan
        loo_abs_gt_1p96 = np.nan
    return (
        parameters,
        bool(result.success),
        str(result.message),
        int(output_indices.size),
        loo_score_sd,
        loo_abs_gt_1p96,
    )


def fit_gaussian_process(
    reference_data: np.ndarray,
    demographics: np.ndarray,
    min_reference_subjects: int | None = None,
    max_optimization_outputs: int = 128,
    max_optimization_iterations: int = 35,
    kernel_name: str = "rbf_ard",
) -> GaussianProcessFit:
    """Fit independent output GPs using shared demographic kernel parameters."""

    reference_data = np.asarray(reference_data, dtype=float)
    demographics = np.asarray(demographics, dtype=float)
    if reference_data.ndim != 2:
        raise ValueError("reference_data must have shape (subjects, outputs)")
    if demographics.ndim != 2 or demographics.shape[0] != reference_data.shape[0]:
        raise ValueError("demographics must have shape (subjects, covariates)")
    if not np.isfinite(demographics).all():
        raise ValueError("demographics must be finite before GP fitting")
    kernel_name = normalize_gaussian_process_kernel(kernel_name)

    standardized_demographics, demo_center, demo_scale = _standardize_demographics(
        demographics
    )
    parameters, success, message, optimization_vertices, loo_sd, loo_rate = (
        _optimize_kernel(
            standardized_demographics,
            reference_data,
            max_optimization_outputs,
            max_optimization_iterations,
            kernel_name,
        )
    )
    signal, noise, length_scales, shape_parameter = _unpack_kernel_parameters(
        parameters,
        standardized_demographics.shape[1],
        kernel_name,
    )

    minimum = (
        int(min_reference_subjects)
        if min_reference_subjects is not None
        else max(demographics.shape[1] + 3, 8)
    )
    finite = np.isfinite(reference_data)
    reference_counts = finite.sum(axis=0).astype(int)
    patterns: dict[bytes, tuple[np.ndarray, list[int]]] = {}
    constant_vertices = []
    constant_values = []
    for vertex in range(reference_data.shape[1]):
        subject_mask = finite[:, vertex]
        if int(subject_mask.sum()) < minimum:
            continue
        values = reference_data[subject_mask, vertex]
        scale = float(np.std(values))
        if not np.isfinite(scale) or scale <= EPS:
            constant_vertices.append(vertex)
            constant_values.append(float(np.mean(values)))
            continue
        key = np.packbits(subject_mask).tobytes()
        if key not in patterns:
            patterns[key] = (subject_mask, [])
        patterns[key][1].append(vertex)

    groups = []
    for subject_mask, vertices_list in patterns.values():
        vertices = np.asarray(vertices_list, dtype=int)
        subjects = np.flatnonzero(subject_mask)
        values = reference_data[np.ix_(subjects, vertices)]
        center = np.mean(values, axis=0)
        scale = np.std(values, axis=0)
        standardized = (values - center) / scale
        _, cholesky = _factor_kernel(
            standardized_demographics[subjects],
            signal,
            noise,
            length_scales,
            kernel_name,
            shape_parameter,
        )
        alpha = linalg.cho_solve(
            (cholesky, True),
            standardized,
            check_finite=False,
        )
        groups.append(
            GaussianProcessGroup(
                vertices=vertices,
                subjects=subjects,
                center=center.astype(np.float32),
                scale=scale.astype(np.float32),
                alpha=alpha.astype(np.float32),
                cholesky=cholesky,
            )
        )

    return GaussianProcessFit(
        demographic_center=demo_center,
        demographic_scale=demo_scale,
        standardized_demographics=standardized_demographics,
        kernel_name=kernel_name,
        kernel_shape_parameter=shape_parameter,
        signal_variance=float(signal),
        noise_variance=float(noise),
        length_scales=length_scales,
        groups=groups,
        constant_vertices=np.asarray(constant_vertices, dtype=int),
        constant_values=np.asarray(constant_values, dtype=float),
        reference_counts=reference_counts,
        optimization_success=success,
        optimization_message=message,
        optimization_vertices=optimization_vertices,
        loo_score_sd=loo_sd,
        loo_abs_gt_1p96=loo_rate,
    )


def predict_gaussian_process(
    fit: GaussianProcessFit,
    patient_data: np.ndarray,
    patient_demographics: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return predictive deviation scores and predictive variances."""

    patient_data = np.asarray(patient_data, dtype=float).reshape(-1)
    patient_demographics = np.asarray(patient_demographics, dtype=float).reshape(-1)
    if patient_demographics.size != fit.demographic_center.size:
        raise ValueError("Patient demographics do not match the fitted GP covariates")
    standardized_patient = (
        patient_demographics - fit.demographic_center
    ) / fit.demographic_scale
    scores = np.full(patient_data.size, np.nan, dtype=float)
    predictive_variances = np.full(patient_data.size, np.nan, dtype=float)

    for group in fit.groups:
        training_demographics = fit.standardized_demographics[group.subjects]
        cross_covariance = _kernel(
            standardized_patient[None, :],
            training_demographics,
            fit.signal_variance,
            fit.length_scales,
            fit.kernel_name,
            fit.kernel_shape_parameter,
        ).reshape(-1)
        predicted_standardized = cross_covariance @ group.alpha
        solved = linalg.solve_triangular(
            group.cholesky,
            cross_covariance,
            lower=True,
            check_finite=False,
        )
        variance = max(
            fit.signal_variance
            + fit.noise_variance
            - float(solved @ solved),
            1e-8,
        )
        observed = patient_data[group.vertices]
        observed_standardized = (observed - group.center) / group.scale
        valid = np.isfinite(observed_standardized)
        group_scores = np.full(group.vertices.size, np.nan, dtype=float)
        group_scores[valid] = (
            observed_standardized[valid] - predicted_standardized[valid]
        ) / np.sqrt(variance)
        scores[group.vertices] = group_scores
        predictive_variances[group.vertices] = variance

    if fit.constant_vertices.size:
        observed = patient_data[fit.constant_vertices]
        matches = np.isfinite(observed) & np.isclose(
            observed,
            fit.constant_values,
            rtol=1e-6,
            atol=1e-8,
        )
        scores[fit.constant_vertices[matches]] = 0.0
        predictive_variances[fit.constant_vertices[matches]] = 0.0
    return scores, predictive_variances


def leave_one_out_gaussian_process(fit: GaussianProcessFit) -> np.ndarray:
    """Return analytic leave-one-out predictive scores for all fitted outputs."""

    n_subjects = fit.standardized_demographics.shape[0]
    n_outputs = fit.reference_counts.size
    scores = np.full((n_subjects, n_outputs), np.nan, dtype=np.float32)
    for group in fit.groups:
        inverse = linalg.cho_solve(
            (group.cholesky, True),
            np.eye(group.subjects.size),
            check_finite=False,
        )
        inverse_diagonal = np.maximum(np.diag(inverse), EPS)
        group_scores = group.alpha / np.sqrt(inverse_diagonal[:, None])
        scores[np.ix_(group.subjects, group.vertices)] = group_scores.astype(
            np.float32
        )
    if fit.constant_vertices.size:
        scores[:, fit.constant_vertices] = 0.0
    return scores

def _prepare_gp_feature_arrays(
    reference_data: np.ndarray,
    patient_data: np.ndarray,
    wscore_preprocessing: str,
    blur_depth_model: str,
    intensity_depth_model: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Apply the existing feature reductions before GP normative fitting."""

    reference_data = np.asarray(reference_data, dtype=float).copy()
    patient_data = np.asarray(patient_data, dtype=float).copy()
    metadata: dict[str, object] = {
        "input_was_multidepth": reference_data.ndim > 2,
        "blur_component_vertices": None,
        "intensity_component_vertices": None,
        "intensity_component_count": 0,
    }

    if reference_data.ndim > 2:
        if patient_data.ndim <= 1:
            raise ValueError(
                "Patient multi-depth data must contain the control depth dimension"
            )
        if reference_data.shape[2] != patient_data.shape[1]:
            raise ValueError(
                "Control and patient data must contain the same number of depths"
            )
        if intensity_depth_model == "multisurface_median_abs_dominant":
            if reference_data.shape[2] != 4:
                raise ValueError(
                    "multisurface_median_abs_dominant requires four depths"
                )
            metadata["intensity_component_vertices"] = reference_data.shape[1]
            metadata["intensity_component_count"] = reference_data.shape[2]
            reference_data = np.transpose(reference_data, (0, 2, 1)).reshape(
                reference_data.shape[0], -1
            )
            patient_data = patient_data.T.reshape(-1)
        elif blur_depth_model == "mean":
            reference_data = np.mean(reference_data, axis=2)
            patient_data = np.mean(patient_data, axis=1)
        elif blur_depth_model == "gradient_flattening":
            if reference_data.shape[2] != 2:
                raise ValueError(
                    "gradient_flattening expects gray-white and white-SWM1 components"
                )
            metadata["blur_component_vertices"] = reference_data.shape[1]
            reference_data = np.transpose(reference_data, (0, 2, 1)).reshape(
                reference_data.shape[0], -1
            )
            patient_data = patient_data.T.reshape(-1)
        else:
            depth = np.arange(reference_data.shape[2], dtype=float)
            depth -= np.mean(depth)
            weights = depth / np.sum(depth**2)
            reference_mean = np.mean(reference_data, axis=2)
            reference_slope = np.tensordot(reference_data, weights, axes=([2], [0]))
            patient_mean = np.mean(patient_data, axis=1)
            patient_slope = np.tensordot(patient_data, weights, axes=([1], [0]))
            metadata["blur_component_vertices"] = reference_mean.shape[1]
            reference_data = np.concatenate(
                [reference_mean, reference_slope], axis=1
            )
            patient_data = np.concatenate([patient_mean, patient_slope])
    elif patient_data.ndim > 1:
        patient_data = np.mean(patient_data, axis=1)

    def spatial_normalize(values: np.ndarray, robust: bool) -> np.ndarray:
        one_dimensional = values.ndim == 1
        matrix = values[None, :] if one_dimensional else values
        if robust:
            center = np.nanmedian(matrix, axis=1, keepdims=True)
            scale = 1.482602218505602 * np.nanmedian(
                np.abs(matrix - center), axis=1, keepdims=True
            )
        else:
            center = np.nanmean(matrix, axis=1, keepdims=True)
            scale = np.nanstd(matrix, axis=1, keepdims=True)
        scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
        normalized = (matrix - center) / scale
        return normalized[0] if one_dimensional else normalized

    if wscore_preprocessing != "none":
        robust = wscore_preprocessing == "spatial_robust_z"
        component_vertices = (
            metadata["intensity_component_vertices"]
            or metadata["blur_component_vertices"]
        )
        if component_vertices is None:
            reference_data = spatial_normalize(reference_data, robust)
            patient_data = spatial_normalize(patient_data, robust)
        else:
            component_count = (
                metadata["intensity_component_count"]
                if metadata["intensity_component_vertices"] is not None
                else 2
            )
            for component in range(int(component_count)):
                component_slice = slice(
                    component * int(component_vertices),
                    (component + 1) * int(component_vertices),
                )
                reference_data[:, component_slice] = spatial_normalize(
                    reference_data[:, component_slice], robust
                )
                patient_data[component_slice] = spatial_normalize(
                    patient_data[component_slice], robust
                )
    return reference_data, patient_data, metadata


def calculate_gaussian_process_maps(
    reference_data,
    patient_data,
    demographics_ref,
    demographics_pat,
    output_file,
    normative_columns,
    verbose,
    gaussian_process_kernel,
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
    """Fit a demographic GP and save predictive deviation scores as GIFTI.

    When ``prediction_variance_percentile`` is not None, output vertices whose
    GP predictive variance exceeds that within-map percentile are masked to NaN
    (excluded from downstream detection). This is the stage-4 "GP + demographics
    + prediction-uncertainty filtering" rung of the greedy optimizer: it suppresses
    scores where the normative model is least confident.
    """

    if reference_vertex_covariates is not None or patient_vertex_covariates is not None:
        raise ValueError(
            "gaussian_process currently models demographic covariates only; "
            "disable use_curvature_covariates"
        )
    cache_reference_data = reference_data
    reference_data, patient_data, feature_metadata = _prepare_gp_feature_arrays(
        reference_data,
        patient_data,
        wscore_preprocessing,
        blur_depth_model,
        intensity_depth_model,
    )

    demo_ref = demographics_ref[normative_columns].to_numpy(dtype=float)
    demo_pat = demographics_pat[normative_columns]
    demo_pat = np.asarray(demo_pat, dtype=float).reshape(-1)
    demographic_mask = np.isfinite(demo_ref).all(axis=1)
    demo_ref = demo_ref[demographic_mask]
    reference_data = reference_data[demographic_mask]
    if not np.isfinite(demo_pat).all():
        raise ValueError("Patient demographics must be complete for GP scoring")
    if demo_ref.shape[0] < max(len(normative_columns) + 3, 8):
        raise ValueError("Too few controls with complete demographics for GP scoring")

    signature = (
        id(cache_reference_data),
        id(demographics_ref),
        tuple(normative_columns),
        wscore_preprocessing,
        blur_depth_model,
        intensity_depth_model,
        gaussian_process_kernel,
    )
    cache_key = "gaussian_process_surface_current"
    cached = wscore_fit_cache.get(cache_key) if wscore_fit_cache is not None else None
    if cached is not None and cached.get("signature") == signature:
        fit = cached["fit"]
    else:
        fit = fit_gaussian_process(
            reference_data,
            demo_ref,
            min_reference_subjects=min_reference_subjects,
            kernel_name=gaussian_process_kernel,
        )
        cached = {"signature": signature, "fit": fit}
        if wscore_fit_cache is not None:
            wscore_fit_cache.clear()
            wscore_fit_cache[cache_key] = cached

    scores, predictive_variances = predict_gaussian_process(
        fit,
        patient_data,
        demo_pat,
    )
    intensity_vertices = feature_metadata["intensity_component_vertices"]
    blur_vertices = feature_metadata["blur_component_vertices"]
    dominant_depth_indices = None
    dominant_depth_file = None
    dominant_null_threshold = None
    dominant_calibration_factor = 1.0
    # Per-vertex prediction uncertainty for the optional GP filter = the
    # LEAVE-ONE-OUT residual dispersion (genuine model reliability). At each vertex
    # it is the mean squared LOO standardized score over controls: ~1 where the GP
    # predicts its own held-out controls well, and >1 where it predicts them poorly
    # (least trustworthy). This replaces the previous control-SD^2 proxy, which was
    # NOT model confidence -- it just flagged the highest-amplitude cortex, i.e.
    # exactly the gray-white-blurring regions where FCD signal lives, so masking by
    # it discarded the target signal. LOO dispersion is amplitude-independent
    # (scores are already standardized) and high precisely where the normative fit
    # is unreliable and false-positive prone.
    loo_scores_all = leave_one_out_gaussian_process(fit)   # (n_controls, n_outputs)
    with np.errstate(invalid="ignore"):
        loo_uncertainty = np.nanmean(np.square(np.asarray(loo_scores_all, dtype=float)), axis=0)
    loo_uncertainty = np.where(np.isfinite(loo_uncertainty), loo_uncertainty, np.inf)
    if fit.constant_vertices.size:
        loo_uncertainty[fit.constant_vertices] = np.inf   # unusable vertices -> maskable
    pre_reduction_variance = loo_uncertainty
    # Reduced to the same output-vertex layout as ``scores`` (default: already
    # output layout for plain regional maps).
    output_variance = pre_reduction_variance

    if intensity_vertices is not None:
        intensity_vertices = int(intensity_vertices)
        component_count = int(feature_metadata["intensity_component_count"])
        if "dominant_calibration_factor" not in cached:
            loo_scores = leave_one_out_gaussian_process(fit).reshape(
                demo_ref.shape[0],
                component_count,
                intensity_vertices,
            )
            control_dominant = np.nanmax(np.abs(loo_scores), axis=1)
            finite_null = control_dominant[np.isfinite(control_dominant)]
            dominant_null_threshold = (
                float(np.percentile(finite_null, 95.0))
                if finite_null.size
                else 1.96
            )
            dominant_calibration_factor = max(
                dominant_null_threshold / 1.96,
                1e-8,
            )
            cached["dominant_null_threshold"] = dominant_null_threshold
            cached["dominant_calibration_factor"] = dominant_calibration_factor
        else:
            dominant_null_threshold = cached["dominant_null_threshold"]
            dominant_calibration_factor = cached["dominant_calibration_factor"]
        component_scores = scores.reshape(component_count, intensity_vertices)
        dominant_depth_indices = np.argmax(np.abs(component_scores), axis=0)
        scores = np.take_along_axis(
            component_scores,
            dominant_depth_indices[None, :],
            axis=0,
        )[0] / dominant_calibration_factor
        output_variance = np.take_along_axis(
            pre_reduction_variance.reshape(component_count, intensity_vertices),
            dominant_depth_indices[None, :],
            axis=0,
        )[0]
    elif blur_vertices is not None:
        blur_vertices = int(blur_vertices)
        components = scores.reshape(2, blur_vertices)
        if blur_depth_model == "gradient_flattening":
            scores = -np.mean(components, axis=0)
        else:
            magnitude = np.sqrt(np.mean(components**2, axis=0))
            scores = np.copysign(magnitude, components[0])
        # A blur vertex is uncertain if either depth-gradient component is.
        output_variance = np.nanmax(
            pre_reduction_variance.reshape(2, blur_vertices),
            axis=0,
        )

    output_vertices = scores.size
    smoothing_applied = bool(
        wscore_surface_smoothing_iterations > 0 and output_vertices == 32492
    )
    if smoothing_applied:
        scores = surface_smoother(
            scores,
            wscore_surface_smoothing_iterations,
        )

    # Stage-4 GP prediction-uncertainty filter: mask out vertices whose GP
    # predictive variance exceeds the requested within-map percentile so that
    # low-confidence vertices do not contribute to lesion detection.
    prediction_variance_threshold = None
    masked_uncertain_vertices = 0
    if prediction_variance_percentile is not None and output_variance.size == scores.size:
        finite_variance = output_variance[np.isfinite(output_variance)]
        if finite_variance.size:
            prediction_variance_threshold = float(
                np.percentile(finite_variance, float(prediction_variance_percentile))
            )
            high_uncertainty = np.isfinite(output_variance) & (
                output_variance > prediction_variance_threshold
            )
            scores = np.asarray(scores, dtype=float).copy()
            scores[high_uncertainty] = np.nan
            masked_uncertain_vertices = int(high_uncertainty.sum())
            if verbose:
                print(
                    f"    GP uncertainty filter: masked {masked_uncertain_vertices}/"
                    f"{scores.size} vertices with predictive variance > p"
                    f"{prediction_variance_percentile:g} ({prediction_variance_threshold:.4g})"
                )

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    data_array = nib.gifti.GiftiDataArray(
        data=np.asarray(scores, dtype=np.float32),
        intent="NIFTI_INTENT_NORMAL",
    )
    nib.save(nib.gifti.GiftiImage(darrays=[data_array]), output_file)

    if dominant_depth_indices is not None:
        dominant_depth_file = output_file.replace(
            ".func.gii",
            "_dominant-depth.shape.gii",
        )
        depth_array = nib.gifti.GiftiDataArray(
            data=dominant_depth_indices.astype(np.float32),
            intent="NIFTI_INTENT_NORMAL",
        )
        nib.save(nib.gifti.GiftiImage(darrays=[depth_array]), dominant_depth_file)

    finite_counts = fit.reference_counts[fit.reference_counts > 0]
    model_info = {
        "method": "gaussian_process",
        "score": "predictive_deviation_z",
        "wscore_covariate_model": wscore_covariate_model,
        "wscore_preprocessing": wscore_preprocessing,
        "blur_depth_model": (
            blur_depth_model
            if feature_metadata["input_was_multidepth"]
            and intensity_vertices is None
            else None
        ),
        "intensity_depth_model": intensity_depth_model,
        "wscore_surface_smoothing_iterations": int(
            wscore_surface_smoothing_iterations
        ),
        "wscore_surface_smoothing_applied": smoothing_applied,
        "uses_prediction_uncertainty": True,
        "uses_vertex_covariates": False,
        "reference_count": int(demo_ref.shape[0]),
        "reference_count_min": int(np.min(finite_counts)),
        "reference_count_max": int(np.max(finite_counts)),
        "dominant_depth_file": dominant_depth_file,
        "intensity_dominant_uncalibrated_control_p95": dominant_null_threshold,
        "intensity_dominant_calibration_factor": dominant_calibration_factor,
        "prediction_variance_percentile": prediction_variance_percentile,
        "prediction_variance_threshold": prediction_variance_threshold,
        "masked_uncertain_vertices": masked_uncertain_vertices,
        **fit.metadata(list(normative_columns)),
    }
    model_file = output_file.replace(".func.gii", "_model.json")
    with open(model_file, "w", encoding="utf-8") as stream:
        json.dump(model_info, stream, indent=2)

    if verbose:
        print(
            "    Saved Gaussian-process predictive deviation map: "
            f"{output_file}"
        )
        print(
            f"      Kernel: {GP_KERNEL_LABELS[fit.kernel_name]}; "
            f"LOO SD={fit.loo_score_sd:.3f}, "
            f"LOO |z|>1.96={fit.loo_abs_gt_1p96:.3%}"
        )

    finite_prediction = predictive_variances[np.isfinite(predictive_variances)]
    return {
        "wscore_file": output_file,
        "mean_wscore": float(np.nanmean(scores)),
        "std_wscore": float(np.nanstd(scores)),
        "normative_data": None,
        "reference_count": int(demo_ref.shape[0]),
        "reference_count_min": int(np.min(finite_counts)),
        "reference_count_max": int(np.max(finite_counts)),
        "vertex_covariate_count": 0,
        "uses_prediction_uncertainty": True,
        "wscore_distribution": "gaussian_process_predictive",
        "wscore_preprocessing": wscore_preprocessing,
        "wscore_covariate_model": wscore_covariate_model,
        "wscore_surface_smoothing_iterations": int(
            wscore_surface_smoothing_iterations
        ),
        "wscore_surface_smoothing_applied": smoothing_applied,
        "blur_depth_model": model_info["blur_depth_model"],
        "intensity_depth_model": intensity_depth_model,
        "dominant_depth_file": dominant_depth_file,
        "dominant_depth_indices": dominant_depth_indices,
        "predictive_variance_min": (
            float(np.min(finite_prediction)) if finite_prediction.size else np.nan
        ),
        "predictive_variance_max": (
            float(np.max(finite_prediction)) if finite_prediction.size else np.nan
        ),
        "prediction_variance_percentile": prediction_variance_percentile,
        "prediction_variance_threshold": prediction_variance_threshold,
        "masked_uncertain_vertices": masked_uncertain_vertices,
        "model_file": model_file,
    }


def calculate_gaussian_process_csv(
    reference_data,
    patient_data,
    demographics_ref,
    demographics_pat,
    output_file,
    normative_columns,
    verbose,
    gaussian_process_kernel,
    wscore_covariate_model,
    wscore_preprocessing,
    wscore_fit_cache,
):
    """Fit demographic GPs to subcortical structures and save deviation scores."""

    structures = [
        column
        for column in reference_data.columns
        if column not in {"structure", "SubjID"}
    ]
    reference_values = reference_data[structures].to_numpy(dtype=float)
    patient_values = patient_data[structures].iloc[0].to_numpy(dtype=float)

    if wscore_preprocessing != "none":
        robust = wscore_preprocessing == "spatial_robust_z"

        def normalize(values: np.ndarray) -> np.ndarray:
            one_dimensional = values.ndim == 1
            matrix = values[None, :] if one_dimensional else values
            if robust:
                center = np.nanmedian(matrix, axis=1, keepdims=True)
                scale = 1.482602218505602 * np.nanmedian(
                    np.abs(matrix - center), axis=1, keepdims=True
                )
            else:
                center = np.nanmean(matrix, axis=1, keepdims=True)
                scale = np.nanstd(matrix, axis=1, keepdims=True)
            scale = np.where(np.isfinite(scale) & (scale > EPS), scale, 1.0)
            result = (matrix - center) / scale
            return result[0] if one_dimensional else result

        reference_values = normalize(reference_values)
        patient_values = normalize(patient_values)

    demo_ref = demographics_ref[normative_columns].to_numpy(dtype=float)
    demo_pat = np.asarray(
        demographics_pat[normative_columns],
        dtype=float,
    ).reshape(-1)
    demographic_mask = np.isfinite(demo_ref).all(axis=1)
    demo_ref = demo_ref[demographic_mask]
    reference_values = reference_values[demographic_mask]

    signature = (
        id(reference_data),
        id(demographics_ref),
        tuple(normative_columns),
        wscore_preprocessing,
        gaussian_process_kernel,
    )
    cache_key = "gaussian_process_csv_current"
    cached = wscore_fit_cache.get(cache_key) if wscore_fit_cache is not None else None
    if cached is not None and cached.get("signature") == signature:
        fit = cached["fit"]
    else:
        fit = fit_gaussian_process(
            reference_values,
            demo_ref,
            kernel_name=gaussian_process_kernel,
        )
        cached = {"signature": signature, "fit": fit}
        if wscore_fit_cache is not None:
            wscore_fit_cache.clear()
            wscore_fit_cache[cache_key] = cached

    scores, predictive_variances = predict_gaussian_process(
        fit,
        patient_values,
        demo_pat,
    )
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    pd.DataFrame([{name: scores[index] for index, name in enumerate(structures)}]).to_csv(
        output_file,
        index=False,
    )
    structure_stats = {
        name: {
            "count": int(fit.reference_counts[index]),
            "predictive_variance": float(predictive_variances[index]),
            "kernel": GP_KERNEL_LABELS[fit.kernel_name],
        }
        for index, name in enumerate(structures)
    }
    if verbose:
        print(
            "    Saved subcortical Gaussian-process deviation CSV: "
            f"{output_file}"
        )
    return {
        "wscore_file": output_file,
        "wscores": dict(zip(structures, scores)),
        "structure_stats": structure_stats,
        "uses_prediction_uncertainty": True,
        "wscore_distribution": "gaussian_process_predictive",
        "wscore_preprocessing": wscore_preprocessing,
        "wscore_covariate_model": wscore_covariate_model,
        "effective_normative_columns": list(normative_columns),
        "gp_metadata": fit.metadata(list(normative_columns)),
    }

