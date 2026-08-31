import json
import tempfile
import unittest
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from zbrains.analysis import (
    _bilateral_multisurface_median_abs_scale,
    _bilateral_blur_gradient_flattening_components,
    _blur_gradient_flattening_components,
    _filter_reference_controls_to_subjects,
    _fit_shash_quantile_parameters,
    _normalize_wscore_covariate_model,
    _normalize_wscore_distribution,
    _normalize_wscore_preprocessing,
    _normalize_intensity_depth_model,
    _normalize_blur_depth_model,
    _spatial_zscore,
    _smooth_fslr32k_scores,
    _shash_to_normal_score,
    _white_swm1_direction_cosine,
    calculate_wscore_csv,
    calculate_wscore_maps,
)


def shash_from_normal(z, location, scale, skew, tail_weight):
    return location + scale * np.sinh((np.arcsinh(z) + skew) / tail_weight)


class TestShashWscore(unittest.TestCase):
    def test_distribution_validation_and_aliases(self):
        self.assertEqual(_normalize_wscore_distribution("normal"), "gaussian")
        self.assertEqual(_normalize_wscore_distribution("sinh-arcsinh"), "shash")
        with self.assertRaises(ValueError):
            _normalize_wscore_distribution("not-a-distribution")
        self.assertEqual(_normalize_wscore_distribution("mad"), "gaussian_mad")
        self.assertEqual(_normalize_wscore_distribution("winsorized"), "gaussian_winsor")
        self.assertEqual(_normalize_wscore_distribution("rank-normal"), "empirical")
        self.assertEqual(_normalize_wscore_preprocessing("spatial-z"), "spatial_robust_z")
        self.assertEqual(_normalize_wscore_preprocessing("zscore"), "spatial_zscore")
        self.assertEqual(
            _normalize_wscore_preprocessing("robust-zscore"), "spatial_robust_z"
        )

    def test_spatial_zscore_standardizes_each_subject_map(self):
        values = np.array([[1.0, 2.0, 3.0], [10.0, 14.0, 18.0]])
        standardized = _spatial_zscore(values)
        np.testing.assert_allclose(np.mean(standardized, axis=1), 0.0, atol=1e-12)
        np.testing.assert_allclose(np.std(standardized, axis=1), 1.0, atol=1e-12)

    def test_white_swm1_direction_cosine_is_scale_invariant(self):
        depths = np.array(
            [
                [2.0, 4.0, 1.0, 2.0],
                [1.0, 2.0, 3.0, 4.0],
                [0.5, 1.0, 0.5, 0.25],
            ]
        )
        observed = _white_swm1_direction_cosine(depths, floor_percentile=0.0)
        expected = (depths[:, 1] - depths[:, 2]) / np.sqrt(
            np.sum(depths**2, axis=1)
        )
        np.testing.assert_allclose(observed, expected)
        np.testing.assert_allclose(
            _white_swm1_direction_cosine(7.0 * depths, floor_percentile=0.0),
            observed,
        )
        self.assertEqual(
            _normalize_intensity_depth_model("depth-cosine"),
            "white_swm1_direction_cosine",
        )
        with self.assertRaises(ValueError):
            _white_swm1_direction_cosine(np.ones((3, 3)))

    def test_multisurface_median_abs_dominant_score_and_calibration(self):
        rng = np.random.default_rng(79)
        left = rng.normal(size=(6, 4))
        right = rng.normal(size=(6, 4))
        scaled_left, scaled_right = _bilateral_multisurface_median_abs_scale(
            left, right
        )
        scaled_left_2, scaled_right_2 = _bilateral_multisurface_median_abs_scale(
            9.0 * left, 9.0 * right
        )
        np.testing.assert_allclose(scaled_left, scaled_left_2)
        np.testing.assert_allclose(scaled_right, scaled_right_2)

        n_controls, n_vertices = 36, 6
        ages = np.linspace(20.0, 70.0, n_controls)
        reference = rng.normal(size=(n_controls, n_vertices, 4))
        reference += 0.01 * ages[:, None, None]
        patient_age = 45.0
        patient = rng.normal(size=(n_vertices, 4)) + 0.01 * patient_age
        design = np.column_stack([np.ones(n_controls), ages])
        patient_design = np.array([1.0, patient_age])
        flat_reference = reference.reshape(n_controls, -1)
        coefficients = np.linalg.lstsq(design, flat_reference, rcond=None)[0]
        residuals = flat_reference - design @ coefficients
        residual_sd = np.std(residuals, axis=0)
        control_scores = (residuals / residual_sd).reshape(
            n_controls, n_vertices, 4
        )
        calibration = (
            np.percentile(np.max(np.abs(control_scores), axis=2), 95.0) / 1.96
        )
        patient_scores = (
            (patient.reshape(-1) - patient_design @ coefficients) / residual_sd
        ).reshape(n_vertices, 4)
        dominant = np.argmax(np.abs(patient_scores), axis=1)
        expected = np.take_along_axis(
            patient_scores, dominant[:, None], axis=1
        )[:, 0] / calibration

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "intensity.func.gii")
            result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"age": ages}),
                pd.Series({"age": patient_age}),
                output_file,
                normative_columns=["age"],
                verbose=False,
                intensity_depth_model="multisurface_median_abs_dominant",
            )
            observed = np.asarray(nib.load(output_file).darrays[0].data)
            observed_depth = np.asarray(
                nib.load(result["dominant_depth_file"]).darrays[0].data
            )
            metadata = json.loads(
                Path(output_file.replace(".func.gii", "_model.json")).read_text()
            )

        np.testing.assert_allclose(observed, expected, rtol=1e-5, atol=1e-6)
        np.testing.assert_array_equal(observed_depth, dominant)
        self.assertEqual(
            metadata["intensity_depth_model"],
            "multisurface_median_abs_dominant",
        )

    def test_blur_mean_slope_rms_matches_independent_component_models(self):
        rng = np.random.default_rng(77)
        n_controls, n_vertices, n_depths = 48, 7, 4
        ages = np.linspace(20.0, 70.0, n_controls)
        reference = rng.normal(size=(n_controls, n_vertices, n_depths))
        reference += 0.01 * ages[:, None, None]
        patient_age = 44.0
        patient = rng.normal(size=(n_vertices, n_depths)) + 0.01 * patient_age

        depth = np.arange(n_depths, dtype=float)
        depth -= depth.mean()
        slope_weights = depth / np.sum(depth**2)
        reference_components = (
            np.mean(reference, axis=2),
            reference @ slope_weights,
        )
        patient_components = (
            np.mean(patient, axis=1),
            patient @ slope_weights,
        )
        design = np.column_stack([np.ones(n_controls), ages])
        patient_design = np.array([1.0, patient_age])
        expected_components = []
        for control_component, patient_component in zip(
            reference_components, patient_components
        ):
            coefficients = np.linalg.lstsq(design, control_component, rcond=None)[0]
            residuals = control_component - design @ coefficients
            expected_components.append(
                (patient_component - patient_design @ coefficients)
                / np.std(residuals, axis=0)
            )
        expected_components = np.asarray(expected_components)
        expected = np.copysign(
            np.sqrt(np.mean(expected_components**2, axis=0)),
            expected_components[0],
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "blur.func.gii")
            result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"age": ages}),
                pd.Series({"age": patient_age}),
                output_file,
                normative_columns=["age"],
                verbose=False,
                blur_depth_model="mean_slope_rms",
            )
            observed = np.asarray(nib.load(output_file).darrays[0].data)
            metadata = json.loads(
                Path(output_file.replace(".func.gii", "_model.json")).read_text()
            )

        np.testing.assert_allclose(observed, expected, rtol=1e-5, atol=1e-6)
        self.assertEqual(result["blur_depth_model"], "mean_slope_rms")
        self.assertEqual(metadata["blur_depth_model"], "mean_slope_rms")
        self.assertEqual(observed.shape, (n_vertices,))

    def test_gradient_flattening_components_and_directional_wscore(self):
        depths = np.array(
            [
                [2.0, 6.0, 4.0, 3.0],
                [1.0, 4.0, 2.0, 1.0],
                [3.0, 5.0, 4.5, 4.0],
            ]
        )
        distance = np.array([2.0, 1.5, 1.0])
        safe_distance = np.maximum(distance, np.percentile(distance, 1.0))
        expected_gradients = np.column_stack(
            [
                np.abs(depths[:, 1] - depths[:, 0]) / safe_distance,
                np.abs(depths[:, 2] - depths[:, 1]),
            ]
        )
        np.testing.assert_allclose(
            _blur_gradient_flattening_components(depths, distance),
            expected_gradients,
        )
        right_depths = depths + np.array([0.5, 1.0, -0.5, 0.25])
        right_distance = np.array([1.5, 2.0, 1.25])
        scaled_components = _bilateral_blur_gradient_flattening_components(
            depths,
            right_depths,
            distance,
            right_distance,
        )
        rescaled_components = _bilateral_blur_gradient_flattening_components(
            9.0 * depths,
            9.0 * right_depths,
            distance,
            right_distance,
        )
        np.testing.assert_allclose(
            scaled_components[0], rescaled_components[0]
        )
        np.testing.assert_allclose(
            scaled_components[1], rescaled_components[1]
        )
        self.assertEqual(
            _normalize_blur_depth_model("gw-swm1-gradient-flattening"),
            "gradient_flattening",
        )

        rng = np.random.default_rng(177)
        n_controls, n_vertices = 40, 6
        ages = np.linspace(20.0, 70.0, n_controls)
        reference = rng.normal(size=(n_controls, n_vertices, 2))
        reference += 0.01 * ages[:, None, None]
        patient_age = 43.0
        patient = rng.normal(size=(n_vertices, 2)) + 0.01 * patient_age
        design = np.column_stack([np.ones(n_controls), ages])
        patient_design = np.array([1.0, patient_age])
        component_scores = []
        for component in range(2):
            controls = reference[:, :, component]
            coefficients = np.linalg.lstsq(design, controls, rcond=None)[0]
            residuals = controls - design @ coefficients
            component_scores.append(
                (patient[:, component] - patient_design @ coefficients)
                / np.std(residuals, axis=0)
            )
        expected = -np.mean(component_scores, axis=0)

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "gradient_flattening.func.gii")
            result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"age": ages}),
                pd.Series({"age": patient_age}),
                output_file,
                normative_columns=["age"],
                verbose=False,
                blur_depth_model="gradient_flattening",
            )
            observed = np.asarray(nib.load(output_file).darrays[0].data)
            metadata = json.loads(
                Path(output_file.replace(".func.gii", "_model.json")).read_text()
            )

        np.testing.assert_allclose(observed, expected, rtol=1e-5, atol=1e-6)
        self.assertEqual(result["blur_depth_model"], "gradient_flattening")
        self.assertEqual(
            metadata["blur_depth_components"],
            [
                "distance_corrected_gray_white_gradient_magnitude",
                "white_to_swm1_gradient_magnitude",
            ],
        )
        self.assertIn("-(gray_white_gradient_wscore", metadata["blur_depth_combination"])
        self.assertEqual(
            metadata["blur_depth_self_normalization"],
            "bilateral_median_absolute_intensity_across_all_four_depths",
        )

    def test_blur_depth_models_support_vertex_covariates_and_legacy_mean(self):
        rng = np.random.default_rng(78)
        n_controls, n_vertices = 24, 5
        ages = np.linspace(25.0, 65.0, n_controls)
        covariates = rng.normal(size=(n_controls, n_vertices))
        patient_covariates = rng.normal(size=n_vertices)
        reference = rng.normal(size=(n_controls, n_vertices, 4))
        reference += 0.2 * covariates[:, :, None]
        patient = rng.normal(size=(n_vertices, 4))
        patient += 0.2 * patient_covariates[:, None]

        with tempfile.TemporaryDirectory() as tmp_dir:
            for model in ("mean_slope_rms", "mean"):
                output_file = str(Path(tmp_dir) / f"{model}.func.gii")
                result = calculate_wscore_maps(
                    reference,
                    patient,
                    pd.DataFrame({"age": ages}),
                    pd.Series({"age": 40.0}),
                    output_file,
                    normative_columns=["age"],
                    verbose=False,
                    blur_depth_model=model,
                    reference_vertex_covariates=covariates,
                    patient_vertex_covariates=patient_covariates,
                )
                scores = np.asarray(nib.load(output_file).darrays[0].data)
                self.assertEqual(scores.shape, (n_vertices,))
                self.assertTrue(np.isfinite(scores).all())
                self.assertEqual(result["blur_depth_model"], model)

    def test_robust_residual_models_and_spatial_preprocessing(self):
        rng = np.random.default_rng(101)
        n_controls = 80
        ages = rng.normal(size=n_controls)
        reference = 0.2 * ages[:, None] + rng.normal(size=(n_controls, 12))
        reference[0, :3] += 25.0
        patient = 0.2 * 0.4 + rng.normal(size=12)

        with tempfile.TemporaryDirectory() as tmp_dir:
            for distribution in (
                "gaussian_mad",
                "gaussian_winsor",
                "student_t",
                "empirical",
            ):
                output_file = str(Path(tmp_dir) / f"{distribution}.func.gii")
                result = calculate_wscore_maps(
                    reference,
                    patient,
                    pd.DataFrame({"age": ages}),
                    pd.Series({"age": 0.4}),
                    output_file,
                    normative_columns=["age"],
                    verbose=False,
                    wscore_distribution=distribution,
                    wscore_preprocessing="spatial_robust_z",
                )
                scores = nib.load(output_file).darrays[0].data
                self.assertTrue(np.isfinite(scores).all())
                self.assertEqual(result["wscore_distribution"], distribution)
                self.assertEqual(result["wscore_preprocessing"], "spatial_robust_z")

            standard_output = str(Path(tmp_dir) / "spatial_zscore.func.gii")
            standard_result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"age": ages}),
                pd.Series({"age": 0.4}),
                standard_output,
                normative_columns=["age"],
                verbose=False,
                wscore_preprocessing="spatial_zscore",
            )
            self.assertTrue(
                np.isfinite(nib.load(standard_output).darrays[0].data).all()
            )
            self.assertEqual(
                standard_result["wscore_preprocessing"], "spatial_zscore"
            )

    def test_fslr_surface_smoothing_preserves_constant_scores(self):
        scores = np.full(32492, 2.5)
        smoothed = _smooth_fslr32k_scores(scores, 3)
        np.testing.assert_allclose(smoothed, scores, atol=1e-12)

    def test_quadratic_age_covariate_removes_curved_age_effect(self):
        rng = np.random.default_rng(202)
        ages = np.linspace(18.0, 70.0, 100)
        age_z = (ages - ages.mean()) / ages.std()
        reference = (3.0 + 0.2 * ages + 1.5 * age_z**2)[:, None]
        reference += rng.normal(scale=0.05, size=reference.shape)
        patient_age = 55.0
        patient_age_z = (patient_age - ages.mean()) / ages.std()
        patient = np.array([3.0 + 0.2 * patient_age + 1.5 * patient_age_z**2])

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "quadratic.func.gii")
            result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"AGE": ages}),
                pd.Series({"AGE": patient_age}),
                output_file,
                normative_columns=["AGE"],
                verbose=False,
                wscore_covariate_model="quadratic_age",
            )
            score = float(nib.load(output_file).darrays[0].data[0])

        self.assertAlmostEqual(score, 0.0, delta=0.5)
        self.assertEqual(result["wscore_covariate_model"], "quadratic_age")

    def test_regional_qc_mask_keeps_asymmetry_inputs_aligned(self):
        subjects = [("sub-01", "ses-01"), ("sub-02", "ses-01"), ("sub-03", "ses-01")]
        data = np.arange(12, dtype=float).reshape(3, 4)
        demographics = pd.DataFrame({"age": [20, 30, 40]})
        covariates = np.arange(24, dtype=float).reshape(3, 4, 2)

        filtered_data, filtered_subjects, filtered_demo, filtered_covariates = (
            _filter_reference_controls_to_subjects(
                data,
                subjects,
                {subjects[0], subjects[2]},
                reference_demographics_df=demographics,
                reference_vertex_covariates=covariates,
                verbose=False,
            )
        )

        self.assertEqual(filtered_subjects, [subjects[0], subjects[2]])
        np.testing.assert_array_equal(filtered_data, data[[0, 2]])
        np.testing.assert_array_equal(filtered_demo["age"].to_numpy(), [20, 40])
        np.testing.assert_array_equal(filtered_covariates, covariates[[0, 2]])

    def test_forward_and_inverse_parameterization(self):
        parameters = (1.2, 2.3, 0.7, 0.75)
        for expected_score in (-2.0, 0.0, 2.0):
            value = shash_from_normal(expected_score, *parameters)
            observed_score = _shash_to_normal_score(value, *parameters)
            self.assertAlmostEqual(observed_score, expected_score, places=10)

    def test_quantile_fit_recovers_latent_scores(self):
        rng = np.random.default_rng(4)
        parameters = (1.2, 2.3, 0.7, 0.75)
        residuals = shash_from_normal(rng.normal(size=5000), *parameters)
        fitted, success = _fit_shash_quantile_parameters(residuals)

        self.assertTrue(success[0])
        for expected_score in (-2.0, 0.0, 2.0):
            value = shash_from_normal(expected_score, *parameters)
            observed_score = _shash_to_normal_score(value, *fitted[0])
            self.assertAlmostEqual(observed_score, expected_score, delta=0.15)

    def test_surface_wscore_uses_shash_latent_normal_scale(self):
        rng = np.random.default_rng(18)
        n_controls = 2000
        ages = rng.normal(size=n_controls)
        design = np.column_stack([np.ones(n_controls), ages])
        coefficients = np.array([[2.0, -1.0, 0.5], [0.4, 0.2, -0.3]])
        parameters = np.array(
            [
                [0.0, 1.0, 0.7, 0.75],
                [1.0, 2.0, -0.5, 1.35],
                [-0.5, 0.8, 0.2, 0.60],
            ]
        )
        residuals = np.column_stack(
            [
                shash_from_normal(rng.normal(size=n_controls), *vertex_parameters)
                for vertex_parameters in parameters
            ]
        )
        reference_data = design @ coefficients + residuals

        patient_age = 0.3
        expected_scores = np.array([0.0, 1.5, -2.0])
        patient_residuals = np.array(
            [
                shash_from_normal(score, *vertex_parameters)
                for score, vertex_parameters in zip(expected_scores, parameters)
            ]
        )
        patient_data = np.array([1.0, patient_age]) @ coefficients + patient_residuals

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "shash.func.gii")
            result = calculate_wscore_maps(
                reference_data,
                patient_data,
                pd.DataFrame({"age": ages}),
                pd.Series({"age": patient_age}),
                output_file,
                normative_columns=["age"],
                verbose=False,
                wscore_distribution="shash",
            )
            observed_scores = nib.load(output_file).darrays[0].data
            with open(result["model_file"]) as model_stream:
                model_info = json.load(model_stream)

        np.testing.assert_allclose(observed_scores, expected_scores, atol=0.25)
        self.assertEqual(result["wscore_distribution"], "shash")
        self.assertTrue(result["shash_fit_success"].all())
        self.assertEqual(model_info["wscore_distribution"], "shash")
        self.assertEqual(model_info["shash_fit_method"], "quantile_matching_10_25_50_75_90")

    def test_subcortical_csv_supports_shash(self):
        rng = np.random.default_rng(31)
        n_controls = 1500
        ages = rng.normal(size=n_controls)
        parameters = (0.5, 1.4, -0.6, 0.8)
        control_values = (
            3.0
            + 0.25 * ages
            + shash_from_normal(rng.normal(size=n_controls), *parameters)
        )
        patient_age = -0.4
        expected_score = 1.75
        patient_value = (
            3.0
            + 0.25 * patient_age
            + shash_from_normal(expected_score, *parameters)
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "shash.csv")
            result = calculate_wscore_csv(
                pd.DataFrame({"thalamus": control_values}),
                pd.DataFrame({"thalamus": [patient_value]}),
                pd.DataFrame({"age": ages}),
                pd.Series({"age": patient_age}),
                output_file,
                normative_columns=["age"],
                verbose=False,
                wscore_distribution="shash",
            )
            saved_score = pd.read_csv(output_file).loc[0, "thalamus"]

        self.assertAlmostEqual(saved_score, expected_score, delta=0.25)
        self.assertEqual(result["wscore_distribution"], "shash")
        self.assertTrue(result["structure_stats"]["thalamus"]["shash_fit_success"])
    def test_gaussian_process_covariate_aliases(self):
        self.assertEqual(
            _normalize_wscore_covariate_model("gp"),
            "gaussian_process",
        )
        self.assertEqual(
            _normalize_wscore_covariate_model("gpr"),
            "gaussian_process",
        )
        self.assertEqual(
            _normalize_wscore_covariate_model("gp_matern32"),
            "gaussian_process_matern32",
        )
        self.assertEqual(
            _normalize_wscore_covariate_model("gp_matern52"),
            "gaussian_process_matern52",
        )
        self.assertEqual(
            _normalize_wscore_covariate_model("gp_rq"),
            "gaussian_process_rational_quadratic",
        )

    def test_gaussian_process_kernel_variants(self):
        rng = np.random.default_rng(705)
        ages = np.linspace(20.0, 72.0, 40)
        sex = np.tile([0.0, 1.0], 20)
        expectation = np.sin(ages / 10.0) + 0.3 * sex
        reference = np.column_stack(
            [
                expectation + rng.normal(scale=0.12, size=ages.size),
                -0.6 * expectation + rng.normal(scale=0.10, size=ages.size),
            ]
        )
        patient_age = 46.0
        patient_sex = 1.0
        patient_expectation = np.sin(patient_age / 10.0) + 0.3 * patient_sex
        patient = np.array([patient_expectation, -0.6 * patient_expectation])
        variants = {
            "gaussian_process_rbf_isotropic": "isotropic_RBF_plus_white_noise",
            "gaussian_process_matern32": "ARD_Matern32_plus_white_noise",
            "gaussian_process_matern52": "ARD_Matern52_plus_white_noise",
            "gaussian_process_rational_quadratic": (
                "ARD_rational_quadratic_plus_white_noise"
            ),
        }

        with tempfile.TemporaryDirectory() as tmp_dir:
            for model, expected_kernel in variants.items():
                with self.subTest(model=model):
                    output_file = str(Path(tmp_dir) / f"{model}.func.gii")
                    result = calculate_wscore_maps(
                        reference,
                        patient,
                        pd.DataFrame({"AGE": ages, "SEX": sex}),
                        pd.Series({"AGE": patient_age, "SEX": patient_sex}),
                        output_file,
                        normative_columns=["AGE", "SEX"],
                        verbose=False,
                        wscore_covariate_model=model,
                    )
                    scores = nib.load(output_file).darrays[0].data
                    metadata = json.loads(
                        Path(
                            output_file.replace(".func.gii", "_model.json")
                        ).read_text()
                    )
                    self.assertTrue(np.isfinite(scores).all())
                    self.assertEqual(result["wscore_covariate_model"], model)
                    self.assertEqual(metadata["wscore_covariate_model"], model)
                    self.assertEqual(metadata["kernel"], expected_kernel)

    def test_gaussian_process_scores_nonlinear_demographic_expectation(self):
        rng = np.random.default_rng(704)
        ages = np.linspace(18.0, 75.0, 72)
        sex = np.tile([0.0, 1.0], 36)
        base = (
            2.0
            + np.sin(ages / 9.0)
            + 0.35 * sex
            + 0.25 * sex * np.cos(ages / 11.0)
        )
        reference = np.column_stack(
            [
                base + rng.normal(scale=0.08, size=ages.size),
                0.5 * base + rng.normal(scale=0.06, size=ages.size),
                -base + rng.normal(scale=0.10, size=ages.size),
            ]
        )
        patient_age = 47.0
        patient_sex = 1.0
        expected = (
            2.0
            + np.sin(patient_age / 9.0)
            + 0.35 * patient_sex
            + 0.25 * patient_sex * np.cos(patient_age / 11.0)
        )
        patient = np.array([expected, 0.5 * expected, -expected])

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = str(Path(tmp_dir) / "gp.func.gii")
            result = calculate_wscore_maps(
                reference,
                patient,
                pd.DataFrame({"AGE": ages, "SEX": sex}),
                pd.Series({"AGE": patient_age, "SEX": patient_sex}),
                output_file,
                normative_columns=["AGE", "SEX"],
                verbose=False,
                wscore_covariate_model="gaussian_process",
            )
            scores = nib.load(output_file).darrays[0].data
            metadata = json.loads(
                Path(output_file.replace(".func.gii", "_model.json")).read_text()
            )

        self.assertTrue(np.isfinite(scores).all())
        self.assertLess(float(np.max(np.abs(scores))), 1.5)
        self.assertEqual(result["wscore_covariate_model"], "gaussian_process")
        self.assertTrue(result["uses_prediction_uncertainty"])
        self.assertEqual(metadata["method"], "gaussian_process")
        self.assertEqual(metadata["kernel"], "ARD_RBF_plus_white_noise")



if __name__ == "__main__":
    unittest.main()
