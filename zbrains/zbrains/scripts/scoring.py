#!/usr/bin/env python
import os
import numpy as np
import nibabel as nib
import pandas as pd

# Snakemake variables
input_files = snakemake.input
output_files = snakemake.output
params = snakemake.params

# Helper: get single or first value from snakemake input/output/params
get = lambda x, k: x[k] if isinstance(x, dict) else x[0] if isinstance(x, (list, tuple)) else x

# Assign variables from snakemake
reference_data_files = input_files['reference_data'] if 'reference_data' in input_files else input_files[0]
patient_data_file = input_files['data_file'] if 'data_file' in input_files else input_files[1]
demographics_csv = input_files['demographics'] if 'demographics' in input_files else params.get('demographics_csv', None)
output_file = output_files['score_file'] if 'score_file' in output_files else output_files[0]

structure = params.get('structure', 'cortex')
method = params.get('method', 'z')
patient_id = params.get('patient_id', None)
normative_columns = params.get('normative_columns', ["age", "sex"])
verbose = params.get('verbose', True)
resolution = params.get('resolution', None)

# --- Scoring functions (as before) ---
def calculate_wscore_maps(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=["age", "sex"], verbose=True):
    if len(reference_data) == 0:
        raise ValueError("No reference data provided")
    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate W-scores")
    if demographics_ref.shape[0] != reference_data.shape[0]:
        raise ValueError("Demographics and reference data must have same number of subjects")
    if len(reference_data.shape) > 2:
        reference_data = np.mean(reference_data, axis=2)
    if len(patient_data.shape) > 1:
        patient_data = np.mean(patient_data, axis=1)
    n_vertices = reference_data.shape[1]
    n_predictors = len(normative_columns)
    normative_data = np.zeros((n_vertices, n_predictors + 2))
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[mask]
    ref_data_clean = reference_data[mask]
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    for i in range(n_vertices):
        vertex_data = ref_data_clean[:, i]
        if np.all(vertex_data == 0):
            normative_data[i, :] = 0
            normative_data[i, -1] = 1
            continue
        coefficients = np.linalg.lstsq(X_ref_with_intercept, vertex_data, rcond=None)[0]
        predicted = X_ref_with_intercept @ coefficients
        residuals = vertex_data - predicted
        std_residuals = np.std(residuals)
        normative_data[i, :-1] = coefficients
        normative_data[i, -1] = std_residuals
    wscores = np.zeros(n_vertices)
    for i in range(n_vertices):
        coefficients = normative_data[i, :-1]
        std_residuals = normative_data[i, -1]
        expected = coefficients[0]
        for j in range(n_predictors):
            expected += coefficients[j+1] * X_pat[j]
        wscores[i] = (patient_data[i] - expected) / std_residuals
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    wscore_data = nib.gifti.gifti.GiftiDataArray(
        data=wscores.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    wscore_gii = nib.gifti.GiftiImage(darrays=[wscore_data])
    nib.save(wscore_gii, output_file)
    if verbose:
        print(f"    Saved w-score map: {output_file}")

def calculate_zscore_maps(reference_data, patient_data, output_file, verbose=True):
    if len(reference_data) == 0:
        raise ValueError("No reference data provided")
    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate z-scores")
    ref_mean = np.mean(reference_data, axis=0)
    ref_std = np.std(reference_data, axis=0)
    ref_std[ref_std == 0] = 1.0
    zscore = (patient_data - ref_mean) / ref_std
    if len(zscore.shape) > 1:
        mean_zscore = np.mean(zscore, axis=1)
        zscore_to_save = mean_zscore
    else:
        zscore_to_save = zscore
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    zscore_data = nib.gifti.gifti.GiftiDataArray(
        data=zscore_to_save.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    zscore_gii = nib.gifti.GiftiImage(darrays=[zscore_data])
    nib.save(zscore_gii, output_file)
    if verbose:
        print(f"    Saved z-score map: {output_file}")

def calculate_wscore_csv(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=["age", "sex"], verbose=True):
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[mask]
    ref_data_clean = reference_data.loc[mask]
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    patient_wscores = {}
    structure_stats = {}
    for structure in structures:
        if structure in reference_data.columns and structure in patient_data.columns:
            ref_values = ref_data_clean[structure].values
            pat_value = patient_data[structure].iloc[0]
            if np.all(ref_values == 0):
                patient_wscores[structure] = 0
                structure_stats[structure] = {
                    'intercept': 0, 'coefficients': [0] * len(normative_columns),
                    'std_residuals': 1, 'count': len(ref_values)
                }
                continue
            coefficients = np.linalg.lstsq(X_ref_with_intercept, ref_values, rcond=None)[0]
            predicted = X_ref_with_intercept @ coefficients
            residuals = ref_values - predicted
            std_residuals = np.std(residuals)
            if std_residuals == 0:
                std_residuals = 1
            expected = coefficients[0] + np.sum([coefficients[j+1] * X_pat[j] for j in range(len(normative_columns))])
            wscore = (pat_value - expected) / std_residuals
            patient_wscores[structure] = wscore
            structure_stats[structure] = {
                'intercept': coefficients[0],
                'coefficients': coefficients[1:].tolist(),
                'std_residuals': std_residuals,
                'count': len(ref_values)
            }
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    wscore_df = pd.DataFrame([patient_wscores])
    wscore_df.to_csv(output_file, index=False)
    if verbose:
        print(f"    Saved w-score CSV: {output_file}")

def calculate_zscore_csv(reference_data, patient_data, output_file, verbose=True):
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)
    structure_stats = {}
    for structure in structures:
        if structure in reference_data.columns:
            ref_values = reference_data[structure].values
            structure_stats[structure] = {
                'mean': np.mean(ref_values),
                'std': np.std(ref_values) if np.std(ref_values) > 0 else 1.0,
                'count': len(ref_values)
            }
    patient_zscores = {}
    for structure in structures:
        if structure in patient_data.columns and structure in structure_stats:
            pat_value = patient_data[structure].iloc[0]
            ref_mean = structure_stats[structure]['mean']
            ref_std = structure_stats[structure]['std']
            zscore = (pat_value - ref_mean) / ref_std
            patient_zscores[structure] = zscore
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    zscore_df = pd.DataFrame([patient_zscores])
    zscore_df.to_csv(output_file, index=False)
    if verbose:
        print(f"    Saved z-score CSV: {output_file}")

# --- Main logic ---
# Load demographics
if demographics_csv is None:
    raise ValueError("demographics_csv must be provided via snakemake.input or snakemake.params")
demographics = pd.read_csv(demographics_csv)
if patient_id is None:
    raise ValueError("patient_id must be provided via snakemake.params")
if patient_id not in demographics["ID"].values:
    raise ValueError(f"Patient ID {patient_id} not found in demographics CSV")
demographics_ref = demographics
demographics_pat = demographics[demographics["ID"] == patient_id].iloc[0]

# Load reference data
if structure == "subcortical":
    reference_data = [pd.read_csv(f) for f in reference_data_files]
    reference_data = pd.concat(reference_data, ignore_index=True)
    patient_data = pd.read_csv(patient_data_file)
    if method == "w":
        calculate_wscore_csv(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=normative_columns, verbose=verbose)
    else:
        calculate_zscore_csv(reference_data, patient_data, output_file, verbose=verbose)
elif structure in ["cortex", "hippocampus"]:
    def load_gifti_or_npy(path):
        if path.endswith(".npy"):
            return np.load(path)
        else:
            img = nib.load(path)
            return np.stack([d.data for d in img.darrays], axis=0) if len(img.darrays) > 1 else img.darrays[0].data[np.newaxis, :]
    reference_data = np.stack([load_gifti_or_npy(f) for f in reference_data_files])
    reference_data = np.squeeze(reference_data)
    patient_data = load_gifti_or_npy(patient_data_file)
    if method == "w":
        calculate_wscore_maps(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=normative_columns, verbose=verbose)
    else:
        calculate_zscore_maps(reference_data, patient_data, output_file, verbose=verbose)
else:
    print(f"Unsupported structure: {structure}")
    exit(1) 