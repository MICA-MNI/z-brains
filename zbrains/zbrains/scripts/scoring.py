#!/usr/bin/env python3
"""
scoring.py: Compute z- or w-scores for patient surface/subcortical maps against a reference cohort.

Inputs (snakemake.input):
  - data_file: patient .func.gii map
  - reference_data: list of reference .func.gii maps
  - demographics_ref: CSV of reference subjects with columns ['ID', ... normative_columns]
  - demographics_pat: CSV or single-row CSV for the patient subject
Params (snakemake.params):
  - method: 'z' or 'w'
  - normative_columns (optional): list of column names in demographics for regression (default ['age','sex'])
Outputs (snakemake.output):
  - out .func.gii file with per-vertex/structure scores
"""
import os
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.gifti import GiftiImage, GiftiDataArray


def load_gifti(path):
    """Load a GIFTI .func.gii and return a 1D or 2D numpy array."""
    img = nib.load(path)
    arrays = [da.data for da in img.darrays]
    return np.stack(arrays, axis=0) if len(arrays) > 1 else arrays[0]


def calculate_zscore(ref_mat, pat_vec):
    """Compute z-score: (x - mean) / std per feature."""
    mean = np.mean(ref_mat, axis=0)
    std = np.std(ref_mat, axis=0)
    std[std == 0] = 1.0
    return (pat_vec - mean) / std


def calculate_wscore(ref_mat, pat_vec, X_ref, X_pat):
    """Compute w-score via linear regression residuals and patient predictors."""
    n_ref = ref_mat.shape[0]
    # build design with intercept
    X = np.hstack([np.ones((n_ref, 1)), X_ref.astype(float)])
    coefs = np.zeros((ref_mat.shape[1], X.shape[1]))
    stds = np.ones(ref_mat.shape[1], dtype=float)
    for i in range(ref_mat.shape[1]):
        y = ref_mat[:, i]
        if not np.all(y == 0):
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            res = y - X.dot(coef)
            std = np.std(res)
            coefs[i, :] = coef
            stds[i] = std if std > 0 else 1.0
    # expected values for patient
    intercepts = coefs[:, 0]
    slopes = coefs[:, 1:]
    expected = intercepts + slopes.dot(X_pat.astype(float))
    return (pat_vec - expected) / stds


def main():
    # Snakemake variables
    refs = snakemake.input.reference_data
    patient_map = snakemake.input.data_file
    demo_ref_csv = snakemake.input.demographics_ref
    demo_pat_csv = snakemake.input.demographics_pat
    method = snakemake.params.method
    normative_cols = snakemake.params.get('normative_columns', ['age', 'sex'])
    out_gii = snakemake.output[0]

    # Load reference cohort data
    ref_list = [load_gifti(f) for f in refs]
    ref_mat = np.stack(ref_list, axis=0)
    if ref_mat.ndim > 2:
        # collapse trailing dimensions
        ref_mat = ref_mat.mean(axis=tuple(range(2, ref_mat.ndim)))

    # Load patient data
    pat_data = load_gifti(patient_map)
    if pat_data.ndim > 1:
        pat_data = pat_data.mean(axis=0)

    # Compute scores
    if method.lower().startswith('w'):
        # Load separate demographics
        df_ref = pd.read_csv(demo_ref_csv)
        df_pat = pd.read_csv(demo_pat_csv)
        if df_pat.shape[0] != 1:
            raise ValueError("Patient demographics CSV must contain exactly one subject row.")
        X_ref = df_ref[normative_cols].values
        X_pat = df_pat[normative_cols].values.flatten()
        scores = calculate_wscore(ref_mat, pat_data, X_ref, X_pat)
    else:
        scores = calculate_zscore(ref_mat, pat_data)

    # Save output .func.gii
    os.makedirs(os.path.dirname(out_gii), exist_ok=True)
    da = GiftiDataArray(scores.astype(np.float32))
    nib.save(GiftiImage(darrays=[da]), out_gii)
    print(f"Saved {method}-score map to {out_gii}")


if __name__ == '__main__':
    main()
