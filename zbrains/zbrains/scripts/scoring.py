#!/usr/bin/env python3
"""
scoring.py: Compute z‐ or w‐scores for a patient .func.gii against
a pickled reference‐matrix (one file per structure/feature/hemi).

Snakemake inputs:
  snakemake.input.data_file     – patient .func.gii
  snakemake.input.ref_mat       – reference numpy array (.pkl)
  snakemake.input.demographics_ref – CSV of reference cohort
  snakemake.input.demographics_pat – single‐row CSV of patient

Snakemake params:
  snakemake.params.method            – 'z' or 'w'
  snakemake.params.normative_columns – list of columns, e.g. ['AGE','SEX']

Snakemake output:
  snakemake.output[0] – scored .func.gii
"""
import os
import pickle
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.gifti import GiftiImage, GiftiDataArray

def load_gifti(path):
    """Load a GIFTI .func.gii and return a 1D vector."""
    img = nib.load(path)
    arrays = [da.data for da in img.darrays]
    arr = np.stack(arrays, axis=0) if len(arrays) > 1 else arrays[0]
    # if it’s multi‐dimensional beyond [n_samples, n_vertices], collapse extras
    if arr.ndim > 1:
        arr = arr.mean(axis=0)
    return arr

def calculate_zscore(ref_mat, pat_vec):
    """(x - mean) / std per feature."""
    mu  = np.mean(ref_mat, axis=0)
    sig = np.std(ref_mat, axis=0)
    sig[sig == 0] = 1.0
    return (pat_vec - mu) / sig

def calculate_wscore(ref_mat, pat_vec, X_ref, X_pat):
    """
    Standardized residual from linear regression:
      (observed - predicted) / residual_std
    Assumes X_ref: n_ref × n_preds, X_pat: n_preds,
    and ref_mat: n_ref × n_features.
    """
    n_ref, n_feat = ref_mat.shape

    # sanity checks
    if X_ref.shape[0] != n_ref:
        raise ValueError(f"X_ref rows ({X_ref.shape[0]}) ≠ ref_mat rows ({n_ref})")
    if X_pat.ndim != 1 or X_pat.shape[0] != X_ref.shape[1]:
        raise ValueError(f"X_pat must be length {X_ref.shape[1]}; got {X_pat.shape}")

    # design matrix with intercept
    X = np.column_stack([np.ones(n_ref, dtype=float), X_ref.astype(float)])

    coefs = np.zeros((n_feat, X.shape[1]), dtype=float)
    stds  = np.ones(n_feat, dtype=float)

    for i in range(n_feat):
        y = ref_mat[:, i].astype(float)
        if np.all(y == 0):
            continue
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res      = y - X.dot(coef)
        sigma    = np.std(res)
        coefs[i, :] = coef
        stds[i]     = sigma if sigma > 0 else 1.0

    expected = coefs[:, 0] + coefs[:, 1:].dot(X_pat.astype(float))
    return (pat_vec.astype(float) - expected) / stds

def main():
    # --- snakemake variables ---
    data_file       = str(snakemake.input.data_file)
    ref_mat_pkl     = str(snakemake.input.ref_mat)
    demo_ref_csv    = str(snakemake.input.demographics_ref)
    demo_pat_csv    = str(snakemake.input.demographics_pat)
    method          = snakemake.params.method
    normative_cols  = snakemake.params.normative_columns
    out_path        = str(snakemake.output[0])

    # 1) load pickled reference matrix
    with open(ref_mat_pkl, "rb") as f:
        ref_mat = pickle.load(f)
    # collapse any extra dims beyond (n_ref × n_vertices)
    if ref_mat.ndim > 2:
        ref_mat = ref_mat.mean(axis=tuple(range(2, ref_mat.ndim)))

    # 2) load patient data
    pat_vec = load_gifti(data_file)

    # 3) choose scoring method
    if method.lower().startswith("w"):
        # demographics → design
        df_ref = pd.read_csv(demo_ref_csv)
        df_pat = pd.read_csv(demo_pat_csv)
        if df_pat.shape[0] != 1:
            raise ValueError("Patient demographics CSV must contain exactly one row")

        # unify case
        df_ref.columns = df_ref.columns.str.upper()
        df_pat.columns = df_pat.columns.str.upper()
        normative_cols = [c.upper() for c in normative_cols]

        # simple map for SEX column (if present)
        if "SEX" in normative_cols:
            df_ref["SEX"] = df_ref["SEX"].map({"M": 1.0, "F": 0.0})
            df_pat["SEX"] = df_pat["SEX"].map({"M": 1.0, "F": 0.0})

        X_ref = df_ref[normative_cols].values
        X_pat = df_pat[normative_cols].values.flatten()
        scores = calculate_wscore(ref_mat, pat_vec, X_ref, X_pat)

    else:
        scores = calculate_zscore(ref_mat, pat_vec)

    # 4) save out .func.gii
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    da = GiftiDataArray(scores.astype(np.float32))
    nib.save(GiftiImage(darrays=[da]), out_path)
    print(f"Saved {method}-score map to {out_path}")

if __name__ == "__main__":
    main()
