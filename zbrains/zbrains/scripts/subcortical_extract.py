#!/usr/bin/env python3
"""
Script: subcortical_extract.py
Extracts subcortical metrics (volume, thickness, FA, ADC, T1map) from FreeSurfer segmentation and scalar maps,
then saves the results as a GIFTI (.func.gii) file, compatible with the `subcortical_import` Snakemake rule.
"""
import os
import numpy as np
import nibabel as nib
from nibabel.gifti import GiftiImage, GiftiDataArray
import pandas as pd

# Constants defining subcortical labels and corresponding output names
SUBCORTICAL_LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
SUBCORTICAL_STRUCTURES = [
    "Laccumb", "Lamyg", "Lcaud", "Lhippo", "Lpal", "Lput", "Lthal",
    "Raccumb", "Ramyg", "Rcaud", "Rhippo", "Rpal", "Rput", "Rthal",
]


def extract_stats(stats_file, labels, feature):
    """
    Parse FreeSurfer aseg.stats file to extract subcortical volume or thickness for given labels.
    If feature is 'thickness', uses the volume column as a proxy.
    """
    # Read headers line
    with open(stats_file) as fh:
        for line in fh:
            if line.startswith('# ColHeaders'):
                headers = line.strip().split()[2:]
                break
        else:
            raise ValueError("Could not find 'ColHeaders' line in stats file")
    # Read data lines into records
    records = []
    with open(stats_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            tokens = line.strip().split()
            if len(tokens) < len(headers):
                continue
            records.append(dict(zip(headers, tokens)))
    df = pd.DataFrame(records)
    # Identify ID column
    id_col = next((col for col in ['SegId', 'StructId', 'Label'] if col in df.columns), None)
    if id_col is None:
        raise ValueError("No ID column found in stats file")
    df[id_col] = pd.to_numeric(df[id_col], errors='coerce')
    # Determine metric column: use volume for thickness
    key = 'volume' if feature.lower() == 'thickness' else feature.lower()
    candidates = [h for h in headers if key in h.lower()]
    if not candidates:
        raise ValueError(f"No stats column matching '{key}' in headers: {headers}")
    metric_col = candidates[0]
    df[metric_col] = pd.to_numeric(df[metric_col], errors='coerce')
    # Map each label to its metric value
    values = []
    for lbl in labels:
        row = df[df[id_col] == lbl]
        values.append(float(row.iloc[0][metric_col]) if not row.empty else np.nan)
    return values


def main():
    # Snakemake-provided inputs, params, and outputs
    seg_file   = snakemake.input.seg_file
    stats_file = snakemake.input.aseg_stats
    feature    = snakemake.params.feature
    volumemap  = snakemake.params.volumemap

    # Load segmentation
    seg_img = nib.load(seg_file)
    seg = np.asarray(seg_img.dataobj)

    # Extract metrics
    if feature in ('volume', 'thickness'):
        if stats_file is None:
            raise ValueError(f"Stats file is required for feature '{feature}' but got None")
        values = extract_stats(stats_file, SUBCORTICAL_LABELS, feature)
    else:
        if volumemap is None:
            raise ValueError(f"Volumetric map required for feature '{feature}' but got None")
        img = nib.load(volumemap)
        data = np.asarray(img.dataobj)
        values = []
        for lbl in SUBCORTICAL_LABELS:
            mask = seg == lbl
            values.append(np.nanmean(data[mask]) if mask.any() else np.nan)

    # Create output directory if needed
    out_path = snakemake.output[0]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    # Create GIFTI (.func.gii)
    da = GiftiDataArray(np.array(values, dtype=np.float32))
    gifti_img = GiftiImage(darrays=[da])
    nib.save(gifti_img, out_path)


if __name__ == '__main__':
    main()
