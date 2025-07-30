#!/usr/bin/env python3
"""
clinical_reports.py: Generate clinical PDF reports from score .func.gii maps.

Snakemake integration:
  input:
    - score_files: list of score .func.gii paths (patient and reference)
    - demographics (optional): CSV with columns ['ID','age','sex']
  params:
    - report_params: dict with keys {
        analyses (list of str),
        threshold (float),
        threshold_alpha (float),
        label_ctx (str),
        label_hip (str),
        color_bar (str),
        cmap (str),
        resolution (dict),
        smoothings (dict)
      }
    - subject: str
    - session: str or None
    - approach: 'zscore' or 'wscore'
    - output_dir: path to write PDF
  output:
    - report: clinical_report.pdf
"""
import os
import itertools
import uuid
import yaml
import pandas as pd
import numpy as np
import nibabel as nib
from xhtml2pdf import pisa
from pathlib import Path
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from scipy.spatial.transform import Rotation

# Helper imports for snakemake
snk = globals().get('snakemake', None)
if snk is None:
    raise RuntimeError('This script must be run via Snakemake')

# ----- Snakemake variables -----
score_files = snk.input.score_files
demographics_csv = snk.input.get('demographics', None)
params = snk.params
report_params = params.report_params
subject = params.get('subject')
session = params.get('session', None)
approach = params.get('approach', 'zscore')
output_pdf = snk.output[0]

# Config params
analyses = report_params.get('analyses', ['regional','asymmetry'])
threshold = report_params.get('threshold', 1.96)
threshold_alpha = report_params.get('threshold_alpha', 0.3)
label_ctx = report_params.get('label_ctx', 'midthickness')
label_hip = report_params.get('label_hip', 'midthickness')
color_bar = report_params.get('color_bar', 'bottom')
cmap = report_params.get('cmap', 'cmo.balance')
resolution = report_params.get('resolution', {})
smoothings = report_params.get('smoothings', {})

# Read demographics if provided
age = sex = None
if demographics_csv and subject:
    df_demo = pd.read_csv(demographics_csv)
    row = df_demo[df_demo['ID'] == subject]
    if len(row) == 1:
        age = row['age'].iloc[0] if 'age' in row else None
        sex = row['sex'].iloc[0] if 'sex' in row else None

# HTML report builder
html = []
# Header
html.append(f'<h1>Clinical Report: Subject {subject}</h1>')
html.append(f'<p>Session: {session or "N/A"} | Age: {age or "n/a"} | Sex: {sex or "n/a"}</p>')

# Iterate analyses and features
# Derive features and structures from score_files naming
features = sorted(set([Path(f).stem.split('_feature-')[1].split('_')[0] for f in score_files]))
structures = sorted(set([Path(f).stem.split('_structure-')[1].split('_')[0] for f in score_files]))

for analysis, feature in itertools.product(analyses, features):
    html.append(f'<h2>{analysis.capitalize()} Analysis - Feature: {feature}</h2>')
    for struct in structures:
        html.append(f'<h3>{struct.capitalize()}</h3>')
                # Collect relevant score maps (matching feature, structure, and approach)
        files = [
            f for f in score_files
            if f'feature-{feature}' in f
            and f'structure-{struct}' in f
            and f'_score-{approach}' in f
        ]
        if not files:
            html.append(f'<p style="color:orange;">No score maps for struct={struct}, feature={feature}, approach={approach}</p>')
            continue
        # Load and average arrays
        mats = []
        for f in files:
            img = nib.load(f)
            arrs = [da.data for da in img.darrays]
            mat = np.stack(arrs,axis=0).mean(axis=0) if len(arrs)>1 else arrs[0]
            mats.append(mat)
        data = np.stack(mats,axis=0).mean(axis=0)
        # Threshold
        data[np.abs(data) < threshold] *= threshold_alpha
        # Plot maps
        # For cortex: load standard surfaces
        if struct == 'cortex':
            lh = read_surface(os.path.join(os.environ['DATA_PATH'], f'fsLR.{label_ctx}.lh.surf.gii'))
            rh = read_surface(os.path.join(os.environ['DATA_PATH'], f'fsLR.{label_ctx}.rh.surf.gii'))
            # plot and save PNG
            out_png = Path(output_pdf).with_suffix(f'.{struct}.{analysis}.{feature}.png')
            plot_hemispheres(
                lh, rh, array_name=data,
                cmap=cmap, color_bar=color_bar,
                filename=str(out_png)
            )
            html.append(f'<img src="{out_png.name}" width="600">')
        else:
            # simple table for subcortical
            # treat data as vector
            html.append('<table><tr><th>Region</th><th>Score</th></tr>')
            labels = SUBCORTICAL_STRUCTURES if struct=='subcortical' else []
            for lbl,val in zip(labels, data):
                html.append(f'<tr><td>{lbl}</td><td>{val:.2f}</td></tr>')
            html.append('</table>')

# Convert HTML to PDF
html_str = ''.join(html)
with open(output_pdf, 'w+b') as res:
    pisa.CreatePDF(html_str, dest=res)
print(f'Clinical report written to {output_pdf}')
