#!/usr/bin/env python3
"""
generate_reports.py: Generate clinical PDF reports using BrainSpace and HippoMaps.
"""
import os
import itertools
import pandas as pd
import numpy as np
import nibabel as nib
from xhtml2pdf import pisa
from pathlib import Path, Path as _Path
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
import hippomaps as hm
import logging

# Enable logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Resource directory for surfaces
_res_dir = (_Path(__file__).resolve().parent.parent / "resources").resolve()

# Main entry
snk = globals().get('snakemake')
if not snk:
    raise RuntimeError('This script must be run via Snakemake')

# Snakemake inputs & parameters
score_files   = snk.input.score_files
demographics  = snk.input.demographics_pat
report_params = snk.params.report_params
subject       = snk.params.subject
session       = snk.params.session
output_pdf    = Path(snk.output[0])

# Report settings
analyses        = report_params.get('analyses', ['regional','asymmetry'])
threshold       = report_params.get('threshold', 1.96)
threshold_alpha = report_params.get('threshold_alpha', 0.3)
cmap            = report_params.get('cmap', 'cmo.balance')
color_bar       = report_params.get('color_bar', 'bottom')
res_dict        = report_params.get('resolution', {'cortex':'32k','hippocampus':'0p5mm'})

# Read demographics
df = pd.read_csv(demographics) if demographics else pd.DataFrame()
row = df[(df.ID==subject)&(df.SES==session)]
age = row.AGE.values[0] if len(row)==1 else None
sex = row.SEX.values[0] if len(row)==1 else None

# Plot cortex data
def plot_cortex(data_l, data_r, out_png):
    logger.info(f"Generating cortex plot: {out_png}")
    res = res_dict.get('cortex', '32k')
    surf_l = read_surface(str(_res_dir/f'fsLR-{res}.L.inflated.surf.gii'), itype='gii')
    surf_r = read_surface(str(_res_dir/f'fsLR-{res}.R.inflated.surf.gii'), itype='gii')
    arr = np.hstack([data_l, data_r])
    plot_hemispheres(
        surf_l, surf_r,
        array_name=arr,
        cmap=cmap,
        color_bar=color_bar,
        share='both',
        screenshot=True,
        interactive=False,
        filename=str(out_png)
    )
    return out_png

# Placeholder for subcortical
def plot_subcortical_placeholder(out_png):
    logger.warning('Subcortical plotting not supported; placeholder used')
    return out_png

# Plot hippocampus using HippoMaps surfplot_canonical_foldunfold
def plot_hippocampus(data_l, data_r, out_png):
    logger.info(f"Generating hippocampus plot: {out_png}")
    # build cdata array of shape V x 2 x 1
    cdata = np.stack([data_l, data_r], axis=1)[..., np.newaxis]
    # plot folded and unfolded hippocampus & dentate
    hm.plotting.surfplot_canonical_foldunfold(
        cdata,
        hemis=['L', 'R'],
        labels=['hipp'],
        unfoldAPrescale=False,
        tighten_cwindow=False,
        resourcesdir=str(_res_dir),
        size=[350, 300],
        screenshot=True,
        filename=str(out_png)
    )
    return out_png

# Build HTML
html = [f'<h1>Clinical Report: {subject} ({session}), Age: {age}, Sex: {sex}</h1>']
features = sorted({Path(f).stem.split('_feature-')[1].split('_')[0] for f in score_files})
structures = ['cortex','subcortical','hippocampus']

# Loop analyses/features
for ana, feat in itertools.product(analyses, features):
    html.append(f'<h2>{ana.capitalize()} - Feature: {feat}</h2>')
    for struct in structures:
        html.append(f'<h3>{struct.capitalize()}</h3>')
        fl = [f for f in score_files if f'feature-{feat}' in f and f'structure-{struct}' in f]
        if not fl:
            html.append('<p style="color:orange;">Missing data</p>')
            continue
        if struct == 'cortex':
            lh = [x for x in fl if 'hemi-L' in x]
            rh = [x for x in fl if 'hemi-R' in x]
            dl = np.stack([nib.load(x).darrays[0].data for x in lh]).mean(0)
            dr = np.stack([nib.load(x).darrays[0].data for x in rh]).mean(0)
            png = plot_cortex(dl, dr, output_pdf.with_suffix(f'.cortex.{ana}.{feat}.png'))
            html.append(f'<p><img src="{png.name}"/></p>')
        elif struct == 'subcortical':
            png = plot_subcortical_placeholder(output_pdf.with_suffix(f'.subcortical.{ana}.{feat}.png'))
            html.append(f'<p style="color:gray;">Subcortical plot unavailable</p>')
        else:
            lh = [x for x in fl if 'hemi-L' in x]
            rh = [x for x in fl if 'hemi-R' in x]
            dl = np.stack([nib.load(x).darrays[0].data for x in lh]).mean(0)
            dr = np.stack([nib.load(x).darrays[0].data for x in rh]).mean(0)
            png = plot_hippocampus(dl, dr, output_pdf.with_suffix(f'.hippocampus.{ana}.{feat}.png'))
            html.append(f'<p><img src="{png.name}"/></p>')

# Write PDF
with open(output_pdf, 'w+b') as f:
    pisa.CreatePDF(''.join(html), dest=f)

print(f'Written {output_pdf}')
