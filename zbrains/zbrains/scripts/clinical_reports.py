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

def main():
    _res_dir = (_Path(__file__).resolve().parent.parent / "resources").resolve()
    snk = globals().get('snakemake')
    if not snk:
        raise RuntimeError('This script must be run via Snakemake')

    score_files    = snk.input.score_files
    demographics   = snk.input.demographics_pat
    report_params  = snk.params.report_params
    subject        = snk.params.subject
    session        = snk.params.session
    output_pdf     = Path(snk.output[0])

    res_dict  = {'cortex':'32k','hippocampus':'0p5mm'}

    df  = pd.read_csv(demographics) if demographics else pd.DataFrame()
    row = df[(df.ID == subject) & (df.SES == session)]
    age = row.AGE.values[0] if len(row) == 1 else None
    sex = row.SEX.values[0] if len(row) == 1 else None

    def apply_threshold(arr):
        return np.clip(arr, -report_params['threshold'], report_params['threshold'])

    def plot_cortex(data_l, data_r, out_png):
        res = res_dict['cortex']
        surf_l = read_surface(str(_res_dir / f'fsLR-{res}.L.inflated.surf.gii'), itype='gii')
        surf_r = read_surface(str(_res_dir / f'fsLR-{res}.R.inflated.surf.gii'), itype='gii')
        arr = apply_threshold(np.hstack([data_l, data_r]))
        plot_hemispheres(
            surf_l, surf_r,
            array_name=arr,
            cmap=report_params['cmap'],
            share='both',
            screenshot=True,
            interactive=False,
            filename=str(out_png),
            size=(800,250),
            color_bar=report_params['color_bar'],
        )
        return out_png

    def plot_subcortical(x, out_png):
        x = apply_threshold(x)
        x_pad = np.insert(x, len(x)//2, np.nan)  # left/right ventricle
        x_pad = np.append(x_pad, np.nan)
        sctx_l = read_surface(str(_res_dir / 'sctx.L.surf.gii'), itype='gii')
        sctx_r = read_surface(str(_res_dir / 'sctx.R.surf.gii'), itype='gii')
        n_vert = [867,1419,3012,3784,1446,4003,3726,7653,838,1457,3208,3742,1373,3871,3699,7180]
        parc_dat = np.repeat(x_pad, n_vert)
        plot_hemispheres(
            sctx_l, sctx_r,
            array_name=parc_dat,
            cmap=report_params['cmap'],
            share='both',
            screenshot=True,
            interactive=False,
            size=(500,200),
            color_bar=report_params['color_bar'],
            filename=str(out_png)
        )
        return out_png

    def plot_hippocampus(data_l, data_r, out_png):
        cdata = np.stack([data_l, data_r], axis=1)[..., np.newaxis]
        cdata = apply_threshold(cdata)
        hm.plotting.surfplot_canonical_foldunfold(
            cdata,
            hemis=['L', 'R'],
            labels=['hipp'],
            den=res_dict['hippocampus'],
            unfoldAPrescale=False,
            size=[200, 300],
            screenshot=True,
            filename=str(out_png),
            cmap=report_params['cmap'],
            color_bar=report_params['color_bar'],
        )
        return out_png

    banner_path = _res_dir / 'zbrains_banner.png'
    html = [
        f'<p><img src="{banner_path}" style="width:100%;"/></p>',
        f'<h1>Clinical Report: {subject} ({session}), Age: {age}, Sex: {sex}</h1>'
    ]
    features   = sorted({Path(f).stem.split('_feature-')[1].split('_')[0] for f in score_files})
    structures = ['cortex', 'subcortical', 'hippocampus']

    for feat in features:
        html.append(f'<h2>Feature: {feat}</h2>')
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
                # apply threshold
                png = plot_cortex(dl, dr, output_pdf.with_suffix(f'.cortex.{feat}.png'))

            elif struct == 'subcortical':
                dna = nib.load(fl[0]).darrays[0].data
                # apply threshold
                png = plot_subcortical(dna, output_pdf.with_suffix(f'.subcortical.{feat}.png'))

            else:
                lh = [x for x in fl if 'hemi-L' in x]
                rh = [x for x in fl if 'hemi-R' in x]
                dl = np.stack([nib.load(x).darrays[0].data for x in lh]).mean(0)
                dr = np.stack([nib.load(x).darrays[0].data for x in rh]).mean(0)
                # apply threshold
                png = plot_hippocampus(dl, dr, output_pdf.with_suffix(f'.hippocampus.{feat}.png'))

            html.append(
                f'<div style="text-align:center;">'
                  f'<img src="{png.name}" />'
                f'</div>'
            )
        html.append('<div style="page-break-after:always;"></div>')

    with open(output_pdf, 'w+b') as f:
        pisa.CreatePDF(''.join(html), dest=f)

    print(f'Written {output_pdf}')

    report_stem = output_pdf.stem
    for png_file in output_pdf.parent.glob(f"{report_stem}*.png"):
        try:
            png_file.unlink()
        except Exception as e:
            print(f"Could not remove temporary file {png_file}: {e}")

if __name__ == '__main__':
    main()
