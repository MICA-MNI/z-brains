# Z-Brains v2.0.0

**Python workflow for normative neuroimaging analysis**

Z-Brains is a Python toolkit for analysing structural and diffusion MRI data, building normative control cohorts, and comparing patient scans through weighted scoring. The pipeline is meant to be scripted directly in Python so you can tailor data paths, covariates, and processing steps in a single place.

---

## Current z-brains directory for MICA Lab:

"/host/verges/tank/data/BIDS_MICs/derivatives/zbrains_2.0/"

## Install

```bash
git clone <repository-url>
cd z-brains
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install .
```

Verify the package is importable:

```bash
python -c "import zbrains; print('Z-Brains ready')"
```

---

## Quick start (5 minutes)

1. **Activate your environment**

   ```bash
   source .venv/bin/activate
   ```

2. **Open the template script**  
   Edit [`example.py`](example.py) and set:
   - `OUTPUT_DIR` to where results should be written.
   - `PIPELINE_DIRS` to the Micapipe, HippUnfold, and FreeSurfer derivative roots.
   - `FEATURES`, `CORTEX`, `HIPPOCAMPUS`, `SUBCORTICAL`, or `reprocess` if needed.

3. **Run the end-to-end pipeline**

   ```bash
   python example.py
   ```

This script handles environment creation, demographics loading, control validation, patient processing, statistical comparison, and report generation.

---

## Adapting the workflow

Key pieces in [`example.py`](example.py):

- [`zbrains.environment.zbenv`](example.py) configures threads and Connectome Workbench.
- [`zbrains.dataset.demographics`](example.py) loads CSVs and declares normative covariates.
- [`zbrains.dataset.zbdataset`](example.py) encapsulates control or patient cohorts.

Typical customizations:

```python
PIPELINE_SETTINGS = dict(features=["FA", "thickness"], cortical_smoothing=5, hippocampal_smoothing=3)
patient_dataset.analyze(output_directory=OUTPUT_DIR, reference=control_dataset, method="wscore")
```

W-score modeling can be configured independently at several stages:

```python
patient_dataset.analyze(
    output_directory=OUTPUT_DIR,
    reference=control_dataset,
    method="wscore",
    wscore_distribution="gaussian_mad",  # also gaussian_winsor, student_t, empirical, shash
    wscore_preprocessing="none",          # or spatial_zscore / spatial_robust_z
    wscore_covariate_model="quadratic_age",
    wscore_surface_smoothing_iterations=16,
    blur_depth_model="mean_slope_rms",   # use "mean" for legacy depth averaging
    intensity_depth_model="raw",         # or "white_swm1_direction_cosine"
)
```

`spatial_zscore` normalizes each subject map by its spatial mean/SD, while
`spatial_robust_z` uses its median/MAD. Both are applied after image-intensity
normalization and before regression. Surface smoothing is applied only to fsLR-32k cortical score maps;
other surface densities and subcortical CSVs are left unsmoothed. Keep a held-out
patient set when selecting among these options.

For multi-depth blur features, `mean_slope_rms` fits separate normative models
to the depth mean and equal-spaced depth slope, then combines their W-scores as
a signed RMS magnitude. The output JSON records the selected depth model.

`white_swm1_direction_cosine` derives the fsLR-32k white-surface T1w/FLAIR
analysis from the companion four-depth map as `(white - SWM1) / depth_energy`,
with a within-map denominator floor for stability. Include the corresponding
`T1w*blur` or `FLAIR*blur` feature whenever this toggle is enabled.

`multisurface_median_abs_dominant` scales all four surfaces in both hemispheres
by one within-subject median absolute intensity, fits a signed W-score at every
depth, and returns the depth with the largest absolute abnormality. Its maximum
statistic is control-calibrated back to a 1.96 threshold, and a companion
`.shape.gii` records the selected depth at every vertex.


Use `reprocess = True` to regenerate control derivatives after upstream updates.

---

## T1w/FLAIR normalization

When `T1w`, `T1w*blur`, `FLAIR`, or `FLAIR*blur` are requested, Z-Brains normalizes those scans before any volume-to-surface extraction:

- Each subject is WhiteStripe-normalized in Python using SynthSeg white matter labels (`2`, `41`) as the normal-appearing white matter estimate.
- By default, control datasets are normalized first, then a Python RAVEL model is fitted from the WhiteStripe-normalized controls and saved under `${output_directory}/ravel/`.
- SynthSeg CSF labels are used as RAVEL control voxels.
- Patient datasets run WhiteStripe first, then apply the already-fitted control RAVEL model.
- Set `normalization="whitestripe"` in `zbdataset.process(...)` to skip RAVEL and use `maps/{pid}_{sid}_space-nativepro_desc-whitestripe_{T1w,FLAIR}.nii.gz` for extraction. The same setting can be kept in a shared settings dictionary passed to `validate(...)`.
- With the default `normalization="ravel"`, surface, hippocampal, blur, and subcortical extraction use `maps/{pid}_{sid}_space-nativepro_desc-ravel_{T1w,FLAIR}.nii.gz`.

This replaces the previous T1w/FLAIR and FLAIR/T1w ratio-map normalization path.

---

## Expected inputs

Ian-Goodall-Halliwell:  can you analyze #file:dataset.py and #file:processing.py and #file:analysis.py and other files in #codebase and figure out exactly what files are needed to run the processing step, for each feature and each structure. specify the exact file paths

GitHub Copilot: Searched codebase for "can you analyze #file:dataset.py and #file:processing.py and #file:analysis.py and other files in #codebase and figure out exactly what files are needed to run the processing step, for each feature and each structure. specify the exact file paths", 22 results

Here is the exact file inventory required for `zbdataset.process()` derived from dataset.py, processing.py, and analysis.py.

### 1. Structure-wide prerequisites

| Structure | Location root | Required files |
| --- | --- | --- |
| Cortex | `${micapipe_directory}/{participant_id}/{session_id}` | `anat/{pid}_{sid}_space-nativepro_T1w.nii.gz`, `surf/{pid}_{sid}_hemi-{L,R}_surf-fsnative_label-sphere.surf.gii`, `surf/{pid}_{sid}_hemi-{L,R}_space-nativepro_surf-fsnative_label-{pial,midthickness,white}.surf.gii`, `surf/{pid}_{sid}_hemi-{L,R}_space-nativepro_surf-fsLR-32k_label-{midthickness,pial}.surf.gii`, plus FreeSurfer exports `freesurfer_directory/{pid}_{sid}/surf/{lh,rh}.white` and `.../label/{lh,rh}.cortex.label` for medial-wall derivation. |
| Hippocampus | `${hippunfold_directory}/hippunfold/{participant_id}/{session_id}` | `surf/{pid}_{sid}_hemi-{L,R}_space-T1w_den-0p5mm_label-hipp_{midthickness,inner,outer}.surf.gii` and `surf/{pid}_{sid}_hemi-{L,R}_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii`. |
| Subcortical | `${freesurfer_directory}/{participant_id}_{session_id}` + micapipe | `stats/aseg.stats`, `parc/{pid}_{sid}_space-nativepro_T1w_atlas-subcortical.nii.gz`. |

Blur processing additionally needs:
- FreeSurfer parcellation `mri/aparc+aseg.mgz`.
- Native-pro brain volume `anat/{pid}_{sid}_space-nativepro_T1w_brain.nii.gz`.
- Surfaces `surf/{pid}_{sid}_hemi-{L,R}_space-nativepro_surf-fsnative_label-{white,midthickness}.surf.gii`.
- Output Laplace file `structural/{pid}_{sid}-laplace.nii.gz` and generated shift surfaces `structural/{pid}_{sid}_{hemi}_sfwm-{1.0,2.0}mm.surf.gii`.

### 2. Feature-specific inputs (cortex + blur)

| Feature token | Micapipe inputs required (per hemi + label) |
| --- | --- |
| `FA` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-{midthickness,white}_FA.func.gii`, `maps/{pid}_{sid}_space-nativepro_model-DTI_map-FA.nii.gz`. |
| `ADC` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-{midthickness,white}_ADC.func.gii`, `maps/{pid}_{sid}_space-nativepro_model-DTI_map-ADC.nii.gz`. |
| `thickness` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-thickness.func.gii`. |
| `SA` | `surf/{pid}_{sid}_hemi-{L,R}_space-nativepro_surf-fsnative_label-{midthickness,white}.surf.gii`. |
| `T1w` | `anat/{pid}_{sid}_space-nativepro_T1w.nii.gz`; normalized derivatives are written to `maps/{pid}_{sid}_space-nativepro_desc-{whitestripe,ravel}_T1w.nii.gz`. |
| `FLAIR` | `maps/{pid}_{sid}_space-nativepro_map-flair.nii.gz`; normalized derivatives are written to `maps/{pid}_{sid}_space-nativepro_desc-{whitestripe,ravel}_FLAIR.nii.gz`. |
| `qT1` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-{midthickness,white}_T1map.func.gii`, `maps/{pid}_{sid}_space-nativepro_map-T1map.nii.gz`. |
| `FLAIR*blur` | Same FLAIR volumetric input above; blur outputs stored as `maps/cortex/{pid}_{sid}_hemi-{hemi}_feature-FLAIR*blur_*`. |
| `qT1*blur` | Same qT1 inputs; outputs `...feature-qT1*blur_*`. |
| `fMRI` | `func/desc-se_task-rest_acq-AP_bold/surf/{pid}_{sid}_surf-fsLR-32k_desc-timeseries_clean.shape.gii`. Generates `rmssd`, `timescales`, `alff`, `falff`. |

### 3. Hippocampal feature inputs

For every non-blur feature allowed in hippocampus (`thickness`, `FLAIR`, `ADC`, `FA`, `qT1`):
- Volume map from micapipe (`space-nativepro_map-*` or `model-DTI_map-*`) except `thickness`, which uses hippunfold `*_label-hipp_thickness.shape.gii`.
- Hippunfold surfaces listed above for `wb_command -ribbon-constrained`.

### 4. Subcortical feature inputs

- Base `feature-volume`: FreeSurfer `aseg.stats` only.
- Other features (`FLAIR`, `ADC`, `FA`, `qT1`): micapipe volumes (`space-nativepro_map-*` or `model-DTI_map-*`) plus segmentation `parc/...atlas-subcortical.nii.gz` to extract region means.

### 5. Outputs/checks consumed later

Validation expects, per subject, directories under `${output_directory}/{pid}/{sid}`:
- `structural/` copies of all surfaces, Laplace file, medial wall labels.
- `maps/cortex`, `maps/hippocampus`, `maps/subcortical` populated with the smoothed feature files following naming patterns described above (fsLR resolutions 32k & 5k, labels `midthickness` & `white`, hippocampal `den-0p5mm`, subcortical CSVs).
- For blur: raw/dist/grad files `..._feature-{feat}*blur_surf-fsnative_desc-{raw,dist,grad}.func.gii` and smoothed `..._surf-fsnative_smooth-{cortical_smoothing}mm.func.gii`.
- For fMRI: `maps/cortex` populated with `rmssd`, `timescales`, `alff`, `falff` metric maps at 32k and 5k resolutions.

This inventory comes directly from the file checks in `zbdataset.add_features()` and `validate()` plus the runtime needs of `apply_blurring`, `apply_cortical_processing`, `apply_hippocampal_processing`, and `apply_subcortical_processing`.

## Outputs

`OUTPUT_DIR/` will contain:

```
├── <sub>/<ses>/maps/...
├── <sub>/<ses>/analysis/...
├── <sub>/<ses>/structural/...
└── logs/
```

Patient reports summarize cortical, hippocampal, and subcortical deviations relative to controls.

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| Missing module | `pip install --force-reinstall .` |
| Workbench not found | Install it or point `connectome_workbench_path` to the binary directory. |
| CSV mismatch | Ensure both demographics files share normative columns and types. |
| Performance/memory constraints | Reduce `features`, smoothing kernels, or thread counts in the script. |

---

## Need help?

Collect the executed command (`python example.py`), stack trace, `git rev-parse HEAD`, demographics headers, and derivative directory listings before opening an issue.
