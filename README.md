# Z-Brains v2.0.0

**Python workflow for normative neuroimaging analysis**

Z-Brains is a Python toolkit for analysing structural and diffusion MRI data, building normative control cohorts, and comparing patient scans through weighted scoring. The pipeline is meant to be scripted directly in Python so you can tailor data paths, covariates, and processing steps in a single place.

---

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

Use `reprocess = True` to regenerate control derivatives after upstream updates.

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
| `FLAIR` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-{midthickness,white}_flair.func.gii`, `maps/{pid}_{sid}_space-nativepro_map-flair.nii.gz`. |
| `qT1` | `maps/{pid}_{sid}_hemi-{hemi}_surf-fsLR-32k_label-{midthickness,white}_T1map.func.gii`, `maps/{pid}_{sid}_space-nativepro_map-T1map.nii.gz`. |
| `FLAIR*blur` | Same FLAIR volumetric/surface files above; blur outputs stored as `maps/cortex/{pid}_{sid}_hemi-{hemi}_feature-FLAIR*blur_*`. |
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