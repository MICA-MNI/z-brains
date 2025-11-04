# Z-Brains v2.0.0

**Python workflow for normative neuroimaging analysis**

Z-Brains is a Python toolkit for analysing structural and diffusion MRI data, building normative control cohorts, and comparing patient scans through weighted scoring. The pipeline is meant to be scripted directly in Python so you can tailor data paths, covariates, and processing steps in a single place.

---

## Install

```bash
git clone <repository-url>
cd z-brains-v2.0
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

| Item | Description |
|------|-------------|
| `data/participants_mics_hc_all.csv` | Control demographics with columns like `participant_id`, `session_id`, `AGE`, `SEX`. |
| `data/participants_mics_px_all.csv` | Patient demographics; extra covariates should mirror the control sheet. |
| Micapipe / HippUnfold / FreeSurfer directories | BIDS-derivative folders that match the subject and session labels in the CSV files. |

---

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