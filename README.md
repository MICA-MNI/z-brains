# zbrains

Neuroimaging processing and analysis workflow (controls vs patients) with automated validation and clinical report generation.

## 1. Installation

```bash
# From repository root
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install .
```

Verify install:

```bash
python -c "import zbrains, sys; print('zbrains imported, Python', sys.version)"
```

If import fails, ensure the package layout is `src/zbrains/` (package name must not be `src`). Reinstall after fixing.

## 2. Required Inputs

The example script [test.py](test.py) expects the following existing derivative directories:

```
/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0
/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0
/data/mica3/BIDS_MICs/derivatives/freesurfer
```

Demographics CSV files (relative to repo root):
```
data/participants_mics_hc_all.csv
data/participants_mics_px_all.csv
```

Each CSV must include the normative columns:
- AGE (int)
- SEX (binary; expected encoding consistent between control & patient)

## 3. Quick Start (Run Example Pipeline)

```bash
source .venv/bin/activate
python test.py
```

What it does:
1. Imports core API:
   - `zbrains.dataset.zbdataset`
   - `zbrains.dataset.demographics`
   - `zbrains.environment.zbenv`
2. Defines feature list:
   ```
   ["FA", "ADC", "thickness", "qT1", "qT1*Blur", "FLAIR", "FLAIR*Blur"]
   ```
3. Builds an execution environment (threads + external tool path).
4. Loads control demographics, constructs control dataset, processes & validates.
5. Loads patient demographics (referencing control for normalization), processes, validates.
6. Runs comparative analysis (`method='wscore'`).
7. Generates patient clinical report.
8. Prints `end` when complete.

Outputs written to:
```
/host/verges/tank/data/ian/zbrains_outputs
```

Ensure this path is writable.

## 4. Customizing

Edit [test.py](test.py):
- Change `features` list to include/remove modalities available in your derivatives.
- Adjust smoothing (`cortical_smoothing`, `hippocampal_smoothing`).
- Point the *directory* arguments to your own derivative locations.
- Switch analysis method in `patient_dataset.analyze(..., method='wscore')` if other methods implemented.

## 5. Headless / Cluster Execution

Set thread counts in `zbenv` to match available cores. For batch runs, wrap invocation:

```bash
python -u test.py > run.log 2>&1 &
tail -f run.log
```

## 6. Expected Directory Structure (Minimal)

```
z-brains-v2.0/
├─ src/
│  └─ zbrains/            (package modules)
├─ data/
│  ├─ participants_mics_hc_all.csv
│  └─ participants_mics_px_all.csv
├─ test.py
└─ README.md
```

## 7. Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| `ModuleNotFoundError: No module named 'zbrains'` | Package not installed or wrong layout | Confirm `src/zbrains/__init__.py`, then `pip install --force-reinstall .` |
| `FileNotFoundError` for derivatives | Paths in script not valid | Update paths in [test.py](test.py) to your environment |
| Demographics normalization errors | Missing / mismatched columns | Ensure both CSVs contain AGE, SEX with consistent dtypes |
| Slow run | Too many smoothing/threads mismatch | Lower smoothing kernels or adjust thread counts in `zbenv` |

Reinstall cleanly if needed:

```bash
pip uninstall -y zbrains
pip install .
```

## 8. Reproducibility Notes

- Record git commit: `git rev-parse HEAD`
- Capture environment: `pip freeze > requirements.lock.txt`
- Store exact command line used.

## 9. Support

Open an issue with:
- Command executed
- Full traceback
- git commit hash
- Snippet from demographics CSV header
- Listing of output directory

## 10. License

BSD-3-Clause (see included license metadata).
