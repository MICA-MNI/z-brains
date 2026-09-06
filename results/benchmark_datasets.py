#!/usr/bin/env python3
"""Shared dataset definitions for the lesion-detection feature benchmarks.

Each benchmark estimates normative parameters from a cohort's controls and
evaluates lesion separation on that cohort's patients. To benchmark across both
the MICs and NOEL cohorts (each contributing both FCD and TLE lesions), every
benchmark iterates over ``DATASETS``: it computes per-dataset results and also
pools the raw lesion-vs-outside vertices across datasets for a combined score.

The processed-map roots below are the ``zbrains_WB`` derivatives that hold the
``maps/cortex`` GIFTIs the benchmarks read. Lesion ground truth lives under
``results/lesiondata/<COHORT>`` (FCD CSVs and TLE / hippocampal cavity GIFTIs).
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

DATASETS = [
    {
        "name": "MICs",
        "derivatives": Path("/host/verges/tank/data/BIDS_MICs/derivatives/zbrains_WB"),
        "micapipe": Path("/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0"),
        "control_csv": REPO_ROOT / "data/participants_mics_hc_all.csv",
        "patient_csv": REPO_ROOT / "data/participants_mics_px_all_test.csv",
        "lesion_dir": REPO_ROOT / "results/lesiondata/MICs",
        "hipp_density": "8k",
    },
    {
        "name": "NOEL",
        "derivatives": Path("/host/verges/tank/data/ian/BIDS_NOEL/derivatives/zbrains_WB"),
        "micapipe": Path("/host/verges/tank/data/ian/BIDS_NOEL/derivatives/micapipe_v0.2.0"),
        "control_csv": REPO_ROOT / "data/participants_noel_hc.csv",
        # TLE (NOEL-TLE_raw.xlsx, sub-PX###) + FCD (NOEL-FCD_raw.xlsx,
        # sub-FCDPX###) NOEL patients: 20 TLE + 45 FCD = 65 cases.
        "patient_csv": REPO_ROOT / "data/participants_noel_all.csv",
        "lesion_dir": REPO_ROOT / "results/lesiondata/NOEL",
        # The optimization is now locked to HippUnfold v2 surfaces.
        "hipp_density": "8k",
    },
]

POOLED_NAME = "POOLED"
