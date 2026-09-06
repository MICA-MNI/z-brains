"""Run the greedy staged pipeline optimization pooled across MICs + NOEL.

This is the entry point for the sequential greedy search described in
``zbrains_staged``: start from the simplest pipeline (no intensity normalization
-> raw nativepro T1w/FLAIR, simplest blur, white-surface intensity sampling, plain
z-score, no demographics), then optimize one stage at a time -- smoothing, self
normalization, intensity sampling, blur, normative model -- keeping whichever
choice maximizes the per-subject macro AUROC for lesion detection, pooled across
both cohorts as a single pipeline.

Cost note: because smoothing changes the processed map filenames, the greedy
search may build separate processed bases for the visited smoothing x
normalization choices. Once a base exists, later candidates reuse it via cheap
re-analysis + AUROC scoring.

Run with the pipeline env (has skfmm/brainspace/ants):

    PYTHONPATH=src /host/verges/tank/data/ian/anaconda3/envs/newzb/bin/python example_staged.py

To split the 5 objectives across TWO machines sharing storage, run the second with
the objective order REVERSED (each resume-skips whatever the other already finished):

    # machine B
    ZBRAINS_REVERSE=1 PYTHONPATH=src .../newzb/bin/python example_staged.py
"""

import os

from zbrains.dataset import zbdataset, demographics
from zbrains.environment import zbenv

from zbrains_staged import (
    Cohort, run_staged_optimization, cross_evaluate_configs, plot_cross_evaluation,
)

# Keep the selective-reuse optimization separate from both the earlier v1 run and
# the first all-V2 attempt. The latter remains intact in ``results/hippunfold_v2``.
# Paths are absolute because the processing pipeline changes its working directory.
RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "results",
    "hippunfold_v2_selective",
)


# Lesion detection features. qT1 is MICs-only; NOEL omits it (handled per cohort).
MICS_FEATURES = ["T1w", "T1w*blur", "FA", "ADC", "thickness", "FLAIR", "FLAIR*blur", "qT1", "qT1*blur"]
NOEL_FEATURES = ["T1w", "T1w*blur", "FA", "ADC", "thickness", "FLAIR", "FLAIR*blur"]

env = zbenv(connectome_workbench_path="/usr/bin", num_threads=64, num_threads_wb=16)


def _build_cohort(
    name,
    score_name,
    output_dir_prefix,
    pipeline_dirs,
    features,
    control_csv,
    patient_csv,
    hippunfold_reuse_mode,
):
    control_demo = demographics(
        control_csv, normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"]
    )
    patient_demo = demographics(
        patient_csv, reference=control_demo,
        normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"],
    )
    common = dict(cortex=True, hippocampus=True, subcortical=True, **pipeline_dirs)
    control_dataset = zbdataset("controls", demographics=control_demo, **common)
    patient_dataset = zbdataset("patients", demographics=patient_demo, **common)
    return Cohort(
        name=name,
        score_name=score_name,
        control_dataset=control_dataset,
        patient_dataset=patient_dataset,
        output_dir_prefix=output_dir_prefix,
        base_pipeline_settings=dict(
            features=features, cortical_smoothing=5, hippocampal_smoothing=2
        ),
        hippunfold_reuse_mode=hippunfold_reuse_mode,
    )


mics = _build_cohort(
    name="MICs",
    score_name="MICs",
    output_dir_prefix="/host/verges/tank/data/BIDS_MICs/derivatives",
    pipeline_dirs={
        "micapipe_directory": "/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0",
        "hippunfold_directory": "/host/verges/tank/data/BIDS_MICs/derivatives/hippunfold_v2.0.0",
        "hippunfold_version": 2,
        "freesurfer_directory": "/data/mica3/BIDS_MICs/derivatives/freesurfer",
        "raw_data_directory": "/data/mica3/BIDS_MICs/rawdata",
    },
    features=MICS_FEATURES,
    control_csv="data/participants_mics_hc_all.csv",
    patient_csv="data/participants_mics_px_all.csv",
    # MICs was already processed with this exact V2 source: retain its original
    # base/output paths and reuse all existing processed maps.
    hippunfold_reuse_mode="all",
)

noel = _build_cohort(
    name="NOEL",
    score_name="NOEL",
    output_dir_prefix="/host/verges/tank/data/ian/BIDS_NOEL/derivatives",
    pipeline_dirs={
        "micapipe_directory": "/host/verges/tank/data/ian/BIDS_NOEL/derivatives/micapipe_v0.2.0",
        "hippunfold_directory": "/host/verges/tank/data/ian/BIDS_NOEL/derivatives/hippunfold_v2.0.0",
        "hippunfold_version": 2,
        "freesurfer_directory": "/host/verges/tank/data/ian/BIDS_NOEL/derivatives/freesurfer",
        "raw_data_directory": "/host/verges/tank/data/ian/BIDS_NOEL/rawdata",
    },
    features=NOEL_FEATURES,
    control_csv="data/participants_noel_hc.csv",
    patient_csv="data/participants_noel_all.csv",  # TLE (sub-PX) + FCD (sub-FCDPX) NOEL patients
    # Build V2-tagged bases, but seed cortex/subcortex/normalization from each
    # matching V1 base and regenerate only HippUnfold-dependent products.
    hippunfold_reuse_mode="nonhippocampal",
)


if __name__ == "__main__":
    # Run the greedy optimization FIVE times, one objective (comparison) after the
    # other, each into its own history so the locked pipelines can be compared.
    # The greedy may reach a DIFFERENT solution under each objective.
    #   1. within_tle       -- lesion |z| vs the patient's OWN non-lesion vertices, TLE only
    #   2. within_fcd       -- same, FCD only
    #   3. within           -- same, TLE and FCD together
    #   4. control          -- lesion |z| vs ALL held-out-control vertices
    #   5. disease_control  -- lesion |z| vs the OTHER disease's patients' non-lesional
    #                          vertices (TLE vs FCD, disease control)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    ANALYSES = [
        ("within_tle", os.path.join(RESULTS_DIR, "greedy_staged_history_within_tle.csv")),
        ("within_fcd", os.path.join(RESULTS_DIR, "greedy_staged_history_within_fcd.csv")),
        ("within", os.path.join(RESULTS_DIR, "greedy_staged_history_within_all.csv")),
        ("control", os.path.join(RESULTS_DIR, "greedy_staged_history_control.csv")),
        ("disease_control", os.path.join(RESULTS_DIR, "greedy_staged_history_disease_control.csv")),
    ]

    # Objective run ORDER -- forward (default) or reverse. To parallelize across TWO
    # machines, run one forward and one with ZBRAINS_REVERSE=1: they share the
    # processed bases and the per-objective history CSVs on shared storage, so each
    # machine resume-SKIPS (near-instantly, no re-processing) any objective the other
    # has already finished. Starting from opposite ends, the two split the heavy work
    # and meet in the middle; both still end with all 5 winners (computed or resumed)
    # so the cross-evaluation runs correctly on either.
    #   machine A:                       python example_staged.py
    #   machine B:  ZBRAINS_REVERSE=1    python example_staged.py
    if os.environ.get("ZBRAINS_REVERSE", "").strip().lower() in ("1", "true", "yes", "on"):
        ANALYSES = list(reversed(ANALYSES))
        print("[example_staged] ZBRAINS_REVERSE set -> running the 5 objectives in REVERSE order")
    else:
        print("[example_staged] running the 5 objectives in FORWARD order "
              "(set ZBRAINS_REVERSE=1 on a second machine to split the work)")

    # Reference-CV scope (see run_staged_optimization): "control_disease" (default)
    # reference-CVs ONLY the control + disease_control objectives (per-fold mean+/-SD);
    # the within_* objectives stay single-shot (cheap). Set "full" to reference-CV
    # every objective (~5x cost on within_*), or "off" for all single-shot.
    CV_SCOPE = "control_disease"

    results = {}
    for mode, history_path in ANALYSES:
        print(f"\n{'#' * 78}\n### GREEDY ANALYSIS: objective_mode = {mode!r}\n{'#' * 78}")
        best_config, _history = run_staged_optimization(
            cohorts=[mics, noel],
            env=env,
            objective_mode=mode,
            history_path=history_path,
            reprocess_controls=False,
            cv_scope=CV_SCOPE,
        )
        results[mode] = best_config

    print(f"\n{'#' * 78}\n### BEST PIPELINE PER OBJECTIVE\n{'#' * 78}")
    for mode, best_config in results.items():
        print(f"\nBest pipeline configuration ({mode}):")
        for key, value in best_config.items():
            print(f"  {key} = {value}")

    # The objectives optimize different things, so the greedy can lock DIFFERENT
    # pipelines -- surface each config key where the winners disagree.
    all_keys = sorted({k for cfg in results.values() for k in cfg})
    print(f"\n{'=' * 78}\nWHERE THE {len(results)} WINNERS DIFFER (blank = same as first):")
    for key in all_keys:
        vals = {mode: cfg.get(key) for mode, cfg in results.items()}
        if len(set(map(repr, vals.values()))) > 1:
            print(f"  {key}:")
            for mode in results:
                print(f"      {mode:16s} = {vals[mode]!r}")

    # CROSS-EVALUATION: score each of the 5 winning pipelines under ALL 5 metrics,
    # average (+/- SD) across metrics, and pick the most robust overall winner.
    print(f"\n{'#' * 78}\n### CROSS-EVALUATION (each winner tested on all 5 metrics)\n{'#' * 78}")
    matrix = cross_evaluate_configs(results, [mics, noel], env, verbose=True)
    matrix_csv = os.path.join(RESULTS_DIR, "greedy_cross_eval_matrix.csv")
    plot_png = os.path.join(RESULTS_DIR, "greedy_cross_eval.png")
    matrix.to_csv(matrix_csv)
    print("\nCross-evaluation AUROC matrix (rows = winner, cols = metric):")
    print(matrix.to_string(float_format=lambda v: f"{v:.4f}"))

    winner = plot_cross_evaluation(matrix, plot_png)
    print(f"\nWrote {matrix_csv} and {plot_png}")
    print(f"\n{'=' * 78}\nFINAL WINNER (highest mean AUROC across all 5 metrics): {winner!r}")
    print(f"  mean AUROC = {matrix.loc[winner, 'mean']:.4f} +/- {matrix.loc[winner, 'std']:.4f}")
    print("  pipeline configuration:")
    for key, value in results[winner].items():
        print(f"    {key} = {value}")
