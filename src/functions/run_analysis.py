import sys
import logging
import argparse
import ast
import numpy as np

from .surface_to_volume import surface_to_volume
from .constants import (
    LIST_ANALYSES,
    Resolution,
)

from .utils_analysis import (
    get_id,
    get_session,
    run_analysis,
    get_bids_id,
    load_demo,
)
from .clinical_reports import generate_clinical_report


################################################################################
# Analysis
################################################################################
def main(
    dataset,
    zbrains_ref,
    demo_ref,
    column_map,
    subject_id,
    session,
    demo,
    zbrains,
    struct,
    feat,
    normative,
    deconfound,
    smooth_ctx,
    smooth_hip,
    threshold,
    approach,
    resolution,
    labels_ctx,
    labels_hip,
    tmp,
    logger,
    micapipename,
    hippunfoldname,
    n_jobs,
    n_jobs_wb,
    workbench_path,
    dicoms,
):
    # Some checks
    # logger = logging.getLogger(tmp)
    if len(zbrains_ref) != len(demo_ref):
        raise ValueError(
            "The number of values provided with --zbrains_ref "
            "and --demo_ref must be the same."
        )

    # Handle columns names -----------------------------------------------------
    # Column mapping
    expected_to_actual = dict(
        participant_id="participant_id",
        session_id="session_id",
        age="age",
        sex="sex",
        site="site",
    )
    unknown_cols = set(column_map.keys()).difference(expected_to_actual.keys())
    if unknown_cols:
        raise ValueError(
            f"Unknown column names: {unknown_cols}. Allowed options "
            f"are: {list(expected_to_actual.keys())}"
        )

    expected_to_actual.update(column_map)
    actual_to_expected = {v: k for k, v in expected_to_actual.items()}

    # Column types
    col_dtypes = {
        actual: (str if expected != "age" else float)
        for actual, expected in actual_to_expected.items()
    }

    # Rename covariates for normative modeling
    cov_normative = normative
    if cov_normative is not None:
        cov_normative = [actual_to_expected.get(col, col) for col in cov_normative]

    # Rename covariates for deconfounding
    cov_deconfound = deconfound
    if cov_deconfound is not None:
        for i, col in enumerate(cov_deconfound):
            prefix = ""
            if col.startswith("-"):
                prefix = "-"
                col = col[1:]
            cov_deconfound[i] = prefix + actual_to_expected.get(col, col)

    # Identifiers --------------------------------------------------------------
    px_id = get_id(subject_id, add_prefix=True)
    px_ses = get_session(session, add_predix=True)
    bids_id = get_bids_id(px_id, px_ses)

    # subject_dir = get_subject_dir(zbrains_path, sid, ses)
    # path_analysis = f'{subject_dir}/{approach_folder}'

    # Load px dataframe --------------------------------------------------------
    px_demo = None
    if demo is not None:
        px_demo = load_demo(demo, rename=actual_to_expected, dtypes=col_dtypes, tmp=tmp)
        px_demo = px_demo.loc[
            (px_demo["participant_id"] == px_id) & (px_demo["session_id"] == px_ses)
        ]

        # If no such row exists, create an empty DataFrame with the same columns
        if px_demo.empty:
            px_demo = None
            msg = f"Cannot find {bids_id} in demographics file.\nFile: {demo}\n"
        elif px_demo.shape[0] != 1:
            msg = f"Provided {bids_id} is not unique in demographics file.\nFile: {demo}\n"
        else:
            msg = None
        if msg and (cov_normative is not None or cov_deconfound is not None):
            raise ValueError(msg)
        elif msg:
            logger.warning(msg)
            px_demo = None  # Don't use

    # # Run analyses -------------------------------------------------------------
    logger.info("\n\nStarting analysis")
    for label in labels_ctx:
        available_features = run_analysis(
            px_sid=px_id,
            px_ses=px_ses,
            cn_zbrains=zbrains_ref,
            cn_demo_paths=demo_ref,
            px_zbrains=zbrains,
            px_demo=px_demo,
            structures=struct,
            features=feat,
            cov_normative=cov_normative,
            cov_deconfound=cov_deconfound,
            smooth_ctx=smooth_ctx,
            smooth_hip=smooth_hip,
            resolutions=resolution,
            labels_ctx=label,
            labels_hip=labels_hip,
            actual_to_expected=actual_to_expected,
            analyses=LIST_ANALYSES,
            approach=approach,
            col_dtypes=col_dtypes,
            tmp=tmp,
            n_jobs=n_jobs,
        )
    # # Generate volumes ----------------------------------------------------------
    # logger.info("\n\nStarting volume generation")
    # surface_to_volume(
    #     dataset,
    #     feat,
    #     LIST_ANALYSES,
    #     struct,
    #     smooth_ctx,
    #     smooth_hip,
    #     zbrains_ref,
    #     px_id,
    #     px_ses,
    #     px_demo,
    #     micapipename,
    #     hippunfoldname,
    #     tmp,
    #     n_jobs=n_jobs,
    #     n_jobs_wb=n_jobs_wb,
    #     workbench_path=workbench_path,
    #     dicoms=dicoms,
    # )

    # Generate report ----------------------------------------------------------
    logger.info("\n\nStarting report generation")

    threshold = abs(threshold)

    res_hip: Resolution = "high"
    res_ctx: Resolution = "high"
    if "high" not in resolution:
        res_ctx = res_hip = "low"

    lab_ctx = labels_ctx[0][0]
    lab_hip = labels_hip[0]
    print(lab_ctx, lab_hip)
    age = None
    sex = None
    if px_demo is not None:
        age = px_demo.iloc[0].get("age", None)
        sex = px_demo.iloc[0].get("sex", None)

    feat_ctx = available_features["cortex"][res_ctx][lab_ctx]
    feat_sctx = available_features["subcortex"]
    feat_hip = available_features["hippocampus"][res_hip][lab_hip]

    feat_report = list(np.union1d(np.union1d(feat_ctx, feat_sctx), feat_hip))
    # feat_report = [feat for feat in feat_report if "blur" not in feat]
    feat_ctx = [feat for feat in feat_ctx if "blur" not in feat]
    multi = None
    for feat in [feat_ctx, feat_sctx, feat_hip]:
        print("multi: ", feat)
        if multi is None:
            multi = feat
        elif len(feat) > 1 and multi is not None and len(feat) > len(multi):
            multi = feat

        # print(feat)
        # if len(feat) > 1 and feat not in feat_report:
        #     feat_report.append(feat)
    if multi is not None and len(multi) > 0:
        feat_report.append(multi)

    print("feat_report={}".format(feat_report))
    # del feat_report[-1]
    tmp = "/tmp" if tmp is None else tmp
    generate_clinical_report(
        zbrains_path=zbrains,
        sid=px_id,
        ses=px_ses,
        age=age,
        sex=sex,
        analyses=LIST_ANALYSES,
        features=feat_report,
        approach=approach,
        threshold=threshold,
        smooth_ctx=smooth_ctx,
        smooth_hip=smooth_hip,
        res_ctx=res_ctx,
        res_hip=res_hip,
        label_ctx=lab_ctx,
        label_hip=lab_hip,
        tmp_dir=tmp,
    )


def run(
    subject_id,
    zbrains,
    demo_ref,
    zbrains_ref,
    session=None,
    demo=None,
    struct=None,
    feat=None,
    normative=None,
    deconfound=None,
    smooth_ctx=None,
    smooth_hip=None,
    threshold=None,
    approach=None,
    resolution=None,
    labels_ctx=None,
    labels_hip=None,
    logfile=None,
    tmp=None,
    verbose=None,
    filter_warnings=None,
    column_map=None,
    n_jobs=None,
    micapipename=None,
    hippunfoldname=None,
    n_jobs_wb=None,
    workbench_path=None,
    dataset=None,
    dicoms=None,
):

    # Logging settings
    if verbose == 0:
        logging_level = logging.ERROR
    elif verbose == 1:
        logging_level = logging.WARNING
    elif verbose == 2:
        logging_level = logging.INFO
    else:
        logging_level = logging.DEBUG

    # Create a logger
    logger = logging.getLogger(tmp)
    logger.setLevel(logging.DEBUG)  # Default level

    # Create a console handler
    console_handler = logging.StreamHandler(sys.stdout)
    fmt = logging.Formatter("%(message)s")
    console_handler.setFormatter(fmt)
    console_handler.setLevel(logging_level)
    if filter_warnings:
        console_handler.addFilter(lambda record: record.levelno != logging.WARNING)
    logger.addHandler(console_handler)

    # Create a file handler - logs everything
    if logfile is not None:
        file_handler = logging.FileHandler(logfile, mode="w")
        fmt = logging.Formatter("%(asctime)s - %(levelname)-10.10s: %(message)s")
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    # Creating a handler to add unhandled exceptions to logger
    def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            # Will call default excepthook
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        # Create a critical level log message with info from the except hook.
        logger.critical(
            "Unhandled exception", exc_info=(exc_type, exc_value, exc_traceback)
        )

    # Assign the excepthook to the handler
    sys.excepthook = handle_unhandled_exception

    main(
        dataset,
        zbrains_ref,
        demo_ref,
        column_map,
        subject_id,
        session,
        demo,
        zbrains,
        struct,
        feat,
        normative,
        deconfound,
        smooth_ctx,
        smooth_hip,
        threshold,
        approach,
        resolution,
        labels_ctx,
        labels_hip,
        tmp,
        logger,
        micapipename,
        hippunfoldname,
        n_jobs,
        n_jobs_wb,
        workbench_path,
        dicoms,
    )


if __name__ == "__main__":
    # Parse the command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject_id", required=True)
    parser.add_argument("--session", required=True)
    parser.add_argument("--zbrains_ref", required=True)
    parser.add_argument("--demo_ref", required=True)
    parser.add_argument("--zbrains", required=True)
    parser.add_argument("--struct", required=True)
    parser.add_argument("--feat", required=True)
    parser.add_argument("--smooth_ctx", required=True)
    parser.add_argument("--smooth_hip", required=True)
    parser.add_argument("--threshold", required=True)
    parser.add_argument("--approach", required=True)
    parser.add_argument("--resolution", required=True)
    parser.add_argument("--labels_ctx", required=True)
    parser.add_argument("--labels_hip", required=True)
    parser.add_argument("--logfile", required=True)
    parser.add_argument("--tmp", required=True)
    parser.add_argument("--verbose", required=True)
    parser.add_argument("--demo", default=None)
    parser.add_argument("--normative", default=None)
    parser.add_argument("--deconfound", default=None)
    parser.add_argument("--column_map", required=True)
    parser.add_argument("--n_jobs", type=int, required=True)
    parser.add_argument("--n_jobs_wb", type=int, required=True)
    parser.add_argument("--micapipe", type=str, required=True)
    parser.add_argument("--hippunfold", type=str, required=True)
    parser.add_argument("--workbench_path", type=str, required=True)
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--dicoms", type=int, required=True)

    # Parse the arguments.
    args = parser.parse_args()
    args.labels_ctx = eval(args.labels_ctx)
    args.struct = args.struct.split("-")
    args.feat = args.feat.split("-")
    args.demo_ref = args.demo_ref.split("-")
    args.zbrains_ref = args.zbrains_ref.split("-")
    args.resolution = args.resolution.split("-")
    print(args.labels_ctx)
    run(
        subject_id=args.subject_id,
        zbrains=args.zbrains,
        demo_ref=args.demo_ref,
        zbrains_ref=args.zbrains_ref,
        session=args.session,
        demo=args.demo,
        struct=args.struct,
        feat=args.feat,
        normative=args.normative,
        deconfound=args.deconfound,
        smooth_ctx=args.smooth_ctx,
        smooth_hip=args.smooth_hip,
        threshold=float(args.threshold),
        approach=args.approach,
        resolution=args.resolution,
        labels_ctx=[args.labels_ctx],
        labels_hip=[args.labels_hip],
        logfile=args.logfile,
        tmp=args.tmp,
        verbose=int(args.verbose),
        filter_warnings=False,
        column_map=ast.literal_eval(args.column_map),
        n_jobs=int(args.n_jobs),
        n_jobs_wb=int(args.n_jobs_wb),
        micapipename=args.micapipe,
        hippunfoldname=args.hippunfold,
        workbench_path=args.workbench_path,
        dataset=args.dataset,
        dicoms=args.dicoms,
    )
