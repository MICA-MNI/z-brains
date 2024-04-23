import sys
import logging
import argparse
from pathlib import Path

import numpy as np

from functions.new_constants import (
    LIST_FEATURES, LIST_STRUCTURES, LIST_APPROACHES, LIST_RESOLUTIONS,
    DEFAULT_SMOOTH_CTX, DEFAULT_SMOOTH_HIP, DEFAULT_THRESHOLD, LIST_ANALYSES,
    Resolution,
)

from functions.utils_analysis import (
    get_id, get_session, run_analysis, get_bids_id, load_demo
)
from functions.clinical_reports import generate_clinical_report




# class MapAction(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         alias_columns = dict(item.split('=') for item in values)
#         setattr(namespace, self.dest, alias_columns)


# cli = argparse.ArgumentParser(description='Regional and  asymmetry analysis')
# cli.add_argument(
#     '--subject_id', metavar='ID', type=str, required=True,
#     help='Patient ID for target subject. Example: \'sub-HC001\'.'
# )
# cli.add_argument(
#     '--session', metavar='SESSION', type=str, default='',
#     help='Session id of target subject. Example: \'ses-002\'.'
# )

# cli.add_argument(
#     '--demo', metavar='PATH', type=Path, default=None,
#     help='CSV/TSV files with demographics for target subjects.'
# )
# cli.add_argument(
#     '--zbrains', metavar='PATH', type=Path, required=True,
#     help='Path to the zbrains derivative folder containing data for target '
#          'subjects.'
# )
# cli.add_argument(
#     '--demo_ref', metavar='PATH', type=Path, nargs='+',
#     required=True, help='CSV/TSV files with demographics for reference '
#                         'subjects. You must provide one CSV file for each '
#                         'zbrains folder in --zbrains_ref.'
# )
# cli.add_argument(
#     '--zbrains_ref', metavar='PATH', type=Path, nargs='+',
#     required=True, help='Paths to zbrains derivatives folders containing data '
#                         'for reference subjects.'
# )
# cli.add_argument(
#     '--struct', type=str, nargs='+', choices=LIST_STRUCTURES,
#     default=LIST_STRUCTURES, help='Structures.'
# )
# cli.add_argument(
#     '--feat', type=str, nargs='+', choices=LIST_FEATURES,
#     default=LIST_FEATURES, help='Features.'
# )
# cli.add_argument(
#     '--normative', type=str, nargs='+', default=None,
#     help='Perform normative modeling based on provided covariates'
# )
# cli.add_argument(
#     '--deconfound', type=str, nargs='+', default=None,
#     help='Perform deconfounding based on provided covariates. If covariates '
#          'include \'site\', deconfounding is performed using ComBat. '
#          'Otherwise, linear regression is used to regress out the effects of '
#          'the covariates from the data. By default, ComBat preserves '
#          'additional covariates. To remove the effects of a covariate, prepend '
#          'with \'-\' (ignored when not using ComBat). '
#          'Example: \'--deconfound site -age -sex group\' to harmonize data '
#          'while preserving the effects of group and removing those of age and '
#          'sex.'
# )
# cli.add_argument(
#     '--smooth_ctx', type=str, default=DEFAULT_SMOOTH_CTX,
#     help='Size of gaussian smoothing kernel for cortex, in millimeters.'
# )
# cli.add_argument(
#     '--smooth_hip', type=str, default=DEFAULT_SMOOTH_HIP,
#     help='Size of gaussian smoothing kernel for hippocampus, in millimeters.'
# )
# cli.add_argument(
#     '--threshold', type=float, default=DEFAULT_THRESHOLD,
#     help='Z-score threshold.'
# )
# cli.add_argument(
#     '--approach', type=str, choices=LIST_APPROACHES,
#     default='zscore', help='Comparison approach.'
# )
# cli.add_argument(
#     '--resolution', type=str, nargs='+', choices=LIST_RESOLUTIONS,
#     default=LIST_RESOLUTIONS,
#     help='Resolutions of cortical and hippocampal surfaces.'
# )
# cli.add_argument(
#     '--labels_ctx', type=str, nargs='+', default=['white'],
#     help='Cortical surfaces used in the volume to surface mapping. '
#          'Example: --labels_ctx white'
# )
# cli.add_argument(
#     '--labels_hip', type=str, nargs='+', default=['midthickness'],
#     help='Hippocampal surfaces used in the volume to surface mapping. '
#          'Example: --labels_hip midthickness'
# )
# cli.add_argument(
#     '--logfile', metavar='PATH', type=Path, default=None,
#     help='Specify the path to the log file'
# )
# cli.add_argument(
#     '--tmp', metavar='PATH', type=Path, default=None,
#     help='Temporary folder'
# )
# cli.add_argument(
#     '--verbose', metavar='LEVEL', type=int, default=-1,
#     help='Verbosity level (see zbrains script)'
# )
# cli.add_argument(
#     '--filter_warnings', default=False, action='store_true',
#     help='Filter warnings messages'
# )
# cli.add_argument(
#     '--column_map', metavar='VAR=name', nargs='+', default=dict(),
#     action=MapAction,
#     help='Map expected column names (ID, SES, AGE, SEX, and SITE) to actual '
#          'column names in CSV files. Example: --map ID=subject_id'
# )

# args = cli.parse_args()


# ################################################################################
# # Logging settings
# ################################################################################
# match args.verbose:
#     case 0:
#         logging_level = logging.ERROR
#     case 1:
#         logging_level = logging.WARNING
#     case 2:
#         logging_level = logging.INFO
#     case _:
#         logging_level = logging.DEBUG


# # Create a logger
# logger = logging.getLogger('analysis_logger')
# logger.setLevel(logging.DEBUG)  # Default level

# # formatter = logging.Formatter('%(asctime)s - %(levelname)-10.10s: %(message)s')

# # Create a console handler
# console_handler = logging.StreamHandler(sys.stdout)
# fmt = logging.Formatter('%(message)s')
# console_handler.setFormatter(fmt)
# console_handler.setLevel(logging_level)
# if args.filter_warnings:
#     console_handler.addFilter(lambda record: record.levelno != logging.WARNING)
# logger.addHandler(console_handler)

# # Create a file handler - logs everything
# if args.logfile is not None:
#     file_handler = logging.FileHandler(args.logfile, mode='w')
#     fmt = logging.Formatter('%(asctime)s - %(levelname)-10.10s: %(message)s')
#     file_handler.setFormatter(fmt)
#     logger.addHandler(file_handler)


# # Creating a handler to add unhandled exceptions to logger
# def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
#     if issubclass(exc_type, KeyboardInterrupt):
#         # Will call default excepthook
#         sys.__excepthook__(exc_type, exc_value, exc_traceback)
#         return
#         # Create a critical level log message with info from the except hook.
#     logger.critical("Unhandled exception",
#                     exc_info=(exc_type, exc_value, exc_traceback))


# # Assign the excepthook to the handler
# sys.excepthook = handle_unhandled_exception


################################################################################
# Analysis
################################################################################
def main(zbrains_ref,demo_ref,column_map,subject_id,session,demo,zbrains,struct,feat,normative,deconfound,smooth_ctx,smooth_hip,threshold,approach,resolution,labels_ctx,labels_hip,tmp,logger):
    # Some checks
    # logger = logging.getLogger(tmp)
    if len(zbrains_ref) != len(demo_ref):
        raise ValueError('The number of values provided with --zbrains_ref '
                         'and --demo_ref must be the same.')

    # Handle columns names -----------------------------------------------------
    # Column mapping
    expected_to_actual = dict(
        participant_id='participant_id', session_id='session_id', age='age',
        sex='sex', site='site'
    )
    unknown_cols = set(column_map.keys()).difference(expected_to_actual.keys())
    if unknown_cols:
        raise ValueError(
            f'Unknown column names: {unknown_cols}. Allowed options '
            f'are: {list(expected_to_actual.keys())}')

    expected_to_actual.update(column_map)
    actual_to_expected = {v: k for k, v in expected_to_actual.items()}

    # Column types
    col_dtypes = {actual: (str if expected != 'age' else float)
                  for actual, expected in actual_to_expected.items()}

    # Rename covariates for normative modeling
    cov_normative = normative
    if cov_normative is not None:
        cov_normative = [actual_to_expected.get(col, col)
                         for col in cov_normative]

    # Rename covariates for deconfounding
    cov_deconfound = deconfound
    if cov_deconfound is not None:
        for i, col in enumerate(cov_deconfound):
            prefix = ''
            if col.startswith('-'):
                prefix = '-'
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
        px_demo = px_demo.loc[(px_demo['participant_id'] == px_id) & (px_demo['session_id'] == px_ses)]

        # If no such row exists, create an empty DataFrame with the same columns
        if px_demo.empty:
            px_demo = None
            msg = f'Cannot find {bids_id} in demographics file.\nFile: {demo}\n'
        elif px_demo.shape[0] != 1:
            msg = f'Provided {bids_id} is not unique in demographics file.\nFile: {demo}\n'
        else:
            msg = None

        if msg and (cov_normative is not None or cov_deconfound is not None):
            raise ValueError(msg)
        elif msg:
            logger.warning(msg)
            px_demo = None  # Don't use

    # Run analyses -------------------------------------------------------------
    logger.info('\n\nStarting analysis')
    available_features = run_analysis(
        px_sid=px_id, px_ses=px_ses, cn_zbrains=zbrains_ref,
        cn_demo_paths=demo_ref, px_zbrains=zbrains, px_demo=px_demo,
        structures=struct, features=feat, cov_normative=cov_normative,
        cov_deconfound=cov_deconfound, smooth_ctx=smooth_ctx,
        smooth_hip=smooth_hip, resolutions=resolution,
        labels_ctx=labels_ctx, labels_hip=labels_hip,
        actual_to_expected=actual_to_expected, analyses=LIST_ANALYSES,
        approach=approach, col_dtypes=col_dtypes,tmp=tmp
        )

    # Generate report ----------------------------------------------------------
    logger.info('\n\nStarting report generation')

    threshold = abs(threshold)

    res_hip: Resolution = 'high'
    res_ctx: Resolution = 'high'
    if 'high' not in resolution:
        res_ctx = res_hip = 'low'

    lab_ctx = labels_ctx[0]
    lab_hip = labels_hip[0]

    age = None
    sex = None
    if px_demo is not None:
        age = px_demo.iloc[0].get('age', None)
        sex = px_demo.iloc[0].get('sex', None)

    feat_ctx = available_features['cortex'][res_ctx][lab_ctx]
    feat_sctx = available_features['subcortex']
    feat_hip = available_features['hippocampus'][res_hip][lab_hip]

    feat_report = list(np.union1d(np.union1d(feat_ctx, feat_sctx), feat_hip))
    multi = None
    for feat in [feat_ctx, feat_sctx, feat_hip]:
        if multi is None:
            multi = feat
        elif len(feat) > 1 and multi is not None and len(feat) > len(multi):
            multi = feat

        # print(feat)
        # if len(feat) > 1 and feat not in feat_report:
        #     feat_report.append(feat)
    if multi is not None and len(multi) > 0:
        feat_report.append(multi)

    print('feat_report={}'.format(feat_report))
    print(feat_report)
    exit()
    tmp = '/tmp' if tmp is None else tmp
    generate_clinical_report(
        zbrains_path=zbrains, sid=px_id, ses=px_ses, age=age, sex=sex,
        analyses=LIST_ANALYSES, features=feat_report, approach=approach,
        threshold=threshold, smooth_ctx=smooth_ctx,
        smooth_hip=smooth_hip, res_ctx=res_ctx, res_hip=res_hip,
        label_ctx=lab_ctx, label_hip=lab_hip, tmp_dir=tmp
    )

def run(subject_id, zbrains, demo_ref, zbrains_ref, session=None, demo=None, struct=None, feat=None, normative=None, deconfound=None, smooth_ctx=None, smooth_hip=None, threshold=None, approach=None, resolution=None, labels_ctx=None, labels_hip=None, logfile=None, tmp=None, verbose=None, filter_warnings=None, column_map=None):
    # # Check if mandatory options are provided
    # if not subject_id:
    #     raise ValueError("Error: --subject_id is required")
    # if not zbrains:
    #     raise ValueError("Error: --zbrains is required")
    # if not demo_ref:
    #     raise ValueError("Error: --demo_ref is required")
    # if not zbrains_ref:
    #     raise ValueError("Error: --zbrains_ref is required")

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
    fmt = logging.Formatter('%(message)s')
    console_handler.setFormatter(fmt)
    console_handler.setLevel(logging_level)
    if filter_warnings:
        console_handler.addFilter(lambda record: record.levelno != logging.WARNING)
    logger.addHandler(console_handler)

    # Create a file handler - logs everything
    if logfile is not None:
        file_handler = logging.FileHandler(logfile, mode='w')
        fmt = logging.Formatter('%(asctime)s - %(levelname)-10.10s: %(message)s')
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    # Creating a handler to add unhandled exceptions to logger
    def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            # Will call default excepthook
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        # Create a critical level log message with info from the except hook.
        logger.critical("Unhandled exception",
                        exc_info=(exc_type, exc_value, exc_traceback))

    # Assign the excepthook to the handler
    sys.excepthook = handle_unhandled_exception
    
    main(zbrains_ref,demo_ref,column_map,subject_id,session,demo,zbrains,struct,feat,normative,deconfound,smooth_ctx,smooth_hip,threshold,approach,resolution,labels_ctx,labels_hip,tmp,logger)
if __name__ == '__main__':
    main()
