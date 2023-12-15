import sys
import logging
import argparse
from pathlib import Path

import numpy as np

from utils_analysis import (ALIAS_COLUMNS, LIST_FEATURES, LIST_STRUCTURES, LIST_ANALYSES,
                            LIST_RESOLUTIONS, CONFIG, get_id, get_session, run_analysis)


# Usage
# python run_analysis.py --subject_id $id --session $session --zbrains_dir $out_dir \
#                        --controls $demo --structure cortex hippocampus --feature feat1 feat2 \
#                        --smooth-ctx $smooth_ctx --smooth-hip $smooth_hip \
#                        --resolution all --threshold $threshold --asymmetry


class MapAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        alias_columns = dict(item.split('=') for item in values)
        setattr(namespace, self.dest, alias_columns)


cli = argparse.ArgumentParser(description='Performs regional/asymmetry analysis')
cli.add_argument('--subject_id', metavar='ID', type=str, required=True,
                 help='Patient ID')
cli.add_argument('--session', metavar='SESSION', type=str, default='SINGLE',
                 help='Session id. Use "SINGLE" for single session.')
cli.add_argument('--path_zbrains', metavar='PATH', type=Path, required=True,
                 help="Path to zbrains derivatives folder.")
cli.add_argument('--struct', type=str, nargs='+', choices=LIST_STRUCTURES, required=True,
                 help='Structures.')
cli.add_argument('--feat', type=str, nargs='+', choices=LIST_FEATURES, required=True,
                 help='Features.')
cli.add_argument('--demo_cn', metavar='PATH', type=Path, required=True,
                 help='CSV file with demographics for controls.')
cli.add_argument('--demo_px', metavar='PATH', type=Path, default=None,
                 help='CSV file with demographics for the patient.')
cli.add_argument('--normative', type=str, nargs='+', default=None,
                 help='Perform normative modeling based on provided covariates')
cli.add_argument('--deconfound', type=str, nargs='+', default=None,
                 help='Perform deconfounding based on provided covariates')
cli.add_argument('--smooth_ctx', type=float, required=True,
                 help='Size of gaussian smoothing kernel for cortex, in millimeters.')
cli.add_argument('--smooth_hip', type=float, required=True,
                 help='Size of gaussian smoothing kernel for hippocampus, in millimeters.')
cli.add_argument('--threshold', type=float, default=CONFIG['DEFAULT_THRESHOLD'],
                 help='Z-score threshold.')
cli.add_argument('--analysis', type=str, choices=LIST_ANALYSES, default='norm-z',
                 help='Type of analysis.')
cli.add_argument('--resolution', type=str, nargs='+', choices=LIST_RESOLUTIONS,
                 required=True, help='Resolutions for cortical and hippocampal surfaces.')
cli.add_argument('--labels_ctx', type=str, nargs='+', required=True,
                 help='Names of cortical surfaces used in the volume to surface mapping. '
                      'Example: --labels_ctx white')
cli.add_argument('--labels_hip', type=str, nargs='+', required=True,
                 help='Names of hippocampal surfaces used in the volume to surface mapping. '
                      'Example: --labels_hip midthickness')
cli.add_argument('--logfile', metavar='PATH', type=Path, required=True,
                 help='Specify the path to the log file.')
cli.add_argument('--map', metavar='VAR=name', nargs='+', type=str, default=dict(),
                 action=MapAction, help='Map expected column names (ID, SES, AGE, SEX, and SITE) '
                                        'to actual column names. Example: --map ID=subject_id')


args = cli.parse_args()
px_id = get_id(args.sucject_id, add_prefix=True)
px_ses = get_session(args.session, add_predix=True)


threshold = np.abs(args.threshold)
analysis_type = args.analysis

analysis_folder = CONFIG["FOLDER_NORM_Z"]
if analysis_type == "norm-modelling":
    analysis_folder = CONFIG["FOLDER_NORM_MODEL"]

unknown_cols = list(set(args.map.keys()).difference(ALIAS_COLUMNS.keys()))
if len(unknown_cols) > 0:
    raise ValueError(f'Unknown column names: {unknown_cols}. '
                     f'Allowed options are: {list(ALIAS_COLUMNS.keys())}')
ALIAS_COLUMNS.update(args.map)

reverse_alias_columns = {v: k for k, v in ALIAS_COLUMNS.items()}
cov_deconfound = [reverse_alias_columns.get(col, col) for col in args.deconfound]
cov_normative = [reverse_alias_columns.get(col, col) for col in args.normative]


# ##################################################################################################
# # Logging settings
# ##################################################################################################
logging.basicConfig(filename=args.logfile, encoding='utf-8', filemode='w',
                    format='%(asctime)s - %(levelname)-10.10s: %(message)s', level=logging.DEBUG)

logger = logging.getLogger('run_analysis')
logger_file = logging.getLogger('analysis_logger.file')
logger_both = logging.getLogger('analysis_logger.both')
logger_both.addHandler(logging.StreamHandler(sys.stdout))


# Creating a handler to add unhandled exceptions to logger
def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        # Will call default excepthook
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
        # Create a critical level log message with info from the except hook.
    logger_both.critical("Unhandled exception", exc_info=(exc_type, exc_value, exc_traceback))


# Assign the excepthook to the handler
sys.excepthook = handle_unhandled_exception


##################################################################################################
# Analysis
##################################################################################################
def main():
    for asymmetry in [False, True]:
        run_analysis(args.path_zbrains, analysis_folder=analysis_folder, px_sid=px_id,
                     path_demo_cn=args.demo_cn, structures=args.struct, features=args.feat,
                     resolutions=args.resolution, labels_ctx=args.labels_ctx,
                     labels_hip=args.labels_hip, smooth_ctx=args.smooth_ctx,
                     smooth_hip=args.smooth_hip, px_ses=px_ses, path_demo_px=args.demo_px,
                     asymmetry=asymmetry, col_aliases=reverse_alias_columns,
                     cov_deconfound=cov_deconfound)


if __name__ == '__main__':
    main()
