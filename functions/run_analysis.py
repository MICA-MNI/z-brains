import sys
import logging
import argparse
from pathlib import Path


from utils_analysis import get_id, get_session, run_analysis
from constants import LIST_FEATURES, LIST_STRUCTURES, LIST_ANALYSES, \
    LIST_RESOLUTIONS, DEFAULT_SMOOTH_CTX, DEFAULT_SMOOTH_HIP, \
    DEFAULT_THRESHOLD, FOLDER_NORM_Z, FOLDER_NORM_MODEL


class MapAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        alias_columns = dict(item.split('=') for item in values)
        setattr(namespace, self.dest, alias_columns)


cli = argparse.ArgumentParser(description='Regional and  asymmetry analysis')
cli.add_argument(
    '--subject_id', metavar='ID', type=str, required=True,
    help='Patient ID for target subject. Example: \'sub-HC001\'.'
)
cli.add_argument(
    '--session', metavar='SESSION', type=str, default='',
    help='Session id of target subject. Example: \'ses-002\'.'
)

cli.add_argument(
    '--demo', metavar='PATH', type=Path, default=None,
    help='CSV/TSV files with demographics for target subjects.'
)
cli.add_argument(
    '--zbrains', metavar='PATH', type=Path, required=True,
    help='Path to the zbrains derivative folder containing data for target '
         'subjects.'
)
cli.add_argument(
    '--demo_ref', metavar='PATH', type=Path, nargs='+',
    required=True, help='CSV/TSV files with demographics for reference '
                        'subjects. You must provide one CSV file for each '
                        'zbrains folder in --zbrains_ref.'
)
cli.add_argument(
    '--zbrains_ref', metavar='PATH', type=Path, nargs='+',
    required=True, help='Paths to zbrains derivatives folders containing data '
                        'for reference subjects.'
)
cli.add_argument(
    '--struct', type=str, nargs='+', choices=LIST_STRUCTURES,
    default=LIST_STRUCTURES, help='Structures.'
)
cli.add_argument(
    '--feat', type=str, nargs='+', choices=LIST_FEATURES,
    default=LIST_FEATURES, help='Features.'
)
cli.add_argument(
    '--normative', type=str, nargs='+', default=None,
    help='Perform normative modeling based on provided covariates'
)
cli.add_argument(
    '--deconfound', type=str, nargs='+', default=None,
    help='Perform deconfounding based on provided covariates. If covariates '
         'include \'site\', deconfounding is performed using ComBat. Otherwise, '
         'linear regression is used to regress out the effects of the '
         'covariates from the data. By default, ComBat preserves additional '
         'covariates. To remove the effects of a covariate, prepend with \'-\' '
         '(ignored when not using ComBat). '
         'Example: \'--deconfound site -age -sex group\' to harmonize data while '
         'preserving the effects of group and removing those of age and sex.'
)
cli.add_argument(
    '--smooth_ctx', type=str, default=DEFAULT_SMOOTH_CTX,
    help='Size of gaussian smoothing kernel for cortex, in millimeters.'
)
cli.add_argument(
    '--smooth_hip', type=str, default=DEFAULT_SMOOTH_HIP,
    help='Size of gaussian smoothing kernel for hippocampus, in millimeters.'
)
cli.add_argument(
    '--threshold', type=float, default=DEFAULT_THRESHOLD,
    help='Z-score threshold.'
)
cli.add_argument(
    '--analysis', type=str, choices=LIST_ANALYSES, default='norm-z',
    help='Type of analysis.'
)
cli.add_argument(
    '--resolution', type=str, nargs='+', choices=LIST_RESOLUTIONS,
    default=LIST_RESOLUTIONS,
    help='Resolutions of cortical and hippocampal surfaces.'
)
cli.add_argument(
    '--labels_ctx', type=str, nargs='+', default=['white'],
    help='Cortical surfaces used in the volume to surface mapping. '
         'Example: --labels_ctx white'
)
cli.add_argument(
    '--labels_hip', type=str, nargs='+', default=['midthickness'],
    help='Hippocampal surfaces used in the volume to surface mapping. '
         'Example: --labels_hip midthickness'
)
cli.add_argument(
    '--logfile', metavar='PATH', type=Path, default=None,
    help='Specify the path to the log file'
)
cli.add_argument(
    '--verbose', metavar='LEVEL', type=int, default=-1,
    help='Verbosity level (see zbrains script)'
)
cli.add_argument(
    '--column_map', metavar='VAR=name', nargs='+', default=dict(),
    action=MapAction,
    help='Map expected column names (ID, SES, AGE, SEX, and SITE) to actual '
         'column names in CSV files. Example: --map ID=subject_id'
)

args = cli.parse_args()


################################################################################
# Logging settings
################################################################################
match args.verbose:
    case 0:
        logging_level = logging.ERROR
    case 1:
        logging_level = logging.WARNING
    case 2:
        logging_level = logging.INFO
    case _:
        logging_level = logging.DEBUG


# Create a logger
logger = logging.getLogger('analysis_logger')
logger.setLevel(logging.DEBUG)  # Default level

formatter = logging.Formatter('%(asctime)s - %(levelname)-10.10s: %(message)s')

# Create a console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
console_handler.setLevel(logging_level)
logger.addHandler(console_handler)

# Create a file handler - logs everything
if args.logfile is not None:
    file_handler = logging.FileHandler(args.logfile, mode='w')
    file_handler.setFormatter(formatter)
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


################################################################################
# Analysis
################################################################################
px_id = get_id(args.subject_id, add_prefix=True)
px_ses = get_session(args.session, add_predix=True)
threshold = abs(args.threshold)

analysis_type = args.analysis
analysis_folder = FOLDER_NORM_Z
if analysis_type == "norm-modelling":
    analysis_folder = FOLDER_NORM_MODEL

expected_to_actual = dict(
    participant_id='participant_id', session_id='session_id', age='age',
    sex='sex', site='site'
)
unknown_cols = set(args.column_map.keys()).difference(expected_to_actual.keys())
if len(unknown_cols) > 0:
    raise ValueError(f'Unknown column names: {unknown_cols}. '
                     f'Allowed options are: {list(expected_to_actual.keys())}')
expected_to_actual.update(args.column_map)

actual_to_expected = {v: k for k, v in expected_to_actual.items()}

cov_normative = args.normative
if cov_normative is not None:
    cov_normative = [actual_to_expected.get(col, col) for col in cov_normative]

cov_deconfound = args.deconfound
if cov_deconfound is not None:
    for i, col in enumerate(cov_deconfound):
        prefix = ''
        if col.startswith('-'):
            prefix = '-'
            col = col[1:]
        cov_deconfound[i] = prefix + actual_to_expected.get(col, col)

if len(args.zbrains_ref) != len(args.demo_ref):
    raise ValueError('The number of values provided with --zbrains_ref '
                     'and --demo_ref must be the same.')


def main():
    for asymmetry in [False, True]:
        run_analysis(
            px_sid=px_id, px_ses=px_ses, dir_cn=args.zbrains_ref,
            demo_cn=args.demo_ref, dir_px=args.zbrains, demo_px=args.demo,
            structures=args.struct, features=args.feat,
            cov_normative=cov_normative, cov_deconfound=cov_deconfound,
            smooth_ctx=args.smooth_ctx, smooth_hip=args.smooth_hip,
            resolutions=args.resolution, labels_ctx=args.labels_ctx,
            labels_hip=args.labels_hip, actual_to_expected=actual_to_expected,
            asymmetry=asymmetry, analysis_folder=analysis_folder
        )

        # Generate report - threshold


if __name__ == '__main__':
    main()
