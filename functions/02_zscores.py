import time

import numpy as np
import nibabel as nib
import glob
import sys
import pandas as pd
import warnings

import logging

# Add arg parser
import argparse
from pathlib import Path

from typing import get_args, Literal

# Usage
# python 02_zscores.py --patient PX004 --session 1 --feature All --controls /home/ob/PD/Data/zbrains/control.csv
# --smooth-ctx 5mm --smooth-hipp 2mm  --patient-dir /home/ob/PD/Data/zbrains/
#
# Using propensity scores - TODO - probably remove?
# python 02_zscores.py --patient PX004 --session 1 --feature All --controls /home/ob/PD/Data/zbrains/control.csv
# --smooth-ctx 5mm --smooth-hipp 2mm  --patient-dir /home/ob/PD/Data/zbrains/ --propensity-scores
# --demographics /home/ob/PD/Data/zbrains/PX004.csv --covariates AGE SEX --n-closest 30


# Some defaults - TODO: should we remove this from here?
DEFAULT_PATH_LIST_CONTROLS = '/data_/mica1/03_projects/jessica/hackathon2023/lists/control.csv'

Structure = Literal['cortex', 'hippocampus', 'subcortex']
LIST_STRUCTURES: list[Structure] = list(get_args(Structure))

Feature = Literal['ADC', 'FA', 'T1map', 'thickness']
LIST_FEATURES: list[Feature] = list(get_args(Feature))

Analysis = Literal['norm-z', 'norm-ps', 'norm-normative']

parser = argparse.ArgumentParser(description='Performs z-scoring')
parser.add_argument('-p', '--patient', metavar='ID', type=str, required=True, help='Patient ID')
parser.add_argument('-s', '--session', type=int, default=1, help='Session')
parser.add_argument('-d', '--patient-dir', metavar='DIR', type=Path, required=True,
                    help="Patient directory. This is the output folder from micapipe")
parser.add_argument('-f', '--feature', type=str, nargs='+', choices=LIST_FEATURES + ['All'], default='all',
                    help='Features maps.')

parser.add_argument('-cl', '--controls', metavar='PATH', type=Path, default=DEFAULT_PATH_LIST_CONTROLS,
                    help='File with the list of control subjects. If not specified, use default list.')

parser.add_argument('--smooth-ctx', type=str, help='Size of gaussian smoothing kernel for cortex.', required=True)
parser.add_argument('--smooth-hipp', type=str, help='Size of gaussian smoothing kernel for hippocampus.', required=True)

parser.add_argument('-ps', '--propensity-scores', type=bool, default=False, help='Use propensity scores',
                    action=argparse.BooleanOptionalAction)
parser.add_argument('--covariates', type=str, nargs='+',
                    help='Covariates used to compute prensity scores. Column names in list of control subjects.')
parser.add_argument('--demographics', metavar='PATH', type=Path, required=False,
                    help='CSV file with the patient demographics. Required for propensity scores.')
parser.add_argument('-n', '--n-closest', type=int, default=10, help='Minimun number of controls needed for z-scoring.')
parser.add_argument('-c', '--caliper', type=float, default=0.2, help='Caliper used to select controls for z-scoring')

args = parser.parse_args()

controls_list = args.controls
px_id = args.patient  # patient id
px_id = f'sub-{px_id}'

featureList = args.feature
if 'All' in featureList:
    featureList = LIST_FEATURES
smoothCtx = args.smooth_ctx
smoothHipp = args.smooth_hipp
outDir = args.patient_dir
px_ses = f'ses-{args.session:>02}'

use_ps = args.propensity_scores
covariates = args.covariates
px_demographics = args.demographics
n_min = args.n_closest
caliper = args.caliper

analysis = 'norm-z'
if use_ps:
    analysis = 'norm-propensity'

##################################################################################################
# Logging settings
##################################################################################################
logging_filename = f'{outDir}/{px_id}/{px_ses}/logs/02_zscores.log'
logging.basicConfig(filename=logging_filename, encoding='utf-8',
                    format='%(asctime)s - %(levelname)-10.10s: %(message)s', level=logging.DEBUG)

logger = logging.getLogger('02_zscores')
logger_file = logging.getLogger('02_zscores.file')
logger_both = logging.getLogger('02_zscores.both')
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
# Read/write functions
##################################################################################################
def _load_one(sid: str, ses: str, feat: Feature, struct: Structure = 'cortex', raise_error: bool = True) \
        -> np.ndarray | pd.DataFrame | None:
    """ Load subject data

    Parameters
    ----------
    sid: str
        Subject id in the form 'sub-XXXXX'
    ses: str
        Session, in the form 'ses-XX'
    feat: str
        Feature name
    struct: str, default='cortex'
        Structure.
    raise_error: bool, default=True
        Raise error if file not found

    Returns
    -------
    data: ndarray of shape (n_hemispheres, n_points) or DataFrame of shape (1, n_points,)
        Subject data. If structure is 'cortex' or 'hippocampus', return ndarray of shape (2, n_vertices).
        Otherwise, return DataFrame (1, n_subcortical_structures) including both hemispheres.
        Return None if no data for one or 2 hemispheres.
    """

    folder = f'{outDir}/{sid}/{ses}/maps/{struct}'
    smooth = smoothCtx if struct == 'cortex' else smoothHipp
    if struct == 'subcortex' and feat == 'thickness':
        feat = 'volume'

    if struct == 'subcortex':
        ipth = f'{folder}/{sid}_{ses}_feature-{feat}.csv'
        try:
            return pd.read_csv(ipth, header=[0], index_col=0)
        except FileNotFoundError as e:
            if raise_error:
                logger_file.error(e, stack_info=True, exc_info=True)

                sys.excepthook = sys.__excepthook__  # skip unhandled exception
                raise

            logger_file.warning(f'File not found: "{ipth}"')
            return None

    else:  # cortex & hippocampus
        x = []
        for hemi in ['L', 'R']:
            ipth = f'{folder}/{sid}_{ses}_hemi-{hemi}_feature-{feat}_smooth-{smooth}.func.gii'
            try:
                x.append(nib.load(ipth).darrays[0].data)
            except FileNotFoundError as e:
                if raise_error:
                    logger_file.error(e, stack_info=True, exc_info=True)

                    sys.excepthook = sys.__excepthook__  # skip unhandled exception
                    raise

                logger_file.warning(f'File not found: "{ipth}"')
                return None

        return np.vstack(x)


def _load_data_cn(feat: Feature, df_controls: pd.DataFrame, struct: Structure = 'cortex') \
        -> tuple[pd.DataFrame | np.ndarray, pd.DataFrame]:
    """

    Parameters
    ----------
    feat: str
        Feature name
    df_controls:
        Data frame with control subjects.
    struct: str, default='cortex'
        Structure.

    Returns
    -------
    data: np.ndarray of shape (n_hemispheres, n_subjects, n_points) or  DataFrame of (n_subjects, n_points)
        Data for CN. DataFrame is used for 'subcortex'.
    """

    missing_subjects = np.ones(df_controls.shape[0], dtype=bool)
    data = []
    for i, (sid, subj) in enumerate(df_controls.iterrows()):
        ses = f'ses-{subj["SES"]:>02}'

        x = _load_one(sid, ses, feat, struct=struct, raise_error=False)
        if x is not None:
            data.append(x)
            missing_subjects[i] = False

    if struct == 'subcortex':
        return pd.concat(data, axis=0), df_controls[~missing_subjects]
    return np.stack(data, axis=1), df_controls[~missing_subjects]


def _load_data_px(feat: Feature, struct: Structure = 'cortex') -> np.ndarray | pd.DataFrame:
    """

    Parameters
    ----------
    feat: str
        Feature name
    struct: str, default='cortex'
        Structure.

    Returns
    -------
    data: DataFrame of shape (1, n_points)
        Patient data.
    """

    return _load_one(px_id, px_ses, feat, struct=struct, raise_error=True)


def _save(x: np.ndarray | pd.DataFrame, sid: str, ses: str, feat: Feature, struct: Structure = 'cortex'):
    """ Save results

    Parameters
    ----------
    x: ndarray of shape (2, n_points) or (n_points,)
        Patient data to save
    sid: str
        Subject id in the form 'sub-XXXXX'
    ses: str
        Session, in the form 'ses-XX'
    feat: str
        Feature name
    struct: str, default='cortex'
        Structure.

    """

    folder = f'{outDir}/{sid}/{ses}/{analysis}/{struct}'
    smooth = smoothCtx if struct == 'cortex' else smoothHipp

    if struct == 'subcortex':
        x.to_csv(f'{folder}/{sid}_{ses}_feature-{feat}.csv')
        return

    for i, hemi in enumerate(['L', 'R']):
        data_array = nib.gifti.GiftiDataArray(data=x[i])  # per hemisphere
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)

        opth = f'{folder}/{sid}_{ses}_hemi-{hemi}_feature-{feat}_smooth-{smooth}.func.gii'
        nib.save(image, opth)


##################################################################################################
# Analysis
##################################################################################################


def run_propensity_scores(data_cn: np.ndarray, demo_cn: pd.DataFrame, demo_px: pd.DataFrame) -> np.ndarray:
    from propensity_scores import estimate_ps, get_matches

    mask = pd.notna(demo_cn).all(axis=1)
    demo_cn = demo_cn[mask]
    data_cn = data_cn[mask] if data_cn.ndim == 2 else data_cn[:, mask]

    n = demo_cn.shape[0]
    y = np.zeros(n * 2 + 1)
    y[n:] = 1

    demo_joint = pd.concat([demo_cn, demo_cn, demo_px], axis=0)

    categorical = [c for c in demo_joint.columns if demo_joint[c].dtype == 'O']
    ps = estimate_ps(demo_joint, y, cat=categorical)
    ps_px = ps[-1]
    ps_cn = ps[:n]

    np.empty_like(mask, dtype=float)
    idx_closest = get_matches(ps_px, ps_cn, caliper=None, n_min=n_min, n_max=None)
    data_cn = data_cn[idx_closest] if data_cn.ndim == 2 else data_cn[:, idx_closest]

    return data_cn


def main():
    logger_file.info(f'Logging call: {sys.argv[0]} {" ".join(sys.argv[1:])}')
    logger_both.debug(f'Logging to {logging_filename}\n')

    # Data frame of control subjects
    orig_df_cn = pd.read_csv(controls_list, header=[0], index_col=0)
    orig_n_cn = orig_df_cn.shape[0]

    # Demographics for patient
    df_px = None
    if use_ps:
        df_px = pd.read_csv(px_demographics, header=[0], index_col=0)

    for feat_name in featureList:
        logger_both.info(f'Feature: {feat_name}')
        for structure in LIST_STRUCTURES:

            try:
                data_cn, df_cn = _load_data_cn(feat_name, df_controls=orig_df_cn, struct=structure)
                data_px = _load_data_px(feat_name, struct=structure)
                # print(df_cn)
            except Exception as e:
                continue

            data_px_as_df = None
            if structure == 'subcortex':
                data_px_as_df = data_px
                data_cn = data_cn.to_numpy()
                data_px = data_px.to_numpy()

                # TODO: some have 15 columns - an additional column for ICV - remove
                data_cn = data_cn[..., :14]
                data_px = data_px[..., :14]
                data_px_as_df = data_px_as_df.iloc[:, :14]

            if use_ps:
                data_cn = run_propensity_scores(data_cn, df_cn[covariates], df_px[covariates])

            logger_both.info(f'\tStructure: {structure:<12} \t[Number of controls: {data_cn.shape[-2]}/{orig_n_cn}]')

            # Analysis
            zscores = data_px - np.nanmean(data_cn, axis=-2)
            zscores /= np.nanstd(data_cn, axis=-2)

            if data_px_as_df is not None:
                zscores = pd.DataFrame(zscores, index=data_px_as_df.index, columns=data_px_as_df.columns)

            _save(zscores, px_id, px_ses, feat_name, struct=structure)

    logger_both.info('Done!\n\n')


if __name__ == '__main__':
    main()
