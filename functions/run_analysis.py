import os
import sys
import logging
import argparse
import time
from pathlib import Path
from functools import reduce
from typing import get_args, Literal

import numpy as np
import pandas as pd
import nibabel as nib


# Usage
# cmd = "python run_analysis.py --subject_id $id --session $session --zbrains_dir $out_dir \
#                               --controls $demo --structure cortex hippocampus --feature feat1 feat2 \
#                               --smooth-ctx $smooth_ctx --smooth-hipp $smooth_hipp \
#                               --resolution all --threshold $threshold --asymmetry"

def _read_config(file_path):
    conf = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if len(line) == 0 or line.startswith('#'):
                continue

            key, value = map(str.strip, line.split('='))
            if value.startswith('('):  # Array
                value = [v.strip('"\'') for v in value.strip('()').split()]
            else:
                value = value.strip('"\'')
            conf[key] = value

    return conf


config_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'config.cfg')
config = _read_config(config_file)


map_feature = dict(thickness='thickness', volume='volume', flair='flair', ADC='ADC', FA='FA', qT1='T1map')

Structure = Literal['cortex', 'hippocampus', 'subcortex']
LIST_STRUCTURES: list[Structure] = list(get_args(Structure))
map_struct_to_folder = dict(cortex=config['FOLDER_CTX'], subcortex=config['FOLDER_SCTX'],
                            hippocampus=config['FOLDER_HIPP'])

Feature = Literal['flair', 'ADC', 'FA', 'qT1', 'thickness']
LIST_FEATURES: list[Feature] = list(get_args(Feature))

Analysis = Literal['norm-z', 'norm-modelling']
LIST_ANALYSES: list[Analysis] = list(get_args(Analysis))

Resolution = Literal['low', 'high']
LIST_RESOLUTIONS: list[Resolution] = list(get_args(Resolution))

parser = argparse.ArgumentParser(description='Performs z-scoring')
parser.add_argument('-id', '--subject_id', metavar='ID', type=str, required=True, help='Patient ID')
parser.add_argument('-s', '--session', type=int, default=1, help='Session')
parser.add_argument('-zd', '--zbrains_dir', metavar='DIR', type=Path, required=True,
                    help="This is the z-brains derivative folder.")
parser.add_argument('-f', '--feature', type=str, nargs='+', choices=LIST_FEATURES + ['all'], default='all',
                    help='Features maps.')
parser.add_argument('-st', '--structure', type=str, nargs='+', choices=LIST_STRUCTURES + ['all'],
                    default='all', help='Structures.')
parser.add_argument('-d', '--demo', metavar='PATH', type=Path, required=True,
                    help='CSV File with the list of control subjects, and potentially other demographic info.')
parser.add_argument('--smooth-ctx', type=str,
                    help='Size of gaussian smoothing kernel for cortex (e.g., 5mm).', required=True)
parser.add_argument('--smooth-hipp', type=str,
                    help='Size of gaussian smoothing kernel for hippocampus (e.g., 2mm).', required=True)
parser.add_argument('--asymmetry', action='store_true', default=False,
                    help='Perform asymmetry analysis instead of regional analysis.')
parser.add_argument('--threshold', type=float, default=config['DEFAULT_THRESHOLD'],
                    help='Z-score threshold.')
parser.add_argument('--analysis', type=str, choices=LIST_ANALYSES, default='norm-z',
                    help='Type of analysis.')
parser.add_argument('-r', '--resolution', type=str, nargs='+', choices=LIST_RESOLUTIONS + ['all'],
                    default='all', help='Resolution for cortical and hippocampal surfaces.')

args = parser.parse_args()

controls_list = args.demo
px_id = args.subject_id  # patient id
px_id = f'sub-{px_id}'

features = args.feature
if 'all' in features:
    features = LIST_FEATURES
structures = args.structure
if 'all' in structures:
    structures = LIST_STRUCTURES
smooth_ctx = args.smooth_ctx
smooth_hipp = args.smooth_hipp
zbrains_dir = args.zbrains_dir
px_ses = f'ses-{args.session:>02}'
do_asymmetry = args.asymmetry
threshold = np.abs(args.threshold)
analysis_type = args.analysis

analysis_folder = config["FOLDER_NORM_Z"]
if analysis_type == "norm-modelling":
    analysis_folder = config["FOLDER_NORM_MODEL"]

resolutions = args.resolution
if 'all' in resolutions:
    resolutions = LIST_RESOLUTIONS

subject_dir = f'{zbrains_dir}/{px_id}/{px_ses}'
pth_analysis = f'{subject_dir}/{analysis_folder}'
pth_logs = f'{subject_dir}/{config["FOLDER_LOGS"]}'

##################################################################################################
# Logging settings
##################################################################################################
logging_filename = f'{pth_logs}/analysis.txt'
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
def mahalanobis_distance(x_cn: np.ndarray, x_px: np.ndarray):
    """Compute mahalanobis distance for multiple vertices

    Parameters
    ----------
    x_cn: ndarray of shape=(n_subjects, n_vertices, n_features)
         Data for controls.
    x_px: ndarray of shape=(n_vertices, n_features)
        Patient data.

    Returns
    -------
    dist: np.ndarray of shape=(n_vertices,)
        Mahalanobis distance.
    """

    n_subjects, n_vertices, _ = x_cn.shape

    dist = np.zeros(n_vertices, dtype=np.float32)
    for i, xp in enumerate(x_px):
        xc = x_cn[:, i]
        mu = xc.mean(axis=0)

        xc = xc - mu
        xp = xp - mu

        cov = xc.T @ xc
        cov /= n_subjects - 1
        try:
            cov_inv = np.linalg.inv(cov)
            dist[i] = np.sqrt(xp.T @ cov_inv @ xp)
        except np.linalg.LinAlgError:
            pass  # set distance to zero (default)

    return dist

    # Faster but cannot handle singular matrix exception
    # x_cn = x_cn.swapaxes(0, 1)
    # n_vertices, n_subjects, n_features = x_cn.shape
    #
    # mu = x_cn.mean(axis=1, keepdims=True)
    # x_cn_centered = x_cn - mu
    #
    # cov = (x_cn_centered.swapaxes(1, 2) @ x_cn_centered)
    # cov /= n_subjects - 1
    # cov_inv = np.linalg.inv(cov)
    #
    # x_px_centered = x_px - mu.squeeze()
    # dist = np.sqrt(x_px_centered[:, None] @ cov_inv @ x_px_centered[..., None]).squeeze()
    # return dist


def _map_resolution(struct: Structure, resolution: Resolution):
    if struct == 'cortex':
        res_suffix = config['HIGH_RESOLUTION_CTX'] if resolution == 'high' else config['LOW_RESOLUTION_CTX']
        return f'surf-fsLR-{res_suffix}'

    res_suffix = config['HIGH_RESOLUTION_HIPP'] if resolution == 'high' else config['LOW_RESOLUTION_HIPP']
    return f'den-{res_suffix}'


def _load_one(sid: str, ses: str, feat: Feature, struct: Structure = 'cortex', resolution: Resolution | None = None,
              raise_error: bool = True) -> np.ndarray | pd.DataFrame | None:
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
    resolution: str, optional
        Resolution. Required for cortex and hippocampus.
    raise_error: bool, default=True
        Raise error if file not found

    Returns
    -------
    data: ndarray of shape (n_hemispheres, n_points) or DataFrame of shape (1, n_points,)
        Subject data. If structure is 'cortex' or 'hippocampus', return ndarray of shape (2, n_vertices).
        Otherwise, return DataFrame (1, n_subcortical_structures) including both hemispheres.
        Return None if no data for one or 2 hemispheres.
    """

    folder = f'{zbrains_dir}/{sid}/{px_ses}/{config["FOLDER_MAPS"]}/{map_struct_to_folder[struct]}'
    smooth = smooth_ctx if struct == 'cortex' else smooth_hipp
    if struct == 'subcortex' and feat == 'thickness':
        feat = 'volume'
    feat = map_feature[feat]

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
            res = _map_resolution(struct, resolution)
            ipth = f'{folder}/{sid}_{ses}_hemi-{hemi}_{res}_feature-{feat}_smooth-{smooth}mm.func.gii'
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


def _load_data_cn(feat: Feature, df_controls: pd.DataFrame, struct: Structure = 'cortex',
                  resolution: Resolution | None = None) -> tuple[pd.DataFrame | np.ndarray, pd.DataFrame]:
    """

    Parameters
    ----------
    feat: str
        Feature name
    df_controls:
        Data frame with control subjects.
    struct: str, default='cortex'
        Structure.
    resolution: str, optional
        Resolution. Required for cortex and hippocampus.

    Returns
    -------
    data: np.ndarray of shape (n_subjects, n_hemispheres, n_points) or DataFrame of (n_subjects, n_points)
        Data for CN. DataFrame is used for 'subcortex'.
    df_controls_available: DataFrame of shape (n_available_subjects, n_points)
        Dataframe only including those rows in 'df_controls' with available data in the maps folder
    """

    missing_subjects = np.ones(df_controls.shape[0], dtype=bool)
    data = []
    for i, (sid, subj) in enumerate(df_controls.iterrows()):
        ses = f'ses-{subj["SES"]:>02}'

        x = _load_one(sid, ses, feat, struct=struct, resolution=resolution, raise_error=False)
        if x is not None:
            data.append(x)
            missing_subjects[i] = False

    if struct == 'subcortex':
        return pd.concat(data, axis=0), df_controls[~missing_subjects]
    return np.stack(data, axis=0), df_controls[~missing_subjects]


def _load_data_px(feat: Feature, struct: Structure = 'cortex', resolution: Resolution | None = None) \
        -> np.ndarray | pd.DataFrame:
    """

    Parameters
    ----------
    feat: str
        Feature name
    struct: str, default='cortex'
        Structure.
    resolution: str, optional
        Resolution. Required for cortex and hippocampus.

    Returns
    -------
    data: ndarray of shape (n_hemispheres, n_points) or DataFrame of shape (1, n_points,)
        Patient data. If structure is 'cortex' or 'hippocampus', return ndarray of shape (2, n_vertices).
        Otherwise, return DataFrame (1, n_subcortical_structures) including both hemispheres.
    """

    return _load_one(px_id, px_ses, feat, struct=struct, resolution=resolution, raise_error=True)


def _save(x: np.ndarray | pd.DataFrame, sid: str, ses: str, feat: Feature | list[Feature], struct: Structure = 'cortex',
          resolution: Resolution | None = None, thresh: float | None = None):
    """ Save results

    Parameters
    ----------
    x: ndarray of shape (2, n_points) or (n_points,)
        Patient data to save
    sid: str
        Subject id in the form 'sub-XXXXX'
    ses: str
        Session, in the form 'ses-XX'
    feat: str or list[str]
        Feature name or names for Mahalanobis
    struct: str, default='cortex'
        Structure.
    resolution: str, optional
        Resolution. Required for cortex and hippocampus.
    thresh: float, optional
        Z-score threshold.
    """

    folder = f'{pth_analysis}/{map_struct_to_folder[struct]}'
    smooth = smooth_ctx if struct == 'cortex' else smooth_hipp

    # Handle the case when feat is a list of string (used for Mahalanobis)
    is_list = isinstance(feat, list)
    feat = feat if is_list else [feat]
    feat = [map_feature['volume' if k == 'thickness' else k] for k in feat]
    feat = '-'.join(feat) if is_list else feat[0]

    reg_or_asym = 'asymmetry' if do_asymmetry else 'regional'

    if struct == 'subcortex':
        fname_prefix = f'{folder}/{sid}_{ses}_feature-{feat}_analysis-{reg_or_asym}'
        if thresh is None:
            x.to_csv(f'{fname_prefix}.csv')
        else:
            x.to_csv(f'{fname_prefix}_threshold-{thresh}.csv')
        return

    res = _map_resolution(struct, resolution)
    for i, h in enumerate(['L', 'R']):
        if do_asymmetry:
            data_array = nib.gifti.GiftiDataArray(data=x)  # one hemisphere
        else:
            data_array = nib.gifti.GiftiDataArray(data=x[i])  # per hemisphere
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)

        opth_prefix = f'{folder}/{sid}_{ses}_hemi-{h}_{res}_feature-{feat}_smooth-{smooth}mm_analysis-{reg_or_asym}'
        opth = f'{opth_prefix}.func.gii' if thresh is None else f'{opth_prefix}_threshold-{thresh}.func.gii'
        nib.save(image, opth)

        if do_asymmetry:
            break


def _prepare_data(data_cn: np.ndarray | pd.DataFrame, data_px: np.ndarray | pd.DataFrame, make_asymmetric: bool = False):
    df_index, df_cols = None, None
    if isinstance(data_px, pd.DataFrame):
        df_index, df_cols = data_px.index, data_px.columns

        # TODO: some have 15 columns - an additional column for ICV - remove
        data_cn = data_cn.to_numpy()[..., :14].reshape(-1, 2, 7)
        data_px = data_px.to_numpy()[..., :14].reshape(2, 7)
        df_cols = df_cols[:14]

    if make_asymmetric:
        lh_data_cn, rh_data_cn = data_cn[:, 0], data_cn[:, 1]
        den_cn = .5 * (lh_data_cn + rh_data_cn)
        data_cn = np.divide(lh_data_cn - rh_data_cn, den_cn, out=den_cn, where=den_cn > 0)

        lh_data_px, rh_data_px = data_px[0], data_px[1]
        den_px = .5 * (lh_data_px + rh_data_px)
        data_px = np.divide(lh_data_px - rh_data_px, den_px, out=den_px, where=den_px > 0)

        if df_cols is not None:
            df_cols = df_cols[:df_cols.size // 2]

    return data_cn, data_px, df_index, df_cols


def _save_results(x: np.ndarray, feat: Feature | list[Feature], index: list[str] | None = None,
                  columns: list[str] | None = None, **kwargs):

    # Add header and index to dataframe if available
    if columns is not None:
        x = pd.DataFrame(x.reshape(1, -1), index=index, columns=columns)

    _save(x, px_id, px_ses, feat, **kwargs)


##################################################################################################
# Analysis
##################################################################################################
def main():
    logger_file.info(f'Logging call: {sys.argv[0]} {" ".join(sys.argv[1:])}')
    logger_both.debug(f'Logging to {logging_filename}\n')

    # Data frame of control subjects
    orig_df_cn = pd.read_csv(controls_list, header=[0], index_col=0)
    orig_n_cn = orig_df_cn.shape[0]

    struct_res_pair = []
    for st in structures:
        if st == 'subcortex':
            struct_res_pair.append((st, None))
        else:
            struct_res_pair.extend([(st, res) for res in resolutions])

    for struct_name, res in struct_res_pair:
        msg = f' [res = {res:<4}]' if struct_name != 'subcortex' else ''
        logger_both.info(f'\nStructure: {struct_name}{msg}')

        # Collect data for Mahalanobis
        list_df_cn = []  # Because some subjects may not have some features
        list_data_cn = []
        list_data_px = []
        index_df, cols_df = None, None
        available_feat = []

        for feat_name in features:
            # Load data
            try:
                data_cn, df_cn = _load_data_cn(feat_name, df_controls=orig_df_cn, struct=struct_name, resolution=res)
                data_px = _load_data_px(feat_name, struct=struct_name, resolution=res)
            except FileNotFoundError:  # When patient files not available
                logger_both.info(f'\t{feat_name:<15}: \tNo data available.')
                continue

            data_cn, data_px, index_df, cols_df = _prepare_data(data_cn, data_px, make_asymmetric=do_asymmetry)

            # Save data for mahalanobis
            available_feat.append(feat_name)
            list_data_cn.append(data_cn)
            list_df_cn.append(df_cn)
            list_data_px.append(data_px)

            # Analysis: z-scoring
            mask = np.any(data_cn != 0, axis=0)  # ignore all zeros
            data_cn[:, ~mask] = np.finfo(float).eps

            z = data_px - np.nanmean(data_cn, axis=0)
            z = np.divide(z, np.nanstd(data_cn, axis=0), out=z, where=mask)
            z[~mask] = 0

            # thresholding
            thresh_z = np.where((z >= -threshold) & (z <= threshold), 0, z)

            # Save results
            kwds = dict(index=index_df, columns=cols_df, struct=struct_name, resolution=res)
            _save_results(z, feat_name, **kwds)
            _save_results(thresh_z, feat_name, thresh=threshold, **kwds)

            logger_both.info(f'\t{feat_name:<15}: \t[{data_cn.shape[0]}/{orig_n_cn} controls available]')

        # Mahalanobis
        if len(available_feat) < 2:
            continue

        common_ids = reduce(np.intersect1d, [df.index for df in list_df_cn])
        for i, (df, x) in enumerate(zip(list_df_cn, list_data_cn)):
            mask = df.index.isin(common_ids)
            list_data_cn[i] = x[mask].reshape(x[mask].shape[0], -1)
        x_cn = np.stack(list_data_cn, axis=-1)
        x_px = np.stack(list_data_px, axis=-1).reshape(-1, len(list_data_px))

        md = mahalanobis_distance(x_cn, x_px)
        thresh_md = np.where((md >= -threshold) & (md <= threshold), 0, md)  # thresholding

        # Save results
        # feat_name = '&'.join(available_feat)
        kwds = dict(index=index_df, columns=cols_df, struct=struct_name, resolution=res)
        _save_results(md, available_feat, **kwds)
        _save_results(thresh_md, available_feat, thresh=threshold, **kwds)

        logger_both.info(f'\n\t{"Mahalanobis":<15}: \t[{x_cn.shape[0]}/{orig_n_cn} controls available]\n')

    logger_both.info('Done!\n\n')


if __name__ == '__main__':
    main()
