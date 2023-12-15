import itertools
import os
import sys
import logging
from functools import reduce
from typing import get_args, Literal

import numpy as np
import pandas as pd
import nibabel as nib


from deconfounding import CombatModel, RegressOutModel


# Zbrains config file -----------------------------------------------------------------------------
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


CONFIG_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           'config.cfg')
CONFIG = _read_config(CONFIG_FILE)


# Global variables and types  ---------------------------------------------------------------------
map_feature = dict(thickness='thickness', volume='volume', flair='flair', ADC='ADC', FA='FA',
                   qT1='T1map')

Structure = Literal['cortex', 'hippocampus', 'subcortex']
LIST_STRUCTURES: list[Structure] = list(get_args(Structure))

map_struct_to_folder = dict(cortex=CONFIG['FOLDER_CTX'], subcortex=CONFIG['FOLDER_SCTX'],
                            hippocampus=CONFIG['FOLDER_HIP'])

Feature = Literal['flair', 'ADC', 'FA', 'qT1', 'thickness']
LIST_FEATURES: list[Feature] = list(get_args(Feature))

Analysis = Literal['norm-z', 'norm-modelling']
LIST_ANALYSES: list[Analysis] = list(get_args(Analysis))

Resolution = Literal['low', 'high']
LIST_RESOLUTIONS: list[Resolution] = list(get_args(Resolution))


ALIAS_COLUMNS = dict(ID='ID', SES='SES', AGE='AGE', SEX='SEX', SITE='SITE')


# Helpers -----------------------------------------------------------------------------------------
def get_id(sid, add_prefix=True):
    if sid.startswith('sub-'):
        sid = sid[4:]
    if add_prefix:
        sid = f'sub-{sid}'
    return sid


def get_session(session, add_predix=True):
    if pd.isnull(session) or session == '':
        session = 'SINGLE'
    if session.startswith('ses-'):
        session = session[4:]
    if add_predix and session != 'SINGLE':
        session = f'ses-{session}'
    return session


def _get_bids_id(sid, ses):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    if ses == 'SINGLE':
        return sid
    return f'{sid}_{ses}'


def _get_subject_dir(root_pth, sid, ses):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    return f'{root_pth}/{sid}' if ses == 'SINGLE' else f'{root_pth}/{sid}/{ses}'


# Mahalanobis -------------------------------------------------------------------------------------
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


# Read/write functions ----------------------------------------------------------------------------
def _map_resolution(struct: Structure, resolution: Resolution):
    if struct == 'cortex':
        return CONFIG['HIGH_RESOLUTION_CTX' if resolution == 'high' else 'LOW_RESOLUTION_CTX']
    return CONFIG['HIGH_RESOLUTION_HIP' if resolution == 'high' else 'LOW_RESOLUTION_HIP']


def _get_ipath_from_template(struct, **kwargs):
    if struct == 'subcortex':
        ipth = '{root_path}/{bids_id}_feat-{feat}.csv'

    elif struct == 'cortex':
        ipth = ('{root_path}/{bids_id}_hemi-{hemi}_surf-fsLR-{res}_label-{label}_feat-{feat}'
                '_smooth-{smooth}mm.func.gii')

    else:
        ipth = ('{root_path}/{bids_id}_hemi-{hemi}_den-{res}_label-{label}_feat-{feat}'
                '_smooth-{smooth}mm.func.gii')

    if 'res' in kwargs:
        kwargs['res'] = _map_resolution(struct, kwargs['res'])

    return ipth.format(**kwargs)


def _get_opath_from_template(struct, **kwargs):
    if struct == 'subcortex':
        opth = '{root_path}/{bids_id}_feature-{feat}_analysis-{analysis}.csv'

    elif struct == 'cortex':
        opth = ('{root_path}/{bids_id}_hemi-{hemi}_surf-fsLR-{res}_label-{label}_feature-{feat}'
                '_smooth-{smooth}mm_analysis-{analysis}.func.gii')
    else:
        opth = ('{root_path}/{bids_id}_hemi-{hemi}_den-{res}_label-{label}_feature-{feat}'
                '_smooth-{smooth}mm_analysis-{analysis}.func.gii')

    if 'res' in kwargs:
        kwargs['res'] = _map_resolution(struct, kwargs['res'])
    return opth.format(**kwargs)


def _load_one(pth_zbrains: str, *, sid: str, ses: str, struct: Structure, feat: Feature,
              resolution: Resolution | None = None, label: str | None = None,
              smooth: float | None = None, raise_error: bool = True) \
        -> np.ndarray | pd.DataFrame | None:
    """ Load subject data

    Parameters
    ----------
    pth_zbrains:
        Path to the zbrains derivatives folder.
    sid:
        Subject id.
    ses:
        Session identifier.
    struct:
        Structure.
    feat:
        Feature name.
    resolution:
        Resolution. Required when struct='cortex' or struct='hippocampus'.
    label:
        Label indicates the surfaces used in the volume to surface mapping. Required when
        struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or struct='hippocampus'.
    raise_error:
        Raise error if file not found

    Returns
    -------
    x:
        Subject data. If structure is 'cortex' or 'hippocampus', shape ndarray of shape
        (2 * n_vertices_per_hemisphere,). Otherwise, return DataFrame of shape
        (1, 2 * n_subcortical_structures_per_hemisphere).
        None if no data available for at least one hemisphere.
    """

    logger_file = logging.getLogger('analysis_logger.file')
    logger_both = logging.getLogger('analysis_logger.both')

    bids_id = _get_bids_id(sid, ses)
    subject_dir = _get_subject_dir(pth_zbrains, sid, ses)
    if not os.path.isdir(subject_dir):
        logger_both.warning(f"Subject '{bids_id}' zbrains directory does not exist")

    folder = f'{subject_dir}/{CONFIG["FOLDER_MAPS"]}/{map_struct_to_folder[struct]}'

    if struct == 'subcortex' and feat == 'thickness':
        feat = 'volume'

    if struct == 'subcortex':
        ipth = _get_ipath_from_template(struct, root_path=folder, bids_id=bids_id, feat=feat)

        try:
            return pd.read_csv(ipth, header=[0], index_col=0)
        except FileNotFoundError as e:
            if raise_error:
                logger_file.error(e, stack_info=True, exc_info=True)

                sys.excepthook = sys.__excepthook__  # skip unhandled exception
                raise

            logger_file.warning(f'File not found: "{ipth}"')
            return None

    x = []
    for h in ['L', 'R']:
        ipth = _get_ipath_from_template(struct, root_path=folder, bids_id=bids_id, hemi=h,
                                        res=resolution, label=label, feat=feat, smooth=smooth)
        try:
            x.append(nib.load(ipth).darrays[0].data)
        except FileNotFoundError as e:
            if raise_error:
                logger_file.error(e, stack_info=True, exc_info=True)

                sys.excepthook = sys.__excepthook__  # skip unhandled exception
                raise

            logger_file.warning(f'File not found: "{ipth}"')
            return None

        return np.concatenate(x)


def _load_data(pth_zbrains: str, *, df_subjects: pd.DataFrame, struct: Structure, feat: Feature,
               resolution: Resolution | None = None, label: str | None = None, smooth: float) \
        -> tuple[pd.DataFrame | np.ndarray, pd.DataFrame]:
    """ Load data form all subjects in 'df_subjects'.

    Parameters
    ----------
    pth_zbrains:
        Path to the zbrains derivatives folder.
    df_subjects:
        Data frame with subjects. Must contain ID column. SES column is optional.
    struct:
        Structure.
    feat:
        Feature name.
    resolution:
        Resolution. Required when struct='cortex' or struct='hippocampus'.
    label:
        Label indicates the surfaces used in the volume to surface mapping. Required when
        struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or struct='hippocampus'.

    Returns
    -------
    x:
        Data for CN. Return ndarray of shape (n_available_subjects, 2 * n_points_per_hemisphere).
        If struct='subcortex', return DataFrame.
    df_controls_available:
        Dataframe of shape (n_available_subjects, n_cols), only including those rows in
        'df_controls' with available data.
    """

    missing_subjects = np.ones(df_subjects.shape[0], dtype=bool)
    data = []
    for i, row in df_subjects.iterrows():
        sid = row.get('ID')
        ses = row.get('SES', 'SINGLE')

        x = _load_one(pth_zbrains, sid=sid, ses=ses, struct=struct, feat=feat,
                      resolution=resolution, label=label, smooth=smooth, raise_error=False)
        # x = _load_one(sid, ses, feat, struct=struct, resolution=resolution, raise_error=False)
        if x is not None:
            data.append(x)
            missing_subjects[i] = False

    if struct == 'subcortex':
        return pd.concat(data, axis=0), df_subjects[~missing_subjects]
    return np.stack(data, axis=0), df_subjects[~missing_subjects]


def _save(pth_analysis: str, *, x: np.ndarray | pd.DataFrame, sid: str, struct: Structure,
          feat: Feature | list[Feature], ses: str = 'SINGLE', resolution: Resolution | None = None,
          label: str | None = None, smooth: float | None = None, asymmetry=False):
    """ Save results

    Parameters
    ----------
    pth_analysis:
        Path to the analysis folder in the zbrains derivatives folder.
    x:
        Patient data to save. shape ndarray of shape (2, n_points) or (n_points,)
    sid:
        Subject id.
    struct:
        Structure.
    feat:
        One feature or list of features if Mahalanobis.
    ses:
        Session identifier.
    resolution:
        Resolution. Required when struct='cortex' or struct='hippocampus'.
    label:
        Label indicates the surfaces used in the volume to surface mapping. Required when
        struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or struct='hippocampus'.
    asymmetry:
        If true, data is from asymmetry analysis, so only left hemisphere.
    """

    bids_id = _get_bids_id(sid, ses)
    folder = f'{pth_analysis}/{map_struct_to_folder[struct]}'

    # Handle the case when feat is a list of string (used for Mahalanobis)
    is_list = isinstance(feat, list)
    feat = feat if is_list else [feat]
    feat = ['volume' if (k == 'thickness' and struct == 'subcortex') else k for k in feat]
    feat = [map_feature[k] for k in feat]
    feat = '-'.join(feat) if is_list else feat[0]

    analysis = 'asymmetry' if asymmetry else 'regional'

    if struct == 'subcortex':
        opth = _get_opath_from_template(struct, root_path=folder, bids_id=bids_id, feat=feat,
                                        analysis=analysis)
        x.to_csv(opth)
        return

    for i, h in enumerate(['L', 'R']):
        data_array = nib.gifti.GiftiDataArray(data=x if asymmetry else x[i])  # per hemisphere
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)

        opth = _get_opath_from_template(struct, root_path=folder, bids_id=bids_id, hemi=h,
                                        res=resolution, label=label, feat=feat, smooth=smooth,
                                        analysis=analysis)
        nib.save(image, opth)

        if asymmetry:
            break


def _save_results(pth_analysis: str, *, x: np.ndarray | pd.DataFrame, sid: str, struct: Structure,
                  feat: Feature | list[Feature], ses: str = 'SINGLE',
                  resolution: Resolution | None = None, label: str | None = None,
                  smooth: float | None = None, asymmetry=False, index: list[str] | None = None,
                  columns: list[str] | None = None):

    if columns is not None:
        x = pd.DataFrame(x.reshape(1, -1), index=index, columns=columns)

    _save(pth_analysis, x=x, sid=sid, struct=struct, feat=feat, ses=ses, resolution=resolution,
          label=label, smooth=smooth, asymmetry=asymmetry)


# Prepare data ------------------------------------------------------------------------------------
def deconfound(*, covariates: list[str], data_cn: np.ndarray, df_cn: pd.DataFrame,
               data_px: np.ndarray, df_px: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:

    cols = [k for k in covariates if k != 'SITE']
    if 'SITE' in covariates:
        deconfounder = CombatModel(site_key='SITE', remove=cols)
    else:
        deconfounder = RegressOutModel(remove=cols)

    data_cn = deconfounder.fit_transform(data_cn, df_cn)
    data_px = deconfounder.transform(data_px, df_px)

    return data_cn, data_px


def make_asymmetric(*, data_cn: np.ndarray, data_px: np.ndarray):
    data_cn = data_cn.reshape(data_cn.shape[0], 2, -1)
    lh_data_cn, rh_data_cn = data_cn[:, 0], data_cn[:, 1]
    den_cn = .5 * (lh_data_cn + rh_data_cn)
    data_cn = np.divide(lh_data_cn - rh_data_cn, den_cn, out=den_cn, where=den_cn > 0)

    data_px = data_px.reshape(2, -1)
    lh_data_px, rh_data_px = data_px[0], data_px[1]
    den_px = .5 * (lh_data_px + rh_data_px)
    data_px = np.divide(lh_data_px - rh_data_px, den_px, out=den_px, where=den_px > 0)

    return data_cn, data_px


def zscore(*, data_cn: np.ndarray, data_px: np.ndarray):
    mask = np.any(data_cn != 0, axis=0)  # ignore all zeros
    data_cn[:, ~mask] = np.finfo(float).eps

    z = data_px - np.nanmean(data_cn, axis=0)
    z = np.divide(z, np.nanstd(data_cn, axis=0), out=z, where=mask)
    z[~mask] = 0

    return z


##################################################################################################
# Analysis
##################################################################################################
def run_analysis(pth_zbrains: str, *, analysis_folder: str, px_sid: str, path_demo_cn: str,
                 structures: list[Structure], features: list[Feature],
                 resolutions: list[Resolution], labels_ctx: list[str], labels_hip: list[str],
                 smooth_ctx: float, smooth_hip: float, px_ses: str = 'SINGLE',
                 path_demo_px: str | None = None, asymmetry=False,
                 col_aliases: dict[str, str] | None = None,
                 cov_deconfound: list[str] | None = None):

    logger_file = logging.getLogger('analysis_logger.file')
    logger_both = logging.getLogger('analysis_logger.both')

    logger_file.info(f'Logging call: {sys.argv[0]} {" ".join(sys.argv[1:])}')
    logger_both.debug(f'{"Asymmetry" if asymmetry else "Regional"} analysis\n')

    pth_analysis = f'{_get_subject_dir(pth_zbrains, px_sid, px_ses)}/{analysis_folder}'

    # Data frame of control subjects
    # All columns in ALIAS_COLUMNS are categorical except AGE
    dtypes = {v: (str if k != 'AGE' else float) for k, v in ALIAS_COLUMNS.items()}
    orig_df_cn = pd.read_csv(path_demo_cn, header=[0], dtype=dtypes)
    if col_aliases is not None:  # Rename columns to expected names
        orig_df_cn.rename(columns=col_aliases, inplace=True)
    orig_n_cn = orig_df_cn.shape[0]

    df_px = None
    if path_demo_px is not None:
        df_px = pd.read_csv(path_demo_px, header=[0], dtype=dtypes)
        if col_aliases is not None:
            df_px.rename(columns=col_aliases, inplace=True)

    iterables = []
    for st in structures:
        if st == 'subcortex':
            iterables.append((st, None, None))
        elif st == 'cortex':
            iterables += itertools.product([st], resolutions, labels_ctx)
        else:
            iterables += itertools.product([st], resolutions, labels_hip)

    # Main loop
    for struct, resol, label in iterables:
        s = f'[resolution = {resol:<4}\tlabel = {label:<15}]' if struct != 'subcortex' else ''
        logger_both.info(f'\nStructure: {struct} {s}')

        smooth = smooth_ctx if struct == 'cortex' else smooth_hip

        index_df, cols_df = None, None

        # Collect data for Mahalanobis
        list_df_cn = []  # Because some subjects may not have some features
        list_data_cn = []
        list_data_px = []
        feat_mahalanobis = []

        for feat in features:
            # Load data
            try:
                data_cn, df_cn = _load_data(pth_zbrains, df_subjects=orig_df_cn, struct=struct,
                                            feat=feat, resolution=resol, label=label, smooth=smooth)
                data_px = _load_one(pth_zbrains, sid=px_sid, ses=px_ses, struct=struct, feat=feat,
                                    resolution=resol, label=label, smooth=smooth, raise_error=True)
            except FileNotFoundError:  # When patient files not available
                logger_both.info(f'\t{feat:<15}: \tNo data available.')
                continue

            # For subcortex, we have dataframes
            if isinstance(data_px, pd.DataFrame):
                index_df, cols_df = data_px.index, data_px.columns
                data_cn = data_cn.to_numpy()
                data_px = data_px.to_numpy()

            # Deconfounding
            if cov_deconfound is not None and df_px is not None:
                data_cn, data_px = deconfound(covariates=cov_deconfound, data_cn=data_cn,
                                              df_cn=df_cn, data_px=data_px, df_px=df_px)

            # Use asymmetry analysis
            if asymmetry:
                data_cn, data_px = make_asymmetric(data_cn=data_cn, data_px=data_px)

                if cols_df is not None:
                    cols_df = cols_df[:cols_df.size // 2]

            # Store data for mahalanobis
            feat_mahalanobis.append(feat)
            list_data_cn.append(data_cn)
            list_df_cn.append(df_cn)
            list_data_px.append(data_px)

            # Analysis: z-scoring
            z = zscore(data_cn=data_cn, data_px=data_px)
            z = z if asymmetry else z.reshape(2, -1)

            # Save results
            _save_results(pth_analysis, x=z, sid=px_sid, struct=struct, feat=feat, ses=px_ses,
                          resolution=resol, label=label, smooth=smooth, asymmetry=asymmetry,
                          index=index_df, columns=cols_df)

            logger_both.info(f'\t{feat:<15}: \t[{data_cn.shape[0]}/{orig_n_cn} controls available]')

        # Mahalanobis
        if len(feat_mahalanobis) < 2:
            continue

        common_ids = reduce(np.intersect1d, [df.index for df in list_df_cn])
        for i, (df, x) in enumerate(zip(list_df_cn, list_data_cn)):
            mask = df.index.isin(common_ids)
            list_data_cn[i] = x[mask]  # .reshape(x[mask].shape[0], -1)
        x_cn = np.stack(list_data_cn, axis=-1)
        x_px = np.stack(list_data_px, axis=-1)  # .reshape(-1, len(list_data_px))

        md = mahalanobis_distance(x_cn, x_px)
        md = md if asymmetry else md.reshape(2, -1)

        # Save results
        _save_results(pth_analysis, x=md, sid=px_sid, struct=struct, feat=feat_mahalanobis,
                      ses=px_ses, resolution=resol, label=label, smooth=smooth, asymmetry=asymmetry,
                      index=index_df, columns=cols_df)

        logger_both.info(f'\n\t{"Mahalanobis":<15}: \t[{x_cn.shape[0]}/{orig_n_cn} '
                         f'controls available]\n')

    logger_both.info('Done!\n\n')



