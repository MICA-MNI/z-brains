import itertools
import os
import sys
import logging
from functools import reduce
from collections import defaultdict

import numpy as np
import pandas as pd
import nibabel as nib

from functions.deconfounding import CombatModel, RegressOutModel
from functions.constants import Structure, Resolution, Feature, struct_to_folder, \
    map_feature_to_file, HIGH_RESOLUTION_CTX, LOW_RESOLUTION_CTX, HIGH_RESOLUTION_HIP, \
    LOW_RESOLUTION_HIP, FOLDER_MAPS


COLUMN_DATASET = '__zbrains_dataset_identifier__'


logger = logging.getLogger('analysis_logger')


# Helpers -----------------------------------------------------------------------------------------
def get_id(sid, add_prefix=True):
    if sid.startswith('sub-'):
        sid = sid[4:]
    if add_prefix:
        sid = f'sub-{sid}'
    return sid


def get_session(session, add_predix=True):
    if pd.isnull(session) or session == 'n/a':
        session = ''
    if session.startswith('ses-'):
        session = session[4:]
    if add_predix and session != '':
        session = f'ses-{session}'
    return session


def _get_bids_id(sid, ses):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    if ses == '':
        return sid
    return f'{sid}_{ses}'


def _get_subject_dir(root_pth, sid, ses):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    return f'{root_pth}/{sid}' if ses == '' else f'{root_pth}/{sid}/{ses}'


# Mahalanobis -------------------------------------------------------------------------------------
def get_deconfounder(*, covariates: list[str]) -> RegressOutModel | CombatModel:
    """ Build deconfounder based on covariates.

    If covariates include 'site', use ComBat. Otherwise, use RegressOutModel.

    For ComBat, all remaining covariates will be preserved (default), unless
    prefixed with '-', in which case the covariate will be removed (e.g., -age).
    The '-' prefix is ignored for 'site' and when using RegressOutModel.

    Parameters
    ----------
    covariates :
        List of covariates. Covariates prepended with '-' will be removed
        from the data when using ComBat.

    Returns
    -------
    dec :
        Deconfounder object.
    """

    covariates = [s[1:] if s == '-site' else s for s in covariates]

    if 'site' in covariates:
        remove = [k[1:] for k in covariates if k.startswith('-')]
        keep = [k for k in covariates if not k.startswith('-') and k != 'site']
        return CombatModel(site_key='site', keep=keep, remove=remove)

    cols = [s[1:] if s.startswith('-') else s for s in covariates]
    return RegressOutModel(remove=cols)


def compute_asymmetry(x_lh: np.ndarray, x_rh: np.ndarray) -> np.ndarray:
    """ Compute asymmetry.

    Parameters
    ----------
    x_lh
        Left hemisphere data. Shape (n_subjects, n_points) or (n_points,).
    x_rh
        Right hemisphere data. Shape (n_subjects, n_points) or (n_points,).

    Returns
    -------
    s
        Output data. Shape (n_subjects, n_points) or (n_points,).
    """

    den = x_lh + x_rh
    den *= .5
    x = np.divide(x_lh - x_rh, den, out=den, where=den > 0)
    return x


def zscore(x_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    """ Calculate z-scores for the test data based on the training data.

    Parameters
    ----------
    x_train
        Training data. Shape (n_train, n_points).
    x_test
        Test data. Shape (n_test, n_points) or (n_points,).

    Returns
    -------
    z
        Z-scores for the test data. Shape (n_test, n_points) or (n_points,).
    """

    mask = np.any(x_train != 0, axis=0)  # ignore all zeros
    # x_train[:, ~mask] = np.finfo(float).eps

    z = x_test - np.nanmean(x_train, axis=0)
    z = np.divide(z, np.nanstd(x_train, axis=0), out=z, where=mask)
    z[~mask] = 0

    return z


def mahalanobis_distance(x_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    """Compute mahalanobis distance.

    Parameters
    ----------
    x_train
        Training data. Shape (n_train, n_points, n_feat)
    x_test: ndarray of shape=(n_points, n_features)
        Test data. Shape (n_test, n_points, n_feat) or (n_points, n_feat)

    Returns
    -------
    dist
        Mahalanobis distance for the test data. Shape (n_test, n_points) or
        (n_points,).
    """

    n_train = x_train.shape[0]
    mu = x_train.mean(axis=0)
    x_train = x_train - mu

    cov = np.moveaxis(x_train, 0, -1) @ x_train.swapaxes(0, 1)
    cov /= n_train - 1
    cov_inv = np.linalg.inv(cov)

    x_test = x_test - mu
    dist = np.sqrt(x_test[:, None] @ cov_inv @ x_test[..., None]).squeeze()

    # n_train, n_points, n_feat = x_train.shape
    # dist = np.zeros(n_points, dtype=np.float32)
    # for i, xp in enumerate(x_test):
    #     xc = x_train[:, i]
    #     mu = xc.mean(axis=0)
    #
    #     xc = xc - mu
    #     xp = xp - mu
    #
    #     cov = xc.T @ xc
    #     cov /= n_train - 1
    #     try:
    #         cov_inv = np.linalg.inv(cov)
    #         dist[i] = np.sqrt(xp.T @ cov_inv @ xp)
    #     except np.linalg.LinAlgError:
    #         pass  # set distance to zero (default)

    return dist


# Read/write functions ----------------------------------------------------------------------------
def _map_resolution(struct: Structure, resolution: Resolution):
    if struct == 'cortex':
        return HIGH_RESOLUTION_CTX if resolution == 'high' else LOW_RESOLUTION_CTX
    if struct == 'hippocampus':
        return HIGH_RESOLUTION_HIP if resolution == 'high' else LOW_RESOLUTION_HIP
    raise ValueError(f'Mapping resolution for unknown structure: {struct}')


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
              smooth: float | None = None, asymmetry=False, raise_error: bool = True) \
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

    bids_id = _get_bids_id(sid, ses)
    subject_dir = _get_subject_dir(pth_zbrains, sid, ses)
    if not os.path.isdir(subject_dir):
        logger.warning(f"Subject '{bids_id}' zbrains directory does not exist")

    folder = f'{subject_dir}/{FOLDER_MAPS}/{struct_to_folder[struct]}'

    if struct == 'subcortex' and feat == 'thickness':
        feat = 'volume'

    if struct == 'subcortex':
        ipth = _get_ipath_from_template(struct, root_path=folder, bids_id=bids_id, feat=feat)

        try:
            # return pd.read_csv(ipth, header=[0], index_col=0)
            # return pd.read_csv(ipth, header=[0]).iloc[:, 1:]  # ignore index
            x = pd.read_csv(ipth, header=[0]).iloc[:, 1:]
        except FileNotFoundError as e:
            if raise_error:
                logger.error(e, stack_info=True, exc_info=True)

                sys.excepthook = sys.__excepthook__  # skip unhandled exception
                raise

            logger.warning(f'File not found: "{ipth}"')
            return None

        if asymmetry:
            y = x.to_numpy()
            n = x.shape[1] // 2
            x.iloc[:, :n] = compute_asymmetry(y[:, :n], y[:, n:])
            return x.iloc[:, :n]
        return x

    x = []
    for h in ['L', 'R']:
        ipth = _get_ipath_from_template(struct, root_path=folder, bids_id=bids_id, hemi=h,
                                        res=resolution, label=label, feat=feat, smooth=smooth)
        try:
            x.append(nib.load(ipth).darrays[0].data)
        except FileNotFoundError as e:
            if raise_error:
                logger.error(e, stack_info=True, exc_info=True)

                sys.excepthook = sys.__excepthook__  # skip unhandled exception
                raise

            logger.warning(f'File not found: "{ipth}"')
            return None

        if asymmetry:
            return compute_asymmetry(x[0], x[1])
        return np.concatenate(x)


def _load_data(
        pth_zbrains: str | list[str], *,
        df_subjects: pd.DataFrame | list[pd.DataFrame], struct: Structure,
        feat: Feature, resolution: Resolution | None = None,
        label: str | None = None, smooth: float | None = None, asymmetry=False) \
        -> tuple[pd.DataFrame | np.ndarray, pd.DataFrame]:
    """ Load data form all subjects in 'df_subjects'.

    Parameters
    ----------
    pth_zbrains:
        Path to the zbrains derivatives folder.
    df_subjects:
        Data frame with subjects. Must contain participant_id column. session_id col optional.
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

    if isinstance(pth_zbrains, list):
        list_data, list_dfs = [], []
        for i, (pth, df) in enumerate(zip(pth_zbrains, df_subjects)):
            x, df = _load_data(
                pth, df_subjects=df, struct=struct, feat=feat,
                resolution=resolution, label=label, smooth=smooth,
                asymmetry=asymmetry
            )

            list_data.append(x)

            df[COLUMN_DATASET] = f'Dataset{i:>03}'
            list_dfs.append(df)

        common_cols = reduce(np.intersect1d, [df.columns for df in list_dfs])
        list_dfs = [df[common_cols] for df in list_dfs]
        df = pd.concat(list_dfs, axis=0, ignore_index=True)

        if struct == 'subcortex':
            return pd.concat(list_data, axis=0, ignore_index=True), df
        return np.vstack(list_data), df

    missing_subjects = np.ones(df_subjects.shape[0], dtype=bool)
    data = []
    for i, row in df_subjects.iterrows():
        sid = row.get('participant_id')
        ses = row.get('session_id', '')

        x = _load_one(pth_zbrains, sid=sid, ses=ses, struct=struct, feat=feat,
                      resolution=resolution, label=label, smooth=smooth,
                      asymmetry=asymmetry, raise_error=False)
        if x is not None:
            data.append(x)
            missing_subjects[i] = False

    if struct == 'subcortex':
        return pd.concat(data, axis=0, ignore_index=True), df_subjects[~missing_subjects]
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
    folder = f'{pth_analysis}/{struct_to_folder[struct]}'

    # Handle the case when feat is a list of string (used for Mahalanobis)
    is_list = isinstance(feat, list)
    feat = feat if is_list else [feat]
    feat = ['volume' if (k == 'thickness' and struct == 'subcortex') else k for k in feat]
    feat = [map_feature_to_file[k] for k in feat]
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


def load_demo(path: str | list[str], *, rename: dict[str, str] | None = None,
              dtypes: dict | None = None):

    if not (is_list := isinstance(path, list)):
        path = [path]

    list_df = []
    for p in path:
        sep = '\t' if p.endswith('.tsv') else ','
        df = pd.read_csv(p, header=[0], dtype=dtypes, sep=sep)
        if rename is not None:
            df.rename(columns=rename, inplace=True)
        list_df.append(df)

    if not is_list:
        return list_df[0]
    return list_df


def _subject_zscore(
        *, dir_px, px_sid, px_ses, data_cn: np.ndarray, feat: Feature,
        deconfounder: CombatModel | RegressOutModel | None,
        df_px: pd.DataFrame | None, **kwargs):

    # Load patient data
    try:
        data_px = _load_one(dir_px, sid=px_sid, ses=px_ses, raise_error=True,
                            **kwargs)
    except FileNotFoundError:
        logger.warning(f'\t{feat:<15}: \tNo data available.')
        raise

    # For subcortex, we have dataframes
    cols_df = index_df = None
    if is_df := isinstance(data_px, pd.DataFrame):
        index_df, cols_df = data_px.index, data_px.columns
        data_px = data_px.to_numpy()

    # Deconfounding
    if deconfounder:
        data_px = deconfounder.transform(data_px, df_px)

    # Analysis: z-scoring
    z = zscore(data_cn, data_px)
    if is_df:
        z = pd.DataFrame(z.reshape(1, -1), index=index_df, columns=cols_df)

    # Store data for mahalanobis
    return dict(z=z, data_px=data_px, index_df=index_df, cols_df=cols_df)


def _subject_mahalanobis(*, data: defaultdict[str, list]):

    list_df_cn = data['data_df']
    list_data_cn = data['data_cn']
    common_ids = reduce(np.intersect1d, [df.index for df in list_df_cn])
    for i, (df, x) in enumerate(zip(list_df_cn, list_data_cn)):
        mask = df.index.isin(common_ids)
        list_data_cn[i] = x[mask]

    data_cn = np.stack(list_data_cn, axis=-1)
    data_px = np.stack(data['data_px'], axis=-1)

    md = mahalanobis_distance(data_cn, data_px)
    cols_df = data['cols_df'][0]
    index_df = data['index_df'][0]
    if index_df is not None:
        md = pd.DataFrame(md.reshape(1, -1), index=index_df, columns=cols_df)

    return dict(md=md, data_cn=data_cn, data_px=data_px)


def run_analysis(
        *, px_sid: str, px_ses='SINGLE', dir_cn: list[str],
        demo_cn: list[str], dir_px: str, demo_px: str | None = None,
        structures: list[Structure], features: list[Feature],
        cov_normative: list[str] | None = None,
        cov_deconfound: list[str] | None = None, smooth_ctx: float,
        smooth_hip: float, resolutions: list[Resolution], labels_ctx: list[str],
        labels_hip: list[str], actual_to_expected: dict[str, str],
        asymmetry=False, analysis_folder: str
):

    logger.debug(f'Logging call: {sys.argv[0]} {" ".join(sys.argv[1:])}')
    logger.info(f'{"Asymmetry" if asymmetry else "Regional"} analysis\n')

    pth_analysis = f'{_get_subject_dir(dir_px, px_sid, px_ses)}/{analysis_folder}'

    # Load dataframes ----------------------------------------------------------
    # Read participant_id, session_id, sex, site as categorical if available
    dtypes = {actual: (str if expected != 'age' else float)
              for actual, expected in actual_to_expected.items()}

    # Controls
    list_df_cn = load_demo(demo_cn, rename=actual_to_expected, dtypes=dtypes)
    n_cn = sum(len(df) for df in list_df_cn)

    # Patient
    df_px = None
    if demo_px:
        df_px = load_demo(demo_px, rename=actual_to_expected, dtypes=dtypes)

    # Load dataframes ----------------------------------------------------------
    iterables = []
    for st in structures:
        if st == 'subcortex':
            iterables.append((st, None, None))
        elif st == 'cortex':
            iterables += itertools.product([st], resolutions, labels_ctx)
        else:
            iterables += itertools.product([st], resolutions, labels_hip)

    # Main loop ----------------------------------------------------------------
    for struct, resol, label in iterables:
        s = f'[resolution = {resol:<4}\tlabel = {label:<15}]' \
            if struct != 'subcortex' else ''
        logger.info(f'\nStructure: {struct} {s}')

        # Shared kwds
        smooth = smooth_ctx if struct == 'cortex' else smooth_hip
        kwds = dict(struct=struct, resolution=resol, label=label, smooth=smooth,
                    asymmetry=asymmetry)

        data_mahalanobis = defaultdict(list)
        for feat in features:

            # Load control data
            kwds |= dict(feat=feat)
            data_cn, df_cn = _load_data(dir_cn, df_subjects=list_df_cn, **kwds)
            if isinstance(data_cn, pd.DataFrame):
                data_cn = data_cn.to_numpy()

            # Deconfounding
            dec = None
            if cov_deconfound is not None and df_px is not None:
                dec = get_deconfounder(covariates=cov_deconfound)
                data_cn = dec.fit_transform(data_cn, df_cn)

            logger.info(f'\t{feat:<15}: \t[{df_cn.shape[0]}/{n_cn} '
                        f'controls available]')

            # zscore
            res = _subject_zscore(
                dir_px=dir_px, px_sid=px_sid, px_ses=px_ses, data_cn=data_cn,
                deconfounder=dec, df_px=df_px, **kwds)

            # Save results
            z = res['z'] if not asymmetry else res['z'].reshape(2, -1)
            _save(pth_analysis, x=z, sid=px_sid, ses=px_ses, **kwds)

            res |= dict(data_cn=data_cn, df_cn=df_cn, feat=feat)
            for k, v in res.items():
                data_mahalanobis[k].append(v)

        # Mahalanobis
        if len(data_mahalanobis['feat']) < 2:
            continue

        # Analysis: mahalanobis distance
        res = _subject_mahalanobis(data=data_mahalanobis)

        n_available_cn = res['data_cn'].shape[0]
        logger.info(f'\n\t{"Mahalanobis":<15}: \t[{n_available_cn}/{n_cn} '
                    f'controls available]\n')

        # Save results
        md = res['md'] if not asymmetry else res['md'].reshape(2, -1)

        kwds |= dict(feat=data_mahalanobis['feat'])
        _save(pth_analysis, x=md, sid=px_sid, ses=px_ses, **kwds)

    logger.info('Done!\n\n')
