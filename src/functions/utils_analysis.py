import os
import re
import sys
import logging
import itertools
from pathlib import Path
from typing import Union, Optional, List, Dict, Tuple, Any, DefaultDict
from functools import reduce
from collections import defaultdict
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import nibabel as nib
import copy
from .utilities import (
    add_field_to_xml,
    replace_field_in_xml,
    edit_root_attributes_in_xml,
    remove_doctype_from_xml,
)
from .constants import ProcessingException


from .deconfounding import CombatModel, RegressOutModel
from .constants import (
    Structure,
    Analysis,
    Approach,
    Resolution,
    Feature,
    struct_to_folder,
    approach_to_folder,
    map_feature_to_file,
    HIGH_RESOLUTION_CTX,
    LOW_RESOLUTION_CTX,
    HIGH_RESOLUTION_HIP,
    LOW_RESOLUTION_HIP,
    FOLDER_MAPS,
    LIST_FEATURES,
)

COLUMN_DATASET = "__zbrains_dataset_identifier__"

PathType = Union[str, os.PathLike]

# logger = logging.getLogger('analysis_logger')


# Helpers -----------------------------------------------------------------------------------------
def get_id(sid: str, add_prefix=True):
    if sid.startswith("sub-"):
        sid = sid[4:]
    if add_prefix:
        sid = f"sub-{sid}"
    return sid


def get_session(session, add_predix=True):
    if pd.isnull(session) or session == "n/a" or session == "":
        return None
    if session.startswith("ses-"):
        session = session[4:]
    if add_predix:
        session = f"ses-{session}"
    return session


def get_bids_id(sid: str, ses: Union[str, None] = None):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    if ses is None:
        return sid
    return f"{sid}_{ses}"


def get_subject_dir(root_pth: PathType, sid: str, ses: Union[str, None] = None):
    sid = get_id(sid, add_prefix=True)
    ses = get_session(ses, add_predix=True)
    p = f"{root_pth}/{sid}" if ses is None else f"{root_pth}/{sid}/{ses}"
    return Path(p)


# Mahalanobis -------------------------------------------------------------------------------------
def get_deconfounder(
    *, covariates: Optional[List[str]]
) -> Union[RegressOutModel, CombatModel]:
    """Build deconfounder based on covariates.

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

    covariates = [s[1:] if s == "-site" else s for s in covariates]

    if "site" in covariates:
        remove = [k[1:] for k in covariates if k.startswith("-")]
        keep = [k for k in covariates if not k.startswith("-") and k != "site"]
        return CombatModel(site_key="site", keep=keep, remove=remove)

    cols = [s[1:] if s.startswith("-") else s for s in covariates]
    return RegressOutModel(remove=cols)


def compute_asymmetry(x_lh: np.ndarray, x_rh: np.ndarray) -> np.ndarray:
    """Compute asymmetry.

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
    den *= 0.5
    return np.divide(x_lh - x_rh, den, out=den, where=den > 0)


def zscore(x_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    """Calculate z-scores for the test data based on the training data.

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
    z[..., ~mask] = 0
    if len(z.shape) > 1:
        if z.shape[1] != 1:
            z = np.mean(z, axis=1)
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
    try:
        cov_inv = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        raise ProcessingException(
            "Singular matrix, unable to compute Mahalanobis distance. Check your smoothing kernels, they might be too small."
        )

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
def map_resolution(struct: Structure, res: Resolution):
    if struct == "cortex":
        return HIGH_RESOLUTION_CTX if res == "high" else LOW_RESOLUTION_CTX
    if struct == "hippocampus":
        return HIGH_RESOLUTION_HIP if res == "high" else LOW_RESOLUTION_HIP
    raise ValueError(f"Mapping resolution for unknown structure: {struct}")


def get_feature_path_from_template(struct: Structure, **kwargs) -> Path:
    if struct == "subcortex":
        ipth = "{root_path}/{bids_id}_feature-{feat}.csv"

    elif struct == "cortex":
        ipth = (
            "{root_path}/{bids_id}_hemi-{hemi}_surf-fsLR-{res}_"
            "label-{label}_feature-{feat}_smooth-{smooth}mm.func.gii"
        )

    else:
        ipth = (
            "{root_path}/{bids_id}_hemi-{hemi}_den-{res}_"
            "label-{label}_feature-{feat}_smooth-{smooth}mm.func.gii"
        )

    if "res" in kwargs:
        kwargs["res"] = map_resolution(struct, kwargs["res"])

    return Path(ipth.format(**kwargs))


def get_analysis_path_from_template(struct: Structure, **kwargs) -> Path:
    kwargs["feat"] = (
        kwargs["feat"].replace("qT1", "T1map").replace("thickness", "volume")
    )
    # print(kwargs['feat'])
    # if kwargs['feat'] == 'qT1':
    #     kwargs['feat'] = "T1map"
    # if isinstance(kwargs['feat'], list):
    #     kwargs['feat'].replace('qT1', 'T1map')
    if struct == "subcortex":
        opth = "{root_path}/{bids_id}_feature-{feat}_analysis-{analysis}.csv"

    elif struct == "cortex":
        opth = (
            "{root_path}/{bids_id}_hemi-{hemi}_surf-fsLR-{res}_"
            "label-{label}_feature-{feat}_smooth-{smooth}mm_"
            "analysis-{analysis}.func.gii"
        )
    else:
        opth = (
            "{root_path}/{bids_id}_hemi-{hemi}_den-{res}_label-{label}_"
            "feature-{feat}_smooth-{smooth}mm_analysis-{analysis}.func.gii"
        )

    if "res" in kwargs:
        kwargs["res"] = map_resolution(struct, kwargs["res"])
    # if kwargs['feat'] == 'T1map':
    #     kwargs['feat'] = 'qT1'
    # if isinstance(kwargs['feat'], list):
    #     kwargs['feat'].replace('T1map', 'qT1')
    # kwargs['feat'] = kwargs['feat'].replace('T1map', 'qT1')
    return Path(opth.format(**kwargs))


def _load_one(
    pth_zbrains: PathType,
    *,
    sid: str,
    ses: str,
    struct: Structure,
    feat: str,
    resolution: Union[str, None] = None,
    label: Union[str, None] = None,
    smooth: Union[float, None] = None,
    # analysis: Analysis,
    raise_error: bool = True,
    tmp: str = True,
) -> Union[np.ndarray, pd.DataFrame, None]:
    """Load subject data

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
        Label indicates the surfaces used in the volume to surface mapping.
        Required when struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or
         struct='hippocampus'.
    raise_error:
        Raise error if file not found

    Returns
    -------
    x:
        Subject data. If structure is 'cortex' or 'hippocampus', ndarray of
        shape (2 * n_vertices_per_hemisphere,). Otherwise, return DataFrame of
        shape (1, 2 * n_subcortical_structures_per_hemisphere).
        None if no data available for at least one hemisphere.
    """
    logger = logging.getLogger(tmp)
    feat = feat.replace("qT1", "T1map")
    bids_id = get_bids_id(sid, ses)
    subject_dir = get_subject_dir(pth_zbrains, sid, ses)
    if not os.path.isdir(subject_dir):
        logger.debug(f"Subject '{bids_id}' zbrains directory does not exist")

    folder = f"{subject_dir}/{FOLDER_MAPS}/{struct_to_folder[struct]}"

    if struct == "subcortex" and feat == "thickness":
        feat = "volume"

    kwds = dict(root_path=folder, bids_id=bids_id, feat=feat)

    if struct == "subcortex":
        ipth = get_feature_path_from_template(struct, **kwds)

        try:
            x = pd.read_csv(ipth, header=[0], index_col=0)
        except FileNotFoundError:
            if raise_error:
                raise

            logger.debug(f'File not found: "{ipth}"')
            return None

        # if analysis == 'asymmetry':
        #     y = x.to_numpy()
        #     n = x.shape[1] // 2
        #     x.iloc[:, :n] = compute_asymmetry(y[:, :n], y[:, n:])
        #     return x.iloc[:, :n]
        return x

    x = []
    for h in ["L", "R"]:
        ipth = get_feature_path_from_template(
            struct, hemi=h, res=resolution, label=label, smooth=smooth, **kwds
        )
        try:
            nbimage = nib.load(ipth)
            # arrsize = (len(nbimage.darrays[0].data), len(nbimage.darrays))
            nbarrs = np.array([i.data for i in nbimage.darrays]).T
            x.append(nbarrs)
        except FileNotFoundError:
            if raise_error:
                raise

            logger.debug(f'File not found: "{ipth}"')
            return None

    # if analysis == 'asymmetry':
    #     return compute_asymmetry(x[0], x[1])
    return np.concatenate(x)


def _load_data(
    pth_zbrains: Union[str, Optional[List[str]]],
    *,
    struct: Structure,
    feat: Feature,
    df_subjects: Union[pd.DataFrame, Optional[List[pd.DataFrame]]],
    # analysis: Analysis,
    resolution: Union[Resolution, None] = None,
    label: Union[str, None] = None,
    smooth: Union[float, None] = None,
    tmp: str,
) -> Tuple[Union[pd.DataFrame, None, np.ndarray, pd.DataFrame, None]]:
    """Load data form all subjects in 'df_subjects'.

    Parameters
    ----------
    pth_zbrains:
        Path to the zbrains derivatives folder.
    df_subjects:
        Data frame with subjects. Must contain participant_id column.
        session_id col optional.
    struct:
        Structure.
    feat:
        Feature name.
    resolution:
        Resolution. Required when struct='cortex' or struct='hippocampus'.
    label:
        Label indicates the surfaces used in the volume to surface mapping.
        Required when struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or
        struct='hippocampus'.

    Returns
    -------
    x:
        Data for CN. Return ndarray of shape (n_available_subjects,
        2 * n_points_per_hemisphere). If struct='subcortex', return DataFrame.
    df_controls_available:
        Dataframe of shape (n_available_subjects, n_cols), only including those
        rows in 'df_controls' with available data.
    """

    if isinstance(pth_zbrains, list):
        list_data, list_dfs = [], []
        for i, (pth, df) in enumerate(zip(pth_zbrains, df_subjects)):
            x, df = _load_data(
                pth,
                df_subjects=df,
                struct=struct,
                feat=feat,
                resolution=resolution,
                label=label,
                smooth=smooth,
                tmp=tmp,
                # analysis=analysis
            )

            if x is not None:
                list_data.append(x)

                df[COLUMN_DATASET] = f"Dataset{i:>03}"
                list_dfs.append(df)

        if len(list_data) == 0:
            return None, None

        common_cols = reduce(np.intersect1d, [df.columns for df in list_dfs])
        list_dfs = [df[common_cols] for df in list_dfs]
        df = pd.concat(list_dfs, axis=0, ignore_index=True)

        if struct == "subcortex":
            return pd.concat(list_data, axis=0, ignore_index=True), df
        return np.vstack(list_data), df

    missing_subjects = np.ones(df_subjects.shape[0], dtype=bool)
    data = []
    for i, row in df_subjects.iterrows():
        sid = row.get("participant_id")
        ses = row.get("session_id", None)

        x = _load_one(
            pth_zbrains,
            sid=sid,
            ses=ses,
            struct=struct,
            feat=feat,
            resolution=resolution,
            label=label,
            smooth=smooth,
            # analysis=analysis,
            raise_error=False,
            tmp=tmp,
        )
        if x is not None:
            data.append(x)
            missing_subjects[i] = False

    if missing_subjects.all():
        return None, None

    df_subjects = df_subjects[~missing_subjects].copy()
    if struct == "subcortex":
        return pd.concat(data, axis=0, ignore_index=True), df_subjects

    # data = [np.where(x == 0, 0.0000000001, x) for x in data]
    # zeros = np.where(data_mahalanobis['data_cn'][3][1] == 0)[0]
    return np.stack(data, axis=0), df_subjects


def _save(
    pth_analysis: str,
    *,
    x: Union[np.ndarray, pd.DataFrame],
    sid: str,
    struct: Structure,
    feat: Union[Feature, Optional[List[Feature]]],
    ses: str = None,
    resolution: Union[Resolution, None] = None,
    label: Union[str, None] = None,
    smooth: Union[float, None] = None,
    analysis: Analysis,
):
    """Save results

    Parameters
    ----------
    pth_analysis:
        Path to the analysis folder in the zbrains derivatives folder.
    x:
        Patient data to save. shape ndarray of shape (2, n_points)
        or (n_points,)
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
        Label indicates the surfaces used in the volume to surface mapping.
        Required when struct='cortex' or struct='hippocampus'.
    smooth:
        Size of gaussian smoothing kernel. Required when struct='cortex' or
        struct='hippocampus'.
    analysis:
        If 'asymmetry', only save left hemisphere.
    """

    bids_id = get_bids_id(sid, ses)
    folder = f"{pth_analysis}/{struct_to_folder[struct]}"

    # Handle the case when feat is a list of string (used for Mahalanobis)
    is_list = isinstance(feat, list)
    feat = feat if is_list else [feat]
    feat = [
        "volume" if (k == "thickness" and struct == "subcortex") else k for k in feat
    ]
    feat = [map_feature_to_file[k] for k in feat]
    if is_list:
        feat = [z for z in feat if "blur" not in z]
        feat = "-".join(feat)
    else:
        feat = feat[0]

    kwds = dict(root_path=folder, bids_id=bids_id, feat=feat, analysis=analysis)
    if struct == "subcortex":
        opth = get_analysis_path_from_template(struct, **kwds)
        x.to_csv(opth)
        return

    for i, h in enumerate(["L", "R"]):
        data = x if analysis == "asymmetry" else x[i]

        metadata = nib.gifti.GiftiMetaData(
            {
                "AnatomicalStructurePrimary": (
                    "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT"
                ),
            }
        )
        header = """<GIFTI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
            xsi:noNamespaceSchemaLocation="http://brainvis.wustl.edu/caret6/xml_schemas/GIFTI_Caret.xsd"
            Version="1"
            NumberOfDataArrays="1">
                """
        data_array = nib.gifti.GiftiDataArray(
            data=data, meta=metadata
        )  # per hemisphere
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)

        opth = get_analysis_path_from_template(
            struct, hemi=h, res=resolution, label=label, smooth=smooth, **kwds
        )
        nib.save(image, opth)
        # edit_root_attributes_in_xml(
        #     opth,
        #     {
        #         "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
        #         "xsi:noNamespaceSchemaLocation": "http://brainvis.wustl.edu/caret6/xml_schemas/GIFTI_Caret.xsd",
        #         "Version": "1",
        #         "NumberOfDataArrays": "1",
        #     },
        # )
        # remove_doctype_from_xml(opth)

        if analysis == "asymmetry":
            break


def load_demo(
    path: Union[PathType, Optional[List[PathType]]],
    *,
    rename: Union[Optional[Dict[str, str]], None] = None,
    dtypes: Union[Optional[Dict[str, type]], None] = None,
    tmp: str,
):
    logger = logging.getLogger(tmp)
    is_list = isinstance(path, list)
    if not is_list:
        path = [path]
    path = [Path(p) for p in path]

    list_df = []
    for p in path:
        sep = "\t" if p.suffix == ".tsv" else ","
        df = pd.read_csv(p, header=[0], dtype=dtypes, sep=sep)
        if rename is not None:
            df.rename(columns=rename, inplace=True)

            if "participant_id" in df:
                pids = df["participant_id"].tolist()
                for v in pids:
                    if not re.match(r"^sub-.+", v):
                        msg = (
                            f'Participant ID must have the form "sub-XXXX".'
                            f"\nCheck demographics file: {p}"
                        )
                        logger.error(msg)
                        exit(1)

            if "session_id" in df:
                pids = df["session_id"].tolist()
                for v in pids:
                    if not re.match(r"^ses-.+", v):
                        msg = (
                            f'Session ID must have the form "ses-XXXX".'
                            f"\nCheck demographics file: {p}"
                        )
                        logger.error(msg)
                        exit(1)

        list_df.append(df)

    return list_df[0] if not is_list else list_df


def load_px_demo(
    *,
    sid: str,
    ses: Union[str, None] = None,
    demo_px: PathType,
    actual_to_expected: Dict[str, str],
    col_dtypes: Optional[Dict[str, Any]] = None,
    tmp: str,
):

    df_px = load_demo(demo_px, rename=actual_to_expected, dtypes=col_dtypes, tmp=tmp)
    mask_px = df_px["participant_id"] == sid
    if ses is not None:
        mask_px &= df_px["session_id"] == ses
    return df_px[mask_px]


def _subject_zscore(
    *,
    data_cn: np.ndarray,
    data_px: np.ndarray,
    index_df=None,
    cols_df=None,
    analyses: Optional[List[Analysis]],
):

    res = {}

    if "regional" in analyses:
        z = zscore(data_cn, data_px)
        if index_df is not None:
            z = pd.DataFrame(z.reshape(1, -1), index=index_df, columns=cols_df)

        res["regional"] = z

    if "asymmetry" in analyses:
        if data_cn.shape[-1] == 1:
            xh_cn = data_cn.reshape(2, data_cn.shape[0], -1)
        elif len(data_cn.shape) > 2:
            xh_cn = data_cn.reshape(2, data_cn.shape[0], -1, data_cn.shape[-1])
        else:
            xh_cn = data_cn.reshape(2, data_cn.shape[0], -1)
        data_cn_asym = compute_asymmetry(xh_cn[0], xh_cn[1])

        if data_cn.shape[-1] == 1:
            xh_px = data_px.reshape(2, -1)
        elif len(data_px.shape) > 1:
            xh_px = data_px.reshape(2, -1, data_cn.shape[-1])
        else:
            xh_px = data_px.reshape(2, -1)
        data_px_asym = compute_asymmetry(xh_px[0], xh_px[1])

        cols_df = None if cols_df is None else cols_df[: xh_px.shape[1]]
        res["asymmetry"] = _subject_zscore(
            data_cn=data_cn_asym,
            data_px=data_px_asym,
            index_df=index_df,
            cols_df=cols_df,
            analyses=["regional"],
        )["regional"]

    return res


def fixblur(data, px=False):
    output = []
    for arr in data:
        if len(arr.shape) > 1 if px else len(arr.shape) > 2:
            if arr.shape[-1] > 1:
                for i in range(arr.shape[-1]):
                    if len(arr.shape) == 3:
                        output.append(arr[:, :, i])
                    elif len(arr.shape) == 2:
                        output.append(arr[:, i])
            else:
                output.append(arr.squeeze())
        else:
            output.append(arr.squeeze())
    return output


def _subject_mahalanobis(
    *, data: DefaultDict[str, list], analyses: Optional[List[Analysis]]
):
    list_df_cn = data["df_cn"]
    list_data_cn = []
    common_ids = reduce(np.intersect1d, [df.index for df in list_df_cn])
    for i, (df, x) in enumerate(zip(list_df_cn, data["data_cn"])):
        mask = df.index.isin(common_ids)
        list_data_cn.append(x[mask])
        # list_df_cn[i] = df[mask]

    cols_df = data["cols_df"][0]
    index_df = data["index_df"][0]

    res = {}
    if "regional" in analyses:
        data_cn = fixblur(list_data_cn)
        data_cn = np.stack(data_cn, axis=-1)
        data_px = fixblur(data["data_px"], px=True)
        data_px = np.stack(data_px, axis=-1)

        md = mahalanobis_distance(data_cn, data_px)

        if index_df is not None:
            md = pd.DataFrame(md.reshape(1, -1), index=index_df, columns=cols_df)

        res["regional"] = dict(md=md, data_cn=data_cn, data_px=data_px)

    if "asymmetry" in analyses:
        n = len(data["data_px"])
        list_data_cn = [None] * n
        list_data_px = [None] * n
        list_cols_df = [None] * n
        for i, x_cn in enumerate(data["data_cn"]):
            ogdim = x_cn.shape
            xh_cn = x_cn.reshape(2, x_cn.shape[0], -1)
            list_data_cn[i] = compute_asymmetry(xh_cn[0], xh_cn[1])
            if ogdim[-1] > 1 and len(ogdim) > 2:
                list_data_cn[i] = list_data_cn[i].reshape(
                    ogdim[0], ogdim[1] // 2, ogdim[2]
                )
            x_px = data["data_px"][i]
            ogdim = x_px.shape
            xh_px = x_px.reshape(2, -1)
            list_data_px[i] = compute_asymmetry(xh_px[0], xh_px[1])
            if ogdim[-1] > 1 and len(ogdim) > 1:
                list_data_px[i] = list_data_px[i].reshape(ogdim[0] // 2, ogdim[1])
            if i > len(data["cols_df"]) - 1:
                list_cols_df[i] = data["data_px"][i][: xh_px.shape[1]]
            elif data["cols_df"][i] is None:
                list_cols_df[i] = data["data_px"][i][: xh_px.shape[1]]
            else:
                list_cols_df[i] = data["cols_df"][i][: xh_px.shape[1]]

        data["data_cn"] = list_data_cn
        data["data_px"] = list_data_px
        data["cols_df"] = list_cols_df

        res["asymmetry"] = _subject_mahalanobis(data=data, analyses=["regional"])[
            "regional"
        ]

    return res


def process_feature(
    feat,
    kwds,
    cn_zbrains,
    list_df_cn,
    tmp,
    cov_deconfound,
    px_demo,
    px_zbrains,
    px_sid,
    px_ses,
    analyses,
    n_cn,
    struct,
    pth_analysis,
    logs,
):
    data_mahalanobis = defaultdict(list)
    log = logs[feat]
    kwds = copy.deepcopy(kwds)
    # Load control data
    kwds.update(dict(feat=feat))
    data_cn, df_cn = _load_data(cn_zbrains, df_subjects=list_df_cn, tmp=tmp, **kwds)
    if data_cn is None:
        log["warning"] = f"\t{feat:<15}: \tNo data available for reference subjects."
        return None, log

    if isinstance(data_cn, pd.DataFrame):
        data_cn = data_cn.to_numpy()

    # Deconfounding
    dec = None
    if cov_deconfound is not None and px_demo is not None:
        dec = get_deconfounder(covariates=cov_deconfound)
        data_cn = dec.fit_transform(data_cn, df_cn)

    # Load patient data
    data_px = _load_one(
        px_zbrains, sid=px_sid, ses=px_ses, raise_error=False, tmp=tmp, **kwds
    )
    if data_px is None:
        log["warning"] = f"\t{feat:<15}: \tNo data available for target subject."
        return None, log

    # For subcortex, we have dataframes
    cols_df = index_df = None
    is_df = isinstance(data_px, pd.DataFrame)
    if is_df:
        index_df, cols_df = data_px.index, data_px.columns
        data_px = data_px.to_numpy().ravel()

    # Deconfounding
    if dec:
        df_px = px_demo.to_frame().T
        data_px = dec.transform(data_px.reshape(1, -1), df_px)[0]

    # Analysis: zscoring
    res = _subject_zscore(
        data_cn=data_cn,
        data_px=data_px,
        index_df=index_df,
        cols_df=cols_df,
        analyses=analyses,
    )

    log["info"] = (
        f"\t{feat:<15}: \t[{df_cn.shape[0]}/{n_cn} reference subjects available]"
    )

    # Save results
    for analysis in analyses:
        z = res[analysis]
        if analysis == "regional" and struct != "subcortex":
            z = z.reshape(2, -1)
        _save(pth_analysis, x=z, sid=px_sid, ses=px_ses, analysis=analysis, **kwds)

    # store for mahalanobis
    res = dict(
        data_cn=data_cn,
        data_px=data_px,
        df_cn=df_cn,
        feat=feat,
        index_df=index_df,
        cols_df=cols_df,
    )

    for k, v in res.items():
        data_mahalanobis[k].append(v)
    return data_mahalanobis, log


def run_analysis(
    *,
    px_sid: str,
    px_ses: str = None,
    cn_zbrains: Optional[List[PathType]],
    cn_demo_paths: Optional[List[PathType]],
    px_zbrains: PathType,
    px_demo: Union[pd.Series, None] = None,
    structures: Optional[List[Structure]],
    features: Optional[List[Feature]],
    cov_normative: Union[Optional[List[str]], None] = None,
    cov_deconfound: Union[Optional[List[str]], None] = None,
    smooth_ctx: float,
    smooth_hip: float,
    resolutions: Optional[List[Resolution]],
    labels_ctx: Optional[List[str]],
    labels_hip: Optional[List[str]],
    actual_to_expected: Dict[str, str],
    analyses: Optional[List[Analysis]],
    approach: Approach,
    col_dtypes: Union[Optional[Dict[str, type]], None] = None,
    tmp: str,
    n_jobs: int,
):
    """

    Parameters
    ----------
    px_sid:
        Patient ID
    px_ses:
        Patient Ses
    cn_zbrains:
        Names of zbrains folders in for the controls datasets (zbrains, abrains_v03)
    cn_demo_paths:
        CSVs with demographics
    px_zbrains:
        Name of zbrains folder for patient
    px_demo:
        CSV of patient
    structures:
        list of struct
    features:
        list of feat
    cov_normative:
        list de covariates
    cov_deconfound
        list de covariates
    smooth_ctx
    smooth_hip
    resolutions
    labels_ctx
    labels_hip
    actual_to_expected
    analyses
    approach
    col_dtypes

    Returns
    -------

    """
    logger = logging.getLogger(tmp)

    approach_folder = approach_to_folder[approach]

    logger.debug(f'Logging call: {sys.argv[0]} {" ".join(sys.argv[1:])}')
    # logger.info(f'{analysis.capitalize()} analysis\n')

    pth_analysis = (
        f"{get_subject_dir(px_zbrains, px_sid, px_ses)}/" f"{approach_folder}"
    )

    # Load dataframes ----------------------------------------------------------
    list_df_cn = load_demo(
        cn_demo_paths, rename=actual_to_expected, dtypes=col_dtypes, tmp=tmp
    )
    n_cn = sum(len(df) for df in list_df_cn)

    # Iterables ----------------------------------------------------------------
    iterables = []
    for st in structures:
        if st == "subcortex":
            iterables.append((st, None, None))
        elif st == "cortex":
            iterables += itertools.product([st], resolutions, labels_ctx)
        else:
            iterables += itertools.product([st], resolutions, labels_hip)

    # Available features -------------------------------------------------------
    available_features = {k: defaultdict(dict) for k in structures}
    if "subcortex" in structures:
        available_features["subcortex"] = []

    # Main loop ----------------------------------------------------------------
    for struct, resol, label in iterables:
        print(struct, resol, label)
        s = (
            f"[resolution = {resol:<4}\tlabel = {label:<15}]"
            if struct != "subcortex"
            else ""
        )
        logger.info(f"\nStructure: {struct} {s}")

        # Shared kwds
        smooth = smooth_ctx if struct == "cortex" else smooth_hip
        kwds = dict(
            struct=struct,
            resolution=resol,
            label=label,
            smooth=smooth,
            # analysis=analysis
        )
        log_features = {x: {"info": None, "warning": None} for x in features}
        results = Parallel(n_jobs=n_jobs)(
            delayed(process_feature)(
                feat,
                kwds,
                cn_zbrains,
                list_df_cn,
                tmp,
                cov_deconfound,
                px_demo,
                px_zbrains,
                px_sid,
                px_ses,
                analyses,
                n_cn,
                struct,
                pth_analysis,
                log_features,
            )
            for feat in features
        )
        mahalanobis_dicts, logs = zip(*results)
        for log in logs:
            if log["info"]:
                logger.info(log["info"])
            if log["warning"]:
                logger.warning(log["warning"])

        data_mahalanobis = defaultdict(list)

        # Iterate over the dictionaries and the keys.
        for d in mahalanobis_dicts:

            if d is None:
                continue
            for key, value in d.items():
                data_mahalanobis[key].extend(value)
        # store available features
        if struct == "subcortex":
            available_features[struct] = data_mahalanobis["feat"]
        else:
            available_features[struct][resol][label] = data_mahalanobis["feat"]

        # Mahalanobis
        if len(data_mahalanobis["feat"]) < 2:
            continue

        # Analysis: mahalanobis distance
        res = _subject_mahalanobis(data=data_mahalanobis, analyses=analyses)

        # Save results
        kwds.update({"feat": data_mahalanobis["feat"]})

        n_available_cn = 0
        for analysis in analyses:
            data = res[analysis]
            md = data["md"]
            n_available_cn = data["data_cn"].shape[0]

            if analysis == "regional" and struct != "subcortex":
                md = md.reshape(2, -1)

            _save(pth_analysis, x=md, sid=px_sid, ses=px_ses, analysis=analysis, **kwds)

        # n_available_cn = res['data_cn'].shape[0]
        logger.info(
            f'\n\t{"Mahalanobis":<15}: \t[{n_available_cn}/{n_cn} '
            f"controls available]\n"
        )

    logger.info("Done!\n\n")

    return available_features
