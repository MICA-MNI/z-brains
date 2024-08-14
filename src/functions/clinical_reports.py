"""
Generate the clinical reports for ZBRAINS

zbrains MUST exist on the global environment

USAGE:
    
    outdir='<PATH to output directory>/derivatives/z-brains'
    zbrains='<PATH to repository>/z-brains'

    # Run the function per norm ['z', 'normative'] to create one PDF per norm
    clinical_reports(norm='z', 
                     sub='sub-PX069', 
                     ses='ses-01', 
                     outdir=outdir, 
                     cmap='cmo.balance', 
                     cmap_asymmetry='PRGn',
                     Range=(-3,3),
                     thr=1.96, 
                     thr_alpha=0.3, 
                     color_bar='left', 
                     den='0p5', 
                     tmp='/tmp', 
                     sex='F', 
                     age='40',
                     analysis=['asymmetry', 'regional'], 
                     feature=['ADC', 'FA', 'T1map', 'flair', 'thickness', 'volume']
                     )

"""

import glob
import uuid
import logging
import itertools
from pathlib import Path
import os
from typing import List, Union, Tuple, Dict
import cmocean
import numpy as np
import pandas as pd
import nibabel as nib
from xhtml2pdf import pisa
import platform
from brainspace.datasets import load_mask
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting.surface_plotting import plot_surf
from scipy.spatial.transform import Rotation

if (
    "DISPLAY" not in os.environ or not os.environ["DISPLAY"]
) and platform.system() != "Windows":
    from pyvirtualdisplay import Display


from .constants import (
    LIST_ANALYSES,
    Analysis,
    Approach,
    Resolution,
    Structure,
    Feature,
    struct_to_folder,
    approach_to_folder,
)
from .utils_analysis import (
    get_bids_id,
    map_resolution,
    get_analysis_path_from_template,
    get_subject_dir,
    PathType,
)

cmaps = cmocean.cm.cmap_d


DATA_PATH = Path(__file__).resolve().parent.parent / "data"


# logger = logging.getLogger('analysis_logger')


def plot_hemisphere_lh(
    surf_lh,
    array_name=None,
    color_bar=False,
    color_range=None,
    label_text=None,
    layout_style="row",
    cmap="viridis",
    nan_color=(0, 0, 0, 1),
    zoom=1,
    background=(1, 1, 1),
    size=(400, 400),
    interactive=True,
    embed_nb=False,
    screenshot=False,
    filename=None,
    scale=(1, 1),
    transparent_bg=True,
    **kwargs,
):
    """Plot left hemisphere in lateral and medial views.

    Parameters
    ----------
    surf_lh : vtkPolyData or BSPolyData
        Left hemisphere.
    array_name : str, list of str, ndarray or list of ndarray, optional
        Name of point data array to plot. If ndarray, the array is split for
        the left hemisphere. If list, plot one row per array.
        If None, plot surfaces without any array data. Default is None.
    color_bar : bool, optional
        Plot color bar for each array (row). Default is False.
    color_range : {'sym'}, tuple or sequence.
        Range for each array name. If 'sym', uses a symmetric range. Only used
        if array has positive and negative values. Default is None.
    label_text : dict[str, array-like], optional
        Label text for column/row. Possible keys are {'left', 'right',
        'top', 'bottom'}, which indicate the location. Default is None.
    layout_style : str
        Layout style for hemispheres. If 'row', layout is a single row
        alternating lateral and medial views, from left to right. If 'grid',
        layout is a 2x2 grid, with lateral views in the top row, medial
        views in the bottom row. Default is 'row'.
    nan_color : tuple
        Color for nan values. Default is (0, 0, 0, 1).
    zoom : float or sequence of float, optional
        Zoom applied to the surfaces in each layout entry.
    background : tuple
        Background color. Default is (1, 1, 1).
    cmap : str, optional
        Color map name (from matplotlib). Default is 'viridis'.
    size : tuple, optional
        Window size. Default is (800, 200).
    interactive : bool, optional
        Whether to enable interaction. Default is True.
    embed_nb : bool, optional
        Whether to embed figure in notebook. Only used if running in a
        notebook. Default is False.
    screenshot : bool, optional
        Take a screenshot instead of rendering. Default is False.
    filename : str, optional
        Filename to save the screenshot. Default is None.
    transparent_bg : bool, optional
        Whether to us a transparent background. Only used if
        ``screenshot==True``. Default is False.
    scale : tuple, optional
        Scale (magnification). Only used if ``screenshot==True``.
        Default is None.
    kwargs : keyword-valued args
        Additional arguments passed to the plotter.

    Returns
    -------
    figure : Ipython Image or None
        Figure to plot. None if using vtk for rendering (i.e.,
        ``embed_nb == False``).

    See Also
    --------
    :func:`build_plotter`
    :func:`plot_surf`

    """
    if color_bar is True:
        color_bar = "right"

    surfs = {"lh": surf_lh}
    layout = ["lh", "lh"]

    if isinstance(array_name, np.ndarray):
        if array_name.ndim == 2:
            array_name = [a for a in array_name]
        elif array_name.ndim == 1:
            array_name = [array_name]

    to_remove = []
    if isinstance(array_name, list):
        layout = [layout] * len(array_name)
        array_name2 = []
        n_pts_lh = surf_lh.n_points
        for an in array_name:
            if isinstance(an, np.ndarray):
                name = surf_lh.append_array(an[:n_pts_lh], at="p")
                array_name2.append(name)
                to_remove.append(name)
            else:
                array_name2.append(an)
        array_name = np.asarray(array_name2)[:, None]

    if layout_style == "grid":

        # create 2x2 grid for each array_name and stack altogether
        n_arrays = len(array_name)
        array_names, layouts = [], []
        for a, l in zip(array_name, layout):
            array_names.append(np.full((2, 1), fill_value=a[0]))
            layouts.append(np.array(l).reshape(1, 2).T.tolist())
        array_name = np.vstack(array_names)
        layout = np.vstack(layouts)

        view = [["lateral"], ["medial"]] * n_arrays
        share = "both"
    else:
        view = ["lateral", "medial"]
        share = "r"

    if isinstance(cmap, list):
        cmap = np.asarray(cmap)[:, None]

    # when embed_nb=True, interactive (with panel) only supports one renderer,
    # here we have at least 4
    if embed_nb:
        interactive = False

    kwds = {"view": view, "share": share}
    kwds.update(kwargs)
    res = plot_surf(
        surfs,
        layout,
        array_name=array_name,
        color_bar=color_bar,
        color_range=color_range,
        label_text=label_text,
        cmap=cmap,
        nan_color=nan_color,
        zoom=zoom,
        background=background,
        size=size,
        interactive=interactive,
        embed_nb=embed_nb,
        screenshot=screenshot,
        filename=filename,
        scale=scale,
        transparent_bg=transparent_bg,
        **kwds,
    )

    # remove arrays added to surfaces if any
    # cannot do it if return_plotter=True
    if not kwargs.get("return_plotter", False):
        surf_lh.remove_array(name=to_remove, at="p")

    return res


# -----------------------------------------------------------------------------
def adjectivize_struct(struct: Structure):
    if struct == "cortex":
        return "Cortical"
    elif struct == "subcortex":
        return "Subcortical"
    elif struct == "hippocampus":
        return "Hippocampal"
    else:
        raise ValueError(f"Unknown structure {struct}")


# -----------------------------------------------------------------------------
# Functions for plotting
def plot_surfs(
    surfaces,
    values: List[np.ndarray],
    views: Union[List[str], None] = None,
    size: Union[int, Tuple[int, int], None] = None,
    zoom: Union[float, List[float]] = 1.75,
    color_bar="bottom",
    share="both",
    color_range=(-2, 2),
    cmap="cmo.balance",
    transparent_bg=False,
    **kwargs,
):
    """
    surfaces = [hip_mid_l, hip_unf_l, hip_unf_r, hip_mid_r]  Can be 1 or more
    views = ['dorsal', 'lateral', 'lateral', 'dorsal'] Can be 1 or more
    """

    # Append values to surfaces
    my_surfs = {}
    array_names = []
    for i, surf in enumerate(surfaces):
        name = f"surf{i + 1}"
        surf.append_array(values[i], name=name)

        my_surfs[name] = surf
        array_names.append(name)

    # Set the layout as the list of keys from my_surfs
    layout = [list(my_surfs.keys())]

    if size is None:
        size = (200 * len(surfaces), 350)

    return plot_surf(
        my_surfs,
        layout,
        array_name=array_names,
        view=views,
        color_bar=color_bar,
        color_range=color_range,
        share=share,
        cmap=cmap,
        zoom=zoom,
        size=size,
        transparent_bg=transparent_bg,
        **kwargs,
    )


# -----------------------------------------------------------------------------
# Functions for report
def convert_html_to_pdf(source_html, output_filename: PathType):
    with open(output_filename, "w+b") as result_file:
        pisa_status = pisa.CreatePDF(source_html, dest=result_file)

    # True on success and False on errors
    return pisa_status.err


def report_header_template(
    *,
    sid: str,
    ses: Union[str, None] = None,
    age: Union[float, None] = None,
    sex: Union[str, None] = None,
    analysis: Analysis,
):

    # if ses is None:
    #     ses = ''

    if age is None:
        age = "n/a"
    if sex is None:
        sex = "n/a"

    style = (
        "margin-bottom:0;font-family:gill sans,sans-serif;text-align:center;"
        "font-size:{fontsize}px;color:#505050;{extra}"
    )

    ses_str = "" if ses is None else f" &nbsp; <b>Session</b>: {ses},"

    # Header
    report_header = (
        # Banner
        f'<img id="top" src="{DATA_PATH}/zbrains_banner.png" alt="zbrains">'
        # Title
        f'<p style="{style.format(fontsize=20, extra="")}">'
        f"Clinical Summary &nbsp; | &nbsp; <b> {analysis.capitalize()} analysis </b> "
        f"</p>"
        # Subject's ID - Session - Basic demographics
        f'<p style="{style.format(fontsize=14, extra="margin-top:-100px;")}">'
        f"<b>Subject</b>: {sid},{ses_str} "
        f"&nbsp; <b>Sex</b>: {sex}, &nbsp; <b>Age</b>: {age}"
        f"</p>"
    )

    return report_header


def report_colors(analysis: Analysis = "regional"):
    report = "<hr>"

    style = (
        "margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;"
        "text-align:center;font-size:14px;color:{color}"
    )

    if analysis != "regional":
        report += (
            f'<p style="{style.format(color="#b31b2c")}">'
            "<b> Red </b> = <b>right</b> MORE THAN <b>left</b> "
            "</p>"
            f'<p style="{style.format(color="#13365d")}">'
            "<b> Blue </b> = <b>left</b> MORE THAN <b>right</b> "
            "</p>"
        )
    else:
        report += (
            f'<p style="{style.format(color="#b31b2c")}">'
            "<b> Red </b> = INCREASED compared to controls "
            "</p>"
            f'<p style="{style.format(color="#13365d")}">'
            "<b> Blue </b> = DECREASED compared to controls "
            "</p>"
        )
    return report


def feature_header_template(feature: Union[Feature, List[Feature]], extra=""):
    style = (
        "border:0px solid #666;padding-top:10px;padding-left:5px;"
        "background-color:#eee;font-family:Helvetica,sans-serif;font-size:14px;"
        "text-align:center;color:#5d5070"
    )

    if isinstance(feature, list):
        info = f'Features: {" & ".join(feature)} {extra}'
    else:
        info = f"Feature: {feature} {extra}"

    return f'<p style="{style}">' f"<b> {info} </b>" "</p>"


def report_1x2_table(fig1: PathType, fig2: PathType, height=250):
    style = (
        "style=padding-top:4px;padding-bottom:4px;padding-left:3px;"
        "padding-right:3px;text-align:center"
    )

    return (
        '<table style="border:1px solid white;width:100%">'
        "<tr>"
        f"<td style={style}>"
        f'<img style="height:{height}px;margin-top:-100px;" src="{fig1}">'
        "</td>"
        f"<td style={style}>"
        f'<img style="height:{height}px;margin-top:-100px;" src="{fig2}">'
        "</td>"
        "</tr>"
        "</table>"
    )


# Utility functions ------------------------------------------------------------
def map_subcortical_vertices(x) -> np.ndarray:
    """
    Taken from the ENIGMA toolbox
    https://github.com/MICA-MNI/ENIGMA/blob/master/enigmatoolbox

    Map one value per subcortical area to surface vertices (author:
    @saratheriver)

    Parameters
    ----------
    x:
        Subcortical values with shape = (16,) in the following order:
        L_accumbens, L_amygdala, L_caudate, L_hippocampus, L_pallidun,
        L_putamen, L_thalamus, L_ventricles,
        R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun,
        R_putamen, R_thalamus, R_ventricles

    zbrains output (the fisrt colum is SubjID):
        'Laccumb', 'Lamyg', 'Lcaud', 'Lhippo', 'Lpal', 'Lput','Lthal',
        'Raccumb', 'Ramyg', 'Rcaud', 'Rhippo', 'Rpal', 'Rput', 'Rthal'

    Returns
    -------
    data :
        Mapped data, shape = (51278,)
    """

    n_vert = [
        867,
        1419,
        3012,
        3784,
        1446,
        4003,
        3726,
        7653,
        838,
        1457,
        3208,
        3742,
        1373,
        3871,
        3699,
        7180,
    ]

    x = np.asarray(x)
    if x.shape == (16,):
        return np.repeat(x, n_vert)
    else:
        raise ValueError("Input data must have 16 values.")


# Load surfaces ----------------------------------------------------------------
def _load_surfaces_ctx(resolution: Resolution = "high"):

    res_ctx = map_resolution("cortex", resolution)

    inf_lh = read_surface(f"{DATA_PATH}/fsLR-{res_ctx}.L.inflated.surf.gii")
    inf_rh = read_surface(f"{DATA_PATH}/fsLR-{res_ctx}.R.inflated.surf.gii")

    mask = load_mask(join=True)
    return inf_lh, inf_rh, mask


def _load_surfaces_hip(res: Resolution = "high"):
    res_hip = map_resolution("hippocampus", res)
    label = "midthickness"

    pth_canonical = (
        f"{DATA_PATH}/tpl-avg_space-canonical_den-{res_hip}"
        f"_label-hipp_{label}.surf.gii"
    )
    pth_unfold = (
        f"{DATA_PATH}/tpl-avg_space-unfold_den-{res_hip}"
        f"_label-hipp_{label}.surf.gii"
    )

    mid_rh = read_surface(pth_canonical)
    unf_rh = read_surface(pth_unfold)
    mid_lh = read_surface(pth_canonical)
    unf_lh = read_surface(pth_unfold)

    # Flip right to left surface
    mid_lh.Points[:, 0] *= -1
    unf_lh.Points[:, 0] *= -1

    # vflip = np.ones(hipp_mid_l.Points.shape)
    # vflip[:, 0] = -1
    # hipp_mid_l.Points = hipp_mid_l.Points * vflip
    # hipp_unf_l.Points = hipp_unf_l.Points * vflip

    # Rotate surfaces because Reinder didn't accept my pull request

    # Rotate around Y-axis 270
    rot_y270 = Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0]))
    unf_rh.Points = rot_y270.apply(unf_rh.Points)

    # Rotate around X-axis 90
    rot_y90 = Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0]))
    unf_lh.Points = rot_y90.apply(unf_lh.Points)

    # Rotate around Z-axis 180
    rot_z = Rotation.from_rotvec(np.pi * np.array([0, 0, 1]))
    unf_rh.Points = rot_z.apply(unf_rh.Points)

    # Right Antero-posterior lateral
    lat_rh = read_surface(pth_canonical)
    lat_rh.Points = rot_y270.apply(lat_rh.Points)

    # Left Antero-posterior lateral
    lat_lh = read_surface(pth_canonical)
    lat_lh.Points = rot_y90.apply(
        mid_lh.Points
    )  # TODO: should it be lat_lh instead of mid_lh?

    return lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh


# Load data --------------------------------------------------------------------
def _load_data_sctx(
    sctx_file: PathType,
    analysis: Analysis,
    threshold: Union[float, None] = None,
    threshold_alpha: float = 0.5,
):

    # Load CSV data
    x = pd.read_csv(sctx_file, header=[0], index_col=0).to_numpy().ravel()
    if threshold is not None:
        x[np.abs(x) < threshold] *= threshold_alpha
    if analysis == "asymmetry":
        print("e")
    # Array of data
    array_16 = np.full(16, np.nan)
    array_16[0:7] = x[0:7]
    if analysis != "asymmetry":
        array_16[8:15] = x[7:]
    else:
        array_16[8:15] = x[0:7]
        # lest is positive, rights is negative
        # array_16[0:7][(array_16[0:7]) < 0] = 0
        array_16[8:15] = 0

    # Map array to vertices of surfaces
    feat_map = map_subcortical_vertices(array_16)

    # Map array to left and right surfaces
    data_lh = feat_map[0:25910]
    data_rh = feat_map[25910:]

    return data_lh, data_rh


def load_data_struct(
    struct: Structure,
    *,
    file_lh: PathType,
    file_rh: Union[PathType, None] = None,
    analysis: Analysis,
    threshold: Union[float, None] = None,
    threshold_alpha: float = 0.5,
):

    if struct == "subcortex":
        return _load_data_sctx(
            file_lh, analysis, threshold=threshold, threshold_alpha=threshold_alpha
        )

    data_lh = nib.load(file_lh).darrays[0].data
    if analysis != "asymmetry":
        data_rh = nib.load(file_rh).darrays[0].data
    else:
        data_rh = data_lh.copy()
        # left is positive, right is negative
        # data_lh[data_lh < 0] = 0
        data_rh[:] = 0

    if threshold is not None:
        data_lh[np.abs(data_lh) < threshold] *= threshold_alpha
        data_rh[np.abs(data_rh) < threshold] *= threshold_alpha

    return data_lh, data_rh


# Generate figures -------------------------------------------------------------
def _make_png_ctx(
    *,
    analysis,
    data_lh: np.ndarray,
    data_rh: np.ndarray,
    out_png: PathType,
    res: Resolution = "high",
    cmap="cmo.balance",
    color_range=(-2, 2),
    color_bar="bottom",
):

    slh, srh, mask = _load_surfaces_ctx(resolution=res)

    # Replace values in f with NaN where mask_32k is False
    feat_map = np.hstack(np.concatenate((data_lh, data_rh), axis=0))
    feat_map[~mask] = np.nan

    kwds = dict(
        cmap=cmap,
        color_bar=color_bar,
        color_range=color_range,
        transparent_bg=False,
        nan_color=(0, 0, 0, 1),
        embed_nb=True,
        interactive=False,
        screenshot=True,
    )

    # Plot the mean FEATURE on fsLR-32k
    out_png = Path(out_png)
    if analysis == "asymmetry":
        plot_hemisphere_lh(
            slh,
            array_name=feat_map,
            layout_style="grid",
            label_text={"left": ["Lateral", "Medial"], "top": [""]},
            share=None,
            zoom=1.25,
            size=(600, 600),
            filename=out_png,
            **kwds,
        )
    else:
        plot_hemispheres(
            slh,
            srh,
            array_name=feat_map,
            layout_style="grid",
            label_text={"left": ["Lateral", "Medial"], "top": ["Left", "Right"]},
            share=None,
            zoom=1.25,
            size=(600, 600),
            filename=out_png,
            **kwds,
        )

    # Plot cortex with different views
    out_png2 = os.path.join(
        os.path.dirname(out_png),
        os.path.splitext(os.path.basename(out_png))[0]
        + "_SI"
        + os.path.splitext(out_png)[1],
    )
    if analysis == "asymmetry":
        kwds["text__textproperty"] = {"fontSize": 30}
        plot_surfs(
            surfaces=[slh, slh],
            values=[data_lh, data_lh],
            views=["dorsal", "ventral"],
            label_text={
                "bottom": [
                    "Left Superior",
                    "Left Inferior",
                ]
            },
            # share='both', zoom=2.95, size=(550, 350),
            share="both",
            zoom=1.25,
            size=(550, 450),
            filename=out_png2,
            **kwds,
        )

    else:
        plot_surfs(
            surfaces=[slh, srh, srh, slh],
            values=[data_lh, data_rh, data_rh, data_lh],
            views=["dorsal", "dorsal", "ventral", "ventral"],
            label_text={
                "bottom": [
                    "Left Superior",
                    "Right Superior",
                    "Right Inferior",
                    "Left Inferior",
                ]
            },
            # share='both', zoom=2.95, size=(550, 350),
            share="both",
            zoom=2.75,
            size=(550, 350),
            filename=out_png2,
            **kwds,
        )

    return report_1x2_table(fig1=out_png, fig2=out_png2, height=300)


def _make_png_sctx(
    *,
    analysis,
    data_lh: np.ndarray,
    data_rh: np.ndarray,
    out_png: PathType,
    cmap="cmo.balance",
    color_range=(-2, 2),
):

    slh = read_surface(f"{DATA_PATH}/sctx.L.surf.gii", itype="gii")
    srh = read_surface(f"{DATA_PATH}/sctx.R.surf.gii", itype="gii")
    # Plot subcortical structures
    if analysis == "asymmetry":
        kwds = dict()
        kwds["text__textproperty"] = {"fontSize": 50}
        plot_surfs(
            surfaces=[slh, slh],
            values=[data_lh, data_lh],
            views=["lateral", "medial"],
            label_text={"left": ["left"]},
            cmap=cmap,
            color_bar="bottom",
            color_range=color_range,
            share="both",
            transparent_bg=True,
            nan_color=(0, 0, 0, 0),
            zoom=1.4,
            size=(900, 250),
            embed_nb=True,
            interactive=False,
            screenshot=True,
            filename=out_png,
            **kwds,
        )

    else:
        plot_surfs(
            surfaces=[slh, slh, srh, srh],
            values=[data_lh, data_lh, data_rh, data_rh],
            views=["lateral", "medial", "lateral", "medial"],
            label_text={"left": ["left"], "right": ["right"]},
            cmap=cmap,
            color_bar="bottom",
            color_range=color_range,
            share="both",
            transparent_bg=True,
            nan_color=(0, 0, 0, 0),
            zoom=1.4,
            size=(900, 250),
            embed_nb=True,
            interactive=False,
            screenshot=True,
            filename=out_png,
        )

    return (
        f'<p style="text-align:center;margin-left=0px;"> '
        f'<a href="{out_png}" target="_blank">'
        f'<img style="height:150px;margin-top:-100px;" src="{out_png}"> '
        f"</a> "
        f"</p>"
    )


def _make_png_hip(
    *,
    analysis,
    data_lh: np.ndarray,
    data_rh: np.ndarray,
    out_png: PathType,
    res: Resolution = "high",
    cmap="cmo.balance",
    color_range=(-2, 2),
    color_bar="bottom",
):

    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip(res=res)
    if analysis == "asymmetry":
        kwds = dict()
        kwds["text__textproperty"] = {"fontSize": 50}
        plot_surfs(
            surfaces=[lat_lh, mid_lh, unf_lh],
            values=[data_lh, data_lh, data_lh],
            views=["dorsal", "dorsal", "lateral"],
            color_bar=color_bar,
            zoom=1.75,
            cmap=cmap,
            color_range=color_range,
            interactive=False,
            screenshot=True,
            filename=out_png,
        )

    else:
        plot_surfs(
            surfaces=[lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh],
            values=[data_lh, data_lh, data_lh, data_rh, data_rh, data_rh],
            views=["dorsal", "dorsal", "lateral", "lateral", "dorsal", "dorsal"],
            color_bar=color_bar,
            zoom=1.75,
            cmap=cmap,
            color_range=color_range,
            interactive=False,
            screenshot=True,
            filename=out_png,
        )

    return (
        f'<p style="text-align:center;margin-left=0px;"> '
        f'<a href="{out_png}" target="_blank">'
        f'<img style="height:175px;margin-top:-100px;" src="{out_png}"> '
        f"</a> "
        f"</p>"
    )


def make_png(
    struct: Structure,
    *,
    feat_lh: np.ndarray,
    feat_rh: Union[np.ndarray, None] = None,
    out_png: PathType,
    res: Union[Resolution, None] = None,
    cmap="cmo.balance",
    color_range=(-2, 2),
    color_bar="bottom",
    analysis=None,
):

    if struct == "cortex":
        return _make_png_ctx(
            analysis=analysis,
            data_lh=feat_lh,
            data_rh=feat_rh,
            out_png=out_png,
            res=res,
            cmap=cmap,
            color_range=color_range,
            color_bar=color_bar,
        )

    if struct == "hippocampus":
        return _make_png_hip(
            analysis=analysis,
            data_lh=feat_lh,
            data_rh=feat_rh,
            out_png=out_png,
            res=res,
            cmap=cmap,
            color_range=color_range,
            color_bar=color_bar,
        )

    return _make_png_sctx(
        analysis=analysis,
        data_lh=feat_lh,
        data_rh=feat_rh,
        out_png=out_png,
        cmap=cmap,
        color_range=color_range,
    )


def make_png_missing(struct: Structure):
    st = adjectivize_struct(struct)
    return (
        f'<p style="margin-bottom:0;margin-top:0;font-family:gill sans,'
        f'sans-serif;text-align:center;font-size:14px;color:#ffb311"> '
        f"<b> [WARNING] </b>{st} file was not found </p>"
    )


def report_struct(
    *,
    struct: Structure,
    path_analysis: PathType,
    sid: str,
    ses: Union[str, None] = None,
    analysis: Analysis,
    approach: Approach,
    feat: Union[Feature, List[Feature]],
    thr: Union[float, None] = None,
    thr_alpha=0.5,
    smooth: Union[float, None] = None,
    res: Union[Resolution, None] = None,
    label: Union[str, None] = None,
    cmap="cmo.balance",
    color_range=(-2, 2),
    color_bar="bottom",
    tmp_dir: PathType = "/tmp",
):
    logger = logging.getLogger(tmp_dir)
    bids_id = get_bids_id(sid, ses)
    root_path = f"{path_analysis}/{struct_to_folder[struct]}"

    if isinstance(feat, list):
        feat = "-".join(feat)

    kwds = dict(root_path=root_path, bids_id=bids_id, feat=feat, analysis=analysis)

    file_rh = None
    if struct != "subcortex":
        kwds.update(dict(feat=feat, smooth=smooth, res=res, label=label, hemi="L"))
        file_lh = get_analysis_path_from_template(struct, **kwds)

        if analysis != "asymmetry":
            kwds.update(dict(hemi="R"))
            file_rh = get_analysis_path_from_template(struct, **kwds)
    else:
        file_lh = get_analysis_path_from_template(struct, **kwds)
    # html
    thr_str = "" if thr is None else f" | threshold={thr}"
    info = ""
    if struct != "subcortex":
        struct_res = map_resolution(struct, res)
        info = f"| {smooth}mm smooth | resolution {struct_res}"

    html = (
        '<p style="margin-bottom:0;margin-top:0;'
        "font-family:gill sans,sans-serif;text-align:left;font-size:14px;"
        'color:#5d5070"> '
        # f'<b>{adjectivize_struct(struct)} {feat}</b> | {info} '
        f"<b>{struct.capitalize()}</b>"
        # f'{analysis} analysis | approach-{approach} {thr_str}'
        f"| {approach} approach {info}"
        "</p>"
    )

    lh_exists = file_lh.exists()
    if not lh_exists or file_rh is not None and not file_rh.exists():

        missing = file_rh if lh_exists else file_lh

        logger.warning(f"{adjectivize_struct(struct)} file was not found:\n{missing}")

        png_block = make_png_missing(struct)
        html += png_block
        return html

    feat_lh, feat_rh = load_data_struct(
        struct,
        file_lh=file_lh,
        file_rh=file_rh,
        analysis=analysis,
        threshold=thr,
        threshold_alpha=thr_alpha,
    )

    if analysis == "asymmetry":
        color_range = (-1.5, 1.5)

    out_png = (
        f"{tmp_dir}/{bids_id}_{struct}_feature-{feat}_"
        f"analysis-{analysis}{thr}_{uuid.uuid4()}.png"
    )
    # if analysis == "asymmetry":
    #     feat_lh = feat_lh + feat_rh
    png_block = make_png(
        struct,
        feat_lh=feat_lh,
        feat_rh=feat_rh,
        out_png=out_png,
        res=res,
        cmap=cmap,
        color_range=color_range,
        color_bar=color_bar,
        analysis=analysis,
    )

    html += png_block
    return html


def generate_clinical_report(
    *,
    zbrains_path: PathType = None,
    sid: str,
    ses: Union[str, None] = None,
    age: Union[float, None] = None,
    sex: Union[str, None] = None,
    analyses: Union[List[Analysis], None] = None,
    features: Union[List[Feature], None] = None,
    approach: Approach = "zscore",
    threshold=1.96,
    threshold_alpha=0.3,
    smooth_ctx: float = 5,
    smooth_hip: float = 2,
    res_ctx: Resolution = "high",
    res_hip: Resolution = "high",
    label_ctx="midthickness",
    label_hip="midthickness",
    color_bar="bottom",
    cmap="cmo.balance",
    color_range=(-2, 2),
    cmap_asymmetry="cmo.balance_r",
    tmp_dir: PathType = "/tmp",
    subject_dir=None,
    output_dir=None,
    tag=None,
):
    """Zbrains: Clinical report generator

      Global Parameters
      -----------------
      zbrains  :  str
      MUST exist and is the PATH to the z-brain repository <path_to/z-brains>

      Parameters
      ----------
      zbrains_path:
          Output directory for processed z-brains derivatives.
      sid
          Participant id
      ses
          OPTIONAL flag that indicates the session name (if omitted will
          manage as SINGLE session)
      approach:
          Comparison approach.
      features:
          Default=None it will run all the features unless a list of
          features is provided
      analyses :
          Default=None it will run all the analyses unless a list of
          analyses is provided
      cmap    : str
          Default='cmo.balance'. Colormap for plotting the regional changes.
      cmap_asymmetry
          Default='PRGn'. Colormap for plotting the regional changes.
      color_range   : tuple
          Default=(-2,2). color_range of values for the regional plot
      color_bar: str
          Default='right' position of the color bar for hippocampal and
          cortical plots
      res_ctx :
          Surface resolution for cortex.
      res_hip :
          Surface resolution for hippocampus.
      label_ctx :
          Cortical surface used when mapping the data from volume to
          surface.
      label_hip :
          Hippocampal surface used when mapping the data from volume to
          surface.
      threshold :
          Default=1.96. Threshold for the maps
      threshold_alpha :
          Default=0.3 Alpha channel for the colormap when thr is applied
          0=white, 1=full color
      smooth_ctx     : int
          Default=5 Smoothing for the cortical data in mm.
      smooth_hip     : int
          default=2 Smoothing for the hippocampal data in mm.
      age     : int
          Defaul='unknown'. OPTIONAL Covariate for the normative model.
      sex:
          Participant's sex. Must be {'F', 'M'} .
      tmp_dir:
          Default='/tmp'. Working directory for temporary files

    NOTE: clinical_reports does not take 'volume' as a feature.
          if it exists it will be removed. To process morphometry only
          'thickness' is allowed

    McGill University, MNI, MICA-lab, Created 24 November 2023
    'Zbrains Hackathon', MICA collective
    https://github.com/MICA-MNI/micapipe
    https://mica-mni.github.io/
    """
    logger = logging.getLogger(tmp_dir)
    approach_folder = approach_to_folder[approach]

    bids_id = get_bids_id(sid, ses)
    if zbrains_path is None:
        if subject_dir is not None:
            path_analysis = f"{subject_dir}/{approach_folder}"
        else:
            print("Either zbrains_path or subject_dir must be provided")
    else:
        "If zbrains_path is provided, the subject_dir is calculated"
        subject_dir = get_subject_dir(zbrains_path, sid, ses)
        path_analysis = f"{subject_dir}/{approach_folder}"

    # List all files of norm
    subses_files = glob.glob(f"{path_analysis}/*/*")

    # If ANALYSIS is empty run all available analysis
    if analyses is None:
        available_analyses = sorted(
            list(
                set(
                    [
                        file.split(".")[0].split("analysis-")[1].split("_")[0]
                        for file in subses_files
                    ]
                )
            )
        )
        analyses: list[str] = [a for a in LIST_ANALYSES if a in available_analyses]

    # If FEATURES is empty run all available features
    if features is None:
        features = sorted(
            list(
                set([file.split("feature-")[1].split("_")[0] for file in subses_files])
            )
        )

    # logger.info(f'{bids_id} features:\n{features}')

    # Remove volume from the features
    if "volume" in features:
        features.remove("volume")

    report = ""
    # for analysis in analyses:

    display_flag = False
    if (
        "DISPLAY" not in os.environ or not os.environ["DISPLAY"]
    ) and platform.system() != "Windows":

        os.environ["PYVIRTUALDISPLAY_DISPLAYFD"] = "0"
        # Display for headless plotting
        dsize = (900, 750)
        display = Display(visible=False, size=dsize)
        display.start()
        display_flag = True

    for analysis, feat, thresh in itertools.product(
        analyses, features, [None, threshold]
    ):

        if isinstance(feat, list):
            feat_sctx = ["volume" if feat == "thickness" else v for v in feat]
        else:
            feat_sctx = "volume" if feat == "thickness" else feat

        if thresh is None:
            logger.info(
                f"Running summary of analysis={analysis}, "
                f"approach={approach}, feature={feat}"
            )

        # -------------------------------------------------------------------------
        # Generate the report Title
        report += report_header_template(
            sid=sid, ses=ses, age=age, sex=sex, analysis=analysis
        )

        extra = "" if thresh is None else f"| Threshold: {thresh}"
        report += feature_header_template(feat, extra=extra)

        kwds = dict(
            path_analysis=path_analysis,
            sid=sid,
            ses=ses,
            analysis=analysis,
            approach=approach,
            thr=thresh,
            thr_alpha=threshold_alpha,
            color_range=color_range,
            cmap=cmap_asymmetry if analysis == "asymmetry" else cmap,
            color_bar=color_bar,
            tmp_dir=tmp_dir,
        )

        report += report_struct(
            struct="cortex",
            feat=feat,
            res=res_ctx,
            label=label_ctx,
            smooth=smooth_ctx,
            **kwds,
        )

        report += report_struct(struct="subcortex", feat=feat_sctx, **kwds)

        report += report_struct(
            struct="hippocampus",
            feat=feat,
            res=res_hip,
            label=label_hip,
            smooth=smooth_hip,
            **kwds,
        )

        report += report_colors(analysis=analysis)

        # page break
        report += '<div style="page-break-after: always;"></div>'

    # Report file name
    if output_dir and tag:
        file_pdf = f"{output_dir}/{tag}.pdf"
        if os.path.exists(file_pdf):
            os.remove(file_pdf)
    else:
        file_pdf = f"{subject_dir}/{bids_id}_approach-{approach}_summary-report.pdf"

    # Create the HTML file
    convert_html_to_pdf(report, file_pdf)

    logger.info(f"Clinical report successfully created: {file_pdf}")
    if display_flag:
        display.stop()
        del display

    if output_dir and tag:
        return os.path.realpath(file_pdf)
