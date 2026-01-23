import io
import re
import glob
import base64
import uuid
import time
import random
import string
import segno
from datetime import date
import matplotlib.pyplot as plt
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
from brainspace.mesh.mesh_creation import build_polydata
from scipy.spatial.transform import Rotation
from zbrains.report_functions import project_to_native_surface
if (
    "DISPLAY" not in os.environ or not os.environ["DISPLAY"]
) and platform.system() != "Windows":
    from pyvirtualdisplay import Display


from zbrains.constants import (
    LIST_ANALYSES,
    Analysis,
    Approach,
    Resolution,
    Structure,
    Feature,
)
from zbrains.utils_analysis import (
    map_resolution,
    get_subject_dir,
    PathType,
)
import copy

cmaps = cmocean.cm.cmap_d


DATA_PATH = Path(__file__).resolve().parent / "data"

def uuid_qr(qr_path="/tmp", scale=5, border=0):
    """Generate a UUID-based QR code for report tracking."""
    # Generate a UUID using "zbrain", timestamp, and a random string
    timestamp = str(int(time.time()))
    rand_str = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    base_string = f"zbrain-{timestamp}-{rand_str}"
    uuid_generated = uuid.uuid5(uuid.NAMESPACE_DNS, base_string)
    
    # Create QR code file path
    qr_path_file = os.path.join(qr_path, f'{uuid_generated}.png')
    
    # Generate QR code using segno
    qr = segno.make(str(uuid_generated))
    
    # Save it as PNG with foreground color #A569BD and white background
    qr.save(
        qr_path_file,
        scale=scale,
        border=border,
        dark='#A569BD',
        light='white'
    )
    return uuid_generated, qr_path_file

def features_table(sid, ses, valid_subjects):
    """Generate an HTML table showing available features for each brain structure."""
    
    # Generate HTML table
    html = """
    <style>
        table.features {
            width: 100%;
            border-collapse: collapse;
            font-family: "Gill Sans", sans-serif;
            font-size: 12px;
        }
        table.features th {
            background-color: #e3d7f4;
            color: #505050;
            padding: 4px;
            text-align: center;
            line-height: 1.2;
        }
        table.features td {
            padding: 4px;
            border: 1px solid #ddd;
            text-align: center;
            line-height: 1.2;
        }
        table.features tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        table.features tr:hover {
            background-color: #e0f7fa;
        }
    </style>

    <table class="features">
        <thead>
            <tr>
                <th>Feature</th>
                <th>Cortex</th>
                <th>Hippocampus</th>
                <th>Subcortex</th>
            </tr>
        </thead>
        <tbody>
    """
    
    # Use HTML entities instead of Unicode characters for better PDF compatibility
    CHECK = "&#10004;"  # HTML entity for checkmark
    CROSS = "&#10008;"  # HTML entity for cross mark
    
    # Current subject tuple
    current_subject = (sid, ses)
    
    # Collect all features from valid_subjects structure
    all_features = set()
    
    # Get features from each feature key in valid_subjects (excluding 'base' and 'structures')
    for key in valid_subjects.keys():
        if key not in ['base', 'structures'] and isinstance(valid_subjects[key], dict):
            all_features.add(key)
    
    # Sort features for consistent display
    sorted_features = sorted(all_features)
    
    for feature in sorted_features:
        # Check if this subject has this feature for each structure
        ctx_check = CROSS
        hip_check = CROSS
        sub_check = CROSS
        
        if feature in valid_subjects:
            feature_data = valid_subjects[feature]
            
            # Check cortex
            if 'structures' in feature_data and 'cortex' in feature_data['structures']:
                if current_subject in feature_data['structures']['cortex']:
                    ctx_check = CHECK
            
            # Check hippocampus
            if 'structures' in feature_data and 'hippocampus' in feature_data['structures']:
                if current_subject in feature_data['structures']['hippocampus']:
                    hip_check = CHECK
            
            # Check subcortical
            if 'structures' in feature_data and 'subcortical' in feature_data['structures']:
                if current_subject in feature_data['structures']['subcortical']:
                    sub_check = CHECK
    
        html += f"""
        <tr>
            <td>{feature}</td>
            <td>{ctx_check}</td>
            <td>{hip_check}</td>
            <td>{sub_check}</td>
        </tr>
        """
    
    html += "</tbody></table>"
    
    return html

def controls_summary_html(data, subject_age=None):
    """Generate HTML with control database statistics and age distribution chart."""
    # Compute statistics
    df = data.data
    
    # Create case-insensitive mapping from lowercase names to actual column names
    lower_normcols = [col.lower() for col in data.normative_columns]
    column_map = {col.lower(): col for col in data.normative_columns}
    
    # Check if age and sex columns exist (case insensitive)
    if 'age' in lower_normcols and 'sex' in lower_normcols:
        age_col = column_map['age']
        sex_col = column_map['sex']
        
        # Get binary encoding dictionary if it exists
        binary_encoding = getattr(data, 'binary_encodings', {}).get(sex_col, None)
        
        # If sex is binary encoded (e.g., {'m':1, 'f':0}), handle appropriately
        is_binary_encoded = binary_encoding is not None
        
        if is_binary_encoded:
            # Create reverse mapping from codes to labels
            sex_codes = {v: k.upper() for k, v in binary_encoding.items()}
            male_code = binary_encoding.get('M', None)  # Default to None if not specified
            female_code = binary_encoding.get('F', None)  # Default to None if not specified

            # Count by binary codes
            Fnumber = (df[sex_col] == female_code).sum()
            Mnumber = (df[sex_col] == male_code).sum()
        else:
            # Assume traditional 'M'/'F' or similar strings
            Fnumber = (df[sex_col] == 'F').sum()
            Mnumber = (df[sex_col] == 'M').sum()
        
        mean_age = df[age_col].mean()
        sd_age = df[age_col].std()
        Ntotal = len(df)
        Fperc = 100 * Fnumber / Ntotal if Ntotal else 0
        Mperc = 100 * Mnumber / Ntotal if Ntotal else 0

        # Bin and count
        bins = np.arange(df[age_col].min() // 5 * 5, df[age_col].max() // 5 * 5 + 10, 5)
        labels = [f"{int(b)}-{int(b + 4)}" for b in bins[:-1]]
        df['age_bin'] = pd.cut(df[age_col], bins=bins, labels=labels, right=True, include_lowest=True)
        
        # Add a text sex column for groupby if binary encoded
        if is_binary_encoded:
            df['sex_label'] = df[sex_col].map(sex_codes)
            group_sex_col = 'sex_label'
        else:
            group_sex_col = sex_col
            
        age_bin_counts = df.groupby(['age_bin', group_sex_col]).size().unstack(fill_value=0)
        
        # Ensure both columns exist
        m_column = 'M'
        f_column = 'F'
        
        if is_binary_encoded:
            m_column = sex_codes.get(male_code, 'M')
            f_column = sex_codes.get(female_code, 'F')
            
        for sex in [m_column, f_column]:
            if sex not in age_bin_counts.columns:
                age_bin_counts[sex] = 0
                
        # Create plot
        fig, ax = plt.subplots(figsize=(5, 2), dpi=150)
        ax.bar(age_bin_counts.index, age_bin_counts[m_column], color="#34495E", label=m_column, width=0.9)
        ax.bar(age_bin_counts.index, age_bin_counts[f_column], bottom=age_bin_counts[m_column], color="#839192", label=f_column, width=0.9)
        
        ax.set_xlabel('Age')
        ax.set_ylabel('Count')
        plt.xticks(rotation=45)
        
        # Add vertical line at subject age if available
        try:
            if subject_age is not None:
                subject_age = float(subject_age)
                subject_bin = pd.cut([subject_age], bins=bins, labels=labels, right=True, include_lowest=True)[0]
                if pd.notna(subject_bin):
                    x_tick_labels = list(age_bin_counts.index)
                    if subject_bin in x_tick_labels:
                        x_position = x_tick_labels.index(subject_bin)
                        ax.axvline(x=x_position, color="#A93226", linestyle='-', linewidth=2)
        except Exception as e:
            pass  # Continue without age line if there's an error
        
        # Custom legend
        handles, labels_ = ax.get_legend_handles_labels()
        label_order = [f_column, m_column]
        ordered_handles = [handles[labels_.index(lbl)] for lbl in label_order if lbl in labels_]
        ax.legend(ordered_handles, label_order)
        
        # Clean up plot
        plt.tight_layout()
        for spine in ['top', 'right', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)
        
        # Save to buffer as base64 PNG
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        img = base64.b64encode(buf.getvalue()).decode('utf-8')
        buf.close()

        # HTML layout
        html = (
            f"""
        <table style="width:100%; border-collapse:collapse; font-family:gill sans, sans-serif; font-size:12px;">
          <tr>
            <td style="vertical-align:top; padding:8px;">
              <table style="border:none; border-collapse:collapse; font-size:12px;">
                <tr>
                  <td style="padding:4px 8px; text-align:left; font-weight:bold;">Age</td>
                  <td style="padding:4px 8px; text-align:left;">{mean_age:.1f} &plusmn; {sd_age:.1f}</td>
                </tr>
                <tr>
                  <td style="padding:4px 8px; text-align:left; font-weight:bold;">Total</td>
                  <td style="padding:4px 8px; text-align:left;">{Ntotal}</td>
                </tr>
                <tr>
                  <td style="padding:4px 8px; text-align:left; font-weight:bold;">Female</td>
                  <td style="padding:4px 8px; text-align:left;">{Fnumber} ({Fperc:.1f}%)</td>
                </tr>
                <tr>
                  <td style="padding:4px 8px; text-align:left; font-weight:bold;">Male</td>
                  <td style="padding:4px 8px; text-align:left;">{Mnumber} ({Mperc:.1f}%)</td>
                </tr>
              </table>
            </td>
            <td style="vertical-align:top; padding:8px;">
              <img src="data:image/png;base64,{img}" style="width:100%; max-width:525px;"/>
            </td>
          </tr>
        </table>
        """
        )
        return html

def pipeline_info(tmp_dir="/tmp"):
    """Generate HTML with pipeline information and QR code."""
    try:
        user = os.getlogin()
    except Exception:
        user = "unknown"
    try:
        workstation = os.uname().nodename
    except Exception:
        workstation = "unknown"

    # Badge and link HTML
    version_badge = (
        '<img src="https://img.shields.io/github/v/tag/MICA-MNI/z-brains" '
        'alt="Version" style="vertical-align:middle;">'
    )
    license_badge = (
        '<img src="https://img.shields.io/badge/license-BSD-brightgreen" '
        'alt="License" style="vertical-align:middle;">'
    )
    license_link = (
        '<a href="https://github.com/MICA-MNI/z-brains/blob/main/LICENSE" '
        'target="_blank" style="text-decoration:none;color:#5d5070;">  </a>'
    )
    github_link = (
        '<a href="https://github.com/MICA-MNI/z-brains" '
        'target="_blank" style="text-decoration:none;color:#5d5070;">github.com/MICA-MNI/z-brains</a>'
    )
    docs_link = (
        '<a href="https://z-brains.readthedocs.io" '
        'target="_blank" style="text-decoration:none;color:#5d5070;">z-brains.readthedocs.io</a>'
    )
    today_date = date.today().strftime("%B %d, %Y")

    # Create a UUID-QR
    uuid_generated, qr_path = uuid_qr(tmp_dir, scale=5, border=0)
    
    # Table HTML with explicit width control (50%) and table-based layout
    info_table = f'''
    <table style="width:50%; font-family:gill sans, sans-serif;">
      <tr>
        <td style="vertical-align:top; width:70%;">
          <table style="border-collapse:collapse; font-size:11px; width:100%;">
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">Version</td>
              <td style="padding:4px 8px;text-align:left;">{version_badge}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">License</td>
              <td style="padding:4px 8px;text-align:left;">{license_badge} {license_link}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">GitHub</td>
              <td style="padding:4px 8px;text-align:left;">{github_link}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">Documentation</td>
              <td style="padding:4px 8px;text-align:left;">{docs_link}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">Date</td>
              <td style="padding:4px 8px;text-align:left;">{today_date}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">User</td>
              <td style="padding:4px 8px;text-align:left;">{user}</td>
            </tr>
            <tr>
              <td style="padding:4px 8px;text-align:left;font-weight:bold;">Workstation</td>
              <td style="padding:4px 8px;text-align:left;">{workstation}</td>
            </tr>
          </table>
        </td>
        <td style="vertical-align:top; width:30%;">
          <img src="{qr_path}" alt="QR Code" style="width:150px; height:auto; border:1px solid #ccc;" />
        </td>
      </tr>
    </table>
    '''

    return info_table

def inflate_surf(orig_lh, orig_rh, ref_lh, ref_rh, W=0.15):
    """
    This functions inflates a surface a determined amount [W{0:1}] to a reference
    
    Parameters
    ----------
        orig_lh : vtkPolyData or VTKObjectWrapper
                    left input original surface
        orig_rh : vtkPolyData or VTKObjectWrapper
                    right input original surface
        ref_lh  : vtkPolyData or VTKObjectWrapper
                    left input reference surface
        ref_rh  : vtkPolyData or VTKObjectWrapper
                    right input reference surface
        W       : float value from 0 to 1
        
    Returns
    -------
        new_lh : vtkPolyData or VTKObjectWrapper
                    left output inflated surface
        new_rh : vtkPolyData or VTKObjectWrapper
                    right output inflated surface
        
    """
    def inflate(orig, ref, Winf):
        # Convert BSPolyData objects to numpy arrays
        inf_coord = copy.copy(orig.points)
        inf_triag = copy.copy(orig.GetCells2D())
        
        # Inflated mean surface
        maxs = np.max(orig.points, axis=0)
        mins = np.min(orig.points, axis=0)
        maxsp = np.max(ref.points, axis=0)
        minsp = np.min(ref.points, axis=0)
        
        for i in range(3):
            inf_coord[:,i] = (((ref.points[:, i] - minsp[i]) / (maxsp[i] - minsp[i]))
                                 * (maxs[i] - mins[i]) + mins[i]) * Winf + orig.points[:, i] * (1 - Winf)
        
        # Create the new surface
        new_surf = build_polydata(inf_coord, cells=inf_triag)
        
        return(new_surf)
    
    new_lh = inflate(orig_lh, ref_lh, Winf=W)
    new_rh = inflate(orig_rh, ref_rh, Winf=W)
    
    return(new_lh, new_rh)

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
    analysis = None,
):

    analysis = analysis + " analysis" if analysis != 'Control statistics' else 'Control statistics'
    # if ses is None:
    #     ses = ''

    if age is None:
        age = "n/a"
    else:
        # Handle both scalar values and pandas Series
        if hasattr(age, 'values'):
            age = age.values[0]
        else:
            age = age  # Already a scalar value
    
    if sex is None:
        sex = "n/a"
    else:
        # Handle both scalar values and pandas Series
        if hasattr(sex, 'values'):
            sex = sex.values[0]
        else:
            sex = sex  # Already a scalar value

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
        f"Clinical Summary &nbsp; | &nbsp; <b> {analysis.capitalize()} </b> "
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


def feature_header_template(
    feature: Union[Feature, List[Feature]] = "",
    extra: str = "",
    align: str = "center",
    size: int = 14,
    info = None
) -> str:
    """
    Generates an HTML header for a feature or list of features, with configurable alignment
    and the option to provide a custom info string.

    Args:
        feature: A single feature or a list of features.
        extra: Additional text to append.
        align: Text alignment ('left', 'center', etc.).
        size: Font size in pixels.
        info: Custom info string to display. If provided, overrides feature/extra.

    Returns:
        HTML string for the feature header.
    """
    style = (
        "border:0px solid #666;padding-top:10px;padding-left:5px;"
        f"background-color:#eee;font-family:gill sans,sans-serif;font-size:{size}px;"
        f"text-align:{align};color:#5d5070"
    )

    if info is not None:
        display_info = info
    elif isinstance(feature, list):
        display_info = f'Features: {" & ".join(feature)} {extra}'
    else:
        display_info = f"Feature: {feature} {extra}"

    return f'<p style="{style}"><b> {display_info} </b></p>'


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

from scipy.interpolate import griddata

def metric_resample(s32k_metric,native_lh_path, native_rh_path):
    # Split the metric in L and R
    s32k_metric_L = s32k_metric[0:32492]
    s32k_metric_R = s32k_metric[32492:64984]
    
    # Initialize lists to store resampled metrics
    s5k_metric_L = []
    s5k_metric_R = []
    
    for hemi in ['L', 'R']:
        # Load the fsLR-32k surface
        s32k_surf = nib.load(f'{DATA_PATH}/fsLR-32k.{hemi}.inflated.surf.gii')
        s32k_coords = s32k_surf.darrays[0].data
        
        # Load the low-resolution surface (s5k)
        s5k_surf = nib.load(native_rh_path if hemi == 'R' else native_lh_path)
        s5k_coords = s5k_surf.darrays[0].data
        
        # Interpolate the metric data from s32k to s5k
        if hemi == 'L':
            s5k_metric_L = griddata(s32k_coords, s32k_metric_L, s5k_coords, method='nearest')
        else:
            s5k_metric_R = griddata(s32k_coords, s32k_metric_R, s5k_coords, method='nearest')

    
    return s5k_metric_L, s5k_metric_R

# Load surfaces ----------------------------------------------------------------
def _load_surfaces_ctx(
    subject_dir=None,
    participant_id=None,
    session_id=None,
):
    """Load cortical surfaces - prioritize native subject surfaces if available"""
    
    # Try to load native subject-specific surfaces if subject info is provided
    if subject_dir and participant_id and session_id:
        native_lh_path = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii"
        )
        native_rh_path = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii"
        )

        if os.path.exists(native_lh_path) and os.path.exists(native_rh_path):
            try:
                print(f"Loading native cortical surfaces for {participant_id}/{session_id}")
                inf_lh = read_surface(native_lh_path)
                inf_rh = read_surface(native_rh_path)
                sphere_lh = read_surface(os.path.join(
                        subject_dir, 
                        "structural", 
                        f"{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii"
                    ))
                sphere_rh = read_surface(os.path.join(
                        subject_dir, 
                        "structural", 
                        f"{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii"
                    ))
                inf_lh, inf_rh = inflate_surf(inf_lh, inf_rh, sphere_lh, sphere_rh, W=0.25)

                mask_L = nib.load(
                    os.path.join(
                        subject_dir, 
                        "structural", 
                        f"{participant_id}_{session_id}_hemi-L_surf-fsnative_label-medialwall.label.gii"
                    )
                ).darrays[0].data
                mask_R = nib.load(
                    os.path.join(
                        subject_dir, 
                        "structural", 
                        f"{participant_id}_{session_id}_hemi-R_surf-fsnative_label-medialwall.label.gii"
                    )
                ).darrays[0].data
                mask_L, mask_R = mask_L.astype(bool), mask_R.astype(bool)
                
                return inf_lh, inf_rh, mask_L, mask_R, native_lh_path, native_rh_path
            except Exception as e:
                print(f"Failed to load native surfaces: {e}.")


def _load_surfaces_hip(
    res: Resolution = "high",
    subject_dir=None,
    participant_id=None,
    session_id=None
):
    """Load hippocampal surfaces - prioritize native subject surfaces if available"""
    
    # Map resolution
    if res == "8k":
        res_hip = "8k"
    else:
        res_hip = map_resolution("hippocampus", res)
        
    label = "midthickness"
    
    # Try to load native subject-specific surfaces if subject info is provided
    if subject_dir and participant_id and session_id:
        # Check for native subject-specific hippocampal surfaces for both hemispheres
        # Left hemisphere native files
        native_canonical_path_lh = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-L_space-T1w_den-{res_hip}_label-hipp_{label}.surf.gii"
        )
        native_unfold_path_lh = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-L_space-unfold_den-{res_hip}_label-hipp_{label}.surf.gii"
        )
        
        # Right hemisphere native files
        native_canonical_path_rh = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-R_space-T1w_den-{res_hip}_label-hipp_{label}.surf.gii"
        )
        native_unfold_path_rh = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-R_space-unfold_den-{res_hip}_label-hipp_{label}.surf.gii"
        )
        
        # Check if both hemispheres' files exist
        lh_exists = os.path.exists(native_canonical_path_lh) and os.path.exists(native_unfold_path_lh)
        rh_exists = os.path.exists(native_canonical_path_rh) and os.path.exists(native_unfold_path_rh)
        
        if lh_exists and rh_exists:
            try:
                print(f"Loading native hippocampal surfaces for both hemispheres: {participant_id}/{session_id}")
                
                # Load surfaces for both hemispheres
                mid_lh = read_surface(native_canonical_path_lh)
                unf_lh = read_surface(native_unfold_path_lh)
                mid_rh = read_surface(native_canonical_path_rh)
                unf_rh = read_surface(native_unfold_path_rh)
                
                # Flip right to left surface
                mid_lh.Points[:, 0] *= -1
                unf_lh.Points[:, 0] *= -1

                # Rotate surfaces
                # Rotate around Y-axis 270 for right hemisphere unfolded
                rot_y270 = Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0]))
                unf_rh.Points = rot_y270.apply(unf_rh.Points)
                
                # Rotate around X-axis 90 for left hemisphere unfolded
                rot_y90 = Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0]))
                unf_lh.Points = rot_y90.apply(unf_lh.Points)
                
                # Rotate around Z-axis 180 for right hemisphere unfolded
                rot_z = Rotation.from_rotvec(np.pi * np.array([0, 0, 1]))
                unf_rh.Points = rot_z.apply(unf_rh.Points)
                
                # Create lateral views from canonical
                lat_rh = read_surface(native_canonical_path_rh)
                lat_rh.Points = rot_y270.apply(lat_rh.Points)
                
                lat_lh = read_surface(native_canonical_path_lh)
                lat_lh.Points = rot_y90.apply(mid_lh.Points)
                
                return lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh
            
            except Exception as e:
                print(f"Failed to load native hippocampal surfaces: {e}")
                print("Falling back to template hippocampal surfaces")
        elif rh_exists:
            # Only right hemisphere exists, mirror it to create left
            try:
                print(f"Loading native right hippocampal surface and mirroring for left: {participant_id}/{session_id}")
                
                # Load right hemisphere surfaces
                mid_rh = read_surface(native_canonical_path_rh)
                unf_rh = read_surface(native_unfold_path_rh)
                
                # Create left hemisphere by copying and mirroring right hemisphere
                mid_lh = mid_rh.copy()
                unf_lh = unf_rh.copy()
                
                # Flip right to left surface (mirror along X axis)
                mid_lh.Points[:, 0] *= -1
                unf_lh.Points[:, 0] *= -1
                
                # Apply rotations as before
                rot_y270 = Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0]))
                unf_rh.Points = rot_y270.apply(unf_rh.Points)
                
                rot_y90 = Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0]))
                unf_lh.Points = rot_y90.apply(unf_lh.Points)
                
                rot_z = Rotation.from_rotvec(np.pi * np.array([0, 0, 1]))
                unf_rh.Points = rot_z.apply(unf_rh.Points)
                
                lat_rh = mid_rh.copy()
                lat_rh.Points = rot_y270.apply(lat_rh.Points)
                
                lat_lh = mid_lh.copy()
                lat_lh.Points = rot_y90.apply(mid_lh.Points)
                
                return lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh
                
            except Exception as e:
                print(f"Failed to load and mirror native hippocampal surfaces: {e}")
                print("Falling back to template hippocampal surfaces")
        elif lh_exists:
            # Only left hemisphere exists, mirror it to create right
            try:
                print(f"Loading native left hippocampal surface and mirroring for right: {participant_id}/{session_id}")
                
                # Load left hemisphere surfaces
                mid_lh = read_surface(native_canonical_path_lh)
                unf_lh = read_surface(native_unfold_path_lh)
                
                # Create right hemisphere by copying and mirroring left hemisphere
                mid_rh = mid_lh.copy()
                unf_rh = unf_lh.copy()
                
                # Flip left to right surface (mirror along X axis)
                mid_rh.Points[:, 0] *= -1
                unf_rh.Points[:, 0] *= -1
                
                # Apply rotations as before
                rot_y270 = Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0]))
                unf_rh.Points = rot_y270.apply(unf_rh.Points)
                
                rot_y90 = Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0]))
                unf_lh.Points = rot_y90.apply(unf_lh.Points)
                
                rot_z = Rotation.from_rotvec(np.pi * np.array([0, 0, 1]))
                unf_rh.Points = rot_z.apply(unf_rh.Points)
                
                lat_rh = mid_rh.copy()
                lat_rh.Points = rot_y270.apply(lat_rh.Points)
                
                lat_lh = mid_lh.copy()
                lat_lh.Points = rot_y90.apply(mid_lh.Points)
                
                return lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh
                
            except Exception as e:
                print(f"Failed to load and mirror native hippocampal surfaces: {e}")
                print("Falling back to template hippocampal surfaces")
    
    # Fall back to template surfaces if native surfaces aren't available or failed to load
    print("Using template hippocampal surfaces")
    
    # Template paths
    pth_canonical = f"{DATA_PATH}/tpl-avg_space-canonical_den-{res_hip}_label-hipp_{label}.surf.gii"
    pth_unfold = f"{DATA_PATH}/tpl-avg_space-unfold_den-{res_hip}_label-hipp_{label}.surf.gii"
    
    # Load surfaces
    try:
        mid_rh = read_surface(pth_canonical)
        unf_rh = read_surface(pth_unfold)
        
        # Create left hemisphere by copying and mirroring right hemisphere
        mid_lh = mid_rh.copy()
        unf_lh = unf_rh.copy()
        
        # Flip right to left surface
        mid_lh.Points[:, 0] *= -1
        unf_lh.Points[:, 0] *= -1
        
        # Rotate surfaces
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
        lat_rh = mid_rh.copy()
        lat_rh.Points = rot_y270.apply(lat_rh.Points)
        
        # Left Antero-posterior lateral
        lat_lh = mid_lh.copy()
        lat_lh.Points = rot_y90.apply(mid_lh.Points)
        
        return lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh
    
    except Exception as e:
        print(f"Error loading hippocampal surfaces: {e}")
        raise


# Load data --------------------------------------------------------------------
def _load_data_sctx(
    sctx_file: PathType,
    analysis: Analysis,
    threshold: Union[float, None] = None,
    threshold_alpha: float = 0.5,
):

    # Load CSV data
    x = pd.read_csv(sctx_file, header=[0], index_col=None).to_numpy().ravel()
    if threshold is not None:
        x[np.abs(x) < threshold] *= threshold_alpha
    # if analysis == "asymmetry":
    #     print("e")
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

    # Fill ventricle slots with 0 (or another default) instead of leaving as NaN
    array_16[7] = 0.0   # Left ventricle
    array_16[15] = 0.0  # Right ventricle

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
    subject_dir=None,
    participant_id=None,
    session_id=None,
    env=None,
    tmp_dir=None
):

    slh, srh, mask_L, mask_R, native_lh_path, native_rh_path = _load_surfaces_ctx(
        subject_dir=subject_dir,
        participant_id=participant_id,
        session_id=session_id,
    )


    feat_map_L, feat_map_R = project_to_native_surface(
        data_lh,
        data_rh, 
        native_lh_path, 
        native_rh_path,
        env,
        DATA_PATH=DATA_PATH,
        subject_dir=subject_dir,
        participant_id=participant_id,
        session_id=session_id,
        tmp_dir=tmp_dir
    ) 
    feat_map_L[mask_L] = np.nan
    feat_map_R[mask_R] = np.nan

    feat_map = np.hstack((feat_map_L, feat_map_R))

    kwds = dict(
        cmap=cmap,
        color_bar=color_bar,
        color_range=color_range,
        transparent_bg=False,
        nan_color=(211, 211, 211, 1),
        embed_nb=True,
        interactive=False,
        screenshot=True,
    )

    # Plot the mean FEATURE on fsLR-32k
    out_png = Path(out_png)
    if analysis == "asymmetry":
        plot_hemisphere_lh(
            slh,
            array_name=feat_map_L,
            layout_style="grid",
            label_text={"left": ["Lateral", "Medial"], "top": [""]},
            share=None,
            zoom=1.1,
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
            zoom=1.1,
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
            values=[feat_map_L, feat_map_L],
            views=["dorsal", "ventral"],
            label_text={
                "bottom": [
                    "Left Superior",
                    "Left Inferior",
                ]
            },
            # share='both', zoom=2.95, size=(550, 350),
            share="both",
            zoom=1.1,
            size=(550, 450),
            filename=out_png2,
            **kwds,
        )

    else:
        plot_surfs(
            surfaces=[slh, srh, srh, slh],
            values=[feat_map_L, feat_map_R, feat_map_R, feat_map_L],
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
            zoom=2.5,
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
            nan_color=(211, 211, 211, 1),
            zoom=1.2,
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
            nan_color=(211, 211, 211, 1),
            zoom=1.2,
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
    subject_dir=None,
    participant_id=None,
    session_id=None
):

    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip(res=res,    
                                                                        subject_dir=subject_dir,
                                                                        participant_id=participant_id,
                                                                        session_id=session_id)
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
    subject_dir=None,
    participant_id=None,
    session_id=None,
    env=None,
    tmp_dir=None
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
            subject_dir=subject_dir,
            participant_id=participant_id,
            session_id=session_id,
            env=env,
            tmp_dir=tmp_dir
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
            subject_dir=subject_dir,
            participant_id=participant_id,
            session_id=session_id,
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
        f'<p style="margin-bottom:0;margin-top:0;font-family=gill sans,'
        f'sans-serif;text-align:center;font-size:14px;color=#ffb311"> '
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
    feature_means=None,
    subject_dir=None,
    env=None,
):
    """
    Generate report section for a specific brain structure adapted for zbdataset structure.
    """
    import logging
    import uuid
    import os
    import numpy as np
    from pathlib import Path
    
    logger = logging.getLogger(tmp_dir)
    bids_id = f"{sid}_{ses}" if ses else sid
    
    # Handle feature name
    if isinstance(feat, list):
        feat = "-".join(feat)
    
    # Map structure names to directory names used by zbdataset
    struct_dir_map = {
        'cortex': 'cortex',
        'subcortex': 'subcortical', 
        'hippocampus': 'hippocampus'
    }
    
    struct_dir = struct_dir_map.get(struct, struct)
    
    # Build file paths based on structure type and zbdataset naming conventions
    file_lh = None
    file_rh = None

    if struct == "subcortex":
        # For subcortical, look for CSV files with analysis type
        file_lh = os.path.join(path_analysis, struct_dir, f"{bids_id}_feature-{feat}_analysis-{analysis}.csv")
    else:
        # For cortical and hippocampal, build surface file paths
        if struct == "cortex":
            # Map resolution for cortical files
            res_mapped = "32k" if res == "high" else "5k"
            
            # Cortical files - always include analysis type
            file_lh = os.path.join(
                path_analysis, struct_dir,
                f"{bids_id}_hemi-L_surf-fsLR-{res_mapped}_label-{label}_feature-{feat}_smooth-{smooth}mm_analysis-{analysis}.func.gii"
            )
            
            if analysis != "asymmetry":
                file_rh = os.path.join(
                    path_analysis, struct_dir,
                    f"{bids_id}_hemi-R_surf-fsLR-{res_mapped}_label-{label}_feature-{feat}_smooth-{smooth}mm_analysis-{analysis}.func.gii"
                )
            # For asymmetry, we only need the left hemisphere file
            
        elif struct == "hippocampus":
            # Map resolution for hippocampal files
            if res == "8k":
                res_mapped = "8k"
                suffix = ""
            else:
                res_mapped = "0p5" if res == "high" else "1p0"
                suffix = "mm"
            
            # Hippocampal files - always include analysis type
            file_lh = os.path.join(
                path_analysis, struct_dir,
                f"{sid}_{ses}_hemi-L_den-{res_mapped}{suffix}_label-hipp_{label}_feature-{feat}_smooth-{smooth}mm_analysis-{analysis}.func.gii"
            )
            
            if analysis != "asymmetry":
                file_rh = os.path.join(
                    path_analysis, struct_dir,
                    f"{sid}_{ses}_hemi-R_den-{res_mapped}{suffix}_label-hipp_{label}_feature-{feat}_smooth-{smooth}mm_analysis-{analysis}.func.gii"
                )
            # For asymmetry, we only need the left hemisphere file
    
    # Generate info string
    thr_str = "" if thr is None else f" | threshold={thr}"
    info = ""
    if struct != "subcortex":
        if struct == "cortex":
            res_display = "32k" if res == "high" else "5k"
        else:
            if res == "8k":
                res_display = "8k"
            else:
                res_display = "0.5mm" if res == "high" else "1.0mm"
        info = f"| {smooth}mm smooth | resolution {res_display} "
    
    # Check if files exist
    lh_exists = file_lh and os.path.exists(file_lh)
    rh_exists = file_rh and os.path.exists(file_rh) if file_rh else True
    
    if not lh_exists or not rh_exists:
        missing_file = file_rh if lh_exists else file_lh
        logger.warning(f"{struct.capitalize()} file was not found: {missing_file}")
        
        png_block = make_png_missing(struct)
        html = (
            '<p style="margin-bottom:0;margin-top:0;'
            "font-family:gill sans,sans-serif;text-align:left;font-size:14px;"
            'color:#5d5070"> '
            f"<b>{struct.capitalize()}</b> | {approach} approach {info}"
            "</p>"
        )
        html += png_block
        return html
    
    # Load and process data
    try:
        feat_lh, feat_rh = load_data_struct(
            struct,
            file_lh=file_lh,
            file_rh=file_rh,
            analysis=analysis,
            threshold=thr,
            threshold_alpha=thr_alpha,
        )
        
        # Calculate means for display
        mean_str = ""
        if feat_lh is not None:
            mean_lh = np.mean(feat_lh)
            mean_str += f"| left mean={mean_lh:.2f} "
        if feat_rh is not None:
            mean_rh = np.mean(feat_rh)
            mean_str += f"| right mean={mean_rh:.2f}"
        
        # Generate HTML info
        html = (
            '<p style="margin-bottom:0;margin-top:0;'
            "font-family:gill sans,sans-serif;text-align:left;font-size:14px;"
            'color:#5d5070"> '
            f"<b>{struct.capitalize()}</b> | {approach} approach {info}{mean_str}"
            "</p>"
        )
        
        # Adjust color range for asymmetry
        if analysis == "asymmetry":
            color_range = (-1.5, 1.5)
        
        # Generate output PNG
        out_png = os.path.join(
            tmp_dir, 
            f"{bids_id}_{struct}_feature-{feat}_analysis-{analysis}{thr}_{uuid.uuid4()}.png"
        )
        
        # Create PNG visualization
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
            subject_dir=subject_dir,
            participant_id=sid,
            session_id=ses,
            env=env,
            tmp_dir=tmp_dir
        )
        
        html += png_block
        return html
        
    except Exception as e:
        logger.error(f"Error processing {struct} data: {e}")
        png_block = make_png_missing(struct)
        html = (
            '<p style="margin-bottom:0;margin-top:0;'
            "font-family:gill sans,sans-serif;text-align:left;font-size:14px;"
            'color:#5d5070"> '
            f"<b>{struct.capitalize()}</b> | {approach} approach {info} | ERROR: {e}"
            "</p>"
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
    approach: Approach = "wscore",
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
    color_range=(-3, 3),
    cmap_asymmetry="cmo.balance_r",
    tmp_dir: PathType = "/tmp",
    subject_dir=None,
    output_dir=None,
    tag=None,
    feature_means=None,
    env=None,
    verbose=True,
    control_data=None,
    valid_subjects=None,
):
    """Zbrains: Clinical report generator adapted for zbdataset

    Parameters
    ----------
    zbrains_path:
        Output directory for processed z-brains derivatives. (Legacy parameter)
    sid:
        Participant id
    ses:
        Session name
    age:
        Subject age
    sex:
        Subject sex
    approach:
        Comparison approach ('zscore' or 'wscore').
    features:
        List of features to include in report
    analyses:
        List of analyses to include in report
    cmap:
        Colormap for plotting the regional changes.
    cmap_asymmetry:
        Colormap for plotting the asymmetry changes.
    color_range:
        Color range of values for the regional plot
    color_bar:
        Position of the color bar for hippocampal and cortical plots
    res_ctx:
        Surface resolution for cortex.
    res_hip:
        Surface resolution for hippocampus.
    label_ctx:
        Cortical surface used when mapping the data from volume to surface.
    label_hip:
        Hippocampal surface used when mapping the data from volume to surface.
    threshold:
        Threshold for the maps
    threshold_alpha:
        Alpha channel for the colormap when threshold is applied
    smooth_ctx:
        Smoothing for the cortical data in mm.
    smooth_hip:
        Smoothing for the hippocampal data in mm.
    tmp_dir:
        Working directory for temporary files
    subject_dir:
        Subject-specific directory containing analysis results
    output_dir:
        Output directory for the report
    tag:
        Tag for the output filename
    verbose:
        Verbose output

    Returns
    -------
    str
        Path to generated PDF report
    """
    import logging
    import itertools
    import os
    import platform
    from pathlib import Path
    
    # Set up logging
    logger = logging.getLogger(tmp_dir)
    
    # Convert approach to folder name
    approach_folder = f"{approach}_maps"
    
    # Get BIDS ID
    bids_id = f"{sid}_{ses}" if ses else sid
    
    # Set up paths - use subject_dir if provided, otherwise use zbrains_path
    if subject_dir is not None:
        path_analysis = os.path.join(subject_dir, approach_folder)
    elif zbrains_path is not None:
        subject_dir = get_subject_dir(zbrains_path, sid, ses)
        path_analysis = f"{subject_dir}/{approach_folder}"
    else:
        raise ValueError("Either subject_dir or zbrains_path must be provided")
    
    if not os.path.exists(path_analysis):
        raise ValueError(f"Analysis results not found at {path_analysis}")
    
    # List all files of approach
    subses_files = []
    for root, dirs, files in os.walk(path_analysis):
        subses_files.extend([os.path.join(root, f) for f in files])
    
    # If ANALYSIS is empty run all available analysis
    if analyses is None:
        available_analyses = sorted(
            list(
                set(
                    [
                        os.path.basename(file).split(".")[0].split("analysis-")[1].split("_")[0]
                        for file in subses_files
                        if "analysis-" in os.path.basename(file)
                    ]
                )
            )
        )
        analyses = [a for a in LIST_ANALYSES if a in available_analyses]
        if not analyses:
            analyses = ['regional']  # Default fallback
    
    # If FEATURES is empty run all available features
    if features is None:
        features = sorted(
            list(
                set([
                    os.path.basename(file).split("feature-")[1].split("_")[0] 
                    for file in subses_files
                    if "feature-" in os.path.basename(file)
                ])
            )
        )
    
    # Remove volume from the features if present
    if "volume" in features:
        features.remove("volume")
    
    # Set up display for headless plotting
    display_flag = False
    if (
        "DISPLAY" not in os.environ or not os.environ["DISPLAY"]
    ) and platform.system() != "Windows":
        try:
            from pyvirtualdisplay import Display
            os.environ["PYVIRTUALDISPLAY_DISPLAYFD"] = "0"
            dsize = (900, 750)
            display = Display(visible=False, size=dsize)
            display.start()
            display_flag = True
        except ImportError:
            if verbose:
                print("Warning: pyvirtualdisplay not available, plotting may fail in headless environment")
    
    # Generate report content
    report = ""
    
    # Generate first page with overview
    try:
        # First page header with just "Clinical Report" as title
        report += report_header_template(
            sid=sid, 
            ses=ses, 
            age=age, 
            sex=sex, 
            analysis='Control statistics'
        )
        
        # Features table section
        report += feature_header_template(align="left", info="Subject feature data", size=14)
        if valid_subjects:
            report += features_table(sid, ses, valid_subjects)
        else:
            # If no zbrains_path provided, just note what features will be in the report
            feature_list = ", ".join(features)
            report += f"""
            <p style="font-family:gill sans,sans-serif;text-align:center;">
            The following features will be analyzed: {feature_list}
            </p>
            """
        
        # Control database section
        report += feature_header_template(align="left", info="Control database", size=14)
        if control_data is not None:
            report += controls_summary_html(control_data, subject_age=age)
        else:
            report += """
            <p style="font-family:gill sans,sans-serif;text-align:center;">
            No control database information available.
            </p>
            """
        
        # Pipeline information section
        report += feature_header_template(align="left", info="Pipeline details", size=12)
        report += pipeline_info(tmp_dir)
        
        # Add page break after first page
        report += '<div style="page-break-after: always;"></div>'
    except Exception as e:
        # If there's an error creating the first page, log it but continue with analysis
        logger.warning(f"Error generating first page: {e}")

    # Function to check if files exist for a given feature/analysis combination
    def check_feature_exists(feat, analysis):
        """Check if required files exist for a feature/analysis combination"""
        feat_sctx = "volume" if feat == "thickness" else feat
        
        # Check cortex files
        res_mapped_ctx = "32k" if res_ctx == "high" else "5k"
        file_ctx_lh = os.path.join(
            path_analysis, "cortex",
            f"{bids_id}_hemi-L_surf-fsLR-{res_mapped_ctx}_label-{label_ctx}_feature-{feat}_smooth-{smooth_ctx}mm_analysis-{analysis}.func.gii"
        )
        
        # Check hippocampus files
        if res_hip == "8k":
            res_mapped_hip = "8k"
            suffix = ""
        else:
            res_mapped_hip = "0p5" if res_hip == "high" else "1p0"
            suffix = "mm"

        file_hip_lh = os.path.join(
            path_analysis, "hippocampus",
            f"{sid}_{ses}_hemi-L_den-{res_mapped_hip}{suffix}_label-hipp_{label_hip}_feature-{feat}_smooth-{smooth_hip}mm_analysis-{analysis}.func.gii"
        )
        
        # Check subcortical files
        file_sctx = os.path.join(
            path_analysis, "subcortical",
            f"{bids_id}_feature-{feat_sctx}_analysis-{analysis}.csv"
        )
        
        # Return True if at least one structure has data for this feature/analysis
        return (os.path.exists(file_ctx_lh) or 
                os.path.exists(file_hip_lh) or 
                os.path.exists(file_sctx))
    
    # Fix fMRI features if present:
    fmri_features = ['rmssd', 'timescales', 'alff', 'falff']
    if 'fMRI' in features:
        features.remove('fMRI')
        features.extend(fmri_features)

    # Filter out features/analyses combinations that don't have any data
    valid_combinations = []
    for analysis, feat in itertools.product(analyses, features):
        if check_feature_exists(feat, analysis):
            valid_combinations.append((analysis, feat))
        elif verbose:
            logger.info(f"Skipping {analysis}/{feat} - no data files found")
    
    if not valid_combinations:
        raise ValueError("No valid feature/analysis combinations found with existing data files")

    try:
        # Process each valid analysis and feature combination
        for analysis, feat in valid_combinations:
            for thresh in [None, threshold]:
                if thresh is None and verbose:
                    logger.info(f"Generating report for analysis={analysis}, approach={approach}, feature={feat}")
                
                # Handle feature name for subcortical (volume instead of thickness)
                if isinstance(feat, list):
                    feat_sctx = ["volume" if f == "thickness" else f for f in feat]
                else:
                    feat_sctx = "volume" if feat == "thickness" else feat
                
                # Generate report header
                report += report_header_template(
                    sid=sid, 
                    ses=ses, 
                    age=age, 
                    sex=sex, 
                    analysis=analysis
                )
                
                # Generate feature header
                extra = "" if thresh is None else f" | Threshold: {thresh}"
                report += feature_header_template(feat, extra=extra)
                
                # Common kwargs for report generation
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
                    subject_dir=subject_dir,
                )
                
                # Generate cortical report section
                report += report_struct(
                    struct="cortex",
                    feat=feat,
                    res=res_ctx,
                    label=label_ctx,
                    smooth=smooth_ctx,
                    feature_means=feature_means,
                    env=env,
                    **kwds,
                )
                
                # Generate subcortical report section
                report += report_struct(
                    struct="subcortex", 
                    feat=feat_sctx, 
                    feature_means=feature_means,
                    **kwds
                )
                
                # Generate hippocampal report section
                report += report_struct(
                    struct="hippocampus",
                    feat=feat,
                    res=res_hip,
                    label=label_hip,
                    smooth=smooth_hip,
                    feature_means=feature_means,
                    **kwds,
                )
                
                # Add color legend
                report += report_colors(analysis=analysis)
                
                # Add page break
                report += '<div style="page-break-after: always;"></div>'
        
        # Generate output filename
        if output_dir and tag:
            file_pdf = os.path.join(output_dir, f"{tag}.pdf")
            if os.path.exists(file_pdf):
                os.remove(file_pdf)
        else:
            file_pdf = os.path.join(subject_dir, f"{bids_id}_approach-{approach}_summary-report.pdf")
        
        # Convert HTML to PDF
        convert_html_to_pdf(report, file_pdf)
        
        if verbose:
            logger.info(f"Clinical report successfully created: {file_pdf}")
        
        return os.path.realpath(file_pdf)
        
    finally:
        # Clean up display
        if display_flag:
            display.stop()
            del display
