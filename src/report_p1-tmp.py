import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
import re
import os
import glob
import base64
from datetime import date
from typing import List, Union, Optional
from pathlib import Path
from xhtml2pdf import pisa

import uuid
import time
import random
import string
import segno

# -----------------------------------------------------------------------------
# global variables
try:
    DATA_PATH = Path(__file__).resolve().parent.parent / "data"
except NameError:
    DATA_PATH = "/host/yeatman/local_raid/rcruces/git_here/z-brains/data"
    
Analysis="Regional"
sid="HC062"
ses="01" 
age="37" 
sex="N"
analysis=Analysis
df = pd.DataFrame({
    'sex': np.random.choice(['F', 'M'], size=500),
    'age': np.random.normal(loc=30, scale=10, size=500).clip(0, 90)
})
PathType = "/tmp"
ZBRAINS_DIR="/data/mica3/BIDS_MICs/derivatives/zbrains_ianjudy_test"

def features_table(sid, ses, ZBRAINS_DIR):
    def extract_features_from_glob(pattern):
        feature_pattern = re.compile(r'feature-([a-zA-Z0-9]+)')
        features = []
    
        for file_path in glob.glob(pattern):
            if os.path.isfile(file_path):
                # Get just the base filename
                base_name = os.path.basename(file_path).split('.')[0]
    
                # Extract the feature using regex
                match = feature_pattern.search(base_name)
                if match:
                    features.append(match.group(1))
    
        return sorted(set(features))
    
    # Function to check both hemispheres
    def both_hemispheres_exist(base_path_fn, feature):
        return all(os.path.isfile(base_path_fn(hemi, feature)) for hemi in ["L", "R"])
    
    # Path generators
    def cortex_path(hemi, feature):
        return f"{ZBRAINS_DIR}/sub-{sid}/ses-{ses}/maps/cortex/sub-{sid}_ses-{ses}_hemi-{hemi}_surf-fsLR-32k_label-midthickness_feature-{feature}_smooth-10mm.func.gii"
    
    def hippocampus_path(hemi, feature):
        return f"{ZBRAINS_DIR}/sub-{sid}/ses-{ses}/maps/hippocampus/sub-{sid}_ses-{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-5mm.func.gii"
    
    def subcortex_path(feature):
        sub_feature = "volume" if feature == "thickness" else feature
        return f"{ZBRAINS_DIR}/sub-{sid}/ses-{ses}/maps/subcortex/sub-{sid}_ses-{ses}_feature-{sub_feature}.csv"
    
    # Generate HTML table
    html = """
    <style>
        table.features {
            width: 100%;
            border-collapse: collapse;
            font-family: "Gill Sans", sans-serif;
            font-size: 14px;
        }
        table.features th {
            background-color: #e3d7f4;
            color: #505050;
            padding: 8px;
            text-align: center;
        }
        table.features td {
            padding: 8px;
            border: 1px solid #ddd;
            text-align: center;
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
    
    # Get all features
    all_features = extract_features_from_glob(f'{ZBRAINS_DIR}/sub-{sid}/ses-{ses}/maps/*/*')
    
    CHECK = "✔️"
    CROSS = "❌"
    
    for feature in sorted(all_features):
        ctx_check = CHECK if both_hemispheres_exist(cortex_path, feature) else CROSS
        hip_check = CHECK if both_hemispheres_exist(hippocampus_path, feature) else CROSS
        sub_check = CHECK if os.path.isfile(subcortex_path(feature)) else CROSS
    
        html += f"""
        <tr>
            <td>{feature}</td>
            <td>{ctx_check}</td>
            <td>{hip_check}</td>
            <td>{sub_check}</td>
        </tr>
        """
    
    html += "</table>"
    
    return(html)

# -----------------------------------------------------------------------------
# Functions
def report_header_template(
    *,
    sid: str,
    ses: Union[str, None] = None,
    age: Union[float, None] = None,
    sex: Union[str, None] = None,
    analysis: Analysis,
    title=None,
):

    # if ses is None:
    #     ses = ''
    if title is None:
        title = f"Clinical Summary &nbsp; | &nbsp; <b> {analysis.capitalize()} analysis </b> "

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

    # Responsive banner and spacing
    report_header = (
        f'<div style="margin-bottom:30px;">'
        f'<img id="top" src="{DATA_PATH}/zbrains_banner.png" alt="zbrains" '
        f'style="display:block;margin-left:auto;margin-right:auto;width:100%;height:auto;margin-bottom:10px;">'
        f'<p style="{style.format(fontsize=20, extra="margin-bottom:10px;")}">'
        f"{title}"
        f"</p>"
        f'<p style="{style.format(fontsize=14, extra="margin-top:-10px;margin-bottom:20px;")}">'
        f"<b>Subject</b>: {sid},{ses_str} "
        f"&nbsp; <b>Sex</b>: {sex}, &nbsp; <b>Age</b>: {age}"
        f"</p>"
        f'</div>'
    )

    return report_header

def feature_header_template(
    feature: Union[str, List[str]] = "",
    extra: str = "",
    align: str = "center",
    size= 14,
    info: Optional[str] = None
) -> str:
    """
    Generates an HTML header for a feature or list of features, with configurable alignment
    and the option to provide a custom info string.

    Args:
        feature (Union[str, List[str]], optional): A single feature or a list of features.
        extra (str, optional): Additional text to append. Defaults to "".
        align (str, optional): Text alignment ('left', 'center', etc.). Defaults to "center".
        info (str, optional): Custom info string to display. If provided, overrides feature/extra.

    Returns:
        str: HTML string for the feature header.
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

def uuid_qr(qr_path="/tmp", scale=5, border=0):
    # Step 1: Generate a UUID using "zbrain", timestamp, and a random string
    timestamp = str(int(time.time()))
    rand_str = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    base_string = f"zbrain-{timestamp}-{rand_str}"
    uuid_generated = uuid.uuid5(uuid.NAMESPACE_DNS, base_string)
    print("Generated UUID:", uuid_generated)
    
    # Step 2: Generate QR code using segno
    qr = segno.make(str(uuid_generated))
    
    # Step 3: Save it as PNG with foreground color #A569BD and white background
    qr.save(
        f'{qr_path}/{uuid_generated}.png',
        scale=scale,
        border=border,
        dark='#A569BD',
        light='white'
    )
    return(uuid_generated)

def pipeline_info():
    """
    Generates an HTML table with project info, styled for tight left alignment.

    Args:
        workstation (str): Workstation name or identifier.
        user (str, optional): Username to display. If None, uses os.getlogin().

    Returns:
        str: HTML string for the info table.
    """
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
    tmp_dir='/home/bic/rcruces/Desktop'
    uuid_generated = uuid_qr(tmp_dir, scale=5, border=0)
    
    # Table HTML
    info_table = f'''
    <div style="display:flex; align-items:flex-start; font-family:gill sans, sans-serif; gap:40px;">
    
      <!-- Left: Info Table -->
      <div>
        <table style="border-collapse:collapse; font-size:11px;">
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
      </div>
    
      <!-- Right: QR Code Image -->
      <div style="flex-shrink:0;">
        <img src="{tmp_dir}/{uuid_generated}.png" alt="QR Code" style="width:150px; height:auto; border:1px solid #ccc;" />
      </div>
    
    </div>
    '''

    return info_table



def controls_summary_html(df):    
    # --- Compute statistics ---
    mean_age = df['age'].mean()
    sd_age = df['age'].std()
    Ntotal = len(df)
    Fnumber = (df['sex'] == 'F').sum()
    Mnumber = (df['sex'] == 'M').sum()
    Fperc = 100 * Fnumber / Ntotal if Ntotal else 0
    Mperc = 100 * Mnumber / Ntotal if Ntotal else 0

    # --- Bin and count ---
    bins = np.arange(df['age'].min() // 5 * 5, df['age'].max() // 5 * 5 + 10, 5)
    labels = [f"{int(b)}-{int(b + 4)}" for b in bins[:-1]]
    df['age_bin'] = pd.cut(df['age'], bins=bins, labels=labels, right=True, include_lowest=True)
    age_bin_counts = df.groupby(['age_bin', 'sex']).size().unstack(fill_value=0)
    
    # Ensure both columns exist
    for sex in ['M', 'F']:
        if sex not in age_bin_counts.columns:
            age_bin_counts[sex] = 0
    age_bin_counts = age_bin_counts[['M', 'F']]
    
    # --- Plot ---
    fig, ax = plt.subplots(figsize=(5, 2), dpi=150)
    ax.bar(age_bin_counts.index, age_bin_counts['M'], color="#34495E", label='M', width=0.9)
    ax.bar(age_bin_counts.index, age_bin_counts['F'], bottom=age_bin_counts['M'], color="#839192", label='F', width=0.9)
    
    ax.set_xlabel('Age')
    ax.set_ylabel('Count')
    plt.xticks(rotation=45)
    
    # --- Vertical line at correct bin position ---
    try:
        subject_age = float(age)
        subject_bin = pd.cut([subject_age], bins=bins, labels=labels, right=True, include_lowest=True)[0]
        if pd.notna(subject_bin):
            # Find the center x-position of the categorical bin
            #bin_index = labels.index(subject_bin)
            x_tick_labels = list(age_bin_counts.index)
            x_position = x_tick_labels.index(subject_bin)
            ax.axvline(x=x_position, color="#A93226", linestyle='-', linewidth=2)  # no label
    except Exception as e:
        print("Error drawing subject age line:", e)
    
    # --- Custom legend: Only F (top) and M ---
    handles, labels_ = ax.get_legend_handles_labels()
    label_order = ['F', 'M']  # only bars, no Subject Age
    ordered_handles = [handles[labels_.index(lbl)] for lbl in label_order if lbl in labels_]
    ax.legend(ordered_handles, label_order)
    
    # --- Clean up ---
    plt.tight_layout()
    for spine in ['top', 'right', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)
    
    # --- Save to buffer as base64 PNG ---
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    img = base64.b64encode(buf.getvalue()).decode('utf-8')
    buf.close()

    # --- HTML layout (no border on table) ---
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

# Functions for report
def convert_html_to_pdf(source_html, output_filename: PathType):
    with open(output_filename, "w+b") as result_file:
        pisa_status = pisa.CreatePDF(source_html, dest=result_file)

    # True on success and False on errors
    return pisa_status.err


# -----------------------------------------------------------------------------
#               REPORT

# Generate HTML header
today_date = date.today().strftime("%B %d, %Y")
report = f"""
<html>
<head>
<style>
@page {{
        @bottom-center {{
            content: "Date: {today_date}";
            font-size: 10px;
            color: #888;
        }}
    }}
    </style>
    </head>
"""

# Generate FIRST page report header
report = report_header_template(
                sid=sid, 
                ses=ses, 
                age=age, 
                sex=sex, 
                analysis=analysis,
                title=f"Clinical Report</b>"
            )

# SUBTITLE section: Control database
report += feature_header_template(align="left", info="Subject feature data", size=14)
report += features_table(sid, ses, ZBRAINS_DIR)
        
# SUBTITLE section: Control database
report += feature_header_template(align="left", info="Control database", size=14)
report += controls_summary_html(df)

# SUBTITLE section: Pipeline information
report += feature_header_template(align="left", info="Pipeline details", size=12)

### Add pipeline details
report += pipeline_info()

# Add page break
report += '<div style="page-break-after: always;"></div>'

# CLOSE
report += "</html>"


# -----------------------------------------------------------------------------
# SAVE HTML
with open("/home/bic/rcruces/Desktop/report_p1.html", "w") as f:
    f.write(report)
    

# -----------------------------------------------------------------------------
# SAVE PDF
convert_html_to_pdf(report, '/home/bic/rcruces/Desktop/report_p1.pdf')



