import os
from .utilities import (
    show_title,
    show_note,
    show_info,
    show_warning,
    add_field_to_xml,
)
import subprocess
from .constants import *
import time
import glob
from typing import Union
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import shutil
from .blurring import compute_blurring


def subcortical_mapping(
    subject_id: str,
    image: Union[str, Path, None] = None,
    seg: Union[str, Path, None] = None,
    vol: Union[str, Path, None] = None,
    output: Union[str, Path, None] = None,
    include_icv: bool = False,
):
    LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
    STRUCTURES = [
        "Laccumb",
        "Lamyg",
        "Lcaud",
        "Lhippo",
        "Lpal",
        "Lput",
        "Lthal",
        "Raccumb",
        "Ramyg",
        "Rcaud",
        "Rhippo",
        "Rpal",
        "Rput",
        "Rthal",
    ]
    if seg and image and vol is None:
        seg = np.asanyarray(nib.load(seg).dataobj)
        img = np.asanyarray(nib.load(image).dataobj)
        cols = STRUCTURES + ["ICV"] if include_icv else STRUCTURES
        df = pd.DataFrame(
            np.nan, columns=cols, index=pd.Index([subject_id], name="SubjID")
        )
        for i, k in enumerate(LABELS):
            df.loc[subject_id, STRUCTURES[i]] = img[seg == k].mean()

    elif seg is None and image is None and vol:
        df_volumes = pd.read_csv(
            vol,
            comment="#",
            index_col=0,
            sep=r"\s+",
            header=None,
            usecols=[1, 3, 4],
            names=["label", "volume", "structure"],
        )
        volumes = df_volumes.loc[LABELS].volume.to_numpy()
        cols = STRUCTURES
        if include_icv:
            icv = float(
                next(
                    (
                        line.split(",")[-2].strip()
                        for line in open(vol)
                        if line.startswith("# Measure EstimatedTotalIntraCranialVol")
                    ),
                    None,
                )  # type: ignore
            )
            volumes = np.r_[volumes, icv]
            cols = STRUCTURES + ["ICV"]

        df = pd.DataFrame(
            volumes[None], columns=cols, index=pd.Index([subject_id], name="SubjID")
        )

    else:
        raise ValueError(
            "Ambiguous inputs. Please provide either 'image' and 'seg' together or 'vol' alone."
        )

    df.to_csv(output)


def map_subcortex(
    bids_id,
    feat,
    subject_surf_dir,
    subject_micapipe_dir,
    subject_output_dir,
    folder_maps,
    folder_sctx,
    script_dir,
    subject_plugin_dir,
):
    if "blur" in feat:
        return
    show_info(f"{bids_id}: Mapping '{feat}' to subcortical structures")
    map_input = {
        "flair": "map-flair",
        "adc": "model-DTI_map-ADC",
        "fa": "model-DTI_map-FA",
        "qt1": "map-T1map",
    }
    # Input & output locations
    aseg_stats_file = os.path.join(subject_surf_dir, "stats", "aseg.stats")
    seg_file = os.path.join(
        subject_micapipe_dir,
        "parc",
        f"{bids_id}_space-nativepro_T1w_atlas-subcortical.nii.gz",
    )
    input_dir = os.path.join(subject_micapipe_dir, "maps")
    # If feat starts with 'plugin-', update input_dir and remove 'plugin_' from feat
    if feat.startswith("plugin_"):
        input_dir = os.path.join(subject_plugin_dir, "maps")
        feat = feat[7:]
    feat_lower = feat.lower()
    if feat_lower in map_input:
        input_feat = map_input[feat_lower]
        map_output = {
            "volume": "volume",
            "flair": "flair",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "T1map",
        }

        output_feat = map_output[feat_lower]
    else:
        input_feat = feat
        output_feat = feat
    output_dir = os.path.join(subject_output_dir, folder_maps, folder_sctx)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{bids_id}_feature-{output_feat}.csv")

    # check that input files exist & if not volume
    if feat != "volume":
        # Mappings from features names to the way they appear in the input and output filenames
        input_file = os.path.join(
            input_dir, f"{bids_id}_space-nativepro_{input_feat}.nii.gz"
        )

        for file in [seg_file, input_file]:
            if not os.path.isfile(file):

                show_warning(
                    f"{bids_id}: cannot map '{feat}' to subcortical structures. Missing file: {file}"
                )
                return
        subcortical_mapping(bids_id, image=input_file, seg=seg_file, output=output_file)
        # subprocess.run(["python", os.path.join(script_dir, "subcortical_mapping.py"), "-id", bids_id, "-i", input_file, "-s", seg_file, "-o", output_file])

    else:
        if not os.path.isfile(aseg_stats_file):
            show_warning(
                f"{bids_id}: cannot map '{feat}' to subcortical structures. Missing file {aseg_stats_file}"
            )
            return
        subcortical_mapping(bids_id, vol=aseg_stats_file, output=output_file)
        # subprocess.run(["python", os.path.join(script_dir, "subcortical_mapping.py"), "-id", bids_id, "-v", aseg_stats_file, "-o", output_file])

    if os.path.isfile(output_file):
        show_note(f"{bids_id}: '{feat}' successfully mapped.")
    else:
        show_warning(f"{bids_id}: could not map '{feat}' to subcortical structures.")


def map_cortex(
    bids_id,
    feat,
    resol,
    label: str,
    fwhm,
    workbench_path: str,
    subject_micapipe_dir: str,
    subject_output_dir: str,
    folder_maps,
    folder_ctx,
    subject_plugin_dir: str,
    tmp_dir,
):
    show_info(
        f"{bids_id}: Mapping '{feat}' to cortex [label={label}, resolution={resol}]"
    )
    root_dir = subject_micapipe_dir
    if feat.startswith("plugin_"):
        root_dir = subject_plugin_dir
        feat = feat[7:]
    # Input & output locations
    surf_dir = os.path.join(root_dir, "surf")
    input_dir = os.path.join(root_dir, "maps")
    output_dir = os.path.join(subject_output_dir, folder_maps, folder_ctx)
    os.makedirs(output_dir, exist_ok=True)

    # Mappings from features names to the way they appear in the input and output filenames
    map_feat = {
        "thickness": "thickness",
        "flair": "flair",
        "adc": "ADC",
        "fa": "FA",
        "qt1": "T1map",
    }

    # feat name in filenames
    feat_lower = feat.lower()
    if feat_lower in map_feat:
        input_feat = map_feat[feat_lower]
        output_feat = map_feat[feat_lower]
    elif "blur" in feat_lower:
        feat_base = feat_lower.replace("_blur", "")
        if feat_base in map_feat:
            input_feat = map_feat[feat_base]
            output_feat = map_feat[feat_base]
        feat_base = input_feat + "_blur"
        output_feat = output_feat + "_blur"
    else:
        input_feat = feat
        output_feat = feat
    n = 0
    for h in ["L", "R"]:
        # Set paths
        surf_file = os.path.join(
            surf_dir,
            f"{bids_id}_hemi-{h}_space-nativepro_surf-fsLR-{resol}_label-{label}.surf.gii",
        )
        prefix = f"{bids_id}_hemi-{h}_surf-fsLR-{resol}"
        if feat == "thickness":
            input_file = os.path.join(
                input_dir, f"{prefix}_label-{input_feat}.func.gii"
            )
        else:
            input_file = os.path.join(
                input_dir, f"{prefix}_label-{label}_{input_feat}.func.gii"
            )
        output_file = os.path.join(
            output_dir,
            f"{prefix}_label-{label}_feature-{output_feat}_smooth-{fwhm}mm.func.gii",
        )

        if "_blur" in feat_lower:
            inter_file = os.path.join(
                output_dir,
                f"{prefix}_label-{label}_feature-{input_feat}",
            )
            wmBoundaryDataPath = f"{input_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-white_{input_feat}.func.gii"
            wmBoundarySurfacePath = f"{surf_dir}/{bids_id}_hemi-{h}_space-nativepro_surf-fsnative_label-white.surf.gii"
            midthicknessDataPath = f"{input_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-midthickness_{input_feat}.func.gii"
            midthicknessSurfacePath = f"{surf_dir}/{bids_id}_hemi-{h}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"

            sfwm1mmDataPath = f"{input_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm1.0mm_{input_feat}.func.gii"
            sfwm1mmSurfacePath = (
                f"{surf_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm1.0mm.surf.gii"
            )
            sfwm2mmDataPath = f"{input_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm2.0mm_{input_feat}.func.gii"
            sfwm2mmSurfacePath = (
                f"{surf_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm2.0mm.surf.gii"
            )
            sfwm3mmDataPath = f"{input_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm3.0mm_{input_feat}.func.gii"
            sfwm3mmSurfacePath = (
                f"{surf_dir}/{bids_id}_hemi-{h}_surf-fsnative_label-swm3.0mm.surf.gii"
            )

            for file in [
                wmBoundaryDataPath,
                wmBoundarySurfacePath,
                midthicknessDataPath,
                midthicknessSurfacePath,
                sfwm1mmDataPath,
                sfwm1mmSurfacePath,
                sfwm2mmDataPath,
                sfwm2mmSurfacePath,
                sfwm3mmDataPath,
                sfwm3mmSurfacePath,
            ]:
                if not os.path.isfile(file):
                    show_warning(
                        f"{bids_id}: cannot map '{feat}' [label={label}, resolution={resol}] to cortex. Missing file: {file}"
                    )
                    return
            output_path = compute_blurring(
                input_dir,
                surf_dir,
                bids_id,
                h,
                inter_file,
                input_feat,
                workbench_path,
                resol,
                fwhm,
                surf_file,
                output_file,
                tmp_dir,
            )
            input_file = output_path
        # Check if file exists
        for file in [surf_file, input_file]:
            if not os.path.isfile(file):
                show_warning(
                    f"{bids_id}: cannot map '{feat}' [label={label}, resolution={resol}] to cortex. Missing file: {file}"
                )
                return

        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                input_file,
                "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT",
            ]
        )
        # Perform mapping
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-metric-smoothing",
                surf_file,
                input_file,
                str(fwhm),
                output_file,
            ]
        )
        if os.path.isfile(output_file):
            n += 1
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                output_file,
                "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT",
            ]
        )
        # add_field_to_xml(
        #     output_file,
        #     "./DataArray/MetaData",
        #     "AnatomicalStructurePrimary",
        #     "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT",
        # )

    if n == 2:
        show_note(
            f"{bids_id}: '{feat}' [label={label}, resolution={resol}] successfully mapped."
        )
    else:
        show_warning(
            f"{bids_id}: could not map '{feat}' [label={label}, resolution={resol}] to cortex."
        )


def map_hippocampus(
    bids_id,
    feat,
    resol,
    label: str,
    workbench_path: str,
    subject_hippunfold_dir: str,
    subject_micapipe_dir: str,
    subject_output_dir: str,
    folder_maps,
    folder_hip,
    fwhm,
    tmp_dir,
    subject_plugin_dir: str,
):
    if "blur" in feat:
        return
    show_info(
        f"{bids_id}: Mapping '{feat}' to hippocampus [label={label}, resolution={resol}]"
    )

    # Input & output locations
    surf_dir = os.path.join(subject_hippunfold_dir, "surf")
    input_dir = os.path.join(subject_micapipe_dir, "maps")
    is_surf = False
    if feat.startswith("plugin_"):
        surf_dir = os.path.join(subject_plugin_dir, "surf")
        input_dir = os.path.join(subject_plugin_dir, "maps")
        feat = feat[7:]
    output_dir = os.path.join(subject_output_dir, folder_maps, folder_hip)
    os.makedirs(output_dir, exist_ok=True)

    # Mappings from features names to the way they appear in the input, intermediate and output filenames
    map_input = {
        "thickness": "thickness",
        "flair": "map-flair",
        "adc": "model-DTI_map-ADC",
        "fa": "model-DTI_map-FA",
        "qt1": "map-T1map",
    }
    # feat name in filenames
    feat_lower = feat.lower()
    if feat_lower in map_input:
        input_feat = map_input[feat_lower]
        map_inter = {
            "thickness": "thickness",
            "flair": "flair",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "T1map",
        }
        inter_feat = map_inter[feat_lower]
        map_output = {
            "thickness": "thickness",
            "flair": "flair",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "T1map",
        }

        output_feat = map_output[feat_lower]
    else:
        input_feat = feat
        inter_feat = feat
        output_feat = feat

    n = 0
    for h in ["L", "R"]:

        # Set paths
        prefix = f"{bids_id}_hemi-{h}"
        surf_file = os.path.join(
            surf_dir, f"{prefix}_space-T1w_den-{resol}_label-hipp_{label}.surf.gii"
        )
        input_file = os.path.join(
            input_dir, f"{bids_id}_space-nativepro_{input_feat}.nii.gz"
        )  # Not used for thickness
        output_file = os.path.join(
            output_dir,
            f"{prefix}_den-{resol}_label-{label}_feature-{output_feat}_smooth-{fwhm}mm.func.gii",
        )

        # Note the logic here is that if a [shape|func].gii exists, use that. Otherwise map a .nii.gz file
        surf_files = glob.glob(
            os.path.join(
                surf_dir, f"{prefix}_space-T1w_den-{resol}_label-hipp_{feat}.*.gii"
            )
        )
        input_files = glob.glob(
            os.path.join(
                input_dir, f"{prefix}_space-T1w_den-{resol}_label-hipp_{feat}.*.gii"
            )
        )

        if surf_files:
            inter_file = surf_files[0]
            is_surf = True
        elif input_files:
            inter_file = input_files[0]
            is_surf = True
        else:
            inter_file = os.path.join(
                tmp_dir,
                f"{prefix}_space-T1w_desc-{inter_feat}_den-{resol}_feature-hipp_{label}.func.gii",
            )
            is_surf = False

        # Check if file exists
        check_file = inter_file if is_surf else input_file
        for file in [surf_file, check_file]:
            if not os.path.isfile(file):
                show_warning(
                    f"{bids_id}: cannot map '{feat}' [label={label}, resolution={resol}] to hippocampus. Missing file: {file}"
                )
                return

        # Perform mapping
        if not is_surf:
            subprocess.run(
                [
                    os.path.join(workbench_path, "wb_command"),
                    "-volume-to-surface-mapping",
                    input_file,
                    surf_file,
                    inter_file,
                    "-trilinear",
                ]
            )
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                inter_file,
                "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT",
            ]
        )
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-metric-smoothing",
                surf_file,
                inter_file,
                str(fwhm),
                output_file,
            ]
        )

        if os.path.isfile(output_file):
            n += 1
        subprocess.run(
            [
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                output_file,
                "CORTEX_LEFT" if h == "L" else "CORTEX_RIGHT",
            ]
        )
    if n == 2:
        show_note(
            f"{bids_id}: '{feat}' [label={label}, resolution={resol}] successfully mapped."
        )
    else:
        show_warning(
            f"{bids_id}: could not map '{feat}' [label={label}, resolution={resol}] to hippocampus."
        )


def run(
    structure,
    features,
    tmp_dir,
    WORKBENCH_PATH,
    subject_micapipe_dir,
    subject_output_dir,
    folder_maps,
    folder_ctx,
    folder_sctx,
    folder_hip,
    subject_surf_dir,
    subject_hippunfold_dir,
    script_dir,
    BIDS_ID,
    VERBOSE,
    fwhm=None,
    resolutions=[],
    labels=[],
    subject_plugin_dir=None,
):
    if structure != "subcortex":
        if fwhm is None:
            raise ValueError(
                "Error: --fwhm is required when --struct is not 'subcortex'"
            )
        if resolutions is None:
            raise ValueError(
                "Error: --resolution is required when --struct is not 'subcortex'"
            )
        if labels is None:
            raise ValueError(
                "Error: --labels is required when --struct is not 'subcortex'"
            )

    start_time = time.time()

    # Define mapping functions
    map_struct = {
        "cortex": "Cortical",
        "subcortex": "Subcortical",
        "hippocampus": "Hippocampal",
    }

    show_title(f"{map_struct[structure]} feature mapping: {BIDS_ID}")

    # Perform the mapping

    # thickness -> volume
    if structure == "subcortex":
        features = [feat.replace("thickness", "volume") for feat in features]

    # do the mapping
    for feat in features:
        if structure == "cortex":
            for res in resolutions:

                map_cortex(
                    BIDS_ID,
                    feat,
                    res,
                    labels,
                    fwhm,
                    WORKBENCH_PATH,
                    subject_micapipe_dir,
                    subject_output_dir,
                    folder_maps,
                    folder_ctx,
                    subject_plugin_dir,  # type: ignore
                    tmp_dir,
                )
        elif structure == "subcortex":
            map_subcortex(
                BIDS_ID,
                feat,
                subject_surf_dir,
                subject_micapipe_dir,
                subject_output_dir,
                folder_maps,
                folder_sctx,
                script_dir,
                subject_plugin_dir,
            )
        elif structure == "hippocampus":
            for res in resolutions:

                map_hippocampus(
                    BIDS_ID,
                    feat,
                    res,
                    labels,
                    WORKBENCH_PATH,
                    subject_hippunfold_dir,
                    subject_micapipe_dir,
                    subject_output_dir,
                    folder_maps,
                    folder_hip,
                    fwhm,
                    tmp_dir,
                    subject_plugin_dir,  # type: ignore
                )

    # Copy base T1w to output folder
    if not os.path.exists(os.path.join(subject_output_dir, "structural")):
        os.makedirs(os.path.join(subject_output_dir, "structural"))

    shutil.copyfile(
        os.path.join(
            subject_micapipe_dir,
            "anat",
            f"{BIDS_ID}_space-nativepro_T1w_brain.nii.gz",
        ),
        os.path.join(
            subject_output_dir,
            "structural",
            f"{BIDS_ID}_space-nativepro_T1w_brain.nii.gz",
        ),
    )
    shutil.copyfile(
        os.path.join(
            subject_micapipe_dir,
            "parc",
            f"{BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz",
        ),
        os.path.join(
            subject_output_dir,
            "structural",
            f"{BIDS_ID}_space-nativepro_T1w_atlas-subcortical.nii.gz",
        ),
    )
    for hemi in ["L", "R"]:
        shutil.copyfile(
            os.path.join(
                subject_micapipe_dir,
                "surf",
                f"{BIDS_ID}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
            os.path.join(
                subject_output_dir,
                "structural",
                f"{BIDS_ID}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
        )

        for surf in ["white", "pial", "midthickness"]:
            shutil.copyfile(
                os.path.join(
                    subject_micapipe_dir,
                    "surf",
                    f"{BIDS_ID}_hemi-{hemi}_space-nativepro_surf-fsnative_label-{surf}.surf.gii",
                ),
                os.path.join(
                    subject_output_dir,
                    "structural",
                    f"{BIDS_ID}_hemi-{hemi}_space-nativepro_surf-fsnative_label-{surf}.surf.gii",
                ),
            )

        for surf in ["inner", "outer", "midthickness"]:
            shutil.copyfile(
                os.path.join(
                    subject_hippunfold_dir,
                    "surf",
                    f"{BIDS_ID}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_{surf}.surf.gii",
                ),
                os.path.join(
                    subject_output_dir,
                    "structural",
                    f"{BIDS_ID}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_{surf}.surf.gii",
                ),
            )

    # Wrap up
    elapsed = round((time.time() - start_time) / 60, 2)
    print(
        f"{map_struct[structure]} feature mapping for {BIDS_ID} ended in {elapsed} minutes"
    )
