import os
import subprocess
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import nibabel as nib
from .niidcm import convert_nifti_to_dicom
import matplotlib.pyplot as plt
import re
from time import time
import shutil

hemis = ["L", "R"]


def float_array_to_hot(array):
    """
    Converts a floating-point array to a "hot" colormap representation.

    This function takes a floating-point array, normalizes it to the range [0, 1],
    maps the values to colors in the "hot" colormap, and returns the RGB array along with a mask.

    Args:
        array: A numpy array containing floating-point values.

    Returns:
        Tuple containing the RGB array representing the "hot" colormap and a mask array.
    """
    if np.sum(array) == 0:
        dims = array.shape + (3,)
        mask = array != 0.0
        return np.zeros(dims, dtype=np.uint8), mask
    # Normalize the array to the range [0, 1]

    array = np.clip(array, -4, 4)
    mask = array != 0.0
    # Normalize the clipped array to the range [0, 1]
    array = (array - -4) / (4 - -4)

    # Get the "hot" colormap
    cmap = plt.get_cmap("hot")

    # Convert the normalized values to colors in the "hot" colormap
    rgb_array = cmap(array)

    # Convert the colors to integers in the range [0, 255]
    rgb_array = (rgb_array[:, :, :, :3] * 255).astype(float)

    return rgb_array, mask


def float_array_to_grayscale(array):
    """
    Converts a floating-point array to a grayscale image representation.

    This function takes a floating-point array, normalizes it to the range [0, 1],
    converts the values to integers in the range [0, 255], and returns a grayscale image.

    Args:
        array: A numpy array containing floating-point values.

    Returns:
        A numpy array representing the grayscale image.
    """

    # Normalize the array to the range [0, 1]
    array = (array - np.min(array)) / (np.max(array) - np.min(array))

    # Convert the values to integers in the range [0, 255]
    int_array = (array * 255).astype(float)

    return np.stack([int_array, int_array, int_array], axis=-1)


def savevolume(
    rootmicafolder,
    subj,
    ses,
    feature,
    analysis,
    thresh,
    outdir,
    smooth_ctx,
    smooth_hipp,
    rootzbrainfolder,
    tmp_dir=None,
    vol=None,
):
    """
    Saves a volume image with specified features and analysis.

    This function processes and saves a volume image based on the provided parameters,
    including thresholding, color mapping, and saving the resulting image.

    Args:
        rootmicafolder: The root folder containing MICA data.
        subj: Subject identifier.
        ses: Session identifier.
        feature: Specific feature of the volume image.
        analysis: Type of analysis to be performed.
        thresh: Threshold value for the volume image.
        outdir: Output directory to save the processed image.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        tmp_dir: Temporary directory path (default is None).
        vol: Volume data (default is None).

    Returns:
        None
    """

    if vol is None and tmp_dir is not None:
        vol = nib.load(f"{tmp_dir}/temp.nii.gz")
        vol = vol.get_fdata()

    template = nib.load(
        f"{rootzbrainfolder}/structural/{subj}_{ses}_space-nativepro_T1w_brain.nii.gz"
    )
    template_data = template.get_fdata()

    template_data = float_array_to_grayscale(template_data)
    vol_thresh = threshold(vol, thresh)
    vol_thresh, mask_thresh = float_array_to_hot(vol_thresh)

    template_data[mask_thresh] = vol_thresh[mask_thresh]
    template = nib.Nifti1Image(template_data, template.affine, template.header)
    template.to_filename(
        f"{outdir}/full_burned/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}_threshold-{thresh}.nii.gz"
    )

    # vol, _ = float_array_to_hot(vol)
    vol = vol.astype(np.float16)
    vol_nifti = nib.Nifti1Image(vol, template.affine, template.header)
    vol_nifti.to_filename(
        f"{outdir}/full/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}.nii.gz"
    )


def process_cortex(
    feature,
    hemi,
    analysis,
    smooth,
    rootzbrainfolder,
    rootmicafolder,
    outdir,
    subj,
    ses,
    struct,
    workbench_path,
    tmp,
):
    """
    Processes cortex data for a specific feature and analysis.

    This function handles the processing of cortex data based on the provided parameters,
    including resampling, mapping to volume, and saving the processed data.

    Args:
        feature: Specific feature of the cortex data.
        hemi: Hemisphere (L or R).
        analysis: Type of analysis to be performed.
        smooth: Level of smoothing.
        rootzbrainfolder: Root folder for z-brain data.
        rootmicafolder: Root folder for MICA data.
        outdir: Output directory to save the processed data.
        subj: Subject identifier.
        ses: Session identifier.
        struct: Specific structure for processing.
        workbench_path: Path to the workbench command.
        tmp: Temporary directory path.

    Returns:
        None
    """
    if analysis == "asymmetry" and hemi == "R":
        return
    metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_surf-fsLR-32k_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
    metricsphere = f"src/data/fsLR-32k.{hemi}.sphere.reg.surf.gii"
    nativesphere = f"{rootzbrainfolder}/structural/{subj}_{ses}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii"
    boundingpattern = f"{rootzbrainfolder}/structural/{subj}_{ses}_hemi-{hemi}_space-nativepro_surf-fsnative_label-"
    if not os.path.isfile(metricfile):
        print(
            f"{feature} is not available for {subj}_{ses}_{hemi} at {smooth} smoothing in the {struct}, skipping path: {metricfile}"
        )
        return

    outputmetric = (
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_fsnative_temp.func.gii"
    )

    command_struct = [
        os.path.join(workbench_path, "wb_command"),
        "-set-structure",
        metricfile,
        "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
    ]

    # command_struct_native = [
    #     os.path.join(workbench_path, "wb_command"),
    #     "-set-structure",
    #     nativesphere,
    #     "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
    # ]
    command1 = [
        os.path.join(workbench_path, "wb_command"),
        "-metric-resample",
        metricfile,
        metricsphere,
        nativesphere,
        "BARYCENTRIC",
        outputmetric,
    ]

    command_struct_2 = [
        os.path.join(workbench_path, "wb_command"),
        "-set-structure",
        outputmetric,
        "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
    ]

    command2 = [
        os.path.join(workbench_path, "wb_command"),
        "-metric-to-volume-mapping",
        outputmetric,
        f"{boundingpattern}midthickness.surf.gii",
        f"{rootzbrainfolder}/structural/{subj}_{ses}_space-nativepro_T1w_brain.nii.gz",  # the structural image that the metric map is based off
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",  # the output file (mine gets renamed later)
        "-ribbon-constrained",
        f"{boundingpattern}white.surf.gii",  # white surf
        f"{boundingpattern}pial.surf.gii",  # pial surf
    ]

    # Run the commands
    subprocess.run(command_struct)
    # subprocess.run(command_struct_native)

    subprocess.run(command1)
    subprocess.run(command_struct_2)
    subprocess.run(command2)

    os.replace(
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
        f"{outdir}/{subj}_{ses}_hemi-{hemi}_surf-fsLR-32k_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz",
    )


def process_hippocampus(
    feature,
    hemi,
    analysis,
    smooth,
    rootzbrainfolder,
    rootmicafolder,
    outdir,
    subj,
    ses,
    struct,
    micapipename,
    rootfolder,
    workbench_path,
    tmp,
):
    """
    Processes hippocampus data for a specific feature and analysis.

    This function manages the processing of hippocampus data based on the provided parameters,
    including resampling, mapping to volume, and saving the processed data.

    Args:
        feature: Specific feature of the hippocampus data.
        hemi: Hemisphere (L or R).
        analysis: Type of analysis to be performed.
        smooth: Level of smoothing.
        rootzbrainfolder: Root folder for z-brain data.
        rootmicafolder: Root folder for MICA data.
        outdir: Output directory to save the processed data.
        subj: Subject identifier.
        ses: Session identifier.
        struct: Specific structure for processing.
        micapipename: Name of the micapipe.
        rootfolder: Root folder for data.
        workbench_path: Path to the workbench command.
        tmp: Temporary directory path.

    Returns:
        None
    """
    if analysis == "asymmetry" and hemi == "R":
        return
    metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
    boundingpattern = f"{rootzbrainfolder}/structural/{subj}_{ses}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_"
    micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)

    if not os.path.isfile(metricfile):
        print(
            f"{feature} is not available for {subj}_{ses}_{hemi} at {smooth} smoothing in the {struct}, skipping"
        )
        return

    command_struct = [
        os.path.join(workbench_path, "wb_command"),
        "-set-structure",
        metricfile,
        "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
    ]

    command2 = [
        os.path.join(workbench_path, "wb_command"),
        "-metric-to-volume-mapping",
        metricfile,
        f"{boundingpattern}midthickness.surf.gii",
        f"{rootzbrainfolder}/structural/{subj}_{ses}_space-nativepro_T1w_brain.nii.gz",
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
        "-ribbon-constrained",
        f"{boundingpattern}inner.surf.gii",
        f"{boundingpattern}outer.surf.gii",
    ]

    subprocess.run(command_struct)

    subprocess.run(command2)

    os.replace(
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
        f"{outdir}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz",
    )


def process_subcortex(
    feature,
    hemi,
    analysis,
    rootzbrainfolder,
    rootmicafolder,
    outdir,
    subj,
    ses,
    struct,
    tmp,
):
    """
    Processes subcortex data for a specific feature and analysis.

    This function handles the processing of subcortex data based on the provided parameters,
    including matching data to subcortical structures, creating an output atlas, and saving the processed data.

    Args:
        feature: Specific feature of the subcortex data.
        hemi: Hemisphere (L or R).
        analysis: Type of analysis to be performed.
        rootzbrainfolder: Root folder for z-brain data.
        rootmicafolder: Root folder for MICA data.
        outdir: Output directory to save the processed data.
        subj: Subject identifier.
        ses: Session identifier.
        struct: Specific structure for processing.
        tmp: Temporary directory path.

    Returns:
        None
    """
    if analysis == "asymmetry" and hemi == "R":
        return
    if hemi == "R":
        return

    metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_feature-{feature}_analysis-{analysis}.csv"
    if not os.path.isfile(metricfile):
        print(f"{feature} is not available for {subj}_{ses} in the {struct}, skipping")
        return
    STRUCTURES = {
        "Laccumb": 26,
        "Lamyg": 18,
        "Lcaud": 11,
        "Lhippo": 17,
        "Lpal": 13,
        "Lput": 12,
        "Lthal": 10,
        "Raccumb": 58,
        "Ramyg": 54,
        "Rcaud": 50,
        "Rhippo": 53,
        "Rpal": 52,
        "Rput": 51,
        "Rthal": 49,
    }

    atlas = nib.load(
        f"{rootzbrainfolder}/structural/{subj}_{ses}_space-nativepro_T1w_atlas-subcortical.nii.gz"
    )
    atlasdata = atlas.get_fdata()
    metricdata = pd.read_csv(metricfile).to_dict()
    matcheddata = {
        STRUCTURES[x]: metricdata[x][0] for x in metricdata if x in STRUCTURES
    }

    outputatlas = np.zeros_like(atlasdata, dtype=np.float64)
    for key, value in matcheddata.items():
        outputatlas[atlasdata == key] = value

    output_img = nib.Nifti1Image(
        outputatlas, atlas.affine, atlas.header, dtype=np.float64
    )
    output_img.to_filename(
        f"{tmp}/{feature}_{analysis}_{struct}_temp.nii.gz",
    )

    os.replace(
        f"{tmp}/{feature}_{analysis}_{struct}_temp.nii.gz",
        f"{outdir}/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz",
    )


def process(
    feature,
    hemi,
    analysis,
    rootzbrainfolder,
    rootfolder,
    outdir,
    subj,
    ses,
    struct,
    micapipename,
    hippunfoldname,
    smooth_ctx,
    smooth_hipp,
    workbench_path,
    tmp,
):
    """
    Processes different brain structures based on the provided parameters.

    This function orchestrates the processing of various brain structures, such as cortex, hippocampus, and subcortex,
    by calling specific processing functions based on the structure type and input parameters.

    Args:
        feature: Specific feature of the brain structure data.
        hemi: Hemisphere (L or R).
        analysis: Type of analysis to be performed.
        rootzbrainfolder: Root folder for z-brain data.
        rootfolder: Root folder for data.
        outdir: Output directory to save the processed data.
        subj: Subject identifier.
        ses: Session identifier.
        struct: Specific brain structure for processing.
        micapipename: Name of the micapipe.
        hippunfoldname: Name of the hippocampus unfolding.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        workbench_path: Path to the workbench command.
        tmp: Temporary directory path.

    Returns:
        None
    """
    outdir = os.path.join(outdir, struct)
    print(f"Processing structure: {struct}")
    if struct == "cortex":
        subdir = micapipename
        smooth = smooth_ctx
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
        print(f"Processing cortex with smoothing level: {smooth}")
        process_cortex(
            feature,
            hemi,
            analysis,
            smooth,
            rootzbrainfolder,
            rootsubdir,
            outdir,
            subj,
            ses,
            struct,
            workbench_path,
            tmp,
        )
    elif struct == "hippocampus":
        subdir = f"{hippunfoldname}/hippunfold"
        smooth = smooth_hipp
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
        print(f"Processing hippocampus with smoothing level: {smooth}")
        process_hippocampus(
            feature,
            hemi,
            analysis,
            smooth,
            rootzbrainfolder,
            rootsubdir,
            outdir,
            subj,
            ses,
            struct,
            micapipename,
            rootfolder,
            workbench_path,
            tmp,
        )
    elif struct == "subcortex":
        subdir = micapipename
        smooth = None
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
        print("Processing subcortex without smoothing")
        process_subcortex(
            feature,
            hemi,
            analysis,
            rootzbrainfolder,
            rootsubdir,
            outdir,
            subj,
            ses,
            struct,
            tmp,
        )
    print("Processing completed.")


def gluetogether(
    outdir,
    subj,
    ses,
    feature,
    smooth_ctx,
    smooth_hipp,
    analysis,
    rootfolder,
    micapipename,
    tmp,
    rootzbrainfolder,
    thresh,
):
    """
    Combines different brain structure data into a single volume image.

    This function merges cortex, hippocampus, and subcortex data into a unified volume image,
    handling asymmetry cases and saving the combined data to the specified output directory.

    Args:
        outdir: Output directory for saving the combined volume image.
        subj: Subject identifier.
        ses: Session identifier.
        feature: Specific feature of the brain structure data.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        analysis: Type of analysis to be performed.
        rootfolder: Root folder for data.
        micapipename: Name of the micapipe.
        tmp: Temporary directory path.
        thresh: Threshold value for the volume image.

    Returns:
        None
    """
    print(f"Combining brain structure data for subject {subj}, session {ses}.")
    cort = f"{outdir}/cortex/{subj}_{ses}_hemi-L_surf-fsLR-32k_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
    hippo = f"{outdir}/hippocampus/{subj}_{ses}_hemi-L_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth_hipp}_analysis-{analysis}.nii.gz"
    subcort = (
        f"{outdir}/subcortex/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"
    )
    if "blur" in feature:
        if not os.path.isfile(cort):
            print(f"{feature} is not available for {subj}_{ses}, skipping")
            return
        else:
            print(f"Loading file: {cort}")
    else:
        for each in [cort, hippo, subcort]:
            if not os.path.isfile(each):
                print(f"{feature} is not available for {subj}_{ses}, skipping")
                return
            else:
                print(f"Loading file: {each}")

    cortnifti = nib.load(cort)
    if not "blur" in feature:
        hipponifti = nib.load(hippo)
        subcortnifti = nib.load(subcort)
    print("Data loaded, starting to combine structures.")

    cortdata = cortnifti.get_fdata()
    if not "blur" in feature:
        hippodata = hipponifti.get_fdata()
        subcortdata = subcortnifti.get_fdata()
    outputnifti = np.zeros_like(cortdata)
    outputnifti[cortdata != 0] = cortdata[cortdata != 0]
    if not "blur" in feature:
        outputnifti[subcortdata != 0] = subcortdata[subcortdata != 0]
        outputnifti[hippodata != 0] = hippodata[hippodata != 0]

    print("Left hemisphere data combined.")

    if analysis != "asymmetry":
        cort = f"{outdir}/cortex/{subj}_{ses}_hemi-R_surf-fsLR-32k_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
        hippo = f"{outdir}/hippocampus/{subj}_{ses}_hemi-R_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth_hipp}_analysis-{analysis}.nii.gz"

        cortnifti = nib.load(cort)
        if not "blur" in feature:
            hipponifti = nib.load(hippo)

        cortdata = cortnifti.get_fdata()
        if not "blur" in feature:
            hippodata = hipponifti.get_fdata()
        outputnifti[cortdata != 0] = cortdata[cortdata != 0]
        if not "blur" in feature:
            outputnifti[hippodata != 0] = hippodata[hippodata != 0]
        print("Right hemisphere data combined.")

    micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)
    print("Saving combined volume image.")
    savevolume(
        micapipefolder,
        subj,
        ses,
        feature,
        analysis,
        thresh,
        outdir,
        smooth_ctx,
        smooth_hipp,
        rootzbrainfolder,
        vol=outputnifti,
    )
    print("Combined volume image saved.")


def threshold(array, threshold):
    """
    Apply thresholding to an array based on a specified threshold value.

    This function sets array elements to 0 if they fall within the specified threshold range,
    effectively thresholding the array values.

    Args:
        array: Input numpy array to be thresholded.
        threshold: Threshold value for the array.

    Returns:
        Numpy array with thresholding applied.
    """
    array = np.where((array < threshold) & (array > 0), 0, array)
    array = np.where((array > -threshold) & (array < 0), 0, array)
    return array


def dicomify(
    outdir,
    subj,
    ses,
    feature,
    smooth_ctx,
    smooth_hipp,
    analysis,
    tmp,
    thresh,
    px_demo=None,
):
    """
    Convert a NIfTI image to DICOM format for a specific subject and session.

    This function converts a NIfTI image to DICOM format based on the provided parameters,
    creating DICOM files in the specified output directory.

    Args:
        outdir: Output directory for saving the DICOM files.
        subj: Subject identifier.
        ses: Session identifier.
        feature: Specific feature of the image.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        analysis: Type of analysis performed on the image.
        tmp: Temporary directory path.
        px_demo: Path to the participant demographics file (default is None).

    Returns:
        None
    """
    path = f"{outdir}/full_burned/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}_threshold-{thresh}.nii.gz"

    if not os.path.isfile(path):
        print(
            f"File not found at {path}. {feature} and the {analysis} analysis is not available for {subj}_{ses}, skipping"
        )
        return

    outpath = f"{outdir}/DICOM/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}_threshold-{thresh}"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    tempnii = nib.load(path)

    array = tempnii.get_fdata()

    convert_nifti_to_dicom(
        array.astype(np.int16),
        tempnii.header,
        tempnii.affine,
        outpath,
        feature,
        smooth_ctx,
        smooth_hipp,
        analysis,
        px_demo,
    )


def surface_to_volume(
    rootfolder,
    features,
    analyses,
    structs,
    smooth_ctx,
    smooth_hipp,
    zbrainsdir,
    subj,
    ses,
    px_demo,
    micapipename,
    hippunfoldname,
    tmp,
    n_jobs,
    n_jobs_wb,
    workbench_path,
    thresh=3,
    dicoms=None,
):
    """
    Process surface data to generate volumetric images for specified features and analyses.

    This function orchestrates the conversion of surface data to volumetric images,
    handling various features, structures, and analyses based on the provided parameters.

    Args:
        rootfolder: Root folder for data processing.
        features: List of specific features to process.
        analyses: List of analysis types to perform.
        structs: List of brain structures to process.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        zbrainsdir: Directory containing z-brain data.
        subj: Subject identifier.
        ses: Session identifier.
        px_demo: Participant demographics data.
        micapipename: Name of the micapipe.
        hippunfoldname: Name of the hippocampus unfolding.
        tmp: Temporary directory path.
        n_jobs: Number of parallel jobs to run.
        n_jobs_wb: Number of parallel jobs for workbench.
        workbench_path: Path to the workbench command.

    Returns:
        None
    """
    rootfolder = os.path.join(rootfolder, "derivatives")
    zbrainsdir = zbrainsdir[0]
    smooth_ctx = f"{str(smooth_ctx)}mm"
    smooth_hipp = f"{str(smooth_hipp)}mm"

    os.environ["OMP_NUM_THREADS"] = str(n_jobs_wb)
    if isinstance(px_demo, pd.DataFrame):
        px_demo = px_demo[px_demo["participant_id"] == subj]
        px_demo = px_demo[px_demo["session_id"] == ses]
        px_demo = px_demo.reset_index(drop=True)
    rootzbrainfolder = os.path.join(rootfolder, zbrainsdir, subj, ses)
    outdir = os.path.join(rootzbrainfolder, "norm-z-volumetric")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for struct in structs:
        structdir = os.path.join(outdir, struct)
        if not os.path.exists(structdir):
            os.makedirs(structdir)

    features = [f.replace("qT1", "T1map") for f in features]
    available_features = os.listdir(os.path.join(rootzbrainfolder, "maps", "cortex"))
    feats = []
    for feat in available_features:

        # Use re.search() to find the first match of the pattern in the string
        match = re.search(r"feature-([a-zA-Z0-9]+(?:_blur)?)(?:_smooth)?", feat)

        # Extract the matched group
        feature = match[1] if match else None
        feats.append(feature)

    features = list(set(features).intersection(feats))
    if "thickness" in features:
        features[features.index("thickness")] = "volume"
    features = sorted(features, key=str.lower)
    features.append("-".join([x for x in features if "blur" not in x]))
    print("feats: ", features)
    # shutil.copyfile(
    #     os.path.join(
    #         rootfolder,
    #         micapipename,
    #         subj,
    #         ses,
    #         "anat",
    #         f"{subj}_{ses}_space-nativepro_T1w_brain.nii.gz",
    #     ),
    #     os.path.join(outdir, "base_T1w.nii.gz"),
    # )

    Parallel(n_jobs=n_jobs)(
        delayed(process)(
            feature,
            hemi,
            analysis,
            rootzbrainfolder,
            rootfolder,
            outdir,
            subj,
            ses,
            struct,
            micapipename,
            hippunfoldname,
            smooth_ctx,
            smooth_hipp,
            workbench_path,
            tmp,
        )
        for feature in features
        for hemi in hemis
        for analysis in analyses
        for struct in structs
    )

    if not os.path.exists(f"{outdir}/full"):
        os.makedirs(f"{outdir}/full")
    if not os.path.exists(f"{outdir}/full_burned"):
        os.makedirs(f"{outdir}/full_burned")
    Parallel(n_jobs=n_jobs)(
        delayed(gluetogether)(
            outdir,
            subj,
            ses,
            feature,
            smooth_ctx,
            smooth_hipp,
            analysis,
            rootfolder,
            micapipename,
            tmp,
            rootzbrainfolder,
            thresh=thresh,
        )
        for feature in features
        for analysis in analyses
    )
    if dicoms == 1:
        print("Converting to DICOM")
        timepre = time()
        Parallel(n_jobs=n_jobs)(
            delayed(dicomify)(
                outdir,
                subj,
                ses,
                feature,
                smooth_ctx,
                smooth_hipp,
                analysis,
                tmp,
                thresh,
                px_demo=px_demo,
            )
            for feature in features
            for analysis in analyses
        )
        timepost = time() - timepre
        print(f"Time taken to convert to DICOM: {timepost}")


if __name__ == "__main__":
    features = ["ADC", "FA", "flair", "T1map", "volume"]
    # features = ["flair"]
    hemis = ["L", "R"]
    analyses = ["asymmetry", "regional"]
    structs = ["cortex", "hippocampus", "subcortex"]
    smooth_ctx = "10"
    smooth_hipp = "5"
    zbrainsdir = ["Test"]
    subj = "sub-PX001"
    ses = "ses-02"

    px_demo = pd.read_csv("E:/BIDS_MICS_Test/PX_participants.csv")
    micapipename = "micapipe"
    hippunfoldname = "hippunfold"
    n_jobs = 1
    n_jobs_wb = 1
    # Define the commands
    rootfolder = "E:/BIDS_MICS_Test/data"
    tmp = "E:/tmp"
    os.makedirs(tmp, exist_ok=True)
    surface_to_volume(
        rootfolder,
        features,
        analyses,
        structs,
        smooth_ctx,
        smooth_hipp,
        zbrainsdir,
        subj,
        ses,
        px_demo,
        micapipename,
        hippunfoldname,
        tmp,
        n_jobs,
        n_jobs_wb,
        "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64",
    )
