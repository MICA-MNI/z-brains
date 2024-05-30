from ants.registration import apply_transforms
import os
from ants import image_read, image_clone
import subprocess
import numpy as np
import pandas as pd
from .utilities import tempdir
from joblib import Parallel, delayed
import nibabel as nib
from .niidcm import convert_nifti_to_dicom
import matplotlib.pyplot as plt
import re

hemis = ["L", "R"]


def float_array_to_hot(array):
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
    tmp_dir=None,
    vol=None,
):
    if vol is None and tmp_dir is not None:
        vol = nib.load(f"{tmp_dir}/temp.nii.gz")
        vol = vol.get_fdata()

    template = nib.load(
        f"{rootmicafolder}/anat/{subj}_{ses}_space-nativepro_T1w_brain.nii.gz"
    )
    template_data = template.get_fdata()
    template_data = float_array_to_grayscale(template_data)
    print(np.min(vol), np.max(vol))
    vol = threshold(vol, thresh)
    vol, mask = float_array_to_hot(vol)
    template_data[mask] = vol[mask]
    template = nib.Nifti1Image(template_data, template.affine, template.header)
    template.to_filename(
        f"{outdir}/full/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}.nii.gz"
    )


# @memory.cache
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
    micapipename,
    workbench_path,
    tmp,
    thresh=2,
    mni=False,
):

    if analysis == "asymmetry" and hemi == "R":
        return
    metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_surf-fslr-32k_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
    metricsphere = "src/data/templates/fsLR-32k.L.sphere.reg.surf.gii"
    nativesphere = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii"
    boundingpattern = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_space-nativepro_surf-fsnative_label-"
    if not os.path.isfile(metricfile):
        print(
            f"{feature} is not available for {subj}_{ses}_{hemi} at {smooth} smoothing in the {struct}, skipping"
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

    command1 = [
        os.path.join(workbench_path, "wb_command"),
        "-metric-resample",
        metricfile,
        metricsphere,
        nativesphere,
        "BARYCENTRIC",
        outputmetric,
    ]

    command2 = [
        os.path.join(workbench_path, "wb_command"),
        "-metric-to-volume-mapping",
        outputmetric,
        f"{boundingpattern}midthickness.surf.gii",
        (
            f"{rootmicafolder}/maps/{subj}_{ses}_space-nativepro_map-T1map.nii.gz"
            if feature != "flair"
            else f"{rootmicafolder}/maps/{subj}_{ses}_space-nativepro_map-flair.nii.gz"
        ),
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
        "-ribbon-constrained",
        f"{boundingpattern}white.surf.gii",
        f"{boundingpattern}pial.surf.gii",
    ]

    # Run the commands
    subprocess.run(command_struct)

    subprocess.run(command1)

    subprocess.run(command2)
    if mni:
        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp}/temp.nii.gz",
        )

        array = Moving_img.numpy()
        array = threshold(array, thresh)
        Moving_img = image_clone(Moving_img, array)

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            interpolator="linear",
            transformlist=[
                f"{rootmicafolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{rootmicafolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_surf-fslr-32k_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz"
        )
    else:
        os.replace(
            f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_surf-fslr-32k_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz",
        )


# @memory.cache
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
    thresh=2,
    mni=False,
):
    if analysis == "asymmetry" and hemi == "R":
        return
    metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
    boundingpattern = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_"
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
        (
            f"{micapipefolder}/maps/{subj}_{ses}_space-nativepro_map-T1map.nii.gz"
            if feature != "flair"
            else f"{micapipefolder}/maps/{subj}_{ses}_space-nativepro_map-flair.nii.gz"
        ),
        f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
        "-ribbon-constrained",
        f"{boundingpattern}inner.surf.gii",
        f"{boundingpattern}outer.surf.gii",
    ]

    subprocess.run(command_struct)

    subprocess.run(command2)

    if mni:
        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp_dir}/temp.nii.gz",
        )

        array = Moving_img.numpy()
        array = threshold(array, thresh)
        Moving_img = image_clone(Moving_img, array)

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            interpolator="linear",
            transformlist=[
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz"
        )
    else:
        os.replace(
            f"{tmp}/{feature}_{analysis}_{struct}_{smooth}_{hemi}_temp.nii.gz",
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz",
        )


# @memory.cache
def process_subcortex(
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
    tmp,
    thresh=2,
    mni=False,
):

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
        f"{rootmicafolder}/parc/{subj}_{ses}_space-nativepro_T1w_atlas-subcortical.nii.gz"
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

    micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)
    if mni:
        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp_dir}/temp.nii.gz",
        )
        array = Moving_img.numpy()
        array = threshold(array, thresh)
        Moving_img = image_clone(Moving_img, array)

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            interpolator="linear",
            transformlist=[
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"
        )
    else:
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
    thresh=2,
):
    outdir = os.path.join(outdir, struct)
    if struct == "cortex":
        subdir = micapipename
        smooth = smooth_ctx
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
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
            micapipename,
            workbench_path,
            tmp,
            thresh=thresh,
        )
    elif struct == "hippocampus":
        subdir = f"{hippunfoldname}/hippunfold"
        smooth = smooth_hipp
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
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
            thresh=thresh,
        )
    elif struct == "subcortex":
        subdir = micapipename
        smooth = None
        rootsubdir = os.path.join(rootfolder, subdir, subj, ses)
        process_subcortex(
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
            tmp,
            thresh=thresh,
        )


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
    thresh,
):

    cort = f"{outdir}/cortex/{subj}_{ses}_hemi-L_surf-fslr-32k_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
    hippo = f"{outdir}/hippocampus/{subj}_{ses}_hemi-L_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth_hipp}_analysis-{analysis}.nii.gz"
    subcort = (
        f"{outdir}/subcortex/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"
    )
    for each in [cort, hippo, subcort]:
        if not os.path.isfile(each):
            print(f"{feature} is not available for {subj}_{ses}, skipping")
            return
    cortnifti = nib.load(cort)
    hipponifti = nib.load(hippo)
    subcortnifti = nib.load(subcort)

    cortdata = cortnifti.get_fdata()
    hippodata = hipponifti.get_fdata()
    subcortdata = subcortnifti.get_fdata()
    print(np.average(cortdata), np.average(hippodata), np.average(subcortdata))
    outputnifti = np.zeros_like(cortdata)
    outputnifti[cortdata != 0] = cortdata[cortdata != 0]
    outputnifti[subcortdata != 0] = subcortdata[subcortdata != 0]
    outputnifti[hippodata != 0] = hippodata[hippodata != 0]

    if analysis != "asymmetry":
        cort = f"{outdir}/cortex/{subj}_{ses}_hemi-R_surf-fslr-32k_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
        hippo = f"{outdir}/hippocampus/{subj}_{ses}_hemi-R_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth_hipp}_analysis-{analysis}.nii.gz"
        # subcort = f"{outdir}/subcortex/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"

        cortnifti = nib.load(cort)
        hipponifti = nib.load(hippo)
        # subcortnifti = nib.load(subcort)

        cortdata = cortnifti.get_fdata()
        hippodata = hipponifti.get_fdata()
        # subcortdata = subcortnifti.get_fdata()
        # print(np.average(cortdata), np.average(hippodata), np.average(subcortdata))
        outputnifti[cortdata != 0] = cortdata[cortdata != 0]
        # outputnifti[subcortdata != 0] = subcortdata[subcortdata != 0]
        outputnifti[hippodata != 0] = hippodata[hippodata != 0]
        # else:
        #     outputnifti = outputnifti - np.flip(outputnifti, axis=0)

    micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)
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
        vol=outputnifti,
    )


def threshold(array, threshold):
    array = np.where((array < threshold) & (array > 0), 0, array)
    array = np.where((array > -threshold) & (array < 0), 0, array)
    return array


def dicomify(
    outdir, subj, ses, feature, smooth_ctx, smooth_hipp, analysis, tmp, px_demo=None
):
    path = f"{outdir}/full/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}.nii.gz"

    if not os.path.isfile(path):
        print(
            f"{feature} and the {analysis} analysis is not available for {subj}_{ses}, skipping"
        )
        return

    outpath = f"{outdir}/DICOM/{subj}_{ses}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    tempnii = nib.load(path)

    array = tempnii.get_fdata()
    tempnii = nib.Nifti1Image(array.astype(np.int16), tempnii.affine, tempnii.header)
    convert_nifti_to_dicom(
        tempnii, outpath, feature, smooth_ctx, smooth_hipp, analysis, px_demo
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
):

    rootfolder = os.path.join(rootfolder, "derivatives")
    zbrainsdir = zbrainsdir[0]
    smooth_ctx = str(smooth_ctx) + "mm"
    smooth_hipp = str(smooth_hipp) + "mm"

    os.environ["OMP_NUM_THREADS"] = str(n_jobs_wb)
    print(px_demo)
    # px_demo = pd.read_csv(px_demo)
    px_demo = px_demo[px_demo["participant_id"] == subj]
    px_demo = px_demo[px_demo["session_id"] == ses]

    rootzbrainfolder = os.path.join(rootfolder, zbrainsdir, subj, ses)
    outdir = os.path.join(rootzbrainfolder, "norm-z-volumetric")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Move the definition of Fixed_img outside the loop
    # Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

    for struct in structs:
        dir = os.path.join(outdir, struct)
        if not os.path.exists(dir):
            os.makedirs(dir)

    print(features)
    if "qT1" in features:
        features[features.index("qT1")] = "T1map"
    available_features = os.listdir(os.path.join(rootzbrainfolder, "maps", "cortex"))
    feats = []
    for feat in available_features:
        # Use re.search() to find the first match of the pattern in the string
        match = re.search(r"feature-([a-zA-Z0-9]+)", feat)

        # Extract the matched group
        feature = match.group(1) if match else None
        feats.append(feature)

    features = list(set(features).intersection(feats))
    if "thickness" in features:
        features[features.index("thickness")] = "volume"
    features.sort()
    features.append("-".join(features))
    print(features)
    # Use joblib to parallelize the loops
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
    # Use joblib to parallelize the loops
    if not os.path.exists(f"{outdir}/full"):
        os.makedirs(f"{outdir}/full")
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
            thresh=2,
        )
        for feature in features
        for analysis in analyses
    )

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
            px_demo=px_demo,
        )
        for feature in features
        for analysis in analyses
    )


if __name__ == "__main__":
    # features = ["ADC", "FA", "flair", "T1map", "volume", "ADC-FA-flair-T1map-volume"]
    # features = ["flair"]
    hemis = ["L", "R"]
    analyses = ["asymmetry", "regional"]
    structs = ["cortex", "hippocampus", "subcortex"]
    smooth_ctx = "10mm"
    smooth_hipp = "5mm"
    zbrainsdir = "Test"
    subj = "sub-PX001"
    ses = "ses-02"
    px_demo = "E:/BIDS_MICS_Test/PX_participants.csv"
    micapipename = "micapipe"
    hippunfoldname = "hippunfold"
    n_jobs = 4
    n_jobs_wb = 1
    # Define the commands
    rootfolder = "E:/BIDS_MICS_Test/data/derivatives"
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
        n_jobs,
        n_jobs_wb,
    )
