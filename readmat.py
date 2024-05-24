from ants.registration import apply_transforms
import os
from ants import image_read
import subprocess
import numpy as np

import pandas as pd
from src.functions.utilities import tempdir
from joblib import Parallel, delayed
import nibabel as nib
from joblib import Memory

from nii2dcm.run import run_nii2dcm
from niidcm import convert_nifti_to_dicom


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
    surf,
    struct,
    micapipename,
):
    if analysis == "asymmetry" and hemi == "R":
        return

    with tempdir(f"{rootzbrainfolder}", prefix="z_brains_temp.") as tmp_dir:

        metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_surf-{surf}_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
        metricsphere = "src/data/templates/fsLR-32k.L.sphere.reg.surf.gii"
        nativesphere = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii"
        outputmetric = f"{tmp_dir}/fsnative_surf_L.func.gii"
        boundingpattern = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_space-nativepro_surf-fsnative_label-"

        command1 = [
            "wb_command",
            "-metric-resample",
            metricfile,
            metricsphere,
            nativesphere,
            "BARYCENTRIC",
            outputmetric,
        ]

        command2 = [
            "wb_command",
            "-metric-to-volume-mapping",
            outputmetric,
            f"{boundingpattern}midthickness.surf.gii",
            (
                f"{rootmicafolder}/maps/{subj}_{ses}_space-nativepro_map-T1map.nii.gz"
                if feature != "flair"
                else f"{rootmicafolder}/maps/{subj}_{ses}_space-nativepro_map-flair.nii.gz"
            ),
            f"{tmp_dir}/temp.nii.gz",
            "-ribbon-constrained",
            f"{boundingpattern}white.surf.gii",
            f"{boundingpattern}pial.surf.gii",
        ]

        # Run the commands
        subprocess.run(command1)

        subprocess.run(command2)

        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp_dir}/temp.nii.gz",
        )

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            transformlist=[
                f"{rootmicafolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{rootmicafolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_surf-{surf}_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz"
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
    surf,
    struct,
    micapipename,
):
    if analysis == "asymmetry" and hemi == "R":
        return

    with tempdir(f"{rootzbrainfolder}", prefix="z_brains_temp.") as tmp_dir:
        metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
        boundingpattern = f"{rootmicafolder}/surf/{subj}_{ses}_hemi-{hemi}_space-T1w_den-0p5mm_label-hipp_"
        micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)

        command2 = [
            "wb_command",
            "-metric-to-volume-mapping",
            metricfile,
            f"{boundingpattern}midthickness.surf.gii",
            (
                f"{micapipefolder}/maps/{subj}_{ses}_space-nativepro_map-T1map.nii.gz"
                if feature != "flair"
                else f"{micapipefolder}/maps/{subj}_{ses}_space-nativepro_map-flair.nii.gz"
            ),
            f"{tmp_dir}/temp.nii.gz",
            "-ribbon-constrained",
            f"{boundingpattern}inner.surf.gii",
            f"{boundingpattern}outer.surf.gii",
        ]

        subprocess.run(command2)

        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp_dir}/temp.nii.gz",
        )

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            transformlist=[
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_hemi-{hemi}_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.nii.gz"
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
    surf,
    struct,
    micapipename,
):
    if analysis == "asymmetry" and hemi == "R":
        return
    if hemi == "R":
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
    with tempdir(f"{rootzbrainfolder}", prefix="z_brains_temp.") as tmp_dir:
        metricfile = f"{rootzbrainfolder}/norm-z/{struct}/{subj}_{ses}_feature-{feature}_analysis-{analysis}.csv"
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
            f"{tmp_dir}/temp.nii.gz",
        )

        micapipefolder = os.path.join(rootfolder, micapipename, subj, ses)

        Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

        Moving_img = image_read(
            f"{tmp_dir}/temp.nii.gz",
        )

        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
        outp = apply_transforms(
            fixed=Fixed_img,
            moving=Moving_img,
            transformlist=[
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
                f"{micapipefolder}/xfm/{subj}_{ses}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
            ],
            verbose=True,
        )
        outp.to_filename(
            f"{outdir}/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"
        )


def process(
    feature,
    hemi,
    analysis,
    rootzbrainfolder,
    rootmicafolder,
    outdir,
    subj,
    ses,
    surf,
    struct,
    micapipename,
    smooth_ctx,
    smooth_hipp,
):
    outdir = os.path.join(outdir, struct)
    if struct == "cortex":
        subdir = "micapipe"
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
            surf,
            struct,
            micapipename,
        )
    elif struct == "hippocampus":
        subdir = "hippunfold/hippunfold"
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
            surf,
            struct,
            micapipename,
        )
    elif struct == "subcortex":
        subdir = "micapipe"
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
            surf,
            struct,
            micapipename,
        )


def gluetogether(outdir, subj, ses, feature, smooth_ctx, smooth_hipp, analysis):

    cort = f"{outdir}/cortex/{subj}_{ses}_hemi-L_surf-{surf}_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
    hippo = f"{outdir}/hippocampus/{subj}_{ses}_hemi-L_den-0p5mm_label-midthickness_feature-{feature}_smooth-{smooth_hipp}_analysis-{analysis}.nii.gz"
    subcort = (
        f"{outdir}/subcortex/{subj}_{ses}_feature-{feature}_analysis-{analysis}.nii.gz"
    )

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
        cort = f"{outdir}/cortex/{subj}_{ses}_hemi-R_surf-{surf}_label-midthickness_feature-{feature}_smooth-{smooth_ctx}_analysis-{analysis}.nii.gz"
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
    else:
        outputnifti = outputnifti - np.flip(outputnifti, axis=0)

    outputnifti = nib.Nifti1Image(outputnifti, cortnifti.affine, cortnifti.header)
    if not os.path.exists(f"{outdir}/full"):
        os.makedirs(f"{outdir}/full")
    outputnifti.to_filename(
        f"{outdir}/full/{subj}_{ses}_surf-{surf}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}.nii.gz"
    )


def dicomify(
    outdir,
    subj,
    ses,
    feature,
    smooth_ctx,
    smooth_hipp,
    analysis,
):
    path = f"{outdir}/full/{subj}_{ses}_surf-{surf}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}.nii.gz"

    outpath = f"{outdir}/DICOM/{subj}_{ses}_surf-{surf}_label-midthickness_feature-{feature}_smooth-ctx-{smooth_ctx}_smooth-hipp-{smooth_hipp}_analysis-{analysis}"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # convert_nifti_to_dicom(path, outpath)
    with tempdir(f"{outdir}/DICOM", prefix="dicom_temp.") as tmp_dir:
        tempnii = nib.load(path)

        array = tempnii.get_fdata()
        # Normalize the array to the range [0, 1]
        array = (array - np.min(array)) / (np.max(array) - np.min(array))

        # Scale the array to the range [0, sitk.sitkInt16]
        array = array * np.iinfo(np.int16).max

        tempnii = nib.Nifti1Image(
            array.astype(np.int16), tempnii.affine, tempnii.header
        )

        tempnii.to_filename(f"{tmp_dir}/temp.nii.gz")
        run_nii2dcm(
            f"{tmp_dir}/temp.nii.gz", outpath, dicom_type="MR", ref_dicom_file=None
        )


# Define the parameters for the loops
# sourcery skip: avoid-builtin-shadow
os.environ["OMP_NUM_THREADS"] = "4"
features = ["ADC", "FA", "flair", "T1map", "volume", "ADC-FA-flair-T1map-volume"]
# features = ["flair"]
hemis = ["L", "R"]
analyses = ["asymmetry", "regional"]
structs = ["cortex", "hippocampus", "subcortex"]
smooth_ctx = "10mm"
smooth_hipp = "5mm"
zbrainsdir = "Test"
subj = "sub-PX001"
ses = "ses-02"

surf = "fsLR-32k"
micapipename = "micapipe"

# Define the commands
rootfolder = "E:/BIDS_MICS_Test/data/derivatives"


rootzbrainfolder = os.path.join(rootfolder, zbrainsdir, subj, ses)
outdir = os.path.join(rootzbrainfolder, "norm-z-volumetric")
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Move the definition of Fixed_img outside the loop
Fixed_img = image_read("src/data/templates/MNI152_T1_0.8mm_brain.nii.gz")

for struct in structs:
    dir = os.path.join(outdir, struct)
    if not os.path.exists(dir):
        os.makedirs(dir)


# Use joblib to parallelize the loops
# Parallel(n_jobs=4)(
#     delayed(process)(
#         feature,
#         hemi,
#         analysis,
#         rootzbrainfolder,
#         rootfolder,
#         outdir,
#         subj,
#         ses,
#         surf,
#         struct,
#         micapipename,
#         smooth_ctx,
#         smooth_hipp,
#     )
#     for feature in features
#     for hemi in hemis
#     for analysis in analyses
#     for struct in structs
# )
# Use joblib to parallelize the loops
# Parallel(n_jobs=os.cpu_count())(
#     delayed(gluetogether)(
#         outdir,
#         subj,
#         ses,
#         feature,
#         smooth_ctx,
#         smooth_hipp,
#         analysis,
#     )
#     for feature in features
#     for analysis in analyses
# )

Parallel(n_jobs=2)(
    delayed(dicomify)(
        outdir,
        subj,
        ses,
        feature,
        smooth_ctx,
        smooth_hipp,
        analysis,
    )
    for feature in features
    for analysis in analyses
)
