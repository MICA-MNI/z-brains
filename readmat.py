from ants.registration import apply_transforms
import os
from ants import image_read
import subprocess
from src.functions.utilities import tempdir

zbrainsdir = "Test"
subj = "sub-PX001"
ses = "ses-02"
feature = "FA"
smooth = "10mm"
analysis = "asymmetry"
hemi = "L"
surf = "fsLR-32k"

# Define the commands
rootfolder = "E:/BIDS_MICS_Test/data/derivatives"
rootmicafolder = f"{rootfolder}/micapipe/{subj}/{ses}"
rootzbrainfolder = f"{rootfolder}/{zbrainsdir}/{subj}/{ses}"


with tempdir("test", prefix="z_brains_temp.") as tmp_dir:

    metricfile = f"{rootzbrainfolder}/norm-z/{subj}_{ses}_hemi-{hemi}_surf-{surf}_label-midthickness_feature-{feature}_smooth-{smooth}_analysis-{analysis}.func.gii"
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
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/patient_surfs/sub-PX001_ses-02_space-nativepro_map-T1map.nii.gz",
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/Nativepro_L.nii.gz",
        "-ribbon-constrained",
        f"{boundingpattern}white.surf.gii",
        f"{boundingpattern}pial.surf.gii",
    ]

    # Run the commands
    subprocess.run(command1, check=True)
    subprocess.run(command2, check=True)

    Fixed_img = image_read(
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/templates/MNI152_T1_0.8mm_brain.nii.gz"
    )

    Moving_img = image_read(
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/Nativepro_L.nii.gz"
    )

    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
    outp = apply_transforms(
        fixed=Fixed_img,
        moving=Moving_img,
        transformlist=[
            "E:/data/derivatives/micapipe/sub-PX001/ses-02/xfm/sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz",
            "E:/data/derivatives/micapipe/sub-PX001/ses-02/xfm/sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat",
        ],
        verbose=True,
    )
    outp.to_filename(
        "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/Output_L.nii.gz"
    )
