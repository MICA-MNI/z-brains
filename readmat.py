from ants.registration import apply_transforms
import os
from ants import image_read

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
