import scipy

mat = scipy.io.loadmat(
    "E:/data/derivatives/micapipe/sub-PX001/ses-01/xfm/sub-PX001_ses-01_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_0GenericAffine.mat"
)

print(mat)
