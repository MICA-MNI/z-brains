import scipy.io
import nibabel as nib
from nibabel.affines import apply_affine
from nilearn.image import resample_img, resample_to_img
import numpy as np

mat = scipy.io.loadmat(
    "E:/data/derivatives/micapipe/sub-PX001/ses-02/xfm/sub-PX001_ses-02_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_0GenericAffine.mat"
)


def ea_antsmat2mat(afftransform, m_Center):
    # reshape the afftransform to a 3x3 matrix and transpose it
    mat = np.reshape(afftransform[:9], (3, 3)).T
    # append the rest of the afftransform as a new column to the matrix
    mat = np.column_stack((mat, afftransform[9:12]))
    m_Translation = mat[:, 3]
    # append a new row to the matrix
    mat = np.vstack((mat, [0, 0, 0, 1]))
    m_Offset = np.zeros(3)
    for i in range(3):
        m_Offset[i] = m_Translation[i] + m_Center[i]
        for j in range(3):
            m_Offset[i] -= mat[i, j] * m_Center[j]
    mat[:3, 3] = m_Offset
    mat = np.linalg.inv(mat)
    # convert RAS to LPS
    mat = mat * np.array([[1, 1, -1, -1], [1, 1, -1, -1], [-1, -1, 1, 1], [1, 1, 1, 1]])
    return mat


mat_transf = ea_antsmat2mat(mat["AffineTransform_double_3_3"], mat["fixed"])

img = nib.load(
    "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/output_okbrain.nii.gz"
)
img2 = nib.load(
    "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/mni152.nii.gz"
)
# outimg = apply_affine(mat_transf, img.get_fdata())
# clipped_img = nib.Nifti1Image(outimg, img2.affine, img.header)

outimg = resample_img(img=img, target_affine=mat_transf)
outimg.to_filename(
    "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/output_okbrain_resampled.nii.gz"
)
print("done")
