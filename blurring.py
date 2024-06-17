import nibabel as nib
import numpy as np


def load_gifti_data(filepath):
    data = nib.load(filepath)
    return data.darrays[0].data


sub = "sub-PX103"
ses = "ses-01"
surface = "fsLR-32k"
midthicknessDataArr = load_gifti_data(
    f"E:\\data\\derivatives\\micapipe\\{sub}\\{ses}\\maps\\{sub}_{ses}_hemi-L_surf-{surface}_label-midthickness_T1map.func.gii"
)
whiteMatterDataArr = load_gifti_data(
    f"E:\\data\\derivatives\\micapipe\\{sub}\\{ses}\\maps\\{sub}_{ses}_hemi-L_surf-{surface}_label-white_T1map.func.gii"
)
midthicknessSurfaceArr = load_gifti_data(
    f"E:\\data\\derivatives\\micapipe\\{sub}\\{ses}\\surf\\{sub}_{ses}_hemi-L_space-nativepro_surf-{surface}_label-midthickness.surf.gii"
)
whiteMatterSurfaceArr = load_gifti_data(
    f"E:\\data\\derivatives\\micapipe\\{sub}\\{ses}\\surf\\{sub}_{ses}_hemi-L_space-nativepro_surf-{surface}_label-white.surf.gii"
)

euclidianDistanceSquared = (whiteMatterSurfaceArr - midthicknessSurfaceArr) ** 2
euclidianDistanceSummed = np.sum(euclidianDistanceSquared, axis=1)
euclidean_distance = np.sqrt(euclidianDistanceSummed)
# euclidean_distance = np.where(euclidean_distance == 0, 1, euclidean_distance)


wmmthckDifference = midthicknessDataArr - whiteMatterDataArr
blurring = np.multiply(wmmthckDifference, euclidean_distance)
nans = np.sum(np.isnan(blurring))

mask = euclidean_distance == 0
values = wmmthckDifference[mask]

data_array = nib.gifti.gifti.GiftiDataArray(
    data=blurring,
    intent="NIFTI_INTENT_NORMAL",
)

gii = nib.gifti.GiftiImage(darrays=[data_array])
nib.save(gii, f"{sub}_{ses}_multiply.func.gii")

print(values)
print("here")
