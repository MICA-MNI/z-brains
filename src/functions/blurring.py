import nibabel as nib
import numpy as np
from scipy.stats import mode
import subprocess
import os
from .utilities import show_warning


def load_gifti_data(filepath):
    data = nib.load(filepath)
    return data.darrays[0].data


def calcdist(surf1, surf2):

    euclidianDistanceSquared = (surf1 - surf2) ** 2
    euclidianDistanceSummed = np.sum(euclidianDistanceSquared, axis=1)
    return np.sqrt(euclidianDistanceSummed)


def computegrad(data, dists):
    data = np.ediff1d(data)
    data[dists == 0] = 0
    dists[dists == 0] = 1
    return np.divide(data, dists)


def compute_blurring(
    input_dir,
    surf_dir,
    bids_id,
    hemi,
    output_file,
    feat,
    workbench_path,
    resol,
    fwhm,
    surf_file,
    output_file_final,
    tmp_dir,
):

    wmBoundaryDataArr = load_gifti_data(
        f"{input_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-white_{feat}.func.gii"
    )
    wmBoundarySurfaceArr = load_gifti_data(
        f"{surf_dir}/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii"
    )

    modeofboundary = mode(wmBoundaryDataArr, keepdims=True)

    midthicknessDataArr = load_gifti_data(
        f"{input_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-midthickness_{feat}.func.gii"
    )
    midthicknessSurfaceArr = load_gifti_data(
        f"{surf_dir}/{bids_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"
    )

    surfarr = [
        [midthicknessDataArr, midthicknessSurfaceArr],
        [wmBoundaryDataArr, wmBoundarySurfaceArr],
    ]
    for dist in ["1", "2", "3"]:
        whiteMatterDataArr = load_gifti_data(
            f"{input_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-swm{dist}.0mm_{feat}.func.gii"
        )

        whiteMatterSurfaceArr = load_gifti_data(
            f"{surf_dir}/{bids_id}_hemi-{hemi}_surf-fsnative_label-swm{dist}.0mm.surf.gii"
        )
        surfarr.append(
            [
                whiteMatterDataArr,
                whiteMatterSurfaceArr,
            ]
        )

    distances = np.zeros(shape=(len(midthicknessDataArr), len(surfarr) - 1))
    dataArr = np.zeros(shape=(len(midthicknessDataArr), len(surfarr)))
    for e, ds in enumerate(surfarr):
        data, surf = ds
        dataArr[:, e] = np.divide(data, modeofboundary.mode[0])
        if e == len(surfarr) - 1:
            break
        nextdata, nextsurt = surfarr[e + 1]
        distance = calcdist(surf, nextsurt)

        distances[:, e] = distance

    blurring = np.zeros(
        shape=(len(midthicknessDataArr), len(surfarr) - 1), dtype=np.float32
    )
    for i in range(len(dataArr) - 1):
        gradient = computegrad(dataArr[i], distances[i])
        gradient = np.nan_to_num(gradient)
        if gradient[0] == 0:
            gradient = np.zeros_like(gradient)
        blurring[i] = gradient

    # for e, i in enumerate(["midthickness-0mm", "0mm-1mm", "1mm-2mm", "2mm-3mm"]):
    data_array = nib.gifti.gifti.GiftiDataArray(
        data=blurring,
        intent="NIFTI_INTENT_NORMAL",
    )

    gii = nib.gifti.GiftiImage(darrays=[data_array])
    nib.save(
        gii,
        os.path.join(
            tmp_dir,
            f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad.func.gii",
        ),
    )

    all_blurred = []
    # for i in ["midthickness-0mm", "0mm-1mm", "1mm-2mm", "2mm-3mm"]:
    subprocess.run(
        [
            os.path.join(workbench_path, "wb_command"),
            "-metric-resample",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad.func.gii",
            ),
            # "E:\data\derivatives\micapipe\sub-PX103\ses-01\surf\sub-PX103_ses-01_Normal.func.gii",
            os.path.join(
                surf_dir,
                f"{bids_id}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii",
            ),
            f"src/data/fsLR-{resol}.{hemi}.sphere.reg.surf.gii",
            "BARYCENTRIC",
            os.path.join(
                tmp_dir,
                f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad-output.func.gii",
            ),
        ]
    )
    return os.path.join(
        tmp_dir,
        f"{bids_id}-{hemi}-{feat}-{resol}-{fwhm}-surf-fsnative_grad-output.func.gii",
    )


if __name__ == "__main__":
    sub = "sub-PX103"
    ses = "ses-01"
    surface = "fsnative"
    micapipe = "micapipe"
    hemi = "L"
    input_dir = f"E:/data/derivatives/{micapipe}/{sub}/{ses}/maps/"
    surf_dir = f"E:/data/derivatives/{micapipe}/{sub}/{ses}/surf/"
    output_dir = "."
    bids_id = f"{sub}_{ses}"
    compute_blurring(
        input_dir,
        surf_dir,
        bids_id,
        hemi,
        f"{output_dir}/{bids_id}_hemi-{hemi}_blurring.func.gii",
    )
