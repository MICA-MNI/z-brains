import argparse
import nibabel as nib
import numpy as np
import os
import gc
import subprocess
import scipy
import tempfile
from src.sWM import laplace_solver, surface_generator
from src.utils import reshape_distances

def fixmatrix(path, inputmap, outputmap, basemap, BIDS_ID, temppath, wb_path, mat_path):
    mat = scipy.io.loadmat(os.path.join(path, "xfm", f"{BIDS_ID}_{mat_path}.mat"))
    affine_transform = mat["AffineTransform_double_3_3"].flatten()
    fixed = mat["fixed"].flatten()
    temp = np.identity(4)
    for i in range(3):
        temp[i, 3] = affine_transform[9 + i] + fixed[i]
        for j in range(3):
            temp[i, j] = affine_transform[i * 3 + j]
            temp[i, 3] -= temp[i, j] * fixed[j]
    flips = np.identity(4)
    flips[0, 0] = -1
    flips[1, 1] = -1
    m_matrix = np.linalg.inv(flips @ temp @ flips)
    with open(os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"), "w") as f:
        for row in m_matrix:
            f.write(" ".join(map(str, row)) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply depth-dependent blurring to a feature.")
    parser.add_argument("--participant_id", required=True)
    parser.add_argument("--session_id", required=True)
    parser.add_argument("--feature", required=True)
    parser.add_argument("--hemi", required=True)
    parser.add_argument("--output_directory", required=True)
    parser.add_argument("--workbench_path", required=True)
    parser.add_argument("--micapipe_directory", required=True)
    parser.add_argument("--freesurfer_directory", required=True)
    parser.add_argument("--tmp_dir", required=True)
    parser.add_argument("--smoothing_fwhm", type=float, required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--dist", required=True)
    parser.add_argument("--grad", required=True)
    parser.add_argument("--smooth", required=True)
    args = parser.parse_args()

    participant_id = args.participant_id
    session_id = args.session_id
    feature = args.feature
    hemi = args.hemi
    output_directory = args.output_directory
    workbench_path = args.workbench_path
    micapipe_directory = args.micapipe_directory
    freesurfer_directory = args.freesurfer_directory
    tmp_dir = args.tmp_dir
    smoothing_fwhm = args.smoothing_fwhm
    raw_data_path = args.raw
    dist_file_path = args.dist
    grad_file_path = args.grad
    smoothed_data_path = args.smooth

    input_dir = os.path.join(micapipe_directory, participant_id, session_id)
    bids_id = f"{participant_id}_{session_id}"
    freesurfer_path = os.path.join(freesurfer_directory, f"{participant_id}_{session_id}")
    temp_parc_path = os.path.join(tmp_dir, f"{participant_id}_{session_id}_{hemi}_surf-fsnative_label-temp.nii.gz")
    output_path = os.path.join(output_directory, participant_id, session_id, "structural", f"{participant_id}_{session_id}-laplace.nii.gz")
    aparc_path = os.path.join(freesurfer_path, "mri", "aparc+aseg.mgz")

    # Transform FreeSurfer parcellation to native space if needed
    if not os.path.exists(temp_parc_path):
        fixmatrix(
            path=input_dir,
            BIDS_ID=f"{participant_id}_{session_id}",
            temppath=tmp_dir,
            wb_path=workbench_path,
            inputmap=aparc_path,
            outputmap=temp_parc_path,
            basemap=os.path.join(input_dir, "anat", f"{participant_id}_{session_id}_space-nativepro_T1w_brain.nii.gz"),
            mat_path="from-fsnative_to_nativepro_T1w_0GenericAffine"
        )
        laplace_solver.solve_laplace(temp_parc_path, output_path)
        surface_generator.shift_surface(
            os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii"),
            output_path,
            os.path.join(output_directory, participant_id, session_id, "structural", f"{participant_id}_{session_id}_{hemi}_sfwm-"),
            [1.0, 2.0]
        )

    # Standardize feature name for consistent file paths
    feature_lower = feature.lower()
    if feature_lower == "t1map" or feature_lower == "qt1":
        volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz")
        feature_name = "qT1"
    elif feature_lower in ["adc", "fa"]:
        volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_model-DTI_map-{feature_lower}.nii.gz")
        feature_name = feature_lower
    else:
        volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-{feature_lower}.nii.gz")
        feature_name = feature_lower

    # Define temporary file paths
    midthickness_surf = os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii")
    midthickness_func = os.path.join(tmp_dir, f"{participant_id}_{session_id}_hemi-{hemi}_{feature}_midthickness.func.gii")
    white_surf = os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii")
    white_func = os.path.join(tmp_dir, f"{participant_id}_{session_id}_hemi-{hemi}_{feature}_white.func.gii")

    # Map volume to surfaces
    subprocess.run([
        os.path.join(workbench_path, "wb_command"),
        "-volume-to-surface-mapping",
        volumemap,
        midthickness_surf,
        midthickness_func,
        "-trilinear"
    ])
    subprocess.run([
        os.path.join(workbench_path, "wb_command"),
        "-volume-to-surface-mapping",
        volumemap,
        white_surf,
        white_func,
        "-trilinear"
    ])

    pial_data = nib.load(midthickness_func).darrays[0].data.astype(np.float32)
    pial_surface = nib.load(midthickness_surf).darrays[0].data.astype(np.float32)
    white_data = nib.load(white_func).darrays[0].data.astype(np.float32)
    white_surface = nib.load(white_surf).darrays[0].data.astype(np.float32)
    num_vertices = len(pial_data)
    surface_data = np.zeros((num_vertices, 4), dtype=np.float32)
    surface_data[:, 0] = pial_data
    surface_data[:, 1] = white_data
    del pial_data, white_data
    gc.collect()
    distances = np.zeros((num_vertices, 3), dtype=np.float32)
    gradients = np.zeros((num_vertices, 3), dtype=np.float32)
    dist_pial_white = np.sqrt(np.sum((pial_surface - white_surface)**2, axis=1)).astype(np.float32)
    distances[:, 0] = dist_pial_white
    with np.errstate(divide='ignore', invalid='ignore'):
        gradient = np.where(dist_pial_white > 0, (surface_data[:, 1] - surface_data[:, 0]) / dist_pial_white, 0).astype(np.float32)
    gradients[:, 0] = gradient
    del gradient, dist_pial_white
    gc.collect()
    prev_surface = white_surface
    prev_data = surface_data[:, 1].copy()
    for i, surf_dist in enumerate([1.0, 2.0]):
        swm_surf = os.path.join(output_directory, participant_id, session_id, "structural", f"{participant_id}_{session_id}_{hemi}_sfwm-{surf_dist}mm.surf.gii")
        swm_func = os.path.join(tmp_dir, f"{participant_id}_{session_id}_{hemi}_{feature}_sfwm-{surf_dist}mm.func.gii")
        subprocess.run([
            os.path.join(workbench_path, "wb_command"),
            "-volume-to-surface-mapping",
            volumemap,
            swm_surf,
            swm_func,
            "-trilinear"
        ])
        curr_data = nib.load(swm_func).darrays[0].data.astype(np.float32)
        curr_surface = nib.load(swm_surf).darrays[0].data.astype(np.float32)
        surface_data[:, i + 2] = curr_data
        dist = np.sqrt(np.sum((prev_surface - curr_surface)**2, axis=1)).astype(np.float32)
        distances[:, i + 1] = dist
        with np.errstate(divide='ignore', invalid='ignore'):
            grad = np.where(dist > 0, (curr_data - prev_data) / dist, 0).astype(np.float32)
        gradients[:, i + 1] = grad
        prev_surface = curr_surface
        prev_data = curr_data
        del curr_data, curr_surface, dist, grad
        gc.collect()
    del pial_surface, white_surface, prev_surface, prev_data
    gc.collect()
    distances = reshape_distances(distances)
    # Save outputs
    from nibabel import gifti
    data_array = gifti.GiftiDataArray(data=surface_data, intent="NIFTI_INTENT_NORMAL")
    gii_data = gifti.GiftiImage(darrays=[data_array])
    nib.save(gii_data, raw_data_path)
    del data_array, gii_data
    gc.collect()
    dist_array = gifti.GiftiDataArray(data=distances, intent="NIFTI_INTENT_NORMAL")
    gii_dist = gifti.GiftiImage(darrays=[dist_array])
    nib.save(gii_dist, dist_file_path)
    del distances, dist_array, gii_dist
    gc.collect()
    grad_array = gifti.GiftiDataArray(data=gradients, intent="NIFTI_INTENT_NORMAL")
    gii_grad = gifti.GiftiImage(darrays=[grad_array])
    nib.save(gii_grad, grad_file_path)
    del gradients, grad_array, gii_grad
    gc.collect()
    subprocess.run([
        os.path.join(workbench_path, "wb_command"),
        "-metric-smoothing",
        midthickness_surf,
        raw_data_path,
        str(smoothing_fwhm),
        smoothed_data_path
    ], check=False)
    subprocess.run([
        os.path.join(workbench_path, "wb_command"),
        "-set-structure",
        smoothed_data_path,
        "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT"
    ], check=False)
    del surface_data
    gc.collect() 