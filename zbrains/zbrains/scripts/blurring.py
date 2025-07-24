#!/usr/bin/env python3
"""
Standalone blurring pipeline script replicating blurring.smk logic.
- No imports from src/ (self-contained)
- Uses nibabel, numpy, scipy, subprocess, argparse
- For each subject/session/feature:
    1. Laplace solver on aparc+aseg.nii.gz
    2. Shift white surface to 1mm/2mm using Laplace field
    3. Map feature volume to midthickness, white, sfwm surfaces
    4. Calculate distances/gradients
    5. Save raw, dist, grad, smoothed outputs
    6. Apply metric smoothing and set-structure
"""
import os
import argparse
import subprocess
import numpy as np
import nibabel as nib
from scipy.ndimage import binary_dilation, convolve
import skfmm
import tempfile
import gc

def fast_convolve(data, kernel):
    nan_mask = np.isnan(data)
    if np.any(nan_mask):
        data_filled = np.copy(data)
        data_filled[nan_mask] = 0.0
        result = convolve(data_filled, kernel, mode='constant', cval=0.0)
        weights = convolve((~nan_mask).astype(np.float32), kernel, mode='constant', cval=0.0)
        valid_weights = weights > 1e-10
        result[valid_weights] /= weights[valid_weights]
        result[~valid_weights] = np.nan
        return result
    else:
        return convolve(data, kernel, mode='constant', cval=0.0)

def laplace_solver(in_seg, out_laplace, verbose=True):
    # Parameters
    convergence_threshold = 1e-4
    max_iters = 500
    kernelSize = 3
    initial_lr = 0.5
    lr_decay = 0.999
    fg_labels = [41, 2]
    src_labels = np.concatenate((np.arange(1000, 2999), [0]))
    alpha = 0.1
    # Load segmentation
    lbl_nib = nib.load(in_seg)
    lbl = np.asarray(lbl_nib.dataobj, dtype=np.int16)
    fg = np.isin(lbl, fg_labels)
    fg = binary_dilation(fg, iterations=1)
    source = np.isin(lbl, src_labels)
    source = source & ~fg
    sink = ~(fg | source)
    phi = np.ones_like(lbl, dtype=np.float32)
    phi[source] = 0.0
    mask = ~(fg | source)
    phi_masked = np.ma.MaskedArray(phi, mask)
    forward = skfmm.travel_time(phi_masked, np.ones_like(lbl, dtype=np.float32))
    init_coords = forward.data
    finite_mask = np.isfinite(init_coords)
    if np.any(finite_mask):
        min_val = np.min(init_coords[finite_mask])
        max_val = np.max(init_coords[finite_mask])
        if max_val > min_val:
            init_coords = (init_coords - min_val) / (max_val - min_val)
        else:
            init_coords.fill(0)
    init_coords[fg] = 0.0
    min_bbox, max_bbox = get_bounding_box(fg, padding=kernelSize)
    cropped_init = crop(init_coords, min_bbox, max_bbox)
    cropped_fg = crop(fg, min_bbox, max_bbox)
    cropped_source = crop(source, min_bbox, max_bbox)
    cropped_sink = crop(sink, min_bbox, max_bbox)
    coords = laplace(
        cropped_init, cropped_fg, cropped_source, cropped_sink,
        kernelSize=kernelSize, convergence_threshold=convergence_threshold,
        max_iters=max_iters, initial_lr=initial_lr, lr_decay=lr_decay, verbose=verbose
    )
    full_coords = np.zeros_like(init_coords, dtype=np.float32)
    full_coords[min_bbox[0]:max_bbox[0], min_bbox[1]:max_bbox[1], min_bbox[2]:max_bbox[2]] = coords
    full_coords[source] = 0.0
    full_coords[sink] = 1.0
    if alpha > 0:
        full_coords = full_coords * (1 - alpha) + init_coords * alpha
    coords_nib = nib.Nifti1Image(full_coords, lbl_nib.affine, lbl_nib.header)
    coords_nib.set_data_dtype(np.float32)
    nib.save(coords_nib, out_laplace)

def get_bounding_box(mask, padding=1):
    coords = np.where(mask)
    if len(coords[0]) == 0:
        return (0, 0, 0), mask.shape
    min_coords = np.maximum(0, [np.min(c) - padding for c in coords])
    max_coords = np.minimum(mask.shape, [np.max(c) + padding + 1 for c in coords])
    return tuple(min_coords), tuple(max_coords)

def crop(x, min_bbox, max_bbox):
    return x[min_bbox[0]:max_bbox[0], min_bbox[1]:max_bbox[1], min_bbox[2]:max_bbox[2]]

def laplace(init_coords, fg, source, sink, kernelSize=3, convergence_threshold=1e-4, max_iters=500, initial_lr=1, lr_decay=0.995, verbose=False):
    hl = np.ones([kernelSize, kernelSize, kernelSize], dtype=np.float32)
    hl = hl / np.sum(hl)
    coords = np.zeros(init_coords.shape, dtype=np.float32)
    coords[fg] = init_coords[fg].astype(np.float32)
    coords[source] = 0.0
    coords[sink] = 1.0
    update_mask = fg & ~source & ~sink
    n_update = np.sum(update_mask)
    lr = initial_lr
    last_ssd = np.inf
    min_improvement = convergence_threshold * 0.01
    stagnation_count = 0
    coords_new = np.copy(coords)
    for i in range(max_iters):
        upd_coords = fast_convolve(coords, hl)
        coords_new[update_mask] = (coords[update_mask] * (1 - lr) + upd_coords[update_mask] * lr)
        coords_new[source] = 0.0
        coords_new[sink] = 1.0
        diff = coords[update_mask] - coords_new[update_mask]
        ssd = np.sum(diff * diff)
        lr *= lr_decay
        improvement = abs(last_ssd - ssd) / max(last_ssd, 1e-10)
        if improvement < min_improvement:
            stagnation_count += 1
            if stagnation_count >= 5:
                break
        else:
            stagnation_count = 0
        last_ssd = ssd
        if ssd < convergence_threshold:
            break
        coords, coords_new = coords_new, coords
    return coords

def build_adjacency_map(F):
    n_vertices = np.max(F) + 1
    adjacency = [[] for _ in range(n_vertices)]
    for face in F:
        for i in range(3):
            for j in range(3):
                if i != j:
                    adjacency[face[i]].append(face[j])
    adjacency = [np.array(list(set(adj))) for adj in adjacency]
    return adjacency

def avg_neighbours_vectorized(adjacency, cdat, zero_indices):
    result = np.zeros(len(zero_indices))
    for i, vertex_idx in enumerate(zero_indices):
        neighbors = adjacency[vertex_idx]
        if len(neighbors) > 0:
            result[i] = np.nanmean(cdat[neighbors])
        else:
            result[i] = 0.0
    return result

def shift_surface(in_surf, in_laplace, out_surf_prefix, depth_mm=[1, 2]):
    surf = nib.load(in_surf)
    V_original = surf.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0].data.astype(np.float32)
    F = surf.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")[0].data
    laplace = nib.load(in_laplace)
    lp = laplace.get_fdata().astype(np.float32)
    affine = laplace.affine
    xres, yres, zres = affine[0, 0], affine[1, 1], affine[2, 2]
    step_size = 0.1
    dx, dy, dz = np.gradient(lp)
    dx = dx / xres
    dy = dy / yres
    dz = dz / zres
    adjacency = build_adjacency_map(F)
    translation = affine[:3, 3]
    scales = np.array([affine[0, 0], affine[1, 1], affine[2, 2]])
    inv_scales = 1.0 / scales
    V = V_original.copy()
    for depth in depth_mm:
        nsteps = int(depth / step_size)
        V_vox = (V - translation) * inv_scales
        for step in range(nsteps):
            V_int = np.round(V_vox).astype(np.int32)
            V_int[:, 0] = np.clip(V_int[:, 0], 0, dx.shape[0] - 1)
            V_int[:, 1] = np.clip(V_int[:, 1], 0, dx.shape[1] - 1)
            V_int[:, 2] = np.clip(V_int[:, 2], 0, dx.shape[2] - 1)
            stepx = dx[V_int[:, 0], V_int[:, 1], V_int[:, 2]]
            stepy = dy[V_int[:, 0], V_int[:, 1], V_int[:, 2]]
            stepz = dz[V_int[:, 0], V_int[:, 1], V_int[:, 2]]
            zero_mask = (stepx == 0) & (stepy == 0) & (stepz == 0)
            zero_indices = np.where(zero_mask)[0]
            if len(zero_indices) > 0:
                stepx[zero_indices] = avg_neighbours_vectorized(adjacency, stepx, zero_indices)
                stepy[zero_indices] = avg_neighbours_vectorized(adjacency, stepy, zero_indices)
                stepz[zero_indices] = avg_neighbours_vectorized(adjacency, stepz, zero_indices)
            magnitude = np.sqrt(stepx**2 + stepy**2 + stepz**2)
            valid_mask = magnitude > 0
            scale_factor = step_size / magnitude
            stepx[valid_mask] *= scale_factor[valid_mask]
            stepy[valid_mask] *= scale_factor[valid_mask]
            stepz[valid_mask] *= scale_factor[valid_mask]
            V[valid_mask, 0] += stepx[valid_mask]
            V[valid_mask, 1] += stepy[valid_mask]
            V[valid_mask, 2] += stepz[valid_mask]
            V_vox = (V - translation) * inv_scales
        # Save shifted surface
        out_surf = f"{out_surf_prefix}{depth:.1f}mm.surf.gii"
        new_surf = nib.GiftiImage(
            darrays=[
                nib.gifti.GiftiDataArray(data=V, intent="NIFTI_INTENT_POINTSET"),
                nib.gifti.GiftiDataArray(data=F, intent="NIFTI_INTENT_TRIANGLE"),
            ]
        )
        nib.save(new_surf, out_surf)

def calcdist(surf1, surf2):
    euclidianDistanceSquared = (surf1 - surf2) ** 2
    euclidianDistanceSummed = np.sum(euclidianDistanceSquared, axis=1)
    return np.sqrt(euclidianDistanceSummed)

def computegrad(data, dists):
    data = np.ediff1d(data)
    data[dists == 0] = 0
    dists[dists == 0] = 1
    return np.divide(data, dists)

def main():
    parser = argparse.ArgumentParser(description="Standalone blurring pipeline (no src imports)")
    parser.add_argument("--participant_id", required=True)
    parser.add_argument("--session_id", required=True)
    parser.add_argument("--features", nargs="+", required=True)
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--workbench_path", required=True)
    parser.add_argument("--tmp_dir", required=True)
    parser.add_argument("--smoothing_fwhm", type=float, default=5)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # Directory setup
    subject_output_dir = os.path.join(args.output_dir, args.participant_id, args.session_id)
    os.makedirs(subject_output_dir, exist_ok=True)
    swm_dir = os.path.join(subject_output_dir, "maps", "cortex")
    os.makedirs(swm_dir, exist_ok=True)
    struct_dir = os.path.join(subject_output_dir, "structural")
    os.makedirs(struct_dir, exist_ok=True)

    for hemi in ["L", "R"]:
        if args.verbose:
            print(f"Processing hemisphere {hemi} for {args.participant_id}/{args.session_id}...")
        # Laplace solver
        aparc_aseg = os.path.join(args.input_dir, f"{args.participant_id}_{args.session_id}_aparc+aseg.nii.gz")
        laplace_out = os.path.join(struct_dir, f"{args.participant_id}_{args.session_id}_laplace-wm.nii.gz")
        if not os.path.exists(laplace_out):
            laplace_solver(aparc_aseg, laplace_out, verbose=args.verbose)
        # Shift surface
        white_surf = os.path.join(args.input_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_label-white.surf.gii")
        shift_surface(white_surf, laplace_out, os.path.join(struct_dir, f"{args.participant_id}_{args.session_id}_{hemi}_sfwm-"), [1.0, 2.0])
        for feature in args.features:
            if args.verbose:
                print(f"  Processing feature {feature} for {hemi}")
            # Map volume to surfaces
            volumemap = os.path.join(args.input_dir, "maps", f"{args.participant_id}_{args.session_id}_space-nativepro_map-{feature}.nii.gz")
            midthickness_surf = os.path.join(args.input_dir, "surf", f"{args.participant_id}_{args.session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii")
            midthickness_func = os.path.join(args.tmp_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_{feature}_midthickness.func.gii")
            white_func = os.path.join(args.tmp_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_{feature}_white.func.gii")
            subprocess.run([
                os.path.join(args.workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                midthickness_surf,
                midthickness_func,
                "-trilinear"
            ])
            subprocess.run([
                os.path.join(args.workbench_path, "wb_command"),
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
                swm_surf = os.path.join(struct_dir, f"{args.participant_id}_{args.session_id}_{hemi}_sfwm-{surf_dist}mm.surf.gii")
                swm_func = os.path.join(args.tmp_dir, f"{args.participant_id}_{args.session_id}_{hemi}_{feature}_sfwm-{surf_dist}mm.func.gii")
                subprocess.run([
                    os.path.join(args.workbench_path, "wb_command"),
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
            # Save outputs
            raw_data_path = os.path.join(swm_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-raw.func.gii")
            dist_file_path = os.path.join(swm_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-dist.func.gii")
            grad_file_path = os.path.join(swm_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_desc-grad.func.gii")
            smoothed_data_path = os.path.join(swm_dir, f"{args.participant_id}_{args.session_id}_hemi-{hemi}_feature-{feature}-blur_surf-fsnative_smooth-{args.smoothing_fwhm}mm.func.gii")
            data_array = nib.gifti.gifti.GiftiDataArray(data=surface_data, intent="NIFTI_INTENT_NORMAL")
            gii_data = nib.gifti.GiftiImage(darrays=[data_array])
            nib.save(gii_data, raw_data_path)
            del data_array, gii_data
            gc.collect()
            dist_array = nib.gifti.gifti.GiftiDataArray(data=distances, intent="NIFTI_INTENT_NORMAL")
            gii_dist = nib.gifti.GiftiImage(darrays=[dist_array])
            nib.save(gii_dist, dist_file_path)
            del distances, dist_array, gii_dist
            gc.collect()
            grad_array = nib.gifti.gifti.GiftiDataArray(data=gradients, intent="NIFTI_INTENT_NORMAL")
            gii_grad = nib.gifti.GiftiImage(darrays=[grad_array])
            nib.save(gii_grad, grad_file_path)
            del gradients, grad_array, gii_grad
            gc.collect()
            # Smoothing
            subprocess.run([
                os.path.join(args.workbench_path, "wb_command"),
                "-metric-smoothing",
                midthickness_surf,
                raw_data_path,
                str(args.smoothing_fwhm),
                smoothed_data_path
            ], check=False)
            # Set structure
            subprocess.run([
                os.path.join(args.workbench_path, "wb_command"),
                "-set-structure",
                smoothed_data_path,
                "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT"
            ], check=False)
            del surface_data
            gc.collect()
            if args.verbose:
                print(f"    Saved {feature} blur data, distances, gradients, and applied {args.smoothing_fwhm}mm smoothing for {hemi}")
if __name__ == "__main__":
    main() 