import nibabel as nib
import numpy as np
import gc
from src.utils import reshape_distances

# This script is now intended to be run as a Snakemake script, using snakemake.input and snakemake.output

def load_gifti_data(path):
    return nib.load(path).darrays[0].data.astype(np.float32)

def load_gifti_surface(path):
    return nib.load(path).darrays[0].data.astype(np.float32)

def main():
    # Inputs from snakemake
    midthickness_func = snakemake.input["midthickness_func"]
    white_func = snakemake.input["white_func"]
    sfwm_funcs = snakemake.input["sfwm_funcs"]
    midthickness_surf = snakemake.input["midthickness_surf"]
    white_surf = snakemake.input["white_surf"]
    sfwm_surfs = snakemake.input["sfwm_surfs"]
    raw_out = snakemake.output["raw"]
    dist_out = snakemake.output["dist"]
    grad_out = snakemake.output["grad"]

    # Load functional data
    func_files = [midthickness_func, white_func] + list(sfwm_funcs)
    surf_files = [midthickness_surf, white_surf] + list(sfwm_surfs)
    surface_data = [load_gifti_data(f) for f in func_files]
    surface_coords = [load_gifti_surface(f) for f in surf_files]
    num_vertices = len(surface_data[0])
    n_surfs = len(surface_data)
    # Stack data: shape (num_vertices, n_surfs)
    data_stack = np.stack(surface_data, axis=1)
    # Calculate distances and gradients between consecutive surfaces
    distances = np.zeros((num_vertices, n_surfs-1), dtype=np.float32)
    gradients = np.zeros((num_vertices, n_surfs-1), dtype=np.float32)
    for i in range(n_surfs-1):
        dist = np.sqrt(np.sum((surface_coords[i] - surface_coords[i+1])**2, axis=1)).astype(np.float32)
        distances[:, i] = dist
        with np.errstate(divide='ignore', invalid='ignore'):
            grad = np.where(dist > 0, (data_stack[:, i+1] - data_stack[:, i]) / dist, 0).astype(np.float32)
        gradients[:, i] = grad
    # Optionally reshape distances if needed
    distances = reshape_distances(distances)
    # Save outputs as GIFTI
    from nibabel import gifti
    raw_array = gifti.GiftiDataArray(data=data_stack, intent="NIFTI_INTENT_NORMAL")
    nib.save(gifti.GiftiImage(darrays=[raw_array]), raw_out)
    dist_array = gifti.GiftiDataArray(data=distances, intent="NIFTI_INTENT_NORMAL")
    nib.save(gifti.GiftiImage(darrays=[dist_array]), dist_out)
    grad_array = gifti.GiftiDataArray(data=gradients, intent="NIFTI_INTENT_NORMAL")
    nib.save(gifti.GiftiImage(darrays=[grad_array]), grad_out)
    gc.collect()

if __name__ == "__main__":
    main() 