import os
import tempfile
import subprocess
import nibabel as nib
from nibabel.gifti import GiftiImage, GiftiDataArray
import numpy as np
import shutil
import random

def project_to_native_surface(
    data_lh: np.ndarray,
    data_rh: np.ndarray, 
    native_lh_path: str,
    native_rh_path: str,
    env = None,
    resampling_method: str = "ADAP_BARY_AREA",
    DATA_PATH: str = None,
    subject_dir=None,
    participant_id=None,
    session_id=None,
    verbose: bool = False,
    tmp_dir: str = None
):
    """
    Project data from fsLR standard surfaces to native surfaces using Connectome Workbench.
    
    Parameters:
    -----------
    data_lh : numpy.ndarray
        Left hemisphere data on fsLR standard surface
    data_rh : numpy.ndarray
        Right hemisphere data on fsLR standard surface
    fsLR_sphere_lh_path : str
        Path to left hemisphere fsLR sphere surface
    fsLR_sphere_rh_path : str
        Path to right hemisphere fsLR sphere surface
    native_sphere_lh_path : str
        Path to left hemisphere native sphere surface
    native_sphere_rh_path : str
        Path to right hemisphere native sphere surface
    native_lh_path : str
        Path to left hemisphere native surface
    native_rh_path : str
        Path to right hemisphere native surface
    workbench_path : str, optional
        Path to Connectome Workbench binary (default: looks for "wb_command" in PATH)
    resampling_method : str, optional
        Method for resampling: ADAP_BARY_AREA (default), BARYCENTRIC, or CUBIC
    verbose : bool, optional
        Whether to print detailed information about the processing
        
    Returns:
    --------
    tuple
        (resampled_lh_data, resampled_rh_data) - Data resampled to native surfaces
    """
    
    workbench_path = env.connectome_workbench_path
    
    if verbose:
        print(f"Using Workbench: {workbench_path}")
        print(f"Input data shapes - LH: {data_lh.shape}, RH: {data_rh.shape}")
    
    

    temp_lh_metric = os.path.join(tmp_dir, f"fsLR_lh_{''.join([str(random.randint(0, 9)) for _ in range(10)])}.func.gii")
    temp_rh_metric = os.path.join(tmp_dir, f"fsLR_rh_{''.join([str(random.randint(0, 9)) for _ in range(10)])}.func.gii")
    resampled_lh_metric = os.path.join(tmp_dir, f"native_lh_{''.join([str(random.randint(0, 9)) for _ in range(10)])}.func.gii")
    resampled_rh_metric = os.path.join(tmp_dir, f"native_rh_{''.join([str(random.randint(0, 9)) for _ in range(10)])}.func.gii")
    
    # Create GIFTI metric files from data arrays
    for hemi, data, output_file in [
        ("LH", data_lh, temp_lh_metric),
        ("RH", data_rh, temp_rh_metric)
    ]:
        # Create metric GIFTI
        metric_gii = GiftiImage()
        data_array = GiftiDataArray(data=data.astype(np.float32))
        metric_gii.add_gifti_data_array(data_array)
        
        # Save to temporary file
        metric_gii.to_filename(output_file)
        if verbose:
            print(f"Created temporary {hemi} metric file: {output_file}")
    
    if subject_dir and participant_id and session_id:
        native_sphere_lh_path = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii"
        )
        native_sphere_rh_path = os.path.join(
            subject_dir, 
            "structural", 
            f"{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii"
        )
        fsLR_sphere_lh_path = os.path.join(
            DATA_PATH, 
            "fsLR-32k.L.sphere.reg.surf.gii"
        )
        fsLR_sphere_rh_path = os.path.join(
            DATA_PATH, 
            "fsLR-32k.R.sphere.reg.surf.gii"
        )



    # Check if files exist
    for file_path, desc in [
        (fsLR_sphere_lh_path, "fsLR left sphere"),
        (fsLR_sphere_rh_path, "fsLR right sphere"),
        (native_sphere_lh_path, "native left sphere"),
        (native_sphere_rh_path, "native right sphere"),
        (native_lh_path, "native left surface"),
        (native_rh_path, "native right surface"),
        (temp_lh_metric, "left metric"),
        (temp_rh_metric, "right metric")
    ]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing required file: {desc} at {file_path}")
    
    # Get standard midthickness surfaces from fsLR
    fsLR_midthick_lh = fsLR_sphere_lh_path.replace("sphere", "")
    fsLR_midthick_rh = fsLR_sphere_rh_path.replace("sphere", "")
    
    # If standard midthickness surfaces don't exist, use defaults
    if not os.path.exists(fsLR_midthick_lh):
        fsLR_midthick_lh = f"{DATA_PATH}/fsLR-32k.L.surf.gii"
    if not os.path.exists(fsLR_midthick_rh):
        fsLR_midthick_rh = f"{DATA_PATH}/fsLR-32k.R.surf.gii"
    
    # Define resampling commands
    cmds = [
        # Left hemisphere
        [
            os.path.join(workbench_path, "wb_command"), "-metric-resample",
            temp_lh_metric,
            fsLR_sphere_lh_path,
            native_sphere_lh_path,
            resampling_method,
            resampled_lh_metric,
            "-area-surfs",
            fsLR_midthick_lh,
            native_lh_path
        ],
        # Right hemisphere
        [
            os.path.join(workbench_path, "wb_command"), "-metric-resample",
            temp_rh_metric,
            fsLR_sphere_rh_path,
            native_sphere_rh_path,
            resampling_method,
            resampled_rh_metric,
            "-area-surfs",
            fsLR_midthick_rh,
            native_rh_path
        ]
    ]
    
    # Execute resampling for both hemispheres
    for i, (cmd, hemi) in enumerate(zip(cmds, ["Left", "Right"])):
        if verbose:
            print(f"Resampling {hemi} hemisphere data...")
            print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR in {hemi} hemisphere resampling:")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"Failed to resample {hemi} hemisphere: {result.stderr}")
        elif verbose and result.stdout:
            print(f"{hemi} resampling output: {result.stdout}")
    
    # Load resampled data
    if verbose:
        print("Loading resampled data from GIFTI files")
    
    resampled_lh_data = nib.load(resampled_lh_metric).darrays[0].data
    resampled_rh_data = nib.load(resampled_rh_metric).darrays[0].data
    
    if verbose:
        print(f"Output data shapes - LH: {resampled_lh_data.shape}, RH: {resampled_rh_data.shape}")
    
    return resampled_lh_data, resampled_rh_data
        