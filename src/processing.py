import os
import subprocess
import nibabel as nib
import numpy as np
import shutil
from src.sWM import laplace_solver, surface_generator
from src.utils import reshape_distances
import scipy
import tempfile

def _process_single_subject(
    subject,
    features,
    blur_features,
    base_features,
    output_directory,
    cortical_smoothing,
    hippocampal_smoothing,
    env,
    micapipe_directory,
    hippunfold_directory,
    freesurfer_directory,
    cortex,
    hippocampus,
    subcortical,
    verbose=False,  # Reduce verbosity in parallel processing to avoid mixed output
    valid_subjects_dict=None  # Dictionary of valid subjects by feature/structure
):
    """
    Process a single subject. This function is called by each parallel worker.
    
    Returns:
    --------
    bool
        True if processing was successful, False otherwise
    """
    import shutil
    
    participant_id, session_id = subject
    
    # Create session-specific tmp directory
    session_tmp_dir = os.path.join(output_directory, participant_id, session_id, f"tmp_{participant_id}_{session_id}")
    
    try:
        os.makedirs(session_tmp_dir, exist_ok=True)
        
        # Apply blurring to features that need it
        if blur_features:
            # Check if all base features are available for this subject
            all_base_available = True
            if valid_subjects_dict:
                for feature in base_features:
                    if subject not in valid_subjects_dict.get(feature, {}).get('all', []):
                        all_base_available = False
                        break
            
            if all_base_available:
                apply_blurring(
                    participant_id=participant_id,
                    session_id=session_id,
                    features=base_features,
                    output_directory=output_directory,
                    workbench_path=env.connectome_workbench_path,
                    micapipe_directory=micapipe_directory,
                    freesurfer_directory=freesurfer_directory,
                    tmp_dir=session_tmp_dir,
                    verbose=verbose
                )

        # Process cortical features if cortex is enabled
        if cortex and (not valid_subjects_dict or subject in valid_subjects_dict.get('structures', {}).get('cortex', [])):
            # Filter features to only those valid for this subject's cortex
            features_to_process = []
            if valid_subjects_dict:
                for f in features:
                    if subject in valid_subjects_dict.get(f, {}).get('structures', {}).get('cortex', []):
                        features_to_process.append(f)
            else:
                features_to_process = features
            
            if features_to_process:
                apply_cortical_processing(
                    participant_id=participant_id,
                    session_id=session_id,
                    features=features_to_process,
                    output_directory=output_directory,
                    workbench_path=env.connectome_workbench_path,
                    micapipe_directory=micapipe_directory,
                    tmp_dir=session_tmp_dir,
                    cortical_smoothing=cortical_smoothing,
                    resolutions=["32k", "5k"],
                    labels=["midthickness", "white"],
                    verbose=verbose
                )
        
        # If hippocampus is enabled, process hippocampal data
        if hippocampus and hippunfold_directory and (not valid_subjects_dict or subject in valid_subjects_dict.get('structures', {}).get('hippocampus', [])):
            # Filter to non-blur features valid for this subject's hippocampus
            features_to_process = []
            if valid_subjects_dict:
                for f in features:
                    if not f.endswith("-blur") and subject in valid_subjects_dict.get(f, {}).get('structures', {}).get('hippocampus', []):
                        features_to_process.append(f)
            else:
                features_to_process = [f for f in features if not f.endswith("-blur")]
            
            if features_to_process:
                apply_hippocampal_processing(
                    participant_id=participant_id,
                    session_id=session_id,
                    features=features_to_process,
                    output_directory=output_directory,
                    workbench_path=env.connectome_workbench_path,
                    micapipe_directory=micapipe_directory,
                    hippunfold_directory=hippunfold_directory,
                    tmp_dir=session_tmp_dir,
                    smoothing_fwhm=hippocampal_smoothing,
                    verbose=verbose
                )

        # If subcortex is enabled, extract subcortical stats
        if subcortical and freesurfer_directory and (not valid_subjects_dict or subject in valid_subjects_dict.get('structures', {}).get('subcortical', [])):
            # Filter to non-blur features valid for this subject's subcortical structures
            features_to_process = []
            if valid_subjects_dict:
                for f in features:
                    if not f.endswith("-blur") and subject in valid_subjects_dict.get(f, {}).get('structures', {}).get('subcortical', []):
                        features_to_process.append(f)
            else:
                features_to_process = [f for f in features if not f.endswith("-blur")]
            
            if features_to_process:
                apply_subcortical_processing(
                    participant_id=participant_id,
                    session_id=session_id,
                    features=features_to_process,
                    output_directory=output_directory,
                    micapipe_directory=micapipe_directory,
                    freesurfer_directory=freesurfer_directory,
                    verbose=verbose
                )
        
        return True
        
    except Exception as e:
        if verbose:
            print(f"Error processing {participant_id}/{session_id}: {e}")
        return False
        
    finally:
        # Clean up session-specific tmp directory
        if os.path.exists(session_tmp_dir):
            try:
                shutil.rmtree(session_tmp_dir)
            except Exception as e:
                if verbose:
                    print(f"Warning: Could not clean up {session_tmp_dir}: {e}")

def fixmatrix(path, inputmap, outputmap, basemap, BIDS_ID, temppath, wb_path, mat_path):
    """
    Transform a volume using an affine transformation matrix.
    
    Parameters:
    -----------
    path : str
        Path to the directory containing transformation matrices
    inputmap : str
        Path to input volume (can be MGZ or NIfTI)
    outputmap : str
        Path where the transformed volume will be saved
    basemap : str
        Path to reference volume for resampling
    BIDS_ID : str
        Subject ID used in filenames
    temppath : str
        Directory for temporary files
    wb_path : str
        Path to Connectome Workbench executables
    mat_path : str
        Base name of the transformation matrix file
    """
    # Load the .mat file
    mat = scipy.io.loadmat(
        os.path.join(
            path,
            "xfm",
            f"{BIDS_ID}_{mat_path}.mat",
        )
    )

    # Extract variables from the .mat file
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

    print(m_matrix)
    with open(
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
        "w",
    ) as f:
        for row in m_matrix:
            f.write(" ".join(map(str, row)) + "\n")

    # Check if inputmap is an MGZ file
    if inputmap.endswith('.mgz'):
        # Create a temporary NIfTI file for wb_command
        with tempfile.NamedTemporaryFile(suffix='.nii.gz', dir=temppath, delete=False) as tmp:
            temp_nifti = tmp.name
            
        # Load and save as NIfTI (wb_command requires NIfTI)
        img = nib.load(inputmap)
        nib.save(img, temp_nifti)
        
        # Use the temporary file for the wb_command
        input_for_wb = temp_nifti
    else:
        input_for_wb = inputmap

    command3 = [
        os.path.join(wb_path, "wb_command"),
        "-volume-resample",
        input_for_wb,
        basemap,
        "ENCLOSING_VOXEL",
        outputmap,
        "-affine",
        os.path.join(temppath, f"{BIDS_ID}_real_world_affine.txt"),
    ]

    subprocess.run(command3)
    
    # Clean up temporary file if created
    if inputmap.endswith('.mgz') and os.path.exists(temp_nifti):
        os.unlink(temp_nifti)

def apply_blurring(participant_id, session_id, features, output_directory, workbench_path, micapipe_directory, freesurfer_directory, tmp_dir, smoothing_fwhm=5, verbose=True):
    """
    Apply depth-dependent blurring to one or more features for a subject.
    Only generates 1mm and 2mm superficial white matter surfaces.
    Calculates intensity gradients between adjacent surfaces.
    Applies surface smoothing to the blur features directly.
    
    Parameters:
    -----------
    participant_id : str
        Subject ID (e.g., "sub-HC001")
    session_id : str
        Session ID (e.g., "ses-01")
    features : str or list
        Feature(s) to process (e.g., "FLAIR", "qT1" or ["FLAIR", "qT1"])
    output_directory : str
        Directory to store output files
    workbench_path : str
        Path to Connectome Workbench executables
    micapipe_directory : str
        Path to micapipe outputs
    freesurfer_directory : str
        Path to FreeSurfer outputs
    tmp_dir : str
        Directory for temporary files
    smoothing_fwhm : float, default=5
        Full width at half maximum for smoothing (in mm)
    verbose : bool, default=True
        Print detailed information during processing
    """
    # Convert single feature to list if needed
    if isinstance(features, str):
        features = [features]

    # Create directories for outputs
    subject_output_dir = os.path.join(output_directory, participant_id, session_id)
    os.makedirs(subject_output_dir, exist_ok=True)
    swm_dir = os.path.join(subject_output_dir, "maps", "cortex")
    os.makedirs(swm_dir, exist_ok=True)
    struct_dir = os.path.join(subject_output_dir, "structural")
    os.makedirs(struct_dir, exist_ok=True)
    
    # Input directory for the subject
    input_dir = os.path.join(micapipe_directory, participant_id, session_id)
    bids_id = f"{participant_id}_{session_id}"
    
    # Process each hemisphere
    for hemi in ["L", "R"]:
        if verbose:
            print(f"    Processing hemisphere {hemi} for {participant_id}/{session_id}...")
        
        # Define paths
        freesurfer_path = os.path.join(freesurfer_directory, f"{participant_id}_{session_id}")
        temp_parc_path = os.path.join(tmp_dir, f"{participant_id}_{session_id}_{hemi}_surf-fsnative_label-temp.nii.gz")
        output_path = os.path.join(struct_dir, f"{participant_id}_{session_id}-laplace.nii.gz")
        
        # Get FreeSurfer parcellation path (using .mgz directly)
        aparc_path = os.path.join(freesurfer_path, "mri", "aparc+aseg.mgz")
        
        # Transform FreeSurfer parcellation to native space if needed (done once per hemisphere)
        if not os.path.exists(temp_parc_path):
            # Pass the MGZ file directly to fixmatrix
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
            
            # Solve Laplace equation
            laplace_solver.solve_laplace(temp_parc_path, output_path)
            
            # Generate ONLY 1mm and 2mm shifted surfaces
            surface_generator.shift_surface(
                os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii"),
                output_path,
                os.path.join(struct_dir, f"{participant_id}_{session_id}_{hemi}_sfwm-"),
                [1.0, 2.0], 
            )
        
        # Process each feature
        for feature in features:
            if verbose:
                print(f"      Processing feature {feature} for {participant_id}/{session_id}, hemisphere {hemi}...")
            
            # Standardize feature name for consistent file paths
            feature_lower = feature.lower()
            
            # Set up volume map path based on feature type
            if feature_lower == "t1map" or feature_lower == "qt1":
                volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz")
                # Use qT1 as the feature name for consistency
                feature_name = "qT1"  # Changed from using feature_lower
            elif feature_lower in ["adc", "fa"]:
                volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_model-DTI_map-{feature_lower}.nii.gz")
                feature_name = feature_lower
            else:
                volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-{feature_lower}.nii.gz")
                feature_name = feature_lower
            
            # Skip if the volume map doesn't exist
            if not os.path.exists(volumemap):
                if verbose:
                    print(f"      Warning: Volume map not found for feature {feature}, skipping")
                continue
                
            # Generate multi-depth data
            surfarr = []
            
            # Load midthickness surface and map data
            midthickness_surf = os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-midthickness.surf.gii")
            midthickness_func = os.path.join(tmp_dir, f"{participant_id}_{session_id}_hemi-{hemi}_{feature}_midthickness.func.gii")
            
            # Map volume to midthickness surface
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                midthickness_surf,
                midthickness_func,
                "-trilinear"
            ])
            
            # Load data from midthickness
            pial_data = nib.load(midthickness_func).darrays[0].data
            pial_surface = nib.load(midthickness_surf).darrays[0].data
            surfarr.append([pial_data, pial_surface])
            
            # Load white matter boundary surface and map data
            white_surf = os.path.join(input_dir, "surf", f"{participant_id}_{session_id}_hemi-{hemi}_space-nativepro_surf-fsnative_label-white.surf.gii") 
            white_func = os.path.join(tmp_dir, f"{participant_id}_{session_id}_hemi-{hemi}_{feature}_white.func.gii")
            
            # Map volume to white surface
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                white_surf,
                white_func,
                "-trilinear"
            ])
            
            # Load data from white surface
            white_data = nib.load(white_func).darrays[0].data
            white_surface = nib.load(white_surf).darrays[0].data
            surfarr.append([white_data, white_surface])
            
            # Map data to ONLY 1mm and 2mm shifted white matter surfaces
            for surf_dist in [1.0, 2.0]:  # Only process 1mm and 2mm
                swm_surf = os.path.join(struct_dir, f"{participant_id}_{session_id}_{hemi}_sfwm-{surf_dist}mm.surf.gii")
                swm_func = os.path.join(swm_dir, f"{participant_id}_{session_id}_{hemi}_{feature}_sfwm-{surf_dist}mm_metric.func.gii")
                
                # Map data to shifted surface
                subprocess.run([
                    os.path.join(workbench_path, "wb_command"),
                    "-volume-to-surface-mapping",
                    volumemap,
                    swm_surf,
                    swm_func,
                    "-trilinear"
                ])
                
                # Add to surface array
                surfarr.append([
                    nib.load(swm_func).darrays[0].data,
                    nib.load(swm_surf).darrays[0].data
                ])
            
            # Calculate distances between surfaces
            distances = np.zeros(shape=(len(pial_data), len(surfarr) - 1))
            dataArr_nonmode = np.zeros(shape=(len(pial_data), len(surfarr)), dtype=np.float32)
            
            # Calculate gradients between adjacent surfaces (intensity change / distance)
            gradients = np.zeros(shape=(len(pial_data), len(surfarr) - 1), dtype=np.float32)
            
            # Calculate Euclidean distances and intensity values
            for e, ds in enumerate(surfarr):
                data, surf = ds
                dataArr_nonmode[:, e] = data
                
                if e == len(surfarr) - 1:
                    break
                    
                nextdata, nextsurt = surfarr[e + 1]
                
                # Calculate distance between surfaces
                distance = np.sqrt(np.sum((surf - nextsurt) ** 2, axis=1))  # Euclidean distance
                distances[:, e] = distance
                
                # Calculate intensity gradient (change in intensity / distance)
                intensity_change = nextdata - data  # Change in intensity
                
                # Avoid division by zero
                with np.errstate(divide='ignore', invalid='ignore'):
                    gradient = np.where(distance > 0, intensity_change / distance, 0)
                    gradients[:, e] = gradient
            
            # Reshape distances and save outputs with STANDARDIZED NAMING
            distances = reshape_distances(distances)
            
            # Save feature data across surfaces (unsmoothed version) - UPDATED NAMING
            raw_data_path = os.path.join(swm_dir, 
                f"{participant_id}_{session_id}_hemi-{hemi}_feature-{feature_name}-blur_surf-fsnative_desc-raw.func.gii")
            data_non_grad = nib.gifti.gifti.GiftiDataArray(
                data=dataArr_nonmode,
                intent="NIFTI_INTENT_NORMAL"
            )
            gii_non_grad = nib.gifti.GiftiImage(darrays=[data_non_grad])
            nib.save(gii_non_grad, raw_data_path)
            
            # Save distance data - UPDATED NAMING
            data_dist = nib.gifti.gifti.GiftiDataArray(
                data=distances.astype(np.float32),
                intent="NIFTI_INTENT_NORMAL"
            )
            gii_dist = nib.gifti.GiftiImage(darrays=[data_dist])
            nib.save(
                gii_dist,
                os.path.join(swm_dir, f"{participant_id}_{session_id}_hemi-{hemi}_feature-{feature_name}-blur_surf-fsnative_desc-dist.func.gii")
            )
            
            # Save gradient data - UPDATED NAMING
            data_grad = nib.gifti.gifti.GiftiDataArray(
                data=gradients.astype(np.float32),
                intent="NIFTI_INTENT_NORMAL"
            )
            gii_grad = nib.gifti.GiftiImage(darrays=[data_grad])
            nib.save(
                gii_grad,
                os.path.join(swm_dir, f"{participant_id}_{session_id}_hemi-{hemi}_feature-{feature_name}-blur_surf-fsnative_desc-grad.func.gii")
            )
            
            # Apply smoothing directly in native space - UPDATED NAMING
            smoothed_data_path = os.path.join(swm_dir, 
                f"{participant_id}_{session_id}_hemi-{hemi}_feature-{feature_name}-blur_surf-fsnative_smooth-{smoothing_fwhm}mm.func.gii")
            
            # Apply smoothing to midthickness surface since it's better for smoothing
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-metric-smoothing",
                midthickness_surf,
                raw_data_path,
                str(smoothing_fwhm),
                smoothed_data_path
            ], check=False)
            
            # Set structure attribute
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                smoothed_data_path,
                "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
            ], check=False)
            
            if verbose:
                print(f"      Saved {feature} blur data, distances, gradients, and applied {smoothing_fwhm}mm smoothing for {hemi}")

def apply_hippocampal_processing(
    participant_id, 
    session_id, 
    features, 
    output_directory, 
    workbench_path, 
    micapipe_directory, 
    hippunfold_directory, 
    tmp_dir, 
    smoothing_fwhm=2, 
    resolution="0p5mm", 
    verbose=True
):
    """
    Apply processing to hippocampal features for a subject.
    Maps volume features to hippocampal surfaces and applies smoothing.
    
    Parameters:
    -----------
    participant_id : str
        Subject ID (e.g., "sub-HC001")
    session_id : str
        Session ID (e.g., "ses-01")
    features : str or list
        Feature(s) to process (e.g., "FLAIR", "qT1" or ["FLAIR", "qT1"])
    output_directory : str
        Directory to store output files
    workbench_path : str
        Path to Connectome Workbench executables
    micapipe_directory : str
        Path to micapipe outputs
    hippunfold_directory : str
        Path to hippunfold outputs
    tmp_dir : str
        Directory for temporary files
    smoothing_fwhm : float, default=2
        Full width at half maximum for smoothing (in mm)
    resolution : str, default="0p5mm"
        Resolution for hippocampal surfaces (e.g., "0p5mm")
    verbose : bool, default=True
        Print detailed information during processing
    """
    # Convert single feature to list if needed
    if isinstance(features, str):
        features = [features]

    # Create directories for outputs
    subject_output_dir = os.path.join(output_directory, participant_id, session_id)
    os.makedirs(subject_output_dir, exist_ok=True)
    hippocampus_dir = os.path.join(subject_output_dir, "maps", "hippocampus")
    os.makedirs(hippocampus_dir, exist_ok=True)
    
    # Input directories for the subject
    input_dir = os.path.join(micapipe_directory, participant_id, session_id)
    hipp_dir = os.path.join(hippunfold_directory, "hippunfold", participant_id, session_id)
    
    if verbose:
        print(f"  Processing hippocampal data for {participant_id}/{session_id}")
    
    # Validate hippunfold directory exists
    if not os.path.exists(hipp_dir):
        if verbose:
            print(f"  Warning: Hippunfold directory not found at {hipp_dir}")
        return False
    
    # Process each feature
    for feature in features:
        if verbose:
            print(f"    Processing hippocampal feature {feature} for {participant_id}/{session_id}")
        
        # Skip features with "blur" suffix
        if feature.lower().endswith("-blur"):
            continue
        
        # Set up volume map path based on feature type
        if feature.lower() == "thickness":
            intermediate_file = os.path.join(
                hipp_dir, 
                "surf",
                f"{participant_id}_{session_id}_hemi-L_den-{resolution}_label-hipp_midthickness_thickness.func.gii"
            )
            output_feat = "thickness"
        elif feature.lower() == "t1map" or feature.lower() == "qt1":
            volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz")
            output_feat = "qT1"  # Changed from "T1map" to "qT1"
        elif feature.lower() in ["adc", "fa"]:
            volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_model-DTI_map-{feature.upper()}.nii.gz")
            output_feat = feature.upper()
        else:
            volumemap = os.path.join(input_dir, "maps", f"{participant_id}_{session_id}_space-nativepro_map-{feature.lower()}.nii.gz")
            output_feat = feature.lower()
        
        # Skip if volume map doesn't exist
        if not os.path.exists(volumemap):
            if verbose:
                print(f"    Warning: Volume map not found for feature {feature}: {volumemap}")
            continue
        
        # Process each hemisphere
        n_processed = 0
        for hemi in ["L", "R"]:
            # Define paths for hippocampal surfaces
        
            # Set up paths
            surf_file = os.path.join(
                hipp_dir, 
                "surf", 
                f"{participant_id}_{session_id}_hemi-{hemi}_space-T1w_den-{resolution}_label-hipp_midthickness.surf.gii"
            )
            inner_surf_file = os.path.join(
                hipp_dir, 
                "surf", 
                f"{participant_id}_{session_id}_hemi-{hemi}_space-T1w_den-{resolution}_label-hipp_inner.surf.gii"
            )
            outer_surf_file = os.path.join(
                hipp_dir, 
                "surf", 
                f"{participant_id}_{session_id}_hemi-{hemi}_space-T1w_den-{resolution}_label-hipp_outer.surf.gii"
            )
            
            # Check if surface exists
            if not os.path.exists(surf_file):
                if verbose:
                    print(f"    Warning: Hippocampal surface not found: {surf_file}")
                continue
            
            if feature.lower() != "thickness":
                # Define intermediate and output files
                intermediate_file = os.path.join(
                    tmp_dir,
                    f"{participant_id}_{session_id}_hemi-{hemi}_den-{resolution}_label-hipp_midthickness_{feature}.func.gii"
                )
            output_file = os.path.join(
                hippocampus_dir,
                f"{participant_id}_{session_id}_hemi-{hemi}_den-{resolution}_label-hipp_midthickness_feature-{output_feat}_smooth-{smoothing_fwhm}mm.func.gii"
            )
            
            # Map volume to surface
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-volume-to-surface-mapping",
                volumemap,
                surf_file,
                intermediate_file,
                "-ribbon-constrained",
                inner_surf_file,
                outer_surf_file,
            ], check=False)
            
            # Apply smoothing
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-metric-smoothing",
                surf_file,
                intermediate_file,
                str(smoothing_fwhm),
                output_file
            ], check=False)
            
            # Set structure attribute
            subprocess.run([
                os.path.join(workbench_path, "wb_command"),
                "-set-structure",
                output_file,
                "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
            ], check=False)
            
            if os.path.exists(output_file):
                n_processed += 1
        
        if verbose:
            print(f"    Completed hippocampal processing for feature {feature}: {n_processed} surfaces mapped")
    
    if verbose:
        print(f"  Completed hippocampal processing for {participant_id}/{session_id}")
    
    return True

def apply_subcortical_processing(
    participant_id,
    session_id,
    features,
    output_directory,
    micapipe_directory,
    freesurfer_directory,
    verbose=True
):
    """
    Apply subcortical processing to features for a subject.
    Maps volume features to subcortical structures and extracts statistics.
    
    Parameters:
    -----------
    participant_id : str
        Subject ID (e.g., "sub-HC001")
    session_id : str
        Session ID (e.g., "ses-01")
    features : str or list
        Feature(s) to process (e.g., "FLAIR", "qT1" or ["FLAIR", "qT1"])
    output_directory : str
        Directory to store output files
    micapipe_directory : str
        Path to micapipe outputs
    freesurfer_directory : str
        Path to FreeSurfer outputs
    verbose : bool, default=True
        Print detailed information during processing
    """
    import pandas as pd
    import numpy as np
    import nibabel as nib
    
    # Define subcortical structure labels and names
    SUBCORTICAL_LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
    SUBCORTICAL_STRUCTURES = [
        "Laccumb", "Lamyg", "Lcaud", "Lhippo", "Lpal", "Lput", "Lthal",
        "Raccumb", "Ramyg", "Rcaud", "Rhippo", "Rpal", "Rput", "Rthal",
    ]
    
    # Convert single feature to list if needed
    if isinstance(features, str):
        features = [features]
    
    # Create directories for outputs
    subject_output_dir = os.path.join(output_directory, participant_id, session_id)
    os.makedirs(subject_output_dir, exist_ok=True)
    subcortical_dir = os.path.join(subject_output_dir, "maps", "subcortical")
    os.makedirs(subcortical_dir, exist_ok=True)
    
    # Set paths
    bids_id = f"{participant_id}_{session_id}"
    subject_micapipe_dir = os.path.join(micapipe_directory, participant_id, session_id)
    subject_freesurfer_dir = os.path.join(freesurfer_directory, bids_id)
    seg_file = os.path.join(subject_micapipe_dir, "parc", f"{bids_id}_space-nativepro_T1w_atlas-subcortical.nii.gz")
    aseg_stats_file = os.path.join(subject_freesurfer_dir, "stats", "aseg.stats")
    
    if verbose:
        print(f"  Processing subcortical data for {participant_id}/{session_id}")
    
    # Validate necessary files exist
    if not os.path.exists(seg_file):
        if verbose:
            print(f"  Warning: Subcortical segmentation file not found at {seg_file}")
        return False
        
    if not os.path.exists(aseg_stats_file):
        if verbose:
            print(f"  Warning: FreeSurfer aseg.stats file not found at {aseg_stats_file}")
        return False
    
    # Process each feature
    for feature in features:
        if verbose:
            print(f"    Processing subcortical feature {feature} for {participant_id}/{session_id}")
        
        # Skip features with "blur" suffix
        if feature.lower().endswith("-blur"):
            continue
        
        # Set up feature type and output name
        feat_lower = feature.lower()
        if feat_lower == "thickness" or feat_lower == "volume":
            # Volume measurements come from aseg.stats
            output_file = os.path.join(subcortical_dir, f"{bids_id}_feature-volume.csv")
            
            if verbose:
                print(f"    Extracting volume data from FreeSurfer stats")
                
            # Extract volumes from aseg.stats
            try:
                # Read aseg.stats and extract volume data
                df_volumes = pd.read_csv(
                    aseg_stats_file,
                    comment="#",
                    index_col=0,
                    sep=r"\s+",
                    header=None,
                    usecols=[1, 3, 4],
                    names=["label", "volume", "structure"]
                )
                
                # Extract volumes for specified labels
                volumes = df_volumes.loc[SUBCORTICAL_LABELS].volume.to_numpy()
                
                # Create dataframe with subcortical volumes
                df = pd.DataFrame(
                    volumes[None], 
                    columns=SUBCORTICAL_STRUCTURES, 
                    index=pd.Index([bids_id], name="SubjID")
                )
                
                # Save to CSV
                df.to_csv(output_file)
                
                if verbose:
                    print(f"    Saved volume data to {output_file}")
            
            except Exception as e:
                if verbose:
                    print(f"    Error extracting volume data: {e}")
                continue
        else:
            # Feature intensity measurements
            if feat_lower == "t1map" or feat_lower == "qt1":
                input_file = os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-T1map.nii.gz")
                output_feat = "qT1"
            elif feat_lower in ["adc", "fa"]:
                input_file = os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_model-DTI_map-{feat_lower.upper()}.nii.gz")
                output_feat = feat_lower.upper()
            elif feat_lower == "flair":
                input_file = os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-flair.nii.gz")
                output_feat = "flair"
            else:
                input_file = os.path.join(subject_micapipe_dir, "maps", f"{bids_id}_space-nativepro_map-{feat_lower}.nii.gz")
                output_feat = feat_lower
            
            output_file = os.path.join(subcortical_dir, f"{bids_id}_feature-{output_feat}.csv")
            
            # Skip if input file doesn't exist
            if not os.path.exists(input_file):
                if verbose:
                    print(f"    Warning: Input file not found for feature {feature}: {input_file}")
                continue
            
            try:
                # Load segmentation and feature image
                seg = np.asarray(nib.load(seg_file).dataobj)
                img = np.asarray(nib.load(input_file).dataobj)
                
                # Create dataframe for results
                df = pd.DataFrame(
                    np.nan, 
                    columns=SUBCORTICAL_STRUCTURES, 
                    index=pd.Index([bids_id], name="SubjID")
                )
                
                # Extract mean values for each structure
                for i, label in enumerate(SUBCORTICAL_LABELS):
                    region_mask = seg == label
                    if np.any(region_mask):
                        df.loc[bids_id, SUBCORTICAL_STRUCTURES[i]] = np.nanmean(img[region_mask])
                    else:
                        if verbose:
                            print(f"    Warning: No voxels found for structure {SUBCORTICAL_STRUCTURES[i]} (label {label})")
                
                # Save to CSV
                df.to_csv(output_file)
                
                if verbose:
                    print(f"    Saved {feature} subcortical data to {output_file}")
            
            except Exception as e:
                if verbose:
                    print(f"    Error processing {feature} for subcortical structures: {e}")
                continue
    
    if verbose:
        print(f"  Completed subcortical processing for {participant_id}/{session_id}")
    
    return True

def apply_cortical_processing(
    participant_id,
    session_id,
    features,
    output_directory,
    workbench_path,
    micapipe_directory,
    tmp_dir,
    cortical_smoothing=5,
    resolutions=["32k", "5k"],
    labels=["midthickness", "white"],
    verbose=True
):
    """
    Apply processing to cortical features for a subject.
    Maps volume features to cortical surfaces and applies smoothing.
    For blur features, uses pre-smoothed data from apply_blurring.
    
    Parameters:
    -----------
    participant_id : str
        Subject ID (e.g., "sub-HC001")
    session_id : str
        Session ID (e.g., "ses-01")
    features : str or list
        Feature(s) to process (e.g., "FLAIR", "qT1" or ["FLAIR", "qT1"])
    output_directory : str
        Directory to store output files
    workbench_path : str
        Path to Connectome Workbench executables
    micapipe_directory : str
        Path to micapipe outputs
    tmp_dir : str
        Directory for temporary files
    cortical_smoothing : float, default=5
        Full width at half maximum for smoothing (in mm)
    resolutions : list, default=["32k", "5k"]
        Resolutions for cortical surfaces (e.g., ["32k", "5k"])
    labels : list, default=["midthickness", "white"]
        Labels for cortical surfaces (e.g., ["midthickness", "white"])
    verbose : bool, default=True
        Print detailed information during processing
    """
    import numpy as np
    import nibabel as nib
    from src.utils import reshape_distances
    
    # Path to standard sphere templates
    data_dir = "data"
    
    # Convert single feature to list if needed
    if isinstance(features, str):
        features = [features]

    # Create directories for outputs
    subject_output_dir = os.path.join(output_directory, participant_id, session_id)
    os.makedirs(subject_output_dir, exist_ok=True)
    cortex_dir = os.path.join(subject_output_dir, "maps", "cortex")
    os.makedirs(cortex_dir, exist_ok=True)
    struct_dir = os.path.join(subject_output_dir, "structural")
    os.makedirs(struct_dir, exist_ok=True)
    
    # Input directory for the subject
    subject_dir = os.path.join(micapipe_directory, participant_id, session_id)
    surf_dir = os.path.join(subject_dir, "surf")
    input_dir = os.path.join(subject_dir, "maps")
    bids_id = f"{participant_id}_{session_id}"
    
    if verbose:
        print(f"  Processing cortical data for {participant_id}/{session_id}")
    
    # Map features to their file naming conventions
    feature_mapping = {
        "thickness": {"input": "thickness", "output": "thickness"},
        "flair": {"input": "flair", "output": "flair"},
        "adc": {"input": "ADC", "output": "ADC"},
        "fa": {"input": "FA", "output": "FA"},
        "qt1": {"input": "T1map", "output": "qT1"},  # Changed output from "T1map" to "qT1"
        "qt1-blur": {"input": "T1map", "output": "qT1-blur"}  # Added for consistency with blur features
    }
    
    # Process each feature
    for feature in features:
        feat_lower = feature.lower()
        is_blur = feat_lower.endswith("-blur")
        
        # Get base feature name and mapping
        if is_blur:
            base_feat = feat_lower.replace("-blur", "")
            if base_feat in feature_mapping:
                input_feat = feature_mapping[base_feat]["input"]
                output_feat = feature_mapping[base_feat]["output"] + "-blur"
            else:
                input_feat = base_feat
                output_feat = base_feat + "-blur"
        else:
            if feat_lower in feature_mapping:
                input_feat = feature_mapping[feat_lower]["input"]
                output_feat = feature_mapping[feat_lower]["output"]
            else:
                input_feat = feat_lower
                output_feat = feat_lower
        
        if verbose:
            print(f"    Processing cortical feature {feature} ({input_feat} → {output_feat})")
        
        # Process each hemisphere
        for hemi in ["L", "R"]:
            # Use standard sphere from data directory
            sphere_native = os.path.join(surf_dir, f"{bids_id}_hemi-{hemi}_surf-fsnative_label-sphere.surf.gii")
            
            # Skip if native sphere doesn't exist
            if not os.path.exists(sphere_native):
                if verbose:
                    print(f"      Warning: Native sphere not found: {sphere_native}")
                continue
            
            # Check and process each resolution
            for resolution in resolutions:
                # Use standard fsLR sphere from data directory
                sphere_fsLR = os.path.join(data_dir, f"fsLR-{resolution}.{hemi}.sphere.reg.surf.gii")
                
                # Process each label (surface type) - For blur features, only use midthickness
                label_list = ["midthickness"] if is_blur else labels
                for label in label_list:
                    if verbose:
                        print(f"      Processing {feature} for hemi-{hemi}, {resolution}, label-{label}")
                    
                    # Set up surface file path
                    surf_file = os.path.join(
                        surf_dir, 
                        f"{bids_id}_hemi-{hemi}_space-nativepro_surf-fsLR-{resolution}_label-{label}.surf.gii"
                    )
                    
                    # Skip if surface doesn't exist
                    if not os.path.exists(surf_file):
                        if verbose:
                            print(f"      Warning: Surface file not found: {surf_file}")
                        continue
                    
                    # Define output file
                    output_file = os.path.join(
                        cortex_dir,
                        f"{bids_id}_hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{cortical_smoothing}mm.func.gii"
                    )
                    
                    # Handle blur features differently - use existing pre-smoothed blurred data
                    if is_blur:
                        if verbose:
                            print(f"      Using pre-smoothed blur data for {feature}")
                        
                        # Use already smoothed blur files generated by apply_blurring - UPDATED PATH
                        blur_file = os.path.join(
                            cortex_dir, 
                            f"{participant_id}_{session_id}_hemi-{hemi}_feature-{output_feat}_surf-fsnative_smooth-{cortical_smoothing}mm.func.gii"
                        )
                        
                        if not os.path.exists(blur_file):
                            if verbose:
                                print(f"      Warning: Smoothed blur file not found: {blur_file}")
                                print(f"      Ensure apply_blurring has been run first for {base_feat}")
                            continue
                        
                        # Map the already-smoothed blur data to fsLR space (no additional smoothing needed)
                        subprocess.run([
                            os.path.join(workbench_path, "wb_command"),
                            "-metric-resample",
                            blur_file,
                            sphere_native,
                            sphere_fsLR,
                            "BARYCENTRIC",
                            output_file,
                            "-largest"
                        ], check=False)
                        
                    else:
                        # Regular feature processing (non-blur)
                        # Define input file based on feature type
                        if feat_lower == "thickness":
                            input_file = os.path.join(
                                input_dir,
                                f"{bids_id}_hemi-{hemi}_surf-fsLR-{resolution}_label-{input_feat}.func.gii"
                            )
                        else:
                            input_file = os.path.join(
                                input_dir,
                                f"{bids_id}_hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_{input_feat}.func.gii"
                            )
                        
                        # Skip if input file doesn't exist
                        if not os.path.exists(input_file):
                            if verbose:
                                print(f"      Warning: Input file not found: {input_file}")
                            continue
                        
                        # Apply smoothing
                        subprocess.run([
                            os.path.join(workbench_path, "wb_command"),
                            "-metric-smoothing",
                            surf_file,
                            input_file,
                            str(cortical_smoothing),
                            output_file
                        ], check=False)
                    
                    # Set structure attribute
                    if os.path.exists(output_file):
                        subprocess.run([
                            os.path.join(workbench_path, "wb_command"),
                            "-set-structure",
                            output_file,
                            "CORTEX_LEFT" if hemi == "L" else "CORTEX_RIGHT",
                        ], check=False)
                        
                        if verbose:
                            print(f"      Successfully processed {output_file}")
    
    if verbose:
        print(f"  Completed cortical processing for {participant_id}/{session_id}")
    
    return True