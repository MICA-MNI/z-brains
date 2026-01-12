import os
import numpy as np
import nibabel as nib
import pandas as pd

def calculate_wscore_maps(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=['age', 'sex'], verbose=True):
    """
    Calculate W-scores for patient data against reference data using normative modeling and save as GIFTI file.
    
    Parameters:
    -----------
    reference_data : np.ndarray
        Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
    patient_data : np.ndarray
        Patient data with shape (n_vertices,) or (n_vertices, n_depths)
    demographics_ref : pd.DataFrame
        Demographics data for reference subjects with normative columns
    demographics_pat : pd.Series or pd.DataFrame
        Demographics data for patient (single row if DataFrame)
    output_file : str
        Path to save the W-score map
    normative_columns : list, default=['age', 'sex']
        List of demographic columns to use for normative modeling
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing W-score statistics, normative data, and file path
    """
    if len(reference_data) == 0:
        raise ValueError("No reference data provided")
    
    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate W-scores")
    
    if demographics_ref.shape[0] != reference_data.shape[0]:
        raise ValueError("Demographics and reference data must have same number of subjects")
    
    # Handle multi-depth data (blur features) by averaging across depths
    if len(reference_data.shape) > 2:
        reference_data = np.mean(reference_data, axis=2)
    if len(patient_data.shape) > 1:
        patient_data = np.mean(patient_data, axis=1)
    
    n_vertices = reference_data.shape[1]
    
    # Get number of predictor variables
    n_predictors = len(normative_columns)
    
    # Prepare normative data storage - dynamic size based on number of predictors
    # Format: [intercept, coef_1, coef_2, ..., coef_n, std_residuals]
    normative_data = np.zeros((n_vertices, n_predictors + 2))
    
    # Prepare demographics data
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    
    # Convert to numpy arrays
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    
    # Remove subjects with missing demographic data
    mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[mask]
    ref_data_clean = reference_data[mask]
    
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")
    
    # Add intercept term
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    
    # Fit normative model for each vertex
    for i in range(n_vertices):
        vertex_data = ref_data_clean[:, i]
        
        # Skip vertices with all zeros
        if np.all(vertex_data == 0):
            # Initialize with zeros and set std to 1 to avoid division by zero
            normative_data[i, :] = 0
            normative_data[i, -1] = 1
            continue
        
        # Fit linear regression: y = intercept + coef_1*x_1 + coef_2*x_2 + ... + coef_n*x_n
        coefficients = np.linalg.lstsq(X_ref_with_intercept, vertex_data, rcond=None)[0]
        
        # Calculate residuals and their standard deviation
        predicted = X_ref_with_intercept @ coefficients
        residuals = vertex_data - predicted
        std_residuals = np.std(residuals)
        
        # Store all coefficients plus std_residuals
        normative_data[i, :-1] = coefficients
        normative_data[i, -1] = std_residuals
            

    
    # Calculate W-scores for patient
    wscores = np.zeros(n_vertices)
    
    for i in range(n_vertices):
        # Get coefficients and standard deviation
        coefficients = normative_data[i, :-1]  # All coefficients (intercept + predictors)
        std_residuals = normative_data[i, -1]  # Standard deviation of residuals
        
        # Calculate expected value: intercept + sum(coef_j * X_pat[j])
        expected = coefficients[0]  # Start with intercept
        for j in range(n_predictors):
            expected += coefficients[j+1] * X_pat[j]
        
        # Calculate W-score
        wscores[i] = (patient_data[i] - expected) / std_residuals
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save W-score map
    wscore_data = nib.gifti.gifti.GiftiDataArray(
        data=wscores.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    wscore_gii = nib.gifti.GiftiImage(darrays=[wscore_data])
    nib.save(wscore_gii, output_file)
    
    if verbose:
        print(f"    Saved W-score map: {output_file}")
    
    return {
        'wscore_file': output_file,
        'mean_wscore': np.mean(wscores),
        'std_wscore': np.std(wscores),
        'normative_data': normative_data,
        'reference_count': len(X_ref)
    }

def calculate_wscore_csv(reference_data, patient_data, demographics_ref, demographics_pat, output_file, normative_columns=['age', 'sex'], verbose=True):
    """
    Calculate W-scores for subcortical data using normative modeling and save as CSV file.
    
    Parameters:
    -----------
    reference_data : pd.DataFrame
        Combined reference data from all subjects
    patient_data : pd.DataFrame
        Patient data (single subject)
    demographics_ref : pd.DataFrame
        Demographics data for reference subjects
    demographics_pat : pd.Series or pd.DataFrame
        Demographics data for patient
    output_file : str
        Path to save the W-score CSV
    normative_columns : list, default=['age', 'sex']
        List of demographic columns to use for normative modeling
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing W-score statistics and file path
    """
    # Get structure columns (exclude metadata columns)
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)
    
    # Prepare demographics data
    demo_ref = demographics_ref[normative_columns].copy()
    demo_pat = demographics_pat[normative_columns].copy() if isinstance(demographics_pat, pd.DataFrame) else demographics_pat[normative_columns]
    
    # Convert to numpy arrays
    X_ref = demo_ref.values.astype(float)
    X_pat = demo_pat.values.astype(float) if isinstance(demo_pat, pd.DataFrame) else demo_pat.values.astype(float)
    if len(X_pat.shape) > 1:
        X_pat = X_pat.flatten()
    
    # Remove subjects with missing demographic data
    mask = ~np.isnan(X_ref).any(axis=1)
    X_ref = X_ref[mask]
    ref_data_clean = reference_data.loc[mask]
    
    if len(X_ref) == 0:
        raise ValueError("No reference subjects with complete demographic data")
    
    # Add intercept term
    X_ref_with_intercept = np.hstack([np.ones((X_ref.shape[0], 1)), X_ref])
    
    # Calculate W-scores for each structure
    patient_wscores = {}
    structure_stats = {}
    
    for structure in structures:
        if structure in reference_data.columns and structure in patient_data.columns:
            ref_values = ref_data_clean[structure].values
            pat_value = patient_data[structure].iloc[0]
            
            # Skip if all reference values are zero
            if np.all(ref_values == 0):
                patient_wscores[structure] = 0
                structure_stats[structure] = {
                    'intercept': 0, 'coefficients': [0] * len(normative_columns),
                    'std_residuals': 1, 'count': len(ref_values)
                }
                continue
            
            # Fit linear regression
            coefficients = np.linalg.lstsq(X_ref_with_intercept, ref_values, rcond=None)[0]
            
            # Calculate residuals and their standard deviation
            predicted = X_ref_with_intercept @ coefficients
            residuals = ref_values - predicted
            std_residuals = np.std(residuals)
            
            if std_residuals == 0:
                std_residuals = 1  # Avoid division by zero
            
            # Predict expected value for patient
            expected = coefficients[0] + np.sum([coefficients[j+1] * X_pat[j] for j in range(len(normative_columns))])
            
            # Calculate W-score
            wscore = (pat_value - expected) / std_residuals
            patient_wscores[structure] = wscore
            
            structure_stats[structure] = {
                'intercept': coefficients[0],
                'coefficients': coefficients[1:].tolist(),
                'std_residuals': std_residuals,
                'count': len(ref_values)
            }
                
    # Create output directory and save CSV
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    wscore_df = pd.DataFrame([patient_wscores])
    wscore_df.to_csv(output_file, index=False)
    
    if verbose:
        print(f"    Saved subcortical W-score CSV: {output_file}")
    
    return {
        'wscore_file': output_file,
        'wscores': patient_wscores,
        'structure_stats': structure_stats
    }

def calculate_zscore_maps(reference_data, patient_data, output_file, analysis='regional', verbose=True):
    """
    Calculate z-scores for patient data against reference data and save as GIFTI file.
    
    Parameters:
    -----------
    reference_data : np.ndarray
        Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
    patient_data : np.ndarray
        Patient data with shape (n_vertices,) or (n_vertices, n_depths)
    output_file : str
        Path to save the z-score map
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing z-score statistics and file path
    """
    if len(reference_data) == 0:
        raise ValueError("No reference data provided")
    
    if len(reference_data) == 1:
        raise ValueError("Only one subject in reference data, cannot calculate z-scores")
    
    # Calculate reference statistics
    ref_mean = np.mean(reference_data, axis=0)
    ref_std = np.std(reference_data, axis=0)
    
    # Avoid division by zero
    ref_std[ref_std == 0] = 1.0
    
    # Calculate z-score
    zscore = (patient_data - ref_mean) / ref_std
    
    # For multi-depth data (blur features), average across depths
    if len(zscore.shape) > 1:
        mean_zscore = np.mean(zscore, axis=1)  # Average across depths
        zscore_to_save = mean_zscore
    else:
        zscore_to_save = zscore
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save z-score map
    zscore_data = nib.gifti.gifti.GiftiDataArray(
        data=zscore_to_save.astype(np.float32),
        intent="NIFTI_INTENT_NORMAL"
    )
    zscore_gii = nib.gifti.GiftiImage(darrays=[zscore_data])
    nib.save(zscore_gii, output_file)
    
    if verbose:
        print(f"    Saved z-score map: {output_file}")
    
    return {
        'zscore_file': output_file,
        'mean_zscore': np.mean(zscore_to_save),
        'std_zscore': np.std(zscore_to_save),
        'reference_mean': ref_mean,
        'reference_std': ref_std,
        'reference_count': len(reference_data)
    }

def calculate_zscore_csv(reference_data, patient_data, output_file, verbose=True):
    """
    Calculate z-scores for subcortical data and save as CSV file.
    
    Parameters:
    -----------
    reference_data : pd.DataFrame
        Combined reference data from all subjects
    patient_data : pd.DataFrame
        Patient data (single subject)
    output_file : str
        Path to save the z-score CSV
    verbose : bool, default=True
        If True, prints processing information
    
    Returns:
    --------
    dict
        Dictionary containing z-score statistics and file path
    """
    # Get structure columns (exclude metadata columns)
    structures = reference_data.columns.tolist()
    for col in ['structure', 'SubjID']:
        if col in structures:
            structures.remove(col)
    
    # Calculate statistics for each structure
    structure_stats = {}
    for structure in structures:
        if structure in reference_data.columns:
            ref_values = reference_data[structure].values
            structure_stats[structure] = {
                'mean': np.mean(ref_values),
                'std': np.std(ref_values) if np.std(ref_values) > 0 else 1.0,
                'count': len(ref_values)
            }
    
    # Calculate z-scores for patient
    patient_zscores = {}
    for structure in structures:
        if structure in patient_data.columns and structure in structure_stats:
            pat_value = patient_data[structure].iloc[0]
            ref_mean = structure_stats[structure]['mean']
            ref_std = structure_stats[structure]['std']
            
            zscore = (pat_value - ref_mean) / ref_std
            patient_zscores[structure] = zscore
    
    # Create output directory and save CSV
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    zscore_df = pd.DataFrame([patient_zscores])
    zscore_df.to_csv(output_file, index=False)
    
    if verbose:
        print(f"    Saved subcortical z-score CSV: {output_file}")
    
    return {
        'zscore_file': output_file,
        'zscores': patient_zscores,
        'structure_stats': structure_stats
    }

def load_reference_surface_data(reference_subjects, output_directory, file_suffix, analysis='regional', verbose=True):
    """
    Load reference surface data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry maps
    verbose : bool, default=True
        If True, prints loading information
    
    Returns:
    --------
    tuple: (np.ndarray, list)
        - Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        ref_bids_id = f"{ref_pid}_{ref_sid}"
        
        if analysis == 'asymmetry':
            # For asymmetry, we need both left and right hemisphere data
            ref_file_lh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
            )
            ref_file_rh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
            )
            
            if os.path.exists(ref_file_lh) and os.path.exists(ref_file_rh):
                try:
                    ref_img_lh = nib.load(ref_file_lh)
                    ref_img_rh = nib.load(ref_file_rh)
                    
                    # Handle multi-depth data (blur features)
                    if len(ref_img_lh.darrays) > 1:
                        data_lh = np.zeros(shape=(ref_img_lh.darrays[0].data.shape[0], len(ref_img_lh.darrays)))
                        data_rh = np.zeros(shape=(ref_img_rh.darrays[0].data.shape[0], len(ref_img_rh.darrays)))
                        
                        for e, (darray_lh, darray_rh) in enumerate(zip(ref_img_lh.darrays, ref_img_rh.darrays)):
                            data_lh[:, e] = darray_lh.data
                            data_rh[:, e] = darray_rh.data
                        
                        asymmetry_data = compute_asymmetry(data_lh, data_rh)
                        reference_data.append(asymmetry_data)
                        successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                    else:
                        data_lh = ref_img_lh.darrays[0].data
                        data_rh = ref_img_rh.darrays[0].data
                        
                        asymmetry_data = compute_asymmetry(data_lh, data_rh)
                        reference_data.append(asymmetry_data)
                        successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                        
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference files {ref_file_lh} or {ref_file_rh}: {e}")
        else:
            # Regional analysis - load single hemisphere data
            ref_file = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "cortex",
                f"{ref_bids_id}_{file_suffix}"
            )
            
            if os.path.exists(ref_file):
                try:
                    ref_img = nib.load(ref_file)
                    
                    # Handle multi-depth data (blur features)
                    if len(ref_img.darrays) > 1:
                        darrays = np.zeros(shape=(ref_img.darrays[0].data.shape[0], len(ref_img.darrays)))
                        for e, darray in enumerate(ref_img.darrays):
                            darrays[:, e] = darray.data
                        reference_data.append(darrays)
                        successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                    else:
                        reference_data.append(ref_img.darrays[0].data)
                        successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                        
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed cortical reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return np.array(reference_data) if reference_data else np.array([]), successfully_loaded_subjects

def load_reference_hippocampal_data(reference_subjects, output_directory, file_suffix, analysis='regional', verbose=True):
    """
    Load reference hippocampal data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry maps
    verbose : bool, default=True
        If True, prints loading information
    
    Returns:
    --------
    tuple: (np.ndarray, list)
        - Array of reference data with shape (n_subjects, n_vertices) or (n_subjects, n_vertices, n_depths)
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
        For asymmetry analysis, returns left-right asymmetry maps
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        if analysis == 'asymmetry':
            # For asymmetry, we need both left and right hemisphere data
            ref_file_lh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
            )
            ref_file_rh = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
            )
            
            if os.path.exists(ref_file_lh) and os.path.exists(ref_file_rh):
                try:
                    ref_img_lh = nib.load(ref_file_lh)
                    ref_img_rh = nib.load(ref_file_rh)
                    
                    data_lh = ref_img_lh.darrays[0].data
                    data_rh = ref_img_rh.darrays[0].data
                    
                    # Compute asymmetry
                    asymmetry_data = compute_asymmetry(data_lh, data_rh)
                    reference_data.append(asymmetry_data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                    
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference files {ref_file_lh} or {ref_file_rh}: {e}")
        else:
            # Regional analysis - load single hemisphere data
            ref_file = os.path.join(
                output_directory,
                ref_pid, ref_sid, "maps", "hippocampus",
                f"{ref_pid}_{ref_sid}_{file_suffix}"
            )
            
            if os.path.exists(ref_file):
                try:
                    ref_img = nib.load(ref_file)
                    reference_data.append(ref_img.darrays[0].data)
                    successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                except Exception as e:
                    if verbose:
                        print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed hippocampal reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return np.array(reference_data) if reference_data else np.array([]), successfully_loaded_subjects

def load_reference_subcortical_data(reference_subjects, output_directory, file_suffix, analysis='regional', verbose=True):
    """
    Load reference subcortical data from multiple subjects.
    
    Parameters:
    -----------
    reference_subjects : list
        List of (participant_id, session_id) tuples
    output_directory : str
        Base output directory
    file_suffix : str
        File suffix pattern for the data files
    analysis : str, default='regional'
        Analysis type. If 'asymmetry', computes left-right asymmetry for subcortical structures
    verbose : bool, default=True
        If True, prints loading information
    
    Returns:
    --------
    tuple: (pd.DataFrame, list)
        - Combined reference data from all subjects
        - List of (participant_id, session_id) tuples for subjects that successfully loaded
        For asymmetry analysis, returns asymmetry values for paired structures
    """
    reference_data = []
    successfully_loaded_subjects = []
    failed_subjects = []  # Track failures
    
    for ref_pid, ref_sid in reference_subjects:
        ref_bids_id = f"{ref_pid}_{ref_sid}"
        ref_file = os.path.join(
            output_directory,
            ref_pid, ref_sid, "maps", "subcortical",
            f"{ref_bids_id}_{file_suffix}"
        )
        
        if os.path.exists(ref_file):
            try:
                ref_df = pd.read_csv(ref_file)
                
                if analysis == 'asymmetry':
                    # Compute asymmetry for paired structures
                    asymmetry_data = {}
                    for col in ref_df.columns:
                        if col.startswith('L'):
                            right_col = col.replace('L', 'R')
                            if right_col in ref_df.columns:
                                left_val = ref_df[col].iloc[0]
                                right_val = ref_df[right_col].iloc[0]
                                avg_val = (left_val + right_val) / 2
                                asym_val = (left_val - right_val) / avg_val if avg_val > 0 else 0
                                asymmetry_data[col] = asym_val
                    
                    # Convert to DataFrame
                    ref_df = pd.DataFrame([asymmetry_data])
                
                reference_data.append(ref_df)
                successfully_loaded_subjects.append((ref_pid, ref_sid))  # Track success
                
            except Exception as e:
                if verbose:
                    print(f"    Warning: Could not load reference file {ref_file}: {e}")
    
    if verbose and failed_subjects:
        print("    Failed subcortical reference subjects:")
        for pid, sid, reason in failed_subjects:
            print(f"      {pid}/{sid} -> {reason}")
    return (pd.concat(reference_data, ignore_index=True) if reference_data else pd.DataFrame(),
            successfully_loaded_subjects)

def analyze_dataset(dataset, reference, method='zscore', output_directory=None, verbose=True):
    """
    Analyze a dataset by comparing it to a reference dataset using specified method.
    
    Parameters:
    -----------
    dataset : zbdataset
        Dataset to analyze (patient dataset)
    reference : zbdataset
        Reference dataset to compare against (e.g., control dataset)
    method : str, default='zscore'
        Analysis method to use. Options: 'zscore', 'wscore'
    normative_columns : list, optional
        List of demographic columns to use for normative modeling (only used for wscore)
        Default: ['age', 'sex']
    output_directory : str, optional
        Directory where score maps will be saved. If None, uses validation output directory
    verbose : bool, default=True
        If True, prints detailed information about the analysis process
    
    Returns:
    --------
    dict
        Dictionary containing analysis results for each feature and region
    """
    if method not in ['zscore', 'wscore']:
        raise ValueError(f"Method '{method}' not supported. Currently supports 'zscore' and 'wscore'.")
    
    # Check that features match
    if set(dataset.features) != set(reference.features):
        raise ValueError("Feature sets don't match between patient and reference datasets.")
    
    # Check if demographics data is available for wscore
    if method == 'wscore':
        if not hasattr(dataset, 'demographics') or not hasattr(reference, 'demographics'):
            raise ValueError("Demographics data required for both datasets when using wscore method.")
        
        if not hasattr(dataset.demographics, 'normative_columns') or not hasattr(reference.demographics, 'normative_columns'):
            raise ValueError("Normative columns must be specified in dataset and reference demographics for wscore method.")
        
        # Check if required columns exist
        for col in dataset.demographics.normative_columns:
            if col not in dataset.demographics.data.columns:
                raise ValueError(f"Column '{col}' not found in dataset demographics.")
            if col not in reference.demographics.data.columns:
                raise ValueError(f"Column '{col}' not found in reference demographics.")
    
    # Use default output directory if not provided
    if output_directory is None:
        raise ValueError("output_directory must be specified for saving score maps.")
    
    if verbose:
        print(f"Analyzing dataset {dataset.name} against reference {reference.name} using {method} method...")
        print(f"Features to analyze: {', '.join(dataset.features)}")
        if method == 'wscore':
            print(f"Normative columns: {', '.join(dataset.demographics.normative_columns)}")
    
    # Get valid subjects from both datasets
    if not hasattr(dataset, 'valid_subjects'):
        raise ValueError("Dataset does not have valid_subjects attribute. Please validate the dataset first.")
    if not hasattr(reference, 'valid_subjects'):
        raise ValueError("Reference dataset does not have valid_subjects attribute. Please validate the reference dataset first.")
    
    # Print summary of valid subjects
    if verbose:
        print(f"Patient dataset: {len(dataset.valid_subjects['base'])} base subjects")
        print(f"Reference dataset: {len(reference.valid_subjects['base'])} base subjects")
        
        # Print feature-specific counts
        for feature in dataset.features:
            if feature in dataset.valid_subjects:
                pat_count = len(dataset.valid_subjects[feature]['all']) if 'all' in dataset.valid_subjects[feature] else 0
                ref_count = len(reference.valid_subjects[feature]['all']) if feature in reference.valid_subjects and 'all' in reference.valid_subjects[feature] else 0
                print(f"Feature {feature}: {pat_count} patient subjects, {ref_count} reference subjects")
    
    # Feature mapping for different naming conventions
    feature_mapping = {
        "thickness": {"output": "thickness"},
        "flair": {"output": "FLAIR"},
        "adc": {"output": "ADC"},
        "fa": {"output": "FA"},
        "qt1": {"output": "qT1"},  # Changed from "T1map" to "qT1"
        "fmri": {"output": "fMRI"}
    }
    
    # Define analysis parameters
    resolutions = ["32k", "5k"]
    labels = ["midthickness", "white"]
    hemispheres = ["L", "R"]
    
    # Store analysis results
    analysis_results = {
        "cortical": {},
        "hippocampal": {},
        "subcortical": {},
        "blur": {}
    }
    
    # Get demographics data for wscore
    if method == 'wscore':
        # Create mapping from (participant_id, session_id) to demographics
        ref_demo_dict = {}
        for _, row in reference.demographics.data.iterrows():
            key = (row['participant_id'], row.get('session_id', 'ses-01'))
            ref_demo_dict[key] = row
        
        pat_demo_dict = {}
        for _, row in dataset.demographics.data.iterrows():
            key = (row['participant_id'], row.get('session_id', 'ses-01'))
            pat_demo_dict[key] = row
    
    # Get normative columns for W-score calculation
    normative_columns = dataset.demographics.normative_columns if method == 'wscore' else None
    
    # Process each feature
    for feature in dataset.features:
        if verbose:
            print(f"\nAnalyzing feature: {feature}")
        
        feat_lower = feature.lower()
        is_blur = feat_lower.endswith("*blur")
        is_fmri = feat_lower == "fmri"
        # Get base feature name and mapping
        if is_blur:
            base_feat = feat_lower.replace("*blur", "")
            if base_feat in feature_mapping:
                output_feat = feature_mapping[base_feat]["output"] + "*blur"
            else:
                output_feat = base_feat + "*blur"
        else:
            if feat_lower in feature_mapping:
                output_feat = feature_mapping[feat_lower]["output"]
            else:
                output_feat = feat_lower
        
        # 1. CORTICAL ANALYSIS
        if dataset.cortex and reference.cortex and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing cortical data for {feature}...")
            
            # Get valid subjects for this feature and structure
            dataset_cortical_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_cortical_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_cortical_subjects or not reference_cortical_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for cortical {feature}")
                continue
                
            cortical_results = {}
            
            for hemi in hemispheres:
                for resolution in resolutions:
                    for label in labels:
                        for analysis in ['regional', 'asymmetry']:
                            map_key = f"{hemi}_{resolution}_{label}_{analysis}"
                            
                            file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"
                            
                            # Load reference data and track which subjects actually loaded
                            reference_data, successfully_loaded_subjects = load_reference_surface_data(
                                reference_cortical_subjects, output_directory, file_suffix, analysis, verbose
                            )
                            
                            if len(reference_data) == 0:
                                if verbose:
                                    print(f"    Warning: No reference data found for {map_key}")
                                continue
                            
                            # Get demographics for reference subjects (if using wscore)
                            if method == 'wscore':
                                ref_demographics = []
                                valid_ref_subjects = []
                                # Now iterate over subjects that ACTUALLY loaded, not all subjects
                                for ref_pid, ref_sid in successfully_loaded_subjects:
                                    key = (ref_pid, ref_sid)
                                    if key in ref_demo_dict:
                                        ref_demographics.append(ref_demo_dict[key])
                                        valid_ref_subjects.append(key)
                                
                                if len(ref_demographics) == 0:
                                    if verbose:
                                        print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                                    continue
                                
                                ref_demographics_df = pd.DataFrame(ref_demographics)
                                # Filter reference data to match demographics
                                # Now indices will match because we're iterating over successfully_loaded_subjects
                                ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                                reference_data = reference_data[ref_indices]
                            
                            # Process patient data
                            patient_scores = []
                            for pat_pid, pat_sid in dataset_cortical_subjects:
                                pat_bids_id = f"{pat_pid}_{pat_sid}"
                                
                                if analysis == 'asymmetry':
                                    # For asymmetry, we need both left and right hemisphere data
                                    pat_file_lh = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
                                    )
                                    pat_file_rh = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                                    )
                                    
                                    if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                        try:
                                            pat_img_lh = nib.load(pat_file_lh)
                                            pat_img_rh = nib.load(pat_file_rh)
                                            
                                            # Handle multi-depth data (blur features)
                                            if len(pat_img_lh.darrays) > 1:
                                                pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                                pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                                
                                                for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                                    pat_data_lh[:, e] = darray_lh.data
                                                    pat_data_rh[:, e] = darray_rh.data
                                                
                                                # Compute asymmetry
                                                pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                            else:
                                                pat_data_lh = pat_img_lh.darrays[0].data
                                                pat_data_rh = pat_img_rh.darrays[0].data
                                                
                                                # Compute asymmetry
                                                pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                            
                                            # Create score output file with analysis type in filename
                                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                            score_file = os.path.join(
                                                score_dir, 
                                                f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                            )
                                            
                                            # Calculate scores
                                            if method == 'zscore':
                                                score_result = calculate_zscore_maps(
                                                    reference_data, pat_data, score_file, analysis, verbose
                                                )
                                            else:  # wscore
                                                # Get patient demographics
                                                pat_key = (pat_pid, pat_sid)
                                                if pat_key not in pat_demo_dict:
                                                    if verbose:
                                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                    continue
                                                
                                                pat_demographics = pat_demo_dict[pat_key]
                                                score_result = calculate_wscore_maps(
                                                    reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                    score_file, normative_columns, verbose
                                                )
                                            
                                            score_result['subject'] = (pat_pid, pat_sid)
                                            score_result['analysis'] = analysis
                                            patient_scores.append(score_result)
                                            
                                        except Exception as e:
                                            if verbose:
                                                print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                                    else:
                                        if verbose:
                                            print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                                
                                else:
                                    # Regional analysis - load single hemisphere data
                                    pat_file = os.path.join(
                                        output_directory,
                                        pat_pid, pat_sid, "maps", "cortex",
                                        f"{pat_bids_id}_{file_suffix}"
                                    )
                                    
                                    if os.path.exists(pat_file):
                                        try:
                                            pat_img = nib.load(pat_file)
                                            pat_data = pat_img.darrays[0].data
                                            
                                            # Create score output file with analysis type in filename
                                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                            score_file = os.path.join(
                                                score_dir, 
                                                f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                            )
                                            
                                            # Calculate scores
                                            if method == 'zscore':
                                                score_result = calculate_zscore_maps(
                                                    reference_data, pat_data, score_file, analysis, verbose
                                                )
                                            else:  # wscore
                                                # Get patient demographics
                                                pat_key = (pat_pid, pat_sid)
                                                if pat_key not in pat_demo_dict:
                                                    if verbose:
                                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                    continue
                                                
                                                pat_demographics = pat_demo_dict[pat_key]
                                                score_result = calculate_wscore_maps(
                                                    reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                    score_file, normative_columns, verbose
                                                )
                                            
                                            score_result['subject'] = (pat_pid, pat_sid)
                                            score_result['analysis'] = analysis
                                            patient_scores.append(score_result)
                                            
                                        except Exception as e:
                                            if verbose:
                                                print(f"    Warning: Could not process patient file {pat_file}: {e}")
                            
                            cortical_results[map_key] = {
                                f'patient_{method}s': patient_scores
                            }

            analysis_results["cortical"][feature] = cortical_results
        
        # 2. HIPPOCAMPAL ANALYSIS
        if dataset.hippocampus and reference.hippocampus and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing hippocampal data for {feature}...")
            
            # Determine resolution string from dataset
            hipp_res = "8k" if getattr(dataset, "hippunfold_version", 1) == 2 else "0p5mm"
            
            # Get valid subjects for this feature and structure
            dataset_hippocampal_subjects = dataset.valid_subjects[feature]['structures']['hippocampus'] if feature in dataset.valid_subjects else []
            reference_hippocampal_subjects = reference.valid_subjects[feature]['structures']['hippocampus'] if feature in reference.valid_subjects else []
            
            if not dataset_hippocampal_subjects or not reference_hippocampal_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for hippocampal {feature}")
                continue
                
            hippocampal_results = {}
            
            for hemi in hemispheres:
                for analysis in ['regional', 'asymmetry']:
                    map_key = f"{hemi}_hippocampus_{analysis}"  # Include analysis in map_key
                    file_suffix = f"hemi-{hemi}_den-{hipp_res}_label-hipp_midthickness_feature-{output_feat}_smooth-{dataset.hippocampal_smoothing}mm.func.gii"
                    
                    # Load reference data
                    reference_data, successfully_loaded_subjects = load_reference_hippocampal_data(
                        reference_hippocampal_subjects, output_directory, file_suffix, analysis, verbose
                    )
                    
                    if len(reference_data) == 0:
                        if verbose:
                            print(f"    Warning: No reference data found for {map_key}")
                        continue
                    
                    # Get demographics for reference subjects (if using wscore)
                    if method == 'wscore':
                        ref_demographics = []
                        valid_ref_subjects = []
                        for ref_pid, ref_sid in successfully_loaded_subjects:
                            key = (ref_pid, ref_sid)
                            if key in ref_demo_dict:
                                ref_demographics.append(ref_demo_dict[key])
                                valid_ref_subjects.append(key)
                        
                        if len(ref_demographics) == 0:
                            if verbose:
                                print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                            continue
                        
                        ref_demographics_df = pd.DataFrame(ref_demographics)
                        # Filter reference data to match demographics
                        ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                        reference_data = reference_data[ref_indices]
                    
                    # Process patient data
                    patient_scores = []
                    for pat_pid, pat_sid in dataset_hippocampal_subjects:
                        if analysis == 'asymmetry':
                            # For asymmetry, we need both left and right hemisphere data
                            pat_file_lh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"

                            )
                            pat_file_rh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                            )
                            
                            if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                try:
                                    pat_img_lh = nib.load(pat_file_lh)
                                    pat_img_rh = nib.load(pat_file_rh)
                                    
                                    pat_data_lh = pat_img_lh.darrays[0].data
                                    pat_data_rh = pat_img_rh.darrays[0].data
                                    
                                    # Compute asymmetry
                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                    
                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "hippocampus")
                                    score_file = os.path.join(
                                        score_dir, 
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                            else:
                                if verbose:
                                    print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                        
                        else:
                            # Regional analysis - load single hemisphere data
                            pat_file = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "hippocampus",
                                f"{pat_pid}_{pat_sid}_{file_suffix}"
                            )
                            
                            if os.path.exists(pat_file):
                                try:
                                    pat_img = nib.load(pat_file)
                                    pat_data = pat_img.darrays[0].data
                                    
                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "hippocampus")
                                    score_file = os.path.join(
                                        score_dir, 
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient file {pat_file}: {e}")
                    
                    hippocampal_results[map_key] = {
                        f'patient_{method}s': patient_scores
                    }
            
            analysis_results["hippocampal"][feature] = hippocampal_results
        
        # 3. BLUR ANALYSIS
        if is_blur:
            if verbose:
                print(f"  Processing blur depth data for {feature}...")
            
            # For blur features, use cortical subjects from the base feature
            base_feature = feat_lower.replace("*blur", "")
            dataset_blur_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_blur_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_blur_subjects or not reference_blur_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for blur {feature}")
                continue
                
            blur_results = {}
            
            for hemi in hemispheres:
                for analysis in ['regional', 'asymmetry']:
                    map_key = f"{hemi}_blur_depths_{analysis}"  # Include analysis in map_key
                    label = "midthickness"
                    resolution = "32k"
                    
                    file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"
                    
                    # Load reference data
                    reference_data, successfully_loaded_subjects = load_reference_surface_data(
                        reference_blur_subjects, output_directory, file_suffix, analysis, verbose
                    )
                    
                    if len(reference_data) == 0:
                        if verbose:
                            print(f"    Warning: No reference data found for {map_key}")
                        continue
                    
                    if len(reference_data) == 1:
                        if verbose:
                            print(f"    Warning: only one subject in reference data for {map_key}, skipping")
                        continue
                    
                    # Get demographics for reference subjects (if using wscore)
                    if method == 'wscore':
                        ref_demographics = []
                        valid_ref_subjects = []
                        for ref_pid, ref_sid in successfully_loaded_subjects:
                            key = (ref_pid, ref_sid)
                            if key in ref_demo_dict:
                                ref_demographics.append(ref_demo_dict[key])
                                valid_ref_subjects.append(key)
                        
                        if len(ref_demographics) == 0:
                            if verbose:
                                print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                            continue
                        
                        ref_demographics_df = pd.DataFrame(ref_demographics)
                        # Filter reference data to match demographics
                        ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                        reference_data = reference_data[ref_indices]
                    
                    # Process patient data
                    patient_scores = []
                    for pat_pid, pat_sid in dataset_blur_subjects:
                        if analysis == 'asymmetry':
                            # For asymmetry, we need both left and right hemisphere data
                            pat_file_lh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"

                            )
                            pat_file_rh = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                            )
                            
                            if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                try:
                                    pat_img_lh = nib.load(pat_file_lh)
                                    pat_img_rh = nib.load(pat_file_rh)
                                    
                                    # Handle multi-depth data
                                    pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                    pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                    
                                    for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                        pat_data_lh[:, e] = darray_lh.data
                                        pat_data_rh[:, e] = darray_rh.data
                                    
                                    # Compute asymmetry
                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                    
                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                    score_file = os.path.join(
                                        score_dir,
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores (will automatically average across depths)
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    score_result[f'avg_{method}_file'] = score_file
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                            else:
                                if verbose:
                                    print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                        
                        else:
                            # Regional analysis - load single hemisphere data
                            pat_file = os.path.join(
                                output_directory,
                                pat_pid, pat_sid, "maps", "cortex",
                                f"{pat_pid}_{pat_sid}_{file_suffix}"
                            )
                            
                            if os.path.exists(pat_file):
                                try:
                                    pat_img = nib.load(pat_file)
                                    pat_data = np.zeros(shape=(pat_img.darrays[0].data.shape[0], len(pat_img.darrays)))
                                    for e, darray in enumerate(pat_img.darrays):
                                        pat_data[:, e] = darray.data
                                    
                                    # Create score output file with analysis type in filename
                                    score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                    score_file = os.path.join(
                                        score_dir,
                                        f"{pat_pid}_{pat_sid}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                    )
                                    
                                    # Calculate scores (will automatically average across depths)
                                    if method == 'zscore':
                                        score_result = calculate_zscore_maps(
                                            reference_data, pat_data, score_file, analysis, verbose
                                        )
                                    else:  # wscore
                                        # Get patient demographics
                                        pat_key = (pat_pid, pat_sid)
                                        if pat_key not in pat_demo_dict:
                                            if verbose:
                                                print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                            continue
                                        
                                        pat_demographics = pat_demo_dict[pat_key]
                                        score_result = calculate_wscore_maps(
                                            reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                            score_file, normative_columns, verbose
                                        )
                                    
                                    score_result['subject'] = (pat_pid, pat_sid)
                                    score_result['analysis'] = analysis
                                    score_result[f'avg_{method}_file'] = score_file
                                    patient_scores.append(score_result)
                                    
                                except Exception as e:
                                    if verbose:
                                        print(f"    Warning: Could not process patient file {pat_file}: {e}")
                    
                    blur_results[map_key] = {
                        f'patient_{method}s': patient_scores
                    }
            
            analysis_results["blur"][feature] = blur_results
        
        # 4. SUBCORTICAL ANALYSIS
        if dataset.subcortical and reference.subcortical and not is_blur and not is_fmri:
            if verbose:
                print(f"  Processing subcortical data for {feature}...")
            
            # Get valid subjects for this feature and structure
            dataset_subcortical_subjects = dataset.valid_subjects[feature]['structures']['subcortical'] if feature in dataset.valid_subjects else []
            reference_subcortical_subjects = reference.valid_subjects[feature]['structures']['subcortical'] if feature in reference.valid_subjects else []
            
            if not dataset_subcortical_subjects or not reference_subcortical_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for subcortical {feature}")
                continue
            
            subcortical_results = {}
            
            for analysis in ['regional', 'asymmetry']:
                map_key = f"subcortical_{analysis}"  # Include analysis in map_key
                
                if feat_lower in ["thickness", "volume"]:
                    file_suffix = "feature-volume.csv"
                else:
                    file_suffix = f"feature-{output_feat}.csv"
                
                # Load reference data
                reference_data, successfully_loaded_subjects = load_reference_subcortical_data(
                    reference_subcortical_subjects, output_directory, file_suffix, analysis, verbose
                )
                
                if reference_data.empty:
                    if verbose:
                        print(f"    Warning: No reference data found for {map_key}")
                    continue
                
                # Get demographics for reference subjects (if using wscore)
                if method == 'wscore':
                    ref_demographics = []
                    valid_ref_indices = []
                    for i, (ref_pid, ref_sid) in enumerate(successfully_loaded_subjects):
                        key = (ref_pid, ref_sid)
                        if key in ref_demo_dict:
                            ref_demographics.append(ref_demo_dict[key])
                            valid_ref_indices.append(i)
                    
                    if len(ref_demographics) == 0:
                        if verbose:
                            print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                        continue
                    
                    ref_demographics_df = pd.DataFrame(ref_demographics)
                    # Filter reference data to match demographics
                    reference_data = reference_data.iloc[valid_ref_indices]
                
                # Process patient data
                patient_scores = []
                for pat_pid, pat_sid in dataset_subcortical_subjects:
                    pat_bids_id = f"{pat_pid}_{pat_sid}"
                    pat_file = os.path.join(
                        output_directory,
                        pat_pid, pat_sid, "maps", "subcortical",
                        f"{pat_bids_id}_{file_suffix}"
                    )
                    
                    if os.path.exists(pat_file):
                        try:
                            pat_df = pd.read_csv(pat_file)
                            
                            # For asymmetry analysis, compute left-right differences
                            if analysis == 'asymmetry':
                                # Create asymmetry data by computing differences between left and right structures
                                asymmetry_data = {}
                                for col in pat_df.columns:
                                    if col.startswith('L'):
                                        right_col = col.replace('L', 'R')
                                        if right_col in pat_df.columns:
                                            left_val = pat_df[col].iloc[0]
                                            right_val = pat_df[right_col].iloc[0]
                                            # Compute asymmetry: (L - R) / ((L + R) / 2)
                                            avg_val = (left_val + right_val) / 2
                                            asym_val = (left_val - right_val) / avg_val if avg_val > 0 else 0
                                            asymmetry_data[col] = asym_val
                                
                                # Convert to DataFrame
                                pat_df = pd.DataFrame([asymmetry_data])
                            
                            # Create score output file with analysis type in filename
                            score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "subcortical")
                            score_file = os.path.join(
                                score_dir, 
                                f"{pat_bids_id}_{file_suffix.replace('.csv', f'_analysis-{analysis}.csv')}"

                            )
                            
                            # Calculate scores
                            if method == 'zscore':
                                score_result = calculate_zscore_csv(
                                    reference_data, pat_df, score_file, verbose
                                )
                            else:  # wscore
                                # Get patient demographics
                                pat_key = (pat_pid, pat_sid)
                                if pat_key not in pat_demo_dict:
                                    if verbose:
                                        print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                    continue
                                
                                pat_demographics = pat_demo_dict[pat_key]
                                score_result = calculate_wscore_csv(
                                    reference_data, pat_df, ref_demographics_df, pat_demographics, 
                                    score_file, normative_columns, verbose
                                )
                            
                            score_result['subject'] = (pat_pid, pat_sid)
                            score_result['analysis'] = analysis
                            patient_scores.append(score_result)
                            
                        except Exception as e:
                            if verbose:
                                print(f"    Warning: Could not process patient file {pat_file}: {e}")
                
                subcortical_results[map_key] = {
                    f'patient_{method}s': patient_scores
                }
            
            analysis_results["subcortical"][feature] = subcortical_results
    
        # 5. fMRI ANALYSIS
        if is_fmri:
            if verbose:
                print(f"  Processing fMRI data for {feature}...")
            
            # Get valid subjects for this feature
            dataset_fmri_subjects = dataset.valid_subjects[feature]['structures']['cortex'] if feature in dataset.valid_subjects else []
            reference_fmri_subjects = reference.valid_subjects[feature]['structures']['cortex'] if feature in reference.valid_subjects else []
            
            if not dataset_fmri_subjects or not reference_fmri_subjects:
                if verbose:
                    print(f"    Warning: No valid subjects found for fMRI {feature}")
                continue
            
            fmri_features = ['rmssd', 'timescales', 'alff', 'falff']
            
            for fmri_feat in fmri_features:
                if verbose:
                    print(f"    Processing fMRI feature: {fmri_feat}")
                
                fmri_results = {}
                
                for hemi in hemispheres:
                    for resolution in resolutions:
                        for label in ['midthickness']:
                            for analysis in ['regional', 'asymmetry']:
                                map_key = f"{hemi}_{resolution}_{label}_{analysis}"
                                
                                file_suffix = f"hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{fmri_feat}_smooth-{dataset.cortical_smoothing}mm.func.gii"
                                
                                # Load reference data and track which subjects actually loaded
                                reference_data, successfully_loaded_subjects = load_reference_surface_data(
                                    reference_fmri_subjects, output_directory, file_suffix, analysis, verbose
                                )
                                
                                if len(reference_data) == 0:
                                    if verbose:
                                        print(f"    Warning: No reference data found for {map_key}")
                                    continue
                                
                                # Get demographics for reference subjects (if using wscore)
                                if method == 'wscore':
                                    ref_demographics = []
                                    valid_ref_subjects = []
                                    # Now iterate over subjects that ACTUALLY loaded, not all subjects
                                    for ref_pid, ref_sid in successfully_loaded_subjects:
                                        key = (ref_pid, ref_sid)
                                        if key in ref_demo_dict:
                                            ref_demographics.append(ref_demo_dict[key])
                                            valid_ref_subjects.append(key)
                                    
                                    if len(ref_demographics) == 0:
                                        if verbose:
                                            print(f"    Warning: No demographics data found for reference subjects in {map_key}")
                                        continue
                                    
                                    ref_demographics_df = pd.DataFrame(ref_demographics)
                                    # Filter reference data to match demographics
                                    # Now indices will match because we're iterating over successfully_loaded_subjects
                                    ref_indices = [i for i, subj in enumerate(successfully_loaded_subjects) if subj in valid_ref_subjects]
                                    reference_data = reference_data[ref_indices]
                                
                                # Process patient data
                                patient_scores = []
                                for pat_pid, pat_sid in dataset_fmri_subjects:
                                    pat_bids_id = f"{pat_pid}_{pat_sid}"
                                    
                                    if analysis == 'asymmetry':
                                        # For asymmetry, we need both left and right hemisphere data
                                        pat_file_lh = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-L').replace('hemi-R', 'hemi-L')}"
                                        )
                                        pat_file_rh = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix.replace('hemi-L', 'hemi-R').replace('hemi-R', 'hemi-R')}"
                                        )
                                        
                                        if os.path.exists(pat_file_lh) and os.path.exists(pat_file_rh):
                                            try:
                                                pat_img_lh = nib.load(pat_file_lh)
                                                pat_img_rh = nib.load(pat_file_rh)
                                                
                                                # Handle multi-depth data (blur features)
                                                if len(pat_img_lh.darrays) > 1:
                                                    pat_data_lh = np.zeros(shape=(pat_img_lh.darrays[0].data.shape[0], len(pat_img_lh.darrays)))
                                                    pat_data_rh = np.zeros(shape=(pat_img_rh.darrays[0].data.shape[0], len(pat_img_rh.darrays)))
                                                    
                                                    for e, (darray_lh, darray_rh) in enumerate(zip(pat_img_lh.darrays, pat_img_rh.darrays)):
                                                        pat_data_lh[:, e] = darray_lh.data
                                                        pat_data_rh[:, e] = darray_rh.data
                                                    
                                                    # Compute asymmetry
                                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                                else:
                                                    pat_data_lh = pat_img_lh.darrays[0].data
                                                    pat_data_rh = pat_img_rh.darrays[0].data
                                                    
                                                    # Compute asymmetry
                                                    pat_data = compute_asymmetry(pat_data_lh, pat_data_rh)
                                                
                                                # Create score output file with analysis type in filename
                                                score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                                score_file = os.path.join(
                                                    score_dir, 
                                                    f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                                )
                                                
                                                # Calculate scores
                                                if method == 'zscore':
                                                    score_result = calculate_zscore_maps(
                                                        reference_data, pat_data, score_file, analysis, verbose
                                                    )
                                                else:  # wscore
                                                    # Get patient demographics
                                                    pat_key = (pat_pid, pat_sid)
                                                    if pat_key not in pat_demo_dict:
                                                        if verbose:
                                                            print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                        continue
                                                    
                                                    pat_demographics = pat_demo_dict[pat_key]
                                                    score_result = calculate_wscore_maps(
                                                        reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                        score_file, normative_columns, verbose
                                                    )
                                                
                                                score_result['subject'] = (pat_pid, pat_sid)
                                                score_result['analysis'] = analysis
                                                patient_scores.append(score_result)
                                                
                                            except Exception as e:
                                                if verbose:
                                                    print(f"    Warning: Could not process patient files {pat_file_lh} or {pat_file_rh}: {e}")
                                        else:
                                            if verbose:
                                                print(f"    Warning: Missing files for asymmetry analysis: {pat_file_lh} or {pat_file_rh}")
                                    
                                    else:
                                        # Regional analysis - load single hemisphere data
                                        pat_file = os.path.join(
                                            output_directory,
                                            pat_pid, pat_sid, "maps", "cortex",
                                            f"{pat_bids_id}_{file_suffix}"
                                        )
                                        
                                        if os.path.exists(pat_file):
                                            try:
                                                pat_img = nib.load(pat_file)
                                                pat_data = pat_img.darrays[0].data
                                                
                                                # Create score output file with analysis type in filename
                                                score_dir = os.path.join(output_directory, pat_pid, pat_sid, f"{method}_maps", "cortex")
                                                score_file = os.path.join(
                                                    score_dir, 
                                                    f"{pat_bids_id}_{file_suffix.replace('.func.gii', f'_analysis-{analysis}.func.gii')}"
                                                )
                                                
                                                # Calculate scores
                                                if method == 'zscore':
                                                    score_result = calculate_zscore_maps(
                                                        reference_data, pat_data, score_file, analysis, verbose
                                                    )
                                                else:  # wscore
                                                    # Get patient demographics
                                                    pat_key = (pat_pid, pat_sid)
                                                    if pat_key not in pat_demo_dict:
                                                        if verbose:
                                                            print(f"    Warning: No demographics data found for patient {pat_pid}/{pat_sid}")
                                                        continue
                                                    
                                                    pat_demographics = pat_demo_dict[pat_key]
                                                    score_result = calculate_wscore_maps(
                                                        reference_data, pat_data, ref_demographics_df, pat_demographics, 
                                                        score_file, normative_columns, verbose
                                                    )
                                                
                                                score_result['subject'] = (pat_pid, pat_sid)
                                                score_result['analysis'] = analysis
                                                patient_scores.append(score_result)
                                                
                                            except Exception as e:
                                                if verbose:
                                                    print(f"    Warning: Could not process patient file {pat_file}: {e}")
                                
                                fmri_results[map_key] = {
                                    f'patient_{method}s': patient_scores
                                }
                
                analysis_results["cortical"][fmri_feat] = fmri_results

    if verbose:
        print(f"\nAnalysis complete! {method.upper()} maps saved to {method}_maps directories")
        
        # Print summary statistics
        for region_type, region_results in analysis_results.items():
            if region_results:
                print(f"\n{region_type.capitalize()} analysis:")
                for feature, feature_results in region_results.items():
                    if region_type == "subcortical":
                        if f'patient_{method}s' in feature_results:
                            print(f"  {feature}: {len(feature_results[f'patient_{method}s'])} subjects analyzed")
                    else:
                        map_count = sum(len(maps[f'patient_{method}s']) for maps in feature_results.values())
                        print(f"  {feature}: {map_count} {method} maps generated")
    
    return analysis_results

def compute_asymmetry(data_lh, data_rh):
    """Compute asymmetry.

    Parameters
    ----------
    x_lh
        Left hemisphere data. Shape (n_subjects, n_points) or (n_points,).
    x_rh
        Right hemisphere data. Shape (n_subjects, n_points) or (n_points,).

    Returns
    -------
    s
        Output data. Shape (n_subjects, n_points) or (n_points,).
    """

    den = data_lh + data_rh
    den *= 0.5
    return np.divide(data_lh - data_rh, den, out=den, where=den > 0)