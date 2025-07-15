import os
import sys
import datetime
from src.processing import apply_blurring, apply_hippocampal_processing, apply_subcortical_processing, apply_cortical_processing
from src.analysis import analyze_dataset
import shutil
import multiprocessing
from joblib import Parallel, delayed

class LogRedirect:
    """
    Class for redirecting stdout/stderr to a log file while still displaying output.
    """
    def __init__(self, log_file):
        self.log_file = log_file
        self.terminal_stdout = sys.stdout
        
        # Remove old log file if it exists
        if os.path.exists(log_file):
            os.remove(log_file)
            
        # Open in write mode to create a new file
        self.log = open(self.log_file, "w", encoding="utf-8")
        
    def __del__(self):
        if hasattr(self, 'log') and self.log:
            self.log.close()
        
    def write(self, message):
        self.terminal_stdout.write(message)
        self.log.write(message)
        self.log.flush()
        
    def flush(self):
        self.terminal_stdout.flush()
        self.log.flush()

class demographics():
    def __init__(self, csv_file, column_mapping=None, normative_columns=None, normative_dtypes=None, reference=None):
        self.csv_file = csv_file
        self.data = None
        self.column_mapping = {
            "ID": "participant_id",
            "SES": "session_id",
        } if column_mapping is None else column_mapping
        self.normative_columns = normative_columns
        self.normative_dtypes = normative_dtypes
        self.reference = reference
        self.binary_encodings = {}  # Store binary variable encodings
        self.load_data()

    def __repr__(self):
        return f"Demographics(csv_file={self.csv_file})"
    def __str__(self):
        return f"Demographics data loaded from {self.csv_file}"
    
    def load_data(self):
        import pandas as pd
        import numpy as np
        try:
            self.data = pd.read_csv(self.csv_file)
        except Exception as e:
            raise ValueError(f"Failed to load demographics data: {e}")
        if self.data.empty:
            raise ValueError("Demographics data is empty.")
        # Rename columns based on mapping
        self.data.rename(columns=self.column_mapping, inplace=True)
        print(f"Demographics data loaded from {self.csv_file} with {len(self.data)} entries.")
        if not set(self.column_mapping.values()).issubset(self.data.columns):
            raise ValueError("Demographics data does not contain all required columns.")
        
        # Use reference normative columns if available and none are specified
        single_subject = len(self.data) == 1
        using_reference = False
        
        if self.reference is not None:
            if self.normative_columns is None and hasattr(self.reference, 'normative_columns'):
                self.normative_columns = self.reference.normative_columns
                self.normative_dtypes = self.reference.normative_dtypes
                using_reference = True
                print(f"Using normative columns from reference: {self.normative_columns}")
        
        # Validate normative columns and dtypes
        normative_columns = self.normative_columns if self.normative_columns else []
        
        # Check if normative_dtypes is provided when normative_columns is specified
        if normative_columns and self.normative_dtypes is None:
            raise ValueError("normative_dtypes must be provided when normative_columns is specified")
            
        # Check if normative_dtypes has the same length as normative_columns
        if normative_columns and len(normative_columns) != len(self.normative_dtypes):
            raise ValueError(f"normative_dtypes ({len(self.normative_dtypes)}) must have the same length as normative_columns ({len(normative_columns)})")
            
        # Validate each normative column
        for col in normative_columns:
            if col not in self.data.columns:
                raise ValueError(f"Normative column '{col}' not found in demographics data.")
            
            # Check for missing values
            if self.data[col].isnull().any():
                raise ValueError(f"Normative column '{col}' contains missing values.")
        
        # Validate data types if both normative_columns and normative_dtypes are specified
        if normative_columns and self.normative_dtypes:
            for i, (col, dtype) in enumerate(zip(normative_columns, self.normative_dtypes)):
                try:
                    # Validate column values based on specified data type
                    if dtype.lower() == 'int':
                        # Try converting to int and check if values are preserved
                        original_values = self.data[col].values
                        int_values = self.data[col].astype(int).values
                        # Check if conversion to int maintains the original values
                        if not np.allclose(original_values, int_values, rtol=1e-05, atol=1e-08, equal_nan=True):
                            problematic_values = self.data.loc[~np.isclose(original_values, int_values)][col].unique()
                            raise ValueError(f"Column '{col}' has non-integer values: {problematic_values}")
                        # Convert to int
                        self.data[col] = self.data[col].astype(int)
                    
                    elif dtype.lower() == 'float':
                        # Try converting to float
                        self.data[col] = self.data[col].astype(float)
                    
                    elif dtype.lower() == 'binary':
                        # Check if we should use reference encoding
                        if single_subject and self.reference is not None and hasattr(self.reference, 'binary_encodings') and col in self.reference.binary_encodings:
                            # Get encoding from reference
                            ref_value_to_binary = self.reference.binary_encodings[col]
                            
                            # Convert using reference encoding
                            self.data[col] = self.data[col].astype(str).str.lower().map(
                                lambda x: ref_value_to_binary.get(x.lower(), 0)
                            )
                            print(f"Applied reference encoding for binary column '{col}'")
                            
                        else:
                            # Check if the column has exactly 2 unique values
                            unique_values = self.data[col].astype(str).str.lower().unique()
                            if len(unique_values) != 2:
                                raise ValueError(f"Column '{col}' has {len(unique_values)} unique values, expected exactly 2 for binary data type")
                            
                            # Convert to numeric binary values (0/1) - arbitrary but consistent mapping
                            first_value = unique_values[0]
                            value_to_binary = {val.lower(): 1 if val.lower() == first_value.lower() else 0 for val in unique_values}
                            
                            # Store the encoding for future reference
                            self.binary_encodings[col] = value_to_binary
                            
                            # Apply encoding
                            self.data[col] = self.data[col].astype(str).str.lower().map(value_to_binary)
                    
                    else:
                        raise ValueError(f"Unsupported dtype '{dtype}' for column '{col}'. "
                                        f"Supported types: 'int', 'float', 'binary'")
                        
                except Exception as e:
                    raise ValueError(f"Error validating column '{col}' with dtype '{dtype}': {e}")
            
            print(f"Normative columns: {normative_columns} with types {self.normative_dtypes} successfully validated.")
            if using_reference and single_subject:
                print("Used reference encoding for binary variables to ensure consistency in normative modeling.")
        else:
            print(f"Normative columns: {normative_columns} successfully validated.")
    
        return self

class zbdataset():
    def __init__(self, name, demographics : demographics, micapipe_directory, hippunfold_directory=None, freesurfer_directory=None, cortex=True, hippocampus=True, subcortical=True):
        
        self.name = name
        self.demographics = demographics
        self.micapipe_directory = micapipe_directory
        self.hippunfold_directory = hippunfold_directory
        self.freesurfer_directory = freesurfer_directory
        self.features = []
        self.cortex = cortex
        self.hippocampus = hippocampus
        self.subcortical = subcortical
        self.valid_dataset = False
       
    def __repr__(self):
        return f"Dataset(name={self.name})"

    def __str__(self):
        return f"Dataset: {self.name}"
    
    def check_directories(self):
        if self.hippunfold_directory is None:
            print("Warning: hippunfold_directory is not specified, hippocampal data will not be available.")
        if self.freesurfer_directory is None:
            print("Warning: freesurfer_directory is not specified, subcortical data will not be available.")
        print("Checking agreement between demographics data and micapipe directory...")

        if not os.path.exists(self.micapipe_directory):
            raise ValueError(f"Micapipe directory {self.micapipe_directory} does not exist.")
        if self.hippunfold_directory and not os.path.exists(self.hippunfold_directory):
            raise ValueError(f"Hippunfold directory {self.hippunfold_directory} does not exist.")
        if self.freesurfer_directory and not os.path.exists(self.freesurfer_directory):
            raise ValueError(f"Freesurfer directory {self.freesurfer_directory} does not exist.")
        print("All directories exist.")

        # Initialize valid_subjects structure
        self.valid_subjects = {
            'base': [],
            'structures': {
                'cortex': [],
                'hippocampus': [],
                'subcortical': []
            }
        }
        
        # Check if demographics data matches micapipe directory
        missing_micapipe = []
        missing_hippunfold = []
        missing_freesurfer = []
        
        for _, row in self.demographics.data.iterrows():
            participant_id = row['participant_id']
            session_id = row['session_id']
            
            # Track if this participant/session has all required directories
            is_valid = True
            
            # Check for directory in micapipe
            participant_session_path = os.path.join(self.micapipe_directory, f"{participant_id}", f"{session_id}")
            if not os.path.exists(participant_session_path):
                missing_micapipe.append((participant_id, session_id))
                is_valid = False
                
            # Check for directory in hippunfold if specified
            if self.hippunfold_directory:
                hippunfold_path = os.path.join(self.hippunfold_directory, "hippunfold", f"{participant_id}", f"{session_id}")
                if not os.path.exists(hippunfold_path):
                    missing_hippunfold.append((participant_id, session_id))
                    is_valid = False
                    
            # Check for directory in freesurfer if specified
            if self.freesurfer_directory:
                freesurfer_path = os.path.join(self.freesurfer_directory, f"{participant_id}_{session_id}")
                if not os.path.exists(freesurfer_path):
                    missing_freesurfer.append((participant_id, session_id))
                    is_valid = False
            
            # If all required directories exist, add to valid subjects
            if is_valid:
                self.valid_subjects['base'].append((participant_id, session_id))
        
        # Report missing directories
        if missing_micapipe:
            print(f"Warning: {len(missing_micapipe)} participant/session directories not found in micapipe directory:")
            for p, s in missing_micapipe:
                print(f"  - {p}/{s}")
                
        if missing_hippunfold:
            print(f"Warning: {len(missing_hippunfold)} participant/session directories not found in hippunfold directory:")
            for p, s in missing_hippunfold:
                print(f"  - {p}/{s}")
                
        if missing_freesurfer:
            print(f"Warning: {len(missing_freesurfer)} participant/session directories not found in freesurfer directory:")
            for p, s in missing_freesurfer:
                print(f"  - {p}/{s}")
        
        print(f"Found {len(self.valid_subjects['base'])} valid subjects with complete directory structure.")
        return self
    
    def add_features(self, *features, verbose=True):
        """
        Add features to the dataset.
        Features should be specified as strings, e.g., "FA", "ADC", etc.
        
        This function checks for the existence of required files in the micapipe directory
        for each feature and each valid subject.
        
        Parameters:
        -----------
        *features : str
            Names of features to add
        verbose : bool, default=True
            If True, prints detailed information about missing files and features.
            If False, only prints the summary with percentages.
        """
        
        self.check_directories()

        self.features = features
        if verbose:
            print(f"Features {self.features} added to the dataset {self.name}.")
        
        # Define file patterns for each feature
        feature_files = {
            "FA": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_FA.func.gii",
                   "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-FA.nii.gz",],
            "ADC": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_ADC.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-ADC.nii.gz"],
            "thickness": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-thickness.func.gii",],
            "FLAIR": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_flair.func.gii",
                      "maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_T1map.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
            "FLAIR-blur": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_flair.func.gii",
                           "maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1-blur": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_T1map.func.gii",
                         "maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
        }

        # Define required files for different structures
        required_cortical_files = [
            "anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz",
            "anat/{participant_id}_{session_id}_space-nativepro_T1w_brain_mask.nii.gz",
            "surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-white.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-white.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
        ]
        
        required_hippocampal_files = [
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
        ]

        required_subcortical_files = [
            "{participant_id}_{session_id}/stats/aseg.stats",
        ]

        # Track missing files and feature availability
        missing_files = {}
        missing_features = {}
        feature_availability = {feature: 0 for feature in features}
        subjects_with_complete_data = []
        
        # Initialize feature-specific valid subjects with nested structure
        for feature in features:
            self.valid_subjects[feature] = {
                'all': [],  # All structures
                'structures': {
                    'cortex': [],
                    'hippocampus': [],
                    'subcortical': []
                }
            }
        
        # Check for each valid subject from the base list
        for participant_id, session_id in self.valid_subjects['base']:
            subject_missing_files = []
            subject_missing_features = set()
            has_required_cortical = True
            has_required_hippocampal = True
            has_required_subcortical = True
            subject = (participant_id, session_id)
            
            # Check required cortical files
            for file_pattern in required_cortical_files:
                file_path = os.path.join(
                    self.micapipe_directory,
                    participant_id,
                    session_id,
                    file_pattern.format(participant_id=participant_id, session_id=session_id)
                )
                if not os.path.exists(file_path):
                    subject_missing_files.append(file_path)
                    has_required_cortical = False
            
            # Check hippocampal files if enabled
            if self.hippocampus and self.hippunfold_directory:
                for file_pattern in required_hippocampal_files:
                    file_path = os.path.join(
                        self.hippunfold_directory,
                        "hippunfold",
                        participant_id,
                        session_id,
                        file_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    if not os.path.exists(file_path):
                        subject_missing_files.append(file_path)
                        has_required_hippocampal = False
            
            # Check subcortical files if enabled
            if self.subcortical and self.freesurfer_directory:
                for file_pattern in required_subcortical_files:
                    file_path = os.path.join(
                        self.freesurfer_directory,
                        file_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    if not os.path.exists(file_path):
                        subject_missing_files.append(file_path)
                        has_required_subcortical = False
            
            # Update structure-specific lists
            if has_required_cortical:
                self.valid_subjects['structures']['cortex'].append(subject)
            
            if has_required_hippocampal:
                self.valid_subjects['structures']['hippocampus'].append(subject)
            
            if has_required_subcortical:
                self.valid_subjects['structures']['subcortical'].append(subject)
            
            # Check feature-specific files and structure availability
            for feature in features:
                feature_complete = True
                # Keep track of which structures have this feature
                feature_has_cortical = has_required_cortical
                feature_has_hippocampal = has_required_hippocampal
                feature_has_subcortical = has_required_subcortical
                
                if feature not in feature_files:
                    print(f"Warning: Unknown feature '{feature}'. No file checks performed.")
                    continue
                
                for file_pattern in feature_files[feature]:
                    # Handle file patterns with hemisphere placeholders
                    if "{hemi}" in file_pattern:
                        for hemi in ["L", "R"]:
                            # Handle file patterns with surface type placeholders
                            if "{surfacetype}" in file_pattern:
                                for surfacetype in ["midthickness", "white"]:
                                    file_path = os.path.join(
                                        self.micapipe_directory,
                                        participant_id,
                                        session_id,
                                        file_pattern.format(
                                            participant_id=participant_id, 
                                            session_id=session_id,
                                            hemi=hemi,
                                            surfacetype=surfacetype
                                        )
                                    )
                                    if not os.path.exists(file_path):
                                        subject_missing_files.append(file_path)
                                        feature_complete = False
                                        feature_has_cortical = False  # Surface files affect cortical processing
                            else:
                                # Pattern has hemi but not surfacetype
                                file_path = os.path.join(
                                    self.micapipe_directory,
                                    participant_id,
                                    session_id,
                                    file_pattern.format(
                                        participant_id=participant_id, 
                                        session_id=session_id,
                                        hemi=hemi
                                    )
                                )
                                if not os.path.exists(file_path):
                                    subject_missing_files.append(file_path)
                                    feature_complete = False
                                    feature_has_cortical = False
                    else:
                        # Pattern has neither hemi nor surfacetype (volume file)
                        file_path = os.path.join(
                            self.micapipe_directory,
                            participant_id,
                            session_id,
                            file_pattern.format(
                                participant_id=participant_id, 
                                session_id=session_id
                            )
                        )
                        if not os.path.exists(file_path):
                            subject_missing_files.append(file_path)
                            feature_complete = False
                            # Volume files affect all structures
                            feature_has_cortical = False
                            feature_has_hippocampal = False
                            feature_has_subcortical = False
                
                # Update feature availability
                if feature_complete:
                    feature_availability[feature] += 1
                    self.valid_subjects[feature]['all'].append(subject)
                    
                    # Update structure-specific availability for this feature
                    if feature_has_cortical:
                        self.valid_subjects[feature]['structures']['cortex'].append(subject)
                    
                    if feature_has_hippocampal and not feature.lower().endswith('-blur'):
                        # Blur features don't apply to hippocampus
                        self.valid_subjects[feature]['structures']['hippocampus'].append(subject)
                    
                    if feature_has_subcortical and not feature.lower().endswith('-blur'):
                        # Blur features don't apply to subcortical
                        self.valid_subjects[feature]['structures']['subcortical'].append(subject)
                else:
                    subject_missing_features.add(feature)
            
            # Track missing features for this subject
            if subject_missing_features:
                missing_features[subject] = subject_missing_features
                
            # Track missing files and complete subjects
            if subject_missing_files:
                missing_files[subject] = subject_missing_files
            else:
                subjects_with_complete_data.append(subject)
        
        # Store the results
        self.subjects_with_complete_data = subjects_with_complete_data
        self.missing_files = missing_files
        self.missing_features = missing_features
        self.feature_availability = feature_availability
        
        # Report findings
        if verbose:
            if subjects_with_complete_data:
                print(f"Found {len(subjects_with_complete_data)} subjects with complete data for all requested features.")
            
            if len(self.valid_subjects['structures']['cortex']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['cortex'])} subjects are missing required cortical files.")
            
            if self.hippocampus and len(self.valid_subjects['structures']['hippocampus']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['hippocampus'])} subjects are missing required hippocampal files.")
                
            if self.subcortical and len(self.valid_subjects['structures']['subcortical']) < len(self.valid_subjects['base']):
                print(f"Warning: {len(self.valid_subjects['base']) - len(self.valid_subjects['structures']['subcortical'])} subjects are missing required subcortical files.")
        
        # Always print feature availability summary
        print("\nFeature availability summary:")
        total_subjects = len(self.valid_subjects['base'])
        for feature in features:
            if feature in feature_availability:
                avail_count = feature_availability[feature]
                print(f"  {feature}: {avail_count}/{total_subjects} subjects ({avail_count/total_subjects*100:.1f}%)")
                
                # Print structure-specific availability for this feature
                if verbose:
                    cortex_count = len(self.valid_subjects[feature]['structures']['cortex'])
                    print(f"    - cortex: {cortex_count}/{total_subjects} subjects ({cortex_count/total_subjects*100:.1f}%)")
                    
                    if self.hippocampus and not feature.lower().endswith('-blur'):
                        hippo_count = len(self.valid_subjects[feature]['structures']['hippocampus'])
                        print(f"    - hippocampus: {hippo_count}/{total_subjects} subjects ({hippo_count/total_subjects*100:.1f}%)")
                    
                    if self.subcortical and not feature.lower().endswith('-blur'):
                        subcort_count = len(self.valid_subjects[feature]['structures']['subcortical'])
                        print(f"    - subcortical: {subcort_count}/{total_subjects} subjects ({subcort_count/total_subjects*100:.1f}%)")
        
        # Print structure availability summary
        print("\nStructure availability summary:")
        print(f"  cortex: {len(self.valid_subjects['structures']['cortex'])}/{total_subjects} subjects ({len(self.valid_subjects['structures']['cortex'])/total_subjects*100:.1f}%)")
        if self.hippocampus:
            print(f"  hippocampus: {len(self.valid_subjects['structures']['hippocampus'])}/{total_subjects} subjects ({len(self.valid_subjects['structures']['hippocampus'])/total_subjects*100:.1f}%)")
        if self.subcortical:
            print(f"  subcortical: {len(self.valid_subjects['structures']['subcortical'])}/{total_subjects} subjects ({len(self.valid_subjects['structures']['subcortical'])/total_subjects*100:.1f}%)")
        
        # Only print detailed missing features if verbose
        if verbose and missing_features:
            print(f"\nMissing features by subject:")
            for (pid, sid), missing_feats in missing_features.items():
                print(f"  {pid}/{sid}: missing {', '.join(missing_feats)}")
        
        return self

    def process(self, output_directory, features, cortical_smoothing=5, hippocampal_smoothing=2, env=None, verbose=True, n_jobs=None):
        """
        Process the dataset with specified features and smoothing parameters using joblib parallelization.
        Logs all output to a log file in each subject's session directory.
        
        Parameters:
        -----------
        output_directory : str
            Directory to store processed data
        features : list or None, default=None
            List of features to process. If None, uses previously added features.
        cortical_smoothing : int, default=5
            Smoothing parameter for cortical features
        hippocampal_smoothing : int, default=2
            Smoothing parameter for hippocampal features
        env : object
            Environment object containing paths to required tools
        verbose : bool, default=True
            If True, prints detailed processing information
        n_jobs : int, optional
            Number of parallel jobs to run. If None, uses all available CPU cores.
            If 1, runs sequentially.
            
        Returns:
        --------
        self
            The dataset object
        """
        self.cortical_smoothing = cortical_smoothing
        self.hippocampal_smoothing = hippocampal_smoothing
        self.add_features(*features, verbose=verbose)
        
        if verbose:
            print(f"Processing dataset {self.name} with cortical smoothing {cortical_smoothing} and hippocampal smoothing {hippocampal_smoothing}.")
        
        if env is None:
            raise ValueError("Environment (zbenv) must be provided for processing.")

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
            print(f"Created output directory: {output_directory}")
        else:
            print(f"Output directory already exists: {output_directory}")

        # Create main log file
        timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        main_log_file = os.path.join(output_directory, f"processing_{self.name}_{timestamp}.log")
        with open(main_log_file, 'w', encoding='utf-8') as f:
            f.write(f"===== Processing dataset {self.name} =====\n")
            f.write(f"Date/Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output directory: {output_directory}\n")
            f.write(f"Features: {features}\n")
            f.write(f"Cortical smoothing: {cortical_smoothing}\n")
            f.write(f"Hippocampal smoothing: {hippocampal_smoothing}\n")
        print(f"Main processing log started at: {main_log_file}")

        # Determine number of jobs
        if n_jobs is None:
            n_jobs = env.num_threads or multiprocessing.cpu_count()
        
        # Use all base valid subjects
        valid_subjects_to_process = self.valid_subjects['base']
        
        if verbose:
            print(f"Using {n_jobs if n_jobs > 1 else 'sequential'} processing for {len(valid_subjects_to_process)} subjects")

        # Identify blur features
        blur_features = [feature for feature in self.features if feature.endswith("-blur")]
        base_features = [feature.replace("-blur", "") for feature in blur_features]
        
        # Define a function to process a single subject
        def process_single_subject(subject):
            participant_id, session_id = subject
            
            # Create session-specific directories
            session_output_dir = os.path.join(output_directory, participant_id, session_id)
            os.makedirs(session_output_dir, exist_ok=True)
            
            # Create session log file
            log_file = os.path.join(session_output_dir, f"{participant_id}_{session_id}_processing.log")
            
            # Create session-specific tmp directory
            session_tmp_dir = os.path.join(session_output_dir, f"tmp_{participant_id}_{session_id}")
            os.makedirs(session_tmp_dir, exist_ok=True)
            
            # Redirect stdout to log file for this subject
            original_stdout = sys.stdout
            log_redirect = None
            try:
                # Create log file
                log_redirect = LogRedirect(log_file)
                sys.stdout = log_redirect
                
                print(f"Processing subject {participant_id}/{session_id}...")
                print(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"Features: {features}")
                print(f"Cortical smoothing: {cortical_smoothing}mm")
                print(f"Hippocampal smoothing: {hippocampal_smoothing}mm")
                print(f"Log file: {log_file}")
                print("-" * 50)
                
                try:
                    # Copy structural files
                    self._copy_structural_files(participant_id, session_id, output_directory, verbose=verbose)
                    
                    # Apply blurring to features that need it
                    if blur_features:
                        # Check which base features are available for this subject
                        available_features = []
                        for feature in base_features:
                            if subject in self.valid_subjects[feature]['all']:
                                available_features.append(feature)
                        
                        if available_features:  # If ANY features are available, proceed with blurring
                            if verbose:
                                print(f"  Applying additional blur processing for features: {', '.join([f+'-blur' for f in available_features])}")
                            
                            apply_blurring(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=available_features,  # Only process available features
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                freesurfer_directory=self.freesurfer_directory,
                                tmp_dir=session_tmp_dir,
                                smoothing_fwhm=self.cortical_smoothing,
                                verbose=verbose
                            )

                    # Process cortical features if cortex is enabled
                    if self.cortex and subject in self.valid_subjects['structures']['cortex']:
                        # Get valid features for cortex for this subject
                        valid_cortical_features = [f for f in self.features if subject in self.valid_subjects[f]['structures']['cortex']]
                        
                        if valid_cortical_features:
                            if verbose:
                                print(f"  Processing cortical data for features: {', '.join(valid_cortical_features)}")
                            
                            apply_cortical_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_cortical_features,
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                tmp_dir=session_tmp_dir,
                                cortical_smoothing=cortical_smoothing,
                                resolutions=["32k", "5k"],
                                labels=["midthickness", "white"],
                                verbose=verbose
                            )
                
                    # If hippocampus is enabled, process hippocampal data
                    if self.hippocampus and self.hippunfold_directory and subject in self.valid_subjects['structures']['hippocampus']:
                        # Get non-blur features for hippocampus
                        valid_hipp_features = [f for f in self.features 
                                            if not f.endswith("-blur") 
                                            and subject in self.valid_subjects[f]['structures']['hippocampus']]
                        
                        if valid_hipp_features:
                            if verbose:
                                print(f"  Processing hippocampal data for features: {', '.join(valid_hipp_features)}")
                            
                            apply_hippocampal_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_hipp_features,
                                output_directory=output_directory,
                                workbench_path=env.connectome_workbench_path,
                                micapipe_directory=self.micapipe_directory,
                                hippunfold_directory=self.hippunfold_directory,
                                tmp_dir=session_tmp_dir,
                                smoothing_fwhm=hippocampal_smoothing,
                                verbose=verbose
                            )

                    # If subcortex is enabled, extract subcortical stats
                    if self.subcortical and self.freesurfer_directory and subject in self.valid_subjects['structures']['subcortical']:
                        # Get non-blur features for subcortical
                        valid_subcort_features = [f for f in self.features 
                                                if not f.endswith("-blur") 
                                                and subject in self.valid_subjects[f]['structures']['subcortical']]
                        
                        if valid_subcort_features:
                            if verbose:
                                print(f"  Processing subcortical data for features: {', '.join(valid_subcort_features)}")
                            
                            apply_subcortical_processing(
                                participant_id=participant_id,
                                session_id=session_id,
                                features=valid_subcort_features,
                                output_directory=output_directory,
                                micapipe_directory=self.micapipe_directory,
                                freesurfer_directory=self.freesurfer_directory,
                                verbose=verbose
                            )
                    
                    print(f"Completed processing {participant_id}/{session_id} at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                    return (participant_id, session_id, True)
                
                except Exception as e:
                    if verbose:
                        print(f"Error processing {participant_id}/{session_id}: {e}")
                    import traceback
                    traceback.print_exc()
                    return (participant_id, session_id, False)
                
                finally:
                    # Clean up session-specific tmp directory
                    if os.path.exists(session_tmp_dir):
                        shutil.rmtree(session_tmp_dir)
                        if verbose:
                            print(f"  Cleaned up temporary directory: {session_tmp_dir}")
                    
                    print("-" * 50)
                    print(f"Processing finished at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            finally:
                sys.stdout = original_stdout
                # Ensure log is closed
                if log_redirect and hasattr(log_redirect, 'log'):
                    log_redirect.log.close()
            
            return (participant_id, session_id, True)
        
        # Process subjects using joblib (or sequentially if n_jobs=1)
        if n_jobs == 1:
            print(f"Running sequential processing for {len(valid_subjects_to_process)} subjects")
            results = [process_single_subject(subject) for subject in valid_subjects_to_process]
        else:
            print(f"Running parallel processing with {n_jobs} jobs for {len(valid_subjects_to_process)} subjects")
            results = Parallel(n_jobs=n_jobs, verbose=10 if verbose else 0)(
                delayed(process_single_subject)(subject) for subject in valid_subjects_to_process
            )
        
        # Process results
        failed_subjects = [(pid, sid) for pid, sid, success in results if not success]
        
        # Report results
        if verbose:
            successful_count = len(valid_subjects_to_process) - len(failed_subjects)
            print(f"\nProcessing complete: {successful_count}/{len(valid_subjects_to_process)} subjects successful")
            
            if failed_subjects:
                print(f"Failed subjects ({len(failed_subjects)}):")
                for participant_id, session_id in failed_subjects:
                    print(f"  - {participant_id}/{session_id}")

        # Append summary to main log
        with open(main_log_file, 'a', encoding='utf-8') as f:
            f.write("\n===== PROCESSING SUMMARY =====\n")
            f.write(f"Completed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Subjects processed successfully: {len(valid_subjects_to_process) - len(failed_subjects)}/{len(valid_subjects_to_process)}\n")
            if failed_subjects:
                f.write(f"Failed subjects ({len(failed_subjects)}):\n")
                for participant_id, session_id in failed_subjects:
                    f.write(f"  - {participant_id}/{session_id}\n")

        print(f"Dataset {self.name} processed successfully.")
        print(f"Processing summary written to: {main_log_file}")
        
        self.validate(
            features=self.features,
            output_directory=output_directory,
            cortical_smoothing=cortical_smoothing,
            hippocampal_smoothing=hippocampal_smoothing,
            verbose=verbose
        )
        return self


    def validate(self, output_directory, features, cortical_smoothing=5, hippocampal_smoothing=2, verbose=True, error_threshold=0, warning_threshold=25):
        """
        Validate that all expected processed files exist for each subject and feature.
        
        Parameters:
        -----------
        output_directory : str
            Directory where processed data should be stored
        features : list or None, default=None
            List of features to validate. If None, uses previously added features.
        cortical_smoothing : int, default=5
            Smoothing parameter used for cortical features
        hippocampal_smoothing : int, default=2
            Smoothing parameter used for hippocampal features
        verbose : bool, default=True
            If True, prints detailed information about missing files
        error_threshold : int, default=0
            Percentage of complete subjects below which to raise an error (0 means error only if no subjects)
        warning_threshold : int, default=25
            Percentage of complete subjects below which to raise a warning

        Returns:
        --------
        dict
            A dictionary with validation results containing:
            - 'valid_subjects': List of subjects with complete data
            - 'missing_files': Dictionary of missing files per subject
            - 'summary': Overall validation summary
        
        Raises:
        -------
        ValueError
            When there are no valid subjects or percentage is below error_threshold
        Warning
            When the percentage of valid subjects is below warning_threshold
        """
        # Make sure we have the input features
        self.add_features(*features, verbose=verbose)
        
        # Check for output directories instead of input directories
        output_subjects = self.check_output_directories(output_directory, verbose=verbose)
        
        if not output_subjects:
            raise ValueError(f"No subject output directories found in {output_directory}. Nothing to validate.")
        
        if verbose:
            print(f"Validating dataset {self.name} with {len(output_subjects)} subjects and {len(self.features)} features...")
        
        # Identify blur features
        blur_features = [feature for feature in features if feature.endswith("-blur")]
        base_features = [feature.replace("-blur", "") for feature in blur_features]
        
        # Feature mapping for different naming conventions
        feature_mapping = {
            "thickness": {"output": "thickness"},
            "flair": {"output": "flair"},
            "adc": {"output": "ADC"},
            "fa": {"output": "FA"},
            "qt1": {"output": "qT1"}  # Changed from "T1map" to "qT1"
        }
        
        # Track missing files
        missing_files = {}
        valid_subjects = []
        all_files_count = 0
        missing_files_count = 0
        
        # Define surface resolutions and labels
        resolutions = ["32k", "5k"]
        labels = ["midthickness", "white"]
        
        # Process each subject with output directory
        for subject in output_subjects:
            participant_id, session_id = subject
            subject_missing_files = []
            bids_id = f"{participant_id}_{session_id}"
            
            # Define output directories
            subject_output_dir = os.path.join(output_directory, participant_id, session_id)
            cortex_dir = os.path.join(subject_output_dir, "maps", "cortex")
            hippocampus_dir = os.path.join(subject_output_dir, "maps", "hippocampus")
            subcortical_dir = os.path.join(subject_output_dir, "maps", "subcortical")
            
            # 2. Check cortical feature files
            if self.cortex:
                for feature in self.features:
                    feat_lower = feature.lower()
                    is_blur = feat_lower.endswith("-blur")
                    
                    # Get base feature name and mapping
                    if is_blur:
                        base_feat = feat_lower.replace("-blur", "")
                        if base_feat in feature_mapping:
                            output_feat = feature_mapping[base_feat]["output"] + "-blur"
                        else:
                            output_feat = base_feat + "-blur"
                    else:
                        if feat_lower in feature_mapping:
                            output_feat = feature_mapping[feat_lower]["output"]
                        else:
                            output_feat = feat_lower
                    
                    for hemi in ["L", "R"]:
                        for resolution in resolutions:
                            # For blur features, only check midthickness
                            label_list = ["midthickness"] if is_blur else labels
                            for label in label_list:
                                # Check cortical feature file
                                cortical_file = os.path.join(
                                    cortex_dir,
                                    f"{bids_id}_hemi-{hemi}_surf-fsLR-{resolution}_label-{label}_feature-{output_feat}_smooth-{cortical_smoothing}mm.func.gii"
                                )
                                all_files_count += 1
                                if not os.path.exists(cortical_file):
                                    subject_missing_files.append(cortical_file)
                                    missing_files_count += 1
            
            # 3. Check hippocampal feature files
            if self.hippocampus and self.hippunfold_directory:
                for feature in self.features:
                    if not feature.lower().endswith("-blur"):  # Blur features are not processed for hippocampus
                        feat_lower = feature.lower()
                        
                        # Get output feature name based on mapping
                        if feat_lower in feature_mapping:
                            output_feat = feature_mapping[feat_lower]["output"]
                        else:
                            output_feat = feat_lower
                        
                        for hemi in ["L", "R"]:
                            # Check hippocampal feature file
                            hippo_file = os.path.join(
                                hippocampus_dir,
                                f"{participant_id}_{session_id}_hemi-{hemi}_den-0p5mm_label-hipp_midthickness_feature-{output_feat}_smooth-{hippocampal_smoothing}mm.func.gii"
                        )
                        all_files_count += 1
                        if not os.path.exists(hippo_file):
                            subject_missing_files.append(hippo_file)
                            missing_files_count += 1
            
            # 4. Check subcortical feature files
            if self.subcortical and self.freesurfer_directory:
                # Check volume data
                volume_file = os.path.join(subcortical_dir, f"{bids_id}_feature-volume.csv")
                all_files_count += 1
                if not os.path.exists(volume_file):
                    subject_missing_files.append(volume_file)
                    missing_files_count += 1
                
                # Check other feature data
                for feature in self.features:
                    if not feature.lower().endswith("-blur") and feature.lower() not in ["thickness", "volume"]:
                        feat_lower = feature.lower()
                        
                        # Get output feature name based on mapping
                        if feat_lower in feature_mapping:
                            output_feat = feature_mapping[feat_lower]["output"]
                        else:
                            output_feat = feat_lower
                        
                        subcort_file = os.path.join(subcortical_dir, f"{bids_id}_feature-{output_feat}.csv")
                        all_files_count += 1
                        if not os.path.exists(subcort_file):
                            subject_missing_files.append(subcort_file)
                            missing_files_count += 1
            
            # Track missing files and valid subjects
            if subject_missing_files:
                missing_files[(participant_id, session_id)] = subject_missing_files
            else:
                valid_subjects.append((participant_id, session_id))
        
        # Calculate statistics
        total_subjects = len(self.valid_subjects['base'])
        complete_subjects = len(valid_subjects)
        complete_percentage = (complete_subjects / total_subjects) * 100 if total_subjects > 0 else 0
        file_presence_percentage = ((all_files_count - missing_files_count) / all_files_count) * 100 if all_files_count > 0 else 0
        
        # Create summary
        summary = {
            "total_subjects": total_subjects,
            "complete_subjects": complete_subjects,
            "complete_percentage": complete_percentage,
            "total_files_expected": all_files_count,
            "missing_files_count": missing_files_count,
            "file_presence_percentage": file_presence_percentage
        }
        
        # Print summary
        if verbose:
            print(f"\nValidation Summary:")
            print(f"  Total subjects: {total_subjects}")
            print(f"  Subjects with complete data: {complete_subjects} ({complete_percentage:.1f}%)")
            print(f"  Total expected files: {all_files_count}")
            print(f"  Missing files: {missing_files_count} ({100-file_presence_percentage:.1f}%)")
            
            if missing_files:
                print(f"\nSubjects with missing files: {len(missing_files)}")
                if verbose > 1:  # Extra verbosity for detailed missing file listing
                    for (pid, sid), files in missing_files.items():
                        print(f"  {pid}/{sid}: {len(files)} missing files")
                        if verbose > 2:  # Extremely verbose, list all missing files
                            for file in files:
                                print(f"    - {os.path.basename(file)}")
    
        # Error handling based on valid subjects
        if complete_subjects == 0:
            self.valid_dataset = False
            raise ValueError(f"Validation failed: No subjects have complete data. Check processing logs for details.")
        
        # Error if below error threshold (but not zero, which is handled above)
        if error_threshold > 0 and complete_percentage < error_threshold:
            self.valid_dataset = False
            raise ValueError(f"Validation failed: Only {complete_percentage:.1f}% of subjects have complete data (below error threshold of {error_threshold}%).")
        
        # Warning if below warning threshold
        if complete_percentage < warning_threshold:
            import warnings
            warnings.warn(f"Validation warning: Only {complete_percentage:.1f}% of subjects have complete data (below warning threshold of {warning_threshold}%).")

        self.valid_dataset = True
        self.cortical_smoothing = cortical_smoothing
        self.hippocampal_smoothing = hippocampal_smoothing
        return {
            "valid_subjects": valid_subjects,
            "missing_files": missing_files,
            "summary": summary
        }
    
    def check_output_directories(self, output_directory, verbose=True):
    
        """
        Check the output directory structure and identify valid subjects with output data.
        This is used during validation to determine which subjects to validate.
        
        Parameters:
        -----------
        output_directory : str
            Directory where processed data should be stored
        verbose : bool, default=True
            If True, prints detailed information about the directories

        Returns:
        --------
        list
            List of tuples containing (participant_id, session_id) for subjects with output directories
        """
        if not os.path.exists(output_directory):
            if verbose:
                print(f"Output directory {output_directory} does not exist. No subjects to validate.")
            return []
        
        # Find all subject directories in the output directory
        subject_dirs = [d for d in os.listdir(output_directory) 
                    if os.path.isdir(os.path.join(output_directory, d)) and d.startswith('sub-')]
        
        if verbose:
            print(f"Found {len(subject_dirs)} subject directories in output directory.")
        
        valid_output_subjects = []
        missing_cortex = []
        missing_hippocampus = []
        missing_subcortical = []
        
        # Check each subject directory for session directories
        for subject_dir in subject_dirs:
            participant_id = subject_dir
            subject_path = os.path.join(output_directory, subject_dir)
            
            # Find session directories
            session_dirs = [d for d in os.listdir(subject_path) 
                        if os.path.isdir(os.path.join(subject_path, d)) and d.startswith('ses-')]
            
            for session_dir in session_dirs:
                session_id = session_dir
                session_path = os.path.join(subject_path, session_dir)
                
                # Check for required output directories
                has_cortex = os.path.exists(os.path.join(session_path, "maps", "cortex"))
                has_hippocampus = os.path.exists(os.path.join(session_path, "maps", "hippocampus"))
                has_subcortical = os.path.exists(os.path.join(session_path, "maps", "subcortical"))
                
                # Determine if this is a valid subject based on enabled components
                is_valid = True
                
                if self.cortex and not has_cortex:
                    missing_cortex.append((participant_id, session_id))
                    is_valid = False
                    
                if self.hippocampus and not has_hippocampus:
                    missing_hippocampus.append((participant_id, session_id))
                    is_valid = False
                    
                if self.subcortical and not has_subcortical:
                    missing_subcortical.append((participant_id, session_id))
                    is_valid = False
                
                # If all required directories exist, add to valid subjects
                if is_valid:
                    valid_output_subjects.append((participant_id, session_id))
        
        # Report missing directories
        if verbose:
            if missing_cortex and self.cortex:
                print(f"Warning: {len(missing_cortex)} subjects missing cortex output directory:")
                for p, s in missing_cortex[:5]:  # Show only first 5 to avoid excessive output
                    print(f"  - {p}/{s}")
                if len(missing_cortex) > 5:
                    print(f"  - ... and {len(missing_cortex) - 5} more")
                    
            if missing_hippocampus and self.hippocampus:
                print(f"Warning: {len(missing_hippocampus)} subjects missing hippocampus output directory:")
                for p, s in missing_hippocampus[:5]:
                    print(f"  - {p}/{s}")
                if len(missing_hippocampus) > 5:
                    print(f"  - ... and {len(missing_hippocampus) - 5} more")
                    
            if missing_subcortical and self.subcortical:
                print(f"Warning: {len(missing_subcortical)} subjects missing subcortical output directory:")
                for p, s in missing_subcortical[:5]:
                    print(f"  - {p}/{s}")
                if len(missing_subcortical) > 5:
                    print(f"  - ... and {len(missing_subcortical) - 5} more")

        return valid_output_subjects
    def analyze(self, reference, method='zscore', output_directory=None, verbose=True):
        """
        Analyze the dataset by comparing it to a reference dataset using specified method.
        
        Parameters:
        -----------
        reference : zbdataset
            Reference dataset to compare against (e.g., control dataset)
        method : str, default='zscore'
            Analysis method to use. Currently supports 'zscore'
        output_directory : str, optional
            Directory where z-score maps will be saved. If None, uses validation output directory
        verbose : bool, default=True
            If True, prints detailed information about the analysis process
        
        Returns:
        --------
        dict
            Dictionary containing analysis results for each feature and region
        """
        from src.analysis import analyze_dataset

        # Call the analysis function and store results
        results = analyze_dataset(self, reference, method, output_directory, verbose)
        self.analysis_results = results
        
        return results
    
    def clinical_report(self, output_directory=None, approach='wscore', analyses=['regional','asymmetry'], features=None, 
                    threshold=1.96, threshold_alpha=0.3, color_range=(-3, 3), 
                    cmap='cmo.balance', cmap_asymmetry='cmo.balance_r', 
                    color_bar='bottom', tmp_dir=None, verbose=True):
        """
        Generate clinical reports for each subject in the dataset.
        
        Parameters:
        -----------
        output_directory : str, optional
            Directory where processed data is stored. If None, uses the directory from analysis
        approach : str, default='wscore'
            Analysis approach used ('zscore' or 'wscore')
        analyses : list, optional
            List of analyses to include in report. If None, uses ['regional']
        features : list, optional
            List of features to include in report. If None, uses dataset features
        threshold : float, default=1.96
            Threshold for significance in maps
        threshold_alpha : float, default=0.3
            Alpha transparency for thresholded regions
        color_range : tuple, default=(-3, 3)
            Color range for visualization
        cmap : str, default='cmo.balance'
            Colormap for regional analysis
        cmap_asymmetry : str, default='cmo.balance_r'
            Colormap for asymmetry analysis
        color_bar : str, default='bottom'
            Position of color bar
        tmp_dir : str, optional
            Base temporary directory. If None, uses session directory
        verbose : bool, default=True
            If True, prints detailed information
            
        Returns:
        --------
        list
            List of generated PDF report file paths
        """
        from src.clinical_reports import generate_clinical_report
        
        # Check if analysis has been run
        if not hasattr(self, 'analysis_results'):
            raise ValueError("No analysis results found. Please run dataset.analyze() first.")
            
        # Use features from the dataset if not specified
        if features is None:
            features = list(self.features)
            
        # Default analyses
        if analyses is None:
            analyses = ['regional']  # Default to regional analysis
            
        # Validate output directory
        if output_directory is None:
            raise ValueError("output_directory must be specified")
            
        # Don't create a separate reports directory - save directly in subject folders
        generated_reports = []
        
        # Extract subjects from analysis results
        valid_subjects = []
        for region_type, region_results in self.analysis_results.items():
            for feature, feature_results in region_results.items():
                if region_type == "subcortical":
                    if f'patient_{approach}s' in feature_results:
                        for result in feature_results[f'patient_{approach}s']:
                            if 'subject' in result:
                                valid_subjects.append(result['subject'])
                else:
                    for map_key, map_results in feature_results.items():
                        if f'patient_{approach}s' in map_results:
                            for result in map_results[f'patient_{approach}s']:
                                if 'subject' in result:
                                    valid_subjects.append(result['subject'])
        
        # Remove duplicates while preserving order
        valid_subjects = list(dict.fromkeys(valid_subjects))
        
        if not valid_subjects:
            if verbose:
                print("No valid subjects found in analysis results. Using all subjects from dataset.")
            valid_subjects = self.valid_subjects['base']
        
        if verbose:
            print(f"Generating clinical reports for {len(valid_subjects)} subjects...")
            
        # Generate report for each subject
        for participant_id, session_id in valid_subjects:
            try:
                # Get demographics for this subject
                subject_demo = self.demographics.data[
                    (self.demographics.data['participant_id'] == participant_id) &
                    (self.demographics.data['session_id'] == session_id)
                ]
                
                if subject_demo.empty:
                    if verbose:
                        print(f"Warning: No demographics found for {participant_id}/{session_id}, skipping...")
                    continue
                    
                # Extract demographics
                age = subject_demo['AGE'].iloc[0] if 'AGE' in subject_demo.columns else None
                sex = subject_demo['SEX'].iloc[0] if 'SEX' in subject_demo.columns else None
                
                # Convert binary sex encoding back to string if needed
                if sex is not None and isinstance(sex, (int, float)):
                    if hasattr(self.demographics, 'binary_encodings') and 'SEX' in self.demographics.binary_encodings:
                        # Find the original value that maps to this binary code
                        encoding = self.demographics.binary_encodings['SEX']
                        sex = next((k for k, v in encoding.items() if v == sex), sex)
                    else:
                        sex = 'M' if sex == 1 else 'F'  # Default mapping
                
                # Subject directory path
                subject_dir = os.path.join(output_directory, participant_id, session_id)
                
                # Create session-specific temporary directory
                session_tmp_dir = os.path.join(subject_dir, "tmp_clinical_report")
                os.makedirs(session_tmp_dir, exist_ok=True)
                
                # Generate report - save directly in subject directory, not in separate reports folder
                report_path = generate_clinical_report(
                    sid=participant_id,
                    ses=session_id,
                    age=age,
                    sex=sex,
                    analyses=analyses,
                    features=features,
                    approach=approach,
                    threshold=threshold,
                    threshold_alpha=threshold_alpha,
                    color_range=color_range,
                    cmap=cmap,
                    cmap_asymmetry=cmap_asymmetry,
                    color_bar=color_bar,
                    tmp_dir=session_tmp_dir,  # Use session-specific tmp directory
                    subject_dir=subject_dir,
                    output_dir=None,  # Don't use separate output directory
                    tag=f"{participant_id}_{session_id}_{approach}_clinical_report",
                    smooth_ctx=self.cortical_smoothing,
                    smooth_hip=self.hippocampal_smoothing,
                    verbose=verbose
                )
                
                generated_reports.append(report_path)
                
                if verbose:
                    print(f"Generated report for {participant_id}/{session_id}: {report_path}")
                
                # Clean up temporary directory after successful report generation
                import shutil
                if os.path.exists(session_tmp_dir):
                    shutil.rmtree(session_tmp_dir)
                    if verbose:
                        print(f"  Cleaned up temporary directory: {session_tmp_dir}")
                    
            except Exception as e:
                if verbose:
                    print(f"Error generating report for {participant_id}/{session_id}: {e}")
                # Clean up temporary directory even on error
                import shutil
                session_tmp_dir = os.path.join(output_directory, participant_id, session_id, "tmp_clinical_report")
                if os.path.exists(session_tmp_dir):
                    shutil.rmtree(session_tmp_dir)
                continue
        
        if verbose:
            print(f"Generated {len(generated_reports)} clinical reports in subject directories")
            
        return generated_reports
    
    def _copy_structural_files(self, participant_id, session_id, output_directory, verbose=False):
        """
        Copy structural surface files from micapipe directory to output directory.
        
        Parameters:
        -----------
        participant_id : str
            Subject ID
        session_id : str
            Session ID
        output_directory : str
            Base output directory
        verbose : bool
            Whether to print verbose messages
        
        Returns:
        --------
        bool
            Whether the copying was successful
        """
        import shutil
        
        # Define source and target paths
        structural_output_dir = os.path.join(output_directory, participant_id, session_id, "structural")
        os.makedirs(structural_output_dir, exist_ok=True)
        
        # Files to copy
        surface_files = [
            f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
            f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
            f"anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz"
        ]
        
        success = True
        for file_path in surface_files:
            source_file = os.path.join(self.micapipe_directory, participant_id, session_id, file_path)
            target_file = os.path.join(structural_output_dir, os.path.basename(file_path))
            
            if os.path.exists(source_file):
                if verbose:
                    print(f"  Copying structural file: {os.path.basename(source_file)}")
                try:
                    shutil.copy2(source_file, target_file)
                except Exception as e:
                    if verbose:
                        print(f"  Error copying {source_file}: {e}")
                    success = False
            else:
                if verbose:
                    print(f"  Warning: Structural file not found: {source_file}")
                success = False
                
        return success