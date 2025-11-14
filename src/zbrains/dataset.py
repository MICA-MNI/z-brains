import os
import sys
import datetime
from zbrains.processing import apply_blurring, apply_hippocampal_processing, apply_subcortical_processing, apply_cortical_processing
from zbrains.analysis import analyze_dataset
from zbrains.clinical_reports import generate_clinical_report
import shutil
import multiprocessing
from joblib import Parallel, delayed
import nibabel as nib
import numpy as np

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
    def __init__(self, csv_file, column_mapping=None, normative_columns=None, normative_dtypes=None, reference=None, subset=None):
        self.csv_file = csv_file
        self.data = None
        self.column_mapping = {
            "ID": "participant_id",
            "SES": "session_id",
        } if column_mapping is None else column_mapping
        self.normative_columns = normative_columns
        self.normative_dtypes = normative_dtypes
        self.reference = reference
        self.subset = subset
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
        
        if self.subset is not None:
            subset_pairs = []
            for item in self.subset:
                if not isinstance(item, (list, tuple)) or len(item) != 2:
                    raise ValueError(f"Subset entries must be 2-element lists/tuples. Invalid entry: {item}")
                subset_pairs.append(tuple(item))
            subset_set = set(subset_pairs)
            available_pairs = set(zip(self.data['participant_id'], self.data['session_id']))
            missing_pairs = [pair for pair in subset_pairs if pair not in available_pairs]
            if missing_pairs:
                raise ValueError(f"Subset entries not found in demographics data: {missing_pairs}")
            mask = self.data[['participant_id', 'session_id']].apply(tuple, axis=1).isin(subset_set)
            self.data = self.data[mask].reset_index(drop=True)
            if self.data.empty:
                raise ValueError("Subset filtering resulted in empty demographics data.")
            print(f"Filtered demographics to subset with {len(self.data)} entries.")
        
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
        
        # Standardize feature names consistently - apply mapping first
        feature_mapping = {
            "fa": "FA",
            "adc": "ADC", 
            "thickness": "thickness",
            "flair": "FLAIR",
            "qt1": "qT1",
            "qt1*blur": "qT1*blur",
            "flair*blur": "FLAIR*blur",
        }
        
        # Normalize all features to consistent case
        normalized_features = []
        for feature in features:
            # Convert to lowercase for mapping lookup
            feature_lower = str(feature).lower()
            
            # Apply mapping to get standardized case
            if feature_lower in feature_mapping:
                normalized_features.append(feature_mapping[feature_lower])
            else:
                # For unknown features, keep original case but warn user
                normalized_features.append(str(feature))
                if verbose:
                    print(f"Warning: Unknown feature '{feature}' - using as provided")
        
        # Store the normalized features
        self.features = normalized_features
        
        if verbose:
            print(f"Features {self.features} added to the dataset {self.name}.")
        
        # Define file patterns for each feature (using standardized names)
        feature_files = {
            "FA": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_FA.func.gii",
                   "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-FA.nii.gz"],
            "ADC": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_ADC.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_model-DTI_map-ADC.nii.gz"],
            "thickness": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-thickness.func.gii"],
            "FLAIR": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_flair.func.gii",
                      "maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_T1map.func.gii", 
                    "maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
            "FLAIR*blur": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_flair.func.gii",
                           "maps/{participant_id}_{session_id}_space-nativepro_map-flair.nii.gz"],
            "qT1*blur": ["maps/{participant_id}_{session_id}_hemi-{hemi}_surf-fsLR-32k_label-{surfacetype}_T1map.func.gii",
                         "maps/{participant_id}_{session_id}_space-nativepro_map-T1map.nii.gz"],
        }
        
        # Rest of the method remains the same...
        # Define required files for different structures
        required_cortical_files = [
            "anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz",
            # "anat/{participant_id}_{session_id}_space-nativepro_T1w_brain_mask.nii.gz",
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
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
        ]
        
        required_hippocampal_files = [
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii"
        ]

        required_freesurfer_files = [
            "{participant_id}_{session_id}/surf/lh.white",
            "{participant_id}_{session_id}/surf/rh.white",
            "{participant_id}_{session_id}/label/lh.cortex.label",
            "{participant_id}_{session_id}/label/rh.cortex.label"
        ]

        required_subcortical_files = [
            "{participant_id}_{session_id}/stats/aseg.stats",
        ]

        # Track missing files and feature availability
        missing_files = {}
        missing_features = {}
        feature_availability = {feature: 0 for feature in self.features}
        subjects_with_complete_data = []
        
        # Initialize feature-specific valid subjects with nested structure
        for feature in self.features:
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
                    print(f"Missing required cortical file for subject {participant_id}/{session_id}: {file_path}")
                    has_required_cortical = False

            for file_pattern in required_freesurfer_files:
                file_path = os.path.join(
                    self.freesurfer_directory,
                    file_pattern.format(participant_id=participant_id, session_id=session_id)
                )
                if not os.path.exists(file_path):
                    subject_missing_files.append(file_path)
                    print(f"Missing required cortical file for subject {participant_id}/{session_id}: {file_path}")
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
                        print(f"Missing required hippocampal file for subject {participant_id}/{session_id}: {file_path}")
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
                        print(f"Missing required subcortical file for subject {participant_id}/{session_id}: {file_path}")
                        has_required_subcortical = False
            
            # Update structure-specific lists
            if has_required_cortical:
                self.valid_subjects['structures']['cortex'].append(subject)
            
            if has_required_hippocampal:
                self.valid_subjects['structures']['hippocampus'].append(subject)
            
            if has_required_subcortical:
                self.valid_subjects['structures']['subcortical'].append(subject)
            
            # Check feature-specific files and structure availability
            for feature in self.features:
                feature_complete = True
                # Keep track of which structures have this feature
                feature_has_cortical = has_required_cortical
                feature_has_hippocampal = has_required_hippocampal
                feature_has_subcortical = has_required_subcortical
                
                if feature not in feature_files:
                    if verbose:
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
                    
                    if feature_has_hippocampal and not feature.lower().endswith('*blur'):
                        # Blur features don't apply to hippocampus
                        self.valid_subjects[feature]['structures']['hippocampus'].append(subject)
                    
                    if feature_has_subcortical and not feature.lower().endswith('*blur'):
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
        for feature in self.features:
            if feature in feature_availability:
                avail_count = feature_availability[feature]
                print(f"  {feature}: {avail_count}/{total_subjects} subjects ({avail_count/total_subjects*100:.1f}%)")
                
                # Print structure-specific availability for this feature
                if verbose:
                    cortex_count = len(self.valid_subjects[feature]['structures']['cortex'])
                    print(f"    - cortex: {cortex_count}/{total_subjects} subjects ({cortex_count/total_subjects*100:.1f}%)")
                    
                    if self.hippocampus and not feature.lower().endswith('*blur'):
                        hippo_count = len(self.valid_subjects[feature]['structures']['hippocampus'])
                        print(f"    - hippocampus: {hippo_count}/{total_subjects} subjects ({hippo_count/total_subjects*100:.1f}%)")
                    
                    if self.subcortical and not feature.lower().endswith('*blur'):
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
        # Use the same consistent feature mapping
        feature_mapping = {
            "fa": "FA",
            "adc": "ADC",
            "thickness": "thickness",
            "flair": "FLAIR",
            "qt1": "qT1",
            "flair*blur": "FLAIR*blur",
            "qt1*blur": "qT1*blur"
        }
        
        # Normalize features consistently
        normalized_features = []
        for feature in features:
            feature_lower = str(feature).lower()
            if feature_lower in feature_mapping:
                normalized_features.append(feature_mapping[feature_lower])
            else:
                normalized_features.append(str(feature))
                if verbose:
                    print(f"Warning: Unknown feature '{feature}' in process method")
        
        features = normalized_features
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

        # Timestamp for all logs
        timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        
        # NO MAIN LOG FILE CREATION - REMOVED

        # Determine number of jobs
        if n_jobs is None:
            n_jobs = env.num_threads or multiprocessing.cpu_count()
        
        # Use all base valid subjects
        valid_subjects_to_process = self.valid_subjects['base']
        
        if verbose:
            print(f"Using {n_jobs if n_jobs > 1 else 'sequential'} processing for {len(valid_subjects_to_process)} subjects")

        # Identify blur features
        blur_features = [feature for feature in self.features if feature.endswith("*blur")]
        base_features = [feature.replace("*blur", "") for feature in blur_features]
        
        # Define a function to process a single subject
        def process_single_subject(subject):
            participant_id, session_id = subject
            
            # Create session-specific directories
            session_output_dir = os.path.join(output_directory, participant_id, session_id)
            os.makedirs(session_output_dir, exist_ok=True)
            
            # Create session log file - kept in subject directory
            log_file = os.path.join(session_output_dir, f"{participant_id}_{session_id}_processing_{timestamp}.log")
            
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
                
                # Add dataset-level information to each subject log
                print(f"Dataset: {self.name}")
                print(f"Total subjects in dataset: {len(valid_subjects_to_process)}")
                print(f"Processing job: {n_jobs if n_jobs > 1 else 'sequential'}")
                print("-" * 50)
                
                try:
                    # Copy structural files
                    self._copy_structural_files(participant_id, session_id, output_directory, verbose=verbose)
                    self._create_midline_from_freesurfer(participant_id, session_id, output_directory, verbose=verbose)
                    # Apply blurring to features that need it
                    if blur_features:
                        # Check which base features are available for this subject
                        available_features = []
                        for feature in base_features:
                            if subject in self.valid_subjects[feature]['all']:
                                available_features.append(feature)
                        
                        if available_features:  # If ANY features are available, proceed with blurring
                            if verbose:
                                print(f"  Applying additional blur processing for features: {', '.join([f+'*blur' for f in available_features])}")
                            
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
                        # Get nonBlur features for hippocampus
                        valid_hipp_features = [f for f in self.features 
                                            if not f.endswith("*blur") 
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
                        # Get nonBlur features for subcortical
                        valid_subcort_features = [f for f in self.features 
                                                if not f.endswith("*blur") 
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

        # Create a summary file in each subject's directory
        for participant_id, session_id in valid_subjects_to_process:
            session_dir = os.path.join(output_directory, participant_id, session_id)
            summary_file = os.path.join(session_dir, f"processing_summary_{timestamp}.txt")
            
            try:
                with open(summary_file, 'w', encoding='utf-8') as f:
                    f.write(f"===== PROCESSING SUMMARY FOR {participant_id}/{session_id} =====\n")
                    f.write(f"Dataset: {self.name}\n")
                    f.write(f"Completed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Total subjects processed: {len(valid_subjects_to_process)}\n")
                    f.write(f"Subjects processed successfully: {len(valid_subjects_to_process) - len(failed_subjects)}/{len(valid_subjects_to_process)}\n")
                    
                    if failed_subjects:
                        f.write(f"Failed subjects ({len(failed_subjects)}):\n")
                        for pid, sid in failed_subjects:
                            f.write(f"  - {pid}/{sid}\n")
                            
                    # Mark this subject's status
                    is_failed = (participant_id, session_id) in failed_subjects
                    f.write(f"\nThis subject status: {'FAILED' if is_failed else 'SUCCESS'}\n")
            except Exception as e:
                if verbose:
                    print(f"Error writing summary for {participant_id}/{session_id}: {e}")

        print(f"Dataset {self.name} processed. Summary files written to each subject directory.")
        
        self.validate(
            features=self.features,
            output_directory=output_directory,
            cortical_smoothing=cortical_smoothing,
            hippocampal_smoothing=hippocampal_smoothing,
            verbose=verbose
        )
        return self

    def validate(self, output_directory, features=None, cortical_smoothing=5, hippocampal_smoothing=2,
                 verbose=True, error_threshold=0, warning_threshold=25, only_demographics_subjects=True):
        """
        Validate that all expected processed files exist for each subject.

        Returns
        -------
        dict
            {
                "valid_subjects": [(participant_id, session_id), ...],
                "missing_files": { (participant_id, session_id): [paths...] },
                "summary": {...}
            }
        """
        feature_case_mapping = {
            "fa": "FA",
            "adc": "ADC",
            "thickness": "thickness",
            "flair": "FLAIR",
            "qt1": "qT1",
            "t1map": "qT1",
            "flair*blur": "FLAIR*blur",
            "qt1*blur": "qT1*blur"
        }

        if features is None:
            normalized_features = list(self.features)
        else:
            normalized_features = []
            for feat in features:
                feat_lower = str(feat).lower()
                if feat_lower.endswith("*blur"):
                    base = feat_lower.replace("*blur", "")
                    mapped = feature_case_mapping.get(feat_lower, f"{feature_case_mapping.get(base, feat[:-5])}*blur")
                    normalized_features.append(mapped)
                else:
                    normalized_features.append(feature_case_mapping.get(feat_lower, str(feat)))
        if not normalized_features:
            normalized_features = list(self.features)
        self.features = normalized_features

        if not hasattr(self, "valid_subjects") or "base" not in self.valid_subjects:
            self.check_directories()

        output_subjects = self.check_output_directories(
            output_directory,
            only_demographics_subjects=only_demographics_subjects,
            verbose=verbose
        )
        if not output_subjects:
            raise ValueError(f"No subject output directories found in {output_directory}. Nothing to validate.")

        if verbose:
            print(f"Validating dataset {self.name} with {len(output_subjects)} output subjects "
                  f"and {len(self.features)} features...")

        feature_output_mapping = {
            "thickness": "thickness",
            "flair": "FLAIR",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "qT1",
            "t1map": "qT1"
        }

        feature_meta = []
        for feat in self.features:
            feat_lower = feat.lower()
            is_blur = feat_lower.endswith("*blur")
            base_lower = feat_lower.replace("*blur", "") if is_blur else feat_lower
            cortical_token = feat  # already normalised to processing output naming

            if is_blur:
                if base_lower in {"qt1", "t1map"}:
                    blur_token = "qT1"
                elif base_lower == "flair":
                    blur_token = "FLAIR"
                elif base_lower in {"adc", "fa"}:
                    blur_token = base_lower
                else:
                    blur_token = feat[:-5]
            else:
                blur_token = None

            hippo_token = None if is_blur else feature_output_mapping.get(base_lower, feat)
            subcortical_token = None
            if not is_blur and base_lower not in {"thickness", "volume"}:
                subcortical_token = feature_output_mapping.get(base_lower, feat.upper())

            feature_meta.append({
                "original": feat,
                "base_lower": base_lower,
                "is_blur": is_blur,
                "cortical_token": cortical_token,
                "blur_token": blur_token,
                "hippo_token": hippo_token,
                "subcortical_token": subcortical_token,
                "requires_cortex": self.cortex,
                "requires_hippocampus": (
                    self.hippocampus and self.hippunfold_directory and not is_blur and hippo_token is not None
                ),
                "requires_subcortical": (
                    self.subcortical and self.freesurfer_directory and not is_blur and subcortical_token is not None
                )
            })

        resolutions = ["32k", "5k"]
        surface_labels = ["midthickness", "white"]
        hippocampal_resolution = "0p5mm"
        blur_suffixes = [
            "_surf-fsnative_desc-raw.func.gii",
            "_surf-fsnative_desc-dist.func.gii",
            "_surf-fsnative_desc-grad.func.gii"
        ]

        processed_valid_subjects = {
            'base': [],
            'structures': {
                'cortex': [],
                'hippocampus': [],
                'subcortical': []
            }
        }
        feature_processed = {
            meta['original']: {
                'all': [],
                'structures': {
                    'cortex': [],
                    'hippocampus': [],
                    'subcortical': []
                }
            } for meta in feature_meta
        }

        missing_files = {}
        missing_features = {}
        feature_availability = {meta['original']: 0 for meta in feature_meta}
        subjects_with_complete_data = []
        all_files_count = 0
        missing_files_count = 0

        for participant_id, session_id in output_subjects:
            subject_missing = []
            subject_missing_features = set()
            bids_id = f"{participant_id}_{session_id}"
            subject_dir = os.path.join(output_directory, participant_id, session_id)
            cortex_dir = os.path.join(subject_dir, "maps", "cortex")
            hippocampus_dir = os.path.join(subject_dir, "maps", "hippocampus")
            subcortical_dir = os.path.join(subject_dir, "maps", "subcortical")
            structural_dir = os.path.join(subject_dir, "structural")

            subject_structures = {
                'cortex': True if self.cortex else None,
                'hippocampus': True if self.hippocampus and self.hippunfold_directory else None,
                'subcortical': True if self.subcortical and self.freesurfer_directory else None
            }
            subject_feature_structures = {
                meta['original']: {
                    'cortex': True if meta['requires_cortex'] else None,
                    'hippocampus': True if meta['requires_hippocampus'] else None,
                    'subcortical': True if meta['requires_subcortical'] else None
                } for meta in feature_meta
            }

            structural_expected = [
                f"{bids_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii",
                f"{bids_id}_hemi-L_surf-fsnative_label-sphere.surf.gii",
                f"{bids_id}_hemi-R_surf-fsnative_label-sphere.surf.gii",
                f"{bids_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii",
                f"{bids_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
                f"{bids_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii",
                f"{bids_id}_space-nativepro_T1w.nii.gz",
                f"{bids_id}_hemi-L_surf-fsnative_label-medialwall.label.gii",
                f"{bids_id}_hemi-R_surf-fsnative_label-medialwall.label.gii"
            ]
            
            # Only require Laplace file if blur features are present
            has_blur_features = any(feat.lower().endswith("*blur") for feat in self.features)
            if has_blur_features:
                structural_expected.append(f"{participant_id}_{session_id}-laplace.nii.gz")
            
            if self.hippocampus:
                structural_expected.extend([
                    f"{bids_id}_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-R_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
                    f"{bids_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii"
                ])

            for filename in structural_expected:
                path = os.path.join(structural_dir, filename)
                all_files_count += 1
                if not os.path.exists(path):
                    subject_missing.append(path)
                    missing_files_count += 1
                    if "hipp" in filename and subject_structures['hippocampus'] is not None:
                        subject_structures['hippocampus'] = False
                    elif "laplace" in filename or filename.endswith(".surf.gii") or filename.endswith(".label.gii"):
                        if subject_structures['cortex'] is not None:
                            subject_structures['cortex'] = False

            if self.cortex:
                for meta in feature_meta:
                    if not meta['requires_cortex']:
                        continue
                    feature_name = meta['original']
                    for hemi in ["L", "R"]:
                        for res in resolutions:
                            for label in (["midthickness"] if meta['is_blur'] else surface_labels):
                                cortical_path = os.path.join(
                                    cortex_dir,
                                    f"{bids_id}_hemi-{hemi}_surf-fsLR-{res}_label-{label}_feature-"
                                    f"{meta['cortical_token']}_smooth-{cortical_smoothing}mm.func.gii"
                                )
                                all_files_count += 1
                                if not os.path.exists(cortical_path):
                                    subject_missing.append(cortical_path)
                                    missing_files_count += 1
                                    subject_feature_structures[feature_name]['cortex'] = False
                                    subject_missing_features.add(feature_name)
                    if meta['is_blur'] and meta['blur_token']:
                        for hemi in ["L", "R"]:
                            blur_prefix = os.path.join(
                                cortex_dir,
                                f"{participant_id}_{session_id}_hemi-{hemi}_feature-{meta['blur_token']}*blur"
                            )
                            blur_paths = [f"{blur_prefix}{suffix}" for suffix in blur_suffixes]
                            blur_paths.append(f"{blur_prefix}_surf-fsnative_smooth-{cortical_smoothing}mm.func.gii")
                            for blur_path in blur_paths:
                                all_files_count += 1
                                if not os.path.exists(blur_path):
                                    subject_missing.append(blur_path)
                                    missing_files_count += 1
                                    subject_feature_structures[feature_name]['cortex'] = False
                                    subject_missing_features.add(feature_name)

            if self.hippocampus and self.hippunfold_directory:
                for meta in feature_meta:
                    if not meta['requires_hippocampus']:
                        continue
                    feature_name = meta['original']
                    for hemi in ["L", "R"]:
                        hippo_path = os.path.join(
                            hippocampus_dir,
                            f"{participant_id}_{session_id}_hemi-{hemi}_den-{hippocampal_resolution}_label-hipp_"
                            f"midthickness_feature-{meta['hippo_token']}_smooth-{hippocampal_smoothing}mm.func.gii"
                        )
                        all_files_count += 1
                        if not os.path.exists(hippo_path):
                            subject_missing.append(hippo_path)
                            missing_files_count += 1
                            subject_feature_structures[feature_name]['hippocampus'] = False
                            subject_missing_features.add(feature_name)

            if self.subcortical and self.freesurfer_directory:
                volume_file = os.path.join(subcortical_dir, f"{bids_id}_feature-volume.csv")
                all_files_count += 1
                if not os.path.exists(volume_file):
                    subject_missing.append(volume_file)
                    missing_files_count += 1
                    if subject_structures['subcortical'] is not None:
                        subject_structures['subcortical'] = False

                for meta in feature_meta:
                    if not meta['requires_subcortical']:
                        continue
                    feature_name = meta['original']
                    subcort_path = os.path.join(
                        subcortical_dir,
                        f"{bids_id}_feature-{meta['subcortical_token']}.csv"
                    )
                    all_files_count += 1
                    if not os.path.exists(subcort_path):
                        subject_missing.append(subcort_path)
                        missing_files_count += 1
                        subject_feature_structures[feature_name]['subcortical'] = False
                        subject_missing_features.add(feature_name)

            if any(flag is True for flag in subject_structures.values() if flag is not None):
                processed_valid_subjects['base'].append((participant_id, session_id))
                for structure_name, flag in subject_structures.items():
                    if flag is True:
                        processed_valid_subjects['structures'][structure_name].append((participant_id, session_id))

            for meta in feature_meta:
                feature_name = meta['original']
                structure_flags = subject_feature_structures[feature_name]
                relevant_flags = [flag for flag in structure_flags.values() if flag is not None]
                if not relevant_flags:
                    continue
                if any(flag is True for flag in relevant_flags):
                    feature_processed[feature_name]['all'].append((participant_id, session_id))
                    for structure_name, flag in structure_flags.items():
                        if flag is True:
                            feature_processed[feature_name]['structures'][structure_name].append((participant_id, session_id))
                else:
                    subject_missing_features.add(feature_name)

            if subject_missing_features:
                missing_features[(participant_id, session_id)] = subject_missing_features

            if subject_missing:
                missing_files[(participant_id, session_id)] = subject_missing
            else:
                subjects_with_complete_data.append((participant_id, session_id))

        for feature_name, data in feature_processed.items():
            data['all'] = list(dict.fromkeys(data['all']))
            for structure_name in data['structures']:
                data['structures'][structure_name] = list(dict.fromkeys(data['structures'][structure_name]))
            feature_availability[feature_name] = len(data['all'])

        for structure_name in processed_valid_subjects['structures']:
            processed_valid_subjects['structures'][structure_name] = list(dict.fromkeys(processed_valid_subjects['structures'][structure_name]))

        processed_valid_subjects.update(feature_processed)

        self.valid_subjects = processed_valid_subjects
        self.subjects_with_complete_data = subjects_with_complete_data
        self.missing_files = missing_files
        self.missing_features = missing_features
        self.feature_availability = feature_availability

        total_subjects = len(output_subjects)
        complete_subjects = len(subjects_with_complete_data)
        complete_percentage = (complete_subjects / total_subjects) * 100 if total_subjects else 0.0
        file_presence_percentage = ((all_files_count - missing_files_count) / all_files_count) * 100 if all_files_count else 0.0

        summary = {
            "total_subjects": total_subjects,
            "complete_subjects": complete_subjects,
            "complete_percentage": complete_percentage,
            "total_files_expected": all_files_count,
            "missing_files_count": missing_files_count,
            "file_presence_percentage": file_presence_percentage
        }
        print("\nValidation Summary:")
        print(f"  Total subjects (outputs): {total_subjects}")
        print(f"  Subjects with complete data: {complete_subjects} ({complete_percentage:.1f}%)")
        print(f"  Total expected files: {all_files_count}")
        print(f"  Missing files: {missing_files_count} ({100 - file_presence_percentage:.1f}%)")
        if verbose and missing_files:
            print(f"\nSubjects with missing files: {len(missing_files)}")
            for (pid, sid), files in missing_files.items():
                print(f"  {pid}/{sid}: {len(files)} missing")
                for path in files[:10]:
                    print(f"    - {os.path.basename(path)}")
                if len(files) > 10:
                    print(f"    - ... and {len(files) - 10} more")

        if complete_subjects == 0:
            self.valid_dataset = False
            raise ValueError("Validation failed: no subjects have complete processed data.")

        if error_threshold > 0 and complete_percentage < error_threshold:
            self.valid_dataset = False
            raise ValueError(f"Validation failed: {complete_percentage:.1f}% complete "
                             f"(threshold {error_threshold}%).")

        if warning_threshold and complete_percentage < warning_threshold:
            import warnings
            warnings.warn(f"Validation warning: only {complete_percentage:.1f}% of subjects have complete data "
                          f"(warning threshold {warning_threshold}%).")

        self.valid_dataset = True
        self.cortical_smoothing = cortical_smoothing
        self.hippocampal_smoothing = hippocampal_smoothing

        if verbose:
            print("\nFeature availability summary:")
            denom = total_subjects if total_subjects else 1
            for feature_name in self.features:
                avail_count = feature_availability.get(feature_name, 0)
                print(f"  {feature_name}: {avail_count}/{total_subjects} subjects ({(avail_count/denom)*100:.1f}%)")
                if feature_name in feature_processed:
                    for structure_name in ['cortex', 'hippocampus', 'subcortical']:
                        struct_list = feature_processed[feature_name]['structures'][structure_name]
                        if struct_list:
                            print(f"    - {structure_name}: {len(struct_list)}/{total_subjects} subjects ({(len(struct_list)/denom)*100:.1f}%)")

            print("\nStructure availability summary:")
            for structure_name in ['cortex', 'hippocampus', 'subcortical']:
                struct_list = processed_valid_subjects['structures'][structure_name]
                if struct_list:
                    print(f"  {structure_name}: {len(struct_list)}/{total_subjects} subjects ({(len(struct_list)/denom)*100:.1f}%)")

        return {
            "valid_subjects": processed_valid_subjects['base'],
            "missing_files": missing_files,
            "summary": summary
        }
    
    def check_output_directories(self, output_directory, only_demographics_subjects=True, verbose=True):
        """
        Check the output directory structure and identify valid subjects with output data.
        This is used during validation to determine which subjects to validate.
        
        Parameters:
        -----------
        output_directory : str
            Directory where processed data should be stored
        only_demographics_subjects : bool, default=True
            If True, only check directories for subjects in the demographics data
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
        
        # If we're only checking demographics subjects, create a list of valid subject IDs
        valid_demo_subjects = []
        if only_demographics_subjects and self.demographics is not None:
            for _, row in self.demographics.data.iterrows():
                valid_demo_subjects.append((row['participant_id'], row['session_id']))
            
            if verbose:
                print(f"Filtering to {len(valid_demo_subjects)} subjects specified in demographics data.")
        
        valid_output_subjects = []
        missing_cortex = []
        missing_hippocampus = []
        missing_subcortical = []
        skipped_non_demo = []
        
        # Check each subject directory for session directories
        for subject_dir in subject_dirs:
            participant_id = subject_dir
            subject_path = os.path.join(output_directory, subject_dir)
            
            # Find session directories
            session_dirs = [d for d in os.listdir(subject_path) 
                        if os.path.isdir(os.path.join(subject_path, d)) and d.startswith('ses-')]
            
            for session_dir in session_dirs:
                session_id = session_dir
                
                # Skip if only checking demographics subjects and this subject isn't in the demographics
                if only_demographics_subjects and (participant_id, session_id) not in valid_demo_subjects:
                    skipped_non_demo.append((participant_id, session_id))
                    continue
                    
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
            if skipped_non_demo:
                print(f"Skipped {len(skipped_non_demo)} subject directories not found in demographics data.")
                if len(skipped_non_demo) < 6:
                    for p, s in skipped_non_demo:
                        print(f"  - {p}/{s}")
                else:
                    for p, s in skipped_non_demo[:5]:
                        print(f"  - {p}/{s}")
                    print(f"  - ... and {len(skipped_non_demo) - 5} more")
            
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
            Analysis method to use. Currently supports 'zscore'and 'wscore'
        output_directory : str, optional
            Directory where z-score maps will be saved. If None, uses validation output directory
        verbose : bool, default=True
            If True, prints detailed information about the analysis process
        
        Returns:
        --------
        dict
            Dictionary containing analysis results for each feature and region
        """
        
        # Call the analysis function and store results
        results = analyze_dataset(self, reference, method, output_directory, verbose)
        self.analysis_results = results
        self.reference_demographics = reference.demographics
        return results
    
    def clinical_report(self, output_directory=None, approach='wscore', analyses=['regional','asymmetry'], features=None, 
                    threshold=1.96, threshold_alpha=0.3, color_range=(-3, 3), 
                    cmap='cmo.balance', cmap_asymmetry='cmo.balance_r', 
                    color_bar='bottom', tmp_dir=None, env=None, verbose=True):
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

        if env == None:
            raise ValueError("env must be specified to access workbench and other paths")

        # Check if analysis has been run
        if not hasattr(self, 'analysis_results'):
            raise ValueError("No analysis results found. Please run dataset.analyze() first.")
            
        # Use features from the dataset if not specified
        if features is None:
            features = list(self.features)

        feature_mapping = {
            "thickness": "thickness",
            "flair": "FLAIR",
            "adc": "ADC",
            "fa": "FA",
            "qt1": "qT1",
            "qt1*blur": "qT1*blur",
            "flair*blur": "FLAIR*blur"
        }

        features = [feature_mapping.get(f.lower(), f) for f in features]

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
                sex = int(sex)
                # Convert binary sex encoding back to string if needed
                if sex is not None and isinstance(sex, (int, float)):
                    if hasattr(self.reference_demographics, 'binary_encodings') and 'SEX' in self.demographics.binary_encodings:
                        # Find the original value that maps to this binary code
                        print("ref binary encodings:", self.reference_demographics.binary_encodings)
                        encoding = self.reference_demographics.binary_encodings['SEX']
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
                    env=env,
                    verbose=verbose,
                    control_data=self.reference_demographics,
                    valid_subjects=self.valid_subjects,
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
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_surf-fsnative_label-sphere.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_surf-fsnative_label-sphere.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsnative_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"),
            os.path.join(self.micapipe_directory, participant_id, session_id, f"anat/{participant_id}_{session_id}_space-nativepro_T1w.nii.gz"),
        ]
        
        if self.hippocampus:
            surface_files.append(os.path.join(self.hippunfold_directory, "hippunfold", participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii"))
            surface_files.append(os.path.join(self.hippunfold_directory, "hippunfold", participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii"))
            surface_files.append(os.path.join(self.hippunfold_directory, "hippunfold", participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii"))
            surface_files.append(os.path.join(self.hippunfold_directory, "hippunfold", participant_id, session_id, f"surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii"))

        success = True
        for file_path in surface_files:
            source_file = os.path.join(file_path)
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
    
    def _create_midline_from_freesurfer(self, participant_id, session_id, output_directory, verbose=False):
        # Pick the surface you consider "fsnative" (coords/triangles don’t actually matter for the mask)
        lh_surf = os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}/surf/lh.white")
        rh_surf = os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}/surf/rh.white")
        lh_coords, _ = nib.freesurfer.read_geometry(str(lh_surf))
        rh_coords, _ = nib.freesurfer.read_geometry(str(rh_surf))
        n_lh, n_rh = lh_coords.shape[0], rh_coords.shape[0]

        # Read cortex labels (indices of vertices that are true cortex)
        lh_cortex_idx = nib.freesurfer.read_label(os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}", "label", "lh.cortex.label"))
        rh_cortex_idx = nib.freesurfer.read_label(os.path.join(self.freesurfer_directory,f"{participant_id}_{session_id}", "label", "rh.cortex.label"))

        # Boolean masks (True = cortex, False = medial wall)
        lh_is_cortex = np.zeros(n_lh, dtype=bool); lh_is_cortex[lh_cortex_idx] = True
        rh_is_cortex = np.zeros(n_rh, dtype=bool); rh_is_cortex[rh_cortex_idx] = True

        # Medial-wall masks (True = medial wall)
        lh_medial_wall = ~lh_is_cortex
        rh_medial_wall = ~rh_is_cortex

        structural_output_dir = os.path.join(output_directory, participant_id, session_id, "structural")
        os.makedirs(structural_output_dir, exist_ok=True)

        lh_medial_wall_file = os.path.join(structural_output_dir, f"{participant_id}_{session_id}_hemi-L_surf-fsnative_label-medialwall.label.gii")
        rh_medial_wall_file = os.path.join(structural_output_dir, f"{participant_id}_{session_id}_hemi-R_surf-fsnative_label-medialwall.label.gii")

        # Save medial wall labels as gifti files
        lh_medial_wall_gifti = nib.gifti.GiftiImage()
        rh_medial_wall_gifti = nib.gifti.GiftiImage()
        lh_medial_wall_gifti.add_gifti_data_array(nib.gifti.GiftiDataArray(data=lh_medial_wall.astype(np.uint8), intent="NIFTI_INTENT_LABEL"))
        rh_medial_wall_gifti.add_gifti_data_array(nib.gifti.GiftiDataArray(data=rh_medial_wall.astype(np.uint8), intent="NIFTI_INTENT_LABEL"))
        nib.save(lh_medial_wall_gifti, lh_medial_wall_file)
        nib.save(rh_medial_wall_gifti, rh_medial_wall_file)
        if verbose:
            print(f"  Created medial wall labels for {participant_id}/{session_id} at {structural_output_dir}")


