import os
from src.processing import apply_blurring, apply_hippocampal_processing, apply_subcortical_processing, apply_cortical_processing

class demographics():
    def __init__(self, csv_file, column_mapping=None, normative_columns=None):
        self.csv_file = csv_file
        self.data = None
        self.column_mapping = {
            "ID": "participant_id",
            "SES": "session_id",
        } if column_mapping is None else column_mapping
        self.normative_columns = normative_columns
        self.load_data()

    def __repr__(self):
        return f"Demographics(csv_file={self.csv_file})"
    def __str__(self):
        return f"Demographics data loaded from {self.csv_file}"
    
    def load_data(self):
        import pandas as pd
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
        normative_columns = self.normative_columns if self.normative_columns else []
        for col in normative_columns:
            if col not in self.data.columns:
                raise ValueError(f"Normative column '{col}' not found in demographics data.")
        # Check to see if any normative columns are missing for any participant/session
        for col in normative_columns:
            if self.data[col].isnull().any():
                raise ValueError(f"Normative column '{col}' contains missing values.")
        print(f"Normative columns: {normative_columns} sucessfully validated.")
    
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
        self.check_directories()

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

        # Check if demographics data matches micapipe directory
        missing_micapipe = []
        missing_hippunfold = []
        missing_freesurfer = []
        valid_subjects = []  # Track valid subjects with complete data
        
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
                valid_subjects.append((participant_id, session_id))
        
        # Store valid subjects in the class instance
        self.valid_subjects = valid_subjects
        
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
        
        print(f"Found {len(valid_subjects)} valid subjects with complete directory structure.")
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

        # Define required surfaces
        required_files = [
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
        ]
        # If hippunfold is enabled, add hippocampal surfaces
        required_hippocampal_files = [
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_inner.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_outer.surf.gii",
            "surf/{participant_id}_{session_id}_hemi-R_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii",
            "surf/{participant_id}_{session_id}_hemi-L_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii"
        ]

        required_subcortical_files = [
            "{participant_id}_{session_id}/stats/aseg.stats",
        ]

        # Track missing files for each feature and subject
        missing_files = {}
        missing_features = {}  # Track which features are missing for each subject
        feature_availability = {feature: 0 for feature in features}  # Count availability per feature
        subjects_with_complete_data = []
        missing_required_surfaces = []  # Track subjects missing required surfaces
        missing_hippocampal = []  # Track subjects missing hippocampal files
        missing_subcortical = []  # Track subjects missing subcortical files
        
        # Check for each valid subject
        for participant_id, session_id in self.valid_subjects:
            subject_missing_files = []
            subject_missing_features = set()
            has_required_surfaces = True
            has_required_hippocampal = True
            has_required_subcortical = True
            
            # Check required surfaces
            for surface_pattern in required_files:
                surface_path = os.path.join(
                    self.micapipe_directory,
                    participant_id,
                    session_id,
                    surface_pattern.format(participant_id=participant_id, session_id=session_id)
                )
                if not os.path.exists(surface_path):
                    subject_missing_files.append(surface_path)
                    has_required_surfaces = False
            
            # Check hippocampal files if enabled
            if self.hippocampus and self.hippunfold_directory:
                for hipp_pattern in required_hippocampal_files:
                    hipp_path = os.path.join(
                        self.hippunfold_directory,
                        "hippunfold",
                        participant_id,
                        session_id,
                        hipp_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    if not os.path.exists(hipp_path):
                        subject_missing_files.append(hipp_path)
                        has_required_hippocampal = False
            
            # Check subcortical files if enabled
            if self.subcortical and self.freesurfer_directory:
                for subcort_pattern in required_subcortical_files:
                    subcort_path = os.path.join(
                        self.freesurfer_directory,
                        subcort_pattern.format(participant_id=participant_id, session_id=session_id)
                    )
                    if not os.path.exists(subcort_path):
                        subject_missing_files.append(subcort_path)
                        has_required_subcortical = False
            
            # Track missing components
            if not has_required_surfaces:
                missing_required_surfaces.append((participant_id, session_id))
            if not has_required_hippocampal:
                missing_hippocampal.append((participant_id, session_id))
            if not has_required_subcortical:
                missing_subcortical.append((participant_id, session_id))
            
            # Check feature-specific files
            for feature in features:
                feature_complete = True
                
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
                    else:
                        # Pattern has neither hemi nor surfacetype
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
                
                # Track feature availability
                if feature_complete:
                    feature_availability[feature] += 1
                else:
                    subject_missing_features.add(feature)
            
            # Track missing features for this subject
            if subject_missing_features:
                missing_features[(participant_id, session_id)] = subject_missing_features
                
            # Track missing files and complete subjects
            if subject_missing_files:
                missing_files[(participant_id, session_id)] = subject_missing_files
            else:
                subjects_with_complete_data.append((participant_id, session_id))
        
        # Store the results
        self.subjects_with_complete_data = subjects_with_complete_data
        self.missing_files = missing_files
        self.missing_features = missing_features
        self.feature_availability = feature_availability
        
        # Report findings - only show basic info if not verbose
        if verbose:
            if subjects_with_complete_data:
                print(f"Found {len(subjects_with_complete_data)} subjects with complete data for all requested features.")
            
            if missing_required_surfaces:
                print(f"Warning: {len(missing_required_surfaces)} subjects are missing required surface files.")
            
            if self.hippocampus and missing_hippocampal:
                print(f"Warning: {len(missing_hippocampal)} subjects are missing required hippocampal files.")
                
            if self.subcortical and missing_subcortical:
                print(f"Warning: {len(missing_subcortical)} subjects are missing required subcortical files.")
        
        # Always print feature availability summary
        print("\nFeature availability summary:")
        total_subjects = len(self.valid_subjects)
        for feature in features:
            if feature in feature_availability:
                avail_count = feature_availability[feature]
                print(f"  {feature}: {avail_count}/{total_subjects} subjects ({avail_count/total_subjects*100:.1f}%)")
        
        # Only print detailed missing features if verbose
        if verbose and missing_features:
            print(f"\nMissing features by subject:")
            for (pid, sid), missing_feats in missing_features.items():
                print(f"  {pid}/{sid}: missing {', '.join(missing_feats)}")
        
        return self

    # Update the process method to use the blurring functionality
    def process(self, output_directory, cortical_smoothing=5, hippocampal_smoothing=2, env=None, verbose=True):
        """
        Process the dataset with specified smoothing parameters.
        
        Parameters:
        -----------
        output_directory : str
            Directory to store processed data
        cortical_smoothing : int, default=5
            Smoothing parameter for cortical features
        hippocampal_smoothing : int, default=2
            Smoothing parameter for hippocampal features
        env : object
            Environment object containing paths to required tools
        verbose : bool, default=True
            If True, prints detailed processing information
        """
        import shutil
        
        if verbose:
            print(f"Processing dataset {self.name} with cortical smoothing {cortical_smoothing} and hippocampal smoothing {hippocampal_smoothing}.")
        
        if env is None:
            raise ValueError("Environment (zbenv) must be provided for processing.")
        if not self.features:
            raise ValueError("No features have been added to the dataset. Please add features before processing.")
        
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
            print(f"Created output directory: {output_directory}")
        else:
            print(f"Output directory already exists: {output_directory}")

        # Identify blur features
        blur_features = [feature for feature in self.features if feature.endswith("-blur")]
        base_features = [feature.replace("-blur", "") for feature in blur_features]
        
        # Process each subject
        for subject in self.valid_subjects:
            participant_id, session_id = subject
            print(f"Processing subject {participant_id}/{session_id}...")
            
            # Create session-specific tmp directory
            session_tmp_dir = os.path.join(output_directory, f"tmp_{participant_id}_{session_id}")
            os.makedirs(session_tmp_dir, exist_ok=True)
            
            try:
                # Apply blurring to features that need it
                if blur_features:
                    print(f"  Applying additional blur processing for features: {', '.join(blur_features)}")
                    
                    apply_blurring(
                        participant_id=participant_id,
                        session_id=session_id,
                        features=base_features,
                        output_directory=output_directory,
                        workbench_path=env.connectome_workbench_path,
                        micapipe_directory=self.micapipe_directory,
                        freesurfer_directory=self.freesurfer_directory,
                        tmp_dir=session_tmp_dir,
                        verbose=verbose
                    )

                # Process cortical features if cortex is enabled
                if self.cortex:
                    print(f"  Processing cortical data for features: {', '.join(self.features)}")
                    
                    apply_cortical_processing(
                        participant_id=participant_id,
                        session_id=session_id,
                        features=self.features,
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
                if self.hippocampus and self.hippunfold_directory:
                    print(f"  Processing hippocampal data for features: {', '.join(self.features)}")
                    
                    apply_hippocampal_processing(
                        participant_id=participant_id,
                        session_id=session_id,
                        features=self.features,
                        output_directory=output_directory,
                        workbench_path=env.connectome_workbench_path,
                        micapipe_directory=self.micapipe_directory,
                        hippunfold_directory=self.hippunfold_directory,
                        tmp_dir=session_tmp_dir,
                        smoothing_fwhm=hippocampal_smoothing,
                        verbose=verbose
                    )

                # If subcortex is enabled, extract subcortical stats
                if self.subcortical and self.freesurfer_directory:
                    print(f"  Processing subcortical data for features: {', '.join(self.features)}")
                    
                    apply_subcortical_processing(
                        participant_id=participant_id,
                        session_id=session_id,
                        features=self.features,
                        output_directory=output_directory,
                        micapipe_directory=self.micapipe_directory,
                        freesurfer_directory=self.freesurfer_directory,
                        verbose=verbose
                    )

            finally:
                # Clean up session-specific tmp directory
                if os.path.exists(session_tmp_dir):
                    shutil.rmtree(session_tmp_dir)
                    if verbose:
                        print(f"  Cleaned up temporary directory: {session_tmp_dir}")
    
        print(f"Dataset {self.name} processed successfully.")
        return self
