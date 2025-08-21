#!/usr/bin/env python3
"""
Z-Brains Neuroimaging Pipeline

This module provides the core functionality of the Z-Brains neuroimaging
processing and analysis pipeline. It handles the complete workflow from loading
demographic data through processing control and patient datasets, performing
statistical analysis, and generating clinical reports.
"""


import os
import sys
import logging
from typing import List
from zbrains.dataset import zbdataset, demographics
from zbrains.environment import zbenv


def setup_logger(log_file: str = None) -> logging.Logger:
    """
    Set up logging configuration.
    
    Parameters
    ----------
    log_file : str, optional
        Path to log file. If None, no file logging is performed.
        
    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    logger = logging.getLogger('zbrains')
    logger.setLevel(logging.INFO)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Always add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Add file handler only if log_file is provided
    if log_file:
        # Create log directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        logger.info(f"Logging to file: {log_file}")
    
    return logger


def run(
    features: List[str], 
    wb_path: str = "wb_command", 
    num_threads: int = 4, 
    num_threads_wb: int = 2, 
    control_demo_path: str = None, 
    patient_demo_path: str = None,
    normative_columns: List[str] = None, 
    normative_dtypes: List[str] = None,
    micapipe_dir: str = None, 
    hippunfold_dir: str = None, 
    freesurfer_dir: str = None,
    cortex: bool = True, 
    hippocampus: bool = True, 
    subcortical: bool = True,
    output_directory: str = "./zbrains_output",
    cortical_smoothing: float = 10, 
    hippocampal_smoothing: float = 5, 
    verbose: bool = True, 
    force_reprocess: bool = False, 
    method: str = 'wscore',
    log_file: str = None
) -> None:
    """
    Run the complete Z-Brains neuroimaging processing pipeline.
    
    This function orchestrates the entire workflow from preprocessing through
    analysis and reporting. It handles control and patient datasets separately,
    with controls providing the normative model for patient analysis.
    
    Parameters
    ----------
    features : List[str]
        List of neuroimaging features/modalities to process and analyze.
        Common examples: "FA", "ADC", "thickness", "qT1", "FLAIR"
        
    wb_path : str, optional
        Path to Connectome Workbench command (wb_command).
        Default is "wb_command" (assumes it's in PATH)
        
    num_threads : int, optional
        Number of parallel processing threads.
        Default is 4
        
    num_threads_wb : int, optional
        Number of threads for Workbench operations.
        Default is 2
        
    control_demo_path : str, optional
        Path to CSV file with control/healthy subject demographics.
        Required columns depend on normative_columns.
        Must include subject IDs matching derivatives.
        
    patient_demo_path : str, optional
        Path to CSV file with patient demographics.
        Format should match control_demo_path.
        
    normative_columns : List[str], optional
        Column names in demographics files used for normative modeling.
        Default is ["AGE", "SEX"]
        
    normative_dtypes : List[str], optional
        Data types corresponding to normative_columns.
        Options include "int", "float", "binary", "categorical"
        Default is ["int", "binary"]
        
    micapipe_dir : str, optional
        Path to micapipe derivatives directory (BIDS derivatives structure).
        Required if using cortical or subcortical features.
        
    hippunfold_dir : str, optional
        Path to hippunfold derivatives directory.
        Required if using hippocampal features.
        
    freesurfer_dir : str, optional
        Path to FreeSurfer derivatives directory.
        Required for certain features like cortical thickness.
        
    cortex : bool, optional
        Whether to process cortical structures.
        Default is True
        
    hippocampus : bool, optional
        Whether to process hippocampal structures. 
        Default is True
        
    subcortical : bool, optional
        Whether to process subcortical structures.
        Default is True
        
    output_directory : str, optional
        Directory where processed data and reports will be saved.
        Default is "./zbrains_output"
        
    cortical_smoothing : float, optional
        Smoothing kernel size in mm for cortical surfaces.
        Default is 10 mm
        
    hippocampal_smoothing : float, optional
        Smoothing kernel size in mm for hippocampal surfaces.
        Default is 5 mm
        
    verbose : bool, optional
        Print detailed progress information.
        Default is True
        
    force_reprocess : bool, optional
        Force reprocessing of data even if it already exists.
        Default is False
        
    method : str, optional
        Statistical method for patient-control comparison.
        Options: 'wscore' (Weighted score) or 'zscore'
        Default is 'wscore'
        
    log_file : str, optional
        Path to log file. If None, logs to console only.
        
    Returns
    -------
    None
        Results are saved to the specified output_directory
        
    Notes
    -----
    - The pipeline requires neuroimaging derivatives from micapipe, hippunfold,
      and/or FreeSurfer to be available in their respective directories.
    - Control processing must complete successfully before patient analysis.
    - The output directory will be created if it doesn't exist.
    - For weighted score method (wscore), age and sex must be included as normative
      variables.
    """
    
    # Setup logging
    logger = setup_logger(log_file)
    logger.info("Starting Z-Brains neuroimaging pipeline")
    
    # Log pipeline configuration
    logger.info(f"Pipeline configuration:")
    logger.info(f"  Features: {features}")
    logger.info(f"  Method: {method}")
    logger.info(f"  Output directory: {output_directory}")
    logger.info(f"  Threads: {num_threads}")
    logger.info(f"  Workbench threads: {num_threads_wb}")
    logger.info(f"  Cortical smoothing: {cortical_smoothing}mm")
    logger.info(f"  Hippocampal smoothing: {hippocampal_smoothing}mm")
    logger.info(f"  Processing structures: cortex={cortex}, hippocampus={hippocampus}, subcortical={subcortical}")
    
    # Set default values if not provided
    if normative_columns is None:
        normative_columns = ["AGE", "SEX"]
    if normative_dtypes is None:
        normative_dtypes = ["int", "binary"]
    
    logger.info(f"  Normative columns: {normative_columns}")
    logger.info(f"  Normative dtypes: {normative_dtypes}")
    
    # Validate required parameters
    if not features:
        logger.error("No features specified for processing")
        raise ValueError("No features specified for processing")
    
    if cortex or subcortical:
        if not micapipe_dir:
            logger.error("micapipe_dir must be provided when processing cortical or subcortical structures")
            raise ValueError("micapipe_dir must be provided when processing cortical or subcortical structures")
    
    if hippocampus and not hippunfold_dir:
        logger.error("hippunfold_dir must be provided when processing hippocampal structures")
        raise ValueError("hippunfold_dir must be provided when processing hippocampal structures")
    
    if not control_demo_path or not patient_demo_path:
        logger.error("Both control_demo_path and patient_demo_path must be provided")
        raise ValueError("Both control_demo_path and patient_demo_path must be provided")
    
    # Create output directory
    os.makedirs(output_directory, exist_ok=True)
    logger.info(f"Output will be saved to: {output_directory}")
    
    # -------------------------------------
    # Initialize environment
    # -------------------------------------
    logger.info(f"Initializing environment with {num_threads} threads...")
    env = zbenv(
        connectome_workbench_path=wb_path, 
        num_threads=num_threads, 
        num_threads_wb=num_threads_wb
    )
    logger.info("Environment initialized successfully")
    
    # -------------------------------------
    # Process control dataset
    # -------------------------------------
    logger.info(f"Loading control demographics from {control_demo_path}")
    
    # Load control demographics
    try:
        control = demographics(
            control_demo_path, 
            normative_columns=normative_columns, 
            normative_dtypes=normative_dtypes
        )
        logger.info(f"Successfully loaded {len(control.data)} control subjects")
    except Exception as e:
        logger.error(f"Failed to load control demographics: {e}")
        raise
        
    # Set up control dataset object
    logger.info("Setting up control dataset")
    control_dataset = zbdataset(
        "controls", 
        demographics=control,
        micapipe_directory=micapipe_dir, 
        hippunfold_directory=hippunfold_dir,
        freesurfer_directory=freesurfer_dir,
        cortex=cortex,
        hippocampus=hippocampus,
        subcortical=subcortical,
    )

    # Process control data if forced or needed
    if force_reprocess:
        logger.info("Force reprocessing control dataset...")
        control_dataset.process(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            env=env, 
            verbose=verbose
        )
        logger.info("Control dataset processing completed")

    # Validate control dataset
    try:
        logger.info("Validating control dataset...")
        control_dataset.validate(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            verbose=verbose
        )
        if control_dataset.valid_dataset:
            logger.info("Control dataset validation passed")
        else:
            logger.warning("Control dataset validation failed")
    except Exception as e:
        logger.error(f"Control dataset validation failed: {e}")
        logger.info("Attempting to reprocess control dataset...")

    # Reprocess if validation failed
    if not control_dataset.valid_dataset:
        logger.info("Processing control dataset after validation failure")
        control_dataset.process(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            env=env, 
            verbose=verbose
        )
        logger.info("Control dataset reprocessing completed")
        
    # Check control dataset validity again
    if not control_dataset.valid_dataset:
        logger.error("Control dataset still invalid after reprocessing")
        raise RuntimeError("Failed to create valid control dataset")
        
    # -------------------------------------
    # Process patient dataset
    # -------------------------------------
    logger.info(f"Loading patient demographics from {patient_demo_path}")
    
    # Load patient demographics (using control as reference for normalization)
    try:
        patient = demographics(
            patient_demo_path, 
            reference=control, 
            normative_columns=normative_columns, 
            normative_dtypes=normative_dtypes
        )
        logger.info(f"Successfully loaded {len(patient.data)} patient subjects")
    except Exception as e:
        logger.error(f"Failed to load patient demographics: {e}")
        raise
    
    # Create patient dataset object
    logger.info("Setting up patient dataset")
    patient_dataset = zbdataset(
        "patients", 
        demographics=patient,
        micapipe_directory=micapipe_dir,
        hippunfold_directory=hippunfold_dir,
        freesurfer_directory=freesurfer_dir,
        cortex=cortex,
        hippocampus=hippocampus,
        subcortical=subcortical,
    )

    # Process patients if forced or needed
    if force_reprocess:
        logger.info("Force reprocessing patient dataset")
        patient_dataset.process(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            env=env, 
            verbose=verbose
        )
        logger.info("Patient dataset processing completed")

    # Validate patient dataset
    try:
        logger.info("Validating patient dataset")
        patient_dataset.validate(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            verbose=verbose
        )
        if patient_dataset.valid_dataset:
            logger.info("Patient dataset validation passed")
        else:
            logger.warning("Patient dataset validation failed")
    except Exception as e:
        logger.error(f"Patient dataset validation failed: {e}")
        logger.info("Attempting to reprocess patient dataset...")

    # Reprocess if validation failed
    if not patient_dataset.valid_dataset:
        logger.info("Processing patient dataset after validation failure")
        patient_dataset.process(
            output_directory=output_directory, 
            features=features, 
            cortical_smoothing=cortical_smoothing, 
            hippocampal_smoothing=hippocampal_smoothing, 
            env=env, 
            verbose=verbose
        )
        logger.info("Patient dataset reprocessing completed")
        
    # Check patient dataset validity again
    if not patient_dataset.valid_dataset:
        logger.error("Patient dataset still invalid after reprocessing")
        raise RuntimeError("Failed to create valid patient dataset")

    # -------------------------------------
    # Analysis and reporting
    # -------------------------------------
    logger.info(f"Analyzing patient data using method: {method}")
    
    # Compare patients to controls using selected method
    patient_dataset.analyze(
        output_directory=output_directory, 
        reference=control_dataset, 
        method=method
    )
    logger.info("Statistical analysis completed")

    # Generate clinical reports
    logger.info("Generating clinical reports")
    patient_dataset.clinical_report(
        output_directory=output_directory, 
        approach=method,
        features=features,
        env=env,
        verbose=verbose
    )
    logger.info("Clinical reports generated")
    
    logger.info(f"Pipeline complete. Results saved to: {output_directory}")
    
    # Log completion summary
    logger.info("Pipeline execution summary:")
    logger.info(f"  Control subjects processed: {len(control.data)}")
    logger.info(f"  Patient subjects processed: {len(patient.data)}")
    logger.info(f"  Features analyzed: {', '.join(features)}")
    logger.info(f"  Analysis method: {method}")
    logger.info(f"  Output location: {output_directory}")