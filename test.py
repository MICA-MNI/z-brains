#!/usr/bin/env python3
"""
Z-Brains Neuroimaging Pipeline

This script provides a user-friendly interface to run the Z-Brains neuroimaging 
processing and analysis pipeline. It handles the complete workflow from loading 
demographic data through processing control and patient datasets, performing 
statistical analysis, and generating clinical reports.

The pipeline performs these main steps:
1. Initializes the processing environment (threads, external tools)
2. Loads and validates control/patient demographic data 
3. Processes neuroimaging features for control subjects
4. Processes neuroimaging features for patient subjects
5. Performs statistical analysis comparing patients to control distribution
6. Generates clinical reports with visualization of abnormalities

Usage examples:
    # As a Python module
    from zbrains.test import run
    run(features=["thickness", "FA"], output_directory="./output")
    
    # As a command-line tool
    python test.py --features FA thickness qT1 --output-dir ./output
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from typing import List, Union, Dict, Optional

from zbrains.dataset import zbdataset, demographics
from zbrains.environment import zbenv


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
    # Set up logging
    logger = logging.getLogger("zbrains")
    logger.setLevel(logging.INFO if verbose else logging.WARNING)
    
    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(console_handler)
    
    # Add file handler if log_file is specified
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Set default values if not provided
    if normative_columns is None:
        normative_columns = ["AGE", "SEX"]
    if normative_dtypes is None:
        normative_dtypes = ["int", "binary"]
    
    # Validate required parameters
    if not features:
        raise ValueError("No features specified for processing")
    
    if cortex or subcortical:
        if not micapipe_dir:
            raise ValueError("micapipe_dir must be provided when processing cortical or subcortical structures")
    
    if hippocampus and not hippunfold_dir:
        raise ValueError("hippunfold_dir must be provided when processing hippocampal structures")
    
    if not control_demo_path or not patient_demo_path:
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
    except Exception as e:
        logger.warning(f"Control dataset validation failed: {e}")
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
    except Exception as e:
        logger.warning(f"Patient dataset validation failed: {e}")
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

    # Generate clinical reports
    logger.info("Generating clinical reports")
    patient_dataset.clinical_report(
        output_directory=output_directory, 
        approach=method,
        features=features,
        env=env,
        verbose=verbose
    )
    
    logger.info(f"Pipeline complete. Results saved to: {output_directory}")


def parse_args():
    """Parse command-line arguments for the Z-Brains pipeline."""
    parser = argparse.ArgumentParser(
        description="Z-Brains neuroimaging processing and analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        "--features", 
        nargs="+", 
        required=True,
        choices=["FA", "ADC", "thickness", "qT1", "FLAIR", "qT1-blur", "FLAIR-blur"],
        help="List of features/modalities to process (e.g. FA ADC thickness qT1 FLAIR)"
    )
    
    # Demographics files
    demo_group = parser.add_argument_group('Demographics files')
    demo_group.add_argument(
        "--control-demographics", 
        dest="control_demo_path",
        required=True,
        help="Path to control demographics CSV file"
    )
    demo_group.add_argument(
        "--patient-demographics", 
        dest="patient_demo_path",
        required=True,
        help="Path to patient demographics CSV file"
    )
    demo_group.add_argument(
        "--normative-columns", 
        nargs="+", 
        default=["AGE", "SEX"],
        help="Column names in demographics for normative modeling"
    )
    demo_group.add_argument(
        "--normative-dtypes", 
        nargs="+", 
        default=["int", "binary"],
        help="Data types for normative columns (e.g., int, binary, float)"
    )
    
    # Input paths
    paths_group = parser.add_argument_group('Input data paths')
    paths_group.add_argument(
        "--micapipe-dir", 
        required=True,
        help="Path to micapipe derivatives directory"
    )
    paths_group.add_argument(
        "--hippunfold-dir", 
        required=True,
        help="Path to hippunfold derivatives directory"
    )
    paths_group.add_argument(
        "--freesurfer-dir", 
        required=True,
        help="Path to FreeSurfer derivatives directory"
    )
    
    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        "--output-dir", 
        dest="output_directory",
        default="./zbrains_output",
        help="Directory where processed data and reports will be saved"
    )
    output_group.add_argument(
        "--log-file",
        help="Path to log file (if not specified, logs to console only)"
    )
    
    # Processing parameters
    proc_group = parser.add_argument_group('Processing parameters')
    proc_group.add_argument(
        "--cortical-smoothing", 
        type=float, 
        default=10.0,
        help="Smoothing kernel size (mm) for cortical surfaces"
    )
    proc_group.add_argument(
        "--hippocampal-smoothing", 
        type=float, 
        default=5.0,
        help="Smoothing kernel size (mm) for hippocampal surfaces"
    )
    proc_group.add_argument(
        "--method", 
        choices=["wscore", "zscore"], 
        default="wscore",
        help="Statistical method for patient-control comparison"
    )
    proc_group.add_argument(
        "--no-cortex", 
        dest="cortex",
        action="store_false", 
        help="Skip cortical processing"
    )
    proc_group.add_argument(
        "--no-hippocampus", 
        dest="hippocampus",
        action="store_false", 
        help="Skip hippocampal processing"
    )
    proc_group.add_argument(
        "--no-subcortical", 
        dest="subcortical",
        action="store_false", 
        help="Skip subcortical processing"
    )
    
    # Environment parameters
    env_group = parser.add_argument_group('Environment parameters')
    env_group.add_argument(
        "--wb-path", 
        default="wb_command",
        help="Path to Connectome Workbench command (wb_command)"
    )
    env_group.add_argument(
        "--threads", 
        dest="num_threads",
        type=int, 
        default=4,
        help="Number of parallel processing threads"
    )
    env_group.add_argument(
        "--wb-threads", 
        dest="num_threads_wb",
        type=int, 
        default=2,
        help="Number of threads for Workbench operations"
    )
    
    # Other options
    other_group = parser.add_argument_group('Other options')
    other_group.add_argument(
        "--force-reprocess", 
        action="store_true",
        help="Force reprocessing of data even if it already exists"
    )
    other_group.add_argument(
        "--quiet", 
        dest="verbose",
        action="store_false", 
        help="Suppress detailed progress information"
    )
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    
    # Convert args to dict and pass to run function
    arg_dict = vars(args)
    run(**arg_dict)