#!/usr/bin/env python3
"""
Z-Brains: Command-line interface

This module provides the command-line interface for the Z-Brains neuroimaging 
processing and analysis pipeline.

Usage:
    zbrains --features FA thickness --control-demographics controls.csv [options]
"""

import argparse
import sys
from zbrains.run import run


class Colors:
    """ANSI color codes for terminal output."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    # Colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    
    # Background colors
    BG_BLUE = '\033[44m'
    BG_GREEN = '\033[42m'
    BG_YELLOW = '\033[43m'
    BG_RED = '\033[41m'


def colorize(text, color_code):
    """Add color to text if stdout is a terminal."""
    if sys.stdout.isatty():
        return f"{color_code}{text}{Colors.RESET}"
    return text


def show_detailed_help():
    """Display detailed help information about the Z-Brains pipeline."""
    help_text = f"""
{colorize('╔═══════════════════════════════════════════════════════════════════════════════╗', Colors.CYAN + Colors.BOLD)}
{colorize('║', Colors.CYAN + Colors.BOLD)}                           {colorize('Z-BRAINS v2.0.0', Colors.WHITE + Colors.BOLD)}                                    {colorize('║', Colors.CYAN + Colors.BOLD)}
{colorize('║', Colors.CYAN + Colors.BOLD)}                {colorize('Neuroimaging Processing and Analysis Pipeline', Colors.BLUE)}                  {colorize('║', Colors.CYAN + Colors.BOLD)}
{colorize('╚═══════════════════════════════════════════════════════════════════════════════╝', Colors.CYAN + Colors.BOLD)}

{colorize('DESCRIPTION:', Colors.YELLOW + Colors.BOLD)}
    Z-Brains is a comprehensive neuroimaging pipeline for processing and analyzing
    structural and diffusion MRI data. It performs patient-control comparisons using
    normative modeling approaches to identify brain abnormalities.

{colorize('WORKFLOW:', Colors.YELLOW + Colors.BOLD)}
    {colorize('1.', Colors.GREEN)} Load control and patient demographics
    {colorize('2.', Colors.GREEN)} Process control dataset to create normative models
    {colorize('3.', Colors.GREEN)} Process patient dataset using the same preprocessing
    {colorize('4.', Colors.GREEN)} Compare patients to controls using statistical methods (zscore/wscore)
    {colorize('5.', Colors.GREEN)} Generate individual clinical reports with brain maps

{colorize('REQUIRED ARGUMENTS:', Colors.RED + Colors.BOLD)}
    {colorize('--features', Colors.CYAN + Colors.BOLD)} FEATURE [FEATURE ...]
        Neuroimaging features to process and analyze.
        
        Available features:
        {colorize('•', Colors.GREEN)} {colorize('FA', Colors.MAGENTA)}              - Fractional Anisotropy (DTI)
        {colorize('•', Colors.GREEN)} {colorize('ADC', Colors.MAGENTA)}             - Apparent Diffusion Coefficient (DTI)
        {colorize('•', Colors.GREEN)} {colorize('thickness', Colors.MAGENTA)}       - Cortical thickness
        {colorize('•', Colors.GREEN)} {colorize('qT1', Colors.MAGENTA)}             - Quantitative T1 relaxometry
        {colorize('•', Colors.GREEN)} {colorize('FLAIR', Colors.MAGENTA)}           - FLAIR intensity
        {colorize('•', Colors.GREEN)} {colorize('qT1*blur', Colors.MAGENTA)}        - qT1 with additional spatial blurring
        {colorize('•', Colors.GREEN)} {colorize('FLAIR*blur', Colors.MAGENTA)}      - FLAIR with additional spatial blurring
        
        {colorize('Example:', Colors.BLUE)} --features FA ADC thickness

    {colorize('--control-demographics', Colors.CYAN + Colors.BOLD)} PATH
        CSV file containing control/healthy subject demographics.
        
        {colorize('Required columns:', Colors.BLUE)}
        {colorize('•', Colors.GREEN)} {colorize('participant_id', Colors.MAGENTA)}  - Subject identifier (e.g., sub-001)
        {colorize('•', Colors.GREEN)} {colorize('session_id', Colors.MAGENTA)}      - Session identifier (e.g., ses-01)
        {colorize('•', Colors.GREEN)} {colorize('AGE', Colors.MAGENTA)}             - Age in years (integer)
        {colorize('•', Colors.GREEN)} {colorize('SEX', Colors.MAGENTA)}             - Sex (binary: M/F, 0/1, etc.)
        
        {colorize('Optional columns (if using custom normative modeling):', Colors.BLUE)}
        {colorize('•', Colors.GREEN)} {colorize('EDUCATION', Colors.MAGENTA)}       - Years of education
        {colorize('•', Colors.GREEN)} {colorize('HANDEDNESS', Colors.MAGENTA)}      - Handedness score
        {colorize('•', Colors.GREEN)} Any other demographic variables

    {colorize('--patient-demographics', Colors.CYAN + Colors.BOLD)} PATH
        CSV file containing patient demographics.
        Must have the same column structure as control demographics.

    {colorize('--micapipe-dir', Colors.CYAN + Colors.BOLD)} PATH
        Path to micapipe derivatives directory (BIDS format).
        Contains preprocessed cortical and subcortical data.

    {colorize('--hippunfold-dir', Colors.CYAN + Colors.BOLD)} PATH
        Path to hippunfold derivatives directory.
        Required for hippocampal analysis.

    {colorize('--freesurfer-dir', Colors.CYAN + Colors.BOLD)} PATH
        Path to FreeSurfer derivatives directory.
        Required for cortical thickness and subcortical volumes.

{colorize('OPTIONAL ARGUMENTS:', Colors.YELLOW + Colors.BOLD)}
    {colorize('Output Options:', Colors.BLUE + Colors.BOLD)}
    {colorize('--output-dir', Colors.CYAN)} PATH           Output directory (default: ./zbrains_output)
    {colorize('--log-file', Colors.CYAN)} PATH             Log file path (default: console output only)

    {colorize('Processing Parameters:', Colors.BLUE + Colors.BOLD)}
    {colorize('--cortical-smoothing', Colors.CYAN)} INT    Cortical smoothing kernel in mm (default: 10)
    {colorize('--hippocampal-smoothing', Colors.CYAN)} INT Hippocampal smoothing kernel in mm (default: 5)
    {colorize('--method', Colors.CYAN)} {{wscore,zscore}}    Statistical comparison method (default: wscore)
    
    {colorize('Structure Selection:', Colors.BLUE + Colors.BOLD)}
    {colorize('--no-cortex', Colors.CYAN)}                 Skip cortical processing
    {colorize('--no-hippocampus', Colors.CYAN)}           Skip hippocampal processing  
    {colorize('--no-subcortical', Colors.CYAN)}           Skip subcortical processing

    {colorize('Normative Modeling:', Colors.BLUE + Colors.BOLD)}
    {colorize('--normative-columns', Colors.CYAN)} COL [COL ...]    
        Demographic columns for normative modeling (default: AGE SEX)
    {colorize('--normative-dtypes', Colors.CYAN)} TYPE [TYPE ...]   
        Data types for normative columns (default: int binary)
        Options: int, float, binary, categorical

    {colorize('Environment:', Colors.BLUE + Colors.BOLD)}
    {colorize('--wb-path', Colors.CYAN)} PATH              Path to wb_command
    {colorize('--threads', Colors.CYAN)} INT               Number of processing threads (default: 4)
    {colorize('--wb-threads', Colors.CYAN)} INT            Workbench threads (default: 2)

    {colorize('Other:', Colors.BLUE + Colors.BOLD)}
    {colorize('--force-reprocess', Colors.CYAN)}           Force reprocessing even if data exists
    {colorize('--quiet', Colors.CYAN)}                     Suppress detailed output
    {colorize('--version', Colors.CYAN)}                   Show version and exit

{colorize('STATISTICAL METHODS:', Colors.YELLOW + Colors.BOLD)}
    {colorize('wscore', Colors.GREEN + Colors.BOLD)} (Weighted Score):
    {colorize('•', Colors.GREEN)} Uses normative modeling with demographic covariates
    {colorize('•', Colors.GREEN)} Accounts for age, sex, and other variables in control distribution
    {colorize('•', Colors.GREEN)} Provides weighted z-scores based on expected values
    {colorize('•', Colors.GREEN)} {colorize('Recommended for clinical applications', Colors.BLUE)}

    {colorize('zscore', Colors.GREEN + Colors.BOLD)} (Z-Score):
    {colorize('•', Colors.GREEN)} Simple standardization using control mean and standard deviation
    {colorize('•', Colors.GREEN)} Does not account for demographic variables
    {colorize('•', Colors.GREEN)} Faster computation but less precise normalization

{colorize('SUPPORTED BRAIN STRUCTURES:', Colors.YELLOW + Colors.BOLD)}
    {colorize('Cortical:', Colors.GREEN + Colors.BOLD)}
    {colorize('•', Colors.GREEN)} Surface-based analysis on fsLR-32k and fsLR-5k meshes
    {colorize('•', Colors.GREEN)} Midthickness and white matter surfaces
    {colorize('•', Colors.GREEN)} Features: thickness, FA, ADC, qT1, FLAIR, blur variants

    {colorize('Hippocampal:', Colors.GREEN + Colors.BOLD)}
    {colorize('•', Colors.GREEN)} Hippocampus-specific surface analysis
    {colorize('•', Colors.GREEN)} Unfolded hippocampal coordinates
    {colorize('•', Colors.GREEN)} Features: FA, ADC, qT1, FLAIR (no blur variants)

    {colorize('Subcortical:', Colors.GREEN + Colors.BOLD)}
    {colorize('•', Colors.GREEN)} Volume-based analysis of deep brain structures
    {colorize('•', Colors.GREEN)} FreeSurfer segmentation-based
    {colorize('•', Colors.GREEN)} Features: volumes, FA, ADC, qT1, FLAIR (no blur variants)

{colorize('EXAMPLE USAGE:', Colors.YELLOW + Colors.BOLD)}
    {colorize('Basic usage:', Colors.BLUE + Colors.BOLD)}
    {colorize('zbrains', Colors.WHITE + Colors.BOLD)} {colorize('--features', Colors.CYAN)} FA thickness \\
            {colorize('--control-demographics', Colors.CYAN)} controls.csv \\
            {colorize('--patient-demographics', Colors.CYAN)} patients.csv \\
            {colorize('--micapipe-dir', Colors.CYAN)} /data/micapipe \\
            {colorize('--hippunfold-dir', Colors.CYAN)} /data/hippunfold \\
            {colorize('--freesurfer-dir', Colors.CYAN)} /data/freesurfer

    {colorize('Advanced usage with custom parameters:', Colors.BLUE + Colors.BOLD)}
    {colorize('zbrains', Colors.WHITE + Colors.BOLD)} {colorize('--features', Colors.CYAN)} FA ADC thickness qT1 FLAIR \\
            {colorize('--control-demographics', Colors.CYAN)} controls.csv \\
            {colorize('--patient-demographics', Colors.CYAN)} patients.csv \\
            {colorize('--micapipe-dir', Colors.CYAN)} /data/micapipe \\
            {colorize('--hippunfold-dir', Colors.CYAN)} /data/hippunfold \\
            {colorize('--freesurfer-dir', Colors.CYAN)} /data/freesurfer \\
            {colorize('--output-dir', Colors.CYAN)} /results/zbrains \\
            {colorize('--cortical-smoothing', Colors.CYAN)} 10 \\
            {colorize('--hippocampal-smoothing', Colors.CYAN)} 5 \\
            {colorize('--method', Colors.CYAN)} wscore \\
            {colorize('--threads', Colors.CYAN)} 16 \\
            {colorize('--normative-columns', Colors.CYAN)} AGE SEX EDUCATION \\
            {colorize('--normative-dtypes', Colors.CYAN)} int binary int

    {colorize('Skip hippocampal processing:', Colors.BLUE + Colors.BOLD)}
    {colorize('zbrains', Colors.WHITE + Colors.BOLD)} {colorize('--features', Colors.CYAN)} thickness FA \\
            {colorize('--control-demographics', Colors.CYAN)} controls.csv \\
            {colorize('--patient-demographics', Colors.CYAN)} patients.csv \\
            {colorize('--micapipe-dir', Colors.CYAN)} /data/micapipe \\
            {colorize('--freesurfer-dir', Colors.CYAN)} /data/freesurfer \\
            {colorize('--no-hippocampus', Colors.CYAN)}

{colorize('FILE REQUIREMENTS:', Colors.YELLOW + Colors.BOLD)}
    {colorize('Demographics CSV format:', Colors.BLUE + Colors.BOLD)}
    {colorize('participant_id,session_id,AGE,SEX', Colors.DIM)}
    {colorize('sub-001,ses-01,25,M', Colors.DIM)}
    {colorize('sub-002,ses-01,30,F', Colors.DIM)}
    {colorize('sub-003,ses-01,28,M', Colors.DIM)}

    {colorize('Directory structure (BIDS derivatives):', Colors.BLUE + Colors.BOLD)}
    micapipe_dir/
    ├── sub-001/
    │   └── ses-01/
    │       ├── anat/
    │       ├── surf/
    │       └── maps/
    
    hippunfold_dir/
    ├── hippunfold/
    │   └── sub-001/
    │       └── ses-01/
    │           └── surf/

{colorize('OUTPUT STRUCTURE:', Colors.YELLOW + Colors.BOLD)}
    zbrains_output/
    ├── sub-001/
    │   └── ses-01/
    │       ├── maps/
    │       │   ├── cortex/     {colorize('# Processed cortical maps', Colors.DIM)}
    │       │   ├── hippocampus/ {colorize('# Processed hippocampal maps', Colors.DIM)}
    │       │   └── subcortical/ {colorize('# Subcortical statistics', Colors.DIM)}
    │       ├── structural/     {colorize('# Surface files', Colors.DIM)}
    │       ├── analysis/       {colorize('# Statistical comparison results', Colors.DIM)}
    │       └── *.pdf          {colorize('# Clinical report', Colors.DIM)}

{colorize('NOTES:', Colors.YELLOW + Colors.BOLD)}
    {colorize('•', Colors.GREEN)} Control dataset must be processed successfully before patient analysis
    {colorize('•', Colors.GREEN)} Missing data is handled gracefully - subjects with incomplete data are skipped
    {colorize('•', Colors.GREEN)} Processing is parallelized across subjects for efficiency
    {colorize('•', Colors.GREEN)} All intermediate files are preserved for quality control
    {colorize('•', Colors.GREEN)} Clinical reports include brain visualizations and statistical summaries

{colorize('For more information, visit:', Colors.BLUE)} {colorize('https://github.com/MICA-MNI/z-brains', Colors.CYAN + Colors.BOLD)}
"""
    print(help_text)


def parse_args():
    """Parse command-line arguments for the Z-Brains pipeline."""
    
    # Check if no arguments provided
    if len(sys.argv) == 1:
        show_detailed_help()
        sys.exit(0)
    
    parser = argparse.ArgumentParser(
        description="Z-Brains neuroimaging processing and analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False  # We'll add custom help
    )
    
    # Add custom help argument
    parser.add_argument(
        '-h', '--help',
        action='store_true',
        help='Show this help message and exit'
    )
    
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        "--features", 
        nargs="+", 
        required=True,
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
        type=int, 
        default=10,
        help="Smoothing kernel size (mm) for cortical surfaces"
    )
    proc_group.add_argument(
        "--hippocampal-smoothing", 
        type=int, 
        default=5,
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
        required=True,
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
    other_group.add_argument(
        "--version", 
        action="version", 
        version="zbrains 2.0.0",
        help="Show program's version number and exit"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Handle custom help
    if args.help:
        show_detailed_help()
        sys.exit(0)
    
    # Validate that normative-columns and normative-dtypes have the same length
    if args.normative_columns and args.normative_dtypes and len(args.normative_columns) != len(args.normative_dtypes):
        parser.error(f"normative-dtypes ({len(args.normative_dtypes)}) must have the same length as normative-columns ({len(args.normative_columns)})")
    
    return args


def main():
    """Main entry point for the Z-Brains CLI."""
    try:
        args = parse_args()
        
        # Convert args to dict and remove 'help' key if present
        arg_dict = vars(args)
        
        # Remove the 'help' argument since it's handled in parse_args()
        if 'help' in arg_dict:
            del arg_dict['help']
        
        run(**arg_dict)
        
        return 0
    except KeyboardInterrupt:
        print("\nOperation canceled by user", file=sys.stderr)
        return 130
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())