import os
import sys
import subprocess
import tempfile
from pathlib import Path
from functions.utilities import assert_same_size, assert_required,  show_error, assert_exists,show_info, show_note,do_cmd, show_title
import shutil
import glob
import time
import datetime
import atexit
from functions.new_constants import *
from functions import run_proc, run_analysis
import colorama
from functions.environment import setenv
from joblib import Parallel, delayed
import copy
def main(args):

    # Define the ZBRAINS and script_dir variables
    ZBRAINS = Path(os.path.realpath(__file__)).parent
    script_dir = ZBRAINS / 'functions'
    if isinstance(args.column_map, list):
        args.column_map = dict([arg.split('=') for arg in args.column_map]) if args.column_map else None

    # Check for unknown arguments
    unknown_args = vars(args)
    for arg in sys.argv:
        if arg.startswith('--'):
            arg = arg[2:]
            if arg not in unknown_args:
                print(f"Unknown option '--{arg}'")
                sys.exit(1)


    # ------------------------------------------------- Check arguments ------------------------------------------------- #
    VERBOSE=args.verbose
    tasks = args.run
    if "all" in tasks:
        tasks = LIST_TASKS

    # Check options required for processing
    if "proc" in tasks:
        assert_required("--micapipe", args.micapipe, "--micapipe option is required for post-processing.")
        if "hippocampus" in args.struct:
            assert_required("--hippunfold", args.hippunfold, "--hippunfold option is required for post-processing.")

    # Check options required for regional/asymmetry analysis
    if "analysis" in tasks:
        assert_required("--dataset_ref", args.dataset_ref, "--dataset_ref option is required for analysis.")
        assert_required("--zbrains_ref", args.zbrains_ref, "--zbrains_ref option is required for analysis.")
        assert_required("--demo_ref", args.demo_ref, "--demo_ref option is required for analysis.")
        assert_same_size("--dataset_ref", args.dataset_ref, "--demo_ref", args.demo_ref)

        n_datasets = len(args.dataset_ref)
        if n_datasets != len(args.zbrains_ref):
            if len(args.zbrains_ref) > 1:
                assert_same_size("--dataset_ref", args.dataset_ref, "--zbrains_ref", args.zbrains_ref)
            else:
                args.zbrains_ref = [args.zbrains_ref[0]] * n_datasets

        if args.normative or args.deconfound:
            assert_required("--demo", args.demo, "--demo option is required for analysis.")

    # Some defaults
    sid = f"sub-{args.sub.replace('sub-', '')}"
    ses = f"ses-{args.ses.replace('ses-', '')}" if args.ses else None

    structures = args.struct or ["all"]
    if "all" in structures:
        structures = LIST_STRUCTURES
    structures.sort(key=str.lower)  # case-insensitive sort

    features = args.feat or ["all"]
    if "all" in features:
        features = LIST_FEATURES
    features.sort(key=str.lower)  # case-insensitive sort

    resolutions = args.resolution or ["all"]
    if "all" in resolutions:
        resolutions = LIST_RESOLUTIONS

    labels_ctx = args.label_ctx or ["midthickness"]
    labels_hip = args.label_hip or ["midthickness"]

    smooth_ctx = args.smooth_ctx or DEFAULT_SMOOTH_CTX
    smooth_hip = args.smooth_hip or DEFAULT_SMOOTH_HIP
    threshold = args.threshold or DEFAULT_THRESHOLD



    # ------------------------------------------------ Check dependencies ----------------------------------------------- #
    if "proc" in tasks:
        
        
        if os.environ['WORKBENCH_PATH'] is not None:
            binary = os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command')
        else:
            binary = shutil.which('wb_command')
        if binary is None:
            show_error("Workbench not found. Please set the WORKBENCH_PATH environment variable to the location of Workbench binaries.")
            sys.exit(1)
        os.environ['WORKBENCH_PATH'] = os.path.dirname(binary)



    # -------------------------------------------- Check files & directories -------------------------------------------- #
    dataset_path = os.path.realpath(os.path.join(args.dataset, "derivatives"))
    assert_exists(dataset_path)
    px_zbrains_path = os.path.join(dataset_path, args.zbrains)

    # Export BIDS_ID, SUBJECT_OUTPUT_DIR
    BIDS_ID = f"{sid}_{ses}" if ses else sid
    SUBJECT_OUTPUT_DIR = os.path.join(px_zbrains_path, sid, ses) if ses else os.path.join(px_zbrains_path, sid)
    # Sanity checks required for processing
    if "proc" in tasks:

        # Subject's micapipe directory exists
        SUBJECT_MICAPIPE_DIR = os.path.join(dataset_path, args.micapipe, sid, ses) if ses else os.path.join(dataset_path, args.micapipe, sid)
        assert_exists(SUBJECT_MICAPIPE_DIR, f"{BIDS_ID} micapipe directory does not exist.")
        
        # Subject's hippunfold directory exists
        if "hippocampus" in structures:
            SUBJECT_HIPPUNFOLD_DIR = os.path.join(dataset_path, args.hippunfold, "hippunfold", sid, ses) if ses else os.path.join(dataset_path, args.hippunfold, "hippunfold", sid)
            assert_exists(SUBJECT_HIPPUNFOLD_DIR, f"{BIDS_ID} hippunfold directory does not exist.")
        if args.plugin:
            SUBJECT_PLUGIN_DIR = os.path.join(dataset_path, args.plugin, sid, ses or "")
            if not os.path.exists(SUBJECT_PLUGIN_DIR):
                sys.exit(f"{BIDS_ID} plugin directory does not exist.")
        else:
            SUBJECT_PLUGIN_DIR = None
        # Check if subject's freesurfer/fastsurfer directory exists - only needed for subcortex
        if "subcortex" in structures:
            # Set surface directory and check if subject has a surface directory
            subject_micapipe_qc = os.path.join(SUBJECT_MICAPIPE_DIR, "QC")
            Nrecon = len(glob.glob(f"{subject_micapipe_qc}/{BIDS_ID}_module-proc_surf-*.json"))
            if Nrecon < 1:
                show_error(f"{BIDS_ID} doesn't have a module-proc_surf: run -proc_surf")
                sys.exit(1)
            elif Nrecon == 1:
                module_qc = glob.glob(f"{subject_micapipe_qc}/{BIDS_ID}_module-proc_surf-*.json")[0]
                recon = os.path.splitext(module_qc)[0].split('proc_surf-')[1]
            elif Nrecon > 1:
                show_error(f"{BIDS_ID} has been processed with freesurfer and fastsurfer. Not supported yet")
                sys.exit(1)

            # recon is 'freesurfer' or 'fastsurfer'
            SUBJECT_SURF_DIR = os.path.join(dataset_path, recon, BIDS_ID)
            assert_exists(SUBJECT_SURF_DIR, f"{BIDS_ID} {recon} directory does not exist.")



    # ------------------------------------------------------------------------------------------------------------------- #
    # -------------------------------------------- Start processing/analysis -------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------- #

    show_info("zbrains is running with:")
    if "proc" in tasks:

        # Get WorkBench version
        workbench_version = subprocess.check_output([os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command'), '-version']).decode().split()[1]
        show_note("WorkBench...", workbench_version)
        show_note("            ", os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command'))

    if "analysis" in tasks or "subcortex" in structures:
        # Get Python version
        python_version = subprocess.check_output(['python', '--version']).decode().split()[1]
        show_note("python......", python_version)
        show_note("            ", shutil.which('python'))

    show_note("", "")

    # -------------------------------------------------- Initiate timer ------------------------------------------------- #
    start_time = time.time()

    # -------------------------------------------- Create zbrains directory --------------------------------------------- #
    # Create subject output directory structure
    if not os.path.isdir(SUBJECT_OUTPUT_DIR):
        show_info(f"Subject {BIDS_ID} directory doesn't exist, creating...")

    # Folder for logs
    logs_dir = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS)
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir,exist_ok=True)
    # do_cmd(['mkdir', '-m', '770', '-p', logs_dir])

    if "proc" in tasks:
        dirs = [FOLDER_SCTX, FOLDER_CTX, FOLDER_HIP]
        for dir in dirs:
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_MAPS, dir)
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)

    # Folders for regional/asymmetry analysis
    if "analysis" in tasks:
        dirs = [FOLDER_SCTX, FOLDER_CTX, FOLDER_HIP]
        for dir in dirs:
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_NORM_Z, dir)
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_NORM_MODEL, dir)
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)

    # Temporary folder
    tmp_dir = tempfile.mkdtemp(dir=SUBJECT_OUTPUT_DIR, prefix="z_brains_temp.")
    try:
        os.chmod(SUBJECT_OUTPUT_DIR, 0o770)

        # ----------------------------------------------------- Cleanup ----------------------------------------------------- #
        def cleanup(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)

        atexit.register(cleanup, tmp_dir)

        # -------------------------------------------------- Postprocessing ------------------------------------------------- #
        if "proc" in tasks:
            logfile = os.path.join(logs_dir, f"proc_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt")
            with open(logfile, 'w') as f:
                f.write("\n")
            
            for struct in structures:
                # cmd = [os.path.join(script_dir, 'run_proc.sh'), '--struct', struct, '--feat', ' '.join(features), '--tmp', tmp_dir]

                if struct == "hippocampus":
                    list_resolutions = [map_resolution_hip[k] for k in resolutions]
                    temp_labels = labels_hip
                    fwhm = smooth_hip
                elif struct == "cortex":
                    list_resolutions = [map_resolution_ctx[k] for k in resolutions]
                    temp_labels = labels_ctx
                    fwhm = smooth_ctx
                    # list_resolutions = [map_resolution_ctx[k] for k in resolutions]
                    # cmd.extend(['--fwhm', smooth_ctx, '--resolution', ' '.join(list_resolutions), '--labels', ' '.join(labels_ctx)])
                else:
                    list_resolutions = None
                    temp_labels = None
                    fwhm = None
                run_proc.run(structure=struct, features=features, tmp_dir=tmp_dir, WORKBENCH_PATH=os.environ['WORKBENCH_PATH'], subject_micapipe_dir=SUBJECT_MICAPIPE_DIR, subject_output_dir=SUBJECT_OUTPUT_DIR, folder_maps=FOLDER_MAPS, folder_ctx=FOLDER_CTX, folder_sctx=FOLDER_SCTX, folder_hip=FOLDER_HIP, subject_surf_dir=SUBJECT_SURF_DIR, subject_hippunfold_dir=SUBJECT_HIPPUNFOLD_DIR, script_dir=script_dir, BIDS_ID=BIDS_ID, VERBOSE=VERBOSE, fwhm=fwhm, resolutions=list_resolutions, labels=temp_labels,subject_plugin_dir=SUBJECT_PLUGIN_DIR)
                # do_cmd(cmd)

        # ----------------------------------------------------- Analysis ---------------------------------------------------- #
        if "analysis" in tasks:
            logfile = os.path.join(logs_dir, f"analysis_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt")

            ref_zbrains_paths = [os.path.join(args.dataset_ref[idx], 'derivatives', args.zbrains_ref[idx]) for idx in range(len(args.dataset_ref))]
            
            run_analysis.run(
            subject_id=sid,
            session=ses,
            zbrains_ref=ref_zbrains_paths,
            demo_ref=args.demo_ref,
            zbrains=px_zbrains_path,
            struct=structures,
            feat=features,
            smooth_ctx=smooth_ctx,
            smooth_hip=smooth_hip,
            threshold=threshold,
            approach='zscore',
            resolution=resolutions,
            labels_ctx=labels_ctx,
            labels_hip=labels_hip,
            logfile=logfile,
            tmp=tmp_dir,
            verbose=VERBOSE,
            demo=args.demo,
            normative=args.normative, 
            deconfound=args.deconfound, 
            column_map=args.column_map ,
        )
            # print(' '.join(cmd))
            # do_cmd(cmd)

        # Total running time
        elapsed = (time.time() - start_time) / 60
        show_title(f"Total elapsed time for {BIDS_ID}: {elapsed:.2f} minutes")
        try:
            cleanup(tmp_dir)
        except:
            pass
    except Exception as e:
        print(e)
        try:
            cleanup(tmp_dir)
        except:
            pass
if __name__ == '__main__':
    pcolor = colorama.Fore.MAGENTA  # Purple
    rcolor = colorama.Fore.RED  # Red
    gcolor = colorama.Fore.GREEN  # Green
    bcolor = colorama.Fore.BLUE  # Blue
    gray = colorama.Fore.LIGHTBLACK_EX  # Gray
    nc = colorama.Style.RESET_ALL  # No color
    help = f"""
    {pcolor}COMMAND:{nc}
    zbrains.py


    {pcolor}OPTIONS:{nc}
    \t{rcolor}--sub{nc} ID                  : Subject ID. This is the target subject. Example: 'sub-PX001'.
    \t{rcolor}--dataset{nc} path            : Path to the BIDS dataset containing the target subject's data.
                                        Example: '/path/to/BIDSDataset'.
    \t{rcolor}--zbrains{nc} dir             : Name of the zbrains derivative folder in the target BIDS dataset. The
                                        folder will be created if it does not exist. Example: '--zbrains zbrainsdir' for
                                        '/path/to/BIDSDataset/derivatives/zbrainsdir'.
    \t{gcolor}--run{nc} task                : Tasks to perform. Options:
                                        - {bcolor}proc{nc}          : perform post-processing of target subject (default).
                                        - analysis      : perform analysis (regional & asymmetry) and generate clinical
                                                            report for target subject. The analysis is performed wrt the
                                                            reference subjects (those provided with --demo_ref).
                                                            Post-processing of both target and reference subjects must be
                                                            performed beforehand.
                                        - all           : perform all tasks

    \t{gcolor}--ses{nc} [ses]               : Identifier for a session in the target subject data. If omitted, data will
                                        be managed as a single session. Example: 'ses-001'.
    \t{gcolor}--micapipe{nc} [dir]          : Name of the micapipe derivative folder in the target BIDS dataset. Required
                                        only for post-processing. Example: '--micapipe micapipedir' for
                                        '/path/to/BIDSDataset/derivatives/micapipedir'.
    \t{gcolor}--hippunfold{nc} [dir]        : Name of the hippunfold derivative folder in the target BIDS dataset. Required
                                        only for post-processing. Example: '--hippunfold hipdir' for
                                        '/path/to/BIDSDataset/derivatives/hipdir'.
    \t${gcolor}--plugin${nc} [dir]            : Name of a plugin derivative folder in the target BIDS dataset. zbrains can accept
                                    data outside of micapipe and hippunfold as a 'plugin' folder. However, these data MUST
                                    be formatted as BIDS-derivatives exactly as in micapipe and hippunfold. If hippocampal
                                    surface data are present then they will be used but otherwise volumetric data will be
                                    mapped to hippocampal and subcortical surfaces. 
                                    '/path/to/BIDSDataset/derivatives/plugindir'.
    \t{gcolor}--demo{nc} [path]             : CSV/TSV file with demographics for the target subject. Required only for
                                        analysis when provided by --normative or --deconfound. Additionally, the file is
                                        also used to extract the subject's age and sex, if available, to be included in the
                                        clinical report. The file must include one row with target subject ID and session.
                                        Expected column names:
                                        - participant_id: Subject ID. Required.
                                        - session_id    : Session ID. Use 'n/a' if no session available.
                                        - age           : Subject age. Required only when used by --normative or
                                                            --deconfound. If provided, it will also be included in the
                                                            clinical report.
                                        - sex           : Subject sex. Possible values: 'F' or 'M'. Required only when
                                                            used by --normative or --deconfound. If provided, it will
                                                            also be included in the clinical report.
                                        - site          : Acquisition site. Required only when used by --normative or
                                                            --deconfound.
                                        - Other         : Any other columns used by --normative or --deconfound.
                                        Use the --column_map option to indicate if these variables are under different
                                        column names in your file.

    \t{gcolor}--dataset_ref{nc} [path ...]  : Paths to the BIDS datasets containing the reference subjects data. Required
                                        only for analysis. Each dataset must correspond to one file in --demo_ref.
    \t{gcolor}--zbrains_ref{nc} [dir ...]   : Names of the zbrains derivative folder in each of the reference datasets.
                                        Required only for analysis. If only one folder name is provided but there are
                                        multiple datasets, we assume the same name in all datasets.
    \t{gcolor}--demo_ref{nc} [path ...]     : CSV/TSV files with demographics for reference subjects. Required only for
                                        analysis. There must be one file for each reference dataset (--dataset_ref).
                                        Required only for analysis. See --demo for expected column names.

    \t{gcolor}--struct{nc} [structure ...]  : Structures to use in processing and/or analysis. Options:
                                        - cortex        : cortical data
                                        - subcortex     : subcortical data
                                        - hippocampus   : hippocampal data
                                        - {bcolor}all{nc}           : all structures (default)
    \t{gcolor}--feat{nc} [feature ...]      : Features to use in processing and/or analysis. Options:
                                        - ADC           : apparent diffusion coefficient
                                        - FA            : fractional anisotropy
                                        - flair         : FLAIR
                                        - qT1           : quantitative T1
                                        - thickness     : cortical thickness (for subcortex, volume is used)
                                        - {bcolor}all{nc}           : all features (default)
                                        - plugin-*      : when pulling data from a plugin, feature names must be given the 
                                                        'plugin-' prefix (but this is not needed in the actual file name)
    \t{gcolor}--normative{nc} [cov ...]     : Normative modeling based on provided covariates. Covariates must match
                                        columns in --demo and --demo_ref files. Note that --normative expects some
                                        covariates to have specific names (see --column_map).
                                        Example: '--normative site age sex'.
    \t{gcolor}--deconfound{nc} [[-]cov ...] : Deconfounding based on provided covariates. Covariates must match columns in
                                        --demo and --demo_ref CSV files. If the covariates include 'site', deconfounding is
                                        performed using ComBat. Otherwise, linear regression is used to regress out the
                                        effects of the covariates from the data. By default, ComBat preserves additional
                                        covariates. To remove the effects of a covariate, prepend with '-' (ignored when not
                                        using ComBat). Note that --deconfound expects some covariates to have specific names
                                        (see --column_map). Example: '--deconfound site -age -sex group' to harmonize data
                                        while preserving the effects of group and removing those of age and sex.
    \t{gcolor}--resolution{nc} [res ...]    : Surface resolutions to use for cortex and hippocampus. Options:
                                        - low           : 5k cortical & 2mm hippocampal surfaces
                                        - high          : 32k cortical surfaces & 0p5mm hippocampal surfaces
                                        - {bcolor}all{nc}           : all resolutions (default)
    \t{gcolor}--label_ctx{nc} [label]       : Cortical surfaces used in the volume to surface mapping. Options:
                                        - {bcolor}white{nc}         : WM surface (default)
                                        - midthickness  : Midthickness surface
                                        - pial          : Pial surface
                                        - swmD          : Superficial white matter, where D indicates the distance in
                                                            millimeters. Example: --label_ctx swm2
    \t{gcolor}--label_hip{nc} [label]       : Hippocampal surface used in the volume to surface mapping. Options:
                                        - {bcolor}midthickness{nc}  : Midthickness surface (default)
    \t{gcolor}--smooth_ctx{nc} [size]       : Size of gaussian smoothing kernel in mm used for cortical features.
                                        Default is {bcolor}5{nc}.
    \t{gcolor}--smooth_hip{nc} [size]       : Size of gaussian smoothing kernel in mm used for hippocampal features.
                                        Default is {bcolor}2{nc}.
    \t{gcolor}--threshold{nc} [th]          : Threshold for statistical maps used for clinical reports.
                                        Default is {bcolor}1.96{nc}.
    \t{gcolor}--column_map{nc} [VAR=col ...]: Map expected to actual column names in the CSV/TSV files:
                                        - participant_id: Subject ID is assumed to be provided by the 'participant_id'
                                                            column, unless indicated otherwise. For example, if subject ID
                                                            is under the column ‘SubID’, you can indicate this with
                                                            --column_map participant_id=SubID.
                                        - session_id    : Session ID is assumed to be provided by the ‘session_id’ column,
                                                            unless indicated otherwise (e.g., --column_map session_id=ses)
                                        - age           : Age is assumed to be provided by the ‘age’ column, unless
                                                            indicated otherwise (e.g., --column_map age=\"Subject age\")
                                        - sex           : Sex is assumed to be provided by the 'sex' column, unless
                                                            indicated otherwise (e.g., --column_map ses=\"Subject sex\")
                                        - site          : Acquisition site is assumed to be provided by the ‘site’ column,
                                                            unless indicated otherwise (e.g., --column_map site=center)
    \t{gcolor}--init{nc} [path]             : Initialization script that will be sourced before executing the main script.
                                        Useful for setting up environment variables, activating virtual environments, or any
                                        other setup tasks needed before running your script (see DEPENDENCIES below).
    \t{gcolor}--verbose{nc} [level]         : Verbosity level (default is {bcolor}-1{nc}). Levels:
                                        - 0             : Only errors
                                        - 1             : Warning messages and previous levels
                                        - 2             : Information messages and previous levels
                                        - 3             : Command logs and previous levels
                                        - >3 or <0      : All messages
    \t{gcolor}--help{nc}                    : Print help
    \t{gcolor}--version{nc}                 : Print software version



    {pcolor}USAGE:{nc}
        {gray}# Post-processing{nc}
        {pcolor}zbrains.py{nc} {gcolor}--run{nc} proc
                {rcolor}--sub{nc} <participant_id>
                {gcolor}--ses{nc} <session_id>
                {rcolor}--dataset{nc} <path_bids_dataset>
                {rcolor}--zbrains{nc} <zbrains_directory>
                {gcolor}--micapipe{nc} <micapipe_directory>
                {gcolor}--hippunfold{nc} <hipunfold_directory>

        {gray}# Analysis{nc}
        {pcolor}zbrains.py{nc} {gcolor}--run{nc} analysis
                {rcolor}--sub{nc} <participant_id>
                {gcolor}--ses{nc} <session_id>
                {rcolor}--dataset{nc} <participant_dataset>
                {rcolor}--zbrains{nc} <participant_zbrains_dir>
                {gcolor}--demo_ref{nc} <reference_subjects1.csv> <reference_subjects2.csv>
                {gcolor}--dataset_ref{nc} <reference_dataset1> <reference_dataset2>
                {gcolor}--zbrains_ref{nc} <reference_zbrains_dir1> <reference_zbrains_dir2>



    {pcolor}DEPENDENCIES:{nc}
        > workbench   1.4.2   (https://www.humanconnectome.org/software/workbench-command)
        > ANTs        2.3.4   (https://github.com/ANTsX/ANTs)
        > python      3.10    (https://www.python.org)

        To customize binary locations, use the following environment variables:
        - Set ANTSPATH for ANTs
        - Set WORKBENCH_PATH for Workbench

        Control the number of threads:
        - Set ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS for ANTs
        - Set OMP_NUM_THREADS for Workbench

        Example:
        {gray}# Set threads for ANTs{nc}
        $ export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4

        {gray}# Set threads for Workbench{nc}
        $ export OMP_NUM_THREADS=4

        Use the --init option to specify an initialization script for environment variables, Python activation, etc.



    McGill University, MNI, MICA lab, April 2023
    https://github.com/MICA-MNI/micapipe
    https://github.com/MICA-MNI/z-brains
    http://mica-mni.github.io/
    """

    import argparse
    import sys

    # Define the argument parser
    class Parser(argparse.ArgumentParser):

        def print_help(self):
            print(help)
            
    parser = Parser(description='Handle all arguments and create variables')

    # Add arguments to the parser
    parser.add_argument('--sub', type=str,required=True)
    parser.add_argument('--dataset', type=str,required=True)
    parser.add_argument('--zbrains', type=str,required=True)
    parser.add_argument('--ses', type=str, default=None)
    parser.add_argument('--micapipe', type=str, default=None)
    parser.add_argument('--hippunfold', type=str, default=None)
    parser.add_argument('--plugin', type=str, default=None)
    parser.add_argument('--dataset_ref', nargs='*', default=None)
    parser.add_argument('--zbrains_ref', nargs='*', default=None)
    parser.add_argument('--run', nargs='*', default='all')
    parser.add_argument('--struct', nargs='*', default='all')
    parser.add_argument('--feat', nargs='*', default='all')
    parser.add_argument('--demo_ref', nargs='*', default=None)
    parser.add_argument('--demo', type=str, default=None)
    parser.add_argument('--normative', nargs='*', default=None)
    parser.add_argument('--deconfound', nargs='*', default=None)
    parser.add_argument('--resolution', nargs='*', default='all')
    parser.add_argument('--label_ctx', type=str, default=None)
    parser.add_argument('--label_hip', type=str, default=None)
    parser.add_argument('--smooth_ctx', type=str, default=None)
    parser.add_argument('--smooth_hip', type=str, default=None)
    parser.add_argument('--threshold', type=str, default=None)
    parser.add_argument('--column_map', nargs='*', default=None)
    parser.add_argument('--init', type=str, default=None)
    parser.add_argument('--n_jobs', type=int, default=1)
    parser.add_argument('--wb_path', type=str, default="/data/mica1/01_programs/workbench-1.4.2/bin_linux64")
    parser.add_argument('--verbose', type=int, default=-1)
    parser.add_argument('--version', action='version', version='1.0.0')
  
    

    # Parse the arguments
    args = parser.parse_args()
    WORKBENCH_PATH = setenv(args.wb_path)
    os.environ['WORKBENCH_PATH'] = WORKBENCH_PATH
    runs = args.run
    
    procjobs = []
    analysisjobs = []
    
    def check_sub(args, sub, ses=None):
        micapipe_path = os.path.join(args.dataset, 'derivatives', args.micapipe, sub)
        hippunfold_path = os.path.join(args.dataset, 'derivatives', args.hippunfold, 'hippunfold', sub)

        if ses is not None:
            micapipe_path = os.path.join(micapipe_path, ses)
            hippunfold_path = os.path.join(hippunfold_path, ses)

        if os.path.exists(micapipe_path):
            if os.path.exists(hippunfold_path):
                return True
            else:
                if ses is not None:
                    print(f'No hippunfold for {sub}-{ses}, skipping')
                else:
                    print(f'No hippunfold for {sub}, skipping')
                return False
        else:
            if ses is not None:
                print(f'No micapipe for {sub}-{ses}, skipping')
            else:
                print(f'No micapipe for {sub}, skipping')
            return False
        
    if args.sub == 'all':
        subs = ''
        for sub in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe)):
            if check_sub(args,sub):
                subs += f'{sub} '
        args.sub = subs


    if len(args.sub.split(' ')) > 1:
        subs = args.sub.split(' ')
        if args.ses:
            print(f'Will run selected subs using only session {args.ses}')

            if 'proc' in runs:
                for sub in subs:
                    job = copy.copy(args)
                    job.sub, job.ses, job.run = sub, args.ses, 'proc'
                    if check_sub(args,sub,args.ses):
                        procjobs.append(job)

            if 'analysis' in runs:  
                for sub in subs:
                    job = copy.copy(args)
                    job.sub, job.ses, job.run = sub, args.ses, 'analysis'
                    if check_sub(args,sub,args.ses):
                        analysisjobs.append(job)
        else:
            print('Will run selected subs and all sessions')

            if 'proc' in runs: 
                for sub in subs:
                    
                    for ses in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe, sub)):
                        job = copy.copy(args)
                        job.sub, job.ses, job.run = sub, ses, 'proc'
                        if check_sub(args,sub,ses):
                            procjobs.append(job)

            if 'analysis' in runs: 
                for sub in subs:
                    for ses in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe, sub)):
                        job = copy.copy(args)
                        job.sub, job.ses, job.run = sub, ses, 'analysis'
                        if check_sub(args,sub,ses):
                            analysisjobs.append(job)
    else:
        if args.ses:
            print(f'Will run selected subs using only session {args.ses}')

            if 'proc' in runs:
                job = copy.copy(args)
                job.ses, job.run = args.ses, 'proc'
                if check_sub(args,args.sub,args.ses):
                    procjobs.append(job)
                
            if 'analysis' in runs:  
                job = copy.copy(args)
                job.ses, job.run = args.ses, 'analysis'
                if check_sub(args,args.sub,args.ses):
                    analysisjobs.append(job)
        else:
            print('Will run selected subs and all sessions')

            if 'proc' in runs: 
                for ses in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe, args.sub)):
                    job = copy.copy(args)
                    job.ses, job.run = ses, 'proc'
                    if check_sub(args,args.sub,ses):
                        procjobs.append(job)
            if 'analysis' in runs: 
                
                for ses in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe, args.sub)):
                    job = copy.copy(args)
                    job.ses, job.run = ses, 'analysis'
                    if check_sub(args,args.sub,ses):
                        analysisjobs.append(job)

    for job in procjobs:
        print(job.sub, job.ses)
    Parallel(n_jobs=args.n_jobs)(delayed(main)(job) for job in procjobs)
    Parallel(n_jobs=args.n_jobs)(delayed(main)(job) for job in analysisjobs)


