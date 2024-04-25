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
from functions.constants import FOLDER_LOGS, FOLDER_MAPS, FOLDER_NORM_Z, FOLDER_NORM_MODEL, FOLDER_SCTX, FOLDER_CTX, FOLDER_HIP, LIST_TASKS, LIST_STRUCTURES, LIST_FEATURES, LIST_RESOLUTIONS, DEFAULT_SMOOTH_CTX, DEFAULT_SMOOTH_HIP, DEFAULT_THRESHOLD, map_resolution_ctx, map_resolution_hip, ProcessingException
from functions import run_proc, run_analysis
from functions.environment import setenv
from joblib import Parallel, delayed
import copy
from contextlib import contextmanager
import argparse
import sys
from functions.help import help
from functions.cleantemps import delete_temp_folders
import subprocess
import gc
@contextmanager
def tempdir(SUBJECT_OUTPUT_DIR, prefix):
    path = tempfile.mkdtemp(dir=SUBJECT_OUTPUT_DIR, prefix=prefix)
    try:
        yield path
    finally:
        print(f"Cleaning up temp dir {path}")
        try:
            shutil.rmtree(path)
        except IOError:
            sys.stderr.write('Failed to clean up temp dir {}'.format(path))
            
def jobloop(args):
    try:
        main_func(args)
    except Exception as e:
        print(e)
    gc.collect()

def parse_args(args):
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
                # print(f"Unknown option '--{arg}'")
                raise ProcessingException(f"Unknown option '--{arg}'")

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

    return args, ZBRAINS, script_dir, VERBOSE, tasks, sid, ses, structures, features, resolutions, labels_ctx, labels_hip, smooth_ctx, smooth_hip, threshold

def check_workbench_dependency(tasks):
    if "proc" in tasks:
        if os.environ['WORKBENCH_PATH'] is not None:
            binary = os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command')
        else:
            binary = shutil.which('wb_command')
        if binary is None:
            show_error("Workbench not found. Please set the WORKBENCH_PATH environment variable to the location of Workbench binaries.")
            raise ProcessingException("Workbench not found. Please set the WORKBENCH_PATH environment variable to the location of Workbench binaries.")
        os.environ['WORKBENCH_PATH'] = os.path.dirname(binary)

def check_files_and_directories(args, tasks, structures, sid, ses):
    dataset_path = os.path.realpath(os.path.join(args.dataset, "derivatives"))
    assert_exists(dataset_path)
    px_zbrains_path = os.path.join(dataset_path, args.zbrains)

    BIDS_ID = f"{sid}_{ses}" if ses else sid
    SUBJECT_OUTPUT_DIR = os.path.join(px_zbrains_path, sid, ses) if ses else os.path.join(px_zbrains_path, sid)

    if "proc" or "analysis" in tasks:
        SUBJECT_MICAPIPE_DIR = os.path.join(dataset_path, args.micapipe, sid, ses) if ses else os.path.join(dataset_path, args.micapipe, sid)
        # print(SUBJECT_MICAPIPE_DIR)
        assert_exists(SUBJECT_MICAPIPE_DIR, f"{BIDS_ID} micapipe directory does not exist.")

        if "hippocampus" in structures:
            SUBJECT_HIPPUNFOLD_DIR = os.path.join(dataset_path, args.hippunfold, "hippunfold", sid, ses) if ses else os.path.join(dataset_path, args.hippunfold, "hippunfold", sid)
            assert_exists(SUBJECT_HIPPUNFOLD_DIR, f"{BIDS_ID} hippunfold directory does not exist.")

        SUBJECT_PLUGIN_DIR = os.path.join(dataset_path, args.plugin, sid, ses or "") if args.plugin else None
        if SUBJECT_PLUGIN_DIR and not os.path.exists(SUBJECT_PLUGIN_DIR):
            sys.exit(f"{BIDS_ID} plugin directory does not exist.")

        if "subcortex" in structures:
            subject_micapipe_qc = os.path.join(SUBJECT_MICAPIPE_DIR, "QC")
            Nrecon = len(glob.glob(f"{subject_micapipe_qc}/{BIDS_ID}_module-proc_surf-*.json"))
            if Nrecon < 1:
                show_error(f"{BIDS_ID} doesn't have a module-proc_surf: run -proc_surf")
                raise ProcessingException(f"{BIDS_ID} doesn't have a module-proc_surf: run -proc_surf")
            elif Nrecon == 1:
                module_qc = glob.glob(f"{subject_micapipe_qc}/{BIDS_ID}_module-proc_surf-*.json")[0]
                recon = os.path.splitext(module_qc)[0].split('proc_surf-')[1]
            elif Nrecon > 1:
                show_error(f"{BIDS_ID} has been processed with freesurfer and fastsurfer. Not supported yet")
                raise ProcessingException(f"{BIDS_ID} has been processed with freesurfer and fastsurfer. Not supported yet")

            SUBJECT_SURF_DIR = os.path.join(dataset_path, recon, BIDS_ID)
            assert_exists(SUBJECT_SURF_DIR, f"{BIDS_ID} {recon} directory does not exist.")

    return BIDS_ID, SUBJECT_OUTPUT_DIR, SUBJECT_MICAPIPE_DIR, SUBJECT_HIPPUNFOLD_DIR, SUBJECT_PLUGIN_DIR, SUBJECT_SURF_DIR, px_zbrains_path

def create_directories(BIDS_ID, SUBJECT_OUTPUT_DIR, FOLDER_LOGS, FOLDER_MAPS, FOLDER_NORM_Z, FOLDER_NORM_MODEL, tasks):
    # Create subject output directory structure
    if not os.path.isdir(SUBJECT_OUTPUT_DIR):
        show_info(f"Subject {BIDS_ID} directory doesn't exist, creating...")

    # Folder for logs
    logs_dir = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS)
    os.makedirs(logs_dir, exist_ok=True)

    if "proc" in tasks:
        dirs = [FOLDER_SCTX, FOLDER_CTX, FOLDER_HIP]
        for dir in dirs:
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_MAPS, dir)
            os.makedirs(path, exist_ok=True)

    # Folders for regional/asymmetry analysis
    if "analysis" in tasks:
        dirs = [FOLDER_SCTX, FOLDER_CTX, FOLDER_HIP]
        for dir in dirs:
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_NORM_Z, dir)
            os.makedirs(path, exist_ok=True)
            path = os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_NORM_MODEL, dir)
            os.makedirs(path, exist_ok=True)
    
    return logs_dir
            
def main_func(args):
    
    args, ZBRAINS, script_dir, VERBOSE, tasks, sid, ses, structures, features, resolutions, labels_ctx, labels_hip, smooth_ctx, smooth_hip, threshold = parse_args(args)

    check_workbench_dependency(tasks)

    BIDS_ID, SUBJECT_OUTPUT_DIR, SUBJECT_MICAPIPE_DIR, SUBJECT_HIPPUNFOLD_DIR, SUBJECT_PLUGIN_DIR, SUBJECT_SURF_DIR, px_zbrains_path = check_files_and_directories(args, tasks, structures, sid, ses)




    # ------------------------------------------------------------------------------------------------------------------- #
    # -------------------------------------------- Start processing/analysis -------------------------------------------- #
    # ------------------------------------------------------------------------------------------------------------------- #

    # -------------------------------------------------- Initiate timer ------------------------------------------------- #
    start_time = time.time()

    # -------------------------------------------- Create zbrains directory --------------------------------------------- #
    
    logs_dir = create_directories(BIDS_ID, SUBJECT_OUTPUT_DIR, FOLDER_LOGS, FOLDER_MAPS, FOLDER_NORM_Z, FOLDER_NORM_MODEL, tasks)


    # Temporary folder
    with tempdir(SUBJECT_OUTPUT_DIR, prefix="z_brains_temp.") as tmp_dir:
        try:
            os.chmod(SUBJECT_OUTPUT_DIR, 0o770)

            # -------------------------------------------------- Postprocessing ------------------------------------------------- #
            if "proc" in tasks:
                logfile = os.path.join(logs_dir, f"proc_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt")
                with open(logfile, 'w') as f:
                    f.write("\n")
                for struct in structures:
                    if struct == "hippocampus":
                        list_resolutions = [map_resolution_hip[k] for k in resolutions]
                        temp_labels = labels_hip
                        fwhm = smooth_hip
                    elif struct == "cortex":
                        list_resolutions = [map_resolution_ctx[k] for k in resolutions]
                        temp_labels = labels_ctx
                        fwhm = smooth_ctx
                    else:
                        list_resolutions = None
                        temp_labels = None
                        fwhm = None
                    run_proc.run(structure=struct, features=features, tmp_dir=tmp_dir, WORKBENCH_PATH=os.environ['WORKBENCH_PATH'], subject_micapipe_dir=SUBJECT_MICAPIPE_DIR, subject_output_dir=SUBJECT_OUTPUT_DIR, folder_maps=FOLDER_MAPS, folder_ctx=FOLDER_CTX, folder_sctx=FOLDER_SCTX, folder_hip=FOLDER_HIP, subject_surf_dir=SUBJECT_SURF_DIR, subject_hippunfold_dir=SUBJECT_HIPPUNFOLD_DIR, script_dir=script_dir, BIDS_ID=BIDS_ID, VERBOSE=VERBOSE, fwhm=fwhm, resolutions=list_resolutions, labels=temp_labels,subject_plugin_dir=SUBJECT_PLUGIN_DIR)

            # ----------------------------------------------------- Analysis ---------------------------------------------------- #
            if "analysis" in tasks:
                logfile = os.path.join(logs_dir, f"analysis_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt")
                ref_zbrains_paths = [os.path.join(args.dataset_ref[idx], 'derivatives', args.zbrains_ref[idx]) for idx in range(len(args.dataset_ref))]
                args_list = [
                    os.path.join(ZBRAINS, 'functions', 'run_analysis.py'), 
                    '--sub', sid, 
                    '--ses', ses, 
                    '--zbrains', px_zbrains_path, 
                    '--struct', '-'.join(structures),
                    '--feat', '-'.join(features),
                    '--demo_ref', '-'.join(args.demo_ref),
                    '--zbrains_ref', '-'.join(ref_zbrains_paths),
                    '--resolution', '-'.join(resolutions),
                    '--labels_ctx', labels_ctx[0],
                    '--labels_hip', labels_hip[0],
                    '--smooth_ctx', smooth_ctx,
                    '--smooth_hip', smooth_hip,
                    '--threshold', str(threshold),
                    '--approach', 'zscore',
                    '--logfile', logfile,
                    '--tmp', tmp_dir,
                    '--verbose', str(VERBOSE),
                    '--demo', args.demo,
                    '--column_map', str(args.column_map),
                    '--n_jobs', str(args.n_jobs)
                ]

                if args.normative:
                    args_list.extend(['--normative', args.normative])

                if args.deconfound:
                    args_list.extend(['--deconfound', args.deconfound])
                out = subprocess.call(['python', *args_list])

            #     run_analysis.run(
            #     subject_id=sid,
            #     session=ses,
            #     zbrains_ref=ref_zbrains_paths,
            #     demo_ref=args.demo_ref,
            #     zbrains=px_zbrains_path,
            #     struct=structures,
            #     feat=features,
            #     smooth_ctx=smooth_ctx,
            #     smooth_hip=smooth_hip,
            #     threshold=threshold,
            #     approach='zscore',
            #     resolution=resolutions,
            #     labels_ctx=labels_ctx,
            #     labels_hip=labels_hip,
            #     logfile=logfile,
            #     tmp=tmp_dir,
            #     verbose=VERBOSE,
            #     demo=args.demo,
            #     normative=args.normative, 
            #     deconfound=args.deconfound, 
            #     column_map=args.column_map ,
            #     n_jobs=args.n_jobs
            # )
            elapsed = (time.time() - start_time) / 60
            show_title(f"Total elapsed time for {BIDS_ID}: {elapsed:.2f} minutes")
        except Exception as e:
            print(e)

def check_sub(args, sub, ses=None):
    micapipe_path = os.path.join(args.dataset, 'derivatives', args.micapipe, sub)
    hippunfold_path = os.path.join(args.dataset, 'derivatives', args.hippunfold, 'hippunfold', sub)

    if ses is not None:
        micapipe_path = os.path.join(micapipe_path, ses)
        hippunfold_path = os.path.join(hippunfold_path, ses)

    if not os.path.exists(micapipe_path):
        print(f'No micapipe at {micapipe_path} for {sub}-{f"-{ses}" if ses else ""}, skipping')
        return False

    if not os.path.exists(hippunfold_path):
        print(f'No hippunfold at {hippunfold_path} for {sub}{f"-{ses}" if ses else ""}, skipping')
        return False

    return True

def create_jobs(args, subs, ses, run_type):
    jobs = []
    for sub in subs:
        for s in ses or os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe, sub)):
            if check_sub(args, sub, s):
                job = copy.copy(args)
                job.sub, job.ses, job.run = sub, s, run_type
                if args.patient_prefix in job.sub or run_type == 'proc':
                    jobs.append(job)
    return jobs

def main(args):
    # Assume args is parsed from command line arguments
    args.ses = args.ses.split(' ') if args.ses else None
    args.sub = args.sub.split(' ') if args.sub != 'all' else [sub for sub in os.listdir(os.path.join(args.dataset, 'derivatives', args.micapipe)) if check_sub(args, sub)]

    if args.ses and len(args.ses) != len(args.sub):
        print('Number of subs and sessions do not match')
        sys.exit()

    if 'proc' in args.run:
        procjobs = create_jobs(args, args.sub, args.ses, 'proc')
        Parallel(n_jobs=args.n_jobs)(delayed(jobloop)(job) for job in procjobs)
    if 'analysis' in args.run:
        analysisjobs = create_jobs(args, args.sub, args.ses, 'analysis')
        [jobloop(job) for job in analysisjobs]

class Parser(argparse.ArgumentParser):

        def print_help(self):
            print(help)

if __name__ == '__main__':
    

    # Define the argument parser
            
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
    parser.add_argument('--n_jobs_wb', type=int, default=8)
    parser.add_argument('--wb_path', type=str, default="/data/mica1/01_programs/workbench-1.4.2/bin_linux64")
    parser.add_argument('--patient_prefix', type=str, default="PX")
    
    parser.add_argument('--verbose', type=int, default=-1)
    parser.add_argument('--version', action='version', version='1.0.0')
  
    

    # Parse the arguments
    args = parser.parse_args()
    WORKBENCH_PATH = setenv(args.wb_path)
    os.environ['WORKBENCH_PATH'] = WORKBENCH_PATH
    runs = args.run
    os.environ['OMP_NUM_THREADS'] = str(args.n_jobs_wb)
    display_flag = False
    # delete_temp_folders(os.path.join(args.dataset, 'derivatives', args.zbrains))
    show_info("zbrains is running with:")
    if "proc" in args.run:
        # Get WorkBench version
        workbench_version = subprocess.check_output([os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command'), '-version']).decode().split()[1]
        show_note("WorkBench...", workbench_version)
        show_note("            ", os.path.join(os.environ['WORKBENCH_PATH'], 'wb_command'))
    
    # Get Python version
    python_version = subprocess.check_output(['python', '--version']).decode().split()[1]
    show_note("python......", python_version)
    show_note("            ", shutil.which('python'))

    show_note("", "")
    
    main(args)

    


