import argparse
from functions.utilities import assert_same_size, assert_required, show_error, assert_exists, show_warning, assert_columns_in_csv,show_info,do_cmd,submit_job
import os
import signal
import subprocess
from multiprocessing import cpu_count
import time
from functions.new_constants import FOLDER_LOGS
from pathlib import Path
from datetime import datetime
import colorama
import csv
def main(known, unknown):
    known.column_map = dict([arg.split('=') for arg in known.column_map]) if known.column_map else None
    assert_same_size("--dataset", known.dataset, "--demo", known.demo)
    n_datasets = len(known.dataset)

    if n_datasets != len(known.zbrains):
        if len(known.zbrains) > 1:
            assert_same_size("--dataset", known.dataset, "--zbrains", known.zbrains)
        else:
            zbrains_dirs = [known.zbrains[0]] * n_datasets

    if known.run == "proc":
        assert_required("--micapipe", known.micapipe)
        if n_datasets != len(known.micapipe):
            if len(known.micapipe) > 1:
                assert_same_size("--dataset", known.dataset, "--micapipe", known.micapipe)
            else:
                micapipe_dirs = [known.micapipe[0]] * n_datasets

        assert_required("--hippunfold", known.hippunfold)
        if n_datasets != len(known.hippunfold):
            if len(known.hippunfold) > 1:
                assert_same_size("--dataset", known.dataset, "--hippunfold", known.hippunfold)
            else:
                hip_dirs = [known.hippunfold[0]] * n_datasets

    if known.run == "analysis":
        assert_required("--demo_ref", known.demo_ref, "--demo_ref option required for analysis.")
        assert_required("--dataset_ref", known.dataset_ref, "--dataset_ref option required for analysis.")
        assert_same_size("--dataset_ref", known.dataset_ref, "--demo_ref", known.demo_ref)

        n_ref_datasets = len(known.dataset_ref)
        assert_required("--zbrains_ref", known.zbrains_ref, "--zbrains_ref option required for analysis.")
        if n_ref_datasets != len(known.zbrains_ref):
            if len(known.zbrains_ref) > 1:
                assert_same_size("--dataset_ref", known.dataset_ref, "--zbrains_ref", known.zbrains_ref)
            else:
                ref_zbrains_dirs = [known.zbrains_ref[0]] * n_ref_datasets

    scheduler = known.scheduler or "LOCAL" 
    scheduler_options = known.scheduler_options or ""

    if hasattr(unknown, 'sub'):
        show_error("Subject ID should not be provided. Please remove --sub.")
        exit(1)
    if hasattr(unknown, 'ses'):
        show_error("Session should not be provided: Please remove --ses.")
        exit(1)

    # --------------------------------------------------- CSV columns --------------------------------------------------- #
    expected_to_actual = {'participant_id': 'participant_id', 'session_id': 'session_id', 'age': 'age', 'sex': 'sex', 'site': 'site'}
    for key, value in known.column_map.items():
        key = key.strip()
        value = value.strip()
        if key in expected_to_actual:
            expected_to_actual[key] = value

    # Columns needed for clinical report
    columns_report = [expected_to_actual["age"], expected_to_actual["sex"]]

    # Required columns
    required_columns = ['participant_id']
    if known.run == "analysis":
        required_columns.extend(known.normative)
        required_columns.extend(known.deconfound)
    for idx, key in enumerate(required_columns):
        key = key.lstrip('-')  # For deconfound, column names may be prepended with a -
        required_columns[idx] = expected_to_actual.get(key, key)

    # Remove duplicates
    required_columns = list(set(required_columns))

    # Store column indices in each csv file
    col2idx = {}
    for _csv in known.demo + known.demo_ref:
        assert_exists(_csv)
        assert_columns_in_csv(_csv, required_columns)

        delimiter = '\t' if _csv.endswith('.tsv') else ','
        with open(_csv, 'r') as f:
            header = next(csv.reader(f, delimiter=delimiter))
        if known.run == "analysis" and _csv in known.demo:
            for col in columns_report:
                if col not in header:
                    show_warning(f"Column '{col}' not found in {_csv} for clinical report.")
        for idx, col in enumerate(header):
            col2idx[f"{col}"] = idx
    print(f"Debug: col2idx Array: {list(col2idx.keys())}")


    # ------------------------------------------- Setting important variables ------------------------------------------- #
    # Get number of parallel processes for LOCAL execution
    if scheduler == "LOCAL":
        n_jobs = 1  # default is 1 -> run serially
        max_jobs = cpu_count()
        n_jobs = min(n_jobs, max_jobs)

        scheduler_options = ""

    # ----------------------------------------------------- Cleanup ----------------------------------------------------- #
    # Function to clean up and kill child processes
    def cleanup(signum, frame):
        print("Terminating the script. Cleaning up...")

        # Kill all child processes
        subprocess.run(['pkill', '-P', str(os.getpid())])

        print("Cleanup complete. Exiting.")
        exit()

    # Trap signals and call the cleanup function
    signal.signal(signal.SIGINT, cleanup)
    signal.signal(signal.SIGTERM, cleanup)

    def submit_job_with_timer(scheduler, scheduler_options, log_file, cmd,run,sid,ses):
        show_info(f"Starting {run} for {sid}{f'_{ses}' if ses else ''}")

        start_time = time.time()

        # Submit the job to the scheduler
        process = subprocess.Popen([scheduler, scheduler_options, do_cmd, cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        with open(log_file, 'w') as f:
            for line in iter(process.stdout.readline, b''):
                f.write(line.decode('utf-8'))

        process.wait()

        end_time = time.time()

        elapsed_time = round(end_time - start_time, 2)

        # Print the elapsed time
        show_info(f"Job {run} for {sid}{f'_{ses}' if ses else ''} finished in {elapsed_time} seconds", f"Check logs: {log_file}")


    # ---------------------------------- Post-processing/analysis of controls/patients ---------------------------------- #
    if known.run == "proc":
        show_info("Running post-processing")
    else:
        show_info("Running analysis")

    col_sid = col2idx["participant_id"]
    col_ses = col2idx["session_id"]
    ZBRAINS = Path(os.path.realpath(__file__)).parent
    for iter, csv_path in enumerate(known.demo):
        with open(csv_path, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t' if csv_path.endswith('.tsv') else ',')
            next(reader)  # Skip the header
            for row in reader:
                sid = f"sub-{row[col_sid].replace('sub-', '')}"
                ses = row[col_ses]
                ses = "" if ses == "n/a" else f"ses-{ses.replace('ses-', '')}"
                dataset_path = known.dataset[iter]
                zbrains_dir = known.zbrains[iter]

                if known.run == "proc":
                    cmd = f"{ZBRAINS}/zbrains.py --run {known.run} --sub {sid} --dataset {dataset_path} --zbrains {zbrains_dir} \
                            --micapipe {known.micapipe[iter]} --hippunfold {known.hippunfold[iter]} {' '.join(unknown)}"
                else:
                    cmd = f"{ZBRAINS}/zbrains.py --run {known.run} --sub {sid} --dataset {dataset_path} --zbrains {zbrains_dir} \
                            --demo_ref {known.demo_ref} --dataset_ref {known.dataset_ref} \
                            --zbrains_ref {ref_zbrains_dirs} --demo {csv_path} {' '.join(unknown)}"

                if ses:
                    cmd += f" --ses {ses}"

                print(cmd)

                out_dir = os.path.join(dataset_path, "derivatives", zbrains_dir)
                logs_dir = os.path.join(out_dir, sid, ses if ses else "", FOLDER_LOGS)
                os.makedirs(logs_dir, exist_ok=True)

                name = f"{sid}_{ses if ses else ''}_run-{known.run}"
                log_file = os.path.join(logs_dir, f"{name}_{datetime.now().strftime('%d-%m-%Y')}_{scheduler.lower()}.txt")

                if scheduler == "LOCAL" and n_jobs > 1:
                    submit_job_with_timer(name, scheduler, scheduler_options, log_file, cmd)
                    while len(os.popen('jobs -r -p').read().splitlines()) >= n_jobs:
                        time.sleep(1)
                else:
                    if scheduler in ["PBS", "SGE"]:
                        opts = f"-N {name} -cwd -S/bin/bash -v ANTSPATH,WORKBENCH_PATH -j oe -o {log_file}"
                        scheduler_options_ext = f"{scheduler_options} {opts}"
                        submit_job(scheduler, scheduler_options_ext, cmd)
                        # do_cmd("SUBMIT_JOB", scheduler, scheduler_options_ext, cmd)
                    else:
                        submit_job(scheduler, scheduler_options, cmd)
                        # do_cmd("SUBMIT_JOB", scheduler, scheduler_options, cmd)

    if scheduler == "LOCAL" and n_jobs > 1:
        while len(os.popen('jobs -r -p').read().splitlines()) > 0:
            time.sleep(1)

if __name__ == '__main__':
    
    
    
    pcolor = colorama.Fore.MAGENTA  # Purple
    rcolor = colorama.Fore.RED  # Red
    gcolor = colorama.Fore.GREEN  # Green
    bcolor = colorama.Fore.BLUE  # Blue
    gray = colorama.Fore.LIGHTBLACK_EX  # Gray
    nc = colorama.Style.RESET_ALL  # No color
    
    LIST_TASKS = ["proc", "analysis"]
    LIST_SCHEDULERS = ["LOCAL", "SGE", "PBS"]  # SLURM
    
    help = f"""
    {pcolor}COMMAND:{nc}
    zbrains_batch.py --run <run> --demo <csv_subjects> [options...]


   For --run proc: run post-processing for all subjects provided with --demo target_subjects1 [target_subjects2 ...].

   For --run analysis: run analysis for all subjects provided with --demo target_subjects1 [target_subjects2 ...] wrt
   the reference group (i.e., --demo_ref csv_reference_subjects1 [csv_reference_subjects2 ...]). Analysis requires both
   target and reference subjects to be post-processed beforehand.

   Additional options that are not listed below are passed to the zbrains script.


{pcolor}OPTIONS:{nc}
\t{rcolor}--run{nc} run                 : Task to perform. Options:
                                      - proc          : perform post-processing for each target subject (i.e., those
                                                        provided with --demo).
                                      - analysis      : perform analysis (regional & asymmetry) and generate clinical
                                                        report. The analysis is performed for each target subject wrt
                                                        the reference subjects (--demo_ref). Post-processing of both
                                                        target and reference subjects must be performed beforehand.
\t{rcolor}--demo{nc} path [path ...]    : CSV/TSV files with demographics for target subjects. There must be one file
                                    for each target dataset (i.e., those provided with --dataset). These files must
                                    include subject and session (if available) identifiers, and other demographics. See
                                    the zbrains script for more info.
\t{rcolor}--dataset{nc} path [path ...] : Paths to BIDS datasets containing data for the target subjects. Each
                                    dataset must correspond to one file in --demo.
\t{rcolor}--zbrains{nc} dir [dir ...]   : Names of the zbrains derivative folder in each of the target datasets. If
                                    only one folder name is provided but there are multiple target datasets, we assume
                                    the same folder name for all datasets. For post-processing, the folder will be
                                    created if it does not exist.

\t{gcolor}--micapipe{nc} [dir ...]      : Name of the micapipe derivative folder in each of the target datasets. If
                                    only one folder name is provided but there are multiple target datasets, we assume
                                    the same folder name for all datasets. Required only for post-processing.
\t{gcolor}--hippunfold{nc} [dir ...]    : Name of the hippunfold derivative folder in each of the target datasets. If
                                    only one folder name is provided but there are multiple target datasets, we assume
                                    the same folder name for all datasets. Required only for post-processing.
\t{gcolor}--demo_ref{nc} [path ...]     : CSV/TSV files with demographics for reference subjects. Required only for
                                    analysis. There must be one file for each reference dataset (i.e., those provided
                                    with --dataset_ref). These files must include subject and session (if available)
                                    identifiers, and other demographics. See the zbrains script for more info.
\t{gcolor}--dataset_ref{nc} [path ...]  : Paths to BIDS datasets containing data for the target subjects. Each dataset
                                    must correspond to one file in --demo_ref. Required only for analysis.
\t{gcolor}--zbrains_ref{nc} [dir ...]   :Names of the zbrains derivative folder in each of the reference datasets. If
                                    only one folder name is provided but there are multiple reference datasets, we
                                    assume the same folder name for all datasets. Required only for analysis.
\t{gcolor}--scheduler{nc} [scheduler]   : Control for parallel computation. Options:
                                      - {bcolor}LOCAL{nc}       : run locally (default)
                                      - SGE           : SGE qsub
                                      - PBS           : PBS qsub
\t{gcolor}--scheduler_opts{nc} [opts]   : Scheduler options for SGE and PBS. For LOCAL, each subject is processed
                                    serially (default) but you can provide \"-n N\" to run in parallel. Options
                                    should be provided as a single string. Example: \"-q queue -M my_email\".
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
    {pcolor}zbrains_batch.py{nc} {rcolor}--run{nc} proc
                  {rcolor}--demo{nc} <target_subjects1.csv> <target_subjects2.csv>
                  {rcolor}--dataset{nc} <target_dataset1> <target_dataset2>
                  {rcolor}--zbrains{nc} <target_zbrains_dir1> <target_zbrains_dir2>
                  {gcolor}--micapipe{nc} <target_micapipe_dir1> <target_micapipe_dir2>
                  {gcolor}--hippunfold{nc} <target_hippunfold_dir>    {gray}# same folder name for both datasets{nc}
                  {gcolor}--scheduler{nc} LOCAL
                  {gcolor}--scheduler_opts{nc} \"-n 6\"                 {gray}# run in parallel{nc}

    {gray}# Analysis{nc}
    {pcolor}zbrains_batch.py{nc} {rcolor}--run{nc} analysis
                     {rcolor}--demo{nc} <target_subjects1.csv>
                     {rcolor}--dataset{nc} <target_dataset1>
                     {rcolor}--zbrains{nc} <target_zbrains_dir1>
                     {gcolor}--demo_ref{nc} <reference_subjects1.csv> <reference_subjects2.csv>
                     {gcolor}--dataset_ref{nc} <reference_dataset1> <reference_dataset2>
                     {gcolor}--zbrains_ref{nc} <reference_zbrains_dir1> <reference_zbrains_dir2>


{pcolor}DEPENDENCIES:{nc}
    See the zbrains script for more info.


McGill University, MNI, MICA lab, Nov 2023
https://github.com/MICA-MNI/micapipe
https://github.com/MICA-MNI/z-brains
http://mica-mni.github.io/
    """
    class Parser(argparse.ArgumentParser):

        def print_help(self):
            print(help)
            
    parser = Parser()
    

    

    parser.add_argument("--run", required=True)
    parser.add_argument("--demo", nargs='+', required=True)
    parser.add_argument("--dataset", nargs='+', required=True)
    parser.add_argument("--micapipe", nargs='+', default=None)
    parser.add_argument("--hippunfold", nargs='+', default=None)
    parser.add_argument("--zbrains", nargs='+', required=True)
    parser.add_argument("--demo_ref", nargs='+', default=None)
    parser.add_argument("--dataset_ref", nargs='+', default=None)
    parser.add_argument("--zbrains_ref", nargs='+', default=None)
    parser.add_argument("--deconfound", nargs='+', default=None)
    parser.add_argument("--normative", nargs='+', default=None)
    parser.add_argument("--column_map", nargs='+', default=None)
    parser.add_argument("--scheduler", default=None)
    parser.add_argument("--scheduler_options", default=None)
    parser.add_argument("--verbose", action='store_true')

    known, unknown = parser.parse_known_args()
    main(known, unknown)