# FILEPATH: /c:/Users/Ian/Documents/GitHub/z-brains/tests/test_zbrains.py

import os
import struct
from unittest import mock
import pytest
from unittest.mock import patch, MagicMock
from argparse import Namespace
import tempfile
import shutil
import sys

from .pytest_constants import HC_DEMOGRAPHICS, PX_DEMOGRAPHICS, DATASET_DIR

from functions.constants import (
    ProcessingException,
    LIST_FEATURES,
    LIST_RESOLUTIONS,
    DEFAULT_SMOOTH_CTX,
    DEFAULT_SMOOTH_HIP,
    DEFAULT_THRESHOLD,
    FOLDER_NORM_Z,
    FOLDER_NORM_MODEL,
    FOLDER_MAPS,
    FOLDER_LOGS,
)
from zbrains import (
    main,
    create_jobs,
    check_sub,
    main_func,
    create_directories,
    check_files_and_directories,
    check_workbench_dependency,
    parse_args,
    tempdir,
)

# Test parse_args function


def test_tempdir():
    SUBJECT_OUTPUT_DIR = tempfile.gettempdir()
    prefix = "test_"

    # Test that the context manager creates a directory
    with tempdir(SUBJECT_OUTPUT_DIR, prefix) as temp_dir:
        assert os.path.isdir(temp_dir)
        assert temp_dir.startswith(os.path.join(SUBJECT_OUTPUT_DIR, prefix))

    # Test that the context manager cleans up the directory
    assert not os.path.exists(temp_dir)


# Test parse_args function


def test_parse_args():
    # Mock the args
    mock_args = Namespace(
        column_map=None,
        verbose=True,
        run=["proc"],
        micapipe="micapipe",
        struct=["hippocampus"],
        hippunfold="hippunfold",
        dataset_ref=DATASET_DIR,
        zbrains_ref=["zbrains"],
        demo_ref=HC_DEMOGRAPHICS,
        normative=True,
        demo=PX_DEMOGRAPHICS,
        sub="sub-001",
        ses="ses-001",
        feat=["all"],
        resolution=["all"],
        label_ctx=["midthickness"],
        label_hip=["midthickness"],
        smooth_ctx=None,
        smooth_hip=None,
        threshold=None,
        plugin=None,
    )

    # Call the function with the mock args
    result = parse_args(mock_args)
    FEATURES_NOPLUGIN = LIST_FEATURES.copy()
    FEATURES_NOPLUGIN.remove("plugin-*")
    # Check the result
    assert result[0] == mock_args
    assert str(result[1]) == os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    assert result[2].name == "functions"
    assert result[3] == True
    assert result[4] == ["proc"]
    assert result[5] == "sub-001"
    assert result[6] == "ses-001"
    assert result[7] == ["hippocampus"]
    assert result[8] == FEATURES_NOPLUGIN
    assert result[9] == LIST_RESOLUTIONS
    assert result[10] == ["midthickness"]
    assert result[11] == ["midthickness"]
    assert result[12] == DEFAULT_SMOOTH_CTX
    assert result[13] == DEFAULT_SMOOTH_HIP
    assert result[14] == DEFAULT_THRESHOLD


def test_parse_args_missing_required_option():
    # Mock the args
    mock_args = Namespace(
        column_map=None,
        verbose=True,
        run=["proc"],
        micapipe=None,
        struct=["hippocampus"],
        hippunfold="hippunfold",
        dataset_ref=["dataset_ref"],
        zbrains_ref=["zbrains_ref"],
        demo_ref=["demo_ref"],
        normative=True,
        demo="demo",
        sub="sub-001",
        ses="ses-001",
        feat=["all"],
        resolution=["all"],
        label_ctx=["midthickness"],
        label_hip=["midthickness"],
        smooth_ctx=None,
        smooth_hip=None,
        threshold=None,
    )

    # Call the function with the mock args and check for the exception
    with pytest.raises(ProcessingException):
        parse_args(mock_args)


# Test check_workbench_dependency function


@patch("os.environ", {"WORKBENCH_PATH": "/path/to/workbench"})
@patch("os.path.join")
@patch("shutil.which")
@patch("zbrains.show_error")
def test_check_workbench_dependency_workbench_path_set(
    mock_show_error, mock_which, mock_join
):
    mock_join.return_value = "/path/to/workbench/wb_command"
    tasks = ["proc"]
    check_workbench_dependency(tasks)
    mock_join.assert_called_once_with("/path/to/workbench", "wb_command")
    mock_show_error.assert_not_called()


@patch("os.environ", {"WORKBENCH_PATH": None})
@patch("os.path.join")
@patch("shutil.which")
@patch("zbrains.show_error")
def test_check_workbench_dependency_workbench_path_not_set(
    mock_show_error, mock_which, mock_join
):
    mock_which.return_value = "/usr/local/bin/wb_command"
    tasks = ["proc"]
    check_workbench_dependency(tasks)
    mock_which.assert_called_once_with("wb_command")
    mock_show_error.assert_not_called()


@patch("os.environ", {"WORKBENCH_PATH": None})
@patch("os.path.join")
@patch("shutil.which")
@patch("zbrains.show_error")
def test_check_workbench_dependency_wb_command_not_found(
    mock_show_error, mock_which, mock_join
):
    mock_which.return_value = None
    tasks = ["proc"]
    with pytest.raises(ProcessingException):
        check_workbench_dependency(tasks)
    mock_which.assert_called_once_with("wb_command")
    mock_show_error.assert_called_once_with(
        "Workbench not found. Please set the WORKBENCH_PATH environment variable to the location of Workbench binaries."
    )


# Test check_files_and_directories function


@patch("os.path.realpath")
@patch("os.path.join")
@patch("os.path.exists")
@patch("glob.glob")
def test_check_files_and_directories(mock_glob, mock_exists, mock_join, mock_realpath):
    HC_DEMOGRAPHICS = "E:/HC_participants.csv"
    PX_DEMOGRAPHICS = "E:/PX_participants.csv"
    DATASET_DIR = "E:/data/derivatives"

    # Mock the args
    mock_args = MagicMock()
    mock_args.dataset = DATASET_DIR
    mock_args.zbrains = "zbrains"
    mock_args.micapipe = "micapipe"
    mock_args.hippunfold = "hippunfold"
    mock_args.plugin = "plugin"

    # Mock the os and glob functions
    mock_realpath.return_value = DATASET_DIR
    mock_join.side_effect = lambda *args: "/".join(args)
    mock_exists.return_value = True
    mock_glob.return_value = ["/QC/BIDS_ID_module-proc_surf-1.json"]
    sid = "sub-HC002"
    ses = "ses-01"
    # Call the function with the mock args
    result = check_files_and_directories(
        mock_args, ["proc"], ["hippocampus", "subcortex"], sid, ses
    )

    # Check the result
    assert result == (
        f"{sid}_{ses}",
        f"{DATASET_DIR}/zbrains/{sid}/{ses}",
        f"{DATASET_DIR}/micapipe/{sid}/{ses}",
        f"{DATASET_DIR}/hippunfold/hippunfold/{sid}/{ses}",
        f"{DATASET_DIR}/plugin/{sid}/{ses}",
        f"{DATASET_DIR}/1/{sid}_{ses}",
        f"{DATASET_DIR}/zbrains",
    )


@patch("os.path.realpath")
@patch("os.path.join")
@patch("os.path.exists")
@patch("glob.glob")
def test_check_files_and_directories_missing_dir(
    mock_glob, mock_exists, mock_join, mock_realpath
):
    HC_DEMOGRAPHICS = "E:/HC_participants.csv"
    PX_DEMOGRAPHICS = "E:/PX_participants.csv"
    DATASET_DIR = "E:/data/derivatives"

    # Mock the args
    mock_args = MagicMock()
    mock_args.dataset = DATASET_DIR
    mock_args.zbrains = "zbrains"
    mock_args.micapipe = "micapipe"
    mock_args.hippunfold = "hippunfold"
    mock_args.plugin = "plugin"

    # Mock the os and glob functions
    mock_realpath.return_value = DATASET_DIR
    mock_join.side_effect = lambda *args: "/".join(args)
    mock_exists.return_value = False
    mock_glob.return_value = ["/QC/BIDS_ID_module-proc_surf-1.json"]

    sid = "sub-HC002"
    ses = "ses-01"
    # Call the function with the mock args and check for the exception
    with pytest.raises(ProcessingException):
        check_files_and_directories(
            mock_args, ["proc"], ["hippocampus", "subcortex"], sid, ses
        )


@patch("os.path.realpath")
@patch("os.path.join")
@patch("os.path.exists")
@patch("glob.glob")
def test_check_files_and_directories_missing_module(
    mock_glob, mock_exists, mock_join, mock_realpath
):
    # Mock the args
    mock_args = MagicMock()
    mock_args.dataset = "/dataset"
    mock_args.zbrains = "zbrains"
    mock_args.micapipe = "micapipe"
    mock_args.hippunfold = "hippunfold"
    mock_args.plugin = "plugin"

    # Mock the os and glob functions
    mock_realpath.return_value = "/real/dataset"
    mock_join.side_effect = lambda *args: "/".join(args)
    mock_exists.return_value = True
    mock_glob.return_value = []

    # Call the function with the mock args and check for the exception
    with pytest.raises(ProcessingException):
        check_files_and_directories(mock_args, ["proc"], ["subcortex"], "sid", "ses")


# Test create_directories function


@patch("os.makedirs")
@patch("os.path.isdir")
@patch("zbrains.show_info")
def test_create_directories(mock_show_info, mock_isdir, mock_makedirs):
    BIDS_ID = "test"
    SUBJECT_OUTPUT_DIR = "/path/to/output"
    FOLDER_LOGS = "logs"
    FOLDER_MAPS = "maps"
    FOLDER_NORM_Z = "norm_z"
    FOLDER_NORM_MODEL = "norm_model"
    FOLDER_SCTX = "sctx"
    FOLDER_CTX = "ctx"
    FOLDER_HIP = "hip"
    tasks = ["proc", "analysis"]

    # Test case where directory does not exist
    mock_isdir.return_value = False
    result = create_directories(
        BIDS_ID,
        SUBJECT_OUTPUT_DIR,
        FOLDER_LOGS,
        FOLDER_MAPS,
        FOLDER_NORM_Z,
        FOLDER_NORM_MODEL,
        tasks,
    )
    assert result == os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS)
    mock_show_info.assert_called_once_with(
        f"Subject {BIDS_ID} directory doesn't exist, creating..."
    )
    mock_makedirs.assert_any_call(
        os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS), exist_ok=True
    )

    mock_show_info.reset_mock()
    mock_makedirs.reset_mock()

    # Test case where directory exists
    mock_isdir.return_value = True
    result = create_directories(
        BIDS_ID,
        SUBJECT_OUTPUT_DIR,
        FOLDER_LOGS,
        FOLDER_MAPS,
        FOLDER_NORM_Z,
        FOLDER_NORM_MODEL,
        tasks,
    )
    assert result == os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS)
    mock_show_info.assert_not_called()
    mock_makedirs.assert_any_call(
        os.path.join(SUBJECT_OUTPUT_DIR, FOLDER_LOGS), exist_ok=True
    )


# Test main_func function


@patch("zbrains.parse_args")
@patch("zbrains.check_workbench_dependency")
@patch("zbrains.check_files_and_directories")
@patch("zbrains.create_directories")
@patch("zbrains.tempdir")
@patch("zbrains.run_proc.run")
def test_main_func(
    mock_run_proc,
    mock_tempdir,
    mock_create_directories,
    mock_check_files_and_directories,
    mock_check_workbench_dependency,
    mock_parse_args,
    create_directory,
):
    mock_args = MagicMock()
    mock_args.tasks = ["proc", "analysis"]
    mock_args.structures = ["hippocampus", "cortex"]
    mock_args.features = ["ADC", "FA"]
    mock_args.resolutions = ["high", "low"]
    mock_args.labels_hip = ["white", "midthickness"]
    mock_args.labels_ctx = ["midthickness"]
    mock_args.smooth_hip = 1.0
    mock_args.smooth_ctx = 2.0
    mock_args.threshold = 0.5
    mock_args.WORKBENCH_PATH = "/path/to/workbench"

    mock_parse_args.return_value = (
        mock_args,
        "ZBRAINS",
        "/path/to/script",
        True,
        mock_args.tasks,
        "sid",
        "ses",
        mock_args.structures,
        mock_args.features,
        mock_args.resolutions,
        mock_args.labels_ctx,
        mock_args.labels_hip,
        mock_args.smooth_ctx,
        mock_args.smooth_hip,
        mock_args.threshold,
    )

    mock_check_files_and_directories.return_value = (
        "BIDS_ID",
        "/path/to/SUBJECT_OUTPUT_DIR",
        "/path/to/SUBJECT_MICAPIPE_DIR",
        "/path/to/SUBJECT_HIPPUNFOLD_DIR",
        "/path/to/SUBJECT_PLUGIN_DIR",
        "/path/to/SUBJECT_SURF_DIR",
        "/path/to/px_zbrains_path",
    )

    mock_create_directories.return_value = "/path/to/logs_dir"

    mock_tempdir.return_value.__enter__.return_value = "/path/to/tmp_dir"
    with patch.dict("os.environ", {"WORKBENCH_PATH": "/path/to/workbench"}):
        main_func(mock_args)

    mock_parse_args.assert_called_once_with(mock_args)
    mock_check_workbench_dependency.assert_called_once_with(mock_args.tasks)
    mock_check_files_and_directories.assert_called_once_with(
        mock_args, mock_args.tasks, mock_args.structures, "sid", "ses"
    )
    mock_create_directories.assert_called_once_with(
        "BIDS_ID",
        "/path/to/SUBJECT_OUTPUT_DIR",
        FOLDER_LOGS,
        FOLDER_MAPS,
        FOLDER_NORM_Z,
        FOLDER_NORM_MODEL,
        mock_args.tasks,
    )
    mock_tempdir.assert_called_once_with(
        "/path/to/SUBJECT_OUTPUT_DIR", prefix="z_brains_temp."
    )
    mock_run_proc.assert_called()


@pytest.fixture
def create_directory():
    os.makedirs("/path/to/SUBJECT_OUTPUT_DIR", exist_ok=True)
    os.makedirs("/path/to/logs_dir", exist_ok=True)
    yield
    os.rmdir("/path/to/SUBJECT_OUTPUT_DIR")


# Test check_sub function


@patch("os.path.exists")
@patch("builtins.print")
def test_check_sub(mock_print, mock_exists):
    mock_args = MagicMock()
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.hippunfold = "hippunfold"
    sub = "sub-001"
    ses = "ses-01"

    # Test case where directories exist
    mock_exists.return_value = True
    result = check_sub(mock_args, sub, ses)
    assert result == True
    assert mock_print.call_count == 0

    # Test case where micapipe directory does not exist
    mock_exists.side_effect = [False, True]
    result = check_sub(mock_args, sub, ses)
    assert result == False
    mock_print.assert_called_once_with(
        "No micapipe at /path/to/dataset\\derivatives\\micapipe\\sub-001\\ses-01 for sub-001-ses-01, skipping"
    )

    mock_print.reset_mock()

    # Test case where hippunfold directory does not exist
    mock_exists.side_effect = [True, False]
    result = check_sub(mock_args, sub, ses)
    assert result == False
    mock_print.assert_called_once_with(
        "No hippunfold at /path/to/dataset\\derivatives\\hippunfold\\hippunfold\\sub-001\\ses-01 for sub-001-ses-01, skipping"
    )


## Test create_jobs function


@patch("zbrains.os.listdir")
@patch("zbrains.check_sub")
@patch("copy.copy")
def test_create_jobs(mock_copy, mock_check_sub, mock_listdir):
    mock_args = MagicMock()
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.patient_prefix = "patient"
    subs = ["sub1", "sub2"]
    ses = ["ses1", "ses2"]
    run_type = "proc"

    mock_listdir.return_value = ["ses1", "ses2"]
    mock_check_sub.return_value = True
    mock_copy.return_value = mock_args

    jobs = create_jobs(mock_args, subs, ses, run_type)

    assert len(jobs) == 2
    mock_copy.assert_called()
    mock_check_sub.assert_called()
    mock_listdir.assert_not_called()


@patch("zbrains.os.listdir")
@patch("zbrains.check_sub")
@patch("copy.copy")
def test_create_jobs_no_ses(mock_copy, mock_check_sub, mock_listdir):
    mock_args = MagicMock()
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.patient_prefix = "patient"
    subs = ["sub1", "sub2"]
    ses = None
    run_type = "proc"

    mock_listdir.return_value = ["ses1", "ses2"]
    mock_check_sub.return_value = True
    mock_copy.return_value = mock_args

    jobs = create_jobs(mock_args, subs, ses, run_type)

    assert len(jobs) == 4
    mock_copy.assert_called()
    mock_check_sub.assert_called()
    mock_listdir.assert_called()


@patch("zbrains.os.listdir")
@patch("zbrains.check_sub")
@patch("copy.copy")
def test_create_jobs_check_sub_false(mock_copy, mock_check_sub, mock_listdir):
    mock_args = MagicMock()
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.patient_prefix = "patient"
    subs = ["sub1", "sub2"]
    ses = ["ses1", "ses2"]
    run_type = "proc"

    mock_listdir.return_value = ["ses1", "ses2"]
    mock_check_sub.return_value = False
    mock_copy.return_value = mock_args

    jobs = create_jobs(mock_args, subs, ses, run_type)

    assert len(jobs) == 0
    mock_copy.assert_not_called()
    mock_check_sub.assert_called()
    mock_listdir.assert_not_called()


## Test main function
@patch("zbrains.delete_temp_folders")
@patch("zbrains.show_info")
@patch("zbrains.show_note")
@patch("zbrains.check_sub")
@patch("zbrains.create_jobs")
@patch("zbrains._jobloop")
def test_main(
    mock__jobloop,
    mock_create_jobs,
    mock_check_sub,
    mock_show_note,
    mock_show_info,
    mock_delete_temp_folders,
):
    mock_args = MagicMock()
    mock_args.wb_path = "/path/to/workbench"
    mock_args.n_jobs_wb = 4
    mock_args.delete_temps = False
    mock_args.run = ["proc"]
    mock_args.ses = None
    mock_args.sub = "all"
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.n_jobs = 2
    mock_args.zbrains = "zbrains"

    mock_check_sub.return_value = True
    mock_create_jobs.return_value = ["job1", "job2"]
    with pytest.raises(FileNotFoundError):
        main(mock_args)


@patch("zbrains.delete_temp_folders")
def test_main_delete_temps(
    mock_delete_temp_folders,
):
    mock_args = MagicMock()
    mock_args.wb_path = "/path/to/workbench"
    mock_args.delete_temps = True
    mock_args.dataset = "/path/to/dataset"
    mock_args.zbrains = "zbrains"
    with pytest.raises(SystemExit):
        main(mock_args)

    mock_delete_temp_folders.assert_called_once_with(
        os.path.join(mock_args.dataset, "derivatives", mock_args.zbrains)
    )


@patch("zbrains.show_info")
@patch("zbrains.show_note")
@patch("zbrains.check_sub")
@patch("zbrains.create_jobs")
@patch("zbrains._jobloop")
def test_main_mismatch_subs_ses(
    mock__jobloop,
    mock_create_jobs,
    mock_check_sub,
    mock_show_note,
    mock_show_info,
):
    mock_args = MagicMock()
    mock_args.wb_path = "/path/to/workbench"
    mock_args.n_jobs_wb = 4
    mock_args.delete_temps = False
    mock_args.run = ["proc"]
    mock_args.ses = "ses1 ses2"
    mock_args.sub = "sub1"
    mock_args.dataset = "/path/to/dataset"
    mock_args.micapipe = "micapipe"
    mock_args.n_jobs = 2
    mock_args.zbrains = "zbrains"

    with pytest.raises(SystemExit):
        main(mock_args)


def test_main_fullrun():
    rootdir = "E:/BIDS_MICS_Test"
    DATASET_DIR = f"{rootdir}/data"
    HC_DEMOGRAPHICS = f"{rootdir}/HC_participants.csv"
    PX_DEMOGRAPHICS = f"{rootdir}/PX_participants.csv"
    with tempfile.TemporaryDirectory() as tmpdir:

        WBPATH = (
            "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64"
        )

        sub = "sub-HC005 sub-PX001"
        ses = "all"
        mock_args = Namespace(
            wb_path=WBPATH,
            sub=sub,
            ses=ses,
            run=["proc", "analysis"],
            dataset=DATASET_DIR,
            zbrains=tmpdir,
            micapipe="micapipe",
            hippunfold="hippunfold",
            plugin=None,
            column_map="participant_id=ID session_id=SES",
            verbose=3,
            struct=["all"],
            dataset_ref=DATASET_DIR,
            zbrains_ref=tmpdir,
            demo_ref=HC_DEMOGRAPHICS,
            normative=None,
            deconfound=None,
            demo=PX_DEMOGRAPHICS,
            feat=["all"],
            resolution=["all"],
            label_ctx=["midthickness"],
            label_hip=["midthickness"],
            smooth_ctx=5,
            smooth_hip=2,
            delete_temps=False,
            patient_prefix="PX",
            threshold=None,
            n_jobs_wb=4,
            n_jobs=4,
            test=True,
        )

        try:
            main(mock_args)
        except Exception as e:
            pytest.fail(f"main raised an exception: {e}")
        print_directory_tree(tmpdir)
        # print(len(os.listdir(os.path.join(tmpdir, "ses-02"))))
        assert "sub-PX001_ses-02_approach-zscore_summary-report.pdf" in os.listdir(
            os.path.join(tmpdir, "sub-PX001", "ses-02")
        )
        assert "logs" in os.listdir(os.path.join(tmpdir, "sub-PX001", "ses-02"))
        assert "maps" in os.listdir(os.path.join(tmpdir, "sub-PX001", "ses-02"))
        assert "norm-normative" in os.listdir(
            os.path.join(tmpdir, "sub-PX001", "ses-02")
        )
        assert "norm-z" in os.listdir(os.path.join(tmpdir, "sub-PX001", "ses-02"))
        assert (
            len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "maps", "cortex")
                )
            )
            + len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "maps", "subcortex")
                )
            )
            + len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "maps", "hippocampus")
                )
            )
        ) == 45
        assert (
            len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "norm-z", "cortex")
                )
            )
            + len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "norm-z", "subcortex")
                )
            )
            + len(
                os.listdir(
                    os.path.join(tmpdir, "sub-PX001", "ses-02", "norm-z", "hippocampus")
                )
            )
        ) == 84


def print_directory_tree(path):
    for root, dirs, files in os.walk(path):
        level = root.replace(path, "").count(os.sep)
        indent = " " * 4 * (level)
        print("{}{}/".format(indent, os.path.basename(root)))
        subindent = " " * 4 * (level + 1)
        for f in files:
            print("{}{}".format(subindent, f))
