import os
import pytest
from unittest.mock import patch, mock_open
from functions.run_proc import (
    map_hippocampus,
    map_cortex,
    map_subcortex,
    subcortical_mapping,
)
import pandas as pd
import numpy as np


@patch("numpy.asanyarray")
@patch("nib.load")
@patch("pandas.DataFrame")
@patch("pandas.read_csv")
@patch("builtins.open")
def test_subcortical_mapping(
    mock_open,
    mock_read_csv,
    mock_DataFrame,
    mock_nib_load,
    mock_asanyarray,
):
    # Setup
    mock_asanyarray.return_value = np.array([1, 2, 3])
    mock_nib_load.return_value = mock_open(read_data="data")
    mock_DataFrame.return_value = pd.DataFrame()
    mock_read_csv.return_value = pd.DataFrame(
        {
            "label": [26, 18, 11],
            "volume": [1, 2, 3],
            "structure": ["Laccumb", "Lamyg", "Lcaud"],
        }
    )

    # Test when seg and image are provided and vol is None
    subcortical_mapping(
        "subject_id",
        "image",
        "seg",
        None,
        "output",
        False,
    )
    # Assertions
    mock_asanyarray.assert_called()
    mock_nib_load.assert_called()
    mock_DataFrame.assert_called()

    # Test when seg and image are None and vol is provided
    subcortical_mapping(
        "subject_id",
        None,
        None,
        "vol",
        "output",
        False,
    )
    # Assertions
    mock_read_csv.assert_called()
    mock_DataFrame.assert_called()

    # Test when seg, image and vol are None
    try:
        subcortical_mapping(
            "subject_id",
            None,
            None,
            None,
            "output",
            False,
        )
    except ValueError as e:
        assert (
            str(e)
            == "Ambiguous inputs. Please provide either 'image' and 'seg' together or 'vol' alone."
        )

    # Test when seg, image and vol are provided
    try:
        subcortical_mapping(
            "subject_id",
            "image",
            "seg",
            "vol",
            "output",
            False,
        )
    except ValueError as e:
        assert (
            str(e)
            == "Ambiguous inputs. Please provide either 'image' and 'seg' together or 'vol' alone."
        )


# Test map_subcortex function


@patch("os.path.isfile")
@patch("os.makedirs")
@patch("functions.run_proc.subcortical_mapping")
@patch("functions.utilities.show_warning")
@patch("functions.utilities.show_note")
def test_map_subcortex(
    mock_show_note,
    mock_show_warning,
    mock_subcortical_mapping,
    mock_os_makedirs,
    mock_os_path_isfile,
):
    # Setup
    mock_os_path_isfile.return_value = True
    mock_os_makedirs.return_value = None
    mock_subcortical_mapping.return_value = None

    # Test when feat starts with "plugin-"
    map_subcortex(
        "bids_id",
        "plugin-test",
        "subject_surf_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_sctx",
        "script_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subcortical_mapping.assert_called()

    # Test when feat is in map_input
    map_subcortex(
        "bids_id",
        "flair",
        "subject_surf_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_sctx",
        "script_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subcortical_mapping.assert_called()

    # Test when feat is not in map_input
    map_subcortex(
        "bids_id",
        "not_in_map_input",
        "subject_surf_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_sctx",
        "script_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subcortical_mapping.assert_called()

    # Test when required files do not exist
    mock_os_path_isfile.return_value = False
    map_subcortex(
        "bids_id",
        "feat",
        "subject_surf_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_sctx",
        "script_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subcortical_mapping.assert_not_called()
    mock_show_warning.assert_called()

    # Test when output file is not created successfully
    mock_os_path_isfile.side_effect = [True, True, False]
    map_subcortex(
        "bids_id",
        "feat",
        "subject_surf_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_sctx",
        "script_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subcortical_mapping.assert_called()
    mock_show_warning.assert_called()


# Test map_hippocampus function


@patch("os.path.isfile")
@patch("os.makedirs")
@patch("subprocess.run")
@patch("glob.glob")
def test_map_hippocampus(
    mock_glob, mock_subprocess_run, mock_os_makedirs, mock_os_path_isfile
):
    # Setup
    mock_os_path_isfile.return_value = True
    mock_subprocess_run.return_value = None
    mock_os_makedirs.return_value = None
    mock_glob.return_value = ["file1"]

    # Test when feat starts with "plugin-"
    map_hippocampus(
        "bids_id",
        "plugin-test",
        "resol",
        "label",
        "workbench_path",
        "subject_hippunfold_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_hip",
        "fwhm",
        "tmp_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is in map_input
    map_hippocampus(
        "bids_id",
        "thickness",
        "resol",
        "label",
        "workbench_path",
        "subject_hippunfold_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_hip",
        "fwhm",
        "tmp_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is not in map_input
    map_hippocampus(
        "bids_id",
        "not_in_map_input",
        "resol",
        "label",
        "workbench_path",
        "subject_hippunfold_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_hip",
        "fwhm",
        "tmp_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when required files do not exist
    mock_os_path_isfile.return_value = False
    map_hippocampus(
        "bids_id",
        "feat",
        "resol",
        "label",
        "workbench_path",
        "subject_hippunfold_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_hip",
        "fwhm",
        "tmp_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_not_called()

    # Test when output file is not created successfully
    mock_os_path_isfile.side_effect = [True, True, False]
    map_hippocampus(
        "bids_id",
        "feat",
        "resol",
        "label",
        "workbench_path",
        "subject_hippunfold_dir",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_hip",
        "fwhm",
        "tmp_dir",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()


# Test map_cortex function
@patch("os.path.isfile")
@patch("os.makedirs")
@patch("subprocess.run")
def test_map_cortex(mock_subprocess_run, mock_os_makedirs, mock_os_path_isfile):
    # Setup
    mock_os_path_isfile.return_value = True
    mock_subprocess_run.return_value = None
    mock_os_makedirs.return_value = None

    # Test when feat starts with "plugin-"
    map_cortex(
        "bids_id",
        "plugin-test",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is in map_feat
    map_cortex(
        "bids_id",
        "thickness",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is not in map_feat
    map_cortex(
        "bids_id",
        "not_in_map_feat",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is "thickness"
    map_cortex(
        "bids_id",
        "thickness",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when feat is not "thickness"
    map_cortex(
        "bids_id",
        "not_thickness",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()

    # Test when required files do not exist
    mock_os_path_isfile.return_value = False
    map_cortex(
        "bids_id",
        "feat",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_not_called()

    # Test when output file is not created successfully
    mock_os_path_isfile.side_effect = [True, True, False]
    map_cortex(
        "bids_id",
        "feat",
        "resol",
        "label",
        "fwhm",
        "workbench_path",
        "subject_micapipe_dir",
        "subject_output_dir",
        "folder_maps",
        "folder_ctx",
        "subject_plugin_dir",
    )
    # Assertions
    mock_os_makedirs.assert_called()
    mock_subprocess_run.assert_called()
