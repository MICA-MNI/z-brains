from pathlib import Path
import os
import sys

import pytest

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from src.functions.utils_analysis import (
    get_session,
    get_bids_id,
    get_subject_dir,
    get_id,
    zscore,
    mahalanobis_distance,
    _load_data,
    _subject_zscore,
    _save,
    get_deconfounder,
    run_analysis,
)
import pandas as pd
import numpy as np
from numpy.testing import assert_array_almost_equal
from src.functions.utils_analysis import process_feature
from unittest.mock import MagicMock, patch
from src.functions.constants import ProcessingException


@patch("src.functions.utils_analysis.load_demo")
@patch("src.functions.utils_analysis.get_subject_dir")
@patch("src.functions.utils_analysis.Parallel")
@patch(
    "builtins.zip",
    return_value=([{"Test": "value"}], [{"info": "test", "warning": None}]),
)
def test_run_analysis(mock_zip, mock_parallel, mock_get_subject_dir, mock_load_demo):
    # Setup
    mock_get_subject_dir.return_value = "/path/to/subject_dir"
    mock_load_demo.return_value = [pd.DataFrame()]
    mock_parallel.return_value = MagicMock()

    # Call function
    run_analysis(
        px_sid="test_px_sid",
        px_ses="test_px_ses",
        cn_zbrains=["/path/to/cn_zbrains"],
        cn_demo_paths=["/path/to/cn_demo_paths"],
        px_zbrains="/path/to/px_zbrains",
        px_demo=pd.Series(),
        structures=["subcortex"],
        features=["ADC", "FA"],
        cov_normative=["cov1", "cov2"],
        cov_deconfound=["cov3", "cov4"],
        smooth_ctx=1.0,
        smooth_hip=1.0,
        resolutions=["high", "low"],
        labels_ctx=["label1", "label2"],
        labels_hip=["label3", "label4"],
        actual_to_expected={"actual": "expected"},
        analyses=["asymmetry", "regional"],
        approach="zscore",
        col_dtypes={"col1": str, "col2": int},
        tmp="tmp",
        n_jobs=1,
    )

    # Asserts
    mock_get_subject_dir.assert_called_once_with(
        "/path/to/px_zbrains", "test_px_sid", "test_px_ses"
    )
    mock_load_demo.assert_called_once_with(
        ["/path/to/cn_demo_paths"],
        rename={"actual": "expected"},
        dtypes={"col1": str, "col2": int},
        tmp="tmp",
    )
    mock_parallel.assert_called_once_with(n_jobs=1)


# Test process_feature


@patch("src.functions.utils_analysis._subject_zscore")
@patch("src.functions.utils_analysis._save")
@patch("src.functions.utils_analysis._load_one", return_value=None)
@patch(
    "src.functions.utils_analysis._load_data",
    return_value=(None, pd.DataFrame()),
)
@patch("src.functions.utils_analysis.get_deconfounder")
def test_process_feature_no_data_cn(
    mock_get_deconfounder,
    mock_load_data,
    mock_load_one,
    mock_subject_zscore,
    mock__save,
):
    log = process_feature(
        "feat",
        {},
        _load_data,
        [],
        "tmp",
        None,
        None,
        None,
        "px_sid",
        "px_ses",
        [],
        0,
        "struct",
        "pth_analysis",
        {"feat": {"warning": None, "info": None}},
    )
    assert (
        log[1]["warning"]
        == "\tfeat           : \tNo data available for reference subjects."
    )


@patch("src.functions.utils_analysis._subject_zscore")
@patch("src.functions.utils_analysis._save")
@patch("src.functions.utils_analysis._load_one", return_value=None)
@patch(
    "src.functions.utils_analysis._load_data",
    side_effect=[(np.array([1, 2, 3]), pd.DataFrame())],
)
@patch("src.functions.utils_analysis.get_deconfounder")
def test_process_feature_no_data_px(
    mock_get_deconfounder,
    mock_load_data,
    mock_load_one,
    mock_subject_zscore,
    mock__save,
):

    log = process_feature(
        "feat",
        {},
        None,
        [],
        "tmp",
        None,
        None,
        None,
        "px_sid",
        "px_ses",
        [],
        0,
        "struct",
        "pth_analysis",
        {"feat": {"warning": None, "info": None}},
    )
    print(log[1])
    assert (
        log[1]["warning"]
        == "\tfeat           : \tNo data available for target subject."
    )


@patch("src.functions.utils_analysis._subject_zscore")
@patch("src.functions.utils_analysis._save")
@patch("src.functions.utils_analysis._load_one")
@patch(
    "src.functions.utils_analysis._load_data",
    side_effect=[(np.array([1, 2, 3]), pd.DataFrame()), np.array([4, 5, 6])],
)
@patch("src.functions.utils_analysis.get_deconfounder")
def test_process_feature_with_data(
    mock_get_deconfounder,
    mock_load_data,
    mock_load_one,
    mock_subject_zscore,
    mock__save,
):
    mock_subject_zscore.return_value = {"analysis": np.array([7, 8, 9])}
    log = process_feature(
        "feat",
        {},
        _load_data,
        [],
        "tmp",
        None,
        None,
        _subject_zscore,
        "px_sid",
        "px_ses",
        ["analysis"],
        0,
        "struct",
        "pth_analysis",
        {"feat": {"warning": None, "info": None}},
    )
    assert log[1]["info"] == "\tfeat           : \t[0/0 reference subjects available]"
    mock__save.assert_called_once()


@patch("src.functions.utils_analysis._subject_zscore")
@patch("src.functions.utils_analysis._save")
@patch("src.functions.utils_analysis._load_one")
@patch(
    "src.functions.utils_analysis._load_data",
    side_effect=[(np.array([1, 2, 3]), pd.DataFrame()), np.array([4, 5, 6])],
)
@patch("src.functions.utils_analysis.get_deconfounder")
def test_process_feature_with_deconfounding(
    mock_get_deconfounder,
    mock_load_data,
    mock_load_one,
    mock_subject_zscore,
    mock__save,
):
    mock_load_one.return_value = MagicMock()
    _subject_zscore = MagicMock(return_value={"analysis": np.array([7, 8, 9])})
    _save = MagicMock()

    dec = MagicMock()
    mock_get_deconfounder.return_value = dec
    dec.fit_transform.return_value = np.array([10, 11, 12])
    dec.transform.return_value = np.array([13, 14, 15])
    mock_subject_zscore.return_value = np.array([-3, 0, 3])
    log = process_feature(
        "ADC",
        {},
        mock_load_data,
        [],
        "tmp",
        ["cov"],
        pd.Series(),
        _subject_zscore,
        "px_sid",
        "px_ses",
        ["regional"],
        0,
        "struct",
        "pth_analysis",
        {"ADC": {"warning": None, "info": None}},
    )
    assert log[1]["info"] == "\tADC            : \t[0/0 reference subjects available]"
    mock__save.assert_called_once()
    mock_get_deconfounder.assert_called_once_with(covariates=["cov"])
    dec.fit_transform.assert_called_once()
    dec.transform.assert_called_once()


# Test mahalanobis_distance


def test_mahalanobis_distance_2d():
    x_train = np.array([[10, 10, 9], [1, 2, 3], [1, 1, 1]])
    x_test = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    print(x_train.mean(axis=0))
    result = mahalanobis_distance(x_train, x_test)
    assert result.shape[0] == x_test.shape[0]


def test_mahalanobis_distance_with_zeros():
    x_train = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    x_test = np.array([[2, 3, 4], [2, 3, 4], [2, 3, 4]])
    with pytest.raises(ProcessingException):
        mahalanobis_distance(x_train, x_test)


def test_mahalanobis_distance_with_nan():
    x_train = np.array([[1, np.nan, 5], [1, np.nan, 5], [1, np.nan, 5]])
    x_test = np.array([[2, 3, 4], [2, 3, 4], [2, 3, 4]])
    result = mahalanobis_distance(x_train, x_test)
    assert result.shape[0] == x_test.shape[0]


# Test zscore


def test_zscore_1d():
    x_train = np.array([1, 2, 3, 4, 5])
    x_test = np.array([2, 3, 4, 5, 6])
    expected = np.array([-0.707107, 0.0, 0.707107, 1.414214, 2.12132])
    result = zscore(x_train, x_test)
    assert_array_almost_equal(result, expected)


def test_zscore_2d():
    x_train = np.array([[1, 2, 3], [4, 5, 6]])
    x_test = np.array([[2, 3, 4], [5, 6, 7]])
    expected = np.array(
        [[-0.333333, -0.333333, -0.333333], [1.666666, 1.666666, 1.666666]]
    )
    result = zscore(x_train, x_test)
    assert_array_almost_equal(result, expected)


def test_zscore_with_zeros():
    x_train = np.array([0, 0, 0, 0, 0])
    x_test = np.array([2, 3, 4, 5, 6])
    expected = np.array([0, 0, 0, 0, 0])
    result = zscore(x_train, x_test)
    assert_array_almost_equal(result, expected)


def test_zscore_with_nan():
    x_train = np.array([-3, 1, np.nan, 1, 3])
    x_test = np.array([1, 1, np.nan, 1, 1])
    expected = np.array([0.229416, 0.229416, np.nan, 0.229416, 0.229416])
    result = zscore(x_train, x_test)
    assert_array_almost_equal(result, expected)


# Test get_subject_dir
def test_get_subject_dir_with_session():
    assert get_subject_dir("/root/path", "1234", "1") == Path(
        "/root/path/sub-1234/ses-1"
    )
    assert get_subject_dir("/root/path", "sub-1234", "ses-1") == Path(
        "/root/path/sub-1234/ses-1"
    )


def test_get_subject_dir_without_session():
    assert get_subject_dir("/root/path", "1234") == Path("/root/path/sub-1234")
    assert get_subject_dir("/root/path", "sub-1234") == Path("/root/path/sub-1234")


def test_get_subject_dir_with_none_session():
    assert get_subject_dir("/root/path", "1234", None) == Path("/root/path/sub-1234")
    assert get_subject_dir("/root/path", "sub-1234", None) == Path(
        "/root/path/sub-1234"
    )


def test_get_subject_dir_with_empty_session():
    assert get_subject_dir("/root/path", "1234", "") == Path("/root/path/sub-1234")
    assert get_subject_dir("/root/path", "sub-1234", "") == Path("/root/path/sub-1234")


def test_get_bids_id_with_session():
    assert get_bids_id("1234", "1") == "sub-1234_ses-1"
    assert get_bids_id("sub-1234", "ses-1") == "sub-1234_ses-1"


def test_get_bids_id_without_session():
    assert get_bids_id("1234") == "sub-1234"
    assert get_bids_id("sub-1234") == "sub-1234"


def test_get_bids_id_with_none_session():
    assert get_bids_id("1234", None) == "sub-1234"
    assert get_bids_id("sub-1234", None) == "sub-1234"


# Test get_session


def test_get_session_with_prefix():
    assert get_session("ses-1234") == "ses-1234"
    assert get_session("1234") == "ses-1234"
    assert get_session(pd.NA) == None
    assert get_session("n/a") == None
    assert get_session("") == None


def test_get_session_without_prefix():
    assert get_session("ses-1234", add_predix=False) == "1234"
    assert get_session("1234", add_predix=False) == "1234"
    assert get_session(pd.NA, add_predix=False) == None
    assert get_session("n/a", add_predix=False) == None
    assert get_session("", add_predix=False) == None


# Test get_id


def test_get_id_with_prefix():
    assert get_id("sub-1234") == "sub-1234"
    assert get_id("1234") == "sub-1234"


def test_get_id_without_prefix():
    assert get_id("sub-1234", add_prefix=False) == "1234"
    assert get_id("1234", add_prefix=False) == "1234"
