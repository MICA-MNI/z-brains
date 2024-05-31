# FILEPATH: /c:/Users/Ian/Documents/GitHub/z-brains/tests/test_clinical_reports.py
from unittest import mock
import pytest
import os
import sys
import pandas as pd

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
from brainspace.vtk_interface.wrappers.data_object import BSPolyData

from src.functions.clinical_reports import (
    adjectivize_struct,
    plot_surfs,
    convert_html_to_pdf,
    report_header_template,
    report_colors,
    feature_header_template,
    report_1x2_table,
    map_subcortical_vertices,
    DATA_PATH,
    _load_surfaces_ctx,
    _load_surfaces_hip,
    _load_data_sctx,
    load_data_struct,
    report_struct,
    generate_clinical_report,
)
import numpy as np
from unittest.mock import patch, MagicMock, mock_open, PropertyMock
from xhtml2pdf import pisa
from scipy.spatial.transform import Rotation
from pathlib import Path
import tempfile

# Test generate_clinical_report


@patch("src.functions.clinical_reports.get_bids_id")
@patch("src.functions.clinical_reports.get_subject_dir")
@patch("glob.glob")
def test_generate_clinical_report(mock_glob, mock_get_subject_dir, mock_get_bids_id):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Mocking the dependencies
        mock_get_bids_id.return_value = "dummy_bids_id"
        mock_get_subject_dir.return_value = tmpdir
        mock_glob.return_value = [
            "analysis-asymmetrey-feature-ADC-dummy_file1_dummy",
            "analysis-asymemtrey-feature-ADC-dummy_file2_dummy",
        ]

        # Test when analyses and features are None
        generate_clinical_report(
            zbrains_path="/dummy/zbrains/path",
            sid="dummy_sid",
            ses="dummy_ses",
            approach="zscore",
            tmp_dir="/dummy/tmp_dir",
        )
        mock_get_bids_id.assert_called_once_with("dummy_sid", "dummy_ses")
        mock_get_subject_dir.assert_called_once_with(
            "/dummy/zbrains/path", "dummy_sid", "dummy_ses"
        )
        mock_glob.assert_called_once_with(f"{tmpdir}/norm-z/*/*")

        # Test when analyses and features are not None
        mock_get_bids_id.reset_mock()
        mock_get_subject_dir.reset_mock()
        mock_get_subject_dir.return_value = tmpdir
        mock_glob.reset_mock()
        generate_clinical_report(
            zbrains_path="/dummy/zbrains/path",
            sid="dummy_sid",
            ses="dummy_ses",
            approach="zscore",
            analyses=["regional", "asymmetry"],
            features=["ADC", "FA"],
            tmp_dir="/dummy/tmp_dir",
        )
        mock_get_bids_id.assert_called_with("dummy_sid", "dummy_ses")
        mock_get_subject_dir.assert_called_once_with(
            "/dummy/zbrains/path", "dummy_sid", "dummy_ses"
        )


# Test report_struct


@patch("src.functions.clinical_reports.get_analysis_path_from_template")
@patch("src.functions.utils_analysis.get_bids_id")
@patch("src.functions.clinical_reports.make_png_missing")
@patch("src.functions.clinical_reports.make_png")
@patch("src.functions.clinical_reports.load_data_struct")
@patch("logging.getLogger")
def test_report_struct(
    mock_getLogger,
    mock_load_data_struct,
    mock_make_png,
    mock_make_png_missing,
    mock_get_bids_id,
    mock_get_analysis_path_from_template,
):
    # Mocking the dependencies
    mock_get_bids_id.return_value = "dummy_bids_id"
    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = False
    mock_get_analysis_path_from_template.return_value = mock_path
    mock_load_data_struct.return_value = (None, None)
    mock_make_png.return_value = "dummy_png_block"
    mock_make_png_missing.return_value = "dummy_png_missing_block"
    mock_logger = MagicMock()
    mock_getLogger.return_value = mock_logger

    # Test when file does not exist

    html = report_struct(
        struct="cortex",
        path_analysis="/dummy/analysis/path",
        sid="dummy_sid",
        analysis="asymmetry",
        approach="zscore",
        feat="ADC",
        tmp_dir="/dummy/tmp_dir",
    )
    mock_logger.warning.assert_called_once()
    assert html[-23:] == "dummy_png_missing_block"

    # Test when file exists and analysis is not "asymmetry"
    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True
    mock_get_analysis_path_from_template.return_value = mock_path

    mock_logger.warning.reset_mock()
    html = report_struct(
        struct="cortex",
        path_analysis="/dummy/analysis/path",
        sid="dummy_sid",
        analysis="regional",
        approach="zscore",
        feat="ADC",
        tmp_dir="/dummy/tmp_dir",
    )
    mock_logger.warning.assert_not_called()
    assert html[-15:] == "dummy_png_block"

    # Test when file exists and analysis is "asymmetry"
    html = report_struct(
        struct="cortex",
        path_analysis="/dummy/analysis/path",
        sid="dummy_sid",
        analysis="asymmetry",
        approach="zscore",
        feat="ADC",
        tmp_dir="/dummy/tmp_dir",
    )
    assert html[-15:] == "dummy_png_block"


# Test load_data_struct


@patch("nibabel.load")
@patch("src.functions.clinical_reports._load_data_sctx")
def test_load_data_struct_subcortex(mock_load_data_sctx, mock_nib_load):
    mock_load_data_sctx.return_value = (
        np.array([1, 2, 3], dtype=np.float64),
        np.array([4, 5, 6], dtype=np.float64),
    )
    data_lh, data_rh = load_data_struct(
        "subcortex",
        file_lh="dummy_path",
        analysis="regional",
        threshold=0.5,
        threshold_alpha=0.5,
    )
    mock_load_data_sctx.assert_called_once_with(
        "dummy_path", "regional", threshold=0.5, threshold_alpha=0.5
    )
    np.testing.assert_array_equal(data_lh, np.array([1, 2, 3]))
    np.testing.assert_array_equal(data_rh, np.array([4, 5, 6]))


@patch("nibabel.load")
def test_load_data_struct_cortex_asymmetry(mock_nib_load):
    mock_data = MagicMock()
    mock_data.data = np.array([-1, 0, 1], dtype=np.float64)
    mock_nib_load.return_value.darrays = [mock_data]
    data_lh, data_rh = load_data_struct(
        "cortex",
        file_lh="dummy_path",
        file_rh="dummy_path",
        analysis="asymmetry",
        threshold=0.5,
        threshold_alpha=0.5,
    )
    np.testing.assert_array_equal(data_lh, np.array([0, 0, 1]))
    np.testing.assert_array_equal(data_rh, np.array([-1, 0, 0]))


@patch("nibabel.load")
def test_load_data_struct_cortex_regional(mock_nib_load):
    mock_data = MagicMock()
    mock_data.data = np.array([-1, 0, 1], dtype=np.float64)
    mock_nib_load.return_value.darrays = [mock_data]
    data_lh, data_rh = load_data_struct(
        "cortex",
        file_lh="dummy_path",
        file_rh="dummy_path",
        analysis="regional",
        threshold=0.5,
        threshold_alpha=0.5,
    )
    np.testing.assert_array_equal(data_lh, np.array([-1, 0, 1]))
    np.testing.assert_array_equal(data_rh, np.array([-1, 0, 1]))


# Test _load_data_sctx


@patch("src.functions.clinical_reports.map_subcortical_vertices")
@patch("pandas.read_csv")
def test_load_data_sctx(mock_read_csv, mock_map_subcortical_vertices):
    # Mocking the read_csv function
    mock_read_csv.return_value = pd.DataFrame(
        {
            "column_name": [
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                1.0,
                1.1,
                1.2,
                1.3,
                1.4,
            ]
        }
    )
    # Mocking the map_subcortical_vertices function
    mock_map_subcortical_vertices.return_value = np.array(list(range(51820)))

    # Test with analysis != "asymmetry" and threshold is not None
    data_lh, data_rh = _load_data_sctx("dummy_path", "regional", 0.5, 0.5)
    mock_read_csv.assert_called_once_with("dummy_path", header=[0], index_col=0)
    np.testing.assert_array_equal(data_lh, np.array(list(range(25910))))
    np.testing.assert_array_equal(data_rh, np.array(list(range(25910, 51820))))

    # Test with analysis == "asymmetry" and threshold is not None
    data_lh, data_rh = _load_data_sctx("dummy_path", "asymmetry", 0.5, 0.5)
    np.testing.assert_array_equal(data_lh, np.array(list(range(25910))))
    np.testing.assert_array_equal(data_rh, np.array(list(range(25910, 51820))))

    # Test with analysis != "asymmetry" and threshold is None
    data_lh, data_rh = _load_data_sctx("dummy_path", "regional", None, 0.5)
    np.testing.assert_array_equal(data_lh, np.array(list(range(25910))))
    np.testing.assert_array_equal(data_rh, np.array(list(range(25910, 51820))))

    # Test with analysis == "asymmetry" and threshold is None
    data_lh, data_rh = _load_data_sctx("dummy_path", "asymmetry", None, 0.5)
    np.testing.assert_array_equal(data_lh, np.array(list(range(25910))))
    np.testing.assert_array_equal(data_rh, np.array(list(range(25910, 51820))))


# Test _load_surfaces_hip
def test_load_surfaces_hip_high_resolution(var):
    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip("high")

    assert isinstance(lat_lh, BSPolyData)
    assert isinstance(mid_lh, BSPolyData)
    assert isinstance(unf_lh, BSPolyData)
    assert isinstance(unf_rh, BSPolyData)
    assert isinstance(mid_rh, BSPolyData)
    assert isinstance(lat_rh, BSPolyData)
    assert len(lat_lh.Points) == 419
    assert len(mid_lh.Points) == 419
    assert len(unf_lh.Points) == 419
    assert len(unf_rh.Points) == 419
    assert len(mid_rh.Points) == 419
    assert len(lat_rh.Points) == 419


def test_load_surfaces_hip_low_resolution(var):

    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip("low")

    assert isinstance(lat_lh, BSPolyData)
    assert isinstance(mid_lh, BSPolyData)
    assert isinstance(unf_lh, BSPolyData)
    assert isinstance(unf_rh, BSPolyData)
    assert isinstance(mid_rh, BSPolyData)
    assert isinstance(lat_rh, BSPolyData)
    assert len(lat_lh.Points) == 7262
    assert len(mid_lh.Points) == 7262
    assert len(unf_lh.Points) == 7262
    assert len(unf_rh.Points) == 7262
    assert len(mid_rh.Points) == 7262
    assert len(lat_rh.Points) == 7262


# Test _load_surfaces_ctx


@pytest.fixture
def var(mocker):
    return mocker.patch(
        "src.functions.clinical_reports.DATA_PATH", new="src/data", autospec=False
    )


def test_load_surfaces_ctx_high_resolution(var):

    inf_lh, inf_rh, mask = _load_surfaces_ctx("high")
    assert isinstance(inf_lh, BSPolyData)
    assert isinstance(inf_rh, BSPolyData)
    assert len(mask) == 64984


# Test map_subcortical_vertices


def test_map_subcortical_vertices_correct_input():
    x = np.arange(16)
    result = map_subcortical_vertices(x)
    assert result.shape == (51278,)
    assert np.array_equal(result[:867], np.zeros(867))
    assert np.array_equal(result[867:2286], np.ones(1419))
    assert np.array_equal(result[-7180:], np.full(7180, 15))


def test_map_subcortical_vertices_incorrect_input_shape():
    x = np.arange(15)
    with pytest.raises(ValueError) as excinfo:
        map_subcortical_vertices(x)
    assert str(excinfo.value) == "Input data must have 16 values."


# Test feature_header_template


def test_feature_header_template_single_feature_no_extra():
    result = feature_header_template("ADC")
    assert (
        result
        == '<p style="border:0px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica,sans-serif;font-size:14px;text-align:center;color:#5d5070"><b> Feature: ADC  </b></p>'
    )


def test_feature_header_template_single_feature_with_extra():
    result = feature_header_template("ADC", extra="Extra")
    assert (
        result
        == '<p style="border:0px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica,sans-serif;font-size:14px;text-align:center;color:#5d5070"><b> Feature: ADC Extra </b></p>'
    )


def test_feature_header_template_multiple_features_no_extra():
    result = feature_header_template(["ADC", "FA"])
    assert (
        result
        == '<p style="border:0px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica,sans-serif;font-size:14px;text-align:center;color:#5d5070"><b> Features: ADC & FA  </b></p>'
    )


def test_feature_header_template_multiple_features_with_extra():
    result = feature_header_template(["ADC", "FA"], extra="Extra")
    assert (
        result
        == '<p style="border:0px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica,sans-serif;font-size:14px;text-align:center;color:#5d5070"><b> Features: ADC & FA Extra </b></p>'
    )


# Test report_colors


def test_report_colors_regional():
    result = report_colors("regional")
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#b31b2c"><b> Red </b> = INCREASED compared to controls </p>'
        in result
    )
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#13365d"><b> Blue </b> = DECREASED compared to controls </p>'
        in result
    )


def test_report_colors_not_regional():
    result = report_colors("asymmetry")
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#7532a8"><b> Purple </b> = <b>right</b> MORE THAN <b>left</b> </p>'
        in result
    )
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#32a852"><b> Green </b> = <b>left</b> MORE THAN <b>right</b> </p>'
        in result
    )


def test_report_colors_default():
    result = report_colors()
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#b31b2c"><b> Red </b> = INCREASED compared to controls </p>'
        in result
    )
    assert (
        '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#13365d"><b> Blue </b> = DECREASED compared to controls </p>'
        in result
    )


# Test report_header_template


def test_report_header_template_no_optional_params():
    result = report_header_template(sid="123", analysis="asymmetry")
    assert "<b>Subject</b>: 123," in result
    assert "<b>Sex</b>: n/a" in result
    assert "<b>Age</b>: n/a" in result


def test_report_header_template_with_ses():
    result = report_header_template(sid="123", ses="001", analysis="regional")
    assert "<b>Subject</b>: 123," in result
    assert "<b>Session</b>: 001," in result


def test_report_header_template_with_age():
    result = report_header_template(sid="123", age=30.5, analysis="regional")
    assert "<b>Age</b>: 30.5" in result


def test_report_header_template_with_sex():
    result = report_header_template(sid="123", sex="M", analysis="asymmetry")
    print(result)
    assert "<b>Sex</b>: M," in result


def test_report_header_template_with_all_params():
    result = report_header_template(
        sid="123", ses="001", age=30.5, sex="M", analysis="asymmetry"
    )
    assert "<b>Subject</b>: 123," in result
    assert "<b>Session</b>: 001," in result
    assert "<b>Sex</b>: M," in result
    assert "<b>Age</b>: 30.5" in result


# Test convert_html_to_pdf


def test_convert_html_to_pdf_success():
    with patch("builtins.open", mock_open()) as mocked_file, patch.object(
        pisa, "CreatePDF", return_value=MockPisaStatus(err=0)
    ) as mocked_pisa:
        assert convert_html_to_pdf("<html></html>", "output.pdf") == 0
        mocked_file.assert_called_once_with("output.pdf", "w+b")
        mocked_pisa.assert_called_once_with(
            "<html></html>", dest=mocked_file.return_value
        )


def test_convert_html_to_pdf_failure():
    with patch("builtins.open", mock_open()) as mocked_file, patch.object(
        pisa, "CreatePDF", return_value=MockPisaStatus(err=1)
    ) as mocked_pisa:
        assert convert_html_to_pdf("<html></html>", "output.pdf") == 1
        mocked_file.assert_called_once_with("output.pdf", "w+b")
        mocked_pisa.assert_called_once_with(
            "<html></html>", dest=mocked_file.return_value
        )


class MockPisaStatus:
    def __init__(self, err):
        self.err = err


# Test plot_surfs


@pytest.fixture
def mock_surface():
    mock = MagicMock()
    mock.append_array = MagicMock()
    return mock


# Test adjectivize_struct
def test_adjectivize_struct_cortex():
    assert adjectivize_struct("cortex") == "Cortical"


def test_adjectivize_struct_subcortex():
    assert adjectivize_struct("subcortex") == "Subcortical"


def test_adjectivize_struct_hippocampus():
    assert adjectivize_struct("hippocampus") == "Hippocampal"


def test_adjectivize_struct_invalid():
    with pytest.raises(ValueError):
        adjectivize_struct("invalid_structure")  # type: ignore
