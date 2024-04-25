# FILEPATH: /c:/Users/Ian/Documents/GitHub/z-brains/tests/test_clinical_reports.py
import pytest
import os
import sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from functions.clinical_reports import (
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
from unittest.mock import patch, MagicMock, mock_open
from xhtml2pdf import pisa
from scipy.spatial.transform import Rotation
from pathlib import Path

# Test generate_clinical_report


@patch("functions.utils_analysis.get_bids_id")
@patch("functions.utils_analysis.get_subject_dir")
@patch("glob.glob")
def test_generate_clinical_report(mock_glob, mock_get_subject_dir, mock_get_bids_id):
    # Mocking the dependencies
    mock_get_bids_id.return_value = "dummy_bids_id"
    mock_get_subject_dir.return_value = "/dummy/subject/dir"
    mock_glob.return_value = ["dummy_file1", "dummy_file2"]

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
    mock_glob.assert_called_once_with("/dummy/subject/dir/zscore/*/*")

    # Test when analyses and features are not None
    mock_get_bids_id.reset_mock()
    mock_get_subject_dir.reset_mock()
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
    mock_get_bids_id.assert_called_once_with("dummy_sid", "dummy_ses")
    mock_get_subject_dir.assert_called_once_with(
        "/dummy/zbrains/path", "dummy_sid", "dummy_ses"
    )
    mock_glob.assert_not_called()


# Test report_struct


@patch("functions.utils_analysis.get_analysis_path_from_template")
@patch("functions.utils_analysis.get_bids_id")
@patch("functions.clinical_reports.make_png_missing")
@patch("functions.clinical_reports.make_png")
@patch("functions.clinical_reports.load_data_struct")
@patch("builtins.logging")
def test_report_struct(
    mock_logging,
    mock_load_data_struct,
    mock_make_png,
    mock_make_png_missing,
    mock_get_bids_id,
    mock_get_analysis_path_from_template,
):
    # Mocking the dependencies
    mock_get_bids_id.return_value = "dummy_bids_id"
    mock_get_analysis_path_from_template.return_value = Path("/dummy/path")
    mock_load_data_struct.return_value = (None, None)
    mock_make_png.return_value = "dummy_png_block"
    mock_make_png_missing.return_value = "dummy_png_missing_block"
    mock_logger = MagicMock()
    mock_logging.getLogger.return_value = mock_logger

    # Test when file does not exist
    mock_get_analysis_path_from_template.return_value.exists.return_value = False
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
    assert html == "dummy_png_missing_block"

    # Test when file exists and analysis is not "asymmetry"
    mock_get_analysis_path_from_template.return_value.exists.return_value = True
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
    assert html == "dummy_png_block"

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
    assert html == "dummy_png_block"


# Test load_data_struct


@patch("nib.load")
@patch("functions.clinical_reports._load_data_sctx")
def test_load_data_struct_subcortex(mock_load_data_sctx, mock_nib_load):
    mock_load_data_sctx.return_value = (np.array([1, 2, 3]), np.array([4, 5, 6]))
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


@patch("nib.load")
def test_load_data_struct_cortex_asymmetry(mock_nib_load):
    mock_data = MagicMock()
    mock_data.data = np.array([-1, 0, 1])
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


@patch("nib.load")
def test_load_data_struct_cortex_regional(mock_nib_load):
    mock_data = MagicMock()
    mock_data.data = np.array([-1, 0, 1])
    mock_nib_load.return_value.darrays = [mock_data]
    data_lh, data_rh = load_data_struct(
        "cortex",
        file_lh="dummy_path",
        file_rh="dummy_path",
        analysis="regional",
        threshold=0.5,
        threshold_alpha=0.5,
    )
    np.testing.assert_array_equal(data_lh, np.array([-0.5, 0, 1]))
    np.testing.assert_array_equal(data_rh, np.array([-0.5, 0, 1]))


# Test _load_data_sctx


@patch("pd.read_csv")
@patch("functions.clinical_reports.map_subcortical_vertices")
def test_load_data_sctx(mock_map_subcortical_vertices, mock_read_csv):
    # Mocking the read_csv function
    mock_read_csv.return_value.to_numpy.return_value.ravel.return_value = np.array(
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
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
@patch("brainspace.mesh.mesh_io.read_surface")
@patch("functions.utils_analysis.map_resolution")
def test_load_surfaces_hip_high_resolution(mock_map_resolution, mock_read_surface):
    mock_map_resolution.return_value = "high"
    mock_read_surface.return_value = MagicMock()
    mock_read_surface.return_value.Points = np.array([[1, 2, 3]])

    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip("high")

    mock_map_resolution.assert_called_once_with("hippocampus", "high")
    assert mock_read_surface.call_count == 6
    mock_read_surface.assert_any_call(
        f"{DATA_PATH}/tpl-avg_space-canonical_den-0p5mm_label-hipp_midthickness.surf.gii"
    )
    mock_read_surface.assert_any_call(
        f"{DATA_PATH}/tpl-avg_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii"
    )

    # Check transformations
    assert np.array_equal(mid_lh.Points, np.array([[-1, 2, 3]]))  # Flipped
    assert np.array_equal(
        unf_rh.Points,
        Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        unf_lh.Points,
        Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[-1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        unf_rh.Points,
        Rotation.from_rotvec(np.pi * np.array([0, 0, 1])).apply(unf_rh.Points),
    )  # Rotated
    assert np.array_equal(
        lat_rh.Points,
        Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        lat_lh.Points,
        Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[-1, 2, 3]])
        ),
    )  # Rotated


@patch("brainspace.mesh.mesh_io.read_surface")
@patch("functions.utils_analysis.map_resolution")
def test_load_surfaces_hip_low_resolution(mock_map_resolution, mock_read_surface):
    mock_map_resolution.return_value = "low"
    mock_read_surface.return_value = MagicMock()
    mock_read_surface.return_value.Points = np.array([[1, 2, 3]])

    lat_lh, mid_lh, unf_lh, unf_rh, mid_rh, lat_rh = _load_surfaces_hip("low")

    mock_map_resolution.assert_called_once_with("hippocampus", "low")
    assert mock_read_surface.call_count == 6
    mock_read_surface.assert_any_call(
        f"{DATA_PATH}/tpl-avg_space-canonical_den-2mm_label-hipp_midthickness.surf.gii"
    )
    mock_read_surface.assert_any_call(
        f"{DATA_PATH}/tpl-avg_space-unfold_den-2mm_label-hipp_midthickness.surf.gii"
    )

    # Check transformations
    assert np.array_equal(mid_lh.Points, np.array([[-1, 2, 3]]))  # Flipped
    assert np.array_equal(
        unf_rh.Points,
        Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        unf_lh.Points,
        Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[-1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        unf_rh.Points,
        Rotation.from_rotvec(np.pi * np.array([0, 0, 1])).apply(unf_rh.Points),
    )  # Rotated
    assert np.array_equal(
        lat_rh.Points,
        Rotation.from_rotvec(3 * np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[1, 2, 3]])
        ),
    )  # Rotated
    assert np.array_equal(
        lat_lh.Points,
        Rotation.from_rotvec(np.pi / 2 * np.array([0, 1, 0])).apply(
            np.array([[-1, 2, 3]])
        ),
    )  # Rotated


# Test _load_surfaces_ctx


@patch("brainspace.mesh.mesh_io.read_surface")
@patch("functions.utils_analysis.map_resolution")
@patch("brainspace.datasets.load_mask")
def test_load_surfaces_ctx_high_resolution(
    mock_load_mask, mock_map_resolution, mock_read_surface
):
    mock_map_resolution.return_value = "high"
    mock_read_surface.return_value = "mock_surface"
    mock_load_mask.return_value = "mock_mask"

    inf_lh, inf_rh, mask = _load_surfaces_ctx("high")

    mock_map_resolution.assert_called_once_with("cortex", "high")
    mock_read_surface.assert_any_call(f"{DATA_PATH}/fsLR-32k.L.inflated.surf.gii")
    mock_read_surface.assert_any_call(f"{DATA_PATH}/fsLR-32k.R.inflated.surf.gii")
    mock_load_mask.assert_called_once_with(join=True)
    assert inf_lh == "mock_surface"
    assert inf_rh == "mock_surface"
    assert mask == "mock_mask"


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


def test_map_subcortical_vertices_incorrect_input_type():
    x = list(range(16))
    with pytest.raises(TypeError):
        map_subcortical_vertices(x)


# Test report_1x2_table


def test_report_1x2_table_default_height():
    fig1 = "path/to/fig1.png"
    fig2 = "path/to/fig2.png"
    result = report_1x2_table(fig1, fig2)
    expected = (
        '<table style="border:1px solid white;width:100%">'
        "<tr>"
        "<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center>"
        '<img style="height:250px;margin-top:-100px;" src="path/to/fig1.png">'
        "</td>"
        "<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center>"
        '<img style="height:250px;margin-top:-100px;" src="path/to/fig2.png">'
        "</td>"
        "</tr>"
        "</table>"
    )
    assert result == expected


def test_report_1x2_table_custom_height():
    fig1 = "path/to/fig1.png"
    fig2 = "path/to/fig2.png"
    result = report_1x2_table(fig1, fig2, height=300)
    expected = (
        '<table style="border:1px solid white;width:100%">'
        "<tr>"
        "<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center>"
        '<img style="height:300px;margin-top:-100px;" src="path/to/fig1.png">'
        "</td>"
        "<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center>"
        '<img style="height:300px;margin-top:-100px;" src="path/to/fig2.png">'
        "</td>"
        "</tr>"
        "</table>"
    )
    assert result == expected


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
    assert "Subject: 123," in result
    assert "Sex: n/a," in result
    assert "Age: n/a" in result
    assert "Test analysis" in result


def test_report_header_template_with_ses():
    result = report_header_template(sid="123", ses="001", analysis="regional")
    assert "Subject: 123, &nbsp; <b>Session</b>: 001," in result


def test_report_header_template_with_age():
    result = report_header_template(sid="123", age=30.5, analysis="regional")
    assert "Age: 30.5" in result


def test_report_header_template_with_sex():
    result = report_header_template(sid="123", sex="M", analysis="asymmetry")
    assert "Sex: M," in result


def test_report_header_template_with_all_params():
    result = report_header_template(
        sid="123", ses="001", age=30.5, sex="M", analysis="asymmetry"
    )
    assert "Subject: 123, &nbsp; <b>Session</b>: 001," in result
    assert "Sex: M," in result
    assert "Age: 30.5" in result
    assert "Test analysis" in result


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


@pytest.fixture
def mock_plot_surf():
    with patch("brainspace.plotting.surface_plotting.plot_surf", autospec=True) as mock:
        yield mock


def test_plot_surfs_empty(mock_surface, mock_plot_surf):
    plot_surfs([], [], views=None, size=None)
    mock_plot_surf.assert_not_called()


def test_plot_surfs_one_item(mock_surface, mock_plot_surf):
    plot_surfs([mock_surface], [np.array([1, 2, 3])], views=["dorsal"], size=None)
    mock_surface.append_array.assert_called_once_with(np.array([1, 2, 3]), name="surf1")
    mock_plot_surf.assert_called_once()


def test_plot_surfs_multiple_items(mock_surface, mock_plot_surf):
    plot_surfs(
        [mock_surface, mock_surface],
        [np.array([1, 2, 3]), np.array([4, 5, 6])],
        views=["dorsal", "lateral"],
        size=None,
    )
    assert mock_surface.append_array.call_count == 2
    mock_plot_surf.assert_called_once()


def test_plot_surfs_size_none(mock_surface, mock_plot_surf):
    plot_surfs([mock_surface], [np.array([1, 2, 3])], views=["dorsal"], size=None)
    mock_plot_surf.assert_called_once_with(
        {"surf1": mock_surface},
        [["surf1"]],
        array_name=["surf1"],
        view=["dorsal"],
        color_bar="bottom",
        color_range=(-2, 2),
        share="both",
        cmap="cmo.balance",
        zoom=1.75,
        size=(200, 350),
        transparent_bg=False,
    )


def test_plot_surfs_size_int(mock_surface, mock_plot_surf):
    plot_surfs([mock_surface], [np.array([1, 2, 3])], views=["dorsal"], size=500)
    mock_plot_surf.assert_called_once_with(
        {"surf1": mock_surface},
        [["surf1"]],
        array_name=["surf1"],
        view=["dorsal"],
        color_bar="bottom",
        color_range=(-2, 2),
        share="both",
        cmap="cmo.balance",
        zoom=1.75,
        size=500,
        transparent_bg=False,
    )


def test_plot_surfs_size_tuple(mock_surface, mock_plot_surf):
    plot_surfs([mock_surface], [np.array([1, 2, 3])], views=["dorsal"], size=(500, 600))
    mock_plot_surf.assert_called_once_with(
        {"surf1": mock_surface},
        [["surf1"]],
        array_name=["surf1"],
        view=["dorsal"],
        color_bar="bottom",
        color_range=(-2, 2),
        share="both",
        cmap="cmo.balance",
        zoom=1.75,
        size=(500, 600),
        transparent_bg=False,
    )


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
