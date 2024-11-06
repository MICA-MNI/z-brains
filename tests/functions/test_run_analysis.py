import pytest
from unittest.mock import patch, MagicMock
from src.functions.run_analysis import main


@patch("src.functions.clinical_reports.generate_clinical_report")
@patch("src.functions.utils_analysis.run_analysis")
@patch("src.functions.utils_analysis.load_demo")
@patch("src.functions.utils_analysis.approach_to_folder")
def test_main_mismatched_zbrains_demo_refs(
    mock_approach_to_folder,
    mock_load_demo,
    mock_run_analysis,
    mock_generate_clinical_reports,
):
    logger = MagicMock()
    with pytest.raises(
        ValueError,
        match="The number of values provided with --zbrains_ref and --demo_ref must be the same.",
    ):
        main(
            dataset=None,
            micapipename=None,
            hippunfoldname=None,
            n_jobs_wb=None,
            workbench_path=None,
            zbrains_ref=["ref1"],
            demo_ref=["ref1", "ref2"],
            column_map={},
            subject_id="sub-001",
            session=None,
            demo=None,
            zbrains=None,
            struct=None,
            feat=None,
            normative=None,
            deconfound=None,
            smooth_ctx=None,
            smooth_hip=None,
            threshold=None,
            approach=None,
            resolution=None,
            labels_ctx=None,
            labels_hip=None,
            tmp=None,
            logger=logger,
            n_jobs=None,
        )


@patch("src.functions.clinical_reports.generate_clinical_report")
@patch("src.functions.utils_analysis.run_analysis")
@patch("src.functions.utils_analysis.load_demo")
@patch("src.functions.utils_analysis.approach_to_folder")
def test_main_unknown_column_names(
    mock_approach_to_folder,
    mock_load_demo,
    mock_run_analysis,
    mock_generate_clinical_reports,
):
    logger = MagicMock()
    with pytest.raises(
        ValueError, match="Unknown column names: {'unknown'}. Allowed options are:"
    ):
        main(
            dataset=None,
            micapipename=None,
            hippunfoldname=None,
            n_jobs_wb=None,
            workbench_path=None,
            zbrains_ref=["ref1"],
            demo_ref=["ref1"],
            column_map={"unknown": "known"},
            subject_id="sub-001",
            session=None,
            demo=None,
            zbrains=None,
            struct=None,
            feat=None,
            normative=None,
            deconfound=None,
            smooth_ctx=None,
            smooth_hip=None,
            threshold=None,
            approach=None,
            resolution=None,
            labels_ctx=None,
            labels_hip=None,
            tmp=None,
            logger=logger,
            n_jobs=None,
        )


@patch("src.functions.clinical_reports.generate_clinical_report")
@patch("src.functions.utils_analysis.run_analysis")
@patch("src.functions.run_analysis.load_demo")
@patch("src.functions.utils_analysis.approach_to_folder")
def test_main_multiple_demo_rows_2(
    mock_approach_to_folder,
    mock_load_demo,
    mock_run_analysis,
    mock_generate_clinical_reports,
):
    mock_run_analysis.return_value = {
        "cortex": {"high": {"midthickness": "ADC"}},
        "subcortex": "ADC",
        "hippocampus": {"high": {"midthickness": "ADC"}},
    }
    mock_generate_clinical_reports.return_value = "clinical_reports"
    mock_approach_to_folder.return_value = "approach_folder"
    return_val = MagicMock()
    new_mock = MagicMock()
    new_mock.empty = True

    return_val.loc.__getitem__.return_value = new_mock

    mock_load_demo.return_value = return_val
    logger = MagicMock()
    demo = MagicMock()
    with pytest.raises(ValueError, match="Cannot find sub-001 in demographics"):
        main(
            dataset=None,
            micapipename=None,
            hippunfoldname=None,
            n_jobs_wb=None,
            workbench_path=None,
            zbrains_ref=["ref1"],
            demo_ref=["ref1"],
            column_map={},
            subject_id="sub-001",
            session=None,
            demo=demo,
            zbrains=None,
            struct=None,
            feat=None,
            normative=MagicMock(),
            deconfound=None,
            smooth_ctx=None,
            smooth_hip=None,
            threshold=1,
            approach=None,
            resolution=["high"],
            labels_ctx=["midthickness"],
            labels_hip=["midthickness"],
            tmp=None,
            logger=logger,
            n_jobs=None,
        )


@patch("src.functions.clinical_reports.generate_clinical_report")
@patch("src.functions.utils_analysis.run_analysis")
@patch("src.functions.run_analysis.load_demo")
@patch("src.functions.utils_analysis.approach_to_folder")
def test_main_multiple_demo_rows(
    mock_approach_to_folder,
    mock_load_demo,
    mock_run_analysis,
    mock_generate_clinical_reports,
):
    mock_run_analysis.return_value = {
        "cortex": {"high": {"midthickness": "ADC"}},
        "subcortex": "ADC",
        "hippocampus": {"high": {"midthickness": "ADC"}},
    }
    mock_generate_clinical_reports.return_value = "clinical_reports"
    mock_approach_to_folder.return_value = "approach_folder"
    return_val = MagicMock()
    new_mock = MagicMock()
    new_mock.shape.__getitem__.return_value = 2
    new_mock.empty = False

    return_val.loc.__getitem__.return_value = new_mock

    mock_load_demo.return_value = return_val
    logger = MagicMock()
    demo = MagicMock()
    with pytest.raises(
        ValueError, match="Provided sub-001 is not unique in demographics file."
    ):
        main(
            dataset=None,
            micapipename=None,
            hippunfoldname=None,
            n_jobs_wb=None,
            workbench_path=None,
            zbrains_ref=["ref1"],
            demo_ref=["ref1"],
            column_map={},
            subject_id="sub-001",
            session=None,
            demo=demo,
            zbrains=None,
            struct=None,
            feat=None,
            normative=MagicMock(),
            deconfound=None,
            smooth_ctx=None,
            smooth_hip=None,
            threshold=1,
            approach=None,
            resolution=["high"],
            labels_ctx=["midthickness"],
            labels_hip=["midthickness"],
            tmp=None,
            logger=logger,
            n_jobs=None,
        )
