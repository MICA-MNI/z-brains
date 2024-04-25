import pytest
import logging
import os
import sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from unittest.mock import patch, MagicMock
from functions.run_analysis import run, main


def test_main_mismatched_zbrains_demo_refs():
    logger = MagicMock()
    with pytest.raises(
        ValueError,
        match="The number of values provided with --zbrains_ref and --demo_ref must be the same.",
    ):
        main(
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


def test_main_unknown_column_names():
    logger = MagicMock()
    with pytest.raises(
        ValueError, match="Unknown column names: {'unknown'}. Allowed options are:"
    ):
        main(
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


def test_main_no_demo_row():
    logger = MagicMock()
    demo = MagicMock()
    demo.empty = True
    main(
        zbrains_ref=["ref1"],
        demo_ref=["ref1"],
        column_map={},
        subject_id="sub-001",
        session=None,
        demo=demo,
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
    logger.warning.assert_called_once()


def test_main_multiple_demo_rows():
    logger = MagicMock()
    demo = MagicMock()
    demo.empty = False
    demo.shape = (2, 2)
    with pytest.raises(
        ValueError, match="Provided sub-001 is not unique in demographics file."
    ):
        main(
            zbrains_ref=["ref1"],
            demo_ref=["ref1"],
            column_map={},
            subject_id="sub-001",
            session=None,
            demo=demo,
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


def test_run():
    # Mock the main function
    with patch("functions.run_analysis.main") as mock_main:
        # Mock the logger and its methods
        mock_logger = MagicMock()
        with patch("logging.getLogger", return_value=mock_logger):
            # Mock the console handler and its methods
            mock_console_handler = MagicMock()
            with patch("logging.StreamHandler", return_value=mock_console_handler):
                # Mock the file handler and its methods
                mock_file_handler = MagicMock()
                with patch("logging.FileHandler", return_value=mock_file_handler):
                    # Call the run function
                    run(
                        subject_id="sub-001",
                        zbrains="zbrains",
                        demo_ref="demo_ref",
                        zbrains_ref="zbrains_ref",
                        verbose=2,
                        filter_warnings=True,
                        logfile="logfile.log",
                    )

    # Check that the main function was called with the correct arguments
    mock_main.assert_called_once_with(
        "zbrains_ref",
        "demo_ref",
        None,
        "sub-001",
        None,
        None,
        "zbrains",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        mock_logger,
        None,
    )

    # Check that the logger was created with the correct name
    logging.getLogger.assert_called_once_with(None)

    # Check that the console handler was created and configured correctly
    mock_console_handler.assert_called_once_with(sys.stdout)
    mock_console_handler.setFormatter.assert_called_once()
    mock_console_handler.setLevel.assert_called_once_with(logging.INFO)
    mock_console_handler.addFilter.assert_called_once()

    # Check that the file handler was created and configured correctly
    mock_console_handler.assert_called_once_with("logfile.log", mode="w")
    mock_file_handler.setFormatter.assert_called_once()

    # Check that the logger was configured correctly
    mock_logger.setLevel.assert_called_once_with(logging.DEBUG)
    assert mock_logger.addHandler.call_count == 2


def test_run_unhandled_exception():
    # Mock the main function to raise an exception
    with patch("functions.run_analysis.main", side_effect=Exception):
        # Mock the logger and its methods
        mock_logger = MagicMock()
        with patch("logging.getLogger", return_value=mock_logger):
            # Mock the console handler and its methods
            mock_console_handler = MagicMock()
            with patch("logging.StreamHandler", return_value=mock_console_handler):
                # Call the run function and check for the exception
                with pytest.raises(Exception):
                    run(
                        subject_id="sub-001",
                        zbrains="zbrains",
                        demo_ref="demo_ref",
                        zbrains_ref="zbrains_ref",
                        verbose=2,
                        filter_warnings=True,
                        logfile="logfile.log",
                    )

    # Check that the logger was called with the correct arguments
    mock_logger.critical.assert_called_once()
