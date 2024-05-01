import os
import shutil
import fnmatch
from unittest import mock
import sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from src.functions.utilities import delete_temp_folders


# Test delete_temp_folders


def test_delete_temp_folders():
    with mock.patch("os.walk") as mock_walk, mock.patch("shutil.rmtree") as mock_rmtree:
        mock_walk.return_value = [
            ("root", ["z_brains_temp.12345678", "other_folder"], ["file1", "file2"]),
            ("root/z_brains_temp.12345678", [], ["file3", "file4"]),
            ("root/other_folder", [], ["file5", "file6"]),
        ]

        delete_temp_folders("root")

        mock_rmtree.assert_called_once_with("root\\z_brains_temp.12345678")


def test_delete_temp_folders_no_match():
    with mock.patch("os.walk") as mock_walk, mock.patch("shutil.rmtree") as mock_rmtree:
        mock_walk.return_value = [
            ("root", ["other_folder"], ["file1", "file2"]),
            ("root/other_folder", [], ["file5", "file6"]),
        ]

        delete_temp_folders("root")

        mock_rmtree.assert_not_called()


def test_delete_temp_folders_multiple_matches():
    with mock.patch("os.walk") as mock_walk, mock.patch("shutil.rmtree") as mock_rmtree:
        mock_walk.return_value = [
            (
                "root",
                ["z_brains_temp.12345678", "z_brains_temp.87654321", "other_folder"],
                ["file1", "file2"],
            ),
            ("root/z_brains_temp.12345678", [], ["file3", "file4"]),
            ("root/z_brains_temp.87654321", [], ["file7", "file8"]),
            ("root/other_folder", [], ["file5", "file6"]),
        ]

        delete_temp_folders("root")

        assert mock_rmtree.call_count == 2
        mock_rmtree.assert_any_call("root\\z_brains_temp.12345678")
        mock_rmtree.assert_any_call("root\\z_brains_temp.87654321")
