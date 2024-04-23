import os
import tempfile
import shutil
from contextlib import contextmanager
from zbrains import tempdir

def test_tempdir_creation():
    with tempdir(tempfile.gettempdir(), 'test_prefix') as temp_dir:
        assert os.path.exists(temp_dir), "Temp directory should exist after creation"
        assert temp_dir.startswith(tempfile.gettempdir()), "Temp directory should be created in the specified directory"
        assert temp_dir.endswith('test_prefix'), "Temp directory should have the specified prefix"

def test_tempdir_cleanup():
    with tempdir(tempfile.gettempdir(), 'test_prefix') as temp_dir:
        pass
    assert not os.path.exists(temp_dir), "Temp directory should be cleaned up after use"

def test_tempdir_cleanup_on_error():
    try:
        with tempdir(tempfile.gettempdir(), 'test_prefix') as temp_dir:
            raise Exception("Test exception")
    except Exception:
        pass
    assert not os.path.exists(temp_dir), "Temp directory should be cleaned up after use even if an error occurs"