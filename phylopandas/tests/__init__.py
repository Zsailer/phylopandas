import pytest
import os, glob

@pytest.fixture
def path_to_dat():
    """Get path to dat folder with test data."""
    # Get path to test directory.
    path_to_test = os.path.dirname(os.path.realpath(__file__))

    # Build path to dat
    path_to_dat = os.path.join(path_to_test, 'dat')
    return path_to_dat

@pytest.fixture()
def clean_dat(path_to_dat):
    yield clean_dat

    # Get files in dat_files
    dat_files = glob.glob(os.path.join(path_to_dat,"*"))

    # Remove files from dat folder that begin with 'test'.
    for datf in dat_files:
        path, f = os.path.split(datf)
        if f[:4] == 'test':
            os.remove(datf)
