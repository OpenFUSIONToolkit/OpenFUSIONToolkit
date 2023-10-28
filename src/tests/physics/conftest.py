import pytest

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: Mark test as too slow for CI testing")
    config.addinivalue_line("markers", "mpi: Mark test as parallel using MPI")
    config.addinivalue_line("markers", "linear: Mark MHD test as linear")
