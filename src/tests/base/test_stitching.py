# SPDX-License-Identifier: LGPL-3.0-only
import os
import sys

import pytest

test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT  # noqa: E402

oft_in = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/
"""


@pytest.mark.coverage
@pytest.mark.parametrize("nthreads", (1, 2))
@pytest.mark.parametrize("nproc", (1, pytest.param(2, marks=pytest.mark.mpi)))
def test_stitching_dot(nthreads, nproc, monkeypatch):
    """Check native dot products with local and distributed ownership.

    @param nthreads Number of OpenMP threads requested for the native executable
    @param nproc Number of MPI processes
    @param monkeypatch Pytest fixture restoring environment and working directory
    """
    monkeypatch.chdir(test_dir)
    monkeypatch.setenv("OMP_NUM_THREADS", str(nthreads))
    monkeypatch.setenv("OMP_DYNAMIC", "FALSE")
    with open('oft.in', 'w') as fid:
        fid.write(oft_in)
    assert run_OFT("./test_stitching", nproc, 20)
