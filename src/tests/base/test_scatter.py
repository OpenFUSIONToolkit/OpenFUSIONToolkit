#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for scatter plans (test_scatter.F90): forward and reverse,
# insert and add, real and complex values, unique and repeated indices, masks,
# local weights and the local part alone, on 1, 2 and 4 processes
#------------------------------------------------------------------------------
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/
"""

# Checks per process: for each of the 4 plans (unique/repeated indices, with/without
# local weights) one uniqueness flag check, then for each direction and mode 4
# communicating cases (mask, nolocal) and the local part, on real and complex values
NCHECK = 4*(1 + 2*2*(4 + 1)*2)


def run_scatter(nprocs):
    os.chdir(test_dir)
    if os.path.isfile('scatter.results'):
        os.remove('scatter.results')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template)
    if not run_OFT("./test_scatter", nprocs, 60):
        return False
    with open('scatter.results', 'r') as fid:
        nfail = int(fid.readline())
        ncheck = int(fid.readline())
    assert ncheck == NCHECK
    return nfail == 0


@pytest.mark.coverage
@pytest.mark.parametrize('nprocs', [pytest.param(1, id='serial'),
                                    pytest.param(2, id='mpi2', marks=pytest.mark.mpi),
                                    pytest.param(4, id='mpi4', marks=pytest.mark.mpi)])
def test_scatter(nprocs):
    '''Scatter plan operations against values computed directly from the plan lists'''
    assert run_scatter(nprocs)
