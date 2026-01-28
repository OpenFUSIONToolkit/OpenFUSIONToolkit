from __future__ import print_function
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

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={0}
 nbase={1}
 grid_order=1
/

&cube_options
 mesh_type=1
/

&test_arpack_options
 order={2}
 minlev={3}
/
"""

# Common setup function and process handling
def arpack_setup(order, parallel = False):
    nbase = 3
    if parallel:
        nlevels = 4
        nproc = 2
    else:
        nlevels = 3
        nproc = 1
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels,nbase,order,nlevels))
    return run_OFT("./test_arpack", nproc, 60)

def validate_result(expected_lams):
    fid = open('arpack.results','r')
    for (i,expected_lam) in enumerate(expected_lams):
        string_vals = fid.readline().split()
        val1 = float(string_vals[0])
        if abs(val1-expected_lam) > 1.E-6:
            print("FAILED: Incorrect eigenvalue from IRLM solver (order = {0})!".format(i+1))
            print("  Expected = {0}".format(expected_lam))
            print("  Actual =   {0}".format(val1))
            return False
        val2 = float(string_vals[1])
        if abs(val2-expected_lam) > 1.E-6:
            print("FAILED: Incorrect eigenvalue from IRAM solver (order = {0})!".format(i+1))
            print("  Expected = {0}".format(expected_lam))
            print("  Actual =   {0}".format(val2))
            return False
    return True

#============================================================================
expected_lams = (1.7109941816427034, 2.1784222152553760, 2.7901970304756087, 3.6590583531008862)
# Test runner for base test case
@pytest.mark.coverage
def test_base():
    assert arpack_setup(4)
    assert validate_result(expected_lams)

# Test runner for MPI test case
@pytest.mark.coverage
def test_mpi():
    assert arpack_setup(4,True)
    assert validate_result(expected_lams)
