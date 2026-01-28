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
 use_petsc={4}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=1
/

&cube_options
 mesh_type={3}
/

&test_hcurl_options
 order={2}
/
"""

# Common setup function and process handling
def hcurl_sop_setup(nbase, nlevels, order, grid_type, petsc):
    petsc_flag=('T' if petsc else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, petsc_flag))
    return run_OFT("./test_hcurl_sop", nproc, 60)

def validate_result():
    fid = open('hcurl.results','r')
    _ = int(fid.readline())
    error_test = float(fid.readline())
    if abs(1.-error_test) > 1.E-6:
        print("FAILED: Field energy incorrect!")
        print("  Expected = {0}".format(1.))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

# Single level test function
def single_level(nlevels, order, mpi=False, grid_type=1, petsc=False):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert hcurl_sop_setup(nbase, nlev, order, grid_type=grid_type, petsc=petsc)
    assert validate_result()

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("grid_type", (1, 2))
def test_r2_p1(grid_type, petsc_flag=False):
    single_level(2, 1, grid_type=grid_type, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1(grid_type, mpi, petsc_flag):
    single_level(3, 1, grid_type=grid_type, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("grid_type", (1, 2))
def test_r2_p2(grid_type, petsc_flag=False):
    single_level(2, 2, grid_type=grid_type, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2(grid_type, mpi, petsc_flag):
    single_level(3, 2, grid_type=grid_type, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("grid_type", (1, 2))
def test_r2_p3(grid_type, petsc_flag=False):
    single_level(2, 3, grid_type=grid_type, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3(grid_type, mpi, petsc_flag):
    single_level(3, 3, grid_type=grid_type, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("grid_type", (1, 2))
def test_r2_p4(grid_type, petsc_flag=False):
    single_level(2, 4, grid_type=grid_type, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4(grid_type, mpi, petsc_flag):
    single_level(3, 4, grid_type=grid_type, mpi=mpi, petsc=petsc_flag)