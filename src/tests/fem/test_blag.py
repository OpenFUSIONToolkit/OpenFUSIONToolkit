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
 use_petsc={6}
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

&lag_op_options
 df_lop={4}
 nu_lop={5}
/

&test_blag_options
 order={2}
/
"""

# Common setup function and process handling
def blagrange_setup(nbase, nlevels, order, grid_type, mg=False, df='', nu='', petsc=False):
    petsc_flag=('T' if petsc else 'F')
    # mg_flag=('T' if mg else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, df, nu, petsc_flag))
    return run_OFT("./test_blag", nproc, 60)

def validate_result(iteration_count,converged_error):
    fid = open('lagrange.results','r')
    its_test = int(fid.readline())
    if iteration_count != None:
        if abs(iteration_count-its_test) >= max(3,0.05*iteration_count):
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            return False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

# Single level test function
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False,petsc=False):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert blagrange_setup(nbase, nlev, order, grid_type, petsc=petsc)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
def test_r2_p1(petsc_flag=False):
    single_level(2, 1, 1, 2.3437500000000000E-002, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 3, 0.14575195312500000, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_r2_p2(petsc_flag=False):
    single_level(2, 2, 3, 0.15562500000000043, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 10, 0.65168062291538276, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_r2_p3(petsc_flag=False):
    single_level(2, 3, 9, 0.36475123309992658, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 38, 1.4704555201530953, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_r2_p4(petsc_flag=False):
    single_level(2, 4, 18, 0.65348841031767191, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 74, 2.6149174992726563, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=1
def test_hex_r2_p1(petsc_flag=False):
    single_level(2, 1, 1, 5.2734375000000000E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 3, 0.17916772959183688, grid_type=2, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_hex_r2_p2(petsc_flag=False):
    single_level(2, 2, 3, 0.15991066475591753, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 10, 0.65263272101701786, grid_type=2, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_hex_r2_p3(petsc_flag=False):
    single_level(2, 3, 6, 0.36786823385360656, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 21, 1.4709839088462502, grid_type=2, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_hex_r2_p4(petsc_flag=False):
    single_level(2, 4, 12, 0.65297513559274001, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 45, 2.6148592894354170, grid_type=2, mpi=mpi, petsc=petsc_flag)
