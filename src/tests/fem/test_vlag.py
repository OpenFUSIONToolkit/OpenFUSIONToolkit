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
 use_petsc={7}
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
 df_lop={5}
 nu_lop={6}
/

&test_lag_options
 order={2}
 mg_test={4}
/
"""

# Common setup function and process handling
def vec_lagrange_setup(nbase, nlevels, order, grid_type, mg='F', df='', nu='', petsc_flag='F'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg, df, nu, petsc_flag))
    return run_OFT("./test_vlag", nproc, 120)

def validate_result(converged_error):
    retval = True
    fid = open('lagrange.results','r')
    _ = int(fid.readline())
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-8:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        retval = False
    fid.close()
    return retval
# Single level test function
def single_level(nlevels,order,exp_error=1.,grid_type=1,mpi=False,petsc_flag='F'):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert vec_lagrange_setup(nbase, nlev, order, grid_type, petsc_flag=petsc_flag)
    assert validate_result(exp_error)
# Multi-level test function
df_string='1.,1.,1.,.826,.645,.491'
nu_string='10,64,8,4,2,1'
def multi_level(nlevels,order,exp_error,grid_type=1,petsc_flag='F'):
    assert vec_lagrange_setup(nlevels, nlevels, order, grid_type, mg='T', df=df_string, nu=nu_string, petsc_flag=petsc_flag)
    assert validate_result(exp_error)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p1(petsc_flag):
    single_level(2, 1, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1_mg(petsc_flag):
    multi_level(3, 1, 3.*(8.0768670061291750E-2), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p2(petsc_flag):
    single_level(2, 2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 3.*(0.63752436752392194), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p3(petsc_flag):
    single_level(2, 3, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3_mg(petsc_flag):
    multi_level(3, 3, 3.*(2.1575632982055652), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p4(petsc_flag):
    single_level(2, 4, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 3.*(5.1147216970434615), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p1(petsc_flag):
    single_level(2, 1, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p1_mg(petsc_flag):
#     multi_level(3, 1, grid_type=2, 3.*(8.0768670061291750E-2), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p2(petsc_flag):
    single_level(2, 2, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p2_mg(petsc_flag):
#     multi_level(3, 2, grid_type=2, 3.*(0.63752436752392194), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p3(petsc_flag):
    single_level(2, 3, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p3_mg(petsc_flag):
#     multi_level(3, 3, grid_type=2, 3.*(2.1575632982055652), petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p4(petsc_flag):
    single_level(2, 4, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p4_mg(petsc_flag):
#     multi_level(3, 4, grid_type=2, 3.*(5.1147216970434615), petsc_flag=petsc_flag)
