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

&hcurl_grad_op_options
 df_mop={5}
 nu_mop={6}
/

&test_hcurl_grad_options
 order={2}
 mg_test={4}
/
"""

# Common setup function and process handling
def hcurl_grad_setup(nbase, nlevels, order, grid_type, mg='F', df='', nu='', petsc_flag='F'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg, df, nu, petsc_flag))
    return run_OFT("./test_hcurl_grad", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('hcurl_grad.results','r')
    its_test = int(fid.readline())
    if iteration_count != None:
        if abs(iteration_count-its_test) >= max(2,.05*iteration_count):
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
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False,petsc_flag='F'):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert hcurl_grad_setup(nbase, nlev, order, grid_type, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
df_string='0.,.39,.36,.217,.172,.142'
nu_string='0,80,16,8,4,2'
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1,petsc_flag='F'):
    assert hcurl_grad_setup(nlevels, nlevels, order, grid_type, mg='T', df=df_string, nu=nu_string, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
iteration_count = None
converged_error = 1.
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p1(petsc_flag):
    single_level(2, 1, iteration_count, converged_error, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, iteration_count, converged_error, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1_mg(petsc_flag):
    multi_level(3, 1, iteration_count, converged_error, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p2(petsc_flag):
    single_level(2, 2, iteration_count, converged_error, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, iteration_count, converged_error, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, iteration_count, converged_error, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p3(petsc_flag):
    single_level(2, 3, iteration_count, converged_error, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, iteration_count, converged_error, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3_mg(petsc_flag):
    multi_level(3, 3, iteration_count, converged_error, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p4(petsc_flag):
    single_level(2, 4, iteration_count, converged_error, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, iteration_count, converged_error, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, iteration_count, converged_error, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=1
iteration_count = None
converged_error = 1.
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p1(petsc_flag):
    single_level(2, 1, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, iteration_count, converged_error, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p1_mg(petsc_flag):
#     multi_level(3, 1, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p2(petsc_flag):
    single_level(2, 2, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, iteration_count, converged_error, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p2_mg(petsc_flag):
#     multi_level(3, 2, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p3(petsc_flag):
    single_level(2, 3, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, iteration_count, converged_error, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_r3_p3_mg(petsc_flag):
#     multi_level(3, 3, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p4(petsc_flag):
    single_level(2, 4, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, iteration_count, converged_error, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p4_mg(petsc_flag):
#     multi_level(3, 4, iteration_count, converged_error, grid_type=2, petsc_flag=petsc_flag)
