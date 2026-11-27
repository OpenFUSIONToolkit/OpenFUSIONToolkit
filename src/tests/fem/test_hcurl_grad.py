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
 zerob_test={8}
 adjoint_test={9}
/
"""

# Common setup function and process handling
def hcurl_grad_setup(nbase, nlevels, order, grid_type, mg=False, df='', nu='', petsc=False, zerob=False, adjoint=False):
    petsc_flag=('T' if petsc else 'F')
    mg_flag=('T' if mg else 'F')
    zerob_flag=('T' if zerob else 'F')
    adjoint_flag=('T' if adjoint else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg_flag, df, nu, petsc_flag, zerob_flag, adjoint_flag))
    return run_OFT("./test_hcurl_grad", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('hcurl_grad.results','r')
    its_test = int(fid.readline())
    if iteration_count is not None:
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

def validate_zerob_result(exp_error):
    with open('hcurl_grad.results','r') as fid:
        its_test = int(fid.readline())
        res_test = float(fid.readline())
        error_test = float(fid.readline())
    if (its_test <= 0) or (not res_test < 1.E-8):
        print("FAILED: Solver did not converge!")
        print("  Iterations = {0}".format(its_test))
        print("  Residual =   {0}".format(res_test))
        return False
    if not abs(error_test-exp_error) < 0.05*exp_error:
        print("FAILED: Projection error incorrect!")
        print("  Expected = {0}".format(exp_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

# Mesh levels for single level tests (MPI cases use one more refinement after decomposition)
def single_level_mesh(nlevels,mpi=False):
    if mpi:
        return nlevels-1, nlevels+1
    else:
        return nlevels, nlevels

# 'zerob' mass matrix projection test function
def zerob_projection(nlevels,order,exp_error,grid_type=1,mpi=False,petsc=False):
    nbase, nlev = single_level_mesh(nlevels, mpi)
    assert hcurl_grad_setup(nbase, nlev, order, grid_type, petsc=petsc, zerob=True)
    assert validate_zerob_result(exp_error)

# Single level test function
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False,petsc=False):
    nbase, nlev = single_level_mesh(nlevels, mpi)
    assert hcurl_grad_setup(nbase, nlev, order, grid_type, petsc=petsc)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1,petsc=False):
    if grid_type == 1:
        df_string='0.0, 0.0, 0.39, 0.36, 0.217, 0.172, 0.142'
        nu_string='0, 0, 80, 16, 8, 4, 2'
    else:
        df_string='0.000, 0.000, 0.000, 0.41, 0.25, 0.24, 0.21'
        nu_string='0, 0, 0, 40, 8, 4, 2'
    assert hcurl_grad_setup(2, nlevels+1, order, grid_type, mg=True, df=df_string, nu=nu_string, petsc=petsc)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
def test_r2_p1(petsc_flag=False):
    single_level(2, 1, 57, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 61, 1.0, mpi=mpi, petsc=petsc_flag)
def test_r3_p1_mg(petsc_flag=False):
    multi_level(3, 1, 9, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_r2_p2(petsc_flag=False):
    single_level(2, 2, 167, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 181, 1.0, mpi=mpi, petsc=petsc_flag)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 29, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_r2_p3(petsc_flag=False):
    single_level(2, 3, 361, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 379, 1.0, mpi=mpi, petsc=petsc_flag)
def test_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 64, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_r2_p4(petsc_flag=False):
    single_level(2, 4, 712, 1.0, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 743, 1.0, mpi=mpi, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 156, 1.0, petsc=petsc_flag)

#============================================================================
# Test runners for NP=1
def test_hex_r2_p1(petsc_flag=False):
    single_level(2, 1, 7, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 7, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
# def test_hex_r3_p1_mg(petsc_flag=False):
#     multi_level(3, 1, 1, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=2
def test_hex_r2_p2(petsc_flag=False):
    single_level(2, 2, 192, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 278, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p2_mg(petsc_flag=False):
    multi_level(3, 2, 55, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=3
def test_hex_r2_p3(petsc_flag=False):
    single_level(2, 3, 456, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 580, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 122, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Test runners for NP=4
def test_hex_r2_p4(petsc_flag=False):
    single_level(2, 4, 903, 1.0, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 1035, 1.0, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p4_mg(petsc_flag=False):
    multi_level(3, 4, 248, 1.0, grid_type=2, petsc=petsc_flag)

#============================================================================
# Projection tests using 'zerob' mass matrix
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
@pytest.mark.parametrize("order,exp_error", ((1, 2.7564E-2), (2, 2.0973E-3), (3, 1.4600E-4)))
def test_zerob_proj(order, exp_error, mpi, petsc_flag):
    """L2 projection of an analytic field with zero tangential boundary component
    on a tet cube using the 'zerob' mass matrix, checked against reference errors."""
    zerob_projection(3, order, exp_error, mpi=mpi, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
@pytest.mark.parametrize("order,exp_error", ((1, 3.6377E-2), (2, 3.6105E-3), (3, 1.1195E-4)))
def test_hex_zerob_proj(order, exp_error, mpi, petsc_flag):
    """Same as test_zerob_proj on a hex cube, exercising the hex vertex/boundary DOF layout."""
    zerob_projection(3, order, exp_error, grid_type=2, mpi=mpi, petsc=petsc_flag)

#============================================================================
# Geometric MG prolongation/restriction adjointness test (<y, P x> == <R y, x>)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_mg_adjoint(mpi, petsc_flag):
    """Verify exact transpose relationship between geometric interp and inject for H(curl)+Grad(H^1)."""
    if mpi:
        assert hcurl_grad_setup(1, 3, 1, 1, petsc=petsc_flag, adjoint=True)
    else:
        assert hcurl_grad_setup(3, 3, 1, 1, petsc=petsc_flag, adjoint=True)
    with open('hcurl_grad.results', 'r') as fid:
        rel_err = float(fid.readline())
    assert rel_err < 1.E-12
