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
 ni=3,3,0
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
def blagrange_setup(nbase, nlevels, order, grid_type, mg='F', df='', nu='', petsc_flag='F'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, df, nu, petsc_flag))
    return run_OFT("./test_lag_2d", nproc, 60)

def validate_result(iteration_count,converged_error):
    fid = open('lagrange.results','r')
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
def single_level(levels,order,exp_its,exp_error,grid_type=1,petsc_flag='F'):
    assert blagrange_setup(levels[0], levels[1], order, grid_type, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p1_base(levels, petsc_flag):
    single_level(levels, 1, 3, 3.0669984042999903E-002, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p1_1ref(levels, petsc_flag):
    single_level(levels, 1, 12, 0.12247778581844551, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p2_base(levels, petsc_flag):
    single_level(levels, 2, 12, 0.12223325061763389, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p2_1ref(levels, petsc_flag):
    single_level(levels, 2, 40, 0.49022558771158881, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p3_base(levels, petsc_flag):
    single_level(levels, 3, 26, 0.27578221875222847, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p3_1ref(levels, petsc_flag):
    single_level(levels, 3, 76, 1.1032233406550747, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p4_base(levels, petsc_flag):
    single_level(levels, 4, 50, 0.49032174402311279, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_tri_p4_1ref(levels, petsc_flag):
    single_level(levels, 4, 132, 1.9612908119621140, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p1_base(levels, petsc_flag):
    single_level(levels, 1, 1, 1.7777777777777774E-002, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p1_1ref(levels, petsc_flag):
    single_level(levels, 1, 6, 6.4027316352716732E-002, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p2_base(levels, petsc_flag):
    single_level(levels, 2, 6, 6.0977378920293084E-002, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p2_1ref(levels, petsc_flag):
    single_level(levels, 2, 21, 0.24507299559374232, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p3_base(levels, petsc_flag):
    single_level(levels, 3, 10, 0.13790992889879716, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p3_1ref(levels, petsc_flag):
    single_level(levels, 3, 40, 0.55161392909934681, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p4_base(levels, petsc_flag):
    single_level(levels, 4, 25, 0.24510302835546499, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_quad_p4_1ref(levels, petsc_flag):
    single_level(levels, 4, 85, 0.98063152633981288, grid_type=2, petsc_flag=petsc_flag)
