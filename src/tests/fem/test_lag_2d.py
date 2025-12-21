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
 ref_per={7},{8},F
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
def blagrange_setup(nbase, nlevels, order, grid_type, mg=False, df='', nu='', petsc=False, xperiodic=False, yperiodic=False):
    petsc_flag=('T' if petsc else 'F')
    # mg_flag=('T' if mg else 'F')
    xper_flag = ('T' if xperiodic else 'F')
    yper_flag = ('T' if yperiodic else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, df, nu, petsc_flag, xper_flag, yper_flag))
    return run_OFT("./test_lag_2d", nproc, 60)

def validate_result(iteration_count,converged_error):
    fid = open('lagrange.results','r')
    its_test = int(fid.readline())
    ret_val = True
    if iteration_count is not None:
        if abs(iteration_count-its_test) >= max(3,0.05*iteration_count):
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            ret_val = False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        ret_val = False
    return ret_val

# Single level test function
def single_level(levels,order,exp_its,exp_error,grid_type=1,petsc=False,xperiodic=False,yperiodic=False):
    assert blagrange_setup(levels[0], levels[1], order, grid_type, petsc=petsc, xperiodic=xperiodic, yperiodic=yperiodic)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p1_base(levels, petsc_flag):
    single_level(levels, 1, 3, 3.0669984042999903E-002, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p1_1ref(levels, petsc_flag):
    single_level(levels, 1, 12, 0.12247778581844551, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p1_periodic(levels, petsc_flag):
    single_level(levels, 1, 7, 0.5928317097479414, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 1, 7, 0.5928317097479414, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p2_base(levels, petsc_flag):
    single_level(levels, 2, 12, 0.12223325061763389, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p2_1ref(levels, petsc_flag):
    single_level(levels, 2, 40, 0.49022558771158881, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p2_periodic(levels, petsc_flag):
    single_level(levels, 2, 19, 2.399992766203704, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 2, 19, 2.399992766203704, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p3_base(levels, petsc_flag):
    single_level(levels, 3, 26, 0.27578221875222847, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p3_1ref(levels, petsc_flag):
    single_level(levels, 3, 76, 1.1032233406550747, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p3_periodic(levels, petsc_flag):
    single_level(levels, 3, 63, 5.399996784979429, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 3, 63, 5.399996784979429, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p4_base(levels, petsc_flag):
    single_level(levels, 4, 50, 0.49032174402311279, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p4_1ref(levels, petsc_flag):
    single_level(levels, 4, 132, 1.9612908119621140, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_tri_p4_periodic(levels, petsc_flag):
    single_level(levels, 4, 113, 9.599998191551755, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 4, 113, 9.599998191551755, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p1_base(levels, petsc_flag):
    single_level(levels, 1, 1, 1.7777777777777774E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p1_1ref(levels, petsc_flag):
    single_level(levels, 1, 6, 6.4027316352716732E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p1_periodic(levels, petsc_flag):
    single_level(levels, 1, 3, 0.2997685185185187, grid_type=2, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 1, 3, 0.2997685185185187, grid_type=2, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p2_base(levels, petsc_flag):
    single_level(levels, 2, 6, 6.0977378920293084E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p2_1ref(levels, petsc_flag):
    single_level(levels, 2, 21, 0.24507299559374232, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p2_periodic(levels, petsc_flag):
    single_level(levels, 2, 12, 1.1999421296296282, grid_type=2, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 2, 12, 1.1999421296296282, grid_type=2, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p3_base(levels, petsc_flag):
    single_level(levels, 3, 10, 0.13790992889879716, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p3_1ref(levels, petsc_flag):
    single_level(levels, 3, 40, 0.55161392909934681, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p3_periodic(levels, petsc_flag):
    single_level(levels, 3, 20, 2.69997427983548, grid_type=2, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 3, 20, 2.69997427983548, grid_type=2, petsc=petsc_flag, yperiodic=True)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("levels", ((1,1), (1,2)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p4_base(levels, petsc_flag):
    single_level(levels, 4, 25, 0.24510302835546499, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p4_1ref(levels, petsc_flag):
    single_level(levels, 4, 85, 0.98063152633981288, grid_type=2, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("levels", ((2,2), (2,3), (1,3)))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_quad_p4_periodic(levels, petsc_flag):
    single_level(levels, 4, 44, 4.799985532407216, grid_type=2, petsc=petsc_flag, xperiodic=True)
    single_level(levels, 4, 44, 4.799985532407216, grid_type=2, petsc=petsc_flag, yperiodic=True)
