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

&h0_op_options
 df_lop={5}
 nu_lop={6}
/

&test_h0_options
 order={2}
 mg_test={4}
/
"""

# Common setup function and process handling
def h0_setup(nbase, nlevels, order, grid_type, mg='F', df='', nu='', petsc_flag='F'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg, df, nu, petsc_flag))
    return run_OFT("./test_h0", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('h0.results','r')
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
    assert h0_setup(nbase, nlev, order, grid_type, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1,petsc_flag='F'):
    if grid_type==1:
        df_string='1.,1.,1.,.554,.440,.375,.346'
        nu_string='10,64,16,8,4,2,1'
    else:
        df_string='0.,0.,1.48,0.42,0.40,0.36,0.36'
        nu_string='0,0,8,2,2,2,2'
    assert h0_setup(nlevels, nlevels, order, grid_type, mg='T', df=df_string, nu=nu_string, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p1(petsc_flag):
    single_level(2, 1, 3, 1.1382046922363371E-002, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 17, 8.0768670061291750E-2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1_mg(petsc_flag):
    multi_level(3, 1, 3, 8.0768670061291750E-2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p2(petsc_flag):
    single_level(2, 2, 20, 4.6015646711738915E-003, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 54, 1.5558865184434981E-002, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 11, 1.5558865184434981E-002, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p3(petsc_flag):
    single_level(2, 3, 72, 2.0060529012794457E-002, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 113, 1.6127657918294260E-002, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3_mg(petsc_flag):
    multi_level(3, 3, 22, 1.6127657918294260E-002, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p4(petsc_flag):
    single_level(2, 4, 152, 4.6856137119590809E-002, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 208, 2.1650972314821715E-002, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 40, 2.1650972314821715E-002, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=5
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p5(petsc_flag):
    single_level(2, 5, 316, 9.3830177411447419E-002, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p5(mpi, petsc_flag):
    single_level(3, 5, 373, 4.1329169495863446E-002, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p5_mg(petsc_flag):
     multi_level(3, 5, 87, 4.1329169495863446E-002, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p1(petsc_flag):
    single_level(2, 1, 1, 8.7890624999999896E-003, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 4, 4.8062570071622540E-002, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p1_mg(petsc_flag):
    multi_level(3, 1, 4, 4.8062570071622540E-002, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p2(petsc_flag):
    single_level(2, 2, 5, 1.4816014657915060, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 30, 0.27417620914726598, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 18, 0.27417620914726598, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p3(petsc_flag):
    single_level(2, 3, 24, 1.9945583361965995, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 123, 0.33381655550161526, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3_mg(petsc_flag):
    multi_level(3, 3, 50, 0.33381655550161526, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p4(petsc_flag):
    single_level(2, 4, 49, 1.7015535251800660, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 250, 0.29892496503486543, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 89, 0.29892496503486543, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=5
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p5(petsc_flag):
    single_level(2, 5, 117, 1.6928011950453470, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p5(mpi, petsc_flag):
    single_level(3, 5, 410, 0.29647336194316470, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p5_mg(petsc_flag):
     multi_level(3, 5, 123, 0.29647336194316470, grid_type=2, petsc_flag=petsc_flag)
