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

&hcurl_op_options
 df_wop={5}
 nu_wop={6}
/

&test_hcurl_options
 order={2}
 mg_test={4}
/
"""

# Common setup function and process handling
def hcurl_setup(nbase, nlevels, order, grid_type, mg='F', df='', nu='', petsc_flag='F'):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg, df, nu, petsc_flag))
    return run_OFT("./test_hcurl", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('hcurl.results','r')
    its_test = int(fid.readline())
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
petsc_flags = ('F','T')
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False, petsc_flag='F'):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert hcurl_setup(nbase, nlev, order, grid_type, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1, petsc_flag='F'):
    if grid_type==1:
        df_string='0.,.65,.65,.374,.323,.295'
        nu_string='0,64,8,4,2,1'
    else:
        df_string='0.,0.,0.74,0.56,0.56,0.56'
        nu_string='0,0,8,2,2,2,2'
    assert hcurl_setup(nlevels, nlevels, order, grid_type, mg='T', df=df_string, nu=nu_string, petsc_flag=petsc_flag)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p1(petsc_flag):
    single_level(2, 1, 29, 1.97622421429438861E-003, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 67, -5.63096078009781378E-004, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p1_mg(petsc_flag):
    multi_level(3, 1, 9, -5.63096078009781378E-004, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p2(petsc_flag):
    single_level(2, 2, 76, 3.66223581984657107E-003, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 124, -8.90637229607081507E-004, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 22, -8.90637229607081507E-004, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p3(petsc_flag):
    single_level(2, 3, 129, 6.23311455636365187E-003, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 159, 4.63583786151764994E-004, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p3_mg(petsc_flag):
    multi_level(3, 3, 43, 4.63583786151764994E-004, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r2_p4(petsc_flag):
    single_level(2, 4, 244, 4.65839958594977104E-003, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 264, -5.54338137661367580E-004, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 93, -5.54338137661367580E-004, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=1
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p1(petsc_flag):
    single_level(2, 1, 1, -2.6101217871994083E-052, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 1, 3.9439369402228879E-036, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
# @pytest.mark.parametrize("petsc_flag", ('F','T'))
# def test_hex_r3_p1_mg(petsc_flag):
#     multi_level(3, 1, 1, 3.9439369402228879E-036, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=2
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p2(petsc_flag):
    single_level(2, 2, 1, 5.2956909820633636E-035, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 72, 2.4107779811264456E-005, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p2_mg(petsc_flag):
    multi_level(3, 2, 30, 2.4107779811264456E-005, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=3
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p3(petsc_flag):
    single_level(2, 3, 72, 1.2988112691509165E-004, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 111, -1.1574898861854135E-004, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p3_mg(petsc_flag):
    multi_level(3, 3, 43, -1.1574898861854135E-004, grid_type=2, petsc_flag=petsc_flag)

#============================================================================
# Test runners for NP=4
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r2_p4(petsc_flag):
    single_level(2, 4, 106, -1.5679728404500014E-004, grid_type=2, petsc_flag=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 185, -1.6139281906695047E-004, grid_type=2, mpi=mpi, petsc_flag=petsc_flag)
@pytest.mark.parametrize("petsc_flag", ('F','T'))
def test_hex_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 65, -1.6139281906695047E-004, grid_type=2, petsc_flag=petsc_flag)
