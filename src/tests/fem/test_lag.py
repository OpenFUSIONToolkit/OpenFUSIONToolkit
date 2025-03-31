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
{12}
 ref_per={8},{9},{10}
 shift=-0.5,-0.5,-0.5
/

&lag_op_options
 df_lop={5}
 nu_lop={6}
/

&test_lag_options
 order={2}
 mg_test={4}
 grnd_dir={11}
/
"""

# Common setup function and process handling
def lagrange_setup(nbase, nlevels, order, grid_type, mg=False, df='', nu='', petsc=False, xperiodic=False, yperiodic=False, zperiodic=False):
    petsc_flag=('T' if petsc else 'F')
    mg_flag=('T' if mg else 'F')
    xper_flag = ('T' if xperiodic else 'F')
    yper_flag = ('T' if yperiodic else 'F')
    zper_flag = ('T' if zperiodic else 'F')
    if xperiodic:
        if yperiodic:
            grnd_dir = 3
        else:
            grnd_dir = 2
    else:
        grnd_dir = 1
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    if (grid_type == 2) and (xperiodic or yperiodic or zperiodic):
        ni_line = ' ni=2,2,2'
    else:
        ni_line = '! ni=2,2,2'
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type, mg_flag, df, nu, petsc_flag, xper_flag, yper_flag, zper_flag, grnd_dir, ni_line))
    return run_OFT("./test_lag", nproc, 120)

def validate_result(iteration_count,converged_error):
    fid = open('lagrange.results','r')
    its_test = int(fid.readline())
    ret_val = True
    if iteration_count != None:
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
def single_level(nlevels,order,exp_its,exp_error,grid_type=1,mpi=False,petsc=False,xperiodic=False,yperiodic=False,zperiodic=False):
    if mpi:
        nbase = nlevels-1
        nlev = nlevels+1
    else:
        nbase = nlevels
        nlev = nlevels
    assert lagrange_setup(nbase, nlev, order, grid_type, petsc=petsc, xperiodic=xperiodic, yperiodic=yperiodic, zperiodic=zperiodic)
    assert validate_result(exp_its, exp_error)
# Multi-level test function
def multi_level(nlevels,order,exp_its,exp_error,grid_type=1,petsc=False):
    if grid_type==1:
        df_string='0.0, 0.0, 1.0, 1.0, 0.82, 0.64, 0.49'
        nu_string='0, 0, 64, 8, 4, 2, 1'
    else:
        df_string='0.0, 0.0, 0.0, 1.48, 1.40, 0.78, 0.40'
        nu_string='0, 0, 0, 32, 4, 2, 1'
    assert lagrange_setup(2, nlevels+1, order, grid_type, mg=True, df=df_string, nu=nu_string, petsc=petsc)
    assert validate_result(exp_its, exp_error)

#============================================================================
# Test runners for NP=1
def test_r2_p1(petsc_flag=False):
    single_level(2, 1, 3, 1.1382046922363371E-002, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 17, 8.0768670061291750E-2, mpi=mpi, petsc=petsc_flag)
def test_r3_p1_mg(petsc_flag=False):
    multi_level(3, 1, 4, 8.0768670061291750E-2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1_periodic1(mpi, petsc_flag):
    single_level(3, 1, 27, 0.370222100762129, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 1, 23, 0.370125278783165, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 1, 27, 0.370222100762129, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1_periodic2(mpi, petsc_flag):
    single_level(3, 1, 30, 3.5062876049795166, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 1, 27, 3.529838658073092, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 1, 30, 3.5062876049795166, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=2
def test_r2_p2(petsc_flag=False):
    single_level(2, 2, 19, 7.8106713098390637E-002, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 47, 0.63752436752392194, mpi=mpi, petsc=petsc_flag)
def test_r3_p2_mg(petsc_flag=False):
    multi_level(3, 2, 9, 0.63752436752392194, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_periodic1(mpi, petsc_flag):
    single_level(3, 2, 77, 3.478845838102087, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 2, 73, 3.4786810541763407, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 2, 77, 3.478845838102087, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_periodic2(mpi, petsc_flag):
    single_level(3, 2, 80, 34.13281249999963, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 2, 76, 34.13281249999963, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 2, 80, 34.13281249999963, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=3
def test_r2_p3(petsc_flag=False):
    single_level(2, 3, 43, 0.26895148657016710, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 87, 2.1575632982055652, mpi=mpi, petsc=petsc_flag)
def test_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 12, 2.1575632982055652, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3_periodic1(mpi, petsc_flag):
    single_level(3, 3, 144, 11.767079839311773, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 3, 136, 11.767079839311773, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 3, 143, 11.767079839311773, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3_periodic2(mpi, petsc_flag):
    single_level(3, 3, 151, 115.19965277776963, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 3, 144, 115.19965277776963, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 3, 151, 115.19965277776963, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=4
def test_r2_p4(petsc_flag=False):
    single_level(2, 4, 79, 0.63940369781706052, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 154, 5.1147216970434615, mpi=mpi, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 21, 5.1147216970434615, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_periodic1(mpi, petsc_flag):
    single_level(3, 4, 260, 27.893884170145597, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 4, 246, 27.893884170145597, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 4, 260, 27.893884170145597, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_periodic2(mpi, petsc_flag):
    single_level(3, 4, 272, 273.0664062499967, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 4, 259, 273.0664062499967, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 4, 272, 273.0664062499967, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=1 Hex grid
def test_hex_r2_p1(petsc_flag=False):
    single_level(2, 1, 1, 8.7890624999999896E-003, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1(mpi, petsc_flag):
    single_level(3, 1, 4, 4.8062570071622540E-002, grid_type=2, mpi=mpi, petsc=petsc_flag)
# def test_hex_r3_p1_mg(petsc_flag=False):
#     multi_level(3, 1, 4, 4.8062570071622540E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1_periodic1(mpi, petsc_flag):
    single_level(3, 1, 10, 0.8938986709682393, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 1, 10, 0.8938986709682393, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 1, 10, 0.8938986709682393, grid_type=2, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p1_periodic2(mpi, petsc_flag):
    single_level(3, 1, 4, 4.2656249999999245, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 1, 4, 4.2656249999999245, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 1, 4, 4.2656249999999245, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=2
def test_hex_r2_p2(petsc_flag=False):
    single_level(2, 2, 4, 3.7489146284692222E-002, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2(mpi, petsc_flag):
    single_level(3, 2, 18, 0.31784865062629630, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p2_mg(petsc_flag=False):
    multi_level(3, 2, 9, 0.31784865062629630, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2_periodic1(mpi, petsc_flag):
    single_level(3, 2, 31, 6.97266891509049, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 2, 31, 6.97266891509049, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 2, 31, 6.97266891509049, grid_type=2, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2_periodic2(mpi, petsc_flag):
    single_level(3, 2, 22, 34.13281249999523, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 2, 22, 34.13281249999523, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 2, 22, 34.13281249999523, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=3
def test_hex_r2_p3(petsc_flag=False):
    single_level(2, 3, 10, 0.13534115933200871, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3(mpi, petsc_flag):
    single_level(3, 3, 39, 1.0791049645224631, grid_type=2, mpi=mpi, petsc=petsc_flag)
def test_hex_r3_p3_mg(petsc_flag=False):
    multi_level(3, 3, 11, 1.0791049645224631, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3_periodic1(mpi, petsc_flag):
    single_level(3, 3, 73, 23.535508841766205, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 3, 73, 23.535508841766205, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 3, 73, 23.535508841766205, grid_type=2, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3_periodic2(mpi, petsc_flag):
    single_level(3, 3, 38, 115.19965277776782, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 3, 38, 115.19965277776782, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 3, 38, 115.19965277776782, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)

#============================================================================
# Test runners for NP=4
def test_hex_r2_p4(petsc_flag=False):
    single_level(2, 4, 26, 0.31862239339433901, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4(mpi, petsc_flag):
    single_level(3, 4, 97, 2.5568335937210307, grid_type=2, mpi=mpi, petsc=petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_mg(petsc_flag):
    multi_level(3, 4, 20, 2.5568335937210307, grid_type=2, petsc=petsc_flag)
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_periodic1(mpi, petsc_flag):
    single_level(3, 4, 188, 55.78760364877135, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True)
    single_level(3, 4, 188, 55.78760364877135, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True)
    single_level(3, 4, 188, 55.78760364877135, grid_type=2, mpi=mpi, petsc=petsc_flag, zperiodic=True)
@pytest.mark.coverage
@pytest.mark.parametrize("mpi", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_periodic2(mpi, petsc_flag):
    single_level(3, 4, 101, 273.06640625001774, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, yperiodic=True)
    single_level(3, 4, 101, 273.06640625001774, grid_type=2, mpi=mpi, petsc=petsc_flag, xperiodic=True, zperiodic=True)
    single_level(3, 4, 101, 273.06640625001774, grid_type=2, mpi=mpi, petsc=petsc_flag, yperiodic=True, zperiodic=True)
