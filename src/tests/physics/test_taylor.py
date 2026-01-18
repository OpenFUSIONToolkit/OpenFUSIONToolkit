from __future__ import print_function
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT, run_command

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={11}
/

&mesh_options
 meshname='sphere'
 cad_type=91
 nlevels={1}
 nbase={0}
 grid_order={10}
 jac_ratio_tol=30.
/

&sphere_options
 mesh_type={3}
/

&hcurl_op_options
 df_wop={6}
 nu_wop={7}
/

&lag_op_options
 df_lop={8}
 nu_lop={9}
/

&test_taylor_options
 order={2}
 nm={4}
 mg_test={5}
/
"""

# Common setup function and process handling
def taylor_setup(nbase,nlevels,order,grid_type=1,mg='F',nm=1,
                 df_wop='',nu_wop='',df_lop='',nu_lop='',rst=False,petsc=False):
    petsc_flag=('T' if petsc else 'F')
    max_order = 2
    if grid_type == 2:
        max_order = 2
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type,
                                       nm, mg, df_wop, nu_wop, df_lop,
                                       nu_lop, min(order,max_order), petsc_flag))
    #
    if not rst:
        run_command("rm -f *.rst")
    return run_OFT("./test_taylor", nproc, 180)

def validate_result(lam,tflux):
    retval = True
    fid = open('taylor.results','r')
    lam_test = [float(x) for x in fid.readline().split()]
    tflux_test = [float(x) for x in fid.readline().split()]
    fid.close()
    if isinstance(lam, list):
        if (len(lam_test) != len(lam)):
            print("FAILED: Number of modes incorrect!")
            print("  Expected = {0}".format(len(lam)))
            print("  Actual =   {0}".format(len(lam_test)))
            return False
        #
        for x,y in zip(lam, lam_test):
            if abs(x-y) > 1.E-6:
                print("FAILED: Eigenvalue incorrect!")
                print("  Expected = {0}".format(x))
                print("  Actual =   {0}".format(y))
                retval = False
        if tflux is not None:
            for x,y in zip(tflux, tflux_test):
                y=abs(y)
                if abs(x-y) > 1.E-6:
                    print("FAILED: Toroidal flux incorrect!")
                    print("  Expected = {0}".format(x))
                    print("  Actual =   {0}".format(y))
                    retval = False
    else:
        if (len(lam_test) > 1):
            print("FAILED: Number of modes incorrect!")
            print("  Expected = {0}".format(1))
            print("  Actual =   {0}".format(len(lam_test)))
            return False
        lam_test=lam_test[0]
        tflux_test=abs(tflux_test[0])
        #
        if abs(lam - lam_test) > 1.E-6:
            print("FAILED: Eigenvalue incorrect! ")
            print("  Expected = {0}".format(lam))
            print("  Actual =   {0}".format(lam_test))
            retval = False
        if tflux is not None:
            if abs(tflux - tflux_test) > 1.E-6:
                print("FAILED: Toroidal flux incorrect! ")
                print("  Expected = {0}".format(tflux))
                print("  Actual =   {0}".format(tflux_test))
                retval = False
    return retval

def run_r2_base(lam, tflux, order, parallel, petsc_flag, grid_type=1):
    if parallel:
        minlev = 2
        nlev = 3
    else:
        minlev = 2
        nlev = 2
    assert taylor_setup(minlev,nlev,order,petsc=petsc_flag,grid_type=grid_type)
    assert validate_result(lam,tflux)
    if parallel:
        assert taylor_setup(minlev,nlev,order,rst=True,petsc=petsc_flag,grid_type=grid_type)
        assert validate_result(lam,tflux)

def run_r3_base(lam, tflux, order, parallel, petsc_flag, grid_type=1):
    if parallel:
        minlev = 2
        nlev = 4
    else:
        minlev = 3
        nlev = 3
    assert taylor_setup(minlev,nlev,order,petsc=petsc_flag,grid_type=grid_type)
    assert validate_result(lam,tflux)
    if parallel:
        assert taylor_setup(minlev,nlev,order,rst=True,petsc=petsc_flag,grid_type=grid_type)
        assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p1(parallel, petsc_flag=False):
    lam = 7.8527035135303347
    tflux = None #0.39122960997664691
    run_r2_base(lam, tflux, 1, parallel, petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p1(parallel, petsc_flag=False):
    lam = 5.245541634758173
    tflux = None #0.4119228695699627
    run_r3_base(lam, tflux, 1, parallel, petsc_flag)
def test_r3_p1_mg(petsc_flag=False):
    lam = 5.245541634758173
    tflux = None #0.4119228695699627
    assert taylor_setup(2,4,1,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
                        df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p2(parallel, petsc_flag=False):
    lam = 4.669057291929678
    tflux = None #0.4349569760329484
    run_r2_base(lam, tflux, 2, parallel, petsc_flag)
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p2(parallel, petsc_flag=False):
    lam = 4.511123080535232
    tflux = None #0.41690001577278524
    run_r3_base(lam, tflux, 2, parallel, petsc_flag)
def test_r3_p2_mg(petsc_flag=False):
    lam = 4.511123080535232
    tflux = None #0.41690001577278524
    assert taylor_setup(2,4,2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p3(parallel, petsc_flag=False):
    lam = 4.528662069895297
    tflux = None #0.41592425711440156
    run_r2_base(lam, tflux, 3, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p3(parallel, petsc_flag=False):
    lam = 4.495335122998029
    tflux = None #0.4146402161332857
    run_r3_base(lam, tflux, 3, parallel, petsc_flag)
def test_r3_p3_mg(petsc_flag=False):
    lam = 4.495335122998029
    tflux = None #0.4146402161332857
    assert taylor_setup(2,4,3,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p4(parallel, petsc_flag=False):
    lam = 4.516835126772638
    tflux = None #6.031771190203547E-6
    run_r2_base(lam, tflux, 4, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p4(parallel, petsc_flag=False):
    lam = 4.4949864203192
    tflux = None #0.06331006366265143
    run_r3_base(lam, tflux, 4, parallel, petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg(petsc_flag):
    lam = 4.4949864203192
    tflux = None #0.06331006366265143
    assert taylor_setup(2,4,4,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1_mg_2mode(petsc_flag):
    lam = [5.245541634758124, 5.287437330321279]
    tflux = None #[0.4119228696818406, 2.1285719296774202E-9]
    assert taylor_setup(2,4,1,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
                        df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_mg_2mode(petsc_flag):
    lam = [4.511123080534722, 4.511989205532844]
    tflux = None #[0.41690001611473354, 4.4170623311047773E-10]
    assert taylor_setup(2,4,2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3_mg_2mode(petsc_flag):
    lam = [4.495335122997833, 4.495357996501883]
    tflux = None #[0.4146402219595743, 1.808984001967127E-5]
    assert taylor_setup(2,4,3,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg_2mode(petsc_flag):
    lam = [4.494986420319032, 4.494986701337575]
    tflux = None #[0.414639953968278, 7.096930579996403E-6]
    assert taylor_setup(2,4,4,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p1(parallel, petsc_flag=False):
    lam = 6.2902752655699627
    tflux = None
    run_r2_base(lam, tflux, 1, parallel, petsc_flag, grid_type=2)
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p1(parallel, petsc_flag=False):
    lam = 4.952188637006833
    tflux = None
    run_r3_base(lam, tflux, 1, parallel, petsc_flag, grid_type=2)
# def test_hex_r3_p1_mg(petsc_flag=False):
#     lam = 4.952188637006833
#     tflux = None
#     assert taylor_setup(2,4,1,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
#                         df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',petsc=petsc_flag)
#     assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p2(parallel, petsc_flag=False):
    lam = 4.628157911049961
    tflux = None
    run_r2_base(lam, tflux, 2, parallel, petsc_flag, grid_type=2)
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p2(parallel, petsc_flag=False):
    lam = 4.503430041172023
    tflux = None
    run_r3_base(lam, tflux, 2, parallel, petsc_flag, grid_type=2)
def test_hex_r3_p2_mg(petsc_flag=False):
    lam = 4.503430041172023
    tflux = None
    assert taylor_setup(2,4,2,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p3(parallel, petsc_flag=False):
    lam = 4.504562497987162
    tflux = None
    run_r2_base(lam, tflux, 3, parallel, petsc_flag, grid_type=2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p3(parallel, petsc_flag=False):
    lam = 4.493760065748465
    tflux = None
    run_r3_base(lam, tflux, 3, parallel, petsc_flag, grid_type=2)
def test_hex_r3_p3_mg(petsc_flag=False):
    lam = 4.493760065748465
    tflux = None
    assert taylor_setup(2,4,3,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# Test runner for base test case
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p4(parallel, petsc_flag=False):
    lam = 4.496860916053858
    tflux = None
    run_r2_base(lam, tflux, 4, parallel, petsc_flag, grid_type=2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p4(parallel, petsc_flag=False):
    lam = 4.493588233334749
    tflux = None
    run_r3_base(lam, tflux, 4, parallel, petsc_flag, grid_type=2)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_mg(petsc_flag):
    lam = 4.493588233334749
    tflux = None
    assert taylor_setup(2,4,4,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',petsc=petsc_flag)
    assert validate_result(lam,tflux)

#============================================================================
# # Test runner for base test case
# @pytest.mark.parametrize("petsc_flag", (True, False))
# def test_hex_r3_p1_mg_2mode(petsc_flag):
#     lam = [5.245541634758124, 5.287437330321279]
#     tflux = [0.4119228696818406, 2.1285719296774202E-9]
#     assert taylor_setup(2,4,1,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
#                         df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',nm=2,petsc=petsc_flag)
#     assert validate_result(lam,tflux)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2_mg_2mode(petsc_flag):
    lam = [4.5034300411720345, 4.5034300411720345]
    tflux = None
    assert taylor_setup(2,4,2,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3_mg_2mode(petsc_flag):
    lam = [4.493760065748454, 4.493760065748786]
    tflux = None
    assert taylor_setup(2,4,3,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_mg_2mode(petsc_flag):
    lam = [4.493588233334745, 4.4935882333347275]
    tflux = None
    assert taylor_setup(2,4,4,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam,tflux)
