#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: Taylor states (force-free eigenmodes, J. B. Taylor, Phys. Rev. Lett. 33, 1139 (1974)) in a sphere
# on tet and hex meshes (test_taylor.F90). The lowest eigenvalue converges with order and refinement to the analytic
# value lambda*a = 4.4934 (Chandrasekhar and Kendall, Astrophys. J. 126, 457 (1957)); the toroidal flux of each mode
# is written but not checked.
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
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
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        # Mesh geometry order is capped at quadratic
        fid.write(oft_in_template.format(nbase, nlevels, order, grid_type,
                                       nm, mg, df_wop, nu_wop, df_lop,
                                       nu_lop, min(order,2), petsc_flag))
    #
    if not rst:
        run_command("rm -f *.rst")
    return run_OFT("./test_taylor", nproc, 180)

def validate_result(lam,tflux=None):
    '''Check the eigenvalue(s) in "taylor.results", and the toroidal flux of each mode if `tflux` is given'''
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

def run_base(lam, order, parallel, petsc_flag, nlevels, grid_type=1):
    '''Taylor state on `nlevels` refinements of the base sphere mesh; in parallel the finest level is distributed over
    two ranks and the run is repeated from the restart file'''
    nbase = 2 if parallel else nlevels
    nlev = nlevels + 1 if parallel else nlevels
    assert taylor_setup(nbase,nlev,order,petsc=petsc_flag,grid_type=grid_type)
    assert validate_result(lam)
    if parallel:
        assert taylor_setup(nbase,nlev,order,rst=True,petsc=petsc_flag,grid_type=grid_type)
        assert validate_result(lam)

@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p1(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the tet mesh with order-1 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 7.8527035135303347
    run_base(lam, 1, parallel, petsc_flag, 2)
@pytest.mark.coverage
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p1(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the tet mesh with order-1 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 5.245541634758173
    run_base(lam, 1, parallel, petsc_flag, 3)
def test_r3_p1_mg(petsc_flag=False):
    '''The three-refinement tet mesh case at order 1 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 5.245541634758173
    assert taylor_setup(2,4,1,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
                        df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p2(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the tet mesh with order-2 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.669057291929678
    run_base(lam, 2, parallel, petsc_flag, 2)
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p2(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the tet mesh with order-2 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.511123080535232
    run_base(lam, 2, parallel, petsc_flag, 3)
def test_r3_p2_mg(petsc_flag=False):
    '''The three-refinement tet mesh case at order 2 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.511123080535232
    assert taylor_setup(2,4,2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p3(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the tet mesh with order-3 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.528662069895297
    run_base(lam, 3, parallel, petsc_flag, 2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p3(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the tet mesh with order-3 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.495335122998029
    run_base(lam, 3, parallel, petsc_flag, 3)
def test_r3_p3_mg(petsc_flag=False):
    '''The three-refinement tet mesh case at order 3 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.495335122998029
    assert taylor_setup(2,4,3,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r2_p4(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the tet mesh with order-4 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.516835126772638
    run_base(lam, 4, parallel, petsc_flag, 2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r3_p4(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the tet mesh with order-4 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.4949864203192
    run_base(lam, 4, parallel, petsc_flag, 3)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg(petsc_flag):
    '''The three-refinement tet mesh case at order 4 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.4949864203192
    assert taylor_setup(2,4,4,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p1_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 1 with the multigrid solver.
    The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 1 with the multigrid solver.'''
    lam = [5.245541634758124, 5.287437330321279]
    assert taylor_setup(2,4,1,mg='T',df_wop='0.0,99.,.68,.64',nu_wop='0,10,64,1',
                        df_lop='0.0,1.,1.,.9',nu_lop='0,10,64,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p2_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 2 with the multigrid solver.
    The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 2 with the multigrid solver.'''
    lam = [4.511123080534722, 4.511989205532844]
    assert taylor_setup(2,4,2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p3_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 3 with the multigrid solver.
    The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 3 with the multigrid solver.'''
    lam = [4.495335122997833, 4.495357996501883]
    assert taylor_setup(2,4,3,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r3_p4_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 4 with the multigrid solver.
    The two lowest Taylor-state eigenvalues on the three-refinement tet mesh at order 4 with the multigrid solver.'''
    lam = [4.494986420319032, 4.494986701337575]
    assert taylor_setup(2,4,4,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)

# Hex mesh cases (P1 multigrid is not exercised on the hex mesh)
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p1(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the hex mesh with order-1 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 6.2902752655699627
    run_base(lam, 1, parallel, petsc_flag, 2, grid_type=2)
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p1(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the hex mesh with order-1 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.952188637006833
    run_base(lam, 1, parallel, petsc_flag, 3, grid_type=2)

@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p2(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the hex mesh with order-2 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.628157911049961
    run_base(lam, 2, parallel, petsc_flag, 2, grid_type=2)
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p2(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the hex mesh with order-2 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.503430041172023
    run_base(lam, 2, parallel, petsc_flag, 3, grid_type=2)
def test_hex_r3_p2_mg(petsc_flag=False):
    '''The three-refinement hex mesh case at order 2 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.503430041172023
    assert taylor_setup(2,4,2,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p3(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the hex mesh with order-3 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.504562497987162
    run_base(lam, 3, parallel, petsc_flag, 2, grid_type=2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p3(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the hex mesh with order-3 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.493760065748465
    run_base(lam, 3, parallel, petsc_flag, 3, grid_type=2)
def test_hex_r3_p3_mg(petsc_flag=False):
    '''The three-refinement hex mesh case at order 3 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.493760065748465
    assert taylor_setup(2,4,3,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r2_p4(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on two refinements of the hex mesh with order-4 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.496860916053858
    run_base(lam, 4, parallel, petsc_flag, 2, grid_type=2)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_hex_r3_p4(parallel, petsc_flag=False):
    '''Lowest Taylor-state eigenvalue in a sphere on three refinements of the hex mesh with order-4 H(curl) elements,
    serially and on two MPI ranks with a restart from file.'''
    lam = 4.493588233334749
    run_base(lam, 4, parallel, petsc_flag, 3, grid_type=2)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_mg(petsc_flag):
    '''The three-refinement hex mesh case at order 4 with the multigrid-preconditioned eigensolver, with per-level
    smoother damping and iteration counts.'''
    lam = 4.493588233334749
    assert taylor_setup(2,4,4,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',petsc=petsc_flag)
    assert validate_result(lam)

@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p2_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement hex mesh at order 2 with the multigrid solver (a
    nearly degenerate pair on the hex mesh).'''
    lam = [4.5034300411720345, 4.5034300411720345]
    assert taylor_setup(2,4,2,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35',nu_wop='0,10,64,2,1',
                        df_lop='0.0,1.,1.,.9,.78',nu_lop='0,10,64,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p3_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement hex mesh at order 3 with the multigrid solver (a
    nearly degenerate pair on the hex mesh).'''
    lam = [4.493760065748454, 4.493760065748786]
    assert taylor_setup(2,4,3,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3',nu_wop='0,10,64,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578',nu_lop='0,10,64,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_r3_p4_mg_2mode(petsc_flag):
    '''The two lowest Taylor-state eigenvalues on the three-refinement hex mesh at order 4 with the multigrid solver (a
    nearly degenerate pair on the hex mesh).'''
    lam = [4.493588233334745, 4.4935882333347275]
    assert taylor_setup(2,4,4,grid_type=2,mg='T',df_wop='0.0,99.,.68,.64,.35,.3,.25',nu_wop='0,10,64,8,4,2,1',
                        df_lop='0.0,1.,1.,.9,.78,.578,.431',nu_lop='0,10,64,8,4,2,1',nm=2,petsc=petsc_flag)
    assert validate_result(lam)
