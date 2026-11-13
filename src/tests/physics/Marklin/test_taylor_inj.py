#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Marklin tests: vacuum field and injector-driven force-free plasma field (lambda = 5) for one flux jump
# plane in a torus from the Fortran API, against regression values of the vacuum energy and field norms
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT, run_command

oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={9}
/

&mesh_options
 meshname='torus'
 cad_type=0
 nlevels={1}
 nbase={0}
 grid_order={8}
 jac_ratio_tol=50.
/

&native_mesh_options
 filename='torus_test.h5'
/

&hcurl_op_options
 df_wop={4}
 nu_wop={5}
/

&lag_op_options
 df_lop={6}
 nu_lop={7}
/

&h1_op_options
 df_lop={6}
 nu_lop={7}
/

&test_taylor_options
 order={2}
 mg_test={3}
/
"""

def taylor_setup(nbase,nlevels,order,mg='F',df_wop='',nu_wop='',df_lop='',nu_lop='',rst=False,petsc=False):
    '''Write `oft.in` and run `test_taylor_inj`, on two MPI ranks when the mesh has distributed levels

    The torus geometry is quadratic for order >= 2. Cached restart files are removed unless `rst` is set.
    '''
    petsc_flag=('T' if petsc else 'F')
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, mg, df_wop, nu_wop,
                                       df_lop, nu_lop, min(order,2), petsc_flag))
    #
    if not rst:
        run_command("rm -f *.rst")
    return run_OFT("./test_taylor_inj", nproc, 180)

def validate_result(energy, mags):
    '''Compare the vacuum energy, the norms of the vacuum field (HVAC), inhomogeneous source field (HCUR) and
    plasma vector potential (GFFA), and the difference between single and grouped injector solves'''
    retval = True
    fid = open('taylor.results','r')
    energy_test = float(fid.readline())
    mag_test = [float(x) for x in fid.readline().split()]
    diff_err = float(fid.readline())
    fid.close()
    #
    if abs((energy - energy_test)/energy) > 1.E-6:
        print("FAILED: Vacuum energy incorrect! ")
        print("  Expected = {0}".format(energy))
        print("  Actual =   {0}".format(energy_test))
        retval = False
    #
    names = ('HVAC', 'HCUR', 'GFFA')
    for (i,mag) in enumerate(mags):
        if abs((mag - mag_test[i])/mag) > 1.E-5:
            print('FAILED: "{0}" magnitude incorrect!'.format(names[i]))
            print("  Expected = {0}".format(mag))
            print("  Actual =   {0}".format(mag_test[i]))
            retval = False
    #
    if abs(diff_err) > 1.E-7:
        print("FAILED: Single calculation does not match group!")
        print("  Error = {0}".format(diff_err))
        retval = False
    return retval

def run_r1_base(energy, mags, order, parallel, petsc_flag):
    '''Base mesh; the parallel run adds a distributed level and repeats the solve from its restart file'''
    if parallel:
        minlev = 1
        nlev = 2
    else:
        minlev = 1
        nlev = 1
    assert taylor_setup(minlev,nlev,order,petsc=petsc_flag)
    assert validate_result(energy, mags)
    if parallel:
        assert taylor_setup(minlev,nlev,order,rst=True,petsc=petsc_flag)
        assert validate_result(energy, mags)

def run_r2_base(energy, mags, order, parallel, petsc_flag):
    '''Refinement level 2; the parallel run adds a distributed level and repeats the solve from its restart file'''
    if parallel:
        minlev = 1
        nlev = 3
    else:
        minlev = 2
        nlev = 2
    assert taylor_setup(minlev,nlev,order,petsc=petsc_flag)
    assert validate_result(energy, mags)
    if parallel:
        assert taylor_setup(minlev,nlev,order,rst=True,petsc=petsc_flag)
        assert validate_result(energy, mags)

#============================================================================
@pytest.mark.parametrize("parallel", (False, True))
def test_r1_p1(parallel, petsc_flag=False):
    '''Linear elements on the base torus mesh; the parallel variant runs on two MPI ranks and repeats the solve
    from the saved restart file.'''
    energy = 7.969787207500343
    mags = (111.78525848354835, 7.3927838909183474E-2, 2.3689570294000877)
    run_r1_base(energy, mags, 1, parallel, petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("parallel", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p1(parallel, petsc_flag):
    '''Linear elements on the torus mesh at refinement level 2; the parallel variant runs on two MPI ranks and
    repeats the solve from the saved restart file (native or PETSc solvers).'''
    energy = 7.584444064306102
    mags = (164.12243228531662, 0.15111708774307858, 7.969982572878115)
    run_r2_base(energy, mags, 1, parallel, petsc_flag)
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p1_mg(petsc_flag):
    '''Linear elements on the torus mesh at refinement level 2, with the multigrid preconditioner on two MPI
    ranks (native or PETSc solvers).'''
    energy = 7.584444064306102
    mags = (164.12243228531662, 0.15111708774307858, 7.969982572878115)
    assert taylor_setup(1,3,1,mg='T',df_wop='0.,0.67,0.63',nu_wop='0,64,1',
                        df_lop='0.,0.96,0.95,0.55',nu_lop='0,64,2,1',petsc=petsc_flag)
    assert validate_result(energy, mags)

#============================================================================
@pytest.mark.parametrize("parallel", (False, True))
def test_r1_p2(parallel, petsc_flag=False):
    '''Quadratic elements on the base torus mesh; the parallel variant runs on two MPI ranks and repeats the
    solve from the saved restart file.'''
    energy = 7.465311073555994
    mags = (99.312891042746003, 0.12514150675409644, 11.854654269477292)
    run_r1_base(energy, mags, 2, parallel, petsc_flag)
@pytest.mark.parametrize("parallel", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p2(parallel, petsc_flag):
    '''Quadratic elements on the torus mesh at refinement level 2; the parallel variant runs on two MPI ranks and
    repeats the solve from the saved restart file (native or PETSc solvers).'''
    energy = 7.4654554644593762
    mags = (159.33067882086129, 0.15585267451452800, 11.630373174047200)
    run_r2_base(energy, mags, 2, parallel, petsc_flag)
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p2_mg(petsc_flag):
    '''Quadratic elements on the torus mesh at refinement level 2, with the multigrid preconditioner on two MPI
    ranks (native or PETSc solvers).'''
    energy = 7.4654554644593762
    mags = (159.33067882086129, 0.15585267451452800, 11.630373174047200)
    assert taylor_setup(1,3,2,mg='T',df_wop='0.,0.67,0.63,0.36',nu_wop='0,64,2,1',
                        df_lop='0.,0.96,0.95,0.55,0.41',nu_lop='0,64,4,2,1',petsc=petsc_flag)
    assert validate_result(energy, mags)

#============================================================================
@pytest.mark.parametrize("parallel", (False, True))
def test_r1_p3(parallel, petsc_flag=False):
    '''Cubic elements on the base torus mesh; the parallel variant runs on two MPI ranks and repeats the solve
    from the saved restart file.'''
    energy = 7.4654401997148261
    mags = (99.963405459273531, 0.13867123084037269, 15.698941734111418)
    run_r1_base(energy, mags, 3, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p3(parallel, petsc_flag):
    '''Cubic elements on the torus mesh at refinement level 2; the parallel variant runs on two MPI ranks and
    repeats the solve from the saved restart file (native or PETSc solvers).'''
    energy = 7.4654749264250251
    mags = (159.56410182996319, 0.15625258194273445, 11.811248986404657)
    run_r2_base(energy, mags, 3, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p3_mg(petsc_flag):
    '''Cubic elements on the torus mesh at refinement level 2, with the multigrid preconditioner on two MPI ranks
    (native or PETSc solvers).'''
    energy = 7.4654749264250251
    mags = (159.56410182996319, 0.15625258194273445, 11.811248986404657)
    assert taylor_setup(1,3,3,mg='T',df_wop='0.,0.67,0.63,0.36,0.3',nu_wop='0,64,4,2,1',
                        df_lop='0.,0.96,0.95,0.55,0.41,0.34',nu_lop='0,64,4,4,2,1',petsc=petsc_flag)
    assert validate_result(energy, mags)

#============================================================================
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
def test_r1_p4(parallel, petsc_flag=False):
    '''Quartic elements on the base torus mesh; the parallel variant runs on two MPI ranks and repeats the solve
    from the saved restart file.'''
    energy = 7.4654670797578868
    mags = (103.02430162981594, 0.14093313194717669, 16.583714239751544)
    run_r1_base(energy, mags, 4, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("parallel", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p4(parallel, petsc_flag):
    '''Quartic elements on the torus mesh at refinement level 2; the parallel variant runs on two MPI ranks and
    repeats the solve from the saved restart file (native or PETSc solvers).'''
    energy = 7.4654796151171006
    mags = (160.85985195059590, 0.15636315584737048, 11.814889822341094)
    run_r2_base(energy, mags, 4, parallel, petsc_flag)
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_r2_p4_mg(petsc_flag):
    '''Quartic elements on the torus mesh at refinement level 2, with the multigrid preconditioner on two MPI
    ranks (native or PETSc solvers).'''
    energy = 7.4654796151171006
    mags = (160.85985195059590, 0.15636315584737048, 11.814889822341094)
    assert taylor_setup(1,3,4,mg='T',df_wop='0.,0.67,0.63,0.36,0.3,0.27',nu_wop='0,64,4,4,2,1',
                        df_lop='0.,0.96,0.95,0.55,0.41,0.34,0.3',nu_lop='0,64,4,4,4,2,1',petsc=petsc_flag)
    assert validate_result(energy, mags)
