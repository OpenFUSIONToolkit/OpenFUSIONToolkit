#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear sound wave in a triply periodic box (xMHD with all fields on Lagrange elements),
# advanced one wave period and compared with the initial wave; also checks the time steps and fields in
# the plot file
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT
from mug_testing import validate_errors, validate_mug_output

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={10}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=1
/

&cube_options
 mesh_type={11}
 ni=1,1,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_sound_options
 order={2}
 minlev={3}
 delta=1.e-3
 linear={5}
 two_temp={12}
/

&xmhd_options
 xmhd_jcb=F
 xmhd_advec=F
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_visc_heat=F
 vbc='all'
 dt={6}
 eta=1.
 k_db=1.
 nu_par=1.e-3
 nsteps={7}
 rst_freq=25
 lin_tol={8}
 nl_tol=1.E-5
 nu_xmhd={4}
 nclean=500
 rst_ind=0
 maxextrap=2
 kappa_par=1.d-6
 kappa_perp=1.d-6
 ittarget=100
 xmhd_mfnk={9}
/

&xmhd_plot_options
 t0=.9e-5
 dt=1.e-5
 rst_start=0
 rst_end={7}
/
"""

def sound_setup(nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,mf=False,petsc=False,hex_mesh=False,two_temp=False):
    os.chdir(test_dir)
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='1.e-6'
    its='100'
    tol='1.E-12'
    mf_flag='F'
    tt_flag=('T' if two_temp else 'F')
    petsc_flag=('T' if petsc else 'F')
    lin_flag=('T' if linear else 'F')
    # Matrix-free Jacobian products are finite-difference approximations, so the linear tolerance is relaxed
    if (not linear) and mf:
        tol='1.E-9'
        mf_flag='T'
    # The level above nbase is the base mesh distributed over the processes (a transfer level, not a refinement),
    # so MPI runs reproduce the serial errors
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nbase,nlevels,order,minlev,nu_mhd,
                                       lin_flag,dt,its,tol,mf_flag,petsc_flag,mesh_type,tt_flag))
    return run_OFT("./test_sound_lag", nproc, 1000)

# Errors are regression values from a reference run; each may exceed its stored value by up to 1%
def validate_result(nerr_exp,terr_exp,verr_exp,steps_exp=11,linear=False,two_temp=False):
    errors_ok = validate_errors('sound.results', (('Density', nerr_exp), ('Temperature', terr_exp), ('Velocity', verr_exp)), 1.01)
    output_ok = validate_mug_output(steps_exp, linear=linear, two_temp=two_temp)
    return errors_ok and output_ok

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave, quadratic elements on a tetrahedral mesh, serial.'''
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, cubic elements
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, quartic elements
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear sound wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave, quadratic elements on a tetrahedral mesh, serial.'''
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quadratic elements
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave, quadratic elements on a hexahedral mesh, serial.'''
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, cubic elements
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave, cubic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, cubic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quartic elements
@pytest.mark.slow
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4(petsc_flag):
    '''Nonlinear sound wave, quartic elements on a hexahedral mesh with p-multigrid, serial.
    Run with native and PETSc linear algebra.'''
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    '''Nonlinear sound wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov, two-temperature, and (default variant only) PETSc runs.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, quadratic elements
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave, quadratic elements on a hexahedral mesh, serial.'''
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave, cubic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, cubic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(two_temp,petsc_flag):
    '''Linearized sound wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers two-temperature and (default variant only) PETSc runs.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_result(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)
