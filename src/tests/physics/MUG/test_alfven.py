#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear shear Alfven wave in a triply periodic box (xMHD, H(curl) magnetic field), advanced
# half a wave period and compared with the analytic (inverted) wave; also checks the time steps and
# fields in the plot file
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
 ni=1,2,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_alfven_options
 order={2}
 minlev={3}
 delta=1.e-3
 linear={5}
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=F
 xmhd_hall=F
 xmhd_adv_temp=F
 xmhd_adv_den=F
 bbc='bc'
 vbc='all'
 dt={6}
 eta=1.e-3
 nu_par=1.e-3
 nsteps={7}
 rst_freq=125
 lin_tol={8}
 nl_tol=1.E-5
 nu_xmhd={4}
 nclean=5000
 rst_ind=0
 maxextrap=2
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

def alfven_setup(nbase,nlevels,order,minlev,nu_mhd='-1',linear=False,
                 mf=False,petsc=False,hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='4.E-7'
    its='250'
    tol='1.E-12'
    mf_flag='F'
    petsc_flag=('T' if petsc else 'F')
    lin_flag=('T' if linear else 'F')
    if not linear:
        if order==4:
            tol='1.E-13'
        if mf:
            mf_flag='T'
    # The level above nbase is the base mesh distributed over the processes (a transfer level, not a refinement),
    # so MPI runs reproduce the serial errors
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, minlev, nu_mhd,
                                       lin_flag, dt, its, tol, mf_flag,
                                       petsc_flag, mesh_type))
    return run_OFT("./test_alfven", nproc, 4000)

# Errors are regression values from a reference run; each may exceed its stored value by up to 1%
def validate_result(berr_exp,verr_exp,steps_exp=11,linear=False):
    errors_ok = validate_errors('alfven.results', (('Magnetic', berr_exp), ('Velocity', verr_exp)), 1.01)
    output_ok = validate_mug_output(steps_exp, linear=linear)
    return errors_ok and output_ok

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.2660712289547556E-002
    verr_exp = 3.2114524936411738E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 4.2660712289547556E-002
    verr_exp = 3.2114524936411738E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, cubic elements
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 6.387170216933607E-003
    verr_exp = 2.3404150169503545E-003
    assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 6.387170216933607E-003
    verr_exp = 2.3404150169503545E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, quartic elements
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 3.8236858815081484E-004
    verr_exp = 2.5485805438810947E-004
    assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 3.8236858815081484E-004
    verr_exp = 2.5485805438810947E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linearized runs, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.2680204579902244E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 4.2680204579902244E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 6.3953902414141325E-003
    verr_exp = 2.3594097996094206E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 6.3953902414141325E-003
    verr_exp = 2.3594097996094206E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 8.595609282107447E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 8.595609282107447E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Nonlinear runs, hexahedral mesh, quadratic elements
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, cubic elements
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, cubic elements on a hexahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quartic elements
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
    #assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
    #assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp)

#============================================================================
# Linearized runs, hexahedral mesh, quadratic elements
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a hexahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    #assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    #assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_result(berr_exp,verr_exp,linear=True)
