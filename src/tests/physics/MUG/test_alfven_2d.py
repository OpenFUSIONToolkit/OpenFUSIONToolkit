#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear shear Alfven wave in a doubly periodic Cartesian slab with the 2D xMHD solver,
# advanced half a wave period and compared with the analytic (inverted) wave
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT
from mug_testing import validate_errors

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={9}
/

&mesh_options
 meshname='test'
 cad_type=92
 nlevels={1}
 nbase={0}
 grid_order=2
 fix_boundary=T
/

&cube_options
 mesh_type={8}
 ni=4,1,0
 rscale=2.0,0.5,0.
 shift=-1.0,-0.25,0.
 ref_per=T,T,F
/

&xmhd_options
  linear={3}
  nl_tol=1.E-8
  cyl=F
  order={2}
  n0=1E19,
  psi0=0.d0,
  velx0=0,
  vely0=0,
  velz0=0,
  t0=1E3,
  by0=0.d0,
  bx0=0.d0,
  bz0=0.d0,
  chi=1.E0
  eta=0.00795774715
  nu=0.59880239
  gamma=1.67
  D_diff=1.E0
  dt={4}
  nsteps={5}
  rst_freq=5
  use_mfnk={6}
  lin_tol={7}
  pm=F
  den_scale=1.E19
/

&xmhd_plot_options
 rst_start=0
 rst_end={5}
/
"""

def alfven_2d_setup(nbase,nlevels,order,
                    linear=False,mf=False,petsc=False,hex_mesh=False):
    mesh_type=1
    if hex_mesh:
        mesh_type=2
    dt='4.0E-7'
    its='250'
    tol='1.E-10'
    mf_flag='F'
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
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, lin_flag,
                                         dt, its, mf_flag, tol, mesh_type, petsc_flag))
    return run_OFT("./test_alfven_2d", nproc, 4000)

# Errors are regression values from a reference run; each may exceed its stored value by up to 5%
def validate_result(verr_exp,berr_exp,steps_exp=11,linear=False):
    return validate_errors('alfven_2d.results', (('Velocity', verr_exp), ('Magnetic', berr_exp)), 1.05)

#============================================================================
# Nonlinear runs, triangular mesh, quadratic elements
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p2(mf,petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a triangular mesh, serial.
    Also covers matrix-free Newton-Krylov.'''
    verr_exp = 3.627799135938326E-2
    berr_exp = 2.7401604608693758E-2
    assert alfven_2d_setup(1,1,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quadratic elements on a triangular mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.627799135938326E-2
    berr_exp = 2.7401604608693758E-2
    assert alfven_2d_setup(1,2,2,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Nonlinear runs, triangular mesh, cubic elements
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p3(mf,petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a triangular mesh, serial.
    Also covers matrix-free Newton-Krylov.'''
    verr_exp = 3.1531577649676595E-2
    berr_exp = 2.7201129077752486E-2
    assert alfven_2d_setup(1,1,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, cubic elements on a triangular mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.1531577649676595E-2
    berr_exp = 2.7201129077752486E-2
    assert alfven_2d_setup(1,2,3,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Nonlinear runs, triangular mesh, quartic elements
@pytest.mark.parametrize("mf", (False, True))
def test_nl_r1_p4(mf,petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a triangular mesh, serial.
    Also covers matrix-free Newton-Krylov.'''
    verr_exp = 3.13872526257328E-2
    berr_exp = 2.719581233696034E-2
    assert alfven_2d_setup(1,1,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quartic elements on a triangular mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    verr_exp = 3.13872526257328E-2
    berr_exp = 2.719581233696034E-2
    assert alfven_2d_setup(1,2,4,mf=mf,petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Linearized runs, triangular mesh, quadratic elements
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a triangular mesh, serial.'''
    verr_exp = 1.811020939389565E-2
    berr_exp = 4.923999966698998E-3
    assert alfven_2d_setup(1,1,2, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a triangular mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    verr_exp = 1.811020939389565E-2
    berr_exp = 4.923999966698998E-3
    assert alfven_2d_setup(1,2,2, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Linearized runs, triangular mesh, cubic elements
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a triangular mesh, serial.'''
    verr_exp = 2.5746106154166067E-3
    berr_exp = 5.148998451109316E-4
    assert alfven_2d_setup(1,1,3, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a triangular mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    verr_exp = 2.5746106154166067E-3
    berr_exp = 5.148998451109316E-4
    assert alfven_2d_setup(1,2,3, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp)

#============================================================================
# Linearized runs, triangular mesh, quartic elements
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a triangular mesh, serial.'''
    verr_exp = 4.997848969169074E-4
    berr_exp = 4.504851399342237E-4
    assert alfven_2d_setup(1,1,4, linear = True, petsc=petsc_flag)
    assert validate_result(verr_exp, berr_exp)

@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a triangular mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    verr_exp = 4.997848969169074E-4
    berr_exp = 4.504851399342237E-4
    assert alfven_2d_setup(1,2,4, petsc=petsc_flag, linear = True)
    assert validate_result(verr_exp, berr_exp)