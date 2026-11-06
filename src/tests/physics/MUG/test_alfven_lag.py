#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: traveling shear Alfven wave in a triply periodic box with the nodal Lagrange xMHD module
# (test_alfven_lag.F90): magnetic field and velocity after half a period against the initial wave
import functools
import pytest
from mug_testing import alfven_setup as _alfven_setup, validate_alfven

alfven_setup = functools.partial(_alfven_setup, 'lagrange')

#============================================================================
# Non-Linear test runners for NP=2
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-2
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    '''The order-2 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-3
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    '''The order-3 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-4
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    '''The order-4 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-2 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''The order-2 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-3 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''The order-3 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-4 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''The order-4 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Hex mesh: at orders 3 and 4 the multigrid hierarchy starts one level higher than on the tet mesh
# Non-Linear test runners for NP=2
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-2
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,petsc_flag):
    '''The order-2 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-3
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    assert alfven_setup(1,1,3,3,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,petsc_flag):
    '''The order-3 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-4
    elements and the nodal Lagrange B representation: magnetic field and velocity errors against the initial wave
    and the plot output.'''
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    assert alfven_setup(1,1,4,3,'0,0,0,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,petsc_flag):
    '''The order-4 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with the native and the PETSc linear algebra backends.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-2 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    '''The order-2 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-3 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    '''The order-3 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-4 elements and the nodal Lagrange B representation: magnetic
    field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    '''The order-4 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
