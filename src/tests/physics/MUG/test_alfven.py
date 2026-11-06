#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: traveling shear Alfven wave in a triply periodic box with the mixed H(curl)/Lagrange xMHD
# module (test_alfven.F90): magnetic field and velocity after half a period against the initial wave
import functools
import pytest
from mug_testing import alfven_setup as _alfven_setup, validate_alfven

alfven_setup = functools.partial(_alfven_setup, 'hcurl')

#============================================================================
# Non-Linear test runners for NP=2
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-2
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 4.2660712289547556E-002
    verr_exp = 3.2114524936411738E-002
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
    berr_exp = 4.2660712289547556E-002
    verr_exp = 3.2114524936411738E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-3
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 6.387170216933607E-003
    verr_exp = 2.3404150169503545E-003
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
    berr_exp = 6.387170216933607E-003
    verr_exp = 2.3404150169503545E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the tet mesh with order-4
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 3.8236858815081484E-004
    verr_exp = 2.5485805438810947E-004
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
    berr_exp = 3.8236858815081484E-004
    verr_exp = 2.5485805438810947E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-2 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 4.2680204579902244E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''The order-2 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 4.2680204579902244E-002
    verr_exp = 3.2104584635227973E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-3 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 6.3953902414141325E-003
    verr_exp = 2.3594097996094206E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''The order-3 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 6.3953902414141325E-003
    verr_exp = 2.3594097996094206E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the tet mesh with order-4 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 8.595609282107447E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''The order-4 tet-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 8.595609282107447E-004
    verr_exp = 3.8358125487992524E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Hex mesh: at orders 3 and 4 the multigrid hierarchy starts one level higher than on the tet mesh
# Non-Linear test runners for NP=2
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-2
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
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
    berr_exp = 1.8217184894941060E-002
    verr_exp = 1.8207419422505190E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-3
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
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
    berr_exp = 2.2785287899224150E-003
    verr_exp = 2.2636641474315156E-003
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear shear Alfven wave advanced half a period in a triply periodic box on the hex mesh with order-4
    elements: magnetic field and velocity errors against the initial wave and the plot output.'''
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
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
    berr_exp = 1.8358574384029231E-004
    verr_exp = 1.8591003585664032E-004
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-2 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    '''The order-2 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 1.8214844434637032E-002
    verr_exp = 1.8184066089285304E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-3 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    '''The order-3 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 2.2981509470689828E-003
    verr_exp = 2.2815413576994501E-003
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized shear Alfven wave (the perturbation evolved about the uniform equilibrium) advanced half a period in
    a triply periodic box on the hex mesh with order-4 elements: magnetic field and velocity errors against the
    initial wave and the plot output.'''
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    '''The order-4 hex-mesh linearized case on two MPI ranks with one more refinement level, with the native and the
    PETSc linear algebra backends.'''
    berr_exp = 1.5943398846056970E-003
    verr_exp = 1.3887929641779406E-003
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
