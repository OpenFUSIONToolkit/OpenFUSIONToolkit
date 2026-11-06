#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: traveling sound wave in a triply periodic box with the mixed H(curl)/Lagrange xMHD module
# (test_sound.F90): density, temperature and velocity after one period against the initial wave
import functools
import pytest
from mug_testing import sound_setup as _sound_setup, validate_sound

sound_setup = functools.partial(_sound_setup, 'hcurl')

#============================================================================
# Non-Linear test runners for NP=2
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the tet mesh with order-2 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    '''The order-2 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 5.9119040036454996E-002
    terr_exp = 5.9231719732675436E-002
    verr_exp = 3.6153775372848231E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the tet mesh with order-3 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    '''The order-3 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 3.4138360614025419E-003
    terr_exp = 3.4144365953483555E-003
    verr_exp = 3.3335039251361225E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the tet mesh with order-4 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    '''The order-4 tet-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 2.7069013754645772E-003
    terr_exp = 2.7074485164640862E-003
    verr_exp = 2.7098032859103993E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the tet mesh with order-2 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(two_temp,petsc_flag):
    '''The order-2 tet-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 5.9421723477144917E-002
    terr_exp = 5.9530267790148306E-002
    verr_exp = 3.6808185835839932E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the tet mesh with order-3 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(two_temp,petsc_flag):
    '''The order-3 tet-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 9.1038624886568657E-003
    terr_exp = 9.1039170745524662E-003
    verr_exp = 9.0864241920798547E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the tet mesh with order-4 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(two_temp,petsc_flag):
    '''The order-4 tet-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 8.8610147661516616E-003
    terr_exp = 8.8613658988151399E-003
    verr_exp = 8.8604718664127627E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=2
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the hex mesh with order-2 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,two_temp,petsc_flag):
    '''The order-2 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 2.1360749524494939E-002
    terr_exp = 2.1362460079360893E-002
    verr_exp = 2.1373898556061369E-002
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=3
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the hex mesh with order-3 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,two_temp,petsc_flag):
    '''The order-3 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 4.5268508665425016E-003
    terr_exp = 4.5266729911366405E-003
    verr_exp = 4.5140710240899129E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Non-Linear test runners for NP=4
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear sound wave advanced one period in a triply periodic box on the hex mesh with order-4 elements:
    density, temperature and velocity errors against the initial wave and the plot output.'''
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,two_temp,petsc_flag):
    '''The order-4 hex-mesh nonlinear case on two MPI ranks with one more refinement level, with and without the
    matrix-free Newton-Krylov solver, with one and two temperatures, with the native and the PETSc linear algebra
    backends.'''
    if (mf or two_temp) and petsc_flag:
        pytest.skip()
    nerr_exp = 3.2905133020702569E-003
    terr_exp = 3.2907964458522860E-003
    verr_exp = 3.2936707864742649E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=2
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the hex mesh with order-2 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(two_temp,petsc_flag):
    '''The order-2 hex-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 2.1008321863506766E-002
    terr_exp = 2.1008318352739700E-002
    verr_exp = 2.0998273206287381E-002
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=3
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the hex mesh with order-3 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(two_temp,petsc_flag):
    '''The order-3 hex-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 9.2415214063574797E-003
    terr_exp = 9.2415774088205841E-003
    verr_exp = 9.2421965458169622E-003
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linear test runners for NP=4
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave (the perturbation evolved about the uniform equilibrium) advanced one period in a triply
    periodic box on the hex mesh with order-4 elements: density, temperature and velocity errors against the initial
    wave and the plot output.'''
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("two_temp", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(two_temp,petsc_flag):
    '''The order-4 hex-mesh linearized case on two MPI ranks with one more refinement level, with one and two
    temperatures, with the native and the PETSc linear algebra backends.'''
    if two_temp and petsc_flag:
        pytest.skip()
    nerr_exp = 8.8588190950046782E-003
    terr_exp = 8.8582074491383549E-003
    verr_exp = 8.8591147973098201E-003
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)
