#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear sound wave in a triply periodic box (xMHD, H(curl) magnetic field), advanced one
# wave period and compared with the initial wave; also checks the time steps and fields in the plot file
import functools
import pytest
from mug_testing import sound_setup as _sound_setup, validate_sound

sound_setup = functools.partial(_sound_setup, 'hcurl')

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave, quadratic elements on a tetrahedral mesh, serial.'''
    nerr_exp = 2.003599E-02
    terr_exp = 2.002004E-02
    verr_exp = 1.916745E-02
    assert sound_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 2.003599E-02
    terr_exp = 2.002004E-02
    verr_exp = 1.916745E-02
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, cubic elements
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 3.136258E-03
    terr_exp = 3.139024E-03
    verr_exp = 3.416849E-03
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 3.136258E-03
    terr_exp = 3.139024E-03
    verr_exp = 3.416849E-03
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, quartic elements
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear sound wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 2.679343E-03
    terr_exp = 2.679116E-03
    verr_exp = 2.683346E-03
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 2.679343E-03
    terr_exp = 2.679116E-03
    verr_exp = 2.683346E-03
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave, quadratic elements on a tetrahedral mesh, serial.'''
    nerr_exp = 2.164925E-02
    terr_exp = 2.164065E-02
    verr_exp = 2.058023E-02
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 2.164925E-02
    terr_exp = 2.164065E-02
    verr_exp = 2.058023E-02
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 9.005905E-03
    terr_exp = 9.006802E-03
    verr_exp = 9.107430E-03
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 9.005905E-03
    terr_exp = 9.006802E-03
    verr_exp = 9.107430E-03
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, tetrahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    nerr_exp = 8.884944E-03
    terr_exp = 8.885293E-03
    verr_exp = 8.885729E-03
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 8.884944E-03
    terr_exp = 8.885293E-03
    verr_exp = 8.885729E-03
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quadratic elements
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear sound wave, quadratic elements on a hexahedral mesh, serial.'''
    nerr_exp = 2.124740E-02
    terr_exp = 2.124749E-02
    verr_exp = 2.126210E-02
    assert sound_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 2.124740E-02
    terr_exp = 2.124749E-02
    verr_exp = 2.126210E-02
    assert sound_setup(1,2,2,3,mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, cubic elements
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear sound wave, cubic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 4.525519E-03
    terr_exp = 4.524093E-03
    verr_exp = 4.513400E-03
    assert sound_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 4.525519E-03
    terr_exp = 4.524093E-03
    verr_exp = 4.513400E-03
    assert sound_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quartic elements
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear sound wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 3.312171E-03
    terr_exp = 3.312507E-03
    verr_exp = 3.315488E-03
    assert sound_setup(1,1,4,2,'0,2,2,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp)
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
    nerr_exp = 3.312171E-03
    terr_exp = 3.312507E-03
    verr_exp = 3.315488E-03
    assert sound_setup(1,2,4,3,'0,0,2,2,2',mf=mf,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, quadratic elements
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized sound wave, quadratic elements on a hexahedral mesh, serial.'''
    nerr_exp = 2.101794E-02
    terr_exp = 2.101794E-02
    verr_exp = 2.100566E-02
    assert sound_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 2.101794E-02
    terr_exp = 2.101794E-02
    verr_exp = 2.100566E-02
    assert sound_setup(1,2,2,3,linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized sound wave, cubic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 9.266467E-03
    terr_exp = 9.266508E-03
    verr_exp = 9.267281E-03
    assert sound_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 9.266467E-03
    terr_exp = 9.266508E-03
    verr_exp = 9.267281E-03
    assert sound_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)

#============================================================================
# Linearized runs, hexahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized sound wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    nerr_exp = 8.885322E-03
    terr_exp = 8.884684E-03
    verr_exp = 8.885686E-03
    assert sound_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True)
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
    nerr_exp = 8.885322E-03
    terr_exp = 8.884684E-03
    verr_exp = 8.885686E-03
    assert sound_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,two_temp=two_temp,hex_mesh=True)
    assert validate_sound(nerr_exp,terr_exp,verr_exp,linear=True,two_temp=two_temp)
