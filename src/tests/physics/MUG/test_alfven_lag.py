#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear shear Alfven wave in a triply periodic box (xMHD with all fields on Lagrange
# elements), advanced half a wave period and compared with the analytic (inverted) wave; also checks the
# time steps and fields in the plot file
import functools
import pytest
from mug_testing import alfven_setup as _alfven_setup, validate_alfven

alfven_setup = functools.partial(_alfven_setup, 'lagrange')

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p2_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 4.2378931381112289E-002
    verr_exp = 4.1830283541690064E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, cubic elements
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p3_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.1613939097864013E-003
    verr_exp = 2.0399660410128605E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, quartic elements
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_nl_r1_p4_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.6353342214151861E-004
    verr_exp = 2.4711224138100239E-004
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linearized runs, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 4.2445515350299369E-002
    verr_exp = 4.1907987520465033E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 2.1743513983662793E-003
    verr_exp = 2.0589634472042181E-003
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 3.3133689642593500E-004
    verr_exp = 3.1807433060756504E-004
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Nonlinear runs, hexahedral mesh, quadratic elements
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,1,2,2,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p2_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.8204622112726550E-002
    verr_exp = 1.8204690721409598E-002
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, cubic elements
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p3_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, cubic elements on a hexahedral mesh, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 2.2675346818877107E-003
    verr_exp = 2.2675252946567641E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quartic elements
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    #assert alfven_setup(1,1,4,2,'0,10,4,2',petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("mf", (False, True))
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_nl_r1_p4_mpi(mf,petsc_flag):
    '''Nonlinear Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Also covers matrix-free Newton-Krylov and (default variant only) PETSc runs.'''
    if mf and petsc_flag:
        pytest.skip()
    berr_exp = 1.7790763796394384E-004
    verr_exp = 1.7781738147450092E-004
    #assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linearized runs, hexahedral mesh, quadratic elements
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 1.8185749172894013E-002
    verr_exp = 1.8185749172889051E-002
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    #assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,3,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a hexahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 2.2830097215737651E-003
    verr_exp = 2.2830075204362107E-003
    #assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    #assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,1,4,3,'0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 2.6757952245217109E-004
    verr_exp = 2.6757961586721258E-004
    #assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
