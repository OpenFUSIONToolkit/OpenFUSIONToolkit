#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: linear shear Alfven wave in a triply periodic box (xMHD, H(curl) magnetic field), advanced
# half a wave period and compared with the analytic (inverted) wave; also checks the time steps and
# fields in the plot file
import functools
import pytest
from mug_testing import alfven_setup as _alfven_setup, validate_alfven

alfven_setup = functools.partial(_alfven_setup, 'hcurl')

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements
def test_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.2660712289547556E-002
    verr_exp = 2.124825E-02
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
    berr_exp = 4.2660712289547556E-002
    verr_exp = 2.124825E-02
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, cubic elements
@pytest.mark.slow
def test_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 6.370534E-03
    verr_exp = 2.357204E-03
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
    berr_exp = 6.370534E-03
    verr_exp = 2.357204E-03
    assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, tetrahedral mesh, quartic elements
@pytest.mark.slow
def test_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 4.097590E-04
    verr_exp = 2.147172E-04
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
    berr_exp = 4.097590E-04
    verr_exp = 2.147172E-04
    assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linearized runs, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, serial.'''
    berr_exp = 4.267325E-02
    verr_exp = 2.124325E-02
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a tetrahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 4.267325E-02
    verr_exp = 2.124325E-02
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs with backward Euler time stepping, tetrahedral mesh, quadratic elements
@pytest.mark.linear
def test_lin_r1_p2_euler(petsc_flag=False):
    '''Linearized Alfven wave with first-order backward Euler steps (`xmhd_lin_bdf2=F`), quadratic elements on a
    tetrahedral mesh, serial.'''
    berr_exp = 5.737182E-02
    verr_exp = 4.360428E-02
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,bdf2=False)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p2_euler_mpi(petsc_flag):
    '''Linearized Alfven wave with first-order backward Euler steps (`xmhd_lin_bdf2=F`), quadratic elements on a
    tetrahedral mesh, on two MPI processes. Run with native and PETSc linear algebra.'''
    berr_exp = 5.737182E-02
    verr_exp = 4.360428E-02
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,bdf2=False)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 6.422117E-03
    verr_exp = 2.483779E-03
    assert alfven_setup(1,1,3,2,'0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p3_mpi(petsc_flag):
    '''Linearized Alfven wave, cubic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 6.422117E-03
    verr_exp = 2.483779E-03
    assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, tetrahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, serial.'''
    berr_exp = 8.963321E-04
    verr_exp = 8.254841E-04
    assert alfven_setup(1,1,4,2,'0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_lin_r1_p4_mpi(petsc_flag):
    '''Linearized Alfven wave, quartic elements on a tetrahedral mesh with p-multigrid, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 8.963321E-04
    verr_exp = 8.254841E-04
    assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Runs with geometric multigrid, refined tetrahedral mesh, quadratic elements
def test_nl_r2_p2_gmg():
    '''Nonlinear Alfven wave, quadratic elements on a once-refined tetrahedral mesh, serial, with a multigrid
    preconditioner that includes a linear-element level on the coarse mesh (geometric and p-multigrid). Reference
    errors match a single-level (direct solve) run on the same mesh.'''
    berr_exp = 1.011290E-02
    verr_exp = 4.422328E-03
    assert alfven_setup(2,2,2,1,'2,2,2')
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.slow
@pytest.mark.mpi
def test_nl_r2_p2_gmg_mpi():
    '''Nonlinear Alfven wave, quadratic elements on a once-refined tetrahedral mesh, on two MPI processes, with
    geometric multigrid down to the distributed base mesh.'''
    berr_exp = 1.011290E-02
    verr_exp = 4.422328E-03
    assert alfven_setup(1,3,2,2,'0,0,2,2')
    assert validate_alfven(berr_exp,verr_exp)
@pytest.mark.linear
@pytest.mark.slow
def test_lin_r2_p2_gmg():
    '''Linearized Alfven wave, quadratic elements on a once-refined tetrahedral mesh, serial, with geometric and
    p-multigrid.'''
    berr_exp = 1.013647E-02
    verr_exp = 4.472189E-03
    assert alfven_setup(2,2,2,1,'2,2,2',linear=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.slow
@pytest.mark.mpi
def test_lin_r2_p2_gmg_mpi():
    '''Linearized Alfven wave, quadratic elements on a once-refined tetrahedral mesh, on two MPI processes, with
    geometric multigrid down to the distributed base mesh.'''
    berr_exp = 1.013647E-02
    verr_exp = 4.472189E-03
    assert alfven_setup(1,3,2,2,'0,0,2,2',linear=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Nonlinear runs, hexahedral mesh, quadratic elements
def test_hex_nl_r1_p2(petsc_flag=False):
    '''Nonlinear Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.818669E-02
    verr_exp = 1.818676E-02
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
    berr_exp = 1.818669E-02
    verr_exp = 1.818676E-02
    assert alfven_setup(1,2,2,3,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, cubic elements
@pytest.mark.slow
def test_hex_nl_r1_p3(petsc_flag=False):
    '''Nonlinear Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.282140E-03
    verr_exp = 2.282148E-03
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
    berr_exp = 2.282140E-03
    verr_exp = 2.282148E-03
    #assert alfven_setup(1,2,3,3,'0,0,2,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Nonlinear runs, hexahedral mesh, quartic elements
@pytest.mark.slow
def test_hex_nl_r1_p4(petsc_flag=False):
    '''Nonlinear Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 2.396858E-04
    verr_exp = 2.393954E-04
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
    berr_exp = 2.396858E-04
    verr_exp = 2.393954E-04
    #assert alfven_setup(1,2,4,3,'0,0,10,4,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',mf=mf,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp)

#============================================================================
# Linearized runs, hexahedral mesh, quadratic elements
@pytest.mark.linear
def test_hex_lin_r1_p2(petsc_flag=False):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, serial.'''
    berr_exp = 1.818575E-02
    verr_exp = 1.818575E-02
    assert alfven_setup(1,1,2,2,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
@pytest.mark.linear
@pytest.mark.mpi
@pytest.mark.coverage
@pytest.mark.parametrize("petsc_flag", (True, False))
def test_hex_lin_r1_p2_mpi(petsc_flag):
    '''Linearized Alfven wave, quadratic elements on a hexahedral mesh, on two MPI processes.
    Run with native and PETSc linear algebra.'''
    berr_exp = 1.818575E-02
    verr_exp = 1.818575E-02
    assert alfven_setup(1,2,2,3,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, cubic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p3(petsc_flag=False):
    '''Linearized Alfven wave, cubic elements on a hexahedral mesh, serial.'''
    berr_exp = 2.440052E-03
    verr_exp = 2.440116E-03
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
    berr_exp = 2.440052E-03
    verr_exp = 2.440116E-03
    #assert alfven_setup(1,2,3,3,'0,0,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,3,4,linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)

#============================================================================
# Linearized runs, hexahedral mesh, quartic elements
@pytest.mark.linear
@pytest.mark.slow
def test_hex_lin_r1_p4(petsc_flag=False):
    '''Linearized Alfven wave, quartic elements on a hexahedral mesh with p-multigrid, serial.'''
    berr_exp = 8.323704E-04
    verr_exp = 8.323704E-04
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
    berr_exp = 8.323704E-04
    verr_exp = 8.323704E-04
    #assert alfven_setup(1,2,4,3,'0,0,2,2,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert alfven_setup(1,2,4,4,'0,0,0,0,2',linear=True,petsc=petsc_flag,hex_mesh=True)
    assert validate_alfven(berr_exp,verr_exp,linear=True)
