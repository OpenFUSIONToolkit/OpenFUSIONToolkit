#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for oft_nksolver Eisenstat-Walker adaptive forcing
# (test_nksolver.F90), run with Native and PETSc GMRES, serially and on 2 ranks
#------------------------------------------------------------------------------
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={0}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=3
 nbase={1}
 grid_order=1
/

&cube_options
 mesh_type=1
/
"""


def run_nksolver(petsc=False, nprocs=1):
    os.chdir(test_dir)
    petsc_flag = ('T' if petsc else 'F')
    # Distribute finest level for MPI runs
    nbase = 2 if nprocs > 1 else 3
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(petsc_flag, nbase))
    if not run_OFT("./test_nksolver oft.in", nprocs, 60):
        return False
    with open('nksolver.results', 'r') as fid:
        nfail = int(fid.readline())
        nlits_fixed, lits_fixed = [int(x) for x in fid.readline().split()]
        nlits_ew, lits_ew = [int(x) for x in fid.readline().split()]
        diff_norm = float(fid.readline())
    print("Fixed: nlits={0}, lits={1} | EW: nlits={2}, lits={3} | rel_diff={4:.3e}".format(
        nlits_fixed, lits_fixed, nlits_ew, lits_ew, diff_norm))
    if nfail != 0:
        print("FAILED: {0} nksolver checks failed".format(nfail))
        return False
    if lits_ew >= lits_fixed:
        print("FAILED: EW did not reduce linear iterations ({0} >= {1})".format(lits_ew, lits_fixed))
        return False
    if diff_norm > 1.e-8:
        print("FAILED: Relative solution difference too large ({0})".format(diff_norm))
        return False
    return True


@pytest.mark.parametrize("petsc_flag", (
    pytest.param(False, marks=pytest.mark.coverage),
    True,
))
def test_nksolver_ew(petsc_flag):
    '''Solve a nonlinear reaction-diffusion equation with oft_nksolver using fixed tight
    GMRES tolerances and opt-in Eisenstat-Walker Choice 2 forcing on a single process.
    Verifies identical nonlinear convergence, strictly fewer total linear iterations with
    EW forcing, and restoration of J_inv tolerances on converged and max-iteration exits.'''
    assert run_nksolver(petsc=petsc_flag, nprocs=1)


@pytest.mark.mpi
@pytest.mark.parametrize("petsc_flag", (False, True))
def test_nksolver_ew_mpi(petsc_flag):
    '''Solve a nonlinear reaction-diffusion equation with oft_nksolver using fixed tight
    GMRES tolerances and opt-in Eisenstat-Walker Choice 2 forcing across 2 MPI ranks.
    Verifies identical nonlinear convergence, strictly fewer total linear iterations with
    EW forcing, and restoration of J_inv tolerances on converged and max-iteration exits.'''
    assert run_nksolver(petsc=petsc_flag, nprocs=2)
