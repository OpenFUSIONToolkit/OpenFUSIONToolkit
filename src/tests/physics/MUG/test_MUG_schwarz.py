#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: native overlapping Schwarz preconditioner with blocks grown by layers of cells
# (`oft_in_xmhd_mg_cell.xml`: one layer of point-adjacent cells, LU sub-solves) as the GMRES smoother
# of every multigrid level, on two MPI processes, for both xMHD representations
#
# The periodic runs use the triply periodic sound wave of the other MUG tests (cubic elements,
# p-multigrid on the distributed mesh level) and must reproduce their reference errors, with blocks
# growing across the periodic seams. The runs in a box with walls also cover the transfer level; the
# wave is then not a solution, so their errors are regression values only: they come from runs with
# the block-Jacobi smoother of the other MUG tests (`oft_in_xmhd_mg.xml`, also run here), which the
# Schwarz runs must reproduce.
import pytest
from mug_testing import REPRESENTATIONS, sound_setup, validate_sound

ERRORS = {
    'hcurl': (1.324993E+00, 1.327518E+00, 4.073417E+00),
    'lagrange': (1.324987E+00, 1.327513E+00, 4.073380E+00),
}
# Reference errors of test_sound.py and test_sound_lag.py (test_nl_r1_p3_mpi)
PERIODIC_ERRORS = {
    'hcurl': (3.067332E-03, 3.070297E-03, 3.355874E-03),
    'lagrange': (3.067374E-03, 3.070338E-03, 3.355883E-03),
}

#============================================================================
# Nonlinear runs, triply periodic tetrahedral mesh, cubic elements, multigrid on the distributed level
@pytest.mark.mpi
@pytest.mark.parametrize("representation", REPRESENTATIONS)
@pytest.mark.parametrize("nparts", (1, 2))
def test_nl_r1_p3_periodic_mpi(representation, nparts):
    '''Nonlinear sound wave on the triply periodic box, cubic elements on a tetrahedral mesh with
    p-multigrid, on two MPI processes, with cell-mode Schwarz smoothing (one or two blocks per process).
    Reproduces the reference errors of the block-Jacobi runs of test_sound.py/test_sound_lag.py.'''
    assert sound_setup(representation,1,2,3,3,'0,0,2,2',xml='oft_in_xmhd_mg_cell.xml',nparts=nparts)
    assert validate_sound(*PERIODIC_ERRORS[representation])

#============================================================================
# Nonlinear runs, tetrahedral mesh, quadratic elements, multigrid on the distributed levels
@pytest.mark.mpi
@pytest.mark.parametrize("representation", REPRESENTATIONS)
@pytest.mark.parametrize("pre,nparts", (('block_jacobi', 1), ('cell', 1), ('cell', 2)))
def test_nl_r1_p2_walls_mpi(representation, pre, nparts):
    '''Nonlinear sound wave in a box with walls, quadratic elements on a tetrahedral mesh with p-multigrid,
    on two MPI processes. The smoother preconditioner is block-Jacobi (reference) or cell-mode Schwarz with
    one block per process or two blocks per process (`xmhd_nparts`, the builder's cell parts).'''
    xml = 'oft_in_xmhd_mg_cell.xml' if pre == 'cell' else None
    assert sound_setup(representation,1,2,2,2,'0,2,2',periodic=False,xml=xml,nparts=nparts)
    assert validate_sound(*ERRORS[representation])
