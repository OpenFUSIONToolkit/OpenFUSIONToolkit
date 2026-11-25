#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for the native overlapping Schwarz preconditioner in graph mode
# (test_native_schwarz.F90), run serially and on 2 and 4 processes, on plain, x-periodic
# and fully periodic cubes
#------------------------------------------------------------------------------
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={0}
 nbase=3
 grid_order=1
/

&cube_options
 mesh_type=1
 ref_per={1}
/
"""

def run_schwarz(nprocs, ref_per='F,F,F', timeout=240):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, ref_per))
    if not run_OFT("./test_native_schwarz", nprocs, timeout):
        return None
    with open('schwarz.results', 'r') as fid:
        nfail = int(fid.readline())
        its = [int(val) for val in fid.readline().split()]
    if nfail != 0:
        print("FAILED: {0} Schwarz checks failed".format(nfail))
        return None
    return its

def check_overlap_iterations(its):
    # One subdomain per process: restricted (its[0:3]), additive (its[0], its[3:5]) and
    # weighted additive (its[0], its[8:10]) counts must not increase with overlap; four blocks
    # per process: restricted (its[7]) and weighted additive (its[10]) overlap 1 against no
    # overlap (its[5]). Plain additive combination of overlapping blocks within a process
    # (its[6]) counts overlap rows twice and need not improve on block-Jacobi, so it is only
    # required to converge.
    print("Iterations: ", its)
    assert its[1] <= its[0]
    assert its[2] <= its[1]
    assert its[3] <= its[0]
    assert its[4] <= its[3]
    assert its[7] <= its[5]
    assert its[8] <= its[0]
    assert its[9] <= its[8]
    assert its[10] <= its[5]

@pytest.mark.coverage
def test_schwarz():
    '''Overlapping Schwarz on one process: overlap 0 and local overlap 1 reproduce block-Jacobi
    (without and with loverlap) bit for bit, subdomain matrices match the assembled matrix, and
    GMRES iteration counts do not increase with overlap.'''
    its = run_schwarz(1)
    assert its is not None
    check_overlap_iterations(its)

@pytest.mark.mpi
def test_schwarz_mpi2():
    '''Overlapping Schwarz on 2 processes: subdomains grown across processes (rows held only by
    the other process are transferred) match the independently assembled global matrix, overlap 0
    reproduces block-Jacobi, and restricted and additive iteration counts decrease with overlap.'''
    its = run_schwarz(2)
    assert its is not None
    check_overlap_iterations(its)
    assert its[1] < its[0]

@pytest.mark.mpi
def test_schwarz_mpi4():
    '''Overlapping Schwarz on 4 processes, as for 2 processes.'''
    its = run_schwarz(4)
    assert its is not None
    check_overlap_iterations(its)
    assert its[1] < its[0]

# The 4 process periodic case needs about twice as many GMRES iterations as the non-periodic one and
# takes about 2.5 minutes with 2 threads per process on 4 cores
@pytest.mark.parametrize('nprocs', [pytest.param(1, id='serial'), pytest.param(2, id='mpi2', marks=pytest.mark.mpi),
                                    pytest.param(4, id='mpi4', marks=[pytest.mark.mpi, pytest.mark.slow])])
def test_schwarz_periodic(nprocs):
    '''Overlapping Schwarz on a partially periodic cube (periodic images of seam DOFs on one process):
    the same checks as on the plain cube, with subdomains grown across periodic seams and processes.'''
    its = run_schwarz(nprocs, 'T,F,F', timeout=600)
    assert its is not None
    check_overlap_iterations(its)
    if nprocs > 1:
        assert its[1] < its[0]

@pytest.mark.parametrize('nprocs', [pytest.param(1, id='serial'), pytest.param(2, id='mpi2', marks=pytest.mark.mpi),
                                    pytest.param(4, id='mpi4', marks=pytest.mark.mpi)])
def test_schwarz_fully_periodic(nprocs):
    '''Overlapping subdomains on a fully periodic cube (images in all three directions, up to eight per
    DOF) match the assembled global matrix. The periodic Laplacian is singular, so no solves are checked.'''
    its = run_schwarz(nprocs, 'T,T,T')
    assert its is not None
