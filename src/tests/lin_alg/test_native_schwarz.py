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

def run_schwarz(nprocs, ref_per='F,F,F'):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, ref_per))
    if not run_OFT("./test_native_schwarz", nprocs, 120):
        return None
    with open('schwarz.results', 'r') as fid:
        nfail = int(fid.readline())
        its = [int(val) for val in fid.readline().split()]
    if nfail != 0:
        print("FAILED: {0} Schwarz checks failed".format(nfail))
        return None
    return its

def check_overlap_iterations(its):
    # One subdomain per process: restricted (its[0:3]) and additive (its[0], its[3:5]) counts
    # must not increase with overlap; four blocks per process: restricted overlap 1 (its[7])
    # against no overlap (its[5]). Plain additive combination of overlapping blocks within a
    # process (its[6]) counts overlap rows twice and need not improve on block-Jacobi, so it is
    # only required to converge.
    print("Iterations: ", its)
    assert its[1] <= its[0]
    assert its[2] <= its[1]
    assert its[3] <= its[0]
    assert its[4] <= its[3]
    assert its[7] <= its[5]

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

@pytest.mark.parametrize('nprocs', [pytest.param(1, id='serial'), pytest.param(2, id='mpi2', marks=pytest.mark.mpi)])
def test_schwarz_periodic(nprocs):
    '''Overlapping Schwarz with several blocks per process on a partially periodic cube (periodic images
    on one process): overlap 0 reproduces block-Jacobi and overlap 1 reproduces block-Jacobi with
    loverlap, bit for bit.'''
    its = run_schwarz(nprocs, 'T,F,F')
    assert its is not None
