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
 ref_per={2}
/

&test_cell_schwarz_options
 order={1}
/
"""

# Result layout (one row per number of blocks per process, 1 then 3):
#  0: layers for every block to cover the mesh, 1-2: restricted/additive iterations with them,
#  3-11: cell mode iterations for (point, edge, face) x layers 0-2,
#  12-14: graph overlap 0-2, 15: before update, 16: after update, 17: after delete and
#  re-setup, 18: after update with a new pattern
def run_cell_schwarz(nprocs, nlevels, order, timeout=600, ref_per='F,F,F'):
    os.chdir(test_dir)
    if os.path.isfile('cell_schwarz.results'):
        os.remove('cell_schwarz.results')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, order, ref_per))
    if not run_OFT("./test_cell_schwarz", nprocs, timeout):
        return None
    with open('cell_schwarz.results', 'r') as fid:
        nfail = int(fid.readline())
        its = [[int(val) for val in fid.readline().split()] for _ in range(2)]
    if nfail != 0:
        print("FAILED: {0} cell Schwarz checks failed".format(nfail))
        return None
    return its

def check_counts(its, expected):
    for ip, (row, ref) in enumerate(zip(its, expected)):
        print("nparts={0}: {1}".format([1, 3][ip], row))
        # Exactness: covering layers and one restricted iteration
        assert row[0] == ref[0]
        assert row[1] == 1
        # Iteration counts within 5% of the expected values
        for val, rval in zip(row[3:15], ref[3:15]):
            assert abs(val-rval) <= max(1, 0.05*rval)
        # Point adjacency: non-increasing with layers; face <= edge <= point overlap per layer
        assert row[4] <= row[3] and row[5] <= row[4]
        for k in range(3):
            assert row[3+k] <= row[6+k] <= row[9+k]
        # One layer of point-adjacent cells is at least as good as one graph layer
        assert row[4] <= row[13]
        # Update matches a fresh setup (checked in the test), re-setup and rebuild reproduce
        assert row[17] == row[15] and row[18] == row[15]

# Expected values (OMP_NUM_THREADS=2; counts may differ by one with the thread count, as the
# reductions and so the final residuals differ in round-off): (nprocs, nlevels, order) -> rows
# for 1 and 3 blocks per process
expected = {
    (1, 3, 3): [[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [3, 1, 1, 14, 7, 4, 14, 8, 5, 14, 10, 8, 40, 13, 7, 7, 7, 7, 7]],
    (2, 4, 3): [[3, 1, 1, 21, 7, 4, 21, 9, 6, 21, 12, 9, 39, 10, 5, 7, 7, 7, 7],
                [4, 1, 1, 22, 9, 6, 22, 10, 8, 22, 13, 11, 47, 39, 39, 9, 9, 9, 9]],
    (4, 4, 2): [[4, 1, 1, 19, 8, 5, 19, 10, 7, 19, 13, 10, 34, 11, 7, 8, 8, 8, 8],
                [5, 1, 1, 19, 9, 6, 19, 11, 8, 19, 13, 11, 39, 35, 34, 9, 9, 9, 9]],
    (2, 5, 2): [[6, 1, 1, 25, 10, 7, 25, 12, 8, 25, 14, 12, 46, 13, 8, 10, 10, 10, 10],
                [8, 1, 1, 25, 12, 9, 25, 14, 11, 25, 16, 14, 56, 47, 48, 12, 12, 12, 12]],
    (4, 5, 2): [[8, 1, 1, 31, 14, 9, 31, 15, 11, 31, 18, 15, 52, 17, 10, 14, 14, 14, 14],
                [10, 1, 1, 33, 14, 10, 33, 15, 12, 33, 19, 16, 62, 56, 56, 14, 14, 14, 14]],
}

# Cube periodic in x (`ref_per=T,F,F`): blocks and graph subdomains grow across the periodic seam
expected_periodic = {
    (1, 3, 3): [[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [4, 1, 1, 35, 14, 11, 35, 16, 11, 35, 19, 16, 64, 21, 12, 14, 14, 14, 14]],
    (2, 4, 3): [[3, 1, 1, 41, 15, 9, 41, 18, 11, 41, 22, 18, 58, 17, 9, 15, 15, 15, 15],
                [4, 1, 1, 43, 16, 11, 43, 19, 14, 43, 25, 19, 76, 61, 55, 16, 16, 16, 16]],
    (4, 4, 2): [[4, 1, 1, 34, 15, 9, 34, 17, 12, 34, 20, 17, 57, 17, 11, 15, 15, 15, 15],
                [5, 1, 1, 38, 16, 11, 38, 18, 14, 38, 23, 18, 66, 59, 57, 16, 16, 16, 16]],
    (2, 5, 2): [[5, 1, 1, 46, 19, 13, 46, 23, 17, 46, 27, 23, 60, 22, 14, 19, 19, 19, 19],
                [9, 1, 1, 52, 22, 14, 52, 26, 18, 52, 31, 26, 82, 66, 63, 22, 22, 22, 22]],
}

@pytest.mark.coverage
def test_serial():
    '''Cell mode Schwarz on a shared mesh (third order): exact blocks converge in one iteration,
    iteration counts with point, edge and face layers against graph overlap, numeric update and
    re-setup.'''
    its = run_cell_schwarz(1, 3, 3)
    assert its is not None
    check_counts(its, expected[(1, 3, 3)])

@pytest.mark.mpi
@pytest.mark.parametrize('nprocs,order', [(2, 3), (4, 2)])
def test_mpi_transfer(nprocs, order):
    '''Cell mode Schwarz on the transfer level: blocks cross process boundaries (third order on 2
    processes, second order on 4 to limit run time).'''
    its = run_cell_schwarz(nprocs, 4, order)
    assert its is not None
    check_counts(its, expected[(nprocs, 4, order)])

@pytest.mark.mpi
def test_mpi_refined():
    '''Cell mode Schwarz one refinement above the transfer level (second order, 2 processes).'''
    its = run_cell_schwarz(2, 5, 2)
    assert its is not None
    check_counts(its, expected[(2, 5, 2)])

@pytest.mark.slow
@pytest.mark.mpi
def test_mpi4_refined():
    '''Cell mode Schwarz one refinement above the transfer level (second order, 4 processes).'''
    its = run_cell_schwarz(4, 5, 2, timeout=2400)
    assert its is not None
    check_counts(its, expected[(4, 5, 2)])

@pytest.mark.parametrize('nprocs,nlevels,order', [pytest.param(1, 3, 3, id='serial'),
                                                  pytest.param(2, 4, 3, id='mpi2', marks=pytest.mark.mpi),
                                                  pytest.param(4, 4, 2, id='mpi4', marks=pytest.mark.mpi),
                                                  pytest.param(2, 5, 2, id='mpi2_refined',
                                                               marks=[pytest.mark.mpi, pytest.mark.slow])])
def test_periodic(nprocs, nlevels, order):
    '''Cell and graph mode Schwarz on a cube periodic in x (periodic images of seam DoFs on one
    process, cells adjacent across the seam): the same checks as on the plain cube.'''
    its = run_cell_schwarz(nprocs, nlevels, order, timeout=1200, ref_per='T,F,F')
    assert its is not None
    check_counts(its, expected_periodic[(nprocs, nlevels, order)])
