#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for the native Schwarz preconditioner in cell mode
# (test_cell_schwarz.F90), run serially and on 2 and 4 processes, on plain and
# x-periodic cubes
#------------------------------------------------------------------------------
import os
import sys
import warnings
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
#  re-setup, 18: after update with a new pattern, 19: weighted additive iterations with the
#  covering layers, 20: weighted additive with one layer of point adjacent cells
# The test also writes a signature of the partition behind each row: hashes of the cell
# partition and of the graph mode blocks within each process, of the DoF colors and of the
# decomposition of the mesh between processes (the rank of each cell), each over the global
# indices of the cells or rows. The cell partition, graph blocks and decomposition come from
# METIS, and so depend on the METIS build. The DoF colors also depend on the tie breaking between
# parts in mesh_local_partition, which uses the compiler's random number generator, so they vary
# with the compiler and the thread count; they can change the covering layers and iteration
# counts as much as the METIS partitions do. The full signature selects the reference row.
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
        sig = [int(val) for val in fid.readline().split()]
    if nfail != 0:
        print("FAILED: {0} cell Schwarz checks failed".format(nfail))
        return None
    return its, [sig[0:4], sig[4:8]]

def check_counts(result, references):
    its, sigs = result
    for ip, (row, sig, refs) in enumerate(zip(its, sigs, references)):
        nparts = [1, 3][ip]
        print("nparts={0}: {1}".format(nparts, row))
        print("nparts={0} partition signature (cells, DoF colors, graph blocks, ranks): {1}".format(nparts, sig))
        # The covering layers and iteration counts depend on the partition: compare them with a
        # reference row only for the partition (signature) it was recorded with
        matches = [(src, ref) for src, ref_sig, ref in refs if ref_sig == sig]
        if matches:
            src, ref = matches[0]
            print("nparts={0}: reference row from {1}".format(nparts, src))
            # Covering layers exact, iteration counts within 5% of the reference values
            assert row[0] == ref[0]
            for val, rval in zip(row[3:15], ref[3:15]):
                assert abs(val-rval) <= max(1, 0.05*rval)
            assert abs(row[20]-ref[20]) <= max(1, 0.05*ref[20])
        else:
            accepted = [ref_sig for _, ref_sig, _ in refs]
            print("=" * 80)
            print("NOTE: nparts={0}: partition signature {1} has no reference row (accepted: {2}); "
                  "checking partition-independent properties only".format(nparts, sig, accepted))
            print("NOTE:   row: {0}".format(row))
            for src, ref_sig, ref in refs:
                print("NOTE:   reference {0} ({1}): {2}".format(ref_sig, src, ref))
            print("=" * 80)
            # Also as a warning: pytest shows its warnings summary for passing tests, whose
            # output is otherwise captured
            warnings.warn(UserWarning("nparts={0}: no reference row for partition signature {1} (cells, DoF "
                                      "colors, graph blocks, ranks), reference table not applied; row {2}; "
                                      "accepted signatures {3}".format(nparts, sig, row, accepted)))
        # Exactness: one restricted iteration with the covering layers
        assert row[1] == 1
        # Point adjacency: non-increasing with layers; face <= edge <= point overlap per layer
        assert row[4] <= row[3] and row[5] <= row[4]
        for k in range(3):
            assert row[3+k] <= row[6+k] <= row[9+k]
        # One layer of point-adjacent cells is at least as good as one graph layer
        assert row[4] <= row[13]
        # Update matches a fresh setup (checked in the test), re-setup and rebuild reproduce
        assert row[17] == row[15] and row[18] == row[15]
        # Weighted additive: exact with covering blocks (also checked in the test); convergence
        # with one layer is checked in the test
        assert row[19] == 1

# Reference rows: (nprocs, nlevels, order) -> references for 1 and 3 blocks per process, each
# (source, partition signature [cells, DoF colors, graph blocks, ranks], row), recorded with
# OMP_NUM_THREADS=2 (counts may differ by one with the thread count, as the reductions and so the
# final residuals differ in round-off, and the DoF colors depend on the thread count). CI builds
# METIS 5.1.0 with build_libs.py (GKRAND) and partitions differently from the Ubuntu system METIS
# 5.1.0 (built without GKRAND); gfortran and Intel ifx CI agree on the METIS partitions but not
# always on the DoF colors. The single block serial rows (no partition within the process) are
# the same on all of these; the slow configurations, which CI does not run, have the system
# METIS reference only.
CI = 'gfortran and Intel ifx CI (build_libs.py METIS 5.1.0, GKRAND; jobs 139762574767, 139762574553)'
GF = 'gfortran CI (GCC 12, MPICH, build_libs.py METIS 5.1.0, GKRAND; job 139762574767)'
IFX = 'Intel ifx CI (Intel 2025.2, MPICH, build_libs.py METIS 5.1.0, GKRAND; job 139762574553)'
SYS = 'gfortran with the Ubuntu system METIS 5.1.0 (no GKRAND)'

expected = {
    (1, 3, 3): [
        [  # 1 block per process
            (SYS, [1369331334, 402523094, 0, 1369331334],
             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        ],
        [  # 3 blocks per process
            (GF, [1857716918, 1372209359, 1301423679, 1369331334],
             [4, 1, 1, 13, 7, 4, 13, 9, 6, 13, 10, 9, 40, 12, 7, 7, 7, 7, 7, 1, 11]),
            (IFX, [1857716918, 27372501, 1301423679, 1369331334],
             [4, 1, 1, 13, 7, 4, 13, 9, 6, 13, 10, 9, 40, 12, 7, 7, 7, 7, 7, 1, 11]),
            (SYS, [864203196, 1564514225, 434180095, 1369331334],
             [3, 1, 1, 14, 7, 4, 14, 8, 5, 14, 10, 8, 40, 13, 7, 7, 7, 7, 7, 1, 11]),
        ],
    ],
    (2, 4, 3): [
        [  # 1 block per process
            (CI, [1369331334, 182997468, 0, 1845466281],
             [3, 1, 1, 21, 7, 4, 21, 10, 6, 21, 12, 10, 39, 11, 5, 7, 7, 7, 7, 1, 7]),
            (SYS, [1369331334, 2061508305, 0, 1767070778],
             [3, 1, 1, 21, 7, 4, 21, 9, 6, 21, 12, 9, 39, 10, 5, 7, 7, 7, 7, 1, 7]),
        ],
        [  # 3 blocks per process
            (GF, [1671197774, 1549941152, 434080154, 1845466281],
             [4, 1, 1, 22, 9, 6, 22, 11, 8, 22, 14, 11, 47, 40, 40, 9, 9, 9, 9, 1, 13]),
            (IFX, [1671197774, 1697198293, 434080154, 1845466281],
             [4, 1, 1, 22, 9, 6, 22, 11, 8, 22, 14, 11, 47, 40, 40, 9, 9, 9, 9, 1, 13]),
            (SYS, [1004189096, 445042440, 1018507348, 1767070778],
             [4, 1, 1, 22, 9, 6, 22, 10, 8, 22, 13, 11, 47, 39, 39, 9, 9, 9, 9, 1, 13]),
        ],
    ],
    (4, 4, 2): [
        [  # 1 block per process
            (CI, [1369331334, 1391291925, 0, 1045781703],
             [4, 1, 1, 19, 9, 5, 19, 11, 7, 19, 14, 10, 33, 12, 7, 9, 9, 9, 9, 1, 12]),
            (SYS, [1369331334, 1518341197, 0, 1643773048],
             [4, 1, 1, 19, 8, 5, 19, 10, 7, 19, 13, 10, 34, 11, 7, 8, 8, 8, 8, 1, 11]),
        ],
        [  # 3 blocks per process
            (GF, [927193452, 1730224380, 992815505, 1045781703],
             [5, 1, 1, 20, 9, 7, 20, 11, 8, 20, 14, 11, 40, 33, 33, 9, 9, 9, 9, 1, 13]),
            (IFX, [927193452, 1189447766, 992815505, 1045781703],
             [5, 1, 1, 19, 9, 7, 19, 11, 8, 19, 14, 10, 40, 33, 33, 9, 9, 9, 9, 1, 13]),
            (SYS, [802654272, 2067542128, 2008392473, 1643773048],
             [5, 1, 1, 19, 9, 6, 19, 11, 8, 19, 13, 11, 39, 35, 34, 9, 9, 9, 9, 1, 12]),
        ],
    ],
    (2, 5, 2): [
        [  # 1 block per process
            (CI, [699839912, 1997069520, 0, 583289318],
             [6, 1, 1, 25, 11, 7, 25, 13, 9, 25, 15, 13, 45, 14, 8, 11, 11, 11, 11, 1, 11]),
            (SYS, [699839912, 871655633, 0, 2015186007],
             [6, 1, 1, 25, 10, 7, 25, 12, 8, 25, 14, 12, 46, 13, 8, 10, 10, 10, 10, 1, 11]),
        ],
        [  # 3 blocks per process
            (GF, [1957863219, 9236773, 35482575, 583289318],
             [9, 1, 1, 26, 12, 9, 26, 14, 11, 26, 16, 14, 54, 45, 46, 12, 12, 12, 12, 1, 16]),
            (IFX, [1957863219, 1450056773, 35482575, 583289318],
             [9, 1, 1, 26, 13, 9, 26, 14, 11, 26, 16, 14, 54, 45, 46, 13, 13, 13, 13, 1, 16]),
            (SYS, [599442501, 2016617584, 817245016, 2015186007],
             [8, 1, 1, 25, 12, 9, 25, 14, 11, 25, 16, 14, 56, 47, 48, 12, 12, 12, 12, 1, 15]),
        ],
    ],
    (4, 5, 2): [
        [  # 1 block per process
            (SYS, [699839912, 403039421, 0, 1673343115],
             [8, 1, 1, 31, 14, 9, 31, 15, 11, 31, 18, 15, 52, 17, 10, 14, 14, 14, 14, 1, 15]),
        ],
        [  # 3 blocks per process
            (SYS, [1936146723, 1360605893, 1953782997, 1673343115],
             [10, 1, 1, 33, 14, 10, 33, 15, 12, 33, 19, 16, 62, 56, 56, 14, 14, 14, 14, 1, 18]),
        ],
    ],
}

# Cube periodic in x (`ref_per=T,F,F`): blocks and graph subdomains grow across the periodic seam
# (only x: with `ref_per=T,T,F` this coarse cube aborts in the mesh setup, in
# mesh_global_igrnd, before any Schwarz code runs; full periodicity is covered by
# test_native_schwarz.py and physics/MUG/test_MUG_schwarz.py)
expected_periodic = {
    (1, 3, 3): [
        [  # 1 block per process
            (SYS, [1145285914, 1921186034, 0, 1145285914],
             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        ],
        [  # 3 blocks per process
            (GF, [1612309209, 1623915770, 1987179523, 1145285914],
             [3, 1, 1, 33, 14, 11, 33, 16, 12, 33, 19, 16, 69, 20, 11, 14, 14, 14, 14, 1, 18]),
            (IFX, [1612309209, 1918727846, 1987179523, 1145285914],
             [4, 1, 1, 35, 14, 11, 35, 16, 12, 35, 20, 16, 69, 20, 11, 14, 14, 14, 14, 1, 19]),
            (SYS, [133795671, 356876798, 1648869163, 1145285914],
             [4, 1, 1, 35, 14, 11, 35, 16, 11, 35, 19, 16, 64, 21, 12, 14, 14, 14, 14, 1, 18]),
        ],
    ],
    (2, 4, 3): [
        [  # 1 block per process
            (CI, [1145285914, 280817119, 0, 367133601],
             [3, 1, 1, 21, 9, 5, 21, 11, 7, 21, 13, 11, 56, 17, 9, 9, 9, 9, 9, 1, 9]),
            (SYS, [1145285914, 948201768, 0, 996153002],
             [3, 1, 1, 41, 15, 9, 41, 18, 11, 41, 22, 18, 58, 17, 9, 15, 15, 15, 15, 1, 15]),
        ],
        [  # 3 blocks per process
            (GF, [370348049, 98060633, 2075970286, 367133601],
             [4, 1, 1, 43, 16, 12, 43, 19, 14, 43, 24, 19, 75, 57, 59, 16, 16, 16, 16, 1, 21]),
            (IFX, [370348049, 1636848999, 2075970286, 367133601],
             [4, 1, 1, 43, 16, 12, 43, 19, 14, 43, 24, 19, 75, 57, 59, 16, 16, 16, 16, 1, 21]),
            (SYS, [885882329, 1600064468, 834950068, 996153002],
             [4, 1, 1, 43, 16, 11, 43, 19, 14, 43, 25, 19, 76, 61, 55, 16, 16, 16, 16, 1, 22]),
        ],
    ],
    (4, 4, 2): [
        [  # 1 block per process
            (CI, [1145285914, 272438151, 0, 284663625],
             [4, 1, 1, 36, 15, 11, 36, 17, 13, 36, 21, 17, 56, 18, 12, 15, 15, 15, 15, 1, 18]),
            (SYS, [1145285914, 208192849, 0, 1956384897],
             [4, 1, 1, 34, 15, 9, 34, 17, 12, 34, 20, 17, 57, 17, 11, 15, 15, 15, 15, 1, 17]),
        ],
        [  # 3 blocks per process
            (GF, [1337460933, 1825505932, 1256372477, 284663625],
             [5, 1, 1, 38, 16, 11, 38, 18, 14, 38, 23, 18, 63, 57, 56, 16, 16, 16, 16, 1, 19]),
            (IFX, [1337460933, 1851534198, 1256372477, 284663625],
             [5, 1, 1, 38, 16, 12, 38, 18, 14, 38, 23, 18, 63, 57, 56, 16, 16, 16, 16, 1, 19]),
            (SYS, [1720874087, 119415708, 1493531299, 1956384897],
             [5, 1, 1, 38, 16, 11, 38, 18, 14, 38, 23, 18, 66, 59, 57, 16, 16, 16, 16, 1, 20]),
        ],
    ],
    (2, 5, 2): [
        [  # 1 block per process
            (SYS, [355298977, 1415024166, 0, 495323792],
             [5, 1, 1, 46, 19, 13, 46, 23, 17, 46, 27, 23, 60, 22, 14, 19, 19, 19, 19, 1, 20]),
        ],
        [  # 3 blocks per process
            (SYS, [484386452, 1974824886, 1957902133, 495323792],
             [9, 1, 1, 52, 22, 14, 52, 26, 18, 52, 31, 26, 82, 66, 63, 22, 22, 22, 22, 1, 24]),
        ],
    ],
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
