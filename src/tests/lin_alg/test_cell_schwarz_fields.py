#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for field based blocks (nlocal=-1 and field groups nlocal<-1) of the native
# Schwarz preconditioner in cell mode (test_cell_schwarz_fields.F90), on a coupled three field
# system, run serially and on 2 and 4 processes
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
/

&test_cell_schwarz_fields_options
 order={1}
/
"""

# Result layout (one row per number of blocks per process, 1 then 3):
#  0: layers for every block to cover the mesh, 1-2: restricted iterations with them (one group
#  per field, groups 1,1,2), 3: weighted additive with them (one group per field),
#  4-6: no field groups, layers 0-2, 7-9: one group per field, layers 0-2, 10-12: groups 1,1,2,
#  layers 0-2, 13: groups 1,1,1 of 2 (empty group), one layer, 14: before update, 15: after
#  update, 16: after delete and re-setup, 17: after update with a new pattern (one group per
#  field, one layer), 18: weighted additive and 19: additive with one group per field and one
#  layer, 20: number of blocks with one group per field
def run_cell_schwarz_fields(nprocs, nlevels, order, timeout=600):
    os.chdir(test_dir)
    if os.path.isfile('cell_schwarz_fields.results'):
        os.remove('cell_schwarz_fields.results')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, order))
    if not run_OFT("./test_cell_schwarz_fields", nprocs, timeout):
        return None
    with open('cell_schwarz_fields.results', 'r') as fid:
        nfail = int(fid.readline())
        its = [[int(val) for val in fid.readline().split()] for _ in range(2)]
    if nfail != 0:
        print("FAILED: {0} cell Schwarz field block checks failed".format(nfail))
        return None
    return its

# Iteration counts depend on the METIS partition and the DoF color tie break. Properties that hold
# by construction are checked exactly; the expected trends (dropping coupling between fields, or
# removing layers, does not help) only within a slack, as GMRES counts are not monotonic in general
# (exact reference counts are in test_cell_schwarz.py)
SLACK = 2

def check_counts(its):
    for ip, row in enumerate(its):
        nparts = [1, 3][ip]
        print("nparts={0}: {1}".format(nparts, row))
        # Covering blocks: block-Jacobi by field group, at most 3 iterations (also checked in the test)
        assert max(row[1:4]) <= 3
        # Trends: no groups <~ groups 1,1,2 <~ one group per field on the same cells, and no
        # increase with layers, each within the slack
        for k in range(3):
            assert row[4+k] <= row[10+k] + SLACK and row[10+k] <= row[7+k] + SLACK
        for base in (4, 7, 10):
            assert row[base+1] <= row[base] + SLACK and row[base+2] <= row[base+1] + SLACK
        # Empty field group: same as no field groups (also checked in the test)
        assert row[13] == row[5]
        # Update matches a fresh setup (checked in the test), re-setup and rebuild reproduce,
        # with the same blocks as the iteration count runs
        assert row[16] == row[14] and row[17] == row[14] and row[14] == row[8]
        # One block per field and builder block
        assert row[20] == 3*nparts

@pytest.mark.coverage
def test_serial():
    '''Field based cell mode Schwarz on a shared mesh: covering field blocks converge in at most 3
    iterations, iteration counts against blocks without field groups, empty groups, update and re-setup.'''
    its = run_cell_schwarz_fields(1, 3, 2)
    assert its is not None
    check_counts(its)

@pytest.mark.mpi
@pytest.mark.parametrize('nprocs', [2, 4])
def test_mpi_transfer(nprocs):
    '''Field based cell mode Schwarz on the transfer level: blocks cross process boundaries.'''
    its = run_cell_schwarz_fields(nprocs, 4, 2)
    assert its is not None
    check_counts(its)

@pytest.mark.mpi
def test_mpi_refined():
    '''Field based cell mode Schwarz one refinement above the transfer level (2 processes).'''
    its = run_cell_schwarz_fields(2, 5, 2, timeout=1200)
    assert its is not None
    check_counts(its)
