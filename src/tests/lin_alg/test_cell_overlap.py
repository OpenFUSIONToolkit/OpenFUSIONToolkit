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
 nbase={1}
/

&cube_options
 mesh_type=1
/
"""

def run_cell_overlap(nprocs, nbase, nlevels):
    os.chdir(test_dir)
    if os.path.isfile('cell_overlap.results'):
        os.remove('cell_overlap.results')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase))
    if not run_OFT("./test_cell_overlap", nprocs, 300):
        return False
    with open('cell_overlap.results', 'r') as fid:
        nfail, nchecks = [int(val) for val in fid.readline().split()]
    print("Checks: {0}, failed: {1}".format(nchecks, nfail))
    return (nfail == 0) and (nchecks > 0)

def test_serial():
    '''Cell overlap builder on a shared mesh (all levels replicated, no halo): blocks match the brute
    force layering for 0-2 layers, point/edge/face adjacency, 1 and 3 blocks, single and composite FE.'''
    assert run_cell_overlap(1, 3, 3)

@pytest.mark.mpi
@pytest.mark.parametrize('nprocs', [2, 4])
def test_mpi_transfer(nprocs):
    '''Cell overlap builder on the transfer level (no refinement above it): halo cells across
    processes match the brute force layering exactly.'''
    assert run_cell_overlap(nprocs, 3, 4)

@pytest.mark.mpi
@pytest.mark.parametrize('nprocs', [2, 4])
def test_mpi_refined(nprocs):
    '''Cell overlap builder one refinement above the transfer level: halo candidates from the base
    mesh bound are sufficient and blocks match the brute force layering exactly.'''
    assert run_cell_overlap(nprocs, 3, 5)
