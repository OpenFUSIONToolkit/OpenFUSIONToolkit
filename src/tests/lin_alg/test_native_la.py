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
 grid_order=1
/

&cube_options
 mesh_type=1
/
"""

def run_native_la(nprocs):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    nbase = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase))
    if not run_OFT("./test_native_la", nprocs, 60):
        return False
    with open('native_la.results', 'r') as fid:
        nfail = int(fid.readline())
    if nfail != 0:
        print("FAILED: {0} native LA checks failed".format(nfail))
        return False
    return True

@pytest.mark.coverage
def test_native_la():
    '''Native transpose products (adjoint identity on a non-symmetric operator), vector norms,
    complex products and mixed real/complex dot products on a single process.'''
    assert run_native_la(1)

@pytest.mark.mpi
def test_native_la_mpi2():
    '''Native transpose products and vector norms with interface DOFs shared between 2 processes.'''
    assert run_native_la(2)

@pytest.mark.mpi
def test_native_la_mpi4():
    '''Native transpose products and vector norms with interface DOFs shared between 4 processes.'''
    assert run_native_la(4)
