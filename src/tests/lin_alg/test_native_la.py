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
 ref_per={2}
/
"""

def run_native_la(nprocs, ref_per='F,F,F'):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    nbase = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, nbase, ref_per))
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
    complex products and mixed real/complex dot products on a single process. Also borrowed
    (zero-copy) vector access in each access mode against get_local/restore_local, copy, maxpy
    (bitwise equal to successive add calls) and mdot (against individual dot products) for real
    and complex vectors, boundary condition operators that skip synchronization, complex stitching
    and accumulate access against get_local/restore_local(add=.TRUE.).'''
    assert run_native_la(1)

@pytest.mark.mpi
def test_native_la_mpi2():
    '''Native transpose products, vector norms, borrowed access with and without synchronization,
    copy, maxpy, mdot (real and complex), boundary condition operators, complex stitching and
    accumulate access with interface DOFs shared between 2 processes.'''
    assert run_native_la(2)

@pytest.mark.mpi
def test_native_la_mpi4():
    '''Native transpose products, vector norms, borrowed access with and without synchronization,
    copy, maxpy, mdot (real and complex), boundary condition operators, complex stitching and
    accumulate access with interface DOFs shared between 4 processes.'''
    assert run_native_la(4)

# Periodic meshes hold several images of a shared entry on one process
periodic_cases = [pytest.param('T,T,T', id='full'), pytest.param('T,F,F', id='partial')]

@pytest.mark.parametrize('ref_per', periodic_cases)
def test_native_la_periodic(ref_per):
    '''All native LA checks on a periodic mesh (`ref_per` gives the periodic directions) on a single
    process, including bitwise consistent summation stitching of periodic images.'''
    assert run_native_la(1, ref_per)

@pytest.mark.mpi
@pytest.mark.parametrize('ref_per', periodic_cases)
def test_native_la_periodic_mpi2(ref_per):
    '''All native LA checks on a periodic mesh with 2 processes, where entries can have several images
    on one process and copies on the other.'''
    assert run_native_la(2, ref_per)

@pytest.mark.mpi
@pytest.mark.parametrize('ref_per', periodic_cases)
def test_native_la_periodic_mpi4(ref_per):
    '''All native LA checks on a periodic mesh with 4 processes.'''
    assert run_native_la(4, ref_per)
