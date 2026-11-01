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
 meshname='test'
 cad_type={3}
 nlevels={1}
 nbase={0}
 grid_order={2}
 jac_ratio_tol=30.
/

&sphere_options
 mesh_type={4}
/

&test_mapping_options
 order={5}
/
"""

# Common setup function and process handling
def mapping_setup(nbase, nlevels, order, cad_type, exec_suffix, grid_type=1, elem_order=1, catch_warning=True):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, order, cad_type, grid_type, elem_order))
    return run_OFT("./test_mapping_" + exec_suffix, nproc, 60)

#
def check_result(filename):
    retval = True
    fid = open(filename,'r')
    fail_count = int(fid.readline())
    fid.close()
    if fail_count > 0:
        print("FAILED: # of errors = {0}".format(fail_count))
        retval = False
    return retval

#============================================================================
# Test runner for find_cell test on linear elements
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
def test_find_lin(grid_type):
    assert mapping_setup(1,1,1,91,'find',grid_type=grid_type)
    assert check_result('mapping_find.results')

#============================================================================
# Test runner for find_cell test on quadratic elements
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
def test_find_quad(grid_type):
    assert mapping_setup(1,1,2,91,'find',grid_type=grid_type)
    assert check_result('mapping_find.results')

#============================================================================
# Test runner for jacobian tests on linear elements
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("elem_order", (2, 3, 4))
def test_jac_lin(grid_type, elem_order):
    assert mapping_setup(1,1,1,91,'jac',grid_type=grid_type,elem_order=elem_order,catch_warning=False)
    assert check_result('mapping_jac.results')

#============================================================================
# Test runner for jacobian tests on quadratic elements
@pytest.mark.coverage
@pytest.mark.parametrize("grid_type", (1, 2))
@pytest.mark.parametrize("elem_order", (2, 3, 4))
def test_jac_quad(grid_type, elem_order):
    assert mapping_setup(1,1,2,91,'jac',grid_type=grid_type,elem_order=elem_order,catch_warning=False)
    assert check_result('mapping_jac.results')

#============================================================================
# Input file for 2D jacobian tests (circle or packed/stretched square triangle mesh)
oft_in_2d_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='test'
 cad_type={0}
 nlevels={1}
 nbase={1}
 grid_order=1
/

&sphere_options
 mesh_type=1
/

&cube_options
 mesh_type=1
 ni=2,3,0
 rscale=1.5,1.0,1.0
 packing=1.3,0.8,1.0
/

&test_mapping_options
 order={2}
/
"""

#============================================================================
# Test runner for 2D (triangle) jacobian tests on linear elements
# (the 2D Hessian does not support curved triangles or quads)
@pytest.mark.coverage
@pytest.mark.parametrize("cad_type", (91, 92))
@pytest.mark.parametrize("elem_order", (2, 3, 4))
def test_jac_2d_lin(cad_type, elem_order):
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_2d_template.format(cad_type, 2, elem_order))
    assert run_OFT("./test_mapping_jac_2d", 1, 60)
    assert check_result('mapping_jac.results')
