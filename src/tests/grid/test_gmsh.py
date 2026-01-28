from __future__ import print_function
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
 meshname='cylinder'
 cad_type={3}
 nlevels={1}
 nbase={0}
 grid_order={2}
/

&gmsh_options
 filename='cyl.mesh'
 order=2
/

&native_mesh_options
 filename='cyl_gmsh.h5'
/
"""

# Common setup function and process handling
def gmsh_setup(nbase, nlevels, grid_order=1, cad_type=3):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, grid_order, cad_type))
    return run_OFT("./test_gmsh", nproc, 60)

#
def check_result(volume_test, area_test):
    retval = True
    fid = open('gmsh.results','r')
    volume = float(fid.readline())
    if abs(volume-volume_test) > 1.E-4:
        print("FAILED: Incorrect domain volume!")
        print("  Expected = {0}".format(volume_test))
        print("  Actual =   {0}".format(volume))
        retval = False
    area = float(fid.readline())
    if abs(area-area_test) > 1.E-4:
        print("FAILED: Incorrect domain surface area!")
        print("  Expected = {0}".format(area_test))
        print("  Actual =   {0}".format(area))
        retval = False
    return retval

#============================================================================
# Test runners for basic Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 3))
def test_base(top_lev,cad_type):
    volume_gmsh = 3.0913
    area_gmsh = 12.3812
    assert gmsh_setup(1,top_lev,cad_type=cad_type)
    assert check_result(volume_gmsh, area_gmsh)

#============================================================================
# Test runner for quadratic Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 3))
def test_quad(top_lev,cad_type):
    volume_gmsh = 3.1415
    area_gmsh = 12.5660
    assert gmsh_setup(1,top_lev,grid_order=2,cad_type=cad_type)
    assert check_result(volume_gmsh, area_gmsh)

#============================================================================
# Test runner for single refinement Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
@pytest.mark.parametrize("cad_type", (0, 3))
def test_1ref(top_lev,cad_type):
    volume_gmsh = 3.1290
    area_gmsh = 12.5198
    minlev = 4 - top_lev
    assert gmsh_setup(minlev,top_lev,cad_type=cad_type)
    assert check_result(volume_gmsh, area_gmsh)
