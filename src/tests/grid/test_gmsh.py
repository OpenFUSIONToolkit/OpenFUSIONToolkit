import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT, run_command

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cylinder'
 cad_type=0
 nlevels={1}
 nbase={0}
 grid_order={2}
/

&native_mesh_options
 filename='{3}.h5'
/
"""

# Common setup function and process handling
def gmsh_setup(nbase, nlevels, mesh_name, grid_order=1):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, grid_order, mesh_name))
    # Run mesh conversion script
    convert_cmd = ["OFT_convert_gmsh.py", "--in_file={0}.msh".format(mesh_name)]
    outs, errs, errcode = run_command(" ".join(convert_cmd))
    if errcode != 0:
        print("FAILED: OFT_convert_gmsh.py exited with non-zero error code!")
        print("========== STD OUTPUT ==========")
        print(outs.decode())
        print("========== ERR OUTPUT ==========")
        print(errs.decode())
        print("ERRCODE = {0}".format(errcode))
        print("========== END OUTPUT ==========")
        return False
    return run_OFT("./test_gmsh", nproc, 60)

#
def check_result(volume_test, area_test):
    retval = True
    fid = open('gmsh.results','r')
    volume = float(fid.readline())
    if volume != volume_test:
        print("FAILED: Incorrect domain volume!")
        print("  Expected = {0}".format(volume_test))
        print("  Actual =   {0}".format(volume))
        retval = False
    area = float(fid.readline())
    if area != area_test:
        print("FAILED: Incorrect domain surface area!")
        print("  Expected = {0}".format(area_test))
        print("  Actual =   {0}".format(area))
        retval = False
    return retval

#============================================================================
# Test runners for basic Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_base(top_lev):
    volume_gmsh = pytest.approx(3.079621, abs=1.E-4)
    area_gmsh = pytest.approx(12.376435, abs=1.E-4)
    assert gmsh_setup(1,top_lev,'cyl_gmsh')
    assert check_result(volume_gmsh, area_gmsh)

#============================================================================
# Test runner for quadratic Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_quad(top_lev):
    volume_gmsh = pytest.approx(3.141495, abs=1.E-4)
    area_gmsh = pytest.approx(12.565964, abs=1.E-4)
    assert gmsh_setup(1,top_lev,'cyl_gmsh',grid_order=2)
    assert check_result(volume_gmsh, area_gmsh)

#============================================================================
# Test runner for single refinement Cylinder mesh
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_1ref(top_lev):
    volume_gmsh = pytest.approx(3.126017, abs=1.E-4)
    area_gmsh = pytest.approx(12.518663, abs=1.E-4)
    minlev = 4 - top_lev
    assert gmsh_setup(minlev,top_lev,'cyl_gmsh')
    assert check_result(volume_gmsh, area_gmsh)
