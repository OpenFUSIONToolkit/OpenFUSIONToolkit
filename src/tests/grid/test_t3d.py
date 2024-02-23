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
 cad_type=1
 nlevels={1}
 nbase={0}
 grid_order={3}
/

&t3d_options
 filename='{2}.t3d'
 inpname='{2}.inp'
 reflect='xy'
 zstretch={4}
/
"""

# Common setup function and process handling
def t3d_setup(nbase, nlevels, prefix, order=1, zstretch=1.):
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, prefix, order, zstretch))
    return run_OFT("./test_t3d", nproc, 60)

# Validate results against expected values
def check_result(volume_test, area_test):
    retval = True
    fid = open('t3d.results','r')
    volume = float(fid.readline())
    area = float(fid.readline())
    fid.close()
    if abs(volume-volume_test) > 1.E-4:
        print("FAILED: Incorrect domain volume!")
        print("  Expected = {0}".format(volume_test))
        print("  Actual =   {0}".format(volume))
        retval = False
    if abs(area-area_test) > 1.E-4:
        print("FAILED: Incorrect domain surface area!")
        print("  Expected = {0}".format(area_test))
        print("  Actual =   {0}".format(area))
        retval = False
    return retval

#============================================================================
# Test runners for Cylindrical mesh
@pytest.mark.parametrize("top_lev", (1, 2))
def test_base(top_lev):
    volume_t3d = 3.1266
    area_t3d = 12.5187
    assert t3d_setup(1,top_lev,'cyl')
    assert check_result(volume_t3d, area_t3d)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_quad(top_lev):
    volume_t3d = 3.1416
    area_t3d = 12.5663
    assert t3d_setup(1,top_lev,'cyl',2)
    assert check_result(volume_t3d, area_t3d)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_1ref(top_lev):
    volume_t3d = 3.1378
    area_t3d = 12.5544
    minlev = 4 - top_lev
    assert t3d_setup(minlev,top_lev,'cyl')
    assert check_result(volume_t3d, area_t3d)

#============================================================================
# Test runners for stretched Cylindrical mesh
@pytest.mark.parametrize("top_lev", (1, 2))
def test_stretch_base(top_lev):
    volume_t3d = 2*3.1266
    area_t3d =    18.7916
    assert t3d_setup(1,top_lev,'cyl',zstretch=2.)
    assert check_result(volume_t3d, area_t3d)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_stretch_quad(top_lev):
    volume_t3d = 2*3.1416
    area_t3d =    18.8495
    assert t3d_setup(1,top_lev,'cyl',2,zstretch=2.)
    assert check_result(volume_t3d, area_t3d)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_stretch_1ref(top_lev):
    volume_t3d = 2*3.1378
    area_t3d =    18.8350
    minlev = 4 - top_lev
    assert t3d_setup(minlev,top_lev,'cyl',zstretch=2.)
    assert check_result(volume_t3d, area_t3d)

