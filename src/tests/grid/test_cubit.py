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
 meshname='sphere'
 cad_type={9}
 nlevels={1}
 nbase={0}
 grid_order={4}
/

&cubit_test_options
 test_surf={8}
 cad_type={9}
/

&cubit_options
 filename='{2}.g'
 {3}
 lf_file=T
 reflect={5}
 per_ns={6}
 zstretch={7}
/

&native_mesh_options
 filename='{2}.h5'
 reflect={5}
 ref_periodic={10}
/
"""

# Common setup function and process handling
def cubit_setup(nbase, nlevels, prefix, inp_prefix=None, grid_order=1,
                reflect='F', per_ns=-1, zstretch=1., test_2d='F', cad_type=2):
    if inp_prefix == None:
        prefix2 = "! inpname='none'"
    else:
        prefix2 = "inpname='{0}.3dm'".format(inp_prefix)
    nproc = 1
    if nbase != nlevels:
        nproc = 2
    ref_periodic = 'F'
    if per_ns > 0:
        ref_periodic = 'T'
    #
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, prefix, prefix2,
                                       grid_order, reflect, per_ns, zstretch, test_2d, cad_type, ref_periodic))
    return run_OFT("./test_cubit", nproc, 60)

# Validate results against expected values
def check_result(volume_test, area_test, nbp_test=None, nbe_test=None, tol=1.E-4):
    retval = True
    fid = open('cubit.results','r')
    volume = float(fid.readline())
    if abs(volume-volume_test) > tol:
        print("FAILED: Incorrect domain volume!")
        print("  Expected = {0}".format(volume_test))
        print("  Actual =   {0}".format(volume))
        retval = False
    area = float(fid.readline())
    if abs(area-area_test) > tol:
        print("FAILED: Incorrect domain surface area!")
        print("  Expected = {0}".format(area_test))
        print("  Actual =   {0}".format(area))
        retval = False
    if nbp_test is not None:
        nbp = int(fid.readline())
        if abs(nbp-nbp_test) > tol:
            print("FAILED: Incorrect number of boundary points!")
            print("  Expected = {0}".format(nbp_test))
            print("  Actual =   {0}".format(nbp))
            retval = False
    if nbe_test is not None:
        nbe = int(fid.readline())
        if abs(nbe-nbe_test) > tol:
            print("FAILED: Incorrect number of boundary edges!")
            print("  Expected = {0}".format(nbe_test))
            print("  Actual =   {0}".format(nbe))
            retval = False
    return retval

#============================================================================
# Test runners for basic Cubit mesh
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_base(top_lev,cad_type):
    volume_cubit = 3.9246
    area_cubit = 12.1320
    assert cubit_setup(1,top_lev,'sphere_tet4_test',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_quad(top_lev):
    volume_cubit = 4.18656
    area_cubit = 12.56197
    assert cubit_setup(1,top_lev,'sphere_tet4_test','sphere_test',grid_order=2)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_1ref(top_lev):
    volume_cubit = 4.11948
    area_cubit = 12.4519
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'sphere_tet4_test','sphere_test')
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_hex_base(top_lev,cad_type):
    volume_cubit = 3.91742
    area_cubit = 12.15673
    assert cubit_setup(1,top_lev,'sphere_hex8_test',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for high order input meshes
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_tet10_quad(top_lev,cad_type):
    volume_cubit = 4.18669
    area_cubit = 12.56224
    assert cubit_setup(1,top_lev,'sphere_tet10_test',grid_order=2,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
@pytest.mark.parametrize("cad_type", (0, 2))
@pytest.mark.parametrize("grid_order", (1, 2))
def test_tet10_1ref(top_lev,cad_type,grid_order):
    volume_cubit = [4.12061, 4.18669]
    area_cubit = [12.4552, 12.5622]
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'sphere_tet10_test',grid_order=grid_order,cad_type=cad_type)
    assert check_result(volume_cubit[grid_order-1], area_cubit[grid_order-1])
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_hex27_quad(top_lev,cad_type):
    volume_cubit = 4.18831
    area_cubit = 12.56542
    assert cubit_setup(1,top_lev,'sphere_hex27_test',grid_order=2,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
@pytest.mark.parametrize("cad_type", (0, 2))
@pytest.mark.parametrize("grid_order", (1, 2))
def test_hex27_1ref(top_lev,cad_type,grid_order):
    volume_cubit = [4.11911, 4.18831]
    area_cubit = [12.4621, 12.5654]
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'sphere_hex27_test',grid_order=grid_order,cad_type=cad_type)
    assert check_result(volume_cubit[grid_order-1], area_cubit[grid_order-1])

#============================================================================
# Test runners for cut meshes
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_cut_base(top_lev,cad_type):
    volume_cubit = 3.918291
    area_cubit = 12.12685
    assert cubit_setup(1,top_lev,'sphere_cut_test',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_cut_1ref(top_lev):
    volume_cubit = 4.12584
    area_cubit = 12.46292
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'sphere_cut_test','sphere_cut_test')
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for reflection
@pytest.mark.parametrize("top_lev", (1, 2))
def test_reflect_tri_base(top_lev,cad_type=0):
    volume_cubit = 0.0
    area_cubit =  2.0
    assert cubit_setup(1,top_lev,'ref_tri3_test',reflect='T',cad_type=cad_type,test_2d='T')
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.parametrize("top_lev", (1, 2))
def test_reflect_quad_base(top_lev,cad_type=0):
    volume_cubit = 0.0
    area_cubit =  2.0
    assert cubit_setup(1,top_lev,'ref_quad4_test',reflect='T',cad_type=cad_type,test_2d='T')
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_reflect_base(top_lev,cad_type):
    volume_cubit = 3.02044
    area_cubit =  12.26295
    assert cubit_setup(1,top_lev,'ref_tet4_test',reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_reflect_hex_base(top_lev,cad_type):
    volume_cubit = 3.00000
    area_cubit =  12.21166
    assert cubit_setup(1,top_lev,'ref_hex8_test',reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_reflect_1ref(top_lev):
    volume_cubit = 3.11110
    area_cubit =  12.49011
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'ref_tet4_test','ref_test',reflect='T')
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_reflect_quad(top_lev):
    volume_cubit = 3.14123
    area_cubit =  12.56531
    assert cubit_setup(1,top_lev,'ref_tet10_test','ref_test',grid_order=2,reflect='T')
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for periodic reflection
@pytest.mark.parametrize("top_lev", (1, 2))
def test_perreflect_tri_base(top_lev,cad_type=0):
    volume_cubit = 0.0
    area_cubit =   2.0
    assert cubit_setup(1,top_lev,'ref_tri3_test',reflect='T',per_ns=1,cad_type=cad_type,test_2d='T')
    assert check_result(volume_cubit, area_cubit, 80, 80)
@pytest.mark.parametrize("top_lev", (1, 2))
def test_perreflect_quad_base(top_lev,cad_type=0):
    volume_cubit = 0.0
    area_cubit =   2.0
    assert cubit_setup(1,top_lev,'ref_quad4_test',reflect='T',per_ns=1,cad_type=cad_type,test_2d='T')
    assert check_result(volume_cubit, area_cubit, 80, 80)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_perreflect_tet_base(top_lev,cad_type):
    volume_cubit = 3.02044
    area_cubit =   6.22207
    assert cubit_setup(1,top_lev,'ref_tet4_test',reflect='T',per_ns=1,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, 26, 78)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_perreflect_hex_base(top_lev,cad_type):
    volume_cubit = 3.00000
    area_cubit =   6.21166
    assert cubit_setup(1,top_lev,'ref_hex8_test',reflect='T',per_ns=1,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, 24, 48)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_perreflect_tet_1ref(top_lev):
    volume_cubit = 3.11110
    area_cubit =   6.26791
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'ref_tet4_test','ref_test',reflect='T',per_ns=1)
    assert check_result(volume_cubit, area_cubit, 68, 167)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_perreflect_tet_quad(top_lev):
    volume_cubit = 3.14123
    area_cubit =   6.28283
    assert cubit_setup(1,top_lev,'ref_tet10_test','ref_test',grid_order=2,reflect='T',per_ns=1)
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for high order reflected meshes
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_tet10reflect_quad(top_lev,cad_type):
    volume_cubit = 3.14123
    area_cubit =  12.56531
    assert cubit_setup(1,top_lev,'ref_tet10_test',grid_order=2,reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_tet10reflect_1ref(top_lev,cad_type):
    volume_cubit = 3.11110
    area_cubit =  12.48994
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'ref_tet10_test',reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_hex27reflect_quad(top_lev,cad_type):
    volume_cubit = 3.14110
    area_cubit =  12.56491
    assert cubit_setup(1,top_lev,'ref_hex27_test',grid_order=2,reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_hex27reflect_1ref(top_lev,cad_type):
    volume_cubit = 3.10583
    area_cubit =  12.47691
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'ref_hex27_test',reflect='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for stretched meshes
@pytest.mark.parametrize("top_lev", (1, 2))
def test_stretch_base(top_lev):
    volume_cubit = 3.02044
    area_cubit =  12.26295
    assert cubit_setup(1,top_lev,'ref_tet4_test',zstretch=2.)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (2, 3))
def test_stretch_1ref(top_lev):
    volume_cubit = 3.11110
    area_cubit =  12.49011
    minlev = 4 - top_lev
    assert cubit_setup(minlev,top_lev,'ref_tet4_test','ref_test',zstretch=2.)
    assert check_result(volume_cubit, area_cubit)
@pytest.mark.coverage
@pytest.mark.parametrize("top_lev", (1, 2))
def test_stretch_quad(top_lev):
    volume_cubit = 3.14123
    area_cubit =  12.56531
    assert cubit_setup(1,top_lev,'ref_tet10_test','ref_test',grid_order=2,zstretch=2.)
    assert check_result(volume_cubit, area_cubit)

#============================================================================
# Test runners for surface meshes
@pytest.mark.parametrize("mesh_type", ("tri3", "tri6", "quad4", "quad9"))
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_surf_circle_base(mesh_type,top_lev,cad_type):
    volume_cubit = 0.0
    area_cubit = 3.1403
    assert cubit_setup(1,top_lev,'circle_{0}_test'.format(mesh_type),test_2d='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)
@pytest.mark.coverage
@pytest.mark.parametrize("mesh_type", ("tri3", "tri6", "quad4", "quad9"))
@pytest.mark.parametrize("levels", ((2, 2), (2, 3), (1, 3)))
def test_surf_circle_1ref(mesh_type,levels):
    volume_cubit = 0.0
    area_cubit = 3.1403
    assert cubit_setup(levels[0],levels[1],'circle_{0}_test'.format(mesh_type),test_2d='T')
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)
@pytest.mark.coverage
@pytest.mark.parametrize("mesh_type", ("tri6", "quad9"))
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_surf_circle_quad(mesh_type,top_lev,cad_type):
    volume_cubit = 0.0
    area_cubit = 3.14160
    assert cubit_setup(1,top_lev,'circle_{0}_test'.format(mesh_type),test_2d='T',grid_order=2,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)
@pytest.mark.parametrize("mesh_type", ("tri3", "quad4"))
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_surf_sphere_base(mesh_type,top_lev,cad_type):
    volume_cubit = 0.0
    area_cubit = 12.5595
    assert cubit_setup(1,top_lev,'sphere_{0}_test'.format(mesh_type),test_2d='T',cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)
@pytest.mark.coverage
@pytest.mark.parametrize("mesh_type", ("tri3", "quad4"))
@pytest.mark.parametrize("levels", ((2, 2), (2, 3), (1, 3)))
def test_surf_sphere_1ref(mesh_type,levels):
    volume_cubit = 0.0
    area_cubit = 12.5595
    assert cubit_setup(levels[0],levels[1],'sphere_{0}_test'.format(mesh_type),test_2d='T')
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)
@pytest.mark.coverage
@pytest.mark.parametrize("mesh_type", ("tri6", "quad9"))
@pytest.mark.parametrize("top_lev", (1, 2))
@pytest.mark.parametrize("cad_type", (0, 2))
def test_surf_sphere_quad(mesh_type,top_lev,cad_type):
    volume_cubit = 0.0
    area_cubit = 12.5664
    assert cubit_setup(1,top_lev,'sphere_{0}_test'.format(mesh_type),test_2d='T',grid_order=2,cad_type=cad_type)
    assert check_result(volume_cubit, area_cubit, tol=1.E-3)