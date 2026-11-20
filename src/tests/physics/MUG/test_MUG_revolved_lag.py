#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: Z-pinch equilibrium on a rotationally periodic (revolved) wedge about the x-axis (xMHD with all fields
# on Lagrange elements). The azimuthal field and quadratic pressure are exactly represented, so the state must
# remain in equilibrium, which requires the vector fields (B and V) to be rotated correctly across the periodic
# boundaries. Run with the magnetic and velocity boundary conditions applied on the curved outer wall.
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..', 'fem')))
from oft_testing import run_OFT
from test_vlag_periodic import write_mesh, x_wedge

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='revolved'
 cad_type=0
 nlevels={nlevels}
 nbase={nbase}
 grid_order=1
/

&native_mesh_options
 filename='revolved_tet.h5'
 reflect=T
 ref_periodic=T
/

&test_revolved_options
 order=2
 minlev={minlev}
 c=0.1
 n0=1.E19
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=T
 xmhd_adv_den=F
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_upwind=F
 bbc='{bbc}'
 vbc='{vbc}'
 tbc='d'
 dt=2.E-7
 eta=1.E-3
 nu_visc=1.E-1
 kappa_par=1.E-6,1.E-6
 kappa_perp=1.E-6,1.E-6
 nsteps=25
 rst_freq=1000
 lin_tol=1.E-12
 nl_tol=1.E-8
 nu_xmhd={nu_mhd}
 maxextrap=2
 ittarget=100
/
"""

# Allowed RMS velocity (relative to the wall Alfven speed) and relative change in B, well above the round-off level
# (~1E-8) reached when the vector fields are rotated correctly across the periodic boundaries
vel_tol = 1.E-6
b_tol = 1.E-6


def revolved_case(bbc, vbc, nbase=1, nlevels=1, minlev=2, nu_mhd='-1', nproc=1):
    '''Write the input and mesh for the Z-pinch equilibrium and run test_MUG_revolved_lag

    The FE levels run from `minlev` to the mesh levels plus one (quadratic elements); the standard multigrid
    preconditioner (`oft_in_xmhd_mg.xml`: GMRES smoother with block-Jacobi/LU on each level) is used if there is
    more than one.
    '''
    os.chdir(test_dir)
    write_mesh('revolved_tet.h5', x_wedge, 'tet')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(bbc=bbc, vbc=vbc, nbase=nbase, nlevels=nlevels, minlev=minlev, nu_mhd=nu_mhd))
    if os.path.exists('revolved.results'):
        os.remove('revolved.results')
    xml_file = 'oft_in_xmhd.xml'
    if nlevels+1 > minlev:
        xml_file = 'oft_in_xmhd_mg.xml'
    return run_OFT("./test_MUG_revolved_lag oft.in {0}".format(xml_file), nproc, 600)


def validate_result():
    '''Check that the flow and the change in B remain at round-off level'''
    with open('revolved.results', 'r') as fid:
        vel_err = float(fid.readline())
        b_err = float(fid.readline())
    print("RMS V/v_A = {0:.4E}, |dB|/|B| = {1:.4E}".format(vel_err, b_err))
    retval = True
    if not vel_err < vel_tol:
        print("FAILED: Spurious flow generated! ({0})".format(vel_err))
        retval = False
    if not b_err < b_tol:
        print("FAILED: Magnetic field changed! ({0})".format(b_err))
        retval = False
    return retval


@pytest.mark.parametrize("bbc,vbc", (('bc', 'all'), ('bc', 'norm'), ('ic', 'tang')))
def test_pinch(bbc, vbc):
    '''Z-pinch equilibrium on a revolved tetrahedral wedge (quadratic elements, single level) with the given
    magnetic and velocity boundary conditions.'''
    assert revolved_case(bbc, vbc)
    assert validate_result()


def test_pinch_mg():
    '''Z-pinch equilibrium on a refined revolved tetrahedral wedge with a multigrid preconditioner, which uses
    the interpolation operators for vector fields on revolved meshes (geometric and polynomial levels).'''
    assert revolved_case('bc', 'all', nbase=2, nlevels=2, minlev=1, nu_mhd='2,2,2')
    assert validate_result()


@pytest.mark.mpi
def test_pinch_mpi():
    '''Z-pinch equilibrium on a revolved tetrahedral wedge distributed over two MPI processes.'''
    assert revolved_case('bc', 'all', nbase=1, nlevels=2, minlev=3, nproc=2)
    assert validate_result()
