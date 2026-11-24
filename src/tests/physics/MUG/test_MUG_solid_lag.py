#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: solid (resistive wall) regions in xMHD with all fields on Lagrange elements. A unit box is split
# into a plasma region (x < 1/2) and a solid region (x > 1/2, XML <region> with type=2), and a velocity
# perturbation in the plasma is advanced with each velocity boundary condition. The velocity must remain exactly
# zero at nodes in the solid region, and the solve (block-Jacobi/LU preconditioner) must converge, which requires
# the velocity, density, and temperature to be fixed at all nodes in the solid region.
import os
import sys
import itertools
import numpy
import h5py
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='solid'
 cad_type=0
 nlevels={nlevels}
 nbase=1
 grid_order=1
/

&native_mesh_options
 filename='solid_box.h5'
/

&test_solid_options
 order=2
 b0=0.1
 n0=1.E19
 t0=10.
 v0=1.E2
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=T
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_upwind=F
 bbc='bc'
 vbc='{vbc}'
 nbc='{nbc}'
 tbc='{tbc}'
 dt=2.E-7
 eta=1.E-2
 nu_visc=1.E1
 d_dens=1.E1
 kappa_par=1.E1,1.E1
 kappa_perp=1.E1,1.E1
 nsteps=10
 rst_freq=1000
 lin_tol=1.E-9
 nl_tol=1.E-6
 maxextrap=2
 ittarget=40
/
"""

# Maximum velocity allowed at nodes in the solid region (relative to the perturbation amplitude)
solid_tol = 1.E-12

NCELLS = 4


def write_box_mesh(filename):
    '''Write a tetrahedral mesh of the unit cube with region 2 for x > 1/2 (each hex split into 6 Kuhn tets)'''
    def node_id(i, j, k):
        return (i*(NCELLS+1) + j)*(NCELLS+1) + k
    r = numpy.array([[i, j, k] for i in range(NCELLS+1) for j in range(NCELLS+1) for k in range(NCELLS+1)],
                    dtype=float)/NCELLS
    cells = []
    for i, j, k in itertools.product(range(NCELLS), repeat=3):
        for perm in itertools.permutations(range(3)):
            corner = [i, j, k]
            tet = [node_id(*corner)]
            for axis in perm:
                corner[axis] += 1
                tet.append(node_id(*corner))
            cells.append(tet)
    lc = numpy.array(cells)
    edges = [r[lc[:, m]]-r[lc[:, 0]] for m in (1, 2, 3)]
    vol = numpy.einsum('ij,ij->i', numpy.cross(edges[0], edges[1]), edges[2])
    lc[vol < 0] = lc[vol < 0][:, [0, 1, 3, 2]]
    reg = numpy.where(r[lc].mean(axis=1)[:, 0] > 0.5, 2, 1)
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')


def solid_case(vbc, nbc, tbc, nlevels=1, nproc=1):
    '''Write the input and mesh for the solid-region case and run test_MUG_solid_lag

    Quadratic elements on the base mesh (`nlevels=1`), or on one refinement of it (`nlevels=2`, needed for MPI).
    '''
    os.chdir(test_dir)
    write_box_mesh('solid_box.h5')
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(vbc=vbc, nbc=nbc, tbc=tbc, nlevels=nlevels))
    if os.path.exists('solid.results'):
        os.remove('solid.results')
    return run_OFT("./test_MUG_solid_lag oft.in oft_in_xmhd_solid.xml", nproc, 300)


def validate_result():
    '''Check that the velocity is zero in the solid region and that the plasma flow is still present'''
    with open('solid.results', 'r') as fid:
        vmax_solid = float(fid.readline())
        vrms = float(fid.readline())
    print("max |V| (solid)/v0 = {0:.4E}, RMS V/v0 = {1:.4E}".format(vmax_solid, vrms))
    retval = True
    if not vmax_solid < solid_tol:
        print("FAILED: Velocity in solid region! ({0})".format(vmax_solid))
        retval = False
    if not (1.E-2 < vrms < 1.0):
        print("FAILED: Unexpected plasma flow amplitude! ({0})".format(vrms))
        retval = False
    return retval


@pytest.mark.parametrize("vbc,nbc,tbc", (('all', 'd', 'd'), ('norm', 'n', 'n'), ('tang', 'n', 'd'), ('none', 'd', 'n')))
def test_solid(vbc, nbc, tbc):
    '''Velocity perturbation next to a solid region with the given velocity, density, and temperature boundary
    conditions.'''
    assert solid_case(vbc, nbc, tbc)
    assert validate_result()


@pytest.mark.mpi
def test_solid_mpi():
    '''Solid-region case with vbc='norm' and Neumann density/temperature BCs distributed over two MPI processes.'''
    assert solid_case('norm', 'n', 'n', nlevels=2, nproc=2)
    assert validate_result()
