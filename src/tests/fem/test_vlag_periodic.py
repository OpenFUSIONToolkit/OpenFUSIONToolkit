#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for periodic native meshes (planar and revolved) and vector
# Lagrange elements on them. Small structured tet and hex meshes are generated here
# and written in the native HDF5 format, then reflected about z = 0 by OFT to form
# the periodic domain.
import os
import sys
import itertools
import numpy
import h5py
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
 meshname='periodic'
 cad_type=0
 nlevels={1}
 nbase={0}
 grid_order=1
/

&native_mesh_options
 filename='{2}'
 reflect=T
 ref_periodic=T
/

&test_lag_options
 order={3}
/
"""

# Grid resolution: nu along the axis direction, nv radially (or across), nw around (or along z)
NU, NV, NW = 4, 4, 2
# Half-angle of the generated wedge; the reflected (periodic) domain spans twice this
WEDGE_ANGLE = numpy.pi/6.0
# Inner radius of the wedges: tet wedges include the axis, hexes cannot collapse onto it
R_MIN = {'tet': 0.0, 'hex': 0.5}
# Hex corner ordering in OFT (inodesp_base in hexmesh_type.F90) as (u,v,w) offsets
HEX_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def structured_mesh(nu, nv, nw, mapping, cell_type):
    '''Mesh a structured grid on the unit cube with tets or hexes and map it to physical space

    For tets, each hex is split into the 6 Kuhn tets sharing its main diagonal, which gives a
    conforming mesh. Points that map to the same location (eg. a collapsed axis) are merged and
    the resulting degenerate cells are dropped. Cells are oriented to have positive volume.
    Returns points, 0-based cells, and the indices of points on the w = 1 face.
    '''
    def node_id(i, j, k):
        return (i*(nv+1) + j)*(nw+1) + k
    uvw = numpy.array([[i/nu, j/nv, k/nw] for i in range(nu+1) for j in range(nv+1) for k in range(nw+1)])
    r_full = numpy.array([mapping(*pt) for pt in uvw])
    cells = []
    for i, j, k in itertools.product(range(nu), range(nv), range(nw)):
        if cell_type == 'hex':
            cells.append([node_id(i+di, j+dj, k+dk) for di, dj, dk in HEX_CORNERS])
            continue
        for perm in itertools.permutations(range(3)):
            corner = [i, j, k]
            tet = [node_id(*corner)]
            for axis in perm:
                corner[axis] += 1
                tet.append(node_id(*corner))
            cells.append(tet)
    # Merge coincident points
    _, unique_ind, inverse = numpy.unique(numpy.round(r_full, 10), axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    order = numpy.argsort(unique_ind)
    new_index = numpy.empty_like(order)
    new_index[order] = numpy.arange(order.size)
    r = r_full[unique_ind[order]]
    lc = new_index[inverse[numpy.array(cells)]]
    lc = lc[[len(set(cell)) == lc.shape[1] for cell in lc]]
    # Orientation from the corner spanned by the first edges along u, v, and w
    if cell_type == 'hex':
        edges = [r[lc[:, 1]]-r[lc[:, 0]], r[lc[:, 3]]-r[lc[:, 0]], r[lc[:, 4]]-r[lc[:, 0]]]
        mirror = [3, 2, 1, 0, 7, 6, 5, 4]  # Reverse the v direction
    else:
        edges = [r[lc[:, 1]]-r[lc[:, 0]], r[lc[:, 2]]-r[lc[:, 0]], r[lc[:, 3]]-r[lc[:, 0]]]
        mirror = [0, 1, 3, 2]
    vol = numpy.einsum('ij,ij->i', numpy.cross(edges[0], edges[1]), edges[2])
    assert numpy.all(numpy.abs(vol) > 1.E-8)
    lc[vol < 0] = lc[vol < 0][:, mirror]
    per_nodes = numpy.unique(new_index[inverse[numpy.isclose(uvw[:, 2], 1.0)]])
    return r, lc, per_nodes


def x_wedge(u, v, w, r_min):
    '''Wedge about the x-axis, with the w = 0 face on z = 0'''
    theta = WEDGE_ANGLE*w
    r = r_min + (1.0-r_min)*v
    return [u, r*numpy.cos(theta), r*numpy.sin(theta)]


def y_wedge(u, v, w, r_min):
    '''Wedge about the y-axis, with the w = 0 face on z = 0'''
    theta = WEDGE_ANGLE*w
    r = r_min + (1.0-r_min)*v
    return [r*numpy.cos(theta), u, r*numpy.sin(theta)]


def z_slab(u, v, w, r_min):
    '''Unit square slab extending from z = 0 to z = 0.5'''
    return [u, v, 0.5*w]


def write_mesh(filename, mapping, cell_type):
    r, lc, per_nodes = structured_mesh(NU, NV, NW, lambda u, v, w: mapping(u, v, w, R_MIN[cell_type]), cell_type)
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=numpy.ones((lc.shape[0],)), dtype='i4')
        h5_file.create_dataset('mesh/periodicity/nodes', data=per_nodes+1, dtype='i4')


def periodic_setup(mesh_type, cell_type, order=1, nbase=1, nlevels=1, nproc=1, expect_abort=None):
    mappings = {'x_wedge': x_wedge, 'y_wedge': y_wedge, 'z_slab': z_slab}
    os.chdir(test_dir)
    mesh_file = 'periodic_{0}_{1}.h5'.format(mesh_type, cell_type)
    write_mesh(mesh_file, mappings[mesh_type], cell_type)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, mesh_file, order))
    if os.path.exists('periodic.results'):
        os.remove('periodic.results')
    return run_OFT("./test_vlag_periodic", nproc, 60, expect_abort=expect_abort)


def refine_counts(counts, cell_type):
    '''Global (point, edge, face, cell) counts after one uniform refinement (each cell split into 8)'''
    np, ne, nf, nc = counts
    if cell_type == 'hex':
        return [np + ne + nf + nc, 2*ne + 4*nf + 6*nc, 4*nf + 12*nc, 8*nc]
    return [np + ne, 2*ne + 3*nf + nc, 4*nf + 8*nc, 8*nc]


def check_results(revolved, nper, cell_type, nbase=1, nlevels=1):
    with open('periodic.results', 'r') as fid:
        revolved_test = [int(val) == 1 for val in fid.readline().split()]
        nper_test = [int(val) for val in fid.readline().split()]
        local_counts = [int(val) for val in fid.readline().split()]
        base_counts = [int(val) for val in fid.readline().split()]
        fine_counts = [int(val) for val in fid.readline().split()]
        err_good = float(fid.readline())
        err_bad = float(fid.readline())
    retval = True
    if revolved_test != [revolved, revolved]:
        print("FAILED: Incorrect periodicity type (revolved on base, finest level = {0})".format(revolved_test))
        retval = False
    if nper_test != nper:
        print("FAILED: Incorrect periodic (point, edge, face) counts")
        print("  Expected = {0}".format(nper))
        print("  Actual =   {0}".format(nper_test))
        retval = False
    # Periodic entities must be merged with their parents in the base level's global indexing
    if [local_counts[i] - base_counts[i] for i in range(3)] != nper_test:
        print("FAILED: Periodic entities not merged in base level global indexing")
        print("  Local (point, edge, face) counts =  {0}".format(local_counts))
        print("  Global (point, edge, face) counts = {0}".format(base_counts[:3]))
        retval = False
    # Refinement and distribution must keep periodic entities merged on the finest level. Levels up to
    # nbase are refined, level nbase+1 is the distributed copy of level nbase, and later levels are refined
    nrefine = (nbase - 1) + max(nlevels - nbase - 1, 0)
    fine_expected = base_counts
    for _ in range(nrefine):
        fine_expected = refine_counts(fine_expected, cell_type)
    if fine_counts != fine_expected:
        print("FAILED: Incorrect global (point, edge, face, cell) counts on the finest level")
        print("  Expected = {0}".format(fine_expected))
        print("  Actual =   {0}".format(fine_counts))
        retval = False
    # A compatible linear field is in the FE space, so only solver error remains
    if err_good > 1.E-8:
        print("FAILED: Compatible field not reproduced, error = {0}".format(err_good))
        retval = False
    # An incompatible field must not be reproduced if periodic DOFs are shared
    if err_bad < 1.E-2:
        print("FAILED: Incompatible field reproduced, periodic DOFs not shared? error = {0}".format(err_bad))
        retval = False
    return retval


# Periodic (point, edge, face) counts: entities on the copied (w = 1) face of the generated
# mesh after reflection. For the tet wedge the axis points and the NU edges between them lie
# on the reflection plane, so they are not copied and are not periodic, but other edges and
# faces touching the axis are.
EXPECTED_COUNTS = {
    ('x_wedge', 'tet'): [(NU+1)*NV, NU*(NV+1) + (NU+1)*NV + NU*NV - NU, 2*NU*NV],
    ('z_slab', 'tet'): [(NU+1)*(NV+1), NU*(NV+1) + (NU+1)*NV + NU*NV, 2*NU*NV],
    ('x_wedge', 'hex'): [(NU+1)*(NV+1), NU*(NV+1) + (NU+1)*NV, NU*NV],
    ('z_slab', 'hex'): [(NU+1)*(NV+1), NU*(NV+1) + (NU+1)*NV, NU*NV]
}


@pytest.mark.coverage
@pytest.mark.parametrize("cell_type", ('tet', 'hex'))
@pytest.mark.parametrize("order", (1, 2))
def test_revolved(order, cell_type):
    '''Revolved periodicity about the x-axis (native_set_periodic and Lagrange vector rotation)

    A wedge about the x-axis (including the axis for tets) is reflected to give two periodic
    faces related by rotation. Checks that the domain is flagged as revolved, the periodic
    point, edge, and face counts, and that a rotation-compatible linear vector field is
    reproduced by projection onto vector Lagrange elements, which requires rotating vectors
    across the periodic boundary, while a non-compatible field is not.
    '''
    assert periodic_setup('x_wedge', cell_type, order=order)
    assert check_results(True, EXPECTED_COUNTS[('x_wedge', cell_type)], cell_type)


@pytest.mark.parametrize("cell_type", ('tet', 'hex'))
def test_revolved_refined(cell_type):
    '''Revolved periodicity is carried through uniform refinement

    Same checks as `test_revolved`, with the projection done on the refined level. The periodic counts
    are for the base level; on the refined level, periodic entities are only visible through the global
    indexing, so its global counts are checked against those expected from refining the base level.
    '''
    assert periodic_setup('x_wedge', cell_type, nbase=2, nlevels=2)
    assert check_results(True, EXPECTED_COUNTS[('x_wedge', cell_type)], cell_type, nbase=2, nlevels=2)


@pytest.mark.mpi
@pytest.mark.parametrize("cell_type", ('tet', 'hex'))
def test_revolved_mpi(cell_type):
    '''Revolved periodicity on a distributed mesh, refined after distribution

    Same checks as `test_revolved_refined`, with the projection done on the distributed, refined level.
    '''
    assert periodic_setup('x_wedge', cell_type, nbase=1, nlevels=3, nproc=2)
    assert check_results(True, EXPECTED_COUNTS[('x_wedge', cell_type)], cell_type, nbase=1, nlevels=3)


@pytest.mark.coverage
@pytest.mark.parametrize("cell_type", ('tet', 'hex'))
def test_planar(cell_type):
    '''Planar periodicity in z (native_set_periodic)

    A slab is reflected about z = 0 to give periodic faces normal to z. Checks that the
    domain is not flagged as revolved, the periodic point, edge, and face counts, and that a
    z-independent linear vector field is reproduced by projection, while a z-dependent field
    is not.
    '''
    assert periodic_setup('z_slab', cell_type)
    assert check_results(False, EXPECTED_COUNTS[('z_slab', cell_type)], cell_type)


@pytest.mark.coverage
@pytest.mark.parametrize("cell_type", ('tet', 'hex'))
def test_revolved_bad_axis(cell_type):
    '''Revolved periodicity about the y-axis is rejected with an error by native_set_periodic'''
    assert periodic_setup('y_wedge', cell_type, expect_abort='Revolved periodic planes must contain the x-axis')
    assert not os.path.exists('periodic.results')
