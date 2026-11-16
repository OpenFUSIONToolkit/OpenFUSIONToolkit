#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Regression tests for periodic native meshes (planar and revolved) and vector
# Lagrange elements on them. Small structured tet meshes are generated here and
# written in the native HDF5 format, then reflected about z = 0 by OFT to form
# the periodic domain.
import os
import sys
import itertools
import numpy
import h5py
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


def structured_tets(nu, nv, nw, mapping):
    '''Tetrahedralize a structured grid on the unit cube and map it to physical space

    Each hex is split into the 6 Kuhn tets sharing its main diagonal, which gives a conforming
    mesh. Points that map to the same location (eg. a collapsed axis) are merged and the
    resulting degenerate tets are dropped. Returns points, 0-based cells, and the indices of
    points on the w = 1 face.
    '''
    def node_id(i, j, k):
        return (i*(nv+1) + j)*(nw+1) + k
    uvw = numpy.array([[i/nu, j/nv, k/nw] for i in range(nu+1) for j in range(nv+1) for k in range(nw+1)])
    r_full = numpy.array([mapping(*pt) for pt in uvw])
    cells = []
    for i, j, k in itertools.product(range(nu), range(nv), range(nw)):
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
    # Drop collapsed cells and orient the rest with positive volume
    lc = lc[[len(set(cell)) == 4 for cell in lc]]
    vol = numpy.einsum('ij,ij->i', numpy.cross(r[lc[:, 1]]-r[lc[:, 0]], r[lc[:, 2]]-r[lc[:, 0]]), r[lc[:, 3]]-r[lc[:, 0]])
    assert numpy.all(numpy.abs(vol) > 1.E-8)
    flip = vol < 0
    lc[numpy.ix_(flip, [2, 3])] = lc[numpy.ix_(flip, [3, 2])]
    per_nodes = numpy.unique(new_index[inverse[numpy.isclose(uvw[:, 2], 1.0)]])
    return r, lc, per_nodes


def x_wedge(u, v, w):
    '''Wedge about the x-axis, including the axis, with the w = 0 face on z = 0'''
    theta = WEDGE_ANGLE*w
    return [u, v*numpy.cos(theta), v*numpy.sin(theta)]


def y_wedge(u, v, w):
    '''Wedge about the y-axis, with the w = 0 face on z = 0'''
    theta = WEDGE_ANGLE*w
    return [v*numpy.cos(theta), u, v*numpy.sin(theta)]


def z_slab(u, v, w):
    '''Unit square slab extending from z = 0 to z = 0.5'''
    return [u, v, 0.5*w]


def write_mesh(filename, mapping):
    r, lc, per_nodes = structured_tets(NU, NV, NW, mapping)
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=numpy.ones((lc.shape[0],)), dtype='i4')
        h5_file.create_dataset('mesh/periodicity/nodes', data=per_nodes+1, dtype='i4')


def periodic_setup(mesh_type, order=1, nbase=1, nlevels=1, nproc=1, expect_abort=False):
    mappings = {'x_wedge': x_wedge, 'y_wedge': y_wedge, 'z_slab': z_slab}
    os.chdir(test_dir)
    mesh_file = 'periodic_{0}.h5'.format(mesh_type)
    write_mesh(mesh_file, mappings[mesh_type])
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nbase, nlevels, mesh_file, order))
    if os.path.exists('periodic.results'):
        os.remove('periodic.results')
    if expect_abort:
        # oft_abort reports on stderr, which run_OFT does not return
        outs, errs, errcode = run_command("./test_vlag_periodic")
        print(outs.decode())
        print(errs.decode())
        return errcode, errs.decode()
    return run_OFT("./test_vlag_periodic", nproc, 60)


def check_results(revolved, nper):
    with open('periodic.results', 'r') as fid:
        revolved_test = int(fid.readline()) == 1
        nper_test = [int(val) for val in fid.readline().split()]
        err_good = float(fid.readline())
        err_bad = float(fid.readline())
    retval = True
    if revolved_test != revolved:
        print("FAILED: Incorrect periodicity type (revolved = {0})".format(revolved_test))
        retval = False
    if nper_test != nper:
        print("FAILED: Incorrect periodic (point, edge, face) counts")
        print("  Expected = {0}".format(nper))
        print("  Actual =   {0}".format(nper_test))
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


# Copied face (w = 1) of the generated mesh, after reflection: the axis points and the NU
# edges between them lie on the reflection plane, so they are not copied and are not
# periodic, but other edges and faces touching the axis are.
X_WEDGE_COUNTS = [(NU+1)*NV, NU*(NV+1) + (NU+1)*NV + NU*NV - NU, 2*NU*NV]
Z_SLAB_COUNTS = [(NU+1)*(NV+1), NU*(NV+1) + (NU+1)*NV + NU*NV, 2*NU*NV]


@pytest.mark.coverage
@pytest.mark.parametrize("order", (1, 2))
def test_revolved(order):
    '''Revolved periodicity about the x-axis (native_set_periodic and Lagrange vector rotation)

    A wedge containing the x-axis is reflected to give two periodic faces related by
    rotation. Checks that the domain is flagged as revolved, the periodic point, edge, and
    face counts, and that a rotation-compatible linear vector field is reproduced by
    projection onto vector Lagrange elements, which requires rotating vectors across the
    periodic boundary, while a non-compatible field is not.
    '''
    assert periodic_setup('x_wedge', order=order)
    assert check_results(True, X_WEDGE_COUNTS)


def test_revolved_refined():
    '''Revolved periodicity is carried through uniform refinement (same checks as `test_revolved`)'''
    assert periodic_setup('x_wedge', nbase=2, nlevels=2)
    assert check_results(True, X_WEDGE_COUNTS)


@pytest.mark.mpi
def test_revolved_mpi():
    '''Revolved periodicity on a distributed mesh (same checks as `test_revolved`)'''
    assert periodic_setup('x_wedge', nbase=1, nlevels=2, nproc=2)
    assert check_results(True, X_WEDGE_COUNTS)


@pytest.mark.coverage
def test_planar():
    '''Planar periodicity in z (native_set_periodic)

    A slab is reflected about z = 0 to give periodic faces normal to z. Checks that the
    domain is not flagged as revolved, the periodic point, edge, and face counts, and that a
    z-independent linear vector field is reproduced by projection, while a z-dependent field
    is not.
    '''
    assert periodic_setup('z_slab')
    assert check_results(False, Z_SLAB_COUNTS)


@pytest.mark.coverage
def test_revolved_bad_axis():
    '''Revolved periodicity about the y-axis is rejected with an error by native_set_periodic'''
    errcode, std_err = periodic_setup('y_wedge', expect_abort=True)
    assert errcode != 0
    assert std_err.find('Revolved periodic planes must contain the x-axis') > -1
    assert not os.path.exists('periodic.results')
