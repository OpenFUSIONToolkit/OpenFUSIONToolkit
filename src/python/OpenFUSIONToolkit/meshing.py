#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! General OFT utilities for mesh generation and manipulation

@authors Chris Hansen
@date July 2026
@ingroup doxy_oft_python
'''
import numpy
import h5py
from .util import oft_warning


def write_native_mesh(filename, mesh_type, r, lc, reg, nodesets=[], sidesets=[], ho_info=None, periodic_info=None, block_attrs=None, block_names=None):
    r'''Create a native HDF5 mesh file for OFT from the given mesh information

    @param filename Filename for mesh file
    @param mesh_type Mesh type string ( "tri", "quad", "tet", "hex")
    @param r Points list [np,3]
    @param lc Cell list [nc,3] (1-based)
    @param reg Region list [nc]
    @param nodesets List of node sets
    @param sidesets List of side sets
    @param ho_info High-order grid information
    @param periodic_info Information for mesh periodicity
    @param block_attrs List of block attributes
    @param block_names List of block names
    '''
    print()
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5_file:
        # Write out basic mesh information
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        dset = h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
        dset.attrs["TYPE"] = mesh_type.encode('ascii')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
        # Write out high-order mesh information (nodes and indexing information)
        if ho_info is not None:
            h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
            h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')
            if ho_info[2] is not None:
                h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2], dtype='i4')
        # Write block names
        if block_names is not None:
            max_len = max(map(len, block_names))
            h5_file.create_dataset('mesh/reg_attr/BLOCK_NAMES', data=numpy.array(block_names, dtype=f"S{max_len}"))
        # Write block attributes
        if len(block_attrs) > 0:
            h5_file.create_dataset('mesh/reg_attr/NUM_ATTR', data=[len(block_attrs),], dtype='i4')
            for i, block_attr in enumerate(block_attrs):
                h5_file.create_dataset('mesh/reg_attr/ATTR{0:04d}'.format(i+1), data=block_attr, dtype='f8')
        # Write nodesets
        if len(nodesets) > 0:
            h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
            for i, node_set in enumerate(nodesets):
                h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
        # Write sidesets (2D entity blocks)
        if len(sidesets) > 0:
            h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(sidesets),], dtype='i4')
            for i, side_set in enumerate(sidesets):
                h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')
        # Write flag for periodic nodes following mesh reflection
        if periodic_info is not None:
            h5_file.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')


def read_native_mesh(filename):
    r'''Read mesh information from a native OFT HDF5 mesh file

    @param filename Filename for mesh file

    @returns Dictionary containing the mesh information
    '''
    with h5py.File(filename, 'r') as h5_file:
        # Read basic mesh information
        mesh = h5_file['mesh']
        mesh_type = mesh['LC'].attrs['TYPE']
        if isinstance(mesh_type, bytes):
            mesh_type = mesh_type.decode('ascii')
        mesh_info = {
            'mesh_type': mesh_type,
            'r': mesh['R'][()],
            'lc': mesh['LC'][()],
            'reg': mesh['REG'][()],
        }
        # Read high-order mesh information (nodes and indexing information)
        if 'ho_info' in mesh:
            ho_group = mesh['ho_info']
            ho_info = [ho_group['R'][()], ho_group['LE'][()], None]
            if 'LF' in ho_group:
                ho_info[2] = ho_group['LF'][()]
            mesh_info['ho_info'] = tuple(ho_info)
        if 'reg_attr' in mesh:
            reg_attr = mesh['reg_attr']
            # Read block names
            if 'BLOCK_NAMES' in reg_attr:
                mesh_info['block_names'] = [name.decode('utf-8') if isinstance(name, bytes) else str(name)
                                            for name in reg_attr['BLOCK_NAMES'][()]]
            # Read block attributes
            if 'NUM_ATTR' in reg_attr:
                mesh_info['block_attrs'] = [reg_attr['ATTR{0:04d}'.format(i+1)][()]
                                            for i in range(reg_attr['NUM_ATTR'][0])]
        # Read nodesets
        if 'NUM_NODESETS' in mesh:
            mesh_info['nodesets'] = [mesh['NODESET{0:04d}'.format(i+1)][()]
                                     for i in range(mesh['NUM_NODESETS'][0])]
        # Read sidesets (2D entity blocks)
        if 'NUM_SIDESETS' in mesh:
            mesh_info['sidesets'] = [mesh['SIDESET{0:04d}'.format(i+1)][()]
                                     for i in range(mesh['NUM_SIDESETS'][0])]
        # Read flag for periodic nodes following mesh reflection
        if 'periodicity' in mesh:
            periodicity = mesh['periodicity']
            if 'nodes' in periodicity:
                mesh_info['periodic_info'] = periodicity['nodes'][()]
    return mesh_info


def convert_mesh_to_pyvista(mesh_type, r, lc):
    '''! Convert mesh to pyvista representation

    @param mesh_type Mesh integer type (10=line, 21=tri, 23=quad, 31=tet, 33=hex)
    @param r Points list [np,3]
    @param lc Cell list [nc,3]

    @returns `pyvista.UnstructuredGrid` object for grid
    '''
    import pyvista
    if mesh_type == 31:
        celltype = pyvista.CellType.TETRA
        ncv = 4
    # elif mesh_type == 32:
    #     celltype = pyvista.CellType.QUADRATIC_TETRA
    #     ncv = 10
    elif mesh_type == 33:
        celltype = pyvista.CellType.HEXAHEDRON
        ncv = 8
    elif mesh_type == 21:
        celltype = pyvista.CellType.TRIANGLE
        ncv = 3
    # elif mesh_type == 22:
    #     celltype = pyvista.CellType.QUADRATIC_TRIANGLE
    #     ncv = 6
    elif mesh_type == 23:
        celltype = pyvista.CellType.QUAD
        ncv = 4
    elif mesh_type == 10:
        celltype = pyvista.CellType.LINE
        ncv = 2
    else:
        raise ValueError("Unsupported mesh type: {0}".format(mesh_type))
    celltypes = numpy.array([celltype for _ in range(lc.shape[0])], dtype=numpy.int8)
    cells = numpy.insert(lc, [0,], ncv, axis=1)
    return pyvista.UnstructuredGrid(cells, celltypes, r)
