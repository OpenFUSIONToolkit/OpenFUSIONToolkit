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


def write_native_mesh(filename, mesh_type, r, lc, reg, nodesets=None, sidesets=None, ho_info=None, periodic_info=None, reg_attrs=None, reg_names=None):
    r'''Create a native HDF5 mesh file for OFT from the given mesh information

    @param filename Filename for mesh file
    @param mesh_type Mesh type string ( "tri", "quad", "tet", "hex")
    @param r Points list [np,3]
    @param lc Cell list [nc,np_per_cell] (0-based)
    @param reg Region list [nc]
    @param nodesets List of node sets (0-based)
    @param sidesets List of side sets (0-based)
    @param ho_info High-order grid information
    @param periodic_info Information for mesh periodicity
    @param reg_attrs List of region attributes
    @param block_names List of block names
    '''
    print()
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5_file:
        # Write out basic mesh information
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        dset = h5_file.create_dataset('mesh/LC', data=lc+1, dtype='i4') # Convert to 1-based indexing for storage
        dset.attrs["TYPE"] = mesh_type.lower().encode('ascii')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
        # Write out high-order mesh information (nodes and indexing information)
        if ho_info is not None:
            h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
            h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1]+1, dtype='i4') # Convert to 1-based indexing for storage
            if ho_info[2] is not None:
                h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2]+1, dtype='i4') # Convert to 1-based indexing for storage
        # Write region names
        if reg_names is not None:
            max_len = max(map(len, reg_names))
            h5_file.create_dataset('mesh/REG_NAMES', data=numpy.array(reg_names, dtype=f"S{max_len}"))
        # Write region attributes
        if (reg_attrs is not None) and (len(reg_attrs) > 0):
            h5_file.create_dataset('mesh/reg_attrs/NUM_ATTR', data=[len(reg_attrs),], dtype='i4')
            for i, reg_attr in enumerate(reg_attrs):
                h5_file.create_dataset('mesh/reg_attrs/ATTR{0:04d}'.format(i+1), data=reg_attr, dtype='f8')
        # Write nodesets
        if (nodesets is not None) and (len(nodesets) > 0):
            h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
            for i, node_set in enumerate(nodesets):
                h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set+1, dtype='i4') # Convert to 1-based indexing for storage
        # Write sidesets (2D entity blocks)
        if (sidesets is not None) and (len(sidesets) > 0):
            h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(sidesets),], dtype='i4')
            for i, side_set in enumerate(sidesets):
                h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set+1, dtype='i4') # Convert to 1-based indexing for storage
        # Write flag for periodic nodes following mesh reflection
        if periodic_info is not None:
            h5_file.create_dataset('mesh/periodicity/NODES', data=periodic_info+1, dtype='i4') # Convert to 1-based indexing for storage


def read_native_mesh(filename):
    r'''Read mesh information from a native OFT HDF5 mesh file

    @param filename Filename for mesh file

    @returns Dictionary containing the mesh information
    '''
    with h5py.File(filename, 'r') as h5_file:
        # Read basic mesh information
        mesh = h5_file['mesh']
        lc = mesh['LC'][()] - 1 # Convert to 0-based indexing
        if 'TYPE' in mesh['LC'].attrs:
            mesh_type = mesh['LC'].attrs['TYPE']
            if isinstance(mesh_type, bytes):
                mesh_type = mesh_type.decode('ascii')
        else:
            if lc.shape[1] == 2:
                mesh_type = 'line'
            elif lc.shape[1] == 3:
                mesh_type = 'tri'
            elif lc.shape[1] == 4:
                print('Warning: Mesh type not specified in file, assuming tetrahedral mesh with 4-node cells')
                mesh_type = 'tet'
            elif lc.shape[1] == 8:
                mesh_type = 'hex'
        mesh_info = {
            'type': mesh_type,
            'r': mesh['R'][()],
            'lc': mesh['LC'][()] - 1, # Convert to 0-based indexing
            'reg': mesh['REG'][()],
        }
        # Read high-order mesh information (nodes and indexing information)
        if 'ho_info' in mesh:
            ho_group = mesh['ho_info']
            ho_info = [ho_group['R'][()], ho_group['LE'][()] - 1, None] # Convert to 0-based indexing
            if 'LF' in ho_group:
                ho_info[2] = ho_group['LF'][()] - 1 # Convert to 0-based indexing
            mesh_info['ho_info'] = tuple(ho_info)
        # Read region names
        if 'REG_NAMES' in mesh:
            mesh_info['reg_names'] = [name.decode('utf-8') if isinstance(name, bytes) else str(name)
                                        for name in mesh['REG_NAMES'][()]]
        # Read region attributes
        if 'reg_attrs' in mesh:
            reg_attrs = mesh['reg_attrs']
            if 'NUM_ATTR' in reg_attrs:
                mesh_info['reg_attrs'] = [reg_attrs['ATTR{0:04d}'.format(i+1)][()]
                                            for i in range(reg_attrs['NUM_ATTR'][0])]
                # Drop if all attributes are empty (h5py.Empty), to avoid unnecessary data in the output
                if all([isinstance(val, h5py.Empty) for val in mesh_info['reg_attrs']]):
                    del mesh_info['reg_attrs']
        # Read nodesets
        if 'NUM_NODESETS' in mesh:
            mesh_info['nodesets'] = [mesh['NODESET{0:04d}'.format(i+1)][()] - 1 # Convert to 0-based indexing
                                     for i in range(mesh['NUM_NODESETS'][0])]
        # Read sidesets (2D entity blocks)
        if 'NUM_SIDESETS' in mesh:
            mesh_info['sidesets'] = [mesh['SIDESET{0:04d}'.format(i+1)][()] - 1 # Convert to 0-based indexing
                                     for i in range(mesh['NUM_SIDESETS'][0])]
        # Read flag for periodic nodes following mesh reflection
        if 'periodicity' in mesh:
            periodicity = mesh['periodicity']
            if 'NODES' in periodicity:
                mesh_info['periodic_info'] = periodicity['NODES'][()] - 1 # Convert to 0-based indexing
    return mesh_info


def convert_mesh_to_pyvista(mesh_type, r, lc):
    '''! Convert mesh to pyvista representation

    @param mesh_type Mesh integer type (10=line, 21=tri, 23=quad, 31=tet, 33=hex)
    @param r Points list [np,3]
    @param lc Cell list [nc,3]

    @returns `pyvista.UnstructuredGrid` object for grid
    '''
    import pyvista
    if isinstance(mesh_type, str):
        mesh_type = mesh_type.lower()
        if mesh_type == 'line':
            mesh_type = 10
        elif mesh_type == 'tri':
            mesh_type = 21
        elif mesh_type == 'quad':
            mesh_type = 23
        elif mesh_type == 'tet':
            mesh_type = 31
        elif mesh_type == 'hex':
            mesh_type = 33
        else:
            raise ValueError("Unsupported mesh type: {0}".format(mesh_type))
    #
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
