#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python cli for Cubit to native Open FUSION Toolkit mesh conversion

@authors Chris Hansen
@date July 2026
@ingroup doxy_oft_python
'''
import os
import numpy as np
from .meshing import write_native_mesh

tri_ed_map = np.array([
    [0,1],
    [1,2],
    [2,0]
])
tet_ed_map = np.array([
    [0,1],
    [1,2],
    [0,2],
    [0,3],
    [1,3],
    [2,3]
])
quad_ed_map = np.array([
    [0,1],
    [1,2],
    [2,3],
    [3,0]
])
hex_ed_map = np.array([
    [0,1],
    [1,2],
    [2,3],
    [3,0],
    [0,4],
    [1,5],
    [2,6],
    [3,7],
    [4,5],
    [5,6],
    [6,7],
    [4,7]
])
hex_face_map = np.array([
    [0,1,2,3],
    [4,5,6,7],
    [0,3,7,4],
    [1,2,6,5],
    [0,1,5,4],
    [2,3,7,6]
])
element_type_map = {
    'TRI3': {'dim': 2, 'ncp': 3, 'ncp_lin': 3, 'type': 'TRI_p1', 'ed_map': tri_ed_map},
    'TRI6': {'dim': 2, 'ncp': 6, 'ncp_lin': 3, 'type': 'TRI_p2', 'ed_map': tri_ed_map},
    'QUAD4': {'dim': 2, 'ncp': 4, 'ncp_lin': 4, 'type': 'QUAD_p1', 'ed_map': quad_ed_map},
    'QUAD9': {'dim': 2, 'ncp': 9, 'ncp_lin': 4, 'type': 'QUAD_p2', 'ed_map': quad_ed_map, 'cellpoint': True},
    'TETRA4': {'dim': 3, 'ncp': 4, 'ncp_lin': 4, 'type': 'TET_p1', 'ed_map': tet_ed_map},
    'TETRA10': {'dim': 3, 'ncp': 10, 'ncp_lin': 4, 'type': 'TET_p2', 'ed_map': tet_ed_map},
    'HEX8': {'dim': 3, 'ncp': 8, 'ncp_lin': 8, 'type': 'HEX_p1', 'ed_map': hex_ed_map, 'face_map': hex_face_map},
    'HEX27': {'dim': 3, 'ncp': 27, 'ncp_lin': 8, 'type': 'HEX_p2', 'ed_map': hex_ed_map, 'face_map': hex_face_map, 'cellpoint': True}
}
element_type_map['TRI'] = element_type_map['TRI3']
element_type_map['TETRA'] = element_type_map['TETRA4']
element_type_map['QUAD'] = element_type_map['QUAD4']
element_type_map['HEX'] = element_type_map['HEX8']

def read_mesh(filename, ignore_attrs):
    import netCDF4
    print()
    print("Reading mesh: {0}".format(filename))
    ncdf_file = netCDF4.Dataset(filename, "r")
    # Read points
    if 'coord' in ncdf_file.variables:
        r = np.asarray(ncdf_file.variables['coord']).transpose().copy()
    else:
        r = []
        r.append(np.asarray(ncdf_file.variables['coordx']))
        r.append(np.asarray(ncdf_file.variables['coordy']))
        if 'coordz' in ncdf_file.variables:
            rz = np.asarray(ncdf_file.variables['coordz'])
            if np.max(abs(rz)) > 0.0:
                r.append(rz)
        r = np.transpose(np.vstack(r)).copy()
    # Read regions
    lc = []
    node_sets = []
    side_sets = []
    block_types = []
    block_attrs = []
    max_logical_dim = 0
    for varname, variable in ncdf_file.variables.items():
        if varname.startswith('connect'):
            for attrname in variable.ncattrs():
                if attrname.startswith('elem_type'):
                    block_type = getattr(variable, attrname)
                    block_type_info = element_type_map.get(block_type,{'type': block_type, 'dim': -1})
                    max_logical_dim = max(max_logical_dim,block_type_info['dim'])
                    block_types.append(block_type_info)
            lc_tmp = np.asarray(variable)
            lc.append(lc_tmp)
            # Retrieve attributes if present
            attr_name = 'attrib' + varname[7:]
            if attr_name in ncdf_file.variables:
                if np.min(ncdf_file.variables[attr_name],axis=0) != np.max(ncdf_file.variables[attr_name],axis=0):
                    raise ValueError("Attributes must have the same value within a block.")
                block_attrs.append(np.asarray(ncdf_file.variables[attr_name][0,:]))
            else:
                block_attrs.append(None)
        elif varname.startswith('node_ns'):
            node_sets.append(np.asarray(variable)-1) # Convert to 0-based indexing
        elif varname.startswith('elem_ss'):
            side_sets.append(np.asarray(variable)-1) # Convert to 0-based indexing
    # Read block names if present
    if 'eb_names' in ncdf_file.variables:
        block_names = ["".join(block_name.compressed().astype(str)) for block_name in ncdf_file.variables['eb_names']]
    else:
        block_names = None
    # Remove lower level geometry
    keep_inds = []
    reg = []
    for i, block_type in enumerate(block_types):
        if block_type['dim'] < max_logical_dim:
            print("  Note: Removing block {0} of type {1}".format(i+1,block_type['type']))
            continue
        keep_inds.append(i)
        nReg = len(keep_inds)
        reg.append(nReg*np.ones((lc[i].shape[0],), dtype=np.int32))
    lc = [lc[i] for i in keep_inds]
    block_types = [block_types[i] for i in keep_inds]
    block_attrs = [block_attrs[i] for i in keep_inds]
    block_names = [block_names[i] for i in keep_inds] if block_names is not None else None
    lc = np.vstack(lc)
    reg = np.hstack(reg)
    mesh_order = 1
    ncp_lin = lc.shape[1]
    mesh_type = block_types[0]['type']
    for block_type in block_types:
        if block_type['type'] != mesh_type:
            raise ValueError("Found blocks of different type {0} {1}".format(mesh_type,block_type['type']))
    if mesh_type[-1] == '2':
        mesh_order = 2
        ncp_lin = block_types[0]['ncp_lin']
    np_orig = r.shape[0]
    # Handle linear elements
    reindex_flag = np.zeros((r.shape[0]+1,), dtype=np.int32)
    reindex_flag[lc[:,:ncp_lin].flatten()] = 1
    r_new = r[reindex_flag[1:] == 1]
    for i, nodeset in enumerate(node_sets):
        node_sets[i] = np.array([node for node in nodeset if reindex_flag[node] == 1])
    rindexed_pts = np.cumsum(reindex_flag)
    lc_new = rindexed_pts[lc[:,:ncp_lin]] - 1 # Convert back to 0-based indexing
    print('LC CHK',lc_new.min(axis=None))
    node_sets = [rindexed_pts[nodeset] for nodeset in node_sets]
    # Build high-order information
    if mesh_order > 1: # Handle high-order if present
        ed_map = block_types[0]['ed_map']
        nce = ed_map.shape[0]
        reindex_flag_ho = np.zeros((r.shape[0]+1,), dtype=np.int32)
        reindex_flag_ho[lc[:,ncp_lin:ncp_lin+nce].flatten()] = 1
        r_ho = r[reindex_flag_ho[1:] == 1].copy()
        rindexed_pts_ho = np.cumsum(reindex_flag_ho)
        ne = np.max(rindexed_pts_ho)
        le_ho = np.zeros((ne,2), dtype=np.int32)
        for i in range(lc.shape[0]):
            for j in range(nce):
                le_ho[rindexed_pts_ho[lc[i,ncp_lin+j]]-1,:] = lc_new[i,ed_map[j,:]]
        face_map = block_types[0].get('face_map',None)
        if face_map is not None:
            ncf = face_map.shape[0]
            reindex_flag_ho = np.zeros((r.shape[0]+1,), dtype=np.int32)
            reindex_flag_ho[lc[:,ncp_lin+1+nce:ncp_lin+1+nce+ncf].flatten()] = 1
            r_ho = np.vstack((r_ho,r[reindex_flag_ho[1:] == 1].copy()))
            rindexed_pts_ho = np.cumsum(reindex_flag_ho)
            nf = np.max(rindexed_pts_ho)
            lf_ho = np.zeros((nf,4), dtype=np.int32)
            for i in range(lc.shape[0]):
                for j in range(ncf):
                    lf_ho[rindexed_pts_ho[lc[i,ncp_lin+1+nce+j]]-1,:] = lc_new[i,face_map[j,:]]
        else:
            lf_ho = None
        if block_types[0].get('cellpoint',False):
            reindex_flag_ho = np.zeros((r.shape[0]+1,), dtype=np.int32)
            reindex_flag_ho[lc[:,ncp_lin+nce]] = 1
            r_ho = np.vstack((r_ho,r[reindex_flag_ho[1:] == 1].copy()))
        ho_info = (r_ho, le_ho, lf_ho)
    else:
        r_ho = np.zeros((0,r_new.shape[1]))
        ho_info = None
    # Check attributes
    if ignore_attrs:
        block_attrs = []
    else:
        nattrs = 0
        for block_attr in block_attrs:
            if block_attr is not None:
                nattrs = block_attr.shape[0]
        for block_attr in block_attrs:
            if block_attr is None:
                if nattrs > 0:
                    raise ValueError("Attributes specified for only some blocks.")
            elif nattrs != block_attr.shape[0]:
                raise ValueError("Attribute size must be the same on all blocks.")
        if nattrs > 0:
            block_attrs = [block_attr for block_attr in np.array(block_attrs).transpose()]
        else:
            block_attrs = []
    #
    print("  Mesh type: {0}".format(mesh_type))
    print("  Dimension: {0}D".format(r_new.shape[1]))
    print("  # of points       = {0} ({1})".format(r_new.shape[0],r_ho.shape[0]))
    print("  # of cells        = {0}".format(lc_new.shape[0]))
    print("  # of regions      = {0}".format(nReg))
    print("  # of nSets        = {0}".format(len(node_sets)))
    print("  # of sSets        = {0}".format(len(side_sets)))
    if len(block_attrs) > 0:
        print("  # of attributes   = {0}".format(len(block_attrs)))
    np_total = r_new.shape[0]+r_ho.shape[0]
    if np_total > np_orig:
        raise ValueError("One or more points referenced by both linear and high-order nodes.")
    if np_total < np_orig:
        print("""
Note: {0} points were not referenced by cells.
This may be normal or could indicate an error""".format(np_total-np_orig))
    #
    return mesh_type, r_new, lc_new, reg, node_sets, side_sets, ho_info, block_attrs, block_names


def convert_cubit_to_native(in_file, out_file=None, periodic_nodeset=None, ignore_attr=False):
    '''! Convert Cubit (exodus) mesh to native Open FUSION Toolkit mesh

    @param in_file Input mesh file
    @param out_file Output mesh file (optional, default is `in_file` with `.h5` extension)
    @param periodic_nodeset Index of periodic nodeset (optional)
    @param ignore_attr Ignore block attributes (optional)
    '''
    out_file = out_file
    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + ".h5"

    # Read input Exodus file
    mesh_type, r, lc, reg, node_sets, side_sets, ho_info, block_attrs, block_names = read_mesh(in_file, ignore_attr)

    # Map periodicity information
    periodic_info = None
    if options.periodic_nodeset is not None:
        if options.periodic_nodeset > len(node_sets):
            raise ValueError("Periodic nodeset ({0}) is out of bounds ({1})".format(options.periodic_nodeset, len(node_sets)))
        periodic_info = node_sets.pop(options.periodic_nodeset-1)

    # Write output file
    write_native_mesh(out_file, mesh_type.split('_')[0], r, lc, reg, nodesets=node_sets, sidesets=side_sets,
                    ho_info=ho_info, periodic_info=periodic_info, reg_attrs=block_attrs, reg_names=block_names)


def script_entry():
    '''! Command line interface for Cubit (exodus) to native Open FUSION Toolkit mesh conversion
    options:
      -h, --help            show this help message and exit
      --in_file IN_FILE     Input mesh file
      --out_file OUT_FILE   Ouput mesh file
      --periodic_nodeset PERIODIC_NODESET
                            Index of perioidic nodeset
      --ignore_attr         Ignore block attributes
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.description = "Convert a Cubit (exodus) mesh file to native Open FUSION Toolkit mesh format"
    parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
    parser.add_argument("--periodic_nodeset", type=int, default=None, help="Index of perioidic nodeset")
    parser.add_argument("--ignore_attr", default=False, action="store_true", help="Ignore block attributes")
    options = parser.parse_args()

    convert_cubit_to_native(options.in_file, out_file=options.out_file, periodic_nodeset=options.periodic_nodeset, ignore_attr=options.ignore_attr)
