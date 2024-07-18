from __future__ import print_function
import argparse
import os
import numpy as np
import h5py
import netCDF4

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

def read_mesh(filename):
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
        elif varname.startswith('node_ns'):
            node_sets.append(np.asarray(variable))
        elif varname.startswith('elem_ss'):
            side_sets.append(np.asarray(variable))
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
    lc_new = rindexed_pts[lc[:,:ncp_lin]]
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
    print("  Mesh type    = {0}".format(mesh_type))
    print("  Dimension    = {0}".format(r_new.shape[1]))
    print("  # of points  = {0} ({1})".format(r_new.shape[0],r_ho.shape[0]))
    print("  # of cells   = {0}".format(lc_new.shape[0]))
    print("  # of regions = {0}".format(nReg))
    print("  # of nSets   = {0}".format(len(node_sets)))
    print("  # of sSets   = {0}".format(len(side_sets)))
    np_total = r_new.shape[0]+r_ho.shape[0]
    if np_total > np_orig:
        raise ValueError("One or more points referenced by both linear and high-order nodes.")
    if np_total < np_orig:
        print("""
Note: {0} points were not referenced by cells.
This may be normal or could indicate an error""".format(np_total-np_orig))
    #
    return r_new, lc_new, reg, node_sets, side_sets, ho_info

def write_file(filename, r, lc, reg, node_sets=[], side_sets=[], ho_info=None, periodic_info=None):
    print()
    print("Saving mesh: {0}".format(filename))
    h5_file = h5py.File(filename, 'w')
    h5_file.create_dataset('mesh/R', data=r, dtype='f8')
    h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
    h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
    if ho_info is not None:
        h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
        h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')
        if ho_info[2] is not None:
            h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2], dtype='i4')
    if len(node_sets) > 0:
        h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(node_sets),], dtype='i4')
        for i, node_set in enumerate(node_sets):
            h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
    if len(side_sets) > 0:
        h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(side_sets),], dtype='i4')
        for i, side_set in enumerate(side_sets):
            h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')
    if periodic_info is not None:
        h5_file.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')


parser = argparse.ArgumentParser()
parser.description = "Pre-processing script for mesh files"
parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
parser.add_argument("--periodic_nodeset", type=int, default=None, help="Index of perioidic nodeset")
options = parser.parse_args()

out_file = options.out_file
if out_file is None:
    out_file = os.path.splitext(options.in_file)[0] + ".h5"

r, lc, reg, node_sets, side_sets, ho_info = read_mesh(options.in_file)

periodic_info = None
if options.periodic_nodeset is not None:
    if options.periodic_nodeset > len(node_sets):
        raise ValueError("Periodic nodeset ({0}) is out of bounds ({1})".format(options.periodic_nodeset, len(node_sets)))
    periodic_info = node_sets.pop(options.periodic_nodeset-1)

write_file(out_file, r, lc, reg, node_sets, side_sets, ho_info, periodic_info)