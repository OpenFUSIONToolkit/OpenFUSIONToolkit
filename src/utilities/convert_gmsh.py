from __future__ import print_function
import argparse
import os
import numpy as np
import h5py

ed_map_tri = np.array([
    [0,1],
    [1,2],
    [2,0]
])

ed_map_tet = np.array([
    [0,1],
    [1,2],
    [0,2],
    [0,3],
    [1,3],
    [2,3]
])


def read_legacy(fid):
    fid.readline() # Should be "Dimension"
    mesh_dim = int(fid.readline())
    # Read in vertices
    fid.readline() # Should be "Vertices"
    mesh_np = int(fid.readline())
    r = np.zeros((mesh_np,mesh_dim))
    for i in range(mesh_np):
        line_split = fid.readline().split()
        r[i,:] = [float(val) for val in line_split[:mesh_dim]]
    # Read in edges
    edge_type = fid.readline().strip() # Should be "Edges" or "EdgesP2"
    if edge_type == "Edges":
        nep = 2
        nfp = 3
        ncp = 4
    elif edge_type == "EdgesP2":
        nep = 3
        nfp = 6
        ncp = 10
    else:
        raise ValueError("Invalid edge element type")
    mesh_ne = int(fid.readline())
    le = np.zeros((mesh_np,nep), dtype=np.int32)
    for i in range(mesh_ne):
        line_split = fid.readline().split()
        le[i,:] = [int(val) for val in line_split[:nep]]
    # Read in faces
    fid.readline() # Should be "Triangles" or "TrianglesP2"
    mesh_nf = int(fid.readline())
    lf = np.zeros((mesh_nf,nfp), dtype=np.int32)
    lf_reg = np.zeros((mesh_nf,), dtype=np.int32)
    for i in range(mesh_nf):
        line_split = fid.readline().split()
        lf[i,:] = [int(val) for val in line_split[:nfp]]
        lf_reg[i] = int(line_split[nfp])
    # Read in cells
    line = fid.readline() # Should be "Tetrahedra" or "TetrahedraP2" or "End"
    if line.strip() == "End": # No tets, this is a surface mesh
        ncp_lin = 3
        mesh_nc = mesh_nf
        lc = lf
        lc_reg = lf_reg
    else:
        ncp_lin = 4
        mesh_nc = int(fid.readline())
        lc = np.zeros((mesh_nc,ncp), dtype=np.int32)
        lc_reg = np.zeros((mesh_nc,), dtype=np.int32)
        for i in range(mesh_nc):
            line_split = fid.readline().split()
            lc[i,:] = [int(val) for val in line_split[:ncp]]
            lc_reg[i] = int(line_split[ncp])
        fid.readline() # Should be "End"
    return r, lc, lc_reg, ncp_lin


def read_new(fid):
    def check_tag_line(tag):
        line = fid.readline().strip()
        tag = "$" + tag
        if line.lower() != tag.lower():
            raise ValueError('Expected line tag "{0}" not found "{1}"'.format(tag,line))
        return line
    fid.readline() # Mesh format numbers
    check_tag_line("EndMeshFormat")
    # Read in vertices
    check_tag_line("Nodes")
    mesh_np = int(fid.readline())
    mesh_dim = 3
    r = np.zeros((mesh_np,mesh_dim))
    for i in range(mesh_np):
        line_split = fid.readline().split()
        r[i,:] = [float(val) for val in line_split[1:mesh_dim+1]]
    check_tag_line("EndNodes")
    # Read in cells
    check_tag_line("Elements")
    mesh_nelems = int(fid.readline())
    lf = []
    lc = []
    for i in range(mesh_nelems):
        line_vals = [int(val) for val in fid.readline().split()]
        if line_vals[1] in (2,9):
            lf.append(line_vals[5:])
        elif line_vals[1] in (4,11):
            lc.append(line_vals[5:])
    check_tag_line("EndElements")
    if len(lc) == 0:
        lc = np.array(lf, dtype=np.int32)
        ncp_lin = 3
    else:
        lc = np.array(lc, dtype=np.int32)
        ncp_lin = 4
    lc_reg = np.ones((lc.shape[0],), dtype=np.int32)
    return r, lc, lc_reg, ncp_lin


def read_mesh(filename):
    print()
    print("Reading mesh: {0}".format(filename))
    with open(filename,'r') as fid:
        mesh_format_line = fid.readline()
        if mesh_format_line.strip() == '$MeshFormat':
            r, lc, lc_reg, ncp_lin = read_new(fid)
        else:
            r, lc, lc_reg, ncp_lin = read_legacy(fid)
    #
    mesh_np = r.shape[0]
    mesh_nc = lc.shape[0]
    nReg = lc_reg.max()
    mesh_order = 1
    if lc.shape[1] != ncp_lin:
        mesh_order = 2
    if ncp_lin == 3:
        mesh_dim = 2
        ed_map = ed_map_tri
    elif ncp_lin == 4:
        mesh_dim = 3
        ed_map = ed_map_tet
    # Reindex points and remove unreferenced points
    if mesh_order > 1: # Handle high-order if present
        reindex_flag_ho = np.zeros((mesh_np+1,), dtype=np.int32)
        reindex_flag_ho[lc[:,ncp_lin:].flatten()] = 1
        r_ho = r[reindex_flag_ho[1:] == 1].copy()
        rindexed_pts_ho = np.cumsum(reindex_flag_ho)
    else:
        r_ho = np.zeros((0,3))
    # Handle linear elements
    reindex_flag = np.zeros((mesh_np+1,), dtype=np.int32)
    reindex_flag[lc[:,:ncp_lin].flatten()] = 1
    r_new = r[reindex_flag[1:] == 1]
    rindexed_pts = np.cumsum(reindex_flag)
    lc_new = rindexed_pts[lc[:,:ncp_lin]]
    # Build high-order information
    if mesh_order > 1: # Handle high-order if present
        le_ho = np.zeros((r_ho.shape[0],2), dtype=np.int32)
        nce = ed_map.shape[0]
        for i in range(lc.shape[0]):
            for j in range(nce):
                le_ho[rindexed_pts_ho[lc[i,ncp_lin+j]]-1,:] = lc_new[i,ed_map[j,:]]
        ho_info = (r_ho, le_ho)
    else:
        ho_info = None
    print("  Dimension    = {0}".format(mesh_dim))
    print("  Mesh order   = {0}".format(mesh_order))
    print("  # of points  = {0} ({1})".format(r_new.shape[0],r_ho.shape[0]))
    print("  # of cells   = {0}".format(mesh_nc))
    print("  # of regions = {0}".format(nReg))
    #
    return r_new, lc_new, lc_reg, ho_info

def write_file(filename, r, lc, reg, ho_info=None):
    print()
    print("Saving mesh: {0}".format(filename))
    h5_file = h5py.File(filename, 'w')
    h5_file.create_dataset('mesh/R', data=r, dtype='f8')
    h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
    h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
    if ho_info is not None:
        h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
        h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')


parser = argparse.ArgumentParser()
parser.description = "Pre-processing script for Gmsh mesh files"
parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
options = parser.parse_args()

out_file = options.out_file
if out_file is None:
    out_file = os.path.splitext(options.in_file)[0] + ".h5"

r, lc, reg, ho_info = read_mesh(options.in_file)
write_file(out_file, r, lc, reg, ho_info)