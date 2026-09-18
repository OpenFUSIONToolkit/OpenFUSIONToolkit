#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python cli for GMSH to native Open FUSION Toolkit mesh conversion

See @ref convert_gmsh

@authors Chris Hansen
@date July 2026
@ingroup doxy_oft_python
'''
import os
import numpy as np
from .meshing import write_native_mesh

ed_map_tri = np.array([
    [0,1],
    [1,2],
    [2,0]
])

## Quadratic tetrahedron edge ordering used by the MEDIT (`.mesh`) format
ed_map_tet_medit = np.array([
    [0,1],
    [1,2],
    [0,2],
    [0,3],
    [1,3],
    [2,3]
])

## Quadratic tetrahedron edge ordering used by the GMSH (`.msh`) formats
ed_map_tet_msh = np.array([
    [0,1],
    [1,2],
    [0,2],
    [0,3],
    [2,3],
    [1,3]
])


def check_tag_line(fid, tag):
    '''! Read the next line and check that it matches the expected GMSH section tag

    @param fid Open file handle
    @param tag Expected section tag (without the leading "$")
    '''
    line = fid.readline().strip()
    tag = "$" + tag
    if line.lower() != tag.lower():
        raise ValueError('Expected line tag "{0}" not found "{1}"'.format(tag,line))
    return line


def read_node_block(fid, end_tag):
    '''! Read a GMSH node block, whose section tag has already been consumed

    @param fid Open file handle
    @param end_tag Section tag terminating the node block (without the leading "$")
    @result Point list [np,3]
    '''
    mesh_np = int(fid.readline())
    mesh_dim = 3
    r = np.zeros((mesh_np,mesh_dim))
    for i in range(mesh_np):
        line_split = fid.readline().split()
        if int(line_split[0]) != i+1:
            raise ValueError('Unexpected point index {0} (expected {1}), contiguous numbering is required'.format(int(line_split[0]),i+1))
        r[i,:] = [float(val) for val in line_split[1:mesh_dim+1]]
    check_tag_line(fid, end_tag)
    return r


def read_element_block(fid, end_tag):
    '''! Read a GMSH element block, whose section tag has already been consumed

    Triangles are only used as cells if no tetrahedra are present (surface mesh).

    @param fid Open file handle
    @param end_tag Section tag terminating the element block (without the leading "$")
    @result Cell list, region list, and number of linear points per cell
    '''
    mesh_nelems = int(fid.readline())
    lf = []
    lc = []
    for i in range(mesh_nelems):
        line_vals = [int(val) for val in fid.readline().split()]
        if line_vals[1] in (2,9):
            lf.append(line_vals[5:])
        elif line_vals[1] in (4,11):
            lc.append(line_vals[5:])
    check_tag_line(fid, end_tag)
    if len(lc) == 0:
        lc = np.array(lf, dtype=np.int32)
        ncp_lin = 3
    else:
        lc = np.array(lc, dtype=np.int32)
        ncp_lin = 4
    reg = np.ones((lc.shape[0],), dtype=np.int32)
    return lc, reg, ncp_lin


def read_medit(fid):
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
        reg = lf_reg
    else:
        ncp_lin = 4
        mesh_nc = int(fid.readline())
        lc = np.zeros((mesh_nc,ncp), dtype=np.int32)
        reg = np.zeros((mesh_nc,), dtype=np.int32)
        for i in range(mesh_nc):
            line_split = fid.readline().split()
            lc[i,:] = [int(val) for val in line_split[:ncp]]
            reg[i] = int(line_split[ncp])
        fid.readline() # Should be "End"
    return r, lc, reg, ncp_lin


def read_msh1(fid):
    '''! Read a GMSH mesh in the version 1 (`$NOD`/`$ELM`) format

    The leading "$NOD" tag is consumed by @ref read_mesh before this is called.

    @param fid Open file handle
    '''
    # Read in vertices
    r = read_node_block(fid, "ENDNOD")
    # Read in cells
    check_tag_line(fid, "ELM")
    lc, reg, ncp_lin = read_element_block(fid, "ENDELM")
    return r, lc, reg, ncp_lin


def read_msh2(fid):
    '''! Read a GMSH mesh in the version 2 (`$Nodes`/`$Elements`) format

    The leading "$MeshFormat" tag is consumed by @ref read_mesh before this is called.

    @param fid Open file handle
    '''
    fid.readline() # Mesh format numbers
    check_tag_line(fid, "EndMeshFormat")
    # Read in vertices
    check_tag_line(fid, "Nodes")
    r = read_node_block(fid, "EndNodes")
    # Read in cells
    check_tag_line(fid, "Elements")
    lc, reg, ncp_lin = read_element_block(fid, "EndElements")
    return r, lc, reg, ncp_lin


def read_mesh(filename):
    print()
    print("Reading mesh: {0}".format(filename))
    with open(filename,'r') as fid:
        mesh_format_line = fid.readline().strip()
        if mesh_format_line == '$MeshFormat': # GMSH version 2 format
            r, lc, reg, ncp_lin = read_msh2(fid)
            ed_map_tet = ed_map_tet_msh
        elif mesh_format_line == '$NOD': # GMSH version 1 format
            r, lc, reg, ncp_lin = read_msh1(fid)
            ed_map_tet = ed_map_tet_msh
        else: # MEDIT format
            r, lc, reg, ncp_lin = read_medit(fid)
            ed_map_tet = ed_map_tet_medit
    #
    mesh_np = r.shape[0]
    mesh_nc = lc.shape[0]
    nReg = reg.max()
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
    lc_new = rindexed_pts[lc[:,:ncp_lin]] - 1 # Convert back to 0-based indexing
    # Build high-order information
    if mesh_order > 1: # Handle high-order if present
        le_ho = np.zeros((r_ho.shape[0],2), dtype=np.int32)
        nce = ed_map.shape[0]
        for i in range(lc.shape[0]):
            for j in range(nce):
                le_ho[rindexed_pts_ho[lc[i,ncp_lin+j]]-1,:] = lc_new[i,ed_map[j,:]]
        ho_info = (r_ho, le_ho, None)
    else:
        ho_info = None
    mesh_type = 'TRI' if mesh_dim == 2 else 'TET'
    if mesh_order == 1:
        mesh_type += '_P1'
    else:
        mesh_type += '_P2'
    print("  Mesh type: {0}".format(mesh_type))
    print("  Dimension: {0}D".format(mesh_dim))
    print("  # of points  = {0} ({1})".format(r_new.shape[0],r_ho.shape[0]))
    print("  # of cells   = {0}".format(mesh_nc))
    print("  # of regions = {0}".format(nReg))
    #
    return mesh_type, r_new, lc_new, reg, ho_info


def convert_gmsh_to_native(in_file, out_file=None):
    '''! Convert GMSH mesh to native Open FUSION Toolkit mesh

    @param in_file Input mesh file
    @param out_file Output mesh file (optional, default is `in_file` with `.h5` extension)
    '''
    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + ".h5"

    # Read input GMSH file
    mesh_type, r, lc, reg, ho_info = read_mesh(in_file)

    # Write output file
    write_native_mesh(out_file, mesh_type.split('_')[0], r, lc, reg, ho_info=ho_info)


def script_entry():
    '''! Command line interface for GMSH to native Open FUSION Toolkit mesh conversion

    See @ref @convert_gmsh
    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.description = "Convert a GMSH mesh file to native Open FUSION Toolkit mesh format"
    parser.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    parser.add_argument("--out_file", type=str, default=None, help="Ouput mesh file")
    options = parser.parse_args()

    convert_gmsh_to_native(options.in_file, options.out_file)


## @page convert_gmsh `OFT_convert_gmsh`: GMSH to native Open FUSION Toolkit mesh conversion
#
# @section convert_gmsh_desc Description and options
# This script converts a GMSH mesh file to native Open FUSION Toolkit mesh format.
#
#```shell
# usage: OFT_convert_gmsh [-h] --in_file IN_FILE [--out_file OUT_FILE]
#
# Convert a GMSH mesh file to native Open FUSION Toolkit mesh format
#
# options:
#   -h, --help           show this help message and exit
#   --in_file IN_FILE    Input mesh file
#   --out_file OUT_FILE  Ouput mesh file
#```