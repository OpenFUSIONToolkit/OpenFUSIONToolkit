#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#------------------------------------------------------------------------------
#
# Script for handling creation of arbitrary order tessellations of the unit
# Tetrahedron and unit Triangle. Used with plotting methods in the Open FUSION Toolkit (OFT) to map
# high order Lagrange nodes to a linear representation by subdividing elements.
#
#------------------------------------------------------------------------------
from __future__ import print_function
try:
    import numpy
except:
    print("================================")
    print("ERROR: NUMPY is required to run this script")
    print("================================")
    raise
try:
    import scipy.version
    if int(scipy.version.version.split('.')[0]) <= 0 and int(scipy.version.version.split('.')[1]) < 9:
        raise ImportError
except:
    print("================================")
    print("ERROR: SCIPY v0.9 or greater is required to run this script")
    print("================================")
    raise
else:
    import scipy.spatial
#----------------------------------------------------------------
# Format NUMPY array to FORTRAN syntax
#----------------------------------------------------------------
def tess_fort_array(a,mask,name='data',char_per_line=40):
    dims = a.shape
    nskip = mask.sum()
    tmp_string = '{0}({2},{1}) = RESHAPE((/'.format(name,dims[0]-nskip,dims[1])
    for i in range(dims[0]):
        if mask[i] == 1:
            continue
        for j in range(dims[1]):
            tmp_string = tmp_string + str(a[i,j]) + ','
    tmp_string = tmp_string[:-1] + '/),(/{1},{0}/))'.format(dims[0]-nskip,dims[1])
    #
    out_string = ''
    ccount = 0
    for i in range(len(tmp_string)):
        ccount = ccount + 1
        out_string = out_string + tmp_string[i]
        if ccount > char_per_line and tmp_string[i] == ',':
            out_string = out_string + '\n'
            ccount = 0
    return out_string
#----------------------------------------------------------------
# Check for zero volume cells
#----------------------------------------------------------------
def check_tets(triang):
    dims = triang.simplices.shape
    mask = numpy.zeros(dims[0],numpy.int)
    for (i,cell) in enumerate(triang.simplices):
        e1 = triang.points[cell[3],:] - triang.points[cell[0],:]
        e2 = triang.points[cell[3],:] - triang.points[cell[1],:]
        e3 = triang.points[cell[3],:] - triang.points[cell[2],:]
        vol = abs(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1]) + e1[1]*(e2[2]*e3[0]-e2[0]*e3[2]) + e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))
        if vol < 1.e-10:
            mask[i]=1
    return mask
#----------------------------------------------------------------
# Triangle tesselation method
#----------------------------------------------------------------
def tesselate_tri(order=1,ed_dofs=[],fc_dofs=[]):
    pts_base = numpy.array([[0., 0.], [1., 0.], [.5, 1.]], numpy.float)
    ed_base = numpy.array([[2, 1], [0, 2], [1, 0]], numpy.int)
    #
    ne_dofs = len(ed_dofs)
    nf_dofs = len(fc_dofs)
    #
    npts = 3 + ne_dofs*3 + nf_dofs
    pts_out = numpy.zeros((npts,2), numpy.float)
    #
    for i in range(3):
        pts_out[i,:] = pts_base[i,:]
    #
    for i in range(ne_dofs):
        node = ed_dofs[i]
        for j in range(3):
            pt1 = pts_base[ed_base[j,0],:]*node[0]
            pt2 = pts_base[ed_base[j,1],:]*node[1]
            pts_out[j+i*3+3,:] = pt1 + pt2
    #
    for i in range(nf_dofs):
        node = fc_dofs[i]
        pt1 = pts_base[0,:]*node[0]
        pt2 = pts_base[1,:]*node[1]
        pt3 = pts_base[2,:]*node[2]
        pts_out[i+3*ne_dofs+3,:] = pt1 + pt2 + pt3
    #
    triang = scipy.spatial.Delaunay(pts_out)
    dims = triang.simplices.shape
    mask = numpy.zeros(dims[0],numpy.int)
    return triang, mask
#----------------------------------------------------------------
# Tetrahedron tesselation method
#----------------------------------------------------------------
def tesselate_tet(order=1,ed_dofs=[],fc_dofs=[],c_dofs=[]):
    pts_base = numpy.array([[0., 0., 0.], [1., 0., 0.], [.5, 1., 0.], [.5, .5, 1.]], numpy.float)
    ed_base = numpy.array([[0, 3], [1, 3], [2, 3], [1, 2], [2, 0], [0, 1]], numpy.int)
    fc_base = numpy.array([[1, 2, 3], [2, 0, 3], [0, 1, 3], [0, 1, 2]], numpy.int)
    #
    ne_dofs = len(ed_dofs)
    nf_dofs = len(fc_dofs)
    nc_dofs = len(c_dofs)
    #
    npts = 4 + ne_dofs*6 + nf_dofs*4 + nc_dofs
    pts_out = numpy.zeros((npts,3), numpy.float)
    #
    for i in range(4):
        pts_out[i,:] = pts_base[i,:]
    #
    for i in range(ne_dofs):
        node = ed_dofs[i]
        for j in range(6):
            pt1 = pts_base[ed_base[j,0],:]*node[0]
            pt2 = pts_base[ed_base[j,1],:]*node[1]
            pts_out[j+i*6+4,:] = pt1 + pt2
    #
    for i in range(nf_dofs):
        node = fc_dofs[i]
        for j in range(4):
            pt1 = pts_base[fc_base[j,0],:]*node[0]
            pt2 = pts_base[fc_base[j,1],:]*node[1]
            pt3 = pts_base[fc_base[j,2],:]*node[2]
            pts_out[j+i*4+6*ne_dofs+4,:] = pt1 + pt2 + pt3
    #
    for i in range(nc_dofs):
        node = c_dofs[i]
        pt1 = pts_base[0,:]
        for j in range(4):
            pt1 = pt1 + pts_base[j,:]*node[j]
        pts_out[i+nf_dofs*4+6*ne_dofs+4,:] = pt1
    #
    triang = scipy.spatial.Delaunay(pts_out)
    mask = check_tets(triang)
    return triang, mask
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    #----------------------------------------------------------------
    # Help Message
    #----------------------------------------------------------------
    import textwrap
    help_message = textwrap.dedent('''\
    Generate tessellation matrices for OFT
    ==========================================
    Tessellation is generated for:
    - Base triangle [(0.,0.), (1.,0.), (.5,1.)]
    - Base tetrahedron [(0.,0.,0.), (1.,0.,0.), (.5,1.,0.), (.5,.5,1.)]

    Tessellation is produced using a Delaunay triangulation of node points on the
    base elements using the SCIPY [https://www.scipy.org] "spatial" module.
    ==========================================
    ''')
    #----------------------------------------------------------------
    # Parse command line inputs
    #----------------------------------------------------------------
    nlevels=1
    linelength=30
    import argparse
    parser = argparse.ArgumentParser(description=help_message,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--nlevels', help='Desired tessellation order.', default=1)
    parser.add_argument("-l", "--linelength", help='Maximum FORTRAN line length.', default=30)
    args=parser.parse_args()
    nlevels=int(args.nlevels)
    linelength=int(args.linelength)
    #----------------------------------------------------------------
    # Import OFT Lagrange representation
    #----------------------------------------------------------------
    try:
        from basis_functions import lagrange
    except ImportError:
        print("Could not import Lagrange representation!!!")
        raise
    else:
        lag_rep = lagrange.lagrange_interp(nlevels)
    #----------------------------------------------------------------
    # Compute tesselation of triangle element
    #---------------------------------------------------------------
    tess, mask = tesselate_tri(nlevels,lag_rep.nodes_edge,lag_rep.nodes_face)
    print("===========START TRI TESSALLATION MATRIX===========\n")
    print(tess_fort_array(tess.vertices,mask,'tess',linelength))
    print("\n============END TRI TESSALLATION MATRIX============")
    #----------------------------------------------------------------
    # Compute tesselation of tetrahedron element
    #---------------------------------------------------------------
    tess, mask = tesselate_tet(nlevels,lag_rep.nodes_edge,lag_rep.nodes_face,lag_rep.nodes_cell)
    print("===========START TET TESSALLATION MATRIX===========\n")
    print(tess_fort_array(tess.vertices,mask,'tess',linelength))
    print("\n============END TET TESSALLATION MATRIX============\n")
