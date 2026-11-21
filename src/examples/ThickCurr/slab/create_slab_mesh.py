#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Create the native mesh for the ThickCurr slab example: a unit cube divided into
# structured tetrahedra, with a conducting slab (region 2) in 0.4 < z < 0.6 and
# vacuum (region 1) elsewhere.
import numpy
from OpenFUSIONToolkit.meshing import write_native_mesh

nx, ny, nz = 6, 6, 15  # Divisions in each direction (z grid lines align with slab faces)
slab = (0.4, 0.6)      # Extent of conducting slab in z

# Build grid points
x = numpy.linspace(0.0, 1.0, nx+1)
y = numpy.linspace(0.0, 1.0, ny+1)
z = numpy.linspace(0.0, 1.0, nz+1)
r = numpy.array([[x[i], y[j], z[k]] for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)])


def pind(i, j, k):
    return i + j*(nx+1) + k*(nx+1)*(ny+1)


# Split each hexahedral block into 6 tetrahedra (Kuhn triangulation, conforming between blocks)
perms = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
lc = []
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            for perm in perms:
                ind = [i, j, k]
                cell = [pind(*ind)]
                for axis in perm:
                    ind[axis] += 1
                    cell.append(pind(*ind))
                lc.append(cell)
lc = numpy.array(lc, dtype=numpy.int32)
# Ensure positive orientation
for cell in lc:
    if numpy.linalg.det(r[cell[1:], :] - r[cell[0], :]) < 0.0:
        cell[[2, 3]] = cell[[3, 2]]

# Set regions from cell centroids
reg = numpy.ones((lc.shape[0],), dtype=numpy.int32)
zc = r[lc, 2].mean(axis=1)
reg[(zc > slab[0]) & (zc < slab[1])] = 2

write_native_mesh('thickcurr_slab.h5', 'tet', r, lc, reg)
