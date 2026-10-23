#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Build surface meshes used to test ThinCurr hole and closure selection (test_ThinCurr_homology.py)
#
# Run from this directory to regenerate the "homology_*.msh" files.
import numpy as np
import gmsh


def start(name, size):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.model.add(name)


def finish(filename):
    gmsh.model.mesh.renumberNodes()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(filename)
    gmsh.finalize()


def surface_mesh(name, size, build):
    '''Mesh the surfaces created by `build` (OpenCASCADE kernel)'''
    start(name, size)
    build()
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    finish(name + '.msh')


def discrete_mesh(name, r, lc):
    '''Write a mesh given explicitly as vertices and triangles (for surfaces that cannot be built from CAD)'''
    start(name, 1.0)
    tag = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, tag, list(range(1, r.shape[0]+1)), r.flatten().tolist())
    gmsh.model.mesh.addElementsByType(tag, 2, [], (np.asarray(lc)+1).flatten().tolist())
    finish(name + '.msh')


# Torus: 2 holes, 1 closure
surface_mesh('homology_torus', 0.12, lambda: gmsh.model.occ.addTorus(0.0, 0.0, 0.0, 1.0, 0.3))


# Torus with a port cut through its wall: 2 holes, no closures
def torus_port():
    torus = gmsh.model.occ.addTorus(0.0, 0.0, 0.0, 1.0, 0.3)
    gmsh.model.occ.synchronize()
    surfaces = gmsh.model.getBoundary([(3, torus)], oriented=False)
    gmsh.model.occ.remove([(3, torus)]) # Keep only the surface
    port = gmsh.model.occ.addCylinder(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.15) # Away from the torus seam
    gmsh.model.occ.cut(surfaces, [(3, port)])
surface_mesh('homology_torus_port', 0.1, torus_port)


# Plate with 4 openings: 4 holes (boundary), no closures
def plate():
    rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
    disks = [(2, gmsh.model.occ.addDisk(x, y, 0.0, 0.1, 0.1)) for x in (0.3, 0.7) for y in (0.3, 0.7)]
    gmsh.model.occ.cut([(2, rect)], disks)
surface_mesh('homology_plate', 0.05, plate)


# Surface of a slab with 4 through holes (closed surface of genus 4): 8 holes, 1 closure
def genus4():
    box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 0.2)
    cyls = [(3, gmsh.model.occ.addCylinder(x, y, 0.0, 0.0, 0.0, 0.2, 0.12)) for x in (0.3, 0.7) for y in (0.3, 0.7)]
    gmsh.model.occ.cut([(3, box)], cyls)
surface_mesh('homology_genus4', 0.06, genus4)


# Two separate tori: 4 holes, 2 closures
def two_tori():
    gmsh.model.occ.addTorus(0.0, 0.0, 0.0, 1.0, 0.3)
    gmsh.model.occ.addTorus(3.0, 0.0, 0.0, 1.0, 0.3)
surface_mesh('homology_two_tori', 0.15, two_tori)


# Three sheets sharing an edge (non-manifold)
def tee():
    sheets = [gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0), gmsh.model.occ.addRectangle(-1.0, 0.0, 0.0, 1.0, 1.0)]
    vert = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
    gmsh.model.occ.rotate([(2, vert)], 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -np.pi/2.0)
    gmsh.model.occ.fragment([(2, s) for s in sheets], [(2, vert)])
surface_mesh('homology_tee', 0.25, tee)


# Two sheets touching at a single vertex
def touching():
    sheets = [gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0), gmsh.model.occ.addRectangle(-1.0, -1.0, 0.0, 1.0, 1.0)]
    gmsh.model.occ.fragment([(2, sheets[0])], [(2, sheets[1])])
surface_mesh('homology_touching', 0.25, touching)


# Mobius strip (non-orientable)
n = 24
u, v = np.meshgrid(np.linspace(0.0, 2.0*np.pi, n, endpoint=False), np.r_[-0.2, 0.2], indexing='ij')
r = np.c_[((1.0+v*np.cos(u/2))*np.cos(u)).ravel(), ((1.0+v*np.cos(u/2))*np.sin(u)).ravel(), (v*np.sin(u/2)).ravel()]
def idx(i, j):
    return (i % n)*2 + j
lc = [[idx(i,0), idx(i+1,0), idx(i+1,1)] for i in range(n-1)] + [[idx(i,0), idx(i+1,1), idx(i,1)] for i in range(n-1)]
lc += [[idx(n-1,0), idx(0,1), idx(0,0)], [idx(n-1,0), idx(0,0), idx(n-1,1)]] # Joined with a half twist
discrete_mesh('homology_mobius', r, lc)
