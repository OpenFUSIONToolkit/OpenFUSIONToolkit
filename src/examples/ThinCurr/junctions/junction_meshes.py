'''Generate the gmsh meshes used in the ThinCurr junctions example (requires the gmsh Python module)

  thincurr_ex-tplate.msh  T-shaped plate: plate halves A (x < 0, region 1) and B (x > 0, region 2) and a
                          stem C (x = 0, region 3) meeting along one edge
  thincurr_ex-ring.msh    Tube (radius 0.5 m) in an upper (region 1) and lower (region 2) half with a flange
                          (region 3, outer radius 0.8 m) around its middle: a closed junction

Each sheet is a gmsh physical surface, written in the version 2.2 format.
'''
import gmsh

H = 0.05


def write(name, groups):
    for tag, surfs in groups.items():
        gmsh.model.addPhysicalGroup(2, surfs, tag)
    gmsh.option.setNumber('Mesh.MeshSizeMin', H)
    gmsh.option.setNumber('Mesh.MeshSizeMax', H)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.write(name)


gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)

# T-plate: three rectangles sharing the edge x = z = 0, y in [0, 1]
gmsh.model.add('tplate')
geo = gmsh.model.geo
p0, p1 = geo.addPoint(0.0, 0.0, 0.0, H), geo.addPoint(0.0, 1.0, 0.0, H)
edge = geo.addLine(p0, p1)
sheets = {}
for tag, (q0, q1) in {1: ((-0.5, 0.0, 0.0), (-0.5, 1.0, 0.0)), 2: ((0.5, 0.0, 0.0), (0.5, 1.0, 0.0)),
                      3: ((0.0, 0.0, 0.5), (0.0, 1.0, 0.5))}.items():
    a, b = geo.addPoint(*q0, H), geo.addPoint(*q1, H)
    loop = geo.addCurveLoop([edge, geo.addLine(p1, b), geo.addLine(b, a), geo.addLine(a, p0)])
    sheets[tag] = [geo.addPlaneSurface([loop])]
geo.synchronize()
write('thincurr_ex-tplate.msh', sheets)

# Tube with a flange around its middle
gmsh.model.add('ring')
occ = gmsh.model.occ
wires = [occ.addWire([occ.addCircle(0.0, 0.0, z, 0.5)]) for z in (-0.5, 0.0, 0.5)]
surfs = occ.addThruSections(wires[:2], makeSolid=False, makeRuled=True)
surfs += occ.addThruSections(wires[1:], makeSolid=False, makeRuled=True)
flange, _ = occ.cut([(2, occ.addDisk(0.0, 0.0, 0.0, 0.8, 0.8))], [(2, occ.addDisk(0.0, 0.0, 0.0, 0.5, 0.5))])
occ.fragment(surfs + flange, [])
occ.synchronize()
groups = {1: [], 2: [], 3: []}
for dim, tag in gmsh.model.getEntities(2):
    z = occ.getCenterOfMass(dim, tag)[2]
    groups[3 if abs(z) < 1.0E-6 else (1 if z > 0.0 else 2)].append(tag)
write('thincurr_ex-ring.msh', groups)
gmsh.finalize()
