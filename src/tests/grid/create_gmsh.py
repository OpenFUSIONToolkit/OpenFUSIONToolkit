import gmsh

# Build cylinder volume mesh using GMSH
gmsh.initialize()
gmsh.model.add("OFT_cylinder")
cylinder = gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)
gmsh.model.occ.synchronize()

boundary_points = gmsh.model.getBoundary([(3, cylinder)], combined=False, oriented=False, recursive=True)
gmsh.model.mesh.setSize(boundary_points, 0.4)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.model.mesh.generate(3)
gmsh.model.mesh.optimize()

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("cyl_gmsh_3D.msh")
# Save the same mesh in the legacy MSH v1 format to test that reader
gmsh.option.setNumber("Mesh.MshFileVersion", 1.0)
gmsh.write("cyl_gmsh_3D_v1.msh")
# Save the same mesh in the MEDIT format to test that reader
gmsh.write("cyl_gmsh_3D.mesh")
gmsh.finalize()

# Build cylinder surface mesh using GMSH
gmsh.initialize()
gmsh.model.add("OFT_cylinder")
cylinder = gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)
gmsh.model.occ.synchronize()

boundary_points = gmsh.model.getBoundary([(3, cylinder)], combined=False, oriented=False, recursive=True)
gmsh.model.mesh.setSize(boundary_points, 0.4)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.optimize()

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("cyl_gmsh_2D.msh")
gmsh.finalize()