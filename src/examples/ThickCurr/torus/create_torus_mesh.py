#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Create the native mesh for the ThickCurr torus example using GMSH: a solid conducting
# torus (region 2) at the center of a spherical vacuum domain (region 1). A disk in the
# z=0 plane is embedded in the hole of the torus so that the flux through the hole can be
# measured with a flux probe. A second-order (curved) mesh is generated to accurately
# represent the curved surfaces.
#
# Requires the `gmsh` python package (eg. `pip install gmsh`)
import gmsh
from OpenFUSIONToolkit._gmsh_cli import convert_gmsh_to_native

R0 = 1.0          # Major radius of torus [m]
a = 0.5           # Minor radius of torus [m]
R_sphere = 3.0    # Radius of spherical domain [m]
r_probe = 0.45    # Radius of embedded flux probe disk [m]
h_torus = 0.3     # Mesh size in and near the torus [m]
h_far = 1.5       # Mesh size near the outer boundary [m]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("thickcurr_torus")
occ = gmsh.model.occ
# Create geometry: torus centered at the origin with axis along z, inside a sphere
sphere = occ.addSphere(0.0, 0.0, 0.0, R_sphere)
torus = occ.addTorus(0.0, 0.0, 0.0, R0, a)
disk = occ.addDisk(0.0, 0.0, 0.0, r_probe, r_probe)
out, out_map = occ.fragment([(3, sphere)], [(3, torus)])
occ.synchronize()
volumes = [tag for (dim, tag) in out if dim == 3]
# Identify torus volume (smallest) and vacuum volume
vol_sizes = {tag: occ.getMass(3, tag) for tag in volumes}
torus_vol = min(vol_sizes, key=vol_sizes.get)
vac_vol = max(vol_sizes, key=vol_sizes.get)
# Embed probe disk in the vacuum region
gmsh.model.mesh.embed(2, [disk], 3, vac_vol)
gmsh.model.addPhysicalGroup(3, [vac_vol], 1)  # Vacuum region
gmsh.model.addPhysicalGroup(3, [torus_vol], 2)  # Torus region
# Set mesh size, refined near the torus
dist = gmsh.model.mesh.field.add("Distance")
torus_surfs = [tag for (dim, tag) in gmsh.model.getBoundary([(3, torus_vol)], oriented=False)]
gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", torus_surfs)
thresh = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresh, "InField", dist)
gmsh.model.mesh.field.setNumber(thresh, "SizeMin", h_torus)
gmsh.model.mesh.field.setNumber(thresh, "SizeMax", h_far)
gmsh.model.mesh.field.setNumber(thresh, "DistMin", 0.25)
gmsh.model.mesh.field.setNumber(thresh, "DistMax", 1.75)
gmsh.model.mesh.field.setAsBackgroundMesh(thresh)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
# Generate second-order mesh
gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("thickcurr_torus.msh")
gmsh.finalize()

# Convert to native format using physical groups as regions
convert_gmsh_to_native("thickcurr_torus.msh", "thickcurr_torus.h5", physical_regions=True)
