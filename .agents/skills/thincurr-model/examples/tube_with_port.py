#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Build a ThinCurr model of an open cylindrical shell with a rectangular port,
using gmsh for the geometry and an XML file for resistivity and coils.

Pipeline: gmsh (OCC) surface -> .msh -> native .h5 -> holes (OFT_ThinCurr_holes)
-> oft_in XML. Writes into --output-dir:
  tube_with_port.msh, tube_with_port.h5          (intermediate, no holes)
  tube_with_port-homology.h5                     (the model mesh)
  tube_with_port.xml                             (resistivity and coils)

Usage:
  python tube_with_port.py [--output-dir DIR] [--h 0.05]
"""
import argparse
import os
import sys

import numpy

# Geometry [m]
RADIUS = 0.5        # Tube radius
LENGTH = 1.0        # Tube length (centered on Z=0)
PORT_WIDTH = 0.3    # Port extent in Y
PORT_HEIGHT = 0.2   # Port extent in Z
# Material: 1 mm stainless steel
ETA_VOL = 7.0E-7    # Volumetric resistivity [Ohm-m]
THICKNESS = 1.0E-3  # Wall thickness [m]


def build_gmsh_mesh(msh_file, h):
    """Mesh the tube surface (minus the port) with gmsh and write a v2.2 .msh file.

    Args:
        msh_file: Output gmsh file.
        h: Target edge length [m].
    """
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add('tube_with_port')
    occ = gmsh.model.occ
    cyl = occ.addCylinder(0.0, 0.0, -LENGTH/2.0, 0.0, 0.0, LENGTH, RADIUS)
    occ.synchronize()
    # Keep only the lateral surface of the cylinder volume (drop the planar end caps)
    bounds = gmsh.model.getBoundary([(3, cyl)], oriented=False)
    lateral = [dimtag for dimtag in bounds if gmsh.model.getType(*dimtag) != 'Plane']
    occ.remove([(3, cyl)])
    occ.remove([dimtag for dimtag in bounds if dimtag not in lateral])
    # Cut the port with a box straddling the wall at +X
    port = occ.addBox(RADIUS-0.2, -PORT_WIDTH/2.0, -PORT_HEIGHT/2.0, 0.4, PORT_WIDTH, PORT_HEIGHT)
    sheet, _ = occ.cut(lateral, [(3, port)])
    occ.synchronize()
    # One physical surface per region; gmsh only writes elements in physical groups
    gmsh.model.addPhysicalGroup(2, [tag for _, tag in sheet], 1)
    gmsh.option.setNumber('Mesh.MeshSizeMin', h)
    gmsh.option.setNumber('Mesh.MeshSizeMax', h)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.write(msh_file)
    gmsh.finalize()


def main():
    parser = argparse.ArgumentParser(description='Build the tube-with-port ThinCurr model')
    parser.add_argument('--output-dir', default='.', help='Output directory (default: .)')
    parser.add_argument('--h', type=float, default=0.05, help='Target edge length [m] (default: 0.05)')
    args = parser.parse_args()

    root_path = os.getenv('OFT_ROOTPATH')
    if root_path is not None:
        sys.path.append(os.path.join(root_path, 'python'))
    from OpenFUSIONToolkit._gmsh_cli import convert_gmsh_to_native
    from OpenFUSIONToolkit.io import write_oft_xml
    from OpenFUSIONToolkit.ThinCurr._homology_cli import compute_homology
    from OpenFUSIONToolkit.ThinCurr.coils import ThinCurr_Icoil, ThinCurr_XML

    os.makedirs(args.output_dir, exist_ok=True)
    base = os.path.join(args.output_dir, 'tube_with_port')

    # 1. Geometry and surface mesh
    build_gmsh_mesh(base+'.msh', args.h)
    # 2. Native mesh, keeping gmsh physical surfaces as regions
    convert_gmsh_to_native(base+'.msh', physical_regions=True)
    # 3. Holes and closures (writes <base>-homology.h5)
    compute_homology(base+'.h5')

    # 4. Resistivity (one value per region) and a coil pair for later driven runs
    xml = ThinCurr_XML()
    xml.set_eta_vol([ETA_VOL])
    xml.set_thickness([THICKNESS])
    coil = ThinCurr_Icoil(name='VF')
    coil.add_subcoil(RZ=numpy.r_[0.7, 0.3], scale=1.0)
    coil.add_subcoil(RZ=numpy.r_[0.7, -0.3], scale=1.0)
    xml.add_Icoil(coil)
    write_oft_xml([xml], base+'.xml')
    print('\nModel: {0}-homology.h5 + {0}.xml'.format(base))
    return 0


if __name__ == '__main__':
    sys.exit(main())
