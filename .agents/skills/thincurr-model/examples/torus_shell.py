#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Build a ThinCurr model of a closed, shaped toroidal shell from a structured
(theta, phi) grid, storing resistivity in the mesh file (no XML needed).

Pipeline: analytic surface -> build_triangles_from_grid() -> write_ThinCurr_mesh()
(with resistivity) -> holes and closure (OFT_ThinCurr_holes). Writes into --output-dir:
  torus_shell.h5            (intermediate, no holes)
  torus_shell-homology.h5   (the model mesh: 2 holes, 1 closure)

Usage:
  python torus_shell.py [--output-dir DIR] [--ntheta 32] [--nphi 64]
"""
import argparse
import os
import sys

import numpy

# Shell shape [m]: Miller-like cross-section
R0 = 1.0
A = 0.4
KAPPA = 1.5
DELTA = 0.3
# Material: 2 cm thick copper-chromium-zirconium
ETA_VOL = 2.0E-8    # Volumetric resistivity [Ohm-m]
THICKNESS = 2.0E-2  # Wall thickness [m]


def main():
    parser = argparse.ArgumentParser(description='Build the shaped torus-shell ThinCurr model')
    parser.add_argument('--output-dir', default='.', help='Output directory (default: .)')
    parser.add_argument('--ntheta', type=int, default=32, help='Poloidal grid points (default: 32)')
    parser.add_argument('--nphi', type=int, default=64, help='Toroidal grid points (default: 64)')
    args = parser.parse_args()

    root_path = os.getenv('OFT_ROOTPATH')
    if root_path is not None:
        sys.path.append(os.path.join(root_path, 'python'))
    from OpenFUSIONToolkit.ThinCurr._homology_cli import compute_homology
    from OpenFUSIONToolkit.ThinCurr.meshing import build_triangles_from_grid, write_ThinCurr_mesh

    os.makedirs(args.output_dir, exist_ok=True)
    base = os.path.join(args.output_dir, 'torus_shell')

    # 1. Structured grid [nphi, ntheta, 3]; neither angle repeats its end point (both wrap)
    theta = numpy.linspace(0.0, 2.0*numpy.pi, args.ntheta, endpoint=False)
    phi = numpy.linspace(0.0, 2.0*numpy.pi, args.nphi, endpoint=False)
    r_pol = R0 + A*numpy.cos(theta + numpy.arcsin(DELTA)*numpy.sin(theta))
    z_pol = KAPPA*A*numpy.sin(theta)
    grid = numpy.zeros((args.nphi, args.ntheta, 3))
    for i, phi_val in enumerate(phi):
        grid[i, :, 0] = r_pol*numpy.cos(phi_val)
        grid[i, :, 1] = r_pol*numpy.sin(phi_val)
        grid[i, :, 2] = z_pol
    r, lc, _, _ = build_triangles_from_grid(grid, wrap_n=True, wrap_m=True)

    # 2. Native mesh with resistivity stored in the file (one value per region)
    reg = numpy.ones((lc.shape[0],), dtype=numpy.int32)
    write_ThinCurr_mesh(base+'.h5', r, lc, reg, eta_vol=numpy.r_[ETA_VOL], thickness=numpy.r_[THICKNESS])

    # 3. Holes (poloidal + toroidal) and the closure of the closed surface
    compute_homology(base+'.h5')
    print('\nModel: {0}-homology.h5'.format(base))
    return 0


if __name__ == '__main__':
    sys.exit(main())
