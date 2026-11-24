Magnetohydrodynamics on Unstructured Grids (MUG): 3D linear/nonlinear extended MHD simulation package     {#doc_mhd_main}
================

[TOC]

\image html mug_logo.png width=25%

MUG solves the linearized or full nonlinear MHD equations in three dimensions using finite elements on
unstructured grids of tetrahedra or hexahedra.

MUG has two different implementations:
 - The primary implementation (\ref xmhd) using H(Curl) and Lagrange basis functions
 - An alternative implementation (\ref xmhd_lag) using only Lagrange basis functions and divergence cleaning

Both implementations support the same equations, options, and time integration (see \ref doc_mug_xmhd_eqns_lag
for the differences, which relate to the magnetic field representation).

\section doc_mhd_main_ex MUG Examples
The following examples illustrate usage of MUG to run time-dependent resistive and Hall MHD simulations. 

 - \subpage doc_mug_sph_ex1
 - \subpage doc_mug_sph_ex2
 - \subpage doc_mug_recon_ex
 - \subpage doc_mug_gyro_grav_ex
 - \subpage doc_mug_gyro_grav_lag_ex

\section doc_mhd_main_eqns MHD System
MUG solves the single- or two-temperature extended MHD equations: continuity, momentum (with optional
viscosity and gyroviscosity), induction with a resistive/Hall/electron-inertia Ohm's law, and
temperature evolution with anisotropic thermal conduction. The complete system advanced by the
\ref xmhd "primary implementation" is documented in:

 - \subpage doc_mug_xmhd_eqns

The same system is also advanced by the \ref xmhd_lag "Lagrange only" variant, see \ref doc_mug_xmhd_eqns_lag.

For 2D slab and axisymmetric simulations see \ref doc_mug2d_main "MUG 2D".
