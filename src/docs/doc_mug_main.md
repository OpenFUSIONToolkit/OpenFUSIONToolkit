Magnetohydrodynamics on Unstructured Grids (MUG): 3D linear/nonlinear extended MHD simulation package     {#doc_mhd_main}
================

[TOC]

MUG solves the linearized or full nonlinear MHD equations in three dimensions using finite elements on
unstructured grids of tetrahedra or hexahedra.

MUG has two different implementations:
 - The primary implementation (\ref xmhd) using H(Curl) and Lagrange basis functions
 - An alternative implementation (\ref xmhd_lag) using only Lagrange basis functions and divergence cleaning

\warning Not all capability is implemented in the \ref xmhd_lag "Lagrange only" variant.

\section doc_mhd_main_ex MUG Examples
The following examples illustrate usage of MUG to run time-dependent resistive and Hall MHD simulations. 

 - \subpage doc_mug_sph_ex1
 - \subpage doc_mug_sph_ex2
 - \subpage doc_mug_recon_ex

\section doc_mhd_main_eqns MHD System
Below is the complete system of equations supported by the \ref xmhd "primary implementation". Most
features are also supported by the \ref xmhd_lag "Lagrange only" variant.

\f[ \frac{\partial n}{\partial t} +  \nabla \cdot \left( n u \right) =
D \nabla^2 n \f]

\f[ \rho \left[ \frac{\partial u}{\partial t} +  u \cdot \nabla u \right]
= J \times B - \nabla nk(T_e + T_i) - \nabla \cdot \Pi \f]

\f[ \frac{\partial B}{\partial t} = - \nabla \times \left( -u \times B +
\eta J + \frac{1}{ne} \left( J \times B - \nabla nkT_e \right) + \frac{m_e}{ne^2}
\frac{\partial J}{\partial t} \right) \f]

\f[ \frac{n}{\gamma-1} \left[ \frac{\partial T_i}{\partial t} + u \cdot \nabla T_i \right] =
-nkT_i \nabla \cdot u - \nabla \cdot q_i  + Q_i \f]

For two temperature
\f[ \frac{n}{\gamma-1} \left[ \frac{\partial T_e}{\partial t} + u_e \cdot \nabla T_e \right] =
-nkT_e \nabla \cdot u_e - \nabla \cdot q_e  + Q_e \f]

\f[ u_e = u - \frac{J}{ne} \f]

**Thermal transport**

\f[ q_s = - n \left[ \chi_{\parallel,s} \hat{b} \hat{b} + \chi_{\perp,s} \left( I
- \hat{b} \hat{b} \right) \right] \cdot \nabla T_s \f]

With `xmhd_brag=T` and single temperature (\f$ T_e = T_i \f$)
\f[ \chi_{\parallel,i} = \chi_{\parallel,e}, \chi_{\perp,i} = min(\chi_{\perp,i},\chi_{\parallel,i}) \f]

**Heat sources**

For single temperature (\f$ T_e = T_i \f$)
\f[ Q_i = \left[ \eta J^2 - \left( \nabla u \right)^T : \Pi \right] / 2 \f]

For two temperature
\f[ Q_i = - \left( \nabla u \right)^T : \Pi, Q_e = \eta J^2 \f]

**Viscosity**

\f[ W = \left( \nabla u + (\nabla u)^T - \frac{2}{3} I \nabla \cdot u \right) \f]

With `visc_type='kin'`
\f[ \Pi = - \nu \nabla u \f]

With `visc_type='iso'`
\f[ \Pi = - \nu W \f]

With `visc_type='ani'`
\f[ \Pi = - \left[ \nu_{\parallel} \hat{b} \hat{b}
+ \nu_{\perp} \left( I - \hat{b} \hat{b} \right) \right] \cdot W \f]

where \f$ \mathcal{P} \f$, \f$ \mathcal{S} \f$, and \f$ \mathcal{C} \f$ are axisymmetric domains corresponding to the plasma, passive conducting structures (eg. vacuum vessels), coils respectively.