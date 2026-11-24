MUG 2D: MHD equations     {#doc_mug_xmhd_2d_eqns}
================

[TOC]

This page documents the system of equations advanced by the nonlinear time integration in
\ref doc_mug2d_main "MUG 2D" (\ref xmhd_2d, driven by \ref xmhd_2d::run_simulation "run_simulation").
The complete system is defined by the nonlinear residual function \ref xmhd_2d::nlfun_apply
"nlfun_apply". The Jacobian built in \ref xmhd_2d::build_approx_jacobian "build_approx_jacobian"
(through \ref xmhd_2d::build_approx_jacobian_comp "build_approx_jacobian_comp" or
\ref xmhd_2d::build_approx_jacobian_incomp "build_approx_jacobian_incomp") is used as the preconditioner, and also as the Newton operator unless the matrix-free option
(`mfnk`) is enabled. It is the exact derivative of the residual, so it is not described separately. The 2D linearized model
(\ref xmhd_2d::run_lin_simulation "run_lin_simulation") is defined by this Jacobian and is summarized
in \ref doc_mug_xmhd_2d_eqns_linear.

The 2D implementation solves a reduced model compared to the
\ref doc_mug_xmhd_eqns "3D implementation": single-fluid resistive MHD with a single temperature,
isotropic thermal conduction, and no Hall, heating, upwinding, floor, or hyper-diffusion terms. It
also supports an incompressible model. All fields use scalar Lagrange elements of the same order (the
incompressible pressure uses one order lower).

Unlike the 3D implementation, the 2D model has no input-file namelist. Settings are components of the
\ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim" object and are set by the driving program. They are listed
in \ref doc_mug2d_main_params.

\section doc_mug_xmhd_2d_eqns_geom Geometry

The 2D mesh coordinates \f$ (x_1, x_2) \f$ are the physical in-plane coordinates \f$ (x, z) \f$. The
out-of-plane direction, \f$ y \f$, is ignorable (\f$ \partial / \partial y = 0 \f$). Two geometries
are supported:
 - **Planar** (default): a slab, uniform in \f$ y \f$, in Cartesian coordinates \f$ (x, y, z) \f$.
 - **Cylindrical** (`cyl_flag=T`): an axisymmetric domain in cylindrical coordinates
   \f$ (R, \phi, Z) \f$, where \f$ x_1 = R \f$, \f$ x_2 = Z \f$, and \f$ y \f$ is the toroidal angle \f$ \phi \f$.

Both \f$ (x, y, z) \f$ and \f$ (R, \phi, Z) \f$ are right-handed.

Vectors are stored in the order \f$ (x, y, z) \f$, or \f$ (R, \phi, Z) \f$, so the second component is out of
plane. In plot files vectors are written in mesh order, \f$ (x, z, y) \f$, and in cylindrical geometry
the plotted magnetic field is \f$ R B \f$ (`B*R`). The plotted flux `psi` and field are the totals, including
\f$ \psi_{vac} \f$ and \f$ B_0 \f$ (defined below). These are read from the step-0 restart file (`xmhd2d_00000.rst`),
where the run drivers save them. If they are not found there, a warning is printed and they are not included. In cylindrical geometry the weak forms are integrated with the volume element \f$ R \, dR \, dZ \f$,
which per radian of \f$ \phi \f$ is the true volume element, so all operators below are the true cylindrical operators.
In particular, for an axisymmetric scalar \f$ f \f$ and vector \f$ V \f$,

\f[ \nabla f = \partial_R f \, \hat{R} + \partial_Z f \, \hat{Z}, \quad
\nabla \cdot V = \partial_R V_R + \frac{V_R}{R} + \partial_Z V_Z, \quad
\Delta^* \psi = R \, \partial_R \left( \frac{1}{R} \partial_R \psi \right) + \partial_Z^2 \psi \f]

Every factor \f$ 1/R \f$ is evaluated as \f$ 1/(R + \epsilon_R) \f$, with \f$ \epsilon_R \f$ = `gs_epsilon` \f$ = 10^{-12} \f$ m.

\section doc_mug_xmhd_2d_eqns_conv Fields, units, and notation

| Field | Index | Units | Notes |
|-------|-------|-------|-------|
| \f$ n \f$ | 1 | m^-3 | Stored scaled by `den_scale` |
| \f$ V_x \f$, \f$ V_y \f$, \f$ V_z \f$ | 2-4 | m/s | \f$ V_y \f$ is out of plane (\f$ V_\phi \f$ in cylindrical geometry) |
| \f$ T \f$ | 5 | eV | Compressible model only |
| \f$ p \f$ | 5 | Pa | Incompressible model only. Uses Lagrange elements of one order lower |
| \f$ \tilde{\psi} \f$ | 6 | Wb/m (planar), Wb/rad (cylindrical) | Evolved poloidal flux |
| \f$ B_y \f$ | 7 | T (planar), T m (cylindrical) | In cylindrical geometry this is \f$ F = R B_\phi \f$ |

The magnetic field is

\f[ B = \nabla \psi \times \hat{y} + B_y \hat{y} + B_0 \quad \text{(planar)}, \qquad
B = \frac{1}{R} \left( \nabla \psi \times \hat{\phi} + F \hat{\phi} \right) + B_0 \quad \text{(cylindrical)} \f]

so that \f$ B_x = -\partial_z \psi \f$ and \f$ B_z = \partial_x \psi \f$ (planar), or
\f$ B_R = -\partial_Z \psi / R \f$, \f$ B_Z = \partial_R \psi / R \f$, and \f$ B_\phi = F / R \f$ (cylindrical).
Here:
 - \f$ \psi = \tilde{\psi} + \psi_{vac} \f$ is the total poloidal flux, where \f$ \tilde{\psi} \f$ is the evolved flux. \f$ \psi_{vac} \f$ is an optional vacuum flux (`psi_vac`, zero
   by default). It is not evolved, but is set by the driver, eg. from coil currents in the
   MUG-TokaMaker coupling (\ref mugtok_td). See \ref doc_mug_xmhd_2d_eqns_ind.
 - \f$ B_0 \f$ is a uniform background field (`B_0`, zero by default). Its in-plane part acts like the vacuum flux
   \f$ \psi_0 = B_{0,z} x - B_{0,x} z \f$ (planar) or \f$ \psi_0 = B_{0,Z} R^2 / 2 \f$ (cylindrical) added to
   \f$ \psi_{vac} \f$, and its out-of-plane part like a uniform \f$ B_{0,y} \f$ added to the initial \f$ B_y \f$. On a
   non-periodic domain these are equivalent. On a periodic domain \f$ \psi_0 \f$ is not periodic, so a uniform
   in-plane field can only be set with `B_0` (as in the Alfven wave regression test). In
   cylindrical geometry only a \f$ Z \f$ component is allowed, because a uniform field with \f$ R \f$ or \f$ \phi \f$
   components is not divergence- and curl-free, and a nonzero `B_0(1)` or `B_0(2)` aborts the run.

The following conventions are used throughout:
 - \f$ J_\mu \equiv \nabla \times B \f$, so the physical current density is \f$ J = J_\mu / \mu_0 \f$.
 - Temperatures are in eV and \f$ k = e \f$.
 - Resistivities, viscosities, and diffusivities are given in diffusivity units (m^2/s). In particular,
   \f$ \eta \f$ is the physical resistivity divided by \f$ \mu_0 \f$.
 - \f$ (\nabla V)_{ij} = \partial_i V_j \f$ and \f$ (\nabla \cdot A)_j = \sum_i \partial_i A_{ij} \f$, with the
   curvature terms of the cylindrical tensor operators included.
 - \f$ \phi \f$ (when not an angle) and \f$ \phi_p \f$ denote test functions in the main and pressure Lagrange spaces.

\subsection doc_mug_xmhd_2d_eqns_params Material properties

The material properties \f$ m_i \f$ (`m_i`), \f$ D \f$ (`D_diff`), \f$ \nu \f$ (`nu`), \f$ \chi \f$ (`chi`),
\f$ \gamma \f$ (`gamma`), \f$ \eta_1 \f$ (`eta(:,1)`, in-plane currents) and \f$ \eta_2 \f$ (`eta(:,2)`, out-of-plane
current) are given per mesh region and are uniform within a region. Their defaults and the checks made on them
are listed in \ref doc_mug2d_main_params.

Only the product \f$ m_i n \f$ (the mass density) enters the momentum equation, so \f$ n \f$ may be set to 1
and \f$ m_i \f$ to the mass density, as is done in the MUG-TokaMaker coupling.

Regions with `ignore_rmask=T` are skipped, so no equations are assembled in them. Their nodes must be
held fixed with Dirichlet conditions (see \ref doc_mug_xmhd_2d_eqns_bc).

\section doc_mug_xmhd_2d_eqns_cont Continuity

\f[ \frac{\partial n}{\partial t} + \nabla \cdot \left( n V \right) = \nabla \cdot \left( D \nabla n \right) \f]

 - The diffusion term is integrated by parts with the boundary term dropped (no diffusive flux where
   \f$ n \f$ is not fixed).
 - With `incomp=T`, \f$ n \f$ is not evolved. By default every density node is held fixed (see
   \ref doc_mug_xmhd_2d_eqns_bc).

\section doc_mug_xmhd_2d_eqns_mom Momentum

\subsection doc_mug_xmhd_2d_eqns_mom_comp Compressible model

\f[ \frac{\partial V}{\partial t} + V \cdot \nabla V = \frac{J_\mu \times B}{\mu_0 m_i n}
- \frac{\nabla p_h}{m_i n} + \frac{1}{n} \nabla \cdot \left( n \nu \nabla V \right) - g \hat{z} \f]

 - The pressure is \f$ p = 2 n k T \f$: electrons and ions are assumed to have the same temperature, \f$ T \f$.
   The force uses \f$ p_h = \sum_j 2 k n_j T_j \phi_j \f$, the interpolant of the nodal values of \f$ p \f$, rather than
   the product of the \f$ n \f$ and \f$ T \f$ fields (see \ref doc_mug_xmhd_2d_eqns_pres).
 - \f$ g \f$ = `g_accel` is a uniform gravitational acceleration along \f$ -\hat{z} \f$ (\f$ -\hat{Z} \f$ in cylindrical
   geometry), ie. along the negative second mesh coordinate. It is a body force per unit mass, so it does not depend on
   the solution and adds nothing to the Jacobian.
 - The Lorentz force is written as the divergence of the Maxwell stress,
   \f$ J_\mu \times B = \nabla \cdot \left( B B - \frac{1}{2} |B|^2 I \right) \f$. This is exact because
   \f$ \nabla \cdot B = 0 \f$ by construction. In the weak form it is integrated by parts onto
   \f$ \phi / (m_i n) \f$, with the boundary term dropped.
 - The viscosity is the kinematic ("Laplacian") form \f$ \nabla V \f$, not the rate-of-strain tensor, which is
   equivalent to `visc_type='kin'` in the 3D implementation. The dynamic viscosity is
   \f$ m_i n \nu \f$, so it is proportional to the local density. In contrast, the 3D implementation uses
   \f$ m_i n_0 \nu \f$. The term is integrated by parts with the boundary term dropped, so it gives a stress-free
   (zero normal derivative) condition wherever the velocity is not fixed. In cylindrical geometry \f$ \nabla \cdot \nabla V \f$ is the vector
   Laplacian, including the \f$ -V_R / R^2 \f$ and \f$ -V_\phi / R^2 \f$ terms.
 - In cylindrical geometry, \f$ V \cdot \nabla V \f$ includes the centrifugal (\f$ -V_\phi^2 / R \f$, \f$ R \f$ component)
   and Coriolis (\f$ V_R V_\phi / R \f$, \f$ \phi \f$ component) terms. The divergence of the Maxwell stress
   \f$ \mathcal{T} = B B - \frac{1}{2} |B|^2 I \f$ includes the curvature terms \f$ -\mathcal{T}_{\phi\phi} / R =
   -\left( B_\phi^2 - B_R^2 - B_Z^2 \right) / (2R) \f$ (\f$ R \f$ component) and \f$ \mathcal{T}_{R\phi} / R = B_R B_\phi / R \f$
   (\f$ \phi \f$ component).

\subsection doc_mug_xmhd_2d_eqns_mom_incomp Incompressible model

With `incomp=T`, the temperature equation is replaced by a divergence constraint, and field 5 is the pressure:

\f[ \frac{\partial V}{\partial t} + V \cdot \nabla V = \frac{J_\mu \times B}{\mu_0 m_i n}
- \frac{\nabla p}{m_i n} + \frac{1}{n} \nabla \cdot \left( n \nu \nabla V \right) - g \hat{z}, \qquad \nabla \cdot V = 0 \f]

 - The pressure uses Lagrange elements one order lower than the other fields (Taylor-Hood elements). It acts
   as a Lagrange multiplier.
 - The Lorentz force and pressure gradient are integrated by parts onto \f$ \phi / (m_i n) \f$, but the
   \f$ \nabla (1/n) \f$ terms are dropped, ie. the weak form is
   \f$ \int \left( \nabla \phi \cdot \mathcal{T} / \mu_0 - p \nabla \phi \right) / (m_i n) \f$ for the test function
   \f$ \phi \f$ of each velocity component, with \f$ \mathcal{T} \f$ the Maxwell stress (plus \f$ -\phi \, p / (R \, m_i n) \f$ in
   the \f$ R \f$ component in cylindrical geometry). This is exact only for uniform \f$ n \f$, which is the normal use
   (\f$ n \f$ is not evolved). The \f$ \nabla n \f$ part of the viscous term is still included.
 - The constraint is applied at the new time level only, as \f$ \int \phi_p \nabla \cdot V^{n+1} = 0 \f$ (see
   \ref doc_mug_xmhd_2d_eqns_time). Any divergence in the initial velocity is removed in the first step.
 - By default the pressure is fixed (to set its gauge) at the first global pressure node. This is only valid for a
   single connected MHD region. For multiple regions the driver must set the pressure boundary conditions.

\subsection doc_mug_xmhd_2d_eqns_pres Pressure discretization (compressible model)

With \f$ n_h = \sum_j n_j \phi_j \f$ and \f$ T_h = \sum_j T_j \phi_j \f$, the pressure force is evaluated from

\f[ p_h = \sum_j 2 k n_j T_j \phi_j, \qquad
\frac{\nabla p_h}{m_i n_h} = \frac{2 k}{m_i n_h} \sum_j n_j T_j \nabla \phi_j \f]

instead of the pointwise product \f$ 2 k \nabla (n_h T_h) / (m_i n_h) = 2 k \left( \nabla T_h + T_h \nabla n_h / n_h \right) / m_i \f$.
This is the group finite element formulation of Fletcher (\ref Fletcher1983 "1983", \ref Fletcher1984 "1984"), also called product
approximation (\ref Christie1981 "Christie et al. 1981") or interpolated coefficients (\ref Douglas1975 "Douglas and Dupont 1975"): a nonlinear group (here the product
\f$ n T \f$) is expanded in the same basis as the unknowns, using its nodal values. It is widely used for compressible
flow, for example in algebraic flux correction schemes (\ref Kuzmin2012 "Kuzmin et al. 2012"). The Jacobian uses the exact derivatives of
\f$ \nabla p_h / (m_i n_h) \f$ with respect to \f$ n_j \f$ and \f$ T_j \f$.

**Why it is used.** With all fields of the same order, the pointwise form is unstable for any stratified state (non-uniform
\f$ n \f$ and \f$ T \f$), even without gravity and at uniform pressure. Linear analysis of the discretization (and the 2D
code) shows a non-oscillatory spurious mode with a growth rate of order \f$ c_s / L_n \f$, where \f$ L_n \f$ is the density
scale length. It does not decrease with mesh refinement or a higher element order, and implicit time stepping does not damp
it. For typical laboratory parameters it is one to two orders of magnitude faster than buoyancy-driven instabilities, such
as the Rayleigh-Taylor instability. With \f$ p_h \f$, which lies in the finite element space, the mode is absent, and the
Rayleigh-Taylor growth rate of an exponential stratification agrees with the incompressible analytic rate to the
expected compressibility and dissipation corrections (\f$ < 0.2\% \f$; see `test_rt_2d`).

**What is unchanged.** For a uniform background (\f$ n_0 \f$, \f$ T_0 \f$ constant) the linearized pressure
\f$ 2 k (n_0 \delta T + T_0 \delta n) \f$ already lies in the finite element space, so linear waves about uniform states are
identical in both forms. The continuity and temperature equations, the magnetic force (the Maxwell stress is still formed
pointwise from \f$ B \f$), and the incompressible model (where \f$ p \f$ is itself a field) are not affected.

**Limitations.**
 - The scheme is not exactly well-balanced: when \f$ n \f$ and \f$ T \f$ are projected from a continuous hydrostatic
   equilibrium, \f$ \nabla p_h / (m_i n_h) + g \hat{z} \f$ is zero only to the accuracy of the interpolation. The
   remaining static imbalance drives small motions. Because hydrostatic balance is a small difference of large terms when
   \f$ m_i g L_n \ll k T \f$, the relative pressure error is amplified by about \f$ p / (m_i n g L_n) \f$. Schemes that preserve
   discrete hydrostatic states exactly need a dedicated (eg. hydrostatic) reconstruction
   (\ref Botta2004 "Botta et al. 2004", \ref Kaeppeli2014 "Käppeli and Mishra 2014"), which is not implemented.
 - For the same reason, particle diffusion (\f$ D \f$) and thermal conduction (\f$ \chi \f$) slowly move a stratified
   equilibrium away from hydrostatic balance, which drives flows. Small values (or zero) should be used when the
   equilibrium must stay at rest.
 - Only the thermal pressure uses the group form.
 - Dirichlet conditions on \f$ n \f$ and \f$ T \f$ at walls drive a similar spurious instability in stratified states
   (see \ref doc_mug_xmhd_2d_eqns_bc).

**References**
 - \anchor Fletcher1983 C. A. J. Fletcher, "The group finite element formulation", Comput. Methods Appl. Mech. Eng. 37,
   225-244 (1983), [doi:10.1016/0045-7825(83)90122-6](https://doi.org/10.1016/0045-7825(83)90122-6)
 - \anchor Fletcher1984 C. A. J. Fletcher, *Computational Galerkin Methods* (Springer, New York, 1984),
   [doi:10.1007/978-3-642-85949-6](https://doi.org/10.1007/978-3-642-85949-6)
 - \anchor Douglas1975 J. Douglas Jr. and T. Dupont, "The effect of interpolating the coefficients in nonlinear parabolic
   Galerkin procedures", Math. Comp. 29, 360-389 (1975),
   [doi:10.1090/S0025-5718-1975-0502033-7](https://doi.org/10.1090/S0025-5718-1975-0502033-7)
 - \anchor Christie1981 I. Christie, D. F. Griffiths, A. R. Mitchell, and J. M. Sanz-Serna, "Product approximation for
   non-linear problems in the finite element method", IMA J. Numer. Anal. 1, 253-266 (1981),
   [doi:10.1093/imanum/1.3.253](https://doi.org/10.1093/imanum/1.3.253)
 - \anchor Kuzmin2012 D. Kuzmin, M. Möller, and M. Gurris, "Algebraic flux correction II", in *Flux-Corrected Transport:
   Principles, Algorithms, and Applications*, 2nd ed. (Springer, Dordrecht, 2012), pp. 193-238,
   [doi:10.1007/978-94-007-4038-9_7](https://doi.org/10.1007/978-94-007-4038-9_7)
 - \anchor Botta2004 N. Botta, R. Klein, S. Langenberg, and S. Lützenkirchen, "Well balanced finite volume methods for
   nearly hydrostatic flows", J. Comput. Phys. 196, 539-565 (2004),
   [doi:10.1016/j.jcp.2003.11.008](https://doi.org/10.1016/j.jcp.2003.11.008)
 - \anchor Kaeppeli2014 R. Käppeli and S. Mishra, "Well-balanced schemes for the Euler equations with gravitation",
   J. Comput. Phys. 259, 199-219 (2014), [doi:10.1016/j.jcp.2013.11.026](https://doi.org/10.1016/j.jcp.2013.11.026)

\section doc_mug_xmhd_2d_eqns_temp Temperature

The temperature equation is solved only in the compressible model:

\f[ \frac{\partial T}{\partial t} + V \cdot \nabla T + (\gamma - 1) T \nabla \cdot V =
\frac{\gamma - 1}{n} \nabla \cdot \left( n \chi \nabla T \right) \f]

In energy form this is
\f$ \frac{n k}{\gamma - 1} \frac{d T}{d t} = -n k T \nabla \cdot V + \nabla \cdot \left( n k \chi \nabla T \right) \f$.

 - The thermal conduction is isotropic, with a constant diffusivity \f$ \chi \f$ in each region. The term is integrated
   by parts with the boundary term dropped.
 - The equation is solved after division by \f$ \gamma - 1 \f$, so the time-derivative term is \f$ \int \phi \, T / (\gamma - 1) \f$.
 - There is no Ohmic or viscous heating, and no separate electron temperature. The energy removed by
   resistivity and viscosity is not returned to the plasma.

\section doc_mug_xmhd_2d_eqns_ind Induction

The field is advanced with the resistive Ohm's law
\f$ E = -V \times B + \eta_2 J_{\mu,y} \hat{y} + \eta_1 J_{\mu,\text{pol}} \f$, where \f$ \eta_2 \f$ multiplies
the out-of-plane current and \f$ \eta_1 \f$ the in-plane current.

**Planar**

\f[ \frac{\partial \psi}{\partial t} + V \cdot \nabla \psi + \left( B_0 \times V \right)_y = \eta_2 \nabla^2 \tilde{\psi} \f]

\f[ \frac{\partial B_y}{\partial t} + \nabla \cdot \left( B_y V \right) = B_\text{pol} \cdot \nabla V_y
+ B_0 \cdot \nabla V_y - B_{0,y} \nabla \cdot V + \nabla \cdot \left( \eta_1 \nabla B_y \right) \f]

where \f$ B_\text{pol} = \nabla \psi \times \hat{y} \f$. In the implementation the stretching term appears as
\f$ -\left( \nabla \psi \times \nabla V_y \right)_y \f$, which is equal to it.

**Cylindrical**

\f[ \frac{\partial \psi}{\partial t} + V \cdot \nabla \psi + R \left( B_0 \times V \right)_\phi = \eta_2 \Delta^* \tilde{\psi} \f]

\f[ \frac{\partial F}{\partial t} + R^2 \nabla \cdot \left( \frac{F V}{R^2} \right) =
R^2 \left( B_\text{pol} + B_0 \right) \cdot \nabla \left( \frac{V_\phi}{R} \right)
+ R^2 \nabla \cdot \left( \frac{\eta_1}{R^2} \nabla F \right) \f]

where \f$ B_\text{pol} = \nabla \psi \times \hat{\phi} / R \f$, and
\f$ R^2 \nabla \cdot \left( F V / R^2 \right) = V \cdot \nabla F + F \left( \partial_R V_R + \partial_Z V_Z \right) - F V_R / R \f$.
The stretching term is \f$ R^2 B \cdot \nabla ( V_\phi / R ) = R B \cdot \nabla V_\phi - B_R V_\phi \f$, so only differential
rotation generates toroidal field. In the implementation the \f$ R B_\text{pol} \cdot \nabla V_\phi \f$ part appears as
\f$ -\left( \nabla \psi \times \nabla V_\phi \right)_\phi \f$.

Notes:
 - The resistive term uses only \f$ \tilde{\psi} \f$, not \f$ \psi \f$, so \f$ \psi_{vac} \f$ is assumed to be a vacuum field
   (\f$ \nabla^2 \psi_{vac} = 0 \f$, or \f$ \Delta^* \psi_{vac} = 0 \f$) that carries no current in the plasma. The time
   derivative and advection act on \f$ \psi \f$. If \f$ \psi_{vac} \f$ is changed between time steps (eg. by
   driven coil currents), the change acts as an applied loop voltage.
 - The \f$ \psi \f$ equation is solved after division by \f$ \eta_2 \f$ (and by \f$ R \f$ in cylindrical geometry), so the
   resistive term has unit weight and the time-derivative term is \f$ \int \phi \, \psi / \eta_2 \f$. This form is
   correct when \f$ \eta_2 \f$ jumps between regions (\f$ \nabla \tilde{\psi} \f$ is weakly continuous, ie. no surface current),
   but \f$ \eta_2 \f$ must be positive in every region. \f$ \eta_1 \f$ appears inside the divergence, which is also correct
   for piecewise-constant values.
 - Both resistive terms are integrated by parts with the boundary term dropped.

\section doc_mug_xmhd_2d_eqns_bc Boundary conditions

Each field has a node flag array (`n_bc`, `velx_bc`, `vely_bc`, `velz_bc`, `T_bc`, `psi_bc`, `by_bc`). Flagged
nodes are held at their initial values: their rows of the residual and Jacobian are replaced by the identity.
These are applied after the weak forms are assembled. If a flag array is not set before
\ref xmhd_2d::setup "setup", the defaults are:
 - All boundary nodes of the mesh, for every field. Every field is fixed at its initial boundary value, including
   no-slip for \f$ V \f$.
 - With `incomp=T`, every density node, and a single pressure node (see \ref doc_mug_xmhd_2d_eqns_mom_incomp).

Where a field is not fixed, the boundary terms dropped in the weak forms above act as natural boundary
conditions: no diffusive particle, heat, or magnetic flux, and no viscous or Maxwell stress.

\note In stratified equilibria (non-uniform \f$ n \f$ and \f$ T \f$) do not fix \f$ n \f$ and \f$ T \f$ at walls. Where the
normal velocity is zero, \f$ \partial n / \partial t = -\nabla \cdot (n V) \f$ is still not zero in general, so holding
\f$ n \f$ and \f$ T \f$ fixed over-constrains the continuity and temperature equations. This drives a spurious instability with a
growth rate of order \f$ c_s / L_n \f$. Set `n_bc` and `T_bc` to `.FALSE.` on such walls (natural, zero-flux
conditions). For a free-slip wall, fix only the normal and out-of-plane velocity components.

\section doc_mug_xmhd_2d_eqns_time Time discretization

Write the semi-discrete system as \f$ M \partial u / \partial t + R(u) = 0 \f$, where \f$ M \f$ includes the row weights
noted above (\f$ 1/(\gamma-1) \f$ for \f$ T \f$, \f$ 1/\eta_2 \f$ for \f$ \psi \f$, and the geometric factors in cylindrical geometry). Each time
step solves

\f[ M \left( u^{n+1} - u^n \right) + \Delta t \left[ \theta R(u^{n+1}) + (1 - \theta) R(u^n) \right] = 0 \f]

with Newton's method. Here \f$ \theta = 1 \f$ (backward Euler) by default and \f$ \theta = 1/2 \f$ (Crank-Nicolson) if
`timestep_cn=T`. Unlike the 3D implementation, all terms are treated with the same \f$ \theta \f$. In the
implementation, the right-hand side is formed by evaluating the residual with a time step of 0 (backward Euler)
or \f$ -\Delta t/2 \f$ (Crank-Nicolson), and the nonlinear solve then uses \f$ \Delta t \f$ or \f$ \Delta t/2 \f$.

Additional details:
 - The incompressibility constraint has no time derivative, so it is omitted from the right-hand-side evaluation
   (any evaluation with a time step \f$ \le 0 \f$), giving \f$ \int \phi_p \nabla \cdot V^{n+1} = 0 \f$.
 - With `mfnk=F` (the default) the Jacobian is used directly as the Newton operator. It is rebuilt from the
   current solution at the start of each step and every 4 Newton iterations. With `mfnk=T` a matrix-free
   finite-difference operator is used, and the Jacobian is used only as the preconditioner. It is then rebuilt at the
   start of every `pre_freq`-th step, and after the time step changes by more than 10%.
 - The time step is adapted after each step to give about `ittarget` linear iterations per step. It never exceeds the
   initial time step, and the run aborts if it falls below \f$ 10^{-3} \f$ times the initial value.

\section doc_mug_xmhd_2d_eqns_linear Linear runs

\ref xmhd_2d::run_lin_simulation "run_lin_simulation" advances a perturbation \f$ \delta u \f$ about a fixed
equilibrium \f$ u_0 \f$ (`u0`). It is available only in planar geometry (see \ref doc_mug2d_main_linear for how to set
up a linear run). The operator \f$ L \f$ is the Jacobian of \f$ R \f$ evaluated at \f$ u_0 \f$, built once before
the first step (and again before the second step with BDF2). The time step is fixed, and each step is a
single GMRES solve of

\f[ \left( M + \Delta t L \right) \delta u^{n+1} = M \delta u^n \f]

(backward Euler). With `timestep_cn=T`, the second-order backward difference formula (BDF2) is used after a first
backward Euler step:
\f$ \left( M + \frac{2}{3} \Delta t L \right) \delta u^{n+1} = \frac{1}{3} M \left( 4 \delta u^n - \delta u^{n-1} \right) \f$.

\f$ L \f$ is the exact linearization of the equations above, including the response to a density perturbation through
the \f$ 1/n \f$ factors of the force terms. For an equilibrium in force balance with gravity this response contains the buoyancy
force, \f$ -g \, \delta n / n_0 \, \hat{z} \f$ (from \f$ +\delta n \, \nabla p_0 / (m_i n_0^2) \f$ with
\f$ \nabla p_0 = -m_i n_0 g \hat{z} \f$). Gravity itself is constant and does not appear in \f$ L \f$. \f$ u_0 \f$ need not be an
exact steady state: \f$ R(u_0) \f$ is not included, so a non-equilibrium \f$ u_0 \f$ only changes the operator about which the
perturbation evolves.

In the history file (`xmhd_2d.hist`) of a linear run, both energies are those of the perturbation, so their growth
rates are twice the growth rate of the perturbation:
 - `ven` is the kinetic energy \f$ \int \frac{1}{2} m_i n_0 |\delta V|^2 \f$, using the equilibrium density.
 - `men` is \f$ \int |\delta B|^2 \f$ (the magnetic energy times \f$ 2 \mu_0 \f$), where \f$ \delta B \f$ is formed from the
   perturbed \f$ \psi \f$ and \f$ B_y \f$ only. It excludes \f$ B_0 \f$, \f$ \psi_v \f$ and the equilibrium field.

In nonlinear runs `ven` and `men` are the total kinetic and magnetic energies, with \f$ B_0 \f$ and \f$ \psi_v \f$ included in
\f$ B \f$. In both cases they are integrated over the domain with the volume element \f$ R \, dR \, dZ \f$ (per radian) in
cylindrical geometry, or per unit length in \f$ y \f$ in planar geometry.
