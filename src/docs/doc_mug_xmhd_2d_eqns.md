MUG: MHD equations (2D Lagrange implementation)     {#doc_mug_xmhd_2d_eqns}
================

[TOC]

This page documents the system of equations advanced by the nonlinear time integration in the
2D MUG implementation (\ref xmhd_2d, driven by \ref xmhd_2d::run_simulation "run_simulation").
The complete system is defined by the nonlinear residual function \ref xmhd_2d::nlfun_apply
"nlfun_apply". The Jacobian built in \ref xmhd_2d::build_approx_jacobian "build_approx_jacobian"
is used as the preconditioner, and also as the Newton operator unless the matrix-free option
(`mfnk`) is enabled. It is constructed as the derivative of the residual, so it is not described
separately. The 2D linearized model
(\ref xmhd_2d::run_lin_simulation "run_lin_simulation") is defined by this Jacobian and is summarized
in \ref doc_mug_xmhd_2d_eqns_linear.

The 2D implementation solves a reduced model compared to the
\ref doc_mug_xmhd_eqns "3D implementation": single-fluid resistive MHD with a single temperature,
isotropic thermal conduction, and no Hall, heating, upwinding, floor, or hyper-diffusion terms. It
also supports an incompressible model. All fields use scalar Lagrange elements of the same order (the
incompressible pressure uses one order lower).

Unlike the 3D implementation, the 2D model has no input-file namelist. Settings are components of the
\ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim" object and are set by the driving program. Material
properties are given per mesh region (see \ref doc_mug_xmhd_2d_eqns_params).

\section doc_mug_xmhd_2d_eqns_geom Geometry

The 2D mesh coordinates \f$ (x_1, x_2) \f$ are the physical in-plane coordinates \f$ (x, z) \f$. The
out-of-plane direction, \f$ y \f$, is ignorable (\f$ \partial / \partial y = 0 \f$). Two geometries
are supported:
 - **Planar** (default): a slab, uniform in \f$ y \f$, in Cartesian coordinates \f$ (x, y, z) \f$.
 - **Cylindrical** (`cyl_flag=T`): an axisymmetric domain in cylindrical coordinates
   \f$ (R, \phi, Z) \f$, where \f$ x_1 = R \f$, \f$ x_2 = Z \f$, and \f$ y \f$ is the toroidal angle \f$ \phi \f$.

Vectors are stored in the order \f$ (x, y, z) \f$, or \f$ (R, \phi, Z) \f$, so the second component is out of
plane. In plot files vectors are written in mesh order, \f$ (x, z, y) \f$, and in cylindrical geometry
the plotted magnetic field is \f$ R B \f$ (`B*R`). In cylindrical geometry the weak forms are integrated with the volume element \f$ R \, dR \, dZ \f$,
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
| \f$ \psi \f$ | 6 | Wb/m (planar), Wb/rad (cylindrical) | Poloidal flux |
| \f$ B_y \f$ | 7 | T (planar), T m (cylindrical) | In cylindrical geometry this is \f$ F = R B_\phi \f$ |

The magnetic field is

\f[ B = \nabla \psi_t \times \hat{y} + B_y \hat{y} + B_0 \quad \text{(planar)}, \qquad
B = \frac{1}{R} \left( \nabla \psi_t \times \hat{\phi} + F \hat{\phi} \right) + B_0 \quad \text{(cylindrical)} \f]

so that \f$ B_x = -\partial_z \psi_t \f$ and \f$ B_z = \partial_x \psi_t \f$ (planar), or
\f$ B_R = -\partial_Z \psi_t / R \f$, \f$ B_Z = \partial_R \psi_t / R \f$, and \f$ B_\phi = F / R \f$ (cylindrical).
Here:
 - \f$ \psi_t = \psi + \psi_v \f$ is the total poloidal flux. \f$ \psi_v \f$ is an optional vacuum flux (`psi_vac`, zero
   by default). It is not evolved, but is set by the driver, eg. from coil currents in the
   MUG-TokaMaker coupling (\ref mugtok_td). See \ref doc_mug_xmhd_2d_eqns_ind.
 - \f$ B_0 \f$ is a uniform background field (`B_0`, zero by default). See the warnings in
   \ref doc_mug_xmhd_2d_eqns_ind.

The following conventions are used throughout:
 - \f$ J_\mu \equiv \nabla \times B \f$, so the physical current density is \f$ J = J_\mu / \mu_0 \f$.
 - Temperatures are in eV and \f$ k = e \f$.
 - Resistivities, viscosities, and diffusivities are given in diffusivity units (m^2/s). In particular,
   \f$ \eta \f$ is the physical resistivity divided by \f$ \mu_0 \f$.
 - \f$ (\nabla V)_{ij} = \partial_i V_j \f$ and \f$ (\nabla \cdot A)_j = \sum_i \partial_i A_{ij} \f$, with the
   curvature terms of the cylindrical tensor operators included.
 - \f$ \phi \f$ (when not an angle) and \f$ \phi_p \f$ denote test functions in the main and pressure Lagrange spaces.

\subsection doc_mug_xmhd_2d_eqns_params Parameters

Material properties are arrays with one entry per mesh region, and are taken as uniform within a region.
All of them except \f$ m_i \f$ default to -1 and must be set before \ref xmhd_2d::setup "setup" is called.
No check is made that they have been set.

| Symbol | Component | Description |
|--------|-----------|-------------|
| \f$ m_i \f$ | `m_i` | Ion mass [kg] (default: proton mass) |
| \f$ D \f$ | `D_diff` | Particle diffusivity [m^2/s] |
| \f$ \nu \f$ | `nu` | Kinematic viscosity [m^2/s] |
| \f$ \chi \f$ | `chi` | Thermal diffusivity [m^2/s] |
| \f$ \gamma \f$ | `gamma` | Adiabatic index |
| \f$ \eta_1 \f$ | `eta(:,1)` | Resistivity for in-plane (poloidal) currents [m^2/s] |
| \f$ \eta_2 \f$ | `eta(:,2)` | Resistivity for out-of-plane (toroidal) current [m^2/s]. Must be positive |

Only the product \f$ m_i n \f$ (the mass density) enters the momentum equation, so \f$ n \f$ may be set to 1
and \f$ m_i \f$ to the mass density, as is done in the MUG-TokaMaker coupling.

Regions with `ignore_rmask=T` are skipped, so no equations are assembled in them. Their nodes must be
held fixed with Dirichlet conditions (see \ref doc_mug_xmhd_2d_eqns_bc).

Other settings of \ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim" referenced below are `cyl_flag`
(cylindrical geometry), `incomp` (incompressible model), `timestep_cn` (time-centered stepping),
`mfnk` (matrix-free Newton-Krylov), `den_scale` (density scale, default \f$ 10^{19} \f$), and `B_0`.

\section doc_mug_xmhd_2d_eqns_cont Continuity

\f[ \frac{\partial n}{\partial t} + \nabla \cdot \left( n V \right) = \nabla \cdot \left( D \nabla n \right) \f]

 - The diffusion term is integrated by parts with the boundary term dropped (no diffusive flux where
   \f$ n \f$ is not fixed).
 - With `incomp=T`, \f$ n \f$ is not evolved. By default every density node is held fixed (see
   \ref doc_mug_xmhd_2d_eqns_bc).

\section doc_mug_xmhd_2d_eqns_mom Momentum

\subsection doc_mug_xmhd_2d_eqns_mom_comp Compressible model

\f[ \frac{\partial V}{\partial t} + V \cdot \nabla V = \frac{J_\mu \times B}{\mu_0 m_i n}
- \frac{\nabla \left( 2 n k T \right)}{m_i n} + \frac{1}{n} \nabla \cdot \left( n \nu \nabla V \right) \f]

 - The pressure is \f$ p = 2 n k T \f$: electrons and ions are assumed to have the same temperature, \f$ T \f$.
 - The Lorentz force is written as the divergence of the Maxwell stress,
   \f$ J_\mu \times B = \nabla \cdot \left( B B - \frac{1}{2} |B|^2 I \right) \f$. This is exact because
   \f$ \nabla \cdot B = 0 \f$ by construction. In the weak form it is integrated by parts onto
   \f$ \phi / (m_i n) \f$, with the boundary term dropped.
 - The viscosity is the kinematic form \f$ \nabla V \f$, not the rate-of-strain tensor. The dynamic viscosity is
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
- \frac{\nabla p}{m_i n} + \frac{1}{n} \nabla \cdot \left( n \nu \nabla V \right), \qquad \nabla \cdot V = 0 \f]

 - The pressure uses Lagrange elements one order lower than the other fields (Taylor-Hood elements). It acts
   as a Lagrange multiplier.
 - The Lorentz force and pressure gradient are integrated by parts onto \f$ \phi / (m_i n) \f$, but the
   \f$ \nabla (1/n) \f$ terms are dropped, ie. the weak form is
   \f$ \int \left( \nabla \phi \cdot \mathcal{T} / \mu_0 - p \nabla \phi \right) / (m_i n) \f$ for the test function
   \f$ \phi \f$ of each velocity component, with \f$ \mathcal{T} \f$ the Maxwell stress (plus \f$ -\phi \, p / (R \, m_i n) \f$ in
   the \f$ R \f$ component in cylindrical geometry). This is exact only for uniform \f$ n \f$, which is the normal use
   (\f$ n \f$ is not evolved). The \f$ \nabla n \f$ part of the viscous term is still included.
 - The constraint is applied as \f$ \int \phi_p \nabla \cdot V^{n+1} = \int \phi_p \nabla \cdot V^n \f$ (see
   \ref doc_mug_xmhd_2d_eqns_time). The divergence of the initial velocity is therefore preserved, not removed,
   so the initial velocity must be divergence-free.
 - By default the pressure is fixed (to set its gauge) at the first global pressure node. This is only valid for a
   single connected MHD region. For multiple regions the driver must set the pressure boundary conditions.

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

\f[ \frac{\partial \psi_t}{\partial t} + V \cdot \nabla \psi_t + \left( B_0 \times V \right)_y = \eta_2 \nabla^2 \psi \f]

\f[ \frac{\partial B_y}{\partial t} + \nabla \cdot \left( B_y V \right) = B_\text{pol} \cdot \nabla V_y
+ \nabla \cdot \left( \eta_1 \nabla B_y \right) \f]

where \f$ B_\text{pol} = \nabla \psi_t \times \hat{y} \f$. In the implementation the stretching term appears as
\f$ -\left( \nabla \psi_t \times \nabla V_y \right)_y \f$, which is equal to it.

**Cylindrical**

\f[ \frac{\partial \psi_t}{\partial t} + V \cdot \nabla \psi_t + \left( B_0 \times V \right)_\phi = \eta_2 \Delta^* \psi \f]

\f[ \frac{\partial F}{\partial t} + R^2 \nabla \cdot \left( \frac{F V}{R^2} \right) = R B_\text{pol} \cdot \nabla V_\phi
+ R^2 \nabla \cdot \left( \frac{\eta_1}{R^2} \nabla F \right) \f]

where \f$ B_\text{pol} = \nabla \psi_t \times \hat{\phi} / R \f$, and
\f$ R^2 \nabla \cdot \left( F V / R^2 \right) = V \cdot \nabla F + F \left( \partial_R V_R + \partial_Z V_Z \right) - F V_R / R \f$.

Notes:
 - The resistive term uses only \f$ \psi \f$, not \f$ \psi_t \f$, so \f$ \psi_v \f$ is assumed to be a vacuum field
   (\f$ \nabla^2 \psi_v = 0 \f$, or \f$ \Delta^* \psi_v = 0 \f$) that carries no current in the plasma. The time
   derivative and advection act on \f$ \psi_t \f$. If \f$ \psi_v \f$ is changed between time steps (eg. by
   driven coil currents), the change acts as an applied loop voltage.
 - The \f$ \psi \f$ equation is solved after division by \f$ \eta_2 \f$ (and by \f$ R \f$ in cylindrical geometry), so the
   resistive term has unit weight and the time-derivative term is \f$ \int \phi \, \psi_t / \eta_2 \f$. This form is
   correct when \f$ \eta_2 \f$ jumps between regions (\f$ \nabla \psi \f$ is weakly continuous, ie. no surface current),
   but \f$ \eta_2 \f$ must be positive in every region. \f$ \eta_1 \f$ appears inside the divergence, which is also correct
   for piecewise-constant values.
 - Both resistive terms are integrated by parts with the boundary term dropped.

\warning In cylindrical geometry the stretching term should be \f$ R^2 B_\text{pol} \cdot \nabla ( V_\phi / R ) =
R B_\text{pol} \cdot \nabla V_\phi - B_R V_\phi \f$. The \f$ -B_R V_\phi \f$ part is not included in either the residual or the Jacobian,
so a rigid toroidal rotation (\f$ V_\phi = \Omega R \f$) generates toroidal field when \f$ B_R \neq 0 \f$. This
affects runs that evolve \f$ V_\phi \f$.

\warning The background field \f$ B_0 \f$ is not treated consistently. It is included in the Lorentz force
and in the \f$ \psi \f$ equation (where its in-plane part acts like the flux \f$ \psi_0 = B_{0,z} x - B_{0,x} z \f$),
but it is absent from the \f$ B_y \f$ (\f$ F \f$) equation, which is missing the terms \f$ B_0 \cdot \nabla V_y - B_{0,y} \nabla \cdot V \f$.
In cylindrical geometry a uniform \f$ B_0 \f$ with \f$ R \f$ or \f$ \phi \f$ components is not a valid
(divergence- and curl-free) field, and the \f$ (B_0 \times V)_\phi \f$ term is missing a factor of \f$ R \f$. \f$ B_0 \f$ is
also not included in the plotted field. A uniform in-plane field is better represented
through the initial \f$ \psi \f$, and a uniform guide field through the initial \f$ B_y \f$.

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
 - The incompressibility constraint has no time derivative or \f$ \Delta t \f$ factor, so it appears unchanged on
   both sides, giving \f$ \int \phi_p \nabla \cdot V^{n+1} = \int \phi_p \nabla \cdot V^n \f$.
 - With `mfnk=F` (the default) the Jacobian is used directly as the Newton operator. It is rebuilt from the
   current solution at the start of each step and every 4 Newton iterations. With `mfnk=T` a matrix-free
   finite-difference operator is used, and the Jacobian is used only as the preconditioner. It is then rebuilt at the
   start of every `pre_freq`-th step, and after the time step changes by more than 10%.
 - The time step is adapted after each step to give about `ittarget` linear iterations per step. It never exceeds the
   initial time step, and the run aborts if it falls below \f$ 10^{-3} \f$ times the initial value.

\section doc_mug_xmhd_2d_eqns_linear Linear runs

\ref xmhd_2d::run_lin_simulation "run_lin_simulation" advances a perturbation \f$ \delta u \f$ about a fixed
equilibrium \f$ u_0 \f$ (`u0`). It is available only in planar geometry. The operator \f$ L \f$ is the Jacobian of
\f$ R \f$ evaluated at \f$ u_0 \f$, and each step solves

\f[ \left( M + \Delta t L \right) \delta u^{n+1} = M \delta u^n \f]

(backward Euler). With `timestep_cn=T`, the second-order backward difference formula (BDF2) is used after a first
backward Euler step:
\f$ \left( M + \frac{2}{3} \Delta t L \right) \delta u^{n+1} = \frac{1}{3} M \left( 4 \delta u^n - \delta u^{n-1} \right) \f$.

The Jacobian is constructed as the linearization of the equations above, except for the response to a density perturbation
\f$ \delta n \f$. There, the \f$ 1/n \f$ dependence of the force and heat-flux terms is evaluated with the steady-state equations
for \f$ u_0 \f$ in place of the force and conduction terms themselves. For momentum,
\f[ \delta \left[ \frac{1}{m_i n} \left( \frac{J_\mu \times B}{\mu_0} - \nabla p \right) \right]_{\delta n}
\approx -\frac{\delta n}{n_0} \left( V_0 \cdot \nabla V_0 \right) - \frac{2 k \nabla \left( T_0 \delta n \right)}{m_i n_0} \f]
and the viscous terms in \f$ \delta n \f$ are dropped. A similar substitution is used for the conduction term in the temperature equation.
The equilibrium must therefore be a steady state of the ideal and thermal parts of the equations.

\note These substitutions are exact for a static (\f$ V_0 = 0 \f$), uniform-temperature equilibrium in force
balance (\f$ J_\mu \times B = \mu_0 \nabla p \f$ at \f$ u_0 \f$).
