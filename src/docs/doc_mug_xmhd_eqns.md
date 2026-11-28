MUG: Extended MHD equations (H(Curl) + Lagrange implementation)     {#doc_mug_xmhd_eqns}
================

[TOC]

This page documents the system of equations advanced by the nonlinear time integration in the
primary MUG implementation (\ref xmhd, driven by \ref xmhd::xmhd_run "xmhd_run"). The Lagrange-only
variant (\ref xmhd_lag) advances the same system, with the differences listed in
\ref doc_mug_xmhd_eqns_lag. The complete
system is defined by the nonlinear residual function \ref xmhd::xmhd_errmatrix_apply
"xmhd_errmatrix_apply". The Jacobian built in \ref xmhd::xmhd_build_ops "xmhd_build_ops" is used only
for preconditioning, so it contains only an approximate subset of the terms below and is not
described here. The linearized model (\ref xmhd::xmhd_lin_run "xmhd_lin_run") is also not covered.

Runtime options referenced below are set in the `xmhd_options` group of the input file (see
\ref xmhd::xmhd_read_settings "xmhd_read_settings"), unless noted otherwise.

\section doc_mug_xmhd_eqns_conv Fields, units, and notation

| Field | Basis | Units | Notes |
|-------|-------|-------|-------|
| \f$ B \f$   | H(Curl) + \f$ \nabla H^1 \f$ | T | |
| \f$ V \f$   | Vector Lagrange | m/s | Stored scaled by `vel_scale` |
| \f$ n \f$   | Lagrange | m^-3 | Stored scaled by `den_scale` |
| \f$ T_i \f$ | Lagrange | eV | Ion temperature |
| \f$ T_e \f$ | Lagrange | eV | Electron temperature (two-temperature model only, enabled by passing `Te` in the initial fields to \ref xmhd::xmhd_run "xmhd_run") |
| \f$ N_2 \f$ | Lagrange | m^-5 | Hyper-diffusion auxiliary field (`d2_dens>0` only) |
| \f$ J_2 \f$ | H(Curl) | T/m^2 | Hyper-resistivity auxiliary field (`eta_hyper>0` only) |

The following conventions are used throughout:
 - \f$ J_\mu \equiv \nabla \times B \f$, so the physical current density is \f$ J = J_\mu / \mu_0 \f$.
 - Temperatures are in eV and \f$ k = e \f$ (`k_boltz = elec_charge`).
 - \f$ m_i = \mu_{ion} m_p \f$ (`mu_ion`) and \f$ m_e' = f_{m_e} m_e \f$, where \f$ f_{m_e} \f$ is the artificial electron
   mass enhancement `me_factor`.
 - \f$ n_0 \f$ is the density scale `den_scale`, which defaults to the RMS value of the initial density.
 - Resistivities, viscosities, and diffusivities are given in diffusivity units (m^2/s). In particular,
   \f$ \eta \f$ is the physical resistivity divided by \f$ \mu_0 \f$.
 - \f$ \hat{b} = B / \left( |B| + \epsilon_B \right) \f$, where \f$ \epsilon_B \f$ is `xmhd_B_eps`.
 - \f$ (\nabla V)_{ij} = \partial_i V_j \f$, \f$ A:C = \sum_{ij} A_{ij} C_{ij} \f$, and
   \f$ (\nabla \cdot A)_j = \sum_i \partial_i A_{ij} \f$.
 - \f$ W = \nabla V + (\nabla V)^T - \frac{2}{3} I \nabla \cdot V \f$ is the rate of strain tensor.
 - \f$ V^n \f$, \f$ n^n \f$, and \f$ J_\mu^n \f$ denote values at the start of the current time step.
 - \f$ \textbf{r} \f$ is the position vector, which is the argument of the optional spatial shaping functions.

The ion temperature \f$ T_i \f$ is always evolved. In the single-temperature model the electron temperature is
not evolved separately, but is a fixed multiple of the ion temperature, \f$ T_e = f_{T_e} T_i \f$, where
\f$ f_{T_e} \f$ is `te_factor`. The two-temperature model, which evolves \f$ T_e \f$ separately, is used when an
initial electron temperature is supplied to \ref xmhd::xmhd_run "xmhd_run".

The momentum and temperature equations are evolved in non-conservative form. Each is divided by
\f$ m_i n \f$ (momentum) or \f$ n k / (\gamma - 1) \f$ (temperature), so the time derivative of every field
has an unweighted mass matrix.

\section doc_mug_xmhd_eqns_cont Continuity

\f[ \frac{\partial n}{\partial t} + \nabla \cdot \left( n V \right) =
\nabla \cdot \left( D \nabla n \right) - \nabla \cdot \left( D_h \nabla N_2 \right)
+ \nabla \cdot \left[ w_n V^n \left( V \cdot \nabla n \right) \right] + S_n \f]

\f[ N_2 = \nabla^2 n \f]

 - \f$ D \f$ = `d_dens` \f$ \times \f$ `diff_shaping`\f$ (\textbf{r}) \f$
 - \f$ D_h \f$ = `d2_dens`. The hyper-diffusion term and the \f$ N_2 \f$ equation are included only if
   `d2_dens>0`. The \f$ N_2 \f$ equation is solved in weak form
   (\f$ \int \phi N_2 + \int \nabla \phi \cdot \nabla n = 0 \f$) with natural boundary conditions, so for uniform
   \f$ D_h \f$ the term is \f$ -D_h \nabla^4 n \f$.
 - The upwinding term (\f$ w_n \f$) is included if `xmhd_upwind=T`. See \ref doc_mug_xmhd_eqns_upwind.
 - \f$ S_n \f$ is the density floor source. See \ref doc_mug_xmhd_eqns_floors.
 - The equation is advanced only if `xmhd_adv_den=T`.

\section doc_mug_xmhd_eqns_mom Momentum

\f[ \frac{\partial V}{\partial t} + V \cdot \nabla V = \frac{J_\mu \times B}{\mu_0 m_i n}
- \frac{\nabla \left[ n k \left( T_i + T_e \right) \right]}{m_i n} - g \hat{z}
+ \frac{1}{n} \nabla \cdot \left( n \mathcal{V} \right)
+ \nabla \cdot \left[ w_V V^n \left( V \cdot \nabla V \right) \right] \f]

where \f$ \mathcal{V} = -\Pi / (m_i n) \f$ is the viscous stress per unit mass, and the last term is the
upwinding contribution, in which \f$ V^n \left( V \cdot \nabla V \right) \f$ is a dyad (see
\ref doc_mug_xmhd_eqns_upwind).

 - The advection term is included if `xmhd_advec=T` and the Lorentz force if `xmhd_jcb=T`. `xmhd_jcb`
   is forced to `F` if the magnetic field is not advanced.
 - The pressure gradient is split into \f$ k (T_i + T_e) \nabla n \f$, which is included only if
   `xmhd_adv_den=T`, and \f$ n k \nabla (T_i + T_e) \f$, which is included only if `xmhd_adv_temp=T`.
 - \f$ g \f$ = `g_accel`, so gravity acts in the \f$ -\hat{z} \f$ direction for `g_accel>0`. Gravity is disabled
   by default (`g_accel=0`). This is a public module variable, not a namelist option.
 - The upwinding term (\f$ w_V \f$) is included only if `xmhd_upwind=T` and `xmhd_advec=T`.

\subsection doc_mug_xmhd_eqns_visc Viscosity

The viscous stress is written as \f$ \frac{1}{n} \nabla \cdot \left( n \mathcal{V} \right) \f$, where
\f$ \mathcal{V} = -\Pi / (m_i n) \f$ and \f$ \Pi = \Pi_f + \Pi_\parallel + \Pi_{gv} \f$. Each part is included
independently and the parts are summed. In the weak form the divergence is integrated by parts onto the
test function, so the residual (written as in \ref doc_mug_xmhd_eqns_time) contains
\f$ \int \partial_i \phi \, \mathcal{V}_{ik} \f$ for velocity component \f$ k \f$. The boundary term is dropped, which
gives a free-slip condition wherever the velocity is not constrained.

The fluid and parallel parts use constant dynamic viscosities \f$ \mu = m_i n_0 \nu \f$, which is why
they carry the factor \f$ n_0 / n \f$ below. The coefficients are
 - \f$ \nu_f \f$ = `nu_visc` \f$ \times \f$ `visc_shaping`\f$ (\textbf{r}) \f$
 - \f$ \nu_\parallel \f$ = `nu_par` \f$ \times \f$ `visc_shaping`\f$ (\textbf{r}) \f$

`nu_visc` must be set in the input file (use `nu_visc=0` for no fluid viscosity). The legacy name `nu_perp` is
still accepted, with a deprecation warning.

**Fluid viscosity** (\f$ \Pi_f \f$)

The fluid viscosity does not depend on the magnetic field direction. Its form is set by `visc_type`.

With `visc_type='kin'`
\f[ \mathcal{V}_f = \nu_f \frac{n_0}{n} \nabla V \f]

With `visc_type='iso'`
\f[ \mathcal{V}_f = \nu_f \frac{n_0}{n} W, \quad \text{i.e.} \quad \Pi_f = -\mu_f W \f]

\note Anisotropic viscosity is obtained by setting `nu_par` in addition to `nu_visc` with either type.

**Parallel viscosity** (\f$ \Pi_\parallel \f$)

If `nu_par>0`
\f[ \mathcal{V}_\parallel = \frac{\nu_\parallel}{2} \frac{n_0}{n} \left( \hat{b} \cdot W \cdot \hat{b} \right)
\left( 3 \hat{b}\hat{b} - I \right), \quad \text{i.e.} \quad
\Pi_\parallel = -\frac{3}{2} \mu_\parallel \left( \hat{b} \cdot W \cdot \hat{b} \right) \left( \hat{b}\hat{b} - \frac{1}{3} I \right) \f]

This is Braginskii's parallel stress, with \f$ \mu_\parallel = m_i n_0 \nu_\parallel \f$ playing the role of
\f$ \eta_0 \f$ (\f$ \Pi_{zz} = -\mu_\parallel W_{zz} \f$ for \f$ \hat{b} = \hat{z} \f$). It adds to the fluid
stress rather than replacing its parallel part, so with `visc_type='iso'` the total parallel
coefficient is \f$ \mu_\parallel + \mu_f \f$.

**Gyroviscosity** (\f$ \Pi_{gv} \f$)

If `nu_gyro>0`
\f[ \mathcal{V}_{gv} = -\frac{\nu_g k T_i}{4 e |B|} \left[ \hat{b} \times W \cdot \left( I + 3 \hat{b}\hat{b} \right)
- \left( I + 3 \hat{b}\hat{b} \right) \cdot W \times \hat{b} \right] \f]

where \f$ \nu_g \f$ = `nu_gyro`, \f$ (\hat{b} \times A)_{ij} = \epsilon_{ikl} \hat{b}_k A_{lj} \f$, and
\f$ (A \times \hat{b})_{ij} = \epsilon_{jkl} A_{ik} \hat{b}_l \f$. This is Braginskii's ion gyroviscous stress,

\f[ \Pi_{gv} = \nu_g \frac{p_i}{4 \Omega_i} \left[ \hat{b} \times W \cdot \left( I + 3 \hat{b}\hat{b} \right)
- \left( I + 3 \hat{b}\hat{b} \right) \cdot W \times \hat{b} \right], \quad
p_i = n k T_i, \quad \Omega_i = \frac{e |B|}{m_i} \f]

with the same sign as Braginskii for a positive ion. The stress is symmetric, and it does
no work (\f$ \Pi_{gv} : \nabla V = 0 \f$). The ion charge is taken to be \f$ e \f$. This is the same form of the
gyroviscous stress used in NIMROD (C. R. Sovinec and J. R. King, J. Comput. Phys. 229, 5803 (2010), Eq. 7).

\warning Gyroviscosity is not dissipative, but it is grouped with the dissipative terms in the time
discretization (see \ref doc_mug_xmhd_eqns_time). With the default `xmhd_diss_centered=F` it is therefore
treated with backward Euler, which adds numerical damping to the waves it supports.

\note The \f$ \nabla n \cdot \mathcal{V} / n \f$ part of the viscous term is included only if `xmhd_adv_den=T`.
Otherwise the term reduces to \f$ \nabla \cdot \mathcal{V} \f$.

\section doc_mug_xmhd_eqns_ind Induction and Ohm's law

\f[ \frac{\partial B}{\partial t} = - \nabla \times E \f]

\f[ E = -V \times B + \eta J_\mu - \eta_h \nabla \times J_2
+ \frac{1}{e n} \left[ \frac{J_\mu \times B}{\mu_0} - \nabla \left( n k T_e \right) \right]
+ \frac{m_e'}{\mu_0 e^2 n} \frac{\partial J_\mu}{\partial t} \f]

\f[ J_2 = -\nabla \times J_\mu \f]

 - The Hall, electron pressure, and electron inertia terms are included only if `xmhd_hall=T`.
 - Hyper-resistivity is included only if `eta_hyper>0`, with \f$ \eta_h \f$ = `eta_hyper`. The \f$ J_2 \f$ equation
   is solved in weak form (\f$ \int \psi \cdot J_2 + \int \nabla \times \psi \cdot J_\mu = 0 \f$) with natural boundary
   conditions, so the term is \f$ \eta_h \nabla \times \nabla \times J_\mu = -\eta_h \nabla^2 J_\mu \f$.
 - The \f$ \nabla H^1 \f$ component of \f$ B \f$ has only a mass term
   (\f$ \int \nabla \chi \cdot \partial B / \partial t = 0 \f$), because \f$ \nabla \times \nabla \chi = 0 \f$.
 - The equation is advanced only if `xmhd_adv_b=T`. If it is not, `xmhd_jcb=F` and `eta_hyper<0` are forced.

**Resistivity**

\f[ \eta = \begin{cases}
\eta_0 \, \eta_{reg} \, f_\eta(\textbf{r}) & \text{if } T_{\eta} \leq 0 \\
\eta_0 \left( T_{\eta} / T_e \right)^{3/2} f_\eta(\textbf{r}) & \text{if } T_{\eta} > 0
\end{cases} \f]

where \f$ \eta_0 \f$ = `eta`, \f$ T_\eta \f$ = `eta_temp`, \f$ \eta_{reg} \f$ is the per-region `eta` value from the XML
`<region>` blocks (1 if not specified), and \f$ f_\eta \f$ is the optional `res_shaping` function. If the
`res_profile` function is set, \f$ \eta_0 \f$ `res_profile`\f$ (T_e) \f$ is used in place of
\f$ \eta_0 \left( T_{\eta} / T_e \right)^{3/2} \f$.

\warning The per-region scale factor \f$ \eta_{reg} \f$ is dropped when the temperature-dependent
resistivity is used (`eta_temp>0`).

\subsection doc_mug_xmhd_eqns_solid Solid (resistive wall) regions

XML `<region>` blocks with `type=2` (the default if `type` is omitted) are solid conducting regions.
In these cells only resistive diffusion of \f$ B \f$ is solved:

\f[ \frac{\partial B}{\partial t} = -\nabla \times \left( \eta_w J_\mu \right),
\quad \eta_w = \eta_{reg} \, f_\eta(\textbf{r}) \f]

Here \f$ \eta_{reg} \f$ is the absolute wall resistivity, not a scale factor on `eta`. The velocity, density,
and temperature equations are not assembled in solid cells. The velocity is held fixed at all nodes in solid
regions for any `vbc`. Nodes in solid regions are treated as boundary nodes for density and temperature when
`nbc='d'` and `tbc='d'`, respectively (\ref xmhd_lag holds them fixed for any `nbc` and `tbc`, see
\ref doc_mug_xmhd_eqns_lag).

\section doc_mug_xmhd_eqns_temp Temperature

The thermal conduction tensor for species \f$ s \f$ is
\f$ \chi_s = \chi_{\parallel,s} \hat{b}\hat{b} + \chi_{\perp,s} \left( I - \hat{b}\hat{b} \right) \f$.

\subsection doc_mug_xmhd_eqns_temp1 Single temperature

\f[ \begin{aligned}
\frac{\partial T_i}{\partial t} + V \cdot \nabla T_i + (\gamma - 1) T_i \nabla \cdot V = {}&
\frac{\gamma - 1}{n} \nabla \cdot \left( n \chi \cdot \nabla T_i \right)
+ \frac{\gamma - 1}{(1 + f_{T_e}) n k} \left( Q_\eta + Q_\nu \right) \\
&+ \nabla \cdot \left[ w_T V^n \left( V \cdot \nabla T_i \right) \right] + S_T
\end{aligned} \f]

In energy form this is
\f$ \frac{n k}{\gamma - 1} \frac{d T_i}{d t} = -n k T_i \nabla \cdot V - \nabla \cdot q + \frac{Q_\eta + Q_\nu}{1 + f_{T_e}} \f$,
with \f$ q = -n k \chi \cdot \nabla T_i \f$. Equivalently, the total thermal energy \f$ (1 + f_{T_e}) n k T_i / (\gamma - 1) \f$
receives all of the heating, and the conductive flux is \f$ (1 + f_{T_e}) q \f$.

\subsection doc_mug_xmhd_eqns_temp2 Two temperatures

\f[ \begin{aligned}
\frac{\partial T_i}{\partial t} + V \cdot \nabla T_i + (\gamma - 1) T_i \nabla \cdot V = {}&
\frac{\gamma - 1}{n} \nabla \cdot \left( n \chi_i \cdot \nabla T_i \right)
+ \frac{\gamma - 1}{n k} Q_\nu - \frac{T_i - T_e}{\tau_{eq}} \\
&+ \nabla \cdot \left[ w_i V^n \left( V \cdot \nabla T_i \right) \right] + S_{T_i}
\end{aligned} \f]

\f[ \begin{aligned}
\frac{\partial T_e}{\partial t} + V_e \cdot \nabla T_e + (\gamma - 1) T_e \nabla \cdot V_e = {}&
\frac{\gamma - 1}{n} \nabla \cdot \left( n \chi_e \cdot \nabla T_e \right)
+ \frac{\gamma - 1}{n k} Q_\eta + \frac{T_i - T_e}{\tau_{eq}} \\
&+ \nabla \cdot \left[ w_e V_e^n \left( V_e \cdot \nabla T_e \right) \right] + S_{T_e}
\end{aligned} \f]

\f[ V_e = V - \frac{J_\mu}{\mu_0 e n}, \quad
\nabla \cdot V_e = \nabla \cdot V + \frac{J_\mu \cdot \nabla n}{\mu_0 e n^2} \f]

 - Thermal equilibration is included only if `xmhd_therm_equil=T`, with
   \f$ \tau_{eq} = 3.125 \times 10^{14} \mu_{ion} T_e^{3/2} / (n \ln \Lambda) \f$ and
   \f$ \ln \Lambda = 24 - \ln \left( \sqrt{n / 10^6} / T_e \right) \f$ (NRL formulary).
 - The temperature equations are advanced only if `xmhd_adv_temp=T`.

\subsection doc_mug_xmhd_eqns_heat Heating

 - Ohmic heating, if `xmhd_ohmic=T`:
   - \f$ Q_\eta = \eta |J_\mu|^2 / \mu_0 \f$, which is the physical \f$ \eta_{phys} J^2 \f$.
 - Viscous heating, if `xmhd_visc_heat=T`:
   - \f$ Q_\nu = -\Pi : \nabla V = Q_f + Q_\parallel \f$
   - \f$ Q_f = \mu_f \nabla V : \nabla V \f$ (`kin`) or \f$ Q_f = \mu_f \nabla V : W \f$ (`iso`)
   - \f$ Q_\parallel = \frac{3}{4} \mu_\parallel \left( \hat{b} \cdot W \cdot \hat{b} \right)^2 \f$ if `nu_par>0`
   - Here \f$ \mu = m_i n_0 \nu \f$ as in \ref doc_mug_xmhd_eqns_visc. Gyroviscosity does no work, so it has no
     heating term.

\note The viscous heating is always time-centered, while the viscous stress is in \f$ R_d \f$ (see
\ref doc_mug_xmhd_eqns_time). With the default `xmhd_diss_centered=F` the two are evaluated at different
time levels, so the kinetic energy removed by viscosity and the heat added are equal only to
\f$ O(\Delta t) \f$.

\subsection doc_mug_xmhd_eqns_cond Thermal conduction coefficients

Below, \f$ \kappa_{\parallel,s} \f$ and \f$ \kappa_{\perp,s} \f$ are `kappa_par(s)` and `kappa_perp(s)`, with
\f$ s = 1 \f$ for ions (or the single temperature) and \f$ s = 2 \f$ for electrons. All coefficients are also
multiplied by the optional `cond_shaping`\f$ (\textbf{r}) \f$ function.

With `xmhd_brag=F` the conductivities are constant:
\f[ \chi_{\parallel,s} = \kappa_{\parallel,s} \frac{n_0}{n}, \quad \chi_{\perp,s} = \kappa_{\perp,s} \frac{n_0}{n} \f]
so that \f$ n k \chi_s = n_0 k \kappa_s \f$.

With `xmhd_brag=T`, \f$ \kappa \f$ are dimensionless multipliers on the Braginskii diffusivities
\f$ \chi^{B} \f$ (\ref mhd_utils, with \f$ \ln \Lambda = 11 \f$ and the physical electron mass). In each
species, \f$ \chi^{B}_\perp \f$ is capped at \f$ \chi^{B}_\parallel \f$.
 - Two temperatures:
   - \f$ \chi_{\parallel,i} = \kappa_{\parallel,1} \chi^{B}_{\parallel,i}(T_i) \f$
   - \f$ \chi_{\perp,i} = \kappa_{\perp,1} \chi^{B}_{\perp,i}(T_i) \f$
   - \f$ \chi_{\parallel,e} = \kappa_{\parallel,2} \chi^{B}_{\parallel,e}(T_e) \f$
   - \f$ \chi_{\perp,e} = \max \left( \chi_{e,min}, \kappa_{\perp,2} \chi^{B}_{\perp,e}(T_e) \right) \f$
 - Single temperature:
   - \f$ \chi_{\parallel} = \kappa_{\parallel,1} \chi^{B}_{\parallel,e}(T_i) \f$
   - \f$ \chi_{\perp} = \kappa_{\perp,1} \max \left( \chi^{B}_{\perp,i}(T_i), \chi^{B}_{\perp,e}(T_i) \right) \f$

\note \f$ \chi_{e,min} \f$ is the compile-time constant `chi_perp_e_min` in \ref xmhd (m^2/s), which is 0 by default
(no limit). If needed, a lower limit on the electron perpendicular diffusivity can be set by changing it.

\section doc_mug_xmhd_eqns_upwind Upwinding

If `xmhd_upwind=T`, streamline diffusion terms of the form
\f$ \nabla \cdot \left[ w \, v_{ad} \left( v \cdot \nabla f \right) \right] \f$ are added to the advected fields,
where \f$ v_{ad} \f$ is the advecting velocity at the start of the time step. This is an inconsistent
streamline-upwind Petrov-Galerkin (SUPG) method: the modified test functions
\f$ \phi + w \, v_{ad} \cdot \nabla \phi \f$ are applied to the advection terms only, so the scheme adds
diffusion along the streamlines. The weight is the optimal value for 1D advection-diffusion from
A. N. Brooks and T. J. R. Hughes, Comput. Methods Appl. Mech. Eng. 32, 199 (1982):

\f[ w = \frac{\beta h_e}{2 |v_{ad}|}, \quad \beta = \coth(Pe) - \frac{1}{Pe}, \quad
Pe = \frac{|v_{ad}| h_e}{2 \nu_{up}} \f]

The element size \f$ h_e \f$ is the node spacing of a regular cell with the same volume \f$ V_{cell} \f$,
\f$ h_e = (6 \sqrt{2} V_{cell})^{1/3} / p \f$ for tetrahedra and \f$ h_e = V_{cell}^{1/3} / p \f$ for hexahedra,
where \f$ p \f$ is the Lagrange element order.

| Equation | \f$ f \f$ | \f$ v_{ad} \f$ | \f$ v \f$ | \f$ \nu_{up} \f$ |
|----------|-----------|----------------|-----------|------------------|
| Continuity | \f$ n \f$ | \f$ V^n \f$ | \f$ V \f$ | \f$ D \f$ |
| Momentum | \f$ V \f$ | \f$ V^n \f$ | \f$ V \f$ | \f$ \nu_f n_0 / n \f$ |
| Ion temperature | \f$ T_i \f$ | \f$ V^n \f$ | \f$ V \f$ | \f$ (\gamma - 1) \chi_{\perp} \f$ or \f$ (\gamma - 1) \chi_{\perp,i} \f$ |
| Electron temperature | \f$ T_e \f$ | \f$ V_e^n \f$ | \f$ V_e \f$ | \f$ (\gamma - 1) \chi_{\perp,e} \f$ |

 - Each term is included only when its equation is advanced, and the momentum term only if `xmhd_advec=T`.
 - The diffusivities of the upwinded fields (`nu_visc`, `d_dens`, and `kappa_perp`) must be positive, and the run
   aborts otherwise. Where a shaping function makes \f$ \nu_{up} \f$ zero locally, full upwinding
   (\f$ \beta = 1 \f$) is used.
 - Upwinding is always disabled in linearized runs (\ref xmhd::xmhd_lin_run "xmhd_lin_run").

\section doc_mug_xmhd_eqns_floors Density and temperature floors

If `den_floor>0` (`temp_floor>0`), the density (temperatures) used to evaluate all coefficients is replaced
by \f$ \max(n, n_{floor}) \f$ (\f$ \max(T, T_{floor}) \f$). The unmodified value is still used in the time derivative.

Cells where a field falls below its floor also receive a source. At each quadrature point the
relative deficit \f$ (f_{floor} - f) / f_{floor} \f$ is computed, and its maximum over the cell, \f$ d_c \f$, is recorded.
At the start of each step the cell source amplitude \f$ s_c \f$ is first halved, and set to zero if the result is below 0.01.
Then, in cells where \f$ d_c > 0 \f$, it is updated as \f$ s_c \leftarrow \max \left( s_c, \min(2, d_c) \right) \f$. The second
update is repeated, without halving, after a failed nonlinear solve. Within cell \f$ c \f$ this gives

\f[ S_f = \frac{s_c f_{floor}}{2 \Delta t} \f]

that is, an increase of \f$ s_c f_{floor} / 2 \f$ per step.

\section doc_mug_xmhd_eqns_time Time discretization

Write the semi-discrete system as \f$ M \partial u / \partial t + R_c(u) + R_d(u) = 0 \f$, where \f$ R_d \f$ contains
the dissipative terms listed below and \f$ R_c \f$ contains everything else. Each time step solves

\f[ M \left( u^{n+1} - u^n \right) + \frac{\Delta t}{2} \left[ R_c(u^{n+1}) + R_c(u^n) \right]
+ \Delta t \left[ \theta R_d(u^{n+1}) + (1 - \theta) R_d(u^n) \right] = 0 \f]

with Newton's method. Here \f$ \theta = 1 \f$ (backward Euler) by default and \f$ \theta = 1/2 \f$ (Crank-Nicolson) if
`xmhd_diss_centered=T`. In the implementation, the right-hand side is formed by evaluating the
residual with a time step of \f$ -\Delta t/2 \f$, which omits \f$ R_d \f$ unless `xmhd_diss_centered=T`. The
nonlinear solve then uses \f$ +\Delta t/2 \f$, with \f$ R_d \f$ doubled unless `xmhd_diss_centered=T`.

Terms in \f$ R_d \f$:
 - Hyper-resistivity
 - Viscosity (isotropic, parallel, and gyroviscosity) and velocity upwinding. Gyroviscosity is not
   dissipative but is grouped here (see \ref doc_mug_xmhd_eqns_visc).
 - Thermal conduction and temperature upwinding
 - Density hyper-diffusion

All other terms, including resistivity, density diffusion, density upwinding, and Ohmic and viscous
heating, are always time-centered.

Additional details:
 - The upwinding velocities (\f$ V^n \f$, \f$ V_e^n \f$) are taken from the start of the step.
 - The electron inertia term is discretized as
   \f$ \frac{m_e'}{\mu_0 e^2} \frac{1}{2} \left( \frac{1}{n^{n+1}} + \frac{1}{n^n} \right) \frac{J_\mu^{n+1} - J_\mu^n}{\Delta t} \f$.
 - The floor sources \f$ S_f \f$ are applied as described in \ref doc_mug_xmhd_eqns_floors.
 - Rows of the auxiliary equations (\f$ N_2 \f$, \f$ J_2 \f$) are scaled by `n2_scale` and `j2_scale`. They are
   assembled only for the implicit solve, and the right-hand side uses the auxiliary fields from the
   previous step.

\section doc_mug_xmhd_eqns_lag Lagrange-only variant

The Lagrange-only variant (\ref xmhd_lag) solves the system above, with the same options, time
discretization, and boundary conditions, but represents the magnetic field using vector Lagrange
elements (the same basis as \f$ V \f$). The differences are:
 - As the nodal magnetic field is not divergence free by construction, a divergence cleaning term is
   added to the induction equation (including solid regions)
   \f[ \frac{\partial B}{\partial t} = - \nabla \times E + \kappa_{db} \nabla \nabla \cdot B \f]
   where \f$ \kappa_{db} \f$ is set by `k_db` (m^2/s, default 1). The term is integrated by parts with the
   boundary term dropped and is always time-centered. The periodic divergence cleaning of \ref xmhd is not
   used (`nclean` is ignored), and the divergence error \f$ \int (\nabla \cdot B)^2 dV \f$ is always computed
   and reported (`xmhd_monitor_div` is ignored).
 - The magnetic field boundary condition (`bbc`) is applied to the nodal field: `bbc='bc'` holds the normal
   component of \f$ B \f$ fixed, and `bbc='ic'` holds the tangential components fixed, on the boundary. On the
   axis of rotationally periodic (revolved) meshes the components normal to the axis are held fixed, as for
   \f$ V \f$.
 - Density and temperature are held fixed at all nodes in solid regions for any `nbc` and `tbc`, as for the
   velocity.
 - The hyper-resistivity terms (`eta_hyper>0`), which use an auxiliary H(Curl) field \f$ J_2 \f$ with the
   Lagrange magnetic field, have not been verified against \ref xmhd.
