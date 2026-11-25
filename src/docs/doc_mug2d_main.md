MUG 2D: Slab and axisymmetric linear/nonlinear MHD simulation package     {#doc_mug2d_main}
================

[TOC]

\image html mug_logo.png width=25%

MUG 2D (\ref xmhd_2d) solves the visco-resistive MHD equations in two dimensions, in either Cartesian slab
or cylindrical (axisymmetric) geometry. It uses scalar Lagrange finite elements on unstructured
triangular grids. The out-of-plane direction (\f$ y \f$ in a slab, \f$ \phi \f$ in cylindrical) is
ignorable (\f$ \partial / \partial y = 0 \f$). Both compressible and incompressible models are available. Each time step is advanced fully
implicitly with a Newton-Krylov solver. A linearized mode (Cartesian only) advances perturbations about a fixed equilibrium.

For 3D simulations see \ref doc_mhd_main "MUG".

\section doc_mug2d_main_ex MUG 2D Examples
The following examples illustrate how to use MUG 2D.

 - \subpage doc_mug2d_harris_ex
 - \subpage doc_mug2d_rt_ex

Additional verification cases are in the regression tests in `src/tests/physics/MUG/`: `test_alfven_2d.F90`
(Alfven wave), `test_sound_2d.F90` (sound wave), `test_taylor_green_2d.F90` (incompressible Taylor-Green vortex) and
`test_rt_2d.F90` (gravity and Rayleigh-Taylor growth rates).
`test_jacobian_2d.F90` checks that the approximate Jacobian matches the nonlinear residual.

\section doc_mug2d_main_eqns MHD System
MUG 2D solves single-temperature (\f$ T_e = T_i \f$) visco-resistive MHD: continuity with particle diffusion,
momentum with kinematic viscosity, temperature with isotropic thermal conduction, and induction with an
anisotropic resistivity. With `incomp=T` the density is held fixed and the temperature equation is replaced
by the incompressibility constraint, with the pressure solved for in its place. The magnetic field is represented
by a poloidal flux \f$ \psi \f$ and an out-of-plane component (\f$ B_y \f$, or \f$ F = R B_\phi \f$ in cylindrical
geometry), plus an optional fixed vacuum flux (`psi_vac`) and uniform background field (`B_0`).
The geometry, field representation and complete system are documented in:

 - \subpage doc_mug_xmhd_2d_eqns

\section doc_mug2d_main_params Run parameters
Simulation parameters are members of \ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim". Allocate the per-region
arrays with size `mesh%nreg` before calling `setup`, or assign them afterwards. Arrays left unset are filled with
the placeholder value -1, which is not physical. Set every coefficient used by your model:
`run_simulation` and `run_lin_simulation` abort if, in any region that is not ignored, `m_i` or `eta(:,2)` is not
positive, `nu` or `eta(:,1)` is negative, or (compressible model only) `D_diff` or `chi` is negative or `gamma` is
not greater than 1.

| Member | Default | Description |
|--------|---------|-------------|
| `cyl_flag` | `F` | Use cylindrical (axisymmetric) coordinates |
| `incomp` | `F` | Solve the incompressible system |
| `linear` | `F` | Linear run (use `run_lin_simulation`) |
| `timestep_cn` | `F` | Second-order time advance (CN nonlinear, BDF2 linear) |
| `mfnk` | `F` | Use matrix-free Newton-Krylov |
| `dt` | -1 | Initial (and maximum) time step [s] |
| `nsteps` | -1 | Number of time steps |
| `rst_freq` | 10 | Restart file frequency (steps) |
| `pre_freq` | 1 | Preconditioner update frequency with `mfnk=T` (steps) |
| `ittarget` | 40 | Target linear iterations per step for adaptive \f$ \Delta t \f$ |
| `lin_tol`, `nl_tol` | \f$ 10^{-8} \f$, \f$ 10^{-5} \f$ | Absolute linear and nonlinear tolerances |
| `B_0(3)` | 0 | Background field [T] (only `B_0(3)` if `cyl_flag=T`). Equivalent to `psi_vac` and the initial \f$ B_y \f$ on a non-periodic domain, see \ref doc_mug_xmhd_2d_eqns_conv |
| `den_scale` | \f$ 10^{19} \f$ | Density normalization |
| `g_accel` | 0 | Gravitational acceleration along \f$ -z \f$ (the second mesh coordinate) [m/s\f$^2\f$] |
| `nu(nreg)` | -1 | Kinematic viscosity \f$ \nu \f$ [m\f$^2\f$/s] |
| `chi(nreg)` | -1 | Thermal diffusivity \f$ \chi \f$ [m\f$^2\f$/s] |
| `D_diff(nreg)` | -1 | Particle diffusivity \f$ D \f$ [m\f$^2\f$/s] |
| `eta(nreg,2)` | -1 | Magnetic diffusivity [m\f$^2\f$/s]: column 1 for \f$ B_y \f$/\f$ F \f$, column 2 for \f$ \psi \f$ |
| `gamma(nreg)` | -1 | Adiabatic index |
| `m_i(nreg)` | \f$ m_p \f$ | Ion mass [kg] |
| `ignore_rmask(nreg)` | `F` | Skip these regions when assembling the equations |
| `psi_vac` | 0 | Vacuum flux \f$ \psi_{vac} \f$ (`oft_blagrange` vector) |

\section doc_mug2d_main_linear Linear mode
With `linear=T`, the equilibrium is set in `u0` and the initial perturbation in `u`, and `run_lin_simulation` is called.
The equilibrium is written to `xmhd2d_eq.rst`. Linear mode is currently Cartesian only. The linearized operator and
time advance are described in \ref doc_mug_xmhd_2d_eqns_linear.
