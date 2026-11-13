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

Additional verification cases are in the regression tests: `src/tests/physics/test_alfven_2d.F90`
(Alfvén wave), `test_sound_2d.F90` (sound wave) and `test_taylor_green_2d.F90` (incompressible Taylor-Green vortex).

\section doc_mug2d_main_geom Geometry and field representation
The mesh is defined in the \f$ (x,z) \f$ plane for a slab, or the \f$ (R,Z) \f$ plane for cylindrical
geometry (`cyl_flag=T`). The first mesh coordinate is \f$ x \f$ or \f$ R \f$, and the second is \f$ Z \f$. Both
\f$ (x,y,z) \f$ and \f$ (R,\phi,Z) \f$ are right-handed.

The magnetic field is represented by a poloidal flux \f$ \psi \f$ and an out-of-plane component.

Slab:
\f[ \mathbf{B} = \nabla \psi \times \hat{y} + B_y \hat{y} + \mathbf{B}_0 \f]

Cylindrical, with \f$ F = R B_\phi \f$:
\f[ \mathbf{B} = \frac{1}{R} \nabla \psi \times \hat{\phi} + \frac{F}{R} \hat{\phi} + B_{0,Z} \hat{Z} \f]

where
 - \f$ \mathbf{B}_0 \f$ (`B_0`) is a uniform, static background field. In cylindrical geometry only
   \f$ B_{0,Z} \f$ is allowed, since \f$ B_{0,R} \hat{R} \f$ and \f$ B_{0,\phi} \hat{\phi} \f$ are not curl and
   divergence-free. A nonzero `B_0(1)` or `B_0(2)` aborts the run.
 - \f$ \psi = \tilde{\psi} + \psi_{vac} \f$. Here \f$ \tilde{\psi} \f$ is the evolved flux, and \f$ \psi_{vac} \f$
   (`psi_vac`) is an optional fixed vacuum flux, e.g. from external coils, that must satisfy
   \f$ \nabla^2 \psi_{vac} = 0 \f$ (slab) or \f$ \Delta^* \psi_{vac} = 0 \f$ (cylindrical). It contributes to
   \f$ \mathbf{B} \f$ and to advection, but not to the current or to resistive diffusion.

\section doc_mug2d_main_eqns MHD System
In all of the equations below:
 - \f$ \rho = m_i n \f$ and \f$ \mathbf{J} = \nabla \times \mathbf{B} / \mu_0 \f$.
 - \f$ \mathbf{B} \f$ is the total field, including \f$ \mathbf{B}_0 \f$ and \f$ \psi_{vac} \f$.
 - \f$ D \f$, \f$ \nu \f$, \f$ \chi \f$ and \f$ \eta \f$ are diffusivities in m\f$^2\f$/s.
 - \f$ T \f$ is in eV, \f$ k \f$ (`k_boltz`) converts it to Joules, and \f$ T_e = T_i = T \f$.

**Compressible** (`incomp=F`)

\f[ \frac{\partial n}{\partial t} + \nabla \cdot \left( n \mathbf{u} \right) = D \nabla^2 n \f]

\f[ \rho \left[ \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} \right]
= \mathbf{J} \times \mathbf{B} - \nabla \left( 2 n k T \right) + \nabla \cdot \left( \rho \nu \nabla \mathbf{u} \right) \f]

\f[ \frac{1}{\gamma-1} \left[ \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right] =
- T \nabla \cdot \mathbf{u} + \frac{\chi}{n} \nabla \cdot \left( n \nabla T \right) \f]

where \f$ \gamma \f$ is the adiabatic index.

\f[ \frac{\partial \mathbf{B}}{\partial t} = \nabla \times \left( \mathbf{u} \times \mathbf{B} \right)
- \nabla \times \left( \eta \nabla \times \mathbf{B} \right) \f]

**Incompressible** (`incomp=T`)

The density is uniform and held fixed, so \f$ \rho \f$ is constant. The temperature equation is
replaced by the incompressibility constraint, and the pressure \f$ p \f$ [Pa] is solved for in place of \f$ T \f$. The
pressure uses Lagrange elements one order lower than the other fields (Taylor-Hood pairing).

\f[ \rho \left[ \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} \right]
= \mathbf{J} \times \mathbf{B} - \nabla p + \rho \nu \nabla^2 \mathbf{u}, \qquad \nabla \cdot \mathbf{u} = 0 \f]

with the same induction equation as above.

**Viscosity**

Viscosity uses the kinematic "Laplacian" form \f$ \Pi = -\rho \nu \nabla \mathbf{u} \f$, which is equivalent to
`visc_type='kin'` in 3D MUG.

\subsection doc_mug2d_main_eqns_ind 2D induction equations
The induction equation is solved as two scalar equations, one for \f$ \psi \f$ and one for the out-of-plane field.
The resistive term uses `eta(:,2)` in the \f$ \psi \f$ equation, because it acts on the out-of-plane current. It uses
`eta(:,1)` in the \f$ B_y \f$/\f$ F \f$ equation, because it acts on the in-plane current.

Slab:
\f[ \frac{\partial \psi}{\partial t} + \mathbf{u} \cdot \nabla \psi
- \left( \mathbf{u} \times \mathbf{B}_0 \right) \cdot \hat{y} = \eta_2 \nabla^2 \tilde{\psi} \f]

\f[ \frac{\partial B_y}{\partial t} + \nabla \cdot \left( B_y \mathbf{u} \right) =
\left( \nabla u_y \times \nabla \psi \right) \cdot \hat{y} + \mathbf{B}_0 \cdot \nabla u_y
- B_{0,y} \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_z}{\partial z} \right)
+ \nabla \cdot \left( \eta_1 \nabla B_y \right) \f]

Cylindrical, with \f$ \Delta^* = R \frac{\partial}{\partial R} \left( \frac{1}{R} \frac{\partial}{\partial R} \right) + \frac{\partial^2}{\partial Z^2} \f$:
\f[ \frac{\partial \psi}{\partial t} + \mathbf{u} \cdot \nabla \psi + R u_R B_{0,Z} = \eta_2 \Delta^* \tilde{\psi} \f]

\f[ \frac{\partial F}{\partial t} = R^2 \left[ \mathbf{B}_p \cdot \nabla \left( \frac{u_\phi}{R} \right)
- \nabla \cdot \left( \frac{F \mathbf{u}}{R^2} \right) \right]
+ R \nabla \cdot \left( \frac{\eta_1}{R} \nabla F \right) \f]

where \f$ \mathbf{B}_p = \frac{1}{R} \nabla \psi \times \hat{\phi} + B_{0,Z} \hat{Z} \f$ is the total poloidal field.

\section doc_mug2d_main_params Run parameters
Simulation parameters are members of \ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim". Allocate the per-region
arrays with size `mesh%nreg` before calling `setup`, or assign them afterwards. Arrays left unset are filled with
the placeholder value -1, which is not physical. Set every coefficient used by your model.

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
| `B_0(3)` | 0 | Background field [T] (only `B_0(3)` if `cyl_flag=T`) |
| `den_scale` | \f$ 10^{19} \f$ | Density normalization |
| `k_boltz` | \f$ e \f$ | Converts \f$ T \f$ to J (temperature in eV) |
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
The time-advance operator is the Jacobian evaluated at `u0`. It is built once and each step is a single GMRES solve with a fixed time step. The equilibrium is written to `xmhd2d_eq.rst`.
Linear mode is currently Cartesian only.
