# Gyroviscous stabilization of the gravitational instability (planned)

**Status: not implemented.** This folder holds notes for a planned MUG example and
non-regression test of Braginskii gyroviscosity (`nu_gyro`). It is not yet added to
`src/examples/MUG/CMakeLists.txt`.

Reference: N. M. Ferraro and S. C. Jardin, "Finite element implementation of Braginskii's
gyroviscous stress with application to the gravitational instability", Phys. Plasmas
**13**, 092101 (2006), <https://doi.org/10.1063/1.2236277>. Equation and figure numbers
below refer to this paper.

## Why this case

The fast regression tests (`src/tests/physics/MUG/test_MUG_visc.py` and `test_MUG_gyro.py`) check the
viscous stresses in a uniform plasma. This case also exercises:

- the full nonlinear gyroviscous stress with gradients in density, pressure, and |B|,
  including the $\nabla n \cdot \mathcal{V} / n$ part of the momentum equation
- the strong, well-characterized stabilizing effect of gyroviscosity on a growth rate,
  even where $k \rho_i \ll 1$
- the interaction of gyroviscosity with the Hall term
- "gyroviscous cancellation" (Sec. III C), where the perturbation is advected by
  $-v_*$ with the ion magnetization velocity $v_*$ (Eq. 19), not the diamagnetic drift

## Physics setup (paper, Sec. III B)

- Isothermal plasma with $p = p_i$ (no electron pressure) in uniform gravity
  $\vec{g} = -g \hat{y}$, and the magnetic field $I \hat{z}$ along the ignorable
  direction.
- Equilibrium (Eq. 13): $n = n_0 e^{y/L_n}$, $p = T n$, and
  $I^2 = I_0^2 - 2 (g L_n + \Gamma T)(n - n_0)$. The field gradient, not only the
  pressure gradient, supports the equilibrium.
- Domain $L_x \times L_y$, periodic in x and conducting, no-slip walls at
  $y = \pm L_y/2$. Perturbation $\sim e^{i(k_x x + k_y y)}$ with $k_x = 2\pi/L_x$
  and $k_y = \pi/L_y$, amplitude $\epsilon \sim 10^{-6}$ in density.
- Linear runs use $\eta = \mu = 0$ and $\theta = 1/2$ (time-centered gyroviscosity).
- Growth rate from Eq. 16, $\gamma = \frac{1}{2} \partial_t \ln \int A^2(x, y=0, t) dx$,
  with A the density or stream-function perturbation, once it has become constant.
- Units: density $\bar{n}$, field $\bar{B}$, time $\Omega_i^{-1} = m_i/(e \bar{B})$,
  length $d_i = c/\omega_{pi}$, so velocities are in units of $v_A$.

Target cases (Fig. 4, with $n_0 = 1$):

| Case | $L_x$ | $L_y$ | $L_n$ | g | $I_0$ | Result |
|------|-------|-------|-------|---|-------|--------|
| Low-β | $2\pi/10$ | 1 | 100 | 0.005 | 100 | γ ≈ 0.0067 at β → 0, stable above β ≈ 0.0027 |
| High-β | $2\pi/0.05$ | $\pi/0.01$ | $10^5$ | 0.0005 | 1 | γ ≈ 5.3e-4 at β → 0, stable above β ≈ 16 |

β is varied through the temperature. Compare with the full dispersion relation (Eq. 14),
its $k \rho_i \ll 1$ limit (Eq. 15), and Roberts and Taylor (Eq. 11, low β only). The
dimensionless parameters are defined in Eqs. 12 and around Eq. 9:
$\bar{\nu} = \rho_i^2 \Omega_i k / 2 v_A$, with $\rho_i^2 = p_0^{(i)} / n_0 B^2$ in these units.

## Mapping to MUG

1. **Geometry.** MUG gravity acts only along z (`g_accel`), so rotate the setup: density
   gradient and gravity along z, field along y (or x), wavevector along x. Use a thin 3D slab,
   periodic in x and y and bounded in z (eg. the `cube` mesh with `ref_per=T,T,F` and one or
   two cells across the field direction). With field along y, the field-aligned direction
   is the ignorable one, as in the paper.
2. **Gravity sign.** The implemented acceleration is $+g\hat{z}$, despite the source
   comment (see the warning in `doc_mug_xmhd_eqns.md`). Downward gravity needs `g_accel<0`.
   `g_accel` is a public module variable, not a namelist option, so set it in the driver.
3. **Pressure model.** Single temperature with `te_factor=0`, so that $p = p_i$ and the
   Hall-term electron pressure vanishes. Check that no enabled option divides by $T_e$ (eg. `eta_temp`,
   `xmhd_brag`, or equilibration), or use a small `te_factor` instead.
4. **Equilibrium.** Re-derive force balance for MUG's equations,
   $\frac{d}{dz} \left( n k T + \frac{B^2}{2 \mu_0} \right) = -m_i n g$, rather than copying
   Eq. 13c. Direct force balance with $p = T n$ gives $(g L_n + T)$ in place of the paper's
   $(g L_n + \Gamma T)$, so check the paper's convention before use. The projected
   equilibrium is not exactly in discrete balance, so check the residual flow from an
   unperturbed run and keep it well below the perturbation.
5. **Dissipation.** Set `d_dens`, `kappa_par`, `kappa_perp`, `eta`, and `nu_perp` to zero or
   negligible values. The default `d_dens=1` would diffuse the equilibrium density profile.
6. **Gyroviscosity.** `nu_gyro=1` with `xmhd_diss_centered=T`. Gyroviscosity is grouped with
   the dissipative terms, so the default backward-Euler treatment would damp the mode.
7. **Hall term.** Run with `xmhd_hall=T` to compare with Eqs. 14 and 15, and with
   `xmhd_hall=F` for comparison with the no-Hall contours (Fig. 3c).
8. **Boundary conditions.** `vbc='all'` for no-slip walls and `bbc='bc'` for conducting walls.
9. **Diagnostics.** Fit the growth rate from the perturbed kinetic energy or the density
   perturbation (equilibrium subtracted) along the mid-plane, following Eq. 16. A reference run
   with `nu_gyro<0` and `xmhd_hall=F` should recover the ideal growth rate.

## Simpler cases from the same paper

These are cheaper and could come first:

- **Homogeneous normal modes** (Sec. III A, Eqs. 9 and 10, Figs. 1 and 2). Fast, Alfvén, and
  slow wave phase velocities in a uniform plasma with Hall and full gyroviscosity, for
  oblique $\vec{k}$. Parameters: $p_0 = 0.25$, $p_0^{(i)} = 0.125$, $B_{x0} = 0.8$,
  $I_0 = 0.6$, $\Gamma = 5/3$, $\mu = \eta = 10^{-3}$, $\epsilon = 10^{-4}$.
- **Compressional gyroviscous force** (Sec. III C, Eqs. 20 and 21, Fig. 7). A perturbation
  $U_1 = \epsilon \cos(k_x x)$ in a uniform plasma, where $\nabla \cdot \Pi = p_i \nabla(\nabla^2 U) / 2 I_0$
  is the only compressional force. Parameters: $p_0 = p_0^{(i)} = 0.1$, $I_0 = 2$,
  $k_x = 2\pi$, $\epsilon = 10^{-4}$, $\mu = \eta = 0$.
