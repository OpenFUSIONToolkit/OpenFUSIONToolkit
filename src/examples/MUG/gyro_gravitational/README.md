# Gyroviscous stabilization of the gravitational instability

MUG example reproducing the low-β case of N. M. Ferraro and S. C. Jardin, "Finite element
implementation of Braginskii's gyroviscous stress with application to the gravitational instability",
Phys. Plasmas **13**, 092101 (2006), <https://doi.org/10.1063/1.2236277>. The setup, units, and
results are documented in the example page (`MUG_gyro_grav.F90`, `doc_mug_gyro_grav_ex`).

Run the scan of growth rate against β from a directory containing the built `MUG_gyro_grav`,
`oft.in`, and `oft_in.xml`:

```bash
python run_scan.py --plot
```

The growth rate is fit over the exponential growth phase of the kinetic energy: the longest part of the
run where d ln(E_k)/dt is constant to within ±10% (`--fit-tol`). This excludes the initial transient and
the late blow-up or saturation, and the window used is printed for each case. `--fit-frac 0.3` restores the
previous fit over the last 30% of each run, and `--no-run` refits the results of a previous scan. The
perturbation amplitude can be changed with `--eps` (default: `eps` in `oft.in`).

## Lagrange variant

`MUG_gyro_grav_lag` (`MUG_gyro_grav_lag.F90`, `doc_mug_gyro_grav_lag_ex`) runs the same case with the
Lagrange-only xMHD module (`xmhd_lag`):

```bash
python run_scan.py --exe ./MUG_gyro_grav_lag --plot
```

The initial magnetic field must be projected onto the vector Lagrange space with a relative tolerance. The
plasma is supported almost entirely by magnetic pressure, so a relative error of 1e-7 in B is already a
force imbalance much larger than the seeded perturbation (`eps=1e-8`). The right-hand side of that
projection scales with the cell volume (about 1e-6 here), so an absolute tolerance of 1e-12 stops the solve
at that level of error. That gave an initial kinetic energy of about 1e-15, independent of `eps`, which
buried the gravitational mode for β ≥ 1e-3. The example now converges the projection to a relative
residual of 1e-14, as does the H(Curl) example (`MUG_gyro_grav.F90`). With `eps=0` the kinetic energy after
the first step is then about 3e-24 for the Lagrange variant and 2e-26 for H(Curl).

## Open items

- **Hall term.** With `xmhd_hall=T` a faster mode (γ ≈ 0.045 Ω_i) overtakes the gravitational mode.
  It also appears without gravity and does not depend on the time step or element order. Until it is
  understood, the example runs without the Hall term.
- **Gyroviscous mode in a non-uniform equilibrium.** Without the Hall term, a second mode grows at
  γ ≈ 0.007 Ω_i (4×1×4 mesh) for every β. At 8×1×8 it grows at γ ≈ 0.019 Ω_i, above the largest possible
  buoyancy growth rate √(g/Ln) ≈ 0.0071 Ω_i. It persists with g = 0, needs gyroviscosity (with `nu_gyro<0`
  the growth stays at the ideal rate), and disappears for a homogeneous plasma. It grows from round-off: with
  the initial field projection converged and `eps=0` (kinetic energy after the first step ≈ 2e-26) it reaches
  the same amplitude as with `eps=1e-8`, so a better initial force balance does not suppress it. It overtakes
  the gravitational mode late in runs with β ≳ 2.5e-3 (at β = 3e-3, stable by Eq. 15, the fit returns
  γ ≈ 0.0067). Separately, some runs fail late: at β = 1e-3 the time step collapses once the kinetic energy
  reaches about 1e-9, with either projection tolerance. The Hall-term mode above may share its cause.
- **Linearized model.** `xmhd_lin_run` cannot be used for this case: its operator drops the density
  dependence of J×B/n, which provides the buoyancy drive when the equilibrium is supported by a field
  gradient.
- **High-β case** (paper, Fig. 4 right: Lx = 2π/0.05, Ly = π/0.01, Ln = 1e5, g = 5e-4, I0 = 1), where
  compressibility and electromagnetic effects matter.
- **Simpler cases from the same paper:** homogeneous normal modes with Hall and gyroviscosity
  (Sec. III A, Figs. 1 and 2) and the compressional gyroviscous force (Sec. III C, Eqs. 20 and 21, Fig. 7).
