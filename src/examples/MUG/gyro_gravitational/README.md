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

## Open items

- **Hall term.** With `xmhd_hall=T` a faster mode (γ ≈ 0.045 Ω_i) overtakes the gravitational mode.
  It also appears without gravity and does not depend on the time step or element order. Until it is
  understood, the example runs without the Hall term.
- **Gyroviscous mode in a non-uniform equilibrium.** Without the Hall term, a mode seeded by the
  discrete equilibrium's force imbalance grows at γ ≈ 0.0066 Ω_i (4×1×4 mesh) for every β. At 8×1×8 it
  grows at γ ≈ 0.019 Ω_i, above the largest possible buoyancy growth rate √(g/Ln) ≈ 0.0071 Ω_i. It
  persists with g = 0, needs gyroviscosity (with `nu_gyro<0` the growth stays at the ideal rate), and disappears for a
  homogeneous plasma. It hides the gravitational mode for β ≳ 2.5e-3, and at lower β the runs fail late
  (the time step collapses) once the kinetic energy reaches about 1e-9. The Hall-term mode above may
  share its cause.
- **Linearized model.** `xmhd_lin_run` cannot be used for this case: its operator drops the density
  dependence of J×B/n, which provides the buoyancy drive when the equilibrium is supported by a field
  gradient.
- **High-β case** (paper, Fig. 4 right: Lx = 2π/0.05, Ly = π/0.01, Ln = 1e5, g = 5e-4, I0 = 1), where
  compressibility and electromagnetic effects matter.
- **Simpler cases from the same paper:** homogeneous normal modes with Hall and gyroviscosity
  (Sec. III A, Figs. 1 and 2) and the compressional gyroviscous force (Sec. III C, Eqs. 20 and 21, Fig. 7).
