# Plan: 2D MUG Rayleigh-Taylor example (with and without magnetic field)

Status: Session A done (see "Session A report" at the end; corrections from Session A are marked
**[Session A]** below). Session B (the example) is next.

## 1. What the 2D solver does today (`src/physics/xmhd_2d.F90`)

- Fields: `n, vx, vy, vz, T, psi, by` on one Lagrange space. The mesh is a 2D surface mesh in
  (mesh x, mesh y). The solver maps them to **model (x, z)**, so **mesh y is the vertical z**
  (`basis_grads(3,:)=basis_grads(2,:)`), `coords(2)` is z, and **y is the ignorable direction**.
- Field: `B = grad(psi) x y_hat + by y_hat + B_0`, so `Bx = -dpsi/dz`, `Bz = dpsi/dx`, `By = by`.
- Momentum is solved per unit mass: `dv/dt = -v.grad v - grad(p)/(m n) + JxB/(m n) + visc`, with `p = 2 n e T`.
- Linear mode (`run_lin_simulation`) advances `J(u0) du^{n+1} = M du^n`, where `J` is
  `build_approx_jacobian` evaluated at the equilibrium `u0`.

What blocks an RT example:

1. **No gravity term.** Nothing like `g_accel` from the 3D `xmhd` module exists in 2D.
2. **The linear operator drops buoyancy.** In Cartesian linear mode the `(vel,n)` Jacobian block
   (`IF(linear)` branch, ~line 1399) differentiates `grad(p)` but freezes the `1/n` factor. It
   also omits the `n`-derivative of `JxB/(m n)`. The missing term, `-dn*F0/(m n0^2)` with
   `F0 = -grad p0 + J0xB0 = m n0 g z_hat`, is exactly the buoyancy drive (`-g dn/n0 z_hat`), so linear RT
   would come out stable. (This is the same limitation the 3D `gyro_gravitational` README reports.)
   The `(T,n)` linear branch (~line 1543) also differs from the exact derivative when `chi*grad(T0) != 0`.
   The non-linear branch of both blocks already contains the exact derivatives of the pressure and
   magnetic terms.
3. **Linear kinetic-energy diagnostic is wrong.** `mfun_apply` (line 673) computes
   `0.5 m n v^2` using `n` from the *perturbation* vector, so in linear runs `ven ~ dn*dv^2`.
   The growth rate is measured from `ven`, so this has to use the equilibrium density from `u0`.

**[Session A]** Two further blockers were found and are fixed or worked around:

4. **Spurious instability of stratified equilibria (fixed in the solver).** The momentum equation evaluated
   `grad(p)/(m n)` with `p = 2 k n T` pointwise from the separate equal-order n and T expansions. In any
   stratified state (even with `g = 0` and uniform pressure) this has a non-oscillatory spurious mode that grows at
   about `c_s/L` (~4e5 /s for the defaults, 30x the RT rate), does not converge away with resolution, and is not
   damped by implicit time stepping. The force now uses the gradient of the nodal pressure interpolant
   `p_h = sum_j 2 k n_j T_j phi_j` (a "group" FE form), which removes the mode: a 1D-in-z/Fourier-in-x model of the
   discretization then gives the RT rate to 0.07% (the compressibility correction) at every resolution tried.
5. **Dirichlet n and T on the walls (avoid in the example).** With `vz = 0` on a wall, `dn/dt = -n dvz/dz` is not
   zero there, so pinning n and T over-constrains the transport equations. In a stratified state this drives another
   spurious mode at ~`c_s/L`, with either pressure form. Use natural (zero-flux) conditions for n and T on the walls.

Minor, noted but out of scope: `mfun_apply` adds `psi_vac` to the perturbation's psi in linear mode.
The viscous `n`-derivative terms in the non-linear `(vel,n)` branch have a spurious term and a factor of 2.
They only matter when `v0 != 0`, but fix them while touching the block.

## 2. Physics problem

Domain: `x in [-Lx/2, Lx/2]` periodic, `z in [-Lz/2, Lz/2]` between rigid, perfectly conducting walls.
Gravity is `-g z_hat`. The unstable stratification is **exponential** (heavy plasma on top), chosen
because it has an exact incompressible growth rate, both with and without field, and has no sharp
interface to resolve:

- Density: `n0(z) = n_0 exp(z/L)`, with `L > 0` unstable and `L < 0` stable.
- In-plane field along the periodic direction, with constant Alfven speed:
  `Bx(z) = B_0 exp(z/(2L))`, `v_A = B_0/sqrt(mu0 m_i n_0)`, and stored as
  `psi0(z) = -2 L B_0 exp(z/(2L))` (field 6, not `B_0`, so it is advected and compressed consistently).
- Optional uniform out-of-plane field `By = by0` (field 7), with `k.B = 0`, to show that only the field
  along k stabilizes.
- Pressure from total force balance `d/dz(p + Bx^2/2mu0) = -m_i n g`:
  `p(z) = p_0 - (m_i n_0 g L + B_0^2/(2 mu0)) (exp(z/L) - 1)`, `p_0 = 2 e n_0 T_0`, `T(z) = p/(2 e n0(z))`.
  The driver must abort if `p <= 0` anywhere.

**Theory** (incompressible, ideal; Chandrasekhar 1961 §92/§97). For perturbations `~ exp(gamma t) cos(kx)`
with `k = 2 pi/Lx`, vertical mode `m`, the eigenfunction is `w(z) = exp(-z/2L) cos(m pi z/Lz)` (m odd) and

    gamma^2 = k^2 [ (g/L) / (k^2 + m^2 pi^2/Lz^2 + 1/(4L^2))  -  v_A^2 ]

A field along x stabilizes every mode once `v_A^2 > (g/L)/(k^2 + pi^2/Lz^2 + 1/(4L^2))`. A field along y
leaves `gamma` unchanged (apart from compressibility). Departures from this formula come from
compressibility (`g L/c_s^2`, `v_A^2/c_s^2`), viscosity and resistivity (`nu k^2`, `eta k^2`),
time discretisation (`gamma dt`), and resolution. The defaults below keep each of these to about 1% or less.

**Proposed defaults** (SI, T in eV): `Lx = Lz = 1 m`, `L = 0.5 m` (density ratio e^2 ≈ 7.4),
`n_0 = 1e19`, `T_0 = 100 eV` (T from 272 eV at the bottom to 37 eV at the top), `g = 1e8 m/s^2`, `gamma_ad = 5/3`.
Then `gamma_hydro = 1.25e4 /s`, `v_A,crit = 1993 m/s` (`B_0,crit = 2.89e-4 T`),
`gamma(B_0 = 0.5 B_c) = 1.08e4 /s`, `gamma(0.8 B_c) = 7.5e3 /s`, and `g/(c_s^2/L) < 0.5%`.
Suggested dissipation: `nu = eta = 1 m^2/s` (`nu k^2/gamma ≈ 0.3%`). `D_diff` and `chi` should be small
(≈ 1 m^2/s): they make the background drift slowly (rate `D/L^2 ≈ 4 /s`), which the linear model does not
see. Session B checks whether the non-linear late phase needs more. **[Session A]** The drift matters more than this
suggests: hydrostatic balance is a small difference (`rho g L/p ≈ 0.26%` here), so a fractional pressure drift `eps`
gives an acceleration of about `eps p/(rho L) ≈ 380 eps g`. With `D_diff = chi = 1` the non-linear run with no seed
reaches a kinetic energy of ~1e-6 J/m (v ~ 10 m/s) by `t = 1e-4 s`; with `D_diff = chi = 0` it stays at ~5e-9 J/m on
the coarse test mesh. Use `D_diff = chi = 0` (or much less than 1) for growth-rate measurements.

Seed: the m = 1 incompressible eigenmode in velocity only (no sound waves, the same seed works for both models):

    vz = V1 exp(-z/2L) cos(pi z/Lz) cos(kx)
    vx = (V1/k) exp(-z/2L) [cos(pi z/Lz)/(2L) + (pi/Lz) sin(pi z/Lz)] sin(kx)

Because `vx != 0` at the walls, use **free-slip**: Dirichlet `vz` (and `vy`) on the walls, none on `vx`. The
natural boundary condition is consistent, because `Bx Bz = 0` on a conducting wall. `psi` and `by` keep the
default Dirichlet wall conditions (`psi` fixed on the wall gives `dBz = 0` there). **[Session A]** `n` and `T` must
**not** be Dirichlet on the walls (blocker 5 in section 1): set `n_bc` and `T_bc` to all `.FALSE.` (zero-flux walls).

Growth rate: `gamma = (1/2) d ln(ven)/dt` from `xmhd_2d.hist`. With backward Euler the measured value
is biased by about `+gamma dt/2`, so use `timestep_cn=T` (BDF2 in linear mode, Crank-Nicolson in non-linear mode)
and/or `gamma dt ≲ 0.02`. **[Session A]** A velocity-only seed starts with zero displacement, so it is half growing and
half decaying mode: `v ∝ cosh(gamma t)` and the instantaneous rate is `gamma tanh(gamma t)` (0.96 gamma at
`gamma t = 2`). Fit only after `gamma t ≳ 3`, or fit `ven ∝ cosh^2(gamma t)`. The linear solver needs the
block-Jacobi/LU preconditioner (`oft_in_taylor_green.xml` style): the default diagonal one stalls. `lin_tol` is
absolute, so keep the perturbation amplitude moderate in long linear runs.

## 3. Session A: solver changes and regression test

All changes are in `src/physics/xmhd_2d.F90` plus a new test. Build and run with the `oft-build` and `oft-test` skills.

A1. **Gravity.** Add `REAL(r8) :: g_accel = 0.d0` to `oft_xmhd_2d_sim`, the acceleration in m/s^2 along `-z`
(mesh y), matching the 3D convention (`xmhd::g_accel`). In `nlfun_apply` add
`res_loc(jr,4) += self%dt*basis_vals(jr)*g_accel*int_factor` (times `coords(1)` in the cylindrical branch).
Gravity is per unit mass, so there is no Jacobian or `mfun` term. Document it in the type and the module header.

A2. **Exact linearization.** In `build_approx_jacobian`, remove the Cartesian `IF(linear)` special
cases in the `(vel,n)` and `(T,n)` blocks so that linear mode uses the exact derivatives (the current
non-linear branch). Correct the viscous `n`-derivative terms to the exact derivative of
`+nu grad(phi).grad(v) - nu phi grad(n).grad(v)/n`. Verify with a temporary finite-difference check
(`J du` against `(F(u+eps du)-F(u))/eps` at a non-uniform state with all fields non-zero). Do not commit the check.

A3. **Linear diagnostics.** In `mfun_apply`, compute `ven` with the equilibrium density from `parent_sim%u0`
(field 1). Leave the other diagnostics as they are.

A4. **Regression test** `src/tests/physics/MUG/test_rt_2d.{F90,py}`, registered in that directory's
`CMakeLists.txt` and modelled on `test_alfven_2d`. It must run in about a minute or less on a coarse mesh:
- (a) A stable stratification (`L < 0`) with `eps = 0` stays at rest (kinetic energy at the noise floor).
  This checks the sign and magnitude of gravity.
- (b) Linear growth rate for `B_0 = 0` and `B_0 = 0.5 B_c` within about 3-5% of the formula.
- (c) A short non-linear run whose growth rate matches the linear one to about 2%.

A5. Re-run all existing 2D MUG tests (`test_sound_2d`, `test_alfven_2d`, `test_taylor_green_2d`,
including the linear and MPI variants). The uniform equilibria they use should not change with A2.

Deliverable: commits on the dev branch, with a short report of the measured growth rates against theory and the test timing.

## 4. Session B: the example (starts from Session A's branch)

New directory `src/examples/MUG/rayleigh_taylor_2d/`, following the 3D examples and `harris_sheet`:

- `MUG_rt_2d.F90`: doxygen page `{#doc_mug2d_rt_ex}` with Problem setup / Theory / Code walk-through /
  Input files / Results sections, like `MUG_gyro_grav.F90`.
  - A helper module with the profile functions (`n`, `T`, `psi`, and the `vx`/`vz` seed, using `pt(2)` as z).
  - Namelist `rt_2d_options`: `order, linear, plot_run, pm, n_0, T_0, L_n, g, B_0, by0, V1 (seed),
    nu, eta, chi, D_diff, gamma, dt, nsteps, rst_freq, use_mfnk, timestep_cn, ittarget, lin_tol, nl_tol, den_scale`.
    The domain size comes from the mesh (as in `MUG_gyro_grav`), and the driver checks the mesh is centred in z.
  - Print the theoretical `gamma`, `B_c`, and the compressibility parameter at start-up.
  - Linear: `u0` is the equilibrium and `u` the seed. Non-linear: `u = u0 + seed`. Then `run_lin_simulation`
    or `run_simulation`. `plot_run=T` calls `xmhd_2d_plot`.
  - Use free-slip wall flags for `vx` and natural (all `.FALSE.`) flags for `n` and `T` (section 2), set after
    `mhd_sim%setup`. **[Session A]** `src/tests/physics/MUG/test_rt_2d.F90` already has the profile functions, seed,
    boundary flags and the `oft_in_taylor_green.xml` preconditioner; start from it.
- `oft.in`: the default non-linear hydro case on the cube mesh (`mesh_type=1`, `ref_per=T,F,F`,
  starting resolution about 32x32 at order 2-3).
- `oft_in.xml`: GMRES with a block-Jacobi/LU preconditioner under `<xmhd2d>`, as in `oft_in_taylor_green.xml`.
- `run_scan.py`: runs linear (and optionally non-linear) cases over `B_0/B_c in {0, 0.25, 0.5, 0.75, 0.9, 1.1}`
  plus a `by0` case, fits `gamma` from `xmhd_2d.hist`, and compares with the formula (table and `--plot`),
  as in `gyro_gravitational/run_scan.py`.
- `CMakeLists.txt`. Add the directory to `src/examples/MUG/CMakeLists.txt`.
- Docs: add a "2D examples" list to `src/docs/doc_mug_main.md` containing the existing (unlisted) Harris sheet page
  and the new page. Images of the non-linear density evolution (B = 0 against B ≈ 0.5 B_c, showing the stabilized
  short scales and the rising bubbles and falling spikes) can go in `src/docs/images`.

Validation to report: a linear scan table (measured against theory); linear vs non-linear growth agreement; the
non-linear noise floor with `V1 = 0` (the discrete imbalance of the projected equilibrium has to stay well below
the seed); and a non-linear run carried into saturation without negative density or temperature, with whatever
dissipation that needs.

## 5. Open choices (defaults chosen as above unless you say otherwise)

1. **Profile.** The exponential profile gives an exact benchmark with field. A classic tanh interface looks more
   like the textbook mushrooms but has no closed-form growth rate with field (only a numerical eigen-solve).
   It could be added later as `profile='tanh'` for non-linear visuals.
2. **Field options.** Constant-`v_A` in-plane field (analytic) plus an optional uniform `by0`. A uniform in-plane
   field via `B_0` was rejected: the static `B_0` is not compressed, and the theory is not closed-form.
3. **Branching.** Session A commits to `claude/practical-fermi-j3e4is`. Session B starts from that branch after A
   finishes and also commits there. One PR at the end.

## Session A report

Branch `claude/practical-fermi-j3e4is`. All solver changes are in `src/physics/xmhd_2d.F90`.

### What changed

| Commit | Change |
|---|---|
| MUG: Add gravity to the 2D xMHD solver | A1. `g_accel` [m/s^2] along -z (second mesh coordinate), applied per unit mass in the momentum residual of both branches. It has no Jacobian or `mfun` term. |
| MUG: Use the exact Jacobian for linear 2D xMHD runs | A2. Removed the Cartesian `IF(linear)` cases in `(vel,n)` and `(T,n)`. Fixed the `(vel,n)` viscous terms (one spurious term, one factor 2). Also fixed three cylindrical inconsistencies found by the FD check: `(vx,n)`/`(vy,n)` hoop terms (spurious `1/R`, missing 1/2, spurious viscous term) and the sign of the `v_R by/R^2` term in `(by,vx)`. |
| MUG: Use the equilibrium density for kinetic energy in linear 2D runs | A3. `ven` uses `n0` from `u0` in linear runs. |
| MUG: Use an interpolated pressure in the 2D xMHD momentum equation | New (blocker 4 in section 1). The pressure force is `grad(p_h)/(m n)`, with `p_h = sum_j 2k n_j T_j phi_j`, in the residual and Jacobian. The incompressible branch is unchanged. |
| MUG: Add 2D Rayleigh-Taylor regression test | A4. `src/tests/physics/MUG/test_rt_2d.{F90,py}`. |

### FD Jacobian check (temporary program, not committed)

Setup: 4x4 non-periodic triangle mesh, order 3, Cartesian and cylindrical (R in [1,2]). The state is non-uniform
with all 7 fields non-zero (|v| ~ 3e4 m/s, T ~ 100 eV with 30% variation, n with 30% variation, psi and by
non-zero), `B_0 = (0.01, 0.02, 0.03)` (Cartesian), `nu, eta, chi, D_diff` = 5, 2, 3, 4, `g = 1e8`, and dt = 1e-6.
The perturbation was applied one field at a time. The check compared `J du` with the central difference
`(F(u+eps du) - F(u-eps du))/(2 eps)`, eps = 1e-6, and reports `|J du - FD| / max|FD|` per 7x7 block, with no
Dirichlet rows.
- After A2 and the pressure change, every block is at most 1.2e-9 in both geometries.
- Before the cylindrical fixes, `(vx,n)` = 1.6e-2, `(vy,n)` = 6.5e-3 and `(by,vx)` = 3.0e-2. All other blocks were
  already exact once the Cartesian `IF(linear)` cases were removed.

### Growth rates (4x4 cells, order 3, dt = 2e-6, BDF2 linear / CN non-linear, nu = eta = 1, D_diff = chi = 0)

Theory: `gamma = k sqrt((g/L)/(k^2 + pi^2 + 1/(4L^2)) (1 - (B_0/B_c)^2))`. The fit is `(1/2) d ln(ven)/dt` over steps
150-200 (`gamma t` = 3.75-5, where `tanh(gamma t)` = 0.999).

| Case | Measured [1/s] | Theory [1/s] | Ratio |
|---|---|---|---|
| Linear, B_0 = 0 | 1.25051e4 | 1.25229e4 | 0.9986 |
| Linear, B_0 = 0.5 B_c | 1.08288e4 | 1.08451e4 | 0.9985 |
| Linear, B_0 = 0, 2 MPI ranks | 1.25051e4 | 1.25229e4 | 0.9986 |
| Non-linear vs linear, V1 = 10 m/s, steps 60-100 | 1.20134e4 | 1.20063e4 (linear) | 1.0006 |

- Resolution barely matters for the linear rate: 4x4 p3, 8x8 p2 and 8x8 p3 agree to 1e-4.
- On 8x8 p3 the instantaneous rate divided by `gamma tanh(gamma t)` settles at 0.9973.
- The remaining ~0.2% comes from viscosity, resistivity and compressibility.
- The stable linear case (L = -0.5) oscillates at the gravity-wave period (about 1.25e4 rad/s), with no growth.

### Test runtimes (4 cores, OMP_NUM_THREADS=2)

| Test | Time |
|---|---|
| `test_stable_at_rest` | 5 s |
| `test_lin_growth[0.0]` | 1.2 s |
| `test_lin_growth[0.5]` | 1.2 s |
| `test_lin_growth_mpi` | 6.5 s |
| `test_nl_growth` | 13 s |
| **Total** | **27 s** |

### A5: existing 2D MUG tests

`test_sound_2d`, `test_alfven_2d`, `test_taylor_green_2d` and `test_rt_2d`, including the linear and MPI variants:
63 passed and 18 skipped, in 9.4 min. The skips are 12 `Not compiled with PETSc` and 6 matrix-free+PETSc
combinations the tests exclude.
- Linear sound-wave errors are unchanged to all printed digits.
- Non-linear sound-wave errors changed by under 0.1% (p3: 2.1413e-3 against the stored 2.1429e-3).

### For Session B

- **Boundary conditions.**
  - `n_bc = T_bc = .FALSE.` everywhere (natural walls). Dirichlet n/T walls are unstable in any stratified state.
  - Free-slip: `velx_bc = .FALSE.`.
  - `vz`, `vy`, `psi` and `by` keep the default walls.
  - Zero the seed's `vz` on boundary nodes after projection.
  - `setup_bc` still defaults n and T to Dirichlet on `gbe`, so the example must override it.
- **Preconditioner and solver settings.**
  - Use `oft_in_taylor_green.xml` (GMRES(1) with block-Jacobi/LU). The default diagonal preconditioner does not
    converge, linear or non-linear.
  - Set `ittarget` large (the test uses 1000) so the adaptive time step stays at `dt`.
  - `lin_tol`/`nl_tol` are absolute: 1e-10 and 1e-8 work at these amplitudes.
- **Dissipation.**
  - Use `D_diff = chi = 0` (or much less than 1) for clean growth. With 1 m^2/s, the no-seed non-linear floor reaches
    ~8e-7 J/m after 40 steps (8e-5 s) on 4x4 p3.
  - `nu = eta = 1` costs about 0.2% of gamma.
  - The late non-linear phase without D_diff/chi has not been tried.
- **Noise floor** (non-linear, unstable stratification, V1 = 0, D_diff = chi = 0, after 40 steps):

  | Mesh | ven [J/m] |
  |---|---|
  | 4x4 p3 | 5e-9 |
  | 8x8 p3 | 2e-11 |
  | 4x4 p4 | 5e-12 |

  A V1 = 1 m/s seed starts at ven ≈ 2.7e-9 J/m. On 4x4 p3 use V1 ≥ 10 m/s (the test does), or a finer mesh.
- **Seed transient.** The velocity-only seed grows as `cosh(gamma t)`. In `run_scan.py`, fit after `gamma t ≳ 3`
  or fit a `cosh^2` model.
- **Cost.**
  - Linear: ~0.01 s per step on 4x4 p3 and ~0.04 s per step on 8x8 p3.
  - Non-linear (the Jacobian and LU are rebuilt every step): ~0.13 s per step on 4x4 p3 and ~1 s per step on 8x8 p3.
- **`B_c` for a stable profile.** The test computes `B_c` with `ABS(g/L)` so that L < 0 cases do not take the
  square root of a negative number.

### Unresolved

- **Stable-side mode in the model.** The 1D-in-z/Fourier-in-x model of the new discretization shows a weak spurious
  real mode for *stable* stratification. With no dissipation its rate is ~0.44 x the RT scale at k = 2 pi, and it
  scales as sqrt(g). `nu = 1` roughly halves it. It was not visible in the 2D linear stable run (pure gravity-wave
  oscillation), nor in the non-linear stable test (KE ~1.5e-7 of free fall). Watch for it if an example uses a
  stable layer.
- **Harris sheet example.** `src/examples/MUG/harris_sheet` uses stratified density with default (Dirichlet) n/T
  walls, so it may be exposed to blocker 5. Not checked. Its results may also shift slightly with the
  interpolated pressure.
- **Incompressible Jacobian.** The `(vel,n)` Jacobian still uses the compressible formula with the pressure field.
  This is harmless because n is fixed when incompressible, so it was left unchanged.
- **`psi_vac` in linear runs.** `mfun_apply` still adds `psi_vac` to the perturbation's psi in linear runs (noted in
  section 1, out of scope).

## Merge with claude/peaceful-dijkstra-9cjvbz

That branch (742c922, dd22599) added the 2D equations page `src/docs/doc_mug_xmhd_2d_eqns.md`, a finite-difference
Jacobian test (`test_jacobian_2d`), input checks (`check_settings`) and several `xmhd_2d.F90` fixes: B_0 and
cylindrical F terms, the incompressible constraint, and `psi_vac` in linear runs. It was merged as b9a1210, using
their `xmhd_2d.F90` as the base with Session A's changes re-applied:
- gravity
- exact linearization: their approximate `IF(linear)` branches in `(vel,n)`/`(T,n)` are removed, because they assumed
  `F0 = m n v.grad v`, which is wrong with gravity
- linear kinetic energy using u0's density
- the group FE pressure

Terms both branches had fixed keep their version.

- `test_jacobian_2d` now also sets `g_accel`.
- The 2D tests give 72 passed and 18 PETSc skips.
- RT rates are unchanged.
- The equations page documents gravity, the pressure discretization (`doc_mug_xmhd_2d_eqns_pres`, with references),
  the natural n/T wall note and the exact linear operator.

For Session B:
- `check_settings` now aborts unless `nu`, `eta(:,1)`, `D_diff` and `chi` are >= 0, `eta(:,2)` and `m_i` are > 0,
  and `gamma` is > 1. The RT example already satisfies this with `D_diff = chi = 0`.
- `doc_mug_main.md` now has a "2D MHD" section listing the equations page and the Harris sheet example. Put the new
  RT example there rather than in a separate list.
