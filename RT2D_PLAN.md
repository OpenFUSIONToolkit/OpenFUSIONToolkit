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

## Session B report

Branch `claude/practical-fermi-j3e4is`. Built against Ubuntu 24.04 system packages (GCC 13, OpenMPI, OpenBLAS, SuperLU,
UMFPACK; `OFT_BUILD_EXAMPLES=ON`) on 4 cores.

### What was added

| Commit | Change |
|---|---|
| MUG: Honor pre_freq in matrix-free 2D xMHD runs | One-line fix in `xmhd_2d.F90` `run_simulation`. The check that forces a preconditioner update after a time-step change, `ABS(1-(dt-dt_old)/dt_old) > 0.1`, is always true for a constant dt, so the LU preconditioner was refactored every step and `pre_freq` did nothing with `mfnk=T`. It now tests `ABS(dt-dt_old) > 0.1*dt_old`. Runs with `mfnk=F` or `pre_freq=1` are unchanged. |
| MUG: Add 2D Rayleigh-Taylor example | `src/examples/MUG/rayleigh_taylor_2d/`: `MUG_rt_2d.F90` (doc page `doc_mug2d_rt_ex`), `oft.in`, `oft_in.xml`, `run_scan.py`, `plot_rt_2d.py`, `CMakeLists.txt`; registered in `src/examples/MUG/CMakeLists.txt`. `doc_mug_main.md`: a "2D examples" subsection (Harris sheet + RT) under the existing "2D MHD" section. |
| MUG: Tune solver settings in the 2D Rayleigh-Taylor input file | `pre_freq=4`, `lin_tol=1e-6`, `nl_tol=1e-5` in `oft.in`. |
| MUG: Document results of the 2D Rayleigh-Taylor example | First results (exponential profile only), superseded by the tanh redesign below. |
| MUG: Add a tanh transition-layer profile to the 2D RT example (and follow-ups) | `profile='exp'` (linear benchmark, `run_scan.py`) or `'tanh'` (required for nonlinear runs). New `run_nonlinear.py` (linear vs nonlinear growth of the cos(kx) component of vz, `--skip_run` to re-analyse). Default `oft.in`: nonlinear tanh, 48x96 p2, adaptive dt, D = chi = 1. Doc results rewritten; images `MUG_rt2d_ex-scan.png`, `MUG_rt2d_ex-modes.png`, `MUG_rt2d_ex-density.png`. |

The driver:
- Follows `test_rt_2d.F90`: profiles, seed, free-slip `vx`, natural n/T walls, vz zeroed on the walls.
- Takes the domain size from the mesh and aborts unless the mesh is centred in z, or if p <= 0 at a wall.
- Prints gamma, B_c, g L/c_s^2, v_A^2/c_s^2 and the wall temperatures.
- Supports linear and nonlinear runs and `plot_run` (`xmhd_2d_plot`). `B_0` is the in-plane field in T, and `by0` a uniform out-of-plane field.
- Namelist `rt_2d_options` as planned, plus `pre_freq`.

`run_scan.py`:
- Runs `B_0/B_c in {0, 0.25, 0.5, 0.75, 0.9, 1.1}` plus `by0 = 1e-2 T` (35 B_c, beta ~ 8), each to gamma t = 6.
- Fits `E_k = A cosh^2(gamma t)` to the whole history, and checks it with the slope over gamma t > 3.
- Forces `profile='exp'` in a 1 m x 1 m box. Options: `--plot`, `--np`, `--ni`, `--order`, `--dt`, `--nu`, `--eta`, `--V1` (`--nonlinear` was removed with the tanh redesign).

`plot_rt_2d.py`:
- Reads the XDMF/HDF5 plot files with `OpenFUSIONToolkit.io.XDMF_plot_file` and merges the MPI pieces.
- Draws density snapshots (tripcolor), with optional per-run times.
- Prints min(n) and min(T) over each run.

### Linear growth-rate scan (8x8, p3, dt = 2e-6, nu = 1, eta = 0.1, D_diff = chi = 0, V1 = 1; ~45 s total on 1 core)

| B_0/B_c | by0 [T] | Theory [1/s] | cosh^2 fit | Ratio | gamma t > 3 | Ratio |
|---|---|---|---|---|---|---|
| 0 | 0 | 1.2523e4 | 1.2508e4 | 0.9988 | 1.2502e4 | 0.9984 |
| 0.25 | 0 | 1.2125e4 | 1.2111e4 | 0.9988 | 1.2105e4 | 0.9984 |
| 0.5 | 0 | 1.0845e4 | 1.0832e4 | 0.9988 | 1.0827e4 | 0.9984 |
| 0.75 | 0 | 8.2831e3 | 8.2751e3 | 0.9990 | 8.2722e3 | 0.9987 |
| 0.9 | 0 | 5.4586e3 | 5.4649e3 | 1.0012 | 5.4645e3 | 1.0011 |
| 1.1 | 0 | stable | stable (no growth, oscillates) | | | |
| 0 | 1e-2 | 1.2523e4 | 1.2507e4 | 0.9987 | 1.2501e4 | 0.9983 |

- **eta = 1 at 0.9 B_c.** With eta = 1 (the test's value), the 0.9 B_c case is 1.0157.
  - Not resolution: unchanged on 16x16.
  - Not the time step: unchanged at dt = 5e-7.
  - It is resistive slippage, amplified near marginal stability by 1/(1 - b^2): eta = 0.1 gives 1.0012, and nu = 0.1 with eta = 1 gives 1.020. So the scan defaults to eta = 0.1. The doc page explains this.
- **Linear vs nonlinear** (`run_scan.py --nonlinear`, same settings, V1 = 1): the slope over gamma t > 3 is
  - 0.9972 (B_0 = 0) and 0.9982 (0.5 B_c),
  - against 0.9984 linear.
- **Noise floor** (nonlinear, V1 = 0): the seeds start well above it.
  - 8x8 p3: 1.6e-11 J/m after 10 steps, growing at about the RT rate (1e-9 by step 100). The V1 = 1 seed starts at 2.7e-9.
  - Default 32x32 p2: 6e-15 J/m after 50 steps. The V1 = 100 seed starts at 2.7e-5.

### Redesign: exponential profile for the linear benchmark, tanh layer for nonlinear runs

At the user's direction the exponential profile is kept only for the linear benchmark (scan above), and nonlinear runs
require `profile='tanh'`: `n = n_0 (1 + A tanh(z/delta))`, uniform `B_x = B_0`, `p` from force balance, seed
`cos(pi z/L_z) exp(-(z/sigma)^2)`. Defaults: 1 m x 2 m box, A = 0.5, delta = 0.1 m, sigma = 0.2 m, T from 200 eV
(bottom) to 66 eV (top). The driver prints the sharp-interface estimate `gamma^2 = A g k tanh(k L_z/2) - k^2 v_A^2`
(1.772e4 /s, B_c = 4.09e-4 T).

Why the exponential profile was dropped for nonlinear runs (all documented on the doc page):
- The bubble/spikes reach the walls while growing, and the runs fail after wall contact (32x32: ~470 us).
- D_diff or chi > 0 fails on the first Newton step. Investigation (temporary FD program, not committed): the Jacobian
  matches FD and passes a Taylor test; the solve converges with a uniform density, a flat n with D, a flat T with chi,
  Dirichlet n/T walls, smaller dt, or with D acting on n - n0 only. So it is the zero-flux walls being inconsistent with
  the equilibrium gradient (a sub-grid boundary layer must form), not a Jacobian bug.
- Dirichlet n/T walls (8x8 p3, nonlinear, V1 = 1): spurious mode at 3.4e4 /s (vs RT 1.25e4), and with D = 1 the run
  fails at step 124.

Choice of tanh parameters (linear and nonlinear test runs):
- delta = 0.04-0.05 m: grid-scale spurious growth on coarse meshes (1.9e5 /s at 16x32, 5e4 /s at 32x64), and a
  kinetic-energy floor (2e-8 J/m) above a 1 m/s seed. delta = 0.1 m on 48x96 p2: floor 4e-10 J/m (an x-independent
  adjustment of the projected layer, which saturates), against a seed energy of 1.7e-9 J/m. Growth is therefore measured
  from the cos(kx) component of vz, which the adjustment does not contribute to.
- nu = 10 m^2/s damps the short-wavelength layer modes; eta = 1.
- Linear growth rate, B_0 = 0: 1.3922e4 (32x64) vs 1.3909e4 (48x96), 0.09%. 64x128 could not be run: UMFPACK fails
  to factorize (not a machine-memory limit; 11 GB free).

### Linear vs nonlinear (tanh, 48x96 p2, dt = 2e-6, V1 = 1, nu = 10, eta = 1, D = chi = 0)

`run_nonlinear.py`, cos(kx) amplitude of vz, both runs fit over the window where the nonlinear amplitude is 3-30 m/s:

| B_0 [T] | Window [us] | gamma linear | gamma nonlinear | Ratio | Sharp interface |
|---|---|---|---|---|---|
| 0 | 240-400 | 1.4010e4 | 1.4004e4 | 0.9996 | 1.7724e4 |
| 2e-4 (0.49 B_c,sharp) | 310-520 | 1.0633e4 | 1.0632e4 | 0.9999 | 1.5461e4 |

The amplitudes agree to < 0.05% up to a few hundred m/s. `run_nonlinear.py` forces D = chi = 0, since diffusion also
spreads the layer in the nonlinear run, which the linearized model does not include (D = chi = 1 lowers the nonlinear
rate by 0.65%).

### Nonlinear evolution (tanh)

- **Fixed dt fails** at ~710 us (B_0 = 0) and ~990 us (2e-4 T): GMRES iterations climb once the flow is nonlinear
  (velocity ~ 1 cell per step). Rebuilding the preconditioner every step (`pre_freq=1`) fails at the same time
  (712 us). With `ittarget=60` (adaptive dt) the runs continue: dt falls to 2-7e-7 s around roll-up/wall contact.
  The default `oft.in` now uses `ittarget=60`.
- **D = chi = 0 with adaptive dt** (B_0 = 0): n and T stay within their initial range until ~730 us, when the
  Kelvin-Helmholtz roll-ups thin the layer to the grid scale; grid-scale oscillations then fill the domain by 840 us
  (n from 0.5x min to 2.3x max). B_0 = 2e-4 T reached 1.44 ms without a solver failure (not inspected for noise).
- **D = chi = 1 with adaptive dt** (the new default): B_0 = 0 runs to 1.0 ms (bubble at the top wall at ~880 us,
  secondary structure afterwards), n within ~5% of its initial range (one 15% transient at 923 us), T >= 63 eV.
  B_0 = 2e-4 T runs to 1.3 ms (stopped by hand) with n within 1% of its initial range; the field suppresses the KH roll-up.
- Cost on 2 OpenMP threads: ~4-5 s/step with fixed dt, 15-20 s/step after roll-up (runs shared 4 cores). The default
  run (450 steps, to 790 us) takes about an hour. `pre_freq=1` doubles the cost of the linear phase (8.5 vs 4.2 s/step).

### Tests

- `test_rt_2d`: 5/5 passed after the final rebuild (OMP_NUM_THREADS=1). The MPI variant needs
  `OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` when run as root in this sandbox (it fails otherwise).
- Earlier, after merging Session A's branch and with the pre_freq fix: the 2D MUG tests (`test_jacobian_2d`,
  `test_sound_2d`, `test_alfven_2d`, `test_taylor_green_2d`, `test_rt_2d`, including linear and MPI variants) gave
  72 passed and 18 skipped (12 `Not compiled with PETSc`, 6 matrix-free+PETSc). No tolerances or reference values were
  changed. `xmhd_2d.F90` has not changed since then.
- `ruff check` passes on the example scripts.

### Unresolved

- Without particle diffusion/thermal conduction, Galerkin advection of n and T is undamped, so late nonlinear phases
  need D, chi > 0 (or a stabilized advection scheme, not attempted).
- With the exponential profile, D or chi > 0 is incompatible with zero-flux walls (above). A diffusion-consistent
  equilibrium, or diffusion acting on n - n0, would avoid it; not pursued since the exponential profile is now linear only.
- UMFPACK fails to factorize the 64x128 p2 Jacobian; finer meshes need another LU package or MPI.
- Doxygen was not available in the sandbox; the doc page was checked through `generate_doc.py` only.

## SUPG upwinding session report

Branch `claude/mug2d-supg-upwind` (from `claude/practical-fermi-j3e4is`). Built against Ubuntu 24.04 system packages
(GCC 13, OpenMPI, OpenBLAS, SuperLU, UMFPACK; `OFT_BUILD_EXAMPLES=ON`) on 4 cores. Run directories (not committed) were
under `builds/rt_runs/`.

### What changed

| Commit | Change |
|---|---|
| MUG: Add streamline upwinding to the 2D xMHD solver | `oft_xmhd_2d_sim%upwind` (default `.FALSE.`). Inconsistent SUPG streamline diffusion `w (V_ad.grad phi)(V.grad f)` on the advection of n, each component of V, and T (compressible only), added to `nlfun_apply` with exact Jacobian terms in `build_approx_jacobian`. Helpers `upwind_active`, `upwind_weight`, `get_vel_weights`. Off in linear runs (`linear=T`, or inside `run_lin_simulation`, via the private flag `lin_run`). `test_jacobian_2d` extended; driver option `upwind` in the RT example; upwinding section `doc_mug_xmhd_2d_eqns_upwind` in the 2D equations page. |
| MUG: Cache the 2D upwinding element sizes | `he_cell(:)` computed once on first use (`upwind_elem_sizes`), instead of in every residual/Jacobian call. Cuts the per-step overhead from ~16% to ~5%. |
| MUG: Use upwinding in the 2D Rayleigh-Taylor example | `oft.in`: `upwind=T`, `D_diff = chi = 0.3`. New doc subsection `doc_mug2d_rt_ex_upwind`, updated evolution/cost/exp-profile text, new images, regenerated `MUG_rt2d_ex-density.png`. |

### Design choices

- **Same scheme as 3D.** Weight `w = beta h_e/(2|V_ad|)`, `beta = coth(Pe) - 1/Pe`, `Pe = |V_ad| h_e/(2 nu_f)`. The test function
  modification is applied to the advection term only (inconsistent SUPG). Diffusivities: `D` for n, `nu` for V (the 2D
  viscosity is already kinematic, so no `n0/n` factor as in 3D), `(gamma-1) chi` for T. In the T row the term is divided by
  `gamma-1` like the rest of the row.
- **Advecting velocity frozen at the start of the step** (`V_ad` from `self%u`, which holds `u^n` during a step in
  `run_simulation`). This mirrors the 3D `vad` from the previous step. It makes `w` and `V_ad` constants during a step, so the
  terms are linear in V and f and the Jacobian is exact. The CN right-hand side (dt = -dt/2, evaluated at `u^n`) uses the same
  `V_ad`, so the term is time-centred like everything else in 2D.
- **In-plane velocity only** for `V_ad` and `|V_ad|`: `V_y` advects nothing (d/dy = 0), and including it in `|V_ad|` would
  shrink the weight.
- **Element size** for 2D: node spacing of an equilateral triangle (`sqrt(4A/sqrt(3))/p`) or square (`sqrt(A)/p`) with the
  cell's planar area (from the quadrature, so curved cells are handled), also in cylindrical geometry.
- **Fields.** n, V and T only. psi and B_y are not upwinded (as in 3D): streamline diffusion of the flux would be an
  artificial resistivity along the flow, breaking flux freezing. In cylindrical geometry each component is upwinded with
  `V.grad V_k`; the centrifugal and Coriolis terms are not part of it. Incompressible: n and V only (field 5 is pressure).
- **Zero diffusivities allowed.** The 3D checks abort if `nu_visc`, `d_dens` or `kappa_perp` are not positive with
  upwinding. That would forbid the main goal here (D = chi = 0), and the weight is continuous as nu -> 0 (beta -> 1, full
  upwinding), so the 2D solver keeps its existing checks (D, chi, nu >= 0) and uses full upwinding where the diffusivity
  is zero. Documented as a deliberate difference.
- **Default `.FALSE.`** (3D defaults to `.TRUE.`). Turning it on changes every existing nonlinear 2D result, including the
  MUG-TokaMaker coupling (`mugtok_td`) and the Harris sheet example, which were validated without it, and the frozen
  `V_ad` relies on `self%u` holding the start-of-step state, which callers that drive `nlfun` themselves (eg. tests) must
  provide. For the small-amplitude wave tests upwinding is also not needed. The RT example turns it on.

### FD Jacobian check (`test_jacobian_2d`)

- New parameters `upwind` and `no_diff`. With `upwind=T`, `mhd_sim%u` (the advecting state) is set to a second
  pseudo-random velocity, different from the state the Jacobian is taken at, and dt = 1e-3 (instead of 1e-5) so the
  advective terms are not swamped by the mass matrix.
- 19 cases: planar/cylindrical x compressible/incompressible x order 2/3 x upwind off/on, plus upwind with D = chi = 0
  (planar and cylindrical), plus the B_0 abort test. All pass; largest block error 4.2e-7 (tolerance 1e-6).
- Sensitivity check (temporary, not committed): scaling the (n,V), (V,V) and (T,T) upwinding Jacobian terms by 1.01 gives
  block errors of 1e-5 to 3e-3, so each term is tested (with dt = 1e-5 the (V,V) and (T,T) errors stayed near 1e-6, hence
  the larger dt).

### Tests

- 2D MUG tests (`test_jacobian_2d`, `test_sound_2d`, `test_alfven_2d`, `test_taylor_green_2d`, `test_rt_2d`, including linear
  and MPI variants): 82 passed, 18 skipped (12 `Not compiled with PETSc`, 6 matrix-free+PETSc). 72 of these are the previous
  tests, unchanged (upwinding is off by default, so no reference values moved); 10 are the new Jacobian cases.
- Linear runs ignore `upwind`: an 8x8 p3 exponential-profile linear run gives bit-identical `ven` with `upwind=T` and `F`,
  so `run_scan.py` results are unchanged.

### RT nonlinear results (tanh, 48x96 p2, nu = 10, eta = 1, V1 = 1, adaptive dt, `pre_freq=4`)

Growth of the cos(kx) component of v_z in the nonlinear run (fit window 3-30 m/s, as `run_nonlinear.py`):

| Case | gamma [1/s] | vs no upwinding, D = chi = 0 |
|---|---|---|
| B_0 = 0, D = chi = 0 | 1.40044e4 | |
| B_0 = 0, D = chi = 0, upwind | 1.40046e4 | +0.001% |
| B_0 = 0, D = chi = 0.1, upwind | 1.39951e4 | -0.07% |
| B_0 = 0, D = chi = 0.3, upwind | 1.39771e4 | -0.20% |
| B_0 = 0, D = chi = 1, upwind | 1.39145e4 | -0.64% |
| B_0 = 0, D = chi = 1 | 1.39144e4 | -0.64% |
| B_0 = 2e-4 T, D = chi = 0 | 1.06317e4 | |
| B_0 = 2e-4 T, D = chi = 0, upwind | 1.06287e4 | -0.03% |

How far the B_0 = 0 runs stay clean (n relative to its initial range 5e18-1.5e19, over the run up to 1 ms):

| D = chi | Upwind | n leaves +-5% at | n range | T range [eV] | Run to | Wall time to 1 ms (1 core, 4 runs sharing) |
|---|---|---|---|---|---|---|
| 0 | no | 747 us (noise from ~730 us) | 0.91-1.29 | 52-220 | 791 us, stopped (noise everywhere, dt 0.35 us) | |
| 0 | yes | 878 us (noise from ~910 us) | 0.38-1.94 | 33-524 | 988 us, stopped (dt 0.23 us, 136 its) | |
| 0.1 | yes | 912 us | 0.61-1.75 | 40-338 | 974 us, stopped | |
| 0.3 | yes | 931 us | 0.95-1.08 | 65-222 | 1 ms | 5800 s, 547 steps |
| 1 | yes | never | 1.00-1.04 | 66-209 | 1 ms | 5100 s, 523 steps |
| 1 | no | 892 us | 0.85-1.08 | 63-241 | 1 ms | 11200 s, 788 steps |

- Upwinding alone (D = chi = 0) moves the onset of grid-scale noise from ~730 us to ~910 us: it keeps the roll-up and the
  spikes' arrival at the lower wall (~870 us) clean. It then fails where the heavy plasma spreading along the lower wall from
  both sides collides at x = 0 (a stagnation point with a thin upward jet) and in the lower corners, with ripples across the
  wall-parallel layers. This is the expected limitation of streamline-only diffusion (no crosswind/discontinuity-capturing
  term). The D = chi = 1 reference shows the same jet and a 15% transient there, but much smaller.
- Some explicit diffusion is still needed past the wall collision; 0.3 m^2/s with upwinding is cleaner than 1 m^2/s without,
  costs 0.20% of the growth rate (vs 0.64%), and runs ~2x faster to 1 ms. The example's `oft.in` now uses `upwind=T`,
  `D_diff = chi = 0.3`.
- Upwinding also keeps the adaptive dt at 2 us until ~790 us (vs ~690-700 us without), i.e. fewer GMRES iterations in the
  roll-up phase.
- B_0 = 2e-4 T, D = chi = 0, upwind: clean to 2.0 ms (n within 1.019-1.026 of its range at 2 ms, no undershoot; T 66-203.5 eV),
  through the bubble reaching the top wall and the pinch-off of its stem; dt fixed at 2 us throughout (stopped at 2 ms). The
  D = chi = 1 reference without upwinding (also run to 2 ms, dt fixed at 2 us) looks the same apart from a slightly broader
  layer (n within 1.011-1.018 at 2 ms). With the new defaults (upwind, D = chi = 0.3) n stays within 1% to 1.3 ms.
- Exponential profile past wall contact (32x32 p2, V1 = 100, D = chi = 0, nonlinear; driver check temporarily bypassed, not
  committed): without upwinding grid-scale noise fills the domain by 500 us (n 0.29-2.17 of its range at 516 us); with
  upwinding the bulk stays smooth, but the collision line at x = 0 still oscillates (n 0.64-1.73 at 530 us). Both stopped.
- Snapshots: `src/docs/images/MUG_rt2d_ex-upwind.png` (B_0 = 0: no upwind D = chi = 0 / upwind D = chi = 0 / upwind
  D = chi = 0.3 / D = chi = 1), `MUG_rt2d_ex-upwind-field.png` (B_0 = 2e-4 T), and the regenerated
  `MUG_rt2d_ex-density.png` (new defaults).

### Cost

- 48x96 p2, 2 OpenMP threads, first 40 steps (fixed dt): 3.78 s/step with upwinding vs 3.61 s/step without (+5%); GMRES
  iterations 5.25 vs 4.58 per step. Before caching the element sizes it was 4.10 s/step (+16%).
- Over a whole run upwinding is cheaper, because dt stays larger (table above).

### Unresolved

- Streamline upwinding does not control stagnation-point/collision oscillations or crosswind ripples (B_0 = 0 after the
  wall collision, ~920 us; exponential profile after wall contact). A crosswind or discontinuity-capturing term (eg. a
  residual-based shock-capturing diffusivity) would be the next step if D = chi = 0 is needed in such phases.
- Default kept `.FALSE.`; the other 2D users (Harris sheet example, MUG-TokaMaker) were not re-run with upwinding.
- With `mfnk=T`, the GMRES iteration count in the linear phase rises ~15% with upwinding (the preconditioner is rebuilt only
  every `pre_freq` steps while `V_ad` changes every step); not investigated further.
- Doxygen was not available; the example page was checked with `generate_doc.py` only.
