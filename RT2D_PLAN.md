# Plan: 2D MUG Rayleigh-Taylor example (with and without magnetic field)

Status: draft for review. The work is split between two sequential sub-sessions (A: solver, B: example).

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
see. Session B checks whether the non-linear late phase needs more.

Seed: the m = 1 incompressible eigenmode in velocity only (no sound waves, the same seed works for both models):

    vz = V1 exp(-z/2L) cos(pi z/Lz) cos(kx)
    vx = (V1/k) exp(-z/2L) [cos(pi z/Lz)/(2L) + (pi/Lz) sin(pi z/Lz)] sin(kx)

Because `vx != 0` at the walls, use **free-slip**: Dirichlet `vz` (and `vy`) on the walls, none on `vx`. The
natural boundary condition is consistent, because `Bx Bz = 0` on a conducting wall. `psi`, `n`, `T`, and `by` keep the
default Dirichlet wall conditions (`psi` fixed on the wall gives `dBz = 0` there).

Growth rate: `gamma = (1/2) d ln(ven)/dt` from `xmhd_2d.hist`. With backward Euler the measured value
is biased by about `+gamma dt/2`, so use `timestep_cn=T` (BDF2 in linear mode, Crank-Nicolson in non-linear mode)
and/or `gamma dt ≲ 0.02`.

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
  - Use free-slip wall flags for `vx` (section 2), set after `mhd_sim%setup`.
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
