# Phase 0: xMHD baseline profile (quadrature-point caching project)

Baseline timing and memory for the 3D xMHD module (`src/physics/xmhd.F90`) on `main`
(f182988), before any quadrature-point caching. Measured with the harness in this directory
(`run_bench.py`; see `README.md`). Raw data: `results/baseline.{json,csv}`, `results/kernels.json`,
`results/perf_*.json` (perf reports are in the gitignored `runs/` tree; summaries below).

## Hardware and build

| | |
|---|---|
| CPU | Intel Xeon @ 2.10 GHz (cloud VM), 4 cores (`nproc`=4), 15 GB RAM |
| OS | Ubuntu 24.04 container, Linux 6.18 |
| Compiler | GCC/gfortran 13 (Ubuntu), `-O2` Release build with LTO (OFT default flags) |
| Libraries | Ubuntu OpenBLAS, SuperLU 6, ARPACK, HDF5 1.10, METIS; OpenMPI 4.1 (`OFT_USE_MPI=ON`, OpenMP on) |
| Build | `builds/` flat tree following the `oft-build` skill's system-library recipe, `OFT_BUILD_TESTS=ON` |

All timings are wall-clock times. Unless noted, runs use 1 MPI rank with 4 OpenMP threads, and nothing
else was running on the machine. Repeated runs (eg. `tet_p3_mg` vs `tet_p3_mg_mpiref`/`kernels`) agree to within about 5%.

## Benchmark problem

`xmhd_bench` (`src/tests/physics/MUG/xmhd_bench.F90`) projects an analytic nonlinear state: a 0.05 T guide
field plus a 0.01 T ABC perturbation, a 2e4 m/s sheared flow, and n = 1e19 m^-3 and Ti = Te = 20 eV with
10% modulation. Physics options (`&xmhd_options` written by `run_bench.py`):
two-temperature, advection of all fields, upwinding, kinematic viscosity (`nu_visc=10`) and parallel
viscosity (`nu_par=100`), anisotropic conduction (`kappa_par=1000`, `kappa_perp=1`), density diffusion,
Ohmic and viscous heating, `eta=1`, Hall off, `dt=2e-7`. Meshes:

- **tet**: built-in cube generator (`cad_type=92`), triply periodic box 1x1x2 m, `ni` = 8 (12288 tets) or 4 (1536 tets).
- **hex**: same box, hexes, `ni=4` (128 hexes).
- **curved**: `examples/MUG/Spheromak_tilt/cyl_tilt.h5` with `grid_order=2` (16305 quadratic tets; non-periodic, `bbc='bc'`, `vbc='all'`).

p-multigrid ("`_mg`") uses FE `minlev` = the p=2 level on the finest mesh, smoother = GMRES(2) + block-Jacobi.
Iterations per level follow the MUG tests (`0,2,2` for p=3, `0,10,4,2` for p=4).
**Deviation from the plan:** coarsening to p=1 (`minlev=mg_mesh%mgdim`) gives NaNs in the MG
preconditioner, and the MUG tests never use it. Hex p-MG to p=2 is commented out in `test_alfven.py`, so
hex runs single-level. The local solver is **ILU** (`--pre ilu`, as in the Spheromak examples). SuperLU on
the p=3 run case (86k rows, 25M nnz) needed >8 GB and >10 min to factor on this machine, so LU (as the MUG
tests use) is only available as `--pre lu`.

## Results

### Profile mode (`xmhd_profile`: 3 builds x 5 NL applies x 10 matvecs)

| case | cells | B/V/N DOFs | threads x ranks | build total [s] | fine lvl [s] | coarse lvls [s] | NL apply [s] | matvec [s] | Jac nnz (all lvls) | max RSS [MB] |
|---|---|---|---|---|---|---|---|---|---|---|
| tet_p2 | 12288 | 118800/49150/16380 | 4x1 | 1.107 | 1.107 | - | 0.0548 | 0.0229 | 3.77e+07 | 788 |
| tet_p3 | 12288 | 305200/165900/55300 | 4x1 | 6.385 | 6.384 | - | 0.2216 | 0.1165 | 1.86e+08 | 3483 |
| tet_p3_mg | 12288 | 305200/165900/55300 | 4x1 | 8.263 | 6.244 | 2.018 | 0.2169 | 0.1110 | 2.24e+08 | 4190 |
| tet_p4 | 1536 | 78340/49150/16380 | 4x1 | 3.248 | 3.248 | - | 0.0880 | 0.0505 | 8.43e+07 | 1556 |
| tet_p4_mg | 1536 | 78340/49150/16380 | 4x1 | 4.942 | 3.268 | 1.673 | 0.0785 | 0.0504 | 1.14e+08 | 2102 |
| hex_p3 | 128 | 15230/10370/3456 | 4x1 | 6.268 | 6.268 | - | 0.2478 | 0.0244 | 3.22e+07 | 687 |
| curved_p3_mg | 16305 | 419500/233800/77950 | 4x1 | 11.230 | 8.442 | 2.788 | 0.3030 | 0.1576 | 2.80e+08 | 5262 |
| tet_p3_mg_omp1 | 12288 | 305200/165900/55300 | 1x1 | 29.160 | 22.085 | 7.076 | 0.7298 | 0.3896 | 2.24e+08 | 4188 |
| tet_p3_mg_mpiref | 1536 | 38140/20740/6912 | 4x1 | 1.096 | 0.842 | 0.254 | 0.0284 | 0.0155 | 2.99e+07 | 601 |
| tet_p3_mg_mpi2 | 1536 | 38140/20740/6912 | 2x2 | 0.953 | 0.722 | 0.231 | 0.0250 | 0.0179 | 1.48e+07 (rank 0) | 353 (one rank) |
| tet_p3_mg_mpi2_bound | 1536 | 38140/20740/6912 | 2x2 | 1.902 | 1.446 | 0.456 | 0.0462 | 0.0348 | 1.48e+07 (rank 0) | 353 (one rank) |

`build total` is the average `xmhd_set_ops` time over all MG levels. "Jac nnz" is the sum of the
assembled Jacobian nonzeros over the MG levels; memory per level is about 12 B/nnz. `tet_p3_mg_mpi2_bound` used
OpenMPI's default `--bind-to core` for np<=2, which puts both OpenMP threads of a rank on one core.
The harness now passes `--bind-to none` for OpenMPI automatically. Max RSS for MPI runs is one process (GNU time).

Per-level builds (from `XMHD_PROF build`). All levels are on the finest mesh (p-MG only):

| case | level (p) | avg build [s] | rows | nnz |
|---|---|---|---|---|
| tet_p3_mg | 3 (p3) / 2 (p2) | 6.24 / 2.02 | 664k / 230k | 1.86e8 / 3.77e7 |
| tet_p4_mg | 4 (p4) / 3 (p3) / 2 (p2) | 3.27 / 1.23 / 0.44 | 188k / 87k / 31k | 8.4e7 / 2.5e7 / 5.1e6 |
| curved_p3_mg | 3 (p3) / 2 (p2) | 8.44 / 2.79 | 887k / 308k | 2.3e8 / 4.7e7 |
| tet_p3_mg_omp1 | 3 (p3) / 2 (p2) | 22.1 / 7.08 | | |

### Time stepping (`mode='run'`, 5 steps, 1536-tet mesh, ILU)

| case | total [s] | step 1 [s] | steps 2-5 avg [s] | lin its | NL its | set_ops [s] (calls) | NL apply [s] (calls) | build+NL share | max RSS [MB] |
|---|---|---|---|---|---|---|---|---|---|
| tet_p2_run (single level) | 15.6 | 5.7 | 2.4 | 56 | 18 | 1.19 (7) | 0.34 (29) | 10% | 406 |
| tet_p3_mg_run (assembled NK) | 96.2 | 33.0 | 15.3 | 217 | 16 | 6.78 (6) | 0.97 (27) | 8% | 1872 |
| tet_p3_mg_mfnk_run (matrix-free) | 88.4 | 29.6 | 14.2 | 164 | 10 | 5.56 (5) | 7.52 (210) | 15% | 1875 |

Step 1 includes the ILU factorizations. `xmhd_pre%update` itself costs ~0 s
(`XMHD_PROF preupdate`) because the factorization is done lazily on first apply. Perf on the two p=3 runs
attributes **~62% (NK) / ~50% (MFNK)** of samples to the linear solver: the serial ILU triangular solves
(`ilusolver_apply`, 29%), CSR matvecs (`mat_apply_vec`, 20-31%) and ILU setup (`native_matrix_setup_full`,
7%). Another ~12% is OpenMP wait time, mostly idle threads during the serial ILU apply. Jacobian builds are
~7-8% and the NL residual 1% (NK) / ~15% (MFNK, including its share of interp/basis).

### Isolated basis/geometry kernels (`xmhd_prof_kernels`, profile mode)

Same call sequence as one fine-level `xmhd_build_ops` pass (all cells x quad points), with no physics or
assembly. The thread-averaged time is directly comparable to the wall time of one build or NL apply. NL apply
times are the averages from the baseline table.

| case | qp/cell | fine build [s] | NL apply [s] | `mesh%jacobian`+cgops [s] | basis lag / h1 / hcurl [s] | basis total [s] | `log2phys` [s] | basis % fine build | basis % NL | log2phys % NL | jacobian % NL |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tet_p2 | 14 | 1.10 | 0.0548 | 0.0014 | 0.0035 / 0.0065 / 0.0062 | 0.016 | 0.0023 | 1.5% | 30% | 4.3% | 2.6% |
| tet_p3_mg | 35 | 6.47 | 0.2169 | 0.0034 | 0.0159 / 0.0189 / 0.0267 | 0.062 | 0.0062 | 1.0% | 28% | 2.8% | 1.6% |
| tet_p4_mg | 61 | 3.29 | 0.0785 | 0.0007 | 0.0057 / 0.0065 / 0.0117 | 0.024 | 0.0013 | 0.7% | 30% | 1.7% | 0.9% |
| hex_p3 | 512 | 6.00 | 0.2478 | 0.0112 | 0.0303 / 0.0528 / 0.0738 | 0.157 | 0.0022 | 2.6% | 63% | 0.9% | 4.5% |
| curved_p3_mg | 35 | 8.03 | 0.3030 | 0.0218 | 0.0204 / 0.0253 / 0.0354 | 0.081 | 0.0110 | 1.0% | 27% | 3.6% | 7.2% |
| tet_p3_mg_omp1 (1 thr) | 35 | 22.1 | 0.7298 | 0.0127 | 0.0598 / 0.0727 / 0.1011 | 0.234 | 0.0237 | 1.1% | 32% | 3.2% | 1.7% |

The kernel timer calls `system_clock` 5 times per quadrature point, so the numbers are slight overestimates.

### Sampling profiles (perf, `perf_runs.sh`, `perf_summary.py`)

`perf record -F 499` on the serial process, mesh scaled by `--size 0.5` (tet p3/p4: 1536 cells).
Self-time shares by category. Profile-mode runs also contain setup: L2 projections with CG mass solves
show up as "linear solver", and setup sorts/graph building as "other". These shares are therefore
lower bounds for the build/residual categories. OFT is built with LTO, so some small routines are inlined into
`xmhd_build_ops`/`xmhd_errmatrix_apply`. Exact basis costs come from the kernel table above.

Jacobian builds only (`--prof-reps 5,0,0`), 4 threads:

| category | tet_p2 | tet_p3_mg | tet_p4_mg | hex_p3 | curved_p3_mg | tet_p3_mg (1 thread) |
|---|---|---|---|---|---|---|
| `xmhd_build_ops` self (pointwise physics + local matrix outer products) | 26% | 46% | 52% | 77% | 46% | 41% |
| local->global scatter (`mat_add_values` + CSR column `search4`) | 23% | 20% | 14% | 3% | 23% | 19% |
| OpenMP lock/critical/barrier wait (libgomp) | 24% | 9% | 11% | 2% | 3% | - |
| `xmhd_interp_apply` (field reconstruction) | 1.0% | 1.7% | 1.2% | 0.9% | 1.8% | 1.6% |
| basis eval functions (visible) | 1.1% | 1.3% | 1.1% | 3.1% | 0.9% | 1.3% |
| `mesh%jacobian`, `log2phys` | <0.3% | <0.3% | <0.3% | <0.3% | 0.4% | <0.3% |
| setup (projection solves, sorts, getmop) | ~13% | ~15% | ~14% | ~7% | ~19% | ~31% |

NL residual dominated (`--prof-reps 1,100,0`; still includes one build and setup):

| category | tet_p2 | tet_p3 | tet_p4 | hex_p3 | curved_p3_mg | tet_p3_mg (1 thread) |
|---|---|---|---|---|---|---|
| `xmhd_errmatrix_apply` self | 10% | 12% | 8% | 6% | 11% | 6% |
| `xmhd_interp_apply` (gather + field reconstruction) | 12% | 15% | 11% | 8% | 13% | 9% |
| basis eval (lag/h1/hcurl, incl. `lag_poly` helpers for hex) | 10% | 12% | 10% | ~25% | 12% | 8% |
| `log2phys` / `mesh%jacobian` | <0.3% / - | <0.3% / - | 0.3% / - | 0.4% / 1.4% | 0.9% / 2.1% | 0.3% / - |
| OpenMP lock/critical wait (per-cell `!$omp critical(xmhderr_red)` + build locks) | 30% | 12% | 19% | 6% | 4% | - |
| `global_dp_r8` (residual stitching) | 3.6% | 3.2% | 1.8% | - | - | 2.0% |

## Conclusions

1. **Where the time goes depends on the mode.**
   - *Profile mode / per operation:* one fine Jacobian build costs ~20-30 NL residual applies (tet: 1.1 s vs
     0.055 s at p2, 6.4 s vs 0.22 s at p3) and ~55 matvecs.
   - *Time stepping (ILU-preconditioned GMRES with p-MG, Jacobian rebuilt rarely, `xmhd_prefreq=100`):* the
     linear solve dominates (85-92% of wall time). Jacobian builds are ~7% and NL residuals 1% (assembled
     NK) or ~8.5% (MFNK, 210 applies in 5 steps). Matrix-free NK took fewer linear iterations here and was
     8% faster overall. Faster quadrature alone therefore helps time stepping only modestly with this
     preconditioner. It matters more where builds are frequent (`xmhd_prefreq` small, nonlinear
     Jacobian updates, `xmhd_lin_run`), for MFNK, and in later GPU-resident solves.
2. **Coarse p-level builds are expensive:** +32% (tet p3) and +51% (tet p4) on top of the fine build;
   2.0 s vs 6.2 s for p2 under p3. Reconstructing the fine state on the fly (`xmhd_interp` without the
   level cache) is *not* the cause: the fine-basis evaluation it performs costs ~0.06 s per pass (the
   kernel table), ~3% of a coarse build. The cause is that `xmhd_build_ops` always uses the **finest-level
   quadrature** (`quad=>xmhd_ML_hcurl%levels(nlevels)%fe%quad`). The p2 level under p3 is integrated with
   35 points instead of 14, so the p2 build takes 2.02 s vs 1.11 s for a standalone p2 build on the same mesh.
   Lowering the coarse-level quadrature order is a large win for MG builds, but it changes the preconditioner,
   so convergence must be checked.
3. **Basis evaluation is ~1-3% of a Jacobian build but ~30% of an NL residual on tets (63% on hex).**
   `mesh%jacobian` is negligible on straight tets (once per cell) and ~7% of the NL apply on curved cells
   (per quadrature point). Per build, the cost is the O(ndof^2 x nqp) local matrix arithmetic inside
   `xmhd_build_ops` (40-77% of build samples) plus the global scatter. A quadrature-point basis cache
   therefore mostly speeds up the **residual** (and MFNK Krylov iterations). The Jacobian build needs a
   different approach (eg. batching the local outer products into BLAS-3/GEMM over quadrature points,
   which a `(dof, comp, qp, cell)` layout enables).
4. **`log2phys` at every quadrature point** costs 1-4% of an NL apply (3.6% curved, 4.3% p2) and ~0.1% of a
   build. It is worth the cheap Phase 1 fix, but it is not a big lever.
5. **Scatter and locking are significant.**
   - Build: `mat_add_values` + binary column search (`oft_sort::search4`) take ~15-23% of build samples on
     tets, even single-threaded. The cost is the CSR lookups themselves, not only contention. Precomputing
     per-cell CSR positions (a scatter map, cacheable alongside the quadrature data) would remove the
     search entirely. Lock waits add 9-24% with 4 threads; p2 is worst because there is less work per cell.
   - Residual: the per-cell `!$omp critical` scatter plus other waits take 12-30% of NL-dominated samples
     (worst at p2). Replacing it with atomics or coloring (Phase 1) should give 10-25% on the residual.
   - Thread scaling 1->4 threads: build 3.5x, NL apply 3.4x, matvec 3.5x.
6. **Hex cells are ~90x more expensive per cell** than tets at p3: 512 tensor quadrature points vs 35, and
   64/125/88 vs 20/35/29 Lagrange/H1/H(curl) DOFs per cell. Basis evaluation is 63% of the hex residual.
   Hex meshes gain the most from basis caching, and their tensor structure (sum factorization)
   is the long-term fix.
7. **Memory.** RSS is dominated by the assembled Jacobians (12 B/nnz). tet p3 at 12288 cells: 186M nnz fine + 38M
   p2, 4.2 GB RSS with MG. The ILU factors roughly double the matrix storage (run cases: ~1.9 GB at 1536
   cells vs 0.6 GB in profile mode). Estimates for a FULL physical-basis table (doubles per cell per quad
   point = 4·n_lag + 3·n_h1 + 6·n_hcurl, ie. values+gradients, grad(H1), H(curl) values+curls):

   | case | doubles/qp | qp | per cell | total (this mesh) | vs Jacobian storage |
   |---|---|---|---|---|---|
   | tet p2 (12288; 10/20/14 DOFs) | 184 | 14 | 21 KB | 0.25 GB | 0.45 GB |
   | tet p3 (12288; 20/35/29 DOFs) | 359 | 35 | 100 KB | 1.2 GB | 2.2 GB (fine) |
   | tet p4 (1536; 35/56/53 DOFs) | 626 | 61 | 305 KB | 0.47 GB | 1.0 GB (fine) |
   | hex p3 (128; 64/125/88 DOFs) | 1159 | 512 | 4.7 MB | 0.61 GB | 0.39 GB |

   (DOFs per cell are Lagrange / H1 (order p+1) / H(curl) (OFT's gradient-free H(curl); gradients are in H1).)

   FULL mode is the same order as the Jacobian itself, and larger for hex. REFERENCE mode (orientation-indexed
   reference tables plus per-cell `gop(3,4)`/`det`) costs a few MB plus ~100 B/cell (per qp only for curved cells).
   It is the only option that scales, especially for hex and GPUs.

## Notes for Phases 1-2

- **A/B workflow.** `python run_bench.py --build-dir <b> --tag <x> [--xmhd-opt xmhd_qcache=N]`, then
  `python run_bench.py --compare results/baseline.json results/<x>.json`. Correctness checksums
  (`nlop_norm`, `matvec_norm` in profile mode; `final_norm_*` and iteration counts in run mode) should agree to
  round-off. Expect about 1e-12 relative from OpenMP reduction order. Use `--matrix quick` for a 1-minute smoke test.
- **Kernel targets.** `XMHD_PROF kernels_s` gives the time to *evaluate* the basis per pass (tet p3: 0.062 s,
  hex p3: 0.157 s). A cached FULL/REFERENCE lookup should cut that to memory-bandwidth time. The kernel
  routine (`xmhd_prof_kernels`) can be pointed at the cache in Phase 2 for a direct comparison.
- **Coarse levels.** Fine-basis on-the-fly evaluation for coarse p-levels is small (~3% of a coarse
  build). The coarse-build cost is fine-level quadrature with coarse basis matrices. Caching coarse-level
  basis tables at the fine quadrature points would also cover the coarse-level `oft_*_eval_all` calls.
- **Big levers, if speed is the goal beyond caching.** (a) Batched local-matrix formation (GEMM over quad
  points) in `xmhd_build_ops`. (b) Precomputed CSR scatter positions instead of `search4`. (c) Coarse-level
  quadrature order. (d) For time stepping, the preconditioner: serial ILU apply (`nlocal=1`) leaves 3
  of 4 threads idle. `nlocal=-1` or `xmhd_nparts` would change this baseline.
- **Pitfalls.**
  - Tet cube meshes require power-of-2 `ni`.
  - p-MG to p=1 gives NaNs.
  - SuperLU at p>=3 is impractical here.
  - Run cases use the 1536-cell mesh.
  - Profile mode calls `oft_finalize` (the program ends there).
  - Unknown `xmhd_options` entries are fatal namelist errors, so `--xmhd-opt` must match the build.

## Instrumentation added to `xmhd.F90` (no numerical effect)

- Public `xmhd_prof_*` counters: per-level build time/count/order/finemesh flag in `xmhd_set_ops`; totals for
  `xmhd_set_ops`, `xmhd_errmatrix_apply` and `xmhd_pre%update` (in `xmhd_run`). `xmhd_prof_report` prints them,
  with Jacobian rows/nnz per level.
- `debug>=1` prints the build time per level on every `xmhd_set_ops` call.
- Profile mode:
  - repetitions are configurable through `xmhd_prof_reps` (default `10,10,10` = old behavior; 0 skips NL applies/matvecs);
  - prints the `nlop_norm`/`matvec_norm` checksums;
  - runs `xmhd_prof_kernels`;
  - **fixes** the `Build =` line, which printed `matvec_time/10` instead of the build time.
- Verified with `test_alfven.py` and `test_sound.py`: see "Regression tests" below.
