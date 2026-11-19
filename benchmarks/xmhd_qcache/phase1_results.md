# Phase 1: conditional `log2phys` and partition-aware scatter (xMHD quadrature-caching project)

Branch `claude/xmhd-qcache-phase1`, starting from `claude/xmhd-qcache-phase0` (`a3c6a74`). Phase 1 makes two
low-risk changes to `src/physics/xmhd.F90`, plus one optional, safe change to the Jacobian build. All of them
skip work or synchronization that the results don't need. The baseline is the unmodified Phase 0 branch, rebuilt
and rerun on this machine.

## Code changes

1. **Conditional `log2phys`** (`xmhd_build_ops`, `xmhd_errmatrix_apply`). These are the only two cell loops where
   the physical point `pt` is optional:
   - in the Jacobian build, `pt` is used only by the four shaping functions;
   - in the NL residual, it is used by the shaping functions and by the `self%dt<0` diagnostics branch.

   A new helper, `xmhd_shaping_active()`, tests whether any of `res/visc/diff/cond_shaping` is associated. Each
   loop now computes `pt` only when shaping is active (or when `dt<0` in the residual). In `xmhd_massmatrix_apply`
   and `xmhd_diag`, `pt` feeds the always-on energy/flux diagnostics, so those loops are unchanged.
   `xmhd_prof_kernels` (profile-only instrumentation) is also unchanged.
2. **No per-cell critical in the residual-type loops** (`xmhd_errmatrix_apply`, `xmhd_massmatrix_apply`, the only
   two `critical(xmhderr_red)` users).
   - `xmhd_setup_cell_shared` flags each cell that shares a vertex with a cell in another thread partition
     (`mesh%tloc_c`). The flags are cached and rebuilt only if `mesh` changes.
   - Every DOF of a cell lives on its vertices, edges, faces or interior, so an unflagged cell's DOFs are touched
     only by cells of its own partition. Partitions run on one thread at a time (`!$omp do` over `ip`).
   - The scatter goes through `xmhd_scatter_add(vec,inds,vals,atomic,scale)`, which uses plain adds for unflagged
     cells and `!$omp atomic` adds for flagged cells.
   - With one thread there is one partition and no flagged cells, so the summation order is exactly the old one.
3. **Lock-free Jacobian adds for unflagged cells** (`xmhd_build_ops`). `fem_mat_add_local`
   (`src/fem/fem_composite.F90`) now takes `tlocks` as an optional argument; when it is absent, no locks are
   taken.
   - When the Jacobian is an `oft_native_matrix`, unflagged cells call it without locks. A native `add_values`
     writes only the given rows/columns, and those rows belong to DOFs that no other partition touches.
   - Flagged cells, and all non-native matrices (eg. PETSc), keep the per-field locks.
4. **Instrumentation and benchmark support** (no numerical effect):
   - new counter `xmhd_prof_tnlop_loop`, the time in the threaded cell loop of the residual, printed as
     `XMHD_PROF nlop_loop total_s` and parsed by `run_bench.py` as `prof_nlop_loop_total_s`;
   - new `xmhd_bench` option `shaping=T`, which sets all four shaping functions to a smooth position-dependent
     function. No regression test sets a shaping function, so this was the only way to exercise those paths.

**No new user options.** During the A/B runs a temporary `xmhd_scatter_mode` namelist entry selected between the
variants below. It has been removed, and variant 3 is now hard-coded:
1. atomics for all cells, with the Jacobian build locked as before;
2. atomics only for flagged cells, with the Jacobian build locked as before;
3. variant 2 plus lock-free Jacobian adds for unflagged cells.

Flagged-cell fraction with 4 thread partitions: 4570/12288 (37%, tet cube `ni=8`), 5003/16305 (31%, curved cylinder),
984/1536 (64%, tet cube `ni=4`), 88/128 (69%, hex). The vertex-adjacency test is conservative, and the partitions
are small at 4 threads.

## Machine and method

| | |
|---|---|
| CPU | Intel Xeon @ 2.10 GHz (cloud VM), 4 vCPUs, 15 GB |
| Build | Ubuntu 24.04 system libraries (`oft-build` skill, section 5): GCC 13 `-O2` + LTO, OpenMPI 4.1, OpenBLAS, SuperLU, METIS. Two trees: `builds/base` (Phase 0 `a3c6a74` in a worktree) and `builds/` (this branch) |
| Runs | `run_bench.py`, default matrix, ILU, 4 OpenMP threads unless noted |

**This VM is noisy.** Repeat runs of the *same* binary differ by up to ±15% (eg. base `tet_p4` NL apply
0.109 s vs 0.086 s). Even the 1-thread basis kernels, which have no synchronization, vary by ±10%. On this host
the residual also scales only ~2.3x from 1 to 4 threads, against 3.4x on the Phase 0 machine, and it does so with
or without synchronization (see below). The container was also restarted (onto a different host) during the
first baseline run, so the first baseline was discarded. All numbers below come from the same container, with
base and Phase 1 runs **interleaved case by case** to spread drift evenly. Committed result files hold the
**minimum** over the interleaved repetitions (`nreps` in the JSON).

## A/B results

### Profile mode (`results/phase1_base.json` vs `results/phase1.json`, 3 builds x 5 NL applies x 10 matvecs)

Ratio = Phase 1 / base (<1 is faster).

| case | threads x ranks | reps | build base [s] | build P1 [s] | ratio | NL apply base [s] | NL apply P1 [s] | ratio | nlop_norm rel. diff |
|---|---|---|---|---|---|---|---|---|---|
| tet_p2 | 4x1 | 2 | 1.245 | 1.097 | 0.88 | 0.0755 | 0.0601 | 0.80 | 7e-14 |
| tet_p3 | 4x1 | 2 | 6.842 | 6.523 | 0.95 | 0.2405 | 0.2218 | 0.92 | 4e-15 |
| tet_p3_mg | 4x1 | 2 | 8.827 | 8.672 | 0.98 | 0.2418 | 0.2503 | 1.04 | 2e-15 |
| tet_p4 | 4x1 | 2 | 3.482 | 3.422 | 0.98 | 0.0858 | 0.0922 | 1.08 | 2e-11 |
| tet_p4_mg | 4x1 | 2 | 5.577 | 5.201 | 0.93 | 0.1008 | 0.0863 | 0.86 | 1e-15 |
| hex_p3 | 4x1 | 2 | 6.326 | 6.539 | 1.03 | 0.2753 | 0.2766 | 1.00 | 6e-11 |
| curved_p3_mg | 4x1 | 2 | 11.730 | 11.470 | 0.98 | 0.3737 | 0.3277 | 0.88 | 0 |
| tet_p3_mg_omp1 | 1x1 | 1 | 27.780 | 27.740 | 1.00 | 0.7054 | 0.6613 | 0.94 | identical |
| tet_p2 (1 thread) | 1x1 | 1 | 3.490 | 3.477 | 1.00 | 0.1761 | 0.1739 | 0.99 | identical |
| tet_p4_mg (1 thread) | 1x1 | 1 | 16.910 | 16.930 | 1.00 | 0.2507 | 0.2373 | 0.95 | identical |
| tet_p3_mg_mpiref | 4x1 | 1 | 1.297 | 1.435 | 1.11 | 0.0345 | 0.0404 | 1.17 | 5e-11 |
| tet_p3_mg_mpi2 | 2x2 | 1 | 1.325 | 1.209 | 0.91 | 0.0426 | 0.0379 | 0.89 | 9e-13 |

`matvec_norm` behaves the same way as `nlop_norm`: at most 6e-11 relative, and identical at 1 thread. Base-vs-base
repeats on this machine differ by up to 5e-11 (`hex_p3`, `tet_p4`). Note that `matvec_norm` checks the Jacobian
assembled by the lock-free build.

### Focused NL-apply runs (1 build + 100 NL applies, 3 interleaved repetitions; `results/phase1_focus_nl_*.json`)

Times are the per-apply median (min) over the 3 repetitions, 4 threads.

| case | base (critical) | variant 1 (all atomic) | Phase 1 (variant 3) | P1/base (median) |
|---|---|---|---|---|
| tet_p2 | 0.0769 (0.0768) | 0.0767 (0.0732) | 0.0797 (0.0777) | 1.04 |
| tet_p3 | 0.2827 (0.2781) | 0.2777 (0.2607) | 0.2674 (0.2612) | 0.95 |
| tet_p4 | 0.1113 (0.1017) | 0.1062 (0.1057) | 0.1033 (0.0985) | 0.93 |

The loop timer (`xmhd_prof_tnlop_loop`, 100 applies; 40 for the MG cases) shows that the threaded cell loop is
**97-98% of the NL apply**. Serial setup, BCs and stitching take only 2-10 ms per apply. Loop time per apply
(two repetitions at 4 threads):

| case | threads | base | variant 1 | Phase 1 |
|---|---|---|---|---|
| tet_p2 | 1 / 2 / 4 | 0.174 / 0.120 / 0.075, 0.080 | 0.173 / 0.123 / 0.075, 0.073 | 0.148 / 0.117 / 0.080, 0.073 |
| tet_p4 | 1 / 2 / 4 | 0.232 / 0.165 / 0.113, 0.111 | 0.241 / 0.183 / 0.098, 0.108 | 0.231 / 0.160 / 0.107, 0.105 |
| tet_p3_mg | 4 | 0.303, 0.310 | 0.268, 0.293 | 0.254, 0.295 |
| curved_p3_mg | 4 | 0.401, 0.412 | 0.407, 0.403 | 0.419, 0.367 |

(The 1-thread Phase 1 code path for tet_p2 is the base path with the `log2phys` skipped, so the
0.148 vs 0.174 difference is partly noise.)

**Perf (task-clock, 4 threads, `--prof-reps 1,100,0`)** gives a less noisy view of the residual loop itself
(`xmhd_errmatrix_apply._omp_fn` self + `xmhd_scatter_add`):

| case | base CPU-s | Phase 1 CPU-s | change | base `GOMP_critical_name_start` |
|---|---|---|---|---|
| tet_p2 | 9.18 | 7.10 | -23% | 0.37% of samples |
| tet_p4 | 11.0 | 8.66 | -21% | <0.3% |

So the residual loop does ~21-23% less CPU work. Most of the saving is probably cache-line traffic and
serialization inside the old critical section, which perf books under the loop body, not under libgomp. But on
this host, wall time per apply doesn't move outside the noise: -5% to -7% median on p3/p4, +4% on p2. The explicit
lock wait that Phase 0 measured (12-30% of NL samples) does not appear here: `GOMP_critical_name_start` is 0.4%,
and the libgomp spin is ~4.7% in both builds. Thread scaling of the loop is ~2.2-2.3x at 4 threads both with
and without synchronization. The limiting factor on this VM is therefore hardware/hypervisor (shared cores or
memory bandwidth), not synchronization. Expect a larger wall-clock gain on a machine where the critical section
actually contends, as on the Phase 0 host. Phase 1 is never measurably slower.

Variant 1 (atomics everywhere) costs about the same as variant 3 here. Variant 3 is kept because it needs no
atomics for 30-60% of the cells, and its partition/flag structure is the one Phase 2 and GPU coloring need.

### Jacobian build (5 builds, 3 interleaved repetitions; `results/phase1_focus_bd_*.json`)

| case | base median (min) [s] | variant 2 (locks for all cells) | Phase 1 (no locks for unflagged cells) | P1/base |
|---|---|---|---|---|
| tet_p2 | 1.338 (1.312) | 1.389 (1.378) | 1.254 (1.221) | 0.94 |
| tet_p4_mg | 5.969 (5.642) | 6.080 (6.022) | 6.131 (5.816) | 1.03 |

Lock-free adds for unflagged cells give ~6% on the p2 build, where per-cell work is small and lock traffic is
relatively largest. At p4 the build is dominated by local matrix arithmetic, and the change is within noise.
The skipped `log2phys` is ~0.1% of a build (Phase 0), which is not measurable here.

### Time stepping (`mode='run'`, 5 steps; single runs, base and Phase 1 in different sessions)

Periodic cube (`cad_type=92`, triply periodic), 1536 tets:

| case | step 2-5 avg base [s] | Phase 1 [s] | lin its base/P1 | NL its base/P1 | final_norm_B / V / Ti rel. diff |
|---|---|---|---|---|---|
| tet_p2_run | 1.87 | 2.14 | 56/56 | 18/18 | 5e-12 / 2e-12 / 1e-12 |
| tet_p3_mg_run (NK) | 15.77 | 16.33 | 217/217 | 16/16 | 1e-10 / 1e-12 / 4e-13 |
| tet_p3_mg_mfnk_run | 14.96 | 15.61 | 164/164 | 10/10 | 2e-10 / 1e-11 / 1e-12 |

Iteration counts are identical, and final norms agree to round-off. That includes the Jacobian assembled
lock-free on a periodic mesh, where periodic copies stay distinct local DOFs until stitching. Run-mode time is
85-92% linear solver (Phase 0), so these step times are noise-dominated and not affected by Phase 1.

## Correctness

- **1 thread: bit-identical to base** for every checksum compared: `nlop_norm`, `matvec_norm`, and the run-mode
  final norms and iterations. Cases: tet_p2, tet_p4_mg, tet_p3_mg_omp1, quick run cases.
- **Shaping paths** (`--bench-opt shaping=T`, all four shaping functions set, so `log2phys` is evaluated):
  - 1 thread: base and Phase 1 are bit-identical for tet_p2, tet_p4 and curved_p3_mg (profile mode), and for
    the `quick_tet_p3_mg_run` final norms and iterations;
  - 4 threads: they agree to <=6e-11 (profile) and <=2e-9 (run-mode `final_norm_B` after 2 steps, same
    iterations).
- **4 threads, no shaping:** all checksums agree to <=6e-11, within the base-vs-base spread (up to 5e-11). Run
  mode on the periodic mesh has identical iteration counts.

## Regression tests

From `builds/build_release/tests`, `OMP_NUM_THREADS=2`, with `run_test.sh`:

| test | Phase 1 |
|---|---|
| alfven `test_nl_r1_p2` | PASSED |
| alfven `test_nl_r1_p2_mpi[False-False]` (NK) / `[False-True]` (MFNK) | PASSED / PASSED |
| alfven `test_lin_r1_p2` | PASSED |
| alfven `test_hex_nl_r1_p2` | PASSED |
| sound `test_nl_r1_p2` | PASSED |
| sound `test_nl_r1_p2_mpi[False-False-False]` / `[False-True-False]` | PASSED / PASSED |
| sound `test_lin_r1_p2` | PASSED |
| sound `test_nl_r1_p3` (p-MG) | PASSED |
| `test_MUG_visc.py` (5 tests) | all PASSED |
| `test_MUG_gyro.py` (6 tests) | all PASSED |
| alfven `test_nl_r1_p3` (p-MG) | FAILED (known; same failure on base, see below) |
| alfven `test_lin_r1_p3` (p-MG) | FAILED (known; same failure on base) |

Totals: 21 passed and 2 failed in 7 minutes. Both failures are the known p-MG environment failures from Phase 0:
"NL Solver failed to converge" (`test_nl_r1_p3`) and "Solver did not converge" (`test_lin_r1_p3`). The base build
(`builds/base`, Phase 0 `a3c6a74` plus only the loop timer and driver option) fails both in the same way on this
machine.

## What's left (NL residual, 4 threads, Phase 1 build, perf self time)

`--prof-reps 1,100,0`, tet_p2 / tet_p4 (these runs still include setup and one build):
- `xmhd_interp_apply` (field reconstruction): 16% / 12%
- `xmhd_errmatrix_apply` loop body (pointwise physics + local vector accumulation): 13% / 8%
- basis evaluation (`h1_geval`, `hcurl_eval/ceval`, `lag_geval`): ~14% / ~12%
- `global_dp_r8` stitching: 4.3% / 1.8%
- libgomp spin/barrier: ~4.7% / ~4.4%
- `xmhd_scatter_add`: 1.3% / 0.4%

So the residual is now dominated by field reconstruction and basis evaluation. That is exactly what Phase 2's
quadrature cache targets.

## Recommendations for Phase 2

1. **Reuse the partition flags.** `xmhd_cell_shared` (cells touching another thread partition's DOFs) is the hook
   for both CPU and GPU scatter.
   - It belongs next to the quadrature data, in `src/fem`, as a per-mesh object rather than an `xmhd` module
     variable. The setup routine is self-contained and mesh-only, so it can move as is.
   - On GPUs, replace it with a cell coloring, or keep atomics: they cost the same as variant 3 here.
   - With 4 partitions, 30-70% of cells are flagged. A finer partition (more parts than threads, eg.
     `oft_env%nparts`), or the face/edge-based test, would cut that if atomic cost ever shows up.
2. **Physical points for shaping.** `xmhd_shaping_active()` plus the `need_pt` pattern is the guard. The cache
   should store physical qp coordinates (or the shaping *values*) only when shaping is active. The bench driver's
   `shaping=T` option provides a test for it.
3. **Build scatter.** The remaining build cost is local matrix arithmetic plus CSR column search (`search4` ~4-6%
   of NL-mode samples, ~15-23% of build samples in Phase 0, even single-threaded).
   - Lock removal helped only ~6% at p2. **Precomputed CSR positions** (`loc_cache`, already supported by
     `mat_add_values` but only for blocks with `mloc%ind`) are the bigger lever; storing a per-cell position map
     alongside the qcache would remove both the search and most of the lock hold time for flagged cells.
   - Note that the locks in `fem_mat_add_local` are held across the column search. Even before a position map
     exists, doing the search outside the lock would shorten lock hold times for flagged cells.
4. **Benchmarking.** Wall-clock A/B on shared cloud VMs needs interleaved repetitions. Prefer the loop timer
   (`prof_nlop_loop_total_s`) and perf task-clock CPU seconds over single profile-mode averages. Effects below
   ~10% are not resolvable on this host from 2-3 repetitions.
