# xMHD quadrature-cache benchmark harness

Benchmark harness for the xMHD (`src/physics/xmhd.F90`) quadrature-point caching work.
Phase 0 used it to record the baseline (`phase0_baseline.md`, `results/baseline.json`).
Later phases rerun it unchanged to get A/B comparisons.

## Contents

| File | Purpose |
|---|---|
| `run_bench.py` | Writes the inputs, runs the case matrix, parses the output, writes `results/<tag>.json` + `.csv`, and compares two result files |
| `../../src/tests/physics/MUG/xmhd_bench.F90` | Benchmark driver (Fortran), built as `<build_dir>/tests/physics/MUG/xmhd_bench` |
| `results/` | Committed result files (`baseline.*`, `perf_*.txt`) |
| `runs/` | Per-case run directories with inputs and `output.log` (gitignored) |
| `phase0_baseline.md` | Phase 0 report |
| `phase1_results.md` | Phase 1 report (conditional `log2phys`, partition-aware scatter) |
| `run_qcache_bench.py` | Phase 2 kernel benchmark driver for the quadrature cache (`src/fem/fem_quad_cache.F90`); runs `<build_dir>/tests/fem/fem_qcache_bench` on the Phase 0 meshes |
| `../../src/tests/fem/fem_qcache_bench.F90` | Phase 2 kernel benchmark (direct evaluation vs FULL/REFERENCE cache, basis/reconstruction/residual kernels) |
| `phase2_results.md` | Phase 2 report (quadrature-point cache module, validation, kernel A/B, memory) |
| `phase3_results.md` | Phase 3 report (xMHD NL residual on the cache: design, equivalence matrix, A/B) |
| `qcheck_summary.py` | Phase 3: summarizes the residual equivalence checks (`test_MUG_qcache.py -v -s` log) |
| `../../src/tests/physics/MUG/test_MUG_qcache.py` | Phase 3/4 equivalence tests (`xmhd_bench` in `mode='qcheck'`): NL residual implementations, and (Phase 4) Jacobians of every MG level |
| `phase4_results.md` | Phase 4 report (Jacobian on the quadrature cache, CSR position map, default `xmhd_qcache=2`) |
| `jac_compare.py` | Phase 4: runs the equivalence configurations with `jac_dump` and compares Jacobian dumps entrywise, including across builds (eg. against a build of the original Jacobian) |
| `phase4_ab.sh`, `phase4_tables.py` | Phase 4 A/B runs (original vs new Jacobian, cache/direct, position map on/off, BLAS/hand GEMM, experimental coarse quadrature) and their tables |

The driver lives in the OFT source tree, not here, so the normal CMake build compiles and links it
against the in-tree libraries (it is built whenever `OFT_BUILD_TESTS=ON`; no pytest wrapper uses it).
Building it from this directory would mean copying OFT's link line.

### Driver

`xmhd_bench` builds the FE spaces on any mesh `multigrid_construct` supports. It projects an analytic
nonlinear state onto them: a guide field plus a divergence-free ABC perturbation, a sheared
divergence-free flow, and modulated density, Ti and Te. Then it calls either
- `xmhd_run(..., profile_only=.TRUE.)` (`mode='profile'`), which runs `xmhd_profile`: `prof_reps(1)` Jacobian
  builds (`xmhd_set_ops`), each followed by `prof_reps(2)` NL residual applies, each followed by `prof_reps(3)` Jacobian
  matvecs and 100 dot products; or
- `xmhd_run(...)` (`mode='run'`), for `nsteps` time steps, then prints the profiling counters and the norms of the final fields.

Physics is set in the standard `&xmhd_options` group. The defaults in `run_bench.py` turn on two-temperature,
advection of everything, upwinding, kinematic and parallel viscosity, anisotropic conduction, Ohmic and viscous
heating, and resistivity. Hall is off.

### Instrumentation in `xmhd.F90`

Profiling counters only (public module variables `xmhd_prof_*`); they do not change results:
- per FE level: Jacobian build time and count in `xmhd_set_ops`, with whether the level is on the finest mesh
  (so it interpolates the full fine solution on the fly) and its polynomial order; printed per call under `debug>=1`;
- total time and count for `xmhd_set_ops`, `xmhd_errmatrix_apply` (NL residual), and `xmhd_pre%update` in `xmhd_run`;
- `xmhd_prof_reps` sets the repetition counts of profile mode (the default `10,10,10` keeps the original behavior);
- `xmhd_prof_report` prints all counters as `XMHD_PROF ...` lines, together with rows and nnz of every Jacobian
  level. Profile mode calls it, and so does the driver after a run;
- profile mode also prints `XMHD_PROF nlop_norm` / `matvec_norm`, the 2-norms of `A(u)` and `J u` at the
  initial state. Use these as correctness checksums for A/B runs.
- profile mode's `Build =` line used to print `matvec_time/10`; it now prints the average build time.

## Running

```bash
# Build OFT with tests (see the oft-build skill), then:
cd benchmarks/xmhd_qcache
python run_bench.py --build-dir <build_dir> --tag baseline        # full matrix (~15-20 min on 4 cores)
python run_bench.py --build-dir <build_dir> --matrix quick --tag smoke   # smoke test (seconds)
python run_bench.py --list                                          # list cases
python run_bench.py --build-dir <build_dir> --cases 'tet_p3.*' --tag x   # subset (regex per comma item)
```

Useful options:
- `--threads N`: OpenMP threads for cases that don't fix their own (default: all cores). `tet_p3_mg_omp1` always uses 1 thread.
- `--size F`: scale the cube mesh resolution, eg. 1.5 on a bigger machine.
- `--pre {ilu,lu}`: local solver in the block-Jacobi preconditioner/MG smoother. The default `ilu` follows the
  Spheromak examples. `lu` (SuperLU, as in the MUG tests) needs >8 GB and >10 min to factor the p=3 run case on 4 cores.
- `--perf`: wrap serial cases in `perf record` and write `perf_report.txt` (self time) and
  `perf_children.txt` next to `output.log`. Pass `--perf-exe` if `perf` is not on `PATH` (eg.
  `/usr/lib/linux-tools/*/perf` on Ubuntu with a non-matching kernel).
- `--debug 1`: OFT debug output, which includes per-level Jacobian build times for every call.
- `--bench-opt jac_blas=F` (Phase 4+): form local Jacobian matrices with the hand-written loop instead of BLAS `dgemm`.
- `--xmhd-opt xmhd_jac_posmap=0` (Phase 4+): assemble the Jacobian by searching the CSR column lists (no position map).
- `--bench-opt shaping=T`: set all four xMHD transport shaping functions (`res/visc/diff/cond_shaping`) to a
  smooth position-dependent function, to exercise the shaping code paths (Phase 1+).
- `--mpiexec CMD`: MPI launcher for the `nprocs=2` case.

Results go to `results/<tag>.json` (all fields, plus `levels` and per-step `steps` arrays) and to
`results/<tag>.csv` (flat main metrics). Key fields:

| Field | Meaning |
|---|---|
| `build_s` | Average `xmhd_set_ops` time (all levels) in profile mode |
| `build_fine_avg_s` / `build_coarse_avg_s` | Average build time of the finest level / sum over coarser p-levels |
| `levels[]` | Per level: order, finemesh, ncalls, total/avg build time, Jacobian rows/nnz |
| `nlop_s`, `matvec_s`, `dot_s` | Average time per NL residual apply, assembled Jacobian matvec, dot product |
| `nlop_norm`, `matvec_norm` | Checksums (must match between A/B runs up to round-off) |
| `prof_nlop_loop_total_s` | Time in the threaded cell loop of the NL residual (profile/run counters; Phase 1+) |
| `steps[]`, `step_time_avg_s`, `step_time_avg_excl1_s`, `lits_total`, `nlits_total` | Run mode per-step data (wall time, linear/NL iterations) |
| `prof_setops_*`, `prof_nlop_*`, `prof_preupdate_*` | Counters accumulated over the whole run (run mode) |
| `final_norm_B/V/Ti` | Norms of the final fields (run mode checksums) |
| `levels[].qc_builds`, `levels[].posmap_mb`, `posmap_mem_mb`, `posmap_setup_s`, `jcache_coarse_mem_mb` | Phase 4+: Jacobian builds on the quadrature cache per level, CSR position map memory/setup, memory of the coarse-level Jacobian caches |
| `max_rss_mb`, `wall_s` | From `/usr/bin/time -v` (with MPI, `max_rss_mb` covers the launcher process tree as GNU time reports it) |

## Phase 2 kernel benchmark

```bash
python run_qcache_bench.py --build-dir <build_dir> --tag phase2_kernels_t4 --threads 4
python run_qcache_bench.py --build-dir <build_dir> --tag phase2_kernels_t1 --threads 1 --cases tet_p3,hex_p3
python run_qcache_bench.py --report results/phase2_kernels_t4.json    # print tables again
```

Each case times one pass over all cells and quadrature points for six variants (direct evaluation,
FULL/REFERENCE per-qp kernels, REFERENCE canonical frame, FULL/REFERENCE cell-batched kernels) and
three kernels (basis only, basis + xMHD field reconstruction, reconstruction + test-function
accumulation), interleaved over `--nreps` repetitions. Results are merged into `results/<tag>.json`
after every case; run directories go to `runs/qcache_<tag>/`.

## A/B comparisons in later phases

1. Build the branch under test. Use a separate build tree, or rebuild the same one between runs.
2. Run the same matrix with a different `--tag`, passing any new `xmhd_options` through `--xmhd-opt`, eg.
   ```bash
   python run_bench.py --build-dir <build_dir> --tag qcache_off  --xmhd-opt xmhd_qcache=0
   python run_bench.py --build-dir <build_dir> --tag qcache_full --xmhd-opt xmhd_qcache=1
   python run_bench.py --build-dir <build_dir> --tag qcache_ref  --xmhd-opt xmhd_qcache=2
   ```
   `--xmhd-opt` adds raw lines to `&xmhd_options` for every case and can be repeated. `--bench-opt` does the same for
   `&xmhd_bench_options`. On `main` an unknown `xmhd_options` entry is a namelist read error, so only pass
   options the build under test understands.
3. Compare:
   ```bash
   python run_bench.py --compare results/baseline.json results/qcache_full.json
   ```
   This prints A/B time ratios for the build, NL apply, matvec, step time and RSS of each case, the relative
   differences of the checksums (`nlop_norm`, `matvec_norm`, `final_norm_*`), and the iteration counts. With
   OpenMP, results differ at round-off level between runs, because of the order of the reductions/critical
   sections, so expect relative differences around 1e-12–1e-10, not exact equality.
4. Run on an otherwise idle machine and record `nproc`/CPU, which are stored in `meta` in the JSON. Compare
   numbers only from the same machine.
