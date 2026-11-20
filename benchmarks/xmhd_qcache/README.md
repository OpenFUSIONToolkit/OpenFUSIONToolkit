# xMHD benchmark harness

Timing harness for the 3D xMHD (MUG) module (`src/physics/xmhd.F90`): Jacobian builds, the nonlinear
residual, Jacobian matvecs and full time steps, on a fixed matrix of meshes, polynomial orders, multigrid,
threading and MPI configurations. It was used to develop and validate the quadrature-point cache
(`xmhd_qcache`) and the CSR position map (`xmhd_jac_posmap`); the design, validation and the reference results
are summarized in `RESULTS.md`.

## Contents

| File | Purpose |
|---|---|
| `run_bench.py` | Writes the inputs, runs the case matrix, parses the output, writes `results/<tag>.json` + `.csv`, and compares two result files |
| `RESULTS.md` | Design and validation summary, and the reference A/B results (4-core Xeon VM) |
| `../../src/tests/physics/MUG/xmhd_bench.F90` | Benchmark driver (Fortran), built as `<build_dir>/tests/physics/MUG/xmhd_bench` whenever `OFT_BUILD_TESTS=ON` |
| `results/`, `runs/` | Output (gitignored): result files, and per-case run directories with inputs and `output.log` |

The driver lives in the OFT source tree so that the normal CMake build compiles and links it. It is not a
regression test by itself, but `src/tests/physics/MUG/test_MUG_qcache.py` uses it (`mode='qcheck'`) for the
equivalence tests of the residual and Jacobian implementations.

### Driver

`xmhd_bench` builds the FE spaces on any mesh `multigrid_construct` supports. It projects an analytic
nonlinear state onto them: a guide field plus a divergence-free ABC perturbation, a sheared divergence-free
flow, and modulated density, Ti and Te. Then it calls
- `xmhd_run(..., profile_only=.TRUE.)` (`mode='profile'`): `prof_reps(1)` Jacobian builds (`xmhd_set_ops`), each
  followed by `prof_reps(2)` NL residual applies, each followed by `prof_reps(3)` Jacobian matvecs and 100 dot
  products; or
- `xmhd_run(...)` (`mode='run'`): `nsteps` time steps, then it prints the profiling counters and the norms of the
  final fields; or
- `xmhd_run(..., qcache_check=.TRUE.)` (`mode='qcheck'`): the equivalence check used by the tests.

Physics is set in the standard `&xmhd_options` group. The defaults in `run_bench.py` turn on two temperatures,
advection of everything, upwinding, kinematic and parallel viscosity, anisotropic conduction, Ohmic and viscous
heating, and resistivity. Hall is off.

Meshes: the built-in periodic cube generator (tets with 8 or 4 boxes per side, hexes with 4), and the curved
quadratic tet cylinder of `src/examples/MUG/Spheromak_tilt/cyl_tilt.h5`. p-multigrid cases coarsen to p=2 on the
finest mesh with GMRES(2) + block-Jacobi smoothing. The local solver is ILU by default (`--pre lu` for SuperLU).

### Profiling output of `xmhd.F90`

Profile mode and the driver print `XMHD_PROF ...` lines (`xmhd_prof_report`), which `run_bench.py` parses. The
counters only measure; they do not change results:
- per level: Jacobian build count and time, polynomial order, whether it is on the finest mesh, rows and
  non-zeros, the number of builds on the quadrature cache, position-map memory and setup time;
- totals for `xmhd_set_ops`, `xmhd_errmatrix_apply` (NL residual) and preconditioner updates;
- quadrature-cache setup time and memory (`XMHD_PROF qcache`, with the `xmhd_qcache` value) and coarse-level
  cache memory;
- in profile mode, `XMHD_PROF nlop_norm` / `matvec_norm`: the 2-norms of `A(u)` and `J u` at the initial
  state. Use them as checksums when comparing runs (they agree to round-off between options and builds).

`xmhd_prof_reps` (public module variable, set from `prof_reps` in `&xmhd_bench_options`) sets the repetition
counts of profile mode; the default `10,10,10` is the behavior of profile mode without the driver.

## Running

```bash
# Build OFT with tests (see the oft-build skill), then:
cd benchmarks/xmhd_qcache
python run_bench.py --list                                               # list cases
python run_bench.py --build-dir <build_dir> --matrix quick --tag smoke   # smoke test (~1 min)
python run_bench.py --build-dir <build_dir> --tag cache                  # default matrix (~15-20 min on 4 cores)
python run_bench.py --build-dir <build_dir> --cases 'tet_p3.*' --tag x   # subset (regex per comma item)
```

Useful options:
- `--threads N`: OpenMP threads for cases that don't fix their own (default: all cores). `tet_p3_mg_omp1` always
  uses 1 thread and `tet_p3_mg_mpi2` uses 2 ranks with half the cores each.
- `--size F`: scale the cube mesh resolution, eg. 1.5 on a bigger machine.
- `--prof-reps B,N,M`: override the profile-mode repetitions, eg. `2,1,0` (Jacobian builds only) or `1,10,0` (NL
  residual only).
- `--xmhd-opt LINE` / `--bench-opt LINE`: add raw lines to `&xmhd_options` / `&xmhd_bench_options` for every case
  (repeatable), eg. `--xmhd-opt xmhd_qcache=0`, `--bench-opt shaping=T` (sets all four transport shaping
  functions to exercise those code paths).
- `--pre {ilu,lu}`: local solver of the block-Jacobi preconditioner/MG smoother. `lu` needs >8 GB and >10 min to
  factor the p=3 run case on 4 cores.
- `--petsc`: use the PETSc linear algebra backend (`use_petsc=T`; OFT must be built with PETSc). The same XML
  preconditioners are used (block-Jacobi becomes PETSc's additive Schwarz without overlap).
- `--perf` (with `--perf-exe` if `perf` is not on `PATH`): wrap serial cases in `perf record` and write
  `perf_report.txt` (self time) and `perf_children.txt` next to `output.log`.
- `--debug 1`: OFT debug output, which includes per-level Jacobian build times for every call.
- `--mpiexec CMD`: MPI launcher for the 2-rank case (OpenMPI rank binding is disabled automatically).
- `--exe PATH`: run a specific `xmhd_bench` binary (eg. a copy, so that rebuilding does not disturb a run).

Results go to `results/<tag>.json` (all fields, plus `levels` and per-step `steps` arrays) and to
`results/<tag>.csv` (flat main metrics). The JSON `meta` block records the CPU, core count, git SHA and build
type. Key fields:

| Field | Meaning |
|---|---|
| `build_s` | Average `xmhd_set_ops` time (all levels) in profile mode |
| `build_fine_avg_s` / `build_coarse_avg_s` | Average build time of the finest level / sum over coarser p-levels |
| `levels[]` | Per level: order, finemesh, ncalls, total/avg build time, Jacobian rows/nnz, `qc_builds`, `posmap_mb`, `posmap_setup_s` |
| `nlop_s`, `matvec_s`, `dot_s` | Average time per NL residual apply, assembled Jacobian matvec, dot product |
| `nlop_norm`, `matvec_norm` | Checksums (must match between A/B runs up to round-off) |
| `qcache_mode`, `qcache_setup_s`, `qcache_mem_mb` | `xmhd_qcache` value, cache setup time, finest-level cache memory |
| `posmap_mem_mb`, `posmap_setup_s`, `jcache_coarse_mem_mb` | CSR position map memory/setup (all levels), coarse-level cache memory |
| `steps[]`, `step_time_avg_s`, `step_time_avg_excl1_s`, `lits_total`, `nlits_total` | Run mode per-step data (wall time, linear/NL iterations) |
| `prof_setops_*`, `prof_nlop_*`, `prof_preupdate_*` | Counters accumulated over the whole run |
| `final_norm_B/V/Ti` | Norms of the final fields (run mode checksums) |
| `max_rss_mb`, `wall_s` | From GNU `/usr/bin/time -v` (`max_rss_mb` also from BSD/macOS `time -l`; omitted if neither is available, eg. install the `time` package on Debian/Ubuntu). With MPI, `max_rss_mb` covers the launcher process tree as reported by `time` |

## Testing on another platform

1. Build OFT with `OFT_BUILD_TESTS=ON` in Release mode. Check the BLAS: the Jacobian calls `dgemm` from inside
   OpenMP regions, so use an OpenMP BLAS (OFT's `build_libs.py` OpenBLAS, MKL, or the OpenMP variant of a system
   OpenBLAS). With a pthreads OpenBLAS, export `OPENBLAS_NUM_THREADS=1` before running (`run_bench.py` passes the
   environment through). Record which BLAS was used.
2. Run the correctness tests first:
   ```bash
   cd <build_dir>/tests
   OMP_NUM_THREADS=2 ./run_test.sh physics/MUG/test_MUG_qcache.py fem/test_fem_qcache.py
   ```
3. Run the matrix with the default options and with the reference implementation, interleaved if the machine
   is noisy (repeat each a few times with different tags and take medians):
   ```bash
   python run_bench.py --build-dir <build_dir> --tag cache_r1
   python run_bench.py --build-dir <build_dir> --tag direct_r1 --xmhd-opt xmhd_qcache=0 --xmhd-opt xmhd_jac_posmap=0
   ```
   To separate the two optimizations, also run `--xmhd-opt xmhd_qcache=0` alone (direct evaluation with the
   position map) or `--xmhd-opt xmhd_jac_posmap=0` alone (cache with column-search assembly). For single-thread
   numbers use `--threads 1`, and `--prof-reps 2,1,0` (builds) or `--prof-reps 1,10,0` (residual) to shorten the
   profile cases.
4. Compare:
   ```bash
   python run_bench.py --compare results/direct_r1.json results/cache_r1.json
   ```
   This prints A/B ratios (A/B > 1 means B is faster) for the build (all levels, fine, coarse), NL apply, matvec,
   step time and RSS of each case, the relative differences of the checksums (`nlop_norm`, `matvec_norm`,
   `final_norm_*`) and whether the per-step iteration counts match. With OpenMP, results differ at round-off level
   between runs (reduction and atomic order), so expect relative checksum differences of 1e-12-1e-9, not exact
   equality; the iteration counts should match exactly.
5. Compare numbers only from the same machine, run on an otherwise idle machine, and report the CPU, core
   count, compiler, BLAS and MPI with the results. `RESULTS.md` has the reference numbers.
