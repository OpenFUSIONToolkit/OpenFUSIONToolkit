---
name: oft-test
description: Run, select, and write Open FUSION Toolkit (OFT) regression tests (pytest driving Fortran test programs and the TokaMaker/ThinCurr/Marklin Python APIs). Use when verifying a change, running the test suite or a single test, listing tests, debugging a failing or timing-out test, or adding a new test under src/tests.
---

# Testing OFT

Tests use **pytest** and live in `src/tests/`, in subdirectories mirroring the source
layout: `base/`, `grid/`, `lin_alg/`, `fem/`, `physics/`. Each holds `test_*.py` pytest
files, the `test_*.F90` Fortran programs they drive, and data files (`.h5`, `.g`,
`.inp`, …).

## Before running anything

1. **Tests run from the build tree, never the source tree.** The compiled Fortran test executables, copies of the test scripts and data, and the Python package all live under `<build_dir>`. Find `<build_dir>` with the `oft-build` skill ("Locate the build tree") — do not assume `builds/build_release`, which is only the CI/ephemeral layout. Also check that the build tree compiles *this* checkout (see "Check which source tree it builds" in `oft-build`), or you will be testing different code.
2. **Rebuild first.** Run `make -j<N>` in `<build_dir>` after any source change, including edits to test `.py` files or data: they are copied into the build tree by `configure_file`, so the build-tree copy is stale until `make` runs. `make install` is not needed for tests.
3. **Python environment.** `make test` and `run_test.sh` call the `pytest` recorded at configure time (`OFT_PYTEST_EXE`) by absolute path, so they work without activating anything as long as that environment still exists. They also add that interpreter's directory to `PATH`, but after the build-tree `python/` and the MPI `bin/` directory: with a system MPI in `/usr/bin`, a bare `python` inside a test is the system interpreter, so tests must launch Python with `sys.executable`. If it is gone or broken, or you run `python`/`pytest` directly, you need the user's environment — ask which to use if not already known; in CI/ephemeral environments only, `source setup_env.sh`.
4. **Always set `OMP_NUM_THREADS`.** Use `1` if only one or two cores are available, otherwise `2`, however many cores there are. Leaving it unset lets OpenMP use every core in every test process, which oversubscribes the machine (especially with MPI tests) and makes timings and timeouts unreliable. Set it in the same command that runs the tests, since agent shells may not keep exported variables between calls:

   ```bash
   OMP_NUM_THREADS=$( [ "$(getconf _NPROCESSORS_ONLN)" -le 2 ] && echo 1 || echo 2 ) make test
   ```

   In a Linux container, `getconf` and `nproc` report the CPUs the process may run on, not a CPU quota. Check the quota with `cat /sys/fs/cgroup/cpu.max` (cgroup v2: `max` means no limit, otherwise `<quota> <period>`) or `cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us /sys/fs/cgroup/cpu/cpu.cfs_period_us` (cgroup v1: a quota of `-1` means no limit). If quota ÷ period is 2 or less, use `1`.

   The only exception is running individual tests specifically to study OpenMP thread scaling or performance, where you set the thread count deliberately.

## Running tests

From `<build_dir>`, with `OMP_NUM_THREADS` set as above (shown as `2` in the examples):

```bash
OMP_NUM_THREADS=2 make test            # pytest -m "not slow" base grid lin_alg fem physics
OMP_NUM_THREADS=2 make test_full       # all tests, including slow ones
OMP_NUM_THREADS=2 make test_examples   # example notebooks/programs (needs OFT_BUILD_EXAMPLES or OFT_BUILD_PYTHON)
```

`make test`/`test_full` write `<build_dir>/tests/OFT.junit.xml`. `make test` alone takes
roughly 50 minutes on a 4-core sandbox, and `test_full` longer — run them in the
background, and prefer targeted runs while iterating.

Targeted runs use `run_test.sh`, generated into `<build_dir>/tests/`. It sets `PATH`
(build-tree `python/`, the configured interpreter, and the MPI that OFT was built
against), `OFT_HAVE_MPI`, and `OFT_DEBUG_TEST`, then passes all arguments to pytest:

```bash
cd <build_dir>/tests
OMP_NUM_THREADS=2 ./run_test.sh physics/test_TokaMaker.py                 # one file
OMP_NUM_THREADS=2 ./run_test.sh physics/test_TokaMaker.py -k "test_name"  # matching tests
./run_test.sh physics/test_TokaMaker.py --collect-only                   # list tests (runs nothing)
OMP_NUM_THREADS=2 ./run_test.sh grid -m "not slow"                        # a folder, skipping slow tests
```

To call one test function directly, optionally with arguments `test_args`, from the
build-tree directory containing `test_file.py`:

```bash
OMP_NUM_THREADS=2 python -c "from test_file import test_name; test_name(*test_args)"
```

This bypasses `run_test.sh`, so use the user's Python environment and make sure the
`PATH`/`OFT_HAVE_MPI` that `run_test.sh` would set are not needed (Fortran-driven and
MPI tests need them).

## Markers and environment

- `@pytest.mark.slow` — excluded from `make test` and default CI runs; included by `make test_full`.
- `@pytest.mark.mpi` — parallel tests needing MPI (`OFT_HAVE_MPI=1`, set by `run_test.sh` for MPI builds). The skip itself happens in `run_OFT()`, which calls `pytest.skip` when `nproc > 1` and `OFT_HAVE_MPI` is not `1`, so non-MPI builds report these as skipped rather than failed.
- `@pytest.mark.coverage` — selected by coverage builds (`OFT_COVERAGE`).
- `@pytest.mark.linear` — linear MHD tests.
- Test timeouts are multiplied by 4× when `OFT_DEBUG_TEST=1` (Debug builds).

### Skips

Large skip counts are normal when optional libraries are missing. Before reading skips
as a broken build, check the reason (`-rs` prints a summary). The common ones come from
`run_OFT()` in `oft_testing.py`:

- `Not compiled with PETSc` — the Fortran program reported that PETSc is unavailable. Without PETSc this covers a few hundred cases across `fem/`, `lin_alg/test_solver_xml.py`, and the MUG/Marklin physics tests.
- `Configuration not supported` — the program printed `SKIP TEST`, eg. `lin_alg/test_native_bjacobi.py` cases whose direct solver (SuperLU, SuperLU-DIST, MUMPS, UMFPACK) this build lacks, and `lin_alg/test_arpack.py` without ARPACK.
- `Not compiled with MPI` — a `nproc > 1` case in a non-MPI build.
- A bare `Skipped` — a parameter combination the test itself excludes (eg. matrix-free together with PETSc in the MUG tests).

## Debugging failures

- Fortran-driven tests use `run_OFT()` from `oft_testing.py`, which launches the executable (under `mpirun` when `nproc > 1`) and kills it on timeout. Re-run the failing test alone with `-x -s` to see its output, or run the Fortran executable by hand in the test's build-tree directory.
- A test that only fails under MPI is often picking up the wrong `mpirun`; run it through `run_test.sh`, which puts the MPI OFT was built against first on `PATH`.
- If results look unaffected by your change, confirm `make` was run and the build tree points at this checkout.

## Writing tests

- Add the test file **and** every data file it reads to the subdirectory's `CMakeLists.txt` (`configure_file(<file> ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)`), and new Fortran programs with `oft_add_test(test_foo.F90)`; otherwise they never reach the build tree.
- Use `run_OFT(command, nproc, timeout)` from `oft_testing.py` to run Fortran executables.
- Only one `OFT_env` may exist per process. Physics Python tests (TokaMaker, ThinCurr, Marklin) run each case through an `mp_run()` helper that starts a `multiprocessing.Process` and returns results through a `multiprocessing.Queue` — follow the pattern in the existing test file rather than creating `OFT_env` at module level or twice in one process.
- Mark cases too slow for routine CI as `@pytest.mark.slow`, and parallel cases as `@pytest.mark.mpi`.
- Test functions do not need Doxygen-style documentation (see "Test and example code" in `AGENTS.md`); a short comment saying what is covered is welcome.
