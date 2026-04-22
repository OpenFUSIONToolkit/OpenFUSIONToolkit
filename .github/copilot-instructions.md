# Copilot Instructions for OpenFUSIONToolkit

## Project Overview

The Open FUSION Toolkit (OFT) is a scientific computing suite for plasma and fusion research. It provides finite element methods on unstructured 2D/3D meshes for MHD equilibrium, stability, and time-dependent simulations.

The codebase is a hybrid Fortran/C/C++/Python project:
- **Core computational code**: Fortran 90 (`.F90` files) with C/C++ bridge files
- **Python interface**: `ctypes`-based wrappers calling into compiled Fortran shared libraries
- **Build system**: CMake (requires ≥ 3.27), orchestrated by `src/utilities/build_libs.py`

### Component Tools

| Tool | Purpose | Key Fortran Source | Python Module |
|---|---|---|---|
| **TokaMaker** | Axisymmetric Grad-Shafranov MHD equilibria | `src/physics/grad_shaf*.F90` | `OpenFUSIONToolkit.TokaMaker` |
| **ThinCurr** | Thin-wall eddy current modeling | `src/physics/thin_wall*.F90` | `OpenFUSIONToolkit.ThinCurr` |
| **Marklin** | 3D force-free ideal MHD equilibria | `src/physics/taylor.F90` | `OpenFUSIONToolkit.Marklin` |
| **MUG** | Time-dependent extended MHD | `src/physics/xmhd*.F90` | *(Fortran executables only)* |

## Repository Layout

```
src/
├── base/           # Core runtime (I/O, XML, sorting, stitching)
├── grid/           # Mesh types, readers (Cubit, Gmsh, native, T3D), multigrid
├── fem/            # Finite element basis functions and operators (H1, Hcurl, Lagrange)
├── lin_alg/        # Linear algebra (native solvers, PETSc, ARPACK, SuperLU, UMFPACK)
├── physics/        # Physics modules (Grad-Shafranov, thin-wall, Taylor, xMHD, tracing)
├── bin/            # Standalone Fortran executables
├── python/
│   ├── OpenFUSIONToolkit/         # Python package (ctypes wrappers)
│   │   ├── TokaMaker/            # TokaMaker Python API
│   │   ├── ThinCurr/             # ThinCurr Python API
│   │   └── Marklin/              # Marklin Python API
│   └── wrappers/                 # Fortran-side C-interop wrapper subroutines
├── ext_libs/       # Bundled 3rd-party sources (triangle, minpack, bvls, dlsode)
├── tests/          # Regression tests (pytest-driven, Fortran + Python)
├── examples/       # Jupyter notebook examples per tool
├── utilities/      # Build scripts, code generators, helper tools
├── docs/           # Doxygen documentation sources
├── cmake/          # CMake find-modules for external dependencies
├── include/        # C/Fortran header files
└── CMakeLists.txt  # Top-level CMake configuration
```

## Build System

OFT uses a two-stage build process:

### Stage 1: Build external libraries

```bash
mkdir builds && cd builds
source ../setup_env.sh  # activates Python venv
python ../src/utilities/build_libs.py \
  --nthread=4 --build_umfpack=1 --build_superlu=1 \
  --build_arpack=1 --oft_build_tests=1 --build_mpich=1
```

This downloads and compiles dependencies (OpenBLAS, HDF5, METIS, etc.) and generates `config_cmake.sh` — a shell script containing the full CMake invocation with all paths.

### Stage 2: Configure, build, and install OFT

```bash
# Still in builds/
bash config_cmake.sh        # runs cmake, creates builds/build_release/
cd build_release
make                         # compile
make install                 # install to builds/install_release/
```

### Key CMake options (set via build_libs.py flags)

- `OFT_BUILD_TESTS` — build test executables (enable with `--oft_build_tests=1`)
- `OFT_BUILD_PYTHON` — build Python wrappers (default ON)
- `OFT_USE_MPI` — enable MPI parallelism (set by `--build_mpich=1` or `--build_openmpi=1`)
- `OFT_BUILD_DOCS` — build Doxygen documentation
- `OFT_BUILD_EXAMPLES` — build example programs

### Important Environment Notes

- The project uses a Python virtual environment at `oft_venv/`. **Always** `source setup_env.sh` before running builds or tests.
- The `copilot-setup-steps.yml` workflow mirrors the Ubuntu 24.04 GCC 14 + OpenMP CI configuration and pre-builds external libraries in a cached `builds/` directory. After this workflow runs, the agent environment has:
  - Compilers: `gcc-14`, `g++-14`, `gfortran-14`
  - Python venv with: `pytest`, `numpy`, `scipy`, `h5py`, `matplotlib`, `xarray`
  - Pre-built external libraries in `builds/`
  - OFT compiled and installed in `builds/install_release/`

## Testing

Tests use **pytest** and are organized under `src/tests/` in subdirectories matching the source layout: `base/`, `grid/`, `lin_alg/`, `fem/`, `physics/`.

### Test structure

Each test subdirectory has:
- `test_*.py` — pytest test files
- `test_*.F90` — corresponding Fortran test programs (compiled during build)
- Various data files (`.h5`, `.g`, `.inp`, etc.)

### Running tests

From `builds/build_release/`:

```bash
source ../../setup_env.sh
make test                    # runs: pytest -m "not slow" base grid lin_alg fem physics
make test_full               # runs all tests including slow ones
make test_examples           # runs example notebook tests
```

Or run individual test files:

```bash
cd builds/build_release/tests
../../run_test.sh physics/test_TokaMaker.py -k "test_name"
```

To list tests for a specific file or folder:

```bash
cd builds/build_release/tests
../../run_test.sh physics/test_TokaMaker.py --collect-only
```

### Test conventions

- Tests marked `@pytest.mark.slow` are excluded from default CI runs
- Tests marked `@pytest.mark.mpi` require MPI (`OFT_HAVE_MPI=1`)
- Tests marked `@pytest.mark.coverage` are for code coverage runs
- `oft_testing.py` provides `run_OFT()` helper for running Fortran executables with timeout
- Physics Python tests (TokaMaker, ThinCurr) use `multiprocessing.Process` to isolate OFT runtime (only one `OFT_env` instance per process)
- Test timeout is multiplied by 4× when `OFT_DEBUG_TEST=1`

## Linting

### Python linting

Python code is linted with **ruff**. Configuration is in `src/python/pyproject.toml`:

```bash
cd src/python && ruff check
```

Rules: Pyflakes (`F`) + pycodestyle (`E`) with ignores for `E722`, `F403`, `F405`. Target: Python 3.7.

### Fortran stack checking

A custom lint checks Fortran debug stack entries:

```bash
cd src && python utilities/generate_stack.py -l
```

This validates that all `SUBROUTINE`/`FUNCTION` entries have matching debug stack annotations. Run from the `src/` directory.

## CI Workflows

| Workflow | File | Trigger | Purpose |
|---|---|---|---|
| **CI Build** | `ci_build.yaml` | push to main, PRs | Full matrix build (GCC, Intel, macOS) × (OpenMP, MPICH, OpenMPI) |
| **Lint** | `lint.yaml` | push to main, PRs | ruff check + Fortran stack check |
| **Coverage** | `cov_build.yaml` | push to main | Build with `--coverage`, upload to Codecov |
| **CD Nightly** | `cd_nightly.yaml` | schedule | Nightly package builds |
| **Copilot Setup** | `copilot-setup-steps.yml` | manual/PR | Agent environment setup |

### CI configuration used in copilot-setup-steps

- **OS**: Ubuntu 24.04
- **Compilers**: `gcc-14` / `g++-14` / `gfortran-14`
- **Parallel**: OpenMP only (no MPI)
- **Libraries**: OpenBLAS, UMFPACK, SuperLU, ARPACK, HDF5, METIS

## Coding Conventions

### Fortran

- Free-form Fortran 90+ (`.F90` extension, preprocessed)
- SPDX license header: `! SPDX-License-Identifier: LGPL-3.0-only`
- Doxygen-style comments with `!>` and `!!` markers
- Module names typically prefixed with `oft_` (e.g., `oft_gs`, `oft_la_base`)
- Line length: unlimited (`-ffree-line-length-none`), but keep reasonable
- New functions/subroutines must have Doxygen documentation

### Python

- SPDX license header comment block at top of each file
- Doxygen-style docstrings with `@param`, `@result`, `@authors`, `@date`
- Target Python 3.7+ (for OMFIT compatibility)
- Use `ctypes` for Fortran interop; wrapper patterns in `_interface.py` files
- All Python wrappers live under `src/python/OpenFUSIONToolkit/`

### Pull Requests

- PRs for a specific tool should be titled: `ToolName: description` (e.g., `TokaMaker: Fix boundary condition`)
- Keep changes focused on a single feature; secondary changes should be minimal
- Comment on whether APIs or input files change
- Run regression tests before submitting

## Common Pitfalls

1. **Single OFT_env instance**: The Python `OFT_env` class enforces a singleton. Tests use `multiprocessing.Process` to work around this. Never create two `OFT_env` instances in the same process.

2. **Source setup_env.sh**: Always source this before any build/test commands. It activates the Python venv. In CI, this file is generated during prerequisites setup.

3. **Build from `builds/` directory**: The `build_libs.py` script must run from the `builds/` directory. It creates `config_cmake.sh` there. CMake then creates `build_release/` and `install_release/` inside `builds/`.

4. **Tests run from build tree**: Tests must be run from `builds/build_release/tests/` (or via `make test` from `builds/build_release/`), not from the source tree, because compiled Fortran test executables are in the build tree.

5. **ext_libs/ is bundled third-party code**: Do not modify files in `src/ext_libs/`. These are upstream sources (triangle, minpack, bvls, dlsode).

6. **Fortran/Python interop**: The Python package calls compiled Fortran via `ctypes`. The Fortran-side wrappers are in `src/python/wrappers/` and use `ISO_C_BINDING`. Changes to Fortran function signatures require corresponding updates to both `wrappers/*_f.F90` and `python/OpenFUSIONToolkit/*/_interface.py`.

7. **CMake template files**: Some files use `@VARIABLE@` CMake substitution (e.g., `__init__.py`, `pyproject.toml.in`, `run_test.sh.in`). Edit the `.in` template, not the generated file.
