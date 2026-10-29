# AI Agent Instructions for OpenFUSIONToolkit

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
├── lin_alg/        # Linear algebra (native solvers, PETSc, ARPACK, SuperLU, UMFPACK, STRUMPACK)
├── physics/        # Physics modules (Grad-Shafranov, thin-wall, Taylor, xMHD, tracing)
├── bin/            # Standalone Fortran executables
├── python/
│   ├── OpenFUSIONToolkit/        # Python package (ctypes wrappers)
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

OFT uses a two-stage build process. The build tree is created by the user and its
location is **not** fixed by the project — read "Locating the build tree" before
running any build or test command.

### Locating the build tree

These instructions use `<build_dir>` for the configured build directory and
`<config_dir>` for its parent (where `build_libs.py` was run, holding
`config_cmake.sh`, the compiled external libraries, and `install_release/`).
**Do not assume `builds/build_release`.**

- **User environments** (the normal case): users commonly keep several build trees side by side, one per toolchain/parallelism combination, nested under `builds/`. For example `builds/GCC/MPI/build_release` and `builds/LLVM/build_release`, with external libraries and `config_cmake.sh` in `builds/GCC/MPI/` and `builds/LLVM/` respectively.
- **CI/ephemeral environments**: a single flat tree — `<config_dir>` is `builds/` and `<build_dir>` is `builds/build_release`, installing to `builds/install_release`.

To find candidate build trees, look for `config_cmake.sh`:

```bash
find builds -maxdepth 4 -name config_cmake.sh
```

If more than one is found, **ask the user which to use** rather than guessing — they
may differ in compiler, MPI configuration, or build type, and some may be stale.

### Stage 1: Build external libraries

Run from the directory that will become `<config_dir>`, with the user's Python
environment active (see "Python environment"):

```bash
mkdir -p builds/GCC/MPI && cd builds/GCC/MPI   # one directory per configuration
python <repo_root>/src/utilities/build_libs.py \
  --nthread=4 --build_umfpack=1 --build_superlu=1 \
  --build_arpack=1 --oft_build_tests=1 --build_mpich=1
```

This downloads and compiles dependencies (OpenBLAS, HDF5, METIS, etc.) and generates `config_cmake.sh` — a shell script containing the full CMake invocation with all paths.

### Stage 2: Configure, build, and install OFT

```bash
# Still in <config_dir>, beside config_cmake.sh
bash config_cmake.sh        # runs cmake, creates <config_dir>/build_release/
cd build_release
make                         # compile
make install                 # install to <config_dir>/install_release/
```

### Key CMake options (set via build_libs.py flags)

- `OFT_BUILD_TESTS` — build test executables (enable with `--oft_build_tests=1`)
- `OFT_BUILD_PYTHON` — build Python wrappers (default ON)
- `OFT_USE_MPI` — enable MPI parallelism (set by `--build_mpich=1` or `--build_openmpi=1`)
- `OFT_BUILD_DOCS` — build Doxygen documentation
- `OFT_BUILD_EXAMPLES` — build example programs

### Python environment

Building and testing OFT needs a Python environment providing `pytest`, `numpy`,
`scipy`, `h5py`, `netcdf4`, and `matplotlib`.

- **User environments** (the normal case): there is no project-managed virtual environment. `oft_venv/` and `setup_env.sh` are **not** created by `build_libs.py`, are gitignored, and will not exist in a normal checkout. If the user has not already specified which interpreter or environment to use, **ask them** — do not guess, and do not fall back to the system `python3`, which frequently lacks `numpy`/`h5py`. Carry the answer forward for the rest of the session.
- **CI/ephemeral environments only**: the workflows create a venv at `oft_venv/` and write `source oft_venv/bin/activate` into `setup_env.sh` at the workspace root. `source setup_env.sh` is valid there, and only there.

To use an OFT install from an arbitrary Python environment, either set
`OFT_ROOTPATH=<config_dir>/install_release` (the Python examples append
`$OFT_ROOTPATH/python` to `sys.path`) or put `<config_dir>/install_release/python`
on `PYTHONPATH`.

### Ephemeral agent environments

The `copilot-setup-steps.yml` workflow mirrors the Ubuntu 24.04 GCC 14 + OpenMP CI configuration and pre-builds external libraries in a cached `builds/` directory. After this workflow runs, the agent environment has:
  - Compilers: `gcc-14`, `g++-14`, `gfortran-14`
  - Python venv at `oft_venv/` with: `pytest`, `numpy`, `scipy`, `h5py`, `netcdf4`, `matplotlib`
  - Pre-built external libraries in `builds/`
  - OFT compiled and installed in `builds/install_release/`

This is the flat-layout case: `<config_dir>` is `builds/` and `<build_dir>` is `builds/build_release`.

### Building against system libraries (ephemeral environments, Ubuntu 24.04)

When no pre-built `builds/` directory is available (eg. a fresh cloud sandbox, or any environment where `build_libs.py` cannot reach the network to download dependencies), building against distribution packages is much faster than running `build_libs.py`. This skips Stage 1, writing the CMake invocation into `config_cmake.sh` directly so the tree still matches "Locating the build tree". On Ubuntu 24.04 the following provides MPI (OpenMPI), ARPACK, and SuperLU in addition to the required LibXML2, HDF5, BLAS/LAPACK, and METIS. It uses the flat layout and creates the `oft_venv/`/`setup_env.sh` pair described in "Python environment", so it should not be used in a user checkout without asking.

```bash
apt-get install -y --no-install-recommends gfortran libopenmpi-dev openmpi-bin \
  libarpack2-dev libparpack2-dev libsuperlu-dev libhdf5-dev libxml2-dev \
  libopenblas-dev liblapack-dev libmetis-dev

# From the repository root: Python venv for tests
python3 -m venv oft_venv
echo "source $PWD/oft_venv/bin/activate" > setup_env.sh
source setup_env.sh
python -m pip install pytest numpy scipy h5py netcdf4 matplotlib

# Write config_cmake.sh into <config_dir>, mirroring what build_libs.py generates
mkdir -p builds && cd builds
cat > config_cmake.sh <<'EOF'
ROOT_PATH=$(pwd)
BUILD_DIR=$ROOT_PATH/build_release
INSTALL_DIR=$ROOT_PATH/install_release

# Activate the venv so the configured pytest is the one in it, however the script is run
source $ROOT_PATH/../setup_env.sh

rm -rf $BUILD_DIR
mkdir $BUILD_DIR && cd $BUILD_DIR

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DOFT_USE_MPI:BOOL=TRUE -DOFT_USE_OpenMP:BOOL=TRUE -DOFT_BUILD_TESTS:BOOL=TRUE \
  -DBLA_VENDOR=OpenBLAS \
  -DOFT_ARPACK_ROOT=/usr -DOFT_SUPERLU_ROOT=/usr \
  -DOFT_PYTEST_EXE=$(which pytest) \
  $ROOT_PATH/../src
EOF

bash config_cmake.sh
cd build_release && make -j$(nproc)
```

Notes for this configuration:
- Optional libraries are only searched for when the corresponding `OFT_<LIB>_ROOT` is set, so set it to `/usr` for system packages. SuperLU (`libsuperlu-dev`) is the baseline and is included above; SuperLU-DIST and UMFPACK are not. Install `libsuperlu-dist-dev` or `libsuitesparse-dev` and add `-DOFT_SUPERLU_DIST_ROOT=/usr` or `-DOFT_UMFPACK_ROOT=/usr` only if a specific test requires them.
- Use the default `gcc`/`g++`/`gfortran` (GCC 13) so the compilers match the ones Ubuntu used to build the OpenMPI Fortran modules.
- Ubuntu ships HDF5 1.10, which lacks the Fortran file-space-strategy API (added in 1.14.3). OFT builds fine but persistent free-space tracking in HDF5 output files is disabled, and a warning is printed when it is requested.
- OpenMPI refuses to run as root by default. In a root sandbox, `export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` before running tests.
- `config_cmake.sh` sources `setup_env.sh` itself, so re-running it later resolves `pytest` from the venv rather than whatever happens to be on `PATH`.

### Python packaging information

Python packaging is handled by a dedicated script template at `src/python/make_package.sh.in` and associated `src/python/pyproject.toml.in`, which are configured/installed during the CMAKE build process if the `OFT_PACKAGE_PYTHON` option is set. Dependencies, scripts, and other information are available in this file a supersede the lint-only file `src/python/pyproject.toml`. 

## Testing

Tests use **pytest** and are organized under `src/tests/` in subdirectories matching the source layout: `base/`, `grid/`, `lin_alg/`, `fem/`, `physics/`.

### Test structure

Each test subdirectory has:
- `test_*.py` — pytest test files
- `test_*.F90` — corresponding Fortran test programs (compiled during build)
- Various data files (`.h5`, `.g`, `.inp`, etc.)

### Running tests

With the user's Python environment active (see "Python environment"), from `<build_dir>`:

```bash
make test                    # runs: pytest -m "not slow" base grid lin_alg fem physics
make test_full               # runs all tests including slow ones
make test_examples           # runs example notebook tests
```

Or run individual test files. `run_test.sh` is generated into `<build_dir>/tests/`,
so invoke it from there:

```bash
cd <build_dir>/tests
./run_test.sh physics/test_TokaMaker.py -k "test_name"
```

To list tests for a specific file or folder:

```bash
cd <build_dir>/tests
./run_test.sh physics/test_TokaMaker.py --collect-only
```

To run a specific test `test_name`, optionally with arguments (`test_args`), within a specific file `test_file.py` (current working director must contain `test_file.py`):

```bash
python -c "from test_file import test_name; test_name(*test_args)"
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

Python code is linted with **ruff**. Configuration is in `src/python/pyproject.toml` (superseded by `src/python/pyproject.toml.in` for python packaging):

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

## CI Workflows (utilizing GitHub actions)

| Workflow | File | Trigger | Purpose |
|---|---|---|---|
| **CI Build** | `ci_build.yaml` | push to main, PRs | Full matrix build (GCC, Intel, macOS) × (OpenMP, MPICH, OpenMPI) |
| **Lint** | `lint.yaml` | push to main, PRs | ruff check + Fortran stack check |
| **Coverage** | `cov_build.yaml` | push to main | Build with `--coverage`, upload to Codecov |
| **CD Build and Deployment** | `cd_combined.yaml` | push to main, pypi_release, PRs | Release and nightly package/PyPI builds (deploy from `pypi_release`) |
| **Container CD Build and Deployment** | `container_cd.yaml` | push to main, container_release, PRs | Container release builds (deploy from `container_release`) |
| **GitHub Pages Build** | `website.yml` | push to main, gh-pages, PRs | Documentation/website build (deploy from `gh-pages`) |
| **Copilot Agent Setup** | `copilot-setup-steps.yml` | manual/PR | Copilot agent environment setup |
| **Copilot Review Setup** | `copilot-code-review.yml` | manual/PR | Copilot code review agent environment setup |

### Representative build environment

A representative build environment that can be used for creating a testing sandbox as used in `copilot-setup-steps.yml` is as follows:
- **OS**: Ubuntu 24.04
- **Compilers**: `gcc-14` / `g++-14` / `gfortran-14`
- **Parallel**: OpenMP + MPICH (MPI enabled)
- **Libraries**: OpenBLAS, UMFPACK, SuperLU, ARPACK, HDF5, METIS

## Coding Conventions

### Fortran

- Free-form Fortran 90+ (`.F90` extension, preprocessed)
- SPDX license header: `! SPDX-License-Identifier: LGPL-3.0-only`
- Doxygen-style comments with `!>` and `!!` markers
- Module names typically prefixed with `oft_` (e.g., `oft_gs`, `oft_la_base`)
- Line length: unlimited (`-ffree-line-length-none`), but keep reasonable
- New functions/subroutines must have Doxygen documentation — this applies to code Doxygen actually processes; see "Test and example code" below
- All fortran keywords, statements, and built-in functions should be in UPPERCASE 
- New functions should use `RESULT` style definition instead of leading variable style
- `INTEGER` and `REAL` definitions should explicitly provide the `KIND`/precision using one of `i4`/`i8` or `r4`/`r8` (defined in `src/base/local.F90`) respectively
- Two space indentation should be used, with the exception of subroutine/function bodies

### Python

- SPDX license header comment block at top of each file
- New functions/subroutines must have Doxygen documentation using docstrings with `@`-style commands (eg. `@param`, `@result`) — this applies to code Doxygen actually processes; see "Test and example code" below
- Target Python 3.7+ (for OMFIT compatibility)
- Use `ctypes` for Fortran interop; wrapper patterns in `_interface.py` files
- All Python wrappers live under `src/python/OpenFUSIONToolkit/`
- Balance readability and efficiency; For example, avoid overly long or nested list comprehensions and inline if statements
- Four space indentation should be used

### Test and example code

Doxygen's `INPUT` list in `src/docs/Doxyfile.in` covers `src/docs`, `base`, `grid`,
`lin_alg`, `fem`, `physics`, `bin`, and `python`. Neither `src/tests/` nor
`src/examples/` is scanned directly, but examples still reach the documentation
indirectly: the `build_ex_docs` target runs `src/utilities/generate_doc.py`, which
reads `examples/*.F90`, `examples/*/*/*.F90`, and `examples/*/*/*.ipynb` and writes
`docs/generated/doc_*.md`, which Doxygen then picks up through its `docs` entry.

- **Example sources are documentation.** Fortran examples under `src/examples/` carry their own `!!` comment format, including Doxygen page anchors and linking commands (eg. `!!API Example 1 ... {#doc_api_ex1}`, `\section`, `\ref`), which `generate_doc.py` parses into published pages. The same applies to the example notebooks. Do not strip or "correct" this markup, and follow the surrounding file when extending it.
- **Everything else under `src/tests/` and `src/examples/`** is not parsed into the API documentation, so its functions do not need Doxygen-style documentation — no `!>`/`!!` function blocks in Fortran, no `@`-style commands (`@param`, `@result`, `@brief`) in Python docstrings. Plain comments and plain docstrings are fine.
- Doxygen **linking** commands (`\section`, `\ref`, …) remain valid anywhere and may be used where they help, even in files whose functions need no documentation blocks.
- Function-level documentation blocks in test code are **optional**.
- Brief documentation is still helpful and welcome: a short line saying what a test covers or what an example demonstrates aids the next reader.

When reviewing changes, do not flag missing documentation on test functions, and do
not request Doxygen-style function documentation for files outside the paths
Doxygen processes.

### Pull Requests

- PRs for a specific tool should be titled: `ToolName: description` (e.g., `TokaMaker: Fix boundary condition`)
- Keep changes focused on a single feature; secondary changes should be minimal
- Comment on whether APIs or input files change
- Run regression tests before submitting

## Common Pitfalls

1. **Single OFT_env instance**: The Python `OFT_env` class enforces a singleton. Tests use `multiprocessing.Process` to work around this. Never create two `OFT_env` instances in the same process.

2. **Tests run from build tree**: Tests must be run from `<build_dir>/tests/` (or via `make test` from `<build_dir>`), not from the source tree, because compiled Fortran test executables are in the build tree. See "Locating the build tree" — `<build_dir>` is only `builds/build_release` in CI/ephemeral environments.

3. **ext_libs/ is bundled third-party code**: Do not modify files in `src/ext_libs/`. These are upstream sources (triangle, minpack, bvls, dlsode).

4. **Fortran/Python interop**: The Python package calls compiled Fortran via `ctypes`. The Fortran-side wrappers are in `src/python/wrappers/` and use `ISO_C_BINDING`. Changes to Fortran function signatures require corresponding updates to both `wrappers/*_f.F90` and `python/OpenFUSIONToolkit/*/_interface.py`.

5. **CMake template files**: Some files use `@VARIABLE@` CMake substitution (e.g., `__init__.py`, `pyproject.toml.in`, `run_test.sh.in`). Edit the `.in` template, not the generated file.
