---
name: oft-build
description: Configure, compile, rebuild, and install the Open FUSION Toolkit (OFT) from source. Use when OFT needs to be (re)built after changes to Fortran, C/C++, or Python wrapper code, when locating or setting up a build tree, when building external libraries with build_libs.py, or when building against Ubuntu system packages in a fresh sandbox.
---

# Building and installing OFT

OFT uses a two-stage build: `src/utilities/build_libs.py` compiles third-party
dependencies and writes `config_cmake.sh`, which runs CMake to create the OFT build
tree. Throughout, `<config_dir>` is the directory holding `config_cmake.sh` (and the
compiled external libraries), and `<build_dir>` is `<config_dir>/build_release`.

## 1. Locate the build tree

The build tree's location is chosen by the user and is **not** fixed by the project.
**Do not assume `builds/build_release`.**

- **User environments** (the normal case): several trees usually sit side by side, one per toolchain/parallelism combination, eg. `builds/GCC/MPI/` and `builds/LLVM/`, each with its own `config_cmake.sh`, external libraries, `build_release/`, and `install_release/`.
- **CI/ephemeral environments**: a single flat tree — `<config_dir>` is `builds/`, `<build_dir>` is `builds/build_release`, installing to `builds/install_release`.

```bash
find builds -maxdepth 4 -name config_cmake.sh
```

If more than one is found, **ask the user which to use** — they may differ in compiler,
MPI configuration, or build type, and some may be stale.

### Check which source tree it builds

The last line of `config_cmake.sh` is the **absolute** path of the `src/` directory it
compiles (also `CMAKE_HOME_DIRECTORY` in `<build_dir>/CMakeCache.txt`). In a git
worktree, or any second checkout, an existing build tree usually points at the *main*
checkout, so building it will not compile your changes. To build the current checkout
without disturbing the user's tree, make a sibling configuration that reuses the same
external libraries:

```bash
mkdir -p <new_config_dir>
sed '$ s|[^ ].*|<this_checkout>/src|' <config_dir>/config_cmake.sh > <new_config_dir>/config_cmake.sh
cd <new_config_dir> && bash config_cmake.sh
```

`config_cmake.sh` uses `$(pwd)` for the build and install paths, and absolute paths for
the external libraries, so the copy is self-contained. Put `<new_config_dir>` somewhere
gitignored (eg. under `builds/`) and tell the user where it is.

## 2. Python environment

Building (and testing) needs a Python environment with `pytest`, `numpy`, `scipy`,
`h5py`, `netcdf4`, and `matplotlib`. CMake records the interpreter and `pytest` it finds
at configure time, so the right environment must be active when `config_cmake.sh` runs.

- **User environments**: there is no project-managed virtual environment. `oft_venv/` and `setup_env.sh` are **not** created by `build_libs.py`, are gitignored, and will not exist in a normal checkout. If the user has not already said which interpreter or environment to use, **ask them** — do not guess, and do not fall back to the system `python3`, which frequently lacks `numpy`/`h5py`. Carry the answer forward for the rest of the session.
- **CI/ephemeral environments only**: the workflows create a venv at `oft_venv/` and write `source oft_venv/bin/activate` into `setup_env.sh` at the workspace root. `source setup_env.sh` is valid there, and only there.

## 3. Rebuild after source changes (the common case)

```bash
cd <build_dir>
make -j<N>          # compile; re-runs CMake automatically if CMake inputs changed
make install        # only if something uses <config_dir>/install_release (see below)
```

- **Do not re-run `config_cmake.sh` just to rebuild.** It starts with `rm -rf build_release`, discarding the whole build tree and forcing a full recompile. Only re-run it to change CMake options or after the source path or external libraries change.
- Regression tests run from the build tree (`<build_dir>/tests`, with `<build_dir>/python` for the Python package), so `make` alone is enough before running tests — see the `oft-test` skill.
- `make install` is needed when using OFT from an arbitrary Python environment, via `OFT_ROOTPATH=<config_dir>/install_release` (the Python examples append `$OFT_ROOTPATH/python` to `sys.path`) or `PYTHONPATH=<config_dir>/install_release/python`.
- Changing a Fortran function called from Python requires matching edits to `src/python/wrappers/*_f.F90` and `src/python/OpenFUSIONToolkit/*/_interface.py`; both are rebuilt/copied by `make`.
- Files using `@VARIABLE@` substitution (`__init__.py`, `pyproject.toml.in`, `run_test.sh.in`, …) are generated at configure time — edit the template in `src/`, then `make`.

## 4. Set up a new configuration from scratch

### Stage 1: external libraries

Run from the directory that will become `<config_dir>`, with the Python environment
active. This downloads and compiles OpenBLAS, HDF5, METIS, etc., which takes a while —
run it in the background.

```bash
mkdir -p builds/GCC/MPI && cd builds/GCC/MPI   # one directory per configuration
python <repo_root>/src/utilities/build_libs.py \
  --nthread=4 --build_umfpack=1 --build_superlu=1 \
  --build_arpack=1 --oft_build_tests=1 --build_mpich=1
```

`python src/utilities/build_libs.py --help` lists every option. Commonly used ones:

| Flag | Effect |
|---|---|
| `--oft_build_tests=1` | Build regression test executables (`OFT_BUILD_TESTS`; default off) |
| `--oft_build_python={0,1}` | Build Python wrappers (`OFT_BUILD_PYTHON`; default on) |
| `--build_mpich=1` / `--build_openmpi=1` | Build MPI and enable `OFT_USE_MPI` |
| `--oft_use_openmp={0,1}` | OpenMP support (default on) |
| `--oft_build_debug=1` | Debug build (runtime checks on; tests get 4× timeouts) |
| `--oft_build_examples=1` | Build examples (`OFT_BUILD_EXAMPLES`) |
| `--oft_build_docs=1` | Build Doxygen documentation (`OFT_BUILD_DOCS`) |
| `--build_arpack=1`, `--build_superlu=1`, `--build_superlu_dist=1`, `--build_umfpack=1`, `--build_petsc=1` | Optional solver libraries |

### Stage 2: configure, build, install

```bash
# In <config_dir>, beside config_cmake.sh
bash config_cmake.sh        # runs cmake, creates <config_dir>/build_release/
cd build_release
make -j<N>
make install                 # installs to <config_dir>/install_release/
```

## 5. Ephemeral agent environments

`.github/workflows/copilot-setup-steps.yml` mirrors the Ubuntu 24.04 GCC 14 + OpenMP CI
configuration. After it runs, the environment has `gcc-14`/`g++-14`/`gfortran-14`, a
venv at `oft_venv/` (activated by `source setup_env.sh`), external libraries pre-built
in `builds/`, and OFT installed in `builds/install_release/` — the flat layout.

### Building against system libraries (Ubuntu 24.04)

When no pre-built `builds/` directory exists (eg. a fresh cloud sandbox, or one where
`build_libs.py` cannot download dependencies), building against distribution packages is
much faster. This skips Stage 1 and writes `config_cmake.sh` by hand so the tree still
matches the layout above. It uses the flat layout and creates the `oft_venv/` +
`setup_env.sh` pair, so **do not use it in a user checkout without asking**.

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

## 6. Python packaging

Packaging builds use `src/python/make_package.sh.in` and `src/python/pyproject.toml.in`,
configured during the CMake build for packaging builds (`OFT_PACKAGE_BUILD` and
`OFT_PACKAGE_PYTHON`, ie. `build_libs.py --oft_package`, with `--oft_package_python`
defaulting to on). Dependencies and scripts listed there supersede the lint-only
`src/python/pyproject.toml`.
