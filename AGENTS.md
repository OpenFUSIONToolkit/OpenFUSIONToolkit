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
├── lin_alg/        # Linear algebra (native solvers, PETSc, ARPACK, SuperLU, UMFPACK)
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

## Skills

Detailed, task-specific procedures live in agent skills under `.agents/skills/`
(symlinked to `.claude/skills/`). Read the relevant `SKILL.md` before doing the task:

| Skill | Use it to |
|---|---|
| [`oft-build`](.agents/skills/oft-build/SKILL.md) | Locate or set up a build tree, build external libraries, configure/compile/install OFT, build against system libraries in a sandbox, or change Python packaging (`make_package.sh.in`/`pyproject.toml.in`, which supersede the lint-only `src/python/pyproject.toml`) |
| [`oft-test`](.agents/skills/oft-test/SKILL.md) | Run, select, debug, or write regression tests |

## Build System

OFT uses a two-stage build: `src/utilities/build_libs.py` builds external libraries and
generates `config_cmake.sh`, which runs CMake (≥ 3.27) to create the build tree; then
`make` and `make install`. Full instructions are in the `oft-build` skill.

### Locating the build tree

The build tree is created by the user and its location is **not** fixed by the
project. These instructions use `<build_dir>` for the configured build directory and
`<config_dir>` for its parent (holding `config_cmake.sh`, the compiled external
libraries, and `install_release/`; `build_debug`/`install_debug` for Debug configurations). **Do not assume `builds/build_release`.**

- **User environments** (the normal case): several build trees side by side, one per toolchain/parallelism combination, eg. `builds/GCC/MPI/build_release` and `builds/LLVM/build_release`.
- **CI/ephemeral environments**: a single flat tree — `<config_dir>` is `builds/`, `<build_dir>` is `builds/build_release`.

Find candidates with `find builds -maxdepth 4 -name config_cmake.sh`. If more than one
is found, **ask the user which to use** rather than guessing. If none is found (eg. in a git
worktree, which never has `builds/`), follow the `oft-build` skill rather than starting a
new build unasked.

### Python environment

Building and testing need a Python environment with `pytest`, `numpy`, `scipy`, `h5py`,
`netcdf4`, and `matplotlib`.

- **User environments**: there is no project-managed virtual environment; `oft_venv/` and `setup_env.sh` are gitignored and will not exist in a normal checkout. If either is present anyway, it was made by hand: do not source it without asking, and check that the venv it activates still exists. If the user has not already specified which interpreter or environment to use, **ask them** — do not fall back to the system `python3`. Carry the answer forward for the rest of the session.
- **CI/ephemeral environments only**: the workflows create `oft_venv/` and a `setup_env.sh` that activates it. `source setup_env.sh` is valid there, and only there.

## Testing

Tests use **pytest** and live in `src/tests/` in subdirectories matching the source
layout (`base/`, `grid/`, `lin_alg/`, `fem/`, `physics/`), pairing `test_*.py` files with
the `test_*.F90` programs they drive. They run from `<build_dir>` (`make test`, or
from `<build_dir>/tests`, `./run_test.sh <file> -k <name>` for targeted runs; the file
path is resolved against the current directory), never from the source tree. See the `oft-test` skill for running, markers, debugging, and writing tests.

Always set `OMP_NUM_THREADS` when running tests: `1` if only one or two cores are
available, otherwise `2`, regardless of core count. The only exception is running
individual tests specifically to study OpenMP thread scaling or performance.

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
| **CI Build** | `ci_build.yaml` | push to main, ci_debug; PRs — only when `src/**` (excluding `src/docs/**`) or the workflow file changes | Full matrix build (GCC, Intel, macOS) x (OpenMP, MPICH, OpenMPI) |
| **Lint** | `lint.yaml` | push to main, PRs | ruff check + Fortran stack check |
| **Coverage** | `cov_build.yaml` | push to main, cov_debug; PRs — only when `src/**` (excluding `src/docs/**`) or the codecov/workflow files change | Build with `--coverage`, upload to Codecov |
| **CD Build and Deployment** | `cd_combined.yaml` | push to main, pypi_release; manual | Release and nightly package/PyPI builds (deploy from `pypi_release`) |
| **Container CD Build and Deployment** | `container_cd.yaml` | push to container_release; manual | Container release builds (deploy from `container_release`) |
| **GitHub Pages Build** | `website.yaml` | push to main, gh-pages, PRs | Documentation/website build (deploy from `gh-pages`) |
| **Cleanup build caches** | `prune_cache.yaml` | PR closed; manual (with PR number) | Delete the GitHub Actions caches created by a PR's branch |
| **Copilot Agent Setup** | `copilot-setup-steps.yml` | manual; PRs changing the workflow file or `src/utilities/build_libs.py` | Copilot agent environment setup |
| **Copilot Review Setup** | `copilot-code-review.yml` | manual; PRs changing the workflow file or `src/utilities/build_libs.py` | Copilot code review agent environment setup |

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

6. **Build trees compile an absolute source path**: `config_cmake.sh` ends with the absolute path of the `src/` it builds, so in a git worktree or second checkout an existing build tree usually compiles the *main* checkout, not your changes. See "Check which source tree it builds" in the `oft-build` skill.
