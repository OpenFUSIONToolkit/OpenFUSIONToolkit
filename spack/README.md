# Spack Build Instructions for OFT

## Platform notes
On all platforms it is recommended to use the newest release of the [Spack Packages repo](https://github.com/spack/spack-packages), which is presently
`v2026.02`. This can be done by adding or modifying the user's `~/.spack/repos.yaml` file (note plural in naming) to include the following:
```shell
repos:
  builtin:
    branch: releases/v2026.02
```
After updating the file you then need to run `spack repo update`. For more information see [the Spack documentation](https://spack.readthedocs.io/en/latest/repositories.html#configuring-repositories-with-repos-yaml)

### Linux

#### Find external packages/utilities
By default, Spack takes a "build everything" approach. However, it is often useful to use existing system
installations of common tools, which are only used and the build process and are not linked to OFT libraries
or binaries. On Linux the following packages may be helpful to use from the base system.

Find system python (if desired)
```shell
spack external find python
```

Find system cmake (if desired)
```shell
spack external find cmake
```

Find system utilities
```shell
spack external find gmake
spack external find pkgconf
spack external find m4
spack external find perl
spack external find curl
spack external find gettext
```

### macOS

#### Find external packages/utilities
By default, Spack takes a "build everything" approach. However, it is often useful to use existing system
installations of common tools, which are only used and the build process and are not linked to OFT libraries
or binaries. On macOS the following packages may be helpful to use from the base system.

Find system python (if desired)
```shell
spack external find python
```

Find system cmake (if desired)
```shell
spack external find cmake
```

Find system utilities (useful to reduce build time)
```shell
spack external find gmake
spack external find perl
spack external find gettext   # Installed with HomeBrew
spack external find hwloc     # Installed with HomeBrew, if using MPI
```
**Note:** On macOS `m4`, `automake` and `autoconf` are often not the correct version for tools, causing build
issues. Therefore it is not recommended to use system/HomeBrew version of these packages.

## Development Environment (recommended)
The recommended way to build OFT with Spack is to create a "development environment" that enables
Spack to detect changes to the codebase and rebuild as needed. In this way the `spack concretize`
and `spack install` steps now replace (roughly) the `build_libs.py` and `make`/`make install` steps
in the normal build process.

### Environment setup

#### Load Spack (if not already loaded)
```shell
. /path/to/spack/share/spack/setup-env.sh
```

#### Find compilers (if needed)
```shell
spack compiler find
```

#### Create development environment
Replace `/path/to/oft` with the path to the root of the OFT git repo (folder should contain `src`)
```shell
spack env create -d /path/to/oft
```

#### Activate environment
```shell
spack env activate /path/to/oft
```

#### Add local repo to environment
```shell
spack repo add /path/to/oft/spack/repo
```

#### Mark OFT repo for development
```shell
spack develop --path /path/to/oft openfusiontoolkit
```

#### Add OFT to install environment
Now add OFT to the environment, being sure to modify the spec options (`+spec`/`~spec`) and
compilers (`gcc@15.2.0`) accordingly. Please note, that at this point OFT requires the same
compiler vendor to be used for both C/C++ and Fortran compilers.

**Recommended/Example OpenMP build**
```shell
# Example for SuperLU
spack add openfusiontoolkit+arpack+superlu %c,cxx,fortran=gcc@15.2.0 ^openblas~dynamic_dispatch
```

**Recommended/Example MPI build**
```shell
# Example for built MPI and SuperLU
spack add openfusiontoolkit+arpack+mpi+superlu %c,cxx,fortran=gcc@15.2.0 ^openblas~dynamic_dispatch ^mpich
```

**Recommended/Example PETSc build**
```shell
# Example for built MPI and PETSc (w/ SuperLU-DIST and MUMPS)
spack add openfusiontoolkit+arpack+mpi+petsc %c,cxx,fortran=gcc@15.2.0 ^openblas~dynamic_dispatch ^mpich ^petsc~examples+superlu-dist+mumps
```

### Environment usage

#### Load Spack (if not already loaded)
```shell
. /path/to/spack/share/spack/setup-env.sh
```

#### Activate environment (if not already activated)
```shell
spack env activate /path/to/oft
```

#### Concretize
```shell
spack concretize -f
```

#### Install
```shell
spack install
```

#### Set PYTHON environment variables
```shell
spack load openfusiontoolkit
```

or manually

```shell
OFT_ROOTPATH=$(spack location -i openfusiontoolkit)
PYTHONPATH=$(spack location -i openfusiontoolkit)/python:$PYTHONPATH
```

### To run tests
Make sure to add `+tests` to the spec list when adding `openfusiontoolkit` and that you have an
active python environment with `pytest` available.
```shell
spack install --keep-stage openfusiontoolkit
cd $(spack location -b openfusiontoolkit)
make test
```

Or to run tests silently using Spack's built-in testing
```shell
spack install --test=root openfusiontoolkit
# TODO: Add capability to save test result file
```