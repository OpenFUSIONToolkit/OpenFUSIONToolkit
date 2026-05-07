# Spack Build Instructions for OFT

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

#### Find external packages/utilities
By default, Spack takes a "build everything" approach. However, it is often useful to use existing system
installations of common tools, which are only used and the build process and are not linked to OFT libraries
or binaries.

Find system python (if desired)
```shell
spack external find python
```

Find system cmake (if desired)
```shell
spack external find cmake
```

Find system utilities (useful on macOS or to reduce build time)
```shell
spack external find gmake
spack external find pkgconf
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

#### Mark OFT repo for development
```shell
spack develop --path /path/to/oft openfusiontoolkit
```

#### Add local repo to environment
```shell
spack repo add /path/to/oft/spack/repo
```

#### Add OFT to install environment
Now add OFT to the environment, being sure to modify the spec options (`+spec`/`~spec`) and
compilers (`gcc@15.2.0`) accordingly. Please note, that at this point OFT requires the same
compiler vendor to be used for both C/C++ and Fortran compilers.
```shell
spack add openfusiontoolkit+superlu %c,cxx=gcc@15.2.0 %fortran=gcc@15.2.0 # Example with GCC 15.2
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