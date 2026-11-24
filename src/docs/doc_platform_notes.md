Platform and configuration notes  {#doc_platform_notes}
===========

[TOC]

This page collects notes on building and running the Open FUSION Toolkit with specific compilers,
platforms, and library configurations. General build instructions are maintained on the
[GitHub wiki](https://github.com/openfusiontoolkit/OpenFUSIONToolkit/wiki).

\section doc_platform_notes_intel Intel oneAPI compilers (icx/icpx/ifx)

\subsection doc_platform_notes_intel_stack Stack size limits

`ifx` places large temporary arrays on the stack by default. With the default Linux stack limit
(typically 8 MB) larger runs can fail with a segmentation fault (signal 11), often during setup
before any output indicates a problem. For example, the MUG `Spheromak_heating` example crashes
at the default limit. Remove the limit in the shell (or batch script) that launches OFT:

~~~{.sh}
ulimit -s unlimited
mpiexec -n 4 ./MUG_sph_heat oft.in oft_in.xml
~~~

Notes:
 - The limit must be raised in the environment that starts each MPI process. For multi-node runs,
   set it in the batch script or shell startup files on every node, and check that your MPI launcher
   propagates it.
 - `ulimit -s` only applies to the main thread of each process. When running with multiple OpenMP
   threads, the stack size of the other threads is set separately with `OMP_STACKSIZE`
   (eg. `export OMP_STACKSIZE=512M`).
 - Alternatively, OFT can be compiled with `-heap-arrays` added to the Fortran flags, which moves
   temporaries to the heap at some performance cost.

\subsection doc_platform_notes_intel_lto Link-time optimization

With `ifx` 2025.1, the link-time optimization (LTO) link of the Python wrapper library
(`liboftpy.so`) crashes inside the compiler's linker (`ld.lld`) with "Broken module found".
Disable LTO for Intel builds by passing `--oft_use_lto=0` to `build_libs.py`
(`-DOFT_USE_LTO:BOOL=FALSE` in `config_cmake.sh`).

\subsection doc_platform_notes_intel_env Compiler environment and MKL

When building the external libraries with `build_libs.py` (eg. `--build_mpich=1 --use_mkl`),
load only the compiler and MKL components of oneAPI, rather than the full `setvars.sh`:

~~~{.sh}
source /opt/intel/oneapi/compiler/latest/env/vars.sh
source /opt/intel/oneapi/mkl/latest/env/vars.sh
export CC=icx CXX=icpx FC=ifx
~~~

The full `setvars.sh` also adds Intel MPI and its bundled libfabric to the compiler search paths,
and MPICH picks that libfabric up at configure time, which can cause the MPICH build to fail.
The same environment (in particular `MKLROOT`) must be loaded when running `config_cmake.sh`, so
CMake can find MKL, and when running OFT, so the Intel runtime and MKL libraries can be found.

When MKL is used for BLAS/LAPACK, the PARDISO direct solver (`mkl` package) is enabled automatically.

\section doc_platform_notes_mpich MPICH and system libfabric

MPICH's configure step adds the static dependencies of an external libfabric
(`pkg-config --static --libs libfabric`) to its link line. Some Linux distributions list
libraries in `Libs.private` of `libfabric.pc` (eg. `-lpsm2 -lpsm_infinipath`) without installing
the corresponding development files, so the MPICH build in `build_libs.py` fails with errors like
`cannot find -lpsm2`. Either install the missing development packages, or put a copy of
`libfabric.pc` with the `Libs.private` line removed in a directory listed first in
`PKG_CONFIG_PATH` while building the libraries.

\section doc_platform_notes_umfpack UMFPACK and large factorizations

By default OFT uses the 32-bit integer UMFPACK interface (`umfpack_di_*`), which limits the size of
UMFPACK's internal workspace to roughly 16 GB per factorization. Larger factorizations fail with
"UMFPACK: Out of Memory in factorization", even when sufficient memory is available. To use the
64-bit integer interface (`umfpack_dl_*`), pass `--umfpack_long` together with `--build_umfpack=1`
to `build_libs.py` (`-DOFT_UMFPACK_LONG:BOOL=TRUE` in `config_cmake.sh`). The option applies only to
the UMFPACK library built by `build_libs.py` for OFT's native solvers; it is not available in PETSc
builds (`--build_petsc=1`), and `build_libs.py` stops with an error if either condition is not met.
