# Auto-Generated on Thu Mar 12 11:52:40 2026
# using library build at /home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs
# on machine runnervm46oaq
# settings: --oblas_dynamic_arch --build_umfpack=1 --build_superlu=1 --no_dl_progress --nthread=3 --build_arpack=1 --oft_build_tests=1

# Setup build and install paths
ROOT_PATH=$(pwd)
BUILD_DIR=$ROOT_PATH/build_release
INSTALL_DIR=$ROOT_PATH/install_release

# Create fresh build directory
rm -rf $BUILD_DIR
mkdir $BUILD_DIR && cd $BUILD_DIR

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DOFT_BUILD_TESTS:BOOL=TRUE \
  -DOFT_PY_KERNEL:STRING=python3 \
  -DOFT_BUILD_EXAMPLES:BOOL=FALSE \
  -DOFT_BUILD_PYTHON:BOOL=TRUE \
  -DOFT_BUILD_DOCS:BOOL=FALSE \
  -DOFT_USE_OpenMP:BOOL=TRUE \
  -DOFT_PACKAGE_BUILD:BOOL=FALSE \
  -DOFT_PACKAGE_NIGHTLY:BOOL=TRUE \
  -DOFT_COVERAGE:BOOL=FALSE \
  -DOFT_DEBUG_STACK:BOOL=FALSE \
  -DOFT_PROFILING:BOOL=FALSE \
  -DOFT_TOKAMAKER_LEGACY:BOOL=FALSE \
  -DOFT_THINCURR_LEGACY:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH=gcc \
  -DCMAKE_CXX_COMPILER:FILEPATH=g++ \
  -DCMAKE_Fortran_COMPILER:FILEPATH=gfortran \
  -DCMAKE_Fortran_FLAGS:STRING="-fallow-argument-mismatch" \
  -DOFT_USE_MPI:BOOL=FALSE \
  -DOFT_METIS_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/metis-5_1_0 \
  -DHDF5_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6 \
  -DBLAS_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/OpenBLAS-0_3_30 \
  -DLAPACK_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/OpenBLAS-0_3_30 \
  -DBLA_VENDOR:STRING=OpenBLAS \
  -DOFT_ARPACK_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/arpack-ng-3_9_1 \
  -DLIBXML2_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2 \
  -DOFT_FoX_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/fox-4_1_2 \
  -DOFT_SUPERLU_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/superlu-7_0_0 \
  -DOFT_UMFPACK_ROOT:PATH=/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/UMFPACK-6_3_5 \
  /home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src
