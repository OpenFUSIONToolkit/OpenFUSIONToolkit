# Auto-Generated on Mon Jan 19 22:31:22 2026
# using library build at /home/andrew-maris/oftExternalLibs
# on machine ADMLenovo
# settings: --nthread=2 --build_umfpack=1 --build_arpack=1 --oblas_threads

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
  -DOFT_BUILD_TESTS:BOOL=FALSE \
  -DOFT_BUILD_EXAMPLES:BOOL=FALSE \
  -DOFT_BUILD_PYTHON:BOOL=TRUE \
  -DOFT_BUILD_DOCS:BOOL=FALSE \
  -DOFT_USE_OpenMP:BOOL=TRUE \
  -DOFT_PACKAGE_BUILD:BOOL=FALSE \
  -DOFT_PACKAGE_NIGHTLY:BOOL=TRUE \
  -DOFT_COVERAGE:BOOL=FALSE \
  -DOFT_DEBUG_STACK:BOOL=FALSE \
  -DOFT_PROFILING:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH=gcc \
  -DCMAKE_CXX_COMPILER:FILEPATH=g++ \
  -DCMAKE_Fortran_COMPILER:FILEPATH=gfortran \
  -DCMAKE_Fortran_FLAGS:STRING="-fallow-argument-mismatch" \
  -DOFT_USE_MPI:BOOL=FALSE \
  -DOFT_METIS_ROOT:PATH=/home/andrew-maris/oftExternalLibs/metis-5_1_0 \
  -DHDF5_ROOT:PATH=/home/andrew-maris/oftExternalLibs/hdf5-1_14_6 \
  -DBLAS_ROOT:PATH=/home/andrew-maris/oftExternalLibs/OpenBLAS-0_3_30 \
  -DLAPACK_ROOT:PATH=/home/andrew-maris/oftExternalLibs/OpenBLAS-0_3_30 \
  -DBLA_VENDOR:STRING=OpenBLAS \
  -DOFT_ARPACK_ROOT:PATH=/home/andrew-maris/oftExternalLibs/arpack-ng-3_9_1 \
  -DOFT_FoX_ROOT:PATH=/home/andrew-maris/oftExternalLibs/fox-4_1_2 \
  -DOFT_UMFPACK_ROOT:PATH=/home/andrew-maris/oftExternalLibs/UMFPACK-6_3_5 \
  /home/andrew-maris/repos/OpenFUSIONToolkit/src
