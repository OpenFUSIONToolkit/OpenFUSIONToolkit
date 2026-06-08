build_base=$(pwd)
rm -f $build_base/build_tmp.stat
export CC=gcc-15
export FC=gfortran-15
make clean; echo $? >> $build_base/build_tmp.stat
make NO_CBLAS=1 NO_LAPACKE=1 NO_SHARED=1 USE_THREAD=0 USE_LOCKING=1 FCOMMON_OPT="-frecursive -fPIC" NO_AVX=1 NO_AVX2=1 MAKE_NB_JOBS=2; echo $? >> $build_base/build_tmp.stat
make NO_CBLAS=1 NO_LAPACKE=1 NO_SHARED=1 USE_THREAD=0 USE_LOCKING=1 FCOMMON_OPT="-frecursive -fPIC" NO_AVX=1 NO_AVX2=1 NO_PARALLEL_MAKE=1 PREFIX=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/OpenBLAS-0_3_30 install; echo $? >> $build_base/build_tmp.stat