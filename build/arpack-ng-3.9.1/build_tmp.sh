build_base=$(pwd)
rm -f $build_base/build_tmp.stat
rm -rf build; echo $? >> $build_base/build_tmp.stat
mkdir build; echo $? >> $build_base/build_tmp.stat
cd build; echo $? >> $build_base/build_tmp.stat
export CC=gcc-15
export F77=gfortran-15
export FC=gfortran-15
export FFLAGS="-fallow-argument-mismatch"
cmake -DCMAKE_INSTALL_PREFIX:PATH=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/arpack-ng-3_9_1 -DCMAKE_INSTALL_LIBDIR=lib -DEXAMPLES=OFF -DTESTS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF -DBLAS_LIBRARIES:PATH=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/OpenBLAS-0_3_30/lib/libopenblas.a -DLAPACK_LIBRARIES:PATH=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/OpenBLAS-0_3_30/lib/libopenblas.a -DMPI=OFF ..; echo $? >> $build_base/build_tmp.stat
make -j2; echo $? >> $build_base/build_tmp.stat
make install; echo $? >> $build_base/build_tmp.stat