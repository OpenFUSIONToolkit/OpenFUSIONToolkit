build_base=$(pwd)
rm -f $build_base/build_tmp.stat
rm -rf build; echo $? >> $build_base/build_tmp.stat
mkdir build; echo $? >> $build_base/build_tmp.stat
cd build; echo $? >> $build_base/build_tmp.stat
export CC=gcc-15
export FC=gfortran-15
../configure --prefix=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/hdf5-1_14_6 --enable-fortran --enable-hl=no --enable-tests=no --enable-shared=yes --enable-static=no --with-pic --enable-hl=no; echo $? >> $build_base/build_tmp.stat
make -j2; echo $? >> $build_base/build_tmp.stat
make install; echo $? >> $build_base/build_tmp.stat