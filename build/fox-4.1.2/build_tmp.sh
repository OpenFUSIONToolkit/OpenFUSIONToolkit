build_base=$(pwd)
rm -f $build_base/build_tmp.stat
make distclean; echo $? >> $build_base/build_tmp.stat
export CC=gcc-15
export FC=gfortran-15
export CFLAGS=-fPIC
export FCFLAGS=-fPIC
export GFORTRAN_UNBUFFERED_ALL=Y
./configure --prefix=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/fox-4_1_2 --enable-dom ; echo $? >> $build_base/build_tmp.stat
make -j2; echo $? >> $build_base/build_tmp.stat
make install; echo $? >> $build_base/build_tmp.stat