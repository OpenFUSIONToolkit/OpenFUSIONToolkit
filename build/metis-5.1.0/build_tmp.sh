build_base=$(pwd)
rm -f $build_base/build_tmp.stat
GKLIB_PATH=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/build/metis-5.1.0/GKlib; echo $? >> $build_base/build_tmp.stat
rm -rf build; echo $? >> $build_base/build_tmp.stat
mkdir -p build; echo $? >> $build_base/build_tmp.stat
cd build; echo $? >> $build_base/build_tmp.stat
cmake -DCMAKE_INSTALL_PREFIX=/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/metis-5_1_0 -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON -DGKRAND:BOOL=ON -DGKLIB_PATH=$GKLIB_PATH -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .. ; echo $? >> $build_base/build_tmp.stat
make -j2; echo $? >> $build_base/build_tmp.stat
make install; echo $? >> $build_base/build_tmp.stat