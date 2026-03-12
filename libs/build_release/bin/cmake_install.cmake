# Install script for directory: /home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/bin

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/install_release")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check"
         RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/bin/oft_mesh_check")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check"
         OLD_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib:"
         NEW_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace"
         RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/bin/oft_trace")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace"
         OLD_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib:"
         NEW_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_trace")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare"
         RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/bin/oft_poincare")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare"
         OLD_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib:"
         NEW_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/oft_poincare")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs"
         RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/bin/marklin_eigs")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs"
         OLD_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib:"
         NEW_RPATH "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/libxml2-v2_15_2/lib:/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/hdf5-1_14_6/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/marklin_eigs")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  
    foreach(python_bin ThinCurr_compute_holes.py;plot_tokamaker_psi.py;plot_mug_hist.py;convert_gmsh.py;convert_cubit.py;convert_hist.py;build_xdmf-legacy.py;build_xdmf.py)
      file(RELATIVE_PATH REL_PATH ${CMAKE_INSTALL_PREFIX}/bin ${CMAKE_INSTALL_PREFIX}/python/${python_bin})
      file(CREATE_LINK ${REL_PATH} ${CMAKE_INSTALL_PREFIX}/bin/${python_bin} SYMBOLIC)
    endforeach()
    
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/bin/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
