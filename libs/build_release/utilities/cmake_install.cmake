# Install script for directory: /home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE FILE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/oft_mpl.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE FILE FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/tokamaker_fit.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/build_xdmf.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/build_xdmf-legacy.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/convert_hist.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/convert_cubit.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/convert_gmsh.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/scripts/plot_mug_hist.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/scripts/plot_tokamaker_psi.py")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "app" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python" TYPE PROGRAM FILES "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/src/utilities/scripts/ThinCurr_compute_holes.py")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/OpenFUSIONToolkit/OpenFUSIONToolkit/libs/build_release/utilities/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
