# Install script for directory: /Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/UMFPACK-6_3_5")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/libsuitesparseconfig.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuitesparseconfig.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuitesparseconfig.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuitesparseconfig.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/suitesparse" TYPE FILE FILES "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/SuiteSparse_config.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse" TYPE FILE FILES
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparseBLAS.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparseBLAS32.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparseBLAS64.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparseLAPACK.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparsePolicy.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparseReport.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/cmake_modules/SuiteSparse__thread.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config/SuiteSparse_configTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config/SuiteSparse_configTargets.cmake"
         "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/CMakeFiles/Export/03046fdaf624384fd2d3f958f353cef7/SuiteSparse_configTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config/SuiteSparse_configTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config/SuiteSparse_configTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config" TYPE FILE FILES "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/CMakeFiles/Export/03046fdaf624384fd2d3f958f353cef7/SuiteSparse_configTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config" TYPE FILE FILES "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/CMakeFiles/Export/03046fdaf624384fd2d3f958f353cef7/SuiteSparse_configTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/SuiteSparse_config" TYPE FILE FILES
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/SuiteSparse_configConfig.cmake"
    "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/SuiteSparse_configConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/SuiteSparse_config.pc")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/taliaangles/Documents/GitHub/OpenFUSIONToolkit/oft_libs/build/SuiteSparse-7.10.1/SuiteSparse_config/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
