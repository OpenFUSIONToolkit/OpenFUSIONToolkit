########################
# Open FUSION Toolkit
########################
cmake_minimum_required( VERSION 3.12 )
project( OpenFUSIONToolkit LANGUAGES NONE VERSION "1.0.0" )
set( OFT_VER_SUFFIX "-beta5" CACHE INTERNAL "Version suffix" )
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake" )

option(OFT_BUILD_TESTS "Build tests?" OFF)
option(OFT_BUILD_EXAMPLES "Build Examples?" OFF)
option(OFT_USE_MPI "Build with MPI support?" OFF)
option(OFT_MPI_HEADER "Use MPI header interface instead of F08?" OFF)
option(OFT_HDF5_HL "Search for and link with HDF5 HL interface?" OFF)
option(OFT_USE_OpenMP "Build with OpenMP support?" ON)
option(OFT_DEBUG_STACK "Build with internal debugging stack?" OFF)
option(OFT_PROFILING "Build with internal profiling enabled?" OFF)
option(OFT_BUILD_PYTHON "Build experimental Python wrappers?" OFF)
option(OFT_BUILD_DOCS "Build documentation?" OFF)
option(OFT_DOCS_ONLY "Build documentation only?" OFF)
option(OFT_PACKAGE_BUILD "Perform a destributable package?" OFF)
option(OFT_PACKAGE_NIGHTLY "Perform a nightly build and include git hash in package filename?" ON)
option(OFT_COVERAGE "Build with code coverage information?" OFF)
option(OFT_TOKAMAKER_LEGACY "Build legacy TokaMaker executables?" OFF)
set(OFT_MPI_PLEN "4" CACHE STRING "Size of strings for MPI processor index")
set(OFT_PATH_SLEN "200" CACHE STRING "Size of path strings")
set(OFT_ERROR_SLEN "200" CACHE STRING "Size of error strings")
set(OFT_SLEN "80" CACHE STRING "General string length")

# Internal variables
set(OFT_SHARED_LIBS CACHE INTERNAL "package_shared_libs")
set(OFT_PYTHON_BINS CACHE INTERNAL "package_python_binaries")
if(OFT_PACKAGE_NIGHTLY)
  execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_REV_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set( OFT_VER_SUFFIX "-${GIT_REV_HASH}" CACHE INTERNAL "Version suffix" )
endif()

if(OFT_BUILD_DOCS OR OFT_PACKAGE_BUILD OR OFT_DEBUG_STACK)
  find_package(Python REQUIRED)
endif()

########################
# Build documentation
########################
if(OFT_BUILD_DOCS)
  find_package(Doxygen REQUIRED)

  file(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/docs/generated)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs/generated)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/mainpage.in ${CMAKE_CURRENT_SOURCE_DIR}/docs/generated/mainpage.md @ONLY)

  add_custom_target( build_ex_docs
    COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/utilities/generate_doc.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building documentation from examples"
    VERBATIM 
  )

  add_custom_target( docs ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Doxygen"
    VERBATIM
  )
  add_dependencies(docs build_ex_docs)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/doc/Documentation.html"
    "<meta http-equiv=\"REFRESH\" content=\"0;URL=html/index.html\">"
  )

  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc"
    DESTINATION .
    COMPONENT docs
  )

  if(OFT_DOCS_ONLY)
    return()
  endif()
endif()

########################
# Compiler and flag setup
########################
enable_language( C )
enable_language( CXX )
enable_language( Fortran )
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  if(OFT_COVERAGE)
    message( FATAL_ERROR "Coverage only supported with GCC compilers" )
  endif()
  set( CMAKE_C_FLAGS_DEBUG "-g -O0" )
  set( CMAKE_CXX_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} )
  set( CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} )
  set( CMAKE_C_FLAGS_RELEASE "" )
  set( CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} )
  set( CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} )
  if(NOT OFT_PACKAGE_BUILD)
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check bounds,pointers,shape,uninit" )
  endif()
  # Needed for new oneAPI Intel compilers
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "IntelLLVM")
    set( EXT_LIBS "-lstdc++" )
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if(OFT_COVERAGE)
    set( CMAKE_C_FLAGS_DEBUG "-Og --coverage -fprofile-update=atomic -g -Wuninitialized" ) # -fprofile-arcs -ftest-coverage
    set( CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_DEBUG} )
    add_definitions( -DOFT_COVERAGE )
    set(CMAKE_C_OUTPUT_EXTENSION_REPLACE ON)
    set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
    set(CMAKE_Fortran_OUTPUT_EXTENSION_REPLACE ON)
  else()
    set( CMAKE_C_FLAGS_DEBUG "-g -O0 -Wuninitialized" )
    set( CMAKE_C_FLAGS_RELEASE "-O2" )
  endif()
  set( CMAKE_CXX_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} )
  set( CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Werror=line-truncation -Wconversion" )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Werror=line-truncation -Wconversion" )
  if(NOT OFT_PACKAGE_BUILD)
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all" )
  endif()
  set( CMAKE_Fortran_FLAGS "-ffree-line-length-none ${CMAKE_Fortran_FLAGS}" )
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set( CMAKE_Fortran_FLAGS "-fallow-argument-mismatch ${CMAKE_Fortran_FLAGS}" )
  endif()
endif()
# Internal profiling
if( OFT_PROFILING )
  add_definitions( -DOFT_STACK -DOFT_PROFILE )
endif()
# Internal debugging stack
if( OFT_DEBUG_STACK )
  add_definitions( -DOFT_STACK )
endif()
# OpenMP
find_package( OpenMP )
# set( EXT_LIBS "${OpenMP_Fortran_LIBRARIES}")
if( OFT_USE_OpenMP )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" )
endif()
# Use fPIC when building shared library
if(OFT_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
# Use build runtime paths for install unless packaging
if(OFT_PACKAGE_BUILD)
  if(NOT APPLE)
    set(CMAKE_INSTALL_RPATH $ORIGIN)
  endif()
else()
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

########################
# Required library detection and setup
########################
# MPI
if( OFT_USE_MPI )
  find_package( MPI COMPONENTS C Fortran REQUIRED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  set( CMAKE_C_COMPILER ${MPI_C_COMPILER} )
  add_definitions( -DHAVE_MPI )
  if( MPI_Fortran_HAVE_F08_MODULE )
    if( NOT OFT_MPI_HEADER )
      add_definitions( -DOFT_MPI_F08 )
      message(STATUS "MPI Fortran 2008 interface found!")
    else()
      message(STATUS "Using header-based MPI interface")
    endif()
  else()
    message(STATUS "MPI Fortran 2008 interface not found! Falling back to header-based")
  endif()
endif()
# BLAS/LAPACK
find_package( LAPACK REQUIRED )
list(REMOVE_DUPLICATES LAPACK_LIBRARIES) # Sometimes BLAS and LAPACK are in the same library
set( EXT_LIBS ${LAPACK_LIBRARIES} ${EXT_LIBS} )
if( BLA_VENDOR MATCHES "^Intel" )
  add_definitions( -DHAVE_MKL )
endif()
# ZLIB
find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )
set( OFT_ZLIB_LIBS ${ZLIB_LIBRARIES} )
# LIBAEC
find_package( libaec REQUIRED )
include_directories( ${LIBAEC_INCLUDE_DIRS} )
set( OFT_LIBAEC_LIBS ${LIBAEC_LIBRARIES} )

# HDF5
if( OFT_HDF5_HL )
  find_package( HDF5 COMPONENTS C Fortran HL REQUIRED )
  set( OFT_HDF5_LIBS ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} )
else()
  find_package( HDF5 COMPONENTS C Fortran REQUIRED )
  set( OFT_HDF5_LIBS ${HDF5_LIBRARIES} )
endif()
include_directories( ${HDF5_INCLUDE_DIRS} )
if( DEFINED HDF5_C_LIBRARY_z )
  set( OFT_HDF5_LIBS ${OFT_HDF5_LIBS} ${HDF5_C_LIBRARY_z} )
endif()
list(REVERSE OFT_HDF5_LIBS)
list(REMOVE_DUPLICATES OFT_HDF5_LIBS) # HDF5 often provides several duplicates
list(REVERSE OFT_HDF5_LIBS)
set( EXT_LIBS ${OFT_HDF5_LIBS} ${EXT_LIBS} )
if( HDF5_IS_PARALLEL )
  if( MPI_FOUND )
    add_definitions( -DHAVE_PHDF5 )
  else()
    message( FATAL_ERROR "Parallel HDF5 found, but MPI not available" )
  endif()
endif()
# PETSc
if( OFT_PETSc_ROOT ) # Only search if PETSc root is specified
  find_package( OFT_PETSc )
endif()
if( NOT OFT_PETSc_FOUND )
  # METIS
  find_package( OFT_METIS REQUIRED )
  include_directories( ${OFT_METIS_INCLUDE_DIRS} )
  set( EXT_LIBS ${OFT_METIS_LIBRARIES} ${EXT_LIBS} )
endif()

########################
# Optional library detection and setup
########################
# ARPACK
if( OFT_USE_MPI )
  find_package( OFT_ARPACK COMPONENTS MPI )
else()
  find_package( OFT_ARPACK )
endif()
if ( OFT_ARPACK_FOUND )
  set( EXT_LIBS ${OFT_ARPACK_LIBRARIES} ${EXT_LIBS} )
  add_definitions( -DHAVE_ARPACK )
endif()
# NETCDF
find_package( OFT_NETCDF COMPONENTS Fortran )
if ( OFT_NETCDF_FOUND )
  include_directories( ${OFT_NETCDF_INCLUDE_DIRS} )
  set( EXT_LIBS ${OFT_NETCDF_LIBRARIES} ${EXT_LIBS} )
  add_definitions( -DHAVE_NCDF )
endif()
# FoX XML
find_package( OFT_FoX )
if ( OFT_FoX_FOUND )
  include_directories( ${OFT_FoX_INCLUDE_DIRS} )
  set( EXT_LIBS ${OFT_FoX_LIBRARIES} ${EXT_LIBS} )
  add_definitions( -DHAVE_XML )
endif()
# OpenNURBS
find_package( OFT_OpenNURBS )
if ( OFT_OpenNURBS_FOUND )
  include_directories( ${OFT_OpenNURBS_INCLUDE_DIRS} )
  set( EXT_LIBS ${OFT_OpenNURBS_LIBRARIES} ${EXT_LIBS} )
  add_definitions( -DHAVE_ONURBS )
endif()
# PETSc
if( OFT_PETSc_FOUND )
  add_definitions( -DHAVE_PETSC )
  add_definitions( -DPETSC_VERSION_MAJOR=${OFT_PETSc_VER_MAJOR} )
  add_definitions( -DPETSC_VERSION_MINOR=${OFT_PETSc_VER_MINOR} )
  include_directories( ${OFT_PETSc_INCLUDE_DIRS} )
  set( EXT_LIBS ${OFT_PETSc_LIBRARIES} ${EXT_LIBS} )
else() # Check for packages also provided by PETSc
  find_package( OFT_SUPERLU )
  find_package( OFT_SUPERLU_DIST )
  find_package( OFT_UMFPACK )
endif()
# SuperLU
if( OFT_SUPERLU_FOUND )
  add_definitions( -DHAVE_SUPERLU )
  if( OFT_SUPERLU_VER_MAJOR GREATER 4 )
    add_definitions( -DSUPERLU_VER_MAJOR=${OFT_SUPERLU_VER_MAJOR} )
  endif()
  if( NOT OFT_PETSc_FOUND )
    include_directories( ${OFT_SUPERLU_INCLUDE_DIRS} )
    set( EXT_LIBS ${OFT_SUPERLU_LIBRARIES} ${EXT_LIBS} )
  endif()
endif()
# SuperLU-DIST
if( OFT_SUPERLU_DIST_FOUND )
  add_definitions( -DHAVE_SUPERLU_DIST )
  if( NOT OFT_PETSc_FOUND )
    include_directories( ${OFT_SUPERLU_DIST_INCLUDE_DIRS} )
    set( EXT_LIBS ${OFT_SUPERLU_DIST_LIBRARIES} ${EXT_LIBS} )
  endif()
endif()
# UMFPACK
if( OFT_UMFPACK_FOUND )
  add_definitions( -DHAVE_UMFPACK )
  if( NOT OFT_PETSc_FOUND )
    include_directories( ${OFT_UMFPACK_INCLUDE_DIRS} )
    set( EXT_LIBS ${OFT_UMFPACK_LIBRARIES} ${EXT_LIBS} )
  endif()
endif()
# MUMPS
if( OFT_MUMPS_FOUND )
  add_definitions( -DHAVE_MUMPS )
endif()

########################
# Setup component libraries
########################
# Setup component flags
if( OFT_TOKAMAKER_LEGACY )
  add_definitions( -DOFT_TOKAMAKER_LEGACY )
endif()
# Include files
add_subdirectory( include )
include_directories( ${CMAKE_BINARY_DIR}/include )
# 3rd party libraries
add_subdirectory( ext_libs )
add_dependencies( oftext oftheaders )
include_directories( ${CMAKE_BINARY_DIR}/ext_libs )
# Base runtime
add_subdirectory( base )
add_dependencies( oftbase oftext )
include_directories( ${CMAKE_BINARY_DIR}/base )
# Grid
add_subdirectory( grid )
add_dependencies( oftgrid oftbase )
include_directories( ${CMAKE_BINARY_DIR}/grid )
# Linear algebra
add_subdirectory( lin_alg )
add_dependencies( oftla oftgrid )
include_directories( ${CMAKE_BINARY_DIR}/lin_alg )
# Finite elements
add_subdirectory( fem )
add_dependencies( oftfem oftla )
include_directories( ${CMAKE_BINARY_DIR}/fem )
# Physics
add_subdirectory( physics )
add_dependencies( oftphysics oftfem )
include_directories( ${CMAKE_BINARY_DIR}/physics )
set( OFT_LIBS oftphysics oftfem oftla oftgrid oftbase oftext )
# Utilities
add_subdirectory( utilities )
# Build master lib for packaging to reduce size
if(OFT_PACKAGE_BUILD)
  if (Python_VERSION_MAJOR VERSION_LESS 3)
    message( FATAL_ERROR "Python 3 required for packaging" )
  endif ()
  file(WRITE null.c "")
  add_library( oft SHARED null.c )
  set_target_properties( oft PROPERTIES LINKER_LANGUAGE Fortran )
  set( EXT_STATIC_LIBS "" )
  set( EXT_SHARED_LIBS "" )
  foreach(LIB IN LISTS EXT_LIBS)
    if( LIB MATCHES ".a$" )
      set( EXT_STATIC_LIBS ${LIB} ${EXT_STATIC_LIBS} )
      #message(STATUS "Static: ${LIB}")
    else()
      set( EXT_SHARED_LIBS ${LIB} ${EXT_SHARED_LIBS} )
      #message(STATUS "Shared: ${LIB}")
    endif()
  endforeach()
  if (WIN32)
    target_link_libraries( oft PRIVATE ${OFT_LIBS} ${EXT_STATIC_LIBS} )
    set_target_properties( oft PROPERTIES
        LINK_FLAGS "/WHOLEARCHIVE"
    )
  elseif (APPLE)
    foreach(LIB IN LISTS OFT_LIBS EXT_STATIC_LIBS)
      target_link_libraries( oft PRIVATE "-Wl,-force_load" ${LIB} )
    endforeach()
    target_link_libraries( oft PRIVATE ${EXT_SHARED_LIBS} )
    set( EXT_SHARED_LIBS "" )
  else ()
    target_link_libraries( oft PRIVATE "-Wl,-whole-archive" ${OFT_LIBS} ${EXT_STATIC_LIBS} "-Wl,-no-whole-archive" )
    # set_target_properties(oft PROPERTIES
    #     LINK_FLAGS "-Wl,--allow-multiple-definition"
    # )
  endif ()
  install(
    TARGETS oft
    DESTINATION bin
    COMPONENT app
  )
  set( OFT_SHARED_LIBS oft ${OFT_SHARED_LIBS} CACHE INTERNAL "package_shared_libs" )
  set( OFT_LIBS oft )
  set( EXT_LIBS ${EXT_SHARED_LIBS} )
endif()

########################
# Python interface
########################
if( OFT_BUILD_PYTHON )
  # Helper function
  function(oft_add_shared)
    set(OBJS ${ARGV})
    list(POP_FRONT OBJS)
    add_library( ${ARGV0} SHARED ${OBJS} )
    set_target_properties( ${ARGV0} PROPERTIES LINKER_LANGUAGE Fortran )
    target_link_libraries( ${ARGV0} ${OFT_LIBS} ${EXT_LIBS} )
    install(TARGETS ${ARGV0}
      DESTINATION bin
      COMPONENT app)
    set(OFT_SHARED_LIBS ${ARGV0} ${OFT_SHARED_LIBS} CACHE INTERNAL "package_shared_libs")
  endfunction(oft_add_shared)
  # Library directories
  add_subdirectory( python/wrappers )
  add_subdirectory( python/OpenFUSIONToolkit )
endif()

########################
# Build executables
########################
# Helper function
function(oft_add_exe)
  string( REPLACE ".F90" "" name ${ARGV0} )
  add_executable( ${name} ${ARGV0} )
  set_target_properties( ${name} PROPERTIES LINKER_LANGUAGE Fortran )
  target_link_libraries( ${name} ${OFT_LIBS} ${EXT_LIBS} )
  install(TARGETS ${name} COMPONENT app)
endfunction()
# Example directories
add_subdirectory( examples )
# Testing
if( OFT_BUILD_TESTS )
  add_subdirectory( tests )
endif()
# Binary directories
add_subdirectory( bin )

########################
# Packaging
########################
if(OFT_PACKAGE_BUILD)
  include(InstallRequiredSystemLibraries)
  set(CPACK_STRIP_FILES      ON)
  set(CPACK_PACKAGE_NAME     ${CMAKE_PROJECT_NAME} )
  set(CPACK_PACKAGE_VERSION  ${CMAKE_PROJECT_VERSION}${OFT_VER_SUFFIX})
  set(CPACK_PACKAGE_CONTACT  "Chris Hansen")
  set(CPACK_PACKAGE_VENDOR   "Chris Hansen")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Open FUSION Toolkit for plasma equilibrium and dynamic simulation.")
  set(CPACK_PACKAGE_DESCRIPTION ${CPACK_PACKAGE_DESCRIPTION_SUMMARY})

  set(CPACK_RESOURCE_FILE_LICENSE     ${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE)
  set(CPACK_RESOURCE_FILE_README      ${CMAKE_CURRENT_SOURCE_DIR}/../README.md)

  if(APPLE)
    set(CPACK_SYSTEM_NAME   "MacOS" )
    set(CPACK_GENERATOR     "productbuild" )
    set(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/OpenFUSIONToolkit")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/fix_signatures.sh.in ${CMAKE_CURRENT_BINARY_DIR}/fix_signatures.sh @ONLY)
    set(CPACK_POSTFLIGHT_APP_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/fix_signatures.sh)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE ${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../README.md ${CMAKE_CURRENT_BINARY_DIR}/README.txt COPYONLY)
    set(CPACK_RESOURCE_FILE_WELCOME     ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MACOS_WELCOME.txt)
    set(CPACK_RESOURCE_FILE_LICENSE     ${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt)
    set(CPACK_RESOURCE_FILE_README      ${CMAKE_CURRENT_BINARY_DIR}/README.txt)
  elseif(UNIX AND NOT APPLE)
    set(CPACK_SYSTEM_NAME   "Linux" )
    set(CPACK_GENERATOR       "TGZ" )
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE ${CMAKE_CURRENT_BINARY_DIR}/LICENSE COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../README.md ${CMAKE_CURRENT_BINARY_DIR}/README.md COPYONLY)
  endif()
  set(CPACK_PACKAGE_FILE_NAME OpenFUSIONToolkit_v${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_HOST_SYSTEM_PROCESSOR})

  include(CPack)
  cpack_add_component( app REQUIRED )
  cpack_add_component( docs )
  cpack_add_component( examples )
endif()
