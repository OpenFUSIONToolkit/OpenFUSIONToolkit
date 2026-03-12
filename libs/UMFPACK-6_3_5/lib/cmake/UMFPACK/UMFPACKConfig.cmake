#-------------------------------------------------------------------------------
# SuiteSparse/UMFPACK/cmake_modules/UMFPACKConfig.cmake
#-------------------------------------------------------------------------------

# The following copyright and license applies to just this file only, not to
# the library itself:
# UMFPACKConfig.cmake, Copyright (c) 2023-2024, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# Finds the UMFPACK include file and compiled library.
# The following targets are defined:
#   SuiteSparse::UMFPACK           - for the shared library (if available)
#   SuiteSparse::UMFPACK_static    - for the static library (if available)

# For backward compatibility the following variables are set:

# UMFPACK_INCLUDE_DIR - where to find umfpack.h
# UMFPACK_LIBRARY     - dynamic UMFPACK library
# UMFPACK_STATIC      - static UMFPACK library
# UMFPACK_LIBRARIES   - libraries when using UMFPACK
# UMFPACK_FOUND       - true if UMFPACK found

# Set ``CMAKE_MODULE_PATH`` to the parent folder where this module file is
# installed.

#-------------------------------------------------------------------------------


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was UMFPACKConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set ( UMFPACK_DATE "Sept 23, 2024" )
set ( UMFPACK_VERSION_MAJOR 6 )
set ( UMFPACK_VERSION_MINOR 3 )
set ( UMFPACK_VERSION_PATCH 5 )
set ( UMFPACK_VERSION "6.3.5" )

# Check for dependent targets
include ( CMakeFindDependencyMacro )

# Look for SuiteSparse_config, AMD and potentially CHOLMOD targets
if ( OFF )
    if ( NOT TARGET SuiteSparse::SuiteSparseConfig )
        # First check in a common build tree
        find_dependency ( SuiteSparse_config 7.10
            PATHS ${CMAKE_SOURCE_DIR}/../SuiteSparse_config/build NO_DEFAULT_PATH )
        # Then, check in the currently active CMAKE_MODULE_PATH
        if ( NOT SuiteSparse_config_FOUND )
            find_dependency ( SuiteSparse_config 7.10 )
        endif ( )
    endif ( )

    if ( NOT TARGET SuiteSparse::AMD )
        # First check in a common build tree
        find_dependency ( AMD 3.3
            PATHS ${CMAKE_SOURCE_DIR}/../AMD/build NO_DEFAULT_PATH )
        # Then, check in the currently active CMAKE_MODULE_PATH
        if ( NOT AMD_FOUND )
            find_dependency ( AMD 3.3 )
        endif ( )
    endif ( )

    if ( OFF )
        if ( NOT TARGET SuiteSparse::CHOLMOD )
            # First check in a common build tree
            find_dependency ( CHOLMOD .
                PATHS ${CMAKE_SOURCE_DIR}/../CHOLMOD/build NO_DEFAULT_PATH )
            # Then, check in the currently active CMAKE_MODULE_PATH
            if ( NOT CHOLMOD_FOUND )
                find_dependency ( CHOLMOD . )
            endif ( )
        endif ( )
    endif ( )

else ( )
    if ( NOT TARGET SuiteSparse::SuiteSparseConfig )
        find_dependency ( SuiteSparse_config 7.10 )
    endif ( )
    if ( NOT TARGET SuiteSparse::AMD )
        find_dependency ( AMD 3.3 )
    endif ( )
    if ( OFF AND NOT TARGET SuiteSparse::CHOLMOD )
        find_dependency ( CHOLMOD . )
    endif ( )
endif ( )

# FIXME: Also check for BLAS libraries here?

if ( NOT SuiteSparse_config_FOUND OR NOT AMD_FOUND 
     OR ( OFF AND NOT CHOLMOD_FOUND ) )
    set ( UMFPACK_FOUND OFF )
    return ( )
endif ( )


# Import target
include ( ${CMAKE_CURRENT_LIST_DIR}/UMFPACKTargets.cmake )

# The following is only for backward compatibility with FindUMFPACK.

set ( _target_shared SuiteSparse::UMFPACK )
set ( _target_static SuiteSparse::UMFPACK_static )
set ( _var_prefix "UMFPACK" )

if ( NOT OFF AND NOT TARGET ${_target_shared} )
    # make sure there is always an import target without suffix )
    add_library ( ${_target_shared} ALIAS ${_target_static} )
endif ( )

get_target_property ( ${_var_prefix}_INCLUDE_DIR ${_target_shared} INTERFACE_INCLUDE_DIRECTORIES )
if ( ${_var_prefix}_INCLUDE_DIR )
    # First item in SuiteSparse targets contains the "main" header directory.
    list ( GET ${_var_prefix}_INCLUDE_DIR 0 ${_var_prefix}_INCLUDE_DIR )
endif ( )
get_target_property ( ${_var_prefix}_LIBRARY ${_target_shared} IMPORTED_IMPLIB )
if ( NOT ${_var_prefix}_LIBRARY )
    get_target_property ( _library_chk ${_target_shared} IMPORTED_LOCATION )
    if ( EXISTS ${_library_chk} )
        set ( ${_var_prefix}_LIBRARY ${_library_chk} )
    endif ( )
endif ( )
if ( TARGET ${_target_static} )
    get_target_property ( ${_var_prefix}_STATIC ${_target_static} IMPORTED_LOCATION )
endif ( )

# Check for most common build types
set ( _config_types "Debug" "Release" "RelWithDebInfo" "MinSizeRel" "None" )

get_property ( _isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG )
if ( _isMultiConfig )
    # For multi-configuration generators (e.g., Visual Studio), prefer those
    # configurations.
    list ( PREPEND _config_types ${CMAKE_CONFIGURATION_TYPES} )
else ( )
    # For single-configuration generators, prefer the current configuration.
    list ( PREPEND _config_types ${CMAKE_BUILD_TYPE} )
endif ( )

list ( REMOVE_DUPLICATES _config_types )

foreach ( _config ${_config_types} )
    string ( TOUPPER ${_config} _uc_config )
    if ( NOT ${_var_prefix}_LIBRARY )
        get_target_property ( _library_chk ${_target_shared}
            IMPORTED_IMPLIB_${_uc_config} )
        if ( EXISTS ${_library_chk} )
            set ( ${_var_prefix}_LIBRARY ${_library_chk} )
        endif ( )
    endif ( )
    if ( NOT ${_var_prefix}_LIBRARY )
        get_target_property ( _library_chk ${_target_shared}
            IMPORTED_LOCATION_${_uc_config} )
        if ( EXISTS ${_library_chk} )
            set ( ${_var_prefix}_LIBRARY ${_library_chk} )
        endif ( )
    endif ( )
    if ( TARGET ${_target_static} AND NOT ${_var_prefix}_STATIC )
        get_target_property ( _library_chk ${_target_static}
            IMPORTED_LOCATION_${_uc_config} )
        if ( EXISTS ${_library_chk} )
            set ( ${_var_prefix}_STATIC ${_library_chk} )
        endif ( )
    endif ( )
endforeach ( )

set ( UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} )

macro ( suitesparse_check_exist _var _files )
  # ignore generator expressions
  string ( GENEX_STRIP "${_files}" _files2 )

  foreach ( _file ${_files2} )
    if ( NOT EXISTS "${_file}" )
      message ( FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist!" )
    endif ( )
  endforeach ()
endmacro ( )

suitesparse_check_exist ( UMFPACK_INCLUDE_DIR ${UMFPACK_INCLUDE_DIR} )
suitesparse_check_exist ( UMFPACK_LIBRARY ${UMFPACK_LIBRARY} )

message ( STATUS "UMFPACK version: ${UMFPACK_VERSION}" )
message ( STATUS "UMFPACK include: ${UMFPACK_INCLUDE_DIR}" )
message ( STATUS "UMFPACK library: ${UMFPACK_LIBRARY}" )
message ( STATUS "UMFPACK static:  ${UMFPACK_STATIC}" )
