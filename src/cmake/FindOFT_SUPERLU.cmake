# FindOFT_SUPERLU.cmake
#
# Finds the SuperLU library for use with the Open FUSION Toolkit
# as installed by version "4.3+" of SuperLU.
#
# This will define the following variables
#
#    OFT_SUPERLU_FOUND
#    OFT_SUPERLU_INCLUDE_DIRS
#    OFT_SUPERLU_LIBRARIES
#
# Author: Chris Hansen

find_library(SUPERLU_LIBRARY
  NAMES superlu superlu_4.3
)

find_path(SUPERLU_INCLUDE_DIR
  NAMES slu_ddefs.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_SUPERLU
  REQUIRED_VARS
    SUPERLU_INCLUDE_DIR
    SUPERLU_LIBRARY
)

if(OFT_SUPERLU_FOUND)
  set(OFT_SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
  set(OFT_SUPERLU_LIBRARIES ${SUPERLU_LIBRARY})
  file(READ "${SUPERLU_INCLUDE_DIR}/slu_util.h" ver_file)
  string(REGEX MATCH "SUPERLU_MAJOR_VERSION[ ]*([0-9]+)" _ ${ver_file})
  set(SUPERLU_VER_MAJOR ${CMAKE_MATCH_1})
  if(SUPERLU_VER_MAJOR)
    set(OFT_SUPERLU_VER_MAJOR ${SUPERLU_VER_MAJOR})
  else()
    message( STATUS "Found SuperLU but could not determine version" )
  endif()
endif()