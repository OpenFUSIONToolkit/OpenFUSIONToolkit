# FindOFT_SUPERLU_DIST.cmake
#
# Finds the SuperLU-DIST library for use with the OpenFUSIONToolkit
# as installed by version "3.3+" of SuperLU.
#
# This will define the following variables
#
#    OFT_SUPERLU_DIST_FOUND
#    OFT_SUPERLU_DIST_INCLUDE_DIRS
#    OFT_SUPERLU_DIST_LIBRARIES
#
# Author: Chris Hansen

find_library(SUPERLU_DIST_LIBRARY
  NAMES superlu_dist superlu_dist_4.1 superlu_dist_3.3
)

find_path(SUPERLU_DIST_INCLUDE_DIR
  NAMES superlu_ddefs.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_SUPERLU_DIST
  REQUIRED_VARS
    SUPERLU_DIST_INCLUDE_DIR
    SUPERLU_DIST_LIBRARY
)

if(OFT_SUPERLU_DIST_FOUND)
  set(OFT_SUPERLU_DIST_INCLUDE_DIRS ${SUPERLU_DIST_INCLUDE_DIR})
  set(OFT_SUPERLU_DIST_LIBRARIES ${SUPERLU_DIST_LIBRARY})
endif()