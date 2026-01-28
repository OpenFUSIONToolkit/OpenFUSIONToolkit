# FindOFT_METIS.cmake
#
# Finds the metis library for use with the Open FUSION Toolkit
# as installed by version "5.1" of metis.
#
# This will define the following variables
#
#    OFT_METIS_FOUND
#    OFT_METIS_INCLUDE_DIRS
#    OFT_METIS_LIBRARIES
#
# Author: Chris Hansen

find_library(METIS_LIBRARY
  NAMES metis
)

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_METIS
  REQUIRED_VARS
    METIS_LIBRARY
    METIS_INCLUDE_DIR
)

if(OFT_METIS_FOUND)
  set(OFT_METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
  set(OFT_METIS_LIBRARIES ${METIS_LIBRARY})
endif()