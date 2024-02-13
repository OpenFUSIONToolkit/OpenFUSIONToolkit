# FindOFT_UMFPACK.cmake
#
# Finds the SuperLU library for use with the Open FUSION Toolkit
# as installed by version "4.3+" of SuperLU.
#
# This will define the following variables
#
#    OFT_UMFPACK_FOUND
#    OFT_UMFPACK_INCLUDE_DIRS
#    OFT_UMFPACK_LIBRARIES
#
# Author: Chris Hansen

find_library(UMFPACK_BASE_LIBRARY
  NAMES umfpack
)

find_library(UMFPACK_AMD_LIBRARY
  NAMES amd
)

find_library(UMFPACK_SSC_LIBRARY
  NAMES suitesparseconfig
)

find_path(UMFPACK_INCLUDE_DIR
  NAMES umfpack.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_UMFPACK
  REQUIRED_VARS
    UMFPACK_INCLUDE_DIR
    UMFPACK_BASE_LIBRARY
    UMFPACK_AMD_LIBRARY
    UMFPACK_SSC_LIBRARY
)

if(OFT_UMFPACK_FOUND)
  set(OFT_UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})
  set(OFT_UMFPACK_LIBRARIES ${UMFPACK_BASE_LIBRARY} ${UMFPACK_AMD_LIBRARY}
    ${UMFPACK_SSC_LIBRARY})
endif()