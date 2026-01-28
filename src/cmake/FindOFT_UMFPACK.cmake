# FindOFT_UMFPACK.cmake
#
# Finds the UMFPACK library for use with the Open FUSION Toolkit
# as installed by version "7.0+" of SuiteSparse.
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

# find_library(UMFPACK_CHOLMOD_LIBRARY
#   NAMES cholmod
# )

find_library(UMFPACK_SSC_LIBRARY
  NAMES suitesparseconfig
)

find_path(UMFPACK_INCLUDE_DIR
  NAMES umfpack.h
  PATH_SUFFIXES suitesparse
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