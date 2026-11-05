# FindOFT_STRUMPACK.cmake
#
# Finds the STRUMPACK library for use with the Open FUSION Toolkit
# as installed by version "8.0+" of STRUMPACK. Only builds without
# MPI and GPU support are accepted (local factorization only).
#
# This will define the following variables
#
#    OFT_STRUMPACK_FOUND
#    OFT_STRUMPACK_INCLUDE_DIRS
#    OFT_STRUMPACK_LIBRARIES
#
# Author: Chris Hansen

find_library(STRUMPACK_LIBRARY
  NAMES strumpack
)

find_path(STRUMPACK_INCLUDE_DIR
  NAMES StrumpackSparseSolver.hpp
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_STRUMPACK
  REQUIRED_VARS
    STRUMPACK_INCLUDE_DIR
    STRUMPACK_LIBRARY
)

if(OFT_STRUMPACK_FOUND)
  set(OFT_STRUMPACK_INCLUDE_DIRS ${STRUMPACK_INCLUDE_DIR})
  set(OFT_STRUMPACK_LIBRARIES ${STRUMPACK_LIBRARY})
  # Check for unsupported features
  if(EXISTS "${STRUMPACK_INCLUDE_DIR}/StrumpackConfig.h")
    file(READ "${STRUMPACK_INCLUDE_DIR}/StrumpackConfig.h" config_file)
    foreach(feature MPI CUDA HIP SYCL)
      if(config_file MATCHES "\n#define STRUMPACK_USE_${feature}[ \n]")
        message( FATAL_ERROR "STRUMPACK built with ${feature} support, which is not supported by OFT" )
      endif()
    endforeach()
  else()
    message( STATUS "Found STRUMPACK but could not find \"StrumpackConfig.h\" to check configuration" )
  endif()
endif()
