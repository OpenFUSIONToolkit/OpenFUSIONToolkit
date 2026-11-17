# FindOFT_STRUMPACK.cmake
#
# Finds the STRUMPACK library for use with the Open FUSION Toolkit
# as installed by version "8.0+" of STRUMPACK. Only builds without
# MPI and GPU support are accepted (local factorization only).
#
# Only the STRUMPACK library itself is linked, relying on OFT's own
# OpenMP, METIS, and BLAS/LAPACK. Builds that require other third party
# libraries (eg. ZFP, SCOTCH) are rejected. The CMake package installed
# by STRUMPACK is not used as it re-enables languages, which conflicts
# with OFT's use of the MPI compiler wrappers.
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
  # Check for unsupported features and dependencies not linked by OFT
  if(EXISTS "${STRUMPACK_INCLUDE_DIR}/StrumpackConfig.h")
    file(READ "${STRUMPACK_INCLUDE_DIR}/StrumpackConfig.h" config_file)
    foreach(feature MPI CUDA HIP SYCL)
      if(config_file MATCHES "\n#define STRUMPACK_USE_${feature}[ \n]")
        message( FATAL_ERROR "STRUMPACK built with ${feature} support, which is not supported by OFT" )
      endif()
    endforeach()
    foreach(feature SCOTCH PTSCOTCH PARMETIS ZFP SZ3 BPACK SLATE_SCALAPACK COMBBLAS MAGMA KBLAS PAPI MATLAB)
      if(config_file MATCHES "\n#define STRUMPACK_USE_${feature}[ \n]")
        message( FATAL_ERROR "STRUMPACK built with ${feature} support, which is not supported by OFT" )
      endif()
    endforeach()
    if((config_file MATCHES "\n#define STRUMPACK_USE_OPENMP[ \n]") AND (NOT OFT_USE_OpenMP))
      message( FATAL_ERROR "STRUMPACK built with OpenMP support, which requires OFT_USE_OpenMP" )
    endif()
  else()
    message( STATUS "Found STRUMPACK but could not find \"StrumpackConfig.h\" to check configuration" )
  endif()
endif()
