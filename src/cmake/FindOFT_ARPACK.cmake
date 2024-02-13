# FindOFT_ARPACK.cmake
#
# Finds the arpack-ng library for use with the Open FUSION Toolkit
# as installed by version "3.5" of arpack-ng.
#
# This will define the following variables
#
#    OFT_ARPACK_FOUND
#    OFT_ARPACK_LIBRARIES
#
# Author: Chris Hansen
set(FIND_PARALLEL OFF)
foreach(component IN LISTS OFT_ARPACK_FIND_COMPONENTS)
  if(component STREQUAL "MPI") # Search for parallel version?
    set(FIND_PARALLEL ON)
  else()
    message(FATAL_ERROR "${component} is not a valid ARPACK component.")
  endif()
endforeach()

find_library(ARPACK_SERIAL_LIBRARY
  NAMES arpack
)
set(ARPACK_LIBRARIES ${ARPACK_SERIAL_LIBRARY})

if(FIND_PARALLEL)
  find_library(ARPACK_PARALLEL_LIBRARY
    NAMES parpack
  )
  list(PREPEND ARPACK_LIBRARIES ${ARPACK_PARALLEL_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_ARPACK
  REQUIRED_VARS ARPACK_LIBRARIES
)

if(OFT_ARPACK_FOUND)
  set(OFT_ARPACK_LIBRARIES ${ARPACK_LIBRARIES})
endif()
