# FindOFT_NETCDF.cmake
#
# Finds the NetCDF library for use with the Open FUSION Toolkit
# as installed by version "4.6.2" of NetCDF and
# version "4.4.4" NetCDF-Fortran.
#
# This will define the following variables
#
#    OFT_NETCDF_FOUND
#    OFT_NETCDF_INCLUDE_DIRS
#    OFT_NETCDF_LIBRARIES
#
# Author: Chris Hansen
set(FIND_FORTRAN OFF)
foreach(component IN LISTS OFT_NETCDF_FIND_COMPONENTS)
  if(component STREQUAL "Fortran") # Search for Fortran interface?
    set(FIND_FORTRAN ON)
  else()
    message(FATAL_ERROR "${component} is not a valid NETCDF component.")
  endif()
endforeach()

# Main library (C)
find_library(NETCDF_LIBRARY
  NAMES netcdf
)
find_path(NETCDF_INCLUDE_DIR
  NAMES netcdf.h
)
set(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
set(NETCDF_INCLUDES ${NETCDF_INCLUDE_DIR})

# Fortran interface
if(FIND_FORTRAN)
  find_library(NETCDF_F90_LIBRARY
    NAMES netcdff
  )
  find_path(NETCDF_F90_MODULE_DIR
    NAMES netcdf.mod
  )
  list(PREPEND NETCDF_LIBRARIES ${NETCDF_F90_LIBRARY})
  list(APPEND NETCDF_INCLUDES ${NETCDF_F90_MODULE_DIR})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_NETCDF
  REQUIRED_VARS
    NETCDF_LIBRARIES
    NETCDF_INCLUDES
)

if(OFT_NETCDF_FOUND)
  set(OFT_NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDES})
  set(OFT_NETCDF_LIBRARIES ${NETCDF_LIBRARIES})
endif()