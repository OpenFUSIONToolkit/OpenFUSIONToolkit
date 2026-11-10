# FindOFT_CVODE.cmake
#
# Finds the serial CVODE library from SUNDIALS (version 6.0+) for use with
# the Open FUSION Toolkit.
#
# This will define the following variables
#
#    OFT_CVODE_FOUND
#    OFT_CVODE_INCLUDE_DIRS
#    OFT_CVODE_LIBRARIES
#    OFT_CVODE_VER_MAJOR
#
# Author: Chris Hansen

find_library(CVODE_LIBRARY
  NAMES sundials_cvode
)

find_library(CVODE_NVECSERIAL_LIBRARY
  NAMES sundials_nvecserial
)

# Functional iteration (non-stiff) nonlinear solver
find_library(CVODE_FIXEDPOINT_LIBRARY
  NAMES sundials_sunnonlinsolfixedpoint
)

# Required by SUNDIALS 7.0+ (older versions bundle this in the above libraries)
find_library(CVODE_CORE_LIBRARY
  NAMES sundials_core
)

find_path(CVODE_INCLUDE_DIR
  NAMES cvode/cvode.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_CVODE
  REQUIRED_VARS
    CVODE_INCLUDE_DIR
    CVODE_LIBRARY
    CVODE_NVECSERIAL_LIBRARY
    CVODE_FIXEDPOINT_LIBRARY
)

if(OFT_CVODE_FOUND)
  set(OFT_CVODE_INCLUDE_DIRS ${CVODE_INCLUDE_DIR})
  set(OFT_CVODE_LIBRARIES ${CVODE_LIBRARY} ${CVODE_FIXEDPOINT_LIBRARY} ${CVODE_NVECSERIAL_LIBRARY})
  if(CVODE_CORE_LIBRARY)
    list(APPEND OFT_CVODE_LIBRARIES ${CVODE_CORE_LIBRARY})
  endif()
  if(UNIX)
    list(APPEND OFT_CVODE_LIBRARIES m)
  endif()
  # Determine major version for API differences
  set(OFT_CVODE_VER_MAJOR "")
  find_file(CVODE_CONFIG_HEADER
    NAMES sundials_config.h
    PATHS ${CVODE_INCLUDE_DIR}
    PATH_SUFFIXES sundials
    NO_DEFAULT_PATH
  )
  if(CVODE_CONFIG_HEADER)
    file(READ "${CVODE_CONFIG_HEADER}" ver_file)
    string(REGEX MATCH "SUNDIALS_VERSION_MAJOR[ ]*([0-9]+)" _ "${ver_file}")
    set(OFT_CVODE_VER_MAJOR ${CMAKE_MATCH_1})
  endif()
  if(NOT OFT_CVODE_VER_MAJOR)
    message(FATAL_ERROR "Found CVODE but could not determine SUNDIALS version")
  endif()
  if(OFT_CVODE_VER_MAJOR LESS 6)
    message(FATAL_ERROR "CVODE from SUNDIALS 6.0 or newer is required (found ${OFT_CVODE_VER_MAJOR}.x)")
  endif()
endif()
