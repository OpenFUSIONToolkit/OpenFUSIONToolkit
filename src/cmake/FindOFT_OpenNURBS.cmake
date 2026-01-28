# FindOFT_OpenNURBS.cmake
#
# Finds the OpenNURBS library for use with the Open FUSION Toolkit
# as installed by version "5.0" of OpenNURBS.
#
# This will define the following variables
#
#    OFT_OpenNURBS_FOUND
#    OFT_OpenNURBS_INCLUDE_DIR
#    OFT_OpenNURBS_LIBRARIES
#
# Author: Chris Hansen

find_library(OpenNURBS_LIBRARY
  NAMES openNURBS
)

find_path(OpenNURBS_INCLUDE_DIR
  NAMES opennurbs.h
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_OpenNURBS
  REQUIRED_VARS
    OpenNURBS_INCLUDE_DIR
    OpenNURBS_LIBRARY
)

if(OFT_OpenNURBS_FOUND)
  set(OFT_OpenNURBS_INCLUDE_DIRS ${OpenNURBS_INCLUDE_DIR})
  set(OFT_OpenNURBS_LIBRARIES ${OpenNURBS_LIBRARY})
endif()