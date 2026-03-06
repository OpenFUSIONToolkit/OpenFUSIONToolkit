# FindOFT_LibXml2.cmake
#
# Finds the LibXml2 library for use with the Open FUSION Toolkit.
#
# This will define the following variables
#
#    OFT_LibXml2_FOUND
#    OFT_LibXml2_INCLUDE_DIRS
#    OFT_LibXml2_LIBRARIES
#
# Author: Chris Hansen

find_package(LibXml2)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_LibXml2
  REQUIRED_VARS
    LIBXML2_LIBRARIES
    LIBXML2_INCLUDE_DIR
)

if(OFT_LibXml2_FOUND)
  set(OFT_LibXml2_INCLUDE_DIRS ${LIBXML2_INCLUDE_DIR})
  set(OFT_LibXml2_LIBRARIES ${LIBXML2_LIBRARIES})
endif()
