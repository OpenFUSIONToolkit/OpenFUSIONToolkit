# FindOFT_FoX.cmake
#
# Finds the FoX library for use with the Open FUSION Toolkit
# as installed by version "4.1" of FoX.
#
# This will define the following variables
#
#    OFT_FoX_FOUND
#    OFT_FoX_INCLUDE_DIR
#    OFT_FoX_LIBRARIES
#
# Author: Chris Hansen

find_library(FoX_DOM_LIBRARY
  NAMES FoX_dom
)
find_library(FoX_SAX_LIBRARY
  NAMES FoX_sax
)
find_library(FoX_FSYS_LIBRARY
  NAMES FoX_fsys
)
find_library(FoX_WXML_LIBRARY
  NAMES FoX_wxml
)
find_library(FoX_UTILS_LIBRARY
  NAMES FoX_utils
)
find_library(FoX_COMMON_LIBRARY
  NAMES FoX_common
)

find_path(FoX_MODULE_DIR
  NAMES fox_dom.mod
  PATH_SUFFIXES "finclude"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OFT_FoX
  REQUIRED_VARS
    FoX_MODULE_DIR
    FoX_DOM_LIBRARY
    FoX_SAX_LIBRARY
    FoX_FSYS_LIBRARY
    FoX_WXML_LIBRARY
    FoX_UTILS_LIBRARY
    FoX_COMMON_LIBRARY
)

if(OFT_FoX_FOUND)
  set(OFT_FoX_INCLUDE_DIRS ${FoX_MODULE_DIR})
  set(OFT_FoX_LIBRARIES ${FoX_DOM_LIBRARY} ${FoX_SAX_LIBRARY}
    ${FoX_FSYS_LIBRARY} ${FoX_WXML_LIBRARY} ${FoX_UTILS_LIBRARY} ${FoX_COMMON_LIBRARY})
endif()