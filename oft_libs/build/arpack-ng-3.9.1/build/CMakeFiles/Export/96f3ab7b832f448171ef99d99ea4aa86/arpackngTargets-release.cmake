#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "arpack" for configuration "Release"
set_property(TARGET arpack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(arpack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libarpack.a"
  )

list(APPEND _cmake_import_check_targets arpack )
list(APPEND _cmake_import_check_files_for_arpack "${_IMPORT_PREFIX}/lib/libarpack.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
