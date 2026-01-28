# Add test executables
if( OFT_BUILD_TESTS )
  oft_add_test( test_io.F90 )
  configure_file( test_io.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
endif()