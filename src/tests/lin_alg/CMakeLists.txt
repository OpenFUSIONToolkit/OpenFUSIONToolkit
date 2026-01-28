# Add test executables
if( OFT_BUILD_TESTS )
  oft_add_test( test_arpack.F90 )
  configure_file( test_arpack.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_native_bjacobi.F90 )
  configure_file( test_native_bjacobi.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_solver_xml.F90 )
  configure_file( test_solver_xml.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
endif()