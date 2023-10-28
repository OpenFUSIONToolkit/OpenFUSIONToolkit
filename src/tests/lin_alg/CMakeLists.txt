# Add test executables
if( OFT_BUILD_TESTS )
  oft_add_test( test_arpack.F90 )
  file( COPY test_arpack.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_native_bjacobi.F90 )
  file( COPY test_native_bjacobi.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_solver_xml.F90 )
  file( COPY test_solver_xml.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()