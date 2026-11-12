if( OFT_BUILD_TESTS )

  # Tests for each physics tool
  add_subdirectory( TokaMaker )
  add_subdirectory( ThinCurr )
  add_subdirectory( Marklin )
  add_subdirectory( MUG )

  # Axisymmetric Green's function tests (shared by TokaMaker and ThinCurr)
  oft_add_test( test_axi_green.F90 )
  configure_file( test_axi_green.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )

  # Combined tests (Marklin, ThinCurr, and TokaMaker in one process)
  configure_file( test_python_comb.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( Marklin/cyl_Marklin.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( ThinCurr/tw_test-plate.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
endif()
