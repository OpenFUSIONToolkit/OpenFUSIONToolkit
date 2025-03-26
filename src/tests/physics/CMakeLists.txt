if( OFT_BUILD_TESTS )

  # Marklin tests
  oft_add_test( test_taylor.F90 )
  configure_file( test_taylor.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( test_Marklin.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( cyl_Marklin.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_taylor_inj.F90 )
  configure_file( test_taylor_inj.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( torus_test.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )

  # MUG tests
  add_library( test_phys_helpers test_phys_helpers.F90 )
  add_dependencies( test_phys_helpers oftphysics )
  oft_add_test( test_sound.F90 test_phys_helpers )
  configure_file( test_sound.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_sound_lag.F90 test_phys_helpers )
  configure_file( test_sound_lag.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_alfven.F90 test_phys_helpers )
  configure_file( test_alfven.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  oft_add_test( test_alfven_lag.F90 test_phys_helpers )
  configure_file( test_alfven_lag.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( cube.inp ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( cube.t3d ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )

  # ThinCurr tests
  configure_file( test_ThinCurr.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( tw_test-cyl.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( tw_test-cyl_hr.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( tw_test-passive.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( tw_test-plate.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( tw_test-torus.h5 ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )

  # TokaMaker tests
  configure_file( test_TokaMaker.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( ITER_geom.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( ITER_test.eqdsk ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( ITER_test.ifile ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  configure_file( LTX_geom.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )

  # Combined tests
  configure_file( test_python_comb.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
endif()