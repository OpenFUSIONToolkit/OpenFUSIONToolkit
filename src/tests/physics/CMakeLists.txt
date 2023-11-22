if( OFT_BUILD_TESTS )
  # Common
  file( COPY conftest.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # Marklin tests
  oft_add_test( test_taylor.F90 )
  file( COPY test_taylor.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_taylor_inj.F90 )
  file( COPY test_taylor_inj.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY torus_test.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # MUG tests
  add_library( test_phys_helpers test_phys_helpers.F90 )
  add_dependencies( test_phys_helpers oftphysics )
  oft_add_test( test_sound.F90 test_phys_helpers )
  file( COPY test_sound.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_sound_lag.F90 test_phys_helpers )
  file( COPY test_sound_lag.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_alfven.F90 test_phys_helpers )
  file( COPY test_alfven.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  oft_add_test( test_alfven_lag.F90 test_phys_helpers )
  file( COPY test_alfven_lag.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY cube.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY cube.t3d DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # ThinCurr tests
  file( COPY test_thin_wall.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY tw_test-cyl.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY tw_test-passive.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY tw_test-plate.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY tw_test-torus.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # TokaMaker tests
  file( COPY test_TokaMaker.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY ITER_geom.json DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()