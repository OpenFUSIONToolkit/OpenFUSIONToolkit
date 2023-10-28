if( OFT_BUILD_TESTS )
  # Quadrature tests
  oft_add_test( test_quad.F90 )
  file( COPY test_quad.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY quad_1d.tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY quad_2d.tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY quad_3d.tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # Mapping tests
  oft_add_test( test_mapping_find.F90 )
  oft_add_test( test_mapping_jac.F90 )
  file( COPY test_mapping.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY mapping_jac.tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # T3D tests
  oft_add_test( test_t3d.F90 )
  file( COPY test_t3d.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY cyl.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY cyl.t3d DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # Cubit tests
  oft_add_test( test_cubit.F90 )
  file( COPY test_cubit.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY ref_test.3dm DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY ref_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY ref_test2.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_test.3dm DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_test2.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_test_cut.3dm DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_test_cut.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY circle_tri3_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY circle_tri6_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY circle_quad4_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY circle_quad9_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_tri3_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_tri6_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_quad4_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY sphere_quad9_test.in.e DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  # GMSH tests
  oft_add_test( test_gmsh.F90 )
  file( COPY test_gmsh.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  file( COPY cyl.mesh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()