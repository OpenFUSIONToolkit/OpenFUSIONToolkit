install(DIRECTORY slab DESTINATION examples/ThickCurr COMPONENT examples)
install(DIRECTORY torus DESTINATION examples/ThickCurr COMPONENT examples)

if( OFT_BUILD_TESTS )
  file(COPY slab DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY torus DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  # Create test driver
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh @ONLY )
  file(
    CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
  )
endif()
