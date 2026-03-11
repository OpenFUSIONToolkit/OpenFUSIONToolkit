install(DIRECTORY cylinder DESTINATION examples/Marklin COMPONENT examples)
install(DIRECTORY NCSX DESTINATION examples/Marklin COMPONENT examples)

if( OFT_BUILD_TESTS )
  file(COPY cylinder DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY NCSX DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  # Create test driver
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh @ONLY )
  file(
    CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
  )
endif()