ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/git_info.h always_rebuild
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/create_git_header.sh > ${CMAKE_CURRENT_BINARY_DIR}/git_info.h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
configure_file( local.h ${CMAKE_CURRENT_BINARY_DIR} @ONLY )

if( OFT_DEBUG_STACK )
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/stack_defs.h always_rebuild
    COMMAND Python::Interpreter ${CMAKE_SOURCE_DIR}/utilities/generate_stack.py -o ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  ADD_CUSTOM_TARGET(oftheaders DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/stack_defs.h ${CMAKE_CURRENT_BINARY_DIR}/local.h ${CMAKE_CURRENT_BINARY_DIR}/git_info.h)
  add_custom_target(stack_clean
    COMMAND Python::Interpreter ${CMAKE_SOURCE_DIR}/utilities/generate_stack.py -c
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  )
else()
  ADD_CUSTOM_TARGET(oftheaders DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/local.h ${CMAKE_CURRENT_BINARY_DIR}/git_info.h)
endif()
