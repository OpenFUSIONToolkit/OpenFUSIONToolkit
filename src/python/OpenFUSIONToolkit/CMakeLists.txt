add_subdirectory( TokaMaker )
add_subdirectory( ThinCurr )
add_subdirectory( Marklin )

file(GLOB MY_SOURCES "*.py" )
foreach(SOURCE_FILE IN LISTS MY_SOURCES)
  configure_file( ${SOURCE_FILE} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  install(PROGRAMS ${SOURCE_FILE} DESTINATION python/OpenFUSIONToolkit COMPONENT app)
endforeach()