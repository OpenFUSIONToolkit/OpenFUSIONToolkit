# Helper function
function(oft_add_test)
  string( REPLACE ".F90" "" name ${ARGV0} )
  add_executable( ${name} ${ARGV0} )
  set_target_properties( ${name} PROPERTIES LINKER_LANGUAGE Fortran )
  if( DEFINED ARGV1 )
    set(OBJS ${ARGV})
    list(POP_FRONT OBJS)
    target_link_libraries( ${name} ${OBJS} ${OFT_LIBS} ${EXT_LIBS} )
  else()
    target_link_libraries( ${name} ${OFT_LIBS} ${EXT_LIBS} )
  endif()
endfunction(oft_add_test)

# Add sub-directories and copy testing script
file( COPY oft_testing.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_subdirectory( grid )
add_subdirectory( lin_alg )
add_subdirectory( fem )
add_subdirectory( physics )

# Get MPI path and prepend if necessary
set( OFT_TEST_PATH "$ENV{PATH}" )
set( OFT_MPI_TEST "0" )
set( OFT_DEBUG_TEST "0" )
if( OFT_USE_MPI )
    get_filename_component(MPI_BIN_PATH ${MPIEXEC_EXECUTABLE} DIRECTORY)
    set( OFT_TEST_PATH "${MPI_BIN_PATH}:$ENV{PATH}" )
    set( OFT_MPI_TEST "1" )
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set( OFT_DEBUG_TEST "1" )
    # set( OFT_TEST_ENV "${OFT_TEST_ENV};OFT_DEBUG_TEST=1")
endif()
# Add testing command to Makefile
add_custom_target(test)
add_custom_command(
    TARGET  test
    COMMAND PATH=${OFT_TEST_PATH} OFT_HAVE_MPI=${OFT_MPI_TEST} OFT_DEBUG_TEST=${OFT_DEBUG_TEST} ${OFT_PYTEST_NAME} -m \"not slow\" grid lin_alg fem physics
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
# Add "full" testing command to Makefile
add_custom_target(test_full)
add_custom_command(
    TARGET  test_full
    COMMAND PATH=${OFT_TEST_PATH} OFT_HAVE_MPI=${OFT_MPI_TEST} OFT_DEBUG_TEST=${OFT_DEBUG_TEST} ${OFT_PYTEST_NAME} grid lin_alg fem physics
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)
