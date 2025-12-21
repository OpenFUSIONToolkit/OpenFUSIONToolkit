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

# Find pytest
if( NOT OFT_PYTEST_EXE )
  find_program( OFT_PYTEST_EXE NAMES pytest py.test pytest-3 py.test-3 )
endif()
if( OFT_PYTEST_EXE )
  message( STATUS "Found pytest: ${OFT_PYTEST_EXE}" )
else()
  message( FATAL_ERROR "pytest not found but required for automated tests" )
endif()

# Add sub-directories and copy testing script
configure_file( oft_testing.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
configure_file( conftest.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
add_subdirectory( base )
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
    POST_BUILD
    COMMAND PATH=${OFT_TEST_PATH} OFT_HAVE_MPI=${OFT_MPI_TEST} OFT_DEBUG_TEST=${OFT_DEBUG_TEST} ${OFT_PYTEST_EXE} --junitxml=OFT.junit.xml -m \"not slow\" base grid lin_alg fem physics
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    POST_BUILD
)
# Add "full" testing command to Makefile
add_custom_target(test_full)
add_custom_command(
    TARGET  test_full
    POST_BUILD
    COMMAND PATH=${OFT_TEST_PATH} OFT_HAVE_MPI=${OFT_MPI_TEST} OFT_DEBUG_TEST=${OFT_DEBUG_TEST} ${OFT_PYTEST_EXE} --junitxml=OFT.junit.xml base grid lin_alg fem physics
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    POST_BUILD
)
# Add "full" testing command to Makefile
if( OFT_COVERAGE )
  configure_file( .coveragerc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY )
  add_custom_target(test_cov)
  add_custom_command(
      TARGET  test_cov
      POST_BUILD
      COMMAND PATH=${OFT_TEST_PATH} COVERAGE_PROCESS_START=.coveragerc OFT_HAVE_MPI=${OFT_MPI_TEST} OFT_DEBUG_TEST=1 coverage run -m pytest -m \"coverage\" base grid lin_alg fem physics
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
      POST_BUILD
  )
endif()
# Create test driver
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/run_test.sh @ONLY )
file(
  CHMOD ${CMAKE_CURRENT_BINARY_DIR}/run_test.sh
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
)
