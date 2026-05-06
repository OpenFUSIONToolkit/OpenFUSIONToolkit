# Add example executables
if( OFT_BUILD_EXAMPLES )
    add_subdirectory( TokaMaker )
    add_subdirectory( ThinCurr )
    add_subdirectory( MUG )
    add_subdirectory( Marklin )
    #
    oft_add_exe( example1.F90 )
    oft_add_exe( example2.F90 )
    oft_add_exe( example3.F90 )
    oft_add_exe( example4.F90 )
else()
    if( OFT_BUILD_PYTHON )
        add_subdirectory( TokaMaker )
        add_subdirectory( ThinCurr )
        add_subdirectory( Marklin )
    endif()
endif()

if( OFT_BUILD_TESTS AND (OFT_BUILD_EXAMPLES OR OFT_BUILD_PYTHON) )
    # Create test driver
    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/test_examples.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh @ONLY )
    file(
        CHMOD ${CMAKE_CURRENT_BINARY_DIR}/test_examples.sh
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ
    )
endif()