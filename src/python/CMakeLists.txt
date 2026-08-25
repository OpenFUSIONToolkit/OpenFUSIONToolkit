add_subdirectory( OpenFUSIONToolkit )
add_subdirectory( wrappers )

# Python scripts
set( OFT_PY_BINS
  OFT_build_xdmf.py
  OFT_convert_hist.py
  OFT_convert_cubit.py
  OFT_convert_gmsh.py
  OFT_ThinCurr_holes.py
  OFT_ThinCurr_mesh_tool.py
)

set( PY_BIN_OUTPUTS "" )
foreach( path ${OFT_PY_BINS} )
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${path}"
    COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/create_stub.py ${path}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${path}" DESTINATION python COMPONENT app)
  list(APPEND PY_BIN_OUTPUTS "${CMAKE_CURRENT_BINARY_DIR}/${path}")
  set( OFT_PYTHON_BINS ${path} ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
endforeach()
add_custom_target(python_scripts ALL DEPENDS ${PY_BIN_OUTPUTS})