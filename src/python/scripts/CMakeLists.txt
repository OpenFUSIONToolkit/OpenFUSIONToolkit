# Python scripts
set( OFT_PY_BINS
  OFT_build_xdmf.py
  OFT_convert_hist.py
  OFT_convert_cubit.py
  OFT_convert_gmsh.py
  OFT_ThinCurr_holes.py
)

foreach( path ${OFT_PY_BINS} )
  install(PROGRAMS ${path} DESTINATION python COMPONENT app)
  get_filename_component(file "${path}" NAME)
  set( OFT_PYTHON_BINS ${file} ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
  configure_file( ${path} ${CMAKE_BINARY_DIR}/python/${file} COPYONLY )
endforeach()