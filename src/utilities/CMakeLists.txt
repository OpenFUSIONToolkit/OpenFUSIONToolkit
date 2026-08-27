# Python modules
set( OFT_PY_HELPER_MODS
  oft_mpl.py
)
foreach( file ${OFT_PY_HELPER_MODS} )
  configure_file( ${file} ${CMAKE_BINARY_DIR}/python/${file} COPYONLY )
  install(FILES ${file} DESTINATION python COMPONENT app)
endforeach()

# Python scripts
set( OFT_PY_BINS
  scripts/plot_mug_hist.py
  build_xdmf-legacy.py
)

foreach( path ${OFT_PY_BINS} )
  install(PROGRAMS ${path} DESTINATION python COMPONENT app)
  get_filename_component(file "${path}" NAME)
  set( OFT_PYTHON_BINS ${file} ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
  configure_file( ${path} ${CMAKE_BINARY_DIR}/python/${file} COPYONLY )
endforeach()