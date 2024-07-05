# Python modules
set( OFT_PY_HELPER_MODS
  oft_mpl.py
  tokamaker_fit.py
)
foreach( file ${OFT_PY_HELPER_MODS} )
  install(FILES ${file} DESTINATION python COMPONENT app)
endforeach()

# Scripts
install(PROGRAMS build_xdmf.py DESTINATION python COMPONENT app)
set( OFT_PYTHON_BINS build_xdmf.py ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
install(PROGRAMS convert_hist.py DESTINATION python COMPONENT app)
set( OFT_PYTHON_BINS convert_hist.py ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
install(PROGRAMS convert_cubit.py DESTINATION python COMPONENT app)
set( OFT_PYTHON_BINS convert_cubit.py ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
install(PROGRAMS scripts/plot_mug_hist.py DESTINATION python COMPONENT app)
set( OFT_PYTHON_BINS plot_mug_hist.py ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )
install(PROGRAMS scripts/plot_tokamaker_psi.py DESTINATION python COMPONENT app)
set( OFT_PYTHON_BINS plot_tokamaker_psi.py ${OFT_PYTHON_BINS} CACHE INTERNAL "package_python_binaries" )