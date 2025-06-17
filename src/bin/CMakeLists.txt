set( OFT_EXES
  oft_mesh_check.F90
  oft_trace.F90
  oft_poincare.F90
)

set( MARKLIN_EXES
  marklin_eigs.F90
)

set( TOKAMAKER_EXES
  tokamaker_gs.F90
  tokamaker_fit.F90
  tokamaker_wall.F90
)

set( THINCURR_EXES
  thincurr_eig.F90
  thincurr_fr.F90
  thincurr_td.F90
  thincurr_from_mode.F90
  thincurr_coupling.F90
)

set( BIN_SRCS ${OFT_EXES} ${MARKLIN_EXES} ${THINCURR_EXES} )
if( OFT_TOKAMAKER_LEGACY )
  set( BIN_SRCS ${BIN_SRCS} ${TOKAMAKER_EXES} )
endif()

foreach( file ${BIN_SRCS} )
  oft_add_exe( ${file} )
endforeach()

# Create symlinks for Python binaries
install(
    CODE
    "
    foreach(python_bin ${OFT_PYTHON_BINS})
      file(RELATIVE_PATH REL_PATH \${CMAKE_INSTALL_PREFIX}/bin \${CMAKE_INSTALL_PREFIX}/python/\${python_bin})
      file(CREATE_LINK \${REL_PATH} \${CMAKE_INSTALL_PREFIX}/bin/\${python_bin} SYMBOLIC)
    endforeach()
    "
    COMPONENT app
)

if(OFT_PACKAGE_BUILD)
  include(OFTPatchPackage)
endif()