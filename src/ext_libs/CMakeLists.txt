set( QUADPACK_OBJS
  dqagse.f
  dqelg.f
  dqk21.f
  dqpsrt.f
)

add_library( oftext
  dlsode.f
  lmdif.f
  bvls.f90
  ${QUADPACK_OBJS}
)

# Local copy of triangle
if( OFT_BUILD_PYTHON )
  add_library(oft_triangle SHARED triangle.c)
  target_compile_definitions(oft_triangle PRIVATE -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -DREAL=double -DVOID=void)
  target_compile_options(oft_triangle PRIVATE -O1)
  install(
    TARGETS oft_triangle
    DESTINATION bin
    COMPONENT app
  )
endif()