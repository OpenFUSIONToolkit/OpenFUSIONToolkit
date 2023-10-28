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
