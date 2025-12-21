set( OBJS
  mhd_utils.F90
  oft_scalar_inits.F90
  oft_vector_inits.F90
  diagnostic.F90
  taylor.F90
  xmhd.F90
  xmhd_lag.F90
  tracing.F90
  gs_eq.F90
  tracing_2d.F90
  axi_green.F90
  grad_shaf.F90
  grad_shaf_profiles.F90
  grad_shaf_prof_phys.F90
  grad_shaf_util.F90
  grad_shaf_fit.F90
  grad_shaf_td.F90
  thin_wall.F90
  thin_wall_solvers.F90
  thin_wall_hodlr.F90
)

add_library( oftphysics ${OBJS} )