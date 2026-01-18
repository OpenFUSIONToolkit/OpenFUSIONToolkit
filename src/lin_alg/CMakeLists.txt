set( LA_OBJS
  lin_alg_base.F90
  native_la.F90
  petsc_la.F90
  deriv_matrices.F90
)
set( SOLVER_OBJS
  solver_base.F90
  native_solvers.F90
  petsc_solvers.F90
  superlu_bridge.c
  superlu_dist_bridge.c
  umfpack_bridge.c
  mkl_bridge.c
  oft_arpack.F90
  oft_lu.F90
)
set( MGMT_OBJS
  lin_alg_utils.F90
  solver_utils.F90
)

add_library( oftla
  ${LA_OBJS}
  ${SOLVER_OBJS}
  ${MGMT_OBJS}
)