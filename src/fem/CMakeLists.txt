set( FEM_OBJS
  fem_base.F90
  fem_composite.F90
  fem_utils.F90
  lagrange_basis.F90
  lagrange_operators.F90
  blag_operators.F90
  h1_basis.F90
  h1_operators.F90
  hcurl_basis.F90
  hcurl_operators.F90
  hcurl_grad_operators.F90
)

add_library( oftfem ${FEM_OBJS} )