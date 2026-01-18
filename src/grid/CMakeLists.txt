set( UTIL_OBJS
  bezier_cad.F90
  c_nurbs.cxx
  nurbs_cad.F90
  metis_bridge.c
  quadrature.F90
  lagrange_poly.F90
)
set( TM_OBJS
  mesh_type.F90
  mesh_local.F90
  mesh_local_util.F90
  mesh_global.F90
  mesh_global_util.F90
  trimesh_type.F90
  trimesh_tessellation.F90
  tetmesh_type.F90
  tetmesh_tessellation.F90
  hexmesh_type.F90
  quadmesh_type.F90
)
set( TM_CAD_OBJS
  mesh_native.F90
  mesh_cube.F90
  mesh_sphere.F90
  mesh_t3d.F90
  mesh_cubit.F90
  mesh_gmsh.F90
)
set( MG_OBJS
  multigrid.F90
  multigrid_build.F90
)

add_library( oftgrid
  ${BASE_OBJS}
  ${UTIL_OBJS}
  ${TM_OBJS}
  ${TM_CAD_OBJS}
  ${MG_OBJS}
)