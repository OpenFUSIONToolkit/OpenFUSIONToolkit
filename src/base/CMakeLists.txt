set( BASE_OBJS
  oft_local_c.c
  oft_sort_c.cxx
  spline.f
  oft_local.F90
  oft_base.F90
  oft_sort.F90
  oft_io.F90
  oft_stitching.F90
)

add_library( oftbase
  ${BASE_OBJS}
)
