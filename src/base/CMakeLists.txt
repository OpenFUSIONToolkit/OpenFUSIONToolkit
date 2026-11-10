set( BASE_OBJS
  oft_local_c.c
  oft_xml_c.c
  oft_sort_c.cxx
  cvode_bridge.c
  spline.f
  oft_local.F90
  oft_base.F90
  oft_sort.F90
  oft_io.F90
  oft_stitching.F90
  oft_xml.F90
  oft_cvode.F90
)

add_library( oftbase
  ${BASE_OBJS}
)
