set( OBJS
    oft_base_f.F90
    tokamaker_f.F90
    thincurr_f.F90
    marklin_f.F90
)

oft_add_shared(oftpy ${OBJS})
