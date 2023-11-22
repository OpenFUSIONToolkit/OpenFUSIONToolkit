install(PROGRAMS util.py DESTINATION python/OpenFUSIONToolkit COMPONENT app)
install(PROGRAMS TokaMaker.py DESTINATION python/OpenFUSIONToolkit COMPONENT app)
install(PROGRAMS Marklin.py DESTINATION python/OpenFUSIONToolkit COMPONENT app)

configure_file( util.py ${CMAKE_CURRENT_BINARY_DIR}/util.py COPYONLY)
configure_file( TokaMaker.py ${CMAKE_CURRENT_BINARY_DIR}/TokaMaker.py COPYONLY)
configure_file( Marklin.py ${CMAKE_CURRENT_BINARY_DIR}/Marklin.py COPYONLY)