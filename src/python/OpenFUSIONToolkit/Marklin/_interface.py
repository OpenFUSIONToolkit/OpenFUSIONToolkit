'''! Python interface for Marklin force-free ideal MHD equilibrium functionality

@authors Chris Hansen
@date September 2023
@ingroup doxy_oft_python
'''
import numpy
from ..util import *

## @cond
# Marklin setup function (mesh and such) order,nmodes,minlev,error_str
marklin_compute = ctypes_subroutine(oftpy_lib.marklin_compute,
    [c_int, c_int, c_int, c_bool, ctypes_numpy_array(numpy.float64,1), c_char_p])

# (basepath,error_str
marklin_setup_io = ctypes_subroutine(oftpy_lib.marklin_setup_io,
    [c_char_p, c_char_p])

#
marklin_save_visit = ctypes_subroutine(oftpy_lib.marklin_save_visit,
    [c_void_p, c_int, c_char_p, c_char_p])

#
marklin_get_aint = ctypes_subroutine(oftpy_lib.marklin_get_aint,
    [c_int, c_void_ptr_ptr, c_bool, c_char_p])

#
marklin_get_bint = ctypes_subroutine(oftpy_lib.marklin_get_bint,
    [c_int, c_void_ptr_ptr, c_char_p])

#
marklin_apply_int = ctypes_subroutine(oftpy_lib.marklin_apply_int,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_int_ptr ,ctypes_numpy_array(numpy.float64,1)])
## @endcond