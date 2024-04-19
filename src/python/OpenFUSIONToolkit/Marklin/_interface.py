'''! Python interface for Marklin force-free ideal MHD equilibrium functionality

@authors Chris Hansen
@date September 2023
@ingroup doxy_oft_python
'''
import numpy
from ..util import *

## @cond
# Marklin setup function (mesh and such) order,nmodes,minlev,error_str
marklin_compute_eig = ctypes_subroutine(oftpy_lib.marklin_compute_eigs,
    [c_int, c_int, c_int, c_bool, c_char_p])

# Marklin setup function (mesh and such) order,minlev,nh,hcpc,hcpv,save_rst,error_str
marklin_compute_vac = ctypes_subroutine(oftpy_lib.marklin_compute_vac,
    [c_int, c_int, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool, c_char_p])

# int_obj,int_type,k_perp,error_str
marklin_compute_pardiff = ctypes_subroutine(oftpy_lib.marklin_compute_pardiff,
    [c_void_p, c_int, c_double, c_char_p])

#
marklin_save_visit = ctypes_subroutine(oftpy_lib.marklin_save_visit,
    [c_char_p])

# hmode_facs,int_obj,error_str
marklin_get_aint = ctypes_subroutine(oftpy_lib.marklin_get_aint,
    [ctypes_numpy_array(numpy.float64,1), c_void_ptr_ptr, c_char_p])

# hmode_facs,vac_facs,int_obj,error_str
marklin_get_bint = ctypes_subroutine(oftpy_lib.marklin_get_bint,
    [ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_void_ptr_ptr, c_char_p])

#
marklin_apply_int = ctypes_subroutine(oftpy_lib.marklin_apply_int,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_int_ptr ,ctypes_numpy_array(numpy.float64,1)])
## @endcond