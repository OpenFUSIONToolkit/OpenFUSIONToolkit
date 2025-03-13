'''! Python interface for Marklin force-free ideal MHD equilibrium functionality

@authors Chris Hansen
@date September 2023
@ingroup doxy_oft_python
'''
import numpy
from .._interface import *

## @cond
# marklin_setup(marklin_ptr,mesh_ptr,order,minlev,error_str)
marklin_setup = ctypes_subroutine(oftpy_lib.marklin_setup,
    [c_void_ptr_ptr, c_void_p, c_int, c_int, c_char_p])

# marklin_compute(marklin_ptr,nmodes,save_rst,eig_vals,error_str)
marklin_compute = ctypes_subroutine(oftpy_lib.marklin_compute,
    [c_void_p, c_int, c_bool, ctypes_numpy_array(numpy.float64,1), c_char_p])

# (basepath,error_str
marklin_setup_io = ctypes_subroutine(oftpy_lib.marklin_setup_io,
    [c_void_p, c_char_p, c_char_p])

# Marklin setup function (mesh and such) order,minlev,nh,hcpc,hcpv,save_rst,error_str
marklin_compute_vac = ctypes_subroutine(oftpy_lib.marklin_compute_vac,
    [c_int, c_int, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool, c_char_p])

# int_obj,int_type,k_perp,error_str
marklin_compute_pardiff = ctypes_subroutine(oftpy_lib.marklin_compute_pardiff,
    [c_void_p, c_int, c_double, c_char_p])

#
marklin_save_visit = ctypes_subroutine(oftpy_lib.marklin_save_visit,
    [c_void_p, c_void_p, c_int, c_char_p, c_char_p])

# hmode_facs,int_obj,error_str
marklin_get_aint = ctypes_subroutine(oftpy_lib.marklin_get_aint,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_void_ptr_ptr, c_bool, c_char_p])

# hmode_facs,vac_facs,int_obj,error_str
marklin_get_bint = ctypes_subroutine(oftpy_lib.marklin_get_bint,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_void_ptr_ptr, c_char_p])

#
marklin_apply_int = ctypes_subroutine(oftpy_lib.marklin_apply_int,
    [c_void_p, c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_int_ptr ,ctypes_numpy_array(numpy.float64,1)])
## @endcond