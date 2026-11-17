#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Fortran interface definitions for ThickCurr

@authors Chris Hansen
@date November 2026
@ingroup doxy_oft_python
'''
import numpy
from .._interface import *

## @cond
# thickcurr_f_setup(tc_ptr,mesh_ptr,order,minlev,neta,eta,xml_ptr,error_str)
thickcurr_f_setup = ctypes_subroutine(oftpy_lib.thickcurr_f_setup,
    [c_void_ptr_ptr, c_void_p, c_int, c_int, c_int, ctypes_numpy_array(numpy.float64,1), c_void_p, c_char_p])

# thickcurr_f_destroy(tc_ptr,error_str)
thickcurr_f_destroy = ctypes_subroutine(oftpy_lib.thickcurr_f_destroy,
    [c_void_p, c_char_p])

# thickcurr_f_set_options(tc_ptr,dt,theta,lin_tol,nsteps,rst_freq,error_str)
thickcurr_f_set_options = ctypes_subroutine(oftpy_lib.thickcurr_f_set_options,
    [c_void_p, c_double, c_double, c_double, c_int, c_int, c_char_p])

# thickcurr_f_get_qpts(tc_ptr,npts,pts,error_str)
thickcurr_f_get_qpts = ctypes_subroutine(oftpy_lib.thickcurr_f_get_qpts,
    [c_void_p, c_int_ptr, c_void_p, c_char_p])

# thickcurr_f_set_field(tc_ptr,vals,zero_bdry,energy,error_str)
thickcurr_f_set_field = ctypes_subroutine(oftpy_lib.thickcurr_f_set_field,
    [c_void_p, c_void_p, c_bool, c_double_ptr, c_char_p])

# thickcurr_f_get_energy(tc_ptr,energy,error_str)
thickcurr_f_get_energy = ctypes_subroutine(oftpy_lib.thickcurr_f_get_energy,
    [c_void_p, c_double_ptr, c_char_p])

# thickcurr_f_eigs(tc_ptr,nev,tau,shift,save_plot,error_str)
thickcurr_f_eigs = ctypes_subroutine(oftpy_lib.thickcurr_f_eigs,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_bool, c_char_p])

# thickcurr_f_run(tc_ptr,sensor_ptr,t_final,error_str)
thickcurr_f_run = ctypes_subroutine(oftpy_lib.thickcurr_f_run,
    [c_void_p, c_void_p, c_double_ptr, c_char_p])

# thickcurr_f_plot(tc_ptr,sensor_ptr,rst_start,rst_end,error_str)
thickcurr_f_plot = ctypes_subroutine(oftpy_lib.thickcurr_f_plot,
    [c_void_p, c_void_p, c_int, c_int, c_char_p])

# thickcurr_f_setup_sensors(tc_ptr,nsensors,names,npa,ptsa,wa,idxa,npb,ptsb,wb,idxb,filename,sensor_ptr,error_str)
thickcurr_f_setup_sensors = ctypes_subroutine(oftpy_lib.thickcurr_f_setup_sensors,
    [c_void_p, c_int, c_char_p, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2),
     ctypes_numpy_array(numpy.int32,1), c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2),
     ctypes_numpy_array(numpy.int32,1), c_char_p, c_void_ptr_ptr, c_char_p])

# thickcurr_f_eval_sensors(tc_ptr,sensor_ptr,vals,error_str)
thickcurr_f_eval_sensors = ctypes_subroutine(oftpy_lib.thickcurr_f_eval_sensors,
    [c_void_p, c_void_p, ctypes_numpy_array(numpy.float64,1), c_char_p])

# thickcurr_f_destroy_sensors(sensor_ptr)
thickcurr_f_destroy_sensors = ctypes_subroutine(oftpy_lib.thickcurr_f_destroy_sensors,
    [c_void_p])
## @endcond
