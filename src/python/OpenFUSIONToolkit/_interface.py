'''! Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''

## @file util.py
#
# Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces
import os
import platform
import ctypes
from ctypes import c_bool, c_int, c_double, c_char_p, c_void_p, create_string_buffer
import numpy
from numpy import float64, int32

# Helper datatypes
## ctypes logical (bool) pointer alias 
c_bool_ptr = ctypes.POINTER(c_bool)
## ctypes logical (bool) double pointer alias
c_bool_ptr_ptr = ctypes.POINTER(c_bool_ptr)
## ctypes 32-bit integer (int) pointer alias
c_int_ptr = ctypes.POINTER(c_int)
## ctypes 32-bit integer (int) double pointer alias
c_int_ptr_ptr = ctypes.POINTER(c_int_ptr)
## ctypes 64-bit floating point (double) pointer type
c_double_ptr = ctypes.POINTER(c_double)
## ctypes 64-bit floating point (double) double pointer
c_double_ptr_ptr = ctypes.POINTER(c_double_ptr)
## ctypes void double pointer
c_void_ptr_ptr = ctypes.POINTER(c_void_p)
## ctypes struct alias
c_struct = ctypes.Structure

def ctypes_numpy_array(type,ndim):
    '''! Create a ctypes argument object for a given numpy array type and dimension

    @param type NumPy c-compatible datatype (https://numpy.org/doc/stable/user/basics.types.html)
    @param ndim Number of dimensions in array
    '''
    return numpy.ctypeslib.ndpointer(dtype=type,ndim=ndim,flags='C_CONTIGUOUS')

def ctypes_subroutine(function, argtypes=None, restype=None):
    '''! Create a ctypes object for a FORTRAN subroutine (no return value)

    @param function Function object from ctypes library
    @param argtypes List of types for each argument (optional)
    '''
    tmp_fun = function
    tmp_fun.restype = restype
    if argtypes is not None:
        tmp_fun.argtypes = argtypes
    return tmp_fun

# Common parameters
## Vacuum magnetic permeability
mu0 = numpy.pi*4.E-7
## Electron charge
eC = 1.60217663e-19

## @cond
root_path = os.path.realpath(os.path.dirname(__file__))
if platform.system() == 'Linux':
    lib_suffix = '.so'
elif platform.system() == 'Darwin':
    lib_suffix = '.dylib'
else:
    raise SystemError('Unsupported platform type')
oftpy_lib = ctypes.CDLL(os.path.join(root_path,'..','..','bin','liboftpy'+lib_suffix))
oft_triangle_lib = ctypes.CDLL(os.path.join(root_path,'..','..','bin','liboft_triangle'+lib_suffix))

# Global init function
oft_init = ctypes_subroutine(oftpy_lib.oftpy_init,
    [c_int, c_char_p, ctypes_numpy_array(int32,1)])

# oftpy_set_debug(debug_level)
oftpy_set_debug = ctypes_subroutine(oftpy_lib.oftpy_set_debug,
    [c_int])

# oftpy_set_nthreads(nthreads)
oftpy_set_nthreads = ctypes_subroutine(oftpy_lib.oftpy_set_nthreads,
    [c_int])

# oftpy_load_xml(xml_file,oft_node_ptr)
oftpy_load_xml = ctypes_subroutine(oftpy_lib.oftpy_load_xml,
    [c_char_p, c_void_ptr_ptr])

# Set mesh in memory: (ndim,np,r_loc,npc,nc,lc_loc,reg_loc)
oft_setup_smesh = ctypes_subroutine(oftpy_lib.oft_setup_smesh,
    [c_int,c_int, ctypes_numpy_array(float64,2) ,c_int, c_int, ctypes_numpy_array(int32,2), ctypes_numpy_array(int32,1), c_int_ptr])

# Set mesh in memory: (ndim,np,r_loc,npc,nc,lc_loc,reg_loc)
oft_setup_vmesh = ctypes_subroutine(oftpy_lib.oft_setup_vmesh,
    [c_int,c_int, ctypes_numpy_array(float64,2) ,c_int, c_int, ctypes_numpy_array(int32,2), ctypes_numpy_array(int32,1), c_int_ptr])

# Dump coverage information if needed
oftpy_dump_cov = ctypes_subroutine(oftpy_lib.dump_cov)
## @endcond