'''! Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''

## @file util.py
#
# Helper interfaces for Open FUSION Toolkit (OFT) Python interfaces
import sys
import os
import platform
import ctypes
from ctypes import c_bool, c_int, c_double, c_char_p, c_void_p, create_string_buffer
import subprocess
import numpy
from numpy import float64, int32
import h5py

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


def build_XDMF(path='.',repeat_static=False,pretty=False):
    '''! Build XDMF plot metadata files 

    @param path Folder to build XDMF files in (must include "dump.dat" file)
    @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
    @param pretty Use pretty printing (indentation) in XDMF files?
    '''
    cmd = [
        "{0}".format(sys.executable),
        "{0}".format(os.path.join(os.path.dirname(__file__),'..','build_xdmf.py'))
    ]
    if repeat_static:
        cmd.append("--repeat_static")
    if pretty:
        cmd.append("--pretty")
    subprocess.run(cmd,cwd=path)


def write_native_mesh(filename, r, lc, reg, nodesets=[], sidesets=[], ho_info=None, periodic_info=None):
    r'''Create a native HDF5 mesh file for OFT from the given mesh information

    @param filename Filename for mesh file
    @param r Points list [np,3]
    @param lc Cell list [nc,3] (1-based)
    @param reg Region list [nc]
    @param nodesets List of node sets
    @param sidesets List of side sets
    @param ho_info High-order grid information
    @param periodic_info Information for mesh periodicity
    '''
    print()
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5_file:
        h5_file.create_dataset('mesh/R', data=r, dtype='f8')
        h5_file.create_dataset('mesh/LC', data=lc, dtype='i4')
        h5_file.create_dataset('mesh/REG', data=reg, dtype='i4')
        if len(nodesets) > 0:
            h5_file.create_dataset('mesh/NUM_NODESETS', data=[len(nodesets),], dtype='i4')
            for i, node_set in enumerate(nodesets):
                h5_file.create_dataset('mesh/NODESET{0:04d}'.format(i+1), data=node_set, dtype='i4')
        if len(sidesets) > 0:
            h5_file.create_dataset('mesh/NUM_SIDESETS', data=[len(sidesets),], dtype='i4')
            for i, side_set in enumerate(sidesets):
                h5_file.create_dataset('mesh/SIDESET{0:04d}'.format(i+1), data=side_set, dtype='i4')
        if ho_info is not None:
            h5_file.create_dataset('mesh/ho_info/R', data=ho_info[0], dtype='f8')
            h5_file.create_dataset('mesh/ho_info/LE', data=ho_info[1], dtype='i4')
            if ho_info[2] is not None:
                h5_file.create_dataset('mesh/ho_info/LF', data=ho_info[2], dtype='i4')
        if periodic_info is not None:
            h5_file.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')