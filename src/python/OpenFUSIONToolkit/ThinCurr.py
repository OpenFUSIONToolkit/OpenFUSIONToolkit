'''! Python interface for ThinCurr thin-wall eddy current functionality

@authors Chris Hansen
@date March 2024
@ingroup doxy_oft_python
'''

##@file ThinCurr.py
#
# Python interface for ThinCurr thin-wall eddy current functionality
import numpy
from .util import *

## @cond
# ThinCurr setup function (load mesh and setup model) (mesh_file,np,r_loc,nc,lc_loc,reg_loc,pmap_loc,tw_ptr,size,error_str)
thincurr_setup = ctypes_subroutine(oftpy_lib.thincurr_setup,
    [c_char_p, c_int, ctypes_numpy_array(float64,2), c_int, ctypes_numpy_array(int32,2), ctypes_numpy_array(int32,1), 
     ctypes_numpy_array(int32,1), c_void_p, ctypes_numpy_array(int32,1), c_char_p])

# ThinCurr setup plotting (tw_ptr,basepath,save_debug,error_str
thincurr_setup_io = ctypes_subroutine(oftpy_lib.thincurr_setup_io,
    [c_void_p, c_char_p, c_bool, c_char_p])

# ThinCurr save current potential field (tw_ptr,vals,fieldname)
thincurr_save_field = ctypes_subroutine(oftpy_lib.thincurr_save_field,
    [c_void_p, ctypes_numpy_array(float64,1), c_char_p])

# ThinCurr save scalar field (tw_ptr,vals,fieldname)
thincurr_save_scalar = ctypes_subroutine(oftpy_lib.thincurr_save_scalar,
    [c_void_p, ctypes_numpy_array(float64,1), c_char_p])

# Compute mutual coupling between models thincurr_cross_coupling(tw_ptr1,tw_ptr2,Mmat,error_str)
thincurr_cross_coupling = ctypes_subroutine(oftpy_lib.thincurr_cross_coupling,
    [c_void_p, c_void_p, ctypes_numpy_array(float64,2), c_char_p, c_char_p])

# Compute model inductance matrix thincurr_curr_Lmat(tw_ptr,Lmat,error_str)
thincurr_curr_Lmat = ctypes_subroutine(oftpy_lib.thincurr_Lmat,
    [c_void_p, c_double_ptr_ptr, c_char_p, c_char_p])

# Compute model resistivity matrix thincurr_curr_Rmat(tw_ptr,Rmat,error_str)
thincurr_curr_Rmat = ctypes_subroutine(oftpy_lib.thincurr_Rmat,
    [c_void_p, ctypes_numpy_array(float64,2), c_char_p])

# Compute current regularization matrix thincurr_curr_regmat(tw_ptr,Rmat,error_str)
thincurr_curr_regmat = ctypes_subroutine(oftpy_lib.thincurr_curr_regmat,
    [c_void_p, ctypes_numpy_array(float64,2), c_char_p])

## @endcond

oft_in_template = """&runtime_options
 debug={DEBUG_LEVEL}
/

&mesh_options
 meshname='none'
 {MESH_TYPE}
/

{MESH_DEF}
"""

class ThinCurr():
    '''! ThinCurr thin-wall eddy current model class'''
    def __init__(self,debug_level=0,nthreads=2):
        '''! Initialize ThinCurr object

        @param debug_level Level of debug printing (0-3)
        '''
        ## Input file settings
        self._psin_dict = {'DEBUG_LEVEL': debug_level, 'MESH_TYPE': '', 'MESH_DEF': ''}
        self._update_psin()
        oft_init(c_int(nthreads))
        ## Thin-wall model object
        self.tw_obj = c_void_p()
        ## Number of regions in mesh
        self.nregs = -1
        ## Number of points in mesh
        self.np = -1
        ## Number of edges in mesh
        self.ne = -1
        ## Number of cells in mesh
        self.nc = -1
        ## Number of vertices that are active in the thin-wall model
        self.np_active = -1
        ## Number of hole elements in the thin-wall model
        self.nholes = -1
        ## Number of V-coil elements in the thin-wall model
        self.n_vcoils = -1
        ## Total number of active DOFs in the thin-wall model
        self.nelems = -1
        ## Number of V-coil elements in the thin-wall model
        self.n_icoils = -1
        ## Model inductance matrix
        self.Lmat = None
        ## Model resistance matrix
        self.Rmat = None
        ## Mesh vertices [np,3] (last column should be all zeros)
        self.r = None
        ## Mesh triangles [nc,3]
        self.lc = None
        ## Mesh regions [nc]
        self.reg = None

    def _update_psin(self):
        '''! Update input file (`oftpyin`) with current settings'''
        with open('oftpyin','w+') as fid:
            fid.write(oft_in_template.format(**self._psin_dict))

    def setup_model(self,r=None,lc=None,reg=None,mesh_file=None,pmap=None):
        '''! Needs Docs
        '''
        if self.nregs != -1:
            raise ValueError('Mesh already setup, delete or create new instance for new model')
        nregs = c_int()
        if mesh_file is not None:
            idummy = c_int(-1)
            rfake = numpy.ones((1,1),dtype=numpy.float64)
            lcfake = numpy.ones((1,1),dtype=numpy.int32)
            regfake = numpy.ones((1,),dtype=numpy.int32)
            if pmap is None:
                pmap = -numpy.ones((1,),dtype=numpy.int32)
            else:
                pmap = numpy.ascontiguousarray(pmap, dtype=numpy.int32)
            sizes = numpy.zeros((8,),dtype=numpy.int32)
            filename = c_char_p(mesh_file.encode())
            cstring = c_char_p(b""*200)
            thincurr_setup(filename,idummy,rfake,idummy,lcfake,regfake,pmap,ctypes.byref(self.tw_obj),sizes,cstring)
            if cstring.value != b'':
                raise Exception(cstring.value)
            self.np = sizes[0]
            self.ne = sizes[1]
            self.nc = sizes[2]
            self.np_active = sizes[3]
            self.nholes = sizes[4]
            self.n_vcoils = sizes[5]
            self.nelems = sizes[6]
            self.n_icoils = sizes[7]
        elif r is not None:
            raise ValueError('Specifying mesh values not yet supported')
            # r = numpy.ascontiguousarray(r, dtype=numpy.float64)
            # lc = numpy.ascontiguousarray(lc, dtype=numpy.int32)
            # np = c_int(r.shape[0])
            # nc = c_int(lc.shape[0])
            # if reg is None:
            #     reg = numpy.ones((nc.value,),dtype=numpy.int32)
            # else:
            #     reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            # if pmap is None:
            #     pmap = -numpy.ones((1,),dtype=numpy.int32)
            # else:
            #     pmap = numpy.ascontiguousarray(pmap, dtype=numpy.float64)
            # sizes = numpy.zeros((8,),dtype=numpy.int32)
            # filename = c_char_p(b"")
            # cstring = c_char_p(b""*200)
            # thincurr_setup(filename,np,r,nc,lc,reg,ctypes.byref(self.tw_obj),sizes,cstring)
            # if cstring.value != b'':
            #     raise Exception(cstring.value)
            # self.np = sizes[0]
            # self.ne = sizes[1]
            # self.nc = sizes[2]
            # self.np_active = sizes[3]
            # self.nholes = sizes[4]
            # self.n_vcoils = sizes[5]
            # self.nelems = sizes[6]
            # self.n_icoils = sizes[7]
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
    
    def setup_io(self,basepath=None,save_debug=False):
        '''! Needs Docs
        '''
        # tw_ptr,basepath,save_debug,error_str
        if basepath is None:
            basepath = c_char_p(b'')
        else:
            if basepath[-1] != '/':
                raise ValueError('"basepath" must end in "/"')
            basepath = c_char_p(basepath.encode())
        cstring = c_char_p(b""*200)
        thincurr_setup_io(self.tw_obj,basepath,c_bool(save_debug),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
    
    def save_current(self,potential,tag):
        '''! Needs Docs
        '''
        # if potential.shape[0] != self.np:
        #     raise ValueError('Incorrect shape of "psi", should be [np]')
        potential = numpy.ascontiguousarray(potential, dtype=numpy.float64)
        cstring = c_char_p(tag.encode())
        thincurr_save_field(self.tw_obj,potential,cstring)
    
    def save_scalar(self,field,tag):
        '''! Needs Docs
        '''
        # if field.shape[0] != self.np:
        #     raise ValueError('Incorrect shape of "psi", should be [np]')
        field = numpy.ascontiguousarray(field, dtype=numpy.float64)
        ctag = c_char_p(tag.encode())
        thincurr_save_scalar(self.tw_obj,field,ctag)
    
    def compute_Lmat(self,cache_file=None):
        '''! Needs Docs
        '''
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        Lmat_loc = c_double_ptr()
        error_string = c_char_p(b""*200)
        thincurr_curr_Lmat(self.tw_obj,ctypes.byref(Lmat_loc),cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self.Lmat = numpy.ctypeslib.as_array(Lmat_loc,shape=(self.nelems,self.nelems))
    
    def compute_Rmat(self):
        '''! Needs Docs
        '''
        self.Rmat = numpy.zeros((self.nelems,self.nelems), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        thincurr_curr_Rmat(self.tw_obj,self.Rmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def cross_coupling(self,model2,cache_file=None):
        '''! Needs Docs
        '''
        Mmat = numpy.zeros((self.nelems,model2.nelems), dtype=numpy.float64)
        if cache_file is None:
            cache_string = c_char_p(b"")
        else:
            cache_string = c_char_p(cache_file.encode())
        error_string = c_char_p(b""*200)
        thincurr_cross_coupling(self.tw_obj,model2.tw_obj,Mmat,cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Mmat
    
    def get_regmat(self):
        '''! Needs Docs
        '''
        Rmat = numpy.zeros((self.nelems,3*self.nc), dtype=numpy.float64)
        error_string = c_char_p(b""*200)
        thincurr_curr_regmat(self.tw_obj,Rmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Rmat