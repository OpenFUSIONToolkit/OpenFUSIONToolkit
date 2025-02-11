'''! Python interface for Marklin force-free ideal MHD equilibrium functionality

@authors Chris Hansen
@date September 2023
@ingroup doxy_oft_python
'''
import numpy
from ._interface import *
from ..util import *

oft_in_template = """&runtime_options
 debug={DEBUG_LEVEL}
/

&mesh_options
 meshname='none'
 grid_order={GRID_ORDER}
 {MESH_TYPE}
/

{MESH_DEF}
"""

class Marklin_field_interpolator():
    '''! Interpolation class for force-free eigenstate vector fields'''
    def __init__(self,int_obj,int_type,dim,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param int_obj Address of FORTRAN interpolation class
        @param int_type Interpolation type (1: vector potential; 2: magnetic field)
        @param dim Dimension of vector field
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        self.cell = c_int(-1)
        self.int_type = int_type
        self.dim = dim
        self.val = numpy.zeros((self.dim,), dtype=numpy.float64)
        self.int_obj = int_obj
        self.fbary_tol = fbary_tol
    
    def __del__(self):
        '''Destroy underlying interpolation object'''
        pt_eval = numpy.zeros((3,), dtype=numpy.float64)
        marklin_apply_int(self.int_obj,-self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.val)

    def eval(self,pt):
        '''! Evaluate field at a given location

        @param pt Location for evaluation [3]
        @result Field at evaluation point [self.dim]
        '''
        marklin_apply_int(self.int_obj,self.int_type,pt,self.fbary_tol,ctypes.byref(self.cell),self.val)
        return self.val

class Marklin():
    '''! Marklin force-free equilibrium solver class'''
    def __init__(self,debug_level=0,nthreads=2,grid_order=1):
        '''! Initialize Marklin object

        @param debug_level Level of debug printing (0-3)
        '''
        ## Input file settings
        self._psin_dict = {'DEBUG_LEVEL': debug_level, 'MESH_TYPE': '', 'GRID_ORDER': grid_order, 'MESH_DEF': ''}
        self._update_psin()
        oft_init(c_int(nthreads))
        ## Number of regions in mesh
        self.nregs = -1
        ## Number of points in mesh
        self.np = -1
        ## Number of cells in mesh
        self.nc = -1
        ## Mesh vertices [np,3] (last column should be all zeros)
        self.r = None
        ## Mesh triangles [nc,3]
        self.lc = None
        ## Mesh regions [nc]
        self.reg = None
        ## Number of modes
        self._nm = -1
        ## Eigenvalues
        self.eig_vals = None
        ## I/O basepath for plotting/XDMF output
        self._io_basepath = "."

    def _update_psin(self):
        '''! Update input file (`oftpyin`) with current settings'''
        with open('oftpyin','w+') as fid:
            fid.write(oft_in_template.format(**self._psin_dict))

    def setup_mesh(self,r=None,lc=None,reg=None,mesh_file=None):
        '''! Setup mesh for Marklin force-free equilibrium calculations

        A mesh should be specified by passing "r", "lc", and optionally "reg" or using a "mesh_file".

        @param r Mesh point list [np,3]
        @param lc Mesh cell list [nc,4] (base one)
        @param reg Mesh region list [nc] (base one)
        @param mesh_file Filename containing mesh to load (native format only)
        '''
        if self.nregs != -1:
            raise ValueError('Mesh already setup, must call "reset" before loading new mesh')
        nregs = c_int()
        if mesh_file is not None:
            ndim = c_int(-1)
            rfake = numpy.ones((1,1),dtype=numpy.float64)
            lcfake = numpy.ones((1,1),dtype=numpy.int32)
            regfake = numpy.ones((1,),dtype=numpy.int32)
            self._psin_dict['MESH_TYPE'] = 'cad_type=0'
            self._psin_dict['MESH_DEF'] = """&native_mesh_options
 filename='{0}'
/""".format(mesh_file)
            self._update_psin()
            oft_setup_vmesh(ndim,ndim,rfake,ndim,ndim,lcfake,regfake,ctypes.byref(nregs))
        elif r is not None:
            r = numpy.ascontiguousarray(r, dtype=numpy.float64)
            lc = numpy.ascontiguousarray(lc, dtype=numpy.int32)
            ndim = c_int(r.shape[1])
            np = c_int(r.shape[0])
            npc = c_int(lc.shape[1])
            nc = c_int(lc.shape[0])
            if reg is None:
                reg = numpy.ones((nc.value,),dtype=numpy.int32)
            else:
                reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            oft_setup_vmesh(ndim,np,r,npc,nc,lc+1,reg,ctypes.byref(nregs))
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
    
    def setup_io(self,basepath=None):
        '''! Setup XDMF+HDF5 I/O for 3D visualization

        @param basepath Path to root directory to use for I/O
        '''
        if basepath is None:
            basepath_c = c_char_p(b'')
            self._io_basepath = "."
        else:
            if basepath[-1] != '/':
                basepath += '/'
            self._io_basepath = basepath[:-1]
            basepath_c = c_char_p(basepath.encode())
        cstring = create_string_buffer(b"",200)
        marklin_setup_io(basepath_c,cstring)
        if cstring.value != b'':
            raise Exception(cstring.value.decode())
    
    def build_XDMF(self,repeat_static=False,pretty=False):
        '''! Build XDMF plot metadata files for model

        @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
        @param pretty Use pretty printing (indentation) in XDMF files?
        '''
        build_XDMF(path=self._io_basepath,repeat_static=repeat_static,pretty=pretty)

    def compute(self,nmodes=1,order=2,minlev=-1,save_rst=True):
        r'''! Compute force-free eigenmodes

        @param nmodes Number of eigenmodes to compute
        @param order Order of FE representation
        @param minlev Minimum level for multigrid solve
        @param save_rst Save restart files? 
        '''
        if self._nm != -1:
            raise ValueError('Eigenstates already computed')
        #
        eig_vals = numpy.zeros((nmodes,),dtype=numpy.float64)
        cstring = create_string_buffer(b"",200)
        marklin_compute(order,nmodes,minlev,save_rst,eig_vals,cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        self._nm = nmodes
        self.eig_vals = eig_vals

    def save_field(self,field,tag):
        '''! Save field to XDMF files for VisIt/Paraview
        
        @param field Field interpolation object
        @param tag Name for field in XDMF files
        '''
        #
        cstring = create_string_buffer(b"",200)
        marklin_save_visit(field.int_obj,field.int_type,c_char_p(tag.encode()),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)

    def get_ainterp(self,imode,bn_gauge=False):
        r'''! Create field interpolator for vector potential

        @param imode Index of eigenstate
        @param bn_gauge Use B-field gauge (A_t = 0 @ wall)?
        @result Field interpolation object
        '''
        if imode > self._nm:
            raise ValueError("Requested mode number exceeds number of available modes")
        int_obj = c_void_p()
        cstring = create_string_buffer(b"",200)
        marklin_get_aint(imode,ctypes.byref(int_obj),bn_gauge,cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        return Marklin_field_interpolator(int_obj,1,3)

    def get_binterp(self,imode):
        r'''! Create field interpolator for magnetic field

        @param imode Index of eigenstate
        @result Field interpolation object
        '''
        if imode > self._nm:
            raise ValueError("Requested mode number exceeds number of available modes")
        int_obj = c_void_p()
        cstring = create_string_buffer(b"",200)
        marklin_get_bint(imode,ctypes.byref(int_obj),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        return Marklin_field_interpolator(int_obj,2,3)