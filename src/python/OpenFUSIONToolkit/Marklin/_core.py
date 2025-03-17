'''! Python interface for Marklin force-free ideal MHD equilibrium functionality

@authors Chris Hansen
@date September 2023
@ingroup doxy_oft_python
'''
import warnings
import numpy
from ._interface import *
from ..io import build_XDMF


class Marklin_field_interpolator():
    '''! Interpolation class for force-free eigenstate vector fields'''
    def __init__(self,marklin_obj,int_ptr,int_type,dim,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param marklin_obj Marklin instance for interpolator
        @param int_ptr Address of FORTRAN interpolation class
        @param int_type Interpolation type (1: vector potential; 2: magnetic field)
        @param dim Dimension of vector field
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        self.cell = c_int(-1)
        self.int_type = int_type
        self.dim = dim
        self.val = numpy.zeros((self.dim,), dtype=numpy.float64)
        self.marklin_obj = marklin_obj
        self._int_ptr = int_ptr
        self.fbary_tol = fbary_tol
    
    def __del__(self):
        '''Destroy underlying interpolation object'''
        pt_eval = numpy.zeros((3,), dtype=numpy.float64)
        marklin_apply_int(self.marklin_obj._marklin_ptr,self._int_ptr,-self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.val)

    def eval(self,pt):
        '''! Evaluate field at a given location

        @param pt Location for evaluation [3]
        @result Field at evaluation point [self.dim]
        '''
        marklin_apply_int(self.marklin_obj._marklin_ptr,self._int_ptr,self.int_type,pt,self.fbary_tol,ctypes.byref(self.cell),self.val)
        return self.val


class Marklin():
    '''! Marklin force-free equilibrium solver class'''
    def __init__(self,OFT_env):
        '''! Initialize Marklin object

        @param OFT_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        '''
        # Create OFT execution environment
        self._oft_env = OFT_env
        ## Internal Marklin solver object
        self._marklin_ptr = c_void_p()
        ## Internal mesh object
        self._mesh_ptr = c_void_p()
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
        self.nm = 0
        ## Needs docs
        self.nh = 0
        ## Needs docs
        self.hcpc = None
        ## Needs docs
        self.hcpv = None
        ## Eigenvalues
        self.eig_vals = None
        ## I/O basepath for plotting/XDMF output
        self._io_basepath = "."
    
    def __del__(self):
        if not self._marklin_ptr:
            return # Nothing to do
        error_string = self._oft_env.get_c_errorbuff()
        marklin_destroy(self._marklin_ptr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def setup_mesh(self,r=None,lc=None,reg=None,mesh_file=None,grid_order=1):
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
            self._oft_env.oft_in_groups['mesh_options'] = {
                'cad_type': "0",
                'grid_order': '{0:d}'.format(grid_order)
            }
            self._oft_env.oft_in_groups['native_mesh_options'] = {'filename': '"{0}"'.format(mesh_file)}
            self._oft_env.update_oft_in()
            oft_setup_vmesh(ndim,ndim,rfake,ndim,ndim,lcfake,regfake,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
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
            oft_setup_vmesh(ndim,np,r,npc,nc,lc+1,reg,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
    
    def setup(self,order=2,minlev=-1):
        '''! Needs docs
        '''
        error_string = self._oft_env.get_c_errorbuff()
        marklin_setup(ctypes.byref(self._marklin_ptr),self._mesh_ptr,order,minlev,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def setup_io(self,basepath=None):
        '''! Setup XDMF+HDF5 I/O for 3D visualization

        @param basepath Path to root directory to use for I/O
        '''
        if basepath is None:
            basepath_c = self._oft_env.path2c('')
            self._io_basepath = "."
        else:
            if basepath[-1] != '/':
                basepath += '/'
            self._io_basepath = basepath[:-1]
            basepath_c = self._oft_env.path2c(basepath)
        error_string = self._oft_env.get_c_errorbuff()
        marklin_setup_io(self._marklin_ptr,basepath_c,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value.decode())
    
    def build_XDMF(self,repeat_static=False,pretty=False):
        '''! Build XDMF plot metadata files for model

        @param repeat_static Repeat static fields (0-th timestep) in all timesteps?
        @param pretty Use pretty printing (indentation) in XDMF files?
        '''
        return build_XDMF(path=self._io_basepath,repeat_static=repeat_static,pretty=pretty)

    def compute_eig(self,nmodes=1,cache_file=None,save_rst=False):
        r'''! Compute force-free eigenmodes

        @param nmodes Number of eigenmodes to compute
        @param cache_file Path to cache file to store/load modes
        @param save_rst Save restart files? (deprecated)
        '''
        if self.nm != 0:
            raise ValueError('Eigenstates already computed')
        if save_rst:
            warnings.warn("Argument `save_rst` is deprecated, use `cache_file` instead", DeprecationWarning)
            cache_file = 'oft_Marklin.rst'
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        #
        eig_vals = numpy.zeros((nmodes,),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        marklin_compute_eig(self._marklin_ptr,nmodes,eig_vals,cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self.nm = nmodes
        self.eig_vals = eig_vals
    
    def compute_vac(self,nh,hcpc,hcpv,cache_file=None):
        r'''! Compute vacuum field with specified fluxes through jump planes

        @param nh Number of jump planes
        @param hcpc Plane specification points
        @param hcpv Plane specification vectors
        @param cache_file Path to cache file to store/load modes
        '''
        if hcpc.shape[0] != nh:
            raise ValueError('Inconsistent sizes for "hcpc[0]" != {0}'.format(nh))
        if hcpc.shape[1] != 3:
            raise ValueError('Inconsistent sizes for "hcpc[0]" != {0}'.format(3))
        if hcpv.shape[0] != nh:
            raise ValueError('Inconsistent sizes for "hcpv[0]" != {0}'.format(nh))
        if hcpv.shape[1] != 3:
            raise ValueError('Inconsistent sizes for "hcpv[0]" != {0}'.format(3))
        if cache_file is None:
            cache_string = self._oft_env.path2c("")
        else:
            cache_string = self._oft_env.path2c(cache_file)
        error_string = self._oft_env.get_c_errorbuff()
        marklin_compute_vac(self._marklin_ptr,nh,hcpc,hcpv,cache_string,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self.nh = nh
        self.hcpc = hcpc
        self.hcpv = hcpv
    
    def compute_par_diff(self,interpolator,k_perp):
        r'''! Compute parallel diffusion with specified vector field and diffusivity

        @param interpolator Interpolator defining vector field
        @param k_perp Perpendicular diffusion value (k_par = 1.0)
        '''
        cstring = c_char_p(b""*200)
        marklin_compute_pardiff(self._marklin_ptr,interpolator._int_ptr,interpolator.int_type,k_perp,cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)

    def save_field(self,field,tag):
        '''! Save field to XDMF files for VisIt/Paraview
        
        @param field Field interpolation object
        @param tag Name for field in XDMF files
        '''
        #
        error_string = self._oft_env.get_c_errorbuff()
        ctag = self._oft_env.path2c(tag)
        marklin_save_visit(self._marklin_ptr,field._int_ptr,field.int_type,ctag,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def get_ainterp(self,hmode_facs,bn_gauge=False):
        r'''! Create field interpolator for vector potential

        @param imode Index of eigenstate
        @param bn_gauge Use B-field gauge (A_t = 0 @ wall)?
        @result Field interpolation object
        '''
        if hmode_facs.shape[0] > self.nm:
            raise ValueError("Requested mode number exceeds number of available modes")
        interpolation_ptr = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        marklin_get_aint(self._marklin_ptr,hmode_facs,ctypes.byref(interpolation_ptr),bn_gauge,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Marklin_field_interpolator(self,interpolation_ptr,1,3)

    def get_binterp(self,hmode_facs=None,vac_facs=None):
        r'''! Create field interpolator for magnetic field

        @param imode Index of eigenstate
        @result Field interpolation object
        '''
        if hmode_facs is None:
            if vac_facs is None:
                raise ValueError('"hmode_facs" or "vac_facs" must be specified.')
            hmode_facs = numpy.zeros((self.nm,), dtype=numpy.float64)
        else:
            if hmode_facs.shape[0] > self.nm:
                raise ValueError("Requested mode number exceeds number of available modes")
        if vac_facs is None:
            vac_facs = numpy.zeros((self.nh,), dtype=numpy.float64)
        interpolation_ptr = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        marklin_get_bint(self._marklin_ptr,hmode_facs,vac_facs,ctypes.byref(interpolation_ptr),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Marklin_field_interpolator(self,interpolation_ptr,2,3)