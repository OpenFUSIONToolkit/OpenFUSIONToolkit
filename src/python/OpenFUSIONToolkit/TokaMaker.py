'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''

##@file oftpy.py
#
# Python interface for TokaMaker Grad-Shafranov functionality
import ctypes
import json
import numpy
from .util import *

# Settings STRUCT
class tokamaker_settings_struct(c_struct):
    r'''! TokaMaker settings structure

     - `pm` Print 'performance' information (eg. iteration count) during run?
     - `free_boundary` Perform free-boundary calculation?
     - `has_plasma` Include plasma effects in calculation, vacuum otherwise?
     - `limited_only` Do not search for X-points when determining LCFS?
     - `maxits` Maximum NL iteration count for G-S solver
     - `mode` Parallel current source formulation used (0 -> define \f$F'\f$, 1 -> define \f$F*F'\f$)
     - `urf` Under-relaxation factor for NL fixed-point iteration
     - `nl_tol` Convergence tolerance for NL solver
     - `rmin` Minimum magnetic axis major radius, used to catch 'lost' equilibria
     - `lim_zmax` Maximum vertical range for limiter points, can be used to exclude complex diverter regions
     - `limiter_file` File containing additional limiter points not included in mesh (default: 'none')
    '''
    _fields_ = [("pm", c_bool),
                ("free_boundary", c_bool),
                ("has_plasma", c_bool),
                ("limited_only", c_bool),
                ("maxits", c_int),
                ("mode", c_int),
                ("urf", c_double),
                ("nl_tol", c_double),
                ("rmin", c_double),
                ("lim_zmax", c_double),
                ("limiter_file", ctypes.c_char*80)]


## @cond
# Set mesh in memory: (ndim,np,r_loc,npc,nc,lc_loc,reg_loc)
tokamaker_alloc = ctypes_subroutine(oftpy_lib.tokamaker_alloc,
    [c_void_ptr_ptr])

tokamaker_setup_regions = ctypes_subroutine(oftpy_lib.tokamaker_setup_regions,
    [c_char_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(int32,1), ctypes_numpy_array(float64,2), c_int])

tokamaker_eval_green = ctypes_subroutine(oftpy_lib.tokamaker_eval_green,
    [c_int, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_double, c_double, ctypes_numpy_array(float64,1)])

# G-S setup function (mesh and such)
tokamaker_setup = ctypes_subroutine(oftpy_lib.tokamaker_setup,
    [c_int, c_bool, c_int_ptr, c_char_p])

# G-S setup function (mesh and such)
tokamaker_reset = ctypes_subroutine(oftpy_lib.tokamaker_reset,
    [c_char_p])

# G-S settings function
tokamaker_set_settings = ctypes_subroutine(oftpy_lib.tokamaker_set_settings,
    [ctypes.POINTER(tokamaker_settings_struct)])

# G-S init function (r0,z0,a,kappa,delta)
tokamaker_init_psi = ctypes_subroutine(oftpy_lib.tokamaker_init_psi,
    [c_double, c_double, c_double, c_double, c_double, c_int_ptr])

# G-S load flux functions (f_file,f_offset,p_file)
tokamaker_load_profiles = ctypes_subroutine(oftpy_lib.tokamaker_load_profiles,
    [c_char_p, c_double, c_char_p, c_char_p, c_char_p])

# G-S run function
tokamaker_run = ctypes_subroutine(oftpy_lib.tokamaker_run, 
    [c_bool, c_int_ptr])

# G-S info function
tokamaker_analyze = ctypes_subroutine(oftpy_lib.tokamaker_analyze)

# G-S time-dependent run function
tokamaker_setup_td = ctypes_subroutine(oftpy_lib.tokamaker_setup_td,
    [c_double, c_double, c_double, c_bool])

# G-S time-dependent run function
tokamaker_eig_td = ctypes_subroutine(oftpy_lib.tokamaker_eig_td,
    [c_double, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool, c_bool])

# G-S time-dependent run function
tokamaker_step_td = ctypes_subroutine(oftpy_lib.tokamaker_step_td,
    [c_double_ptr, c_double_ptr, c_int_ptr, c_int_ptr, c_int_ptr])

# G-S time-dependent run function
tokamaker_eig_wall = ctypes_subroutine(oftpy_lib.tokamaker_eig_wall,
    [c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool])

# G-S mesh retrieval
tokamaker_get_mesh = ctypes_subroutine(oftpy_lib.tokamaker_get_mesh,
    [c_int_ptr, c_double_ptr_ptr, c_int_ptr, c_int_ptr_ptr, c_int_ptr_ptr])

# G-S flux retrieval
tokamaker_get_psi = ctypes_subroutine(oftpy_lib.tokamaker_get_psi,
    [ctypes_numpy_array(numpy.float64,1), c_double_ptr, c_double_ptr])

#
tokamaker_set_psi = ctypes_subroutine(oftpy_lib.tokamaker_set_psi,
    [ctypes_numpy_array(numpy.float64,1)])

#
tokamaker_set_psi_dt = ctypes_subroutine(oftpy_lib.tokamaker_set_psi_dt,
    [ctypes_numpy_array(numpy.float64,1), c_double])

#
tokamaker_get_field_eval = ctypes_subroutine(oftpy_lib.tokamaker_get_field_eval,
    [c_int, c_void_ptr_ptr, c_char_p])

#
tokamaker_apply_field_eval = ctypes_subroutine(oftpy_lib.tokamaker_apply_field_eval,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_int_ptr, c_int, ctypes_numpy_array(numpy.float64,1)])

#
tokamaker_get_coil_currents = ctypes_subroutine(oftpy_lib.tokamaker_get_coil_currents,
    [ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1)])

#
tokamaker_get_coil_Lmat = ctypes_subroutine(oftpy_lib.tokamaker_get_coil_Lmat,
    [ctypes_numpy_array(numpy.float64,2)])

#
tokamaker_get_refs = ctypes_subroutine(oftpy_lib.tokamaker_get_refs, # (o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm)
    [c_double_ptr_ptr, c_double_ptr_ptr, c_double_ptr_ptr, c_bool_ptr_ptr, c_double_ptr_ptr, c_double_ptr_ptr,  c_double_ptr_ptr])

#
tokamaker_trace_surf = ctypes_subroutine(oftpy_lib.tokamaker_trace_surf,  # (psi_surf,points,npoints)
    [c_double, c_double_ptr_ptr, c_int_ptr])

#
tokamaker_get_q = ctypes_subroutine(oftpy_lib.tokamaker_get_q, # (npsi,psi_q,qvals,ravgs,dl,rbounds,zbounds)
    [c_int,ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,2),
     c_double_ptr, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2)])

#
tokamaker_sauter_fc = ctypes_subroutine(oftpy_lib.tokamaker_sauter_fc, # (npsi,psi_saut,fc,r_avgs,modb_avgs)
    [c_int,ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,2),
     ctypes_numpy_array(numpy.float64,2)])

#
tokamaker_get_globals = ctypes_subroutine(oftpy_lib.tokamaker_get_globals, # (Itor,centroid,vol,pvol,dflux,tflux,li)
    [c_double_ptr, ctypes_numpy_array(numpy.float64,1), c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr])

#
tokamaker_gs_calc_vloop = ctypes_subroutine(oftpy_lib.tokamaker_gs_calc_vloop, # (V_loop)
    [c_double_ptr])

#
tokamaker_get_profs = ctypes_subroutine(oftpy_lib.tokamaker_get_profs, # (npsi,psi_in,f,fp,p,pp)
    [c_int,ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), 
     ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1)])

#
tokamaker_set_coil_currents = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_currents,
    [ctypes_numpy_array(numpy.float64,1)])

# G-S set global targets (ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target)
tokamaker_set_targets = ctypes_subroutine(oftpy_lib.tokamaker_set_targets,
    [c_double, c_double, c_double, c_double, c_double, c_double])

# G-S set isoflux targets
tokamaker_set_isoflux = ctypes_subroutine(oftpy_lib.tokamaker_set_isoflux,
    [ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), c_int, c_double])

# G-S set saddle (X-point) targets
tokamaker_set_saddles = ctypes_subroutine(oftpy_lib.tokamaker_set_saddles,
    [ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), c_int])

# G-S set coil regularization matrix
tokamaker_set_coil_regmat = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_regmat,
    [c_int, ctypes_numpy_array(numpy.float64,2),ctypes_numpy_array(numpy.float64,1),ctypes_numpy_array(numpy.float64,1)])

# G-S set coil regularization matrix
tokamaker_set_coil_bounds = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_bounds,
    [ctypes_numpy_array(numpy.float64,2)])

tokamaker_set_coil_vsc = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_vsc,
    [ctypes_numpy_array(numpy.float64,1)])

tokamaker_get_vfixed = ctypes_subroutine(oftpy_lib.tokamaker_get_vfixed, #(npts,pts,fluxes)
    [c_int_ptr, c_double_ptr_ptr, c_double_ptr_ptr])

tokamaker_get_limiter = ctypes_subroutine(oftpy_lib.tokamaker_get_limiter, #(np,r_loc)
    [c_int_ptr,c_double_ptr_ptr])

tokamaker_save_eqdsk = ctypes_subroutine(oftpy_lib.tokamaker_save_eqdsk, #(filename,nr,nz,rbounds,zbounds,run_info,psi_pad)
    [c_char_p, c_int, c_int, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_char_p, c_double])
## @endcond

def tokamaker_default_settings():
    '''! Initialize settings object with default values

    @result tokamaker_settings_struct object
    '''
    settings = tokamaker_settings_struct()
    settings.pm = True
    settings.free_boundary = True
    settings.has_plasma = True
    settings.limited_only = False
    settings.maxits = 40
    settings.mode = 1
    settings.urf = 0.2
    settings.nl_tol = 1.E-6
    settings.rmin = 0.0
    settings.lim_zmax = 1.E99
    settings.limiter_file = 'none'.encode()
    return settings

def create_isoflux(npts, r0, z0, a, kappa, delta, kappaL=None, deltaL=None):
    r'''! Create isoflux points using simple analytic form

    @param npts Number of points to sample (evenly spaced in \f$\theta\f$)
    @param r0 Major radial position for magnetic axis
    @param z0 Vertical position for magnetic axis
    @param a Minor radius
    @param kappa Elongation (upper only if kappaL is set)
    @param delta Triangularity (upper only if deltaL is set)
    @param kappaL Lower elongation (default: kappa)
    @param deltaL Lower triangularity (default: delta)
    @result Point list [npts,2]
    '''
    kappaU = kappa
    deltaU = delta
    if kappaL is None:
        kappaL = kappaU
    if deltaL is None:
        deltaL = deltaU
    x0 = numpy.r_[r0, z0]
    isoflux_pts = numpy.zeros((npts,2))
    for i in range(npts):
        theta = i*2.0*numpy.pi/npts
        delta = ((deltaU + deltaL) + (deltaU - deltaL)*numpy.sin(theta))/2
        kappa = ((kappaU + kappaL) + (kappaU - kappaL)*numpy.sin(theta))/2
        isoflux_pts[i,:] = x0 + a*numpy.r_[numpy.cos(theta+numpy.arcsin(delta)*numpy.sin(theta)),kappa*numpy.sin(theta)]
    return isoflux_pts

def read_eqdsk(filename):
    '''! Read gEQDSK file

    @param filename Path to gEQDSK file
    @result Dictionary containing gEQDSK information
    '''
    def read_1d(fid, n):
        j = 0
        output = numpy.zeros((n,))
        for i in range(n):
            if j == 0:
                line = fid.readline()
            output[i] = line[j:j+16]
            j += 16
            if j == 16*5:
                j = 0
        return output

    def read_2d(fid, n, m):
        j = 0
        output = numpy.zeros((n, m))
        for k in range(n):
            for i in range(m):
                if j == 0:
                    line = fid.readline()
                output[k, i] = line[j:j+16]
                j += 16
                if j == 16*5:
                    j = 0
        return output
    # Read-in data
    eqdsk_obj = {}
    with open(filename, 'r') as fid:
        # Get sizes
        line = fid.readline()
        eqdsk_obj['case'] = line[:48]
        split_line = line[48:].split()
        eqdsk_obj['nr'] = int(split_line[-2])
        eqdsk_obj['nz'] = int(split_line[-1])
        # Read header content
        line_keys = [['rdim',  'zdim',  'rcentr',  'rleft',  'zmid'],
                     ['raxis', 'zaxis', 'psimag', 'psibry', 'bcentr'],
                     ['ip',    'skip',  'skip',   'skip',   'skip'],
                     ['skip',  'skip',  'skip',   'skip',   'skip']]
        for i in range(4):
            line = fid.readline()
            for j in range(5):
                if line_keys[i][j] == 'skip':
                    continue
                line_seg = line[j*16:(j+1)*16]
                eqdsk_obj[line_keys[i][j]] = float(line_seg)
        # Read flux profiles
        keys = ['fpol', 'pres', 'ffprim', 'pprime']
        for key in keys:
            eqdsk_obj[key] = read_1d(fid, eqdsk_obj['nr'])
        # Read PSI grid
        eqdsk_obj['psirz'] = read_2d(fid, eqdsk_obj['nz'],
                                        eqdsk_obj['nr'])
        # Read q-profile
        eqdsk_obj['qpsi'] = read_1d(fid, eqdsk_obj['nr'])
        # Read limiter count
        line = fid.readline()
        eqdsk_obj['nbbs'] = int(line.split()[0])
        eqdsk_obj['nlim'] = int(line.split()[1])
        # Read outer flux surface
        eqdsk_obj['rzout'] = read_2d(fid, eqdsk_obj['nbbs'], 2)
        # Read limiting corners
        eqdsk_obj['rzlim'] = read_2d(fid, eqdsk_obj['nlim'], 2)
    return eqdsk_obj

def create_prof_file(self, filename, profile_dict, name):
    '''! Create profile input file to be read by load_profiles()

    @param filename Name of input file, see options in set_profiles()
    @param profile_dict Dictionary object containing profile values ['y'] and sampled locations 
    in normalized Psi ['x']
    @param filename Name of input quantity, see options in set_profiles()
    '''
    file_lines = [profile_dict['type']]
    if profile_dict['type'] == 'flat':
        pass
    elif profile_dict['type'] == 'linterp':
        x = profile_dict.get('x',None)
        if x is None:
            raise KeyError('No array "x" for piecewise linear profile.')
        else:
            x = numpy.array(x.copy())
        y = profile_dict.get('y',None)
        if y is None:
            raise KeyError('No array "y" for piecewise linear profile.')
        else:
            y = numpy.array(y.copy())
        if numpy.min(numpy.diff(x)) < 0.0:
            raise ValueError("psi values in {0} profile must be monotonically increasing".format(name))
        if (x[0] < 0.0) or (x[-1] > 1.0):
            raise ValueError("Invalid psi values in {0} profile ({1}, {2})".format(name, x[0], x[-1]))
        if self.psi_convention == 0:
            x = 1.0 - x
            sort_inds = x.argsort()
            x = x[sort_inds]
            y = y[sort_inds]
        elif self.psi_convention == 1:
            pass
        else:
            raise ValueError('Unknown convention type, must be 0 (tokamak) or 1 (spheromak)')
        file_lines += [
            "{0} {1}".format(x.shape[0]-1, y[0]),
            "{0}".format(" ".join(["{0}".format(val) for val in x[1:]])),
            "{0}".format(" ".join(["{0}".format(val) for val in y[1:]]))
        ]
    else:
        raise KeyError('Invalid profile type ("flat", "linterp")')
    with open(filename, 'w+') as fid:
        fid.write("\n".join(file_lines))

oft_in_template = """&runtime_options
 debug={DEBUG_LEVEL}
/

&mesh_options
 meshname='none'
 {MESH_TYPE}
/

{MESH_DEF}
"""


class TokaMaker_field_interpolator():
    '''! Interpolation class for Grad-Shafranov fields'''
    def __init__(self,int_obj,int_type,dim,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param int_obj Address of FORTRAN interpolation class
        @param int_type Interpolation type (see @ref TokaMaker.TokaMaker.get_field_eval "get_field_eval")
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
        tokamaker_apply_field_eval(self.int_obj,-self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.dim,self.val)

    def eval(self,pt):
        '''! Evaluate field at a given location

        @param pt Location for evaluation [2]
        @result Field at evaluation point [self.dim]
        '''
        pt_eval = numpy.zeros((3,), dtype=numpy.float64)
        pt_eval[:2] = pt
        tokamaker_apply_field_eval(self.int_obj,self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.dim,self.val)
        return self.val


class TokaMaker():
    '''! TokaMaker G-S solver class'''
    def __init__(self,debug_level=0,nthreads=2):
        '''! Initialize TokaMaker object

        @param debug_level Level of debug printing (0-3)
        '''
        ## Input file settings
        self._oft_in_dict = {'DEBUG_LEVEL': debug_level, 'MESH_TYPE': '', 'MESH_DEF': ''}
        self._update_oft_in()
        oft_init(c_int(nthreads))
        ## Internal Grad-Shafranov object (@ref psi_grad_shaf.gs_eq "gs_eq")
        self.gs_obj = c_void_p()
        tokamaker_alloc(ctypes.byref(self.gs_obj))
        ## General settings object
        self.settings = tokamaker_default_settings()
        ## Conductor definition dictionary
        self._cond_dict = {}
        ## Coil definition dictionary
        self._coil_dict = {}
        ## Coil set definitions, including sub-coils
        self.coil_sets = {}
        ## Vacuum F value
        self._F0 = 0.0
        ## Plasma current target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._Ip_target=c_double(-1.0)
        ## Plasma current target ratio I_p(FF') / I_p(P') (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._Ip_ratio_target=c_double(-1.E99)
        ## Axis pressure target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._pax_target=c_double(-1.0)
        ## Stored energy target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._estore_target=c_double(-1.0)
        ## R0 target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._R0_target=c_double(-1.0)
        ## V0 target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._V0_target=c_double(-1.E99)
        ## F*F' normalization value [1] (use @ref TokaMaker.TokaMaker.alam "alam" property)
        self._alam = None
        ## Pressure normalization value [1] (use @ref TokaMaker.TokaMaker.pnorm "pnorm" property)
        self._pnorm = None
        ## Location of O-point (magnetic axis) [2]
        self.o_point = None
        ## Limiting point (limter or active X-point) [2]
        self.lim_point = None
        ## Location of X-points [20,2]
        self.x_points = None
        ## Diverted (limited) flag [1] (use @ref TokaMaker.TokaMaker.diverted "diverted" property)
        self._diverted = None
        ## Bounding values for \f$\psi\f$ (\f$\psi_a\f$,\f$\psi_0\f$) [2]
        self.psi_bounds = None
        ## Normalized flux convention (0 -> tokamak, 1 -> spheromak)
        self.psi_convention = 0
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
        ## Number of vacuum regions in mesh
        self.nvac = 0
        ## Limiting contour
        self.lim_contour = None
    
    def _update_oft_in(self):
        '''! Update input file (`oftpyin`) with current settings'''
        with open('oftpyin','w+') as fid:
            fid.write(oft_in_template.format(**self._oft_in_dict))

    def reset(self):
        '''! Reset G-S object to enable loading a new mesh and coil configuration'''
        cstring = c_char_p(b""*200)
        tokamaker_reset(cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        self.nregs = -1
        self.np = -1
        self._oft_in_dict['MESH_TYPE'] = ''
        self._oft_in_dict['MESH_DEF'] = ''
        # Reset defaults
        self.settings = tokamaker_default_settings()
        self._cond_dict = {}
        self._coil_dict = {}
        self._F0 = 0.0
        self._Ip_target=c_double(-1.0)
        self._Ip_ratio_target=c_double(-1.E99)
        self._pax_target=c_double(-1.0)
        self._estore_target=c_double(-1.0)
        self._R0_target=c_double(-1.0)
        self._V0_target=c_double(-1.E99)
        self._alam = None
        self._pnorm = None
        self.o_point = None
        self.lim_point = None
        self.x_points = None
        self._diverted = None
        self.psi_bounds = None
        self.nregs = -1
        self.np = -1
        self.nc = -1
        self.r = None
        self.lc = None
        self.reg = None
        self.nvac = 0
        self.lim_contour = None

    def setup_mesh(self,r=None,lc=None,reg=None,mesh_file=None):
        '''! Setup mesh for static and time-dependent G-S calculations

        A mesh should be specified by passing "r", "lc", and optionally "reg" or using a "mesh_file".
        When a region is specified the following ordering should apply:
          - 1: Plasma region
          - 2: Vacuum/air regions
          - 3+: Conducting regions and coils

        @param r Mesh point list [np,2]
        @param lc Mesh cell list [nc,3] (base one)
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
            self._oft_in_dict['MESH_TYPE'] = 'cad_type=0'
            self._oft_in_dict['MESH_DEF'] = """&native_mesh_options
 filename='{0}'
/""".format(mesh_file)
            self._update_oft_in()
            oft_setup_smesh(ndim,ndim,rfake,ndim,ndim,lcfake,regfake,ctypes.byref(nregs))
        elif r is not None:
            ndim = c_int(r.shape[1])
            np = c_int(r.shape[0])
            npc = c_int(lc.shape[1])
            nc = c_int(lc.shape[0])
            if reg is None:
                reg = numpy.ones((nc.value,),dtype=numpy.int32)
            oft_setup_smesh(ndim,np,r,npc,nc,lc+1,reg,ctypes.byref(nregs))
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        self.nregs = nregs.value
    
    def setup_regions(self,cond_dict={},coil_dict={}):
        '''! Define mesh regions (coils and conductors)

        @param cond_dict Dictionary specifying conducting regions
        '''
        self._cond_dict = cond_dict
        self._coil_dict = coil_dict
        xpoint_mask = numpy.ones((self.nregs,),dtype=numpy.int32)
        xpoint_mask[0] = 1
        eta_vals = -2.0*numpy.ones((self.nregs,),dtype=numpy.float64)
        eta_vals[0] = -1.0
        for key in cond_dict:
            eta_vals[cond_dict[key]['reg_id']-1] = cond_dict[key]['eta']/mu0
            xpoint_mask[cond_dict[key]['reg_id']-1] = int(cond_dict[key].get('allow_xpoints',False))
        nCoils = 0
        self.coil_sets = {}
        for key in coil_dict:
            xpoint_mask[coil_dict[key]['reg_id']-1] = int(coil_dict[key].get('allow_xpoints',False))
            eta_vals[coil_dict[key]['reg_id']-1] = -1.0
            coil_set = coil_dict[key].get('coil_set',key)
            if coil_set not in self.coil_sets:
                self.coil_sets[coil_set] = {
                    'id': nCoils,
                    'sub_coils': []
                }
                nCoils += 1
            self.coil_sets[coil_set]['sub_coils'].append(coil_dict[key])
        # Mark vacuum regions
        self.nvac = 0
        for i in range(self.nregs):
            if eta_vals[i] < -1.5:
                eta_vals[i] = 1.E10
                self.nvac += 1 
        coil_nturns = numpy.zeros((nCoils, self.nregs))
        for key in self.coil_sets:
            for sub_coil in self.coil_sets[key]['sub_coils']:
                coil_nturns[self.coil_sets[key]['id'],sub_coil['reg_id']-1] = sub_coil.get('nturns',1.0)
        cstring = c_char_p('none'.encode())
        tokamaker_setup_regions(cstring,eta_vals,xpoint_mask,coil_nturns,nCoils)

    def eval_green(self,x,xc):
        r'''! Evaluate Green's function for a toroidal filament

        @param x Observation point [2]
        @param xc Coil location [:,2]
        @result \f$\psi(x)\f$ due to a coil with unit current [A] at xc
        '''
        n = x.shape[0]
        vals = numpy.zeros((n,),dtype=numpy.float64)
        r = x[:,0].copy()
        z = x[:,1].copy()
        tokamaker_eval_green(c_int(n),r,z,
            c_double(xc[0]),c_double(xc[1]),vals)
        return vals*mu0
    
    def setup(self,order=2,F0=0.0,full_domain=False):
        r'''! Setup G-S solver

        @param order Order of FE representation to use
        @param F0 Vacuum \f$F(\psi)\f$ value (B0*R0)
        '''
        if self.np != -1:
            raise ValueError('G-S instance already setup')
        self.update_settings()
        #
        ncoils = c_int()
        cstring = c_char_p(b""*200)
        # filename = c_char_p(input_filename.encode())
        tokamaker_setup(order,full_domain,ctypes.byref(ncoils),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        ## Number of coils in mesh
        self.ncoils = ncoils.value
        ## Isoflux constraint points (use @ref TokaMaker.TokaMaker.set_isoflux "set_isoflux")
        self._isoflux = None
        ## Saddle constraint points (use @ref TokaMaker.TokaMaker.set_saddles "set_saddles")
        self._saddles = None
        # Get references to internal variables
        o_loc = c_double_ptr()
        lim_loc = c_double_ptr()
        x_loc = c_double_ptr()
        div_flag_loc = c_bool_ptr()
        bounds_loc = c_double_ptr()
        alam_loc = c_double_ptr()
        pnorm_loc = c_double_ptr()
        tokamaker_get_refs(ctypes.byref(o_loc),ctypes.byref(lim_loc),ctypes.byref(x_loc),ctypes.byref(div_flag_loc),
                    ctypes.byref(bounds_loc),ctypes.byref(alam_loc),ctypes.byref(pnorm_loc))
        self.o_point = numpy.ctypeslib.as_array(o_loc,shape=(2,))
        self.lim_point = numpy.ctypeslib.as_array(lim_loc,shape=(2,))
        self.x_points = numpy.ctypeslib.as_array(x_loc,shape=(20, 2))
        self._diverted = numpy.ctypeslib.as_array(div_flag_loc,shape=(1,))
        self.psi_bounds = numpy.ctypeslib.as_array(bounds_loc,shape=(2,))
        self._alam = numpy.ctypeslib.as_array(alam_loc,shape=(1,))
        self._pnorm = numpy.ctypeslib.as_array(pnorm_loc,shape=(1,))
        # Set default targets
        self.alam = 0.1
        self.pnorm = 0.1
        default_prof={
            'type': 'linterp',
            'x': numpy.array([0.0,0.8,1.0]),
            'y': numpy.array([2.0,1.0,0.0])
        }
        self.set_profiles(ffp_prof=default_prof, foffset=F0, pp_prof=default_prof)
        # Get limiter contour
        npts = c_int()
        r_loc = c_double_ptr()
        tokamaker_get_limiter(ctypes.byref(npts),ctypes.byref(r_loc))
        self.lim_contour = numpy.ctypeslib.as_array(r_loc,shape=(npts.value, 2))
        # Get plotting mesh
        np_loc = c_int()
        nc_loc = c_int()
        r_loc = c_double_ptr()
        lc_loc = c_int_ptr()
        reg_loc = c_int_ptr()
        tokamaker_get_mesh(ctypes.byref(np_loc),ctypes.byref(r_loc),ctypes.byref(nc_loc),ctypes.byref(lc_loc),ctypes.byref(reg_loc))
        ## Number of points in mesh
        self.np = np_loc.value
        ## Number of cells in mesh
        self.nc = nc_loc.value
        ## Mesh vertices [np,3] (last column should be all zeros)
        self.r = numpy.ctypeslib.as_array(r_loc,shape=(self.np, 3))
        ## Mesh triangles [nc,3] 
        self.lc = numpy.ctypeslib.as_array(lc_loc,shape=(self.nc, 3))
        ## Mesh regions [nc] 
        self.reg = numpy.ctypeslib.as_array(reg_loc,shape=(self.nc,))

    @property
    def alam(self):
        if self._alam is not None:
            return self._alam[0]
        else:
            return None
    
    @alam.setter
    def alam(self,value):
        if self._alam is not None:
            self._alam[0] = value
        else:
            raise ValueError('Class must be initialized to set "alam"')
    
    @property
    def pnorm(self):
        if self._pnorm is not None:
            return self._pnorm[0]
        else:
            return None
    
    @pnorm.setter
    def pnorm(self,value):
        if self._pnorm is not None:
            self._pnorm[0] = value
        else:
            raise ValueError('Class must be initialized to set "pnorm"')
    
    @property
    def diverted(self):
        if self._diverted is not None:
            return self._diverted[0]
        else:
            return None
    
    def set_coil_reg(self,reg_mat,reg_targets=None,reg_weights=None):
        '''! Set regularization matrix for coils when isoflux and/or saddle constraints are used

        Can be used to enforce "soft" constraints on coil currents. For hard constraints see
        @ref TokaMaker.TokaMaker.set_coil_bounds "set_coil_bounds".

        @param reg_mat Regularization matrix [nregularize,ncoils+1]
        @param reg_targets Regularization targets [nregularize] (default: 0)
        @param reg_weights Weights for regularization terms [nregularize] (default: 1)
        '''
        if reg_mat.shape[1] != self.ncoils+1:
            raise ValueError('Incorrect shape of "reg_mat", should be [nregularize,ncoils+1]')
        nregularize = reg_mat.shape[0]
        if reg_targets is None:
            reg_targets = numpy.zeros((nregularize,), dtype=numpy.float64)
        if reg_weights is None:
            reg_weights = numpy.ones((nregularize,), dtype=numpy.float64)
        if reg_targets.shape[0] != nregularize:
            raise ValueError('Incorrect shape of "reg_targets", should be [nregularize]')
        if reg_weights.shape[0] != nregularize:
            raise ValueError('Incorrect shape of "reg_weights", should be [nregularize]')
        
        tokamaker_set_coil_regmat(nregularize,numpy.copy(reg_mat.transpose(), order='C'), reg_targets, reg_weights)

    def set_coil_bounds(self,coil_bounds):
        '''! Set hard constraints on coil currents

        Can be used with or without regularization terms (see
        @ref TokaMaker.TokaMaker.set_coil_reg "set_coil_reg").

        @param coil_bounds Minimum and maximum allowable coil currents [ncoils+1,2]
        '''
        if (coil_bounds.shape[0] != self.ncoils+1) or (coil_bounds.shape[1] != 2):
            raise ValueError('Incorrect shape of "coil_bounds", should be [ncoils+1,2]')
        bounds = numpy.copy(coil_bounds, order='C')
        tokamaker_set_coil_bounds(bounds)

    def set_coil_vsc(self,coil_gains):
        '''! Define a vertical stability coil set from one or more coils

        @param coil_gains Gains for each coil (absolute scale is arbitrary)
        '''
        if coil_gains.shape[0] != self.ncoils:
            raise ValueError('Incorrect shape of "coil_gains", should be [ncoils]')
        tokamaker_set_coil_vsc(coil_gains)

    def init_psi(self, r0=-1.0, z0=0.0, a=0.0, kappa=0.0, delta=0.0):
        r'''! Initialize \f$\psi\f$ using uniform current distributions

        If r0>0 then a uniform current density inside a surface bounded by
        a curve of the form defined in @ref oftpy.create_isoflux is used.
        If r0<0 then a uniform current density over the entire plasma region is used.

        @param r0 Major radial position for flux surface-based approach
        @param z0 Vertical position for flux surface-based approach
        @param a Minor radius for flux surface-based approach
        @param kappa Elongation for flux surface-based approach
        @param delta Triangularity for flux surface-based approach
        '''
        error_flag = c_int()
        tokamaker_init_psi(c_double(r0),c_double(z0),c_double(a),c_double(kappa),c_double(delta),ctypes.byref(error_flag))
        return error_flag.value

    def load_profiles(self, f_file='f_prof.in', foffset=None, p_file='p_prof.in', eta_file='eta_prof.in', f_NI_file='f_NI_prof.in'):
        r'''! Load flux function profiles (\f$F*F'\f$ and \f$P'\f$) from files

        @param f_file File containing \f$F*F'\f$ (or \f$F'\f$ if `mode=0`) definition
        @param foffset Value of \f$F0=R0*B0\f$
        @param p_file File containing \f$P'\f$ definition
        @param eta_file File containing $\eta$ definition
        @param f_NI_file File containing non-inductive \f$F*F'\f$ definition
        '''
        if foffset is not None:
            self._F0 = foffset
        tokamaker_load_profiles(c_char_p(f_file.encode()),c_double(self._F0),c_char_p(p_file.encode()),c_char_p(eta_file.encode()),c_char_p(f_NI_file.encode()))

    def set_profiles(self, ffp_prof=None, foffset=None, pp_prof=None, ffp_NI_prof=None):
        r'''! Set flux function profiles (\f$F*F'\f$ and \f$P'\f$) using a piecewise linear definition

        @param ffp_prof Dictionary object containing FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param foffset Value of \f$F0=R0*B0\f$
        @param pp_prof Dictionary object containing P' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param ffp_NI_prof Dictionary object containing non-inductive FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        '''
        ffp_file = 'none'
        if ffp_prof is not None:
            ffp_file = 'tokamaker_f.prof'
            create_prof_file(self, ffp_file, ffp_prof, "F*F'")
        pp_file = 'none'
        if pp_prof is not None:
            pp_file = 'tokamaker_p.prof'
            create_prof_file(self, pp_file, pp_prof, "P'")
        eta_file = 'none'
        ffp_NI_file = 'none'
        if ffp_NI_prof is not None:
            ffp_NI_file = 'tokamaker_ffp_NI.prof'
            create_prof_file(self, ffp_NI_file, ffp_NI_prof, "ffp_NI")
        if foffset is not None:
            self._F0 = foffset
        self.load_profiles(ffp_file,foffset,pp_file,eta_file,ffp_NI_file)

    def set_resistivity(self, eta_prof=None):
        r'''! Set flux function profile $\eta$ using a piecewise linear definition

        Arrays should have the form array[i,:] = (\f$\hat{\psi}_i\f$, \f$f(\hat{\psi}_i)\f$) and span
        \f$\hat{\psi}_i = [0,1]\f$.

        @param eta_prof Values defining $\eta$ [:,2]
        '''
        ffp_file = 'none'
        pp_file = 'none'
        eta_file = 'none'
        if eta_prof is not None:
            eta_file = 'tokamaker_eta.prof'
            create_prof_file(self, eta_file, eta_prof, "eta'")
        ffp_NI_file = 'none'
        self.load_profiles(ffp_file,None,pp_file,eta_file,ffp_NI_file)

    def solve(self, vacuum=False):
        '''! Solve G-S equation with specified constraints, profiles, etc.'''
        error_flag = c_int()
        tokamaker_run(c_bool(vacuum),ctypes.byref(error_flag))
        return error_flag.value

    def get_stats(self,lcfs_pad=0.01):
        '''! Get information (Ip, q, kappa, etc.) about current G-S equilbirium

        @param lcfs_pad Padding at LCFS for boundary calculations
        @result Dictionary of equilibrium parameters
        '''
        _,qvals,_,dl,rbounds,zbounds = self.get_q(numpy.r_[1.0-lcfs_pad,0.95,0.02]) # Given backward so last point is LCFS (for dl)
        Ip,centroid,vol,pvol,dflux,tflux,li = self.get_globals()
        _,_,_,p,_ = self.get_profiles(numpy.r_[0.001])
        if self.diverted:
            for i in range(self.x_points.shape[0]):
                if self.x_points[i,0] < 0.0:
                    break
                x_active = self.x_points[i,:]
            if x_active[1] < zbounds[0,1]:
                zbounds[0,:] = x_active
            elif x_active[1] > zbounds[1,1]:
                zbounds[1,:] = x_active
        #
        eq_stats = {
            'Ip': Ip,
            'Ip_centroid': centroid,
            'kappa': (zbounds[1,1]-zbounds[0,1])/(rbounds[1,0]-rbounds[0,0]),
            'kappaU': (zbounds[1,1]-self.o_point[1])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'kappaL': (self.o_point[1]-zbounds[0,1])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'delta': ((rbounds[1,0]+rbounds[0,0])/2.0-(zbounds[1,0]+zbounds[0,0])/2.0)*2.0/(rbounds[1,0]-rbounds[0,0]),
            'deltaU': ((rbounds[1,0]+rbounds[0,0])/2.0-zbounds[1,0])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'deltaL': ((rbounds[1,0]+rbounds[0,0])/2.0-zbounds[0,0])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'vol': vol,
            'q_0': qvals[2],
            'q_95': qvals[1],
            'P_ax': p[0],
            'W_MHD': pvol*1.5,
            'beta_pol': 100.0*(2.0*pvol*mu0/vol)/numpy.power(Ip*mu0/dl,2),
            'dflux': dflux,
            'tflux': tflux,
            'l_i': li
        }
        if self._F0 > 0.0:
            eq_stats['beta_tor'] = 100.0*(2.0*pvol*mu0/vol)/(numpy.power(self._F0/centroid[0],2))
        return eq_stats

    def print_info(self,lcfs_pad=0.01):
        '''! Print information (Ip, q, etc.) about current G-S equilbirium'''
        eq_stats = self.get_stats(lcfs_pad=lcfs_pad)
        print("Equilibrium Statistics:")
        if self.diverted:
            print("  Topology                =   Diverted")
        else:
            print("  Topology                =   Limited")
        print("  Toroidal Current [A]    =   {0:11.4E}".format(eq_stats['Ip']))
        print("  Current Centroid [m]    =   {0:6.3F} {1:6.3F}".format(*eq_stats['Ip_centroid']))
        print("  Magnetic Axis [m]       =   {0:6.3F} {1:6.3F}".format(*self.o_point))
        print("  Elongation              =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['kappa'],eq_stats['kappaU'],eq_stats['kappaL']))
        print("  Triangularity           =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['delta'],eq_stats['deltaU'],eq_stats['deltaL']))
        print("  Plasma Volume [m^3]     =   {0:6.3F}".format(eq_stats['vol']))
        print("  q_0, q_95               =   {0:6.3F} {1:6.3F}".format(eq_stats['q_0'],eq_stats['q_95']))
        print("  Peak Pressure [Pa]      =   {0:11.4E}".format(eq_stats['P_ax']))
        print("  Stored Energy [J]       =   {0:11.4E}".format(eq_stats['W_MHD']))
        print("  <Beta_pol> [%]          =   {0:7.4F}".format(eq_stats['beta_pol']))
        if 'beta_tor' in eq_stats:
            print("  <Beta_tor> [%]          =   {0:7.4F}".format(eq_stats['beta_tor']))
        print("  Diamagnetic flux [Wb]   =   {0:11.4E}".format(eq_stats['dflux']))
        print("  Toroidal flux [Wb]      =   {0:11.4E}".format(eq_stats['tflux']))
        print("  l_i                     =   {0:7.4F}".format(eq_stats['l_i']))
    
    def set_isoflux(self,isoflux,weights=None,grad_wt_lim=-1.0):
        r'''! Set isoflux constraint points (all points lie on a flux surface)

        To constraint points more uniformly in space additional weighting based on
        the gradient of $\psi$ at each point can also be included by setting
        `grad_wt_lim>0`. When set the actual weight will be
        $w_i * min(grad_wt_lim,|\nabla \psi|_{max} / |\nabla \psi|_i)$

        @param isoflux List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        @param grad_wt_lim Limit on gradient-based weighting (negative to disable)
        '''
        if isoflux is None:
            tokamaker_set_isoflux(numpy.zeros((1,1)),numpy.zeros((1,)),0,grad_wt_lim)
            self._isoflux = None
        else:
            if weights is None:
                weights = numpy.ones((isoflux.shape[0],), dtype=numpy.float64)
            if weights.shape[0] != isoflux.shape[0]:
                raise ValueError('Shape of "weights" does not match first dimension of "isoflux"')
            tokamaker_set_isoflux(isoflux,weights,isoflux.shape[0],grad_wt_lim)
            self._isoflux = isoflux.copy()
    
    def set_saddles(self,saddles,weights=None):
        '''! Set saddle constraint points (poloidal field should vanish at each point)

        @param saddles List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        if saddles is None:
            tokamaker_set_saddles(numpy.zeros((1,1)),numpy.zeros((1,)),0)
            self._saddles = None
        else:
            if weights is None:
                weights = numpy.ones((saddles.shape[0],), dtype=numpy.float64)
            if weights.shape[0] != saddles.shape[0]:
                raise ValueError('Shape of "weights" does not match first dimension of "saddles"')
            tokamaker_set_saddles(saddles,weights,saddles.shape[0])
            self._saddles = saddles.copy()
    
    def set_targets(self,Ip=None,Ip_ratio=None,pax=None,estore=None,R0=None,V0=None,retain_previous=False):
        r'''! Set global target values

        Once set, values are retained until they are explicitly set to their respective disabled
        values (see below). By default, all targets are disabled so this function should be called
        at least once to set "sane" values for `alam` and `pnorm`.

        @param alam Scale factor for \f$F*F'\f$ term (disabled if `Ip` is set)
        @param pnorm Scale factor for \f$P'\f$ term (disabled if `pax`, `estore`, or `R0` are set)
        @param Ip Target plasma current [A] (disabled if <0)
        @param Ip_ratio Amplitude of net plasma current contribution from FF' compared to P' (disabled if <-1.E98)
        @param pax Target axis pressure [Pa] (disabled if <0 or if `estore` is set)
        @param estore Target sotred energy [J] (disabled if <0)
        @param R0 Target major radius for magnetic axis (disabled if <0 or if `pax` or `estore` are set)
        @param V0 Target vertical position for magnetic axis (disabled if <-1.E98)
        @param retain_previous Keep previously set targets unless explicitly updated? (default: False)
        '''
        # Reset all targets unless specified
        if not retain_previous:
            self._Ip_target.value = -1.E99
            self._estore_target.value = -1.0
            self._pax_target.value = -1.0
            self._Ip_ratio_target.value = -1.E99
            self._R0_target.value = -1.0
            self._V0_target.value = -1.E99
        # Set new targets
        if Ip is not None:
            self._Ip_target.value=Ip
        if estore is not None:
            self._estore_target.value=estore
        if pax is not None:
            self._pax_target.value=pax
        if Ip_ratio is not None:
            self._Ip_ratio_target.value=Ip_ratio
        if R0 is not None:
            self._R0_target.value=R0
        if V0 is not None:
            self._V0_target.value=V0
        tokamaker_set_targets(self._Ip_target,self._Ip_ratio_target,self._pax_target,self._estore_target,self._R0_target,self._V0_target)

    def get_psi(self,normalized=True):
        r'''! Get poloidal flux values on node points

        @param normalized Normalize (and offset) poloidal flux
        @result \f$\hat{\psi} = \frac{\psi-\psi_0}{\psi_a-\psi_0} \f$ or \f$\psi\f$
        '''
        psi = numpy.zeros((self.np,),dtype=numpy.float64)
        psi_lim = c_double()
        psi_max = c_double()
        tokamaker_get_psi(psi,ctypes.byref(psi_lim),ctypes.byref(psi_max))
        if normalized:
            psi = (psi-psi_lim.value)/(psi_max.value-psi_lim.value)
            if self.psi_convention == 0:
                psi = 1.0 - psi
        return psi

    def set_psi(self,psi):
        '''! Set poloidal flux values on node points

        @param psi Poloidal flux values (should not be normalized!)
        '''
        if psi.shape[0] != self.np:
            raise ValueError('Incorrect shape of "psi", should be [np]')
        tokamaker_set_psi(psi)
    
    def set_psi_dt(self,psi0,dt):
        '''! Set reference poloidal flux and time step for eddy currents in .solve()

        @param psi0 Reference poloidal flux at t-dt (unnormalized)
        @param dt Time since reference poloidal flux
        '''
        if psi0.shape[0] != self.np:
            raise ValueError('Incorrect shape of "psi0", should be [np]')
        tokamaker_set_psi_dt(psi0,c_double(dt))
    
    def get_field_eval(self,field_type):
        r'''! Create field interpolator for vector potential

        @param imode Index of eigenstate
        @result Field interpolation object
        '''
        #
        mode_map = {'B': 1, 'PSI': 2, 'F': 3, 'P': 4}
        imode = mode_map.get(field_type.upper())
        if imode is None:
            raise ValueError('Invalid field type ("B", "psi", "F", "P")')
        #
        int_obj = c_void_p()
        cstring = c_char_p(b""*200)
        tokamaker_get_field_eval(imode,ctypes.byref(int_obj),cstring)
        if cstring.value != b'':
            raise Exception(cstring.value)
        field_dim = 1
        if imode == 1:
            field_dim = 3
        return TokaMaker_field_interpolator(int_obj,imode,field_dim)
    
    def get_coil_currents(self):
        '''! Get currents in each coil [A] and coil region [A-turns]

        @result Coil currents [ncoils], Coil currents by region [nregs]
        '''
        currents = numpy.zeros((self.ncoils,),dtype=numpy.float64)
        currents_reg = numpy.zeros((self.nregs,),dtype=numpy.float64)
        tokamaker_get_coil_currents(currents, currents_reg)
        return currents, currents_reg

    def get_coil_Lmat(self):
        r'''! Get mutual inductance matrix between coils

        @note This is the inductance in terms of A-turns. To get in terms of
        current in a single of the \f$n\f$ windings you must multiply by \f$n_i*n_j\f$.

        @result L[ncoils+1,ncoils+1]
        '''
        Lmat = numpy.zeros((self.ncoils+1,self.ncoils+1),dtype=numpy.float64)
        tokamaker_get_coil_Lmat(Lmat)
        return Lmat
    
    def trace_surf(self,psi):
        r'''! Trace surface for a given poloidal flux

        @param psi Flux surface to trace \f$\hat{\psi}\f$
        @result \f$r(\hat{\psi})\f$
        '''
        if self.psi_convention == 0:
            psi = 1.0-psi
        npoints = c_int()
        points_loc = c_double_ptr()
        tokamaker_trace_surf(c_double(psi),ctypes.byref(points_loc),ctypes.byref(npoints))
        if npoints.value > 0:
            return numpy.ctypeslib.as_array(points_loc,shape=(npoints.value, 2))
        else:
            return None
    
    def get_q(self,psi=None,psi_pad=0.02,npsi=50):
        r'''! Get q-profile at specified or uniformly spaced points

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @result \f$\hat{\psi}\f$, \f$q(\hat{\psi})\f$, \f$[<R>,<1/R>]\f$, length of last surface,
        [r(R_min),r(R_max)], [r(z_min),r(z_max)]
        '''
        if psi is None:
            psi = numpy.linspace(psi_pad,1.0-psi_pad,npsi,dtype=numpy.float64)
            if self.psi_convention == 0:
                psi = numpy.flip(psi).copy()
                psi_save = 1.0 - psi
        else:
            if self.psi_convention == 0:
                psi_save = psi.copy()
                psi = 1.0-psi
        qvals = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        ravgs = numpy.zeros((2,psi.shape[0]), dtype=numpy.float64)
        dl = c_double()
        rbounds = numpy.zeros((2,2),dtype=numpy.float64)
        zbounds = numpy.zeros((2,2),dtype=numpy.float64)
        tokamaker_get_q(psi.shape[0],psi,qvals,ravgs,ctypes.byref(dl),rbounds,zbounds)
        if self.psi_convention == 0:
            return psi_save,qvals,ravgs,dl.value,rbounds,zbounds
        else:
            return psi,qvals,ravgs,dl.value,rbounds,zbounds

    def sauter_fc(self,psi=None,psi_pad=0.02,npsi=50):
        r'''! Evaluate Sauter trapped particle fractions at specified or uniformly spaced points

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @result \f$ f_c \f$, [\f$<R>,<1/R>,<a>\f$], [\f$<|B|>,<|B|^2>\f$]
        '''
        if psi is None:
            psi = numpy.linspace(psi_pad,1.0-psi_pad,npsi,dtype=numpy.float64)
            if self.psi_convention == 0:
                psi = numpy.flip(psi).copy()
                psi_save = 1.0 - psi
        else:
            if self.psi_convention == 0:
                psi_save = psi.copy()
                psi = 1.0-psi
        fc = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        r_avgs = numpy.zeros((3,psi.shape[0]), dtype=numpy.float64)
        modb_avgs = numpy.zeros((2,psi.shape[0]), dtype=numpy.float64)
        tokamaker_sauter_fc(psi.shape[0],psi,fc,r_avgs,modb_avgs)
        if self.psi_convention == 0:
            return psi_save,fc,r_avgs,modb_avgs
        else:
            return psi,fc,r_avgs,modb_avgs

    def get_globals(self):
        r'''! Get global plasma parameters

        @result Ip, [R_Ip, Z_Ip], \f$\int dV\f$, \f$\int P dV\f$, diamagnetic flux,
        enclosed toroidal flux
        '''
        Ip = c_double()
        centroid = numpy.zeros((2,),dtype=numpy.float64)
        vol = c_double()
        pvol = c_double()
        dflux = c_double()
        tflux = c_double()
        Li = c_double()
        tokamaker_get_globals(ctypes.byref(Ip),centroid,ctypes.byref(vol),ctypes.byref(pvol),
            ctypes.byref(dflux),ctypes.byref(tflux),ctypes.byref(Li))
        return Ip.value, centroid, vol.value, pvol.value, dflux.value, tflux.value, Li.value

    def calc_loopvoltage(self):
        r'''! Get plasma loop voltage

        @param eta Dictionary object containing resistivity profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param ffp_NI Dictionary object containing non-inductive FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @result Vloop [Volts]
        '''
        V_loop = c_double()

        tokamaker_gs_calc_vloop(ctypes.byref(V_loop))

        if V_loop.value < 0.:
            raise ValueError('eta array not specified')
        else:
            return V_loop.value

    def get_profiles(self,psi=None,psi_pad=1.E-8,npsi=50):
        r'''! Get G-S source profiles

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @result \f$\hat{\psi}\f$, \f$F(\hat{\psi})\f$, \f$F'(\hat{\psi})\f$,
        \f$P(\hat{\psi})\f$, \f$P'(\hat{\psi})\f$
        '''
        if psi is None:
            psi = numpy.linspace(psi_pad,1.0-psi_pad,npsi,dtype=numpy.float64)
            if self.psi_convention == 0:
                psi = numpy.flip(psi).copy()
                psi_save = 1.0 - psi
        else:
            if self.psi_convention == 0:
                psi_save = psi.copy()
                psi = 1.0-psi
        f = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        fp = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        p = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        pp = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        tokamaker_get_profs(psi.shape[0],psi,f,fp,p,pp)
        if self.psi_convention == 0:
            return psi_save,f,fp,p/mu0,pp/mu0
        else:
            return psi,f,fp,p/mu0,pp/mu0
    
    def get_xpoints(self):
        '''! Get X-points

        @result X-points, is diverted?
        '''
        if self.x_points[0,0] < 0.0:
            return None, False
        else:
            for i in range(self.x_points.shape[0]):
                if self.x_points[i,0] < 0.0:
                    break
            return self.x_points[:i,:], self.diverted
    
    def set_coil_currents(self, currents):
        '''! Set coil currents

        @param currents Current in each coil [A]
        '''
        if currents.shape[0] != self.ncoils:
            raise ValueError('Incorrect shape of "currents", should be [ncoils]')
        tokamaker_set_coil_currents(currents)

    def update_settings(self):
        '''! Update settings after changes to values in python'''
        tokamaker_set_settings(ctypes.byref(self.settings))
    
    def plot_machine(self,fig,ax,vacuum_color='whitesmoke',cond_color='gray',limiter_color='k',
                     coil_color='gray',coil_colormap=None,coil_symmap=False,coil_scale=1.0,coil_clabel=r'$I_C$ [A]'):
        '''! Plot machine geometry

        @param fig Figure to add to
        @param ax Axis to add to
        @param vacuum_color Color to shade vacuum region (None to disable)
        @param cond_color Color for conducting regions (None to disable)
        @param limiter_color Color for limiter contour (None to disable)
        @param coil_color Color for coil regions (None to disable)
        @param coil_colormap Colormap for coil current values
        @param coil_symmap Make coil current colorscale symmetric
        @param coil_scale Scale for coil currents when plotting
        @param coil_clabel Label for coil current colorbar (None to disable colorbar)
        '''
        mask_vals = numpy.ones((self.np,))
        # Shade vacuum region
        if vacuum_color is not None:
            mask = numpy.logical_and(self.reg > 1, self.reg <= self.nvac+1)
            if mask.sum() > 0.0:
                ax.tricontourf(self.r[:,0], self.r[:,1], self.lc[mask,:], mask_vals, colors=vacuum_color)
        # Shade coils
        if coil_colormap is not None:
            _, region_currents = self.get_coil_currents()
            mesh_currents = numpy.zeros((self.lc.shape[0],))
            for i in range(self.ncoils):
                mesh_currents = region_currents[self.reg-1]
            mask = (abs(mesh_currents) > 0.0)
            if mask.sum() > 0.0:
                mesh_currents *= coil_scale
                if coil_symmap:
                    max_curr = abs(mesh_currents).max()
                    clf = ax.tripcolor(self.r[:,0], self.r[:,1], self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap, vmin=-max_curr, vmax=max_curr)
                else:
                    clf = ax.tripcolor(self.r[:,0], self.r[:,1], self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap)
                if coil_clabel is not None:
                    fig.colorbar(clf,ax=ax,label=coil_clabel)
        else:
            for _, coil_reg in self._coil_dict.items():
                mask_tmp = (self.reg == coil_reg['reg_id'])
                ax.tricontourf(self.r[:,0], self.r[:,1], self.lc[mask_tmp,:], mask_vals, colors=coil_color, alpha=1)
        # Shade conductors
        for _, cond_reg in self._cond_dict.items():
            mask_tmp = (self.reg == cond_reg['reg_id'])
            ax.tricontourf(self.r[:,0], self.r[:,1], self.lc[mask_tmp,:], mask_vals, colors=cond_color, alpha=1)
        # Show limiter
        if limiter_color and (self.lim_contour is not None):
            ax.plot(self.lim_contour[:,0],self.lim_contour[:,1],color=limiter_color)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')

    def plot_constraints(self,fig,ax,isoflux_color='tab:red',isoflux_marker='+',saddle_color='tab:green',saddle_marker='x'):
        '''! Plot geometry constraints

        @param fig Figure to add to
        @param ax Axis to add to
        @param isoflux_color Color of isoflux points (None to disable)
        @param saddle_color Color of saddle points (None to disable)
        '''
        # Plot isoflux constraints
        if (isoflux_color is not None) and (self._isoflux is not None):
            ax.plot(self._isoflux[:,0],self._isoflux[:,1],color=isoflux_color,marker=isoflux_marker,linestyle='none')
        # Plot saddle constraints
        if (saddle_color is not None) and (self._saddles is not None):
            ax.plot(self._saddles[:,0],self._saddles[:,1],color=saddle_color,marker=saddle_marker,linestyle='none')

    def plot_psi(self,fig,ax,psi=None,normalized=True,plasma_color=None,plasma_nlevels=8,plasma_levels=None,plasma_colormap=None,
                 vacuum_color='darkgray',vacuum_nlevels=8,vacuum_levels=None,vacuum_colormap=None,
                 xpoint_color='k',xpoint_marker='x',opoint_color='k',opoint_marker='*'):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param psi Flux values to plot (otherwise `self.get_psi()` is called)
        @param normalized Retreive normalized flux, or assume normalized psi if passed as argument
        @param plasma_color Color for plasma contours
        @param plasma_nlevels Number of plasma contours
        @param plasma_levels Explicit levels for plasma contours
        @param plasma_colormap Colormap for plasma contours (cannot be specified with `plasma_color`)
        @param vacuum_color Color for plasma contours
        @param vacuum_nlevels Number of plasma contours
        @param vacuum_levels Explicit levels for plasma contours (cannot be specified with `vacuum_color`)
        @param vacuum_colormap Colormap for plasma contours
        @param xpoint_color Color for X-point markers (None to disable)
        @param xpoint_marker Colormap for plasma contours
        @param opoint_color Colormap for plasma contours (None to disable)
        @param opoint_marker Colormap for plasma contours
        '''
        # Plot poloidal flux
        if psi is None:
            psi = self.get_psi(normalized)
        if normalized and (self.psi_convention == 0):
            psi = 1.0-psi
        if plasma_levels is None:
            if normalized:
                plasma_levels = numpy.linspace(0.0,1.0,plasma_nlevels)
            else:
                plasma_levels = numpy.linspace(psi.min(),psi.max(),plasma_nlevels)
        else:
            if normalized:
                if self.psi_convention == 0:
                    plasma_levels = sorted(1.0-numpy.array(plasma_levels))
                else:
                    plasma_levels = sorted(numpy.array(plasma_levels))
        if vacuum_levels is None:
            if normalized:
                vacuum_levels1 = numpy.zeros((0,))
                vacuum_levels2 = numpy.zeros((0,))
                if psi.min() < -0.1:
                    vacuum_levels1 = numpy.linspace(psi.min(),0.0,vacuum_nlevels,endpoint=False)
                if psi.max() > 1.1:
                    vacuum_levels2 = numpy.linspace(1.0,psi.max(),vacuum_nlevels,endpoint=False)
                vacuum_levels = numpy.hstack((vacuum_levels1,vacuum_levels2))
        else:
            if normalized:
                if self.psi_convention == 0:
                    vacuum_levels = sorted(1.0-numpy.array(vacuum_levels))
                else:
                    vacuum_levels = sorted(numpy.array(vacuum_levels))
        if (plasma_color is None) and (plasma_colormap is None):
            plasma_colormap='viridis'
        if vacuum_levels is not None:
            ax.tricontour(self.r[:,0],self.r[:,1],self.lc,psi,levels=vacuum_levels,colors=vacuum_color,cmap=vacuum_colormap)
        if plasma_levels is not None:
            ax.tricontour(self.r[:,0],self.r[:,1],self.lc,psi,levels=plasma_levels,colors=plasma_color,cmap=plasma_colormap)

        # Plot saddle points
        if xpoint_color is not None:
            x_points, _ = self.get_xpoints()
            if x_points is not None:
                ax.plot(x_points[:,0], x_points[:,1], color=xpoint_color, marker=xpoint_marker, linestyle='none')
        if (opoint_color is not None) and (self.o_point[0] > 0.0):
            ax.plot(self.o_point[0], self.o_point[1], color=opoint_color, marker=opoint_marker)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
    
    def plot_eddy(self,fig,ax,dpsi_dt=None,nlevels=40,colormap='jet',clabel=r'$J_w$ [$A/m^2$]'):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param dpsi_dt dPsi/dt corresponding to eddy currents (eg. from time-dependent simulation)
        @param nlevels Number contour lines used for shading
        @param colormap Colormap to use for shadings
        @param clabel Label for colorbar (None to disable colorbar)
        '''
        # Apply 1/R scale (avoiding divide by zero)
        dpsi_dt = dpsi_dt.copy()
        dpsi_dt[self.r[:,0]>0.0] /= self.r[self.r[:,0]>0.0,0]
        # Loop over conducting regions and get mask/fields
        mesh_currents = numpy.zeros((self.lc.shape[0],))
        mask = numpy.zeros((self.lc.shape[0],), dtype=numpy.int32)
        for name, cond_reg in self._cond_dict.items():
            eta = cond_reg.get('eta',-1.0)
            if eta > 0:
                mask_tmp = (self.reg == cond_reg['reg_id'])
                field_tmp = dpsi_dt/eta
                mesh_currents[mask_tmp] = numpy.sum(field_tmp[self.lc[mask_tmp,:]],axis=1)/3.0
                mask = numpy.logical_or(mask,mask_tmp)
        clf = ax.tripcolor(self.r[:,0],self.r[:,1],self.lc[mask],mesh_currents[mask],cmap=colormap)
        if clabel is not None:
            cb = fig.colorbar(clf,ax=ax)
            cb.set_label(clabel)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')

    def get_vfixed(self):
        '''! Get required vacuum flux values to balance fixed boundary equilibrium

        @result sampling points [:,2], flux values [:]
        '''
        npts = c_int()
        pts_loc = c_double_ptr()
        flux_loc = c_double_ptr()
        tokamaker_get_vfixed(ctypes.byref(npts),ctypes.byref(pts_loc),ctypes.byref(flux_loc))
        return numpy.ctypeslib.as_array(pts_loc,shape=(npts.value, 2)), \
            numpy.ctypeslib.as_array(flux_loc,shape=(npts.value,))

    def save_eqdsk(self,filename,nr=65,nz=65,rbounds=None,zbounds=None,run_info='',lcfs_pad=0.01):
        '''! Save current equilibrium to gEQDSK format

        @param filename Filename to save equilibrium to
        @param nr Number of radial sampling points
        @param nz Number of vertical sampling points
        @param rbounds Extents of grid in R
        @param zbounds Extents of grid in Z
        @param run_info Run information for EQDSK file (maximum of 36 characters)
        @param lcfs_pad Padding in normalized flux at LCFS
        '''
        cfilename = c_char_p(filename.encode())
        if len(run_info) > 36:
            raise ValueError('"run_info" cannot be longer than 36 characters')
        crun_info = c_char_p(run_info.encode())
        if rbounds is None:
            rbounds = numpy.r_[self.lim_contour[:,0].min(), self.lim_contour[:,0].max()]
            dr = rbounds[1]-rbounds[0]
            rbounds += numpy.r_[-1.0,1.0]*dr*0.05
        if zbounds is None:
            zbounds = numpy.r_[self.lim_contour[:,1].min(), self.lim_contour[:,1].max()]
            dr = zbounds[1]-zbounds[0]
            zbounds += numpy.r_[-1.0,1.0]*dr*0.05
        tokamaker_save_eqdsk(cfilename,c_int(nr),c_int(nz),rbounds,zbounds,crun_info,c_double(lcfs_pad))

    def eig_wall(self,neigs=4,pm=False):
        '''! Compute eigenvalues (1 / Tau_L/R) for conducting structures

        @param neigs Number of eigenvalues to compute
        @param pm Print solver statistics and raw eigenvalues?
        @result eigenvalues[neigs], eigenvectors[neigs,:]
        '''
        eig_vals = numpy.zeros((neigs,2),dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.np),dtype=numpy.float64)
        tokamaker_eig_wall(c_int(neigs),eig_vals,eig_vecs,pm)
        return eig_vals, eig_vecs

    def eig_td(self,omega=-1.E4,neigs=4,include_bounds=True,pm=False):
        '''! Compute eigenvalues for the linearized time-dependent system

        @param omega Growth rate localization point (eigenvalues closest to this value will be found)
        @param neigs Number of eigenvalues to compute
        @param include_bounds Include bounding flux terms for constant normalized profiles?
        @param pm Print solver statistics and raw eigenvalues?
        @result eigenvalues[neigs], eigenvectors[neigs,:]
        '''
        eig_vals = numpy.zeros((neigs,2),dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.np),dtype=numpy.float64)
        tokamaker_eig_td(c_double(omega),c_int(neigs),eig_vals,eig_vecs,c_bool(include_bounds),pm)
        return eig_vals, eig_vecs

    def setup_td(self,dt,lin_tol,nl_tol,pre_plasma=False):
        '''! Setup the time-dependent G-S solver

        @param dt Starting time step
        @param lin_tol Tolerance for linear solver
        @param nl_tol Tolerance for non-linear solver
        @param pre_plasma Use plasma contributions in preconditioner (default: False)
        '''
        tokamaker_setup_td(c_double(dt),c_double(lin_tol),c_double(nl_tol),c_bool(pre_plasma))
    
    def step_td(self,time,dt):
        '''! Compute eigenvalues for the time-dependent system

        @param time Growth rate enhancement point (should be approximately expected value)
        @param dt Number of eigenvalues to compute
        @result new time, new dt, # of NL iterations, # of linear iterations, # of retries
        '''
        dt = c_double(dt)
        time = c_double(time)
        nl_its = c_int()
        lin_its = c_int()
        nretry = c_int()
        tokamaker_step_td(ctypes.byref(time),ctypes.byref(dt),ctypes.byref(nl_its),ctypes.byref(lin_its),ctypes.byref(nretry))
        return time.value, dt.value, nl_its.value, lin_its.value, nretry.value


class gs_Domain:
    '''! Grad-Sharvanov domain definitions for TokaMaker with Triangle library meshing'''
    def __init__(self,rextent=None,zextents=[None,None],rpad=1.2,zpad=[1.2,1.2],json_filename=None):
        '''! Create a new Grad-Shafranov domain

        @param region_list List of @ref oftpy.Region objects that define mesh
        @param merge_thresh Distance threshold for merging nearby points
        '''
        if json_filename is not None:
            with open(json_filename, 'r') as fid:
                input_dict = json.load(fid)
            self.rextent = input_dict['rextent']
            self.zextents = input_dict['zextents']
            self.rpad = input_dict['rpad']
            self.zpad = input_dict['zpad']
            self.rmax = input_dict['rmax']
            self.zmin = input_dict['zmin']
            self.zmax = input_dict['zmax']
            self.boundary_reg = input_dict['boundary_reg']
            self.reg_type_counts = input_dict['reg_type_counts']
            self.region_info = input_dict['region_info']
            self.regions = []
            for region in input_dict['regions']:
                region.append(Region(load_dict=region))
        else:
            self.rextent = rextent
            self.zextents = zextents
            self.rpad = None
            if self.rextent is None:
                self.rpad = rpad
            self.zpad = [None,None]
            if self.zextents[0] is None:
                self.zpad[0] = zpad[0]
            if self.zextents[1] is None:
                self.zpad[1] = zpad[1]
            self.rmax = 0.0
            self.zmin = 0.0
            self.zmax = 0.0
            self.boundary_reg = None
            self.regions = []
            self.reg_type_counts = {
                'plasma': 0,
                'vacuum': 0,
                'boundary': 0,
                'conductor': 0,
                'coil': 0,
            }
            self.region_info = {}
            self._extra_reg_defs = []
    
    def define_region(self,name,dx,reg_type,eta=None,nTurns=None,coil_set=None,allow_xpoints=False):
        '''! Define a new region and its properties (geometry is given in a separate call)

        @param name Name of region
        @param dx Target mesh size for region
        @param reg_type Type of region ("plasma", "vacuum", "boundary", "conductor", or "coil")
        @param eta Resistivity for "conductor" regions (raises error if region is other type)
        @param nTurns Number of turns for "cooil" regions (raises error if region is other type)
        @param allow_xpoints Allow X-points in this region (for non-plasma regions only)
        '''
        if (dx is None) or (dx < 0.0):
            raise ValueError('"dx" must have a non-negative value')
        name = name.upper()
        if (name in self.region_info):
            raise KeyError('Region already exists!')
        next_id = -1
        if reg_type == 'plasma':
            next_id = 1
            allow_xpoints = True
        elif reg_type == 'vacuum':
            pass
        elif reg_type == 'boundary':
            self.boundary_reg = name
        elif reg_type in ('conductor', 'coil'):
            pass
        else:
            raise ValueError("Unknown region type")
        if next_id < 0:
            next_id = len(self.region_info) - self.reg_type_counts['plasma'] + 2
        self.reg_type_counts[reg_type] += 1
        self.region_info[name] = {
            'id': next_id,
            'dx': dx,
            'count': 0,
            'type': reg_type,
            'allow_xpoints': allow_xpoints
        }
        if eta is not None:
            if reg_type != 'conductor':
                raise ValueError('Resistivity specification only valid for "conductor" regions')
            else:
                self.region_info[name]['eta'] = eta
        else:
            if reg_type == 'conductor':
                raise ValueError('Resistivity not specified for "conductor" region')
        if nTurns is not None:
            if reg_type != 'coil':
                raise ValueError('"nTurns" specification only valid for "coil" regions')
            else:
                self.region_info[name]['nturns'] = nTurns
        if coil_set is not None:
            if reg_type != 'coil':
                raise ValueError('"coil_set" specification only valid for "coil" regions')
            else:
                self.region_info[name]['coil_set'] = coil_set
        

    def add_annulus(self,inner_countour,inner_name,outer_contour,annulus_name,parent_name=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None):
        '''! Add annular geometry defining region boundaries to the mesh

        @param inner_contour Curve defining inner boundary
        @param inner_name Name of region enclosed by the inner boundary
        @param outer_contour Curve defining outer boundary
        @param annulus_name Name of annular region between inner and outer boundaries
        @param parent_name Name of region beyond the outer boundary
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        '''
        inner_countour = numpy.array(inner_countour)
        outer_contour = numpy.array(outer_contour)
        inner_name = inner_name.upper()
        annulus_name = annulus_name.upper()
        inner_reg = self.region_info.get(inner_name, None)
        annulus_reg = self.region_info.get(annulus_name, None)
        if inner_reg is None:
            raise KeyError('Region "{0}" not defined'.format(inner_name))
        else:
            inner_dx = inner_reg.get('dx', None)
            inner_dx_curve = inner_dx
            if inner_dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(inner_name))
            if inner_countour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "inner_countour"')
        if annulus_reg is None:
            raise KeyError('Region "{0}" not defined'.format(annulus_name))
        else:
            annulus_dx = annulus_reg.get('dx', None)
            outer_dx_curve = annulus_dx
            if annulus_dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(annulus_name))
            else:
                inner_dx_curve = min(inner_dx_curve,annulus_dx)
            if outer_contour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "outer_contour"')
        if parent_name is not None:
            parent_name = parent_name.upper()
            parent_reg = self.region_info.get(parent_name, None)
            if parent_reg is None:
                raise KeyError('Region "{0}" not defined'.format(parent_name))
            else:
                parent_dx = parent_reg.get('dx', None)
                if parent_dx is None:
                    raise ValueError('Resolution for region "{0}" not defined'.format(parent_name))
                else:
                    outer_dx_curve = min(outer_dx_curve,parent_dx)
        # Add inner region
        maxes = inner_countour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,inner_countour[:,1].min())
        self.regions.append(Region(inner_countour,inner_dx,inner_dx_curve,angle_tol,sliver_tol,small_thresh,inner_reg["id"]))
        inner_reg["count"] += 1
        # Add outer region
        maxes = outer_contour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,outer_contour[:,1].min())
        self.regions.append(Region(outer_contour,annulus_dx,outer_dx_curve,angle_tol,sliver_tol,small_thresh,annulus_reg["id"]))
        annulus_reg["count"] += 1
    
    def add_polygon(self,contour,name,parent_name=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None):
        '''! Add polygon geometry defining region boundaries to the mesh

        @param contour Curve defining polygon
        @param name Name of region enclosed by the polygon
        @param parent_name Name of region outside the polygon
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        '''
        contour = numpy.array(contour)
        name = name.upper()
        reg = self.region_info.get(name, None)
        if reg is None:
            raise KeyError('Region "{0}" not defined'.format(name))
        else:
            dx = reg.get('dx', None)
            dx_curve = dx
            if dx is None:
                raise ValueError('Resolution for region "{0}" not defined'.format(name))
            if contour[:,0].min() < 0.0:
                raise ValueError('Negative radial value detected in "inner_countour"')
        if parent_name is not None:
            parent_name = parent_name.upper()
            parent_reg = self.region_info.get(parent_name, None)
            if parent_reg is None:
                raise KeyError('Region "{0}" not defined'.format(parent_name))
            else:
                parent_dx = parent_reg.get('dx', None)
                if parent_dx is None:
                    raise ValueError('Resolution for region "{0}" not defined'.format(parent_name))
                else:
                    dx = min(dx,parent_dx)
        # Add region
        maxes = contour.max(axis=0)
        self.rmax = max(self.rmax,maxes[0])
        self.zmax = max(self.zmax,maxes[1])
        self.zmin = min(self.zmin,contour[:,1].min())
        self.regions.append(Region(contour,dx,dx_curve,angle_tol,sliver_tol,small_thresh,reg["id"]))
        reg["count"] += 1

    def add_rectangle(self,rc,zc,w,h,name,parent_name=None, rot=None):
        '''! Add rectangular geometry defining region boundaries to the mesh

        @param rc Radial center of rectangle
        @param zc Vertical center of rectangle
        @param w Width of rectangle (radial direction)
        @param h Height of the rectangle (vertical direction)
        @param name Name of region enclosed by the polygon
        @param parent_name Name of region outside the polygon
        @param rot Rotation of rectangle (degrees)
        '''
        contour = numpy.asarray([
            [-w/2.0, -h/2.0],
            [+w/2.0, -h/2.0],
            [+w/2.0, +h/2.0],
            [-w/2.0, +h/2.0]
        ])

        if rot is not None:
            rot = numpy.deg2rad(rot)
            rotmat = numpy.asarray([numpy.cos(rot), -numpy.sin(rot), numpy.sin(rot), numpy.cos(rot)]).reshape((2,2))

            contour = numpy.dot(contour,rotmat.T)
        
        contour[:,0] += rc
        contour[:,1] += zc

        self.add_polygon(contour,name,parent_name)
    
    def add_enclosed(self,in_point,name):
        name = name.upper()
        reg = self.region_info.get(name, None)
        if reg is None:
            raise KeyError('Region "{0}" not defined'.format(name))
        else:
            dx = reg['dx']
            id = reg['id']
            self._extra_reg_defs.append([in_point[0], in_point[1], id, dx*dx/2.0])
            reg["count"] += 1
    
    def get_coils(self):
        '''! Get dictionary describing coil regions in domain

        @result Dictionary of coil regions and attributes
        '''
        coil_list = {}
        coil_id = 0
        for key in self.region_info:
            if self.region_info[key]['type'] == 'coil':
                coil_list[key] = {
                    'reg_id': self.region_info[key]['id'],
                    'coil_id': coil_id,
                    'nturns': self.region_info[key].get('nturns',1),
                    'coil_set': self.region_info[key].get('coil_set',key),
                    'allow_xpoints': self.region_info[key].get('allow_xpoints',False)
                }
                coil_id += 1
        return coil_list

    def get_conductors(self):
        '''! Get dictionary describing conducting regions in domain

        @result Dictionary of conducting regions and attributes
        '''
        cond_list = {}
        cond_id = 0
        for key in self.region_info:
            if self.region_info[key]['type'] == 'conductor':
                cond_list[key] = {
                    'reg_id': self.region_info[key]['id'],
                    'cond_id': cond_id,
                    'eta': self.region_info[key]['eta'],
                    'allow_xpoints': self.region_info[key].get('allow_xpoints',False)
                }
                cond_id += 1
        return cond_list
    
    def build_mesh(self,debug=False,merge_thresh=1.E-4,require_boundary=True,setup_only=False):
        '''! Build mesh for specified domains

        @result Meshed representation (pts[np,2], tris[nc,3], regions[nc])
        '''
        # Check for single plasma region
        if self.reg_type_counts['plasma'] > 1:
            raise ValueError('More than one plasma region specified')
        elif self.reg_type_counts['plasma'] == 0:
            raise ValueError('No plasma region specified')
        else:
            # Make sure a boundary exists if we have regions other than plasma
            if ((self.reg_type_counts['vacuum'] > 0) or (self.reg_type_counts['coil'] > 0) or (self.reg_type_counts['conductor'] > 0)) and require_boundary:
                if self.boundary_reg is None:
                    raise ValueError('No boundary region specified')
                # Check or set extents
                if self.rextent is None:
                    self.rextent = self.rpad*self.rmax
                else:
                    if self.rmax > self.rextent:
                        raise ValueError('User specified "rextent" does not enclose all regions')
                if self.zextents[0] is None:
                    self.zextents[0] = self.zpad[0]*self.zmin
                else:
                    if self.zmin < self.zextents[0]:
                        raise ValueError('User specified "zextents[0]" does not enclose all regions')
                if self.zextents[1] is None:
                    self.zextents[1] = self.zpad[1]*self.zmax
                    if self.zmax > self.zextents[1]:
                        raise ValueError('User specified "zextents[1]" does not enclose all regions')
                # Create boundary region
                vac_dx = self.region_info[self.boundary_reg]['dx']
                if vac_dx is None:
                    raise ValueError('Resolution for region "vacuum" not defined')
                vac_contour = numpy.asarray([
                    [0.0,           self.zextents[0]],
                    [self.rextent, self.zextents[0]],
                    [self.rextent, self.zextents[1]],
                    [0.0,           self.zextents[1]]
                ])
                self.regions.append(Region(vac_contour,vac_dx,vac_dx,id=self.region_info[self.boundary_reg]['id']))
                self.region_info[self.boundary_reg]['count'] += 1
        # Check for undefined regions
        for key in self.region_info:
            if self.region_info[key]['count'] == 0:
                raise KeyError('Region "{0}" defined but never created'.format(key))
        # Re-index regions
        reg_reorder = [-1 for _ in self.region_info]
        reg_reorder[0] = 1
        reg_id = 1
        for reg_type in ('boundary', 'vacuum','conductor','coil'):
            for key in self.region_info:
                if self.region_info[key]['type'] == reg_type:
                    reg_id += 1
                    reg_reorder[self.region_info[key]['id']-1] = reg_id
                    self.region_info[key]['id'] = reg_id
        for region in self.regions:
            region._id = reg_reorder[region._id-1]
        for point_def in self._extra_reg_defs:
            point_def[2] = reg_reorder[point_def[2]-1]
        # Generate mesh
        self.mesh = Mesh(self.regions,debug=debug,extra_reg_defs=self._extra_reg_defs,merge_thresh=merge_thresh)
        if setup_only:
            return None, None, None
        return self.mesh.get_mesh()
    
    def save_json(self,filename):
        '''! Create a JSON file containing a description of the mesh 

        @param filename Path to create JSON file
        '''
        output_dict = {
            'rextent': self.rextent,
            'zextents': self.zextents,
            'rpad': self.rpad,
            'zpad': self.zpad,
            'rmax': self.rmax,
            'zmin': self.zmin,
            'zmax': self.zmax,
            'boundary_reg': self.boundary_reg,
            'reg_type_counts': self.reg_type_counts,
            'regions': [],
            'region_info': self.region_info

        }
        for region in self.regions:
            output_dict['regions'].append(region.get_dict())
        with open(filename, 'w+') as fid:
            fid.write(json.dumps(output_dict))


def save_gs_mesh(pts,tris,regions,coil_dict,cond_dict,filename,use_hdf5=True):
        '''! Save G-S mesh to file in HDF5 format

        @param pts[np,2] Vertex list
        @param tris[nc,3] Cell list
        @param regions[nc] Region list
        @param coil_dict Coil region dictionary
        @param cond_dict Conducting region dictionary
        @param filename Path to create HDF5 mesh file
        '''
        if use_hdf5:
            import h5py
            coil_json = json.dumps(coil_dict)
            cond_json = json.dumps(cond_dict)
            with h5py.File(filename, 'w') as h5_file:
                h5_file.create_dataset('mesh/r', data=pts, dtype='f8')
                h5_file.create_dataset('mesh/lc', data=tris, dtype='i4')
                h5_file.create_dataset('mesh/reg', data=regions, dtype='i4')
                string_datatype = h5py.string_dtype('ascii')
                h5_file.create_dataset('mesh/coil_dict', data=coil_json, dtype=string_datatype)
                h5_file.create_dataset('mesh/cond_dict', data=cond_json, dtype=string_datatype)
        else:
            with open(filename, 'w+') as fid:
                fid.write(json.dumps({
                    'mesh': {
                        'r': pts.tolist(),
                        'lc': tris.tolist(),
                        'reg': regions.tolist(),
                        'coil_dict': coil_dict,
                        'cond_dict': cond_dict
                    }
                }))


def load_gs_mesh(filename,use_hdf5=True):
        '''! Load G-S mesh to file in HDF5 format

        @param filename Path to HDF5 mesh file
        @result pts[np,2], tris[nc,3], regions[nc], coil_dict, cond_dict
        '''
        if use_hdf5:
            import h5py
            with h5py.File(filename, 'r') as h5_file:
                pts = numpy.asarray(h5_file['mesh/r'])
                tris = numpy.asarray(h5_file['mesh/lc'])
                regions = numpy.asarray(h5_file['mesh/reg'])
                coil_dict = json.loads(h5_file['mesh/coil_dict'][()])
                cond_dict = json.loads(h5_file['mesh/cond_dict'][()])
        else:
            with open(filename, 'r') as fid:
                input_dict = json.load(fid)
            pts = numpy.asarray(input_dict['mesh']['r'])
            tris = numpy.asarray(input_dict['mesh']['lc'])
            regions = numpy.asarray(input_dict['mesh']['reg'])
            coil_dict = input_dict['mesh']['coil_dict']
            cond_dict = input_dict['mesh']['cond_dict']
        return pts, tris, regions, coil_dict, cond_dict


class Mesh:
    '''! Mesh builder class for triangle library'''
    def __init__(self,region_list,merge_thresh=1.E-4,debug=False,extra_reg_defs=[]):
        '''! Initialize Mesh builder object

        @param region_list List of @ref oftpy.Region objects that define mesh
        @param merge_thresh Distance threshold for merging nearby points
        '''
        self._merge_thresh = merge_thresh
        self._unique_points = []
        self._reg_seg_map = []
        self._segments = []
        print('Assembling regions:')
        # Build list of unique points
        for ireg, region in enumerate(region_list):
            local_seg_map = []
            reg_pt_map = [i+len(self._unique_points) for i in range(region._points.shape[0])]
            ilocal = len(self._unique_points)-1
            for tmp_pts in region._segments:
                tmp_pt_map = [reg_pt_map[i] for i in tmp_pts]
                for ipt, reg_id in enumerate(tmp_pts):
                    reg_pt = region._points[reg_id,:]
                    if reg_pt_map[reg_id] > len(self._unique_points)-1:
                        for jpt, unique_pt in enumerate(self._unique_points):
                            if numpy.linalg.norm(reg_pt-unique_pt) < merge_thresh:
                                if debug:
                                    print('  Merging points:',ireg,reg_pt_map[reg_id],jpt,reg_pt,unique_pt)
                                tmp_pt_map[ipt] = jpt
                                reg_pt_map[reg_id] = jpt
                                break
                        else:
                            ilocal += 1
                            tmp_pt_map[ipt] = ilocal
                            reg_pt_map[reg_id] = ilocal
                            self._unique_points.append(reg_pt)
                for iseg, segment in enumerate(self._segments):
                    nOverlap = 0
                    for test_pt in segment[0]:
                        if test_pt in tmp_pt_map:
                            nOverlap += 1
                    if (nOverlap > 1) and (len(tmp_pts) == len(segment[0])): # Full overlap
                        # Look forward
                        for i,test_pt in enumerate(segment[0]):
                            if tmp_pt_map[i] != test_pt:
                                break
                        else: # Matched segment
                            if debug:
                                print('  Merging curve segments:',ireg,iseg)
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            local_seg_map.append(iseg)
                            break
                        # Look backward
                        for i,test_pt in enumerate(segment[0]):
                            if tmp_pt_map[-i-1] != test_pt:
                                break
                        else: # Matched segment
                            if debug:
                                print('  Merging curve segments:',ireg,iseg)
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            local_seg_map.append(-iseg)
                            break
                    elif  (nOverlap > 1): # Partial match
                        if debug:
                            print('  Merging partially overlapping curve segments:',ireg,iseg)
                        if len(tmp_pts) < len(segment[0]):
                            overlap_start = len(segment[0])
                            overlap_end = -1
                            first_pt = -1
                            for i, test_pt in enumerate(segment[0]):
                                if test_pt in tmp_pt_map:
                                    overlap_start = min(overlap_start,i)
                                    overlap_end = max(overlap_end,i)
                                    if first_pt < 0 :
                                        first_pt = tmp_pt_map.index(test_pt)
                            segment_split = [iseg, -1, -1]
                            if overlap_start > 0:
                                self._segments.append([segment[0][:overlap_start+1], segment[1], segment[2]])
                                segment_split[1] = len(self._segments)-1
                            if overlap_end < len(segment[0])-1:
                                self._segments.append([segment[0][overlap_end:], segment[1], segment[2]])
                                segment_split[2] = len(self._segments)-1
                            segment[0] = segment[0][overlap_start:overlap_end+1]
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            #
                            for reg_seg_map in self._reg_seg_map:
                                try:
                                    ifound = reg_seg_map.index(iseg)
                                    if segment_split[1] >= 0:
                                        reg_seg_map.insert(ifound,segment_split[1])
                                        ifound += 1
                                    if segment_split[2] >= 0:
                                        reg_seg_map.insert(ifound+1,segment_split[2])
                                except ValueError:
                                    ifound = -1
                                    pass
                                try:
                                    ifound = reg_seg_map.index(-iseg)
                                    if ifound >= 0:
                                        if segment_split[2] >= 0:
                                            reg_seg_map.insert(ifound,-segment_split[2])
                                            ifound += 1
                                        if segment_split[1] >= 0:
                                            reg_seg_map.insert(ifound+1,-segment_split[1])
                                except ValueError:
                                    pass
                            #
                            if first_pt == overlap_start:
                                local_seg_map.append(iseg)
                                break
                            else:
                                local_seg_map.append(-iseg)
                                break
                        else:
                            overlap_start = len(tmp_pts)
                            overlap_end = -1
                            first_pt = -1
                            for i, test_pt in enumerate(tmp_pt_map):
                                if test_pt in segment[0]:
                                    overlap_start = min(overlap_start,i)
                                    overlap_end = max(overlap_end,i)
                                    if first_pt < 0 :
                                        first_pt = segment[0].index(test_pt)
                            if overlap_start > 0:
                                self._segments.append([tmp_pt_map[:overlap_start+1], region._dx_curve, region._small_thresh])
                            if overlap_end < len(segment[0])-1:
                                self._segments.append([tmp_pt_map[overlap_end:], region._dx_curve, region._small_thresh])
                            segment[1] = min(segment[1],region._dx_curve)
                            segment[2] = min(segment[2],region._small_thresh)
                            #
                            if first_pt == overlap_start:
                                local_seg_map.append(iseg)
                                break
                            else:
                                local_seg_map.append(-iseg)
                                break
                else:
                    self._segments.append([tmp_pt_map, region._dx_curve, region._small_thresh])
                    local_seg_map.append(len(self._segments)-1)
            self._reg_seg_map.append(local_seg_map)
        # Resample segments
        pts_out = self._unique_points
        self._unique_points = numpy.array(self._unique_points)
        self._resampled_segments = []
        for segment in self._segments:
            dx = segment[1]
            pts_tmp = self._unique_points[segment[0],:]
            seg_tmp = [segment[0][0],]
            dl = numpy.zeros(len(segment[0]))
            for i in range(len(segment[0])-1):
                dl[i+1] = dl[i] + numpy.linalg.norm(pts_tmp[i+1,:]-pts_tmp[i,:])
            if int(dl[-1]/dx) >= 2:
                for dl_samp in numpy.linspace(0.0,dl[-1],int(dl[-1]/dx)+1)[1:-1]:
                    pts_out.append([
                        numpy.interp(dl_samp,dl,pts_tmp[:,0]),
                        numpy.interp(dl_samp,dl,pts_tmp[:,1])
                    ])
                    seg_tmp.append(len(pts_out)-1)
            elif dl[-1] < segment[2]:
                print("  Warning: small feature (dl={0:.2E}) detected at point {1} ({2}, {3})".format(dl[-1], segment[0][0], *pts_tmp[0,:]))
            seg_tmp.append(segment[0][-1])
            self._resampled_segments.append(seg_tmp)
        # Reindex points
        reindex = -numpy.ones((len(pts_out),),dtype=numpy.int32)
        pt_count = 0
        for segment in self._resampled_segments:
            for i in range(len(segment)-1):
                if reindex[segment[i]] < 0:
                    pt_count += 1
                    reindex[segment[i]] = pt_count-1
                if reindex[segment[i+1]] < 0:
                    pt_count += 1
                    reindex[segment[i+1]] = pt_count-1
                segment[i]=reindex[segment[i]]
            segment[-1]=reindex[segment[-1]]
        self._resampled_points = numpy.zeros((pt_count,2))
        for i in range(len(pts_out)):
            if reindex[i] >= 0:
                self._resampled_points[reindex[i],:] = pts_out[i]
        #
        print('  # of unique points    = {0}'.format(self._resampled_points.shape[0]))
        print('  # of unique segments  = {0}'.format(len(self._segments)))
        # Build list of region definitions
        self._reg_defs = extra_reg_defs
        for ireg, region in enumerate(region_list):
            imin = -1
            dmin = 1.E99
            dmax = -1.0
            reg_points = []
            for segment in self._reg_seg_map[ireg]:
                if segment < 0:
                    tmp_pt_map = [entry for entry in reversed(self._resampled_segments[-segment])]
                else:
                    tmp_pt_map = self._resampled_segments[segment]
                for ipt, unique_id in enumerate(tmp_pt_map):
                    if ipt == len(tmp_pt_map)-1:
                        continue
                    reg_pt = self._resampled_points[unique_id]
                    reg_points.append(reg_pt)
                    dmin_loc = 1.E99
                    for jpt, unique_pt in enumerate(self._resampled_points):
                        if unique_id == jpt:
                            continue
                        dmin_loc = min(dmin_loc,numpy.linalg.norm(reg_pt-unique_pt))
                    dmin = min(dmin,dmin_loc)
                    if (dmin_loc > dmax):
                        dmax = dmin_loc
                        imin = len(reg_points)-1
            region._resampled_points = numpy.asarray(reg_points)
            in_pt = region.get_in_point(imin,dmin)
            self._reg_defs.append([in_pt[0], in_pt[1], region._id, region._dx_vol*region._dx_vol/2.0])
    
    def get_mesh(self):
        '''! Generate mesh using triangle

        @result pts[np,2], tris[nc,3], regions[nc] 
        '''
        resampled_flat = []
        for segment in self._resampled_segments:
            resampled_flat += [[segment[i], segment[i+1]] for i in range(len(segment)-1)]
        alpha = dict(vertices=self._resampled_points, segments=resampled_flat, regions=self._reg_defs)
        try:
            import triangle as tr
        except ImportError:
            print('Meshing requires "triangle" python library')
            return None
        print('Generating mesh:')
        beta = tr.triangulate(alpha,'pqaA')
        regions = beta['triangle_attributes'].astype('int32').ravel()
        print('  # of points  = {0}'.format(len(beta['vertices'])))
        print('  # of cells   = {0}'.format(len(beta['triangles'])))
        print('  # of regions = {0}'.format(regions.max()))
        return beta['vertices'], beta['triangles'], regions


class Region:
    '''! Region class for @ref oftpy.Mesh class'''
    def __init__(self,points,dx=None,dx_curve=None,angle_tol=30.0,sliver_tol=120.0,small_thresh=None,id=0,load_dict=None):
        '''! Create Region object from a closed bounding curve

        @param points List of points forming region boundary [:,2]
        @param dx Target mesh resolution inside region
        @param angle_tol Corner tolerance used when resampling curve at desired resolution
        @param sliver_tol Tolerance used for "sliver" region warnings
        @param small_thresh Tolerance used for "small" curve warnings (default: dx/2)
        @param id Region id number (default: 0)
        '''
        if load_dict is not None:
            self._points = load_dict['points']
            self._segments = load_dict['segments']
            self._id = load_dict['id']
            self._dx_curve = load_dict['dx_curve']
            self._dx_vol = load_dict['dx_vol']
        else:
            if numpy.any(points[:,0] < 0):
                raise ValueError("Point with negative radial position detected!")
            self._points = points
            self._id = id
            if dx is None:
                raise ValueError('No target mesh size set')
            else:
                self._dx_vol = dx
                if dx_curve is not None:
                    self._dx_curve = dx_curve
                else:
                    self._dx_curve = dx
            if numpy.linalg.norm(self._points[0,:]-self._points[-1,:]) < dx/1.E3:
                self._points = self._points[:-1,:]
            nv = self._points.shape[0]
            # Detect corner points and split contour
            keep_tol = numpy.cos(numpy.pi*angle_tol/180)
            sliver_tol = numpy.cos(numpy.pi*sliver_tol/180)
            keep_points = [0]
            for i in range(nv):
                if (i == nv-1):
                    tang = self._points[0,:] - self._points[i,:]
                else:
                    tang = self._points[i+1,:] - self._points[i,:]
                tang_norm = numpy.linalg.norm(tang)
                if tang_norm < self._dx_curve/1.E3:
                    print("  Warning: repeated points detected at point {0} ({1}, {2})".format(i, *self._points[i,:]))
                    continue
                tang /= tang_norm
                if i > 0:
                    angle = numpy.dot(tang,tangp)
                    if angle < keep_tol:
                        keep_points.append(i)
                        if angle < sliver_tol:
                            print("  Warning: sliver (angle={0:.1F}) detected at point {1} ({2}, {3})".format(180.-numpy.arccos(angle)*180./numpy.pi, i, *self._points[i,:]))
                tangp = tang
            keep_points.append(nv+1)
            # Index segments
            k=1
            self._segments = []
            seg_tmp = []
            for i in range(nv):
                if i >= keep_points[k]:
                    self._segments.append(seg_tmp + [i,])
                    seg_tmp = [i,]
                    k += 1
                else:
                    seg_tmp.append(i)
            self._segments.append(seg_tmp + [0,])
        # Get small feature threshold
        if small_thresh is None:
            small_thresh = self._dx_curve/2.0
        self._small_thresh = small_thresh
        self._resampled_points = None
    
    def get_resampled_points(self):
        '''! Get resampled points for bounding curve

        @result Point list [np,2]
        '''
        # Resample curves to approximate target dx
        if self._resampled_points is None:
            pts_out = []
            for segment in self._segments:
                pts_tmp = self._points[segment,:]
                dl = numpy.zeros(len(segment))
                for i in range(len(segment)-1):
                    dl[i+1] = dl[i] + numpy.linalg.norm(pts_tmp[i+1,:]-pts_tmp[i,:])
                pts_out.append(pts_tmp[0,:])
                if int(dl[-1]/self._dx_curve) >= 2:
                    for dl_samp in numpy.linspace(0.0,dl[-1],int(dl[-1]/self._dx_curve)+1)[1:-1]:
                        pts_out.append([
                            numpy.interp(dl_samp,dl,pts_tmp[:,0]),
                            numpy.interp(dl_samp,dl,pts_tmp[:,1])
                        ])
                elif dl[-1] < self._small_thresh:
                    print("  Warning: small feature (dl={0:.2E}) detected at point {1} ({2}, {3})".format(dl[-1], segment[0], *pts_tmp[0,:]))
            self._resampled_points = numpy.asarray(pts_out)
        return self._resampled_points

    def check_in_poly(self,pt):
        '''! Check if point is inside region (polygon-based approach)

        @param pt Point to check [2]
        @result Is pt inside the region?
        '''
        ncuts = 0
        for j in range(self._resampled_points.shape[0]-1):
            if (pt[1]-self._resampled_points[j,1])*(pt[1]-self._resampled_points[j+1,1]) <= 0.0:
                if (pt[1]-self._resampled_points[j,1]) == 0.0:
                    continue
                xInter = (self._resampled_points[j+1,0]-self._resampled_points[j,0])* \
                    (pt[1]-self._resampled_points[j,1])/(self._resampled_points[j+1,1]-self._resampled_points[j,1]) \
                    + self._resampled_points[j,0]
                if pt[0] <= xInter:
                    ncuts +=1
        if (pt[1]-self._resampled_points[-1,1])*(pt[1]-self._resampled_points[0,1]) <= 0.0:
            if (pt[1]-self._resampled_points[-1,1]) == 0.0:
                return (ncuts % 2 == 1)
            xInter = (self._resampled_points[0,0]-self._resampled_points[-1,0])* \
                (pt[1]-self._resampled_points[-1,1])/(self._resampled_points[0,1]-self._resampled_points[-1,1]) \
                + self._resampled_points[-1,0]
            if pt[0] <= xInter:
                ncuts +=1
        return (ncuts % 2 == 1)

    def get_in_point(self,i,dx):
        '''! Get a suitable point for defining the "inside" of the region for triangle 

        @param i Index of desired nearest boundary point
        @param dx Offset distance from bounding curve
        @result Point used to define region [2]
        '''
        dx = min(dx,self._dx_curve)/4.0
        if i==self._resampled_points.shape[0]-1:
            that2 = self._resampled_points[0,:] - self._resampled_points[i,:]
        else:
            that2 = self._resampled_points[i+1,:] - self._resampled_points[i,:]
        that1 = self._resampled_points[i,:] - self._resampled_points[i-1,:]
        that1 /= numpy.linalg.norm(that1)
        that2 /= numpy.linalg.norm(that2)
        # nhat = dx*numpy.r_[-that1[1], that1[0]]
        # dx2 = dx-(nhat[1]*that2[0]-nhat[0]*that2[1])
        # nhat += dx2*numpy.r_[-that2[1],that2[0]]
        nhat = dx*(numpy.r_[-that1[1], that1[0]] + numpy.r_[-that2[1], that2[0]])
        pt_out = self._resampled_points[i,:] + nhat
        # Check if inside contour
        if not self.check_in_poly(pt_out):
            pt_out = self._resampled_points[i,:] - nhat
        return pt_out
    
    def get_segments(self):
        segments = []
        for i in range(len(self._segments)):
            segments.append(self._points[self._segments[i],:])
        return segments
    
    def plot_segments(self,fig,ax):
        '''! Plot boundary curve
        
        @param fig Figure to add curves to
        @param ax Axis to add curves to
        '''
        for i in range(len(self._segments)):
            ax.plot(self._points[self._segments[i],0],self._points[self._segments[i],1])
    
    def get_json(self):
        return {
            'points': self._points,
            'segments': self._segments,
            'id': self._id,
            'dx_curve': self._dx_curve,
            'dx_vol': self._dx_vol
        }
