#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Core definitions for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import collections
import ctypes
import numpy
from ._interface import *


def tokamaker_default_settings(oft_env):
    '''! Initialize settings object with default values

    @param oft_env Current runtime environment
    @result tokamaker_settings_struct object
    '''
    settings = tokamaker_settings_struct()
    settings.pm = True
    settings.free_boundary = True
    settings.has_plasma = True
    settings.limited_only = False
    settings.dipole_mode = False
    settings.maxits = 40
    settings.mode = 1
    settings.urf = 0.2
    settings.nl_tol = 1.E-6
    settings.rmin = 0.0
    settings.lim_zmax = 1.E99
    settings.limiter_file = oft_env.path2c('none')
    return settings


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


class TokaMaker():
    '''! TokaMaker G-S solver class'''
    def __init__(self,OFT_env):
        '''! Initialize TokaMaker object

        @param OFT_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        '''
        # Create OFT execution environment
        self._oft_env = OFT_env
        ## Internal Grad-Shafranov object (@ref psi_grad_shaf.gs_eq "gs_eq")
        self._tMaker_ptr = c_void_p()
        ## Internal mesh object
        self._mesh_ptr = c_void_p()
        ## General settings object
        self.settings = tokamaker_default_settings(self._oft_env)
        ## Conductor definition dictionary
        self._cond_dict = {}
        ## Vacuum definition dictionary
        self._vac_dict = {}
        ## Coil definition dictionary
        self._coil_dict = {}
        ## Coil set definitions, including sub-coils
        self.coil_sets = {}
        ## Virtual coils, if present (currently only `'#VSC'`)
        self._virtual_coils = {'#VSC': -1}
        ## Coil set names in order of id number
        self.coil_set_names = []
        ## Distribution coils, only (currently) saved for plotting utility
        self.dist_coils = {}
        ## Vacuum F value
        self._F0 = 0.0
        ## Plasma current target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._Ip_target=c_double(self._oft_env.float_disable_flag)
        ## Plasma current target ratio I_p(FF') / I_p(P') (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._Ip_ratio_target=c_double(self._oft_env.float_disable_flag)
        ## Axis pressure target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._pax_target=c_double(self._oft_env.float_disable_flag)
        ## Stored energy target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._estore_target=c_double(self._oft_env.float_disable_flag)
        ## R0 target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._R0_target=c_double(self._oft_env.float_disable_flag)
        ## V0 target value (use @ref TokaMaker.TokaMaker.set_targets "set_targets")
        self._V0_target=c_double(self._oft_env.float_disable_flag)
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
    
    def __del__(self):
        '''! Free Fortran-side objects by calling `reset()` before object is deleted or GC'd'''
        self.reset()

    def reset(self):
        '''! Reset G-S object to enable loading a new mesh and coil configuration'''
        if not self._tMaker_ptr:
            return # Nothing to do
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_destroy(self._tMaker_ptr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self.nregs = -1
        self.np = -1
        # Reset defaults
        self._tMaker_ptr = c_void_p()
        self._mesh_ptr = c_void_p()
        self.settings = tokamaker_default_settings(self._oft_env)
        self._cond_dict = {}
        self._vac_dict = {}
        self._coil_dict = {}
        self.coil_sets = {}
        self._virtual_coils = {}
        self._F0 = 0.0
        self._Ip_target=c_double(self._oft_env.float_disable_flag)
        self._Ip_ratio_target=c_double(self._oft_env.float_disable_flag)
        self._pax_target=c_double(self._oft_env.float_disable_flag)
        self._estore_target=c_double(self._oft_env.float_disable_flag)
        self._R0_target=c_double(self._oft_env.float_disable_flag)
        self._V0_target=c_double(self._oft_env.float_disable_flag)
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
            self._oft_env.oft_in_groups['mesh_options'] = {'cad_type': "0"}
            self._oft_env.oft_in_groups['native_mesh_options'] = {'filename': '"{0}"'.format(mesh_file)}
            self._oft_env.update_oft_in()
            oft_setup_smesh(ndim,ndim,rfake,ndim,ndim,lcfake,regfake,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
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
                if reg.min() <= 0:
                    raise ValueError('Invalid "reg" array, values must be >= 0')
                reg = numpy.ascontiguousarray(reg, dtype=numpy.int32)
            oft_setup_smesh(ndim,np,r,npc,nc,lc+1,reg,ctypes.byref(nregs),ctypes.byref(self._mesh_ptr))
        else:
            raise ValueError('Mesh filename (native format) or mesh values required')
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_alloc(ctypes.byref(self._tMaker_ptr),self._mesh_ptr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self.update_settings()
        self.nregs = nregs.value
    
    def setup_regions(self,cond_dict={},coil_dict={}):
        '''! Define mesh regions (coils and conductors)

        @param cond_dict Dictionary specifying conducting regions
        '''
        xpoint_mask = numpy.zeros((self.nregs,),dtype=numpy.int32)
        xpoint_mask[0] = 1
        eta_vals = -2.0*numpy.ones((self.nregs,),dtype=numpy.float64)
        eta_vals[0] = -1.0
        contig_flag = numpy.ones((self.nregs,),dtype=numpy.int32)
        # Process conductors and vacuum regions
        self._vac_dict = {}
        for key in cond_dict:
            if 'vac_id' in cond_dict[key]:
                self._vac_dict[key] = cond_dict[key]
            else:
                eta_vals[cond_dict[key]['reg_id']-1] = cond_dict[key]['eta']/mu0
                if cond_dict[key].get('noncontinuous',False):
                    contig_flag[cond_dict[key]['reg_id']-1] = 0
            xpoint_mask[cond_dict[key]['reg_id']-1] = int(cond_dict[key].get('allow_xpoints',False))
            if cond_dict[key].get('inner_limiter',False):
                contig_flag[cond_dict[key]['reg_id']-1] = -1
        # Remove vacuum regions
        for key in self._vac_dict:
            del cond_dict[key]
        self._cond_dict = cond_dict
        # Process coils
        nCoils = 0
        self.coil_sets = {}
        for key in coil_dict:
            xpoint_mask[coil_dict[key]['reg_id']-1] = int(coil_dict[key].get('allow_xpoints',False))
            eta_vals[coil_dict[key]['reg_id']-1] = -1.0
            coil_set = coil_dict[key].get('coil_set',key)
            if coil_set not in self.coil_sets:
                self.coil_sets[coil_set] = {
                    'id': nCoils,
                    'net_turns': 0.0,
                    'sub_coils': []
                }
                nCoils += 1
            self.coil_sets[coil_set]['sub_coils'].append(coil_dict[key])
            self.coil_sets[coil_set]['net_turns'] += coil_dict[key].get('nturns',1.0)
        self._coil_dict = coil_dict
        # Mark vacuum regions
        self.nvac = 0
        for i in range(self.nregs):
            if eta_vals[i] < -1.5:
                eta_vals[i] = 1.E10
                self.nvac += 1 
        self.coil_set_names = ['' for _ in range(nCoils)]
        coil_nturns = numpy.zeros((nCoils, self.nregs))
        for key in self.coil_sets:
            self.coil_set_names[self.coil_sets[key]['id']] = key
            for sub_coil in self.coil_sets[key]['sub_coils']:
                coil_nturns[self.coil_sets[key]['id'],sub_coil['reg_id']-1] = sub_coil.get('nturns',1.0)
        cstring = self._oft_env.path2c('none')
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_setup_regions(self._tMaker_ptr,cstring,eta_vals,contig_flag,xpoint_mask,coil_nturns,nCoils,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
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
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_setup(self._tMaker_ptr,order,full_domain,ctypes.byref(ncoils),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        ## Number of coils in mesh
        self.ncoils = ncoils.value
        ## Isoflux constraint points (use @ref TokaMaker.TokaMaker.set_isoflux "set_isoflux")
        self._isoflux_targets = None
        ## Flux constraint points (use @ref TokaMaker.TokaMaker.set_isoflux "set_flux")
        self._flux_targets = None
        ## Saddle constraint points (use @ref TokaMaker.TokaMaker.set_saddles "set_saddles")
        self._saddle_targets = None
        # Get references to internal variables
        o_loc = c_double_ptr()
        lim_loc = c_double_ptr()
        x_loc = c_double_ptr()
        div_flag_loc = c_bool_ptr()
        bounds_loc = c_double_ptr()
        alam_loc = c_double_ptr()
        pnorm_loc = c_double_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_refs(self._tMaker_ptr,ctypes.byref(o_loc),ctypes.byref(lim_loc),ctypes.byref(x_loc),ctypes.byref(div_flag_loc),
                    ctypes.byref(bounds_loc),ctypes.byref(alam_loc),ctypes.byref(pnorm_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
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
        nloops = c_int()
        loop_ptr = c_int_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_limiter(self._tMaker_ptr,ctypes.byref(npts),ctypes.byref(r_loc),ctypes.byref(nloops),ctypes.byref(loop_ptr),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        loop_ptr = numpy.ctypeslib.as_array(loop_ptr,shape=(nloops.value+1,))
        self.lim_pts = numpy.ctypeslib.as_array(r_loc,shape=(npts.value, 2))
        self.lim_contours = []
        for i in range(nloops.value):
            lim_contour = numpy.vstack((self.lim_pts[loop_ptr[i]-1:loop_ptr[i+1]-1,:],self.lim_pts[loop_ptr[i]-1,:]))
            self.lim_contours.append(lim_contour)
        self.lim_contour = numpy.zeros((0,2))
        for lim_countour in self.lim_contours:
            if lim_countour.shape[0] > self.lim_contour.shape[0]:
                self.lim_contour = lim_countour
        # Get plotting mesh
        np_loc = c_int()
        nc_loc = c_int()
        r_loc = c_double_ptr()
        lc_loc = c_int_ptr()
        reg_loc = c_int_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_mesh(self._tMaker_ptr,ctypes.byref(np_loc),ctypes.byref(r_loc),ctypes.byref(nc_loc),ctypes.byref(lc_loc),ctypes.byref(reg_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
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
        
    def abspsi_to_normalized(self,psi_in):
        r'''! Convert unnormalized \f$ \psi \f$ values to normalized \f$ \hat{\psi} \f$ values
        
        @param psi_in Input \f$ \psi \f$ values
        @returns Normalized \f$ \hat{\psi} \f$ values
        '''
        if self.psi_convention == 0:
            return (psi_in-self.psi_bounds[1])/(self.psi_bounds[0]-self.psi_bounds[1])
        else:
            return (psi_in-self.psi_bounds[0])/(self.psi_bounds[1]-self.psi_bounds[0])
    
    def psinorm_to_absolute(self,psi_in):
        r'''! Convert normalized \f$ \hat{\psi} \f$ values to unnormalized values \f$ \psi \f$
        
        @param psi_in Input \f$ \hat{\psi} \f$ values
        @returns Unnormalized \f$ \psi \f$ values
        '''
        if self.psi_convention == 0:
            return psi_in*(self.psi_bounds[0]-self.psi_bounds[1]) + self.psi_bounds[1]
        else:
            return psi_in*(self.psi_bounds[1]-self.psi_bounds[0]) + self.psi_bounds[0]
        
    def coil_reg_term(self,coffs,target=0.0,weight=1.0):
        r'''! Define coil current regularization term for the form \f$ target = \Sigma_i \alpha_i I_i \f$
        to be used in @ref set_coil_reg.

        @param coffs Dictionary of coefficients \f$ \alpha \f$ (zero for unspecified coils)
        @param target Regularization target (default: 0.0)
        @param weight Weight for regularization term (default: 1.0)
        '''
        coil_reg_term = collections.namedtuple('coil_reg_term', ['coffs', 'target', 'weight'])
        for key in coffs:
            if (key not in self.coil_sets) and (key not in self._virtual_coils):
                raise KeyError('Unknown coil "{0}"'.format(key))
        return coil_reg_term(coffs,target,weight)
    
    def set_coil_reg(self,reg_mat=None,reg_targets=None,reg_weights=None,reg_terms=None):
        '''! Set regularization matrix for coil currents when isoflux and/or saddle constraints are used

        Can be used to enforce "soft" constraints on coil currents. For hard constraints see
        @ref TokaMaker.TokaMaker.set_coil_bounds "set_coil_bounds".

        @param reg_mat Regularization matrix [nregularize,ncoils+1]
        @param reg_targets Regularization targets [nregularize] (default: 0)
        @param reg_weights Weights for regularization terms [nregularize] (default: 1)
        @param reg_terms List of regularization terms created with @ref coil_reg_term
        '''
        if reg_terms is not None:
            if reg_mat is not None:
                raise ValueError('"reg_terms" and "reg_mat" cannot be specified simultaneously')
            if reg_targets is not None:
                raise ValueError('"reg_terms" and "reg_targets" cannot be specified simultaneously')
            if reg_weights is not None:
                raise ValueError('"reg_terms" and "reg_weights" cannot be specified simultaneously')
            nregularize = len(reg_terms)
            reg_mat = numpy.zeros((self.ncoils+1,nregularize), dtype=numpy.float64)
            reg_targets = numpy.zeros((nregularize,), dtype=numpy.float64)
            reg_weights = numpy.ones((nregularize,), dtype=numpy.float64)
            for i, reg_term in enumerate(reg_terms):
                reg_targets[i] = reg_term.target
                reg_weights[i] = reg_term.weight
                for key, value in reg_term.coffs.items():
                    if key in self.coil_sets:
                        reg_mat[self.coil_sets[key]['id'],i] = value
                    elif key in self._virtual_coils:
                        reg_mat[self._virtual_coils[key],i] = value
                    else:
                        raise KeyError('Unknown coil "{0}"'.format(key))
        elif reg_mat is not None:
            if reg_mat.shape[1] != self.ncoils+1:
                raise IndexError('Incorrect shape of "reg_mat", should be [nregularize,ncoils+1]')
            nregularize = reg_mat.shape[0]
            reg_mat = numpy.ascontiguousarray(reg_mat.transpose(), dtype=numpy.float64)
            if reg_targets is None:
                reg_targets = numpy.zeros((nregularize,), dtype=numpy.float64)
            if reg_weights is None:
                reg_weights = numpy.ones((nregularize,), dtype=numpy.float64)
            if reg_targets.shape[0] != nregularize:
                raise IndexError('Incorrect shape of "reg_targets", should be [nregularize]')
            if reg_weights.shape[0] != nregularize:
                raise IndexError('Incorrect shape of "reg_weights", should be [nregularize]')
        else:
            raise ValueError('Either "reg_terms" or "reg_mat" is required')
        # Ensure VSC is constrained
        if (self._virtual_coils.get('#VSC',-1) < 0) and ((abs(reg_mat[-1,:])).max() < 1.E-8):
            new_row = numpy.zeros((self.ncoils+1,), dtype=numpy.float64)
            new_row[-1] = 1.0
            reg_mat = numpy.hstack((reg_mat,new_row.reshape([self.ncoils+1,1])))
            reg_targets = numpy.append(reg_targets, 0.0)
            reg_weights = numpy.append(reg_weights, 1.0)
            nregularize += 1

        reg_targets = numpy.ascontiguousarray(reg_targets, dtype=numpy.float64)
        reg_weights = numpy.ascontiguousarray(reg_weights, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_regmat(self._tMaker_ptr,nregularize,reg_mat,reg_targets,reg_weights,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def set_coil_bounds(self,coil_bounds=None):
        '''! Set hard constraints on coil currents

        Can be used with or without regularization terms (see
        @ref TokaMaker.TokaMaker.set_coil_reg "set_coil_reg").

        @param coil_bounds Minimum and maximum allowable coil currents [ncoils+1,2]
        '''
        bounds_array = numpy.zeros((self.ncoils+1,2), dtype=numpy.float64)
        bounds_array[:,0] = -1.E98; bounds_array[:,1] = 1.E98
        if coil_bounds is not None:
            for coil_key, coil_bound in coil_bounds.items():
                if coil_key in self.coil_sets:
                    bounds_array[self.coil_sets[coil_key]['id'],:] = coil_bound
                elif coil_key in self._virtual_coils:
                    bounds_array[self._virtual_coils[coil_key],:] = coil_bound
                else:
                    raise KeyError('Unknown coil "{0}"'.format(coil_key))
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_bounds(self._tMaker_ptr,bounds_array,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def set_coil_vsc(self,coil_gains):
        '''! Define a vertical stability coil set from one or more coils

        @param coil_gains Gains for each coil (absolute scale is arbitrary)
        '''
        gains_array = numpy.zeros((self.ncoils,), dtype=numpy.float64)
        for coil_key, coil_gain in coil_gains.items():
            gains_array[self.coil_sets[coil_key]['id']] = coil_gain
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_vsc(self._tMaker_ptr,gains_array,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def init_psi(self, r0=-1.0, z0=0.0, a=0.0, kappa=0.0, delta=0.0, curr_source=None):
        r'''! Initialize \f$\psi\f$ using uniform current distributions

        If r0>0 then a uniform current density inside a surface bounded by
        a curve of the form defined in @ref oftpy.create_isoflux is used.
        If r0<0 then a uniform current density over the entire plasma region is used.

        @param r0 Major radial position for flux surface-based approach
        @param z0 Vertical position for flux surface-based approach
        @param a Minor radius for flux surface-based approach
        @param kappa Elongation for flux surface-based approach
        @param delta Triangularity for flux surface-based approach
        @param curr_source Current source for arbitrary current distribution
        '''
        curr_ptr = None
        if curr_source is not None:
            if curr_source.shape[0] != self.np:
                raise IndexError('Incorrect shape of "curr_source", should be [np]')
            curr_source = numpy.ascontiguousarray(curr_source, dtype=numpy.float64)
            curr_ptr = curr_source.ctypes.data_as(c_double_ptr)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_init_psi(self._tMaker_ptr,c_double(r0),c_double(z0),c_double(a),c_double(kappa),c_double(delta),curr_ptr,error_string)
        if error_string.value != b'':
            raise ValueError("Error in initialization: {0}".format(error_string.value.decode()))

    def load_profiles(self, f_file='none', foffset=None, p_file='none', eta_file='none', f_NI_file='none'):
        r'''! Load flux function profiles (\f$F*F'\f$ and \f$P'\f$) from files

        @param f_file File containing \f$F*F'\f$ (or \f$F'\f$ if `mode=0`) definition
        @param foffset Value of \f$F0=R0*B0\f$
        @param p_file File containing \f$P'\f$ definition
        @param eta_file File containing $\eta$ definition
        @param f_NI_file File containing non-inductive \f$F*F'\f$ definition
        '''
        if foffset is not None:
            self._F0 = foffset
        f_file_c = self._oft_env.path2c(f_file)
        p_file_c = self._oft_env.path2c(p_file)
        eta_file_c = self._oft_env.path2c(eta_file)
        f_NI_file_c = self._oft_env.path2c(f_NI_file)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_load_profiles(self._tMaker_ptr,f_file_c,c_double(self._F0),p_file_c,eta_file_c,f_NI_file_c,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def set_profiles(self, ffp_prof=None, foffset=None, pp_prof=None, ffp_NI_prof=None, keep_files=False):
        r'''! Set flux function profiles (\f$F*F'\f$ and \f$P'\f$) using a piecewise linear definition

        @param ffp_prof Dictionary object containing FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param foffset Value of \f$F0=R0*B0\f$
        @param pp_prof Dictionary object containing P' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param ffp_NI_prof Dictionary object containing non-inductive FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param keep_files Retain temporary profile files
        '''
        delete_files = []
        ffp_file = 'none'
        if ffp_prof is not None:
            ffp_file = self._oft_env.unique_tmpfile('tokamaker_f.prof')
            create_prof_file(self, ffp_file, ffp_prof, "F*F'")
            delete_files.append(ffp_file)
        pp_file = 'none'
        if pp_prof is not None:
            pp_file = self._oft_env.unique_tmpfile('tokamaker_p.prof')
            create_prof_file(self, pp_file, pp_prof, "P'")
            delete_files.append(pp_file)
        eta_file = 'none'
        ffp_NI_file = 'none'
        if ffp_NI_prof is not None:
            ffp_NI_file = self._oft_env.unique_tmpfile('tokamaker_ffp_NI.prof')
            create_prof_file(self, ffp_NI_file, ffp_NI_prof, "ffp_NI")
            delete_files.append(ffp_NI_file)
        if foffset is not None:
            self._F0 = foffset
        self.load_profiles(ffp_file,foffset,pp_file,eta_file,ffp_NI_file)
        if not keep_files:
            for file in delete_files:
                try:
                    os.remove(file)
                except:
                    print('Warning: unable to delete temporary file "{0}"'.format(file))

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
        '''! Solve G-S equation with specified constraints, profiles, etc.
        
        @param vacuum Perform vacuum solve? Plasma-related targets (eg. `Ip`) will be ignored.
        '''
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_solve(self._tMaker_ptr,c_bool(vacuum),error_string)
        if error_string.value != b'':
            raise ValueError("Error in solve: {0}".format(error_string.value.decode()))
    
    def vac_solve(self,psi=None,rhs_source=None):
        '''! Solve for vacuum solution (no plasma), with present coil currents
        and optional other currents

        @note If isoflux, flux, or saddle constraints are desired use @ref solve instead.
        
        @param psi Boundary values for vacuum solve
        @param rhs_source Current source (optional)
        '''
        if psi is None:
            psi = numpy.zeros((self.np,),dtype=numpy.float64)
        else:
            if psi.shape[0] != self.np:
                raise IndexError('Incorrect shape of "psi", should be [np]')
            psi = numpy.ascontiguousarray(psi, dtype=numpy.float64)
        rhs_ptr = None
        if rhs_source is not None:
            if rhs_source.shape[0] != self.np:
                raise IndexError('Incorrect shape of "rhs_source", should be [np]')
            rhs_source = numpy.ascontiguousarray(rhs_source, dtype=numpy.float64)
            rhs_ptr = rhs_source.ctypes.data_as(c_double_ptr)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_vac_solve(self._tMaker_ptr,psi,rhs_ptr,error_string)
        if error_string.value != b'':
            raise ValueError("Error in solve: {0}".format(error_string.value.decode()))
        return psi

    def get_stats(self,lcfs_pad=None,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Get information (Ip, q, kappa, etc.) about current G-S equilbirium

        See eq. 1 for `li_normalization='std'` and eq 2. for `li_normalization='iter'`
        in [Jackson et al.](https://dx.doi.org/10.1088/0029-5515/48/12/125002)

        @param lcfs_pad Padding at LCFS for boundary calculations (default: 1.0 for limited; 0.99 for diverted)
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        @result Dictionary of equilibrium parameters
        '''
        if lcfs_pad is None:
            lcfs_pad = 0.0
            if self.diverted or (not self.settings.free_boundary):
                lcfs_pad = 0.01
        _,qvals,_,dl,rbounds,zbounds = self.get_q(numpy.r_[1.0-lcfs_pad,0.95,0.02],compute_geo=True) # Given backward so last point is LCFS (for dl)
        Ip,centroid,vol,pvol,dflux,tflux,Bp_vol = self.get_globals()
        if beta_Ip is not None:
            Ip = beta_Ip
        p_psi = numpy.linspace(0.0,1.0,100)
        p_psi[0] = 0.001
        _,_,_,p,_ = self.get_profiles(p_psi)
        if self.diverted:
            x_points, _ = self.get_xpoints()
            x_active = x_points[-1,:]
            if x_active[1] < zbounds[0,1]:
                zbounds[0,:] = x_active
                # Find first X-point on opposite side
                for i in range(x_points.shape[0]-1):
                    if x_points[-(i+1),1] > zbounds[1,1]:
                        zbounds[1,:] = x_points[-(i+1),:]
                        break
            elif x_active[1] > zbounds[1,1]:
                zbounds[1,:] = x_active
                # Find first X-point on opposite side
                for i in range(x_points.shape[0]-1):
                    if x_points[-(i+1),1] < zbounds[0,1]:
                        zbounds[0,:] = x_points[-(i+1),:]
                        break
        # Compute normalized inductance
        if li_normalization.lower() == 'std':
            li = (Bp_vol/vol)/numpy.power(mu0*Ip/dl,2)
        elif li_normalization.lower() == 'iter':
            li = 2.0*Bp_vol/(numpy.power(mu0*Ip,2)*self.o_point[0])
        else:
            raise ValueError('Invalid "li_normalization"')
        # Compute geometric major/minor radius
        if geom_type == 'mid':
            rlcfs = self.trace_surf(1.0-lcfs_pad)
            rlcfs = rlcfs[rlcfs[:,0]<self.o_point[0],:]
            iLFS = 0
            iHFS = (numpy.abs(rlcfs[:,1]-self.o_point[1])).argmin()
            R_geo = (rlcfs[iLFS,0]+rlcfs[iHFS,0])/2.0
            a_geo = (rlcfs[iLFS,0]-rlcfs[iHFS,0])/2.0
        elif geom_type == 'max':
            R_geo = (rbounds[1,0]+rbounds[0,0])/2.0
            a_geo = (rbounds[1,0]-rbounds[0,0])/2.0
        else:
            raise ValueError('Invalid "geom_type"')
        # Build dictionary
        eq_stats = {
            'Ip': Ip,
            'Ip_centroid': centroid,
            'kappa': (zbounds[1,1]-zbounds[0,1])/(rbounds[1,0]-rbounds[0,0]),
            'kappaU': (zbounds[1,1]-self.o_point[1])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'kappaL': (self.o_point[1]-zbounds[0,1])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'delta': ((rbounds[1,0]+rbounds[0,0])/2.0-(zbounds[1,0]+zbounds[0,0])/2.0)*2.0/(rbounds[1,0]-rbounds[0,0]),
            'deltaU': ((rbounds[1,0]+rbounds[0,0])/2.0-zbounds[1,0])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'deltaL': ((rbounds[1,0]+rbounds[0,0])/2.0-zbounds[0,0])*2.0/(rbounds[1,0]-rbounds[0,0]),
            'R_geo': R_geo,
            'a_geo': a_geo,
            'vol': vol,
            'q_0': qvals[2],
            'q_95': qvals[1],
            'P_ax': p[0],
            'P_max': p.max(),
            'W_MHD': pvol*1.5,
            'beta_pol': 100.0*(2.0*pvol*mu0/vol)/numpy.power(Ip*mu0/dl,2),
            'dflux': dflux,
            'tflux': tflux,
            'l_i': li
        }
        if self._F0 > 0.0:
            eq_stats['beta_tor'] = 100.0*(2.0*pvol*mu0/vol)/(numpy.power(self._F0/R_geo,2))
            eq_stats['beta_n'] = eq_stats['beta_tor']*eq_stats['a_geo']*(self._F0/R_geo)/(Ip/1.E6)
        return eq_stats

    def print_info(self,lcfs_pad=0.01,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Print information (Ip, q, etc.) about current G-S equilbirium
        
        @param lcfs_pad Padding at LCFS for boundary calculations
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        '''
        eq_stats = self.get_stats(lcfs_pad=lcfs_pad,li_normalization=li_normalization,geom_type=geom_type,beta_Ip=beta_Ip)
        print("Equilibrium Statistics:")
        if self.diverted:
            print("  Topology                =   Diverted")
        else:
            print("  Topology                =   Limited")
        print("  Toroidal Current [A]    =   {0:11.4E}".format(eq_stats['Ip']))
        print("  Current Centroid [m]    =   {0:6.3F} {1:6.3F}".format(*eq_stats['Ip_centroid']))
        if self.settings.dipole_mode:
            print("  Inner limiter [m]       =   {0:6.3F} {1:6.3F}".format(*self.o_point))
        else:
            print("  Magnetic Axis [m]       =   {0:6.3F} {1:6.3F}".format(*self.o_point))
        print("  Elongation              =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['kappa'],eq_stats['kappaU'],eq_stats['kappaL']))
        print("  Triangularity           =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['delta'],eq_stats['deltaU'],eq_stats['deltaL']))
        print("  Plasma Volume [m^3]     =   {0:6.3F}".format(eq_stats['vol']))
        if not self.settings.dipole_mode:
            print("  q_0, q_95               =   {0:6.3F} {1:6.3F}".format(eq_stats['q_0'],eq_stats['q_95']))
        if self.settings.dipole_mode:
            print("  Peak Pressure [Pa]      =   {0:11.4E}".format(eq_stats['P_max']))
        else:
            print("  Plasma Pressure [Pa]    =   Axis: {0:11.4E}, Peak: {1:11.4E}".format(eq_stats['P_ax'], eq_stats['P_max']))
        print("  Stored Energy [J]       =   {0:11.4E}".format(eq_stats['W_MHD']))
        print("  <Beta_pol> [%]          =   {0:7.4F}".format(eq_stats['beta_pol']))
        if 'beta_tor' in eq_stats:
            print("  <Beta_tor> [%]          =   {0:7.4F}".format(eq_stats['beta_tor']))
        if 'beta_n' in eq_stats:
            print("  <Beta_n>   [%]          =   {0:7.4F}".format(eq_stats['beta_n']))
        print("  Diamagnetic flux [Wb]   =   {0:11.4E}".format(eq_stats['dflux']))
        print("  Toroidal flux [Wb]      =   {0:11.4E}".format(eq_stats['tflux']))
        print("  l_i                     =   {0:7.4F}".format(eq_stats['l_i']))
    
    def set_isoflux(self,isoflux,weights=None,grad_wt_lim=-1.0,ref_points=None):
        r'''! Set isoflux constraint points (all points lie on a flux surface)

        To constraint points more uniformly in space additional weighting based on
        the gradient of $\psi$ at each point can also be included by setting
        `grad_wt_lim>0`. When set the actual weight will be
        $w_i * min(grad_wt_lim,|\nabla \psi|_{max} / |\nabla \psi|_i)$

        @param isoflux List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        @param grad_wt_lim Limit on gradient-based weighting (negative to disable)
        @param ref_points Reference points for each isoflux point [:,2] (default: `isoflux[0,:]` is used for all points)
        '''
        if isoflux is None:
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_isoflux(self._tMaker_ptr,numpy.zeros((1,1)),numpy.zeros((1,1)),numpy.zeros((1,)),0,grad_wt_lim,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._isoflux_targets = None
        else:
            if ref_points is None:
                ref_points = numpy.zeros((isoflux.shape[0]-1,2), dtype=numpy.float64)
                ref_points[:,0] = isoflux[0,0]; ref_points[:,1] = isoflux[0,1]
                isoflux = isoflux[1:,:]
                if weights is not None:
                    weights = weights[1:]
            if ref_points.shape[0] != isoflux.shape[0]:
                raise ValueError('Shape of "ref_points" does not match first dimension of "isoflux"')
            if weights is None:
                weights = numpy.ones((isoflux.shape[0],), dtype=numpy.float64)
            if weights.shape[0] != isoflux.shape[0]:
                raise ValueError('Shape of "weights" does not match first dimension of "isoflux"')
            isoflux = numpy.ascontiguousarray(isoflux, dtype=numpy.float64)
            weights = numpy.ascontiguousarray(weights, dtype=numpy.float64)
            ref_points = numpy.ascontiguousarray(ref_points, dtype=numpy.float64)
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_isoflux(self._tMaker_ptr,isoflux,ref_points,weights,isoflux.shape[0],grad_wt_lim,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._isoflux_targets = isoflux.copy()
    
    def set_flux(self,locations,targets,weights=None): #,grad_wt_lim=-1.0):
        r'''! Set explicit flux constraint points \f$ \psi(x_i) \f$

        @param locations List of points defining constraints [:,2]
        @param targets Target \f$ \psi \f$ value at each point [:]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        if locations is None:
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_flux(self._tMaker_ptr,numpy.zeros((1,1)),numpy.zeros((1,)),numpy.zeros((1,)),0,-1.0,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._flux_targets = None
        else:
            if targets.shape[0] != locations.shape[0]:
                raise ValueError('Shape of "targets" does not match first dimension of "locations"')
            if weights is None:
                weights = numpy.ones((locations.shape[0],), dtype=numpy.float64)
            if weights.shape[0] != locations.shape[0]:
                raise ValueError('Shape of "weights" does not match first dimension of "locations"')
            locations = numpy.ascontiguousarray(locations, dtype=numpy.float64)
            targets = numpy.ascontiguousarray(targets, dtype=numpy.float64)
            weights = numpy.ascontiguousarray(weights, dtype=numpy.float64)
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_flux(self._tMaker_ptr,locations,targets,weights,locations.shape[0],-1.0,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._flux_targets = (locations.copy(), targets.copy())
    
    def set_saddles(self,saddles,weights=None):
        '''! Set saddle constraint points (poloidal field should vanish at each point)

        @param saddles List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        if saddles is None:
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_saddles(self._tMaker_ptr,numpy.zeros((1,1)),numpy.zeros((1,)),0,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._saddle_targets = None
        else:
            if weights is None:
                weights = numpy.ones((saddles.shape[0],), dtype=numpy.float64)
            if weights.shape[0] != saddles.shape[0]:
                raise ValueError('Shape of "weights" does not match first dimension of "saddles"')
            saddles = numpy.ascontiguousarray(saddles, dtype=numpy.float64)
            weights = numpy.ascontiguousarray(weights, dtype=numpy.float64)
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_saddles(self._tMaker_ptr,saddles,weights,saddles.shape[0],error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._saddle_targets = saddles.copy()
    
    def set_targets(self,Ip=None,Ip_ratio=None,pax=None,estore=None,R0=None,V0=None,retain_previous=False):
        r'''! Set global target values

        @note Values that are not specified are reset to their defaults on each call unless `retain_previous=True`.

        @param alam Scale factor for \f$F*F'\f$ term (disabled if `Ip` is set)
        @param pnorm Scale factor for \f$P'\f$ term (disabled if `pax`, `estore`, or `R0` are set)
        @param Ip Target plasma current [A] (disabled if `OFT_env.float_disable_flag`)
        @param Ip_ratio Amplitude of net plasma current contribution from FF' compared to P' (disabled if `OFT_env.float_disable_flag`)
        @param pax Target axis pressure [Pa] (disabled if `OFT_env.float_disable_flag` or if `estore` is set)
        @param estore Target sotred energy [J] (disabled if `OFT_env.float_disable_flag`)
        @param R0 Target major radius for magnetic axis (disabled if `OFT_env.float_disable_flag` or if `pax` or `estore` are set)
        @param V0 Target vertical position for magnetic axis (disabled if `OFT_env.float_disable_flag`)
        @param retain_previous Keep previously set targets unless explicitly updated? (default: False)
        '''
        # Reset all targets unless specified
        if not retain_previous:
            self._Ip_target.value = self._oft_env.float_disable_flag
            self._estore_target.value = self._oft_env.float_disable_flag
            self._pax_target.value = self._oft_env.float_disable_flag
            self._Ip_ratio_target.value = self._oft_env.float_disable_flag
            self._R0_target.value = self._oft_env.float_disable_flag
            self._V0_target.value = self._oft_env.float_disable_flag
        # Set new targets
        if Ip is not None:
            if (Ip <= 0.0) and (not self._oft_env.float_is_disabled(Ip)):
                raise ValueError("`Ip_target` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._Ip_target.value=Ip
        if estore is not None:
            if (estore <= 0.0) and (not self._oft_env.float_is_disabled(estore)):
                raise ValueError("`estore` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._estore_target.value=estore
        if pax is not None:
            if (pax <= 0.0) and (not self._oft_env.float_is_disabled(pax)):
                raise ValueError("`pax` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._pax_target.value=pax
        if Ip_ratio is not None:
            self._Ip_ratio_target.value=Ip_ratio
        if R0 is not None:
            if (R0 <= 0.0) and (not self._oft_env.float_is_disabled(R0)):
                raise ValueError("`R0` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._R0_target.value=R0
        if V0 is not None:
            self._V0_target.value=V0
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_targets(self._tMaker_ptr,self._Ip_target,self._Ip_ratio_target,self._pax_target,self._estore_target,self._R0_target,self._V0_target,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def get_targets(self):
        r'''! Get global target values

        @result Dictionary of global target values
        '''
        # Get targets
        target_dict = {}
        if (not self._oft_env.float_is_disabled(self._Ip_target.value)):
            target_dict['Ip'] = self._Ip_target.value
        if (not self._oft_env.float_is_disabled(self._estore_target.value)):
            target_dict['estore'] = self._estore_target.value
        if (not self._oft_env.float_is_disabled(self._pax_target.value)):
            target_dict['pax'] = self._pax_target.value
        if (not self._oft_env.float_is_disabled(self._Ip_ratio_target.value)):
            target_dict['Ip_ratio'] = self._Ip_ratio_target.value
        if (not self._oft_env.float_is_disabled(self._R0_target.value)):
            target_dict['R0'] = self._R0_target.value
        if (not self._oft_env.float_is_disabled(self._V0_target.value)):
            target_dict['V0'] = self._V0_target.value
        return target_dict

    def get_delstar_curr(self,psi):
        r'''! Get toroidal current density from \f$ \psi \f$ through \f$ \Delta^{*} \f$ operator
 
        @param psi \f$ \psi \f$ corresponding to desired current density
        @result \f$ J_{\phi} = \textrm{M}^{-1} \Delta^{*} \psi \f$
        '''
        curr = numpy.copy(psi)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_dels_curr(self._tMaker_ptr,curr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return curr/mu0
    
    def get_jtor_plasma(self):
        r'''! Get plasma toroidal current density for current equilibrium
 
        @result \f$ J_{\phi} \f$ by evalutating RHS source terms
        '''
        curr = numpy.zeros((self.np,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_jtor(self._tMaker_ptr,curr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return curr/mu0

    def get_psi(self,normalized=True):
        r'''! Get poloidal flux values on node points

        @param normalized Normalize (and offset) poloidal flux
        @result \f$\hat{\psi} = \frac{\psi-\psi_0}{\psi_a-\psi_0} \f$ or \f$\psi\f$
        '''
        psi = numpy.zeros((self.np,),dtype=numpy.float64)
        psi_lim = c_double()
        psi_max = c_double()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_psi(self._tMaker_ptr,psi,ctypes.byref(psi_lim),ctypes.byref(psi_max),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        if normalized:
            psi = (psi-psi_lim.value)/(psi_max.value-psi_lim.value)
            if self.psi_convention == 0:
                psi = 1.0 - psi
        return psi

    def set_psi(self,psi,update_bounds=False):
        '''! Set poloidal flux values on node points

        @param psi Poloidal flux values (should not be normalized!)
        @param update_bounds Update plasma bounds by determining new limiting points
        '''
        if psi.shape[0] != self.np:
            raise IndexError('Incorrect shape of "psi", should be [np]')
        psi = numpy.ascontiguousarray(psi, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_psi(self._tMaker_ptr,psi,c_bool(update_bounds),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def set_psi_dt(self,psi0,dt):
        '''! Set reference poloidal flux and time step for eddy currents in .solve()

        @param psi0 Reference poloidal flux at t-dt (unnormalized)
        @param dt Time since reference poloidal flux
        '''
        if psi0.shape[0] != self.np:
            raise IndexError('Incorrect shape of "psi0", should be [np]')
        psi0 = numpy.ascontiguousarray(psi0, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_psi_dt(self._tMaker_ptr,psi0,c_double(dt),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def get_field_eval(self,field_type):
        r'''! Create field interpolator for vector potential

        @param field_type Field to interpolate, must be one of ("B", "psi", "F", "P", "dPSI", "dBr", "dBt", or "dBz")
        @result Field interpolation object
        '''
        #
        mode_map = {'B': 1, 'PSI': 2, 'F': 3, 'P': 4, 'DPSI': 5, 'DBR': 6, 'DBT': 7, 'DBZ': 8}
        imode = mode_map.get(field_type.upper())
        if imode is None:
            raise ValueError('Invalid field type ("B", "psi", "F", "P", "dPSI", "dBr", "dBt", "dBz")')
        #
        int_obj = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_field_eval(self._tMaker_ptr,imode,ctypes.byref(int_obj),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        field_dim = 1
        if imode == 1:
            field_dim = 3
        elif imode >= 5:
            field_dim = 2
        return TokaMaker_field_interpolator(self._tMaker_ptr,int_obj,imode,field_dim)
    
    def get_coil_currents(self):
        '''! Get currents in each coil [A] and coil region [A-turns]

        @result Coil currents [ncoils], Coil currents by region [nregs]
        '''
        currents = numpy.zeros((self.ncoils,),dtype=numpy.float64)
        currents_reg = numpy.zeros((self.nregs,),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_coil_currents(self._tMaker_ptr,currents,currents_reg,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        current_dict = {}
        for coil_key, coil_set in self.coil_sets.items():
            current_dict[coil_key] = currents[coil_set['id']]
        return current_dict, currents_reg

    def get_coil_Lmat(self):
        r'''! Get mutual inductance matrix between coils

        @note This is the inductance in terms of A-turns. To get in terms of
        current in a single of the \f$n\f$ windings you must multiply by \f$n_i*n_j\f$.

        @result L[ncoils+1,ncoils+1]
        '''
        Lmat = numpy.zeros((self.ncoils+1,self.ncoils+1),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_coil_Lmat(self._tMaker_ptr,Lmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
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
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_trace_surf(self._tMaker_ptr,c_double(psi),ctypes.byref(points_loc),ctypes.byref(npoints),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        if npoints.value > 0:
            return numpy.ctypeslib.as_array(points_loc,shape=(npoints.value, 2))
        else:
            return None
    
    def get_q(self,psi=None,psi_pad=0.02,npsi=50,compute_geo=False):
        r'''! Get q-profile at specified or uniformly spaced points

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @param compute_geo Compute geometric values for LCFS
        @result \f$\hat{\psi}\f$, \f$q(\hat{\psi})\f$, \f$[<R>,<1/R>,dV/dPsi]\f$, length of last surface,
        [r(R_min),r(R_max)], [r(z_min),r(z_max)]
        '''
        if psi is None:
            psi = numpy.linspace(psi_pad,1.0-psi_pad,npsi,dtype=numpy.float64)
            if self.psi_convention == 0:
                psi = numpy.ascontiguousarray(numpy.flip(psi), dtype=numpy.float64)
                psi_save = 1.0-psi
        else:
            if self.psi_convention == 0:
                psi_save = numpy.copy(psi)
                psi = numpy.ascontiguousarray(1.0-psi, dtype=numpy.float64)
        qvals = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        ravgs = numpy.zeros((3,psi.shape[0]), dtype=numpy.float64)
        if compute_geo:
            dl = c_double(1.0)
        else:
            dl = c_double(-1.0)
        rbounds = numpy.zeros((2,2),dtype=numpy.float64)
        zbounds = numpy.zeros((2,2),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_q(self._tMaker_ptr,psi.shape[0],psi,qvals,ravgs,ctypes.byref(dl),rbounds,zbounds,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        if self.psi_convention == 0:
            if compute_geo:
                return psi_save,qvals,ravgs,dl.value,rbounds,zbounds
            else:
                return psi_save,qvals,ravgs,None,None,None
        else:
            if compute_geo:
                return psi,qvals,ravgs,dl.value,rbounds,zbounds
            else:
                return psi,qvals,ravgs,None,None,None

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
                psi = numpy.ascontiguousarray(numpy.flip(psi), dtype=numpy.float64)
                psi_save = 1.0 - psi
        else:
            if self.psi_convention == 0:
                psi_save = numpy.copy(psi)
                psi = numpy.ascontiguousarray(1.0-psi, dtype=numpy.float64)
        fc = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        r_avgs = numpy.zeros((3,psi.shape[0]), dtype=numpy.float64)
        modb_avgs = numpy.zeros((2,psi.shape[0]), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_sauter_fc(self._tMaker_ptr,psi.shape[0],psi,fc,r_avgs,modb_avgs,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
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
        Bp_vol = c_double()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_globals(self._tMaker_ptr,ctypes.byref(Ip),centroid,ctypes.byref(vol),ctypes.byref(pvol),
            ctypes.byref(dflux),ctypes.byref(tflux),ctypes.byref(Bp_vol),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Ip.value, centroid, vol.value, pvol.value, dflux.value, tflux.value, Bp_vol.value

    def calc_loopvoltage(self):
        r'''! Get plasma loop voltage

        @param eta Dictionary object containing resistivity profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @param ffp_NI Dictionary object containing non-inductive FF' profile ['y'] and sampled locations 
        in normalized Psi ['x']
        @result Vloop [Volts]
        '''
        V_loop = c_double()
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_gs_calc_vloop(self._tMaker_ptr,ctypes.byref(V_loop),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        #
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
                psi = numpy.ascontiguousarray(numpy.flip(psi), dtype=numpy.float64)
                psi_save = 1.0 - psi
        else:
            if self.psi_convention == 0:
                psi_save = numpy.copy(psi)
                psi = numpy.ascontiguousarray(1.0-psi, dtype=numpy.float64)
        f = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        fp = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        p = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        pp = numpy.zeros((psi.shape[0],), dtype=numpy.float64)
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_profs(self._tMaker_ptr,psi.shape[0],psi,f,fp,p,pp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        #
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
    
    def set_coil_currents(self, currents=None):
        '''! Set coil currents

        @param currents Current in each coil [A]
        '''
        current_array = numpy.zeros((self.ncoils,), dtype=numpy.float64)
        if currents is not None:
            for coil_key, coil_current in currents.items():
                current_array[self.coil_sets[coil_key]['id']] = coil_current
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_currents(self._tMaker_ptr,current_array,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def update_settings(self):
        '''! Update settings after changes to values in python'''
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_settings(self._tMaker_ptr,ctypes.byref(self.settings),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def set_dipole_a(self,a_exp):
        r'''! Update anisotropy exponent `a` when dipole mode is used
        
        @param a_exp New value for `a` exponent
        '''
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_dipole_a(self._tMaker_ptr,a_exp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def area_integral(self,field,reg_mask=-1):
        r'''! Compute area integral of field over a specified region

        @param field Field to integrate [np,]
        @param reg_mask ID of region for integration (negative for whole mesh)
        @result \f$ \int f dA \f$
        '''
        result = c_double(0.0)
        field = numpy.ascontiguousarray(field, dtype=numpy.float64)
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_area_int(self._tMaker_ptr,field,c_int(reg_mask),ctypes.byref(result),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return result.value
    
    def flux_integral(self,psi_vals,field_vals):
        r'''! Compute area integral of field over a specified region

        @param field Field to integrate [np,]
        @param reg_mask ID of region for integration (negative for whole mesh)
        @result \f$ \int f dA \f$
        '''
        if psi_vals.shape[0] != field_vals.shape[0]:
            raise ValueError('"psi_vals" and "field_vals" must be the same length')
        if self.psi_convention == 0:
            psi_vals = numpy.flip(1.0-psi_vals)
            field_vals = numpy.flip(field_vals)
        result = c_double(0.0)
        psi_vals = numpy.ascontiguousarray(psi_vals, dtype=numpy.float64)
        field_vals = numpy.ascontiguousarray(field_vals, dtype=numpy.float64)
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_flux_int(self._tMaker_ptr,psi_vals,field_vals,c_int(psi_vals.shape[0]),ctypes.byref(result),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return result.value

    def plot_machine(self,fig,ax,vacuum_color='whitesmoke',cond_color='gray',limiter_color='k',
                     coil_color='gray',coil_colormap=None,coil_symmap=False,coil_scale=1.0,coil_clabel=r'$I_C$ [A]',colorbar=None):
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
        @param colorbar Colorbar instance to overwrite (None to add)
        @result Colorbar instance for coil colors or None
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
            if self.ncoils >0:
                mesh_currents = region_currents[self.reg-1]
            # Adjust current in coils with non-uniform distribution
            if len(self.dist_coils)>0:
                for _, coil_obj in self.coil_sets.items():
                    if (coil_id:=coil_obj["id"]) in self.dist_coils.keys():
                        for sub_coil in coil_obj["sub_coils"]:
                            mask = (self.reg==sub_coil["reg_id"])
                            face_currents = numpy.mean(self.dist_coils[coil_id][self.lc],axis=1)
                            mesh_currents[mask] *= face_currents[mask]
            mask = (abs(mesh_currents) > 0.0)
            if mask.sum() > 0.0:
                mesh_currents *= coil_scale
                if coil_symmap:
                    max_curr = abs(mesh_currents).max()
                    clf = ax.tripcolor(self.r[:,0], self.r[:,1], self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap, vmin=-max_curr, vmax=max_curr)
                else:
                    clf = ax.tripcolor(self.r[:,0], self.r[:,1], self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap)
                if coil_clabel is not None:
                    cax = None
                    if colorbar is not None:
                        cax = colorbar.ax
                    colorbar = fig.colorbar(clf,ax=ax,cax=cax,label=coil_clabel)
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
            for lim_contour in self.lim_contours:
                ax.plot(lim_contour[:,0],lim_contour[:,1],color=limiter_color)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
        return colorbar

    def plot_constraints(self,fig,ax,isoflux_color='tab:red',isoflux_marker='+',saddle_color='tab:green',saddle_marker='x'):
        '''! Plot geometry constraints

        @param fig Figure to add to
        @param ax Axis to add to
        @param isoflux_color Color of isoflux points (None to disable)
        @param saddle_color Color of saddle points (None to disable)
        '''
        # Plot isoflux constraints
        if (isoflux_color is not None) and (self._isoflux_targets is not None):
            ax.plot(self._isoflux_targets[:,0],self._isoflux_targets[:,1],color=isoflux_color,marker=isoflux_marker,linestyle='none')
        # Plot saddle constraints
        if (saddle_color is not None) and (self._saddle_targets is not None):
            ax.plot(self._saddle_targets[:,0],self._saddle_targets[:,1],color=saddle_color,marker=saddle_marker,linestyle='none')

    def plot_psi(self,fig,ax,psi=None,normalized=True,
                 plasma_color=None,plasma_nlevels=8,plasma_levels=None,plasma_colormap=None,plasma_linestyles=None,
                 vacuum_color='darkgray',vacuum_nlevels=8,vacuum_levels=None,vacuum_colormap=None,vacuum_linestyles=None,
                 xpoint_color='k',xpoint_marker='x',xpoint_inactive_alpha=0.5,opoint_color='k',opoint_marker='*'):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param psi Flux values to plot (otherwise `self.get_psi()` is called)
        @param normalized Retreive normalized flux, or assume normalized psi if passed as argument
        @param plasma_color Color for plasma contours
        @param plasma_nlevels Number of plasma contours
        @param plasma_levels Explicit levels for plasma contours
        @param plasma_colormap Colormap for plasma contours (cannot be specified with `plasma_color`)
        @param plasma_linestyles Linestyle for plasma contours
        @param vacuum_color Color for vacuum contours
        @param vacuum_nlevels Number of vacuum contours
        @param vacuum_levels Explicit levels for vacuum contours (cannot be specified with `vacuum_color`)
        @param vacuum_colormap Colormap for vacuum contours
        @param vacuum_linestyles Linestyle for vacuum contours
        @param xpoint_color Color for X-point markers (None to disable)
        @param xpoint_marker Marker style for X-points
        @param xpoint_inactive_alpha Alpha value for inactive X-points
        @param opoint_color Color for O-point markers (None to disable)
        @param opoint_marker Marker style for O-points
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
            ax.tricontour(self.r[:,0],self.r[:,1],self.lc,psi,levels=vacuum_levels,colors=vacuum_color,cmap=vacuum_colormap,linestyles=vacuum_linestyles)
        if plasma_levels is not None:
            ax.tricontour(self.r[:,0],self.r[:,1],self.lc,psi,levels=plasma_levels,colors=plasma_color,cmap=plasma_colormap,linestyles=plasma_linestyles)

        # Plot saddle points
        if xpoint_color is not None:
            x_points, _ = self.get_xpoints()
            if x_points is not None:
                ax.plot(x_points[-1,0], x_points[-1,1], color=xpoint_color, marker=xpoint_marker, linestyle='none')
                ax.plot(x_points[:-1,0], x_points[:-1,1], color=xpoint_color, marker=xpoint_marker, linestyle='none', alpha=xpoint_inactive_alpha)
        if (opoint_color is not None) and (self.o_point[0] > 0.0):
            ax.plot(self.o_point[0], self.o_point[1], color=opoint_color, marker=opoint_marker)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
    
    def get_conductor_currents(self,psi,cell_centered=False):
        r'''! Get toroidal current density in conducting regions for a given \f$ \psi \f$

        @param psi Psi corresponding to field with conductor currents (eg. from time-dependent simulation)
        @param cell_centered Get currents at cell centers
        '''
        curr = self.get_delstar_curr(psi)
        if cell_centered:
            mesh_currents = numpy.zeros((self.lc.shape[0],))
        # Loop over conducting regions and get mask/fields
        mask = numpy.zeros((self.lc.shape[0],), dtype=numpy.int32)
        for _, cond_reg in self._cond_dict.items():
            eta = cond_reg.get('eta',-1.0)
            if eta > 0:
                mask_tmp = (self.reg == cond_reg['reg_id'])
                if cell_centered:
                    mesh_currents[mask_tmp] = numpy.sum(curr[self.lc[mask_tmp,:]],axis=1)/3.0
                mask = numpy.logical_or(mask,mask_tmp)
        if cell_centered:
            return mask, mesh_currents
        else:
            return mask, curr
    
    def get_conductor_source(self,dpsi_dt):
        r'''! Get toroidal current density in conducting regions for a \f$ d \psi / dt \f$ source

        @param dpsi_dt dPsi/dt source eddy currents (eg. from linear stability)
        '''
        # Apply 1/R scale (avoiding divide by zero)
        curr = dpsi_dt.copy()
        curr[self.r[:,0]>0.0] /= self.r[self.r[:,0]>0.0,0]
        # Compute cell areas
        have_noncontinuous = False
        for _, cond_reg in self._cond_dict.items():
            if 'noncontinuous' in cond_reg:
                have_noncontinuous = True
                break
        if have_noncontinuous:
            area = numpy.zeros((self.lc.shape[0],))
            for i in range(self.nc):
                v1 = self.r[self.lc[i,1],:]-self.r[self.lc[i,0],:]
                v2 = self.r[self.lc[i,2],:]-self.r[self.lc[i,0],:]
                area[i] = numpy.linalg.norm(numpy.cross(v1,v2))/2.0
        #
        mesh_currents = numpy.zeros((self.lc.shape[0],))
        # Loop over conducting regions and get mask/fields
        mask = numpy.zeros((self.lc.shape[0],), dtype=numpy.int32)
        for _, cond_reg in self._cond_dict.items():
            eta = cond_reg.get('eta',-1.0)
            if eta > 0:
                mask_tmp = (self.reg == cond_reg['reg_id'])
                field_tmp = -dpsi_dt/eta
                mesh_currents[mask_tmp] = numpy.sum(field_tmp[self.lc[mask_tmp,:]],axis=1)/3.0
                if cond_reg.get('noncontinuous',False):
                    mesh_currents[mask_tmp] -= (mesh_currents[mask_tmp]*area[mask_tmp]).sum()/area[mask_tmp].sum()
                mask = numpy.logical_or(mask,mask_tmp)
        return mask, mesh_currents
    
    def plot_eddy(self,fig,ax,psi=None,dpsi_dt=None,nlevels=40,colormap='jet',clabel=r'$J_w$ [$A/m^2$]',symmap=False):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param psi Psi corresponding to eddy currents (eg. from time-dependent simulation)
        @param dpsi_dt dPsi/dt source eddy currents (eg. from linear stability)
        @param nlevels Number contour lines used for shading (with "psi" only)
        @param colormap Colormap to use for shadings
        @param clabel Label for colorbar (None to disable colorbar)
        @result Colorbar object
        '''
        if psi is not None:
            mask, plot_field = self.get_conductor_currents(psi,cell_centered=(nlevels < 0))
        elif dpsi_dt is not None:
            mask, plot_field = self.get_conductor_source(dpsi_dt)
        if plot_field.shape[0] == self.nc:
            if symmap:
                max_curr = abs(plot_field).max()
                clf = ax.tripcolor(self.r[:,0],self.r[:,1],self.lc[mask,:],plot_field[mask],cmap=colormap,vmin=-max_curr,vmax=max_curr)
            else:
                clf = ax.tripcolor(self.r[:,0],self.r[:,1],self.lc[mask],plot_field[mask],cmap=colormap)
        else:
            if symmap:
                max_curr = abs(plot_field[self.lc[mask,:]]).max(axis=None)
                clf = ax.tricontourf(self.r[:,0],self.r[:,1],self.lc[mask],plot_field,nlevels,cmap=colormap,vmin=-max_curr,vmax=max_curr)
            else:
                clf = ax.tricontourf(self.r[:,0],self.r[:,1],self.lc[mask],plot_field,nlevels,cmap=colormap)
        if clabel is not None:
            cb = fig.colorbar(clf,ax=ax)
            cb.set_label(clabel)
        else:
            cb = None
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
        return cb

    def get_vfixed(self):
        '''! Get required vacuum flux values to balance fixed boundary equilibrium

        @result sampling points [:,2], flux values [:]
        '''
        npts = c_int()
        pts_loc = c_double_ptr()
        flux_loc = c_double_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_vfixed(self._tMaker_ptr,ctypes.byref(npts),ctypes.byref(pts_loc),ctypes.byref(flux_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return numpy.ctypeslib.as_array(pts_loc,shape=(npts.value, 2)), \
            numpy.ctypeslib.as_array(flux_loc,shape=(npts.value,))

    def save_eqdsk(self,filename,nr=65,nz=65,rbounds=None,zbounds=None,run_info='',lcfs_pad=0.01,rcentr=None,truncate_eq=True,limiter_file='',lcfs_pressure=0.0):
        r'''! Save current equilibrium to gEQDSK format

        @param filename Filename to save equilibrium to
        @param nr Number of radial sampling points
        @param nz Number of vertical sampling points
        @param rbounds Extents of grid in R
        @param zbounds Extents of grid in Z
        @param run_info Run information for gEQDSK file (maximum of 40 characters)
        @param lcfs_pad Padding in normalized flux at LCFS
        @param rcentr `RCENTR` value for gEQDSK file (if `None`, geometric axis is used)
        @param truncate_eq Truncate equilibrium at `lcfs_pad`, if `False` \f$ q(\hat{\psi} > 1-pad) = q(1-pad) \f$
        @param limiter_file File containing limiter contour to use instead of TokaMaker limiter
        @param lcfs_pressure Plasma pressure on the LCFS (zero by default)
        '''
        cfilename = self._oft_env.path2c(filename)
        lim_filename = self._oft_env.path2c(limiter_file)
        if len(run_info) > 40:
            raise ValueError('"run_info" cannot be longer than 40 characters')
        crun_info = self._oft_env.path2c(run_info)
        if rbounds is None:
            rbounds = numpy.r_[self.lim_contour[:,0].min(), self.lim_contour[:,0].max()]
            dr = rbounds[1]-rbounds[0]
            rbounds += numpy.r_[-1.0,1.0]*dr*0.05
        if zbounds is None:
            zbounds = numpy.r_[self.lim_contour[:,1].min(), self.lim_contour[:,1].max()]
            dr = zbounds[1]-zbounds[0]
            zbounds += numpy.r_[-1.0,1.0]*dr*0.05
        if rcentr is None:
            rcentr = -1.0
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_save_eqdsk(self._tMaker_ptr,cfilename,c_int(nr),c_int(nz),rbounds,zbounds,crun_info,c_double(lcfs_pad),c_double(rcentr),c_bool(truncate_eq),lim_filename,lcfs_pressure,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def save_ifile(self,filename,npsi=65,ntheta=65,lcfs_pad=0.01,lcfs_pressure=0.0,pack_lcfs=True,single_precision=False):
        r'''! Save current equilibrium to iFile format

        @param filename Filename to save equilibrium to
        @param npsi Number of radial sampling points
        @param ntheta Number of vertical sampling points
        @param lcfs_pad Padding in normalized flux at LCFS
        @param lcfs_pressure Plasma pressure on the LCFS (zero by default)
        @param pack_lcfs Pack toward LCFS with quadraturic sampling?
        @param single_precision Save single precision file? (default: double precision)
        '''
        cfilename = self._oft_env.path2c(filename)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_save_ifile(self._tMaker_ptr,cfilename,npsi,ntheta,lcfs_pad,lcfs_pressure,pack_lcfs,single_precision,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def save_mug(self,filename):
        r'''! Save current equilibrium to MUG transfer format

        @param filename Filename to save equilibrium to
        '''
        cfilename = self._oft_env.path2c(filename)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_save_mug(self._tMaker_ptr,cfilename,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def set_coil_current_dist(self,coil_name,curr_dist):
        '''! Overwrite coil with non-uniform current distribution.

        @param coil_name Name of coil to modify
        @param curr_dist Relative current density [self.np]
        '''
        if curr_dist.shape[0] != self.np:
            raise IndexError('Incorrect shape of "curr_dist", should be [np]')
        if coil_name not in self.coil_sets:
            raise KeyError('Unknown coil "{0}"'.format(coil_name))
        iCoil = self.coil_sets[coil_name]['id']
        self.dist_coils[iCoil] = curr_dist
        curr_dist = numpy.ascontiguousarray(curr_dist, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_current_dist(self._tMaker_ptr,c_int(iCoil+1),curr_dist,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def eig_wall(self,neigs=4,pm=False):
        r'''! Compute eigenvalues (\f$ 1 / \tau_{L/R} \f$) for conducting structures

        @param neigs Number of eigenvalues to compute
        @param pm Print solver statistics and raw eigenvalues?
        @result eigenvalues[neigs], eigenvectors[neigs,:]
        '''
        eig_vals = numpy.zeros((neigs,2),dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.np),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_eig_wall(self._tMaker_ptr,c_int(neigs),eig_vals,eig_vecs,pm,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return eig_vals, eig_vecs

    def eig_td(self,omega=-1.E4,neigs=4,include_bounds=True,pm=False,damping_scale=-1.0):
        '''! Compute eigenvalues for the linearized time-dependent system

        @param omega Growth rate localization point (eigenvalues closest to this value will be found)
        @param neigs Number of eigenvalues to compute
        @param include_bounds Include bounding flux terms for constant normalized profiles?
        @param pm Print solver statistics and raw eigenvalues?
        @param damping_scale Scale factor for damping term to artificially limit growth rate (negative to disable)?
        @result eigenvalues[neigs], eigenvectors[neigs,:]
        '''
        eig_vals = numpy.zeros((neigs,2),dtype=numpy.float64)
        eig_vecs = numpy.zeros((neigs,self.np),dtype=numpy.float64)
        damp_coeff = abs(omega)*damping_scale
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_eig_td(self._tMaker_ptr,c_double(omega),c_int(neigs),eig_vals,eig_vecs,c_bool(include_bounds),c_double(damp_coeff),pm,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return eig_vals, eig_vecs

    def setup_td(self,dt,lin_tol,nl_tol,pre_plasma=False):
        '''! Setup the time-dependent G-S solver

        @param dt Starting time step
        @param lin_tol Tolerance for linear solver
        @param nl_tol Tolerance for non-linear solver
        @param pre_plasma Use plasma contributions in preconditioner (default: False)
        '''
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_setup_td(self._tMaker_ptr,c_double(dt),c_double(lin_tol),c_double(nl_tol),c_bool(pre_plasma),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
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
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_step_td(self._tMaker_ptr,ctypes.byref(time),ctypes.byref(dt),ctypes.byref(nl_its),ctypes.byref(lin_its),ctypes.byref(nretry),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return time.value, dt.value, nl_its.value, lin_its.value, nretry.value


def solve_with_bootstrap(self,ne,Te,ni,Ti,inductive_jtor,Zeff,jBS_scale=1.0,Zis=[1.],max_iterations=6,initialize_eq=True):
    '''! Self-consistently compute bootstrap contribution from H-mode profiles,
    and iterate solution until all functions of Psi converge. 

    @note if using nis and Zis, dnis_dpsi must be specified in sauter_bootstrap() 
    as a list of impurity gradients over Psi. See 
    https://omfit.io/_modules/omfit_classes/utils_fusion.html for more 
    detailed documentation 

    @note if initialize_eq=True, cubic polynomials will be fit to the core of all 
    kinetic profiles in order to flatten the pedestal. This will initialize the G-S 
    solution at an estimated L-mode pressure profile and using the L-mode bootstrap 
    contribution. Initializing the solver in L-mode before raising the pedestal 
    height increases the likelihood that the solver will converge in H-mode.

    @param ne Electron density profile, sampled over psi_norm
    @param Te Electron temperature profile [eV], sampled over psi_norm
    @param ni Ion density profile, sampled over psi_norm
    @param Ti Ion temperature profile [eV], sampled over psi_norm
    @param inductive_jtor Inductive toroidal current, sampled over psi_norm
    @param Zeff Effective Z profile, sampled over psi_norm
    @param scale_jBS Scalar which can scale bootstrap current profile
    @param nis List of impurity density profiles; NOT USED
    @param Zis List of impurity profile atomic numbers; NOT USED. 
    @param max_iterations Maximum number of H-mode mygs.solve() iterations
    @param initialize_eq Initialize equilibrium solve with flattened pedestal. 
    @param return_jBS Return bootstrap profile alongside err_flag
    '''
    try:
        from omfit_classes.utils_fusion import sauter_bootstrap
    except:
        raise ImportError('omfit_classes.utils_fusion not installed')
    
    def ffprime_from_jtor_pprime(jtor, pprime, R_avg, one_over_R_avg):
        r'''! Convert from J_toroidal to FF' using Grad-Shafranov equation

        @param jtor Toroidal current profile
        @param R_avg Flux averaged R, calculated by TokaMaker
        @param one_over_R_avg Flux averaged 1/R, calculated by TokaMaker
        @param pprime dP/dPsi profile
        '''
        ffprime = 2.0*(jtor -  R_avg * (-pprime)) * (mu0 / one_over_R_avg)
        return ffprime

    kBoltz = eC
    pressure = (kBoltz * ne * Te) + (kBoltz * ni * Ti) # 1.602e-19 * [m^-3] * [eV] = [Pa]

    ### Set new pax target
    self.set_targets(pax=pressure[0])

    ### Reconstruct psi_norm and n_psi from input inductive_jtor
    psi_norm = numpy.linspace(0.,1.,len(inductive_jtor))
    n_psi = len(inductive_jtor)

    def profile_iteration(self,pressure,ne,ni,Te,Ti,psi_norm,n_psi,Zeff,inductive_jtor,jBS_scale,Zis,include_jBS=True):

        pprime = numpy.gradient(pressure) / (numpy.gradient(psi_norm) * (self.psi_bounds[1]-self.psi_bounds[0]))

        ### Get final remaining quantities for Sauter from TokaMaker
        psi,f,_,_,_ = self.get_profiles(npsi=n_psi)
        _,fc,r_avgs,_ = self.sauter_fc(npsi=n_psi)
        ft = 1 - fc # Trapped particle fraction on each flux surface
        eps = r_avgs[2] / r_avgs[0] # Inverse aspect ratio
        _,qvals,ravgs,_,_,_ = self.get_q(npsi=n_psi)
        R_avg = ravgs[0]
        one_over_R_avg = ravgs[1]
        
        if include_jBS:
            ### Calculate flux derivatives for Sauter
            dn_e_dpsi = numpy.gradient(ne) / (numpy.gradient(psi_norm) * (self.psi_bounds[1]-self.psi_bounds[0]))
            dT_e_dpsi = numpy.gradient(Te) / (numpy.gradient(psi_norm) * (self.psi_bounds[1]-self.psi_bounds[0]))
            dn_i_dpsi = numpy.gradient(ni) / (numpy.gradient(psi_norm) * (self.psi_bounds[1]-self.psi_bounds[0]))
            dT_i_dpsi = numpy.gradient(Ti) / (numpy.gradient(psi_norm) * (self.psi_bounds[1]-self.psi_bounds[0]))

            ### Solve for bootstrap current profile. See https://omfit.io/_modules/omfit_classes/utils_fusion.html for more detailed documentation 
            j_BS_neo = sauter_bootstrap(
                                    psi_N=psi_norm,
                                    Te=Te,
                                    Ti=Ti,
                                    ne=ne,
                                    p=pressure,
                                    nis=[ni,],
                                    Zis=Zis,
                                    Zeff=Zeff,
                                    gEQDSKs=[None],
                                    R0=0., # not used
                                    device=None,
                                    psi_N_efit=None,
                                    psiraw=psi*(self.psi_bounds[1]-self.psi_bounds[0]) + self.psi_bounds[0],
                                    R=R_avg,
                                    eps=eps, 
                                    q=qvals,
                                    fT=ft,
                                    I_psi=f,
                                    nt=1,
                                    version='neo_2021',
                                    debug_plots=False,
                                    return_units=True,
                                    return_package=False,
                                    charge_number_to_use_in_ion_collisionality='Koh',
                                    charge_number_to_use_in_ion_lnLambda='Zavg',
                                    dT_e_dpsi=dT_e_dpsi,
                                    dT_i_dpsi=dT_i_dpsi,
                                    dn_e_dpsi=dn_e_dpsi,
                                    dnis_dpsi=[dn_i_dpsi,],
                                    )[0]
                
            inductive_jtor[-1] = 0. ### FORCING inductive_jtor TO BE ZERO AT THE EDGE
            j_BS = j_BS_neo*(R_avg / f) ### Convert into [A/m^2]
            j_BS *= jBS_scale ### Scale j_BS by user specified scalar
            j_BS[-1] = 0. ### FORCING j_BS TO BE ZERO AT THE EDGE
            jtor_total = inductive_jtor + j_BS
        else:
            j_BS = None
            inductive_jtor[-1] = 0. ### FORCING inductive_jtor TO BE ZERO AT THE EDGE
            jtor_total = inductive_jtor
        
        ffprime = ffprime_from_jtor_pprime(jtor_total, pprime, R_avg, one_over_R_avg)

        ffp_prof = {
            'type': 'linterp',
            'x': psi_norm,
            'y': ffprime / ffprime[0]
        }

        pp_prof = {
            'type': 'linterp',
            'x': psi_norm,
            'y': pprime / pprime[0]
        }

        return pp_prof, ffp_prof, j_BS

    if initialize_eq:
        x_trimmed = psi_norm.tolist().copy()
        ne_trimmed = ne.tolist().copy()
        Te_trimmed = Te.tolist().copy()
        ni_trimmed = ni.tolist().copy()
        Ti_trimmed = Ti.tolist().copy()

        ### Remove profile values from psi_norm ~0.5 to ~0.99, leaving single value at the edge
        mid_index = int(len(x_trimmed)/2)
        end_index = len(x_trimmed)-1
        del x_trimmed[mid_index:end_index]
        del ne_trimmed[mid_index:end_index]
        del Te_trimmed[mid_index:end_index]
        del ni_trimmed[mid_index:end_index]
        del Ti_trimmed[mid_index:end_index]

        ### Fit cubic polynomials through all core and one edge value
        ne_model = numpy.poly1d(numpy.polyfit(x_trimmed, ne_trimmed, 3))
        Te_model = numpy.poly1d(numpy.polyfit(x_trimmed, Te_trimmed, 3))
        ni_model = numpy.poly1d(numpy.polyfit(x_trimmed, ni_trimmed, 3))
        Ti_model = numpy.poly1d(numpy.polyfit(x_trimmed, Ti_trimmed, 3))

        init_ne = ne_model(psi_norm)
        init_Te = Te_model(psi_norm)
        init_ni = ni_model(psi_norm)
        init_Ti = Ti_model(psi_norm)

        init_pressure = (kBoltz * init_ne * init_Te) + (kBoltz * init_ni * init_Ti)

        ### Initialize equilibirum on L-mode-like P' and inductive j_tor profiles
        print('>>> Initializing equilibrium with pedestal removed:')

        init_pp_prof, init_ffp_prof, j_BS = profile_iteration(self,init_pressure,init_ne,init_ni,init_Te,init_Ti,psi_norm,n_psi,Zeff,inductive_jtor,jBS_scale,Zis,include_jBS=False)

        init_pp_prof['y'][-1] = 0. # Enforce 0.0 at edge
        init_ffp_prof['y'][-1] = 0. # Enforce 0.0 at edge

        init_pp_prof['y'] = numpy.nan_to_num(init_pp_prof['y'])
        init_ffp_prof['y'] = numpy.nan_to_num(init_ffp_prof['y'])

        self.set_profiles(ffp_prof=init_ffp_prof,pp_prof=init_pp_prof)

        try:
            self.solve()
            flag = 0
        except ValueError:
            flag = -1
        print('  Solve flag: ', flag)

    ### Specify original H-mode profiles, iterate on bootstrap contribution until reasonably converged
    n = 0
    flag = -1
    print('>>> Iterating on H-mode equilibrium solution:')
    while n < max_iterations:
        print('> Iteration '+str(n)+':')

        pp_prof, ffp_prof, j_BS = profile_iteration(self,pressure,ne,ni,Te,Ti,psi_norm,n_psi,Zeff,inductive_jtor,jBS_scale,Zis)

        pp_prof['y'][-1] = 0. # Enforce 0.0 at edge
        ffp_prof['y'][-1] = 0. # Enforce 0.0 at edge
    
        pp_prof['y'] = numpy.nan_to_num(pp_prof['y']) # Check for any nan's
        ffp_prof['y'] = numpy.nan_to_num(ffp_prof['y']) # Check for any nan's

        self.set_profiles(ffp_prof=ffp_prof,pp_prof=pp_prof)

        try:
            self.solve()
            flag = 0
        except ValueError:
            flag = -1
        print('  Solve flag: ', flag)

        n += 1
        if (n > 2) and (flag >= 0):
            break
        elif n >= max_iterations:
            raise TypeError('H-mode equilibrium solve did not converge')
    
    return flag, j_BS

