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
import copy
import collections
import ctypes
from warnings import warn
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
    settings.limited_only = False
    settings.dipole_mode = False
    settings.mirror_mode = False
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
    elif (profile_dict['type'] == 'linterp') or (profile_dict['type'] == 'jphi-linterp'):
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
        raise KeyError('Invalid profile type ("flat", "linterp", "jphi-linterp")')
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
        ## Internal Grad-Shafranov object (@ref psi_grad_shaf.gs_factory "gs_factory")
        self._tMaker_ptr = c_void_p()
        ## Internal Grad-Shafranov object (@ref psi_grad_shaf.gs_equil "gs_equil")
        self._tMaker_equil = None
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
        self._virtual_coils = {'#VSC': {'id': -1 ,'facs': {}}}
        ## Voltage coils dictionary. Currently only used for plotting on python side.
        self._vcoils = {}
        ## Coil set names in order of id number
        self.coil_set_names = []
        ## Distribution coils, only (currently) saved for plotting utility
        self.dist_coils = {}
        ## Vacuum F value
        self._F0 = 0.0
        ## Normalized flux convention (0 -> tokamak, 1 -> spheromak)
        self.psi_convention = 0
        ## Number of regions in mesh
        self.nregs = -1
        ## Number of coils in mesh
        self.ncoils = -1
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
        ## Limiting contours (if multiple)
        self.lim_contours = None
        ## Coil self-inductance matrix [ncoils]
        self.Lcoils = None
    
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
        if self._tMaker_equil is not None:
            del self._tMaker_equil
        self.nregs = -1
        self.np = -1
        # Reset defaults
        self._tMaker_ptr = c_void_p()
        self._tMaker_equil = None
        self._mesh_ptr = c_void_p()
        self.settings = tokamaker_default_settings(self._oft_env)
        self._cond_dict = {}
        self._vac_dict = {}
        self._coil_dict = {}
        self._vcoils = {}
        self.coil_sets = {}
        self._virtual_coils = {'#VSC': {'id': -1 ,'facs': {}}}
        self._F0 = 0.0
        self.nregs = -1
        self.ncoils = -1
        self.np = -1
        self.nc = -1
        self.r = None
        self.lc = None
        self.reg = None
        self.nvac = 0
        self.lim_contour = None
        self.lim_contours = None
        self.Lcoils = None

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
            coil_dict[key]['name'] = key
            self.coil_sets[coil_set]['sub_coils'].append(coil_dict[key])
            self.coil_sets[coil_set]['sub_coils'][-1]['name'] = key
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
        Lmat_loc = c_double_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_setup(self._tMaker_ptr,order,full_domain,ctypes.byref(ncoils),ctypes.byref(Lmat_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        # Update vacuum flux
        self._F0 = F0
        # Get coil count and reference to coil self-inductance matrix
        self.ncoils = ncoils.value
        self.Lcoils = numpy.ctypeslib.as_array(Lmat_loc,shape=(self.ncoils,self.ncoils))
        # Create equilibirum object
        self._tMaker_equil = TokaMaker_equilibrium(self)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_equil_set(self._tMaker_ptr,self._tMaker_equil.c_ptr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
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
        lim_pts = numpy.ctypeslib.as_array(r_loc,shape=(npts.value, 2))
        self.lim_contours = []
        for i in range(nloops.value):
            lim_contour = numpy.vstack((lim_pts[loop_ptr[i]-1:loop_ptr[i+1]-1,:],lim_pts[loop_ptr[i]-1,:]))
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
    def c_ptr(self):
        r'''C pointer to Fortran-side TokaMaker object'''
        return self._tMaker_ptr

    @property
    def ffp_scale(self):
        r'''! F*F' scale value'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.ffp_scale
    
    # @cond
    @ffp_scale.setter
    def ffp_scale(self,value):
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        self._tMaker_equil.ffp_scale = value
    # @endcond
    
    @property
    def p_scale(self):
        r'''! Pressure scale value'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.p_scale
    
    # @cond
    @p_scale.setter
    def p_scale(self,value):
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        self._tMaker_equil.p_scale = value
    # @endcond

    @property
    def alam(self):
        r'''! F*F' normalization value
        
        @deprecated Use `ffp_scale` property instead.'''
        return self.ffp_scale
    
    # @cond
    @alam.setter
    def alam(self,value):
        self.ffp_scale = value
    # @endcond
    
    @property
    def pnorm(self):
        r'''! Pressure normalization value
        
        @deprecated Use `p_scale` property instead.'''
        return self.p_scale
    
    # @cond
    @pnorm.setter
    def pnorm(self,value):
        self.p_scale = value
    # @endcond
    
    @property
    def diverted(self):
        r'''! Diverted flag (limited if `False`)'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.diverted
    
    @property
    def o_point(self):
        r'''! Location of O-point (magnetic axis) [2]'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.o_point
    
    @property
    def lim_point(self):
        r'''! Limiting point (limter or active X-point) [2]'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.lim_point
    
    @property
    def psi_bounds(self):
        r'''! Bounding values for \f$\psi\f$ (\f$\psi_a\f$,\f$\psi_0\f$) [2]'''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.psi_bounds
    
    def coil_dict2vec(self,coil_dict=None,keep_virtual=False,default_value=0.0):
        '''! Create coil vector from dictionary of values

        @param coil_dict Input dictionary
        @param keep_virtual Keep virtual coils in vector instead of mapping to component coils
        @param default_value Fill value for unspecified entries
        @returns Ouput vector
        '''
        if coil_dict is None:
            coil_dict = {}
        vector = default_value*numpy.ones((self.ncoils+len(self._virtual_coils),))
        removal_vector = numpy.zeros((self.ncoils+len(self._virtual_coils),))
        for coil_key, value in coil_dict.items():
            if coil_key in self.coil_sets:
                vector[self.coil_sets[coil_key]['id']] += value
                removal_vector[self.coil_sets[coil_key]['id']] = default_value
            elif coil_key in self._virtual_coils:
                if keep_virtual:
                    vector[self._virtual_coils[coil_key]['id']] += value
                    removal_vector[self._virtual_coils[coil_key]['id']] = default_value
                else:
                    for map_key, map_val in self._virtual_coils[coil_key].get('facs',{}).items():
                        vector[self.coil_sets[map_key]['id']] += map_val*value
                        removal_vector[self.coil_sets[map_key]['id']] = default_value
            else:
                raise KeyError('Unknown coil "{0}"'.format(coil_key))
        vector -= removal_vector
        if keep_virtual:
            return vector
        else:
            return vector[:self.ncoils]
    
    def coil_vec2dict(self,coil_vec,always_virtual=False):
        '''! Create coil value dictionary of from vector values

        @param coil_vec Input vector
        @param always_virtual Always include virtual coils even if not present in vector
        @returns Ouput dictionary
        '''
        if (coil_vec.shape[0] != self.ncoils) and (coil_vec.shape[0] != self.ncoils+len(self._virtual_coils)):
            raise ValueError('Input vector has incorrect length, should be {0} or {1}'.format(self.ncoils, self.ncoils+len(self._virtual_coils)))
        coil_dict = {}
        for coil_key in self.coil_sets:
            coil_dict[coil_key] = coil_vec[self.coil_sets[coil_key]['id']]
        if coil_vec.shape[0] > self.ncoils:
            for coil_key in self._virtual_coils:
                coil_dict[coil_key] = coil_vec[self._virtual_coils[coil_key]['id']]
        elif always_virtual:
            for coil_key in self._virtual_coils:
                coil_dict[coil_key] = 0.0
        return coil_dict
        
    def abspsi_to_normalized(self,psi_in):
        r'''! Convert unnormalized \f$ \psi \f$ values to normalized \f$ \hat{\psi} \f$ values
        
        @param psi_in Input \f$ \psi \f$ values
        @returns Normalized \f$ \hat{\psi} \f$ values
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.abspsi_to_normalized(psi_in)
    
    def psinorm_to_absolute(self,psi_in):
        r'''! Convert normalized \f$ \hat{\psi} \f$ values to unnormalized values \f$ \psi \f$
        
        @param psi_in Input \f$ \hat{\psi} \f$ values
        @returns Unnormalized \f$ \psi \f$ values
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.psinorm_to_absolute(psi_in)
        
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
                        reg_mat[self._virtual_coils[key]['id'],i] = value
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
        if (not self._virtual_coils.get('#VSC',{}).get('facs',{})) and ((abs(reg_mat[-1,:])).max() < 1.E-8):
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

        @param coil_bounds Minimum and maximum allowable coil currents (dictionary of form `{coil_name: coil_bound[2]}`)
        '''
        bounds_array = numpy.zeros((self.ncoils+len(self._virtual_coils),2), dtype=numpy.float64)
        bounds_array[:,0] = -1.E98
        bounds_array[:,1] = 1.E98
        if coil_bounds is not None:
            for coil_key, coil_bound in coil_bounds.items():
                if coil_key in self.coil_sets:
                    bounds_array[self.coil_sets[coil_key]['id'],:] = coil_bound
                elif coil_key in self._virtual_coils:
                    bounds_array[self._virtual_coils[coil_key]['id'],:] = coil_bound
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
        gains_array = numpy.ascontiguousarray(self.coil_dict2vec(coil_gains), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_vsc(self._tMaker_ptr,gains_array,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self._virtual_coils['#VSC']['facs'] = coil_gains.copy()
    
    def set_vcoils(self,coil_resistances):
        '''! Set or unset one or more coils as Vcoils by defining their lumped resistances

        @param coil_resistances Lumped coil resistances for Vcoils [Ohms] (dictionary of form `{coil_name: coil_res}`)
        '''
        res_array = numpy.ascontiguousarray(self.coil_dict2vec(coil_resistances,default_value=-1.0), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_vcoil(self._tMaker_ptr,res_array,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self._vcoils = copy.deepcopy(coil_resistances)

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
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.load_profiles(f_file,foffset,p_file,eta_file,f_NI_file)

    def set_profiles(self, ffp_prof=None, foffset=None, pp_prof=None, ffp_NI_prof=None, keep_files=False):
        r'''! Set flux function profiles (\f$F*F'\f$ and \f$P'\f$) using a piecewise linear definition

        @param ffp_prof Dictionary object containing FF' profile ['y'] and sampled locations in normalized Psi ['x']
        @param foffset Value of \f$F0=R0*B0\f$
        @param pp_prof Dictionary object containing P' profile ['y'] and sampled locations in normalized Psi ['x']
        @param ffp_NI_prof Dictionary object containing non-inductive FF' profile ['y'] and sampled locations in normalized Psi ['x']
        @param keep_files Retain temporary profile files
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.set_profiles(ffp_prof,foffset,pp_prof,ffp_NI_prof,keep_files)

    def set_resistivity(self, eta_prof=None):
        r'''! Set flux function profile $\eta$ using a piecewise linear definition

        Arrays should have the form array[i,:] = (\f$\hat{\psi}_i\f$, \f$f(\hat{\psi}_i)\f$) and span
        \f$\hat{\psi}_i = [0,1]\f$.

        @param eta_prof Values defining $\eta$ [:,2]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.set_resistivity(eta_prof)

    def solve(self, vacuum=False):
        '''! Solve G-S equation with specified constraints, profiles, etc.
        
        @param vacuum Perform vacuum solve? Plasma-related targets (eg. `Ip`) will be ignored.
        @result Equilibrium object
        '''
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_solve(self._tMaker_ptr,c_bool(vacuum),error_string)
        if error_string.value != b'':
            raise ValueError("Error in solve: {0}".format(error_string.value.decode()))
        return self.copy_eq()
    
    def vac_solve(self,psi=None,rhs_source=None):
        '''! Solve for vacuum solution (no plasma), with present coil currents
        and optional other currents

        @note If isoflux, flux, or saddle constraints are desired use @ref solve instead.
        
        @param psi Boundary values for vacuum solve
        @param rhs_source Current source [A/m^2] (optional)
        @result Equilibrium object
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
        equil_out = self.copy_eq(skip_targets=True,skip_constraints=True)
        equil_out.set_psi(psi)
        return equil_out

    def get_stats(self,lcfs_pad=None,axis_pad=0.02,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Get information (Ip, q, kappa, etc.) about current G-S equilbirium

        See eq. 1 for `li_normalization='std'` and eq 2. for `li_normalization='iter'`
        in [Jackson et al.](https://dx.doi.org/10.1088/0029-5515/48/12/125002)

        @param lcfs_pad Padding at LCFS for boundary calculations (default: 1.0 for limited; 0.99 for diverted)
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        @result Dictionary of equilibrium parameters
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_stats(lcfs_pad,axis_pad,li_normalization,geom_type,beta_Ip)

    def print_info(self,lcfs_pad=None,axis_pad=0.02,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Print information (Ip, q, etc.) about current G-S equilbirium
        
        @param lcfs_pad Padding at LCFS for boundary calculations (default: 1.0 for limited; 0.99 for diverted)
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.print_info(lcfs_pad,axis_pad,li_normalization,geom_type,beta_Ip)
    
    def set_isoflux(self,isoflux,weights=None,grad_wt_lim=-1.0,ref_points=None):
        r'''! Set isoflux constraint points (all points lie on a flux surface)

        To constraint points more uniformly in space additional weighting based on
        the gradient of $\psi$ at each point can also be included by setting
        `grad_wt_lim>0`. When set the actual weight will be
        $w_i * min(grad_wt_lim,|\nabla \psi|_{max} / |\nabla \psi|_i)$

        @deprecated Use `set_isoflux_constraints` instead.

        @param isoflux List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        @param grad_wt_lim Limit on gradient-based weighting (negative to disable)
        @param ref_points Reference points for each isoflux point [:,2] (default: `isoflux[0,:]` is used for all points)
        '''
        warn(
            "`set_isoflux()` is deprecated, use `set_isoflux_constraints()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        self.set_isoflux_constraints(isoflux,weights,grad_wt_lim,ref_points)
    
    def set_isoflux_constraints(self,isoflux,weights=None,grad_wt_lim=-1.0,ref_points=None):
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
            self._tMaker_equil._isoflux_constraints = None
        else:
            if ref_points is None:
                ref_points = numpy.zeros((isoflux.shape[0]-1,2), dtype=numpy.float64)
                ref_points[:,0] = isoflux[0,0]
                ref_points[:,1] = isoflux[0,1]
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
            self._tMaker_equil._isoflux_constraints = isoflux.copy()
    
    def set_flux(self,locations,targets,weights=None):
        r'''! Set explicit flux constraint points \f$ \psi(x_i) \f$ [Wb/rad]

        @deprecated Use `set_psi_constraints` or `set_flux_constraints` instead.

        @param locations List of points defining constraints [:,2]
        @param targets Target \f$ \psi \f$ value in Wb/rad at each point [:]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        warn(
            "`set_flux()` is deprecated, use `set_psi_constraints()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        self.set_psi_constraints(locations,targets,weights)

    def set_flux_constraints(self,locations,targets,weights=None):
        r'''! Set explicit flux constraint points \f$ \psi(x_i) \f$ [Wb]

        @param locations List of points defining constraints [:,2]
        @param targets Target \f$ \psi \f$ value in Wb at each point [:]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        self.set_psi_constraints(locations,targets/(2.0*numpy.pi),weights)
    
    def set_psi_constraints(self,locations,targets,weights=None):
        r'''! Set explicit flux constraint points \f$ \psi(x_i) \f$ [Wb/rad]

        @param locations List of points defining constraints [:,2]
        @param targets Target \f$ \psi \f$ value in Wb/rad at each point [:]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        if locations is None:
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_flux(self._tMaker_ptr,numpy.zeros((1,1)),numpy.zeros((1,)),numpy.zeros((1,)),0,-1.0,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._tMaker_equil._psi_constraints = None
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
            self._tMaker_equil._psi_constraints = (locations.copy(), targets.copy())
    
    def set_saddles(self,saddles,weights=None):
        '''! Set saddle constraint points (poloidal field should vanish at each point)

        @deprecated Use `set_saddle_constraints` instead.

        @param saddles List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        warn(
            "`set_saddles()` is deprecated, use `set_saddle_constraints()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        self.set_saddle_constraints(saddles,weights)

    def set_saddle_constraints(self,saddles,weights=None):
        '''! Set saddle constraint points (poloidal field should vanish at each point)

        @param saddles List of points defining constraints [:,2]
        @param weights Weight to be applied to each constraint point [:] (default: 1)
        '''
        if saddles is None:
            error_string = self._oft_env.get_c_errorbuff()
            tokamaker_set_saddles(self._tMaker_ptr,numpy.zeros((1,1)),numpy.zeros((1,)),0,error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            self._tMaker_equil._saddle_targets = None
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
            self._tMaker_equil._saddle_targets = saddles.copy()
    
    def set_targets(self,Ip=None,Ip_ratio=None,pax=None,estore=None,R0=None,V0=None,Z0=None,retain_previous=False):
        r'''! Set global target values

        @note Values that are not specified are reset to their defaults on each call unless `retain_previous=True`.

        @param Ip Target plasma current [A]
        @param Ip_ratio Amplitude of net plasma current contribution from FF' compared to P'
        @param pax Target axis pressure [Pa]
        @param estore Target sotred energy [J]
        @param R0 Target major radius for magnetic axis
        @param V0 Target vertical position for magnetic axis
        @param Z0 Target vertical position for magnetic axis
        @param retain_previous Keep previously set targets unless explicitly updated? (default: False)
        '''
        def float_to_c(value):
            if value is None:
                return c_double(self._oft_env.float_disable_flag)
            else:
                return c_double(value)
        # Reset all targets unless specified
        if not retain_previous:
            self._tMaker_equil._Ip_target = None
            self._tMaker_equil._estored_target = None
            self._tMaker_equil._pax_target = None
            self._tMaker_equil._Ip_ratio_target = None
            self._tMaker_equil._R0_target = None
            self._tMaker_equil._Z0_target = None
        # Set new targets
        if Ip is not None:
            if (Ip <= 0.0) and (not self._oft_env.float_is_disabled(Ip)):
                raise ValueError("`Ip_target` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._tMaker_equil._Ip_target = copy.copy(Ip)
        if estore is not None:
            if (estore <= 0.0) and (not self._oft_env.float_is_disabled(estore)):
                raise ValueError("`estore` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._tMaker_equil._estored_target = copy.copy(estore)
        if pax is not None:
            if (pax <= 0.0) and (not self._oft_env.float_is_disabled(pax)):
                raise ValueError("`pax` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._tMaker_equil._pax_target = copy.copy(pax)
        if Ip_ratio is not None:
            self._tMaker_equil._Ip_ratio_target = copy.copy(Ip_ratio)
        if R0 is not None:
            if (R0 <= 0.0) and (not self._oft_env.float_is_disabled(R0)):
                raise ValueError("`R0` must be positive or set to `OFT_env.float_disable_flag` to disable")
            self._tMaker_equil._R0_target = copy.copy(R0)
        if V0 is not None:
            Z0 = V0
        if Z0 is not None:
            self._tMaker_equil._Z0_target = copy.copy(Z0)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_targets(self._tMaker_ptr,
                              float_to_c(self._tMaker_equil._Ip_target),
                              float_to_c(self._tMaker_equil._Ip_ratio_target),
                              float_to_c(self._tMaker_equil._pax_target),
                              float_to_c(self._tMaker_equil._estored_target),
                              float_to_c(self._tMaker_equil._R0_target),
                              float_to_c(self._tMaker_equil._Z0_target),
                              error_string
        )
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def get_targets(self):
        r'''! Get global target values

        @result Dictionary of global target values
        '''
        # Get targets
        target_dict = {}
        if self._tMaker_equil.Ip_target is not None:
            target_dict['Ip'] = self._tMaker_equil.Ip_target
        if self._tMaker_equil.Estored_target is not None:
            target_dict['estore'] = self._tMaker_equil.Estored_target
        if self._tMaker_equil.Pax_target is not None:
            target_dict['pax'] = self._tMaker_equil.Pax_target
        if self._tMaker_equil.Ip_ratio_target is not None:
            target_dict['Ip_ratio'] = self._tMaker_equil.Ip_ratio_target
        if self._tMaker_equil.R0_target is not None:
            target_dict['R0'] = self._tMaker_equil.R0_target
        if self._tMaker_equil.Z0_target is not None:
            target_dict['V0'] = self._tMaker_equil.Z0_target
            target_dict['Z0'] = self._tMaker_equil.Z0_target
        return target_dict

    def get_delstar_curr(self,psi):
        r'''! Get toroidal current density from \f$ \psi \f$ through \f$ \Delta^{*} \f$ operator

        @deprecated Use `calc_delstar_curr` instead.
 
        @param psi \f$ \psi \f$ corresponding to desired current density
        @result \f$ J_{\phi} = \textrm{M}^{-1} \Delta^{*} \psi \f$ [A/m^2]
        '''
        warn(
            "`get_delstar_curr()` is deprecated, use `calc_delstar_curr()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.calc_delstar_curr(psi)
    
    def calc_delstar_curr(self,psi):
        r'''! Get toroidal current density from \f$ \psi \f$ through \f$ \Delta^{*} \f$ operator
 
        @param psi \f$ \psi \f$ corresponding to desired current density
        @result \f$ J_{\phi} = \textrm{M}^{-1} \Delta^{*} \psi \f$ [A/m^2]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_delstar_curr(psi)
    
    def get_jtor_plasma(self):
        r'''! Get plasma toroidal current density for current equilibrium

        @deprecated Use `calc_jtor_plasma` instead.
 
        @result \f$ J_{\phi} \f$ by evalutating RHS source terms
        '''
        warn(
            "`get_jtor_plasma()` is deprecated, use `calc_jtor_plasma()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.calc_jtor_plasma()
    
    def calc_jtor_plasma(self):
        r'''! Get plasma toroidal current density for current equilibrium
 
        @result \f$ J_{\phi} \f$ by evalutating RHS source terms
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_jtor_plasma()

    def copy_eq(self,skip_targets=False,skip_constraints=False):
        '''! Create a copy of the current equilibrium object
        
        @param skip_targets When copying, skip copying target values
        @param skip_constraints When copying, skip copying constraint values
        @result New `TokaMaker_equilibrium` object with copied values
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return TokaMaker_equilibrium(source_eq=self._tMaker_equil,skip_targets=skip_targets,skip_constraints=skip_constraints)

    def get_psi(self,normalized=True):
        r'''! Get poloidal flux values on node points

        @param normalized Normalize (and offset) poloidal flux
        @result \f$\hat{\psi} = \frac{\psi-\psi_0}{\psi_a-\psi_0} \f$ or \f$\psi\f$
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_psi(normalized)

    def set_psi(self,psi,update_bounds=False):
        '''! Set poloidal flux values on node points

        @param psi Poloidal flux values (should not be normalized!)
        @param update_bounds Update plasma bounds by determining new limiting points
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.set_psi(psi,update_bounds)
    
    def set_psi_dt(self,psi0,dt,coil_currents=None,coil_voltages=None):
        '''! Set reference poloidal flux and time step for eddy currents in .solve()

        @param psi0 Reference poloidal flux at t-dt (unnormalized)
        @param dt Time since reference poloidal flux
        @param coil_currents Currents for Vcoils [A] (dictionary of form `{coil_name: coil_curr}`, defaults to current solution)
        @param coil_voltages Voltages for Vcoils [V] (dictionary of form `{coil_name: coil_volt}`)
        '''
        if psi0.shape[0] != self.np:
            raise IndexError('Incorrect shape of "psi0", should be [np]')
        psi0 = numpy.ascontiguousarray(psi0, dtype=numpy.float64)
        if coil_currents is None:
            coil_currents, _ = self.get_coil_currents()
        curr_array = numpy.ascontiguousarray(self.coil_dict2vec(coil_currents,keep_virtual=True), dtype=numpy.float64)
        if coil_voltages is not None:
            volt_array = numpy.ascontiguousarray(self.coil_dict2vec(coil_voltages,keep_virtual=True), dtype=numpy.float64)
        else:
            volt_array = numpy.ascontiguousarray(self.coil_dict2vec(None,keep_virtual=True), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_psi_dt(self._tMaker_ptr,psi0,curr_array,volt_array,c_double(dt),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def get_field_eval(self,field_type):
        r'''! Create field interpolator for vector potential

        @param field_type Field to interpolate, must be one of ("B", "psi", "F", "P", "dPSI", "dBr", "dBt", or "dBz")
        @result Field interpolation object
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_field_eval(field_type)
    
    def get_coil_currents(self):
        '''! Get currents in each coil [A] and coil region [A-turns]

        @result Coil currents [ncoils], Coil currents by region [nregs]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_coil_currents()

    def get_coil_Lmat(self):
        r'''! Get mutual inductance matrix between coils

        @note This is the inductance in terms of A-turns. To get in terms of
        current in a single of the \f$n\f$ windings you must multiply by \f$n_i*n_j\f$.

        @result L[ncoils+1,ncoils+1]
        '''
        Lmat = numpy.zeros((self.ncoils+1,self.ncoils+1),dtype=numpy.float64)
        Lmat[:-1,:-1] = self.Lcoils
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        L_p, M_p_c = self._tMaker_equil.calc_inductance()
        Lmat[-1,-1] = L_p
        Lmat[-1,:-1] = M_p_c
        Lmat[:-1,-1] = M_p_c
        return Lmat
    
    def trace_surf(self,psi):
        r'''! Trace surface for a given poloidal flux

        @param psi Flux surface to trace \f$\hat{\psi}\f$
        @result \f$r(\hat{\psi})\f$
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.trace_surf(psi)
    
    def get_q(self,psi=None,psi_pad=0.02,npsi=50,compute_geo=False):
        r'''! Get q-profile at specified or uniformly spaced points

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @param compute_geo Compute geometric values for LCFS
        @result \f$\hat{\psi}\f$, \f$q(\hat{\psi})\f$, \f$[<R>,<1/R>,dV/dPsi]\f$, length of last surface,
        [r(R_min),r(R_max)], [r(z_min),r(z_max)]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_q(psi,psi_pad,npsi,compute_geo)

    def sauter_fc(self,psi=None,psi_pad=0.02,npsi=50):
        r'''! Evaluate Sauter trapped particle fractions at specified or uniformly spaced points

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @result \f$ f_c \f$, [\f$<R>,<1/R>,<a>\f$], [\f$<|B|>,<|B|^2>\f$]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_sauter_fc(psi,psi_pad,npsi)

    def get_globals(self):
        r'''! Get global plasma parameters

        @result Ip, [R_Ip, Z_Ip], \f$\int dV\f$, \f$\int P dV\f$, diamagnetic flux,
        enclosed toroidal flux
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_globals()

    def calc_loopvoltage(self):
        r'''! Get plasma loop voltage

        @result Vloop [Volts]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_loopvoltage()

    def get_profiles(self,psi=None,psi_pad=1.E-8,npsi=50):
        r'''! Get G-S source profiles

        @param psi Explicit sampling locations in \f$\hat{\psi}\f$
        @param psi_pad End padding (axis and edge) for uniform sampling (ignored if `psi` is not None)
        @param npsi Number of points for uniform sampling (ignored if `psi` is not None)
        @result \f$\hat{\psi}\f$, \f$F(\hat{\psi})\f$, \f$F'(\hat{\psi})\f$,
        \f$P(\hat{\psi})\f$, \f$P'(\hat{\psi})\f$
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_profiles(psi,psi_pad,npsi)
    
    def get_xpoints(self):
        '''! Get X-points

        @result X-points, is diverted?
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.get_xpoints()
    
    def set_coil_currents(self, currents=None):
        '''! Set coil currents

        @param currents Current in each coil [A]
        '''
        if currents is None:
            currents = {}
        current_array = numpy.ascontiguousarray(self.coil_dict2vec(currents), dtype=numpy.float64)
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
    
    def set_dipole_a(self,a_exp=None):
        r'''! Update anisotropy exponent `a` when dipole mode is used, calling with no argument
        will disable pressure anisotropy.
        
        @param a_exp New value for `a` exponent
        '''
        error_string = self._oft_env.get_c_errorbuff()
        if a_exp is None:
            a_exp = -1.0
        tokamaker_set_dipole_a(self._tMaker_ptr,a_exp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        
    def set_mirror_slosh(self,n_exp=None,b_turn=None,z_throat=None):
        r'''! Update anisotropy exponent `a` when dipole mode is used, calling with no arguments
        will disable pressure anisotropy.
        
        @param n_exp New value for `n` exponent
        @param b_turn Relative B-field at ion turning point (b_turn = B/B_0)
        @param z_throat Location of mirror throat
        '''
        error_string = self._oft_env.get_c_errorbuff()
        if (n_exp or b_turn or z_throat) and (not (n_exp and b_turn and z_throat)):
            raise ValueError('All (or no) arguments must be passed.')
        if n_exp is None:
            n_exp = -1.0
            b_turn = -1.0
            z_throat = -1.0
        tokamaker_set_mirror_slosh(self._tMaker_ptr,n_exp,b_turn,z_throat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def area_integral(self,field,reg_mask=-1):
        r'''! Compute area integral of field over a specified region

        @deprecated Use `compute_area_integral` instead.

        @param field Field to integrate [np]
        @param reg_mask ID of region for integration (negative for whole mesh)
        @result \f$ \int f dA \f$
        '''
        warn(
            "`area_integral()` is deprecated, use `compute_area_integral()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.compute_area_integral(field,reg_mask)
    
    def compute_area_integral(self,field,reg_mask=-1):
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
        r'''! Compute area integral of flux function over the plasma

        @deprecated Use `compute_flux_integral` instead.

        @param psi_vals \f$ \hat{\psi} \f$ values defining flux function [:]
        @param field_vals Flux function values at each \f$ \hat{\psi} \f$ value [:]
        @result \f$ \int f dA \f$
        '''
        warn(
            "`flux_integral()` is deprecated, use `compute_flux_integral()` instead. This function will be removed in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.compute_flux_integral(psi_vals,field_vals)
    
    def compute_flux_integral(self,psi_vals,field_vals):
        r'''! Compute area integral of flux function over the plasma

        @param psi_vals \f$ \hat{\psi} \f$ values defining flux function [:]
        @param field_vals Flux function values at each \f$ \hat{\psi} \f$ value [:]
        @result \f$ \int f dA \f$
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.compute_flux_integral(psi_vals,field_vals)

    def plot_machine(self,fig,ax,equilibrium=None,vacuum_color='whitesmoke',cond_color='gray',limiter_color='k',
                     coil_color='gray',coil_colormap=None,coil_symmap=False,coil_scale=1.0,coil_clabel=r'$I_C$ [A]',colorbar=None):
        '''! Plot machine geometry

        @param fig Figure to add to
        @param ax Axis to add to
        @param equilibrium Equilibrium object (if `None`, current equilibrium is used)
        @param vacuum_color Color to shade vacuum region (`None` to disable)
        @param cond_color Color for conducting regions (`None` to disable)
        @param limiter_color Color for limiter contour (`None` to disable)
        @param coil_color Color for coil regions (`None` to disable)
        @param coil_colormap Colormap for coil current values
        @param coil_symmap Make coil current colorscale symmetric
        @param coil_scale Scale for coil currents when plotting
        @param coil_clabel Label for coil current colorbar (`None` to disable colorbar)
        @param colorbar Colorbar instance to overwrite (`None` to add)
        @result Colorbar instance for coil colors or `None`
        '''
        # Get equilibrium object if not set
        if equilibrium is None:
            equilibrium = self._tMaker_equil
        mask_vals = numpy.ones((self.np,))
        if self.settings.mirror_mode:
            r_plot = self.r[:,1]
            z_plot = self.r[:,0]
        else:
            r_plot = self.r[:,0]
            z_plot = self.r[:,1]
        # Shade vacuum region
        if vacuum_color is not None:
            mask = numpy.logical_and(self.reg > 1, self.reg <= self.nvac+1)
            if mask.sum() > 0.0:
                ax.tricontourf(r_plot, z_plot, self.lc[mask,:], mask_vals, colors=vacuum_color)
        # Shade coils
        if coil_colormap is not None:
            _, region_currents = equilibrium.get_coil_currents()
            mesh_currents = numpy.zeros((self.lc.shape[0],))
            if self.ncoils >0:
                mesh_currents = region_currents[self.reg-1]
            # Adjust current in coils with non-uniform distribution
            if len(self.dist_coils)>0:
                for _, coil_obj in self.coil_sets.items():
                    if coil_obj["id"] in self.dist_coils.keys():
                        for sub_coil in coil_obj["sub_coils"]:
                            mask = (self.reg==sub_coil["reg_id"])
                            face_currents = numpy.mean(self.dist_coils[coil_obj["id"]][self.lc],axis=1)
                            mesh_currents[mask] *= face_currents[mask]
            mask = (abs(mesh_currents) > 0.0)
            if mask.sum() > 0.0:
                mesh_currents *= coil_scale
                if coil_symmap:
                    max_curr = abs(mesh_currents).max()
                    clf = ax.tripcolor(r_plot, z_plot, self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap, vmin=-max_curr, vmax=max_curr)
                else:
                    clf = ax.tripcolor(r_plot, z_plot, self.lc[mask,:], mesh_currents[mask], cmap=coil_colormap)
                if coil_clabel is not None:
                    cax = None
                    if colorbar is not None:
                        cax = colorbar.ax
                    colorbar = fig.colorbar(clf,ax=ax,cax=cax,label=coil_clabel)
        else:
            for _, coil_reg in self._coil_dict.items():
                mask_tmp = (self.reg == coil_reg['reg_id'])
                ax.tricontourf(r_plot, z_plot, self.lc[mask_tmp,:], mask_vals, colors=coil_color, alpha=1)
        # Shade conductors
        for _, cond_reg in self._cond_dict.items():
            mask_tmp = (self.reg == cond_reg['reg_id'])
            ax.tricontourf(r_plot, z_plot, self.lc[mask_tmp,:], mask_vals, colors=cond_color, alpha=1)
        # Show limiter
        if limiter_color and (self.lim_contour is not None):
            for lim_contour in self.lim_contours:
                if self.settings.mirror_mode:
                    ax.plot(lim_contour[:,1],lim_contour[:,0],color=limiter_color)
                else:
                    ax.plot(lim_contour[:,0],lim_contour[:,1],color=limiter_color)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
        return colorbar

    def plot_mesh(self,fig,ax,lw=0.5,show_legends=True,col_max=10,split_coil_sets=False,plot_tessellated=False):
        '''! Plot computational mesh and regions

        @param fig Figure to add to (unused)
        @param ax Axes to add to (must be scalar, [2], or [2,2])
        @param lw Width of lines in calls to "triplot()"
        @param show_legends Show legends for plots with more than one region?
        @param col_max Maximum number of entries per column in each legend
        @param split_coil_sets Split coil sets into sub-coils when plotting
        @param plot_tessellated Plot mesh tessellated onto FE node points?
        '''
        if not plot_tessellated:
            if not self._mesh_ptr:
                raise ValueError('`plot_mesh()` can only be called after `setup_mesh()`')
            # Get "raw" mesh
            ndim_c = c_int()
            np_c = c_int()
            npc_c = c_int()
            nc_c = c_int()
            nreg_c = c_int()
            r_loc = c_double_ptr()
            lc_loc = c_int_ptr()
            reg_loc = c_int_ptr()
            error_string = self._oft_env.get_c_errorbuff()
            oft_smesh_get(self._mesh_ptr,ctypes.byref(ndim_c),ctypes.byref(np_c),ctypes.byref(r_loc),ctypes.byref(npc_c),
                          ctypes.byref(nc_c),ctypes.byref(lc_loc),ctypes.byref(reg_loc),ctypes.byref(nreg_c),error_string)
            if error_string.value != b'':
                raise Exception(error_string.value)
            rz_plot = numpy.ctypeslib.as_array(r_loc,shape=(np_c.value, 3))
            lc_plot = numpy.ctypeslib.as_array(lc_loc,shape=(nc_c.value, npc_c.value)) - 1
            reg_plot = numpy.ctypeslib.as_array(reg_loc,shape=(nc_c.value,))
        else:
            if self.r is None:
                raise ValueError('`plot_mesh(plot_tessellated=True)` can only be called after `setup()`')
            rz_plot = self.r
            lc_plot = self.lc
            reg_plot = self.reg
        # Get format type from shape of axis object
        format_type = -1
        try:
            if (ax.shape[0] == 2):
                format_type = 1
                try:
                    if (ax.shape[1] == 2):
                        format_type = 2
                    else:
                        format_type = -1
                except IndexError: # `ax` is a 1D array
                    pass
            else:
                format_type = -1
        except AttributeError: # `ax` is not an array
            format_type = 0
        if format_type < 0:
            raise ValueError("Axes for plotting must be scalar, [2], or [2,2]")
        # Set appropriate axes based on format type
        if format_type == 0:
            plasma_axis = ax
            cond_axis = ax
            coil_axis = ax
            vac_axis = ax
            ax_flat = [ax]
        elif format_type == 1:
            plasma_axis = ax[0]
            vac_axis = ax[0]
            cond_axis = ax[1]
            coil_axis = ax[1]
            ax_flat = ax.flatten()
        else:
            plasma_axis = ax[1,0]
            cond_axis = ax[0,1]
            coil_axis = ax[1,1]
            vac_axis = ax[0,0]
            ax_flat = ax.flatten()
        #
        if self.settings.mirror_mode:
            r_plot = rz_plot[:,1]
            z_plot = rz_plot[:,0]
        else:
            r_plot = rz_plot[:,0]
            z_plot = rz_plot[:,1]
        # Get region count
        nregs = reg_plot.max()
        reg_mark = numpy.zeros((nregs,))
        reg_mark[0] = 1
        # Plot the plasma region
        plasma_axis.triplot(r_plot,z_plot,lc_plot[reg_plot==1,:],lw=lw,label='Plasma')
        # Plot conductor regions
        nCond = 0
        for key, cond in self._cond_dict.items():
            if 'vac_id' in cond:
                continue
            nCond += 1
            reg_mark[cond['reg_id']-1] = 1
            cond_axis.triplot(r_plot,z_plot,lc_plot[reg_plot==cond['reg_id'],:],lw=lw,label=key)
        # Plot coil regions
        coil_colors = {}
        for key, coil in self._coil_dict.items():
            reg_mark[coil['reg_id']-1] = 1
            if split_coil_sets:
                leg_key = key
            else:
                leg_key = coil.get('coil_set',key)
            if leg_key not in coil_colors:
                lines, _ = coil_axis.triplot(r_plot,z_plot,lc_plot[reg_plot==coil['reg_id'],:],lw=lw,label=leg_key)
                coil_colors[leg_key] = lines.get_color()
            else:
                coil_axis.triplot(r_plot,z_plot,lc_plot[reg_plot==coil['reg_id'],:],lw=lw,color=coil_colors[leg_key])
        nCoil = len(coil_colors)
        # Plot the vacuum regions
        nVac = 0
        for i in range(nregs):
            if reg_mark[i] == 0:
                nVac += 1
                vac_axis.triplot(r_plot,z_plot,lc_plot[reg_plot==i+1,:],lw=lw,label='Vacuum_{0}'.format(nVac))
        # Format plots
        for ax_tmp in ax_flat:
            ax_tmp.set_aspect('equal','box')
            if self.settings.mirror_mode:
                ax_tmp.set_xlabel('Z (m)')
                ax_tmp.set_ylabel('R (m)')
            else:
                ax_tmp.set_xlabel('R (m)')
                ax_tmp.set_ylabel('Z (m)')
        if show_legends:
            if format_type == 0:
                ncols = max(1,(1+nCond+nCoil+nVac)//col_max)
                plasma_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncol=ncols)
            elif format_type == 1:
                ncols = max(1,(1+nVac)//col_max)
                plasma_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncol=ncols)
                ncols = max(1,(nCond+nCoil)//col_max)
                cond_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncol=ncols)
            elif format_type == 2:
                ncols = max(1,(nCond)//col_max)
                cond_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncol=ncols)
                ncols = max(1,(nCoil)//col_max)
                coil_axis.legend(bbox_to_anchor=(1.05,0.5), loc='center left', ncol=ncols)
    
    def plot_psi(self,fig,ax,equilibrium=None,psi=None,normalized=True,
        plasma_color=None,plasma_nlevels=8,plasma_levels=None,plasma_colormap=None,plasma_linestyles=None,
        vacuum_color='darkgray',vacuum_nlevels=8,vacuum_levels=None,vacuum_colormap=None,vacuum_linestyles=None,
        xpoint_color='k',xpoint_marker='x',xpoint_inactive_alpha=0.5,opoint_color='k',opoint_marker='*'):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param equilibrium Equilibrium object (if `None`, current equilibrium is used)
        @param psi Flux values to plot (otherwise `equilibrium.get_psi()` is called)
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
        @param xpoint_color Color for X-point markers (`None` to disable)
        @param xpoint_marker Marker style for X-points
        @param xpoint_inactive_alpha Alpha value for inactive X-points
        @param opoint_color Color for O-point markers (`None` to disable)
        @param opoint_marker Marker style for O-points
        '''
        if self.settings.mirror_mode:
            r_plot = self.r[:,1]
            z_plot = self.r[:,0]
        else:
            r_plot = self.r[:,0]
            z_plot = self.r[:,1]
        # Get equilibrium object if not set
        if equilibrium is None:
            equilibrium = self._tMaker_equil
        # Plot poloidal flux
        if psi is None:
            psi = equilibrium.get_psi(normalized)
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
            ax.tricontour(r_plot,z_plot,self.lc,psi,levels=vacuum_levels,colors=vacuum_color,cmap=vacuum_colormap,linestyles=vacuum_linestyles)
        if plasma_levels is not None:
            ax.tricontour(r_plot,z_plot,self.lc,psi,levels=plasma_levels,colors=plasma_color,cmap=plasma_colormap,linestyles=plasma_linestyles)

        # Plot saddle points
        if xpoint_color is not None:
            x_points, diverted = equilibrium.get_xpoints()
            if x_points is not None:
                if diverted:
                    ax.plot(x_points[-1,0], x_points[-1,1], color=xpoint_color, marker=xpoint_marker, linestyle='none')
                    ax.plot(x_points[:-1,0], x_points[:-1,1], color=xpoint_color, marker=xpoint_marker, linestyle='none', alpha=xpoint_inactive_alpha)
                else:
                    ax.plot(x_points[:,0], x_points[:,1], color=xpoint_color, marker=xpoint_marker, linestyle='none', alpha=xpoint_inactive_alpha)
        if (opoint_color is not None) and (equilibrium.o_point[0] > 0.0):
            ax.plot(equilibrium.o_point[0], equilibrium.o_point[1], color=opoint_color, marker=opoint_marker)
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
    
    def plot_constraints(self,fig,ax,equilibrium=None,isoflux_color='tab:red',isoflux_marker='+',saddle_color='tab:green',saddle_marker='x'):
        '''! Plot geometry constraints

        @param fig Figure to add to
        @param ax Axis to add to
        @param equilibrium Equilibrium object (if `None`, current equilibrium is used)
        @param isoflux_color Color of isoflux points (`None` to disable)
        @param saddle_color Color of saddle points (`None` to disable)
        '''
        # Get equilibrium object if not set
        if equilibrium is None:
            equilibrium = self._tMaker_equil
        # Plot isoflux constraints
        if (isoflux_color is not None) and (equilibrium.Isoflux_constraints is not None):
            ax.plot(equilibrium.Isoflux_constraints[:,0],equilibrium.Isoflux_constraints[:,1],color=isoflux_color,marker=isoflux_marker,linestyle='none')
        # Plot saddle constraints
        if (saddle_color is not None) and (equilibrium.Saddle_constraints is not None):
            ax.plot(equilibrium.Saddle_constraints[:,0],equilibrium.Saddle_constraints[:,1],color=saddle_color,marker=saddle_marker,linestyle='none')
    
    def plot_eddy(self,fig,ax,psi=None,equilibrium=None,dpsi_dt=None,nlevels=40,colormap='jet',clabel=r'$J_w$ [$A/m^2$]',symmap=False,include_Vcoils=False):
        r'''! Plot contours of \f$\hat{\psi}\f$

        @param fig Figure to add to
        @param ax Axis to add to
        @param psi Psi corresponding to eddy currents (eg. from time-dependent simulation)
        @param equilibrium Equilibrium object (if `equilibrium=None`, `psi=None` and `dpsi_dt=None`, current equilibrium is used)
        @param dpsi_dt dPsi/dt source eddy currents (eg. from linear stability)
        @param nlevels Number contour lines used for shading (with "psi" only)
        @param colormap Colormap to use for shadings
        @param clabel Label for colorbar (`None` to disable colorbar)
        @param symmap Make colorscale symmetric around zero (only applies if `colormap` is specified)
        @param include_Vcoils Include V-coil regions when plotting eddy currents? (only applies if `dpsi_dt=None` is specified)
        @result Colorbar object
        '''
        if self.settings.mirror_mode:
            r_plot = self.r[:,1]
            z_plot = self.r[:,0]
        else:
            r_plot = self.r[:,0]
            z_plot = self.r[:,1]
        # Get equilibrium object if not set
        if dpsi_dt is not None:
            if (psi is not None) or (equilibrium is not None):
                raise ValueError('Only one of "psi", "equilibrium" or "dpsi_dt" can be specified')
            mask, plot_field = self.get_conductor_source(dpsi_dt)
        else:
            if psi is None:
                if equilibrium is None:
                    equilibrium = self._tMaker_equil
                psi = equilibrium.get_psi(False)
            else:
                if equilibrium is not None:
                    raise ValueError('Only one of "psi", "equilibrium" or "dpsi_dt" can be specified')
        if psi is not None:
            mask, plot_field = self.get_conductor_currents(psi,cell_centered=(nlevels < 0),include_Vcoils=include_Vcoils)
        elif dpsi_dt is not None:
            mask, plot_field = self.get_conductor_source(dpsi_dt)
        if mask.sum() == 0:
            print("Warning: No conducting regions to plot")
            return None
        if plot_field.shape[0] == self.nc:
            if symmap:
                max_curr = abs(plot_field).max()
                clf = ax.tripcolor(r_plot,z_plot,self.lc[mask,:],plot_field[mask],cmap=colormap,vmin=-max_curr,vmax=max_curr)
            else:
                clf = ax.tripcolor(r_plot,z_plot,self.lc[mask],plot_field[mask],cmap=colormap)
        else:
            if symmap:
                max_curr = abs(plot_field[self.lc[mask,:]]).max(axis=None)
                clf = ax.tricontourf(r_plot,z_plot,self.lc[mask],plot_field,nlevels,cmap=colormap,vmin=-max_curr,vmax=max_curr)
            else:
                clf = ax.tricontourf(r_plot,z_plot,self.lc[mask],plot_field,nlevels,cmap=colormap)
        if clabel is not None:
            cb = fig.colorbar(clf,ax=ax)
            cb.set_label(clabel)
        else:
            cb = None
        # Make 1:1 aspect ratio
        ax.set_aspect('equal','box')
        return cb
    
    def get_conductor_currents(self,psi,cell_centered=False,include_Vcoils=False):
        r'''! Get toroidal current density in conducting regions for a given \f$ \psi \f$

        @param psi Psi corresponding to field with conductor currents (eg. from time-dependent simulation)
        @param cell_centered Get currents at cell centers
        @param include_Vcoils Include voltage coils in the calculation?
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_conductor_currents(psi,cell_centered,include_Vcoils)
    
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
                field_tmp = -curr/eta
                mesh_currents[mask_tmp] = numpy.sum(field_tmp[self.lc[mask_tmp,:]],axis=1)/3.0
                if cond_reg.get('noncontinuous',False):
                    mesh_currents[mask_tmp] -= (mesh_currents[mask_tmp]*area[mask_tmp]).sum()/area[mask_tmp].sum()
                mask = numpy.logical_or(mask,mask_tmp)
        return mask, mesh_currents

    def get_vfixed(self):
        '''! Get required vacuum flux values to balance fixed boundary equilibrium

        @result sampling points [:,2], flux values [:]
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        return self._tMaker_equil.calc_vfixed()

    def save_eqdsk(self,filename,nr=65,nz=65,rbounds=None,zbounds=None,run_info='',lcfs_pad=0.01,rcentr=None,truncate_eq=True,limiter_file='',lcfs_pressure=0.0, cocos=7):
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
        @param cocos COCOS version. (Only 2 or 7 supported. `cocos=7` is the default.)
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        self._tMaker_equil.save_eqdsk(filename,nr,nz,rbounds,zbounds,run_info,lcfs_pad,rcentr,truncate_eq,limiter_file,lcfs_pressure,cocos)
    
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
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        self._tMaker_equil.save_ifile(filename,npsi,ntheta,lcfs_pad,lcfs_pressure,pack_lcfs,single_precision)
    
    def save_mug(self,filename):
        r'''! Save current equilibrium to MUG transfer format

        @param filename Filename to save equilibrium to
        '''
        if self._tMaker_equil is None:
            raise ValueError("Equilibrium object is `None`")
        self._tMaker_equil.save_mug(filename)

    def set_coil_current_dist(self,coil_name,curr_dist=None,normalize=False):
        '''! Overwrite coil with non-uniform current distribution.

        @param coil_name Name of coil to modify
        @param curr_dist Relative current density [self.np] (None to disable non-uniform distribution and return to uniform current)
        @param normalize Normalize distribution to have unit current?
        '''
        if coil_name not in self.coil_sets:
            raise KeyError('Unknown coil "{0}"'.format(coil_name))
        iCoil = self.coil_sets[coil_name]['id']
        if curr_dist is None:
            curr_dist = numpy.ones((self.np,), dtype=numpy.float64)
            iCoil_c = c_int(-(iCoil+1))
        else:
            if curr_dist.shape != (self.np,):
                raise ValueError('curr_dist must be the same shape as the number of points in the mesh ({0})'.format(self.np))
            curr_dist = numpy.ascontiguousarray(curr_dist, dtype=numpy.float64)
            iCoil_c = c_int(iCoil+1)
        dist_coil_ptr = c_double_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_coil_current_dist(self._tMaker_ptr,iCoil_c,curr_dist,ctypes.byref(dist_coil_ptr),c_bool(normalize),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        # Update python side coil distribution if successful
        if iCoil_c.value > 0:
            self.dist_coils[iCoil] = numpy.ctypeslib.as_array(dist_coil_ptr,shape=(self.np,))
        else:
            self.dist_coils.pop(iCoil,None)

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

    def step_td(self,time,dt,coil_currents=None,coil_voltages=None):
        '''! Advance time-dependent solution by one time step

        @param time Time at the start of the time step [s]
        @param dt Timestep size [s]
        @param coil_currents Coil currents at the end of the timestep (if `None`, current coil currents are used)
        @param coil_voltages Coil voltages to apply over the timestep (if `None`, zero voltages are assumed)
        @result new time, new dt, # of NL iterations, # of linear iterations, # of retries
        '''
        dt = c_double(dt)
        time = c_double(time)
        nl_its = c_int()
        lin_its = c_int()
        nretry = c_int()
        error_string = self._oft_env.get_c_errorbuff()
        if coil_currents is None:
            coil_currents, _ = self.get_coil_currents()
        coil_currents = numpy.ascontiguousarray(self.coil_dict2vec(coil_currents), dtype=numpy.float64)
        if coil_voltages is None:
            coil_voltages = numpy.zeros((self.ncoils,), dtype=numpy.float64)
        else:
            coil_voltages = numpy.ascontiguousarray(self.coil_dict2vec(coil_voltages), dtype=numpy.float64)
        tokamaker_step_td(self._tMaker_ptr,coil_currents,coil_voltages,ctypes.byref(time),ctypes.byref(dt),
                          ctypes.byref(nl_its),ctypes.byref(lin_its),ctypes.byref(nretry),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return time.value, dt.value, nl_its.value, lin_its.value, nretry.value


class TokaMaker_equilibrium():
    '''! TokaMaker G-S solver class'''
    def __init__(self,TokaMaker_obj=None,source_eq=None,skip_targets=False,skip_constraints=False):
        '''! Initialize TokaMaker object

        @param TokaMaker_obj TokaMaker object (See @ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "TokaMaker")
        @param source_eq Existing equilibrium object to copy
        @param skip_targets When copying from `source_eq`, skip copying target values
        @param skip_constraints When copying from `source_eq`, skip copying constraint values
        '''
        if source_eq is None:
            source_ptr = c_void_p(None)
            if TokaMaker_obj is None:
                raise ValueError('Either "TokaMaker_obj" or "source_eq" must be provided')
        else:
            if TokaMaker_obj is not None:
                raise ValueError('"TokaMaker_obj" and "source_eq" cannot be provided together')
            TokaMaker_obj = source_eq._tMaker
            source_ptr = source_eq._equil_ptr
        ## Reference to parent TokaMaker object
        self._tMaker = TokaMaker_obj
        ## Reference to OpenFUSIONToolkit environment (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        self._oft_env = self._tMaker._oft_env
        ## Internal Grad-Shafranov object (@ref psi_grad_shaf.gs_equil "gs_equil")
        self._equil_ptr = c_void_p()
        if source_eq is None:
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.F0 "F0" property)
            self._F0 = copy.copy(self._tMaker._F0)
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Ip_target "Ip_target" property)
            self._Ip_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Ip_ratio_target "Ip_ratio_target" property)
            self._Ip_ratio_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Pax_target "Pax_target" property)
            self._pax_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Estored_target "Estored_target" property)
            self._estored_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.R0_target "R0_target" property)
            self._R0_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Z0_target "Z0_target" property)
            self._Z0_target = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Isoflux_constraints "Isoflux_constraints" property)
            self._isoflux_constraints = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Psi_constraints "Psi_constraints" property)
            self._psi_constraints = None
            ## Internal value (use @ref TokaMaker.TokaMaker_equilibrium.Saddle_constraints "Saddle_constraints" property)
            self._saddle_targets = None
        else:
            self._F0 = copy.copy(source_eq._F0)
            if not skip_targets:
                self._Ip_target = source_eq._Ip_target
                self._Ip_ratio_target = source_eq._Ip_ratio_target
                self._pax_target = source_eq._pax_target
                self._estored_target = source_eq._estored_target
                self._R0_target = source_eq._R0_target
                self._Z0_target = source_eq._Z0_target
            if not skip_constraints:
                self._isoflux_constraints = source_eq._isoflux_constraints.copy() if source_eq._isoflux_constraints is not None else None
                self._saddle_targets = source_eq._saddle_targets.copy() if source_eq._saddle_targets is not None else None
            self._psi_constraints = (source_eq._psi_constraints[0].copy(), source_eq._psi_constraints[1].copy()) if source_eq._psi_constraints is not None else None
        ## Normalized flux convention (0 -> tokamak, 1 -> spheromak)
        self.psi_convention = self._tMaker.psi_convention
        ## Free or fixed boundary flag
        self.free_boundary = bool(self._tMaker.settings.free_boundary)

        # Create equilibirum object
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_equil_copy(self._tMaker.c_ptr,source_ptr,ctypes.byref(self._equil_ptr),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        
        # Get references to internal variables
        o_loc = c_double_ptr()
        lim_loc = c_double_ptr()
        x_loc = c_double_ptr()
        div_flag_loc = c_bool_ptr()
        bounds_loc = c_double_ptr()
        ffp_scale_loc = c_double_ptr()
        p_scale_loc = c_double_ptr()
        has_plasma_loc = c_bool_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_refs(self._equil_ptr,ctypes.byref(o_loc),ctypes.byref(lim_loc),ctypes.byref(x_loc),ctypes.byref(div_flag_loc),
                    ctypes.byref(bounds_loc),ctypes.byref(ffp_scale_loc),ctypes.byref(p_scale_loc),ctypes.byref(has_plasma_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        ## F*F' scale value [1] (use @ref TokaMaker.TokaMaker.ffp_scale "ffp_scale" property)
        self._ffp_scale = numpy.ctypeslib.as_array(ffp_scale_loc,shape=(1,))
        ## Pressure scale value [1] (use @ref TokaMaker.TokaMaker.p_scale "p_scale" property)
        self._p_scale = numpy.ctypeslib.as_array(p_scale_loc,shape=(1,))
        ## Location of O-point (magnetic axis) [2]
        self._o_point = numpy.ctypeslib.as_array(o_loc,shape=(2,))
        ## Limiting point (limter or active X-point) [2]
        self._lim_point = numpy.ctypeslib.as_array(lim_loc,shape=(2,))
        ## Location of X-points [20,2]
        self._x_points = numpy.ctypeslib.as_array(x_loc,shape=(20, 2))
        ## Diverted (limited) flag [1] (use @ref TokaMaker.TokaMaker.diverted "diverted" property)
        self._diverted = numpy.ctypeslib.as_array(div_flag_loc,shape=(1,))
        ## Bounding values for \f$\psi\f$ (\f$\psi_a\f$,\f$\psi_0\f$) [2]
        self._psi_bounds = numpy.ctypeslib.as_array(bounds_loc,shape=(2,))
        ## Bounding values for \f$\psi\f$ (\f$\psi_a\f$,\f$\psi_0\f$) [2]
        self._has_plasma = numpy.ctypeslib.as_array(has_plasma_loc,shape=(1,))

        if self.Vacuum:
            self.ffp_scale = 0.0
            self.p_scale = 0.0
            self._o_point[:] = [-1.0, 0.0]
            self._lim_point[:] = [-1.0, 0.0]
            for i in range(self._x_points.shape[0]):
                self._x_points[i,:] = [-1.0, 0.0]
            self._diverted[0] = False
            self._psi_bounds[:] = [-1.E99, 1.E99]
        else:
            if source_eq is None:
                self.ffp_scale = 0.1
                self.p_scale = 0.1
                default_prof={
                    'type': 'linterp',
                    'x': numpy.array([0.0,1.0]),
                    'y': numpy.array([1.0,0.0])
                }
                self.set_profiles(ffp_prof=default_prof, pp_prof=default_prof)
    
    def __del__(self):
        '''! Free Fortran-side objects by calling `reset()` before object is deleted or GC'd'''
        if not self._equil_ptr:
            return # Nothing to do
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_equil_destroy(self._equil_ptr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def copy(self):
        r'''! Create a copy of the equilibrium object'''
        return TokaMaker_equilibrium(source_eq=self)
        
    @property
    def c_ptr(self):
        r'''C pointer to Fortran-side equilibrium object'''
        return self._equil_ptr
    
    @property
    def ffp_scale(self):
        r'''! F*F' scale value'''
        return self._ffp_scale[0]
    
    # @cond
    @ffp_scale.setter
    def ffp_scale(self,value):
        self._ffp_scale[0] = value
    # endcond
    
    @property
    def p_scale(self):
        r'''! Pressure scale value'''
        return self._p_scale[0]
    
    # @cond
    @p_scale.setter
    def p_scale(self,value):
        self._p_scale[0] = value
    # endcond
    
    @property
    def alam(self):
        r'''! F*F' scale value
        
        @deprecated Use `ffp_scale` property instead.'''
        return self.ffp_scale
    
    # @cond
    @alam.setter
    def alam(self,value):
        self.ffp_scale = value
    # endcond

    @property
    def pnorm(self):
        r'''! Pressure scale value
        
        @deprecated Use `p_scale` property instead.'''
        return self.p_scale
    
    # @cond
    @pnorm.setter
    def pnorm(self,value):
        self.p_scale = value
    # @endcond
    
    @property
    def diverted(self):
        r'''! Diverted flag (limited if `False`)'''
        return self._diverted[0]
        
    @property
    def o_point(self):
        r'''! Location of O-point (magnetic axis) [2]'''
        return self._o_point
    
    @property
    def lim_point(self):
        r'''! Limiting point (limter or active X-point) [2]'''
        return self._lim_point
    
    @property
    def psi_bounds(self):
        r'''! Bounding values for \f$\psi\f$ (\f$\psi_a\f$,\f$\psi_0\f$) [2]'''
        return self._psi_bounds

    @property
    def Vacuum(self):
        r'''! Equilibrium is vacuum (no plasma contributions)?'''
        return not self._has_plasma[0]

    @property
    def F0(self):
        r'''! Vacuum \f$ F = B_t(R_0) * R_0 \f$'''
        return self._F0

    @property
    def Ip_target(self):
        r'''! Plasma current target'''
        return self._Ip_target
    
    @property
    def Ip_ratio_target(self):
        r'''! Plasma current ratio I_p(FF') / I_p(P') target'''
        return self._Ip_ratio_target

    @property
    def Pax_target(self):
        r'''! Axis pressure target'''
        return self._pax_target
    
    @property
    def Estored_target(self):
        r'''! Stored energy target'''
        return self._estored_target
    
    @property
    def R0_target(self):
        r'''! Magnetic axis radial position target'''
        return self._R0_target.value
    
    @property
    def Z0_target(self):
        r'''! Magnetic axis vertical position target'''
        return self._Z0_target.value
    
    @property
    def Isoflux_constraints(self):
        r'''! Isoflux constraint points'''
        return self._isoflux_constraints

    @property
    def Psi_constraints(self):
        r'''! Flux constraint points and target values'''
        if self._psi_constraints is None:
            return None, None
        return self._psi_constraints[0], self._psi_constraints[1]

    @property
    def Saddle_constraints(self):
        r'''! Saddle constraint points'''
        return self._saddle_targets
    
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
        tokamaker_load_profiles(self.c_ptr,f_file_c,c_double(self._F0),p_file_c,eta_file_c,f_NI_file_c,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def set_profiles(self, ffp_prof=None, foffset=None, pp_prof=None, ffp_NI_prof=None, keep_files=False):
        r'''! Set flux function profiles (\f$F*F'\f$ and \f$P'\f$) using a piecewise linear definition

        @param ffp_prof Dictionary object containing FF' profile ['y'] and sampled locations in normalized Psi ['x']
        @param foffset Value of \f$F0=R0*B0\f$
        @param pp_prof Dictionary object containing P' profile ['y'] and sampled locations in normalized Psi ['x']
        @param ffp_NI_prof Dictionary object containing non-inductive FF' profile ['y'] and sampled locations in normalized Psi ['x']
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
        ffp_NI_file = 'none'
        if ffp_NI_prof is not None:
            ffp_NI_file = self._oft_env.unique_tmpfile('tokamaker_ffp_NI.prof')
            create_prof_file(self, ffp_NI_file, ffp_NI_prof, "ffp_NI")
            delete_files.append(ffp_NI_file)
        self.load_profiles(f_file=ffp_file,foffset=foffset,p_file=pp_file,f_NI_file=ffp_NI_file)
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
        eta_file = 'none'
        if eta_prof is not None:
            eta_file = 'tokamaker_eta.prof'
            create_prof_file(self, eta_file, eta_prof, "eta")
        self.load_profiles(eta_file=eta_file)

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
        
    def get_psi(self,normalized=True):
        r'''! Get poloidal flux values on node points

        @param normalized Normalize (and offset) poloidal flux
        @result \f$\hat{\psi} = \frac{\psi-\psi_0}{\psi_a-\psi_0} \f$ or \f$\psi\f$
        '''
        psi = numpy.zeros((self._tMaker.np,),dtype=numpy.float64)
        psi_lim = c_double()
        psi_max = c_double()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_psi(self._equil_ptr,psi,ctypes.byref(psi_lim),ctypes.byref(psi_max),error_string)
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
        if psi.shape[0] != self._tMaker.np:
            raise IndexError('Incorrect shape of "psi", should be [np]')
        psi = numpy.ascontiguousarray(psi, dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_set_psi(self.c_ptr,psi,c_bool(update_bounds),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def get_coil_currents(self):
        '''! Get currents in each coil [A] and coil region [A-turns]

        @result Coil currents [ncoils], Coil currents by region [nregs]
        '''
        currents = numpy.zeros((self._tMaker.ncoils,),dtype=numpy.float64)
        currents_reg = numpy.zeros((self._tMaker.nregs,),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_coil_currents(self._equil_ptr,currents,currents_reg,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return self._tMaker.coil_vec2dict(currents), currents_reg

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
        tokamaker_get_globals(self._equil_ptr,ctypes.byref(Ip),centroid,ctypes.byref(vol),ctypes.byref(pvol),
            ctypes.byref(dflux),ctypes.byref(tflux),ctypes.byref(Bp_vol),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Ip.value, centroid, vol.value, pvol.value, dflux.value, tflux.value, Bp_vol.value
    
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
        tokamaker_get_profs(self._equil_ptr,psi.shape[0],psi,f,fp,p,pp,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        #
        if self.psi_convention == 0:
            return psi_save,f,fp,p/mu0,pp/mu0
        else:
            return psi,f,fp,p/mu0,pp/mu0
    
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
        tokamaker_get_q(self._equil_ptr,psi.shape[0],psi,qvals,ravgs,ctypes.byref(dl),rbounds,zbounds,error_string)
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
        tokamaker_trace_surf(self._equil_ptr,c_double(psi),ctypes.byref(points_loc),ctypes.byref(npoints),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        if npoints.value > 0:
            return numpy.ctypeslib.as_array(points_loc,shape=(npoints.value, 2))
        else:
            return None
    
    def get_xpoints(self):
        '''! Get X-points

        @result X-points, is diverted?
        '''
        if self._x_points[0,0] < 0.0:
            return None, False
        else:
            for i in range(self._x_points.shape[0]):
                if self._x_points[i,0] < 0.0:
                    break
            return self._x_points[:i,:], self.diverted
    
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
        tokamaker_get_field_eval(self._equil_ptr,imode,ctypes.byref(int_obj),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        field_dim = 1
        if imode == 1:
            field_dim = 3
        elif imode >= 5:
            field_dim = 2
        return TokaMaker_field_interpolator(self._equil_ptr,int_obj,imode,field_dim)

    def get_stats(self,lcfs_pad=None,axis_pad=0.02,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Get information (Ip, q, kappa, etc.) about current G-S equilbirium

        See eq. 1 for `li_normalization='std'` and eq 2. for `li_normalization='iter'`
        in [Jackson et al.](https://dx.doi.org/10.1088/0029-5515/48/12/125002)

        @param lcfs_pad Padding at LCFS for boundary calculations (default: 1.0 for limited; 0.99 for diverted)
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        @result Dictionary of equilibrium parameters
        '''
        Ip,centroid,vol,pvol,dflux,tflux,Bp_vol = self.get_globals()
        if beta_Ip is not None:
            Ip = beta_Ip
        p_psi = numpy.linspace(0.0,1.0,100)
        p_psi[0] = 0.001
        _,_,_,p,_ = self.get_profiles(p_psi)
        # Return reduced information for mirrors
        if self._tMaker.settings.mirror_mode:
            return {
                'Ip': Ip,
                'Ip_centroid': centroid,
                'vol': vol,
                'P_ax': p[0],
                'P_max': p.max(),
                'W_MHD': pvol*1.5,
                'dflux': dflux,
                'tflux': tflux
            }
        # Trace 95% surface and core/edge
        if lcfs_pad is None:
            lcfs_pad = 0.0
            if self.diverted or (not self.free_boundary):
                lcfs_pad = 0.01
        _,qvals,_,dl,rbounds,zbounds = self.get_q(numpy.r_[1.0-lcfs_pad,0.95,axis_pad],compute_geo=True) # Given backward so last point is LCFS (for dl)
        # Get diverted topology information
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

    def print_info(self,lcfs_pad=None,axis_pad=0.02,li_normalization='std',geom_type='max',beta_Ip=None):
        r'''! Print information (Ip, q, etc.) about current G-S equilbirium
        
        @param lcfs_pad Padding at LCFS for boundary calculations (default: 1.0 for limited; 0.99 for diverted)
        @param li_normalization Form of normalized \f$ l_i \f$ ('std', 'ITER')
        @param geom_type Method for computing geometric major/minor radius ('max': Use LCFS extrema, 'mid': Use axis plane extrema)
        @param beta_Ip Override \f$ I_p \f$ used for beta calculations
        '''
        eq_stats = self.get_stats(lcfs_pad=lcfs_pad,axis_pad=axis_pad,li_normalization=li_normalization,geom_type=geom_type,beta_Ip=beta_Ip)
        print("Equilibrium Statistics:")
        if self.diverted:
            print("  Topology                =   Diverted")
        else:
            print("  Topology                =   Limited")
        print("  Toroidal Current [A]    =   {0:11.4E}".format(eq_stats['Ip']))
        print("  Current Centroid [m]    =   {0:6.3F} {1:6.3F}".format(*eq_stats['Ip_centroid']))
        if self._tMaker.settings.dipole_mode:
            print("  Inner limiter [m]       =   {0:6.3F} {1:6.3F}".format(*self.o_point))
        elif not self._tMaker.settings.mirror_mode:
            print("  Magnetic Axis [m]       =   {0:6.3F} {1:6.3F}".format(*self.o_point))
        if not self._tMaker.settings.mirror_mode:
            print("  Elongation              =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['kappa'],eq_stats['kappaU'],eq_stats['kappaL']))
            print("  Triangularity           =   {0:6.3F} (U: {1:6.3F}, L: {2:6.3F})".format(eq_stats['delta'],eq_stats['deltaU'],eq_stats['deltaL']))
        print("  Plasma Volume [m^3]     =   {0:6.3F}".format(eq_stats['vol']))
        if not (self._tMaker.settings.dipole_mode or self._tMaker.settings.mirror_mode):
            print("  q_0, q_95               =   {0:6.3F} {1:6.3F}".format(eq_stats['q_0'],eq_stats['q_95']))
        if self._tMaker.settings.dipole_mode or self._tMaker.settings.mirror_mode:
            print("  Peak Pressure [Pa]      =   {0:11.4E}".format(eq_stats['P_max']))
        else:
            print("  Plasma Pressure [Pa]    =   Axis: {0:11.4E}, Peak: {1:11.4E}".format(eq_stats['P_ax'], eq_stats['P_max']))
        print("  Stored Energy [J]       =   {0:11.4E}".format(eq_stats['W_MHD'])) if 'W_MHD' in eq_stats else None
        print("  <Beta_pol> [%]          =   {0:7.4F}".format(eq_stats['beta_pol'])) if 'beta_pol' in eq_stats else None
        print("  <Beta_tor> [%]          =   {0:7.4F}".format(eq_stats['beta_tor'])) if 'beta_tor' in eq_stats else None
        print("  <Beta_n>   [%]          =   {0:7.4F}".format(eq_stats['beta_n'])) if 'beta_n' in eq_stats else None
        print("  Diamagnetic flux [Wb]   =   {0:11.4E}".format(eq_stats['dflux'])) if 'dflux' in eq_stats else None
        print("  Toroidal flux [Wb]      =   {0:11.4E}".format(eq_stats['tflux'])) if 'tflux' in eq_stats else None
        print("  l_i                     =   {0:7.4F}".format(eq_stats['l_i'])) if 'l_i' in eq_stats else None
    
    def calc_loopvoltage(self):
        r'''! Get plasma loop voltage

        @result Vloop [Volts]
        '''
        V_loop = c_double()
        #
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_gs_calc_vloop(self._equil_ptr,ctypes.byref(V_loop),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        #
        if V_loop.value < 0.:
            raise ValueError('eta array not specified')
        else:
            return V_loop.value
    
    def calc_sauter_fc(self,psi=None,psi_pad=0.02,npsi=50):
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
        tokamaker_sauter_fc(self._equil_ptr,psi.shape[0],psi,fc,r_avgs,modb_avgs,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        if self.psi_convention == 0:
            return psi_save,fc,r_avgs,modb_avgs
        else:
            return psi,fc,r_avgs,modb_avgs
    
    def calc_delstar_curr(self,psi):
        r'''! Get toroidal current density from \f$ \psi \f$ through \f$ \Delta^{*} \f$ operator
 
        @param psi \f$ \psi \f$ corresponding to desired current density
        @result \f$ J_{\phi} = \textrm{M}^{-1} \Delta^{*} \psi \f$ [A/m^2]
        '''
        curr = numpy.copy(psi)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_dels_curr(self._equil_ptr,curr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return curr/mu0
    
    def calc_jtor_plasma(self):
        r'''! Get plasma toroidal current density for current equilibrium
 
        @result \f$ J_{\phi} \f$ by evalutating RHS source terms [A/m^2]
        '''
        curr = numpy.zeros((self._tMaker.np,), dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_jtor(self._equil_ptr,curr,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return curr/mu0

    def calc_conductor_currents(self,psi,cell_centered=False,include_Vcoils=False):
        r'''! Get toroidal current density in conducting regions for a given \f$ \psi \f$

        @param psi Psi corresponding to field with conductor currents (eg. from time-dependent simulation)
        @param cell_centered Get currents at cell centers
        @param include_Vcoils Include voltage coils in the calculation?
        '''
        curr = self.calc_delstar_curr(psi)
        if cell_centered:
            mesh_currents = numpy.zeros((self._tMaker.lc.shape[0],))
        # Loop over conducting regions and get mask/fields
        mask = numpy.zeros((self._tMaker.lc.shape[0],), dtype=numpy.int32)
        for _, cond_reg in self._tMaker._cond_dict.items():
            eta = cond_reg.get('eta',-1.0)
            if eta > 0:
                mask_tmp = (self._tMaker.reg == cond_reg['reg_id'])
                if cell_centered:
                    mesh_currents[mask_tmp] = numpy.sum(curr[self._tMaker.lc[mask_tmp,:]],axis=1)/3.0
                mask = numpy.logical_or(mask,mask_tmp)
        # Treat vcoils as conductors when looking at induced currents
        if include_Vcoils:
            for coil_name, coil_obj in self._tMaker.coil_sets.items():
                if coil_name in self._tMaker._vcoils.keys():
                    for sub_coil in coil_obj["sub_coils"]:
                        mask_tmp = (self._tMaker.reg == sub_coil['reg_id'])
                        if cell_centered:
                            mesh_currents[mask_tmp] = numpy.mean(curr[self._tMaker.lc[mask_tmp]],axis=1)
                        mask = numpy.logical_or(mask,mask_tmp)
        if cell_centered:
            return mask, mesh_currents
        else:
            return mask, curr
    
    def calc_vfixed(self):
        '''! Get required vacuum flux values to balance fixed boundary equilibrium

        @result sampling points [:,2], flux values [:]
        '''
        npts = c_int()
        pts_loc = c_double_ptr()
        flux_loc = c_double_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_vfixed(self._equil_ptr,ctypes.byref(npts),ctypes.byref(pts_loc),ctypes.byref(flux_loc),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return numpy.ctypeslib.as_array(pts_loc,shape=(npts.value, 2)), \
            numpy.ctypeslib.as_array(flux_loc,shape=(npts.value,))
    
    def calc_inductance(self):
        r'''! Get mutual inductance matrix between coils

        @note This is the inductance in terms of A-turns. To get in terms of
        current in a single of the \f$n\f$ windings you must multiply by \f$n_i*n_j\f$.

        @result L[ncoils+1,ncoils+1]
        '''
        Lmat = numpy.zeros((self._tMaker.ncoils+1,),dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_get_plasma_Lmat(self._equil_ptr,Lmat,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return Lmat[-1], Lmat[:-1]

    def compute_area_integral(self,field,reg_mask=-1):
        r'''! Compute area integral of field over a specified region

        @param field Field to integrate [np,]
        @param reg_mask ID of region for integration (negative for whole mesh)
        @result \f$ \int f dA \f$
        '''
        return self._tMaker.compute_area_integral(field,reg_mask)
    
    def compute_flux_integral(self,psi_vals,field_vals):
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
        tokamaker_flux_int(self._equil_ptr,psi_vals,field_vals,c_int(psi_vals.shape[0]),ctypes.byref(result),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return result.value

    def save_eqdsk(self,filename,nr=65,nz=65,rbounds=None,zbounds=None,run_info='',lcfs_pad=0.01,rcentr=None,truncate_eq=True,limiter_file='',lcfs_pressure=0.0, cocos=7):
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
        @param cocos COCOS version. (Only 2 or 7 supported. COCOS=7 is the default.)
        '''
        cfilename = self._oft_env.path2c(filename)
        lim_filename = self._oft_env.path2c(limiter_file)
        if len(run_info) > 40:
            raise ValueError('"run_info" cannot be longer than 40 characters')
        crun_info = self._oft_env.path2c(run_info)
        if rbounds is None:
            rbounds = numpy.r_[self._tMaker.lim_contour[:,0].min(), self._tMaker.lim_contour[:,0].max()]
            dr = rbounds[1]-rbounds[0]
            rbounds += numpy.r_[-1.0,1.0]*dr*0.05
        if zbounds is None:
            zbounds = numpy.r_[self._tMaker.lim_contour[:,1].min(), self._tMaker.lim_contour[:,1].max()]
            dr = zbounds[1]-zbounds[0]
            zbounds += numpy.r_[-1.0,1.0]*dr*0.05
        if rcentr is None:
            rcentr = -1.0
        if cocos not in [2, 7]:
            raise Exception('Unsupported COCOS version. Only supported versions are 2 or 7.')
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_save_eqdsk(self._equil_ptr,cfilename,c_int(nr),c_int(nz),rbounds,zbounds,crun_info,c_double(lcfs_pad),c_double(rcentr),c_bool(truncate_eq),lim_filename,lcfs_pressure,cocos,error_string)
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
        tokamaker_save_ifile(self._equil_ptr,cfilename,npsi,ntheta,lcfs_pad,lcfs_pressure,pack_lcfs,single_precision,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
    
    def save_mug(self,filename):
        r'''! Save current equilibrium to MUG transfer format

        @param filename Filename to save equilibrium to
        '''
        cfilename = self._oft_env.path2c(filename)
        error_string = self._oft_env.get_c_errorbuff()
        tokamaker_save_mug(self._equil_ptr,cfilename,error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)