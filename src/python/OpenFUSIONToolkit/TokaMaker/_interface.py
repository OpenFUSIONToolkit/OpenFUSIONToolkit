'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
from ..util import *


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
    [c_char_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(int32,1), ctypes_numpy_array(int32,1), ctypes_numpy_array(float64,2), c_int])

# tokamaker_eval_green = ctypes_subroutine(oftpy_lib.tokamaker_eval_green,
#     [c_int, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_double, c_double, ctypes_numpy_array(float64,1)])

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

# tokamaker_solve(error_flag)
tokamaker_solve = ctypes_subroutine(oftpy_lib.tokamaker_solve, 
    [c_int_ptr])

# tokamaker_vac_solve(psi_in,error_flag)
tokamaker_vac_solve = ctypes_subroutine(oftpy_lib.tokamaker_vac_solve, 
    [ctypes_numpy_array(float64,1),  c_int_ptr])


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
tokamaker_get_dels_curr = ctypes_subroutine(oftpy_lib.tokamaker_get_dels_curr,
    [ctypes_numpy_array(numpy.float64,1)])

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
tokamaker_get_globals = ctypes_subroutine(oftpy_lib.tokamaker_get_globals, # (Itor,centroid,vol,pvol,dflux,tflux,bp_vol)
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

# tokamaker_set_flux(locations,targets,weights,ntargets,grad_wt_lim)
tokamaker_set_flux = ctypes_subroutine(oftpy_lib.tokamaker_set_flux,
    [ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_int, c_double])

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

tokamaker_save_eqdsk = ctypes_subroutine(oftpy_lib.tokamaker_save_eqdsk, #(filename,nr,nz,rbounds,zbounds,run_info,psi_pad,rcentr,error_str)
    [c_char_p, c_int, c_int, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_char_p, c_double, c_double, c_char_p])
## @endcond


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