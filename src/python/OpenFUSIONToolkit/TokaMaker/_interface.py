#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Fortran interface definitions for TokaMaker

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
from .._interface import *


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
                ("dipole_mode", c_bool),
                ("maxits", c_int),
                ("mode", c_int),
                ("urf", c_double),
                ("nl_tol", c_double),
                ("rmin", c_double),
                ("lim_zmax", c_double),
                ("limiter_file", ctypes.c_char_p)]


## @cond
# tokamaker_alloc(tMaker_ptr,mesh_ptr,error_str)
tokamaker_alloc = ctypes_subroutine(oftpy_lib.tokamaker_alloc,
    [c_void_ptr_ptr, c_void_p, c_char_p])

# tokamaker_setup_regions(tMaker_ptr,coil_file,reg_eta,contig_flag,xpoint_mask,coil_nturns,ncoils,error_str)
tokamaker_setup_regions = ctypes_subroutine(oftpy_lib.tokamaker_setup_regions,
    [c_void_p, c_char_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(int32,1), ctypes_numpy_array(int32,1), ctypes_numpy_array(float64,2), c_int, c_char_p])

# tokamaker_setup(tMaker_ptr,order,full_domain,ncoils,error_str)
tokamaker_setup = ctypes_subroutine(oftpy_lib.tokamaker_setup,
    [c_void_p, c_int, c_bool, c_int_ptr, c_char_p])

# tokamaker_destroy(tMaker_ptr,error_str)
tokamaker_destroy = ctypes_subroutine(oftpy_lib.tokamaker_destroy,
    [c_void_p, c_char_p])

# tokamaker_load_profiles(tMaker_ptr,f_file,f_offset,p_file,eta_file,f_NI_file,error_str)
tokamaker_load_profiles = ctypes_subroutine(oftpy_lib.tokamaker_load_profiles,
    [c_void_p, c_char_p, c_double, c_char_p, c_char_p, c_char_p, c_char_p])

# tokamaker_init_psi(tMaker_ptr,r0,z0,a,kappa,delta,rhs_source,error_str)
tokamaker_init_psi = ctypes_subroutine(oftpy_lib.tokamaker_init_psi,
    [c_void_p, c_double, c_double, c_double, c_double, c_double, c_double_ptr, c_char_p])

# tokamaker_solve(tMaker_ptr,vacuum,error_str)
tokamaker_solve = ctypes_subroutine(oftpy_lib.tokamaker_solve, 
    [c_void_p, c_bool, c_char_p])

# tokamaker_vac_solve(tMaker_ptr,psi_in,rhs_source,error_str)
tokamaker_vac_solve = ctypes_subroutine(oftpy_lib.tokamaker_vac_solve, 
    [c_void_p, ctypes_numpy_array(float64,1), c_double_ptr,  c_char_p])

# tokamaker_setup_td(tMaker_ptr,dt,lin_tol,nl_tol,pre_plasma,error_str)
tokamaker_setup_td = ctypes_subroutine(oftpy_lib.tokamaker_setup_td,
    [c_void_p, c_double, c_double, c_double, c_bool, c_char_p])

# tokamaker_eig_td(tMaker_ptr,omega,neigs,eigs,eig_vecs,include_bounds,eta_plasma,pm,error_str)
tokamaker_eig_td = ctypes_subroutine(oftpy_lib.tokamaker_eig_td,
    [c_void_p, c_double, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool, c_double, c_bool, c_char_p])

# tokamaker_eig_wall(tMaker_ptr,neigs,eigs,eig_vecs,pm,error_str)
tokamaker_eig_wall = ctypes_subroutine(oftpy_lib.tokamaker_eig_wall,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_bool, c_char_p])

# tokamaker_step_td(tMaker_ptr,time,dt,nl_its,lin_its,nretry,error_str)
tokamaker_step_td = ctypes_subroutine(oftpy_lib.tokamaker_step_td,
    [c_void_p, c_double_ptr, c_double_ptr, c_int_ptr, c_int_ptr, c_int_ptr, c_char_p])

# tokamaker_get_mesh(tMaker_ptr,np,r_loc,nc,lc_loc,reg_loc,error_str)
tokamaker_get_mesh = ctypes_subroutine(oftpy_lib.tokamaker_get_mesh,
    [c_void_p, c_int_ptr, c_double_ptr_ptr, c_int_ptr, c_int_ptr_ptr, c_int_ptr_ptr, c_char_p])

# tokamaker_get_limiter(tMaker_ptr,np,r_loc,nloops,loop_ptr,error_str)
tokamaker_get_limiter = ctypes_subroutine(oftpy_lib.tokamaker_get_limiter,
    [c_void_p, c_int_ptr, c_double_ptr_ptr, c_int_ptr, c_int_ptr_ptr, c_char_p])

# tokamaker_get_psi(tMaker_ptr,psi_vals,psi_lim,psi_max,error_str)
tokamaker_get_psi = ctypes_subroutine(oftpy_lib.tokamaker_get_psi,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_double_ptr, c_double_ptr, c_char_p])

# tokamaker_get_dels_curr(tMaker_ptr,psi_vals,error_str)
tokamaker_get_dels_curr = ctypes_subroutine(oftpy_lib.tokamaker_get_dels_curr,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_get_jtor(tMaker_ptr,jtor,error_str)
tokamaker_get_jtor = ctypes_subroutine(oftpy_lib.tokamaker_get_jtor,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_area_int(tMaker_ptr,vec_vals,reg_ind,result,error_str)
tokamaker_area_int = ctypes_subroutine(oftpy_lib.tokamaker_area_int,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_int, c_double_ptr, c_char_p])

# tokamaker_flux_int(tMaker_ptr,psi_vals,field_vals,nvals,result,error_str)
tokamaker_flux_int = ctypes_subroutine(oftpy_lib.tokamaker_flux_int,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_int, c_double_ptr, c_char_p])

# tokamaker_get_coil_currents(tMaker_ptr,currents,reg_currents,error_str)
tokamaker_get_coil_currents = ctypes_subroutine(oftpy_lib.tokamaker_get_coil_currents,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_get_coil_Lmat(tMaker_ptr,Lmat,error_str)
tokamaker_get_coil_Lmat = ctypes_subroutine(oftpy_lib.tokamaker_get_coil_Lmat,
    [c_void_p, ctypes_numpy_array(numpy.float64,2), c_char_p])

# tokamaker_get_refs(tMaker_ptr,o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm,error_str)
tokamaker_get_refs = ctypes_subroutine(oftpy_lib.tokamaker_get_refs,
    [c_void_p, c_double_ptr_ptr, c_double_ptr_ptr, c_double_ptr_ptr, c_bool_ptr_ptr, c_double_ptr_ptr, c_double_ptr_ptr,  c_double_ptr_ptr, c_char_p])

# tokamaker_trace_surf(tMaker_ptr,psi_surf,points,npoints,error_str)
tokamaker_trace_surf = ctypes_subroutine(oftpy_lib.tokamaker_trace_surf,
    [c_void_p, c_double, c_double_ptr_ptr, c_int_ptr, c_char_p])

# tokamaker_get_q(tMaker_ptr,npsi,psi_q,qvals,ravgs,dl,rbounds,zbounds,error_str)
tokamaker_get_q = ctypes_subroutine(oftpy_lib.tokamaker_get_q,
    [c_void_p, c_int,ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,2),
     c_double_ptr, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), c_char_p])

# tokamaker_sauter_fc(tMaker_ptr,npsi,psi_saut,fc,r_avgs,modb_avgs,error_str)
tokamaker_sauter_fc = ctypes_subroutine(oftpy_lib.tokamaker_sauter_fc,
    [c_void_p, c_int,ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,2),
     ctypes_numpy_array(numpy.float64,2), c_char_p])

# tokamaker_get_globals(tMaker_ptr,Itor,centroid,vol,pvol,dflux,tflux,bp_vol,error_str)
tokamaker_get_globals = ctypes_subroutine(oftpy_lib.tokamaker_get_globals,
    [c_void_p, c_double_ptr, ctypes_numpy_array(numpy.float64,1), c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr, c_char_p])

# tokamaker_gs_calc_vloop(tMaker_ptr,vloop,error_str)
tokamaker_gs_calc_vloop = ctypes_subroutine(oftpy_lib.tokamaker_gs_calc_vloop,
    [c_void_p, c_double_ptr, c_char_p])

# tokamaker_get_profs(tMaker_ptr,npsi,psi_in,f,fp,p,pp,error_str)
tokamaker_get_profs = ctypes_subroutine(oftpy_lib.tokamaker_get_profs,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), 
     ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_get_vfixed(tMaker_ptr,npts,pts,fluxes,error_str)
tokamaker_get_vfixed = ctypes_subroutine(oftpy_lib.tokamaker_get_vfixed,
    [c_void_p, c_int_ptr, c_double_ptr_ptr, c_double_ptr_ptr, c_char_p])

# tokamaker_get_field_eval(tMaker_ptr,imode,int_obj,error_str)
tokamaker_get_field_eval = ctypes_subroutine(oftpy_lib.tokamaker_get_field_eval,
    [c_void_p, c_int, c_void_ptr_ptr, c_char_p])

# tokamaker_apply_field_eval(tMaker_ptr,int_obj,int_type,pt,fbary_tol,cell,dim,field)
tokamaker_apply_field_eval = ctypes_subroutine(oftpy_lib.tokamaker_apply_field_eval,
    [c_void_p, c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_double, c_int_ptr, c_int, ctypes_numpy_array(numpy.float64,1)])

# tokamaker_set_psi(tMaker_ptr,psi_vals,update_bounds,error_str)
tokamaker_set_psi = ctypes_subroutine(oftpy_lib.tokamaker_set_psi,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_bool, c_char_p])

# tokamaker_set_psi_dt(tMaker_ptr,psi_vals,dt,error_str)
tokamaker_set_psi_dt = ctypes_subroutine(oftpy_lib.tokamaker_set_psi_dt,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_double, c_char_p])

# tokamaker_set_settings(tMaker_ptr,settings,error_str)
tokamaker_set_settings = ctypes_subroutine(oftpy_lib.tokamaker_set_settings,
    [c_void_p, ctypes.POINTER(tokamaker_settings_struct), c_char_p])

# tokamaker_set_dipole_a(tMaker_ptr,dipole_a,error_str)
tokamaker_set_dipole_a = ctypes_subroutine(oftpy_lib.tokamaker_set_dipole_a,
    [c_void_p, c_double, c_char_p])

# tokamaker_set_targets(tMaker_ptr,ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target,error_str)
tokamaker_set_targets = ctypes_subroutine(oftpy_lib.tokamaker_set_targets,
    [c_void_p, c_double, c_double, c_double, c_double, c_double, c_double, c_char_p])

# tokamaker_set_isoflux(tMaker_ptr,targets,ref_points,weights,ntargets,grad_wt_lim,error_str)
tokamaker_set_isoflux = ctypes_subroutine(oftpy_lib.tokamaker_set_isoflux,
    [c_void_p, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), c_int, c_double, c_char_p])

# tokamaker_set_flux(tMaker_ptr,locations,targets,weights,ntargets,grad_wt_lim,error_str)
tokamaker_set_flux = ctypes_subroutine(oftpy_lib.tokamaker_set_flux,
    [c_void_p, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_int, c_double, c_char_p])

# tokamaker_set_saddles(tMaker_ptr,targets,weights,ntargets,error_str)
tokamaker_set_saddles = ctypes_subroutine(oftpy_lib.tokamaker_set_saddles,
    [c_void_p, ctypes_numpy_array(numpy.float64,2), ctypes_numpy_array(numpy.float64,1), c_int, c_char_p])

# tokamaker_set_coil_currents(tMaker_ptr,currents,error_str)
tokamaker_set_coil_currents = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_currents,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_set_coil_regmat(tMaker_ptr,nregularize,coil_reg_mat,coil_reg_targets,coil_reg_weights,error_str)
tokamaker_set_coil_regmat = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_regmat,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,2),ctypes_numpy_array(numpy.float64,1),ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_set_coil_bounds(tMaker_ptr,coil_bounds,error_str)
tokamaker_set_coil_bounds = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_bounds,
    [c_void_p, ctypes_numpy_array(numpy.float64,2), c_char_p])

# tokamaker_set_coil_vsc(tMaker_ptr,coil_gains,error_str)
tokamaker_set_coil_vsc = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_vsc,
    [c_void_p, ctypes_numpy_array(numpy.float64,1), c_char_p])

# tokamaker_save_eqdsk(tMaker_ptr,filename,nr,nz,rbounds,zbounds,run_info,psi_pad,rcentr,trunc_eq,lim_filename,lcfs_press,error_str)
tokamaker_save_eqdsk = ctypes_subroutine(oftpy_lib.tokamaker_save_eqdsk,
    [c_void_p, c_char_p, c_int, c_int, ctypes_numpy_array(numpy.float64,1), ctypes_numpy_array(numpy.float64,1), c_char_p,
     c_double, c_double, c_bool, c_char_p, c_double, c_char_p])

# tokamaker_save_ifile(tMaker_ptr,filename,npsi,ntheta,psi_pad,lcfs_press,pack_lcfs,single_prec,error_str)
tokamaker_save_ifile = ctypes_subroutine(oftpy_lib.tokamaker_save_ifile,
    [c_void_p, c_char_p, c_int, c_int, c_double, c_double, c_bool, c_bool, c_char_p])

# tokamaker_save_mug(tMaker_ptr,filename,error_str)
tokamaker_save_mug = ctypes_subroutine(oftpy_lib.tokamaker_save_mug,
    [c_void_p, c_char_p, c_char_p])

# tokamaker_set_coil_current_dist(tMaker_ptr,iCoil,curr_dist,error_str)
tokamaker_set_coil_current_dist = ctypes_subroutine(oftpy_lib.tokamaker_set_coil_current_dist,
    [c_void_p, c_int, ctypes_numpy_array(numpy.float64,1), c_char_p])
## @endcond


class TokaMaker_field_interpolator():
    '''! Interpolation class for Grad-Shafranov fields'''
    def __init__(self,tMaker_obj,int_obj,int_type,dim,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param tMaker_obj Address of FORTRAN TokaMaker class
        @param int_obj Address of FORTRAN interpolation class
        @param int_type Interpolation type (see @ref TokaMaker.TokaMaker.get_field_eval "get_field_eval")
        @param dim Dimension of vector field
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        self.cell = c_int(-1)
        self.int_type = int_type
        self.dim_return = dim
        if dim == 2:
            self.dim_eval = 3
        else:
            self.dim_eval = dim
        self.val = numpy.zeros((self.dim_eval,), dtype=numpy.float64)
        self._tMaker_obj = tMaker_obj
        self._int_obj = int_obj
        self.fbary_tol = fbary_tol
    
    def __del__(self):
        '''Destroy underlying interpolation object'''
        pt_eval = numpy.zeros((3,), dtype=numpy.float64)
        tokamaker_apply_field_eval(self._tMaker_obj,self._int_obj,-self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.dim_eval,self.val)

    def eval(self,pt):
        '''! Evaluate field at a given location

        @param pt Location for evaluation [2]
        @result Field at evaluation point [self.dim_return]
        '''
        pt_eval = numpy.zeros((3,), dtype=numpy.float64)
        pt_eval[:2] = pt
        tokamaker_apply_field_eval(self._tMaker_obj,self._int_obj,self.int_type,pt_eval,self.fbary_tol,ctypes.byref(self.cell),self.dim_eval,self.val)
        return self.val[:self.dim_return]

