'''! Python interface for ThinCurr thin-wall eddy current functionality

@authors Chris Hansen
@date March 2024
@ingroup doxy_oft_python
'''
from ..util import *


## @cond
# ThinCurr setup function (load mesh and setup model) (mesh_file,np,r_loc,nc,lc_loc,reg_loc,pmap_loc,tw_ptr,size,error_str)
thincurr_setup = ctypes_subroutine(oftpy_lib.thincurr_setup,
    [c_char_p, c_int, ctypes_numpy_array(float64,2), c_int, ctypes_numpy_array(int32,2), ctypes_numpy_array(int32,1), 
     ctypes_numpy_array(int32,1), c_void_p, ctypes_numpy_array(int32,1), c_char_p])

# ThinCurr setup plotting (tw_ptr,basepath,save_debug,error_str
thincurr_setup_io = ctypes_subroutine(oftpy_lib.thincurr_setup_io,
    [c_void_p, c_char_p, c_bool, c_char_p])

# thincurr_recon_curr(tw_ptr,vals,curr,format)
thincurr_recon_curr = ctypes_subroutine(oftpy_lib.thincurr_recon_curr,
    [c_void_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,2), c_int])

# thincurr_recon_field(tw_ptr,pot,coils,field,hodlr_ptr)
thincurr_recon_field = ctypes_subroutine(oftpy_lib.thincurr_recon_field,
    [c_void_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,2), c_void_p])

# ThinCurr save current potential field (tw_ptr,vals,fieldname)
thincurr_save_field = ctypes_subroutine(oftpy_lib.thincurr_save_field,
    [c_void_p, ctypes_numpy_array(float64,1), c_char_p])

# ThinCurr save scalar field (tw_ptr,vals,fieldname)
thincurr_save_scalar = ctypes_subroutine(oftpy_lib.thincurr_save_scalar,
    [c_void_p, ctypes_numpy_array(float64,1), c_char_p])

# thincurr_scale_va(tw_ptr,vals,div_flag)
thincurr_scale_va = ctypes_subroutine(oftpy_lib.thincurr_scale_va,
    [c_void_p, ctypes_numpy_array(float64,1), c_bool])

# Compute mutual coupling between models thincurr_cross_coupling(tw_ptr1,tw_ptr2,Mmat,error_str)
thincurr_cross_coupling = ctypes_subroutine(oftpy_lib.thincurr_cross_coupling,
    [c_void_p, c_void_p, ctypes_numpy_array(float64,2), c_char_p, c_char_p])

# thincurr_cross_eval(tw_ptr1,tw_ptr2,nrhs,vec1,vec2,error_str)
thincurr_cross_eval = ctypes_subroutine(oftpy_lib.thincurr_cross_eval,
    [c_void_p, c_void_p, c_int, ctypes_numpy_array(float64,2), ctypes_numpy_array(float64,2), c_char_p])

# Compute model inductance matrix thincurr_Lmat(tw_ptr,Lmat,error_str)
thincurr_Lmat = ctypes_subroutine(oftpy_lib.thincurr_Lmat,
    [c_void_p, c_bool, c_void_ptr_ptr, c_char_p, c_char_p])

# thincurr_Bmat(tw_ptr,hodlr_ptr,Bmat_ptr,Bdr_ptr,cache_file,error_str)
thincurr_Bmat = ctypes_subroutine(oftpy_lib.thincurr_Bmat,
    [c_void_p, c_void_p, c_void_ptr_ptr, c_void_ptr_ptr, c_char_p, c_char_p])

# thincurr_Mcoil(tw_ptr,Mc_ptr,cache_file,error_str)
thincurr_Mcoil = ctypes_subroutine(oftpy_lib.thincurr_Mcoil,
    [c_void_p, c_void_ptr_ptr, c_char_p, c_char_p])

# thincurr_Msensor(tw_ptr,sensor_file,Ms_ptr,Msc_ptr,nsensors,cache_file,error_str)
thincurr_Msensor = ctypes_subroutine(oftpy_lib.thincurr_Msensor,
    [c_void_p, c_char_p, c_void_ptr_ptr, c_void_ptr_ptr, c_int_ptr, c_void_p, c_char_p, c_char_p])

# Compute model resistivity matrix thincurr_curr_Rmat(tw_ptr,copy_out,Rmat,error_str)
thincurr_curr_Rmat = ctypes_subroutine(oftpy_lib.thincurr_Rmat,
    [c_void_p, c_bool, ctypes_numpy_array(float64,2), c_char_p])

# Compute current regularization matrix thincurr_curr_regmat(tw_ptr,Rmat,error_str)
thincurr_curr_regmat = ctypes_subroutine(oftpy_lib.thincurr_curr_regmat,
    [c_void_p, ctypes_numpy_array(float64,2), c_char_p])

# thincurr_eigenvalues(tw_ptr,direct,neigs,eig_vals,eig_vec,error_str)
thincurr_eigenvalues = ctypes_subroutine(oftpy_lib.thincurr_eigenvalues,
    [c_void_p, c_bool, c_int, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,2), c_void_p, c_char_p])

# thincurr_freq_response(tw_ptr,direct,fr_limit,freq,fr_driver,hodlr_ptr,error_str)
thincurr_freq_response = ctypes_subroutine(oftpy_lib.thincurr_freq_response,
    [c_void_p, c_bool, c_int, c_double, ctypes_numpy_array(float64,2), c_void_p, c_char_p])

# thincurr_time_domain(tw_ptr,direct,dt,nsteps,cg_tol,timestep_cn,nstatus,nplot,vec_ic,sensor_ptr,ncurr,curr_ptr,nvolt,volt_ptr,hodlr_ptr,error_str)
thincurr_time_domain = ctypes_subroutine(oftpy_lib.thincurr_time_domain,
    [c_void_p, c_bool, c_double, c_int, c_double, c_bool, c_int, c_int, ctypes_numpy_array(float64,1), c_void_p, c_int, ctypes_numpy_array(float64,2), c_int,
     ctypes_numpy_array(float64,2), c_void_p, c_char_p])

# thincurr_time_domain_plot(tw_ptr,compute_B,rebuild_sensors,nsteps,nplot,sensor_ptr,hodlr_ptr,error_str)
thincurr_time_domain_plot = ctypes_subroutine(oftpy_lib.thincurr_time_domain_plot,
    [c_void_p, c_bool, c_bool, c_int, c_int, c_void_p, c_void_p, c_char_p])

# thincurr_reduce_model(tw_ptr,filename,neigs,eig_vec,compute_B,sensor_ptr,hodlr_ptr,error_str)
thincurr_reduce_model = ctypes_subroutine(oftpy_lib.thincurr_reduce_model,
    [c_void_p, c_char_p, c_int, ctypes_numpy_array(float64,2), c_bool, c_void_p, c_void_p, c_char_p])
## @endcond