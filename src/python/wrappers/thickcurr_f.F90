!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file thickcurr_f.F90
!
!> Fortran part of Python wrapper for ThickCurr thick-wall electromagnetic functionality
!!
!! @note As the underlying @ref thick_wall module stores its state at the module
!! level only a single ThickCurr object may be active at a time.
!---------------------------------------------------------------------------------
MODULE thickcurr_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr, c_associated
!---Base
USE oft_base
USE oft_io, ONLY: oft_bin_file
!--Grid
USE multigrid, ONLY: multigrid_mesh, multigrid_reset
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---FE spaces
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup, oft_hcurl_fem
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, &
  oft_hcurl_grad_project, oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob, &
  oft_hcurl_grad_rinterp, oft_hcurl_grad_cinterp
!---Physics
USE diagnostic, ONLY: field_probe
USE thick_wall, ONLY: thickcurr_ML_hcurl, thickcurr_ML_H1grad, thickcurr_ML_hcurl_grad, &
  thickcurr_ML_lagrange, thickcurr_ML_vlagrange, thickcurr_minlev, thickcurr_xml, &
  thickcurr_setup, thickcurr_delete, thickcurr_set_options, thickcurr_eigs, thickcurr_run, &
  thickcurr_plot, thickcurr_diag, oft_thickcurr_probe
!---Wrappers
USE oft_base_f, ONLY: copy_string, copy_string_rev
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> ThickCurr object for Python interface
!---------------------------------------------------------------------------------
TYPE :: thickcurr_py
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL() !< Mesh
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL() !< Lagrange FE space (plotting)
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_vlagrange => NULL() !< Vector Lagrange FE space (plotting)
  TYPE(oft_ml_fem_type), POINTER :: ML_H1 => NULL() !< H^1 FE space
  TYPE(oft_ml_fem_type), POINTER :: ML_H1grad => NULL() !< Grad(H^1) FE space
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL() !< H(Curl) FE space
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_hcurl_grad => NULL() !< H(Curl) + Grad(H^1) FE space
  CLASS(oft_vector), POINTER :: a => NULL() !< Vector potential (initial condition, then final state of last run)
END TYPE thickcurr_py
!---------------------------------------------------------------------------------
!> Interpolation object for values tabulated at the quadrature points used for projection
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: thickcurr_tab_interp
  INTEGER(i4) :: nq = 0 !< Number of quadrature points per cell
  REAL(r8), POINTER, DIMENSION(:,:) :: qpts => NULL() !< Quadrature points in logical coordinates
  REAL(r8), POINTER, DIMENSION(:,:,:) :: vals => NULL() !< Values at each point [3,nq,nc]
CONTAINS
  PROCEDURE :: interp => tab_interp_apply
END TYPE thickcurr_tab_interp
!---------------------------------------------------------------------------------
!> Sensor (probe) object for Python interface
!!
!! Each signal is a weighted sum of the vector potential and magnetic field
!! sampled at sets of points, \f$ s_k = \sum_{i \in k} w_i \cdot A(x_i) +
!! \sum_{j \in k} w_j \cdot B(x_j) \f$, which allows flux loops (\f$ \oint A \cdot dl \f$)
!! and point (Mirnov) sensors (\f$ B \cdot \hat{n} \f$) to be represented.
!---------------------------------------------------------------------------------
TYPE, EXTENDS(oft_thickcurr_probe) :: thickcurr_sensors
  INTEGER(i4) :: nsensors = 0 !< Number of sensor signals
  INTEGER(i4) :: npa = 0 !< Number of vector potential sample points
  INTEGER(i4) :: npb = 0 !< Number of magnetic field sample points
  INTEGER(i4), POINTER, DIMENSION(:) :: idxa => NULL() !< Signal index for each A point
  INTEGER(i4), POINTER, DIMENSION(:) :: idxb => NULL() !< Signal index for each B point
  REAL(r8), POINTER, DIMENSION(:,:) :: wa => NULL() !< Weights for each A point [3,npa]
  REAL(r8), POINTER, DIMENSION(:,:) :: wb => NULL() !< Weights for each B point [3,npb]
  CHARACTER(LEN=40), POINTER, DIMENSION(:) :: names => NULL() !< Sensor names
  CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'floops.hist' !< History file name
  LOGICAL :: hist_ready = .FALSE. !< Has history file been initialized?
  TYPE(field_probe) :: aprobe !< Probe for vector potential
  TYPE(field_probe) :: bprobe !< Probe for magnetic field
  TYPE(oft_hcurl_grad_rinterp) :: ainterp !< Vector potential interpolation object
  TYPE(oft_hcurl_grad_cinterp) :: binterp !< Magnetic field interpolation object
  TYPE(oft_bin_file) :: hist_file !< History file
CONTAINS
  PROCEDURE :: apply => sensors_apply
  PROCEDURE :: flush => sensors_flush
END TYPE thickcurr_sensors
LOGICAL :: tc_active = .FALSE. !< Is a ThickCurr object active?
CONTAINS
!---------------------------------------------------------------------------------
!> Cast C pointer to ThickCurr object
!---------------------------------------------------------------------------------
FUNCTION thickcurr_ccast(tc_cptr,tc_fptr,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: tc_cptr !< C pointer to ThickCurr object
TYPE(thickcurr_py), POINTER, INTENT(out) :: tc_fptr !< Fortran pointer to ThickCurr object
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
IF(PRESENT(error_str))CALL copy_string('',error_str)
IF(.NOT.c_associated(tc_cptr))THEN
  IF(PRESENT(error_str))CALL copy_string('ThickCurr object not associated',error_str)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(tc_cptr,tc_fptr)
success=.TRUE.
END FUNCTION thickcurr_ccast
!---------------------------------------------------------------------------------
!> Setup FE spaces and ThickCurr operators
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_setup(tc_ptr,mesh_ptr,order,minlev,neta,eta_ptr,xml_ptr,error_str) BIND(C,NAME="thickcurr_f_setup")
TYPE(c_ptr), INTENT(out) :: tc_ptr !< Pointer to ThickCurr object
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to mesh object
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< FE order
INTEGER(KIND=c_int), VALUE, INTENT(in) :: minlev !< Lowest MG level (-1: no MG, 0: p-multigrid on finest mesh)
INTEGER(KIND=c_int), VALUE, INTENT(in) :: neta !< Size of resistivity list (0 to use XML)
TYPE(c_ptr), VALUE, INTENT(in) :: eta_ptr !< Resistivity for each region [neta] (<=0 for vacuum)
TYPE(c_ptr), VALUE, INTENT(in) :: xml_ptr !< Pointer to XML root node (optional)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: toplev,lev_min
REAL(r8), POINTER, DIMENSION(:) :: eta_tmp
TYPE(thickcurr_py), POINTER :: self
CALL copy_string('',error_str)
tc_ptr=c_null_ptr
IF(tc_active)THEN
  CALL copy_string('Only one ThickCurr object can be active at a time',error_str)
  RETURN
END IF
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
IF(order<1)THEN
  CALL copy_string('Invalid FE order, must be >= 1',error_str)
  RETURN
END IF
ALLOCATE(self)
CALL c_f_pointer(mesh_ptr,self%ml_mesh)
!---Determine FE levels
toplev=self%ml_mesh%mgmax+order-1
IF(minlev<0)THEN
  lev_min=toplev
ELSE IF(minlev==0)THEN
  lev_min=self%ml_mesh%mgmax ! p-multigrid (linear elements on finest mesh)
ELSE
  lev_min=minlev
  IF(oft_env%nprocs>1)lev_min=MAX(oft_env%nbase+1,lev_min)
END IF
IF(lev_min>toplev)THEN
  CALL copy_string('Invalid "minlev", must be <= highest FE level',error_str)
  DEALLOCATE(self)
  RETURN
END IF
!---Lagrange (for plotting, highest level only)
ALLOCATE(self%ML_lagrange,self%ML_vlagrange)
CALL oft_lag_setup(self%ml_mesh,order,self%ML_lagrange,ML_vlag_obj=self%ML_vlagrange,minlev=toplev)
!---Grad(H^1) subspace
ALLOCATE(self%ML_H1,self%ML_H1grad,self%ML_hcurl,self%ML_hcurl_grad)
CALL oft_h1_setup(self%ml_mesh,order+1,self%ML_H1,minlev=lev_min)
!---H(Curl) subspace
CALL oft_hcurl_setup(self%ml_mesh,order,self%ML_hcurl,minlev=lev_min)
!---Full H(Curl) + Grad(H^1) space
CALL oft_hcurl_grad_setup(self%ML_hcurl,self%ML_H1,self%ML_hcurl_grad,self%ML_H1grad,minlev=lev_min)
IF(lev_min<toplev)THEN
  CALL h1_setup_interp(self%ML_H1)
  CALL hcurl_setup_interp(self%ML_hcurl)
  CALL hcurl_grad_setup_interp(self%ML_hcurl_grad,self%ML_H1,create_full=.TRUE.)
END IF
!---Set FE spaces for ThickCurr
thickcurr_ML_hcurl=>self%ML_hcurl
thickcurr_ML_H1grad=>self%ML_H1grad
thickcurr_ML_hcurl_grad=>self%ML_hcurl_grad
thickcurr_ML_lagrange=>self%ML_lagrange
thickcurr_ML_vlagrange=>self%ML_vlagrange
thickcurr_minlev=lev_min
NULLIFY(thickcurr_xml)
IF(c_associated(xml_ptr))CALL c_f_pointer(xml_ptr,thickcurr_xml)
!---Setup ThickCurr
IF(neta>0)THEN
  CALL c_f_pointer(eta_ptr,eta_tmp,[neta])
  CALL thickcurr_setup(eta_tmp)
ELSE
  CALL thickcurr_setup()
END IF
!---Initial vector potential (zero)
CALL self%ML_hcurl_grad%vec_create(self%a)
CALL self%a%set(0.d0)
tc_active=.TRUE.
tc_ptr=C_LOC(self)
END SUBROUTINE thickcurr_f_setup
!---------------------------------------------------------------------------------
!> Destroy ThickCurr object
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_destroy(tc_ptr,error_str) BIND(C,NAME="thickcurr_f_destroy")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
CALL thickcurr_delete()
IF(ASSOCIATED(self%a))THEN
  CALL self%a%delete()
  DEALLOCATE(self%a)
END IF
IF(ASSOCIATED(self%ML_lagrange))THEN
  CALL self%ML_lagrange%delete()
  DEALLOCATE(self%ML_lagrange)
END IF
IF(ASSOCIATED(self%ML_hcurl_grad))THEN
  CALL self%ML_hcurl_grad%delete()
  DEALLOCATE(self%ML_hcurl_grad)
END IF
IF(ASSOCIATED(self%ML_hcurl))THEN
  CALL self%ML_hcurl%delete()
  DEALLOCATE(self%ML_hcurl)
END IF
IF(ASSOCIATED(self%ML_H1))THEN
  CALL self%ML_H1%delete()
  DEALLOCATE(self%ML_H1)
END IF
IF(ASSOCIATED(self%ml_mesh))THEN
  CALL multigrid_reset(self%ml_mesh)
  DEALLOCATE(self%ml_mesh)
END IF
NULLIFY(thickcurr_ML_hcurl,thickcurr_ML_H1grad,thickcurr_ML_hcurl_grad,thickcurr_ML_lagrange, &
  thickcurr_ML_vlagrange,thickcurr_xml)
DEALLOCATE(self)
tc_active=.FALSE.
END SUBROUTINE thickcurr_f_destroy
!---------------------------------------------------------------------------------
!> Update run options (negative values are ignored)
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_set_options(tc_ptr,dt,theta,lin_tol,nsteps,rst_freq,error_str) BIND(C,NAME="thickcurr_f_set_options")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
REAL(c_double), VALUE, INTENT(in) :: dt !< Time step [s]
REAL(c_double), VALUE, INTENT(in) :: theta !< Time centering parameter
REAL(c_double), VALUE, INTENT(in) :: lin_tol !< Linear solver tolerance
INTEGER(c_int), VALUE, INTENT(in) :: nsteps !< Number of time steps
INTEGER(c_int), VALUE, INTENT(in) :: rst_freq !< Restart (plot) file frequency
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
IF(dt>0.d0)CALL thickcurr_set_options(dt=dt)
IF(theta>0.d0)THEN
  IF(theta<0.5d0.OR.theta>1.d0)THEN
    CALL copy_string('Invalid "theta", must be in [0.5,1]',error_str)
    RETURN
  END IF
  CALL thickcurr_set_options(theta=theta)
END IF
IF(lin_tol>0.d0)CALL thickcurr_set_options(lin_tol=lin_tol)
IF(nsteps>=0)CALL thickcurr_set_options(nsteps=nsteps)
IF(rst_freq>0)CALL thickcurr_set_options(rst_freq=rst_freq)
END SUBROUTINE thickcurr_f_set_options
!---------------------------------------------------------------------------------
!> Get quadrature points used for projection of initial conditions
!!
!! If `pts_ptr` is NULL only the number of points is returned
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_get_qpts(tc_ptr,npts,pts_ptr,error_str) BIND(C,NAME="thickcurr_f_get_qpts")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
INTEGER(c_int), INTENT(out) :: npts !< Number of points on this process
TYPE(c_ptr), VALUE, INTENT(in) :: pts_ptr !< Point locations [3,npts] (output)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: i,m
REAL(r8), POINTER, DIMENSION(:,:) :: pts
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
IF(.NOT.get_hcurl_rep(self,hcurl_rep,error_str))RETURN
npts=hcurl_rep%quad%np*hcurl_rep%mesh%nc
IF(.NOT.c_associated(pts_ptr))RETURN
CALL c_f_pointer(pts_ptr,pts,[3,npts])
!$omp parallel do private(m)
DO i=1,hcurl_rep%mesh%nc
  DO m=1,hcurl_rep%quad%np
    pts(:,(i-1)*hcurl_rep%quad%np+m)=hcurl_rep%mesh%log2phys(i,hcurl_rep%quad%pts(:,m))
  END DO
END DO
END SUBROUTINE thickcurr_f_get_qpts
!---------------------------------------------------------------------------------
!> Set vector potential from values at quadrature points (see @ref thickcurr_f_get_qpts)
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_set_field(tc_ptr,vals_ptr,zero_bdry,energy,error_str) BIND(C,NAME="thickcurr_f_set_field")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
TYPE(c_ptr), VALUE, INTENT(in) :: vals_ptr !< Vector potential at quadrature points [3,npts] (NULL for zero field)
LOGICAL(c_bool), VALUE, INTENT(in) :: zero_bdry !< Zero tangential vector potential on boundary?
REAL(c_double), INTENT(out) :: energy !< Magnetic energy of projected field [J]
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
CLASS(oft_matrix), POINTER :: mop
CLASS(oft_solver), POINTER :: minv
CLASS(oft_vector), POINTER :: v
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
TYPE(oft_hcurl_grad_gzerop) :: gzerop
TYPE(oft_hcurl_grad_zerob) :: zerob
TYPE(thickcurr_tab_interp) :: tab_field
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
IF(.NOT.get_hcurl_rep(self,hcurl_rep,error_str))RETURN
energy=0.d0
CALL self%a%set(0.d0)
IF(.NOT.c_associated(vals_ptr))RETURN
!---Setup tabulated field
tab_field%nq=hcurl_rep%quad%np
tab_field%qpts=>hcurl_rep%quad%pts
tab_field%mesh=>hcurl_rep%mesh
CALL c_f_pointer(vals_ptr,tab_field%vals,[3,tab_field%nq,hcurl_rep%mesh%nc])
!---Mass matrix and solver
NULLIFY(mop)
IF(zero_bdry)THEN
  CALL hcurl_grad_getmop(self%ML_hcurl_grad%current_level,mop,"zerob")
ELSE
  CALL hcurl_grad_getmop(self%ML_hcurl_grad%current_level,mop,"none")
END IF
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-2
CALL create_diag_pre(minv%pre)
!---Project onto H(Curl) + Grad(H^1) basis
CALL self%ML_hcurl_grad%vec_create(v)
CALL oft_hcurl_grad_project(self%ML_hcurl_grad%current_level,tab_field,v)
gzerop%ML_hcurl_grad_rep=>self%ML_hcurl_grad
zerob%ML_hcurl_grad_rep=>self%ML_hcurl_grad
IF(zero_bdry)CALL zerob%apply(v)
CALL gzerop%apply(v)
CALL minv%apply(self%a,v)
CALL thickcurr_diag(self%a,energy)
!---Cleanup
NULLIFY(tab_field%vals,tab_field%qpts)
CALL v%delete()
CALL mop%delete()
CALL minv%pre%delete()
CALL minv%delete()
DEALLOCATE(minv%pre)
DEALLOCATE(v,mop,minv)
END SUBROUTINE thickcurr_f_set_field
!---------------------------------------------------------------------------------
!> Get magnetic energy of current vector potential
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_get_energy(tc_ptr,energy,error_str) BIND(C,NAME="thickcurr_f_get_energy")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
REAL(c_double), INTENT(out) :: energy !< Magnetic energy [J]
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
CALL thickcurr_diag(self%a,energy)
END SUBROUTINE thickcurr_f_get_energy
!---------------------------------------------------------------------------------
!> Compute slowest decaying eigenmodes
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_eigs(tc_ptr,nev,tau_ptr,shift,save_plot,error_str) BIND(C,NAME="thickcurr_f_eigs")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
INTEGER(c_int), VALUE, INTENT(in) :: nev !< Number of eigenvalues to compute
TYPE(c_ptr), VALUE, INTENT(in) :: tau_ptr !< Decay times [nev] (output)
REAL(c_double), VALUE, INTENT(in) :: shift !< Shift for eigenvalue solve (<=0 to use time step)
LOGICAL(c_bool), VALUE, INTENT(in) :: save_plot !< Save eigenmodes for plotting?
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(r8), POINTER, DIMENSION(:) :: tau
LOGICAL :: plot_flag
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
#ifndef HAVE_ARPACK
CALL copy_string('Eigenvalue solves require ARPACK',error_str)
RETURN
#endif
CALL c_f_pointer(tau_ptr,tau,[nev])
plot_flag=save_plot
IF(shift>0.d0)THEN
  CALL thickcurr_eigs(nev,tau,save_plot=plot_flag,shift=shift)
ELSE
  CALL thickcurr_eigs(nev,tau,save_plot=plot_flag)
END IF
END SUBROUTINE thickcurr_f_eigs
!---------------------------------------------------------------------------------
!> Run time-dependent simulation from current vector potential
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_run(tc_ptr,sensor_ptr,t_final,error_str) BIND(C,NAME="thickcurr_f_run")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Pointer to sensor object (optional)
REAL(c_double), INTENT(out) :: t_final !< Time at end of simulation [s]
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(thickcurr_sensors), POINTER :: sensors
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
oft_env%pm=.TRUE.
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr,sensors)
  sensors%hist_ready=.FALSE.
  CALL thickcurr_run(self%a,probes=sensors,t_final=t_final)
ELSE
  CALL thickcurr_run(self%a,t_final=t_final)
END IF
END SUBROUTINE thickcurr_f_run
!---------------------------------------------------------------------------------
!> Generate plot files from restart files
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_plot(tc_ptr,sensor_ptr,rst_start,rst_end,error_str) BIND(C,NAME="thickcurr_f_plot")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Pointer to sensor object (optional)
INTEGER(c_int), VALUE, INTENT(in) :: rst_start !< First restart file to plot
INTEGER(c_int), VALUE, INTENT(in) :: rst_end !< Last restart file to plot
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(thickcurr_sensors), POINTER :: sensors
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr,sensors)
  sensors%hist_ready=.FALSE.
  CALL thickcurr_plot(probes=sensors,rst_range=[rst_start,rst_end])
ELSE
  CALL thickcurr_plot(rst_range=[rst_start,rst_end])
END IF
END SUBROUTINE thickcurr_f_plot
!---------------------------------------------------------------------------------
!> Setup sensor object
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_setup_sensors(tc_ptr,nsensors,names_ptr,npa,ptsa_ptr,wa_ptr,idxa_ptr, &
  npb,ptsb_ptr,wb_ptr,idxb_ptr,filename,sensor_ptr,error_str) BIND(C,NAME="thickcurr_f_setup_sensors")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
INTEGER(c_int), VALUE, INTENT(in) :: nsensors !< Number of sensors
CHARACTER(KIND=c_char), INTENT(in) :: names_ptr(40*nsensors) !< Sensor names (fixed width of 40)
INTEGER(c_int), VALUE, INTENT(in) :: npa !< Number of vector potential sample points
TYPE(c_ptr), VALUE, INTENT(in) :: ptsa_ptr !< Vector potential sample points [3,npa]
TYPE(c_ptr), VALUE, INTENT(in) :: wa_ptr !< Weights for vector potential samples [3,npa]
TYPE(c_ptr), VALUE, INTENT(in) :: idxa_ptr !< Sensor index (1-based) for vector potential samples [npa]
INTEGER(c_int), VALUE, INTENT(in) :: npb !< Number of magnetic field sample points
TYPE(c_ptr), VALUE, INTENT(in) :: ptsb_ptr !< Magnetic field sample points [3,npb]
TYPE(c_ptr), VALUE, INTENT(in) :: wb_ptr !< Weights for magnetic field samples [3,npb]
TYPE(c_ptr), VALUE, INTENT(in) :: idxb_ptr !< Sensor index (1-based) for magnetic field samples [npb]
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< History file name
TYPE(c_ptr), INTENT(out) :: sensor_ptr !< Pointer to sensor object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: i,j,np
INTEGER(i4), POINTER, DIMENSION(:) :: itmp
REAL(r8), POINTER, DIMENSION(:,:) :: rtmp
TYPE(thickcurr_sensors), POINTER :: sensors
TYPE(thickcurr_py), POINTER :: self
sensor_ptr=c_null_ptr
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
ALLOCATE(sensors)
sensors%nsensors=nsensors
ALLOCATE(sensors%names(nsensors))
DO i=1,nsensors
  sensors%names(i)=''
  DO j=1,40
    IF(names_ptr((i-1)*40+j)==c_null_char)EXIT
    sensors%names(i)(j:j)=names_ptr((i-1)*40+j)
  END DO
END DO
CALL copy_string_rev(filename,sensors%filename)
!---Vector potential samples
sensors%npa=npa
IF(npa>0)THEN
  ALLOCATE(sensors%wa(3,npa),sensors%idxa(npa))
  CALL c_f_pointer(wa_ptr,rtmp,[3,npa])
  sensors%wa=rtmp
  CALL c_f_pointer(idxa_ptr,itmp,[npa])
  sensors%idxa=itmp
  CALL c_f_pointer(ptsa_ptr,rtmp,[3,npa])
  sensors%aprobe%mesh=>self%ml_mesh%mesh
  np=npa
  CALL sensors%aprobe%setup(np,pts=rtmp)
  sensors%aprobe%B=>sensors%ainterp
END IF
!---Magnetic field samples
sensors%npb=npb
IF(npb>0)THEN
  ALLOCATE(sensors%wb(3,npb),sensors%idxb(npb))
  CALL c_f_pointer(wb_ptr,rtmp,[3,npb])
  sensors%wb=rtmp
  CALL c_f_pointer(idxb_ptr,itmp,[npb])
  sensors%idxb=itmp
  CALL c_f_pointer(ptsb_ptr,rtmp,[3,npb])
  sensors%bprobe%mesh=>self%ml_mesh%mesh
  np=npb
  CALL sensors%bprobe%setup(np,pts=rtmp)
  sensors%bprobe%B=>sensors%binterp
END IF
sensor_ptr=C_LOC(sensors)
END SUBROUTINE thickcurr_f_setup_sensors
!---------------------------------------------------------------------------------
!> Evaluate sensor signals for current vector potential
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_eval_sensors(tc_ptr,sensor_ptr,vals_ptr,error_str) BIND(C,NAME="thickcurr_f_eval_sensors")
TYPE(c_ptr), VALUE, INTENT(in) :: tc_ptr !< Pointer to ThickCurr object
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Pointer to sensor object
TYPE(c_ptr), VALUE, INTENT(in) :: vals_ptr !< Sensor signals [nsensors] (output)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(r8), POINTER, DIMENSION(:) :: vals
TYPE(thickcurr_sensors), POINTER :: sensors
TYPE(thickcurr_py), POINTER :: self
IF(.NOT.thickcurr_ccast(tc_ptr,self,error_str))RETURN
IF(.NOT.c_associated(sensor_ptr))THEN
  CALL copy_string('Sensor object not associated',error_str)
  RETURN
END IF
CALL c_f_pointer(sensor_ptr,sensors)
CALL c_f_pointer(vals_ptr,vals,[sensors%nsensors])
CALL sensors_eval(sensors,self%a,vals)
END SUBROUTINE thickcurr_f_eval_sensors
!---------------------------------------------------------------------------------
!> Destroy sensor object
!---------------------------------------------------------------------------------
SUBROUTINE thickcurr_f_destroy_sensors(sensor_ptr) BIND(C,NAME="thickcurr_f_destroy_sensors")
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Pointer to sensor object
TYPE(thickcurr_sensors), POINTER :: sensors
IF(.NOT.c_associated(sensor_ptr))RETURN
CALL c_f_pointer(sensor_ptr,sensors)
IF(ASSOCIATED(sensors%names))DEALLOCATE(sensors%names)
IF(ASSOCIATED(sensors%wa))DEALLOCATE(sensors%wa,sensors%idxa)
IF(ASSOCIATED(sensors%wb))DEALLOCATE(sensors%wb,sensors%idxb)
IF(ASSOCIATED(sensors%aprobe%cells))DEALLOCATE(sensors%aprobe%cells,sensors%aprobe%pts,sensors%aprobe%pts_log)
IF(ASSOCIATED(sensors%bprobe%cells))DEALLOCATE(sensors%bprobe%cells,sensors%bprobe%pts,sensors%bprobe%pts_log)
DEALLOCATE(sensors)
END SUBROUTINE thickcurr_f_destroy_sensors
!---------------------------------------------------------------------------------
!> Get H(Curl) FE representation from ThickCurr object
!---------------------------------------------------------------------------------
FUNCTION get_hcurl_rep(self,hcurl_rep,error_str) RESULT(success)
TYPE(thickcurr_py), INTENT(in) :: self !< ThickCurr object
CLASS(oft_hcurl_fem), POINTER, INTENT(out) :: hcurl_rep !< H(Curl) FE representation
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
success=.FALSE.
SELECT TYPE(this=>self%ML_hcurl_grad%current_level%fields(1)%fe)
  CLASS IS(oft_hcurl_fem)
    hcurl_rep=>this
    success=.TRUE.
  CLASS DEFAULT
    CALL copy_string('Invalid H(Curl) FE object',error_str)
END SELECT
END FUNCTION get_hcurl_rep
!---------------------------------------------------------------------------------
!> Evaluate tabulated field at a quadrature point
!---------------------------------------------------------------------------------
SUBROUTINE tab_interp_apply(self,cell,f,gop,val)
CLASS(thickcurr_tab_interp), INTENT(inout) :: self !< Interpolation object
INTEGER(i4), INTENT(in) :: cell !< Cell for interpolation
REAL(r8), INTENT(in) :: f(:) !< Position in cell in logical coordinates
REAL(r8), INTENT(in) :: gop(3,4) !< Logical gradient vectors at f
REAL(r8), INTENT(out) :: val(:) !< Field value at f
INTEGER(i4) :: m
DO m=1,self%nq
  IF(ALL(ABS(f(1:4)-self%qpts(1:4,m))<1.d-14))THEN
    val=self%vals(:,m,cell)
    RETURN
  END IF
END DO
CALL oft_abort('Point is not a tabulated quadrature point','tab_interp_apply',__FILE__)
END SUBROUTINE tab_interp_apply
!---------------------------------------------------------------------------------
!> Evaluate sensor signals
!---------------------------------------------------------------------------------
SUBROUTINE sensors_eval(self,a,signals)
CLASS(thickcurr_sensors), INTENT(inout) :: self !< Sensor object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Vector potential
REAL(r8), INTENT(out) :: signals(:) !< Sensor signals
INTEGER(i4) :: i
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vals
signals=0.d0
IF(self%npa>0)THEN
  self%ainterp%u=>a
  CALL self%ainterp%setup(thickcurr_ML_hcurl_grad%current_level)
  ALLOCATE(vals(3,self%npa))
  CALL self%aprobe%eval(vals)
  DO i=1,self%npa
    signals(self%idxa(i))=signals(self%idxa(i))+DOT_PRODUCT(self%wa(:,i),vals(:,i))
  END DO
  DEALLOCATE(vals)
END IF
IF(self%npb>0)THEN
  self%binterp%u=>a
  CALL self%binterp%setup(thickcurr_ML_hcurl_grad%current_level)
  ALLOCATE(vals(3,self%npb))
  CALL self%bprobe%eval(vals)
  DO i=1,self%npb
    signals(self%idxb(i))=signals(self%idxb(i))+DOT_PRODUCT(self%wb(:,i),vals(:,i))
  END DO
  DEALLOCATE(vals)
END IF
END SUBROUTINE sensors_eval
!---------------------------------------------------------------------------------
!> Sample sensors and save to history file
!---------------------------------------------------------------------------------
SUBROUTINE sensors_apply(self,a,t)
CLASS(thickcurr_sensors), INTENT(inout) :: self !< Sensor object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Vector potential at current time
REAL(r8), INTENT(in) :: t !< Current time
INTEGER(i4) :: i
REAL(r8), ALLOCATABLE, DIMENSION(:) :: signals
!---Setup history file if necessary
IF(.NOT.self%hist_ready)THEN
  IF(oft_env%head_proc)THEN
    CALL self%hist_file%setup(TRIM(self%filename), desc="ThickCurr sensor signals")
    CALL self%hist_file%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,self%nsensors
      CALL self%hist_file%add_field(TRIM(self%names(i)), 'r8')
    END DO
    CALL self%hist_file%write_header
  END IF
  self%hist_ready=.TRUE.
END IF
!---Sample signals
ALLOCATE(signals(self%nsensors))
CALL sensors_eval(self,a,signals)
IF(oft_env%head_proc)THEN
  CALL self%hist_file%open
  CALL self%hist_file%write(data_r8=[t,signals])
  CALL self%hist_file%close
END IF
DEALLOCATE(signals)
END SUBROUTINE sensors_apply
!---------------------------------------------------------------------------------
!> Flush history file (file is closed after each sample)
!---------------------------------------------------------------------------------
SUBROUTINE sensors_flush(self)
CLASS(thickcurr_sensors), INTENT(inout) :: self !< Sensor object
END SUBROUTINE sensors_flush
END MODULE thickcurr_f
