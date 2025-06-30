!------------------------------------------------------------------------------
!> @file tokamaker_f.F90
!
!> Fortran part of Python wrapper for Grad-Shafranov functionality
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_python
!------------------------------------------------------------------------------
MODULE tokamaker_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
  c_f_pointer, c_bool, c_null_ptr, c_associated
USE oft_base
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell
USE multigrid, ONLY: multigrid_mesh, multigrid_reset
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_base, ONLY: oft_afem_type, oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, &
  oft_lag_setup
USE oft_blag_operators, ONLY: oft_lag_brinterp, oft_lag_bginterp, oft_blag_project
USE mhd_utils, ONLY: mu0
USE axi_green, ONLY: green
USE oft_gs, ONLY: gs_eq, gs_save_fields, gs_setup_walls, build_dels, &
  gs_fixed_vflux, gs_get_qprof, gs_trace_surf, gs_b_interp, gs_j_interp, gs_prof_interp, &
  gs_plasma_mutual, gs_source, gs_err_reason, gs_coil_source_distributed, gs_vacuum_solve, &
  gs_coil_mutual, gs_coil_mutual_distributed, gs_project_b, gs_save_mug, gs_update_bounds
#ifdef OFT_TOKAMAKER_LEGACY
USE oft_gs, ONLY: gs_load_regions
#endif
USE oft_gs_util, ONLY: gs_comp_globals, gs_save_eqdsk, gs_save_ifile, gs_profile_load, gs_profile_save, &
  sauter_fc, gs_calc_vloop
USE oft_gs_fit, ONLY: fit_gs, fit_pm
USE oft_gs_td, ONLY: oft_tmaker_td, eig_gs_td
USE oft_gs_mercier, ONLY: create_dipole_b0_prof
USE diagnostic, ONLY: bscal_surf_int
USE oft_base_f, ONLY: copy_string, copy_string_rev, oftpy_init
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
TYPE, BIND(C) :: tokamaker_settings_type
  LOGICAL(KIND=c_bool) :: pm = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: free_boundary = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: has_plasma = .TRUE. !< Needs docs
  LOGICAL(KIND=c_bool) :: limited_only = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: dipole_mode = .FALSE. !< Needs docs
  INTEGER(KIND=c_int) :: maxits = 40 !< Needs docs
  INTEGER(KIND=c_int) :: mode = 1 !< Needs docs
  REAL(KIND=c_double) :: urf = 0.3d0 !< Needs docs
  REAL(KIND=c_double) :: nl_tol = 1.d-6 !< Needs docs
  REAL(KIND=c_double) :: rmin = 0.d0 !< Needs docs
  REAL(KIND=c_double) :: lim_zmax = 1.d99 !< Needs docs
  TYPE(c_ptr) :: limiter_file !< Needs docs
END TYPE tokamaker_settings_type
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
TYPE, BIND(C) :: tokamaker_recon_settings_type
  LOGICAL(KIND=c_bool) :: fitI = .TRUE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitP = .TRUE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitPnorm = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitAlam = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitR0 = .TRUE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitV0 = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitCoils = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fitF0 = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: fixedCentering = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: pm = .FALSE. !< Needs docs
  TYPE(c_ptr) :: infile !< Needs docs
  TYPE(c_ptr) :: outfile !< Needs docs
END TYPE tokamaker_recon_settings_type
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
TYPE :: tokamaker_instance
  INTEGER(i4), POINTER, DIMENSION(:) :: reg_plot => NULL() !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lc_plot => NULL() !< Needs docs
  REAL(r8), POINTER, DIMENSION(:,:) :: r_plot => NULL() !< Needs docs
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL() !< Mesh container
  TYPE(oft_ml_fem_type), POINTER :: ML_oft_blagrange => NULL() !< Finite element container
  TYPE(gs_eq), POINTER :: gs => NULL() !< G-S object
  TYPE(oft_tmaker_td), POINTER :: gs_td => NULL() !< Time-dependent G-S object
END TYPE tokamaker_instance
CONTAINS
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_alloc(tMaker_ptr,mesh_ptr,error_str) BIND(C,NAME="tokamaker_alloc")
TYPE(c_ptr), INTENT(out) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Needs docs
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
CALL copy_string('',error_str)
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
ALLOCATE(tMaker_obj)
ALLOCATE(tMaker_obj%gs,tMaker_obj%gs_td,tMaker_obj%ML_oft_blagrange)
tMaker_ptr=C_LOC(tMaker_obj)
CALL c_f_pointer(mesh_ptr,tMaker_obj%ml_mesh)
END SUBROUTINE tokamaker_alloc
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
FUNCTION tokamaker_ccast(tMaker_cptr,tMaker_obj,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: tMaker_cptr !< C pointer to TokaMaker object
TYPE(tokamaker_instance), POINTER, INTENT(out) :: tMaker_obj
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
!---Clear error flag
IF(PRESENT(error_str))CALL copy_string('',error_str)
IF(.NOT.c_associated(tMaker_cptr))THEN
  IF(PRESENT(error_str))CALL copy_string('TokaMaker object not associated',error_str)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(tMaker_cptr,tMaker_obj)
success=.TRUE.
END FUNCTION tokamaker_ccast
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_eval_green(n,r,z,rc,zc,vals) BIND(C,NAME="tokamaker_eval_green")
INTEGER(c_int), VALUE, INTENT(in) :: n !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: rc !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: zc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: r !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: z !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
INTEGER(4) :: i
REAL(r8), POINTER, DIMENSION(:) :: rtmp,ztmp,vals_tmp
CALL c_f_pointer(r, rtmp, [n])
CALL c_f_pointer(z, ztmp, [n])
CALL c_f_pointer(vals, vals_tmp, [n])
!$omp parallel do if(n>10)
DO i=1,n
  vals_tmp(i)=green(rtmp(i),ztmp(i),rc,zc)
END DO
END SUBROUTINE tokamaker_eval_green
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup_regions(tMaker_ptr,coil_file,reg_eta,contig_flag,xpoint_mask,coil_nturns,ncoils,error_str) BIND(C,NAME="tokamaker_setup_regions")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
CHARACTER(KIND=c_char), INTENT(in) :: coil_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_eta !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: contig_flag !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: xpoint_mask !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_nturns !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ncoils !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
real(r8), POINTER :: eta_tmp(:),nturns_tmp(:,:)
INTEGER(i4), POINTER :: contig_tmp(:)
INTEGER(4) :: i
INTEGER(4), POINTER :: xpoint_tmp(:)
! TYPE(multigrid_mesh), POINTER :: mg_mesh
TYPE(tokamaker_instance), POINTER :: tMaker_obj
CLASS(oft_bmesh), POINTER :: smesh
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
smesh=>tMaker_obj%ml_mesh%smesh
CALL copy_string_rev(coil_file,tMaker_obj%gs%coil_file)
IF(TRIM(tMaker_obj%gs%coil_file)=='none')THEN
  !
  CALL c_f_pointer(xpoint_mask, xpoint_tmp, [smesh%nreg])
  ALLOCATE(tMaker_obj%gs%saddle_rmask(smesh%nreg))
  tMaker_obj%gs%saddle_rmask=LOGICAL(xpoint_tmp==0)
  !
  CALL c_f_pointer(reg_eta, eta_tmp, [smesh%nreg])
  tMaker_obj%gs%ncoil_regs=0
  tMaker_obj%gs%ncond_regs=0
  DO i=2,smesh%nreg
    IF(eta_tmp(i)>0.d0)THEN
      tMaker_obj%gs%ncond_regs=tMaker_obj%gs%ncond_regs+1
    ELSE
      tMaker_obj%gs%ncoil_regs=tMaker_obj%gs%ncoil_regs+1
    END IF
  END DO
  !
  tMaker_obj%gs%ncoils=ncoils
  ALLOCATE(tMaker_obj%gs%coil_vcont(ncoils),tMaker_obj%gs%coil_currs(ncoils))
  tMaker_obj%gs%coil_vcont = 0.d0
  tMaker_obj%gs%coil_currs = 0.d0
  CALL c_f_pointer(coil_nturns, nturns_tmp, [smesh%nreg,ncoils])
  ALLOCATE(tMaker_obj%gs%coil_nturns(smesh%nreg,ncoils))
  tMaker_obj%gs%coil_nturns=0.d0
  tMaker_obj%gs%coil_nturns=nturns_tmp
  !
  CALL c_f_pointer(contig_flag, contig_tmp, [smesh%nreg])
  ALLOCATE(tMaker_obj%gs%cond_regions(tMaker_obj%gs%ncond_regs))
  ALLOCATE(tMaker_obj%gs%coil_regions(tMaker_obj%gs%ncoil_regs))
  tMaker_obj%gs%ncond_regs=0
  tMaker_obj%gs%ncoil_regs=0
  DO i=2,smesh%nreg
    IF(eta_tmp(i)>0.d0)THEN
      tMaker_obj%gs%ncond_regs=tMaker_obj%gs%ncond_regs+1
      IF(eta_tmp(i)<1.d9)THEN
        tMaker_obj%gs%cond_regions(tMaker_obj%gs%ncond_regs)%eta=eta_tmp(i)
      END IF
      tMaker_obj%gs%cond_regions(tMaker_obj%gs%ncond_regs)%id=i
      tMaker_obj%gs%cond_regions(tMaker_obj%gs%ncond_regs)%continuous=(contig_tmp(i)==1)
      IF(tMaker_obj%gs%dipole_mode)tMaker_obj%gs%cond_regions(tMaker_obj%gs%ncond_regs)%inner_limiter=(contig_tmp(i)==-1)
    ELSE
      tMaker_obj%gs%ncoil_regs=tMaker_obj%gs%ncoil_regs+1
      tMaker_obj%gs%coil_regions(tMaker_obj%gs%ncoil_regs)%id=i
    END IF
  END DO
ELSE
#ifdef OFT_TOKAMAKER_LEGACY
  CALL gs_load_regions(tMaker_obj%gs)
#else
  CALL oft_abort("OFT not compiled with legacy TokaMaker support","tokamaker_setup_regions",__FILE__)
#endif
END IF
END SUBROUTINE tokamaker_setup_regions
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_destroy(tMaker_ptr,error_str) BIND(C,NAME="tokamaker_destroy")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,ierr,io_unit,npts,iostat
REAL(8) :: theta
LOGICAL :: file_exists
real(r8), POINTER :: vals_tmp(:)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
!---Destroy objects
IF(ASSOCIATED(tMaker_obj%r_plot))DEALLOCATE(tMaker_obj%r_plot)
IF(ASSOCIATED(tMaker_obj%lc_plot))DEALLOCATE(tMaker_obj%lc_plot)
IF(ASSOCIATED(tMaker_obj%reg_plot))DEALLOCATE(tMaker_obj%reg_plot)
IF(ASSOCIATED(tMaker_obj%gs))THEN
  CALL tMaker_obj%gs%delete()
  DEALLOCATE(tMaker_obj%gs)
END IF
IF(ASSOCIATED(tMaker_obj%ML_oft_blagrange))THEN
  CALL tMaker_obj%ML_oft_blagrange%current_level%delete()
END IF
IF(ASSOCIATED(tMaker_obj%ml_mesh))THEN
  CALL multigrid_reset(tMaker_obj%ml_mesh)
  DEALLOCATE(tMaker_obj%ml_mesh)
END IF
DEALLOCATE(tMaker_obj)
END SUBROUTINE tokamaker_destroy
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup(tMaker_ptr,order,full_domain,ncoils,error_str) BIND(C,NAME="tokamaker_setup")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: full_domain !< Needs docs
INTEGER(KIND=c_int), INTENT(out) :: ncoils !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,ierr,io_unit,npts,iostat
REAL(8) :: theta
LOGICAL :: file_exists
real(r8), POINTER :: vals_tmp(:)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
!------------------------------------------------------------------------------
! Check input files
!------------------------------------------------------------------------------
IF(TRIM(tMaker_obj%gs%coil_file)/='none')THEN
    INQUIRE(EXIST=file_exists,FILE=TRIM(tMaker_obj%gs%coil_file))
    IF(.NOT.file_exists)THEN
        CALL copy_string('Specified "coil_file" cannot be found',error_str)
        RETURN
    END IF
END IF
IF(TRIM(tMaker_obj%gs%limiter_file)/='none')THEN
    INQUIRE(EXIST=file_exists,FILE=TRIM(tMaker_obj%gs%limiter_file))
    IF(.NOT.file_exists)THEN
        CALL copy_string('Specified "limiter_file" cannot be found',error_str)
        RETURN
    END IF
END IF
! IF(TRIM(tMaker_obj%gs%eqdsk_limiter_file)/='none')THEN
!     INQUIRE(EXIST=file_exists,FILE=TRIM(tMaker_obj%gs%eqdsk_limiter_file))
!     IF(.NOT.file_exists)THEN
!         CALL copy_string('Specified "eqdsk_limiter_file" cannot be found',error_str)
!         RETURN
!     END IF
! END IF
!------------------------------------------------------------------------------
! Setup Lagrange Elements
!------------------------------------------------------------------------------
tMaker_obj%ml_mesh%smesh%tess_order=order
CALL oft_lag_setup(tMaker_obj%ml_mesh,order,ML_blag_obj=tMaker_obj%ML_oft_blagrange,minlev=-1)
CALL tMaker_obj%gs%setup(tMaker_obj%ML_oft_blagrange)
!------------------------------------------------------------------------------
! Setup experimental geometry
!------------------------------------------------------------------------------
tMaker_obj%gs%save_visit=.FALSE.
tMaker_obj%gs%full_domain=full_domain
CALL gs_setup_walls(tMaker_obj%gs)
CALL tMaker_obj%gs%load_limiters
CALL tMaker_obj%gs%init()
IF(tMaker_obj%gs%dipole_mode)THEN
  tMaker_obj%gs%dipole_a=0.d0
  CALL create_dipole_b0_prof(tMaker_obj%gs%dipole_B0,64)
END IF
ncoils=tMaker_obj%gs%ncoils
END SUBROUTINE tokamaker_setup
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_load_profiles(tMaker_ptr,f_file,f_offset,p_file,eta_file,f_NI_file,error_str) BIND(C,NAME="tokamaker_load_profiles")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
CHARACTER(KIND=c_char), INTENT(in) :: f_file(OFT_PATH_SLEN) !< F*F' prof.in file
CHARACTER(KIND=c_char), INTENT(in) :: p_file(OFT_PATH_SLEN) !< P' prof.in file
REAL(c_double), VALUE, INTENT(in) :: f_offset !< Vacuum F_0 value (must be > -1E98 to update)
CHARACTER(KIND=c_char), INTENT(in) :: eta_file(OFT_PATH_SLEN) !< Resistivity (eta) prof.in file
CHARACTER(KIND=c_char), INTENT(in) :: f_NI_file(OFT_PATH_SLEN) !< Non-inductive F*F' prof.in file
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
CHARACTER(LEN=OFT_PATH_SLEN) :: tmp_str
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL copy_string_rev(f_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,tMaker_obj%gs%I)
IF(f_offset>-1.d98)tMaker_obj%gs%I%f_offset=f_offset
CALL copy_string_rev(p_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,tMaker_obj%gs%P)
CALL copy_string_rev(eta_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,tMaker_obj%gs%eta)
CALL copy_string_rev(f_NI_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,tMaker_obj%gs%I_NI)
END SUBROUTINE tokamaker_load_profiles
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_init_psi(tMaker_ptr,r0,z0,a,kappa,delta,rhs_source,error_str) BIND(C,NAME="tokamaker_init_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), VALUE, INTENT(in) :: r0 !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: z0 !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: a !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: kappa !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: delta !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: rhs_source !< Current source term (optional)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: ierr
REAL(8), POINTER, DIMENSION(:) :: rhs_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(c_associated(rhs_source))THEN
  CALL c_f_pointer(rhs_source, rhs_tmp, [tMaker_obj%gs%psi%n])
  CALL tMaker_obj%gs%init_psi(ierr,curr_source=rhs_tmp)
ELSE
  CALL tMaker_obj%gs%init_psi(ierr,r0=[r0,z0],a=a,kappa=kappa,delta=delta)
END IF
IF(ierr/=0)CALL copy_string(gs_err_reason(ierr),error_str)
END SUBROUTINE tokamaker_init_psi
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_solve(tMaker_ptr,vacuum,error_str) BIND(C,NAME="tokamaker_solve")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
LOGICAL(c_bool), VALUE, INTENT(in) :: vacuum !< Perform vacuum solve?
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: ierr
LOGICAL :: vac_save
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
tMaker_obj%gs%timing=0.d0
IF(vacuum)THEN
  vac_save=tMaker_obj%gs%has_plasma
  tMaker_obj%gs%has_plasma=.FALSE.
END IF
CALL tMaker_obj%gs%solve(ierr)
IF(vacuum)tMaker_obj%gs%has_plasma=vac_save
IF(ierr/=0)CALL copy_string(gs_err_reason(ierr),error_str)
END SUBROUTINE tokamaker_solve
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_vac_solve(tMaker_ptr,psi_in,rhs_source,error_str) BIND(C,NAME="tokamaker_vac_solve")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: psi_in !< Input: BCs for \f$ \psi \f$, Output: solution
TYPE(c_ptr), VALUE, INTENT(in) :: rhs_source !< Current source term (optional)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: ierr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp,rhs_tmp
CLASS(oft_vector), POINTER :: psi_tmp,rhs_vec
TYPE(oft_lag_brinterp) :: source_field
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
NULLIFY(psi_tmp)
CALL tMaker_obj%gs%psi%new(psi_tmp)
CALL c_f_pointer(psi_in, vals_tmp, [tMaker_obj%gs%psi%n])
CALL psi_tmp%restore_local(vals_tmp)
IF(c_associated(rhs_source))THEN
  NULLIFY(rhs_tmp)
  CALL tMaker_obj%gs%psi%new(rhs_vec)
  CALL c_f_pointer(rhs_source, rhs_tmp, [tMaker_obj%gs%psi%n])
  CALL rhs_vec%restore_local(rhs_tmp)
  source_field%u=>rhs_vec
  CALL source_field%setup(tMaker_obj%gs%fe_rep)
  CALL tMaker_obj%gs%vac_solve(psi_tmp,rhs_source=source_field,ierr=ierr)
  CALL source_field%delete()
  CALL rhs_vec%delete()
  DEALLOCATE(rhs_vec)
ELSE
  CALL tMaker_obj%gs%vac_solve(psi_tmp,ierr=ierr)
END IF
CALL psi_tmp%get_local(vals_tmp)
CALL psi_tmp%delete()
DEALLOCATE(psi_tmp)
END SUBROUTINE tokamaker_vac_solve
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_recon_run(tMaker_ptr,vacuum,settings,error_flag) BIND(C,NAME="tokamaker_recon_run")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
LOGICAL(c_bool), VALUE, INTENT(in) :: vacuum !< Needs docs
TYPE(tokamaker_recon_settings_type), INTENT(in) :: settings !< Needs docs
INTEGER(c_int), INTENT(out) :: error_flag !< Needs docs
LOGICAL :: fitI,fitP,fitPnorm,fitAlam,fitR0,fitV0,fitCoils,fitF0,fixedCentering
CHARACTER(KIND=c_char), POINTER, DIMENSION(:) :: infile_c,outfile_c
CHARACTER(LEN=OFT_PATH_SLEN) :: infile,outfile
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj))THEN
  error_flag=-100
  RETURN
END IF
error_flag=0
IF(vacuum)tMaker_obj%gs%has_plasma=.FALSE.
fitI=settings%fitI
fitP=settings%fitP
fitPnorm=settings%fitPnorm
fitAlam=settings%fitAlam
fitR0=settings%fitR0
fitV0=settings%fitV0
fitCoils=settings%fitCoils
fitF0=settings%fitF0
fixedCentering=settings%fixedCentering
fit_pm=settings%pm
CALL c_f_pointer(settings%infile,infile_c,[OFT_PATH_SLEN])
CALL c_f_pointer(settings%outfile,outfile_c,[OFT_PATH_SLEN])
CALL copy_string_rev(infile_c,infile)
CALL copy_string_rev(outfile_c,outfile)
tMaker_obj%gs%timing=0.d0
CALL fit_gs(tMaker_obj%gs,infile,outfile,fitI,fitP,fitPnorm,&
            fitAlam,fitR0,fitV0,fitCoils,fitF0, &
            fixedCentering)
CALL gs_profile_save(TRIM(outfile)//'_fprof',tMaker_obj%gs%I)
CALL gs_profile_save(TRIM(outfile)//'_pprof',tMaker_obj%gs%P)
tMaker_obj%gs%has_plasma=.TRUE.
END SUBROUTINE tokamaker_recon_run
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup_td(tMaker_ptr,dt,lin_tol,nl_tol,pre_plasma,error_str) BIND(C,NAME="tokamaker_setup_td")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), VALUE, INTENT(in) :: dt !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: lin_tol !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: nl_tol !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pre_plasma !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(ASSOCIATED(tMaker_obj%gs_td))THEN
  CALL tMaker_obj%gs_td%delete()
  DEALLOCATE(tMaker_obj%gs_td)
END IF
ALLOCATE(tMaker_obj%gs_td)
CALL tMaker_obj%gs_td%setup(tMaker_obj%gs,dt,lin_tol,nl_tol,LOGICAL(pre_plasma))
END SUBROUTINE tokamaker_setup_td
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_eig_td(tMaker_ptr,omega,neigs,eigs,eig_vecs,include_bounds,eta_plasma,pm,error_str) BIND(C,NAME="tokamaker_eig_td")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), VALUE, INTENT(in) :: omega !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vecs !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: include_bounds !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: eta_plasma !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pm !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER :: eigs_tmp(:,:),eig_vecs_tmp(:,:)
LOGICAL :: pm_save
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(eigs, eigs_tmp, [2,neigs])
CALL c_f_pointer(eig_vecs, eig_vecs_tmp, [tMaker_obj%gs%psi%n,neigs])
pm_save=oft_env%pm; oft_env%pm=pm
CALL eig_gs_td(tMaker_obj%gs,neigs,eigs_tmp,eig_vecs_tmp,omega,LOGICAL(include_bounds),eta_plasma)
oft_env%pm=pm_save
IF((eigs_tmp(1,1)<-1.d98).AND.(eigs_tmp(2,1)<-1.d98))CALL copy_string('Error in eigenvalue solve',error_str)
END SUBROUTINE tokamaker_eig_td
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_eig_wall(tMaker_ptr,neigs,eigs,eig_vecs,pm,error_str) BIND(C,NAME="tokamaker_eig_wall")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vecs !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pm !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8) :: alam_save,pnorm_save
REAL(8), POINTER :: eigs_tmp(:,:),eig_vecs_tmp(:,:)
LOGICAL :: pm_save
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(eigs, eigs_tmp, [2,neigs])
CALL c_f_pointer(eig_vecs, eig_vecs_tmp, [tMaker_obj%gs%psi%n,neigs])
alam_save=tMaker_obj%gs%alam; tMaker_obj%gs%alam=0.d0
pnorm_save=tMaker_obj%gs%pnorm; tMaker_obj%gs%pnorm=0.d0
pm_save=oft_env%pm; oft_env%pm=pm
CALL eig_gs_td(tMaker_obj%gs,neigs,eigs_tmp,eig_vecs_tmp,0.d0,.FALSE.,-1.d0)
oft_env%pm=pm_save
IF((eigs_tmp(1,1)<-1.d98).AND.(eigs_tmp(2,1)<-1.d98))CALL copy_string('Error in eigenvalue solve',error_str)
tMaker_obj%gs%alam=alam_save
tMaker_obj%gs%pnorm=pnorm_save
END SUBROUTINE tokamaker_eig_wall
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_step_td(tMaker_ptr,time,dt,nl_its,lin_its,nretry,error_str) BIND(C,NAME="tokamaker_step_td")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), INTENT(inout) :: time !< Needs docs
REAL(c_double), INTENT(inout) :: dt !< Needs docs
INTEGER(c_int), INTENT(out) :: nl_its !< Needs docs
INTEGER(c_int), INTENT(out) :: lin_its !< Needs docs
INTEGER(c_int), INTENT(out) :: nretry !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL tMaker_obj%gs_td%step(time,dt,nl_its,lin_its,nretry)
END SUBROUTINE tokamaker_step_td
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_mesh(tMaker_ptr,np,r_loc,nc,lc_loc,reg_loc,error_str) BIND(C,NAME="tokamaker_get_mesh")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), INTENT(out) :: lc_loc !< Needs docs
TYPE(c_ptr), INTENT(out) :: r_loc !< Needs docs
TYPE(c_ptr), INTENT(out) :: reg_loc !< Needs docs
INTEGER(c_int), INTENT(out) :: np !< Needs docs
INTEGER(c_int), INTENT(out) :: nc !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,j,k,id
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL tMaker_obj%gs%mesh%tessellate(tMaker_obj%r_plot, tMaker_obj%lc_plot, tMaker_obj%gs%mesh%tess_order)
np=SIZE(tMaker_obj%r_plot,DIM=2,KIND=c_int)
nc=SIZE(tMaker_obj%lc_plot,DIM=2,KIND=c_int)
r_loc=c_loc(tMaker_obj%r_plot)
lc_loc=c_loc(tMaker_obj%lc_plot)
!
ALLOCATE(tMaker_obj%reg_plot(nc))
k=nc/tMaker_obj%gs%mesh%nc
IF(ASSOCIATED(tMaker_obj%gs%mesh%reg))THEN
  !$omp parallel do private(j,id)
  DO i=1,tMaker_obj%gs%mesh%nc
    id=tMaker_obj%gs%mesh%reg(i)
    DO j=1,k
      tMaker_obj%reg_plot((i-1)*k+j)=id
    END DO
  END DO
ELSE
  tMaker_obj%reg_plot=0
END IF
reg_loc=c_loc(tMaker_obj%reg_plot)
END SUBROUTINE tokamaker_get_mesh
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_limiter(tMaker_ptr,np,r_loc,nloops,loop_ptr,error_str) BIND(C,NAME="tokamaker_get_limiter")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), INTENT(out) :: r_loc !< Needs docs
TYPE(c_ptr), INTENT(out) :: loop_ptr !< Needs docs
INTEGER(c_int), INTENT(out) :: np !< Needs docs
INTEGER(c_int), INTENT(out) :: nloops !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i
REAL(8), POINTER, DIMENSION(:,:) :: r_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
np=tMaker_obj%gs%nlim_con
ALLOCATE(r_tmp(2,tMaker_obj%gs%nlim_con))
r_loc=C_LOC(r_tmp)
DO i=1,tMaker_obj%gs%nlim_con
  r_tmp(:,i)=tMaker_obj%gs%mesh%r(1:2,tMaker_obj%gs%lim_con(i))
END DO
nloops=tMaker_obj%gs%lim_nloops
loop_ptr=C_LOC(tMaker_obj%gs%lim_ptr)
END SUBROUTINE tokamaker_get_limiter
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_psi(tMaker_ptr,psi_vals,psi_lim,psi_max,error_str) BIND(C,NAME="tokamaker_get_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(c_double), INTENT(out) :: psi_lim !< Needs docs
REAL(c_double), INTENT(out) :: psi_max !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(psi_vals, vals_tmp, [tMaker_obj%gs%psi%n])
CALL tMaker_obj%gs%psi%get_local(vals_tmp)
psi_lim = tMaker_obj%gs%plasma_bounds(1)
psi_max = tMaker_obj%gs%plasma_bounds(2)
END SUBROUTINE tokamaker_get_psi
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_dels_curr(tMaker_ptr,psi_vals,error_str) BIND(C,NAME="tokamaker_get_dels_curr")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_solver), POINTER :: minv
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(.NOT.ASSOCIATED(tMaker_obj%gs%dels_full))CALL build_dels(tMaker_obj%gs%dels_full,tMaker_obj%gs,"none")
!
CALL tMaker_obj%gs%psi%new(u)
CALL tMaker_obj%gs%psi%new(v)
CALL c_f_pointer(psi_vals, vals_tmp, [tMaker_obj%gs%psi%n])
CALL u%restore_local(vals_tmp)
!
NULLIFY(minv)
CALL create_cg_solver(minv)
minv%A=>tMaker_obj%gs%mop
minv%its=-2
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
CALL tMaker_obj%gs%dels_full%apply(u,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vals_tmp)
!
CALL u%delete()
CALL v%delete()
CALL minv%pre%delete()
CALL minv%delete()
DEALLOCATE(u,v,minv)
END SUBROUTINE tokamaker_get_dels_curr
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_jtor(tMaker_ptr,jtor,error_str) BIND(C,NAME="tokamaker_get_jtor")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: jtor !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_solver), POINTER :: minv
TYPE(gs_j_interp) :: j_interp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(.NOT.ASSOCIATED(tMaker_obj%gs%dels_full))CALL build_dels(tMaker_obj%gs%dels_full,tMaker_obj%gs,"none")
!
CALL tMaker_obj%gs%psi%new(u)
CALL tMaker_obj%gs%psi%new(v)
!
NULLIFY(minv)
CALL create_cg_solver(minv)
minv%A=>tMaker_obj%gs%mop
minv%its=-2
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
CALL j_interp%setup(tMaker_obj%gs)
CALL oft_blag_project(tMaker_obj%gs%fe_rep,j_interp,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL c_f_pointer(jtor, vals_tmp, [tMaker_obj%gs%psi%n])
CALL u%get_local(vals_tmp)
!
CALL j_interp%delete()
CALL u%delete()
CALL v%delete()
CALL minv%pre%delete()
CALL minv%delete()
DEALLOCATE(u,v,minv)
END SUBROUTINE tokamaker_get_jtor
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_area_int(tMaker_ptr,vec_vals,reg_ind,result,error_str) BIND(C,NAME="tokamaker_area_int")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: vec_vals !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: reg_ind !< Needs docs
REAL(c_double), INTENT(out) :: result !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,m
real(8) :: goptmp(3,3),v,pt(3),valtmp(1)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CLASS(oft_vector), POINTER :: u
TYPE(oft_lag_brinterp) :: field
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
NULLIFY(field%u)
CALL tMaker_obj%gs%psi%new(field%u)
CALL c_f_pointer(vec_vals, vals_tmp, [tMaker_obj%gs%psi%n])
CALL field%u%restore_local(vals_tmp)
CALL field%setup(tMaker_obj%gs%fe_rep)
IF(reg_ind>0)THEN
  result = bscal_surf_int(tMaker_obj%gs%mesh,field,tMaker_obj%gs%fe_rep%quad%order,reg_ind)
ELSE
  result = bscal_surf_int(tMaker_obj%gs%mesh,field,tMaker_obj%gs%fe_rep%quad%order)
END IF
CALL field%u%delete
DEALLOCATE(field%u)
CALL field%delete
END SUBROUTINE tokamaker_area_int
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_flux_int(tMaker_ptr,psi_vals,field_vals,nvals,result,error_str) BIND(C,NAME="tokamaker_flux_int")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: field_vals !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nvals
REAL(c_double), INTENT(out) :: result !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(i4) :: i,m
REAL(r8) :: area,psitmp(1),sgop(3,3)
REAL(8), POINTER, DIMENSION(:) :: psi_tmp,field_tmp
TYPE(gs_prof_interp) :: prof_interp_obj
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
DEBUG_STACK_PUSH
!---Setup
CALL c_f_pointer(psi_vals, psi_tmp, [nvals])
CALL c_f_pointer(field_vals, field_tmp, [nvals])
result=0.d0
prof_interp_obj%mode=4
CALL prof_interp_obj%setup(tMaker_obj%gs)
!$omp parallel do private(m,psitmp,sgop,area) reduction(+:result)
do i=1,tMaker_obj%gs%mesh%nc
  IF(tMaker_obj%gs%mesh%reg(i)/=1)CYCLE
  !---Loop over quadrature points
  do m=1,tMaker_obj%gs%fe_rep%quad%np
    call tMaker_obj%gs%mesh%jacobian(i,tMaker_obj%gs%fe_rep%quad%pts(:,m),sgop,area)
    call prof_interp_obj%interp(i,tMaker_obj%gs%fe_rep%quad%pts(:,m),sgop,psitmp)
    psitmp(1)=linterp(psi_tmp,field_tmp,nvals,psitmp(1),0)
    IF(psitmp(1)>-1.d98)result = result + psitmp(1)*area*tMaker_obj%gs%fe_rep%quad%wts(m)
  end do
end do
!---Global reduction and cleanup
result=oft_mpi_sum(result)
CALL prof_interp_obj%delete()
DEBUG_STACK_POP
END SUBROUTINE tokamaker_flux_int
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_coil_currents(tMaker_ptr,currents,reg_currents,error_str) BIND(C,NAME="tokamaker_get_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: currents !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_currents !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,j
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp,coil_regs
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(reg_currents, coil_regs, [tMaker_obj%gs%mesh%nreg])
CALL c_f_pointer(currents, vals_tmp, [tMaker_obj%gs%ncoils])
vals_tmp=(tMaker_obj%gs%coil_currs + tMaker_obj%gs%coil_vcont*tMaker_obj%gs%vcontrol_val)/mu0
coil_regs = 0.d0
DO j=1,tMaker_obj%gs%ncoil_regs
  DO i=1,tMaker_obj%gs%ncoils
    coil_regs(tMaker_obj%gs%coil_regions(j)%id) = coil_regs(tMaker_obj%gs%coil_regions(j)%id) &
      + vals_tmp(i)*tMaker_obj%gs%coil_nturns(tMaker_obj%gs%coil_regions(j)%id,i)
  END DO
  coil_regs(tMaker_obj%gs%coil_regions(j)%id) = coil_regs(tMaker_obj%gs%coil_regions(j)%id)*tMaker_obj%gs%coil_regions(j)%area
END DO
END SUBROUTINE tokamaker_get_coil_currents
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_coil_Lmat(tMaker_ptr,Lmat,error_str) BIND(C,NAME="tokamaker_get_coil_Lmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), VALUE, INTENT(in) :: Lmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
REAL(8) :: tmp1,tmp2,tmp3,itor
CLASS(oft_vector), POINTER :: rhs,vec1,vec2
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
!---Update plasma row/column
IF(tMaker_obj%gs%has_plasma)THEN
  DO i=1,tMaker_obj%gs%ncoils
    CALL gs_plasma_mutual(tMaker_obj%gs,tMaker_obj%gs%psi_coil(i)%f,tMaker_obj%gs%Lcoils(i,tMaker_obj%gs%ncoils+1),itor)
    tMaker_obj%gs%Lcoils(tMaker_obj%gs%ncoils+1,i)=tMaker_obj%gs%Lcoils(i,tMaker_obj%gs%ncoils+1)
  END DO
  !
  CALL tMaker_obj%gs%psi%new(rhs)
  CALL tMaker_obj%gs%psi%new(vec1)
  CALL tMaker_obj%gs%psi%new(vec2)
  CALL gs_source(tMaker_obj%gs,tMaker_obj%gs%psi,rhs,vec1,vec2,tmp1,tmp2,tmp3)
  CALL vec1%set(0.d0)
  CALL tMaker_obj%gs%lu_solver%apply(vec1,rhs)
  CALL gs_plasma_mutual(tMaker_obj%gs,vec1,tMaker_obj%gs%Lcoils(tMaker_obj%gs%ncoils+1,tMaker_obj%gs%ncoils+1),itor)
  tMaker_obj%gs%Lcoils(tMaker_obj%gs%ncoils+1,tMaker_obj%gs%ncoils+1)=tMaker_obj%gs%Lcoils(tMaker_obj%gs%ncoils+1,tMaker_obj%gs%ncoils+1)/itor
  CALL rhs%delete()
  CALL vec1%delete()
  CALL vec2%delete()
  DEALLOCATE(rhs,vec1,vec2)
ELSE
  tMaker_obj%gs%Lcoils(tMaker_obj%gs%ncoils+1,:)=0.d0
  tMaker_obj%gs%Lcoils(:,tMaker_obj%gs%ncoils+1)=0.d0
END IF
!---Copy out inductance matrix
CALL c_f_pointer(Lmat, vals_tmp, [tMaker_obj%gs%ncoils+1,tMaker_obj%gs%ncoils+1])
vals_tmp=tMaker_obj%gs%Lcoils
END SUBROUTINE tokamaker_get_coil_Lmat
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_refs(tMaker_ptr,o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm,error_str) BIND(C,NAME="tokamaker_get_refs")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
TYPE(c_ptr), INTENT(out) :: o_point !< Needs docs
TYPE(c_ptr), INTENT(out) :: lim_point !< Needs docs
TYPE(c_ptr), INTENT(out) :: x_points !< Needs docs
TYPE(c_ptr), INTENT(out) :: diverted !< Needs docs
TYPE(c_ptr), INTENT(out) :: plasma_bounds !< Needs docs
TYPE(c_ptr), INTENT(out) :: alam !< Needs docs
TYPE(c_ptr), INTENT(out) :: pnorm !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
o_point=c_loc(tMaker_obj%gs%o_point)
lim_point=c_loc(tMaker_obj%gs%lim_point)
x_points=c_loc(tMaker_obj%gs%x_points)
diverted=c_loc(tMaker_obj%gs%diverted)
plasma_bounds=c_loc(tMaker_obj%gs%plasma_bounds)
alam=c_loc(tMaker_obj%gs%alam)
pnorm=c_loc(tMaker_obj%gs%pnorm)
END SUBROUTINE tokamaker_get_refs
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_trace_surf(tMaker_ptr,psi_surf,points,npoints,error_str) BIND(C,NAME="tokamaker_trace_surf")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), VALUE, INTENT(in) :: psi_surf !< Needs docs
TYPE(c_ptr), INTENT(out) ::  points !< Needs docs
INTEGER(c_int), INTENT(out) :: npoints !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:,:) :: pts_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL gs_trace_surf(tMaker_obj%gs,psi_surf,pts_tmp,npoints)
IF(npoints>0)THEN
  points = c_loc(pts_tmp)
ELSE
  points = C_NULL_PTR
END IF
END SUBROUTINE tokamaker_trace_surf
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_q(tMaker_ptr,npsi,psi_q,qvals,ravgs,dl,rbounds,zbounds,error_str) BIND(C,NAME="tokamaker_get_q")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_q(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: qvals(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: ravgs(npsi,3) !< Needs docs
REAL(c_double), INTENT(inout) :: dl !< Needs docs
REAL(c_double), INTENT(out) :: rbounds(2,2) !< Needs docs
REAL(c_double), INTENT(out) :: zbounds(2,2) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(dl>0.d0)THEN
  CALL gs_get_qprof(tMaker_obj%gs,npsi,psi_q,qvals,dl,rbounds,zbounds,ravgs)
ELSE
  CALL gs_get_qprof(tMaker_obj%gs,npsi,psi_q,qvals,ravgs=ravgs)
  dl = -1.d0
  rbounds = 0.d0
  zbounds = 0.d0
END IF
END SUBROUTINE tokamaker_get_q
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_sauter_fc(tMaker_ptr,npsi,psi_saut,fc,r_avgs,modb_avgs,error_str) BIND(C,NAME="tokamaker_sauter_fc")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_saut(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: fc(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: r_avgs(npsi,3) !< Needs docs
REAL(c_double), INTENT(out) :: modb_avgs(npsi,2) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL sauter_fc(tMaker_obj%gs,npsi,psi_saut,fc,r_avgs,modb_avgs)
END SUBROUTINE tokamaker_sauter_fc
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_globals(tMaker_ptr,Itor,centroid,vol,pvol,dflux,tflux,bp_vol,error_str) BIND(C,NAME="tokamaker_get_globals")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), INTENT(out) :: Itor !< Needs docs
REAL(c_double), INTENT(out) :: centroid(2) !< Needs docs
REAL(c_double), INTENT(out) :: vol !< Needs docs
REAL(c_double), INTENT(out) :: pvol !< Needs docs
REAL(c_double), INTENT(out) :: dflux !< Needs docs
REAL(c_double), INTENT(out) :: tflux !< Needs docs
REAL(c_double), INTENT(out) :: bp_vol !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL gs_comp_globals(tMaker_obj%gs,Itor,centroid,vol,pvol,dflux,tflux,bp_vol)
Itor=Itor/mu0
vol=vol*2.d0*pi
pvol=pvol*2.d0*pi/mu0
END SUBROUTINE tokamaker_get_globals
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_gs_calc_vloop(tMaker_ptr,vloop,error_str) BIND(C,NAME="tokamaker_gs_calc_vloop")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
REAL(c_double), INTENT(out) :: vloop
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(.NOT.ASSOCIATED(tMaker_obj%gs%eta))THEN
  vloop=-1.d0
  RETURN
END IF
CALL gs_calc_vloop(tMaker_obj%gs,vloop)
END SUBROUTINE tokamaker_gs_calc_vloop
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_profs(tMaker_ptr,npsi,psi_in,f,fp,p,pp,error_str) BIND(C,NAME="tokamaker_get_profs")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_in(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: f(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: fp(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: p(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: pp(npsi) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i
REAL(8) :: x1,x2,r
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
x1=0.d0; x2=1.d0
IF(tMaker_obj%gs%plasma_bounds(1)>-1.d98)THEN
  x1=tMaker_obj%gs%plasma_bounds(1); x2=tMaker_obj%gs%plasma_bounds(2)
END IF
DO i=1,npsi
  r=psi_in(i)*(x2-x1) + x1
  IF(tMaker_obj%gs%mode==0)THEN
    fp(i)=tMaker_obj%gs%alam*tMaker_obj%gs%I%fp(r)
    f(i)=tMaker_obj%gs%psiscale*tMaker_obj%gs%alam*tMaker_obj%gs%I%f(r) + tMaker_obj%gs%I%f_offset
  ELSE
    f(i)=SQRT(tMaker_obj%gs%psiscale*tMaker_obj%gs%alam*tMaker_obj%gs%I%f(r) + tMaker_obj%gs%I%f_offset**2)
    fp(i)=tMaker_obj%gs%alam*tMaker_obj%gs%I%fp(r)/(2.d0*f(i))
  END IF
  pp(i)=tMaker_obj%gs%psiscale*tMaker_obj%gs%pnorm*tMaker_obj%gs%P%fp(r)
  p(i)=tMaker_obj%gs%psiscale*tMaker_obj%gs%psiscale*tMaker_obj%gs%pnorm*tMaker_obj%gs%P%f(r)
END DO
END SUBROUTINE tokamaker_get_profs
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_vfixed(tMaker_ptr,npts,pts,fluxes,error_str) BIND(C,NAME="tokamaker_get_vfixed")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(c_int), INTENT(out) :: npts !< Needs docs
TYPE(c_ptr), INTENT(out) :: pts !< Needs docs
TYPE(c_ptr), INTENT(out) :: fluxes !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER :: pts_tmp(:,:),fluxes_tmp(:)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL gs_fixed_vflux(tMaker_obj%gs,pts_tmp,fluxes_tmp)
pts=C_LOC(pts_tmp)
fluxes=C_LOC(fluxes_tmp)
npts=SIZE(fluxes_tmp,DIM=1)
END SUBROUTINE tokamaker_get_vfixed
!---------------------------------------------------------------------------------
!> Create an interpolation object for tokamaker fields
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_field_eval(tMaker_ptr,imode,int_obj,error_str) BIND(C,NAME="tokamaker_get_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< Pointer to TokaMaker object
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Field type
TYPE(c_ptr), INTENT(out) :: int_obj !< Pointer to interpolation object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(oft_lag_bginterp), POINTER :: psi_grad_obj
TYPE(gs_prof_interp), POINTER :: prof_interp_obj
TYPE(gs_b_interp), POINTER :: b_interp_obj
TYPE(tokamaker_instance), POINTER :: tMaker_obj
CLASS(oft_vector), POINTER :: tmp1,tmp2,tmp3
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(imode==1)THEN
  ALLOCATE(b_interp_obj)
  b_interp_obj%gs=>tMaker_obj%gs
  CALL b_interp_obj%setup(tMaker_obj%gs)
  int_obj=C_LOC(b_interp_obj)
ELSEIF(imode>=2.AND.imode<=4)THEN
  ALLOCATE(prof_interp_obj)
  prof_interp_obj%gs=>tMaker_obj%gs
  prof_interp_obj%mode=imode-1
  CALL prof_interp_obj%setup(tMaker_obj%gs)
  int_obj=C_LOC(prof_interp_obj)
ELSEIF(imode==5)THEN
  ALLOCATE(psi_grad_obj)
  psi_grad_obj%u=>tMaker_obj%gs%psi
  CALL psi_grad_obj%setup(tMaker_obj%gs%fe_rep)
  int_obj=C_LOC(psi_grad_obj)
ELSEIF(imode>=6.AND.imode<=8)THEN
  ALLOCATE(psi_grad_obj)
  CALL tMaker_obj%gs%psi%new(tmp1)
  CALL tMaker_obj%gs%psi%new(tmp2)
  CALL tMaker_obj%gs%psi%new(tmp3)
  CALL gs_project_b(tMaker_obj%gs,tmp1,tmp2,tmp3)
  IF(imode==6)THEN
    psi_grad_obj%u=>tmp1
    NULLIFY(tmp1)
  ELSE IF(imode==7)THEN
    psi_grad_obj%u=>tmp2
    NULLIFY(tmp2)
  ELSE IF(imode==8)THEN
    psi_grad_obj%u=>tmp3
    NULLIFY(tmp3)
  END IF
  IF(ASSOCIATED(tmp1))THEN
    CALL tmp1%delete()
    DEALLOCATE(tmp1)
  END IF
  IF(ASSOCIATED(tmp2))THEN
    CALL tmp2%delete()
    DEALLOCATE(tmp2)
  END IF
  IF(ASSOCIATED(tmp3))THEN
    CALL tmp3%delete()
    DEALLOCATE(tmp3)
  END IF
  CALL psi_grad_obj%setup(tMaker_obj%gs%fe_rep)
  int_obj=C_LOC(psi_grad_obj)
END IF
END SUBROUTINE tokamaker_get_field_eval
!---------------------------------------------------------------------------------
!> Evaluate a TokaMaker field with an interpolation object created by
!! \ref tokamaker_f::tokamaker_get_field_eval
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_apply_field_eval(tMaker_ptr,int_obj,int_type,pt,fbary_tol,cell,dim,field) BIND(C,NAME="tokamaker_apply_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Pointer to interpolation object
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Field type (negative to destroy)
REAL(c_double), INTENT(in) :: pt(3) !< Location for evaluation [R,Z,0]
REAL(c_double), VALUE, INTENT(in) :: fbary_tol !< Tolerance for physical to logical mapping
INTEGER(c_int), INTENT(inout) :: cell !< Cell containing `pt` (starting guess on input)
INTEGER(c_int), VALUE, INTENT(in) :: dim !< Dimension of field
REAL(c_double), INTENT(out) :: field(dim) !< Field at `pt`
TYPE(oft_lag_bginterp), POINTER :: psi_grad_obj
TYPE(gs_prof_interp), POINTER :: prof_interp_obj
TYPE(gs_b_interp), POINTER :: b_interp_obj
REAL(8) :: f(4),goptmp(3,4),vol,fmin,fmax
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj))CALL oft_abort("TokaMaker object not associated","tokamaker_apply_field_eval",__FILE__)
IF(int_type<0)THEN
  IF(ABS(int_type)==1)THEN
    CALL c_f_pointer(int_obj, b_interp_obj)
    CALL b_interp_obj%delete
  ELSE IF(ABS(int_type)>=2.AND.ABS(int_type)<=4)THEN
    CALL c_f_pointer(int_obj, b_interp_obj)
    CALL b_interp_obj%delete
  ELSE IF(ABS(int_type)>=5)THEN
    CALL c_f_pointer(int_obj, psi_grad_obj)
    IF(ABS(int_type)>=6)THEN
      CALL psi_grad_obj%u%delete
      DEALLOCATE(psi_grad_obj%u)
    END IF
    CALL psi_grad_obj%delete
  END IF
  RETURN
END IF
call bmesh_findcell(tMaker_obj%gs%mesh,cell,pt,f)
IF(cell==0)RETURN
fmin=MINVAL(f); fmax=MAXVAL(f)
IF(( fmax>1.d0+fbary_tol ).OR.( fmin<-fbary_tol ))THEN
  cell=-ABS(cell)
  RETURN
END IF
CALL tMaker_obj%gs%mesh%jacobian(cell,f,goptmp,vol)
IF(int_type==1)THEN
  CALL c_f_pointer(int_obj, b_interp_obj)
  CALL b_interp_obj%interp(cell,f,goptmp,field)
ELSE IF(int_type>=2.AND.int_type<=4)THEN
  CALL c_f_pointer(int_obj, prof_interp_obj)
  CALL prof_interp_obj%interp(cell,f,goptmp,field)
ELSE IF(int_type>=5)THEN
  CALL c_f_pointer(int_obj, psi_grad_obj)
  CALL psi_grad_obj%interp(cell,f,goptmp,field)
END IF
END SUBROUTINE tokamaker_apply_field_eval
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_psi(tMaker_ptr,psi_vals,update_bounds,error_str) BIND(C,NAME="tokamaker_set_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: update_bounds !< Update bounds by determining new limiting points
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(psi_vals, vals_tmp, [tMaker_obj%gs%psi%n])
CALL tMaker_obj%gs%psi%restore_local(vals_tmp)
IF(update_bounds)CALL gs_update_bounds(tMaker_obj%gs)
END SUBROUTINE tokamaker_set_psi
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_psi_dt(tMaker_ptr,psi_vals,dt,error_str) BIND(C,NAME="tokamaker_set_psi_dt")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: dt !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
tMaker_obj%gs%dt=dt
IF(dt>0.d0)THEN
  IF(.NOT.ASSOCIATED(tMaker_obj%gs%psi_dt))CALL tMaker_obj%gs%psi%new(tMaker_obj%gs%psi_dt)
  CALL c_f_pointer(psi_vals, vals_tmp, [tMaker_obj%gs%psi%n])
  CALL tMaker_obj%gs%psi_dt%restore_local(vals_tmp)
ELSE
  IF(ASSOCIATED(tMaker_obj%gs%psi_dt))THEN
    CALL tMaker_obj%gs%psi_dt%delete()
    DEALLOCATE(tMaker_obj%gs%psi_dt)
  END IF
END IF
END SUBROUTINE tokamaker_set_psi_dt
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_settings(tMaker_ptr,settings,error_str) BIND(C,NAME="tokamaker_set_settings")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(tokamaker_settings_type), INTENT(in) :: settings !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
CHARACTER(KIND=c_char), POINTER, DIMENSION(:) :: limfile_c
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
oft_env%pm=settings%pm
tMaker_obj%gs%has_plasma=settings%has_plasma
tMaker_obj%gs%free=settings%free_boundary
tMaker_obj%gs%lim_zmax=settings%lim_zmax
tMaker_obj%gs%rmin=settings%rmin
tMaker_obj%gs%mode=settings%mode
tMaker_obj%gs%urf=settings%urf
tMaker_obj%gs%maxits=settings%maxits
tMaker_obj%gs%nl_tol=settings%nl_tol
tMaker_obj%gs%dipole_mode=settings%dipole_mode
IF(tMaker_obj%gs%dipole_mode)CALL oft_warn("TokaMaker's dipole functionality is experimental, use with caution")
CALL c_f_pointer(settings%limiter_file,limfile_c,[OFT_PATH_SLEN])
CALL copy_string_rev(limfile_c,tMaker_obj%gs%limiter_file)
END SUBROUTINE tokamaker_set_settings
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_dipole_a(tMaker_ptr,dipole_a,error_str) BIND(C,NAME="tokamaker_set_dipole_a")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
REAL(c_double), VALUE, INTENT(in) :: dipole_a !< New value for dipole_a
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
tMaker_obj%gs%dipole_a=dipole_a
END SUBROUTINE tokamaker_set_dipole_a
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_targets(tMaker_ptr,ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target,error_str) BIND(C,NAME="tokamaker_set_targets")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
REAL(c_double), VALUE, INTENT(in) :: ip_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: ip_ratio_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: pax_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: estore_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: R0_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: V0_target !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
tMaker_obj%gs%R0_target=R0_target
tMaker_obj%gs%V0_target=V0_target
tMaker_obj%gs%pax_target=pax_target*mu0
tMaker_obj%gs%estore_target=estore_target*mu0
tMaker_obj%gs%itor_target=ip_target*mu0
tMaker_obj%gs%ip_ratio_target=ip_ratio_target
END SUBROUTINE tokamaker_set_targets
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_isoflux(tMaker_ptr,targets,ref_points,weights,ntargets,grad_wt_lim,error_str) BIND(C,NAME="tokamaker_set_isoflux")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
REAL(c_double), INTENT(in) :: targets(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: ref_points(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: grad_wt_lim !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(ASSOCIATED(tMaker_obj%gs%isoflux_targets))DEALLOCATE(tMaker_obj%gs%isoflux_targets)
tMaker_obj%gs%isoflux_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(tMaker_obj%gs%isoflux_targets(5,tMaker_obj%gs%isoflux_ntargets))
  tMaker_obj%gs%isoflux_targets(1:2,:)=targets
  tMaker_obj%gs%isoflux_targets(3,:)=weights
  tMaker_obj%gs%isoflux_targets(4:5,:)=ref_points
  tMaker_obj%gs%isoflux_grad_wt_lim=1.d0/grad_wt_lim
ELSE
  tMaker_obj%gs%isoflux_grad_wt_lim=-1.d0
END IF
END SUBROUTINE tokamaker_set_isoflux
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_flux(tMaker_ptr,locations,targets,weights,ntargets,grad_wt_lim,error_str) BIND(C,NAME="tokamaker_set_flux")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
REAL(c_double), INTENT(in) :: locations(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: targets(ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: grad_wt_lim !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(ASSOCIATED(tMaker_obj%gs%flux_targets))DEALLOCATE(tMaker_obj%gs%flux_targets)
tMaker_obj%gs%flux_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(tMaker_obj%gs%flux_targets(4,tMaker_obj%gs%flux_ntargets))
  tMaker_obj%gs%flux_targets(1:2,:)=locations
  tMaker_obj%gs%flux_targets(3,:)=targets
  tMaker_obj%gs%flux_targets(4,:)=weights
  ! tMaker_obj%gs%isoflux_grad_wt_lim=1.d0/grad_wt_lim
! ELSE
  ! tMaker_obj%gs%isoflux_grad_wt_lim=-1.d0
END IF
END SUBROUTINE tokamaker_set_flux
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_saddles(tMaker_ptr,targets,weights,ntargets,error_str) BIND(C,NAME="tokamaker_set_saddles")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
REAL(c_double), INTENT(in) :: targets(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(ASSOCIATED(tMaker_obj%gs%saddle_targets))DEALLOCATE(tMaker_obj%gs%saddle_targets)
tMaker_obj%gs%saddle_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(tMaker_obj%gs%saddle_targets(3,tMaker_obj%gs%saddle_ntargets))
  tMaker_obj%gs%saddle_targets(1:2,:)=targets
  tMaker_obj%gs%saddle_targets(3,:)=weights
END IF
END SUBROUTINE tokamaker_set_saddles
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_currents(tMaker_ptr,currents,error_str) BIND(C,NAME="tokamaker_set_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: currents !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
INTEGER(4) :: i
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(currents, vals_tmp, [tMaker_obj%gs%ncoils])
tMaker_obj%gs%coil_currs = vals_tmp*mu0
tMaker_obj%gs%vcontrol_val = 0.d0
END SUBROUTINE tokamaker_set_coil_currents
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_regmat(tMaker_ptr,nregularize,coil_reg_mat,coil_reg_targets,coil_reg_weights,error_str) BIND(C,NAME="tokamaker_set_coil_regmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
INTEGER(c_int), VALUE, INTENT(in) :: nregularize !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_mat !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_targets !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_weights !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
IF(ASSOCIATED(tMaker_obj%gs%coil_reg_mat))DEALLOCATE(tMaker_obj%gs%coil_reg_mat,tMaker_obj%gs%coil_reg_targets)
tMaker_obj%gs%nregularize=nregularize
ALLOCATE(tMaker_obj%gs%coil_reg_mat(tMaker_obj%gs%nregularize,tMaker_obj%gs%ncoils+1))
ALLOCATE(tMaker_obj%gs%coil_reg_targets(tMaker_obj%gs%nregularize))
CALL c_f_pointer(coil_reg_mat, vals_tmp, [tMaker_obj%gs%nregularize,tMaker_obj%gs%ncoils+1])
tMaker_obj%gs%coil_reg_mat=vals_tmp
CALL c_f_pointer(coil_reg_targets, vals_tmp, [tMaker_obj%gs%nregularize,1])
tMaker_obj%gs%coil_reg_targets=vals_tmp(:,1)*mu0
CALL c_f_pointer(coil_reg_weights, vals_tmp, [tMaker_obj%gs%nregularize,1])
DO i=1,tMaker_obj%gs%nregularize
  tMaker_obj%gs%coil_reg_targets(i)=tMaker_obj%gs%coil_reg_targets(i)*vals_tmp(i,1)
  tMaker_obj%gs%coil_reg_mat(i,:)=tMaker_obj%gs%coil_reg_mat(i,:)*vals_tmp(i,1)
END DO
END SUBROUTINE tokamaker_set_coil_regmat
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_bounds(tMaker_ptr,coil_bounds,error_str) BIND(C,NAME="tokamaker_set_coil_bounds")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: coil_bounds !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(coil_bounds, vals_tmp, [2,tMaker_obj%gs%ncoils+1])
IF(.NOT.ASSOCIATED(tMaker_obj%gs%coil_bounds))THEN
  ALLOCATE(tMaker_obj%gs%coil_bounds(2,tMaker_obj%gs%ncoils+1))
  tMaker_obj%gs%coil_bounds(1,:)=-1.d98; tMaker_obj%gs%coil_bounds(2,:)=1.d98
END IF
DO i=1,tMaker_obj%gs%ncoils
  tMaker_obj%gs%coil_bounds([2,1],i)=-vals_tmp(:,i)*mu0
END DO
tMaker_obj%gs%coil_bounds([2,1],tMaker_obj%gs%ncoils+1)=-vals_tmp(:,tMaker_obj%gs%ncoils+1)*mu0
END SUBROUTINE tokamaker_set_coil_bounds
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_vsc(tMaker_ptr,coil_gains,error_str) BIND(C,NAME="tokamaker_set_coil_vsc")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
TYPE(c_ptr), VALUE, INTENT(in) :: coil_gains !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
INTEGER(4) :: i
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(coil_gains, vals_tmp, [tMaker_obj%gs%ncoils])
tMaker_obj%gs%coil_vcont=vals_tmp
END SUBROUTINE tokamaker_set_coil_vsc
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE tokamaker_save_eqdsk(tMaker_ptr,filename,nr,nz,rbounds,zbounds,run_info,psi_pad,rcentr,trunc_eq,lim_filename,lcfs_press,error_str) BIND(C,NAME="tokamaker_save_eqdsk")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: run_info(40) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nr !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nz !< Needs docs
REAL(c_double), INTENT(in) :: rbounds(2) !< Needs docs
REAL(c_double), INTENT(in) :: zbounds(2) !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: psi_pad !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: rcentr !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: trunc_eq !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: lim_filename(OFT_PATH_SLEN) !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: lcfs_press !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
CHARACTER(LEN=40) :: run_info_f
CHARACTER(LEN=OFT_PATH_SLEN) :: filename_tmp,lim_file
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_flag
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL copy_string_rev(run_info,run_info_f)
CALL copy_string_rev(filename,filename_tmp)
CALL copy_string_rev(lim_filename,lim_file)
IF(rcentr>0.d0)THEN
  CALL gs_save_eqdsk(tMaker_obj%gs,filename_tmp,nr,nz,rbounds,zbounds,run_info_f,lim_file,psi_pad, &
    rcentr_in=rcentr,trunc_eq=LOGICAL(trunc_eq),lcfs_press=lcfs_press,error_str=error_flag)
ELSE
  CALL gs_save_eqdsk(tMaker_obj%gs,filename_tmp,nr,nz,rbounds,zbounds,run_info_f,lim_file,psi_pad, &
    trunc_eq=LOGICAL(trunc_eq),lcfs_press=lcfs_press,error_str=error_flag)
END IF
CALL copy_string(TRIM(error_flag),error_str)
END SUBROUTINE tokamaker_save_eqdsk
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_save_ifile(tMaker_ptr,filename,npsi,ntheta,psi_pad,lcfs_press,pack_lcfs,single_prec,error_str) BIND(C,NAME="tokamaker_save_ifile")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntheta !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: psi_pad !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: lcfs_press !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pack_lcfs !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: single_prec !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
CHARACTER(LEN=OFT_PATH_SLEN) :: filename_tmp
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_flag
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL copy_string_rev(filename,filename_tmp)
CALL gs_save_ifile(tMaker_obj%gs,filename_tmp,npsi,ntheta,psi_pad,lcfs_press=lcfs_press, &
  pack_lcfs=LOGICAL(pack_lcfs),single_prec=LOGICAL(single_prec),error_str=error_flag)
CALL copy_string(TRIM(error_flag),error_str)
END SUBROUTINE tokamaker_save_ifile
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_save_mug(tMaker_ptr,filename,error_str) BIND(C,NAME="tokamaker_save_mug")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
CHARACTER(LEN=OFT_PATH_SLEN) :: filename_tmp
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_flag
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL copy_string_rev(filename,filename_tmp)
CALL gs_save_mug(tMaker_obj%gs,filename_tmp)
CALL copy_string(TRIM(error_flag),error_str)
END SUBROUTINE tokamaker_save_mug
!---------------------------------------------------------------------------
!> Overwrites default coil flux contribution to non-uniform current distribution
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_current_dist(tMaker_ptr,iCoil,curr_dist,error_str) BIND(C,NAME="tokamaker_set_coil_current_dist")
TYPE(c_ptr), VALUE, INTENT(in) :: tMaker_ptr !< TokaMaker instance
INTEGER(c_int), VALUE, INTENT(in) :: iCoil
TYPE(c_ptr), VALUE, INTENT(in) :: curr_dist
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
INTEGER(4) :: i
class(oft_vector), pointer :: tmp_vec
TYPE(tokamaker_instance), POINTER :: tMaker_obj
IF(.NOT.tokamaker_ccast(tMaker_ptr,tMaker_obj,error_str))RETURN
CALL c_f_pointer(curr_dist, vals_tmp, [tMaker_obj%gs%psi%n])
! Update coil flux to overwrite old uniform distribution
NULLIFY(tmp_vec)
call tMaker_obj%gs%psi%new(tmp_vec)

CALL gs_coil_source_distributed(tMaker_obj%gs,iCoil,tmp_vec,vals_tmp)

CALL tMaker_obj%gs%zerob_bc%apply(tmp_vec)
CALL gs_vacuum_solve(tMaker_obj%gs,tMaker_obj%gs%psi_coil(iCoil)%f,tmp_vec)
! Update coil mutual inductances
DO i=1,tMaker_obj%gs%ncoils
  CALL gs_coil_mutual(tMaker_obj%gs,i,tMaker_obj%gs%psi_coil(iCoil)%f,tMaker_obj%gs%Lcoils(i,iCoil))
  tMaker_obj%gs%Lcoils(iCoil,i)=tMaker_obj%gs%Lcoils(i,iCoil)
  IF(i==iCoil)THEN
    CALL gs_coil_mutual_distributed(tMaker_obj%gs,i,tMaker_obj%gs%psi_coil(iCoil)%f,vals_tmp,tMaker_obj%gs%Lcoils(i,iCoil))
  ELSE
    CALL gs_coil_mutual(tMaker_obj%gs,i,tMaker_obj%gs%psi_coil(iCoil)%f,tMaker_obj%gs%Lcoils(i,iCoil))
    tMaker_obj%gs%Lcoils(iCoil,i)=tMaker_obj%gs%Lcoils(i,iCoil)
  END IF
END DO
END SUBROUTINE tokamaker_set_coil_current_dist
END MODULE tokamaker_f

