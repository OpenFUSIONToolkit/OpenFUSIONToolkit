!---------------------------------------------------------------------------
!> @file tokamaker_f.F90
!
!> Fortran part of Python wrapper for Grad-Shafranov functionality
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------
MODULE tokamaker_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
  c_f_pointer, c_bool, c_null_ptr
USE oft_base
USE oft_mesh_type, ONLY: smesh, bmesh_findcell
! USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid, ONLY: multigrid_reset
! USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_base, ONLY: oft_afem_type
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, oft_blagrange, &
  oft_lag_setup
USE mhd_utils, ONLY: mu0
USE axi_green, ONLY: green
USE oft_gs, ONLY: gs_eq, gs_save_fields, gs_save_fgrid, gs_setup_walls, build_dels, &
  gs_fixed_vflux, gs_load_regions, gs_get_qprof, gs_trace_surf, gs_b_interp, gs_prof_interp, &
  gs_plasma_mutual, gs_source
USE oft_gs_util, ONLY: gs_save, gs_load, gs_analyze, gs_comp_globals, gs_save_eqdsk, &
  gs_profile_load, sauter_fc, gs_calc_vloop
USE oft_gs_fit, ONLY: fit_gs, fit_pm
USE oft_gs_td, ONLY: oft_tmaker_td, eig_gs_td
USE oft_base_f, ONLY: copy_string, copy_string_rev, oftpy_init
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE, BIND(C) :: tokamaker_settings_type
  LOGICAL(KIND=c_bool) :: pm = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: free_boundary = .FALSE. !< Needs docs
  LOGICAL(KIND=c_bool) :: has_plasma = .TRUE. !< Needs docs
  LOGICAL(KIND=c_bool) :: limited_only = .FALSE. !< Needs docs
  INTEGER(KIND=c_int) :: maxits = 40 !< Needs docs
  INTEGER(KIND=c_int) :: mode = 1 !< Needs docs
  REAL(KIND=c_double) :: urf = 0.3d0 !< Needs docs
  REAL(KIND=c_double) :: nl_tol = 1.d-6 !< Needs docs
  REAL(KIND=c_double) :: rmin = 0.d0 !< Needs docs
  REAL(KIND=c_double) :: lim_zmax = 1.d99 !< Needs docs
  CHARACTER(KIND=c_char) :: limiter_file(80) = 'none' !< Needs docs
END TYPE tokamaker_settings_type
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
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
END TYPE tokamaker_recon_settings_type
!
TYPE(gs_eq), POINTER :: gs_global => NULL() !< Global G-S object
TYPE(oft_tmaker_td), POINTER :: gs_td_global => NULL() !< Global time-dependent object
integer(i4), POINTER :: lc_plot(:,:) => NULL() !< Needs docs
integer(i4), POINTER :: reg_plot(:) => NULL() !< Needs docs
real(r8), POINTER :: r_plot(:,:) => NULL() !< Needs docs
CONTAINS
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_alloc(gs_obj) BIND(C,NAME="tokamaker_alloc")
TYPE(c_ptr), INTENT(out) :: gs_obj !< Needs docs
IF(.NOT.ASSOCIATED(gs_global))ALLOCATE(gs_global)
gs_obj=c_loc(gs_global) ! Do nothing with this for now
END SUBROUTINE tokamaker_alloc
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup_regions(coil_file,reg_eta,contig_flag,xpoint_mask,coil_nturns,ncoils) BIND(C,NAME="tokamaker_setup_regions")
CHARACTER(KIND=c_char), INTENT(in) :: coil_file(80) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_eta !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: contig_flag !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: xpoint_mask !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_nturns !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ncoils !< Needs docs
real(r8), POINTER :: eta_tmp(:),nturns_tmp(:,:)
INTEGER(i4), POINTER :: contig_tmp(:)
INTEGER(4) :: i
INTEGER(4), POINTER :: xpoint_tmp(:)
CALL copy_string_rev(coil_file,gs_global%coil_file)
IF(TRIM(gs_global%coil_file)=='none')THEN
  !
  CALL c_f_pointer(xpoint_mask, xpoint_tmp, [smesh%nreg])
  ALLOCATE(gs_global%saddle_rmask(smesh%nreg))
  gs_global%saddle_rmask=LOGICAL(xpoint_tmp==0)
  !
  CALL c_f_pointer(reg_eta, eta_tmp, [smesh%nreg])
  gs_global%ncoil_regs=0
  gs_global%ncond_regs=0
  DO i=2,smesh%nreg
    IF(eta_tmp(i)>0.d0)THEN
      gs_global%ncond_regs=gs_global%ncond_regs+1
    ELSE
      gs_global%ncoil_regs=gs_global%ncoil_regs+1
    END IF
  END DO
  !
  gs_global%ncoils=ncoils
  ALLOCATE(gs_global%coil_vcont(ncoils),gs_global%coil_currs(ncoils))
  gs_global%coil_vcont = 0.d0
  gs_global%coil_currs = 0.d0
  CALL c_f_pointer(coil_nturns, nturns_tmp, [smesh%nreg,ncoils])
  ALLOCATE(gs_global%coil_nturns(smesh%nreg,ncoils))
  gs_global%coil_nturns=0.d0
  gs_global%coil_nturns=nturns_tmp
  !
  CALL c_f_pointer(contig_flag, contig_tmp, [smesh%nreg])
  ALLOCATE(gs_global%cond_regions(gs_global%ncond_regs))
  ALLOCATE(gs_global%coil_regions(gs_global%ncoil_regs))
  gs_global%ncond_regs=0
  gs_global%ncoil_regs=0
  DO i=2,smesh%nreg
    IF(eta_tmp(i)>0.d0)THEN
      gs_global%ncond_regs=gs_global%ncond_regs+1
      IF(eta_tmp(i)<1.d9)THEN
        gs_global%cond_regions(gs_global%ncond_regs)%eta=eta_tmp(i)
      END IF
      gs_global%cond_regions(gs_global%ncond_regs)%id=i
      gs_global%cond_regions(gs_global%ncond_regs)%continuous=(contig_tmp(i)==1)
    ELSE
      gs_global%ncoil_regs=gs_global%ncoil_regs+1
      gs_global%coil_regions(gs_global%ncoil_regs)%id=i
    END IF
  END DO
ELSE
  CALL gs_load_regions(gs_global)
END IF
END SUBROUTINE tokamaker_setup_regions
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_reset(error_str) BIND(C,NAME="tokamaker_reset")
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
INTEGER(4) :: i,ierr,io_unit,npts,iostat
REAL(8) :: theta
LOGICAL :: file_exists
real(r8), POINTER :: vals_tmp(:)
CHARACTER(LEN=80) :: tmp_str
!---Clear error flag
CALL copy_string('',error_str)
!---
IF(ASSOCIATED(r_plot))DEALLOCATE(r_plot)
IF(ASSOCIATED(lc_plot))DEALLOCATE(lc_plot)
IF(ASSOCIATED(reg_plot))DEALLOCATE(reg_plot)
!---Destroy objects
CALL gs_global%delete()
DEALLOCATE(gs_global)
CALL oft_blagrange%delete()
CALL multigrid_reset
ALLOCATE(gs_global)
END SUBROUTINE tokamaker_reset
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup(order,full_domain,ncoils,error_str) BIND(C,NAME="tokamaker_setup")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: full_domain !< Needs docs
INTEGER(KIND=c_int), INTENT(out) :: ncoils !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
INTEGER(4) :: i,ierr,io_unit,npts,iostat
REAL(8) :: theta
LOGICAL :: file_exists
real(r8), POINTER :: vals_tmp(:)
CHARACTER(LEN=80) :: tmp_str
!---Clear error flag
CALL copy_string('',error_str)
!---------------------------------------------------------------------------
! Check input files
!---------------------------------------------------------------------------
IF(TRIM(gs_global%coil_file)/='none')THEN
    INQUIRE(EXIST=file_exists,FILE=TRIM(gs_global%coil_file))
    IF(.NOT.file_exists)THEN
        CALL copy_string('Specified "coil_file" cannot be found',error_str)
        RETURN
    END IF
END IF
IF(TRIM(gs_global%limiter_file)/='none')THEN
    INQUIRE(EXIST=file_exists,FILE=TRIM(gs_global%limiter_file))
    IF(.NOT.file_exists)THEN
        CALL copy_string('Specified "limiter_file" cannot be found',error_str)
        RETURN
    END IF
END IF
! IF(TRIM(gs_global%eqdsk_limiter_file)/='none')THEN
!     INQUIRE(EXIST=file_exists,FILE=TRIM(gs_global%eqdsk_limiter_file))
!     IF(.NOT.file_exists)THEN
!         CALL copy_string('Specified "eqdsk_limiter_file" cannot be found',error_str)
!         RETURN
!     END IF
! END IF
!---------------------------------------------------------------------------
! Setup Lagrange Elements
!---------------------------------------------------------------------------
! CALL smesh%setup_io(order)
smesh%tess_order=order
CALL oft_lag_setup(order, -1)
!---------------------------------------------------------------------------
! Setup experimental geometry
!---------------------------------------------------------------------------
gs_global%save_visit=.FALSE.
! CALL gs_load_regions(gs_global)
! IF(gs_global%free)THEN
!     CALL gs_global%load_coils
! ELSE
!     CALL gs_load_regions(gs_global)
! END IF
gs_global%full_domain=full_domain
CALL gs_setup_walls(gs_global)
CALL gs_global%load_limiters
CALL gs_global%init()
ncoils=gs_global%ncoils
END SUBROUTINE tokamaker_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_load_profiles(f_file,f_offset,p_file,eta_file,f_NI_file) BIND(C,NAME="tokamaker_load_profiles")
CHARACTER(KIND=c_char), INTENT(in) :: f_file(80) !< F*F' prof.in file
CHARACTER(KIND=c_char), INTENT(in) :: p_file(80) !< P' prof.in file
CHARACTER(KIND=c_char), INTENT(in) :: eta_file(80) !< Resistivity (eta) prof.in file
CHARACTER(KIND=c_char), INTENT(in) :: f_NI_file(80) !< Non-inductive F*F' prof.in file
REAL(c_double), VALUE, INTENT(in) :: f_offset !< Needs docs
CHARACTER(LEN=80) :: tmp_str
CALL copy_string_rev(f_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,gs_global%I)
IF(f_offset>-1.d98)gs_global%I%f_offset=f_offset
CALL copy_string_rev(p_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,gs_global%P)
CALL copy_string_rev(eta_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,gs_global%eta)
CALL copy_string_rev(f_NI_file,tmp_str)
IF(TRIM(tmp_str)/='none')CALL gs_profile_load(tmp_str,gs_global%I_NI)
END SUBROUTINE tokamaker_load_profiles
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_init_psi(r0,z0,a,kappa,delta,ierr) BIND(C,NAME="tokamaker_init_psi")
REAL(c_double), VALUE, INTENT(in) :: r0 !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: z0 !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: a !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: kappa !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: delta !< Needs docs
INTEGER(c_int), INTENT(out) :: ierr !< Needs docs
CALL gs_global%init_psi(ierr,r0=[r0,z0],a=a,kappa=kappa,delta=delta)
END SUBROUTINE tokamaker_init_psi
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_solve(error_flag) BIND(C,NAME="tokamaker_solve")
INTEGER(c_int), INTENT(out) :: error_flag !< Needs docs
CALL gs_global%solve(error_flag)
END SUBROUTINE tokamaker_solve
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_vac_solve(psi_in,error_flag) BIND(C,NAME="tokamaker_vac_solve")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_in !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
INTEGER(c_int), INTENT(out) :: error_flag !< Needs docs
CLASS(oft_vector), POINTER :: psi_tmp
NULLIFY(psi_tmp)
CALL gs_global%psi%new(psi_tmp)
CALL c_f_pointer(psi_in, vals_tmp, [gs_global%psi%n])
CALL psi_tmp%restore_local(vals_tmp)
CALL gs_global%vac_solve(psi_tmp,error_flag)
CALL psi_tmp%get_local(vals_tmp)
CALL psi_tmp%delete()
DEALLOCATE(psi_tmp)
END SUBROUTINE tokamaker_vac_solve
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_analyze() BIND(C,NAME="tokamaker_analyze")
CALL gs_analyze(gs_global)
END SUBROUTINE tokamaker_analyze
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_recon_run(vacuum,settings,error_flag) BIND(C,NAME="tokamaker_recon_run")
LOGICAL(c_bool), VALUE, INTENT(in) :: vacuum !< Needs docs
TYPE(tokamaker_recon_settings_type), INTENT(in) :: settings !< Needs docs
INTEGER(c_int), INTENT(out) :: error_flag !< Needs docs
LOGICAL :: fitI,fitP,fitPnorm,fitAlam,fitR0,fitV0,fitCoils,fitF0,fixedCentering
IF(vacuum)gs_global%has_plasma=.FALSE.
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
CALL fit_gs(gs_global,fitI,fitP,fitPnorm,&
            fitAlam,fitR0,fitV0,fitCoils,fitF0, &
            fixedCentering)
gs_global%has_plasma=.TRUE.
END SUBROUTINE tokamaker_recon_run
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_setup_td(dt,lin_tol,nl_tol,pre_plasma) BIND(C,NAME="tokamaker_setup_td")
REAL(c_double), VALUE, INTENT(in) :: dt !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: lin_tol !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: nl_tol !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pre_plasma !< Needs docs
IF(ASSOCIATED(gs_td_global))THEN
  CALL gs_td_global%delete()
  DEALLOCATE(gs_td_global)
END IF
ALLOCATE(gs_td_global)
CALL gs_td_global%setup(gs_global,dt,lin_tol,nl_tol,LOGICAL(pre_plasma))
END SUBROUTINE tokamaker_setup_td
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_eig_td(omega,neigs,eigs,eig_vecs,include_bounds,pm) BIND(C,NAME="tokamaker_eig_td")
REAL(c_double), VALUE, INTENT(in) :: omega !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vecs !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: include_bounds !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pm !< Needs docs
REAL(8), POINTER :: eigs_tmp(:,:),eig_vecs_tmp(:,:)
LOGICAL :: pm_save
CALL c_f_pointer(eigs, eigs_tmp, [2,neigs])
CALL c_f_pointer(eig_vecs, eig_vecs_tmp, [gs_global%psi%n,neigs])
pm_save=oft_env%pm; oft_env%pm=pm
CALL eig_gs_td(gs_global,neigs,eigs_tmp,eig_vecs_tmp,omega,LOGICAL(include_bounds))
oft_env%pm=pm_save
END SUBROUTINE tokamaker_eig_td
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_eig_wall(neigs,eigs,eig_vecs,pm) BIND(C,NAME="tokamaker_eig_wall")
INTEGER(c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vecs !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: pm !< Needs docs
REAL(8) :: alam_save,pnorm_save
REAL(8), POINTER :: eigs_tmp(:,:),eig_vecs_tmp(:,:)
LOGICAL :: pm_save
CALL c_f_pointer(eigs, eigs_tmp, [2,neigs])
CALL c_f_pointer(eig_vecs, eig_vecs_tmp, [gs_global%psi%n,neigs])
alam_save=gs_global%alam; gs_global%alam=0.d0
pnorm_save=gs_global%pnorm; gs_global%pnorm=0.d0
pm_save=oft_env%pm; oft_env%pm=pm
CALL eig_gs_td(gs_global,neigs,eigs_tmp,eig_vecs_tmp,0.d0,.FALSE.)
oft_env%pm=pm_save
gs_global%alam=alam_save
gs_global%pnorm=pnorm_save
END SUBROUTINE tokamaker_eig_wall
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_step_td(time,dt,nl_its,lin_its,nretry) BIND(C,NAME="tokamaker_step_td")
REAL(c_double), INTENT(inout) :: time !< Needs docs
REAL(c_double), INTENT(inout) :: dt !< Needs docs
INTEGER(c_int), INTENT(out) :: nl_its !< Needs docs
INTEGER(c_int), INTENT(out) :: lin_its !< Needs docs
INTEGER(c_int), INTENT(out) :: nretry !< Needs docs
CALL gs_td_global%step(time,dt,nl_its,lin_its,nretry)
END SUBROUTINE tokamaker_step_td
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_mesh(np,r_loc,nc,lc_loc,reg_loc) BIND(C,NAME="tokamaker_get_mesh")
TYPE(c_ptr), INTENT(out) :: lc_loc !< Needs docs
TYPE(c_ptr), INTENT(out) :: r_loc !< Needs docs
TYPE(c_ptr), INTENT(out) :: reg_loc !< Needs docs
INTEGER(c_int), INTENT(out) :: np !< Needs docs
INTEGER(c_int), INTENT(out) :: nc !< Needs docs
INTEGER(4) :: i,j,k,id
CALL smesh%tessellate(r_plot, lc_plot, smesh%tess_order)
np=SIZE(r_plot,DIM=2,KIND=c_int)
nc=SIZE(lc_plot,DIM=2,KIND=c_int)
r_loc=c_loc(r_plot)
lc_loc=c_loc(lc_plot)
!
ALLOCATE(reg_plot(nc))
k=nc/smesh%nc
IF(ASSOCIATED(smesh%reg))THEN
  !$omp parallel do private(j,id)
  DO i=1,smesh%nc
    id=smesh%reg(i)
    DO j=1,k
      reg_plot((i-1)*k+j)=id
    END DO
  END DO
ELSE
  reg_plot=0
END IF
reg_loc=c_loc(reg_plot)
END SUBROUTINE tokamaker_get_mesh
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_limiter(np,r_loc) BIND(C,NAME="tokamaker_get_limiter")
TYPE(c_ptr), INTENT(out) :: r_loc !< Needs docs
INTEGER(c_int), INTENT(out) :: np !< Needs docs
INTEGER(4) :: i
REAL(8), POINTER, DIMENSION(:,:) :: r_tmp
np=gs_global%nlim_con+1
ALLOCATE(r_tmp(2,gs_global%nlim_con+1))
r_loc=C_LOC(r_tmp)
DO i=1,gs_global%nlim_con
  r_tmp(:,i)=smesh%r(1:2,gs_global%lim_con(i))
END DO
r_tmp(:,gs_global%nlim_con+1)=smesh%r(1:2,gs_global%lim_con(1))
END SUBROUTINE tokamaker_get_limiter
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_psi(psi_vals,psi_lim,psi_max) BIND(C,NAME="tokamaker_get_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(c_double), INTENT(out) :: psi_lim !< Needs docs
REAL(c_double), INTENT(out) :: psi_max !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
CALL gs_global%psi%get_local(vals_tmp)
psi_lim = gs_global%plasma_bounds(1)
psi_max = gs_global%plasma_bounds(2)
END SUBROUTINE tokamaker_get_psi
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_dels_curr(psi_vals) BIND(C,NAME="tokamaker_get_dels_curr")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_solver), POINTER :: minv
IF(.NOT.ASSOCIATED(gs_global%dels_full))CALL build_dels(gs_global%dels_full,gs_global,"none")
!
CALL gs_global%psi%new(u)
CALL gs_global%psi%new(v)
CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
CALL u%restore_local(vals_tmp)
!
NULLIFY(minv)
CALL create_cg_solver(minv)
minv%A=>gs_global%mop
minv%its=-2
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
CALL gs_global%dels_full%apply(u,v)
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_coil_currents(currents,reg_currents) BIND(C,NAME="tokamaker_get_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: currents !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_currents !< Needs docs
INTEGER(4) :: i,j
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp,coil_regs
CALL c_f_pointer(reg_currents, coil_regs, [smesh%nreg])
CALL c_f_pointer(currents, vals_tmp, [gs_global%ncoils])
vals_tmp=(gs_global%coil_currs + gs_global%coil_vcont*gs_global%vcontrol_val)/mu0
coil_regs = 0.d0
DO j=1,gs_global%ncoil_regs
  DO i=1,gs_global%ncoils
    coil_regs(gs_global%coil_regions(j)%id) = coil_regs(gs_global%coil_regions(j)%id) &
      + vals_tmp(i)*gs_global%coil_nturns(gs_global%coil_regions(j)%id,i)
  END DO
  coil_regs(gs_global%coil_regions(j)%id) = coil_regs(gs_global%coil_regions(j)%id)*gs_global%coil_regions(j)%area
END DO
END SUBROUTINE tokamaker_get_coil_currents
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_coil_Lmat(Lmat) BIND(C,NAME="tokamaker_get_coil_Lmat")
TYPE(c_ptr), VALUE, INTENT(in) :: Lmat !< Needs docs
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
REAL(8) :: tmp1,tmp2,tmp3,itor
CLASS(oft_vector), POINTER :: rhs,vec1,vec2
!---Update plasma row/column
IF(gs_global%has_plasma)THEN
  DO i=1,gs_global%ncoils
    CALL gs_plasma_mutual(gs_global,gs_global%psi_coil(i)%f,gs_global%Lcoils(i,gs_global%ncoils+1),itor)
    gs_global%Lcoils(gs_global%ncoils+1,i)=gs_global%Lcoils(i,gs_global%ncoils+1)
  END DO
  !
  CALL gs_global%psi%new(rhs)
  CALL gs_global%psi%new(vec1)
  CALL gs_global%psi%new(vec2)
  CALL gs_source(gs_global,gs_global%psi,rhs,vec1,vec2,tmp1,tmp2,tmp3)
  CALL vec1%set(0.d0)
  CALL gs_global%lu_solver%apply(vec1,rhs)
  CALL gs_plasma_mutual(gs_global,vec1,gs_global%Lcoils(gs_global%ncoils+1,gs_global%ncoils+1),itor)
  gs_global%Lcoils(gs_global%ncoils+1,gs_global%ncoils+1)=gs_global%Lcoils(gs_global%ncoils+1,gs_global%ncoils+1)/itor
  CALL rhs%delete()
  CALL vec1%delete()
  CALL vec2%delete()
  DEALLOCATE(rhs,vec1,vec2)
ELSE
  gs_global%Lcoils(gs_global%ncoils+1,:)=0.d0
  gs_global%Lcoils(:,gs_global%ncoils+1)=0.d0
END IF
!---Copy out inductance matrix
CALL c_f_pointer(Lmat, vals_tmp, [gs_global%ncoils+1,gs_global%ncoils+1])
vals_tmp=gs_global%Lcoils
END SUBROUTINE tokamaker_get_coil_Lmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_refs(o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm) BIND(C,NAME="tokamaker_get_refs")
TYPE(c_ptr), INTENT(out) :: o_point !< Needs docs
TYPE(c_ptr), INTENT(out) :: lim_point !< Needs docs
TYPE(c_ptr), INTENT(out) :: x_points !< Needs docs
TYPE(c_ptr), INTENT(out) :: diverted !< Needs docs
TYPE(c_ptr), INTENT(out) :: plasma_bounds !< Needs docs
TYPE(c_ptr), INTENT(out) :: alam !< Needs docs
TYPE(c_ptr), INTENT(out) :: pnorm !< Needs docs
o_point=c_loc(gs_global%o_point)
lim_point=c_loc(gs_global%lim_point)
x_points=c_loc(gs_global%x_points)
diverted=c_loc(gs_global%diverted)
plasma_bounds=c_loc(gs_global%plasma_bounds)
alam=c_loc(gs_global%alam)
pnorm=c_loc(gs_global%pnorm)
END SUBROUTINE tokamaker_get_refs
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_trace_surf(psi_surf,points,npoints) BIND(C,NAME="tokamaker_trace_surf")
REAL(c_double), VALUE, INTENT(in) :: psi_surf !< Needs docs
TYPE(c_ptr), INTENT(out) ::  points !< Needs docs
INTEGER(c_int), INTENT(out) :: npoints !< Needs docs
REAL(8), POINTER, DIMENSION(:,:) :: pts_tmp
CALL gs_trace_surf(gs_global,psi_surf,pts_tmp,npoints)
IF(npoints>0)THEN
  points = c_loc(pts_tmp)
ELSE
  points = C_NULL_PTR
END IF
END SUBROUTINE tokamaker_trace_surf
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_q(npsi,psi_q,qvals,ravgs,dl,rbounds,zbounds) BIND(C,NAME="tokamaker_get_q")
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_q(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: qvals(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: ravgs(npsi,2) !< Needs docs
REAL(c_double), INTENT(out) :: dl !< Needs docs
REAL(c_double), INTENT(out) :: rbounds(2,2) !< Needs docs
REAL(c_double), INTENT(out) :: zbounds(2,2) !< Needs docs
CALL gs_get_qprof(gs_global,npsi,psi_q,qvals,dl,rbounds,zbounds,ravgs)
END SUBROUTINE tokamaker_get_q
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_sauter_fc(npsi,psi_saut,fc,r_avgs,modb_avgs) BIND(C,NAME="tokamaker_sauter_fc")
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_saut(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: fc(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: r_avgs(npsi,3) !< Needs docs
REAL(c_double), INTENT(out) :: modb_avgs(npsi,2) !< Needs docs
CALL sauter_fc(gs_global,npsi,psi_saut,fc,r_avgs,modb_avgs)
END SUBROUTINE tokamaker_sauter_fc
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_globals(Itor,centroid,vol,pvol,dflux,tflux,bp_vol) BIND(C,NAME="tokamaker_get_globals")
REAL(c_double), INTENT(out) :: Itor !< Needs docs
REAL(c_double), INTENT(out) :: centroid(2) !< Needs docs
REAL(c_double), INTENT(out) :: vol !< Needs docs
REAL(c_double), INTENT(out) :: pvol !< Needs docs
REAL(c_double), INTENT(out) :: dflux !< Needs docs
REAL(c_double), INTENT(out) :: tflux !< Needs docs
REAL(c_double), INTENT(out) :: bp_vol !< Needs docs
CALL gs_comp_globals(gs_global,Itor,centroid,vol,pvol,dflux,tflux,bp_vol)
Itor=Itor/mu0
vol=vol*2.d0*pi
pvol=pvol*2.d0*pi/mu0
END SUBROUTINE tokamaker_get_globals
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_gs_calc_vloop(vloop) BIND(C,NAME="tokamaker_gs_calc_vloop")
REAL(c_double), INTENT(out) :: vloop
IF(.NOT.ASSOCIATED(gs_global%eta))THEN
  vloop=-1.d0
  RETURN
END IF
CALL gs_calc_vloop(gs_global,vloop)
END SUBROUTINE tokamaker_gs_calc_vloop
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_profs(npsi,psi_in,f,fp,p,pp) BIND(C,NAME="tokamaker_get_profs")
INTEGER(c_int), VALUE, INTENT(in) :: npsi !< Needs docs
REAL(c_double), INTENT(in) :: psi_in(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: f(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: fp(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: p(npsi) !< Needs docs
REAL(c_double), INTENT(out) :: pp(npsi) !< Needs docs
INTEGER(4) :: i
REAL(8) :: x1,x2,r
x1=0.d0; x2=1.d0
IF(gs_global%plasma_bounds(1)>-1.d98)THEN
  x1=gs_global%plasma_bounds(1); x2=gs_global%plasma_bounds(2)
END IF
DO i=1,npsi
  r=psi_in(i)*(x2-x1) + x1
  IF(gs_global%mode==0)THEN
    fp(i)=gs_global%alam*gs_global%I%fp(r)
    f(i)=gs_global%psiscale*gs_global%alam*gs_global%I%f(r) + gs_global%I%f_offset
  ELSE
    f(i)=SQRT(gs_global%psiscale*gs_global%alam*gs_global%I%f(r) + gs_global%I%f_offset**2)
    fp(i)=gs_global%alam*gs_global%I%fp(r)/(2.d0*f(i))
  END IF
  pp(i)=gs_global%psiscale*gs_global%pnorm*gs_global%P%fp(r)
  p(i)=gs_global%psiscale*gs_global%psiscale*gs_global%pnorm*gs_global%P%f(r)
END DO
END SUBROUTINE tokamaker_get_profs
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_vfixed(npts,pts,fluxes) BIND(C,NAME="tokamaker_get_vfixed")
INTEGER(c_int), INTENT(out) :: npts !< Needs docs
TYPE(c_ptr), INTENT(out) :: pts !< Needs docs
TYPE(c_ptr), INTENT(out) :: fluxes !< Needs docs
REAL(8), POINTER :: pts_tmp(:,:),fluxes_tmp(:)
CALL gs_fixed_vflux(gs_global,pts_tmp,fluxes_tmp)
pts=C_LOC(pts_tmp)
fluxes=C_LOC(fluxes_tmp)
npts=SIZE(fluxes_tmp,DIM=1)
END SUBROUTINE tokamaker_get_vfixed
!------------------------------------------------------------------------------
!> Create an interpolation object for tokamaker fields
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_get_field_eval(imode,int_obj,error_str) BIND(C,NAME="tokamaker_get_field_eval")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Field type
TYPE(c_ptr), INTENT(out) :: int_obj !< Pointer to interpolation object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Error string (unused)
TYPE(gs_prof_interp), POINTER :: prof_interp_obj
TYPE(gs_b_interp), POINTER :: b_interp_obj
CALL copy_string('',error_str)
IF(imode==1)THEN
  ALLOCATE(b_interp_obj)
  b_interp_obj%gs=>gs_global
  CALL b_interp_obj%setup()
  int_obj=C_LOC(b_interp_obj)
ELSE
  ALLOCATE(prof_interp_obj)
  prof_interp_obj%gs=>gs_global
  prof_interp_obj%mode=imode-1
  CALL prof_interp_obj%setup()
  int_obj=C_LOC(prof_interp_obj)
END IF
END SUBROUTINE tokamaker_get_field_eval
!------------------------------------------------------------------------------
!> Evaluate a TokaMaker field with an interpolation object created by
!! \ref tokamaker_f::tokamaker_get_field_eval
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_apply_field_eval(int_obj,int_type,pt,fbary_tol,cell,dim,field) BIND(C,NAME="tokamaker_apply_field_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Pointer to interpolation object
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Field type (negative to destroy)
REAL(c_double), INTENT(in) :: pt(3) !< Location for evaluation [R,Z,0]
REAL(c_double), VALUE, INTENT(in) :: fbary_tol !< Tolerance for physical to logical mapping
INTEGER(c_int), INTENT(inout) :: cell !< Cell containing `pt` (starting guess on input)
INTEGER(c_int), VALUE, INTENT(in) :: dim !< Dimension of field
REAL(c_double), INTENT(out) :: field(dim) !< Field at `pt`
TYPE(gs_prof_interp), POINTER :: prof_interp_obj
TYPE(gs_b_interp), POINTER :: b_interp_obj
REAL(8) :: f(4),goptmp(3,4),vol,fmin,fmax
IF(int_type<0)THEN
  IF(int_type==-1)THEN
    CALL c_f_pointer(int_obj, b_interp_obj)
    CALL b_interp_obj%delete
  ELSE
    CALL c_f_pointer(int_obj, prof_interp_obj)
    CALL prof_interp_obj%delete
  END IF
  RETURN
END IF
call bmesh_findcell(smesh,cell,pt,f)
IF(cell==0)RETURN
fmin=MINVAL(f); fmax=MAXVAL(f)
IF(( fmax>1.d0+fbary_tol ).OR.( fmin<-fbary_tol ))THEN
  cell=-ABS(cell)
  RETURN
END IF
CALL smesh%jacobian(cell,f,goptmp,vol)
IF(int_type==1)THEN
  CALL c_f_pointer(int_obj, b_interp_obj)
  CALL b_interp_obj%interp(cell,f,goptmp,field)
ELSE
  CALL c_f_pointer(int_obj, prof_interp_obj)
  CALL prof_interp_obj%interp(cell,f,goptmp,field)
END IF
END SUBROUTINE tokamaker_apply_field_eval
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_psi(psi_vals) BIND(C,NAME="tokamaker_set_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
CALL gs_global%psi%restore_local(vals_tmp)
END SUBROUTINE tokamaker_set_psi
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_psi_dt(psi_vals,dt) BIND(C,NAME="tokamaker_set_psi_dt")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: dt !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
gs_global%dt=dt
IF(dt>0.d0)THEN
  IF(.NOT.ASSOCIATED(gs_global%psi_dt))CALL gs_global%psi%new(gs_global%psi_dt)
  CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
  CALL gs_global%psi_dt%restore_local(vals_tmp)
ELSE
  IF(ASSOCIATED(gs_global%psi_dt))THEN
    CALL gs_global%psi_dt%delete()
    DEALLOCATE(gs_global%psi_dt)
  END IF
END IF
END SUBROUTINE tokamaker_set_psi_dt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_settings(settings) BIND(C,NAME="tokamaker_set_settings")
TYPE(tokamaker_settings_type), INTENT(in) :: settings !< Needs docs
oft_env%pm=settings%pm
gs_global%free=settings%free_boundary
gs_global%lim_zmax=settings%lim_zmax
gs_global%rmin=settings%rmin
gs_global%mode=settings%mode
gs_global%urf=settings%urf
gs_global%maxits=settings%maxits
gs_global%nl_tol=settings%nl_tol
gs_global%limited_only=settings%limited_only
CALL copy_string_rev(settings%limiter_file,gs_global%limiter_file)
END SUBROUTINE tokamaker_set_settings
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_targets(ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target) BIND(C,NAME="tokamaker_set_targets")
REAL(c_double), VALUE, INTENT(in) :: ip_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: ip_ratio_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: pax_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: estore_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: R0_target !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: V0_target !< Needs docs
gs_global%R0_target=R0_target
gs_global%V0_target=V0_target
gs_global%pax_target=pax_target*mu0
gs_global%estore_target=estore_target*mu0
gs_global%itor_target=ip_target*mu0
gs_global%ip_ratio_target=ip_ratio_target
END SUBROUTINE tokamaker_set_targets
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_isoflux(targets,weights,ntargets,grad_wt_lim) BIND(C,NAME="tokamaker_set_isoflux")
REAL(c_double), INTENT(in) :: targets(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: grad_wt_lim !< Needs docs
IF(ASSOCIATED(gs_global%isoflux_targets))DEALLOCATE(gs_global%isoflux_targets)
gs_global%isoflux_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(gs_global%isoflux_targets(3,gs_global%isoflux_ntargets))
  gs_global%isoflux_targets(1:2,:)=targets
  gs_global%isoflux_targets(3,:)=weights
  gs_global%isoflux_grad_wt_lim=1.d0/grad_wt_lim
ELSE
  gs_global%isoflux_grad_wt_lim=-1.d0
END IF
END SUBROUTINE tokamaker_set_isoflux
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_flux(locations,targets,weights,ntargets,grad_wt_lim) BIND(C,NAME="tokamaker_set_flux")
REAL(c_double), INTENT(in) :: locations(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: targets(ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: grad_wt_lim !< Needs docs
IF(ASSOCIATED(gs_global%flux_targets))DEALLOCATE(gs_global%flux_targets)
gs_global%flux_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(gs_global%flux_targets(4,gs_global%flux_ntargets))
  gs_global%flux_targets(1:2,:)=locations
  gs_global%flux_targets(3,:)=targets
  gs_global%flux_targets(4,:)=weights
  ! gs_global%isoflux_grad_wt_lim=1.d0/grad_wt_lim
! ELSE
  ! gs_global%isoflux_grad_wt_lim=-1.d0
END IF
END SUBROUTINE tokamaker_set_flux
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_saddles(targets,weights,ntargets) BIND(C,NAME="tokamaker_set_saddles")
REAL(c_double), INTENT(in) :: targets(2,ntargets) !< Needs docs
REAL(c_double), INTENT(in) :: weights(ntargets) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ntargets !< Needs docs
IF(ASSOCIATED(gs_global%saddle_targets))DEALLOCATE(gs_global%saddle_targets)
gs_global%saddle_ntargets=ntargets
IF(ntargets>0)THEN
  ALLOCATE(gs_global%saddle_targets(3,gs_global%saddle_ntargets))
  gs_global%saddle_targets(1:2,:)=targets
  gs_global%saddle_targets(3,:)=weights
END IF
END SUBROUTINE tokamaker_set_saddles
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_currents(currents) BIND(C,NAME="tokamaker_set_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: currents !< Needs docs
INTEGER(4) :: i
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(currents, vals_tmp, [gs_global%ncoils])
gs_global%coil_currs = vals_tmp*mu0
gs_global%vcontrol_val = 0.d0
END SUBROUTINE tokamaker_set_coil_currents
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_regmat(nregularize,coil_reg_mat,coil_reg_targets,coil_reg_weights) BIND(C,NAME="tokamaker_set_coil_regmat")
INTEGER(c_int), VALUE, INTENT(in) :: nregularize !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_mat !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_targets !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_weights !< Needs docs
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
IF(ASSOCIATED(gs_global%coil_reg_mat))DEALLOCATE(gs_global%coil_reg_mat,gs_global%coil_reg_targets)
gs_global%nregularize=nregularize
ALLOCATE(gs_global%coil_reg_mat(gs_global%nregularize,gs_global%ncoils+1))
ALLOCATE(gs_global%coil_reg_targets(gs_global%nregularize))
CALL c_f_pointer(coil_reg_mat, vals_tmp, [gs_global%nregularize,gs_global%ncoils+1])
gs_global%coil_reg_mat=vals_tmp
CALL c_f_pointer(coil_reg_targets, vals_tmp, [gs_global%nregularize,1])
gs_global%coil_reg_targets=vals_tmp(:,1)*mu0
CALL c_f_pointer(coil_reg_weights, vals_tmp, [gs_global%nregularize,1])
DO i=1,gs_global%nregularize
  gs_global%coil_reg_targets(i)=gs_global%coil_reg_targets(i)*vals_tmp(i,1)
  gs_global%coil_reg_mat(i,:)=gs_global%coil_reg_mat(i,:)*vals_tmp(i,1)
END DO
END SUBROUTINE tokamaker_set_coil_regmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_bounds(coil_bounds) BIND(C,NAME="tokamaker_set_coil_bounds")
TYPE(c_ptr), VALUE, INTENT(in) :: coil_bounds !< Needs docs
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
CALL c_f_pointer(coil_bounds, vals_tmp, [2,gs_global%ncoils+1])
IF(.NOT.ASSOCIATED(gs_global%coil_bounds))THEN
  ALLOCATE(gs_global%coil_bounds(2,gs_global%ncoils+1))
  gs_global%coil_bounds(1,:)=-1.d98; gs_global%coil_bounds(2,:)=1.d98
END IF
DO i=1,gs_global%ncoils
  gs_global%coil_bounds([2,1],i)=-vals_tmp(:,i)*mu0
END DO
gs_global%coil_bounds([2,1],gs_global%ncoils+1)=-vals_tmp(:,gs_global%ncoils+1)*mu0
END SUBROUTINE tokamaker_set_coil_bounds
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_set_coil_vsc(coil_gains) BIND(C,NAME="tokamaker_set_coil_vsc")
TYPE(c_ptr), VALUE, INTENT(in) :: coil_gains !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
INTEGER(4) :: i
CALL c_f_pointer(coil_gains, vals_tmp, [gs_global%ncoils])
gs_global%coil_vcont=vals_tmp
END SUBROUTINE tokamaker_set_coil_vsc
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE tokamaker_save_eqdsk(filename,nr,nz,rbounds,zbounds,run_info,psi_pad,rcentr,error_str) BIND(C,NAME="tokamaker_save_eqdsk")
CHARACTER(KIND=c_char), INTENT(in) :: filename(80) !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: run_info(36) !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nr !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nz !< Needs docs
REAL(c_double), INTENT(in) :: rbounds(2) !< Needs docs
REAL(c_double), INTENT(in) :: zbounds(2) !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: psi_pad !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: rcentr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
CHARACTER(LEN=36) :: run_info_f
CHARACTER(LEN=80) :: filename_tmp,lim_file,error_flag
CALL copy_string_rev(run_info,run_info_f)
CALL copy_string_rev(filename,filename_tmp)
lim_file='none'
IF(rcentr>0.d0)THEN
  CALL gs_save_eqdsk(gs_global,filename_tmp,nr,nz,rbounds,zbounds,run_info_f,lim_file,psi_pad,rcentr_in=rcentr,error_str=error_flag)
ELSE
  CALL gs_save_eqdsk(gs_global,filename_tmp,nr,nz,rbounds,zbounds,run_info_f,lim_file,psi_pad,error_str=error_flag)
END IF
CALL copy_string(TRIM(error_flag),error_str)
END SUBROUTINE tokamaker_save_eqdsk
END MODULE tokamaker_f