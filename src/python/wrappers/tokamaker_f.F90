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
USE oft_mesh_type, ONLY: smesh
! USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid, ONLY: multigrid_reset
! USE multigrid_build, ONLY: multigrid_construct_surf
USE fem_base, ONLY: oft_afem_type
USE oft_la_base, ONLY: oft_vector
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, oft_blagrange, &
  oft_lag_setup
USE mhd_utils, ONLY: mu0
USE axi_green, ONLY: green
USE oft_gs, ONLY: gs_eq, gs_save_fields, gs_save_fgrid, gs_setup_walls, gs_save_prof, &
  gs_fixed_vflux, gs_load_regions, gs_get_qprof, gs_trace_surf
USE oft_gs_util, ONLY: gs_save, gs_load, gs_analyze, gs_comp_globals, gs_save_eqdsk, &
  gs_profile_load, sauter_fc, gs_calc_vloop
USE oft_gs_td, ONLY: setup_gs_td, step_gs_td, eig_gs_td
USE oft_base_f, ONLY: copy_string, copy_string_rev, oftpy_init
IMPLICIT NONE
!
TYPE, BIND(C) :: tokamaker_settings_type
  LOGICAL(KIND=c_bool) :: pm = .FALSE.
  LOGICAL(KIND=c_bool) :: free_boundary = .FALSE.
  LOGICAL(KIND=c_bool) :: has_plasma = .TRUE.
  LOGICAL(KIND=c_bool) :: limited_only = .FALSE.
  INTEGER(KIND=c_int) :: maxits = 40
  INTEGER(KIND=c_int) :: mode = 1
  REAL(KIND=c_double) :: urf = 0.3d0
  REAL(KIND=c_double) :: nl_tol = 1.d-6
  REAL(KIND=c_double) :: rmin = 0.d0
  REAL(KIND=c_double) :: lim_zmax = 1.d99
  CHARACTER(KIND=c_char) :: limiter_file(80) = 'none'
END TYPE tokamaker_settings_type
!
TYPE(gs_eq), POINTER :: gs_global => NULL() !< Global G-S object
integer(i4), POINTER :: lc_plot(:,:),reg_plot(:)
real(r8), POINTER :: r_plot(:,:)
CONTAINS
!
SUBROUTINE tokamaker_alloc(gs_obj) BIND(C,NAME="tokamaker_alloc")
TYPE(c_ptr), INTENT(out) :: gs_obj
IF(.NOT.ASSOCIATED(gs_global))ALLOCATE(gs_global)
gs_obj=c_loc(gs_global) ! Do nothing with this for now
END SUBROUTINE tokamaker_alloc
!
SUBROUTINE tokamaker_eval_green(n,r,z,rc,zc,vals) BIND(C,NAME="tokamaker_eval_green")
INTEGER(c_int), VALUE, INTENT(in) :: n
REAL(c_double), VALUE, INTENT(in) :: rc,zc
TYPE(c_ptr), VALUE, INTENT(in) :: r,z,vals
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
!
SUBROUTINE tokamaker_setup_regions(coil_file,reg_eta) BIND(C,NAME="tokamaker_setup_regions")
CHARACTER(KIND=c_char), INTENT(in) :: coil_file(80)
TYPE(c_ptr), VALUE, INTENT(in) :: reg_eta
real(r8), POINTER :: eta_tmp(:)
INTEGER(4) :: i
CALL copy_string_rev(coil_file,gs_global%coil_file)
IF(TRIM(gs_global%coil_file)=='none')THEN
  CALL c_f_pointer(reg_eta, eta_tmp, [smesh%nreg])
  !
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
    ELSE
      gs_global%ncoil_regs=gs_global%ncoil_regs+1
      gs_global%coil_regions(gs_global%ncoil_regs)%id=i
    END IF
  END DO
ELSE
  CALL gs_load_regions(gs_global)
END IF
END SUBROUTINE tokamaker_setup_regions
!
SUBROUTINE tokamaker_reset(error_str) BIND(C,NAME="tokamaker_reset")
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80)
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
!
SUBROUTINE tokamaker_setup(order,full_domain,ncoils,error_str) BIND(C,NAME="tokamaker_setup")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: full_domain
INTEGER(KIND=c_int), INTENT(out) :: ncoils
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80)
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
ncoils=gs_global%ncoil_regs
END SUBROUTINE tokamaker_setup
!
SUBROUTINE tokamaker_load_profiles(f_file,f_offset,p_file,eta_file,f_NI_file) BIND(C,NAME="tokamaker_load_profiles")
CHARACTER(KIND=c_char), INTENT(in) :: f_file(80),p_file(80),eta_file(80),f_NI_file(80)
REAL(c_double), VALUE, INTENT(in) :: f_offset
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
!
SUBROUTINE tokamaker_init_psi(r0,z0,a,kappa,delta,ierr) BIND(C,NAME="tokamaker_init_psi")
REAL(c_double), VALUE, INTENT(in) :: r0,z0,a,kappa,delta
INTEGER(c_int), INTENT(out) :: ierr
CALL gs_global%init_psi(ierr,r0=[r0,z0],a=a,kappa=kappa,delta=delta)
END SUBROUTINE tokamaker_init_psi
!
SUBROUTINE tokamaker_run(vacuum,error_flag) BIND(C,NAME="tokamaker_run")
LOGICAL(c_bool), VALUE, INTENT(in) :: vacuum
INTEGER(c_int), INTENT(out) :: error_flag
IF(vacuum)gs_global%has_plasma=.FALSE.
CALL gs_global%solve(error_flag)
gs_global%has_plasma=.TRUE.
END SUBROUTINE tokamaker_run
!
SUBROUTINE tokamaker_analyze() BIND(C,NAME="tokamaker_analyze")
CALL gs_analyze(gs_global)
END SUBROUTINE tokamaker_analyze
!
SUBROUTINE tokamaker_setup_td(dt,lin_tol,nl_tol) BIND(C,NAME="tokamaker_setup_td")
REAL(c_double), VALUE, INTENT(in) :: dt,lin_tol,nl_tol
CALL setup_gs_td(gs_global,dt,lin_tol,nl_tol)
END SUBROUTINE tokamaker_setup_td
!
SUBROUTINE tokamaker_eig_td(omega,neigs,eigs,eig_vecs,include_bounds,pm) BIND(C,NAME="tokamaker_eig_td")
REAL(c_double), VALUE, INTENT(in) :: omega
INTEGER(c_int), VALUE, INTENT(in) :: neigs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs,eig_vecs
LOGICAL(c_bool), VALUE, INTENT(in) :: include_bounds,pm
REAL(8), POINTER :: eigs_tmp(:,:),eig_vecs_tmp(:,:)
LOGICAL :: pm_save
CALL c_f_pointer(eigs, eigs_tmp, [2,neigs])
CALL c_f_pointer(eig_vecs, eig_vecs_tmp, [gs_global%psi%n,neigs])
pm_save=oft_env%pm; oft_env%pm=pm
CALL eig_gs_td(gs_global,neigs,eigs_tmp,eig_vecs_tmp,omega,LOGICAL(include_bounds))
oft_env%pm=pm_save
END SUBROUTINE tokamaker_eig_td
!
SUBROUTINE tokamaker_eig_wall(neigs,eigs,eig_vecs,pm) BIND(C,NAME="tokamaker_eig_wall")
INTEGER(c_int), VALUE, INTENT(in) :: neigs
TYPE(c_ptr), VALUE, INTENT(in) :: eigs,eig_vecs
LOGICAL(c_bool), VALUE, INTENT(in) :: pm
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
!
SUBROUTINE tokamaker_step_td(time,dt,nl_its,lin_its,nretry) BIND(C,NAME="tokamaker_step_td")
REAL(c_double), INTENT(inout) :: time,dt
INTEGER(c_int), INTENT(out) :: nl_its,lin_its,nretry
CALL step_gs_td(time,dt,nl_its,lin_its,nretry)
END SUBROUTINE tokamaker_step_td
!
SUBROUTINE tokamaker_get_mesh(np,r_loc,nc,lc_loc,reg_loc) BIND(C,NAME="tokamaker_get_mesh")
TYPE(c_ptr), INTENT(out) :: lc_loc,r_loc,reg_loc
INTEGER(c_int), INTENT(out) :: np,nc
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
!
SUBROUTINE tokamaker_get_limiter(np,r_loc) BIND(C,NAME="tokamaker_get_limiter")
TYPE(c_ptr), INTENT(out) :: r_loc
INTEGER(c_int), INTENT(out) :: np
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
!
SUBROUTINE tokamaker_get_psi(psi_vals,psi_lim,psi_max) BIND(C,NAME="tokamaker_get_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals
REAL(c_double), INTENT(out) :: psi_lim,psi_max
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
CALL gs_global%psi%get_local(vals_tmp)
psi_lim = gs_global%plasma_bounds(1)
psi_max = gs_global%plasma_bounds(2)
END SUBROUTINE tokamaker_get_psi
!
SUBROUTINE tokamaker_get_coil_currents(currents,coil_map) BIND(C,NAME="tokamaker_get_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: currents,coil_map
INTEGER(4) :: i
INTEGER(4), POINTER, DIMENSION(:) :: coil_regs
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(coil_map, coil_regs, [gs_global%ncoil_regs])
CALL c_f_pointer(currents, vals_tmp, [gs_global%ncoil_regs])
DO i=1,gs_global%ncoil_regs
  curr = gs_global%coil_regions(i)%curr & 
    + gs_global%coil_regions(i)%vcont_gain*gs_global%vcontrol_val
  vals_tmp(i)=curr*gs_global%coil_regions(i)%area/mu0
  coil_regs(i)=gs_global%coil_regions(i)%id
END DO
END SUBROUTINE tokamaker_get_coil_currents
!
SUBROUTINE tokamaker_get_coil_Lmat(Lmat) BIND(C,NAME="tokamaker_get_coil_Lmat")
TYPE(c_ptr), VALUE, INTENT(in) :: Lmat
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
CALL c_f_pointer(Lmat, vals_tmp, [gs_global%ncoil_regs+1,gs_global%ncoil_regs+1])
vals_tmp=gs_global%Lcoils
END SUBROUTINE tokamaker_get_coil_Lmat
!
SUBROUTINE tokamaker_get_refs(o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm) BIND(C,NAME="tokamaker_get_refs")
TYPE(c_ptr), INTENT(out) :: o_point,lim_point,x_points,diverted,plasma_bounds,alam,pnorm
o_point=c_loc(gs_global%o_point)
lim_point=c_loc(gs_global%lim_point)
x_points=c_loc(gs_global%x_points)
diverted=c_loc(gs_global%diverted)
plasma_bounds=c_loc(gs_global%plasma_bounds)
alam=c_loc(gs_global%alam)
pnorm=c_loc(gs_global%pnorm)
END SUBROUTINE tokamaker_get_refs
!
SUBROUTINE tokamaker_trace_surf(psi_surf,points,npoints) BIND(C,NAME="tokamaker_trace_surf")
REAL(c_double), VALUE, INTENT(in) :: psi_surf
TYPE(c_ptr), INTENT(out) ::  points
INTEGER(c_int), INTENT(out) :: npoints
REAL(8), POINTER, DIMENSION(:,:) :: pts_tmp
CALL gs_trace_surf(gs_global,psi_surf,pts_tmp,npoints)
IF(npoints>0)THEN
  points = c_loc(pts_tmp)
ELSE
  points = C_NULL_PTR
END IF
END SUBROUTINE tokamaker_trace_surf
!
SUBROUTINE tokamaker_get_q(npsi,psi_q,qvals,ravgs,dl,rbounds,zbounds) BIND(C,NAME="tokamaker_get_q")
INTEGER(c_int), VALUE, INTENT(in) :: npsi
REAL(c_double), INTENT(in) :: psi_q(npsi)
REAL(c_double), INTENT(out) :: qvals(npsi),ravgs(npsi,2),dl,rbounds(2,2),zbounds(2,2)
CALL gs_get_qprof(gs_global,npsi,psi_q,qvals,dl,rbounds,zbounds,ravgs)
END SUBROUTINE tokamaker_get_q
!
SUBROUTINE tokamaker_sauter_fc(npsi,psi_saut,fc,r_avgs,modb_avgs) BIND(C,NAME="tokamaker_sauter_fc")
INTEGER(c_int), VALUE, INTENT(in) :: npsi
REAL(c_double), INTENT(in) :: psi_saut(npsi)
REAL(c_double), INTENT(out) :: fc(npsi),r_avgs(npsi,3),modb_avgs(npsi,2)
CALL sauter_fc(gs_global,npsi,psi_saut,fc,r_avgs,modb_avgs)
END SUBROUTINE tokamaker_sauter_fc
!
SUBROUTINE tokamaker_get_globals(Itor,centroid,vol,pvol,dflux,tflux,li) BIND(C,NAME="tokamaker_get_globals")
REAL(c_double), INTENT(out) :: Itor,centroid(2),vol,pvol,dflux,tflux,li
CALL gs_comp_globals(gs_global,Itor,centroid,vol,pvol,dflux,tflux,li)
Itor=Itor/mu0
vol=vol*2.d0*pi
pvol=pvol*2.d0*pi/mu0
li=li*2.d0/(mu0*gs_global%o_point(1))
END SUBROUTINE tokamaker_get_globals
!
SUBROUTINE tokamaker_gs_calc_vloop(vloop) BIND(C,NAME="tokamaker_gs_calc_vloop")
REAL(c_double), INTENT(out) :: vloop
CALL gs_calc_vloop(gs_global,vloop)
END SUBROUTINE tokamaker_gs_calc_vloop
!
SUBROUTINE tokamaker_get_profs(npsi,psi_in,f,fp,p,pp) BIND(C,NAME="tokamaker_get_profs")
INTEGER(c_int), VALUE, INTENT(in) :: npsi
REAL(c_double), INTENT(in) :: psi_in(npsi)
REAL(c_double), INTENT(out) :: f(npsi),fp(npsi),p(npsi),pp(npsi)
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
!
SUBROUTINE tokamaker_get_vfixed(npts,pts,fluxes) BIND(C,NAME="tokamaker_get_vfixed")
INTEGER(c_int), INTENT(out) :: npts
TYPE(c_ptr), INTENT(out) :: pts,fluxes
REAL(8), POINTER :: pts_tmp(:,:),fluxes_tmp(:)
CALL gs_fixed_vflux(gs_global,pts_tmp,fluxes_tmp)
pts=C_LOC(pts_tmp)
fluxes=C_LOC(fluxes_tmp)
npts=SIZE(fluxes_tmp,DIM=1)
END SUBROUTINE tokamaker_get_vfixed
!
SUBROUTINE tokamaker_set_psi(psi_vals) BIND(C,NAME="tokamaker_set_psi")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(psi_vals, vals_tmp, [gs_global%psi%n])
CALL gs_global%psi%restore_local(vals_tmp)
END SUBROUTINE tokamaker_set_psi
!
SUBROUTINE tokamaker_set_psi_dt(psi_vals,dt) BIND(C,NAME="tokamaker_set_psi_dt")
TYPE(c_ptr), VALUE, INTENT(in) :: psi_vals
REAL(c_double), VALUE, INTENT(in) :: dt
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
!
SUBROUTINE tokamaker_set_settings(settings) BIND(C,NAME="tokamaker_set_settings")
TYPE(tokamaker_settings_type), INTENT(in) :: settings
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
!
SUBROUTINE tokamaker_set_targets(ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target) BIND(C,NAME="tokamaker_set_targets")
REAL(c_double), VALUE, INTENT(in) :: ip_target,ip_ratio_target,pax_target,estore_target,R0_target,V0_target
gs_global%R0_target=R0_target
gs_global%V0_target=V0_target
gs_global%pax_target=pax_target*mu0
gs_global%estore_target=estore_target*mu0
gs_global%itor_target=ip_target*mu0
gs_global%ip_ratio_target=ip_ratio_target
END SUBROUTINE tokamaker_set_targets
!
SUBROUTINE tokamaker_set_isoflux(targets,weights,ntargets,grad_wt_lim) BIND(C,NAME="tokamaker_set_isoflux")
REAL(c_double), INTENT(in) :: targets(2,ntargets),weights(ntargets)
INTEGER(c_int), VALUE, INTENT(in) :: ntargets
REAL(c_double), VALUE, INTENT(in) :: grad_wt_lim
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
!
SUBROUTINE tokamaker_set_saddles(targets,weights,ntargets) BIND(C,NAME="tokamaker_set_saddles")
REAL(c_double), INTENT(in) :: targets(2,ntargets),weights(ntargets)
INTEGER(c_int), VALUE, INTENT(in) :: ntargets
IF(ASSOCIATED(gs_global%isoflux_saddles))DEALLOCATE(gs_global%isoflux_saddles)
gs_global%isoflux_nsaddles=ntargets
IF(ntargets>0)THEN
  ALLOCATE(gs_global%isoflux_saddles(3,gs_global%isoflux_nsaddles))
  gs_global%isoflux_saddles(1:2,:)=targets
  gs_global%isoflux_saddles(3,:)=weights
END IF
END SUBROUTINE tokamaker_set_saddles
!
SUBROUTINE tokamaker_set_coil_currents(currents) BIND(C,NAME="tokamaker_set_coil_currents")
TYPE(c_ptr), VALUE, INTENT(in) :: currents
INTEGER(4) :: i
REAL(8) :: curr
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
CALL c_f_pointer(currents, vals_tmp, [gs_global%ncoil_regs])
DO i=1,gs_global%ncoil_regs
  gs_global%coil_regions(i)%curr=vals_tmp(i)*mu0/gs_global%coil_regions(i)%area
END DO
gs_global%vcontrol_val=0.d0
END SUBROUTINE tokamaker_set_coil_currents
!
SUBROUTINE tokamaker_set_coil_regmat(coil_reg_mat,coil_reg_targets,coil_reg_weights) BIND(C,NAME="tokamaker_set_coil_regmat")
TYPE(c_ptr), VALUE, INTENT(in) :: coil_reg_mat,coil_reg_targets,coil_reg_weights
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
IF(.NOT.ASSOCIATED(gs_global%coil_reg_mat))THEN
  ALLOCATE(gs_global%coil_reg_mat(gs_global%ncoil_regs+1,gs_global%ncoil_regs+1))
  ALLOCATE(gs_global%coil_reg_targets(gs_global%ncoil_regs+1))
END IF
CALL c_f_pointer(coil_reg_mat, vals_tmp, [gs_global%ncoil_regs+1,gs_global%ncoil_regs+1])
gs_global%coil_reg_mat=vals_tmp
CALL c_f_pointer(coil_reg_targets, vals_tmp, [gs_global%ncoil_regs+1,1])
gs_global%coil_reg_targets=vals_tmp(:,1)*mu0
CALL c_f_pointer(coil_reg_weights, vals_tmp, [gs_global%ncoil_regs+1,1])
DO i=1,gs_global%ncoil_regs+1
  gs_global%coil_reg_targets(i)=gs_global%coil_reg_targets(i)*vals_tmp(i,1)
  gs_global%coil_reg_mat(i,:)=gs_global%coil_reg_mat(i,:)*vals_tmp(i,1)
END DO
END SUBROUTINE tokamaker_set_coil_regmat
!
SUBROUTINE tokamaker_set_coil_bounds(coil_bounds) BIND(C,NAME="tokamaker_set_coil_bounds")
TYPE(c_ptr), VALUE, INTENT(in) :: coil_bounds
REAL(8), POINTER, DIMENSION(:,:) :: vals_tmp
INTEGER(4) :: i
CALL c_f_pointer(coil_bounds, vals_tmp, [2,gs_global%ncoil_regs+1])
IF(.NOT.ASSOCIATED(gs_global%coil_bounds))THEN
  ALLOCATE(gs_global%coil_bounds(2,gs_global%ncoil_regs+1))
  gs_global%coil_bounds(1,:)=-1.d98; gs_global%coil_bounds(2,:)=1.d98
END IF
DO i=1,gs_global%ncoil_regs
  gs_global%coil_bounds([2,1],i)=-vals_tmp(:,i)*mu0/gs_global%coil_regions(i)%area
END DO
gs_global%coil_bounds([2,1],gs_global%ncoil_regs+1)=-vals_tmp(:,gs_global%ncoil_regs+1)*mu0
END SUBROUTINE tokamaker_set_coil_bounds
!
SUBROUTINE tokamaker_set_coil_vsc(coil_gains) BIND(C,NAME="tokamaker_set_coil_vsc")
TYPE(c_ptr), VALUE, INTENT(in) :: coil_gains
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
INTEGER(4) :: i
CALL c_f_pointer(coil_gains, vals_tmp, [gs_global%ncoil_regs])
DO i=1,gs_global%ncoil_regs
  gs_global%coil_regions(i)%vcont_gain=vals_tmp(i)/gs_global%coil_regions(i)%area
END DO
END SUBROUTINE tokamaker_set_coil_vsc
!
SUBROUTINE tokamaker_save_eqdsk(filename,nr,nz,rbounds,zbounds,run_info,psi_pad) BIND(C,NAME="tokamaker_save_eqdsk")
CHARACTER(KIND=c_char), INTENT(in) :: filename(80),run_info(36)
INTEGER(c_int), VALUE, INTENT(in) :: nr,nz
REAL(c_double), INTENT(in) :: rbounds(2),zbounds(2)
REAL(c_double), VALUE, INTENT(in) :: psi_pad
CHARACTER(LEN=36) :: run_info_f
CHARACTER(LEN=80) :: filename_tmp,lim_file
CALL copy_string_rev(run_info,run_info_f)
CALL copy_string_rev(filename,filename_tmp)
lim_file='none'
CALL gs_save_eqdsk(gs_global,filename_tmp,nr,nz,rbounds,zbounds,run_info_f,lim_file,psi_pad)
END SUBROUTINE tokamaker_save_eqdsk
END MODULE tokamaker_f