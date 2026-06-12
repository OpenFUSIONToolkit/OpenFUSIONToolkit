!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_bootstrap.F90
!
!> Physics-scripts used in the bootstrap calculation inside the Grad-Sharfranov solve
!!
!! @authors Stuart Benjamin (fortranisation of bootstrap.py, Daniel Burgess)
!! @date June 2026
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module grad_shaf_bootstrap
use oft_base
use oft_gs, only: gs_equil, flux_func, gsinv_interp, gs_factory, gs_psi2r, &
 gs_itor_nl
use oft_lag_basis, only: oft_blag_geval
use oft_mesh_type, only: bmesh_findcell
use oft_blag_operators, only: oft_lag_brinterp
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use grad_shaf_prof_phys, only: eval_R_qtmp, build_Ravg_spline, gs_flux_int, &
  jphi_update, jphi_copy, jphi_flux_func
use spline_mod
USE oft_io, ONLY: hdf5_create_group, hdf5_write, hdf5_read, &
  hdf5_field_get_sizes, hdf5_field_exist
use mhd_utils, only: mu0
implicit none
!------------------------------------------------------------------------------
!> Read-in options for the jphi-split-bootstrap current profile update.
!!
!! Must be initialised (via @ref tokamaker_set_boot_ops) before
!! running TokaMaker with a `jphi-split-bootstrap` current profile.  Fields
!! correspond to the optional arguments of @ref calculate_bootstrap in
!! grad_shaf_prof_phys.F90.
!------------------------------------------------------------------------------
TYPE :: boot_ops
  LOGICAL :: initialized = .FALSE. !< Options have been explicitly set?
  LOGICAL :: isolate_edge_jBS = .FALSE. !< Isolate edge bootstrap spike from bulk?
  LOGICAL :: parameterize_jBS = .FALSE. !< Use parametrised skew-normal fit for spike? Overrides `isolate_edge_jBS` if true.
  REAL(r8) :: scale_jBS = 1.0_r8 !< Scaling factor applied to the spike profile (default 1)
  REAL(r8) :: djBS_tol = 1.0e-4_r8 !< Threshold on rel. change in bootstrap current to stop recalculation (increasing solve speed)
  LOGICAL :: diagnose_bs = .FALSE. !< Print alpha/Ip scalars, j_BS stats, and full profile tables each NL iteration
  LOGICAL :: taper_edge_jBS = .TRUE. !< Smooth taper of toroidal current to zero at plasma edge (guards against numerical issues at the separatrix)
  REAL(r8) :: taper_edge_psi0 = 0.999_r8 !< psi_N (standard: 0=axis, 1=LCFS) where edge taper begins
  INTEGER(i4) :: taper_edge_shape = 2 !< Edge taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
END TYPE boot_ops
!------------------------------------------------------------------------------
!> Cached current profiles produced by the last call to @ref jphi_bs_update.
!!
!! All four arrays are allocated/overwritten on every call to jphi_bs_update
!! and remain valid until the next call or until the parent @ref gs_equil is
!! destroyed.  Units are A/m² (before the mu0 normalisation used internally).
!------------------------------------------------------------------------------
TYPE :: boot_profs
  REAL(r8), POINTER, DIMENSION(:) :: psi_n => NULL() !< Normalised psi_N values for these current profiles in OFT convention (0=LCFS, 1=axis); index 0 is the LCFS boundary
  REAL(r8), POINTER, DIMENSION(:) :: j_bs_raw => NULL() !< Raw bootstrap current density output directly from Redl PoP 2021 formula [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: total_j_phi => NULL() !< Total toroidal current density = j_ind_final + j_bs_final [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: j_ind_final => NULL() !< Input jphi, re-scaled & optionally tapered [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: j_bs_final => NULL() !< Bootstrap current density, optionally isolated/parametrised/tapered [A/m²]
END TYPE boot_profs
!------------------------------------------------------------------------------
!> Jphi flux function type for bootstrap current calculation
!------------------------------------------------------------------------------
type, extends(jphi_flux_func) :: jphi_bs_flux_func
  real(8) :: alpha_last = 1.d0 !< Alpha (input Jphi rescaling factor) from previous NL iteration
  real(8) :: rescale_last = 1.d0 !< Damped jphi_rescale from previous NL iteration, to ensure Ip target is met
  logical :: freeze_j_BS = .FALSE. !< Set .TRUE. once j_BS stagnates for 2 steps; skips Sauter call (big speedup)
  logical :: freeze_alpha = .FALSE.   !< Set .TRUE. once dalpha stagnates for 2 steps; skips alpha re-solve (speedup)
  real(8) :: djBS_stol = 1.0e-3_r8 !< RMS tolerance stagnation value: reset no-improve counter if above this
  real(8) :: dalpha_tol(2) = [1.0e-6_r8, 1.0e-3_r8]   !< Hard tolerance; (1) freeze if below this, (2) reset no-improve counter if above this
  integer(4) :: djBS_no_improve = 0   !< Consecutive steps with non-decreasing djBS
  real(8) :: djBS_min = huge(1.0d0)  !< Running minimum djBS seen so far
  integer(4) :: dalpha_no_improve = 0 !< Consecutive steps with non-decreasing dalpha
  real(8) :: dalpha_min = huge(1.0d0) !< Running minimum dalpha seen so far
  real(8), pointer, dimension(:) :: jphi_total_last => NULL() !< Assembled jphi_total from previous NL iteration
  real(8), pointer, dimension(:) :: j_BS_last => NULL() !< j_BS profile from previous NL iteration (for freeze check)
  TYPE(boot_ops) :: boot_ops       !< Python read-in ptions for j_bs_update
  TYPE(boot_profs) :: boot_profs   !< Cached current profiles from the last jphi_bs_update call
contains
  !> Needs docs
  procedure :: save_hdf5 => jphi_bs_save_hdf5
  procedure :: save_txt => jphi_bs_save_txt
  !> Needs docs
  procedure :: load_hdf5 => jphi_bs_load_hdf5
  procedure :: load_txt => jphi_bs_load_txt
  !> Needs docs
  procedure :: copy => jphi_bs_copy
  !> Needs docs
  procedure :: delete => jphi_bs_delete
  !> Update F*F' profile from Jphi and current equilibrium state
  procedure :: update => jphi_bs_update
end type jphi_bs_flux_func
!------------------------------------------------------------------------------
!> Interpolation object for computing Sauter trapped-particle factors
!! (accumulates flux-surface averages of B, R, and bounce integrals)
!------------------------------------------------------------------------------
type, extends(gsinv_interp) :: sauter_interp
  logical :: stage_1 = .FALSE. !< Stage 1 pass (finding Bmax)
  real(8) :: f_surf = 0.d0      !< F surface value (R*Bt)
  real(8) :: bmax = -1.d0       !< Maximum |B| on surface (set in stage 1)
  real(8) :: rmax_surf = -1.d30 !< Maximum R on surface (set in stage 1)
  real(8) :: rmin_surf =  1.d30 !< Minimum R on surface (set in stage 1)
  real(8) :: mag_axis(2) = 0.d0 !< Magnetic axis (R, Z)
contains
  !> Evaluate ODE RHS for Sauter integration
  procedure :: interp => sauter_apply
end type sauter_interp
!------------------------------------------------------------------------------
!> Context type for MINPACK-based fitting of the edge bootstrap profile.
!> Module-level variable `active_edge_jbs` is populated by curve_fit_edge_jbs
!> before invoking lmdif so that edge_jbs_residual can access the target data.
!------------------------------------------------------------------------------
TYPE :: edge_jbs_fit_ctx
  INTEGER(i4) :: n_fit = 0           !< Number of data points in the fit
  REAL(r8), ALLOCATABLE :: psi_fit(:) !< psi_N values at the fitting points
  REAL(r8), ALLOCATABLE :: j_fit(:)   !< Target j_BS values [A/m^2]
  REAL(r8) :: tail_alpha = 1.5_r8    !< Fixed right-side fall-off (Python default)
  REAL(r8) :: lb(7) = 0.0_r8        !< Lower bounds for the 7 parameters
  REAL(r8) :: ub(7) = 1.0_r8        !< Upper bounds for the 7 parameters
END TYPE edge_jbs_fit_ctx
TYPE(edge_jbs_fit_ctx) :: active_edge_jbs !< Module-level context for lmdif callback
contains
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_save_hdf5(self,filename,path)
class(jphi_bs_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
CALL hdf5_write('jphi-split-bootstrap',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%x,filename,path//'/XVALS')
CALL hdf5_write(self%jphi,filename,path//'/YVALS')
CALL hdf5_write(self%j0,filename,path//'/J0')
!---Save boot_ops (only when initialized)
IF(self%boot_ops%initialized)THEN
  CALL hdf5_create_group(filename,path//'/BOOT_OPS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%isolate_edge_jBS),filename,path//'/BOOT_OPS/ISOLATE_EDGE_JBS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%parameterize_jBS),filename,path//'/BOOT_OPS/PARAMETERIZE_JBS')
  CALL hdf5_write(self%boot_ops%scale_jBS,filename,path//'/BOOT_OPS/SCALE_JBS')
  CALL hdf5_write(self%boot_ops%djBS_tol,filename,path//'/BOOT_OPS/DJBS_TOL')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%diagnose_bs),filename,path//'/BOOT_OPS/DIAGNOSE_BS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%taper_edge_jBS),filename,path//'/BOOT_OPS/TAPER_EDGE_JBS')
  CALL hdf5_write(self%boot_ops%taper_edge_psi0,filename,path//'/BOOT_OPS/TAPER_EDGE_PSI0')
  CALL hdf5_write(self%boot_ops%taper_edge_shape,filename,path//'/BOOT_OPS/TAPER_EDGE_SHAPE')
END IF
!---Save cached bootstrap current profiles (only when available)
IF(ASSOCIATED(self%boot_profs%total_j_phi).OR.ASSOCIATED(self%boot_profs%j_bs_raw))THEN
  CALL hdf5_create_group(filename,path//'/BOOT_PROFS')
  IF(ASSOCIATED(self%boot_profs%total_j_phi))THEN
    CALL hdf5_write(self%boot_profs%psi_n,filename,path//'/BOOT_PROFS/PSI_N')
    CALL hdf5_write(self%boot_profs%total_j_phi,filename,path//'/BOOT_PROFS/TOTAL_J_PHI')
    CALL hdf5_write(self%boot_profs%j_ind_final,filename,path//'/BOOT_PROFS/J_IND_FINAL')
    CALL hdf5_write(self%boot_profs%j_bs_final,filename,path//'/BOOT_PROFS/J_BS_FINAL')
    IF(ASSOCIATED(self%boot_profs%j_bs_raw)) &
      CALL hdf5_write(self%boot_profs%j_bs_raw,filename,path//'/BOOT_PROFS/J_BS_RAW')
  END IF
END IF
end subroutine jphi_bs_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_save_txt(self,io_unit)
class(jphi_bs_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'jphi-split-bootstrap'
WRITE(io_unit,*)self%npsi,self%j0
WRITE(io_unit,*)self%x
WRITE(io_unit,*)self%jphi
end subroutine jphi_bs_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_load_hdf5(self,filename,path,success)
class(jphi_bs_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4), allocatable :: dim_sizes(:)
integer(i4) :: npsi
integer(4) :: int_tmp, ndims
real(r8) :: J0
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
ALLOCATE(xvals(npsi),yvals(npsi))
CALL hdf5_read(xvals,filename,path//'/XVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(yvals,filename,path//'/YVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(J0,filename,path//'/J0',success=success)
IF(.NOT.success)RETURN
CALL create_jphi_bs_ff(self,npsi,xvals,yvals,J0) ! Load npsi,xvals,yvals,J0
DEALLOCATE(xvals,yvals)
!---Load boot_ops
self%boot_ops%initialized = .FALSE.
IF(hdf5_field_exist(filename,path//'/BOOT_OPS'))THEN
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/ISOLATE_EDGE_JBS',success=success)
  IF(success) self%boot_ops%isolate_edge_jBS = (int_tmp/=0)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/PARAMETERIZE_JBS',success=success)
  IF(success) self%boot_ops%parameterize_jBS = (int_tmp/=0)
  CALL hdf5_read(self%boot_ops%scale_jBS,filename,path//'/BOOT_OPS/SCALE_JBS',success=success)
  CALL hdf5_read(self%boot_ops%djBS_tol,filename,path//'/BOOT_OPS/DJBS_TOL',success=success)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/DIAGNOSE_BS',success=success)
  IF(success) self%boot_ops%diagnose_bs = (int_tmp/=0)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/TAPER_EDGE_JBS',success=success)
  IF(success) self%boot_ops%taper_edge_jBS = (int_tmp/=0)
  CALL hdf5_read(self%boot_ops%taper_edge_psi0,filename,path//'/BOOT_OPS/TAPER_EDGE_PSI0',success=success)
  CALL hdf5_read(self%boot_ops%taper_edge_shape,filename,path//'/BOOT_OPS/TAPER_EDGE_SHAPE',success=success)
  self%boot_ops%initialized = .TRUE.
END IF
!---Load cached bootstrap current profiles if available
IF(hdf5_field_exist(filename,path//'/BOOT_PROFS'))THEN
  IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/TOTAL_J_PHI'))THEN
    CALL hdf5_field_get_sizes(filename,path//'/BOOT_PROFS/TOTAL_J_PHI',ndims,dim_sizes)
    IF(ASSOCIATED(self%boot_profs%psi_n))DEALLOCATE(self%boot_profs%psi_n)
    ALLOCATE(self%boot_profs%psi_n(0:dim_sizes(1)-1))
    IF(ASSOCIATED(self%boot_profs%total_j_phi))DEALLOCATE(self%boot_profs%total_j_phi)
    ALLOCATE(self%boot_profs%total_j_phi(0:dim_sizes(1)-1))
    IF(ASSOCIATED(self%boot_profs%j_ind_final))DEALLOCATE(self%boot_profs%j_ind_final)
    ALLOCATE(self%boot_profs%j_ind_final(0:dim_sizes(1)-1))
    IF(ASSOCIATED(self%boot_profs%j_bs_final))DEALLOCATE(self%boot_profs%j_bs_final)
    ALLOCATE(self%boot_profs%j_bs_final(0:dim_sizes(1)-1))
    DEALLOCATE(dim_sizes)
    CALL hdf5_read(self%boot_profs%psi_n,filename,path//'/BOOT_PROFS/PSI_N',success=success)
    CALL hdf5_read(self%boot_profs%total_j_phi,filename,path//'/BOOT_PROFS/TOTAL_J_PHI',success=success)
    CALL hdf5_read(self%boot_profs%j_ind_final,filename,path//'/BOOT_PROFS/J_IND_FINAL',success=success)
    CALL hdf5_read(self%boot_profs%j_bs_final,filename,path//'/BOOT_PROFS/J_BS_FINAL',success=success)
    IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/J_BS_RAW'))THEN
      CALL hdf5_field_get_sizes(filename,path//'/BOOT_PROFS/J_BS_RAW',ndims,dim_sizes)
      IF(ASSOCIATED(self%boot_profs%j_bs_raw))DEALLOCATE(self%boot_profs%j_bs_raw)
      ALLOCATE(self%boot_profs%j_bs_raw(0:dim_sizes(1)-1))
      DEALLOCATE(dim_sizes)
      CALL hdf5_read(self%boot_profs%j_bs_raw,filename,path//'/BOOT_PROFS/J_BS_RAW',success=success)
    END IF
  END IF
END IF
end subroutine jphi_bs_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_load_txt(self,io_unit)
class(jphi_bs_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
real(r8) :: J0
real(r8), allocatable :: xvals(:),yvals(:)
READ(io_unit,*)npsi,J0
ALLOCATE(xvals(npsi),yvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)yvals
CALL create_jphi_bs_ff(self,npsi,xvals,yvals,J0)
DEALLOCATE(xvals,yvals)
end subroutine jphi_bs_load_txt
!------------------------------------------------------------------------------
!> Sets values inside jphi_bs_flux_func. Called only after a fresh allocation.
!------------------------------------------------------------------------------
SUBROUTINE create_jphi_bs_ff(func,npsi,psivals,yvals,y0)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
REAL(8), INTENT(in) :: y0
INTEGER(4) :: i,ierr
SELECT TYPE(self=>func)
  TYPE IS(jphi_bs_flux_func)
  !---
  self%npsi=npsi
  self%ndofs=self%npsi
  !---
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  ALLOCATE(self%jphi(self%npsi))
  !---
  self%j0=y0
  self%y0=0.d0
  self%update_on_load = .FALSE. ! Don't update on load to prevent missing kinetic profile errors
  DO i=1,self%npsi
    self%x(i) = psivals(i)
    self%jphi(i) = yvals(i)
    self%yp(i) = psivals(i) ! Dummy initialization
  END DO
  self%yp = self%yp/(SUM(ABS(self%yp))/REAL(self%npsi,8)) ! Consistent (hopefully) normalization
  ierr=self%set_cofs(self%yp)
  IF(oft_debug_print(1))WRITE(*,*)'Jphi bs flux func Created',self%ndofs,self%x,self%j0
class default
  CALL oft_abort('Invalid flux function type in create_jphi_bs_ff','create_jphi_bs_ff',__FILE__)
END SELECT

END SUBROUTINE create_jphi_bs_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_copy(self,new)
class(jphi_bs_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
CALL jphi_copy(self,new)
SELECT TYPE(new)
  CLASS IS(jphi_bs_flux_func)
    new%alpha_last = self%alpha_last
    new%rescale_last = self%rescale_last
    new%freeze_j_BS = self%freeze_j_BS
    new%freeze_alpha = self%freeze_alpha
    new%djBS_stol = self%djBS_stol
    new%dalpha_tol = self%dalpha_tol
    new%djBS_no_improve = self%djBS_no_improve
    new%djBS_min = self%djBS_min
    new%dalpha_no_improve = self%dalpha_no_improve
    new%dalpha_min = self%dalpha_min
    new%boot_ops = self%boot_ops
    IF(ASSOCIATED(self%jphi_total_last)) ALLOCATE(new%jphi_total_last, SOURCE=self%jphi_total_last)
    IF(ASSOCIATED(self%j_BS_last)) ALLOCATE(new%j_BS_last, SOURCE=self%j_BS_last)
    IF(ASSOCIATED(self%boot_profs%psi_n))ALLOCATE(new%boot_profs%psi_n,SOURCE=self%boot_profs%psi_n)
    IF(ASSOCIATED(self%boot_profs%j_bs_raw))ALLOCATE(new%boot_profs%j_bs_raw,SOURCE=self%boot_profs%j_bs_raw)
    IF(ASSOCIATED(self%boot_profs%total_j_phi))ALLOCATE(new%boot_profs%total_j_phi,SOURCE=self%boot_profs%total_j_phi)
    IF(ASSOCIATED(self%boot_profs%j_bs_final))ALLOCATE(new%boot_profs%j_bs_final,SOURCE=self%boot_profs%j_bs_final)
    IF(ASSOCIATED(self%boot_profs%j_ind_final))ALLOCATE(new%boot_profs%j_ind_final,SOURCE=self%boot_profs%j_ind_final)
END SELECT
end subroutine jphi_bs_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_bs_delete(self)
class(jphi_bs_flux_func), intent(inout) :: self
self%j0=0.d0
IF(ASSOCIATED(self%jphi))DEALLOCATE(self%jphi)
IF(ASSOCIATED(self%x))DEALLOCATE(self%x)
IF(ASSOCIATED(self%yp))DEALLOCATE(self%yp)
IF(ASSOCIATED(self%y))DEALLOCATE(self%y)
IF(ASSOCIATED(self%jphi_total_last))DEALLOCATE(self%jphi_total_last)
IF(ASSOCIATED(self%j_BS_last))DEALLOCATE(self%j_BS_last)
!---Destroy cached bootstrap current profiles
IF(ASSOCIATED(self%boot_profs%j_bs_raw))DEALLOCATE(self%boot_profs%j_bs_raw)
IF(ASSOCIATED(self%boot_profs%total_j_phi))DEALLOCATE(self%boot_profs%total_j_phi)
IF(ASSOCIATED(self%boot_profs%j_bs_final))DEALLOCATE(self%boot_profs%j_bs_final)
IF(ASSOCIATED(self%boot_profs%j_ind_final))DEALLOCATE(self%boot_profs%j_ind_final)
IF(ASSOCIATED(self%boot_profs%psi_n))DEALLOCATE(self%boot_profs%psi_n)
end subroutine jphi_bs_delete
!---------------------------------------------------------------------------------
!> Update F*F' profile from inductive Jphi coupled with bootstrap current.
!>
!> Each call (one NL iteration):
!>   1. Build <R>/<1/R> spline on self%x; pre-compute qtmp = <R>*<1/R>.
!>   2. Evaluate bootstrap current j_BS on self%x (or reuse cache if frozen).
!>   3. Apply edge taper to j_BS and jphi_ind component-wise.
!>   4. Compute jphi_rescale to reconcile gs_itor_nl vs gs_flux_int.
!>   5. Solve analytically for alpha: gs_flux_int is linear in alpha, so two
!>      evaluations (alpha=0, alpha=1) give alpha = (Ip_target - Ip_lo)/(Ip_hi - Ip_lo).
!>   6. Assemble jphi_total = alpha*jphi_ind + j_BS; compute F*F' knots.
!>   7. Diagnostics (if diagnose_bs is set).
!---------------------------------------------------------------------------------
SUBROUTINE jphi_bs_update(self, gseq)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: pscale, pprime
REAL(r8), ALLOCATABLE :: qtmp(:)
TYPE(spline_type) :: R_spline
! Bootstrap arrays (on self%x grid)
REAL(r8), ALLOCATABLE :: j_BS(:)
! Optional edge-spike workspace (only allocated when isolate_edge_jBS or parameterize_jBS is set)
REAL(r8), ALLOCATABLE :: j_spike_tmp(:)
REAL(r8), ALLOCATABLE :: j_spike_mask_tmp(:)  !< Raw masked (pre-fit) spike profile
! Working arrays on self%x grid
REAL(r8), ALLOCATABLE :: jphi_total(:)
REAL(r8), ALLOCATABLE :: jphi_ind(:)  !< Tapered copy of self%jphi (= self%jphi when taper off)
! Alpha-solve scalars
REAL(r8) :: alpha, ip_target, ip_ind, ip_result_lo, ip_result_hi, dalpha
! Relative change in bootstrap current for freeze check
REAL(r8) :: djBS
! gs_itor_nl / gs_flux_int reconciliation
REAL(r8) :: itor_nl = 0.0_r8, itor_flint = 0.0_r8, jphi_rescale
CHARACTER(len=256) :: char_buf
!--- First-iteration runs with no bootstrap current
IF(.NOT. gseq%skip_targets) THEN
  CALL jphi_update(self,gseq)
  RETURN
ENDIF
!---
self%plasma_bounds = gseq%plasma_bounds
IF(gseq%mode/=1) &
  CALL oft_abort("Jphi-BS profile requires (F^2)' formulation", &
                 "jphi_bs_update",__FILE__)
IF(gseq%pax_target<0.d0) &
  CALL oft_abort("Jphi-BS profile requires Pax target", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%Te)) &
  CALL oft_abort("Jphi-BS profile requires Te profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%Ti)) &
  CALL oft_abort("Jphi-BS profile requires Ti profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%ne)) &
  CALL oft_abort("Jphi-BS profile requires ne profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%ni)) &
  CALL oft_abort("Jphi-BS profile requires ni profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%Zeff)) &
  CALL oft_abort("Jphi-BS profile requires Zeff profile", &
                 "jphi_bs_update",__FILE__)
!--- 1. Build <R>/<1/R> spline; pre-compute qtmp = <R>*<1/R> on self%x.
!   R_spline stays alive until after the F*F' loop (step 6).
ALLOCATE(qtmp(self%npsi))
CALL build_Ravg_spline(gseq, self%ngeom, R_spline)
CALL eval_R_qtmp(R_spline, self%x, self%npsi, qtmp)
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
!--- 2. Bootstrap current on self%x grid.
ALLOCATE(j_BS(self%npsi))
IF(self%freeze_j_BS .AND. ASSOCIATED(self%j_BS_last)) THEN
  !--- Frozen: reuse cached j_BS.
  j_BS = self%j_BS_last
  djBS = 0.0_r8
ELSE
  !--- Not frozen: run full bootstrap calculation (Sauter).
  IF (self%boot_ops%isolate_edge_jBS .OR. self%boot_ops%parameterize_jBS) THEN
    ALLOCATE(j_spike_tmp(self%npsi), j_spike_mask_tmp(self%npsi))
    CALL calculate_bootstrap(self, gseq, self%npsi, self%x, j_BS, &
        isolate_edge_jBS=self%boot_ops%isolate_edge_jBS, &
        parameterize_jBS=self%boot_ops%parameterize_jBS, &
        scale_jBS=self%boot_ops%scale_jBS, &
        j_spike=j_spike_tmp, j_spike_masked=j_spike_mask_tmp)
    IF (self%boot_ops%diagnose_bs) THEN
      IF (self%boot_ops%parameterize_jBS) THEN
        WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS(bulk)[A/m2]  j_spike[A/m2]   j_spike_masked[A/m2]  jphi[A/m2]'
        DO i = 1, self%npsi
          WRITE(*,'(A,I4,5ES15.5)') '  ', i, self%x(i), j_BS(i), j_spike_tmp(i), j_spike_mask_tmp(i), self%jphi(i)
        END DO
      ELSE
        WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS(bulk)[A/m2]  j_spike[A/m2]   jphi[A/m2]'
        DO i = 1, self%npsi
          WRITE(*,'(A,I4,4ES15.5)') '  ', i, self%x(i), j_BS(i), j_spike_tmp(i), self%jphi(i)
        END DO
      END IF
    END IF
    j_BS = j_spike_tmp
    DEALLOCATE(j_spike_tmp, j_spike_mask_tmp)
  ELSE
    CALL calculate_bootstrap(self, gseq, self%npsi, self%x, j_BS)
    j_BS = j_BS * self%boot_ops%scale_jBS
    IF(self%boot_ops%diagnose_bs)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS[A/m2]      jphi[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,3ES15.5)') '  ', i, self%x(i), j_BS(i), self%jphi(i)
      END DO
    END IF
  END IF
  !   calculate_bootstrap returns j_BS in raw A/m², multiply by mu0.
  j_BS = j_BS * mu0
  !--- 2a. Freeze check: freeze j_BS if RMS change drops below tol, or stagnates.
  IF(ASSOCIATED(self%j_BS_last)) THEN
    djBS = SQRT(SUM((j_BS - self%j_BS_last)**2) / REAL(self%npsi,r8)) / &
           MAX(SQRT(SUM(j_BS**2) / REAL(self%npsi,r8)), 1.0e-30_r8)
    IF(djBS < self%boot_ops%djBS_tol .OR. self%freeze_alpha) THEN
      self%freeze_j_BS = .TRUE.
      IF(oft_env%pm)WRITE(*,*)' Freezing bootstrap solution.'
    ELSE IF(djBS >= self%djBS_min) THEN
      self%djBS_no_improve = self%djBS_no_improve + 1
      IF(self%djBS_no_improve >= 2) THEN
        WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
          'Bootstrap solution convergence stalled,' // &
          ' relative change per nonlinear step = ', djBS, &
          ', above set tolerance (djBS_tol=', self%boot_ops%djBS_tol, ')'
        IF(djBS > self%djBS_stol) THEN
          self%djBS_no_improve = 0 ! Reset counter, give more chances to improve
        ELSE
          self%freeze_j_BS = .TRUE.
          char_buf = TRIM(char_buf) // ' Freezing bootstrap solution.'
        END IF
        IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
      END IF
    ELSE
      self%djBS_no_improve = 0
      self%djBS_min = djBS
    END IF
  ELSE
    djBS = HUGE(1.0_r8)
  END IF
  IF(.NOT.ASSOCIATED(self%j_BS_last)) ALLOCATE(self%j_BS_last(self%npsi))
  self%j_BS_last = j_BS
END IF
!--- 3. Apply edge taper to j_BS and jphi_ind.
!   self%j_BS_last caches the un-tapered j_BS so freeze comparisons track physics.
!   taper_edge_psi0 is in standard convention (0=axis,1=LCFS);
!   threshold in OFT convention (0=LCFS,1=axis) is (1 - taper_edge_psi0).
ALLOCATE(jphi_ind(self%npsi))
jphi_ind = self%jphi
IF (self%boot_ops%taper_edge_jBS) THEN
  CALL apply_edge_taper(self%npsi, self%x, j_BS, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
  CALL apply_edge_taper(self%npsi, self%x, jphi_ind, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
END IF
ALLOCATE(jphi_total(self%npsi))
!--- 4. Reconcile gs_itor_nl vs gs_flux_int.
!   No Ip target: rescale jphi_total so the integrated current matches the
!   FEM solution (gs_itor_nl) rather than the profile quadrature (gs_flux_int).
jphi_rescale = self%rescale_last
IF(ASSOCIATED(self%jphi_total_last) .AND. .NOT. self%freeze_alpha) THEN
  CALL gs_itor_nl(gseq, itor_nl)
  CALL gs_flux_int(gseq, self%x, self%jphi_total_last/qtmp, self%npsi, itor_flint)
  jphi_rescale = (itor_nl/itor_flint + self%rescale_last) / 2.0_r8
  self%rescale_last = jphi_rescale
END IF
!--- 5. Solve analytically for alpha.
!   gs_flux_int is linear in alpha; two evaluations (alpha=0 and alpha=1) give
!   alpha = (Ip_target - Ip_lo) / (Ip_hi - Ip_lo).  Skip once frozen.
ip_target = ABS(gseq%Ip_target)/jphi_rescale
IF(self%freeze_alpha) THEN
  !--- Frozen: reuse last converged alpha.
  alpha = self%alpha_last
  dalpha = 0.0_r8
  ip_result_lo = 0.0_r8
  ip_result_hi = 0.0_r8
ELSE
  !--- Not yet frozen: exact linear solve for alpha.
  jphi_total = j_BS
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, ip_result_lo)
  jphi_total = jphi_ind + j_BS
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, ip_result_hi)
  ip_ind = ip_result_hi - ip_result_lo
  IF(ABS(ip_ind) > 0.0_r8)THEN
    alpha = (ip_target - ip_result_lo) / ip_ind
    IF(alpha < -1.0_r8 .OR. alpha > 10.0_r8) THEN
      WRITE(char_buf,'(A,ES12.4)') '[jphi_bs_update] WARNING: alpha out of expected range [-1,10]: alpha=', alpha
      CALL oft_warn(TRIM(char_buf))
    END IF
  ELSE
    alpha = self%alpha_last
  END IF
  dalpha = ABS(alpha - self%alpha_last)
  !--- Freeze alpha on hard tol OR 2 consecutive steps above running minimum; also freeze bootstrap.
  IF(dalpha < self%dalpha_tol(1)) THEN
    self%freeze_alpha = .TRUE.
    self%freeze_j_BS = .TRUE.
    IF(oft_env%pm)WRITE(*,*)' Freezing inductive current scale.'
  ELSE IF(dalpha >= self%dalpha_min) THEN
    self%dalpha_no_improve = self%dalpha_no_improve + 1
    IF(self%dalpha_no_improve >= 2) THEN
      WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
        'Alpha convergence stalled,' // &
        ' relative change per nonlinear step = ', dalpha, &
        ', above recommended tolerance (dalpha_tol=', self%dalpha_tol(1), ')'
      IF(dalpha > self%dalpha_tol(2)) THEN
        self%dalpha_no_improve = 0 ! Reset counter, give more chances to improve
      ELSE
        self%freeze_alpha = .TRUE.
        self%freeze_j_BS = .TRUE.
        char_buf = TRIM(char_buf) // ' Freezing alpha solution.'
      END IF
      IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
    END IF
  ELSE
    self%dalpha_no_improve = 0
    self%dalpha_min = dalpha
  END IF
END IF
self%alpha_last = alpha
!--- 6. Assemble jphi_total, save profiles
jphi_total = alpha * jphi_ind + j_BS
IF(.NOT.ASSOCIATED(self%jphi_total_last)) ALLOCATE(self%jphi_total_last(self%npsi))
self%jphi_total_last = jphi_total
IF(.NOT.ASSOCIATED(self%boot_profs%total_j_phi))THEN
  ALLOCATE(self%boot_profs%psi_n(0:self%npsi))
  ALLOCATE(self%boot_profs%total_j_phi(0:self%npsi))
  ALLOCATE(self%boot_profs%j_bs_final(0:self%npsi))
  ALLOCATE(self%boot_profs%j_ind_final(0:self%npsi))
END IF
! LCFS boundary (OFT psi=0; self%j0 is jphi_ind at LCFS; j_BS=0 at LCFS)
self%boot_profs%psi_n(0)       = 0.0_r8
self%boot_profs%total_j_phi(0) = alpha*self%j0/mu0
self%boot_profs%j_bs_final(0)  = 0.0_r8
self%boot_profs%j_ind_final(0) = alpha*self%j0/mu0
! Interior knots (OFT psi convention: self%x(1) near LCFS, self%x(npsi) near axis)
self%boot_profs%psi_n(1:)       = self%x
self%boot_profs%total_j_phi(1:) = jphi_total/mu0
self%boot_profs%j_bs_final(1:)  = j_BS/mu0
self%boot_profs%j_ind_final(1:) = alpha * jphi_ind/mu0
!--- Compute updated F*F' profile
IF(ASSOCIATED(gseq%P_ani)) &
  CALL oft_abort('Jphi profiles do not support anisotropic pressure', &
                 'jphi_bs_update',__FILE__)
pscale = gseq%P%f(gseq%plasma_bounds(2))
pscale = gseq%pax_target / pscale
CALL spline_eval(R_spline, 0.d0, 0) ! LCFS point for y0 calculation
pprime = gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(self%j0*alpha - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
DO i = 1, self%npsi
  CALL spline_eval(R_spline, self%x(i), 0)
  pprime = gseq%P%fp(self%x(i)*(gseq%plasma_bounds(2) - &
                                  gseq%plasma_bounds(1)) + &
                                  gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(jphi_total(i) - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
END DO
! Fix F*F' scale (matching is done here instead)
! gseq%skip_targets is already true when jphi_bs_update called
gseq%ffp_scale=1.d0
gseq%p_scale=pscale
!--- 7. Diagnostics.
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_target   = ', ip_target
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_lo= ', ip_result_lo
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_hi= ', ip_result_hi
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] alpha       = ', alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] dalpha      = ', dalpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djBS        = ', djBS
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] bs_frozen   = ', self%freeze_j_BS
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] freeze_alpha= ', self%freeze_alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] j_BS max    = ', MAXVAL(ABS(j_BS))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi max    = ', MAXVAL(ABS(self%jphi))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_rescale= ', jphi_rescale
  !--- Side-by-side Ip comparison: FEM nonlinear solve vs profile flux integral
  CALL gs_itor_nl(gseq, itor_nl)
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, itor_flint)
  WRITE(*,'(A)') '  [jphi_bs_update] --- Ip comparison (current jphi_total) ---'
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(gs_itor_nl)    = ', itor_nl/mu0
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(flux_int/qtmp) = ', itor_flint/mu0
END IF
!--- Clean up
DEALLOCATE(j_BS, jphi_total, jphi_ind, qtmp)
CALL spline_dealloc(R_spline)
i=self%set_cofs(self%yp)
END SUBROUTINE jphi_bs_update
!------------------------------------------------------------------------------
!> Evaluate terms in augmented tracing ODE for computing Sauter factors
!! (see @ref sauter_fc)
!------------------------------------------------------------------------------
subroutine sauter_apply(self,cell,f,gop,val)
class(sauter_interp), intent(inout) :: self !< Interpolation object
integer(4), intent(in) :: cell !< Cell for interpolation
real(8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(8), intent(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
real(8), intent(out) :: val(:) !< Reconstructed field at f [8]
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),pt(3),grad(3)
real(8) :: s,c,Bp2,Bt2,mod_B,bratio
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j)
!---Reconstruct gradient
grad=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  grad=grad+self%uvals(j(jc))*rop
end do
!---Get radial position
pt=self%mesh%log2phys(cell,f)
s=SIN(self%t)
c=COS(self%t)
Bp2 = (grad(1)**2 + grad(2)**2)/pt(1)**2
Bt2 = (self%f_surf/pt(1))**2
mod_B = SQRT(Bp2+Bt2)
!---Position
val(1)=(self%rho*(grad(1)*s-grad(2)*c))/(grad(1)*c+grad(2)*s)
val(2)=pt(1)*SQRT((self%rho**2+val(1)**2)/SUM(grad**2))
!---Geometric factors
val(3)=val(2)*pt(1)         ! <R>
val(4)=val(2)/pt(1)         ! <1/R>
val(5)=val(2)*SQRT((pt(1)-self%mag_axis(1))**2+(pt(2)-self%mag_axis(2))**2)  ! <a>
!---Magnetic field averages
val(6)=val(2)*SQRT(Bp2+Bt2) ! <|B|>
val(7)=val(2)*(Bp2+Bt2)     ! <|B|^2>
IF(self%stage_1)THEN
  self%bmax = MAX(self%bmax,mod_B)
  self%rmax_surf = MAX(self%rmax_surf, pt(1))
  self%rmin_surf = MIN(self%rmin_surf, pt(1))
  val(8) = 0.d0  ! not used in stage 1; zero to prevent NaN accumulation
  val(9) = 0.d0  ! likewise
ELSE
  bratio = MIN(1.d0,mod_B/self%bmax)
  val(8)=val(2)*((1.d0 - SQRT(1.d0 - bratio) * (1.d0 + bratio / 2.d0)) / bratio**2)
  val(9) = val(2)/pt(1)**2   ! q integrand: q = F_surf/(2π) * ∮ val(9)
END IF
deallocate(j)
end subroutine sauter_apply
!------------------------------------------------------------------------------
!> Compute factors required for Sauter bootstrap formula
!------------------------------------------------------------------------------
subroutine sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs,qprof,eps)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(4), intent(in) :: nr !< Number of flux sample points
real(8), intent(in) :: psi_q(nr) !< Location of flux sample points in normalised psi
real(8), intent(out) :: fc(nr) !< Circulating particle fraction \f$ f_c \f$
real(8), intent(out) :: r_avgs(nr,3) !< Flux surface averaged radial coords \f$<R>\f$, \f$<1/R>\f$, \f$<a>\f$
real(8), intent(out) :: modb_avgs(nr,2) !< Flux surface averaged field \f$<|B|>\f$, \f$<|B|^2>\f$
real(8), optional, intent(out) :: qprof(nr) !< Safety factor q on each surface (avoids a separate gs_get_qprof call)
real(8), optional, intent(out) :: eps(nr) !< Local inverse aspect ratio \f$ \varepsilon = (R_{\max}-R_{\min})/(2\langle R \rangle) \f$ on each surface
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,fpol,qpsi,h,h2,hf,ftu,ftl
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
character(len=80) :: warn_str
logical :: edge_trace_failed
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,cell
type(sauter_interp), target :: field
type(gs_factory), pointer :: device
device=>gseq%device
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(device%fe_rep)
!---Find Rmax along Z=zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(device%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_tmp)
  IF( psi_tmp(1) < x1)EXIT
  rmax=pt(1)
END DO
pt_last=[(.1d0*rmax+.9d0*raxis),zaxis,0.d0]
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Axis Position:'
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
  CALL oft_decrease_indent
END IF
call set_tracer(1)
field%u=>gseq%psi
field%mag_axis=gseq%o_point
CALL field%setup(device%fe_rep)
active_tracer%neq=9
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
edge_trace_failed = .FALSE.
do j=1,nr
  psi_surf=psi_q(j)*(x2-x1) + x1
  !--- Axis guard: flux surface degenerates at psi_N=1; prescribe analytically.
  IF(psi_q(j) == 1.d0)THEN
    IF(gseq%mode==0)THEN
      fpol=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
    ELSE
      fpol=SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
        + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
    END IF
    r_avgs(j,1)    = raxis
    r_avgs(j,2)    = 1.d0 / raxis
    r_avgs(j,3)    = 0.d0
    modb_avgs(j,1) = ABS(fpol) / raxis
    modb_avgs(j,2) = (fpol / raxis)**2
    fc(j)          = 1.d0
    IF(PRESENT(qprof))   qprof(j)   = MERGE(qprof(j-1), 0.d0, j > 1)
    IF(PRESENT(eps))     eps(j)     = 0.d0
    CYCLE
  END IF
  IF(gseq%diverted.AND.psi_q(j)<0.02d0)THEN
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  pt=pt_last
  CALL gs_psi2r(gseq,psi_surf,pt,psi_int=psi_int)
  IF(gseq%mode==0)THEN
    fpol=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    fpol=SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  field%f_surf=fpol
  field%bmax=0.d0
  field%rmax_surf = -1.d30
  field%rmin_surf =  1.d30
  field%stage_1=.TRUE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  field%stage_1=.FALSE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  pt_last=pt
  if(active_tracer%status/=1)THEN
    WRITE(*,*)j,pt
    CALL oft_warn("sauter_fc: Trace did not complete")
    IF(j==1) edge_trace_failed = .TRUE.
    CYCLE
  end if
  r_avgs(j,1)=active_tracer%v(3)/active_tracer%v(2)
  r_avgs(j,2)=active_tracer%v(4)/active_tracer%v(2)
  r_avgs(j,3)=active_tracer%v(5)/active_tracer%v(2)
  modb_avgs(j,1)=active_tracer%v(6)/active_tracer%v(2)
  modb_avgs(j,2)=active_tracer%v(7)/active_tracer%v(2)
  h = modb_avgs(j,1)/field%bmax
  h2 = modb_avgs(j,2)/field%bmax**2
  ftu = 1.d0 - h2 / (h**2) * (1.d0 - SQRT(1.d0 - h) * (1.d0 + 0.5d0 * h))
  hf = active_tracer%v(8)/active_tracer%v(2)
  ftl = 1.d0 - h2 * hf
  fc(j) = 1.d0 - (0.75d0 * ftu + 0.25d0 * ftl)
  qpsi = fpol * active_tracer%v(9) / (2*pi)
  IF(PRESENT(qprof)) qprof(j) = qpsi
  IF(PRESENT(eps))   eps(j)   = (field%rmax_surf - field%rmin_surf) / (2.d0 * r_avgs(j,1))
end do
! Edge guard: if psi_q(1) is very close to the separatrix (psi_N < 1e-5) and
! the field-line trace failed there, copy the geometry from the second surface
! as the best available approximation.
IF (nr >= 2 .AND. psi_q(1) < 1.d-5 .AND. edge_trace_failed) THEN
  fc(1)          = fc(2)
  r_avgs(1,:)    = r_avgs(2,:)
  modb_avgs(1,:) = modb_avgs(2,:)
  IF (PRESENT(eps))     eps(1)     = eps(2)
  IF (PRESENT(qprof))   qprof(1)   = 2*qprof(2) ! <this is totally wrong (q diverges at the separatrix), but better than leaving q=0 at the edge for now
  IF (PRESENT(qprof)) THEN
    WRITE(warn_str,'(A,ES10.3)') 'sauter_fc: edge geometry extrapolated from second surface; q unreliable beyond psi_N = ',psi_q(2)
    CALL oft_warn(TRIM(warn_str))
  END IF
END IF
CALL field%delete
CALL active_tracer%delete
CALL psi_int%delete()
end subroutine sauter_fc
!------------------------------------------------------------------------------
!> Apply a smooth edge taper to an array, zeroing it at the plasma edge.
!>
!> @param n          Number of grid points
!> @param psi        Normalised psi grid
!> @param arr        Array to taper (modified in-place)
!> @param psi0       psi value (in the grid's own convention) where taper begins
!> @param shape      Taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
!> @param oft_psi_conv  If .TRUE., psi uses OFT internal convention (0=LCFS, 1=axis)
!>                      and psi0 should be the OFT threshold (= 1 - standard psi0).
!>                      Taper then acts on points where psi <= psi0.
!>                      If .FALSE. (default), standard convention (0=axis, 1=LCFS),
!>                      taper acts on points where psi >= psi0.
!------------------------------------------------------------------------------
SUBROUTINE apply_edge_taper(n, psi, arr, psi0, shape, oft_psi_conv)
INTEGER(i4), INTENT(in)    :: n
REAL(r8),    INTENT(in)    :: psi(n)
REAL(r8),    INTENT(inout) :: arr(n)
REAL(r8),    INTENT(in)    :: psi0
INTEGER(i4), INTENT(in)    :: shape
LOGICAL,     OPTIONAL, INTENT(in) :: oft_psi_conv
!---
INTEGER(i4) :: i
LOGICAL     :: do_oft
REAL(r8)    :: t_taper, w_taper, span
CHARACTER(len=80) :: char_buf
REAL(r8), PARAMETER :: HALF_PI = 1.5707963267948966_r8
do_oft = .FALSE.
IF (PRESENT(oft_psi_conv)) do_oft = oft_psi_conv
! Taper width: distance from threshold to the plasma edge.
! OFT convention: edge at psi=0, so width = psi0.
! Standard convention: edge at psi=1, so width = 1 - psi0.
IF (do_oft) THEN
  span = psi0
ELSE
  span = 1.0_r8 - psi0
END IF
IF (span < 1.0e-6_r8) RETURN
DO i = 1, n
  IF (do_oft) THEN
    ! OFT convention: edge at psi=0, taper region is psi in [0, psi0]
    IF (psi(i) > psi0) CYCLE
    t_taper = (psi0 - psi(i)) / span
  ELSE
    ! Standard convention: edge at psi=1, taper region is psi in [psi0, 1]
    IF (psi(i) < psi0) CYCLE
    t_taper = (psi(i) - psi0) / span
  END IF
  t_taper = MIN(MAX(t_taper, 0.0_r8), 1.0_r8)
  SELECT CASE (shape)
    CASE (1)  ! cos² / Hann
      w_taper = COS(HALF_PI * t_taper)**2
    CASE (2)  ! quintic smoothstep
      w_taper = 1.0_r8 - t_taper**3*(6.0_r8*t_taper**2 - 15.0_r8*t_taper + 10.0_r8)
    CASE (3)  ! cubic power decay
      w_taper = (1.0_r8 - t_taper)**3
    CASE DEFAULT
      WRITE(char_buf,'(A,I0,A)') 'apply_edge_taper: unknown taper shape=', shape, '; no taper applied'
      CALL oft_warn(TRIM(char_buf))
      w_taper = 1.0_r8
  END SELECT
  arr(i) = arr(i) * w_taper
END DO
END SUBROUTINE apply_edge_taper
!------------------------------------------------------------------------------
!> Computes the bootstrap current on a uniform psi_N grid.
!>
!> Translated from Python bootstrap.py: calculate_profiles_and_bootstrap
!> (up to and including the redl_bootstrap call).
!>
!> Fixed conventions: NRL Coulomb logs, Koh ion collisionality model,
!> Redl 2021 jboot1 form with use_sign_q=.TRUE., L34=L31.
!>
!> @param gseq    Equilibrium object (must have Te, Ti, ne, ni, Zeff set)
!> @param n_psi    Number of flux surface samples
!> @param psi_N    Normalised psi grid [0,1], arbitrary spacing
!> @param j_BS Output: average toroidal bootstrap current density [A/m^2] on psi_N grid
!------------------------------------------------------------------------------
SUBROUTINE calculate_bootstrap(gseq, n_psi, psi_N, j_BS, &
                               isolate_edge_jBS, parameterize_jBS, scale_jBS, &
                               j_spike, j_spike_masked)
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(n_psi)  !< Normalised psi grid [0,1], arbitrary spacing
REAL(r8), INTENT(out) :: j_BS(n_psi)
LOGICAL,  OPTIONAL, INTENT(in)  :: isolate_edge_jBS  !< If .TRUE., isolate edge spike from core
LOGICAL,  OPTIONAL, INTENT(in)  :: parameterize_jBS  !< If .TRUE., use parametrised skew-normal fit
REAL(r8), OPTIONAL, INTENT(in)  :: scale_jBS         !< Scaling factor applied to spike profile (default 1)
REAL(r8), OPTIONAL, INTENT(out) :: j_spike(n_psi)        !< Processed spike profile [A/m^2]
REAL(r8), OPTIONAL, INTENT(out) :: j_spike_masked(n_psi) !< Raw masked (pre-fit) spike profile [A/m^2]
!---
INTEGER(i4) :: i
REAL(r8) :: psi_abs(n_psi)
REAL(r8) :: Te(n_psi), Ti(n_psi), ne(n_psi), ni(n_psi), Zeff(n_psi)
REAL(r8) :: pe(n_psi), pi_arr(n_psi)
REAL(r8) :: f(n_psi)        !< I(psi) = R*Bt [T*m]
REAL(r8) :: fc(n_psi)       !< Circulating fraction
REAL(r8) :: ft(n_psi)       !< Trapped particle fraction
REAL(r8) :: eps(n_psi)      !< Inverse aspect ratio
REAL(r8) :: qvals(n_psi)    !< Safety factor
REAL(r8) :: R_avg(n_psi)    !< <R> [m] per flux surface
REAL(r8) :: B_avg(n_psi)    !< <B> [T] per flux surface
REAL(r8) :: r_avgs_saut(n_psi,3)    !< Sauter FSA: <R>, <1/R>, <a>
REAL(r8) :: modb_avgs_saut(n_psi,2) !< Sauter FSA: <|B|>, <|B|^2>
REAL(r8) :: dn_e_dpsi(n_psi), dT_e_dpsi(n_psi)
REAL(r8) :: dn_i_dpsi(n_psi), dT_i_dpsi(n_psi)
REAL(r8) :: ln_le(n_psi), ln_lii(n_psi), Z_lnLam(n_psi)
REAL(r8) :: Zavg(n_psi), Zion(n_psi)
REAL(r8) :: nu_i_star(n_psi), nu_e_star(n_psi)
REAL(r8) :: B_times_Jbs(n_psi)
REAL(r8) :: psi_range, Zdom
REAL(r8), PARAMETER :: EC = 1.602176634e-19_r8
! Locals for optional edge-spike isolation
LOGICAL  :: do_isolate, do_parametrize
REAL(r8) :: scl_jBS
! Workspace for psi_N convention flip (OFT: 0=LCFS,1=axis → standard: 0=axis,1=LCFS)
REAL(r8), ALLOCATABLE :: psi_N_std(:), j_BS_std(:), mask_std(:), param_std(:)
!---
! Compute raw psi values from the caller-supplied psi_N grid
psi_range = gseq%plasma_bounds(2) - gseq%plasma_bounds(1)
psi_abs = gseq%plasma_bounds(1) + psi_N * psi_range
! Evaluate kinetic profiles on the psi grid; convert Te/Ti from keV to eV
DO i = 1, n_psi
  Te(i) = gseq%Te%fp(psi_N(i)) * 1000.0_r8
  Ti(i) = gseq%Ti%fp(psi_N(i)) * 1000.0_r8
  ne(i) = gseq%ne%fp(psi_N(i))
  ni(i) = gseq%ni%fp(psi_N(i))
  Zeff(i) = gseq%Zeff%fp(psi_N(i))
END DO
IF(MAXVAL(Te) > 5.0e5_r8)CALL oft_warn('calculate_bootstrap: max(Te) > 500 keV — profiles should be in keV, not eV')
IF(MAXVAL(Ti) > 5.0e5_r8)CALL oft_warn('calculate_bootstrap: max(Ti) > 500 keV — profiles should be in keV, not eV')
! Get I(psi) = R*Bt = F(psi) profile on the psi_N grid
DO i = 1, n_psi
  IF(gseq%mode==0)THEN
    f(i) = gseq%ffp_scale * gseq%I%f(psi_abs(i)) + gseq%I%f_offset
  ELSE
    f(i) = SQRT(gseq%ffp_scale * gseq%I%f(psi_abs(i)) + gseq%I%f_offset**2) &
           + gseq%I%f_offset * (1.0_r8 - SIGN(1.0_r8, gseq%I%f_offset))
  END IF
END DO
! Get flux-surface geometry: fc, eps, q, and <R> in one tracing pass
CALL sauter_fc(gseq, n_psi, psi_N, fc, r_avgs_saut, modb_avgs_saut, qprof=qvals, eps=eps)
R_avg = r_avgs_saut(:,1)
B_avg = modb_avgs_saut(:,1)
! ===================================================================
ft = 1.0_r8 - fc
! Pressures [Pa]
pe = EC * ne * Te
pi_arr = EC * ni * Ti
! Gradients d/dpsi [Wb^-1] via non-uniform finite differences
! (one-sided at endpoints, central at interior points)
dT_e_dpsi(1)     = (Te(2)     - Te(1))       / (psi_abs(2)     - psi_abs(1))
dT_e_dpsi(n_psi) = (Te(n_psi) - Te(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dT_i_dpsi(1)     = (Ti(2)     - Ti(1))       / (psi_abs(2)     - psi_abs(1))
dT_i_dpsi(n_psi) = (Ti(n_psi) - Ti(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dn_e_dpsi(1)     = (ne(2)     - ne(1))       / (psi_abs(2)     - psi_abs(1))
dn_e_dpsi(n_psi) = (ne(n_psi) - ne(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dn_i_dpsi(1)     = (ni(2)     - ni(1))       / (psi_abs(2)     - psi_abs(1))
dn_i_dpsi(n_psi) = (ni(n_psi) - ni(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
DO i = 2, n_psi-1
  dT_e_dpsi(i) = (Te(i+1) - Te(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dT_i_dpsi(i) = (Ti(i+1) - Ti(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dn_e_dpsi(i) = (ne(i+1) - ne(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dn_i_dpsi(i) = (ni(i+1) - ni(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
END DO
! In the Fortran internal convention psi increases LCFS→axis (opposite to the
! standard Sauter/Redl derivation where psi increases axis→LCFS).  Negate
! all gradients so the bootstrap formula sees the conventional sign.
dT_e_dpsi = -dT_e_dpsi
dT_i_dpsi = -dT_i_dpsi
dn_e_dpsi = -dn_e_dpsi
dn_i_dpsi = -dn_i_dpsi
! Coulomb logarithms (NRL formulary)
! Electron: ne divided by 1e6 to convert m^-3 -> cm^-3
ln_le = 23.5_r8 &
      - LOG(SQRT(ne/1.0e6_r8) * Te**(-1.25_r8)) &
      - SQRT(1.0e-5_r8 + (LOG(Te) - 2.0_r8)**2 / 16.0_r8)
! Ion: constant 30 absorbs the m^-3 -> cm^-3 conversion
Z_lnLam = MAX(ne/ni, 1.0_r8)
ln_lii = 30.0_r8 - LOG(Z_lnLam**3 * SQRT(ni) / Ti**1.5_r8)
ln_le  = MAX(ln_le,  10.0_r8)
ln_lii = MAX(ln_lii, 10.0_r8)
! Ion collisionality: Koh multi-species model
Zdom = 1.0_r8    ! dominant ion charge (deuterium)
Zavg = ne / ni
Zion = (Zdom**2 * Zavg * Zeff)**0.25_r8
nu_i_star = 4.90e-18_r8 * ABS(qvals) * R_avg * ni &
          * Zion**4 * ln_lii / (Ti**2 * eps**1.5_r8)
! Electron collisionality
nu_e_star = 6.921e-18_r8 * ABS(qvals) * R_avg * ne &
          * Zeff * ln_le / (Te**2 * eps**1.5_r8)
! Compute bootstrap current via Redl 2021
CALL redl_bootstrap(n_psi, Te, Ti, ne, ni, pe, pi_arr, Zeff, qvals, eps, ft, f, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    ln_le, ln_lii, nu_e_star, nu_i_star, B_times_Jbs)
! Convert parallel bootstrap to phi component: j_phi = B_times_Jbs * <R> / F
WHERE(ABS(f) > 0.0_r8)
  j_BS = B_times_Jbs / B_avg
ELSEWHERE
  j_BS = 0.0_r8
END WHERE
! Guard NaN (where F -> 0)
WHERE(.NOT.(ABS(j_BS) < 1.0e99_r8)) j_BS = 0.0_r8
! Save raw bootstrap output
IF(.NOT.ASSOCIATED(gseq%boot_profs%j_bs_raw)) ALLOCATE(gseq%boot_profs%j_bs_raw(0:n_psi))
gseq%boot_profs%j_bs_raw(0)  = 0.0_r8  ! j_BS = 0 at LCFS (OFT psi = 0)
gseq%boot_profs%j_bs_raw(1:) = j_BS
IF(gseq%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A)') '  [calculate_bootstrap] geometry & collisionality sample (i=1,mid,n):'
  WRITE(*,'(A,3ES12.4)') '    <R>      : ', r_avgs_saut(1,1), r_avgs_saut(n_psi/2,1), r_avgs_saut(n_psi,1)
  WRITE(*,'(A,3ES12.4)') '    <1/R>    : ', r_avgs_saut(1,2), r_avgs_saut(n_psi/2,2), r_avgs_saut(n_psi,2)
  WRITE(*,'(A,3ES12.4)') '    <B>      : ', B_avg(1), B_avg(n_psi/2), B_avg(n_psi)
  WRITE(*,'(A,3ES12.4)') '    eps      : ', eps(1), eps(n_psi/2), eps(n_psi)
  WRITE(*,'(A,3ES12.4)') '    q        : ', qvals(1), qvals(n_psi/2), qvals(n_psi)
  WRITE(*,'(A,3ES12.4)') '    nu_e_star: ', nu_e_star(1), nu_e_star(n_psi/2), nu_e_star(n_psi)
  WRITE(*,'(A,3ES12.4)') '    j_BS     : ', j_BS(1), j_BS(n_psi/2), j_BS(n_psi)
END IF
!---
! Optionally isolate the edge bootstrap spike and return as j_spike,
! mirroring the isolate_edge_jBS, parameterize_jBS logic in
! bootstrap.py:calculate_profiles_and_bootstrap.
IF (PRESENT(j_spike)) THEN
  do_isolate    = .FALSE.
  do_parametrize = .FALSE.
  scl_jBS       = 1.0_r8
  IF (PRESENT(isolate_edge_jBS)) do_isolate    = isolate_edge_jBS
  IF (PRESENT(parameterize_jBS)) do_parametrize = parameterize_jBS
  IF (PRESENT(scale_jBS))        scl_jBS        = scale_jBS
  IF (do_isolate .OR. do_parametrize) THEN
    ! analyze_bootstrap_edge_spike uses standard psi_N convention (0=axis, 1=LCFS).
    ! OFT internal convention is reversed (0=LCFS, 1=axis).
    ! Flip arrays before calling, then flip outputs back.
    ALLOCATE(psi_N_std(n_psi), j_BS_std(n_psi), mask_std(n_psi))
    IF (do_parametrize) ALLOCATE(param_std(n_psi))
    psi_N_std = 1.0_r8 - psi_N(n_psi:1:-1)
    j_BS_std  = j_BS(n_psi:1:-1)
    IF (do_parametrize) THEN
      CALL analyze_bootstrap_edge_spike(n_psi, psi_N_std, j_BS_std, mask_std, &
                                      param_std, diagnose=gseq%boot_ops%diagnose_bs)
      IF (PRESENT(j_spike))        j_spike        = scl_jBS * param_std(n_psi:1:-1)
      IF (PRESENT(j_spike_masked)) j_spike_masked = scl_jBS * mask_std(n_psi:1:-1)
      DEALLOCATE(param_std)
    ELSE
      CALL analyze_bootstrap_edge_spike(n_psi, psi_N_std, j_BS_std, mask_std)
      IF (PRESENT(j_spike))        j_spike        = scl_jBS * mask_std(n_psi:1:-1)
      IF (PRESENT(j_spike_masked)) j_spike_masked = scl_jBS * mask_std(n_psi:1:-1)
    END IF
    DEALLOCATE(psi_N_std, j_BS_std, mask_std)
  END IF
END IF
END SUBROUTINE calculate_bootstrap
!------------------------------------------------------------------------------
!> Evaluate the skew-normal PDF at a single point.
!>
!> Equivalent to scipy.stats.skewnorm.pdf(x, sk, loc=ctr, scale=scl):
!>   f = (2/scl) * phi((x-ctr)/scl) * Phi(sk*(x-ctr)/scl)
!> where phi and Phi are the standard normal PDF and CDF.
!>
!> @param x   Evaluation point
!> @param sk  Skewness parameter
!> @param ctr Location (mean shift)
!> @param scl Scale (width)
!> @result    PDF value at x
!------------------------------------------------------------------------------
PURE FUNCTION skewnorm_pdf_pt(x, sk, ctr, scl) RESULT(y)
REAL(r8), INTENT(in) :: x, sk, ctr, scl
REAL(r8) :: y
REAL(r8), PARAMETER :: INV_SQRT2PI = 0.39894228040143268_r8  ! 1/sqrt(2*pi)
REAL(r8), PARAMETER :: INV_SQRT2   = 0.70710678118654752_r8  ! 1/sqrt(2)
REAL(r8) :: z, phi_val, cdf_val
z       =  (x - ctr) / scl
phi_val =  INV_SQRT2PI * EXP(-0.5_r8 * z**2)
cdf_val =  0.5_r8 * (1.0_r8 + ERF(sk * z * INV_SQRT2))
y       =  (2.0_r8 / scl) * phi_val * cdf_val
END FUNCTION skewnorm_pdf_pt
!------------------------------------------------------------------------------
!> Compute log(exp(a) + exp(b)) stably (avoids overflow).
!> Mirrors numpy.logaddexp(a, b).
!------------------------------------------------------------------------------
PURE FUNCTION safe_logaddexp(a, b) RESULT(y)
REAL(r8), INTENT(in) :: a, b
REAL(r8) :: y, c
c = MAX(a, b)
y = c + LOG(EXP(a - c) + EXP(b - c))
END FUNCTION safe_logaddexp
!------------------------------------------------------------------------------
!> Locate the peak of the skew-normal shape used in parametrise_edge_jbs
!> using ternary search on the interval
!>   [max(0, center - 3*width), min(1, center + 3*width)]
!>
!> The skew-normal PDF is strictly unimodal for all skewness values (Azzalini
!> 1985), so ternary search is guaranteed to converge.  Each iteration cuts
!> the search interval by a factor of 2/3; 100 iterations reduce an initial
!> width of ~6*width to < 10^{-17} * initial_width (essentially machine
!> precision), while requiring only 200 PDF evaluations — 50x fewer than the
!> previous 10 000-point linear scan.
!>
!> @param center       Gaussian centre (loc)
!> @param width        Gaussian width (scale)
!> @param sk           Skewness parameter
!> @param x_peak       Output: psi_N at the raw-shape peak
!> @param val_peak_raw Output: raw-shape value at the peak
!------------------------------------------------------------------------------
SUBROUTINE find_skewnorm_peak(center, width, sk, x_peak, val_peak_raw)
REAL(r8), INTENT(in)  :: center, width, sk
REAL(r8), INTENT(out) :: x_peak, val_peak_raw
INTEGER(i4), PARAMETER :: MAX_ITER = 100  ! (2/3)^100 * 6*width < 1e-17*width
INTEGER(i4) :: iter
REAL(r8)    :: x_lo, x_hi, x1, x2, f1, f2
x_lo = MAX(0.0_r8, center - 3.0_r8*width)
x_hi = MIN(1.0_r8, center + 3.0_r8*width)
! Guard: if range collapses set a minimal width
IF (x_hi <= x_lo) x_hi = x_lo + 1.0e-6_r8
! Ternary search: at each step evaluate the PDF at the two third-points.
! The unimodality guarantee ensures we can safely discard one third of the
! interval according to which third-point has the lower value.
DO iter = 1, MAX_ITER
  x1 = x_lo + (x_hi - x_lo) / 3.0_r8
  x2 = x_hi - (x_hi - x_lo) / 3.0_r8
  f1 = skewnorm_pdf_pt(x1, sk, center, width)
  f2 = skewnorm_pdf_pt(x2, sk, center, width)
  IF (f1 < f2) THEN
    x_lo = x1   ! peak cannot be in [x_lo, x1]
  ELSE
    x_hi = x2   ! peak cannot be in [x2, x_hi]
  END IF
END DO
x_peak       = 0.5_r8 * (x_lo + x_hi)
val_peak_raw = skewnorm_pdf_pt(x_peak, sk, center, width)
! Avoid division by zero downstream
IF (val_peak_raw < 1.0e-30_r8) val_peak_raw = 1.0e-30_r8
END SUBROUTINE find_skewnorm_peak
!------------------------------------------------------------------------------
!> Evaluate the parametrised skewnorm profile on an arbitrary psi grid.
!>
!> Direct translation of the inner function generate_baseline_prof from
!> Python bootstrap.py:parameterize_edge_jBS.  Constructs the profile by
!> stitching a left-side SoftMax-blended skew-normal spike with a right-side
!> cosine (or cosh) decay at x_peak:
!>
!>   left  (psi <= x_peak):  SoftMax(offset_in, amp*skewnorm/val_peak_raw)
!>   right (psi >  x_peak):  stitch_height * cos(omega*(psi-x_peak))^tail_alpha
!>
!> x_peak and val_peak_raw must be pre-computed via find_skewnorm_peak so
!> the golden-section scan is not repeated for every function evaluation.
!>
!> @param n            Number of psi points
!> @param psi          psi_N grid [0, 1]
!> @param amp          Amplitude of the skew-normal spike
!> @param center       Spike centre in psi_N
!> @param width        Spike width (sigma = FWHM/2.355)
!> @param offset_in    Flat baseline J_BS level left of the spike
!> @param sk           Skewness parameter (a in skewnorm.pdf)
!> @param y_sep        Profile value prescribed at the separatrix (psi_N = 1)
!> @param blend_width  SoftMax blend width (sharpness of left-side stitch)
!> @param tail_alpha   Right-side decay exponent (>= 1)
!> @param x_peak       Peak location of the raw skew-normal (from find_skewnorm_peak)
!> @param val_peak_raw Peak value   of the raw skew-normal (from find_skewnorm_peak)
!> @param profile      Output: profile values on the psi grid
!------------------------------------------------------------------------------
SUBROUTINE eval_baseline_profile(n, psi, amp, center, width, offset_in, sk, &
    y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, profile)
INTEGER(i4), INTENT(in)  :: n
REAL(r8),    INTENT(in)  :: psi(n)
REAL(r8),    INTENT(in)  :: amp, center, width, offset_in, sk
REAL(r8),    INTENT(in)  :: y_sep, blend_width, tail_alpha
REAL(r8),    INTENT(in)  :: x_peak, val_peak_raw
REAL(r8),    INTENT(out) :: profile(n)
!---
INTEGER(i4) :: i
REAL(r8) :: k_smooth, diff_ao, argument, internal_amp
REAL(r8) :: stitch_height, dist_to_edge
REAL(r8) :: target_cos, omega_cos, target_cosh, omega_cosh
REAL(r8) :: spike_val, pL, pR, arg_cos, base_cos
LOGICAL  :: use_cos_decay
! =====================================================================
! Smoothing factor (mirrors Python k_smooth = amp/width * blend_width/4)
! =====================================================================
k_smooth = MAX((amp / width) * (blend_width / 4.0_r8), 1.0e-5_r8)
! =====================================================================
! Left-side parameters (Case A: amp > offset_in; Case B: otherwise)
! =====================================================================
IF (amp > offset_in) THEN
  ! Case A: standard spike
  diff_ao = amp - offset_in
  IF (diff_ao > 1.0e-10_r8) THEN
    argument    = MAX(1.0_r8 - EXP((offset_in - amp) / k_smooth), 1.0e-16_r8)
    internal_amp = amp + k_smooth * LOG(argument)
  ELSE
    internal_amp = amp
  END IF
  stitch_height = amp
ELSE
  ! Case B: dominant offset
  stitch_height = offset_in
END IF
! =====================================================================
! Right-side parameters (cosine decay or cosh rise)
! =====================================================================
dist_to_edge = MAX(1.0_r8 - x_peak, 1.0e-5_r8)
use_cos_decay = (y_sep < stitch_height) .AND. (ABS(stitch_height) > 1.0e-30_r8)
IF (use_cos_decay) THEN
  ! Cosine decay: stitch_height * cos(omega*(psi-x_peak))^tail_alpha
  target_cos = MAX(-1.0_r8, MIN(1.0_r8, (y_sep / stitch_height)**(1.0_r8/tail_alpha)))
  omega_cos  = ACOS(target_cos) / dist_to_edge
ELSE IF (ABS(stitch_height) > 1.0e-30_r8) THEN
  ! Cosh rise (rare: y_sep >= stitch_height)
  target_cosh = (y_sep / stitch_height)**(1.0_r8/tail_alpha)
  omega_cosh  = ACOSH(MAX(target_cosh, 1.0_r8)) / dist_to_edge
ELSE
  omega_cos = 0.0_r8
  omega_cosh = 0.0_r8
END IF
! =====================================================================
! Build full profile: stitch left side and right side at x_peak
! =====================================================================
DO i = 1, n
  ! Left side
  IF (amp > offset_in) THEN
    spike_val = skewnorm_pdf_pt(psi(i), sk, center, width) / val_peak_raw * internal_amp
    pL = safe_logaddexp(offset_in / k_smooth, spike_val / k_smooth) * k_smooth
  ELSE
    pL = offset_in
  END IF
  ! Right side
  IF (use_cos_decay) THEN
    arg_cos = omega_cos * (psi(i) - x_peak)
    base_cos = MAX(COS(arg_cos), 0.0_r8)
    pR = stitch_height * (base_cos**tail_alpha)
  ELSE IF (ABS(stitch_height) > 1.0e-30_r8) THEN
    pR = stitch_height * (COSH(omega_cosh * (psi(i) - x_peak))**tail_alpha)
  ELSE
    pR = 0.0_r8
  END IF
  ! Stitch
  IF (psi(i) <= x_peak) THEN
    profile(i) = pL
  ELSE
    profile(i) = pR
  END IF
  ! Clamp (mirrors: if y_sep >= 0: profile = maximum(profile, 0))
  IF (y_sep >= 0.0_r8) profile(i) = MAX(profile(i), 0.0_r8)
END DO
END SUBROUTINE eval_baseline_profile
!------------------------------------------------------------------------------
!> Fortran equivalent of Python bootstrap.py:parameterize_edge_jBS.
!>
!> Generates the parameterised edge-bootstrap profile on an n-point psi_N
!> grid.  When offset /= 0 an internal 1-D root find (bisection on alpha)
!> ensures that profile(psi(1)) matches the requested offset level, mirroring
!> the root_scalar(brentq) call in the Python.  When offset == 0 the function
!> falls through directly to eval_baseline_profile with a very negative
!> effective offset (reproducing the Python `-1e10` shortcut).
!>
!> Fixed parameter: tail_alpha = 1.5 (Python default, not a fit variable).
!>
!> @param n           Number of psi_N points
!> @param psi         Normalised poloidal flux grid [0, 1]
!> @param amp         Spike amplitude
!> @param center      Spike centre in psi_N
!> @param width       Spike width (sigma of skew-normal)
!> @param offset      Requested flat-core level
!> @param sk          Skewness parameter
!> @param y_sep       Value at the separatrix (psi_N = 1)
!> @param blend_width SoftMax blending width
!> @param tail_alpha  Right-side fall-off exponent
!> @param profile     Output: profile values on the psi grid
!------------------------------------------------------------------------------
SUBROUTINE parametrise_edge_jbs(n, psi, amp, center, width, offset, sk, &
    y_sep, blend_width, tail_alpha, profile)
INTEGER(i4), INTENT(in)  :: n
REAL(r8),    INTENT(in)  :: psi(n)
REAL(r8),    INTENT(in)  :: amp, center, width, offset, sk
REAL(r8),    INTENT(in)  :: y_sep, blend_width, tail_alpha
REAL(r8),    INTENT(out) :: profile(n)
!---
INTEGER(i4) :: iter
REAL(r8)    :: x_peak, val_peak_raw
REAL(r8)    :: a_lo, a_hi, a_mid, f_lo, f_hi, f_mid
REAL(r8)    :: a_optimal
REAL(r8)    :: psi1(1), prof1(1)  ! scratch 1-element arrays for bisection
! =====================================================================
! Pre-compute the raw skew-normal peak (independent of offset/alpha)
! =====================================================================
CALL find_skewnorm_peak(center, width, sk, x_peak, val_peak_raw)
! =====================================================================
! Determine the offset-corrected scaling alpha via bisection.
!
! Objective: f(alpha) = profile_at_psi1(alpha*offset) - offset = 0
! Mirrors Python: root_scalar(obj, bracket=[-1e10, 10*amp], method='brentq')
!
! When offset == 0 skip the root find and use a deeply negative effective
! offset (-1e20) to fully suppress the SoftMax baseline (pure spike, no floor).
! =====================================================================
IF (ABS(offset) < 1.0e-30_r8) THEN
  ! offset = 0: pass a deeply negative effective offset so the SoftMax baseline
  ! is completely suppressed (pure spike, no floor).  -1e20 ensures this holds
  ! even for amp up to ~1e17 A/m^2; -1e10 would become marginal above ~50 MA/m^2.
  CALL eval_baseline_profile(n, psi, amp, center, width, -1.0e20_r8, sk, &
      y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, profile)
  RETURN
END IF
! Bisection bracket: mirrors Python [-1e10, 10*amp]
! Use alpha values that drive alpha*offset to the extremes
a_lo = -1.0e6_r8 / MAX(ABS(offset), 1.0e-30_r8)  ! alpha*offset ≈ -1e6
a_hi = MAX(10.0_r8 * ABS(amp), 1.0_r8) / MAX(ABS(offset), 1.0e-30_r8)
! Evaluate f at the bracket ends
psi1(1) = psi(1)
CALL eval_baseline_profile(1, psi1, amp, center, width, a_lo*offset, sk, &
    y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, prof1)
f_lo = prof1(1) - offset
CALL eval_baseline_profile(1, psi1, amp, center, width, a_hi*offset, sk, &
    y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, prof1)
f_hi = prof1(1) - offset
! If the bracket already straddles the root proceed; otherwise fall back
! to a_optimal = 1 (identity scaling) to avoid NaN propagation.
IF (f_lo * f_hi > 0.0_r8) THEN
  a_optimal = 1.0_r8
ELSE
  ! 60 bisection iterations -> relative bracket width < 1e-18
  DO iter = 1, 60
    a_mid = 0.5_r8 * (a_lo + a_hi)
    CALL eval_baseline_profile(1, psi1, amp, center, width, a_mid*offset, sk, &
        y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, prof1)
    f_mid = prof1(1) - offset
    IF (ABS(f_mid) < 1.0e-6_r8 * MAX(ABS(offset), 1.0e-30_r8)) EXIT
    IF (SIGN(1.0_r8, f_mid) == SIGN(1.0_r8, f_lo)) THEN
      a_lo = a_mid;  f_lo = f_mid
    ELSE
      a_hi = a_mid;  f_hi = f_mid
    END IF
  END DO
  a_optimal = a_mid
END IF
! =====================================================================
! Final evaluation on the full n-point psi grid
! =====================================================================
CALL eval_baseline_profile(n, psi, amp, center, width, a_optimal*offset, sk, &
    y_sep, blend_width, tail_alpha, x_peak, val_peak_raw, profile)
END SUBROUTINE parametrise_edge_jbs
!------------------------------------------------------------------------------
!> MINPACK lmdif residual subroutine for curve_fit_edge_jbs.
!>
!> Computes fvec(i) = parametrise_edge_jbs(psi_fit(i); cofs) - j_fit(i)
!> where cofs(7) = (amp, center, width, offset, sk, y_sep, blend_width) and
!> the fitting data are taken from the module-level active_edge_jbs context.
!>
!> Signature required by lmdif:  SUBROUTINE fcn(m, n, x, fvec, iflag)
!------------------------------------------------------------------------------
SUBROUTINE edge_jbs_residual(m, n, cofs, err, iflag)
INTEGER(4),  INTENT(in)    :: m, n
REAL(8),     INTENT(in)    :: cofs(n)
REAL(8),     INTENT(out)   :: err(m)
INTEGER(4),  INTENT(inout) :: iflag
REAL(r8) :: amp, center, width, offset, sk, y_sep, blend_width
REAL(r8) :: p_bounded(7)
REAL(r8) :: profile(active_edge_jbs%n_fit)
INTEGER(i4) :: k
!---
! Inverse transform: P_bounded = lb + (tanh(P_internal) + 1) * (ub - lb) / 2
DO k = 1, 7
  p_bounded(k) = active_edge_jbs%lb(k) &
               + (TANH(cofs(k)) + 1.0_r8) &
               * (active_edge_jbs%ub(k) - active_edge_jbs%lb(k)) * 0.5_r8
END DO
amp         = p_bounded(1)
center      = p_bounded(2)
width       = p_bounded(3)
offset      = p_bounded(4)
sk          = p_bounded(5)
y_sep       = p_bounded(6)
blend_width = p_bounded(7)
CALL parametrise_edge_jbs(active_edge_jbs%n_fit, active_edge_jbs%psi_fit, &
    amp, center, width, offset, sk, y_sep, blend_width, &
    active_edge_jbs%tail_alpha, profile)
err = profile - active_edge_jbs%j_fit
END SUBROUTINE edge_jbs_residual
!------------------------------------------------------------------------------
!> Levenberg-Marquardt least-squares fit of parametrise_edge_jbs to a set
!> of (psi_fit, j_BS_fit) data points using MINPACK lmdif.
!>
!> Equivalent to the scipy.optimize.curve_fit call in Python's
!> analyze_bootstrap_edge_spike (using the same 7-parameter model but without
!> parameter bounds, which lmdif does not support natively).
!>
!> The 7 parameters and their Python / curve_fit correspondences are:
!>   p(1)  amp         -- spike amplitude
!>   p(2)  center      -- spike centre in psi_N
!>   p(3)  width       -- spike sigma
!>   p(4)  offset      -- flat-core level
!>   p(5)  sk          -- skewness
!>   p(6)  y_sep       -- value at psi_N = 1
!>   p(7)  blend_width -- SoftMax blending width
!>
!> @param n_fit    Number of data points (m == n_fit for lmdif)
!> @param psi_fit  psi_N values of the fitting points
!> @param j_BS_fit Target j_BS values [A/m^2] at psi_fit
!> @param p0       Initial parameter guess (7 elements)
!> @param popt     Output: optimised parameters (7 elements)
!------------------------------------------------------------------------------
SUBROUTINE curve_fit_edge_jbs(n_fit, psi_fit, j_BS_fit, p0, popt, lb_in, ub_in)
INTEGER(i4), INTENT(in)  :: n_fit
REAL(r8),    INTENT(in)  :: psi_fit(n_fit), j_BS_fit(n_fit)
REAL(r8),    INTENT(in)  :: p0(7)
REAL(r8),    INTENT(out) :: popt(7)
REAL(r8),    INTENT(in)  :: lb_in(7), ub_in(7)
!---MINPACK variables (same pattern as gs_psi2pt / circle_interp)
INTEGER(4), PARAMETER :: NPARAMS = 7
REAL(8) :: ftol, xtol, gtol, epsfcn, factor
REAL(8) :: cofs(NPARAMS), error_vec(n_fit)
INTEGER(i4) :: k
REAL(8), ALLOCATABLE :: diag(:), wa1(:), wa2(:), wa3(:), wa4(:), qtf(:)
REAL(8), ALLOCATABLE :: fjac(:,:)
INTEGER(4), ALLOCATABLE :: ipvt(:)
INTEGER(4) :: maxfev, mode, nprint, info, nfev, ldfjac, ncons, ncofs  ! nfev: function eval count (informational)
CHARACTER(len=80) :: char_buf
!---
! Populate module-level context so edge_jbs_residual can access the data
IF (ALLOCATED(active_edge_jbs%psi_fit)) DEALLOCATE(active_edge_jbs%psi_fit)
IF (ALLOCATED(active_edge_jbs%j_fit))   DEALLOCATE(active_edge_jbs%j_fit)
ALLOCATE(active_edge_jbs%psi_fit(n_fit))
ALLOCATE(active_edge_jbs%j_fit(n_fit))
active_edge_jbs%n_fit   = n_fit
active_edge_jbs%psi_fit = psi_fit
active_edge_jbs%j_fit   = j_BS_fit
active_edge_jbs%tail_alpha = 1.5_r8  ! fixed Python default
! Bounds are passed in via the module-level context; declared as extra args below
active_edge_jbs%lb = lb_in
active_edge_jbs%ub = ub_in
! Clamp p0 to lie strictly inside [lb + 1e-6*range, ub - 1e-6*range]
! to avoid atanh(+-1) = +-Inf at the bounds.
DO k = 1, NPARAMS
  cofs(k) = MAX(lb_in(k) + 1.0e-6_r8 * (ub_in(k) - lb_in(k)), &
            MIN(ub_in(k) - 1.0e-6_r8 * (ub_in(k) - lb_in(k)), p0(k)))
END DO
! Forward transform: P_internal = atanh(2*(P_bounded - lb)/(ub - lb) - 1)
DO k = 1, NPARAMS
  cofs(k) = ATANH(2.0_r8 * (cofs(k) - lb_in(k)) / (ub_in(k) - lb_in(k)) - 1.0_r8)
END DO
ncons  = n_fit
ncofs  = NPARAMS
ldfjac = ncons
ALLOCATE(diag(ncofs), fjac(ncons, ncofs))
ALLOCATE(qtf(ncofs), wa1(ncofs), wa2(ncofs))
ALLOCATE(wa3(ncofs), wa4(ncons))
ALLOCATE(ipvt(ncofs))
! MINPACK lmdif settings (mirrors existing usage in grad_shaf.F90)
mode   = 1
factor = 1.0d0
maxfev = 10000     ! matches Python maxfev=10000
ftol   = 1.0d-6
xtol   = 1.0d-6
gtol   = 1.0d-6
epsfcn = 1.0d-6
nprint = 0
CALL lmdif(edge_jbs_residual, ncons, ncofs, cofs, error_vec, &
           ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, &
           info, nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
DEALLOCATE(diag, fjac, qtf, wa1, wa2, wa3, wa4, ipvt)
! info: 1-4 = converged, 0 = improper input, 5 = maxfev exceeded, 6-7 = tolerance too small
IF (info <= 0 .OR. info >= 5) THEN
  WRITE(char_buf,'(A,I0,A,I0)') '[curve_fit_edge_jbs] lmdif did not converge; info=', info, &
    '; nfev=', nfev
  CALL oft_warn(TRIM(char_buf))
END IF
! Inverse transform: recover bounded parameters from internal lmdif solution
DO k = 1, NPARAMS
  popt(k) = lb_in(k) + (TANH(cofs(k)) + 1.0_r8) * (ub_in(k) - lb_in(k)) * 0.5_r8
END DO
END SUBROUTINE curve_fit_edge_jbs
!------------------------------------------------------------------------------
!> Analyse and isolate the edge bootstrap-current spike from a toroidal
!> bootstrap current density profile.
!>
!> Translated from Python bootstrap.py: analyze_bootstrap_edge_spike
!>
!> The routine performs the following steps:
!>   1. Locates the dominant local peak with j_BS > 0 in the edge region
!>      (psi_N > 0.7).  When several peaks exist the one closest to the
!>      separatrix (largest psi_N) is chosen.  A warning is emitted if
!>      that peak is not also the tallest among qualifying peaks.
!>   2. Finds the minimum of j_BS between psi_N = 0.5 and the peak (the
!>      flat-core stitching level, lmin_j_BS).
!>   3. Builds masked_spike: flat at lmin_j_BS for psi_N below the minimum
!>      index, actual j_BS from that index to the separatrix.
!>   4. Smooths the flat/spike junction with a cubic Hermite spline over a
!>      window of width min(0.5 * dist_to_peak, 0.2) centred on the minimum.
!>   5. (Optional) Fits parametrise_edge_jbs to j_BS over the spike region
!>      (psi_N >= psi_N(lmin_idx)) via curve_fit_edge_jbs (MINPACK lmdif).
!>      Only performed when parameterized_spike is present.  An approximate
!>      FWHM is computed first as the initial sigma guess.  The 7 fitted
!>      parameters are used to evaluate parameterized_spike on the full
!>      n-point psi_N grid.
!>
!> @param n                  Number of flux surface samples
!> @param psi_N              Normalised poloidal flux grid [0, 1]
!> @param j_BS           Bootstrap current density j_phi profile [A/m^2]
!> @param masked_spike       Output: isolated edge spike spliced onto flat core [A/m^2]
!> @param parameterized_spike Output: parametrise_edge_jbs fit evaluated on full psi_N grid [A/m^2]
!------------------------------------------------------------------------------
SUBROUTINE analyze_bootstrap_edge_spike(n, psi_N, j_BS, masked_spike, &
                                      parameterized_spike, diagnose)
INTEGER(i4), INTENT(in)  :: n
REAL(r8),    INTENT(in)  :: psi_N(n)              !< Normalised poloidal flux [0=axis, 1=LCFS/plasma edge]
REAL(r8),    INTENT(in)  :: j_BS(n)           !< Bootstrap current density [A/m^2]
REAL(r8),    INTENT(out) :: masked_spike(n)        !< Isolated edge spike spliced onto flat core [A/m^2]
REAL(r8), OPTIONAL, INTENT(out) :: parameterized_spike(n) !< parametrise_edge_jbs fit on full grid [A/m^2]
LOGICAL,  OPTIONAL, INTENT(in)  :: diagnose               !< If .TRUE., print fit parameters and profile table
!---
INTEGER(i4) :: i, peak_idx, lmin_idx, left_idx, right_idx, idx_start, idx_end
INTEGER(i4) :: n_fit
REAL(r8)    :: best_psi, best_height, half_max, dist_tmp
REAL(r8)    :: peak_psi, peak_height, lmin_j_BS, fwhm
REAL(r8)    :: jBS_min_loc, dist_to_peak, blend_width_val, x_start, x_end
REAL(r8)    :: dist_start, dist_end
REAL(r8)    :: y_start, dy_start, y_end, dy_end
REAL(r8)    :: dx_window, t, h00, h10, h01, h11, y_patch
REAL(r8)    :: p0(7), popt(7)
REAL(r8)    :: lb(7), ub(7)
REAL(r8)    :: amp_lo, amp_hi, off_lo, off_hi, ysep_hi, eps_tmp
REAL(r8)    :: amp_fit, center_fit, width_fit, offset_fit, sk_fit, y_sep_fit, bw_fit
REAL(r8), ALLOCATABLE :: psi_fit(:), j_fit(:)
CHARACTER(len=256) :: char_buf
! =====================================================================
! 1. Find the dominant local peak in the edge region (psi_N > 0.7)
!    with j_BS > 0.  Among all qualifying peaks choose the one with
!    the largest psi_N value (closest to the separatrix).
!    Warn if that peak is not also the tallest.
!    Mirrors: peaks[numpy.argmax(psi_edge[peaks])] in Python.
! =====================================================================
peak_idx    = -1
best_psi    = -1.0_r8
best_height = -1.0_r8
DO i = 2, n-1
  IF (psi_N(i) < 0.7_r8)          CYCLE
  IF (j_BS(i) <= 0.0_r8)      CYCLE
  IF (j_BS(i) <= j_BS(i-1)) CYCLE
  IF (j_BS(i) <= j_BS(i+1)) CYCLE
  IF (j_BS(i) > best_height) best_height = j_BS(i)
  IF (psi_N(i) > best_psi) THEN
    peak_idx = i
    best_psi = psi_N(i)
  END IF
END DO
IF (peak_idx >= 1 .AND. j_BS(peak_idx) < best_height) THEN
  WRITE(char_buf,'(A,ES12.4,A,ES12.4)') &
    '[analyze_bootstrap_edge_spike] dominant (rightmost) peak height=', j_BS(peak_idx), &
    ' is not the tallest peak; tallest=', best_height
  CALL oft_warn(TRIM(char_buf))
END IF
IF (peak_idx < 1) THEN
  ! No clear peak found: return flat-zero profiles
  masked_spike = 0.0_r8
  IF (PRESENT(parameterized_spike)) parameterized_spike = 0.0_r8
  RETURN
END IF
peak_psi    = psi_N(peak_idx)
peak_height = j_BS(peak_idx)
! =====================================================================
! 2. Find the minimum of j_BS between psi_N = 0.5 and the peak.
!    This becomes the flat-core stitching level (lmin_j_BS).
!    Mirrors: lmin_j_BS = min(j_bootstrap[(psi_N > 0.5) & (psi_N < peak_psi)])
! =====================================================================
lmin_j_BS = 1.0e30_r8
lmin_idx  = peak_idx      ! fallback: use peak location if no point found
DO i = 1, n
  IF (psi_N(i) > 0.5_r8 .AND. psi_N(i) < peak_psi) THEN
    IF (j_BS(i) < lmin_j_BS) THEN
      lmin_j_BS = j_BS(i)
      lmin_idx  = i
    END IF
  END IF
END DO
! =====================================================================
! 3. Build masked_spike:
!    - Flat at lmin_j_BS for psi_N < psi_N(lmin_idx)
!    - Actual j_BS   for psi_N >= psi_N(lmin_idx)
!    Mirrors: fit_mask = (psi_N >= masked_psi_N[lmin_arg])
! =====================================================================
DO i = 1, n
  IF (psi_N(i) >= psi_N(lmin_idx)) THEN
    masked_spike(i) = j_BS(i)
  ELSE
    masked_spike(i) = lmin_j_BS
  END IF
END DO
! =====================================================================
! 4. Smooth the flat/spike junction with a cubic Hermite spline.
!    Blend window  = min(0.5 * dist_to_peak, 0.2), centred on lmin_loc.
!    Mirrors the Hermite-spline patch in the Python.
! =====================================================================
jBS_min_loc     = psi_N(lmin_idx)
dist_to_peak    = peak_psi - jBS_min_loc
blend_width_val = MIN(0.5_r8 * dist_to_peak, 0.2_r8)
x_start         = jBS_min_loc - 0.5_r8 * blend_width_val
x_end           = jBS_min_loc + 0.5_r8 * blend_width_val
! Find the nearest grid points to x_start and x_end
idx_start  = 1
dist_start = ABS(psi_N(1) - x_start)
DO i = 2, n
  dist_tmp = ABS(psi_N(i) - x_start)
  IF (dist_tmp < dist_start) THEN
    dist_start = dist_tmp
    idx_start  = i
  END IF
END DO
idx_end  = 1
dist_end = ABS(psi_N(1) - x_end)
DO i = 2, n
  dist_tmp = ABS(psi_N(i) - x_end)
  IF (dist_tmp < dist_end) THEN
    dist_end = dist_tmp
    idx_end  = i
  END IF
END DO
! Safety clamp: mirrors Python max(0,...) and min(len-2,...) (1-indexed here)
idx_start = MAX(1,   idx_start)
idx_end   = MIN(n-1, idx_end)   ! keep room for the central-difference at idx_end+1
idx_end   = MAX(2,   idx_end)   ! keep room for the central-difference at idx_end-1
! Ensure idx_end is strictly right of lmin_idx so the central-difference
! stencil for dy_end does not straddle the flat/spike step discontinuity
idx_end   = MAX(idx_end, MIN(lmin_idx + 1, n-1))
! Boundary conditions
! Left boundary: on the flat section -> slope = 0
y_start  = masked_spike(idx_start)
dy_start = 0.0_r8
! Right boundary: on the spike profile -> central-difference slope
y_end  = masked_spike(idx_end)
dy_end = (masked_spike(idx_end+1) - masked_spike(idx_end-1)) &
       / (psi_N(idx_end+1)       - psi_N(idx_end-1))
! Generate the cubic Hermite patch and overwrite the junction window
dx_window = psi_N(idx_end) - psi_N(idx_start)
IF (ABS(dx_window) > 0.0_r8 .AND. idx_end > idx_start) THEN
  DO i = idx_start, idx_end
    t = (psi_N(i) - psi_N(idx_start)) / dx_window
    h00 =  2.0_r8*t**3 - 3.0_r8*t**2 + 1.0_r8   ! weight for y_start
    h10 =         t**3 - 2.0_r8*t**2 + t          ! weight for dy_start (scaled)
    h01 = -2.0_r8*t**3 + 3.0_r8*t**2              ! weight for y_end
    h11 =         t**3 -        t**2               ! weight for dy_end (scaled)
    y_patch = h00*y_start + h10*dx_window*dy_start &
            + h01*y_end   + h11*dx_window*dy_end
    masked_spike(i) = y_patch
  END DO
END IF
! =====================================================================
! 5. Fit parametrise_edge_jbs to j_BS over the edge-spike region.
!    Only performed when parameterized_spike output is requested.
! =====================================================================
IF (PRESENT(parameterized_spike)) THEN
  ! FWHM scan for LM fit initial guess
  half_max  = 0.5_r8 * peak_height
  left_idx  = lmin_idx   ! fallback: spike left boundary if no half-max crossing found before psi=0.7
  DO i = peak_idx-1, 1, -1
    IF (psi_N(i) < 0.7_r8) EXIT
    IF (j_BS(i) <= half_max) THEN
      left_idx = i
      EXIT
    END IF
  END DO
  right_idx = n       ! fallback: spike right boundary if no half-max crossing found before psi=1.0
  DO i = peak_idx+1, n
    IF (j_BS(i) <= half_max) THEN
      right_idx = i
      EXIT
    END IF
  END DO
  fwhm = psi_N(right_idx) - psi_N(left_idx)
  ! Build the (psi, j) fitting arrays over the spike region
  n_fit = COUNT(psi_N >= psi_N(lmin_idx))
  ALLOCATE(psi_fit(n_fit), j_fit(n_fit))
  n_fit = 0
  DO i = 1, n
    IF (psi_N(i) >= psi_N(lmin_idx)) THEN
      n_fit = n_fit + 1
      psi_fit(n_fit) = psi_N(i)
      j_fit(n_fit)   = j_BS(i)
    END IF
  END DO
  ! Bounds for LM fit (mirrors Python bootstrap.py analyze_bootstrap_edge_spike)
  ! amp bounds
  eps_tmp = MAX(1.0e-6_r8, 1.0e-3_r8 * ABS(peak_height))
  IF (peak_height == 0.0_r8) eps_tmp = 1.0e-6_r8
  amp_lo  = MIN(0.9995_r8 * peak_height, peak_height - eps_tmp)
  amp_hi  = MAX(1.0005_r8 * peak_height, peak_height + eps_tmp)
  ! offset bounds
  eps_tmp = MAX(1.0e-6_r8, 1.0e-3_r8 * ABS(lmin_j_BS))
  IF (lmin_j_BS == 0.0_r8) eps_tmp = 1.0e-6_r8
  off_lo  = MIN(0.99_r8 * lmin_j_BS, lmin_j_BS - eps_tmp)
  off_hi  = MAX(1.01_r8 * lmin_j_BS, lmin_j_BS + eps_tmp)
  ! y_sep upper bound: max(2|j_BS(n)|, |j_BS(n)|+1e-6, 1e-6)
  ysep_hi = MAX(2.0_r8 * ABS(j_BS(n)), ABS(j_BS(n)) + 1.0e-6_r8, 1.0e-6_r8)
  ! Assemble bound arrays
  lb(1) = amp_lo;          ub(1) = amp_hi
  lb(2) = 0.8_r8*peak_psi; ub(2) = MIN(1.2_r8*peak_psi, psi_N(n))
  lb(3) = 0.0_r8;          ub(3) = 0.33_r8
  lb(4) = off_lo;           ub(4) = off_hi
  lb(5) = -50.0_r8;         ub(5) = 50.0_r8
  lb(6) = 0.0_r8;           ub(6) = ysep_hi
  lb(7) = 0.001_r8;         ub(7) = 0.2_r8
  ! Initial parameter guess (matches Python p0)                          ! amp
  p0(1) = peak_height                          ! amplitude
  p0(2) = peak_psi                             ! center
  p0(3) = fwhm / 2.355_r8                     ! width (sigma)
  p0(4) = lmin_j_BS                            ! offset
  p0(5) = 1.0_r8                              ! sk
  p0(6) = MAX(0.0_r8, j_BS(n))            ! y_sep
  p0(7) = 0.05_r8                             ! blend_width
  ! Levenberg-Marquardt fit (pass bounds for internal tan-transform)
  CALL curve_fit_edge_jbs(n_fit, psi_fit, j_fit, p0, popt, lb, ub)
  DEALLOCATE(psi_fit, j_fit)
  amp_fit    = popt(1)
  center_fit = popt(2)
  width_fit  = popt(3)
  offset_fit = popt(4)
  sk_fit     = popt(5)
  y_sep_fit  = popt(6)
  bw_fit     = popt(7)
  ! Evaluate fitted profile on the full n-point grid
  CALL parametrise_edge_jbs(n, psi_N, amp_fit, center_fit, width_fit, &
      offset_fit, sk_fit, y_sep_fit, bw_fit, 1.5_r8, parameterized_spike)
  ! Optional verbose output for diagnostics if gseq%boot_ops%diagnose_bs
  IF (PRESENT(diagnose) .AND. diagnose) THEN
    WRITE(*,'(A,7(A,ES12.4))') '  [edge_spike_fit]', &
      ' amp=',    amp_fit,    ' center=', center_fit, ' width=',  width_fit, &
      ' offset=', offset_fit, ' sk=',     sk_fit,     ' y_sep=',  y_sep_fit, &
      ' bw=',     bw_fit
    WRITE(*,'(A)') '  [edge_spike_profile] i  psi_N(std)    parameterized_spike[A/m2]'
    DO i = 1, n
      WRITE(*,'(A,I4,2ES15.5)') '  ', i, psi_N(i), parameterized_spike(i)
    END DO
  END IF
END IF
END SUBROUTINE analyze_bootstrap_edge_spike
!------------------------------------------------------------------------------
!> Redl 2021 bootstrap current formula.
!>
!> Translates Python bootstrap.py:redl_bootstrap with fixed settings:
!>   - use_legacy_L34 = .FALSE.  (L34 = L31)
!>   - use_sign_q     = .TRUE.
!>   - formula_form   = 'jboot1'
!>   - nu_e_star and nu_i_star are passed in pre-computed (no internal fallback)
!>
!> Reference: Redl et al., Phys. Plasmas 28, 022502 (2021)
!>
!> @param n           Number of flux surfaces
!> @param Te          Electron temperature [eV]
!> @param Ti          Ion temperature [eV]
!> @param ne          Electron density [m^-3]
!> @param ni          Ion density [m^-3]
!> @param pe          Electron pressure [Pa]
!> @param pi          Ion pressure [Pa]
!> @param Zeff        Effective charge
!> @param q           Safety factor
!> @param eps         Inverse aspect ratio
!> @param fT          Trapped particle fraction
!> @param I_psi       Toroidal current function F = R*Bt [T*m]
!> @param dT_e_dpsi   d(Te)/d(psi) [eV/Wb]
!> @param dT_i_dpsi   d(Ti)/d(psi) [eV/Wb]
!> @param dn_e_dpsi   d(ne)/d(psi) [m^-3/Wb]
!> @param dn_i_dpsi   d(ni)/d(psi) [m^-3/Wb]
!> @param ln_lambda_e Electron Coulomb logarithm
!> @param ln_lambda_ii Ion Coulomb logarithm
!> @param nu_e_star   Electron collisionality (pre-computed)
!> @param nu_i_star   Ion collisionality (pre-computed)
!> @param avg_j_bootstrap_times_B Output: <j_{bs,parallel}*B> [same units as -I_psi * pe * d/dpsi]
!------------------------------------------------------------------------------
SUBROUTINE redl_bootstrap(n, Te, Ti, ne, ni, pe, pi, Zeff, q, eps, fT, I_psi, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    ln_lambda_e, ln_lambda_ii, nu_e_star, nu_i_star, avg_j_bootstrap_times_B)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: Te(n), Ti(n), ne(n), ni(n)
REAL(r8), INTENT(in) :: pe(n), pi(n), Zeff(n)
REAL(r8), INTENT(in) :: q(n), eps(n), fT(n), I_psi(n)
REAL(r8), INTENT(in) :: dT_e_dpsi(n), dT_i_dpsi(n)
REAL(r8), INTENT(in) :: dn_e_dpsi(n), dn_i_dpsi(n)
REAL(r8), INTENT(in) :: ln_lambda_e(n), ln_lambda_ii(n)
REAL(r8), INTENT(in) :: nu_e_star(n), nu_i_star(n)
REAL(r8), INTENT(out) :: avg_j_bootstrap_times_B(n)
!---
REAL(r8) :: R_pe(n)
REAL(r8) :: ft_31_d1(n), ft_31_d2(n), X31(n), L31(n), L34(n), dZ(n)
REAL(r8) :: dee_2(n), dee_3(n), X32_ee(n), F32_ee(n)
REAL(r8) :: dei_2(n), dei_3(n), X32_ei(n), F32_ei(n), L32(n)
REAL(r8) :: alpha0(n), alpha(n)
REAL(r8) :: dp_dpsi(n), bra1(n), bra2(n), bra3(n)
REAL(r8), PARAMETER :: EC = 1.602176634e-19_r8
!---
R_pe = pe / (pe + pi)
! =====================================================================
! L31 (Redl Eqs. 10-11)
! =====================================================================
ft_31_d1 = (0.67_r8 * (1.0_r8 - 0.7_r8*fT) * SQRT(nu_e_star)) &
         / (0.56_r8 + 0.44_r8*Zeff)
ft_31_d2 = ((0.52_r8 + 0.086_r8*SQRT(nu_e_star)) &
          * (1.0_r8 + 0.87_r8*fT) * nu_e_star) &
         / (1.0_r8 + 1.13_r8*SQRT(MAX(Zeff - 1.0_r8, 0.0_r8)))
X31 = fT / (1.0_r8 + ft_31_d1 + ft_31_d2)
dZ = Zeff**1.2_r8 - 0.71_r8
L31 = (1.0_r8 + 0.15_r8/dZ)*X31 &
    - (0.22_r8/dZ)*X31**2 &
    + (0.01_r8/dZ)*X31**3 &
    + (0.06_r8/dZ)*X31**4
! L34 = L31 (use_legacy_L34 = .FALSE.)
L34 = L31
! =====================================================================
! L32 (Redl Eqs. 12-16)
! =====================================================================
dee_2 = (0.23_r8 * (1.0_r8 - 0.96_r8*fT) * SQRT(nu_e_star)) &
      / SQRT(Zeff)
dee_3 = (0.13_r8 * (1.0_r8 - 0.38_r8*fT) * nu_e_star / (Zeff**2)) &
      * (SQRT(1.0_r8 + 2.0_r8*SQRT(MAX(Zeff - 1.0_r8, 0.0_r8))) &
       + fT**2 * SQRT((0.075_r8 + 0.25_r8*(Zeff - 1.0_r8)**2) * nu_e_star))
X32_ee = fT / (1.0_r8 + dee_2 + dee_3)
F32_ee = ( (0.1_r8 + 0.6_r8*Zeff) &
         / (Zeff*(0.77_r8 + 0.63_r8*(1.0_r8 + (Zeff - 1.0_r8)**1.1_r8))) &
         * (X32_ee - X32_ee**4) &
         + 0.7_r8/(1.0_r8 + 0.2_r8*Zeff) &
         * (X32_ee**2 - X32_ee**4 - 1.2_r8*(X32_ee**3 - X32_ee**4)) &
         + 1.3_r8/(1.0_r8 + 0.5_r8*Zeff) * X32_ee**4 )
dei_2 = (0.87_r8 * (1.0_r8 + 0.39_r8*fT) * SQRT(nu_e_star)) &
      / (1.0_r8 + 2.95_r8*(Zeff - 1.0_r8)**2)
dei_3 = 1.53_r8 * (1.0_r8 - 0.37_r8*fT) * nu_e_star &
      * (2.0_r8 + 0.375_r8*(Zeff - 1.0_r8))
X32_ei = fT / (1.0_r8 + dei_2 + dei_3)
F32_ei = ( -(0.4_r8 + 1.93_r8*Zeff) / (Zeff*(0.8_r8 + 0.6_r8*Zeff)) &
           * (X32_ei - X32_ei**4) &
           + 5.5_r8/(1.5_r8 + 2.0_r8*Zeff) &
           * (X32_ei**2 - X32_ei**4 - 0.8_r8*(X32_ei**3 - X32_ei**4)) &
           - 1.3_r8/(1.0_r8 + 0.5_r8*Zeff) * X32_ei**4 )
L32 = F32_ee + F32_ei
! =====================================================================
! Alpha (Redl Eqs. 20-21)
! =====================================================================
alpha0 = -(0.62_r8 + 0.055_r8*(Zeff - 1.0_r8)) &
        / (0.53_r8 + 0.17_r8*(Zeff - 1.0_r8)) &
        * (1.0_r8 - fT) &
        / (1.0_r8 - (0.31_r8 - 0.065_r8*(Zeff - 1.0_r8))*fT - 0.25_r8*fT**2)
alpha = ((alpha0 + 0.7_r8*Zeff*SQRT(fT)*SQRT(nu_i_star)) &
       / (1.0_r8 + 0.18_r8*SQRT(nu_i_star)) &
       - 0.002_r8*nu_i_star**2*fT**6) &
      / (1.0_r8 + 0.004_r8*nu_i_star**2*fT**6)
! =====================================================================
! Assemble: jboot1 form with use_sign_q = .TRUE.
! =====================================================================
dp_dpsi = (ne*dT_e_dpsi + Te*dn_e_dpsi &
         + ni*dT_i_dpsi + Ti*dn_i_dpsi) * EC
bra1 = L31 * dp_dpsi / pe
bra2 = L32 * dT_e_dpsi / Te
bra3 = L34 * alpha * (1.0_r8 - R_pe) / R_pe * dT_i_dpsi / Ti
avg_j_bootstrap_times_B = -I_psi * pe * (bra1 + bra2 + bra3)
! Apply sign(q)
avg_j_bootstrap_times_B = avg_j_bootstrap_times_B * SIGN(1.0_r8, q)
END SUBROUTINE redl_bootstrap
end module grad_shaf_bootstrap