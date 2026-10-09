!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_bootstrap.F90
!
!> Physics-scripts used in the bootstrap calculation inside the Grad-Sharfranov solve
!!
!! @authors Stuart Benjamin, fortranisation of Daniel Burgess' bootstrap.py
!! @date June 2026
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module grad_shaf_bootstrap
use oft_base
use oft_gs, only: gs_equil, flux_func, gsinv_interp, gs_factory, gs_psi2r, &
 gs_itor_nl, qprof_trace_tol
use oft_gs_profiles, only: torflux_qgeom_backend
use oft_gs_cutcell, only: gs_sauter_cutcell, gs_qgeom_cutcell
use oft_lag_basis, only: oft_blag_geval
use oft_mesh_type, only: bmesh_findcell
use oft_blag_operators, only: oft_lag_brinterp
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use grad_shaf_prof_phys, only: eval_jtor_imas, gs_ravgs, gs_flux_int, &
  gs_flux_cumint, jphi_update, jphi_copy, jphi_flux_func, jphi_psi_nodes
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
  LOGICAL :: taper_edge_jBS = .FALSE. !< Smooth taper of toroidal current (inductive, bootstrap and fixed) to zero at plasma edge (guards against numerical issues at the separatrix)
  REAL(r8) :: taper_edge_psi0 = 0.999_r8 !< psi_N (standard: 0=axis, 1=LCFS) where edge taper begins
  INTEGER(i4) :: taper_edge_shape = 2 !< Edge taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
  REAL(r8) :: saw_q_s = 0.0_r8 !< Sawtooth reset q on axis (q_s); 0 = off (see @ref saw_reset_1d)
  REAL(r8) :: saw_dq = 0.03_r8 !< Reset regions end where q_base crosses saw_q_s + saw_dq
  REAL(r8) :: saw_tol = 1.0e-4_r8 !< Rel. change in j_saw (vs RMS total) below which j_saw is frozen
  REAL(r8) :: saw_relax = 1.0_r8 !< Under-relaxation of the j_saw update (1 = none)
  REAL(r8) :: saw_ramp = 0.01_r8 !< Reset weight ramps over this q deficit below saw_q_s (0 = hard trigger)
  INTEGER(i4) :: saw_rule = 2 !< Reset rule: 1 = fuse (axis to rho_m, FUSE saw_crash!), 2 = local (each dip, two-sided)
END TYPE boot_ops
!------------------------------------------------------------------------------
!> Cached current profiles produced by the last call to @ref jphi_bs_update.
!!
!! All arrays are allocated/overwritten on every call to jphi_bs_update
!! and remain valid until the next call or until the parent @ref gs_equil is
!! destroyed.  Units are A/m² (before the mu0 normalisation used internally).
!------------------------------------------------------------------------------
TYPE :: boot_profs
  REAL(r8), POINTER, DIMENSION(:) :: psi_n => NULL() !< Normalised psi_N values for these current profiles in OFT convention (0=LCFS, 1=axis); index 0 is the LCFS boundary
  REAL(r8), POINTER, DIMENSION(:) :: j_bs_raw => NULL() !< Redl PoP 2021 bootstrap as TokaMaker jphi, <j_BS.B> F<1/R>/<B^2> + P'(<R> - F^2<1/R>/<B^2>), before isolation/scaling [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: jdotb_bs_raw => NULL() !< Redl PoP 2021 <j_BS.B> on the j_bs_raw grid [T A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: total_j_phi => NULL() !< Total toroidal current density = j_ind_final + j_bs_final + jphi_fixed [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: j_ind_final => NULL() !< Input jphi, re-scaled & optionally tapered [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: j_bs_final => NULL() !< Bootstrap current density, optionally isolated/parametrised/tapered [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: jphi_fixed => NULL() !< Fixed (non-rescaled) toroidal current density, optionally tapered [A/m²]
  REAL(r8), POINTER, DIMENSION(:) :: j_saw => NULL() !< Sawtooth current = input jphi_saw + q reset redistribution [A/m²]
  REAL(r8) :: saw_rho_m = 0.0_r8 !< Sawtooth radius (rho_tor_norm): q_s + dq crossing beyond the outermost q_base < q_s (0 = none)
  REAL(r8) :: saw_rho_out = 0.0_r8 !< Outer end of the last reset (rho_tor_norm; 0 = none)
  INTEGER(i4) :: saw_n_dips = 0 !< Regions reset at the last reset (rule 1: separate q_base < saw_q_s regions)
END TYPE boot_profs
!------------------------------------------------------------------------------
!> Jphi flux function type for bootstrap current calculation
!------------------------------------------------------------------------------
type, extends(jphi_flux_func) :: jphi_bs_flux_func
  real(8) :: alpha_last = 1.d0 !< Alpha (input Jphi rescaling factor) from previous NL iteration
  logical :: freeze_j_BS = .FALSE. !< Set .TRUE. once j_BS stagnates for 2 steps; skips Sauter call (big speedup)
  real(8) :: djBS_stol = 1.0e-3_r8 !< RMS tolerance stagnation value: reset no-improve counter if above this
  real(8) :: dalpha_warn_tol = 1.0e-6_r8 !< Relative (|dalpha|/|alpha|): no alpha stall warning below this
  integer(4) :: djBS_no_improve = 0   !< Consecutive steps with non-decreasing djBS
  real(8) :: djBS_min = huge(1.0d0)  !< Running minimum djBS seen so far
  integer(4) :: dalpha_no_improve = 0 !< Consecutive steps with non-decreasing dalpha
  real(8) :: dalpha_last = huge(1.0d0) !< dalpha of the previous update (stall diagnostic)
  real(8), pointer, dimension(:) :: j_BS_last => NULL() !< j_BS profile from previous NL iteration (for freeze check)
  logical :: freeze_saw = .FALSE. !< Set .TRUE. once j_saw converges or stagnates, after j_BS is frozen
  integer(4) :: djsaw_no_improve = 0 !< Consecutive steps with non-decreasing djsaw
  real(8) :: djsaw_min = huge(1.0d0) !< Running minimum djsaw
  real(8), pointer, dimension(:) :: j_saw_last => NULL() !< j_saw from previous NL iteration (mu0*A/m²)
  real(8), pointer, dimension(:) :: jtot_last => NULL() !< Total jphi from previous NL iteration (mu0*A/m²), the current of the equilibrium being traced
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
  real(8) :: bmax = 0.d0        !< Maximum |B| on surface (sampled during the trace)
  real(8) :: rmax_surf = -1.d30 !< Maximum R on surface (sampled during the trace)
  real(8) :: rmin_surf =  1.d30 !< Minimum R on surface (sampled during the trace)
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
REAL(r8) :: sauter_wtime = 0.d0 !< Time of the last @ref sauter_fc call [s] (cut-cell wall time plus tracing time summed over threads)
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
  CALL hdf5_write(self%boot_ops%saw_q_s,filename,path//'/BOOT_OPS/SAW_Q_S')
  CALL hdf5_write(self%boot_ops%saw_dq,filename,path//'/BOOT_OPS/SAW_DQ')
  CALL hdf5_write(self%boot_ops%saw_tol,filename,path//'/BOOT_OPS/SAW_TOL')
  CALL hdf5_write(self%boot_ops%saw_relax,filename,path//'/BOOT_OPS/SAW_RELAX')
  CALL hdf5_write(self%boot_ops%saw_ramp,filename,path//'/BOOT_OPS/SAW_RAMP')
  CALL hdf5_write(self%boot_ops%saw_rule,filename,path//'/BOOT_OPS/SAW_RULE')
END IF
!---Save cached bootstrap current profiles (only when available)
IF(ASSOCIATED(self%boot_profs%total_j_phi).OR.ASSOCIATED(self%boot_profs%j_bs_raw))THEN
  CALL hdf5_create_group(filename,path//'/BOOT_PROFS')
  IF(ASSOCIATED(self%boot_profs%total_j_phi))THEN
    CALL hdf5_write(self%boot_profs%psi_n,filename,path//'/BOOT_PROFS/PSI_N')
    CALL hdf5_write(self%boot_profs%total_j_phi,filename,path//'/BOOT_PROFS/TOTAL_J_PHI')
    CALL hdf5_write(self%boot_profs%j_ind_final,filename,path//'/BOOT_PROFS/J_IND_FINAL')
    CALL hdf5_write(self%boot_profs%j_bs_final,filename,path//'/BOOT_PROFS/J_BS_FINAL')
    IF(ASSOCIATED(self%boot_profs%jphi_fixed)) &
      CALL hdf5_write(self%boot_profs%jphi_fixed,filename,path//'/BOOT_PROFS/JPHI_FIXED')
    IF(ASSOCIATED(self%boot_profs%j_saw)) &
      CALL hdf5_write(self%boot_profs%j_saw,filename,path//'/BOOT_PROFS/J_SAW')
    IF(ASSOCIATED(self%boot_profs%j_bs_raw)) &
      CALL hdf5_write(self%boot_profs%j_bs_raw,filename,path//'/BOOT_PROFS/J_BS_RAW')
    IF(ASSOCIATED(self%boot_profs%jdotb_bs_raw)) &
      CALL hdf5_write(self%boot_profs%jdotb_bs_raw,filename,path//'/BOOT_PROFS/JDOTB_BS_RAW')
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
  IF(hdf5_field_exist(filename,path//'/BOOT_OPS/SAW_Q_S'))THEN
    CALL hdf5_read(self%boot_ops%saw_q_s,filename,path//'/BOOT_OPS/SAW_Q_S',success=success)
    CALL hdf5_read(self%boot_ops%saw_dq,filename,path//'/BOOT_OPS/SAW_DQ',success=success)
    CALL hdf5_read(self%boot_ops%saw_tol,filename,path//'/BOOT_OPS/SAW_TOL',success=success)
    CALL hdf5_read(self%boot_ops%saw_relax,filename,path//'/BOOT_OPS/SAW_RELAX',success=success)
    CALL hdf5_read(self%boot_ops%saw_ramp,filename,path//'/BOOT_OPS/SAW_RAMP',success=success)
    CALL hdf5_read(self%boot_ops%saw_rule,filename,path//'/BOOT_OPS/SAW_RULE',success=success)
  END IF
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
    IF(ASSOCIATED(self%boot_profs%jphi_fixed))DEALLOCATE(self%boot_profs%jphi_fixed)
    IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/JPHI_FIXED'))THEN
      ALLOCATE(self%boot_profs%jphi_fixed(0:SIZE(self%boot_profs%total_j_phi)-1))
      CALL hdf5_read(self%boot_profs%jphi_fixed,filename,path//'/BOOT_PROFS/JPHI_FIXED',success=success)
    END IF
    IF(ASSOCIATED(self%boot_profs%j_saw))DEALLOCATE(self%boot_profs%j_saw)
    IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/J_SAW'))THEN
      ALLOCATE(self%boot_profs%j_saw(0:SIZE(self%boot_profs%total_j_phi)-1))
      CALL hdf5_read(self%boot_profs%j_saw,filename,path//'/BOOT_PROFS/J_SAW',success=success)
    END IF
    IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/J_BS_RAW'))THEN
      CALL hdf5_field_get_sizes(filename,path//'/BOOT_PROFS/J_BS_RAW',ndims,dim_sizes)
      IF(ASSOCIATED(self%boot_profs%j_bs_raw))DEALLOCATE(self%boot_profs%j_bs_raw)
      ALLOCATE(self%boot_profs%j_bs_raw(0:dim_sizes(1)-1))
      DEALLOCATE(dim_sizes)
      CALL hdf5_read(self%boot_profs%j_bs_raw,filename,path//'/BOOT_PROFS/J_BS_RAW',success=success)
    END IF
    IF(ASSOCIATED(self%boot_profs%jdotb_bs_raw))DEALLOCATE(self%boot_profs%jdotb_bs_raw)
    IF(hdf5_field_exist(filename,path//'/BOOT_PROFS/JDOTB_BS_RAW'))THEN
      CALL hdf5_field_get_sizes(filename,path//'/BOOT_PROFS/JDOTB_BS_RAW',ndims,dim_sizes)
      ALLOCATE(self%boot_profs%jdotb_bs_raw(0:dim_sizes(1)-1))
      DEALLOCATE(dim_sizes)
      CALL hdf5_read(self%boot_profs%jdotb_bs_raw,filename,path//'/BOOT_PROFS/JDOTB_BS_RAW',success=success)
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
    new%freeze_j_BS = self%freeze_j_BS
    new%djBS_stol = self%djBS_stol
    new%dalpha_warn_tol = self%dalpha_warn_tol
    new%djBS_no_improve = self%djBS_no_improve
    new%djBS_min = self%djBS_min
    new%dalpha_no_improve = self%dalpha_no_improve
    new%dalpha_last = self%dalpha_last
    new%boot_ops = self%boot_ops
    IF(ASSOCIATED(self%j_BS_last)) ALLOCATE(new%j_BS_last, SOURCE=self%j_BS_last)
    new%freeze_saw = self%freeze_saw
    new%djsaw_no_improve = self%djsaw_no_improve
    new%djsaw_min = self%djsaw_min
    IF(ASSOCIATED(self%j_saw_last)) ALLOCATE(new%j_saw_last, SOURCE=self%j_saw_last)
    IF(ASSOCIATED(self%jtot_last)) ALLOCATE(new%jtot_last, SOURCE=self%jtot_last)
    new%boot_profs%saw_rho_m = self%boot_profs%saw_rho_m
    new%boot_profs%saw_rho_out = self%boot_profs%saw_rho_out
    new%boot_profs%saw_n_dips = self%boot_profs%saw_n_dips
    IF(ASSOCIATED(self%boot_profs%psi_n))ALLOCATE(new%boot_profs%psi_n,SOURCE=self%boot_profs%psi_n)
    IF(ASSOCIATED(self%boot_profs%j_bs_raw))ALLOCATE(new%boot_profs%j_bs_raw,SOURCE=self%boot_profs%j_bs_raw)
    IF(ASSOCIATED(self%boot_profs%jdotb_bs_raw))ALLOCATE(new%boot_profs%jdotb_bs_raw,SOURCE=self%boot_profs%jdotb_bs_raw)
    IF(ASSOCIATED(self%boot_profs%total_j_phi))ALLOCATE(new%boot_profs%total_j_phi,SOURCE=self%boot_profs%total_j_phi)
    IF(ASSOCIATED(self%boot_profs%j_bs_final))ALLOCATE(new%boot_profs%j_bs_final,SOURCE=self%boot_profs%j_bs_final)
    IF(ASSOCIATED(self%boot_profs%j_ind_final))ALLOCATE(new%boot_profs%j_ind_final,SOURCE=self%boot_profs%j_ind_final)
    IF(ASSOCIATED(self%boot_profs%jphi_fixed))ALLOCATE(new%boot_profs%jphi_fixed,SOURCE=self%boot_profs%jphi_fixed)
    IF(ASSOCIATED(self%boot_profs%j_saw))ALLOCATE(new%boot_profs%j_saw,SOURCE=self%boot_profs%j_saw)
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
IF(ASSOCIATED(self%j_BS_last))DEALLOCATE(self%j_BS_last)
IF(ASSOCIATED(self%j_saw_last))DEALLOCATE(self%j_saw_last)
IF(ASSOCIATED(self%jtot_last))DEALLOCATE(self%jtot_last)
!---Destroy cached bootstrap current profiles
IF(ASSOCIATED(self%boot_profs%j_bs_raw))DEALLOCATE(self%boot_profs%j_bs_raw)
IF(ASSOCIATED(self%boot_profs%jdotb_bs_raw))DEALLOCATE(self%boot_profs%jdotb_bs_raw)
IF(ASSOCIATED(self%boot_profs%total_j_phi))DEALLOCATE(self%boot_profs%total_j_phi)
IF(ASSOCIATED(self%boot_profs%j_bs_final))DEALLOCATE(self%boot_profs%j_bs_final)
IF(ASSOCIATED(self%boot_profs%j_ind_final))DEALLOCATE(self%boot_profs%j_ind_final)
IF(ASSOCIATED(self%boot_profs%jphi_fixed))DEALLOCATE(self%boot_profs%jphi_fixed)
IF(ASSOCIATED(self%boot_profs%j_saw))DEALLOCATE(self%boot_profs%j_saw)
IF(ASSOCIATED(self%boot_profs%psi_n))DEALLOCATE(self%boot_profs%psi_n)
end subroutine jphi_bs_delete
!---------------------------------------------------------------------------------
!> Update F*F' profile from inductive Jphi coupled with bootstrap current.
!>
!> Each call (one NL iteration):
!>   1. Pressure scale.
!>   2. Evaluate fixed current jphi_fixed and bootstrap current j_BS on self%x (or reuse cache if frozen);
!>      <R>, <1/R>, <1/R^2> on the nodes from its Sauter pass, else (and on the LCFS) from gs_ravgs.
!>   3. Apply edge taper to j_BS, jphi_ind, jphi_fixed and the input sawtooth current component-wise.
!>   4. Solve analytically for alpha: the exact I_p (gs_flux_int of eval_jtor_imas, eq. A9c) is
!>      affine in alpha, so two evaluations (alpha=0, alpha=1) give
!>      alpha = (Ip_target - Ip_lo)/(Ip_hi - Ip_lo).
!>   4a. Sawtooth reset (saw_q_s > 0): j_saw = jphi_saw + the Ip-neutral current that resets the q of
!>      base = alpha*jphi_ind + j_BS + jphi_fixed + jphi_saw above saw_q_s (@ref saw_redistribute).
!>   5. Assemble jphi_total = alpha*jphi_ind + j_BS + jphi_fixed + j_saw; compute F*F' knots.
!>   6. Diagnostics (if diagnose_bs is set), including wall times of the steps above.
!---------------------------------------------------------------------------------
SUBROUTINE jphi_bs_update(self, gseq)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: pscale, pprime
REAL(r8), ALLOCATABLE :: jtor(:)  !< IMAS-convention current for the exact I_p measure (eval_jtor_imas)
REAL(r8), ALLOCATABLE :: xpsi(:) !< Node locations in normalized poloidal flux
REAL(r8), ALLOCATABLE :: psi_n(:) !< LCFS and node locations [0, xpsi]
REAL(r8), ALLOCATABLE :: ravgs(:,:) !< <R>, <1/R>, <1/R^2> at psi_n
! Bootstrap arrays (on self%x grid)
REAL(r8), ALLOCATABLE :: j_BS(:)
! Optional edge-spike workspace (only allocated when isolate_edge_jBS or parameterize_jBS is set)
REAL(r8), ALLOCATABLE :: j_spike_tmp(:)
REAL(r8), ALLOCATABLE :: j_spike_mask_tmp(:)  !< Raw masked (pre-fit) spike profile
! Working arrays on self%x grid
REAL(r8), ALLOCATABLE :: jphi_total(:)
REAL(r8), ALLOCATABLE :: jphi_ind(:)  !< Tapered copy of self%jphi (= self%jphi when taper off)
REAL(r8), ALLOCATABLE :: jphi_fixed(:)  !< Fixed current from gseq%jphi_fixed (A/m², mu0*A/m² from step 3; 0 if unset)
! Sawtooth current (mu0*A/m² from step 3): input, value entering the alpha solve, final
REAL(r8), ALLOCATABLE :: j_saw_in(:), j_saw_cur(:), j_saw(:), dj_saw(:)
REAL(r8), ALLOCATABLE :: cflds(:,:), cum(:,:)
LOGICAL :: saw_on, do_saw
REAL(r8) :: djsaw
! Alpha-solve scalars
REAL(r8) :: alpha, ip_target, ip_ind, ip_result_lo, ip_result_hi, dalpha
! Relative change in bootstrap current for freeze check
REAL(r8) :: djBS
LOGICAL :: bs_frozen !< j_BS reused from the last update (no Sauter pass)
! Diagnostic I_p comparison
REAL(r8) :: itor_nl, itor_flint
REAL(r8) :: wt(5) !< Wall times [s]: <R> averages, bootstrap, alpha solve, update total, start
CHARACTER(len=256) :: char_buf
!--- First-iteration runs with no bootstrap current
IF(.NOT. gseq%skip_targets) THEN
  CALL jphi_update(self,gseq)
  RETURN
ENDIF
!---
self%plasma_bounds = gseq%plasma_bounds
CALL jphi_psi_nodes(self,xpsi)
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
!--- 1. Pressure scale.
ALLOCATE(jtor(0:self%npsi), psi_n(0:self%npsi), ravgs(0:self%npsi,3))
psi_n = [0.0_r8, xpsi]
wt=0.d0; wt(5)=omp_get_wtime()
saw_on = self%boot_ops%saw_q_s > 0.0_r8
do_saw = saw_on .AND. (.NOT.self%freeze_saw) .AND. ASSOCIATED(self%jtot_last)
djsaw = 0.0_r8
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
IF(ASSOCIATED(gseq%P_ani)) &
  CALL oft_abort('Jphi profiles do not support anisotropic pressure', &
                 'jphi_bs_update',__FILE__)
pscale = gseq%P%f(gseq%plasma_bounds(2))
pscale = gseq%pax_target / pscale
!--- 2. Fixed current [A/m²] and bootstrap current on self%x grid.
ALLOCATE(jphi_fixed(0:self%npsi))
jphi_fixed = 0.0_r8
IF(ASSOCIATED(gseq%jphi_fixed))THEN
  jphi_fixed(0) = gseq%jphi_fixed%fp(0.0_r8)
  DO i = 1, self%npsi
    jphi_fixed(i) = gseq%jphi_fixed%fp(xpsi(i))
  END DO
END IF
ALLOCATE(j_saw_in(0:self%npsi))
j_saw_in = 0.0_r8
IF(ASSOCIATED(gseq%jphi_saw))THEN
  j_saw_in(0) = gseq%jphi_saw%fp(0.0_r8)
  DO i = 1, self%npsi
    j_saw_in(i) = gseq%jphi_saw%fp(xpsi(i))
  END DO
END IF
ALLOCATE(j_BS(0:self%npsi))
bs_frozen = self%freeze_j_BS .AND. ASSOCIATED(self%j_BS_last)
IF(bs_frozen) THEN
  !--- Frozen: reuse cached j_BS.
  j_BS = self%j_BS_last
  djBS = 0.0_r8
ELSE
  !--- Not frozen: run full bootstrap calculation (Sauter).
  wt(2)=omp_get_wtime()
  IF (self%boot_ops%isolate_edge_jBS .OR. self%boot_ops%parameterize_jBS) THEN
    ALLOCATE(j_spike_tmp(0:self%npsi), j_spike_mask_tmp(0:self%npsi))
    CALL calculate_bootstrap(self, gseq, self%npsi, xpsi, j_BS, &
        isolate_edge_jBS=self%boot_ops%isolate_edge_jBS, &
        parameterize_jBS=self%boot_ops%parameterize_jBS, &
        scale_jBS=self%boot_ops%scale_jBS, &
        j_spike=j_spike_tmp, j_spike_masked=j_spike_mask_tmp, ravgs=ravgs(1:,:))
    IF (self%boot_ops%diagnose_bs) THEN
      IF (self%boot_ops%parameterize_jBS) THEN
        WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS(bulk)[A/m2]  j_spike[A/m2]   j_spike_masked[A/m2]  jphi[A/m2]  jphi_fixed[A/m2]'
        DO i = 1, self%npsi
          WRITE(*,'(A,I4,6ES15.5)') '  ', i, xpsi(i), j_BS(i), j_spike_tmp(i), j_spike_mask_tmp(i), self%jphi(i), jphi_fixed(i)
        END DO
      ELSE
        WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS(bulk)[A/m2]  j_spike[A/m2]   jphi[A/m2]  jphi_fixed[A/m2]'
        DO i = 1, self%npsi
          WRITE(*,'(A,I4,5ES15.5)') '  ', i, xpsi(i), j_BS(i), j_spike_tmp(i), self%jphi(i), jphi_fixed(i)
        END DO
      END IF
    END IF
    j_BS = j_spike_tmp
    DEALLOCATE(j_spike_tmp, j_spike_mask_tmp)
  ELSE
    CALL calculate_bootstrap(self, gseq, self%npsi, xpsi, j_BS, ravgs=ravgs(1:,:))
    j_BS = j_BS * self%boot_ops%scale_jBS
    IF(self%boot_ops%diagnose_bs)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS[A/m2]      jphi[A/m2]  jphi_fixed[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,4ES15.5)') '  ', i, xpsi(i), j_BS(i), self%jphi(i), jphi_fixed(i)
      END DO
    END IF
  END IF
  !   calculate_bootstrap returns j_BS in raw A/m², multiply by mu0.
  j_BS = j_BS * mu0
  !--- 2a. Freeze check: freeze j_BS if RMS change drops below tol, or stagnates.
  IF(ASSOCIATED(self%j_BS_last)) THEN
    djBS = SQRT(SUM((j_BS - self%j_BS_last)**2) / REAL(self%npsi+1,r8)) / &
           MAX(SQRT(SUM(j_BS**2) / REAL(self%npsi+1,r8)), 1.0e-30_r8)
    IF(djBS < self%boot_ops%djBS_tol) THEN
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
  IF(.NOT.ASSOCIATED(self%j_BS_last)) ALLOCATE(self%j_BS_last(0:self%npsi))
  self%j_BS_last = j_BS
  wt(2)=omp_get_wtime()-wt(2)
END IF
!--- 2b. <R>, <1/R>, <1/R^2> (I_p measure, F*F' map) not set by the Sauter pass
wt(1)=omp_get_wtime()
IF(bs_frozen)THEN
  CALL gs_ravgs(gseq, self%npsi+1, psi_n, ravgs)
ELSE
  CALL gs_ravgs(gseq, 1, psi_n(0:0), ravgs(0:0,:))
END IF
wt(1)=omp_get_wtime()-wt(1)
!--- 3. Apply edge taper to j_BS, jphi_ind, jphi_fixed and j_saw_in (the last two first converted to mu0*A/m²).
!   self%j_BS_last caches the un-tapered j_BS so freeze comparisons track physics.
!   taper_edge_psi0 is in standard convention (0=axis,1=LCFS);
!   threshold in OFT convention (0=LCFS,1=axis) is (1 - taper_edge_psi0).
ALLOCATE(jphi_ind(0:self%npsi))
jphi_ind = [self%j0, self%jphi]
jphi_fixed = jphi_fixed * mu0
j_saw_in = j_saw_in * mu0
IF (self%boot_ops%taper_edge_jBS) THEN
  CALL apply_edge_taper(self%npsi+1, psi_n, j_BS, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
  CALL apply_edge_taper(self%npsi+1, psi_n, jphi_ind, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
  CALL apply_edge_taper(self%npsi+1, psi_n, jphi_fixed, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
  CALL apply_edge_taper(self%npsi+1, psi_n, j_saw_in, &
                        1.0_r8 - self%boot_ops%taper_edge_psi0, &
                        self%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
END IF
ALLOCATE(jphi_total(0:self%npsi))
ALLOCATE(j_saw_cur(0:self%npsi), j_saw(0:self%npsi))
!   Sawtooth current entering the alpha solve: last iterate (or the input until the first reset)
j_saw_cur = j_saw_in
IF(saw_on .AND. ASSOCIATED(self%j_saw_last))j_saw_cur = self%j_saw_last
!   One quadrature pass for the alpha solve and the reset's enclosed currents: lo (alpha=0 with the
!   input saw), ind (inductive, linear), dsaw (j_saw_cur - input, linear), eq (last total), area.
IF(do_saw)THEN
  ALLOCATE(cflds(self%npsi+1,5), cum(self%npsi+1,5))
  jphi_total = j_BS + jphi_fixed + j_saw_in
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, cflds(:,1))
  jphi_total = jphi_ind + j_BS + jphi_fixed + j_saw_in
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  cflds(:,2) = jtor - cflds(:,1)
  jphi_total = j_BS + jphi_fixed + j_saw_cur
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  cflds(:,3) = jtor - cflds(:,1)
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, self%jtot_last, pscale, cflds(:,4))
  cflds(:,5) = 1.0_r8
  CALL gs_flux_cumint(gseq, psi_n, cflds, self%npsi+1, 5, cum)
END IF
!--- 4. Solve analytically for alpha.
!   gs_flux_int is linear in alpha; two evaluations (alpha=0 and alpha=1) give
!   alpha = (Ip_target - Ip_lo) / (Ip_hi - Ip_lo).  Re-solved every update (cheap): a
!   frozen alpha lets I_p drift while the shape is still settling.
ip_target = ABS(gseq%Ip_target)
wt(3)=omp_get_wtime()
IF(do_saw)THEN
  ip_result_lo = cum(1,1) + cum(1,3)
  ip_result_hi = ip_result_lo + cum(1,2)
ELSE
  jphi_total = j_BS + jphi_fixed + j_saw_cur
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_lo)
  jphi_total = jphi_ind + j_BS + jphi_fixed + j_saw_cur
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_hi)
END IF
ip_ind = ip_result_hi - ip_result_lo
wt(3)=omp_get_wtime()-wt(3)
IF(ABS(ip_ind) > 0.0_r8)THEN
  alpha = (ip_target - ip_result_lo) / ip_ind
  IF(alpha < -1.0_r8 .OR. alpha > 10.0_r8) THEN
    WRITE(char_buf,'(A,ES12.4)') '[jphi_bs_update] WARNING: alpha out of expected range [-1,10]: alpha=', alpha
    CALL oft_warn(TRIM(char_buf))
  END IF
ELSE
  alpha = self%alpha_last
END IF
!--- Relative change: alpha's scale is set by the units of ffp_prof (e.g. ~1e-6 for a jphi in A/m^2)
dalpha = ABS(alpha - self%alpha_last) / MAX(ABS(alpha), TINY(1.0_r8))
!--- Stall diagnostic only (alpha is never frozen): 2 consecutive non-decreasing steps.
IF(dalpha >= self%dalpha_warn_tol .AND. dalpha >= self%dalpha_last) THEN
  self%dalpha_no_improve = self%dalpha_no_improve + 1
  IF(self%dalpha_no_improve >= 2) THEN
    WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
      'Alpha convergence stalled,' // &
      ' relative change per nonlinear step = ', dalpha, &
      ', above recommended tolerance (dalpha_warn_tol=', self%dalpha_warn_tol, ')'
    self%dalpha_no_improve = 0
    IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
  END IF
ELSE
  self%dalpha_no_improve = 0
END IF
self%dalpha_last = dalpha
self%alpha_last = alpha
IF(ip_result_lo > ip_target)THEN
  WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') 'Fixed + bootstrap current (', ip_result_lo/mu0, &
    ' A) exceeds target plasma current (', ip_target/mu0, ' A); inductive current is reversed'
  IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
END IF
!--- 4a. Sawtooth reset of the base current (alpha now known); j_saw freeze check.
IF(do_saw)THEN
  ALLOCATE(dj_saw(0:self%npsi))
  CALL saw_redistribute(self, gseq, psi_n, ravgs, cum(:,1) + alpha*cum(:,2), &
                        cum(:,4), cum(:,5), alpha*jphi_ind + j_BS + jphi_fixed + j_saw_in, dj_saw)
  j_saw = j_saw_in + self%boot_ops%saw_relax*dj_saw + (1.0_r8 - self%boot_ops%saw_relax)*(j_saw_cur - j_saw_in)
  jphi_total = alpha * jphi_ind + j_BS + jphi_fixed + j_saw
  djsaw = SQRT(SUM((j_saw - self%j_saw_last)**2)) / MAX(SQRT(SUM(jphi_total**2)), 1.0e-30_r8)
  ! Freeze only once j_BS has: before that the base current still moves, and an early iterate
  ! without a dip (djsaw = 0) would freeze the reset off for the whole solve
  IF(.NOT.self%freeze_j_BS)THEN
    CONTINUE
  ELSE IF(djsaw < self%boot_ops%saw_tol)THEN
    self%freeze_saw = .TRUE.
    IF(oft_env%pm)WRITE(*,*)' Freezing sawtooth current.'
  ELSE IF(djsaw >= self%djsaw_min)THEN
    self%djsaw_no_improve = self%djsaw_no_improve + 1
    IF(self%djsaw_no_improve >= 2)THEN
      WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
        'Sawtooth current convergence stalled, relative change per nonlinear step = ', djsaw, &
        ', above set tolerance (saw_tol=', self%boot_ops%saw_tol, ')'
      IF(djsaw > self%djBS_stol)THEN
        self%djsaw_no_improve = 0
      ELSE
        self%freeze_saw = .TRUE.
        char_buf = TRIM(char_buf) // ' Freezing sawtooth current.'
      END IF
      IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
    END IF
  ELSE
    self%djsaw_no_improve = 0
    self%djsaw_min = djsaw
  END IF
  DEALLOCATE(dj_saw, cflds, cum)
ELSE
  j_saw = j_saw_cur
END IF
IF(saw_on)THEN
  IF(.NOT.ASSOCIATED(self%j_saw_last))ALLOCATE(self%j_saw_last(0:self%npsi))
  self%j_saw_last = j_saw
END IF
!--- 5. Assemble jphi_total, save profiles
jphi_total = alpha * jphi_ind + j_BS + jphi_fixed + j_saw
IF(saw_on)THEN
  IF(.NOT.ASSOCIATED(self%jtot_last))ALLOCATE(self%jtot_last(0:self%npsi))
  self%jtot_last = jphi_total
END IF
IF(.NOT.ASSOCIATED(self%boot_profs%total_j_phi))THEN
  ALLOCATE(self%boot_profs%psi_n(0:self%npsi))
  ALLOCATE(self%boot_profs%total_j_phi(0:self%npsi))
  ALLOCATE(self%boot_profs%j_bs_final(0:self%npsi))
  ALLOCATE(self%boot_profs%j_ind_final(0:self%npsi))
END IF
IF(.NOT.ASSOCIATED(self%boot_profs%jphi_fixed))ALLOCATE(self%boot_profs%jphi_fixed(0:self%npsi))
IF(.NOT.ASSOCIATED(self%boot_profs%j_saw))ALLOCATE(self%boot_profs%j_saw(0:self%npsi))
self%boot_profs%j_saw       = j_saw/mu0
self%boot_profs%psi_n       = psi_n
self%boot_profs%total_j_phi = jphi_total/mu0
self%boot_profs%j_bs_final  = j_BS/mu0
self%boot_profs%j_ind_final = alpha * jphi_ind/mu0
self%boot_profs%jphi_fixed  = jphi_fixed/mu0
!--- Compute updated F*F' profile
pprime = gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(jphi_total(0) - ravgs(0,1)*pprime*pscale)/ravgs(0,2)
DO i = 1, self%npsi
  pprime = gseq%P%fp(xpsi(i)*(gseq%plasma_bounds(2) - &
                                  gseq%plasma_bounds(1)) + &
                                  gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(jphi_total(i) - ravgs(i,1)*pprime*pscale)/ravgs(i,2)
END DO
! Fix F*F' scale (matching is done here instead)
! gseq%skip_targets is already true when jphi_bs_update called
gseq%ffp_scale=1.d0
gseq%p_scale=pscale
wt(4)=omp_get_wtime()-wt(5)
!--- 6. Diagnostics.
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A,L1,4(A,F9.4))') '  [bs_timing] bs_frozen=', bs_frozen, ' ravg=', wt(1), &
    ' bootstrap=', wt(2), ' alpha=', wt(3), ' update=', wt(4)
  IF(.NOT.bs_frozen)WRITE(*,'(A,F9.4)') '  [bs_timing] sauter_fc time=', sauter_wtime
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_target   = ', ip_target
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_lo= ', ip_result_lo
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_hi= ', ip_result_hi
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] alpha       = ', alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] dalpha      = ', dalpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djBS        = ', djBS
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] bs_frozen   = ', self%freeze_j_BS
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] j_BS max    = ', MAXVAL(ABS(j_BS))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi max    = ', MAXVAL(ABS(self%jphi))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_fixed max = ', MAXVAL(ABS(jphi_fixed))
  IF(saw_on)THEN
    WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] j_saw max   = ', MAXVAL(ABS(j_saw))/mu0
    WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djsaw       = ', djsaw
    WRITE(*,'(A,L1)')     '  [jphi_bs_update] saw_frozen  = ', self%freeze_saw
  END IF
  !--- Side-by-side Ip comparison: FEM nonlinear solve vs profile flux integral
  CALL gs_itor_nl(gseq, itor_nl)
  CALL eval_jtor_imas(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, itor_flint)
  WRITE(*,'(A)') '  [jphi_bs_update] --- Ip comparison (current jphi_total) ---'
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(gs_itor_nl)    = ', itor_nl/mu0
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(flux_int)      = ', itor_flint/mu0
END IF
!--- Clean up
DEALLOCATE(j_BS, jphi_total, jphi_ind, jphi_fixed, jtor, psi_n, ravgs)
DEALLOCATE(j_saw_in, j_saw_cur, j_saw)
i=self%set_cofs(self%yp)
END SUBROUTINE jphi_bs_update
!---------------------------------------------------------------------------------
!> Sawtooth reset current on the nodes `xn` (OFT convention, mu0*A/m², TokaMaker jphi).
!!
!! q_base = q_eq*I_eq/I_base is the q that the base current implies in the equilibrium geometry
!! (q*I is geometric), with q_eq = |F| q/F on the nodes (q/F by cut-cell, @ref gs_qgeom_cutcell).
!! rho = sqrt(Phi_N), Phi from int q_eq dpsi. @ref saw_reset_1d gives the
!! reset on the nodes whose cumulative-area bin holds at least a quarter of the mean bin area (a
!! node-dense axis, e.g. psi_N ~ rho^2, leaves inner bins without quadrature points), interpolated
!! in rho; inside the reset the base jphi between kept nodes is replaced by its interpolant too (weighted
!! by the reset weight wr), so its fine structure is reset with the rest. dI/dA is mapped back to jphi with the linear part of
!! eval_jtor_imas.
!---------------------------------------------------------------------------------
SUBROUTINE saw_redistribute(self, gseq, xn, ravgs, I_base, I_eq, area, j_base, dj)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
REAL(r8), INTENT(in) :: xn(:) !< Nodes, 0 at LCFS, 1 at axis [n]
REAL(r8), INTENT(in) :: ravgs(:,:) !< <R>, <1/R>, <1/R^2> on the nodes [n,3]
REAL(r8), INTENT(in) :: I_base(:) !< Enclosed base current at each node (gs_flux_cumint) [n]
REAL(r8), INTENT(in) :: I_eq(:) !< Enclosed current of the equilibrium [n]
REAL(r8), INTENT(in) :: area(:) !< Enclosed area [n]
REAL(r8), INTENT(in) :: j_base(:) !< Base jphi on the nodes [n]
REAL(r8), INTENT(out) :: dj(:) !< Reset current [n]
INTEGER(i4) :: n, i, k, n_dips, nk
INTEGER(i4), ALLOCATABLE :: kk(:)
REAL(r8), PARAMETER :: psi_ax = 1.0e-4_r8
REAL(r8) :: rho_s, rho_m, rho_out, w, wf, dAmin, alast, psi
REAL(r8), ALLOCATABLE :: q_eq(:), c1(:), rho(:), qb(:), Ib(:), Ar(:), qn(:), djr(:), wr(:), r(:), phi(:), g(:)
n = SIZE(xn)
ALLOCATE(q_eq(n), c1(n), rho(n), qb(n), Ib(n), Ar(n), qn(n), djr(n), wr(n), r(n), phi(n), g(n))
!---q on the nodes, axis-first (index 1 = innermost node). Nodes without a closed surface (g = 0:
!   the LCFS, the axis, a failed surface) take the q of the next node inward (the axis: outward).
CALL gs_qgeom_cutcell(gseq, n, xn, g)
DO i = 1, n
  k = n + 1 - i
  q_eq(i) = 0.0_r8
  IF(g(k) /= 0.0_r8)THEN
    psi = gseq%plasma_bounds(1) + xn(k)*(gseq%plasma_bounds(2) - gseq%plasma_bounds(1))
    q_eq(i) = ABS(g(k))*SQRT(MAX(gseq%ffp_scale*gseq%I%f(psi) + gseq%I%f_offset**2, 0.0_r8))
  END IF
  c1(i) = ravgs(k,3)/ravgs(k,2)**2
  Ib(i) = I_base(k)
  Ar(i) = area(k)
END DO
DO i = 2, n
  IF(q_eq(i) <= 0.0_r8)q_eq(i) = q_eq(i-1)
END DO
DO i = n-1, 1, -1
  IF(q_eq(i) <= 0.0_r8)q_eq(i) = q_eq(i+1)
END DO
!---Nodes within psi_ax of the axis (cut-cell q/F is unreliable on surfaces much smaller than a
!   cell) take the quadratic extrapolation in psi of the first three nodes outside it
nk = COUNT(1.0_r8 - xn > psi_ax)
IF(nk >= 3 .AND. nk < n)THEN
  g(1:3) = xn(nk:nk-2:-1)
  phi(1:3) = q_eq(n-nk+1:n-nk+3)
  DO i = 1, n - nk
    psi = xn(n+1-i)
    q_eq(i) = phi(1)*(psi-g(2))*(psi-g(3))/((g(1)-g(2))*(g(1)-g(3))) &
            + phi(2)*(psi-g(1))*(psi-g(3))/((g(2)-g(1))*(g(2)-g(3))) &
            + phi(3)*(psi-g(1))*(psi-g(2))/((g(3)-g(1))*(g(3)-g(2)))
  END DO
END IF
!---rho_tor_norm from Phi = int q dpsi_N, from the innermost node (its inner piece as a rectangle)
phi(1) = q_eq(1)*(1.0_r8 - xn(n))
DO i = 2, n
  phi(i) = phi(i-1) + 0.5_r8*(q_eq(i) + q_eq(i-1))*(xn(n+2-i) - xn(n+1-i))
END DO
rho = SQRT(phi/phi(n))
!---Nodes kept for the reset: the axis (if it is a node) and every node whose bin adds >= dAmin
ALLOCATE(kk(n))
dAmin = 0.25_r8*Ar(n)/REAL(n-1,r8)
nk = 0
alast = 0.0_r8
IF(1.0_r8 - xn(n) < 1.0e-12_r8)THEN
  nk = 1; kk(1) = 1
END IF
DO i = 2, n
  IF(Ar(i) - alast >= dAmin .OR. (i == n .AND. nk < 2))THEN
    nk = nk + 1; kk(nk) = i; alast = Ar(i)
  ELSE IF(i == n)THEN
    kk(nk) = n
  END IF
END DO
!---q_base on the kept nodes; the axis takes the ratio of the first node off it
DO i = 1, nk
  r(i) = 1.0_r8
  IF(ABS(Ib(kk(i))) > 0.0_r8 .AND. Ar(kk(i)) > 0.0_r8)r(i) = I_eq(n+1-kk(i))/Ib(kk(i))
END DO
IF(Ar(kk(1)) <= 0.0_r8 .AND. nk > 1)r(1) = r(2)
DO i = 1, nk
  qb(i) = q_eq(kk(i))*r(i)
END DO
CALL saw_reset_1d(nk, rho(kk(1:nk)), qb(1:nk), Ib(kk(1:nk)), Ar(kk(1:nk)), c1(kk(1:nk)), &
                  self%boot_ops%saw_q_s, self%boot_ops%saw_dq, self%boot_ops%saw_ramp, &
                  self%boot_ops%saw_rule, qn(1:nk), djr(1:nk), wr(1:nk), rho_s, rho_m, rho_out, n_dips, w)
!---Reset interpolated in rho; inside the reset the base structure between kept nodes is reset too
!   (target total = base + reset is interpolated), weighted by the reset weight (rule 1: w inside
!   rho_m, a sharp cut as in FUSE; rule 2: wr interpolated) so it vanishes with the reset
DO i = 1, nk
  phi(i) = j_base(n+1-kk(i))
END DO
DO i = 1, n
  dj(n + 1 - i) = linterp(rho(kk(1:nk)), djr(1:nk), nk, rho(i), 1)
  IF(self%boot_ops%saw_rule == 1)THEN
    wf = MERGE(w, 0.0_r8, rho(i) < rho_m)
  ELSE
    wf = linterp(rho(kk(1:nk)), wr(1:nk), nk, rho(i), 1)
  END IF
  IF(wf > 0.0_r8)dj(n + 1 - i) = dj(n + 1 - i) + wf*(linterp(rho(kk(1:nk)), phi(1:nk), nk, rho(i), 1) - j_base(n+1-i))
END DO
IF(.NOT.ALL(ABS(dj) < HUGE(1.0_r8)))THEN
  CALL oft_warn('saw_redistribute: non-finite reset current, skipping the reset this update')
  dj = 0.0_r8
END IF
self%boot_profs%saw_rho_m = rho_m
self%boot_profs%saw_rho_out = rho_out
self%boot_profs%saw_n_dips = n_dips
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A,3F8.4,I3,F7.3,2F8.4)') '  [saw] rho_s rho_m rho_out n_dips w q_base(0) q_new(0) = ', &
    rho_s, rho_m, rho_out, n_dips, w, qb(1), qn(1)
END IF
DEALLOCATE(q_eq, c1, rho, qb, Ib, Ar, qn, djr, wr, r, phi, g, kk)
END SUBROUTINE saw_redistribute
!---------------------------------------------------------------------------------
!> Sawtooth q reset on a 1-D profile (axis first).
!!
!! rule 1 (fuse, after FUSE's saw_crash!): from the axis to rho_m, the first crossing of q = q_s + dq
!! beyond the outermost q < q_s point (interpolated, so it moves continuously), q_t = q_s +
!! dq*p(rho/rho_m), p C1 onto q at rho_m (cubic, or x^m for slope m > 3), weighted by
!! w = clamp((q_s - min q)/ramp, 0, 1) (1 for ramp <= 0).
!! rule 2 (local): each maximal q < q_s + dq interval [rho_in, rho_out] (interpolated crossings,
!! rho_in = 0 at the axis) that holds a q < q_s point and has an outer crossing is reset about its
!! q minimum rho_c (parabola-refined, 0 at the axis node): q_s + dq*p((rho-rho_c)/(rho_out-rho_c))
!! outside rho_c, mirrored inside as q_s + D*p((rho_c-rho)/(rho_c-rho_in)), D = q(rho_in) - q_s
!! (max(q(0) - q_s, 0) at the axis; flat if 0), weighted per interval by w_k. Humps inside an interval
!! (between runs, or a run and an end) blend it toward per-run resets split at the hump tops by the hump
!! height (local_reset), so an interval splitting in two moves q_new continuously.
!! q is unchanged elsewhere, so dI vanishes at both ends of every interval.
!! dI = I*(q/q_new - 1), dj = (d dI/dA)/c1. Both rules: rho_s / rho_m = q_s / q_s + dq crossings
!! beyond the outermost q < q_s point; rho_out = outer end of the (outermost) reset; wr = per-node
!! weight (w_k inside interval k, else 0); w = max weight; n_dips = regions reset (rule 1: q < q_s regions).
!---------------------------------------------------------------------------------
SUBROUTINE saw_reset_1d(n, rho, q, I, A, c1, q_s, dq, ramp, rule, q_new, dj, wr, rho_s, rho_m, rho_out, n_dips, w)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: rho(n) !< Radius, ascending from the axis
REAL(r8), INTENT(in) :: q(n) !< Base q
REAL(r8), INTENT(in) :: I(n) !< Enclosed current
REAL(r8), INTENT(in) :: A(n) !< Enclosed area
REAL(r8), INTENT(in) :: c1(n) !< d(jtor)/d(jphi)
REAL(r8), INTENT(in) :: q_s, dq, ramp
INTEGER(i4), INTENT(in) :: rule !< 1 = fuse, 2 = local
REAL(r8), INTENT(out) :: q_new(n), dj(n), wr(n), rho_s, rho_m, rho_out, w
INTEGER(i4), INTENT(out) :: n_dips
INTEGER(i4) :: k, ia, n_low
INTEGER(i4), ALLOCATABLE :: reg_start(:), reg_end(:)
REAL(r8) :: qm
LOGICAL :: ok
REAL(r8), ALLOCATABLE :: dI(:), gq(:), qt(:)
q_new = q; dj = 0.0_r8; wr = 0.0_r8; rho_s = 0.0_r8; rho_m = 0.0_r8; rho_out = 0.0_r8; w = 0.0_r8; n_dips = 0
IF(n < 3)RETURN
!---Separate q < q_s regions [reg_start(k), reg_end(k)]
ALLOCATE(reg_start(n), reg_end(n), gq(n), dI(n), qt(n))
n_low = 0
DO k = 1, n
  IF(q(k) >= q_s)CYCLE
  IF(k == 1)THEN
    n_low = n_low + 1
    reg_start(n_low) = k
  ELSE IF(q(k-1) >= q_s)THEN
    n_low = n_low + 1
    reg_start(n_low) = k
  END IF
  reg_end(n_low) = k
END DO
IF(n_low == 0)RETURN
CALL grad_nonuniform(n, rho, q, gq)
qm = q_s + dq
!---Sawtooth radius (both rules) and rule 1's reset, from the outermost q < q_s region
CALL full_reset(n_low, qt, rho_s, rho_m, ok)
IF(.NOT.ok)THEN
  rho_s = 0.0_r8; rho_m = 0.0_r8
END IF
IF(rule == 1)THEN
  n_dips = n_low
  IF(ok)THEN
    !---Trigger weight from the deepest point inside rho_m
    w = 1.0_r8
    IF(ramp > 0.0_r8)w = MIN(MAX((q_s - MINVAL(q, MASK=(rho < rho_m)))/ramp, 0.0_r8), 1.0_r8)
    q_new = q + w*(qt - q)
    rho_out = rho_m
    WHERE(rho < rho_m)wr = w
  END IF
ELSE
  !---Maximal q < q_s + dq intervals [ia, k] with an outer crossing and a q < q_s point
  ia = 0
  DO k = 1, n
    IF(q(k) >= qm)CYCLE
    IF(ia == 0)ia = k
    IF(k == n)EXIT
    IF(q(k+1) < qm)CYCLE
    IF(MINVAL(q(ia:k)) < q_s)CALL local_reset(ia, k)
    ia = 0
  END DO
END IF
dI = I*(q/q_new - 1.0_r8)
CALL grad_nonuniform(n, A, dI, dj)
dj = dj/c1
DEALLOCATE(reg_start, reg_end, gq, dI, qt)
CONTAINS
!---Full (w = 1) reset from the end of dip region kr; ok = .FALSE. if no mixing radius
SUBROUTINE full_reset(kr, qt_, rs_, rm_, ok_)
INTEGER(i4), INTENT(in) :: kr
REAL(r8), INTENT(out) :: qt_(n), rs_, rm_
LOGICAL, INTENT(out) :: ok_
INTEGER(i4) :: is, jm, kk
REAL(r8) :: qm, slope, m, x, p
qt_ = q; rs_ = 0.0_r8; rm_ = 0.0_r8; ok_ = .FALSE.
is = reg_end(kr)
IF(is >= n)RETURN
rs_ = rho(is) + (q_s - q(is))/(q(is+1) - q(is))*(rho(is+1) - rho(is))
qm = q_s + dq
jm = 0
DO kk = is + 1, n
  IF(q(kk) >= qm)THEN
    jm = kk
    EXIT
  END IF
END DO
IF(jm == 0)RETURN
rm_ = rho(jm-1) + (qm - q(jm-1))/(q(jm) - q(jm-1))*(rho(jm) - rho(jm-1))
!---Slope at rho_m from interpolated nodal gradients (continuous in rho_m)
slope = gq(jm-1) + (rm_ - rho(jm-1))/(rho(jm) - rho(jm-1))*(gq(jm) - gq(jm-1))
m = MAX(rm_*slope/dq, 0.0_r8)
DO kk = 1, jm - 1
  x = rho(kk)/rm_
  IF(m <= 3.0_r8)THEN
    p = (3.0_r8 - 2.0_r8*x)*x**2 + m*(x**3 - x**2)
  ELSE
    p = x**m
  END IF
  qt_(kk) = q_s + dq*p
END DO
ok_ = .TRUE.
END SUBROUTINE full_reset
!---Local reset of interval [ia_, ib_] (q < q_s + dq, ib_ < n). Humps (q maxima) split it into segments:
!   between two q < q_s runs the highest point of the gap, between a run and an interval end its highest
!   interior local maximum. The merged reset (one centre, interval ends) blends into the split one (each
!   run reset between its bounding hump tops, shoulder height h - q_s, zero slope there; no change in a
!   run-free segment) by s = clamp((h - b)/(q_s + dq - b), 0, 1), b = max(q_s, the higher neighbouring
!   minimum), linear in rho between humps: a hump reaching q_s + dq meets the split intervals continuously.
SUBROUTINE local_reset(ia_, ib_)
INTEGER(i4), INTENT(in) :: ia_, ib_
INTEGER(i4) :: kk, j, nr, ih, m, nsg
INTEGER(i4) :: rs(ib_-ia_+1), re(ib_-ia_+1), srun(ib_-ia_+3)
REAL(r8) :: r_in, r_out, r_c, d_in, s_in, s_out, wk, qt_k, sk, dl, dr, sl, sr, qsp, wj
REAL(r8) :: e(0:ib_-ia_+3), hh(ib_-ia_+3), sh(ib_-ia_+3)
!---Ends: crossings of q_s + dq and the slopes there (interpolated nodal gradients)
r_out = rho(ib_) + (qm - q(ib_))/(q(ib_+1) - q(ib_))*(rho(ib_+1) - rho(ib_))
s_out = gq(ib_) + (r_out - rho(ib_))/(rho(ib_+1) - rho(ib_))*(gq(ib_+1) - gq(ib_))
IF(ia_ == 1)THEN
  r_in = 0.0_r8; d_in = MAX(q(1) - q_s, 0.0_r8); s_in = 0.0_r8
ELSE
  r_in = rho(ia_-1) + (qm - q(ia_-1))/(q(ia_) - q(ia_-1))*(rho(ia_) - rho(ia_-1))
  s_in = gq(ia_-1) + (r_in - rho(ia_-1))/(rho(ia_) - rho(ia_-1))*(gq(ia_) - gq(ia_-1))
  d_in = dq
END IF
r_c = centre(ia_, ib_, r_in, r_out)
wk = run_weight(ia_, ib_)
!---q < q_s runs [rs(j), re(j)]
nr = 0
DO kk = ia_, ib_
  IF(q(kk) >= q_s)CYCLE
  IF(kk == ia_)THEN
    nr = nr + 1; rs(nr) = kk
  ELSE IF(q(kk-1) >= q_s)THEN
    nr = nr + 1; rs(nr) = kk
  END IF
  re(nr) = kk
END DO
!---Segments (srun = run, 0 = run-free) and the hump tops e(m) / heights hh(m) / weights sh(m) between them
nsg = 0
ih = edge_hump(MAX(ia_, 2), rs(1) - 1)
IF(ih > 0)THEN
  nsg = nsg + 1; srun(nsg) = 0
  CALL vertex_max(ih, e(nsg), hh(nsg))
  sh(nsg) = hump_weight(hh(nsg), MAX(q_s, MINVAL(q(ia_:ih))))
END IF
DO j = 1, nr
  nsg = nsg + 1; srun(nsg) = j
  IF(j < nr)THEN
    ih = re(j) + MAXLOC(q(re(j)+1:rs(j+1)-1), DIM=1)
    CALL vertex_max(ih, e(nsg), hh(nsg))
    sh(nsg) = hump_weight(hh(nsg), q_s)
  END IF
END DO
ih = edge_hump(re(nr) + 1, ib_)
IF(ih > 0)THEN
  CALL vertex_max(ih, e(nsg), hh(nsg))
  sh(nsg) = hump_weight(hh(nsg), MAX(q_s, MINVAL(q(ih:ib_))))
  nsg = nsg + 1; srun(nsg) = 0
END IF
e(0) = r_in; e(nsg) = r_out
IF(nsg == 1)THEN
  DO kk = ia_, ib_
    qt_k = two_sided(rho(kk), r_in, d_in, s_in, r_c, r_out, dq, s_out)
    q_new(kk) = q(kk) + wk*(qt_k - q(kk))
    wr(kk) = wk
  END DO
ELSE
  DO kk = ia_, ib_
    m = 1 + COUNT(e(1:nsg-1) <= rho(kk))
    !---Split weight: linear in rho between hump tops, constant outside them
    IF(rho(kk) <= e(1))THEN
      sk = sh(1)
    ELSE IF(rho(kk) >= e(nsg-1))THEN
      sk = sh(nsg-1)
    ELSE
      sk = sh(m-1) + (rho(kk) - e(m-1))/(e(m) - e(m-1))*(sh(m) - sh(m-1))
    END IF
    qt_k = two_sided(rho(kk), r_in, d_in, s_in, r_c, r_out, dq, s_out)
    qsp = q(kk); wj = 0.0_r8
    IF(srun(m) > 0)THEN
      j = srun(m)
      dl = d_in; sl = s_in; dr = dq; sr = s_out
      IF(m > 1)THEN
        dl = MAX(hh(m-1) - q_s, 0.0_r8); sl = 0.0_r8
      END IF
      IF(m < nsg)THEN
        dr = MAX(hh(m) - q_s, 0.0_r8); sr = 0.0_r8
      END IF
      qsp = two_sided(rho(kk), e(m-1), dl, sl, centre(rs(j), re(j), e(m-1), e(m)), e(m), dr, sr)
      wj = run_weight(rs(j), re(j))
    END IF
    q_new(kk) = q(kk) + (1.0_r8 - sk)*wk*(qt_k - q(kk)) + sk*wj*(qsp - q(kk))
    wr(kk) = (1.0_r8 - sk)*wk + sk*wj
  END DO
END IF
n_dips = n_dips + 1
w = MAX(w, wk)
rho_out = r_out
END SUBROUTINE local_reset
!---Highest interior local maximum of q in [a_, b_] (0 if none)
INTEGER(i4) FUNCTION edge_hump(a_, b_)
INTEGER(i4), INTENT(in) :: a_, b_
INTEGER(i4) :: kk
edge_hump = 0
DO kk = a_, b_
  IF(q(kk) < q(kk-1) .OR. q(kk) < q(kk+1))CYCLE
  IF(edge_hump == 0)THEN
    edge_hump = kk
  ELSE IF(q(kk) > q(edge_hump))THEN
    edge_hump = kk
  END IF
END DO
END FUNCTION edge_hump
!---Split weight of a hump of height h over the higher neighbouring minimum b
REAL(r8) FUNCTION hump_weight(h, b)
REAL(r8), INTENT(in) :: h, b
hump_weight = 0.0_r8
IF(qm > b)hump_weight = MIN(MAX((h - b)/(qm - b), 0.0_r8), 1.0_r8)
END FUNCTION hump_weight
!---Reset centre in [a_, b_]: q minimum, refined by the parabola through it and its neighbours
!   (0 at the axis node), clamped to [lo, hi]
REAL(r8) FUNCTION centre(a_, b_, lo, hi)
INTEGER(i4), INTENT(in) :: a_, b_
REAL(r8), INTENT(in) :: lo, hi
INTEGER(i4) :: ic
REAL(r8) :: d0, d2, s0, s2, cc
ic = a_ - 1 + MINLOC(q(a_:b_), DIM=1)
centre = 0.0_r8
IF(ic > 1)THEN
  d0 = rho(ic-1) - rho(ic); d2 = rho(ic+1) - rho(ic)
  s0 = (q(ic-1) - q(ic))/d0; s2 = (q(ic+1) - q(ic))/d2
  cc = (s2 - s0)/(d2 - d0)
  centre = rho(ic)
  IF(cc > 0.0_r8)centre = rho(ic) - 0.5_r8*(s0 - cc*d0)/cc
  centre = MIN(MAX(centre, lo), hi)
END IF
END FUNCTION centre
!---Hump top: maximum of the parabola through node i_ and its neighbours (at least q(i_))
SUBROUTINE vertex_max(i_, rv, qv)
INTEGER(i4), INTENT(in) :: i_
REAL(r8), INTENT(out) :: rv, qv
REAL(r8) :: d0, d2, s0, s2, cc, bb, t
rv = rho(i_); qv = q(i_)
d0 = rho(i_-1) - rho(i_); d2 = rho(i_+1) - rho(i_)
s0 = (q(i_-1) - q(i_))/d0; s2 = (q(i_+1) - q(i_))/d2
cc = (s2 - s0)/(d2 - d0)
IF(cc >= 0.0_r8)RETURN
bb = s0 - cc*d0
t = MIN(MAX(-0.5_r8*bb/cc, d0), d2)
rv = rho(i_) + t
qv = MAX(q(i_) + bb*t + cc*t**2, q(i_))
END SUBROUTINE vertex_max
!---Trigger weight of the q < q_s points in [a_, b_]
REAL(r8) FUNCTION run_weight(a_, b_)
INTEGER(i4), INTENT(in) :: a_, b_
run_weight = 1.0_r8
IF(ramp > 0.0_r8)run_weight = MIN(MAX((q_s - MINVAL(q(a_:b_)))/ramp, 0.0_r8), 1.0_r8)
END FUNCTION run_weight
!---Two-sided target at x: q_s at r_c_, shoulders of height d_in_ / d_out_ C1 onto slopes s_in_ /
!   s_out_ at r_in_ / r_out_ (flat at q_s where the height is 0)
REAL(r8) FUNCTION two_sided(x, r_in_, d_in_, s_in_, r_c_, r_out_, d_out_, s_out_)
REAL(r8), INTENT(in) :: x, r_in_, d_in_, s_in_, r_c_, r_out_, d_out_, s_out_
two_sided = q_s
IF(x >= r_c_)THEN
  IF(d_out_ > 0.0_r8 .AND. r_out_ > r_c_)two_sided = q_s + d_out_*pshape((x - r_c_)/(r_out_ - r_c_), &
    MAX((r_out_ - r_c_)*s_out_/d_out_, 0.0_r8))
ELSE IF(d_in_ > 0.0_r8)THEN
  two_sided = q_s + d_in_*pshape((r_c_ - x)/(r_c_ - r_in_), MAX(-(r_c_ - r_in_)*s_in_/d_in_, 0.0_r8))
END IF
END FUNCTION two_sided
!---FUSE's shoulder: cubic p(0) = 0, p(1) = 1, p'(0) = 0, p'(1) = m (m <= 3), else x^m
REAL(r8) FUNCTION pshape(x, m)
REAL(r8), INTENT(in) :: x, m
IF(m <= 3.0_r8)THEN
  pshape = (3.0_r8 - 2.0_r8*x)*x**2 + m*(x**3 - x**2)
ELSE
  pshape = x**m
END IF
END FUNCTION pshape
!---Derivative dy/dx on a non-uniform grid (2nd-order interior, one-sided ends)
SUBROUTINE grad_nonuniform(nn, xx, yy, dy)
INTEGER(i4), INTENT(in) :: nn
REAL(r8), INTENT(in) :: xx(nn), yy(nn)
REAL(r8), INTENT(out) :: dy(nn)
INTEGER(i4) :: kk
REAL(r8) :: h1, h2
dy(1) = (yy(2) - yy(1))/(xx(2) - xx(1))
dy(nn) = (yy(nn) - yy(nn-1))/(xx(nn) - xx(nn-1))
DO kk = 2, nn - 1
  h1 = xx(kk) - xx(kk-1); h2 = xx(kk+1) - xx(kk)
  dy(kk) = (h1**2*yy(kk+1) - h2**2*yy(kk-1) + (h2**2 - h1**2)*yy(kk))/(h1*h2*(h1 + h2))
END DO
END SUBROUTINE grad_nonuniform
END SUBROUTINE saw_reset_1d
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
integer(4) :: j(self%lag_rep%nce)
integer(4) :: jc
real(8) :: rop(3),pt(3),grad(3)
real(8) :: s,c,Bp2,Bt2,mod_B
!---Get dofs
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
!---Magnetic field averages (<|B|^2> first: recorded per step for the trapped fraction)
val(3)=val(2)*(Bp2+Bt2)     ! <|B|^2>
val(4)=val(2)*mod_B         ! <|B|>
!---Geometric factors
val(5)=val(2)*pt(1)         ! <R>
val(6)=val(2)/pt(1)         ! <1/R>
val(7)=val(2)*SQRT((pt(1)-self%mag_axis(1))**2+(pt(2)-self%mag_axis(2))**2)  ! <a>
val(8)=val(2)/pt(1)**2      ! q integrand: q = F_surf/(2π) * ∮ val(8)
!---Surface extrema, sampled at ODE evaluations
self%bmax = MAX(self%bmax,mod_B)
self%rmax_surf = MAX(self%rmax_surf, pt(1))
self%rmin_surf = MIN(self%rmin_surf, pt(1))
end subroutine sauter_apply
!------------------------------------------------------------------------------
!> Compute factors required for Sauter bootstrap formula
!!
!! Surface averages and extrema come from the cut-cell quadrature
!! (@ref oft_gs_cutcell::gs_sauter_cutcell, two passes as the trapped-particle sum
!! needs Bmax). Surfaces it cannot close, or all surfaces when
!! `torflux_qgeom_backend=1` (testing), are traced instead.
!!
!! Each traced surface is traced once. The trapped-particle integral
!! \f$ \langle (1-\sqrt{1-b}(1+b/2))/b^2 \rangle \f$, \f$ b=|B|/B_{max} \f$, needs
!! \f$ B_{max} \f$ of the whole surface, so it is evaluated after the trace by
!! trapezoidal quadrature over the tracer steps. Surfaces are traced in parallel
!! from start points found serially, so results do not depend on the thread count.
!------------------------------------------------------------------------------
subroutine sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs,qprof,eps,invr2)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(4), intent(in) :: nr !< Number of flux sample points
real(8), intent(in) :: psi_q(nr) !< Location of flux sample points in normalised psi
real(8), intent(out) :: fc(nr) !< Circulating particle fraction \f$ f_c \f$
real(8), intent(out) :: r_avgs(nr,3) !< Flux surface averaged radial coords \f$<R>\f$, \f$<1/R>\f$, \f$<a>\f$
real(8), intent(out) :: modb_avgs(nr,2) !< Flux surface averaged field \f$<|B|>\f$, \f$<|B|^2>\f$
real(8), optional, intent(out) :: qprof(nr) !< Safety factor q on each surface (avoids a separate gs_get_qprof call)
real(8), optional, intent(out) :: eps(nr) !< Local inverse aspect ratio \f$ \varepsilon = (R_{\max}-R_{\min})/(2\langle R \rangle) \f$ on each surface
real(8), optional, intent(out) :: invr2(nr) !< \f$<1/R^2>\f$ on each surface
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,h,h2,hf,ftu,ftl,t0
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
real(8), allocatable :: pts(:,:),fpol(:),qtmp(:),epstmp(:),i2tmp(:)
real(8), pointer :: ptout(:,:),vout(:,:)
real(8), allocatable :: acc(:,:),ext(:,:),bref(:)
logical, allocatable :: done(:)
type(oft_lag_brinterp) :: psi_int
real(8), parameter :: tol=1.d-10
integer(4) :: j,cell
type(sauter_interp), pointer :: field
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
ALLOCATE(pts(2,nr),fpol(nr),qtmp(nr),epstmp(nr),i2tmp(nr),done(nr))
DO j=1,nr
  psi_surf=psi_q(j)*(x2-x1) + x1
  IF(gseq%mode==0)THEN
    fpol(j)=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    fpol(j)=SIGN(1.d0,gseq%I%f_offset)*SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2)
  END IF
END DO
!--- Axis guard: flux surface degenerates at psi_N=1; prescribe analytically.
done=(psi_q == 1.d0)
DO j=1,nr
  IF(.NOT.done(j))CYCLE
  r_avgs(j,:)    = [raxis, 1.d0/raxis, 0.d0]
  modb_avgs(j,:) = [ABS(fpol(j))/raxis, (fpol(j)/raxis)**2]
  fc(j)          = 1.d0
  qtmp(j)        = 0.d0 ! Taken from the previous surface below
  epstmp(j)      = 0.d0
  i2tmp(j)       = 1.d0/raxis**2
END DO
!---Cut-cell averages and extrema
sauter_wtime=0.d0
IF(torflux_qgeom_backend/=1)THEN
  t0=omp_get_wtime()
  ALLOCATE(acc(8,nr),ext(3,nr),bref(nr))
  bref=-1.d0
  CALL gs_sauter_cutcell(gseq,nr,psi_q,fpol,bref,acc,ext,ok=done)
  bref=ext(1,:)
  CALL gs_sauter_cutcell(gseq,nr,psi_q,fpol,bref,acc,ext,ok=done)
  DO j=1,nr
    IF((.NOT.done(j)).OR.(psi_q(j)==1.d0))CYCLE
    r_avgs(j,1)=acc(2,j)/acc(1,j)
    r_avgs(j,2)=acc(3,j)/acc(1,j)
    r_avgs(j,3)=acc(4,j)/acc(1,j)
    modb_avgs(j,1)=acc(5,j)/acc(1,j)
    modb_avgs(j,2)=acc(6,j)/acc(1,j)
    h = modb_avgs(j,1)/ext(1,j)
    h2 = modb_avgs(j,2)/ext(1,j)**2
    ftu = 1.d0 - h2 / (h**2) * (1.d0 - SQRT(1.d0 - h) * (1.d0 + 0.5d0 * h))
    hf = acc(8,j)/acc(1,j)
    ftl = 1.d0 - h2 * hf
    fc(j) = 1.d0 - (0.75d0 * ftu + 0.25d0 * ftl)
    qtmp(j) = fpol(j) * acc(7,j) / (2.d0*pi*ABS(x2-x1))
    epstmp(j) = (ext(2,j) + ext(3,j)) / (2.d0 * r_avgs(j,1))
    i2tmp(j) = acc(7,j)/acc(1,j)
  END DO
  done=done.OR.(psi_q == 1.d0)
  DEALLOCATE(acc,ext,bref)
  sauter_wtime=omp_get_wtime()-t0
  IF(oft_debug_print(1).AND.(.NOT.ALL(done)))WRITE(*,'(2A,I0,A)')oft_indent, &
    'sauter_fc: tracing ',COUNT(.NOT.done),' surfaces not closed by the cut-cell pass'
END IF
!---Start points of the remaining surfaces (each from the previous one)
DO j=1,nr
  IF(done(j))CYCLE
  pt=pt_last
  CALL gs_psi2r(gseq,psi_q(j)*(x2-x1) + x1,pt,psi_int=psi_int)
  pts(:,j)=pt(1:2)
  pt_last=pt
END DO
CALL psi_int%delete()
!---Trace the remaining surfaces
IF(.NOT.ALL(done))THEN
call set_tracer(1)
!$omp parallel private(j,field,ptout,vout,h,h2,hf,ftu,ftl,t0) reduction(+:sauter_wtime)
ALLOCATE(field)
field%u=>gseq%psi
field%mag_axis=gseq%o_point
CALL field%setup(device%fe_rep)
active_tracer%neq=8
active_tracer%B=>field
active_tracer%maxsteps=8e4
IF(qprof_trace_tol>0.d0)active_tracer%maxsteps=2e6
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1),vout(3,active_tracer%maxsteps))
!$omp do schedule(dynamic,1)
do j=1,nr
  IF(done(j))CYCLE
  IF(gseq%diverted.AND.psi_q(j)<0.02d0)THEN
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  IF(qprof_trace_tol>0.d0)active_tracer%tol=qprof_trace_tol*MERGE(1.d-2,1.d0,gseq%diverted.AND.psi_q(j)<0.02d0)
  field%f_surf=fpol(j)
  field%bmax=0.d0
  field%rmax_surf = -1.d30
  field%rmin_surf =  1.d30
  t0=omp_get_wtime()
  CALL tracinginv_fs(device%mesh,pts(:,j),ptout,vout)
  sauter_wtime=sauter_wtime+omp_get_wtime()-t0
  if(active_tracer%status/=1)THEN
    WRITE(*,*)j,pts(:,j)
    CALL oft_warn("sauter_fc: Trace did not complete")
    CYCLE
  end if
  r_avgs(j,1)=active_tracer%v(5)/active_tracer%v(2)
  r_avgs(j,2)=active_tracer%v(6)/active_tracer%v(2)
  r_avgs(j,3)=active_tracer%v(7)/active_tracer%v(2)
  modb_avgs(j,1)=active_tracer%v(4)/active_tracer%v(2)
  modb_avgs(j,2)=active_tracer%v(3)/active_tracer%v(2)
  h = modb_avgs(j,1)/field%bmax
  h2 = modb_avgs(j,2)/field%bmax**2
  ftu = 1.d0 - h2 / (h**2) * (1.d0 - SQRT(1.d0 - h) * (1.d0 + 0.5d0 * h))
  hf = sauter_hf(active_tracer%nsteps,ptout(1,1:active_tracer%nsteps+1),vout,field%bmax)
  ftl = 1.d0 - h2 * hf
  fc(j) = 1.d0 - (0.75d0 * ftu + 0.25d0 * ftl)
  qtmp(j) = fpol(j) * active_tracer%v(8) / (2*pi)
  epstmp(j) = (field%rmax_surf - field%rmin_surf) / (2.d0 * r_avgs(j,1))
  i2tmp(j) = active_tracer%v(8)/active_tracer%v(2)
end do
DEALLOCATE(ptout,vout)
CALL active_tracer%delete
CALL field%delete
DEALLOCATE(field)
!$omp end parallel
END IF
DO j=1,nr
  IF(psi_q(j) == 1.d0)qtmp(j)=MERGE(qtmp(MAX(j-1,1)), 0.d0, j > 1)
END DO
IF(PRESENT(qprof))qprof=qtmp
IF(PRESENT(eps))eps=epstmp
IF(PRESENT(invr2))invr2=i2tmp
DEALLOCATE(pts,fpol,qtmp,epstmp,i2tmp,done)
end subroutine sauter_fc
!------------------------------------------------------------------------------
!> Trapped-particle integral \f$ \langle (1-\sqrt{1-b}(1+b/2))/b^2 \rangle \f$ by
!! periodic trapezoidal quadrature over the steps of one surface trace
!------------------------------------------------------------------------------
function sauter_hf(ns,t,vout,bmax) result(hf)
integer(4), intent(in) :: ns !< Number of steps (the last ends at \f$ 2\pi \f$)
real(8), intent(in) :: t(:) !< Angle before each step and at its end [ns+1]
real(8), intent(in) :: vout(:,:) !< ODE RHS at each step; (2) FSA weight, (3) weight times \f$ B^2 \f$ [3,ns]
real(8), intent(in) :: bmax !< Maximum \f$ |B| \f$ on the surface
real(8) :: hf,g(2),w(2),b,wsum
integer(4) :: k,i,kk(2)
hf=0.d0; wsum=0.d0
DO k=1,ns
  kk=[MERGE(k-1,ns,k>1),k]
  DO i=1,2
    w(i)=vout(2,kk(i))
    b=MIN(1.d0,SQRT(vout(3,kk(i))/w(i))/bmax)
    g(i)=w(i)*(1.d0 - SQRT(1.d0 - b)*(1.d0 + b/2.d0))/b**2
  END DO
  hf=hf+(t(k+1)-t(k))*SUM(g)
  wsum=wsum+(t(k+1)-t(k))*SUM(w)
END DO
hf=hf/wsum
end function sauter_hf
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
!> Compute dy/dx with second-order-accurate finite differences on the native grid.
!>
!> Port of `numpy.gradient(y, x, edge_order=2)`: the interior uses the
!> three-point central stencil for the local (possibly non-uniform) spacing,
!> and both endpoints use the second-order one-sided stencil.  The first-order
!> endpoint difference numpy uses by default is badly inaccurate at the
!> magnetic axis and separatrix, where the error propagates straight into
!> on-axis/edge j_BS.
!>
!> This replaces an earlier shape-preserving PCHIP derivative.  On the uniform
!> grids this code uses, the PCHIP endpoint formula reduces algebraically to
!> the second-order one-sided difference, while its monotonicity clipping
!> flattens the pedestal and makes the interior measurably less accurate.
!>
!> The grid is validated rather than repaired: unlike the PCHIP routine it
!> replaces, this does not sort x or average over duplicated abscissae.
!> Duplicated or unsorted flux labels give an undefined derivative, so they
!> abort here instead of being silently smoothed.
!>
!> @param n     Number of input points (>= 3)
!> @param x     Independent variable (e.g. psi); strictly monotonic, either direction
!> @param y     Dependent variable sampled on x
!> @param dydx  Output: dy/dx at each point in x
!------------------------------------------------------------------------------
SUBROUTINE gradient_(n, x, y, dydx)
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: x(n)
REAL(r8), INTENT(in) :: y(n)
REAL(r8), INTENT(out) :: dydx(n)
!---
INTEGER(i4) :: i
REAL(r8) :: hs, hd !< Backward and forward interval widths about the stencil centre
CHARACTER(len=512) :: char_buf
IF(n < 3) &
  CALL oft_abort('fewer than 3 input points; a second-order stencil is undefined', &
    'gradient_', __FILE__)
DO i = 1, n
  IF(.NOT.IEEE_IS_FINITE(x(i)) .OR. .NOT.IEEE_IS_FINITE(y(i)))THEN
    WRITE(char_buf,'(A,I0,A,ES12.5,A,ES12.5)') 'non-finite input at i=', i, &
      ': x=', x(i), ', y=', y(i)
    CALL oft_abort(TRIM(char_buf), 'gradient_', __FILE__)
  END IF
END DO
!--- Strict monotonicity in either direction.  The stencils below are signed,
!   so a descending x needs no special handling, but a repeated or reversed
!   abscissa leaves dy/dx undefined.
DO i = 2, n
  IF((x(i) - x(i-1))*(x(2) - x(1)) > 0.0_r8)CYCLE
  WRITE(char_buf,'(A,I0,A,ES12.5,A,ES12.5,A)') 'x not strictly monotonic at '// &
    'i=', i, ' (x=', x(i-1), ' then ', x(i), '); duplicated or unsorted '// &
    'abscissae give an undefined derivative'
  CALL oft_abort(TRIM(char_buf), 'gradient_', __FILE__)
END DO
!--- Interior: three-point central stencil for the local spacing
DO i = 2, n-1
  hs = x(i) - x(i-1)
  hd = x(i+1) - x(i)
  dydx(i) = -(hd/(hs*(hs + hd)))*y(i-1) &
            + ((hd - hs)/(hs*hd))*y(i) &
            + (hs/(hd*(hs + hd)))*y(i+1)
END DO
!--- Endpoints: second-order one-sided stencils (numpy's edge_order=2)
hs = x(2) - x(1)
hd = x(3) - x(2)
dydx(1) = -((2.0_r8*hs + hd)/(hs*(hs + hd)))*y(1) &
          + ((hs + hd)/(hs*hd))*y(2) &
          - (hs/(hd*(hs + hd)))*y(3)
hs = x(n-1) - x(n-2)
hd = x(n) - x(n-1)
dydx(n) = (hd/(hs*(hs + hd)))*y(n-2) &
          - ((hd + hs)/(hs*hd))*y(n-1) &
          + ((2.0_r8*hd + hs)/(hd*(hs + hd)))*y(n)
!--- Finite input does not guarantee a finite result: the stencil weights
!   overflow for large y on a very fine grid.
DO i = 1, n
  IF(.NOT.IEEE_IS_FINITE(dydx(i)))THEN
    WRITE(char_buf,'(A,I0,A,ES12.5)') 'non-finite derivative at i=', i, &
      ', x=', x(i)
    CALL oft_abort(TRIM(char_buf), 'gradient_', __FILE__)
  END IF
END DO
END SUBROUTINE gradient_
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
!> @param j_BS Output: bootstrap current density as TokaMaker jphi = <j_phi> [A/m^2] on psi_N grid
!------------------------------------------------------------------------------
SUBROUTINE calculate_bootstrap(self, gseq, n_psi, psi_N, j_BS, &
                               isolate_edge_jBS, parameterize_jBS, scale_jBS, &
                               j_spike, j_spike_masked, ravgs)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)  !< Normalised psi grid (0,1], arbitrary spacing. psi_N(1)>0.0
REAL(r8), INTENT(out) :: j_BS(0:n_psi) ! j_BS(0) the LCFS value (extrapolated)
LOGICAL,  OPTIONAL, INTENT(in)  :: isolate_edge_jBS  !< If .TRUE., isolate edge spike from core
LOGICAL,  OPTIONAL, INTENT(in)  :: parameterize_jBS  !< If .TRUE., use parametrised skew-normal fit
REAL(r8), OPTIONAL, INTENT(in)  :: scale_jBS         !< Scaling factor applied to spike profile (default 1)
REAL(r8), OPTIONAL, INTENT(out) :: j_spike(0:n_psi)        !< Processed spike profile [A/m^2]
REAL(r8), OPTIONAL, INTENT(out) :: j_spike_masked(0:n_psi) !< Raw masked (pre-fit) spike profile [A/m^2]
REAL(r8), OPTIONAL, INTENT(out) :: ravgs(n_psi,3) !< <R>, <1/R>, <1/R^2> on psi_N from the Sauter pass
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
REAL(r8) :: invr2(n_psi)            !< Sauter FSA: <1/R^2>
REAL(r8) :: modb_avgs_saut(n_psi,2) !< Sauter FSA: <|B|>, <|B|^2>
REAL(r8) :: dn_e_dpsi(n_psi), dT_e_dpsi(n_psi)
REAL(r8) :: dn_i_dpsi(n_psi), dT_i_dpsi(n_psi)
REAL(r8) :: ln_le(n_psi), ln_lii(n_psi), Z_lnLam(n_psi)
REAL(r8) :: Zavg(n_psi), Zion(n_psi)
REAL(r8) :: nu_i_star(n_psi), nu_e_star(n_psi)
REAL(r8) :: B_times_Jbs(n_psi)
REAL(r8) :: psi_range, Zdom, pscale, pprime
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
    f(i) = SIGN(1.0_r8, gseq%I%f_offset) * SQRT(gseq%ffp_scale * gseq%I%f(psi_abs(i)) + gseq%I%f_offset**2)
  END IF
END DO
! Get flux-surface geometry: fc, eps, q, and <R> in one tracing pass
CALL sauter_fc(gseq, n_psi, psi_N, fc, r_avgs_saut, modb_avgs_saut, qprof=qvals, eps=eps, invr2=invr2)
IF(PRESENT(ravgs))ravgs = RESHAPE([r_avgs_saut(:,1:2), invr2], [n_psi,3])
R_avg = r_avgs_saut(:,1)
B_avg = modb_avgs_saut(:,1)
! ===================================================================
ft = 1.0_r8 - fc
! Pressures [Pa]
pe = EC * ne * Te
pi_arr = EC * ni * Ti
! Gradients d/dpsi [Wb^-1]
CALL gradient_(n_psi, psi_abs, Te, dT_e_dpsi)
CALL gradient_(n_psi, psi_abs, Ti, dT_i_dpsi)
CALL gradient_(n_psi, psi_abs, ne, dn_e_dpsi)
CALL gradient_(n_psi, psi_abs, ni, dn_i_dpsi)
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
! Convert <j_BS.B> to TokaMaker's jphi = <j_phi> = <R>P' + <1/R>FF'/mu0 (exact, see
! doc_tokamaker_current_conventions.md eq. A7): the field-aligned part F<1/R>/<B^2> * <j_BS.B>,
! plus the pressure-driven (diamagnetic + Pfirsch-Schlueter) part P'(<R> - F^2<1/R>/<B^2>),
! which is assigned to the bootstrap (as IMAS includes_bootstrap=true). P' as in jphi_bs_update.
pscale = gseq%pax_target/gseq%P%f(gseq%plasma_bounds(2))
j_BS(0) = 0.0_r8 ! Placeholder until extrap_jBS_boundaries sets the real LCFS value below
DO i = 1, n_psi
  IF(ABS(f(i)) > 0.0_r8 .AND. modb_avgs_saut(i,2) > 0.0_r8)THEN
    pprime = gseq%P%fp(psi_abs(i))*pscale/mu0
    j_BS(i) = B_times_Jbs(i)*f(i)*r_avgs_saut(i,2)/modb_avgs_saut(i,2) &
      + pprime*(r_avgs_saut(i,1) - f(i)**2*r_avgs_saut(i,2)/modb_avgs_saut(i,2))
  ELSE
    j_BS(i) = 0.0_r8
  END IF
END DO
! Guard NaN (where F -> 0)
WHERE(.NOT.(ABS(j_BS) < 1.0e99_r8)) j_BS = 0.0_r8
! Extrapolate to LCFS/axis where q is undefined. Sets LCFS value j_BS(0).
CALL extrap_jBS_boundaries(n_psi, psi_N, j_BS)
! Save raw bootstrap output, and Redl's <j_BS.B> on the same grid
IF(.NOT.ASSOCIATED(self%boot_profs%j_bs_raw)) ALLOCATE(self%boot_profs%j_bs_raw(0:n_psi))
self%boot_profs%j_bs_raw = j_BS
IF(.NOT.ASSOCIATED(self%boot_profs%jdotb_bs_raw)) ALLOCATE(self%boot_profs%jdotb_bs_raw(0:n_psi))
self%boot_profs%jdotb_bs_raw(0) = 0.0_r8
self%boot_profs%jdotb_bs_raw(1:) = B_times_Jbs
WHERE(.NOT.(ABS(self%boot_profs%jdotb_bs_raw) < 1.0e99_r8)) self%boot_profs%jdotb_bs_raw = 0.0_r8
CALL extrap_jBS_boundaries(n_psi, psi_N, self%boot_profs%jdotb_bs_raw)
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A)') '  [calculate_bootstrap] geometry & collisionality sample (i=1,mid,n):'
  WRITE(*,'(A,3ES12.4)') '    <R>      : ', r_avgs_saut(1,1), r_avgs_saut(n_psi/2,1), r_avgs_saut(n_psi,1)
  WRITE(*,'(A,3ES12.4)') '    <1/R>    : ', r_avgs_saut(1,2), r_avgs_saut(n_psi/2,2), r_avgs_saut(n_psi,2)
  WRITE(*,'(A,3ES12.4)') '    <B>      : ', B_avg(1), B_avg(n_psi/2), B_avg(n_psi)
  WRITE(*,'(A,3ES12.4)') '    <B^2>    : ', modb_avgs_saut(1,2), modb_avgs_saut(n_psi/2,2), modb_avgs_saut(n_psi,2)
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
    ALLOCATE(psi_N_std(0:n_psi), j_BS_std(0:n_psi), mask_std(0:n_psi))
    IF (do_parametrize) ALLOCATE(param_std(0:n_psi))
    psi_N_std(0:n_psi-1) = 1.0_r8 - psi_N(n_psi:1:-1)
    psi_N_std(n_psi) = 1.0_r8
    j_BS_std  = j_BS(n_psi:0:-1)
    IF (do_parametrize) THEN
      CALL analyze_bootstrap_edge_spike((n_psi+1), psi_N_std, j_BS_std, mask_std, &
                                      param_std, diagnose=self%boot_ops%diagnose_bs)
      IF (PRESENT(j_spike))        j_spike        = scl_jBS * param_std(n_psi:0:-1)
      IF (PRESENT(j_spike_masked)) j_spike_masked = scl_jBS * mask_std(n_psi:0:-1)
      DEALLOCATE(param_std)
    ELSE
      CALL analyze_bootstrap_edge_spike((n_psi+1), psi_N_std, j_BS_std, mask_std)
      IF (PRESENT(j_spike))        j_spike        = scl_jBS * mask_std(n_psi:0:-1)
      IF (PRESENT(j_spike_masked)) j_spike_masked = scl_jBS * mask_std(n_psi:0:-1)
    END IF
    DEALLOCATE(psi_N_std, j_BS_std, mask_std)
  END IF
END IF
END SUBROUTINE calculate_bootstrap
!------------------------------------------------------------------------------
!> Linearly extrapolate j_BS to LCFS/axis points where q is undefined, using
!> the second-order gradient at the nearest well-defined point.
!------------------------------------------------------------------------------
SUBROUTINE extrap_jBS_boundaries(n_psi, psi_N, j_BS)
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)
REAL(r8), INTENT(inout) :: j_BS(0:n_psi)
INTEGER(i4) :: n_lo, n_hi
REAL(r8) :: djBS_dpsi(n_psi)
! psi_N(0) = 0.0 at the LCFS is implicit.
n_lo = 1
n_hi = MERGE(n_psi-1, n_psi, psi_N(n_psi) == 1.0_r8)
IF(n_hi-n_lo+1 < 3) CALL oft_abort('extrap_jBS_boundaries: too few '// &
  'points with well-defined q to extrapolate j_BS to the LCFS', &
  'extrap_jBS_boundaries', __FILE__)
CALL gradient_(n_hi-n_lo+1, psi_N(n_lo:n_hi), j_BS(n_lo:n_hi), &
  djBS_dpsi(n_lo:n_hi))
j_BS(0) = j_BS(n_lo) + djBS_dpsi(n_lo) * (0.0_r8 - psi_N(n_lo))
IF(psi_N(n_psi) == 1.0_r8)THEN
  n_lo = 1
  n_hi = n_psi-1
  IF(n_hi-n_lo+1 < 3) CALL oft_abort('extrap_jBS_boundaries: too few '// &
    'points with well-defined q to extrapolate j_BS to the axis', &
    'extrap_jBS_boundaries', __FILE__)
  CALL gradient_(n_hi-n_lo+1, psi_N(n_lo:n_hi), j_BS(n_lo:n_hi), &
    djBS_dpsi(n_lo:n_hi))
  j_BS(n_psi) = j_BS(n_hi) + djBS_dpsi(n_hi) * (psi_N(n_psi) - psi_N(n_hi))
ENDIF
END SUBROUTINE extrap_jBS_boundaries
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
SUBROUTINE curve_fit_edge_jbs(n_fit, psi_fit, j_BS_fit, p0, popt, &
    lb_in, ub_in, success)
INTEGER(i4), INTENT(in)  :: n_fit
REAL(r8),    INTENT(in)  :: psi_fit(n_fit), j_BS_fit(n_fit)
REAL(r8),    INTENT(in)  :: p0(7)
REAL(r8),    INTENT(out) :: popt(7)
REAL(r8),    INTENT(in)  :: lb_in(7), ub_in(7)
LOGICAL, OPTIONAL, INTENT(out) :: success  !< lmdif convergence flag
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
  IF (PRESENT(success)) success = .FALSE.
ELSE
  IF (PRESENT(success)) success = .TRUE.
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
!>   1. Locates the highest local peak (whose prominence exceeds 5%)
!>      in psi_N in [0.85,1.0], else the rightmost peak for psi_N > 0.7.  
!>      With no edge region or no qualifying peak the un-isolated profile 
!>      is returned (with a warning) rather than zeros.
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
INTEGER(i4) :: i, k, peak_idx, far_idx, lmin_idx, left_idx, right_idx
INTEGER(i4) :: idx_start, idx_end
INTEGER(i4) :: n_fit
LOGICAL     :: has_edge, fit_ok
REAL(r8)    :: far_height, half_max, dist_tmp
REAL(r8)    :: j_edge_max, min_prominence, prom, left_base, right_base
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
! =====================================================================
! 1. Find the highest local peak (whose prominence exceeds 5%)
!    in psi_N in [0.85,1.0], else the rightmost peak for psi_N > 0.7.  
!    Fall back to the input j_BS if no edge region or no qualifying peak.
! =====================================================================
has_edge   = .FALSE.
j_edge_max = 0.0_r8
DO i = 1, n
  IF (psi_N(i) >= 0.7_r8) THEN
    IF (has_edge) THEN
      j_edge_max = MAX(j_edge_max, j_BS(i))
    ELSE
      j_edge_max = j_BS(i)
      has_edge   = .TRUE.
    END IF
  END IF
END DO
IF (.NOT. has_edge) THEN
  CALL oft_warn('[analyze_bootstrap_edge_spike] no edge '// &
    'region (psi_N >= 0.7) in profile; using '// &
    'bootstrap profile without edge isolation')
  masked_spike = j_BS
  IF (PRESENT(parameterized_spike)) parameterized_spike = j_BS
  RETURN
END IF
IF (j_edge_max > 0.0_r8) THEN
  min_prominence = 0.05_r8 * j_edge_max
ELSE
  min_prominence = 0.0_r8
END IF
peak_idx   = -1        ! rightmost qualifying peak
far_idx    = -1        ! tallest far-edge peak (psi_N > 0.85)
far_height = -1.0_r8
DO i = 2, n-1
  IF (psi_N(i) < 0.7_r8)    CYCLE
  IF (j_BS(i) <= 0.0_r8)    CYCLE
  IF (j_BS(i) <= j_BS(i-1)) CYCLE
  IF (j_BS(i) <= j_BS(i+1)) CYCLE
  ! Prominence = peak minus the higher of the two bases (nearest minima
  ! before a taller point on each side, within the psi_N >= 0.7 slice).
  left_base = j_BS(i)
  DO k = i-1, 1, -1
    IF (psi_N(k) < 0.7_r8)  EXIT
    IF (j_BS(k) > j_BS(i))  EXIT
    left_base = MIN(left_base, j_BS(k))
  END DO
  right_base = j_BS(i)
  DO k = i+1, n
    IF (j_BS(k) > j_BS(i))  EXIT
    right_base = MIN(right_base, j_BS(k))
  END DO
  prom = j_BS(i) - MAX(left_base, right_base)
  IF (prom < min_prominence) CYCLE
  peak_idx = i
  IF (psi_N(i) > 0.85_r8 .AND. j_BS(i) > far_height) THEN
    far_idx    = i
    far_height = j_BS(i)
  END IF
END DO
! Prefer the tallest far-edge peak, else the rightmost.
IF (far_idx >= 1) peak_idx = far_idx
IF (peak_idx < 1) THEN
  CALL oft_warn('[analyze_bootstrap_edge_spike] no clear '// &
    'edge peak found; using bootstrap profile '// &
    'without edge isolation')
  masked_spike = j_BS
  IF (PRESENT(parameterized_spike)) parameterized_spike = j_BS
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
  CALL curve_fit_edge_jbs(n_fit, psi_fit, j_fit, p0, popt, &
    lb, ub, fit_ok)
  DEALLOCATE(psi_fit, j_fit)
  IF (.NOT. fit_ok) THEN
    ! Non-fatal fit failure: fall back to masked_spike.
    CALL oft_warn('[analyze_bootstrap_edge_spike] '// &
      'parametrise_edge_jbs fit failed; '// &
      'falling back to masked_spike')
    parameterized_spike = masked_spike
  ELSE
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
    ! Optional verbose output for diagnostics if boot_ops%diagnose_bs
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