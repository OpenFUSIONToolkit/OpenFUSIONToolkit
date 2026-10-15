!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_bootstrap.F90
!
!> Bootstrap current for the jphi-split-bootstrap Grad-Shafranov profile
!!
!! Fortran port of Daniel Burgess' Python bootstrap solve (TokaMaker/bootstrap.py).
!!
!! @authors Stuart Benjamin, Daniel Burgess
!! @date June 2026
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE grad_shaf_bootstrap
USE oft_base
USE oft_gs, ONLY: gs_equil, flux_func, gsinv_interp, gs_factory, gs_psi2r, &
 gs_itor_nl, qprof_trace_tol, gs_fpol
USE oft_gs_profiles, ONLY: torflux_qgeom_backend, torflux_set_bounds
USE oft_gs_contour, ONLY: gs_sauter_contour, gs_fsa_contour
USE oft_lag_basis, ONLY: oft_blag_geval
USE oft_mesh_type, ONLY: bmesh_findcell
USE oft_blag_operators, ONLY: oft_lag_brinterp
USE tracing_2d, ONLY: set_tracer, active_tracer, tracinginv_fs
USE grad_shaf_prof_phys, ONLY: eval_dIdA, gs_ravgs, gs_ravgs_failed, gs_flux_int, &
  gs_flux_cumint, jphi_update, jphi_copy, jphi_flux_func, jphi_psi_nodes
USE oft_io, ONLY: hdf5_create_group, hdf5_write, hdf5_read, &
  hdf5_field_get_sizes, hdf5_field_exist
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Options of @ref jphi_bs_update (defaults apply unless set via @ref tokamaker_set_boot_ops)
!------------------------------------------------------------------------------
TYPE :: boot_ops
  LOGICAL :: initialized = .FALSE. !< Options have been explicitly set?
  LOGICAL :: isolate_edge_jbs = .FALSE. !< Isolate edge bootstrap spike from bulk?
  REAL(r8) :: scale_jbs = 1.0_r8 !< Scaling factor applied to the field-aligned bootstrap current, not to its p'G part (default 1)
  REAL(r8) :: djbs_tol = 1.0e-4_r8 !< Rel. change in jphi_bs below which it is frozen (not recomputed)
  LOGICAL :: diagnose_bs = .FALSE. !< Print alpha/Ip scalars, bootstrap stats, and full profile tables each NL iteration
  LOGICAL :: taper_edge_jphi = .FALSE. !< Taper all current components (jphi_ind, _bs, _pd, _fixed, _saw) to zero at the LCFS
  REAL(r8) :: taper_edge_psi0 = 0.999_r8 !< psi_N (standard: 0=axis, 1=LCFS) where edge taper begins
  INTEGER(i4) :: taper_edge_shape = 2 !< Edge taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
  REAL(r8) :: saw_q_s = 0.0_r8 !< Sawtooth reset q on axis (q_s); 0 = off (see @ref saw_reset_1d)
  REAL(r8) :: saw_dq = 0.03_r8 !< Reset regions end where q_base crosses saw_q_s + saw_dq
  REAL(r8) :: saw_tol = 1.0e-4_r8 !< Rel. change in jphi_saw (vs RMS total) below which jphi_saw is frozen
  REAL(r8) :: saw_relax = 1.0_r8 !< Under-relaxation of the jphi_saw update (1 = none)
  REAL(r8) :: saw_ramp = 0.01_r8 !< Reset weight ramps over this q deficit below saw_q_s (0 = hard trigger)
  INTEGER(i4) :: saw_rule = 2 !< Reset rule: 1 = fuse (axis to rho_m, FUSE saw_crash!), 2 = local (each dip, two-sided)
END TYPE boot_ops
!------------------------------------------------------------------------------
!> Current profiles of the last @ref jphi_bs_update (diagnostics; the update itself uses its own
!! `jphi_*` arrays). `_raw`: from the last bootstrap evaluation, before scaling, isolation and taper;
!! no suffix: as added into `jphi_total`. TokaMaker jphi [A/m²] on `psi_n` (OFT convention, 0=LCFS,
!! 1=axis; LCFS first); `jdotB_bs_raw` in [T A/m²].
!------------------------------------------------------------------------------
TYPE :: boot_profs
  REAL(r8), ALLOCATABLE :: psi_n(:) !< Normalised psi of the profile nodes, LCFS first
  REAL(r8), ALLOCATABLE :: jdotB_bs_raw(:) !< Redl PoP 2021 <j_BS.B>
  REAL(r8), ALLOCATABLE :: jphi_bs_raw(:) !< Field-aligned bootstrap jphi, <j_BS.B> F<1/R>/<B^2> (with jphi_pd_raw, eq. A7)
  REAL(r8), ALLOCATABLE :: jphi_pd_raw(:) !< Pressure-driven (diamagnetic + Pfirsch-Schlueter) part, P'(<R> - F^2<1/R>/<B^2>)
  REAL(r8), ALLOCATABLE :: jphi_total(:) !< Total = jphi_ind + jphi_bs + jphi_pd + jphi_saw + jphi_fixed
  REAL(r8), ALLOCATABLE :: jphi_bs(:) !< Field-aligned bootstrap jphi, scaled, optionally isolated and tapered
  REAL(r8), ALLOCATABLE :: jphi_pd(:) !< jphi_pd_raw, optionally tapered (not scaled or isolated)
  REAL(r8), ALLOCATABLE :: jphi_ind(:) !< Inductive jphi: input jphi times alpha, optionally tapered
  REAL(r8), ALLOCATABLE :: jphi_saw(:) !< Sawtooth jphi: input jphi_saw (optionally tapered) plus the q reset current
  REAL(r8), ALLOCATABLE :: jphi_fixed(:) !< Fixed (non-rescaled) jphi, optionally tapered
  REAL(r8) :: saw_rho_m = 0.0_r8 !< Sawtooth radius (rho_tor_norm): q_s + dq crossing beyond the outermost q_base < q_s (0 = none)
  REAL(r8) :: saw_rho_out = 0.0_r8 !< Outer end of the last reset (rho_tor_norm; 0 = none)
  INTEGER(i4) :: saw_n_dips = 0 !< Regions reset at the last reset (rule 1: separate q_base < saw_q_s regions)
END TYPE boot_profs
!> Names of the @ref boot_profs arrays, in declaration order (HDF5 datasets are these in upper case)
CHARACTER(LEN=16), PARAMETER :: boot_prof_names(10) = [CHARACTER(LEN=16) :: 'psi_n', 'jdotB_bs_raw', &
  'jphi_bs_raw', 'jphi_pd_raw', 'jphi_total', 'jphi_bs', 'jphi_pd', 'jphi_ind', &
  'jphi_saw', 'jphi_fixed']
REAL(r8), PARAMETER :: djbs_stol = 1.0e-3_r8 !< Stagnation tolerance: a stalled change above this resets the no-improve counter
REAL(r8), PARAMETER :: dalpha_warn_tol = 1.0e-6_r8 !< Relative (|dalpha|/|alpha|): no alpha stall warning below this
INTEGER(i4), PARAMETER :: saw_nsmooth = 8 !< Half-width (kept nodes) of the quadratic fit smoothing q_base in @ref saw_redistribute
!------------------------------------------------------------------------------
!> Jphi flux function type for bootstrap current calculation
!------------------------------------------------------------------------------
TYPE, EXTENDS(jphi_flux_func) :: jphi_bs_flux_func
  REAL(r8) :: alpha_last = 1.d0 !< Alpha (input Jphi rescaling factor) from previous NL iteration
  LOGICAL :: freeze_jbs = .FALSE. !< jphi_bs converged or stagnated: skip the Sauter pass
  INTEGER(i4) :: djbs_no_improve = 0   !< Consecutive steps with non-decreasing djbs
  REAL(r8) :: djbs_min = HUGE(1.0d0)  !< Running minimum djbs seen so far
  INTEGER(i4) :: dalpha_no_improve = 0 !< Consecutive steps with non-decreasing dalpha
  REAL(r8) :: dalpha_last = HUGE(1.0d0) !< dalpha of the previous update (stall diagnostic)
  REAL(r8), ALLOCATABLE :: jphi_bs_last(:) !< Scaled, untapered jphi_bs of the last evaluation (mu0*A/m²)
  REAL(r8), ALLOCATABLE :: jphi_pd_last(:) !< Its pressure-driven part (mu0*A/m²)
  LOGICAL :: freeze_saw = .FALSE. !< Set .TRUE. once jphi_saw converges or stagnates, after jphi_bs is frozen
  INTEGER(i4) :: djsaw_no_improve = 0 !< Consecutive steps with non-decreasing djsaw
  REAL(r8) :: djsaw_min = HUGE(1.0d0) !< Running minimum djsaw
  REAL(r8), ALLOCATABLE :: jphi_saw_last(:) !< jphi_saw of the last update (mu0*A/m²)
  REAL(r8), ALLOCATABLE :: jphi_total_last(:) !< jphi_total of the last update (mu0*A/m²)
  LOGICAL :: ffp_restored = .FALSE. !< F*F' knots were read from a file, so loading keeps them instead of updating
  TYPE(boot_ops) :: boot_ops       !< Options set from Python (@ref tokamaker_set_boot_ops)
  TYPE(boot_profs) :: boot_profs   !< Current profiles of the last update
CONTAINS
  !> Save to HDF5 (profile definition, options and the last profiles)
  PROCEDURE :: save_hdf5 => jphi_bs_save_hdf5
  !> Save the profile definition to a text file
  PROCEDURE :: save_txt => jphi_bs_save_txt
  !> Load from HDF5
  PROCEDURE :: load_hdf5 => jphi_bs_load_hdf5
  !> Load the profile definition from a text file
  PROCEDURE :: load_txt => jphi_bs_load_txt
  !> Copy into a new profile
  PROCEDURE :: copy => jphi_bs_copy
  !> Free the profile's storage
  PROCEDURE :: delete => jphi_bs_delete
  !> Update F*F' profile from Jphi and current equilibrium state
  PROCEDURE :: update => jphi_bs_update
END TYPE jphi_bs_flux_func
!------------------------------------------------------------------------------
!> Interpolation object for computing Sauter trapped-particle factors
!! (accumulates flux-surface averages of B, R, and bounce integrals)
!------------------------------------------------------------------------------
TYPE, EXTENDS(gsinv_interp) :: sauter_interp
  REAL(r8) :: bmax = 0.d0        !< Maximum |B| on surface (sampled during the trace)
  REAL(r8) :: rmax_surf = -1.d30 !< Maximum R on surface (sampled during the trace)
  REAL(r8) :: rmin_surf =  1.d30 !< Minimum R on surface (sampled during the trace)
  REAL(r8) :: mag_axis(2) = 0.d0 !< Magnetic axis (R, Z)
CONTAINS
  !> Evaluate ODE RHS for Sauter integration
  PROCEDURE :: interp => sauter_apply
END TYPE sauter_interp
CONTAINS
!------------------------------------------------------------------------------
!> Save to HDF5: the profile definition, the F*F' knots (so a load restores the profile exactly),
!! the options and the last profiles
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_save_hdf5(self,filename,path)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(LEN=*), INTENT(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('jphi-split-bootstrap',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%x,filename,path//'/XVALS')
CALL hdf5_write(self%jphi,filename,path//'/YVALS')
CALL hdf5_write(self%j0,filename,path//'/J0')
CALL hdf5_write(self%yp,filename,path//'/YPVALS')
CALL hdf5_write(self%y0,filename,path//'/Y0')
IF(self%boot_ops%initialized)THEN
  CALL hdf5_create_group(filename,path//'/BOOT_OPS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%isolate_edge_jbs),filename,path//'/BOOT_OPS/ISOLATE_EDGE_JBS')
  CALL hdf5_write(self%boot_ops%scale_jbs,filename,path//'/BOOT_OPS/SCALE_JBS')
  CALL hdf5_write(self%boot_ops%djbs_tol,filename,path//'/BOOT_OPS/DJBS_TOL')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%diagnose_bs),filename,path//'/BOOT_OPS/DIAGNOSE_BS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%taper_edge_jphi),filename,path//'/BOOT_OPS/TAPER_EDGE_JPHI')
  CALL hdf5_write(self%boot_ops%taper_edge_psi0,filename,path//'/BOOT_OPS/TAPER_EDGE_PSI0')
  CALL hdf5_write(self%boot_ops%taper_edge_shape,filename,path//'/BOOT_OPS/TAPER_EDGE_SHAPE')
  CALL hdf5_write(self%boot_ops%saw_q_s,filename,path//'/BOOT_OPS/SAW_Q_S')
  CALL hdf5_write(self%boot_ops%saw_dq,filename,path//'/BOOT_OPS/SAW_DQ')
  CALL hdf5_write(self%boot_ops%saw_tol,filename,path//'/BOOT_OPS/SAW_TOL')
  CALL hdf5_write(self%boot_ops%saw_relax,filename,path//'/BOOT_OPS/SAW_RELAX')
  CALL hdf5_write(self%boot_ops%saw_ramp,filename,path//'/BOOT_OPS/SAW_RAMP')
  CALL hdf5_write(self%boot_ops%saw_rule,filename,path//'/BOOT_OPS/SAW_RULE')
END IF
CALL boot_profs_save(self%boot_profs,filename,path//'/BOOT_PROFS')
END SUBROUTINE jphi_bs_save_hdf5
!------------------------------------------------------------------------------
!> Save the profile definition to a text file
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_save_txt(self,io_unit)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
INTEGER, INTENT(in) :: io_unit
WRITE(io_unit,*)'jphi-split-bootstrap'
WRITE(io_unit,*)self%npsi,self%j0
WRITE(io_unit,*)self%x
WRITE(io_unit,*)self%jphi
END SUBROUTINE jphi_bs_save_txt
!------------------------------------------------------------------------------
!> Load from HDF5 (see @ref jphi_bs_save_hdf5); options and profiles are optional
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_load_hdf5(self,filename,path,success)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(LEN=*), INTENT(in) :: path
LOGICAL, INTENT(out) :: success
INTEGER(i4) :: npsi
INTEGER(i4) :: int_tmp, ierr
REAL(r8) :: J0
REAL(r8), ALLOCATABLE :: xvals(:),yvals(:)
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
!---F*F' knots (files written before they were saved get the update on load instead)
self%ffp_restored=.FALSE.
IF(hdf5_field_exist(filename,path//'/YPVALS').AND.hdf5_field_exist(filename,path//'/Y0'))THEN
  CALL hdf5_read(self%y0,filename,path//'/Y0',success=success)
  IF(success)CALL hdf5_read(self%yp,filename,path//'/YPVALS',success=success)
  IF(success)THEN
    ierr=self%set_cofs(self%yp)
    self%ffp_restored=.TRUE.
  END IF
END IF
!---Options
self%boot_ops%initialized = .FALSE.
IF(hdf5_field_exist(filename,path//'/BOOT_OPS'))THEN
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/ISOLATE_EDGE_JBS',success=success)
  IF(success) self%boot_ops%isolate_edge_jbs = (int_tmp/=0)
  IF(hdf5_field_exist(filename,path//'/BOOT_OPS/PARAMETERIZE_JBS'))THEN
    CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/PARAMETERIZE_JBS',success=success)
    IF(success.AND.(int_tmp/=0))CALL oft_warn('jphi_bs_load: "parameterize_jbs" is no longer supported, ignored')
  END IF
  CALL hdf5_read(self%boot_ops%scale_jbs,filename,path//'/BOOT_OPS/SCALE_JBS',success=success)
  CALL hdf5_read(self%boot_ops%djbs_tol,filename,path//'/BOOT_OPS/DJBS_TOL',success=success)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/DIAGNOSE_BS',success=success)
  IF(success) self%boot_ops%diagnose_bs = (int_tmp/=0)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/TAPER_EDGE_JPHI',success=success)
  IF(success) self%boot_ops%taper_edge_jphi = (int_tmp/=0)
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
CALL boot_profs_load(self%boot_profs,filename,path//'/BOOT_PROFS')
success=.TRUE. ! BOOT_OPS and BOOT_PROFS are optional
END SUBROUTINE jphi_bs_load_hdf5
!------------------------------------------------------------------------------
!> Load the profile definition from a text file
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_load_txt(self,io_unit)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
INTEGER, INTENT(in) :: io_unit
INTEGER(i4) :: npsi
REAL(r8) :: J0
REAL(r8), ALLOCATABLE :: xvals(:),yvals(:)
READ(io_unit,*)npsi,J0
ALLOCATE(xvals(npsi),yvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)yvals
CALL create_jphi_bs_ff(self,npsi,xvals,yvals,J0)
DEALLOCATE(xvals,yvals)
END SUBROUTINE jphi_bs_load_txt
!------------------------------------------------------------------------------
!> Sets values inside jphi_bs_flux_func. Called only after a fresh allocation.
!------------------------------------------------------------------------------
SUBROUTINE create_jphi_bs_ff(func,npsi,psivals,yvals,y0)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(i4), INTENT(in) :: npsi
REAL(r8), INTENT(in) :: psivals(npsi)
REAL(r8), INTENT(in) :: yvals(npsi)
REAL(r8), INTENT(in) :: y0
INTEGER(i4) :: i,ierr
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
CLASS DEFAULT
  CALL oft_abort('Invalid flux function type in create_jphi_bs_ff','create_jphi_bs_ff',__FILE__)
END SELECT

END SUBROUTINE create_jphi_bs_ff
!------------------------------------------------------------------------------
!> Copy into a new profile (definition, options, solve state and last profiles)
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_copy(self,new)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(flux_func), POINTER, INTENT(inout) :: new
CALL jphi_copy(self,new)
SELECT TYPE(new)
  CLASS IS(jphi_bs_flux_func)
    new%alpha_last = self%alpha_last
    new%freeze_jbs = self%freeze_jbs
    new%djbs_no_improve = self%djbs_no_improve
    new%djbs_min = self%djbs_min
    new%dalpha_no_improve = self%dalpha_no_improve
    new%dalpha_last = self%dalpha_last
    IF(ALLOCATED(self%jphi_bs_last))new%jphi_bs_last = self%jphi_bs_last
    IF(ALLOCATED(self%jphi_pd_last))new%jphi_pd_last = self%jphi_pd_last
    new%freeze_saw = self%freeze_saw
    new%djsaw_no_improve = self%djsaw_no_improve
    new%djsaw_min = self%djsaw_min
    IF(ALLOCATED(self%jphi_saw_last))new%jphi_saw_last = self%jphi_saw_last
    IF(ALLOCATED(self%jphi_total_last))new%jphi_total_last = self%jphi_total_last
    new%boot_ops = self%boot_ops
    new%boot_profs = self%boot_profs
END SELECT
END SUBROUTINE jphi_bs_copy
!------------------------------------------------------------------------------
!> Free the profile's storage
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_delete(self)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
self%j0=0.d0
IF(ASSOCIATED(self%jphi))DEALLOCATE(self%jphi)
IF(ASSOCIATED(self%x))DEALLOCATE(self%x)
IF(ASSOCIATED(self%yp))DEALLOCATE(self%yp)
IF(ASSOCIATED(self%y))DEALLOCATE(self%y)
CALL reset_solve_state(self)
self%boot_profs = boot_profs()
END SUBROUTINE jphi_bs_delete
!------------------------------------------------------------------------------
!> Set the F*F' knots from the equilibrium's current F*F' (`gseq%I` times `ffp_scale`), so a bootstrap
!! profile that replaces it keeps the current state until the next solve. No-op before an equilibrium exists.
!------------------------------------------------------------------------------
SUBROUTINE jphi_bs_carry_ffp(self,gseq)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self !< New bootstrap profile
CLASS(gs_equil), INTENT(inout) :: gseq !< Equilibrium, still holding the profile being replaced
INTEGER(i4) :: i,ierr
REAL(r8) :: b1,b2
REAL(r8), ALLOCATABLE :: xpsi(:)
LOGICAL :: mapped
IF(.NOT.ASSOCIATED(gseq%I))RETURN
b1=gseq%plasma_bounds(1); b2=gseq%plasma_bounds(2)
IF((b1<-1.d98).OR.(b2==b1))RETURN ! No plasma yet
CALL jphi_psi_nodes(self,xpsi,mapped) ! Knots in psi (mapped from toroidal flux if on it)
IF(.NOT.mapped)THEN ! No map yet (eg. replacing a profile on psi_N): build it for the current equilibrium
  CALL gseq%torflux%update(gseq)
  CALL jphi_psi_nodes(self,xpsi,mapped)
END IF
IF(.NOT.mapped)THEN
  CALL oft_warn("No toroidal flux map, F*F' not carried to the new profile")
  RETURN
END IF
self%y0=gseq%ffp_scale*gseq%I%fp(b1)
DO i=1,self%npsi
  self%yp(i)=gseq%ffp_scale*gseq%I%fp(xpsi(i)*(b2-b1)+b1)
END DO
IF(ASSOCIATED(self%phi_wrapper))THEN
  CALL torflux_set_bounds(self%phi_wrapper,gseq%plasma_bounds)
  ierr=self%phi_wrapper%set_cofs(self%yp)
ELSE
  self%plasma_bounds=gseq%plasma_bounds
  ierr=self%set_cofs(self%yp)
END IF
gseq%ffp_scale=1.d0
END SUBROUTINE jphi_bs_carry_ffp
!------------------------------------------------------------------------------
!> Inputs of a jphi-split-bootstrap solve that `gseq` lacks (kinetic profiles, Ip/pax targets), empty if none
!------------------------------------------------------------------------------
FUNCTION jphi_bs_missing_inputs(gseq) RESULT(missing)
CLASS(gs_equil), INTENT(in) :: gseq !< Equilibrium
CHARACTER(LEN=256) :: missing
missing=''
IF(.NOT.ASSOCIATED(gseq%Te))missing=TRIM(missing)//' Te'
IF(.NOT.ASSOCIATED(gseq%Ti))missing=TRIM(missing)//' Ti'
IF(.NOT.ASSOCIATED(gseq%ne))missing=TRIM(missing)//' ne'
IF(.NOT.ASSOCIATED(gseq%ni))missing=TRIM(missing)//' ni'
IF(.NOT.ASSOCIATED(gseq%Zeff))missing=TRIM(missing)//' Zeff'
IF(LEN_TRIM(missing)>0)missing='kinetic profiles (set_kinetic_profiles), missing:'//TRIM(missing)
IF((gseq%Ip_target<0.d0).OR.(gseq%pax_target<0.d0))THEN
  IF(LEN_TRIM(missing)>0)missing=TRIM(missing)//'; '
  missing=TRIM(missing)//'Ip and pax targets (set_targets)'
END IF
END FUNCTION jphi_bs_missing_inputs
!------------------------------------------------------------------------------
!> Forget the convergence state of the last solve (freeze flags, counters, cached currents)
!------------------------------------------------------------------------------
SUBROUTINE reset_solve_state(self)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
self%freeze_jbs = .FALSE.
self%djbs_no_improve = 0
self%djbs_min = HUGE(1.0d0)
self%dalpha_no_improve = 0
self%dalpha_last = HUGE(1.0d0)
IF(ALLOCATED(self%jphi_bs_last))DEALLOCATE(self%jphi_bs_last)
IF(ALLOCATED(self%jphi_pd_last))DEALLOCATE(self%jphi_pd_last)
self%freeze_saw = .FALSE.
self%djsaw_no_improve = 0
self%djsaw_min = HUGE(1.0d0)
IF(ALLOCATED(self%jphi_saw_last))DEALLOCATE(self%jphi_saw_last)
IF(ALLOCATED(self%jphi_total_last))DEALLOCATE(self%jphi_total_last)
CALL saw_scalars_reset(self%boot_profs)
END SUBROUTINE reset_solve_state
!------------------------------------------------------------------------------
!> Zero the sawtooth scalars of @ref boot_profs (no reset)
!------------------------------------------------------------------------------
SUBROUTINE saw_scalars_reset(bp)
TYPE(boot_profs), INTENT(inout) :: bp
bp%saw_rho_m = 0.0_r8
bp%saw_rho_out = 0.0_r8
bp%saw_n_dips = 0
END SUBROUTINE saw_scalars_reset
!------------------------------------------------------------------------------
!> Set (`set`) or point at (`get`) the @ref boot_profs array called `name` (see @ref boot_prof_names);
!! `get` is NULL if that array is unallocated or the name is unknown
!------------------------------------------------------------------------------
SUBROUTINE boot_prof_access(bp,name,get,set)
TYPE(boot_profs), TARGET, INTENT(inout) :: bp
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), POINTER, OPTIONAL, INTENT(out) :: get(:)
REAL(r8), OPTIONAL, INTENT(in) :: set(:)
IF(PRESENT(get))get=>NULL()
SELECT CASE(TRIM(name))
  CASE('psi_n'); CALL access(bp%psi_n)
  CASE('jdotB_bs_raw'); CALL access(bp%jdotB_bs_raw)
  CASE('jphi_bs_raw'); CALL access(bp%jphi_bs_raw)
  CASE('jphi_pd_raw'); CALL access(bp%jphi_pd_raw)
  CASE('jphi_total'); CALL access(bp%jphi_total)
  CASE('jphi_bs'); CALL access(bp%jphi_bs)
  CASE('jphi_pd'); CALL access(bp%jphi_pd)
  CASE('jphi_ind'); CALL access(bp%jphi_ind)
  CASE('jphi_saw'); CALL access(bp%jphi_saw)
  CASE('jphi_fixed'); CALL access(bp%jphi_fixed)
END SELECT
CONTAINS
SUBROUTINE access(a)
REAL(r8), ALLOCATABLE, TARGET, INTENT(inout) :: a(:)
IF(PRESENT(set))a = set
IF(PRESENT(get).AND.ALLOCATED(a))get => a
END SUBROUTINE access
END SUBROUTINE boot_prof_access
!------------------------------------------------------------------------------
!> Write the allocated @ref boot_profs arrays to group `path` (dataset names: upper-case array names)
!------------------------------------------------------------------------------
SUBROUTINE boot_profs_save(bp,filename,path)
TYPE(boot_profs), TARGET, INTENT(inout) :: bp
CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(LEN=*), INTENT(in) :: path
REAL(r8), POINTER :: vals(:)
INTEGER(i4) :: k
LOGICAL :: group_made
group_made = .FALSE.
DO k=1,SIZE(boot_prof_names)
  CALL boot_prof_access(bp,boot_prof_names(k),get=vals)
  IF(.NOT.ASSOCIATED(vals))CYCLE
  IF(.NOT.group_made)CALL hdf5_create_group(filename,path)
  group_made = .TRUE.
  CALL hdf5_write(vals,filename,path//'/'//h5_name(boot_prof_names(k)))
END DO
IF(group_made)THEN
  CALL hdf5_write(bp%saw_rho_m,filename,path//'/SAW_RHO_M')
  CALL hdf5_write(bp%saw_rho_out,filename,path//'/SAW_RHO_OUT')
  CALL hdf5_write(bp%saw_n_dips,filename,path//'/SAW_N_DIPS')
END IF
END SUBROUTINE boot_profs_save
!------------------------------------------------------------------------------
!> Read the @ref boot_profs arrays present in group `path`; the others are left unallocated
!------------------------------------------------------------------------------
SUBROUTINE boot_profs_load(bp,filename,path)
TYPE(boot_profs), TARGET, INTENT(inout) :: bp
CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(LEN=*), INTENT(in) :: path
REAL(r8), ALLOCATABLE :: vals(:)
INTEGER(i4), ALLOCATABLE :: dim_sizes(:)
INTEGER(i4) :: k, ndims
LOGICAL :: success
bp = boot_profs()
IF(.NOT.hdf5_field_exist(filename,path))RETURN
DO k=1,SIZE(boot_prof_names)
  IF(.NOT.hdf5_field_exist(filename,path//'/'//h5_name(boot_prof_names(k))))CYCLE
  CALL hdf5_field_get_sizes(filename,path//'/'//h5_name(boot_prof_names(k)),ndims,dim_sizes)
  ALLOCATE(vals(dim_sizes(1)))
  CALL hdf5_read(vals,filename,path//'/'//h5_name(boot_prof_names(k)),success=success)
  IF(success)CALL boot_prof_access(bp,boot_prof_names(k),set=vals)
  DEALLOCATE(vals,dim_sizes)
END DO
IF(hdf5_field_exist(filename,path//'/SAW_N_DIPS'))THEN
  CALL hdf5_read(bp%saw_rho_m,filename,path//'/SAW_RHO_M',success=success)
  CALL hdf5_read(bp%saw_rho_out,filename,path//'/SAW_RHO_OUT',success=success)
  CALL hdf5_read(bp%saw_n_dips,filename,path//'/SAW_N_DIPS',success=success)
END IF
END SUBROUTINE boot_profs_load
!------------------------------------------------------------------------------
!> HDF5 dataset name of a @ref boot_profs array (its name in upper case)
!------------------------------------------------------------------------------
FUNCTION h5_name(name) RESULT(h5)
CHARACTER(LEN=*), INTENT(in) :: name
CHARACTER(LEN=LEN_TRIM(name)) :: h5
h5 = TRIM(name)
CALL string_to_upper(h5)
END FUNCTION h5_name
!---------------------------------------------------------------------------------
!> Update F*F' profile from the inductive jphi coupled with the bootstrap current.
!!
!! Called once per NL iteration (the first iteration of a solve runs @ref jphi_update and resets the
!! convergence state). Current densities are TokaMaker jphi in mu0*A/m² on the nodes `psi_n`
!! (LCFS and the profile nodes `xpsi`):
!!   1. Pressure scale.
!!   2. jphi_fixed, the input sawtooth current jphi_saw_in and the bootstrap current jphi_bs (here
!!      with its pressure-driven part jphi_pd, which isolation and scale_jbs leave out), optionally
!!      isolated and scaled, or reused from the last evaluation when frozen; <R>, <1/R>, <1/R^2> on the nodes
!!      from its Sauter pass, else (and on the LCFS) from gs_ravgs.
!!   3. Edge taper of jphi_bs, jphi_pd, jphi_ind, jphi_fixed and jphi_saw_in.
!!   4. alpha from the exact I_p (gs_flux_int of eval_dIdA, eq. A9c), which is affine in
!!      alpha: alpha = (Ip_target - Ip_lo)/(Ip_hi - Ip_lo).
!!   4a. Sawtooth reset (saw_q_s > 0): jphi_saw = jphi_saw_in + the Ip-neutral current that resets the
!!      q of the base current alpha*jphi_ind + jphi_bs + jphi_fixed + jphi_saw_in above saw_q_s
!!      (@ref saw_redistribute).
!!   5. jphi_total = alpha*jphi_ind + jphi_bs + jphi_fixed + jphi_saw; boot_profs (jphi_bs without
!!      jphi_pd); F*F' knots.
!!   6. Diagnostics (if diagnose_bs is set).
!---------------------------------------------------------------------------------
SUBROUTINE jphi_bs_update(self, gseq)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4) :: i, ierr
REAL(r8) :: pscale, pprime
REAL(r8), ALLOCATABLE :: xpsi(:) !< Node locations in normalized poloidal flux
REAL(r8), ALLOCATABLE :: psi_n(:) !< Nodes: LCFS (0) and xpsi
REAL(r8), ALLOCATABLE :: ravgs(:,:) !< <R>, <1/R>, <1/R^2> at psi_n
REAL(r8), ALLOCATABLE :: jtor(:)  !< dI/dA for the exact I_p measure (eval_dIdA)
REAL(r8), ALLOCATABLE :: jphi_bs(:), jphi_pd(:), jphi_ind(:), jphi_fixed(:), jphi_total(:)
REAL(r8), ALLOCATABLE :: jphi_bs_spike(:) !< Isolated edge spike [A/m²]
!> Sawtooth current: input, value entering the alpha solve, result, reset increment
REAL(r8), ALLOCATABLE :: jphi_saw_in(:), jphi_saw_cur(:), jphi_saw(:), djphi_saw(:)
REAL(r8), ALLOCATABLE :: cflds(:,:), cum(:,:)
LOGICAL :: saw_on, do_saw
LOGICAL :: bs_frozen !< jphi_bs reused from the last evaluation (no Sauter pass)
REAL(r8) :: alpha, ip_target, ip_ind, ip_result_lo, ip_result_hi, dalpha
REAL(r8) :: djbs !< Relative RMS change of jphi_bs since the last evaluation
REAL(r8) :: djsaw !< Relative change of jphi_saw since the last update
REAL(r8) :: itor_nl, itor_flint
CHARACTER(len=256) :: char_buf
!--- First iteration of a solve: no bootstrap current yet
IF(.NOT. gseq%skip_targets) THEN
  CALL reset_solve_state(self)
  CALL jphi_update(self,gseq)
  self%alpha_last = self%norm_last ! Same scale of jphi_ind as alpha
  RETURN
ENDIF
!---
self%plasma_bounds = gseq%plasma_bounds
CALL jphi_psi_nodes(self,xpsi)
IF(gseq%mode/=1) &
  CALL oft_abort("Jphi-BS profile requires (F^2)' formulation", &
                 "jphi_bs_update",__FILE__)
char_buf = jphi_bs_missing_inputs(gseq) ! Normally caught before the solve (tokamaker_solve)
IF(LEN_TRIM(char_buf) > 0) &
  CALL oft_abort("Jphi-BS profile requires "//TRIM(char_buf), "jphi_bs_update",__FILE__)
!--- 1. Pressure scale.
ALLOCATE(jtor(0:self%npsi), psi_n(0:self%npsi), ravgs(0:self%npsi,3))
psi_n = [0.0_r8, xpsi]
saw_on = self%boot_ops%saw_q_s > 0.0_r8
do_saw = saw_on .AND. (.NOT.self%freeze_saw) .AND. ALLOCATED(self%jphi_total_last)
djsaw = 0.0_r8
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
IF(ASSOCIATED(gseq%P_ani)) &
  CALL oft_abort('Jphi profiles do not support anisotropic pressure', &
                 'jphi_bs_update',__FILE__)
pscale = gseq%P%f(gseq%plasma_bounds(2))
pscale = gseq%pax_target / pscale
!--- 2. Fixed, input sawtooth and bootstrap currents
ALLOCATE(jphi_fixed(0:self%npsi), jphi_saw_in(0:self%npsi), jphi_bs(0:self%npsi), jphi_pd(0:self%npsi))
CALL node_vals(gseq%jphi_fixed, jphi_fixed)
CALL node_vals(gseq%jphi_saw, jphi_saw_in)
bs_frozen = self%freeze_jbs .AND. ALLOCATED(self%jphi_bs_last)
IF(bs_frozen) THEN
  jphi_bs = self%jphi_bs_last
  jphi_pd = self%jphi_pd_last
  djbs = 0.0_r8
ELSE
  CALL calculate_bootstrap(self, gseq, self%npsi, xpsi, jphi_bs, jphi_pd, ierr, ravgs=ravgs(1:,:))
  IF(ierr /= 0)THEN
    gseq%prof_failed=.TRUE. ! The solve ends with error -10
    RETURN
  END IF
  IF(self%boot_ops%isolate_edge_jbs)THEN
    ALLOCATE(jphi_bs_spike(0:self%npsi))
    CALL edge_spike(self%npsi, xpsi, jphi_bs - jphi_pd, jphi_bs_spike) ! field-aligned part only
  END IF
  IF(self%boot_ops%diagnose_bs)THEN
    IF(self%boot_ops%isolate_edge_jbs)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         jphi_bs(bulk)[A/m2]  spike[A/m2]   jphi[A/m2]  jphi_fixed[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,5ES15.5)') '  ', i, xpsi(i), jphi_bs(i) - jphi_pd(i), jphi_bs_spike(i), self%jphi(i), jphi_fixed(i)/mu0
      END DO
    ELSE
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         jphi_bs[A/m2]      jphi[A/m2]  jphi_fixed[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,4ES15.5)') '  ', i, xpsi(i), jphi_bs(i), self%jphi(i), jphi_fixed(i)/mu0
      END DO
    END IF
  END IF
  IF(self%boot_ops%isolate_edge_jbs)THEN
    jphi_bs = jphi_bs_spike + jphi_pd
    DEALLOCATE(jphi_bs_spike)
  END IF
  !--- scale_jbs scales the field-aligned part only, not p'G
  jphi_bs = jphi_pd + self%boot_ops%scale_jbs*(jphi_bs - jphi_pd)
  jphi_bs = jphi_bs * mu0
  jphi_pd = jphi_pd * mu0
  !--- Freeze when the change drops below djbs_tol or stalls
  IF(ALLOCATED(self%jphi_bs_last)) THEN
    djbs = SQRT(SUM((jphi_bs - self%jphi_bs_last)**2) / REAL(self%npsi+1,r8)) / &
           MAX(SQRT(SUM(jphi_bs**2) / REAL(self%npsi+1,r8)), 1.0e-30_r8)
    CALL stall_check(djbs, self%boot_ops%djbs_tol, djbs_stol, self%djbs_min, self%djbs_no_improve, &
                     self%freeze_jbs, 'Bootstrap solution', 'djbs_tol', 'bootstrap solution')
  ELSE
    djbs = HUGE(1.0_r8)
  END IF
  self%jphi_bs_last = jphi_bs
  self%jphi_pd_last = jphi_pd
END IF
!--- <R>, <1/R>, <1/R^2> (I_p measure, F*F' map) not set by the Sauter pass
IF(bs_frozen)THEN
  CALL gs_ravgs(gseq, self%npsi+1, psi_n, ravgs)
ELSE
  CALL gs_ravgs(gseq, 1, psi_n(0:0), ravgs(0:0,:))
END IF
IF(gs_ravgs_failed(gseq, psi_n, ravgs, 'jphi_bs_update'))RETURN
!--- 3. Edge taper (the cache above keeps the untapered jphi_bs, so freeze checks track physics)
ALLOCATE(jphi_ind(0:self%npsi))
jphi_ind = [self%j0, self%jphi]
CALL taper(jphi_bs)
CALL taper(jphi_pd)
CALL taper(jphi_ind)
CALL taper(jphi_fixed)
CALL taper(jphi_saw_in)
ALLOCATE(jphi_total(0:self%npsi), jphi_saw_cur(0:self%npsi), jphi_saw(0:self%npsi))
!   Sawtooth current entering the alpha solve: last iterate (or the input until the first reset)
jphi_saw_cur = jphi_saw_in
IF(saw_on .AND. ALLOCATED(self%jphi_saw_last))jphi_saw_cur = self%jphi_saw_last
!   One quadrature pass for the alpha solve: lo (alpha=0 with the input saw), ind (inductive, linear),
!   dsaw (jphi_saw_cur - input, linear); eq (last total) is for the reset.
IF(do_saw)THEN
  ALLOCATE(cflds(self%npsi+1,4), cum(self%npsi+1,3))
  jphi_total = jphi_bs + jphi_fixed + jphi_saw_in
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, cflds(:,1))
  jphi_total = jphi_ind + jphi_bs + jphi_fixed + jphi_saw_in
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  cflds(:,2) = jtor - cflds(:,1)
  jphi_total = jphi_bs + jphi_fixed + jphi_saw_cur
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  cflds(:,3) = jtor - cflds(:,1)
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, self%jphi_total_last, pscale, cflds(:,4))
  CALL gs_flux_cumint(gseq, psi_n, cflds(:,1:3), self%npsi+1, 3, cum, extrap=1)
END IF
!--- 4. alpha: I_p is affine in alpha, so alpha = (Ip_target - Ip_lo) / (Ip_hi - Ip_lo). Re-solved
!   every update (cheap): a frozen alpha lets I_p drift while the shape is still settling.
ip_target = ABS(gseq%Ip_target)
IF(do_saw)THEN
  ip_result_lo = cum(1,1) + cum(1,3)
  ip_result_hi = ip_result_lo + cum(1,2)
ELSE
  jphi_total = jphi_bs + jphi_fixed + jphi_saw_cur
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_lo, extrap=1)
  jphi_total = jphi_ind + jphi_bs + jphi_fixed + jphi_saw_cur
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_hi, extrap=1)
END IF
ip_ind = ip_result_hi - ip_result_lo
IF(ABS(ip_ind) > 0.0_r8)THEN
  alpha = (ip_target - ip_result_lo) / ip_ind
ELSE
  CALL oft_warn("jphi_bs_update: inductive jphi carries no (finite) net current, keeping the last alpha")
  alpha = self%alpha_last
END IF
!--- Relative change: alpha's scale is set by the units of ffp_prof (e.g. ~1e-6 for a jphi in A/m^2)
dalpha = ABS(alpha - self%alpha_last) / MAX(ABS(alpha), TINY(1.0_r8))
!--- Stall diagnostic only (alpha is never frozen): 2 consecutive non-decreasing steps.
IF(dalpha >= dalpha_warn_tol .AND. dalpha >= self%dalpha_last) THEN
  self%dalpha_no_improve = self%dalpha_no_improve + 1
  IF(self%dalpha_no_improve >= 2) THEN
    WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
      'Alpha convergence stalled,' // &
      ' relative change per nonlinear step = ', dalpha, &
      ', above recommended tolerance (dalpha_warn_tol=', dalpha_warn_tol, ')'
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
!--- 4a. Sawtooth reset of the base current (alpha now known); jphi_saw freeze check.
IF(do_saw)THEN
  ALLOCATE(djphi_saw(0:self%npsi))
  CALL saw_redistribute(self, gseq, psi_n, ravgs, cflds(:,1) + alpha*cflds(:,2), &
                        cflds(:,4), alpha*jphi_ind + jphi_bs + jphi_fixed + jphi_saw_in, djphi_saw)
  jphi_saw = jphi_saw_in + self%boot_ops%saw_relax*djphi_saw + (1.0_r8 - self%boot_ops%saw_relax)*(jphi_saw_cur - jphi_saw_in)
  jphi_total = alpha * jphi_ind + jphi_bs + jphi_fixed + jphi_saw
  djsaw = SQRT(SUM((jphi_saw - self%jphi_saw_last)**2)) / MAX(SQRT(SUM(jphi_total**2)), 1.0e-30_r8)
  ! Freeze only once jphi_bs has: before that the base current still moves, and an early iterate
  ! without a dip (djsaw = 0) would freeze the reset off for the whole solve
  IF(self%freeze_jbs)CALL stall_check(djsaw, self%boot_ops%saw_tol, djbs_stol, self%djsaw_min, &
    self%djsaw_no_improve, self%freeze_saw, 'Sawtooth current', 'saw_tol', 'sawtooth current')
  DEALLOCATE(djphi_saw, cflds, cum)
ELSE
  jphi_saw = jphi_saw_cur
END IF
IF(saw_on)self%jphi_saw_last = jphi_saw
!--- 5. Total current, profiles (A/m²), F*F' knots
jphi_total = alpha * jphi_ind + jphi_bs + jphi_fixed + jphi_saw
IF(saw_on)self%jphi_total_last = jphi_total
self%boot_profs%psi_n = psi_n
self%boot_profs%jphi_total = jphi_total/mu0
self%boot_profs%jphi_bs = (jphi_bs - jphi_pd)/mu0
self%boot_profs%jphi_pd = jphi_pd/mu0
self%boot_profs%jphi_ind = alpha * jphi_ind/mu0
self%boot_profs%jphi_saw = jphi_saw/mu0
self%boot_profs%jphi_fixed = jphi_fixed/mu0
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
!--- 6. Diagnostics.
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_target   = ', ip_target
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_lo= ', ip_result_lo
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_hi= ', ip_result_hi
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] alpha       = ', alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] dalpha      = ', dalpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djbs        = ', djbs
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] bs_frozen   = ', self%freeze_jbs
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_bs max = ', MAXVAL(ABS(jphi_bs))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi max    = ', MAXVAL(ABS(self%jphi))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_fixed max = ', MAXVAL(ABS(jphi_fixed))
  IF(saw_on)THEN
    WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_saw max = ', MAXVAL(ABS(jphi_saw))/mu0
    WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djsaw       = ', djsaw
    WRITE(*,'(A,L1)')     '  [jphi_bs_update] saw_frozen  = ', self%freeze_saw
  END IF
  !--- Side-by-side Ip comparison: FEM nonlinear solve vs profile flux integral
  CALL gs_itor_nl(gseq, itor_nl)
  CALL eval_dIdA(gseq, psi_n, self%npsi+1, ravgs, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, itor_flint, extrap=1)
  WRITE(*,'(A)') '  [jphi_bs_update] --- Ip comparison (current jphi_total) ---'
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(gs_itor_nl)    = ', itor_nl/mu0
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(flux_int)      = ', itor_flint/mu0
END IF
i=self%set_cofs(self%yp)
CONTAINS
!> Flux function `F` on the nodes times mu0 (zero if `F` is not set)
SUBROUTINE node_vals(F, a)
CLASS(flux_func), POINTER, INTENT(in) :: F
REAL(r8), INTENT(out) :: a(0:)
INTEGER(i4) :: k
a = 0.0_r8
IF(.NOT.ASSOCIATED(F))RETURN
DO k = 0, self%npsi
  a(k) = F%fp(psi_n(k))*mu0
END DO
END SUBROUTINE node_vals
!> Edge taper (taper_edge_psi0 is in standard convention, psi_n in OFT convention)
SUBROUTINE taper(a)
REAL(r8), INTENT(inout) :: a(0:)
IF(self%boot_ops%taper_edge_jphi)CALL apply_edge_taper(self%npsi+1, psi_n, a, &
  1.0_r8 - self%boot_ops%taper_edge_psi0, self%boot_ops%taper_edge_shape)
END SUBROUTINE taper
END SUBROUTINE jphi_bs_update
!---------------------------------------------------------------------------------
!> Freeze test of an iterated current: freeze when its relative change `d` drops below `tol`; after
!! two updates without a new minimum, warn, and freeze if `d` <= `stol` (else keep trying)
!---------------------------------------------------------------------------------
SUBROUTINE stall_check(d, tol, stol, dmin, n_no_improve, frozen, what, tol_name, frozen_what)
REAL(r8), INTENT(in) :: d !< Relative change this update
REAL(r8), INTENT(in) :: tol !< Convergence tolerance
REAL(r8), INTENT(in) :: stol !< Stagnation tolerance
REAL(r8), INTENT(inout) :: dmin !< Smallest change so far
INTEGER(i4), INTENT(inout) :: n_no_improve !< Updates since the last new minimum
LOGICAL, INTENT(inout) :: frozen !< Set when frozen
CHARACTER(LEN=*), INTENT(in) :: what !< Name in the warning ("Bootstrap solution")
CHARACTER(LEN=*), INTENT(in) :: tol_name !< Name of `tol` in the warning
CHARACTER(LEN=*), INTENT(in) :: frozen_what !< Name in the freeze message ("bootstrap solution")
CHARACTER(LEN=256) :: char_buf
IF(d < tol)THEN
  frozen = .TRUE.
  IF(oft_env%pm)WRITE(*,*)' Freezing '//frozen_what//'.'
ELSE IF(d >= dmin)THEN
  n_no_improve = n_no_improve + 1
  IF(n_no_improve >= 2)THEN
    WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') what//' convergence stalled,' // &
      ' relative change per nonlinear step = ', d, ', above set tolerance ('//tol_name//'=', tol, ')'
    IF(d > stol)THEN
      n_no_improve = 0 ! Give more chances to improve
    ELSE
      frozen = .TRUE.
      char_buf = TRIM(char_buf) // ' Freezing '//frozen_what//'.'
    END IF
    IF(oft_env%pm)CALL oft_warn(TRIM(char_buf))
  END IF
ELSE
  n_no_improve = 0
  dmin = d
END IF
END SUBROUTINE stall_check
!---------------------------------------------------------------------------------
!> Sawtooth reset current on the nodes `xn` (OFT convention, mu0*A/m², TokaMaker jphi).
!!
!! q_base = q_eq*I_eq/I_base (q*I is geometric) on rho = sqrt(Phi_N); @ref saw_reset_1d resets it on
!! nodes with enough enclosed area, interpolated in rho; dI/dA maps back to jphi via eval_dIdA.
!! Enclosed area and currents are 1-D integrals over the contour \f$ dA/d\hat{\psi} \f$: the reset
!! differentiates them, and binning mesh quadrature points by psi (@ref gs_flux_cumint) leaves
!! percent-level noise in each node interval.
!---------------------------------------------------------------------------------
SUBROUTINE saw_redistribute(self, gseq, xn, ravgs, jt_base, jt_eq, j_base, dj)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
REAL(r8), INTENT(in) :: xn(:) !< Nodes, 0 at LCFS, 1 at axis [n]
REAL(r8), INTENT(in) :: ravgs(:,:) !< <R>, <1/R>, <1/R^2> on the nodes [n,3]
REAL(r8), INTENT(in) :: jt_base(:) !< Base jtor (@ref eval_dIdA) on the nodes [n]
REAL(r8), INTENT(in) :: jt_eq(:) !< jtor of the equilibrium's current on the nodes [n]
REAL(r8), INTENT(in) :: j_base(:) !< Base jphi on the nodes [n]
REAL(r8), INTENT(out) :: dj(:) !< Reset current [n]
INTEGER(i4) :: n, i, k, n_dips, nk, ierr
INTEGER(i4), ALLOCATABLE :: kk(:)
REAL(r8), PARAMETER :: psi_ax = 1.0e-4_r8
REAL(r8) :: rho_s, rho_m, rho_out, w, wf, dAmin, alast, psi, h
REAL(r8), ALLOCATABLE :: q_eq(:), c1(:), rho(:), qb(:), Ib(:), Ie(:), Ar(:), ap(:), qn(:), djr(:), wr(:), r(:), &
  phi(:), g(:), fs(:,:)
LOGICAL, ALLOCATABLE :: ok(:)
n = SIZE(xn)
ALLOCATE(q_eq(n), c1(n), rho(n), qb(n), Ib(n), Ie(n), Ar(n), ap(n), qn(n), djr(n), wr(n), r(n), phi(n), g(n), &
  fs(n,4), ok(n))
!---q and dA/dpsi_N = oint dl/|grad psi_N| = <1/R> dV/dpsi_N / 2pi on the nodes, axis-first (index 1 =
!   innermost node). Nodes without a closed surface (the LCFS, the axis, a failed surface) take the
!   values of the next node inward (the axis: outward).
CALL gs_fsa_contour(gseq, n, xn, g, fs, ok)
DO i = 1, n
  k = n + 1 - i
  q_eq(i) = 0.0_r8; ap(i) = 0.0_r8
  IF(ok(k))THEN
    psi = gseq%plasma_bounds(1) + xn(k)*(gseq%plasma_bounds(2) - gseq%plasma_bounds(1))
    q_eq(i) = ABS(g(k)*gs_fpol(gseq, psi))
    ap(i) = ABS(fs(k,2)*fs(k,4)*(gseq%plasma_bounds(2) - gseq%plasma_bounds(1)))/(2.0_r8*pi)
  END IF
  c1(i) = ravgs(k,3)/ravgs(k,2)**2
END DO
DO i = 2, n
  IF(.NOT.(q_eq(i) > 0.0_r8))q_eq(i) = q_eq(i-1) ! Also catches NaN
  IF(.NOT.(ap(i) > 0.0_r8))ap(i) = ap(i-1)
END DO
DO i = n-1, 1, -1
  IF(.NOT.(q_eq(i) > 0.0_r8))q_eq(i) = q_eq(i+1)
  IF(.NOT.(ap(i) > 0.0_r8))ap(i) = ap(i+1)
END DO
!---Enclosed area and currents (trapezoid in psi_N; the piece inside the innermost node as a rectangle)
h = 1.0_r8 - xn(n)
Ar(1) = h*ap(1); Ib(1) = h*ap(1)*jt_base(n); Ie(1) = h*ap(1)*jt_eq(n)
DO i = 2, n
  k = n + 1 - i
  h = 0.5_r8*(xn(k+1) - xn(k))
  Ar(i) = Ar(i-1) + h*(ap(i) + ap(i-1))
  Ib(i) = Ib(i-1) + h*(ap(i)*jt_base(k) + ap(i-1)*jt_base(k+1))
  Ie(i) = Ie(i-1) + h*(ap(i)*jt_eq(k) + ap(i-1)*jt_eq(k+1))
END DO
!---Nodes within psi_ax of the axis (contour q/F is unreliable on surfaces much smaller than a
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
IF(.NOT.(phi(n) > 0.0_r8))THEN
  CALL oft_warn('saw_redistribute: no closed surfaces for the q profile, skipping the reset this update')
  dj = 0.0_r8
  CALL saw_scalars_reset(self%boot_profs)
  RETURN
END IF
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
  IF(ABS(Ib(kk(i))) > 0.0_r8 .AND. Ar(kk(i)) > 0.0_r8)r(i) = Ie(kk(i))/Ib(kk(i))
END DO
IF(Ar(kk(1)) <= 0.0_r8 .AND. nk > 1)r(1) = r(2)
DO i = 1, nk
  qb(i) = q_eq(kk(i))*r(i)
END DO
!---The reset differentiates q_base, so smooth its grid-scale roughness first (the shoulder itself is analytic)
CALL quad_smooth(nk, rho(kk(1:nk)), qb(1:nk), saw_nsmooth)
CALL saw_reset_1d(nk, rho(kk(1:nk)), qb(1:nk), Ib(kk(1:nk)), Ar(kk(1:nk)), c1(kk(1:nk)), &
                  self%boot_ops%saw_q_s, self%boot_ops%saw_dq, self%boot_ops%saw_ramp, &
                  self%boot_ops%saw_rule, qn(1:nk), djr(1:nk), wr(1:nk), rho_s, rho_m, rho_out, n_dips, w, ierr)
IF(ierr /= 0)THEN
  CALL oft_warn('saw_redistribute: sawtooth reset failed, skipping the reset this update')
  dj = 0.0_r8
  CALL saw_scalars_reset(self%boot_profs)
  RETURN
END IF
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
DEALLOCATE(q_eq, c1, rho, qb, Ib, Ie, Ar, ap, qn, djr, wr, r, phi, g, fs, ok, kk)
END SUBROUTINE saw_redistribute
!---------------------------------------------------------------------------------
!> Smooth y(x) in place by a least-squares quadratic over +-hw points about each point (Savitzky-Golay
!! on a non-uniform grid; exact for quadratics). If x(1) = 0 the points are mirrored about it (y even in x).
!---------------------------------------------------------------------------------
SUBROUTINE quad_smooth(n, x, y, hw)
INTEGER(i4), INTENT(in) :: n, hw
REAL(r8), INTENT(in) :: x(n)
REAL(r8), INTENT(inout) :: y(n)
INTEGER(i4) :: i, j, k, nm
REAL(r8) :: xe(1-n:n), ye(1-n:n), ys(n), s(0:4), t(0:2), d, det
IF(hw < 1 .OR. n < 3)RETURN
xe(1:n) = x; ye(1:n) = y
nm = 0
IF(x(1) == 0.0_r8)nm = MIN(hw, n-1)
DO j = 1, nm
  xe(1-j) = -x(1+j); ye(1-j) = y(1+j)
END DO
DO i = 1, n
  s = 0.0_r8; t = 0.0_r8
  DO j = MAX(i-hw, 1-nm), MIN(i+hw, n)
    d = xe(j) - x(i)
    DO k = 0, 4
      s(k) = s(k) + d**k
    END DO
    DO k = 0, 2
      t(k) = t(k) + ye(j)*d**k
    END DO
  END DO
  !---Constant term of the normal equations (Cramer)
  det = s(0)*(s(2)*s(4) - s(3)**2) - s(1)*(s(1)*s(4) - s(2)*s(3)) + s(2)*(s(1)*s(3) - s(2)**2)
  ys(i) = y(i)
  IF(ABS(det) > 0.0_r8)ys(i) = (t(0)*(s(2)*s(4) - s(3)**2) - s(1)*(t(1)*s(4) - s(3)*t(2)) &
    + s(2)*(t(1)*s(3) - s(2)*t(2)))/det
END DO
y = ys
END SUBROUTINE quad_smooth
!---------------------------------------------------------------------------------
!> Sawtooth q reset on a 1-D profile (axis first).
!!
!! rule 1 (FUSE saw_crash!): q -> q_s + dq*p(rho/rho_m) from the axis to rho_m (the q = q_s + dq crossing
!! beyond the outermost q < q_s point), p C1 onto q at rho_m, weighted by w = clamp((q_s - min q)/ramp, 0, 1).
!! rule 2 (local): each q < q_s + dq interval holding a q < q_s point is reset about its q minimum,
!! weighted per interval (@ref local_reset); q is unchanged elsewhere, so dI = 0 at interval ends.
!! dI = I*(q/q_new - 1), dj = (d dI/dA)/c1.
!---------------------------------------------------------------------------------
SUBROUTINE saw_reset_1d(n, rho, q, I, A, c1, q_s, dq, ramp, rule, q_new, dj, wr, rho_s, rho_m, rho_out, n_dips, w, ierr)
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
INTEGER(i4), INTENT(out) :: ierr !< 0 on success, -1 if a gradient failed (no reset: q_new = q, dj = 0)
INTEGER(i4) :: k, ia, n_low
INTEGER(i4), ALLOCATABLE :: reg_start(:), reg_end(:)
REAL(r8) :: qm
LOGICAL :: ok
REAL(r8), ALLOCATABLE :: dI(:), gq(:), qt(:)
q_new = q; dj = 0.0_r8; wr = 0.0_r8; rho_s = 0.0_r8; rho_m = 0.0_r8; rho_out = 0.0_r8; w = 0.0_r8; n_dips = 0
ierr = 0
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
CALL gradient_(n, rho, q, gq, ierr, edge_order=1)
IF(ierr /= 0)RETURN
qm = q_s + dq
!---Sawtooth radius (both rules) and rule 1's reset, from the outermost q < q_s region
CALL full_reset(n_low, qt, rho_s, rho_m, ok)
IF(.NOT.ok)THEN
  rho_s = 0.0_r8; rho_m = 0.0_r8
END IF
IF(rule == 1)THEN
  IF(ok)THEN
    n_dips = n_low
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
CALL gradient_(n, A, dI, dj, ierr, edge_order=1)
IF(ierr /= 0)THEN
  q_new = q; wr = 0.0_r8
  RETURN
END IF
dj = dj/c1
DEALLOCATE(reg_start, reg_end, gq, dI, qt)
CONTAINS
!---Full (w = 1) reset from the end of dip region kr; ok = .FALSE. if no mixing radius
SUBROUTINE full_reset(kr, qt_, rs_, rm_, ok_)
INTEGER(i4), INTENT(in) :: kr
REAL(r8), INTENT(out) :: qt_(n), rs_, rm_
LOGICAL, INTENT(out) :: ok_
INTEGER(i4) :: is, jm, kk
REAL(r8) :: slope, m
qt_ = q; rs_ = 0.0_r8; rm_ = 0.0_r8; ok_ = .FALSE.
is = reg_end(kr)
IF(is >= n)RETURN
rs_ = rho(is) + (q_s - q(is))/(q(is+1) - q(is))*(rho(is+1) - rho(is))
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
  qt_(kk) = q_s + dq*pshape(rho(kk)/rm_, m)
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
END SUBROUTINE saw_reset_1d
!------------------------------------------------------------------------------
!> Evaluate terms in augmented tracing ODE for computing Sauter factors
!! (see @ref sauter_fc)
!------------------------------------------------------------------------------
SUBROUTINE sauter_apply(self,cell,f,gop,val)
CLASS(sauter_interp), INTENT(inout) :: self !< Interpolation object
INTEGER(i4), INTENT(in) :: cell !< Cell for interpolation
REAL(r8), INTENT(in) :: f(:) !< Position in cell in logical coord [3]
REAL(r8), INTENT(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
REAL(r8), INTENT(out) :: val(:) !< Reconstructed field at f [8]
INTEGER(i4) :: j(self%lag_rep%nce)
INTEGER(i4) :: jc
REAL(r8) :: rop(3),pt(3),grad(3)
REAL(r8) :: s,c,Bp2,Bt2,mod_B
!---Get dofs
CALL self%lag_rep%ncdofs(cell,j)
!---Reconstruct gradient
grad=0.d0
DO jc=1,self%lag_rep%nce
  CALL oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  grad=grad+self%uvals(j(jc))*rop
END DO
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
END SUBROUTINE sauter_apply
!------------------------------------------------------------------------------
!> Compute factors required for Sauter bootstrap formula
!!
!! Surface averages and extrema come from the contour quadrature
!! (@ref oft_gs_contour::gs_sauter_contour, two passes as the trapped-particle sum
!! needs Bmax). Surfaces it cannot close, or all surfaces when
!! `torflux_qgeom_backend=1` (testing), are traced instead.
!!
!! Each traced surface is traced once. The trapped-particle integral
!! \f$ \langle (1-\sqrt{1-b}(1+b/2))/b^2 \rangle \f$, \f$ b=|B|/B_{max} \f$, needs
!! \f$ B_{max} \f$ of the whole surface, so it is evaluated after the trace by
!! trapezoidal quadrature over the tracer steps. Surfaces are traced in parallel
!! from start points found serially, so results do not depend on the thread count.
!------------------------------------------------------------------------------
SUBROUTINE sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs,qprof,eps,invr2,rgeo)
CLASS(gs_equil), INTENT(inout) :: gseq !< G-S object
INTEGER(i4), INTENT(in) :: nr !< Number of flux sample points
REAL(r8), INTENT(in) :: psi_q(nr) !< Location of flux sample points in normalised psi
REAL(r8), INTENT(out) :: fc(nr) !< Circulating particle fraction \f$ f_c \f$
REAL(r8), INTENT(out) :: r_avgs(nr,3) !< Flux surface averaged radial coords \f$<R>\f$, \f$<1/R>\f$, \f$<a>\f$
REAL(r8), INTENT(out) :: modb_avgs(nr,2) !< Flux surface averaged field \f$<|B|>\f$, \f$<|B|^2>\f$
REAL(r8), OPTIONAL, INTENT(out) :: qprof(nr) !< Safety factor q on each surface (avoids a separate gs_get_qprof call)
REAL(r8), OPTIONAL, INTENT(out) :: eps(nr) !< Inverse aspect ratio (R_max-R_min)/(R_max+R_min) on each surface
REAL(r8), OPTIONAL, INTENT(out) :: invr2(nr) !< \f$<1/R^2>\f$ on each surface
REAL(r8), OPTIONAL, INTENT(out) :: rgeo(nr) !< Geometric major radius (R_max+R_min)/2 on each surface
REAL(r8) :: psi_surf,rmax,x1,x2,raxis,zaxis,h,h2,hf,ftu,ftl
REAL(r8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
REAL(r8), ALLOCATABLE :: pts(:,:),fpol(:),qtmp(:),epstmp(:),i2tmp(:),rgtmp(:)
REAL(r8), POINTER :: ptout(:,:),vout(:,:)
REAL(r8), ALLOCATABLE :: acc(:,:),ext(:,:),bref(:)
LOGICAL, ALLOCATABLE :: done(:)
TYPE(oft_lag_brinterp) :: psi_int
REAL(r8), PARAMETER :: tol=1.d-10
INTEGER(i4) :: j,cell
TYPE(sauter_interp), POINTER :: field
TYPE(gs_factory), POINTER :: device
device=>gseq%device
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
!---Surfaces that fail below keep zeros
fc=0.d0; r_avgs=0.d0; modb_avgs=0.d0
IF(PRESENT(qprof))qprof=0.d0
IF(PRESENT(eps))eps=0.d0
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
ALLOCATE(pts(2,nr),fpol(nr),qtmp(nr),epstmp(nr),i2tmp(nr),rgtmp(nr),done(nr))
!---Surfaces that fail below keep zeros
r_avgs=0.d0; modb_avgs=0.d0; fc=0.d0; qtmp=0.d0; epstmp=0.d0; i2tmp=0.d0; rgtmp=0.d0
DO j=1,nr
  psi_surf=psi_q(j)*(x2-x1) + x1
  fpol(j)=gs_fpol(gseq,psi_surf)
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
  rgtmp(j)       = raxis
END DO
!---Contour averages and extrema
IF(torflux_qgeom_backend/=1)THEN
  ALLOCATE(acc(8,nr),ext(3,nr),bref(nr))
  bref=-1.d0
  CALL gs_sauter_contour(gseq,nr,psi_q,fpol,bref,acc,ext,ok=done)
  bref=ext(1,:)
  CALL gs_sauter_contour(gseq,nr,psi_q,fpol,bref,acc,ext,ok=done)
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
    epstmp(j) = (ext(2,j) + ext(3,j)) / (ext(2,j) - ext(3,j)) ! (R_max - R_min)/(R_max + R_min)
    rgtmp(j) = 0.5d0*(ext(2,j) - ext(3,j))
    i2tmp(j) = acc(7,j)/acc(1,j)
  END DO
  done=done.OR.(psi_q == 1.d0)
  DEALLOCATE(acc,ext,bref)
  IF(oft_debug_print(1).AND.(.NOT.ALL(done)))WRITE(*,'(2A,I0,A)')oft_indent, &
    'sauter_fc: tracing ',COUNT(.NOT.done),' surfaces not closed by the contour pass'
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
CALL set_tracer(1)
!$omp parallel private(j,field,ptout,vout,h,h2,hf,ftu,ftl)
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
DO j=1,nr
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
  CALL tracinginv_fs(device%mesh,pts(:,j),ptout,vout)
  IF(active_tracer%status/=1)THEN
    WRITE(*,*)j,pts(:,j)
    CALL oft_warn("sauter_fc: Trace did not complete")
    CYCLE
  END IF
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
  epstmp(j) = (field%rmax_surf - field%rmin_surf) / (field%rmax_surf + field%rmin_surf)
  rgtmp(j) = 0.5d0*(field%rmax_surf + field%rmin_surf)
  i2tmp(j) = active_tracer%v(8)/active_tracer%v(2)
END DO
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
IF(PRESENT(rgeo))rgeo=rgtmp
DEALLOCATE(pts,fpol,qtmp,epstmp,i2tmp,rgtmp,done)
END SUBROUTINE sauter_fc
!------------------------------------------------------------------------------
!> Trapped-particle integral \f$ \langle (1-\sqrt{1-b}(1+b/2))/b^2 \rangle \f$ by
!! periodic trapezoidal quadrature over the steps of one surface trace
!------------------------------------------------------------------------------
FUNCTION sauter_hf(ns,t,vout,bmax) RESULT(hf)
INTEGER(i4), INTENT(in) :: ns !< Number of steps (the last ends at \f$ 2\pi \f$)
REAL(r8), INTENT(in) :: t(:) !< Angle before each step and at its end [ns+1]
REAL(r8), INTENT(in) :: vout(:,:) !< ODE RHS at each step; (2) FSA weight, (3) weight times \f$ B^2 \f$ [3,ns]
REAL(r8), INTENT(in) :: bmax !< Maximum \f$ |B| \f$ on the surface
REAL(r8) :: hf,g(2),w(2),b,wsum
INTEGER(i4) :: k,i,kk(2)
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
END FUNCTION sauter_hf
!------------------------------------------------------------------------------
!> Apply a smooth edge taper to an array, zeroing it at the plasma edge.
!!
!! @param n          Number of grid points
!! @param psi        Normalised psi grid, OFT convention (0=LCFS, 1=axis)
!! @param arr        Array to taper (modified in-place)
!! @param psi0       psi (OFT convention) where the taper begins; it acts on psi <= psi0
!! @param shape      Taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
!------------------------------------------------------------------------------
SUBROUTINE apply_edge_taper(n, psi, arr, psi0, shape)
INTEGER(i4), INTENT(in)    :: n
REAL(r8),    INTENT(in)    :: psi(n)
REAL(r8),    INTENT(inout) :: arr(n)
REAL(r8),    INTENT(in)    :: psi0
INTEGER(i4), INTENT(in)    :: shape
!---
INTEGER(i4) :: i
REAL(r8)    :: t_taper, w_taper, span
CHARACTER(len=80) :: char_buf
REAL(r8), PARAMETER :: HALF_PI = 1.5707963267948966_r8
span = psi0 ! Taper width: the edge is at psi=0
IF (span < 1.0e-6_r8) RETURN
DO i = 1, n
  IF (psi(i) > psi0) CYCLE
  t_taper = (psi0 - psi(i)) / span
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
!> dy/dx on a strictly monotonic, possibly non-uniform grid, as `numpy.gradient(y, x, edge_order=...)`:
!! three-point central interior, one-sided ends of order `edge_order` (1 or 2, default 2). On a
!! non-finite input or result, a repeated or unsorted abscissa, or too few points: warns, ierr = -1.
!------------------------------------------------------------------------------
SUBROUTINE gradient_(n, x, y, dydx, ierr, edge_order)
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
INTEGER(i4), INTENT(in) :: n !< Number of points
REAL(r8), INTENT(in) :: x(n) !< Abscissae
REAL(r8), INTENT(in) :: y(n) !< Values on x
REAL(r8), INTENT(out) :: dydx(n) !< dy/dx on x
INTEGER(i4), INTENT(out) :: ierr !< 0 on success, -1 on failure (dydx = 0)
INTEGER(i4), OPTIONAL, INTENT(in) :: edge_order !< Order of the end stencils (1 or 2, default 2)
!---
INTEGER(i4) :: i, eo
REAL(r8) :: hs, hd !< Backward and forward interval widths about the stencil centre
CHARACTER(len=512) :: char_buf
eo = 2
IF(PRESENT(edge_order))eo = edge_order
IF(n < eo+1)THEN
  CALL fail('too few input points for the end stencils')
  RETURN
END IF
DO i = 1, n
  IF(.NOT.IEEE_IS_FINITE(x(i)) .OR. .NOT.IEEE_IS_FINITE(y(i)))THEN
    WRITE(char_buf,'(A,I0,A,ES12.5,A,ES12.5)') 'non-finite input at i=', i, &
      ': x=', x(i), ', y=', y(i)
    CALL fail(char_buf)
    RETURN
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
  CALL fail(char_buf)
  RETURN
END DO
!--- Interior: three-point central stencil for the local spacing
DO i = 2, n-1
  hs = x(i) - x(i-1)
  hd = x(i+1) - x(i)
  dydx(i) = -(hd/(hs*(hs + hd)))*y(i-1) &
            + ((hd - hs)/(hs*hd))*y(i) &
            + (hs/(hd*(hs + hd)))*y(i+1)
END DO
!--- Endpoints: one-sided stencils of order eo
IF(eo == 1)THEN
  dydx(1) = (y(2) - y(1))/(x(2) - x(1))
  dydx(n) = (y(n) - y(n-1))/(x(n) - x(n-1))
ELSE
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
END IF
!--- Finite input does not guarantee a finite result: the stencil weights
!   overflow for large y on a very fine grid.
DO i = 1, n
  IF(.NOT.IEEE_IS_FINITE(dydx(i)))THEN
    WRITE(char_buf,'(A,I0,A,ES12.5)') 'non-finite derivative at i=', i, &
      ', x=', x(i)
    CALL fail(char_buf)
    RETURN
  END IF
END DO
ierr = 0
CONTAINS
SUBROUTINE fail(msg)
CHARACTER(LEN=*), INTENT(in) :: msg
CALL oft_warn('gradient_: '//TRIM(msg))
dydx = 0.0_r8
ierr = -1
END SUBROUTINE fail
END SUBROUTINE gradient_
!------------------------------------------------------------------------------
!> Bootstrap current (Redl 2021) on the nodes psi_N (bootstrap.py calculate_profiles_and_bootstrap)
!!
!! NRL Coulomb logs, Koh ion collisionality, Redl jboot1 (no sign(q) factor: <j.B> is odd in F), L34=L31.
!!
!! @param gseq    Equilibrium object (must have Te, Ti, ne, ni, Zeff set)
!! @param n_psi   Number of flux surface samples
!! @param psi_N   Normalised psi, OFT convention (0=LCFS, 1=axis), any spacing, psi_N(1)>0
!! @param jphi_bs Output: bootstrap current density as TokaMaker jphi = <j_phi> [A/m^2] on psi_N grid
!! @param jphi_pd Output: its pressure-driven part [A/m^2]
!! @param ierr    Output: 0 on success, -1 if a gradient or the boundary extrapolation failed
!------------------------------------------------------------------------------
SUBROUTINE calculate_bootstrap(self, gseq, n_psi, psi_N, jphi_bs, jphi_pd, ierr, ravgs)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)  !< Normalised psi (OFT convention, 0=LCFS), psi_N(1)>0
REAL(r8), INTENT(out) :: jphi_bs(0:n_psi) !< jphi_bs(0) the LCFS value (extrapolated)
REAL(r8), INTENT(out) :: jphi_pd(0:n_psi) !< Pressure-driven part of jphi_bs
INTEGER(i4), INTENT(out) :: ierr !< 0 on success, -1 on failure (warned)
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
REAL(r8) :: R_geo(n_psi)    !< (R_max+R_min)/2 [m] per flux surface
REAL(r8) :: B_avg(n_psi)    !< <B> [T] per flux surface
REAL(r8) :: r_avgs_saut(n_psi,3)    !< Sauter FSA: <R>, <1/R>, <a>
REAL(r8) :: invr2(n_psi)            !< Sauter FSA: <1/R^2>
REAL(r8) :: modb_avgs_saut(n_psi,2) !< Sauter FSA: <|B|>, <|B|^2>
REAL(r8) :: dn_e_dpsi(n_psi), dT_e_dpsi(n_psi)
REAL(r8) :: dn_i_dpsi(n_psi), dT_i_dpsi(n_psi)
REAL(r8) :: ln_le(n_psi), ln_lii(n_psi), Z_lnLam(n_psi)
REAL(r8) :: Zavg(n_psi), Zion(n_psi)
REAL(r8) :: nu_i_star(n_psi), nu_e_star(n_psi)
REAL(r8) :: jdotB_bs(n_psi)
REAL(r8) :: jdotB(0:n_psi) !< <j_BS.B> with its LCFS/axis values
REAL(r8) :: psi_range, Zdom, pscale, pprime
REAL(r8), PARAMETER :: EC = 1.602176634e-19_r8
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
  f(i) = gs_fpol(gseq, psi_abs(i))
END DO
! Get flux-surface geometry in one pass: fc, eps, q and R_geo (the major radius in nu*, as eps = r/R_geo)
CALL sauter_fc(gseq, n_psi, psi_N, fc, r_avgs_saut, modb_avgs_saut, qprof=qvals, eps=eps, invr2=invr2, rgeo=R_geo)
IF(PRESENT(ravgs))ravgs = RESHAPE([r_avgs_saut(:,1:2), invr2], [n_psi,3])
B_avg = modb_avgs_saut(:,1)
! ===================================================================
ft = 1.0_r8 - fc
! Pressures [Pa]
pe = EC * ne * Te
pi_arr = EC * ni * Ti
! Gradients d/dpsi [Wb^-1]
CALL gradient_(n_psi, psi_abs, Te, dT_e_dpsi, ierr)
IF(ierr == 0)CALL gradient_(n_psi, psi_abs, Ti, dT_i_dpsi, ierr)
IF(ierr == 0)CALL gradient_(n_psi, psi_abs, ne, dn_e_dpsi, ierr)
IF(ierr == 0)CALL gradient_(n_psi, psi_abs, ni, dn_i_dpsi, ierr)
IF(ierr /= 0)RETURN
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
nu_i_star = 4.90e-18_r8 * ABS(qvals) * R_geo * ni &
          * Zion**4 * ln_lii / (Ti**2 * eps**1.5_r8)
! Electron collisionality
nu_e_star = 6.921e-18_r8 * ABS(qvals) * R_geo * ne &
          * Zeff * ln_le / (Te**2 * eps**1.5_r8)
! Compute bootstrap current via Redl 2021
CALL redl_bootstrap(n_psi, Te, Ti, ne, ni, pe, pi_arr, Zeff, eps, ft, f, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    nu_e_star, nu_i_star, jdotB_bs)
! Convert <j_BS.B> to TokaMaker's jphi = <j_phi> = <R>P' + <1/R>FF'/mu0 (exact, see
! doc_tokamaker_current_conventions.md eq. A7): the field-aligned part F<1/R>/<B^2> * <j_BS.B>,
! plus the pressure-driven (diamagnetic + Pfirsch-Schlueter) part P'(<R> - F^2<1/R>/<B^2>),
! which is assigned to the bootstrap (as IMAS includes_bootstrap=true); boot_profs report it apart.
! P' as in jphi_bs_update.
pscale = gseq%pax_target/gseq%P%f(gseq%plasma_bounds(2))
jphi_bs(0) = 0.0_r8 ! LCFS value set by the extrapolation below
jphi_pd(0) = 0.0_r8
DO i = 1, n_psi
  IF(ABS(f(i)) > 0.0_r8 .AND. modb_avgs_saut(i,2) > 0.0_r8)THEN
    pprime = gseq%P%fp(psi_abs(i))*pscale/mu0
    jphi_pd(i) = pprime*(r_avgs_saut(i,1) - f(i)**2*r_avgs_saut(i,2)/modb_avgs_saut(i,2))
    jphi_bs(i) = jdotB_bs(i)*f(i)*r_avgs_saut(i,2)/modb_avgs_saut(i,2) + jphi_pd(i)
  ELSE
    jphi_pd(i) = 0.0_r8
    jphi_bs(i) = 0.0_r8
  END IF
END DO
jdotB = [0.0_r8, jdotB_bs]
CALL finish(jphi_bs)
IF(ierr == 0)CALL finish(jphi_pd)
IF(ierr == 0)CALL finish(jdotB)
IF(ierr /= 0)RETURN
self%boot_profs%jdotB_bs_raw = jdotB
self%boot_profs%jphi_bs_raw = jphi_bs - jphi_pd
self%boot_profs%jphi_pd_raw = jphi_pd
IF(self%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A)') '  [calculate_bootstrap] geometry & collisionality sample (i=1,mid,n):'
  WRITE(*,'(A,3ES12.4)') '    <R>      : ', r_avgs_saut(1,1), r_avgs_saut(n_psi/2,1), r_avgs_saut(n_psi,1)
  WRITE(*,'(A,3ES12.4)') '    <1/R>    : ', r_avgs_saut(1,2), r_avgs_saut(n_psi/2,2), r_avgs_saut(n_psi,2)
  WRITE(*,'(A,3ES12.4)') '    <B>      : ', B_avg(1), B_avg(n_psi/2), B_avg(n_psi)
  WRITE(*,'(A,3ES12.4)') '    <B^2>    : ', modb_avgs_saut(1,2), modb_avgs_saut(n_psi/2,2), modb_avgs_saut(n_psi,2)
  WRITE(*,'(A,3ES12.4)') '    eps      : ', eps(1), eps(n_psi/2), eps(n_psi)
  WRITE(*,'(A,3ES12.4)') '    q        : ', qvals(1), qvals(n_psi/2), qvals(n_psi)
  WRITE(*,'(A,3ES12.4)') '    nu_e_star: ', nu_e_star(1), nu_e_star(n_psi/2), nu_e_star(n_psi)
  WRITE(*,'(A,3ES12.4)') '    jphi_bs  : ', jphi_bs(1), jphi_bs(n_psi/2), jphi_bs(n_psi)
END IF
CONTAINS
!> Zero non-finite values (where F -> 0) and extrapolate to the LCFS/axis, where q is undefined
SUBROUTINE finish(a)
REAL(r8), INTENT(inout) :: a(0:n_psi)
WHERE(.NOT.(ABS(a) < 1.0e99_r8)) a = 0.0_r8
CALL extrap_jbs_boundaries(n_psi, psi_N, a, ierr)
END SUBROUTINE finish
END SUBROUTINE calculate_bootstrap
!------------------------------------------------------------------------------
!> Edge spike of a bootstrap profile (OFT convention), as isolate_edge_jbs in
!! bootstrap.py calculate_profiles_and_bootstrap
!------------------------------------------------------------------------------
SUBROUTINE edge_spike(n_psi, psi_N, jphi_bs, spike)
INTEGER(i4), INTENT(in) :: n_psi !< Number of nodes (without the LCFS)
REAL(r8), INTENT(in) :: psi_N(1:n_psi) !< Normalised psi (OFT convention, 0=LCFS)
REAL(r8), INTENT(in) :: jphi_bs(0:n_psi) !< Bootstrap profile, jphi_bs(0) at the LCFS
REAL(r8), INTENT(out) :: spike(0:n_psi) !< Edge spike
REAL(r8) :: psi_N_std(0:n_psi), jphi_bs_std(0:n_psi), mask_std(0:n_psi)
! analyze_bootstrap_edge_spike uses the standard convention (0=axis, 1=LCFS): flip in and out
psi_N_std(0:n_psi-1) = 1.0_r8 - psi_N(n_psi:1:-1)
psi_N_std(n_psi) = 1.0_r8
jphi_bs_std = jphi_bs(n_psi:0:-1)
CALL analyze_bootstrap_edge_spike(n_psi+1, psi_N_std, jphi_bs_std, mask_std)
spike = mask_std(n_psi:0:-1)
END SUBROUTINE edge_spike
!------------------------------------------------------------------------------
!> Linearly extrapolate a bootstrap profile to the LCFS (and the axis, if a node) where q is
!! undefined, using the second-order gradient at the nearest well-defined point.
!------------------------------------------------------------------------------
SUBROUTINE extrap_jbs_boundaries(n_psi, psi_N, prof, ierr)
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)
REAL(r8), INTENT(inout) :: prof(0:n_psi)
INTEGER(i4), INTENT(out) :: ierr !< 0 on success, -1 if too few points or the gradient failed
INTEGER(i4) :: n_lo, n_hi
REAL(r8) :: dprof_dpsi(n_psi)
! psi_N(0) = 0.0 at the LCFS is implicit; an axis node is excluded from both stencils.
n_lo = 1
n_hi = MERGE(n_psi-1, n_psi, psi_N(n_psi) == 1.0_r8)
ierr = -1
IF(n_hi-n_lo+1 < 3)THEN
  CALL oft_warn('extrap_jbs_boundaries: too few points with well-defined q to extrapolate the '// &
    'bootstrap current to the LCFS and axis')
  RETURN
END IF
CALL gradient_(n_hi-n_lo+1, psi_N(n_lo:n_hi), prof(n_lo:n_hi), &
  dprof_dpsi(n_lo:n_hi), ierr)
IF(ierr /= 0)RETURN
prof(0) = prof(n_lo) + dprof_dpsi(n_lo) * (0.0_r8 - psi_N(n_lo))
IF(psi_N(n_psi) == 1.0_r8) &
  prof(n_psi) = prof(n_hi) + dprof_dpsi(n_hi) * (psi_N(n_psi) - psi_N(n_hi))
END SUBROUTINE extrap_jbs_boundaries
!------------------------------------------------------------------------------
!> Analyse and isolate the edge bootstrap-current spike from a toroidal
!! bootstrap current density profile.
!!
!! Translated from Python bootstrap.py: analyze_bootstrap_edge_spike
!!
!! The routine performs the following steps:
!!   1. Locates the highest local peak (whose prominence exceeds 5%)
!!      in psi_N in [0.85,1.0], else the rightmost peak for psi_N > 0.7.
!!      With no edge region or no qualifying peak the un-isolated profile
!!      is returned (with a warning) rather than zeros.
!!   2. Finds the minimum of jphi_bs between psi_N = 0.5 and the peak (the
!!      flat-core stitching level, lmin_jbs).
!!   3. Builds masked_spike: flat at lmin_jbs for psi_N below the minimum
!!      index, actual jphi_bs from that index to the separatrix.
!!   4. Smooths the flat/spike junction with a cubic Hermite spline over a
!!      window of width min(0.5 * dist_to_peak, 0.2) centred on the minimum.
!!
!! @param n                  Number of flux surface samples
!! @param psi_N              Normalised poloidal flux grid [0, 1]
!! @param jphi_bs           Bootstrap current density j_phi profile [A/m^2]
!! @param masked_spike       Output: isolated edge spike spliced onto flat core [A/m^2]
!------------------------------------------------------------------------------
SUBROUTINE analyze_bootstrap_edge_spike(n, psi_N, jphi_bs, masked_spike)
INTEGER(i4), INTENT(in)  :: n
REAL(r8),    INTENT(in)  :: psi_N(n)              !< Normalised poloidal flux [0=axis, 1=LCFS/plasma edge]
REAL(r8),    INTENT(in)  :: jphi_bs(n)           !< Bootstrap current density [A/m^2]
REAL(r8),    INTENT(out) :: masked_spike(n)        !< Isolated edge spike spliced onto flat core [A/m^2]
!---
INTEGER(i4) :: i, k, peak_idx, far_idx, lmin_idx
INTEGER(i4) :: idx_start, idx_end
LOGICAL     :: has_edge
REAL(r8)    :: far_height, dist_tmp
REAL(r8)    :: j_edge_max, min_prominence, prom, left_base, right_base
REAL(r8)    :: peak_psi, lmin_jbs
REAL(r8)    :: jbs_min_loc, dist_to_peak, blend_width_val, x_start, x_end
REAL(r8)    :: dist_start, dist_end
REAL(r8)    :: y_start, dy_start, y_end, dy_end
REAL(r8)    :: dx_window, t, h00, h10, h01, h11, y_patch
! =====================================================================
! 1. Find the highest local peak (whose prominence exceeds 5%)
!    in psi_N in [0.85,1.0], else the rightmost peak for psi_N > 0.7.
!    Fall back to the input jphi_bs if no edge region or no qualifying peak.
! =====================================================================
has_edge   = .FALSE.
j_edge_max = 0.0_r8
DO i = 1, n
  IF (psi_N(i) >= 0.7_r8) THEN
    IF (has_edge) THEN
      j_edge_max = MAX(j_edge_max, jphi_bs(i))
    ELSE
      j_edge_max = jphi_bs(i)
      has_edge   = .TRUE.
    END IF
  END IF
END DO
IF (.NOT. has_edge) THEN
  CALL oft_warn('[analyze_bootstrap_edge_spike] no edge '// &
    'region (psi_N >= 0.7) in profile; using '// &
    'bootstrap profile without edge isolation')
  masked_spike = jphi_bs
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
  IF (jphi_bs(i) <= 0.0_r8)    CYCLE
  IF (jphi_bs(i) <= jphi_bs(i-1)) CYCLE
  IF (jphi_bs(i) <= jphi_bs(i+1)) CYCLE
  ! Prominence = peak minus the higher of the two bases (nearest minima
  ! before a taller point on each side, within the psi_N >= 0.7 slice).
  left_base = jphi_bs(i)
  DO k = i-1, 1, -1
    IF (psi_N(k) < 0.7_r8)  EXIT
    IF (jphi_bs(k) > jphi_bs(i))  EXIT
    left_base = MIN(left_base, jphi_bs(k))
  END DO
  right_base = jphi_bs(i)
  DO k = i+1, n
    IF (jphi_bs(k) > jphi_bs(i))  EXIT
    right_base = MIN(right_base, jphi_bs(k))
  END DO
  prom = jphi_bs(i) - MAX(left_base, right_base)
  IF (prom < min_prominence) CYCLE
  peak_idx = i
  IF (psi_N(i) > 0.85_r8 .AND. jphi_bs(i) > far_height) THEN
    far_idx    = i
    far_height = jphi_bs(i)
  END IF
END DO
! Prefer the tallest far-edge peak, else the rightmost.
IF (far_idx >= 1) peak_idx = far_idx
IF (peak_idx < 1) THEN
  CALL oft_warn('[analyze_bootstrap_edge_spike] no clear '// &
    'edge peak found; using bootstrap profile '// &
    'without edge isolation')
  masked_spike = jphi_bs
  RETURN
END IF
peak_psi    = psi_N(peak_idx)
! =====================================================================
! 2. Find the minimum of jphi_bs between psi_N = 0.5 and the peak.
!    This becomes the flat-core stitching level (lmin_jbs).
!    Mirrors: lmin_jbs = min(j_bootstrap[(psi_N > 0.5) & (psi_N < peak_psi)])
! =====================================================================
lmin_jbs = 1.0e30_r8
lmin_idx  = peak_idx
DO i = 1, n
  IF (psi_N(i) > 0.5_r8 .AND. psi_N(i) < peak_psi) THEN
    IF (jphi_bs(i) < lmin_jbs) THEN
      lmin_jbs = jphi_bs(i)
      lmin_idx  = i
    END IF
  END IF
END DO
IF (lmin_idx == peak_idx) THEN
  CALL oft_warn('[analyze_bootstrap_edge_spike] no node between psi_N = 0.5 and the '// &
    'edge peak; using bootstrap profile without edge isolation')
  masked_spike = jphi_bs
  RETURN
END IF
! =====================================================================
! 3. Build masked_spike:
!    - Flat at lmin_jbs for psi_N < psi_N(lmin_idx)
!    - Actual jphi_bs   for psi_N >= psi_N(lmin_idx)
!    Mirrors: fit_mask = (psi_N >= masked_psi_N[lmin_arg])
! =====================================================================
DO i = 1, n
  IF (psi_N(i) >= psi_N(lmin_idx)) THEN
    masked_spike(i) = jphi_bs(i)
  ELSE
    masked_spike(i) = lmin_jbs
  END IF
END DO
! =====================================================================
! 4. Smooth the flat/spike junction with a cubic Hermite spline.
!    Blend window  = min(0.5 * dist_to_peak, 0.2), centred on lmin_loc.
!    Mirrors the Hermite-spline patch in the Python.
! =====================================================================
jbs_min_loc     = psi_N(lmin_idx)
dist_to_peak    = peak_psi - jbs_min_loc
blend_width_val = MIN(0.5_r8 * dist_to_peak, 0.2_r8)
x_start         = jbs_min_loc - 0.5_r8 * blend_width_val
x_end           = jbs_min_loc + 0.5_r8 * blend_width_val
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
END SUBROUTINE analyze_bootstrap_edge_spike
!------------------------------------------------------------------------------
!> Redl 2021 bootstrap current formula.
!!
!! Translates Python bootstrap.py:redl_bootstrap with fixed settings:
!!   - use_legacy_L34 = .FALSE.  (L34 = L31)
!!   - no sign(q) factor (<j.B> is odd in F through I_psi)
!!   - formula_form   = 'jboot1'
!!   - nu_e_star and nu_i_star are passed in pre-computed (no internal fallback)
!!
!! Reference: Redl et al., Phys. Plasmas 28, 022502 (2021)
!!
!! @param n           Number of flux surfaces
!! @param Te          Electron temperature [eV]
!! @param Ti          Ion temperature [eV]
!! @param ne          Electron density [m^-3]
!! @param ni          Ion density [m^-3]
!! @param pe          Electron pressure [Pa]
!! @param pi          Ion pressure [Pa]
!! @param Zeff        Effective charge
!! @param eps         Inverse aspect ratio
!! @param fT          Trapped particle fraction
!! @param I_psi       Toroidal current function F = R*Bt [T*m]
!! @param dT_e_dpsi   d(Te)/d(psi) [eV/Wb]
!! @param dT_i_dpsi   d(Ti)/d(psi) [eV/Wb]
!! @param dn_e_dpsi   d(ne)/d(psi) [m^-3/Wb]
!! @param dn_i_dpsi   d(ni)/d(psi) [m^-3/Wb]
!! @param nu_e_star   Electron collisionality (pre-computed)
!! @param nu_i_star   Ion collisionality (pre-computed)
!! @param avg_j_bootstrap_times_B Output: <j_{bs,parallel}*B> [same units as -I_psi * pe * d/dpsi]
!------------------------------------------------------------------------------
SUBROUTINE redl_bootstrap(n, Te, Ti, ne, ni, pe, pi, Zeff, eps, fT, I_psi, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    nu_e_star, nu_i_star, avg_j_bootstrap_times_B)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: Te(n), Ti(n), ne(n), ni(n)
REAL(r8), INTENT(in) :: pe(n), pi(n), Zeff(n)
REAL(r8), INTENT(in) :: eps(n), fT(n), I_psi(n)
REAL(r8), INTENT(in) :: dT_e_dpsi(n), dT_i_dpsi(n)
REAL(r8), INTENT(in) :: dn_e_dpsi(n), dn_i_dpsi(n)
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
         / (Zeff*(0.77_r8 + 0.63_r8*(1.0_r8 + MAX(Zeff - 1.0_r8, 0.0_r8)**1.1_r8))) &
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
! Assemble: jboot1 form
! =====================================================================
dp_dpsi = (ne*dT_e_dpsi + Te*dn_e_dpsi &
         + ni*dT_i_dpsi + Ti*dn_i_dpsi) * EC
bra1 = L31 * dp_dpsi / pe
bra2 = L32 * dT_e_dpsi / Te
bra3 = L34 * alpha * (1.0_r8 - R_pe) / R_pe * dT_i_dpsi / Ti
avg_j_bootstrap_times_B = -I_psi * pe * (bra1 + bra2 + bra3)
END SUBROUTINE redl_bootstrap
END MODULE grad_shaf_bootstrap