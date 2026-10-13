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
 gs_itor_nl, gs_fpol
USE oft_lag_basis, ONLY: oft_blag_geval
USE oft_mesh_type, ONLY: bmesh_findcell
USE oft_blag_operators, ONLY: oft_lag_brinterp
USE tracing_2d, ONLY: set_tracer, active_tracer, tracinginv_fs
USE grad_shaf_prof_phys, ONLY: eval_dIdA, build_Ravg_spline, gs_flux_int, &
  jphi_update, jphi_copy, jphi_flux_func, jphi_psi_nodes
USE spline_mod
USE oft_gs_profiles, ONLY: torflux_set_bounds
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
  LOGICAL :: parameterize_jbs = .FALSE. !< Use parametrised skew-normal fit for spike? Overrides `isolate_edge_jbs` if true.
  REAL(r8) :: scale_jbs = 1.0_r8 !< Scaling factor applied to the bootstrap current (default 1)
  REAL(r8) :: djbs_tol = 1.0e-4_r8 !< Rel. change in jphi_bs below which it is frozen (not recomputed)
  LOGICAL :: diagnose_bs = .FALSE. !< Print alpha/Ip scalars, bootstrap stats, and full profile tables each NL iteration
  LOGICAL :: taper_edge_jbs = .FALSE. !< Taper all current components (jphi_ind, _bs, _pd, _fixed) to zero at the LCFS
  REAL(r8) :: taper_edge_psi0 = 0.999_r8 !< psi_N (standard: 0=axis, 1=LCFS) where edge taper begins
  INTEGER(i4) :: taper_edge_shape = 2 !< Edge taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
END TYPE boot_ops
!------------------------------------------------------------------------------
!> Current profiles of the last @ref jphi_bs_update (diagnostics; the update itself uses its own
!! `jphi_*` arrays). `_raw`: from the last bootstrap evaluation, before scaling, isolation and taper;
!! `_final`: as added into `jphi_total`. TokaMaker jphi [A/m²] on `psi_n` (OFT convention, 0=LCFS,
!! 1=axis; LCFS first); `jdotb_bs_raw` in [T A/m²].
!------------------------------------------------------------------------------
TYPE :: boot_profs
  REAL(r8), ALLOCATABLE :: psi_n(:) !< Normalised psi of the profile nodes, LCFS first
  REAL(r8), ALLOCATABLE :: jdotb_bs_raw(:) !< Redl PoP 2021 <j_BS.B>
  REAL(r8), ALLOCATABLE :: jphi_bs_raw(:) !< Bootstrap jphi, <j_BS.B> F<1/R>/<B^2> + jphi_pd_raw
  REAL(r8), ALLOCATABLE :: jphi_pd_raw(:) !< Pressure-driven (diamagnetic + Pfirsch-Schlueter) part, P'(<R> - F^2<1/R>/<B^2>)
  REAL(r8), ALLOCATABLE :: jphi_total(:) !< Total = jphi_ind_final + jphi_bs_final + jphi_fixed_final
  REAL(r8), ALLOCATABLE :: jphi_bs_final(:) !< Bootstrap jphi, scaled, optionally isolated/parametrised and tapered
  REAL(r8), ALLOCATABLE :: jphi_pd_final(:) !< jphi_pd_raw scaled and optionally tapered (not isolated/parametrised)
  REAL(r8), ALLOCATABLE :: jphi_ind_final(:) !< Inductive jphi: input jphi times alpha, optionally tapered
  REAL(r8), ALLOCATABLE :: jphi_fixed_final(:) !< Fixed (non-rescaled) jphi, optionally tapered
END TYPE boot_profs
!> Names of the @ref boot_profs arrays, in declaration order (HDF5 datasets are these in upper case)
CHARACTER(LEN=16), PARAMETER :: boot_prof_names(9) = [CHARACTER(LEN=16) :: 'psi_n', 'jdotb_bs_raw', &
  'jphi_bs_raw', 'jphi_pd_raw', 'jphi_total', 'jphi_bs_final', 'jphi_pd_final', 'jphi_ind_final', &
  'jphi_fixed_final']
REAL(r8), PARAMETER :: djbs_stol = 1.0e-3_r8 !< Stagnation tolerance: a stalled change above this resets the no-improve counter
REAL(r8), PARAMETER :: dalpha_warn_tol = 1.0e-6_r8 !< Relative (|dalpha|/|alpha|): no alpha stall warning below this
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
  LOGICAL :: stage_1 = .FALSE. !< Stage 1 pass (finding Bmax)
  REAL(r8) :: bmax = -1.d0       !< Maximum |B| on surface (set in stage 1)
  REAL(r8) :: rmax_surf = -1.d30 !< Maximum R on surface (set in stage 1)
  REAL(r8) :: rmin_surf =  1.d30 !< Minimum R on surface (set in stage 1)
  REAL(r8) :: mag_axis(2) = 0.d0 !< Magnetic axis (R, Z)
CONTAINS
  !> Evaluate ODE RHS for Sauter integration
  PROCEDURE :: interp => sauter_apply
END TYPE sauter_interp
!------------------------------------------------------------------------------
!> Context type for MINPACK-based fitting of the edge bootstrap profile.
!! Module-level variable `active_edge_jbs` is populated by curve_fit_edge_jbs
!! before invoking lmdif so that edge_jbs_residual can access the target data.
!------------------------------------------------------------------------------
TYPE :: edge_jbs_fit_ctx
  INTEGER(i4) :: n_fit = 0           !< Number of data points in the fit
  REAL(r8), ALLOCATABLE :: psi_fit(:) !< psi_N values at the fitting points
  REAL(r8), ALLOCATABLE :: j_fit(:)   !< Target jphi_bs values [A/m^2]
  REAL(r8) :: tail_alpha = 1.5_r8    !< Fixed right-side fall-off (Python default)
  REAL(r8) :: lb(7) = 0.0_r8        !< Lower bounds for the 7 parameters
  REAL(r8) :: ub(7) = 1.0_r8        !< Upper bounds for the 7 parameters
END TYPE edge_jbs_fit_ctx
TYPE(edge_jbs_fit_ctx) :: active_edge_jbs !< Module-level context for lmdif callback
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
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%parameterize_jbs),filename,path//'/BOOT_OPS/PARAMETERIZE_JBS')
  CALL hdf5_write(self%boot_ops%scale_jbs,filename,path//'/BOOT_OPS/SCALE_JBS')
  CALL hdf5_write(self%boot_ops%djbs_tol,filename,path//'/BOOT_OPS/DJBS_TOL')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%diagnose_bs),filename,path//'/BOOT_OPS/DIAGNOSE_BS')
  CALL hdf5_write(MERGE(1_i4, 0_i4, self%boot_ops%taper_edge_jbs),filename,path//'/BOOT_OPS/TAPER_EDGE_JBS')
  CALL hdf5_write(self%boot_ops%taper_edge_psi0,filename,path//'/BOOT_OPS/TAPER_EDGE_PSI0')
  CALL hdf5_write(self%boot_ops%taper_edge_shape,filename,path//'/BOOT_OPS/TAPER_EDGE_SHAPE')
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
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/PARAMETERIZE_JBS',success=success)
  IF(success) self%boot_ops%parameterize_jbs = (int_tmp/=0)
  CALL hdf5_read(self%boot_ops%scale_jbs,filename,path//'/BOOT_OPS/SCALE_JBS',success=success)
  CALL hdf5_read(self%boot_ops%djbs_tol,filename,path//'/BOOT_OPS/DJBS_TOL',success=success)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/DIAGNOSE_BS',success=success)
  IF(success) self%boot_ops%diagnose_bs = (int_tmp/=0)
  CALL hdf5_read(int_tmp,filename,path//'/BOOT_OPS/TAPER_EDGE_JBS',success=success)
  IF(success) self%boot_ops%taper_edge_jbs = (int_tmp/=0)
  CALL hdf5_read(self%boot_ops%taper_edge_psi0,filename,path//'/BOOT_OPS/TAPER_EDGE_PSI0',success=success)
  CALL hdf5_read(self%boot_ops%taper_edge_shape,filename,path//'/BOOT_OPS/TAPER_EDGE_SHAPE',success=success)
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
IF(.NOT.ASSOCIATED(gseq%I))RETURN
b1=gseq%plasma_bounds(1); b2=gseq%plasma_bounds(2)
IF((b1<-1.d98).OR.(b2==b1))RETURN ! No plasma yet
CALL jphi_psi_nodes(self,xpsi) ! Knots in psi (mapped from toroidal flux if on it)
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
END SUBROUTINE reset_solve_state
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
  CASE('jdotb_bs_raw'); CALL access(bp%jdotb_bs_raw)
  CASE('jphi_bs_raw'); CALL access(bp%jphi_bs_raw)
  CASE('jphi_pd_raw'); CALL access(bp%jphi_pd_raw)
  CASE('jphi_total'); CALL access(bp%jphi_total)
  CASE('jphi_bs_final'); CALL access(bp%jphi_bs_final)
  CASE('jphi_pd_final'); CALL access(bp%jphi_pd_final)
  CASE('jphi_ind_final'); CALL access(bp%jphi_ind_final)
  CASE('jphi_fixed_final'); CALL access(bp%jphi_fixed_final)
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
!! (LCFS and the profile nodes xpsi):
!!   1. <R>/<1/R>/<1/R^2> spline; pressure scale.
!!   2. jphi_fixed and the bootstrap current jphi_bs (with its pressure-driven part jphi_pd), scaled
!!      and optionally isolated, or reused from the last evaluation when frozen.
!!   3. Edge taper of jphi_bs, jphi_pd, jphi_ind and jphi_fixed.
!!   4. alpha from the exact I_p (gs_flux_int of eval_dIdA, eq. A9c), which is affine in
!!      alpha: alpha = (Ip_target - Ip_lo)/(Ip_hi - Ip_lo).
!!   5. jphi_total = alpha*jphi_ind + jphi_bs + jphi_fixed; boot_profs; F*F' knots.
!!   6. Diagnostics (if diagnose_bs is set).
!---------------------------------------------------------------------------------
SUBROUTINE jphi_bs_update(self, gseq)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: pscale, pprime
REAL(r8), ALLOCATABLE :: xpsi(:) !< Node locations in normalized poloidal flux
TYPE(spline_type) :: R_spline
REAL(r8), ALLOCATABLE :: psi_n(:) !< Nodes: LCFS (0) and xpsi
REAL(r8), ALLOCATABLE :: jtor(:)  !< dI/dA for the exact I_p measure (eval_dIdA)
REAL(r8), ALLOCATABLE :: jphi_bs(:), jphi_pd(:), jphi_ind(:), jphi_fixed(:), jphi_total(:)
REAL(r8), ALLOCATABLE :: jphi_bs_spike(:), jphi_bs_spike_masked(:) !< Isolated / masked edge spike [A/m²]
LOGICAL :: isolate
REAL(r8) :: alpha, ip_target, ip_ind, ip_result_lo, ip_result_hi, dalpha
REAL(r8) :: djbs !< Relative RMS change of jphi_bs since the last evaluation
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
!--- 1. <R>/<1/R>/<1/R^2> spline (I_p measure, F*F' map), kept until the F*F' loop; pressure scale.
ALLOCATE(psi_n(0:self%npsi), jtor(0:self%npsi))
psi_n = [0.0_r8, xpsi]
CALL build_Ravg_spline(gseq, self%ngeom, R_spline)
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
IF(ASSOCIATED(gseq%P_ani)) &
  CALL oft_abort('Jphi profiles do not support anisotropic pressure', &
                 'jphi_bs_update',__FILE__)
pscale = gseq%P%f(gseq%plasma_bounds(2))
pscale = gseq%pax_target / pscale
!--- 2. Fixed and bootstrap currents
ALLOCATE(jphi_fixed(0:self%npsi), jphi_bs(0:self%npsi), jphi_pd(0:self%npsi))
CALL node_vals(gseq%jphi_fixed, jphi_fixed)
IF(self%freeze_jbs .AND. ALLOCATED(self%jphi_bs_last)) THEN
  jphi_bs = self%jphi_bs_last
  jphi_pd = self%jphi_pd_last
  djbs = 0.0_r8
ELSE
  CALL calculate_bootstrap(self, gseq, self%npsi, xpsi, jphi_bs, jphi_pd)
  isolate = self%boot_ops%isolate_edge_jbs .OR. self%boot_ops%parameterize_jbs
  IF(isolate)THEN
    ALLOCATE(jphi_bs_spike(0:self%npsi), jphi_bs_spike_masked(0:self%npsi))
    CALL edge_spike(self%npsi, xpsi, jphi_bs, self%boot_ops%parameterize_jbs, self%boot_ops%diagnose_bs, &
                    jphi_bs_spike, jphi_bs_spike_masked)
  END IF
  IF(self%boot_ops%diagnose_bs)THEN
    IF(isolate)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         jphi_bs(bulk)[A/m2]  spike[A/m2]   spike_masked[A/m2]  jphi[A/m2]  jphi_fixed[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,6ES15.5)') '  ', i, xpsi(i), jphi_bs(i), jphi_bs_spike(i), jphi_bs_spike_masked(i), self%jphi(i), jphi_fixed(i)/mu0
      END DO
    ELSE
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         jphi_bs[A/m2]      jphi[A/m2]  jphi_fixed[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,4ES15.5)') '  ', i, xpsi(i), jphi_bs(i), self%jphi(i), jphi_fixed(i)/mu0
      END DO
    END IF
  END IF
  IF(isolate)THEN
    jphi_bs = jphi_bs_spike
    DEALLOCATE(jphi_bs_spike, jphi_bs_spike_masked)
  END IF
  jphi_bs = jphi_bs * self%boot_ops%scale_jbs
  jphi_pd = jphi_pd * self%boot_ops%scale_jbs
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
!--- 3. Edge taper (the cache above keeps the untapered jphi_bs, so freeze checks track physics)
ALLOCATE(jphi_ind(0:self%npsi))
jphi_ind = [self%j0, self%jphi]
CALL taper(jphi_bs)
CALL taper(jphi_pd)
CALL taper(jphi_ind)
CALL taper(jphi_fixed)
ALLOCATE(jphi_total(0:self%npsi))
!--- 4. alpha: I_p is affine in alpha, so alpha = (Ip_target - Ip_lo) / (Ip_hi - Ip_lo). Re-solved
!   every update (cheap): a frozen alpha lets I_p drift while the shape is still settling.
ip_target = ABS(gseq%Ip_target)
jphi_total = jphi_bs + jphi_fixed
CALL eval_dIdA(gseq, R_spline, psi_n, self%npsi+1, jphi_total, pscale, jtor)
CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_lo, extrap=1)
jphi_total = jphi_ind + jphi_bs + jphi_fixed
CALL eval_dIdA(gseq, R_spline, psi_n, self%npsi+1, jphi_total, pscale, jtor)
CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, ip_result_hi, extrap=1)
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
!--- 5. Total current, profiles (A/m²), F*F' knots
jphi_total = alpha * jphi_ind + jphi_bs + jphi_fixed
self%boot_profs%psi_n = psi_n
self%boot_profs%jphi_total = jphi_total/mu0
self%boot_profs%jphi_bs_final = jphi_bs/mu0
self%boot_profs%jphi_pd_final = jphi_pd/mu0
self%boot_profs%jphi_ind_final = alpha * jphi_ind/mu0
self%boot_profs%jphi_fixed_final = jphi_fixed/mu0
CALL spline_eval(R_spline, 0.d0, 0) ! LCFS point for y0 calculation
pprime = gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(jphi_total(0) - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
DO i = 1, self%npsi
  CALL spline_eval(R_spline, xpsi(i), 0)
  pprime = gseq%P%fp(xpsi(i)*(gseq%plasma_bounds(2) - &
                                  gseq%plasma_bounds(1)) + &
                                  gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(jphi_total(i) - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
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
  !--- Side-by-side Ip comparison: FEM nonlinear solve vs profile flux integral
  CALL gs_itor_nl(gseq, itor_nl)
  CALL eval_dIdA(gseq, R_spline, psi_n, self%npsi+1, jphi_total, pscale, jtor)
  CALL gs_flux_int(gseq, psi_n, jtor, self%npsi+1, itor_flint, extrap=1)
  WRITE(*,'(A)') '  [jphi_bs_update] --- Ip comparison (current jphi_total) ---'
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(gs_itor_nl)    = ', itor_nl/mu0
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(flux_int)      = ', itor_flint/mu0
END IF
CALL spline_dealloc(R_spline)
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
IF(self%boot_ops%taper_edge_jbs)CALL apply_edge_taper(self%npsi+1, psi_n, a, &
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
INTEGER(i4), ALLOCATABLE :: j(:)
INTEGER(i4) :: jc
REAL(r8) :: rop(3),pt(3),grad(3)
REAL(r8) :: s,c,Bp2,Bt2,mod_B,bratio
!---Get dofs
ALLOCATE(j(self%lag_rep%nce))
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
DEALLOCATE(j)
END SUBROUTINE sauter_apply
!------------------------------------------------------------------------------
!> Compute factors required for Sauter bootstrap formula
!------------------------------------------------------------------------------
SUBROUTINE sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs,qprof,eps)
CLASS(gs_equil), INTENT(inout) :: gseq !< G-S object
INTEGER(i4), INTENT(in) :: nr !< Number of flux sample points
REAL(r8), INTENT(in) :: psi_q(nr) !< Location of flux sample points in normalised psi
REAL(r8), INTENT(out) :: fc(nr) !< Circulating particle fraction \f$ f_c \f$
REAL(r8), INTENT(out) :: r_avgs(nr,3) !< Flux surface averaged radial coords \f$<R>\f$, \f$<1/R>\f$, \f$<a>\f$
REAL(r8), INTENT(out) :: modb_avgs(nr,2) !< Flux surface averaged field \f$<|B|>\f$, \f$<|B|^2>\f$
REAL(r8), OPTIONAL, INTENT(out) :: qprof(nr) !< Safety factor q on each surface (avoids a separate gs_get_qprof call)
REAL(r8), OPTIONAL, INTENT(out) :: eps(nr) !< Inverse aspect ratio (R_max-R_min)/(R_max+R_min) on each surface
REAL(r8) :: psi_surf,rmax,x1,x2,raxis,zaxis,fpol,qpsi,h,h2,hf,ftu,ftl
REAL(r8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
TYPE(oft_lag_brinterp) :: psi_int
REAL(r8), PARAMETER :: tol=1.d-10
INTEGER(i4) :: j,cell
TYPE(sauter_interp), TARGET :: field
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
CALL set_tracer(1)
field%u=>gseq%psi
field%mag_axis=gseq%o_point
CALL field%setup(device%fe_rep)
active_tracer%neq=9
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
DO j=1,nr
  psi_surf=psi_q(j)*(x2-x1) + x1
  !--- Axis guard: flux surface degenerates at psi_N=1; prescribe analytically.
  IF(psi_q(j) == 1.d0)THEN
    fpol=gs_fpol(gseq,psi_surf)
    r_avgs(j,1)    = raxis
    r_avgs(j,2)    = 1.d0 / raxis
    r_avgs(j,3)    = 0.d0
    modb_avgs(j,1) = ABS(fpol) / raxis
    modb_avgs(j,2) = (fpol / raxis)**2
    fc(j)          = 1.d0
    IF(PRESENT(qprof).AND.(j>1))qprof(j) = qprof(j-1)
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
  fpol=gs_fpol(gseq,psi_surf)
  field%f_surf=fpol
  field%bmax=0.d0
  field%rmax_surf = -1.d30
  field%rmin_surf =  1.d30
  field%stage_1=.TRUE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  field%stage_1=.FALSE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  pt_last=pt
  IF(active_tracer%status/=1)THEN
    WRITE(*,*)j,pt
    CALL oft_warn("sauter_fc: Trace did not complete")
    CYCLE
  END IF
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
  IF(PRESENT(eps))   eps(j)   = (field%rmax_surf - field%rmin_surf) / (field%rmax_surf + field%rmin_surf)
END DO
CALL field%delete
CALL active_tracer%delete
CALL psi_int%delete()
END SUBROUTINE sauter_fc
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
!> Compute dy/dx with second-order-accurate finite differences on the native grid.
!!
!! Port of `numpy.gradient(y, x, edge_order=2)`: the interior uses the
!! three-point central stencil for the local (possibly non-uniform) spacing,
!! and both endpoints use the second-order one-sided stencil.  The first-order
!! endpoint difference numpy uses by default is badly inaccurate at the
!! magnetic axis and separatrix, where the error propagates straight into
!! on-axis/edge jphi_bs.
!!
!! This replaces an earlier shape-preserving PCHIP derivative.  On the uniform
!! grids this code uses, the PCHIP endpoint formula reduces algebraically to
!! the second-order one-sided difference, while its monotonicity clipping
!! flattens the pedestal and makes the interior measurably less accurate.
!!
!! The grid is validated rather than repaired: unlike the PCHIP routine it
!! replaces, this does not sort x or average over duplicated abscissae.
!! Duplicated or unsorted flux labels give an undefined derivative, so they
!! abort here instead of being silently smoothed.
!!
!! @param n     Number of input points (>= 3)
!! @param x     Independent variable (e.g. psi); strictly monotonic, either direction
!! @param y     Dependent variable sampled on x
!! @param dydx  Output: dy/dx at each point in x
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
!> Bootstrap current (Redl 2021) on the nodes psi_N (bootstrap.py calculate_profiles_and_bootstrap)
!!
!! NRL Coulomb logs, Koh ion collisionality, Redl jboot1 (no sign(q) factor: <j.B> is odd in F), L34=L31.
!!
!! @param gseq    Equilibrium object (must have Te, Ti, ne, ni, Zeff set)
!! @param n_psi   Number of flux surface samples
!! @param psi_N   Normalised psi, OFT convention (0=LCFS, 1=axis), any spacing, psi_N(1)>0
!! @param jphi_bs Output: bootstrap current density as TokaMaker jphi = <j_phi> [A/m^2] on psi_N grid
!! @param jphi_pd Output: its pressure-driven part [A/m^2]
!------------------------------------------------------------------------------
SUBROUTINE calculate_bootstrap(self, gseq, n_psi, psi_N, jphi_bs, jphi_pd)
CLASS(jphi_bs_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)  !< Normalised psi (OFT convention, 0=LCFS), psi_N(1)>0
REAL(r8), INTENT(out) :: jphi_bs(0:n_psi) !< jphi_bs(0) the LCFS value (extrapolated)
REAL(r8), INTENT(out) :: jphi_pd(0:n_psi) !< Pressure-driven part of jphi_bs
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
REAL(r8) :: jdotb_bs(n_psi)
REAL(r8) :: jdotb(0:n_psi) !< <j_BS.B> with its LCFS/axis values
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
! Get flux-surface geometry: fc, eps, q, and <R> in one tracing pass
CALL sauter_fc(gseq, n_psi, psi_N, fc, r_avgs_saut, modb_avgs_saut, qprof=qvals, eps=eps)
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
CALL redl_bootstrap(n_psi, Te, Ti, ne, ni, pe, pi_arr, Zeff, eps, ft, f, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    nu_e_star, nu_i_star, jdotb_bs)
! Convert <j_BS.B> to TokaMaker's jphi = <j_phi> = <R>P' + <1/R>FF'/mu0 (exact, see
! doc_tokamaker_current_conventions.md eq. A7): the field-aligned part F<1/R>/<B^2> * <j_BS.B>,
! plus the pressure-driven (diamagnetic + Pfirsch-Schlueter) part P'(<R> - F^2<1/R>/<B^2>),
! which is assigned to the bootstrap (as IMAS includes_bootstrap=true). P' as in jphi_bs_update.
pscale = gseq%pax_target/gseq%P%f(gseq%plasma_bounds(2))
jphi_bs(0) = 0.0_r8 ! LCFS value set by the extrapolation below
jphi_pd(0) = 0.0_r8
DO i = 1, n_psi
  IF(ABS(f(i)) > 0.0_r8 .AND. modb_avgs_saut(i,2) > 0.0_r8)THEN
    pprime = gseq%P%fp(psi_abs(i))*pscale/mu0
    jphi_pd(i) = pprime*(r_avgs_saut(i,1) - f(i)**2*r_avgs_saut(i,2)/modb_avgs_saut(i,2))
    jphi_bs(i) = jdotb_bs(i)*f(i)*r_avgs_saut(i,2)/modb_avgs_saut(i,2) + jphi_pd(i)
  ELSE
    jphi_pd(i) = 0.0_r8
    jphi_bs(i) = 0.0_r8
  END IF
END DO
jdotb = [0.0_r8, jdotb_bs]
CALL finish(jphi_bs)
CALL finish(jphi_pd)
CALL finish(jdotb)
self%boot_profs%jdotb_bs_raw = jdotb
self%boot_profs%jphi_bs_raw = jphi_bs
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
CALL extrap_jbs_boundaries(n_psi, psi_N, a)
END SUBROUTINE finish
END SUBROUTINE calculate_bootstrap
!------------------------------------------------------------------------------
!> Edge spike of a bootstrap profile (OFT convention), as isolate_edge_jbs / parameterize_jbs in
!! bootstrap.py calculate_profiles_and_bootstrap: the masked spike and, if @p parameterize, its fit
!------------------------------------------------------------------------------
SUBROUTINE edge_spike(n_psi, psi_N, jphi_bs, parameterize, diagnose, spike, spike_masked)
INTEGER(i4), INTENT(in) :: n_psi !< Number of nodes (without the LCFS)
REAL(r8), INTENT(in) :: psi_N(1:n_psi) !< Normalised psi (OFT convention, 0=LCFS)
REAL(r8), INTENT(in) :: jphi_bs(0:n_psi) !< Bootstrap profile, jphi_bs(0) at the LCFS
LOGICAL, INTENT(in) :: parameterize !< Fit the spike (skew-normal)?
LOGICAL, INTENT(in) :: diagnose !< Print the fit parameters?
REAL(r8), INTENT(out) :: spike(0:n_psi) !< Edge spike (the fit if @p parameterize)
REAL(r8), INTENT(out) :: spike_masked(0:n_psi) !< Masked spike, before any fit
REAL(r8) :: psi_N_std(0:n_psi), jphi_bs_std(0:n_psi), mask_std(0:n_psi), param_std(0:n_psi)
! analyze_bootstrap_edge_spike uses the standard convention (0=axis, 1=LCFS): flip in and out
psi_N_std(0:n_psi-1) = 1.0_r8 - psi_N(n_psi:1:-1)
psi_N_std(n_psi) = 1.0_r8
jphi_bs_std = jphi_bs(n_psi:0:-1)
IF(parameterize)THEN
  CALL analyze_bootstrap_edge_spike(n_psi+1, psi_N_std, jphi_bs_std, mask_std, param_std, diagnose=diagnose)
  spike = param_std(n_psi:0:-1)
ELSE
  CALL analyze_bootstrap_edge_spike(n_psi+1, psi_N_std, jphi_bs_std, mask_std)
  spike = mask_std(n_psi:0:-1)
END IF
spike_masked = mask_std(n_psi:0:-1)
END SUBROUTINE edge_spike
!------------------------------------------------------------------------------
!> Linearly extrapolate a bootstrap profile to the LCFS (and the axis, if a node) where q is
!! undefined, using the second-order gradient at the nearest well-defined point.
!------------------------------------------------------------------------------
SUBROUTINE extrap_jbs_boundaries(n_psi, psi_N, prof)
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(1:n_psi)
REAL(r8), INTENT(inout) :: prof(0:n_psi)
INTEGER(i4) :: n_lo, n_hi
REAL(r8) :: dprof_dpsi(n_psi)
! psi_N(0) = 0.0 at the LCFS is implicit.
n_lo = 1
n_hi = MERGE(n_psi-1, n_psi, psi_N(n_psi) == 1.0_r8)
IF(n_hi-n_lo+1 < 3) CALL oft_abort('extrap_jbs_boundaries: too few '// &
  'points with well-defined q to extrapolate the bootstrap current to the LCFS', &
  'extrap_jbs_boundaries', __FILE__)
CALL gradient_(n_hi-n_lo+1, psi_N(n_lo:n_hi), prof(n_lo:n_hi), &
  dprof_dpsi(n_lo:n_hi))
prof(0) = prof(n_lo) + dprof_dpsi(n_lo) * (0.0_r8 - psi_N(n_lo))
IF(psi_N(n_psi) == 1.0_r8)THEN
  n_lo = 1
  n_hi = n_psi-1
  IF(n_hi-n_lo+1 < 3) CALL oft_abort('extrap_jbs_boundaries: too few '// &
    'points with well-defined q to extrapolate the bootstrap current to the axis', &
    'extrap_jbs_boundaries', __FILE__)
  CALL gradient_(n_hi-n_lo+1, psi_N(n_lo:n_hi), prof(n_lo:n_hi), &
    dprof_dpsi(n_lo:n_hi))
  prof(n_psi) = prof(n_hi) + dprof_dpsi(n_hi) * (psi_N(n_psi) - psi_N(n_hi))
ENDIF
END SUBROUTINE extrap_jbs_boundaries
!------------------------------------------------------------------------------
!> Evaluate the skew-normal PDF at a single point.
!!
!! Equivalent to scipy.stats.skewnorm.pdf(x, sk, loc=ctr, scale=scl):
!!   f = (2/scl) * phi((x-ctr)/scl) * Phi(sk*(x-ctr)/scl)
!! where phi and Phi are the standard normal PDF and CDF.
!!
!! @param x   Evaluation point
!! @param sk  Skewness parameter
!! @param ctr Location (mean shift)
!! @param scl Scale (width)
!! @result    PDF value at x
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
!! Mirrors numpy.logaddexp(a, b).
!------------------------------------------------------------------------------
PURE FUNCTION safe_logaddexp(a, b) RESULT(y)
REAL(r8), INTENT(in) :: a, b
REAL(r8) :: y, c
c = MAX(a, b)
y = c + LOG(EXP(a - c) + EXP(b - c))
END FUNCTION safe_logaddexp
!------------------------------------------------------------------------------
!> Locate the peak of the skew-normal shape used in parametrise_edge_jbs
!! using ternary search on the interval
!!   [max(0, center - 3*width), min(1, center + 3*width)]
!!
!! The skew-normal PDF is strictly unimodal for all skewness values (Azzalini
!! 1985), so ternary search is guaranteed to converge.  Each iteration cuts
!! the search interval by a factor of 2/3; 100 iterations reduce an initial
!! width of ~6*width to < 10^{-17} * initial_width (essentially machine
!! precision), while requiring only 200 PDF evaluations — 50x fewer than the
!! previous 10 000-point linear scan.
!!
!! @param center       Gaussian centre (loc)
!! @param width        Gaussian width (scale)
!! @param sk           Skewness parameter
!! @param x_peak       Output: psi_N at the raw-shape peak
!! @param val_peak_raw Output: raw-shape value at the peak
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
!!
!! Direct translation of the inner function generate_baseline_prof from
!! Python bootstrap.py:parameterize_edge_jbs.  Constructs the profile by
!! stitching a left-side SoftMax-blended skew-normal spike with a right-side
!! cosine (or cosh) decay at x_peak:
!!
!!   left  (psi <= x_peak):  SoftMax(offset_in, amp*skewnorm/val_peak_raw)
!!   right (psi >  x_peak):  stitch_height * cos(omega*(psi-x_peak))^tail_alpha
!!
!! x_peak and val_peak_raw must be pre-computed via find_skewnorm_peak so
!! the golden-section scan is not repeated for every function evaluation.
!!
!! @param n            Number of psi points
!! @param psi          psi_N grid [0, 1]
!! @param amp          Amplitude of the skew-normal spike
!! @param center       Spike centre in psi_N
!! @param width        Spike width (sigma = FWHM/2.355)
!! @param offset_in    Flat baseline J_BS level left of the spike
!! @param sk           Skewness parameter (a in skewnorm.pdf)
!! @param y_sep        Profile value prescribed at the separatrix (psi_N = 1)
!! @param blend_width  SoftMax blend width (sharpness of left-side stitch)
!! @param tail_alpha   Right-side decay exponent (>= 1)
!! @param x_peak       Peak location of the raw skew-normal (from find_skewnorm_peak)
!! @param val_peak_raw Peak value   of the raw skew-normal (from find_skewnorm_peak)
!! @param profile      Output: profile values on the psi grid
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
!> Fortran equivalent of Python bootstrap.py:parameterize_edge_jbs.
!!
!! Generates the parameterised edge-bootstrap profile on an n-point psi_N
!! grid.  When offset /= 0 an internal 1-D root find (bisection on alpha)
!! ensures that profile(psi(1)) matches the requested offset level, mirroring
!! the root_scalar(brentq) call in the Python.  When offset == 0 the function
!! falls through directly to eval_baseline_profile with a very negative
!! effective offset (reproducing the Python `-1e10` shortcut).
!!
!! Fixed parameter: tail_alpha = 1.5 (Python default, not a fit variable).
!!
!! @param n           Number of psi_N points
!! @param psi         Normalised poloidal flux grid [0, 1]
!! @param amp         Spike amplitude
!! @param center      Spike centre in psi_N
!! @param width       Spike width (sigma of skew-normal)
!! @param offset      Requested flat-core level
!! @param sk          Skewness parameter
!! @param y_sep       Value at the separatrix (psi_N = 1)
!! @param blend_width SoftMax blending width
!! @param tail_alpha  Right-side fall-off exponent
!! @param profile     Output: profile values on the psi grid
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
!!
!! Computes fvec(i) = parametrise_edge_jbs(psi_fit(i); cofs) - j_fit(i)
!! where cofs(7) = (amp, center, width, offset, sk, y_sep, blend_width) and
!! the fitting data are taken from the module-level active_edge_jbs context.
!!
!! Signature required by lmdif:  SUBROUTINE fcn(m, n, x, fvec, iflag)
!------------------------------------------------------------------------------
SUBROUTINE edge_jbs_residual(m, n, cofs, err, iflag)
INTEGER(i4),  INTENT(in)    :: m, n
REAL(r8),     INTENT(in)    :: cofs(n)
REAL(r8),     INTENT(out)   :: err(m)
INTEGER(i4),  INTENT(inout) :: iflag
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
!! of (psi_fit, jphi_bs_fit) data points using MINPACK lmdif.
!!
!! Equivalent to the scipy.optimize.curve_fit call in Python's
!! analyze_bootstrap_edge_spike (using the same 7-parameter model but without
!! parameter bounds, which lmdif does not support natively).
!!
!! The 7 parameters and their Python / curve_fit correspondences are:
!!   p(1)  amp         -- spike amplitude
!!   p(2)  center      -- spike centre in psi_N
!!   p(3)  width       -- spike sigma
!!   p(4)  offset      -- flat-core level
!!   p(5)  sk          -- skewness
!!   p(6)  y_sep       -- value at psi_N = 1
!!   p(7)  blend_width -- SoftMax blending width
!!
!! @param n_fit    Number of data points (m == n_fit for lmdif)
!! @param psi_fit  psi_N values of the fitting points
!! @param jphi_bs_fit Target jphi_bs values [A/m^2] at psi_fit
!! @param p0       Initial parameter guess (7 elements)
!! @param popt     Output: optimised parameters (7 elements)
!------------------------------------------------------------------------------
SUBROUTINE curve_fit_edge_jbs(n_fit, psi_fit, jphi_bs_fit, p0, popt, &
    lb_in, ub_in, success)
INTEGER(i4), INTENT(in)  :: n_fit
REAL(r8),    INTENT(in)  :: psi_fit(n_fit), jphi_bs_fit(n_fit)
REAL(r8),    INTENT(in)  :: p0(7)
REAL(r8),    INTENT(out) :: popt(7)
REAL(r8),    INTENT(in)  :: lb_in(7), ub_in(7)
LOGICAL, OPTIONAL, INTENT(out) :: success  !< lmdif convergence flag
!---MINPACK variables (same pattern as gs_psi2pt / circle_interp)
INTEGER(i4), PARAMETER :: NPARAMS = 7
REAL(r8) :: ftol, xtol, gtol, epsfcn, factor
REAL(r8) :: cofs(NPARAMS), error_vec(n_fit)
INTEGER(i4) :: k
REAL(r8), ALLOCATABLE :: diag(:), wa1(:), wa2(:), wa3(:), wa4(:), qtf(:)
REAL(r8), ALLOCATABLE :: fjac(:,:)
INTEGER(i4), ALLOCATABLE :: ipvt(:)
INTEGER(i4) :: maxfev, mode, nprint, info, nfev, ldfjac, ncons, ncofs  ! nfev: function eval count (informational)
CHARACTER(len=80) :: char_buf
!---
! Populate module-level context so edge_jbs_residual can access the data
IF (ALLOCATED(active_edge_jbs%psi_fit)) DEALLOCATE(active_edge_jbs%psi_fit)
IF (ALLOCATED(active_edge_jbs%j_fit))   DEALLOCATE(active_edge_jbs%j_fit)
ALLOCATE(active_edge_jbs%psi_fit(n_fit))
ALLOCATE(active_edge_jbs%j_fit(n_fit))
active_edge_jbs%n_fit   = n_fit
active_edge_jbs%psi_fit = psi_fit
active_edge_jbs%j_fit   = jphi_bs_fit
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
!!   5. (Optional) Fits parametrise_edge_jbs to jphi_bs over the spike region
!!      (psi_N >= psi_N(lmin_idx)) via curve_fit_edge_jbs (MINPACK lmdif).
!!      Only performed when parameterized_spike is present.  An approximate
!!      FWHM is computed first as the initial sigma guess.  The 7 fitted
!!      parameters are used to evaluate parameterized_spike on the full
!!      n-point psi_N grid.
!!
!! @param n                  Number of flux surface samples
!! @param psi_N              Normalised poloidal flux grid [0, 1]
!! @param jphi_bs           Bootstrap current density j_phi profile [A/m^2]
!! @param masked_spike       Output: isolated edge spike spliced onto flat core [A/m^2]
!! @param parameterized_spike Output: parametrise_edge_jbs fit evaluated on full psi_N grid [A/m^2]
!------------------------------------------------------------------------------
SUBROUTINE analyze_bootstrap_edge_spike(n, psi_N, jphi_bs, masked_spike, &
                                      parameterized_spike, diagnose)
INTEGER(i4), INTENT(in)  :: n
REAL(r8),    INTENT(in)  :: psi_N(n)              !< Normalised poloidal flux [0=axis, 1=LCFS/plasma edge]
REAL(r8),    INTENT(in)  :: jphi_bs(n)           !< Bootstrap current density [A/m^2]
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
REAL(r8)    :: peak_psi, peak_height, lmin_jbs, fwhm
REAL(r8)    :: jbs_min_loc, dist_to_peak, blend_width_val, x_start, x_end
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
  IF (PRESENT(parameterized_spike)) parameterized_spike = jphi_bs
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
  IF (PRESENT(parameterized_spike)) parameterized_spike = jphi_bs
  RETURN
END IF
peak_psi    = psi_N(peak_idx)
peak_height = jphi_bs(peak_idx)
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
  IF (PRESENT(parameterized_spike)) parameterized_spike = jphi_bs
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
! =====================================================================
! 5. Fit parametrise_edge_jbs to jphi_bs over the edge-spike region.
!    Only performed when parameterized_spike output is requested.
! =====================================================================
IF (PRESENT(parameterized_spike)) THEN
  ! FWHM scan for LM fit initial guess
  half_max  = 0.5_r8 * peak_height
  left_idx  = lmin_idx   ! fallback: spike left boundary if no half-max crossing found before psi=0.7
  DO i = peak_idx-1, 1, -1
    IF (psi_N(i) < 0.7_r8) EXIT
    IF (jphi_bs(i) <= half_max) THEN
      left_idx = i
      EXIT
    END IF
  END DO
  right_idx = n       ! fallback: spike right boundary if no half-max crossing found before psi=1.0
  DO i = peak_idx+1, n
    IF (jphi_bs(i) <= half_max) THEN
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
      j_fit(n_fit)   = jphi_bs(i)
    END IF
  END DO
  ! Bounds for LM fit (mirrors Python bootstrap.py analyze_bootstrap_edge_spike)
  ! amp bounds
  eps_tmp = MAX(1.0e-6_r8, 1.0e-3_r8 * ABS(peak_height))
  IF (peak_height == 0.0_r8) eps_tmp = 1.0e-6_r8
  amp_lo  = MIN(0.99_r8 * peak_height, peak_height - eps_tmp)
  amp_hi  = MAX(1.01_r8 * peak_height, peak_height + eps_tmp)
  ! offset bounds
  eps_tmp = MAX(1.0e-6_r8, 1.0e-3_r8 * ABS(lmin_jbs))
  IF (lmin_jbs == 0.0_r8) eps_tmp = 1.0e-6_r8
  off_lo  = MIN(0.99_r8 * lmin_jbs, lmin_jbs - eps_tmp)
  off_hi  = MAX(1.01_r8 * lmin_jbs, lmin_jbs + eps_tmp)
  ! y_sep upper bound: max(2|jphi_bs(n)|, |jphi_bs(n)|+1e-6, 1e-6)
  ysep_hi = MAX(2.0_r8 * ABS(jphi_bs(n)), ABS(jphi_bs(n)) + 1.0e-6_r8, 1.0e-6_r8)
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
  p0(4) = lmin_jbs                            ! offset
  p0(5) = 1.0_r8                              ! sk
  p0(6) = MAX(0.0_r8, jphi_bs(n))            ! y_sep
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