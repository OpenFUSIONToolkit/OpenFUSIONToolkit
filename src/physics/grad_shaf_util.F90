!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_gs_util.F90
!
!> Grad-Shafranov utility subroutines for TokaMaker
!!
!! @authors Chris Hansen
!! @date September 2017
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE oft_gs_util
USE, INTRINSIC :: iso_c_binding, ONLY: C_LOC
USE oft_base
USE spline_mod
USE oft_io, ONLY: hdf5_create_file, hdf5_create_group, hdf5_write, hdf5_read, &
  hdf5_field_get_sizes, hdf5_field_exist
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_blag_geval
USE oft_blag_operators, ONLY: oft_blag_project, oft_lag_brinterp, oft_lag_bginterp, &
  oft_blag_vproject
USE tracing_2d, ONLY: active_tracer, tracinginv_fs, set_tracer
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_factory, flux_func, gs_dflux, gs_itor_nl, gs_test_bounds, gs_b_interp, &
  gs_get_qprof, gsinv_interp, gs_psi2r, gs_psi2pt, gs_epsilon, gs_update_bounds
USE oft_gs_profiles
USE grad_shaf_prof_phys, ONLY: create_jphi_ff, jphi_flux_func
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
type, extends(gsinv_interp) :: sauter_interp
  logical :: stage_1 = .FALSE.
  real(8) :: f_surf = 0.d0
  real(8) :: bmax = -1.d0
  real(8) :: mag_axis(2) = 0.d0
contains
  !> Reconstruct the gradient of a Lagrange scalar field
  procedure :: interp => sauter_apply
end type sauter_interp
!
CLASS(gs_factory), POINTER :: gs_fit => NULL()
CLASS(flux_func), POINTER :: ff_fit => NULL()
REAL(8), POINTER, DIMENSION(:,:) :: tmpprof => NULL()
CONTAINS
!---------------------------------------------------------------------------------
!> Create flux function object from profile name
!---------------------------------------------------------------------------------
SUBROUTINE gs_profile_alloc(profType,F)
CHARACTER(LEN=*), INTENT(in) :: profType !< Profile type
CLASS(flux_func), POINTER, INTENT(out) :: F !< Flux function object
!---
SELECT CASE(TRIM(profType))
  CASE("zero")
    ALLOCATE(zero_flux_func::F)
  CASE("flat")
    ALLOCATE(flat_flux_func::F)
  CASE("poly")
    ALLOCATE(poly_flux_func::F)
  CASE("linterp")
    ALLOCATE(linterp_flux_func::F)
  CASE("jphi-linterp")
    ALLOCATE(jphi_flux_func::F)
  CASE("wesson")
    ALLOCATE(wesson_flux_func::F)
  CASE DEFAULT
    CALL oft_abort('Invalid profile type.','gs_profile_load',__FILE__)
END SELECT
END SUBROUTINE gs_profile_alloc
!---------------------------------------------------------------------------------
!> Create flux function object from definition file
!---------------------------------------------------------------------------------
SUBROUTINE gs_profile_load(filename,F)
CHARACTER(LEN=*), INTENT(in) :: filename !< File storing function definition
CLASS(flux_func), POINTER, INTENT(out) :: F !< Flux function object
INTEGER(4) :: io_unit
! !---
! REAL(8) :: alpha,beta,sep
! REAL(8), ALLOCATABLE, DIMENSION(:) :: cofs,yvals
! INTEGER(4) :: ncofs,npsi,io_unit
! LOGICAL :: zero_grad
CHARACTER(LEN=15) :: profType
!---
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
READ(io_unit,*)profType
CALL gs_profile_alloc(profType,F)
CALL F%load(io_unit)
CLOSE(io_unit)
! !---
! SELECT CASE(TRIM(profType))
!   CASE("zero")
!     ALLOCATE(zero_flux_func::F)
!   CASE("flat")
!     ALLOCATE(flat_flux_func::F)
!     ! CALL create_flat_f(F)
!   CASE("linear")
!     CALL oft_abort('"linear" profile no longer supported.','gs_profile_load',__FILE__)
!   CASE("poly")
!     ALLOCATE(poly_flux_func::F)
!     ! READ(io_unit,*)ncofs,zero_grad
!     ! ALLOCATE(cofs(ncofs))
!     ! READ(io_unit,*)cofs
!     ! CALL create_poly_ff(F,ncofs,cofs,zero_grad)
!     ! DEALLOCATE(cofs)
!   CASE("linterp")
!     ALLOCATE(linterp_flux_func::F)
!     ! READ(io_unit,*)ncofs,alpha
!     ! ALLOCATE(cofs(ncofs),yvals(ncofs))
!     ! READ(io_unit,*)cofs
!     ! READ(io_unit,*)yvals
!     ! CALL create_linterp_ff(F,ncofs,cofs,yvals,alpha)
!     ! DEALLOCATE(cofs,yvals)
!   CASE("jphi-linterp")
!     ALLOCATE(jphi_flux_func::F)
!     ! READ(io_unit,*)ncofs,alpha
!     ! ALLOCATE(cofs(ncofs),yvals(ncofs))
!     ! READ(io_unit,*)cofs
!     ! READ(io_unit,*)yvals
!     ! CALL create_jphi_ff(F,ncofs,cofs,yvals,alpha)
!     ! DEALLOCATE(cofs,yvals)
!   CASE("idcd")
!     CALL oft_abort('"idcd" profile no longer supported.','gs_profile_load',__FILE__)
!   CASE("step-slant")
!     CALL oft_abort('"step-slant" profile no longer supported.','gs_profile_load',__FILE__)
!   CASE("wesson")
!     ALLOCATE(wesson_flux_func::F)
!     ! READ(io_unit,*)ncofs
!     ! READ(io_unit,*)beta
!     ! CALL create_wesson_ff(F,ncofs,beta)
!   CASE DEFAULT
!     CALL oft_abort('Invalid profile type.','gs_profile_load',__FILE__)
! END SELECT
! CLOSE(io_unit)
END SUBROUTINE gs_profile_load
!---------------------------------------------------------------------------------
!> Save flux function object to definition file
!!
!! @param[in] filename File to store function definition
!! @param[in] f Flux function object
!---------------------------------------------------------------------------------
SUBROUTINE gs_profile_save(filename,F)
CHARACTER(LEN=*), INTENT(in) :: filename !< File to store function definition
CLASS(flux_func), POINTER, INTENT(in) :: F !< Flux function object
INTEGER(4) :: io_unit
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
! !---
! SELECT TYPE(this=>F)
!   TYPE IS(zero_flux_func)
!     WRITE(io_unit,*)"zero"
!   TYPE IS(flat_flux_func)
!     WRITE(io_unit,*)"flat"
!   TYPE IS(poly_flux_func)
!     WRITE(io_unit,*)"poly"
!     ! IF(this%zero_grad)THEN
!     !   WRITE(io_unit,*)this%deg+1,this%zero_grad
!     !   WRITE(io_unit,*)1.d0,this%cofs
!     ! ELSE
!     !   WRITE(io_unit,*)this%deg,this%zero_grad
!     !   WRITE(io_unit,*)this%cofs
!     ! END IF
!   TYPE IS(linterp_flux_func)
!     WRITE(io_unit,*)"linterp"
!     ! WRITE(io_unit,*)this%npsi,this%y0
!     ! WRITE(io_unit,*)this%x
!     ! WRITE(io_unit,*)this%yp
!   TYPE IS(jphi_flux_func)
!     WRITE(io_unit,*)"jphi-linterp"
!     ! WRITE(io_unit,*)this%npsi,this%y0
!     ! WRITE(io_unit,*)this%x
!     ! WRITE(io_unit,*)this%jphi
!   TYPE IS(wesson_flux_func)
!     WRITE(io_unit,*)"wesson"
!     ! WRITE(io_unit,*)this%ncofs
!     ! WRITE(io_unit,*)this%gamma
!   CLASS DEFAULT
!     CALL oft_abort('Invalid profile type.','gs_profile_save',__FILE__)
! END SELECT
CALL F%save(io_unit)
CLOSE(io_unit)
END SUBROUTINE gs_profile_save
!------------------------------------------------------------------------------
!> Compute various global quantities for Grad-Shafranov equilibrium
!------------------------------------------------------------------------------
subroutine gs_comp_globals(self,itor,centroid,vol,pvol,dflux,tflux,bp_vol)
class(gs_equil), intent(inout) :: self !< G-S object
real(8), intent(out) :: itor !< Toroidal current
real(8), intent(out) :: centroid(2) !< Current centroid [2]
real(8), intent(out) :: vol !< Plasma volume
real(8), intent(out) :: pvol !< \f$ \int P dV \f$
real(8), intent(out) :: dflux !< Diamagnetic flux
real(8), intent(out) :: tflux !< Contained toroidal flux
real(8), intent(out) :: bp_vol !< \f$ \int B_p^2 dV \f$
type(oft_lag_brinterp) :: psi_eval
type(oft_lag_bginterp) :: psi_geval
real(8) :: itor_loc,goptmp(3,3),v,psitmp(1),gpsitmp(3)
real(8) :: pt(3),curr_cent(2),Btor,Bpol(2)
integer(4) :: i,m
class(oft_bmesh), pointer :: smesh
type(gs_factory), pointer :: device
!---
device=>self%device
smesh=>device%mesh
psi_eval%u=>self%psi
CALL psi_eval%setup(device%fe_rep)
CALL psi_geval%shared_setup(psi_eval)
!---
itor = 0.d0
centroid = 0.d0
pvol = 0.d0
vol = 0.d0
dflux = 0.d0
tflux = 0.d0
bp_vol = 0.d0
!$omp parallel do private(m,goptmp,v,psitmp,gpsitmp,pt,itor_loc,Btor,Bpol) &
!$omp reduction(+:itor) reduction(+:centroid) reduction(+:pvol) reduction(+:vol) reduction(+:dflux) &
!$omp reduction(+:tflux) reduction(+:bp_vol)
do i=1,smesh%nc
  IF(smesh%reg(i)/=1)CYCLE
  do m=1,device%fe_rep%quad%np
    call smesh%jacobian(i,device%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,device%fe_rep%quad%pts(:,m),goptmp,psitmp)
    IF(psitmp(1)<self%plasma_bounds(1))CYCLE
    pt=smesh%log2phys(i,device%fe_rep%quad%pts(:,m))
    !---Compute Magnetic Field
    IF(gs_test_bounds(self,pt))THEN
      IF(self%mode==0)THEN
        itor_loc = (self%p_scale*pt(1)*self%P%Fp(psitmp(1)) &
        + self%I%Fp(psitmp(1))*((self%ffp_scale**2)*self%I%f(psitmp(1))+self%ffp_scale*self%I%f_offset)/(pt(1)+gs_epsilon))
      ELSE
        itor_loc = (self%p_scale*pt(1)*self%P%Fp(psitmp(1)) &
        + .5d0*self%ffp_scale*self%I%Fp(psitmp(1))/(pt(1)+gs_epsilon))
      END IF
      itor = itor + itor_loc*v*device%fe_rep%quad%wts(m)
      centroid = centroid + itor_loc*pt(1:2)*v*device%fe_rep%quad%wts(m)
      pvol = pvol + (self%p_scale*self%P%F(psitmp(1)))*v*device%fe_rep%quad%wts(m)*pt(1)
      vol = vol + v*device%fe_rep%quad%wts(m)*pt(1)
      !---Compute total toroidal Field
      IF(self%mode==0)THEN
        Btor = (self%ffp_scale*(self%I%F(psitmp(1))) + self%I%f_offset)/(pt(1)+gs_epsilon)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%ffp_scale*self%I%F(psitmp(1)) + self%I%f_offset**2))/(pt(1)+gs_epsilon)
      END IF
      tflux = tflux + Btor*v*device%fe_rep%quad%wts(m)
      !---Compute internal inductance
      call psi_geval%interp(i,device%fe_rep%quad%pts(:,m),goptmp,gpsitmp)
      Bpol = [gpsitmp(1),gpsitmp(2)]/(pt(1)+gs_epsilon)
      bp_vol = bp_vol + SUM(Bpol**2)*v*device%fe_rep%quad%wts(m)*pt(1)
      !---Compute differential toroidal Field
      IF(self%mode==0)THEN
        Btor = self%ffp_scale*(self%I%F(psitmp(1)))/pt(1)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%ffp_scale*self%I%F(psitmp(1)) + self%I%f_offset**2) &
        - self%I%f_offset)/pt(1)
      END IF
      dflux = dflux + Btor*v*device%fe_rep%quad%wts(m)
    END IF
  end do
end do
centroid = centroid/itor
bp_vol=2*pi*bp_vol
!
itor=itor*self%psiscale
pvol=pvol*self%psiscale*self%psiscale
bp_vol=bp_vol*self%psiscale*self%psiscale
dflux=dflux*self%psiscale
tflux=tflux*self%psiscale
CALL psi_eval%delete
CALL psi_geval%delete
end subroutine gs_comp_globals
!------------------------------------------------------------------------------
!> Compute plasma loop voltage
!------------------------------------------------------------------------------
subroutine gs_calc_vloop(self,vloop)
class(gs_equil), intent(inout) :: self !< G-S object
real(8), intent(out) :: vloop !< loop voltage
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
real(8) :: itor_loc ! local toroidal current in integration
real(8) :: itor ! toroidal current
real(8) :: I_NI ! non-inductive F*F'
real(8) :: eta_jsq ! eta*j_NI**2 
real(8) :: goptmp(3,3)
real(8) :: v ! volume
real(8) :: pt(3) ! radial coordinate
real(8) :: psitmp(1) ! magnetic flux coordinate
integer(4) :: i,m
class(oft_bmesh), pointer :: smesh
type(gs_factory), pointer :: device
!---
device=>self%device
smesh=>device%mesh
CALL self%eta%update(self) ! Make sure eta is up to date with current equilibrium
IF(ASSOCIATED(self%I_NI))CALL self%I_NI%update(self) ! Make sure I_NI is up to date with current equilibrium
psi_eval%u=>self%psi
CALL psi_eval%setup(device%fe_rep)
CALL psi_geval%shared_setup(psi_eval)
!---
eta_jsq = 0.d0
I_NI = 0.d0
itor = 0.d0
vloop = 0.d0
!!$omp parallel do private(m,goptmp,v,psitmp,gpsitmp,pt,itor_loc) &
!!$omp reduction(+:itor) reduction(+:vol) &
do i=1,smesh%nc
  IF(smesh%reg(i)/=1)CYCLE
  do m=1,device%fe_rep%quad%np
    call smesh%jacobian(i,device%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,device%fe_rep%quad%pts(:,m),goptmp,psitmp)
    IF(psitmp(1)<self%plasma_bounds(1))CYCLE
    pt=smesh%log2phys(i,device%fe_rep%quad%pts(:,m))
    !---Compute toroidal current itor, and eta*j^2 eta_jsq (numerator of Vloop integral)
    IF(gs_test_bounds(self,pt))THEN
      IF(ASSOCIATED(self%I_NI))I_NI=self%I_NI%Fp(psitmp(1))
      IF(self%mode==0)THEN
        itor_loc = (self%p_scale*pt(1)*self%P%Fp(psitmp(1)) &
          + self%I%Fp(psitmp(1))*((self%ffp_scale**2)*self%I%f(psitmp(1))+self%ffp_scale*self%I%f_offset - I_NI)/(pt(1)+gs_epsilon))
      ELSE
        itor_loc = (self%p_scale*pt(1)*self%P%Fp(psitmp(1)) &
          + (0.5d0*self%ffp_scale*self%I%Fp(psitmp(1)) - I_NI)/(pt(1)+gs_epsilon))
      END IF
      eta_jsq = eta_jsq + (itor_loc**2)*v*device%fe_rep%quad%wts(m)*pt(1)*self%eta%fp(psitmp(1))
      itor = itor + itor_loc*v*device%fe_rep%quad%wts(m)
    END IF
  end do
end do
eta_jsq=eta_jsq*(2*pi/(mu0*mu0))
itor=itor/mu0
!---Vloop = integral(eta_jsq) / itor
vloop=self%psiscale*(eta_jsq/itor)
!
CALL psi_eval%delete
CALL psi_geval%delete
end subroutine gs_calc_vloop
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine gs_save_tokamaker(self,filename)
class(gs_equil), intent(inout) :: self
character(LEN=*), intent(in) :: filename
real(r8), pointer :: vals_tmp(:)
integer(4) :: hash_tmp,io_unit
logical :: pm_save
IF(ASSOCIATED(self%P_ani))CALL oft_abort('Anisotropic pressure profiles not currently supported in tokamaker save.','gs_save_tokamaker',__FILE__)
CALL hdf5_create_file(filename)
! CALL hdf5_create_group(filename,'mesh')
! CALL hdf5_write(self%device%fe_rep%mesh%r,filename,'mesh/R')
! CALL hdf5_write(self%device%fe_rep%mesh%lc,filename,'mesh/LC')
CALL hdf5_create_group(filename,'tokamaker')
CALL hdf5_write(self%device%fe_rep%order,filename,'tokamaker/FE_ORDER')
hash_tmp = oft_simple_hash(C_LOC(self%device%mesh%lc),INT(4*3*self%device%mesh%nc,8))
CALL hdf5_write(hash_tmp,filename,'tokamaker/LC_HASH')
hash_tmp = oft_simple_hash(C_LOC(self%device%mesh%reg),INT(4*self%device%mesh%nc,8))
CALL hdf5_write(hash_tmp,filename,'tokamaker/REG_HASH')
CALL hdf5_write(self%device%ncoils,filename,'tokamaker/NCOILS')
!---
NULLIFY(vals_tmp)
CALL self%psi%get_local(vals_tmp)
CALL hdf5_write(vals_tmp,filename,'tokamaker/PSI')
CALL hdf5_write(self%coil_currs,filename,'tokamaker/COIL_CURRENTS')
!---
CALL hdf5_write(self%ffp_scale,filename,'tokamaker/FFP_SCALE')
CALL hdf5_write(self%I%f_offset,filename,'tokamaker/F0')
CALL hdf5_create_group(filename,'tokamaker/FFP_PROFILE')
CALL self%I%save(filename,'tokamaker/FFP_PROFILE')
CALL hdf5_write(self%p_scale,filename,'tokamaker/P_SCALE')
CALL hdf5_create_group(filename,'tokamaker/PP_PROFILE')
CALL self%P%save(filename,'tokamaker/PP_PROFILE')
IF(ASSOCIATED(self%I_NI))THEN
  CALL hdf5_create_group(filename,'tokamaker/NI_PROFILE')
  CALL self%I_NI%save(filename,'tokamaker/NI_PROFILE')
END IF
IF(ASSOCIATED(self%eta))THEN
  CALL hdf5_create_group(filename,'tokamaker/ETA_PROFILE')
  CALL self%eta%save(filename,'tokamaker/ETA_PROFILE')
END IF
! IF(ASSOCIATED(self%P_ani))THEN
!   CALL hdf5_create_group(filename,'tokamaker/P_ANI')
!   CALL self%P_ani%save(filename,'tokamaker/P_ANI')
! END IF
!---
CALL hdf5_write(self%plasma_bounds,filename,'tokamaker/PSI_BOUNDS')
CALL hdf5_write(self%o_point,filename,'tokamaker/O_POINT')
CALL hdf5_write(self%lim_point,filename,'tokamaker/LIM_POINT')
CALL hdf5_write(self%diverted,filename,'tokamaker/DIVERTED')
!---
IF(self%Itor_target>0.d0)CALL hdf5_write(self%Itor_target,filename,'tokamaker/IP_TARGET')
IF(self%Ip_ratio_target>-1.d98)CALL hdf5_write(self%Ip_ratio_target,filename,'tokamaker/IP_RATIO_TARGET')
IF(self%R0_target>0.d0)CALL hdf5_write(self%R0_target,filename,'tokamaker/R0_TARGET')
IF(self%Z0_target>-1.d98)CALL hdf5_write(self%Z0_target,filename,'tokamaker/Z0_TARGET')
IF(self%pax_target>0.d0)CALL hdf5_write(self%pax_target,filename,'tokamaker/PAX_TARGET')
IF(self%estore_target>0.d0)CALL hdf5_write(self%estore_target,filename,'tokamaker/ESTORE_TARGET')
IF(self%isoflux_ntargets>0)CALL hdf5_write(self%isoflux_targets,filename,'tokamaker/ISOFLUX_TARGETS')
IF(self%flux_ntargets>0)CALL hdf5_write(self%flux_targets,filename,'tokamaker/FLUX_TARGETS')
IF(self%saddle_ntargets>0)CALL hdf5_write(self%saddle_targets,filename,'tokamaker/SADDLE_TARGETS')
end subroutine gs_save_tokamaker
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine gs_load_tokamaker(self,filename,error_string)
class(gs_equil), intent(inout) :: self
character(LEN=*), intent(in) :: filename
CHARACTER(LEN=OFT_ERROR_SLEN), intent(out) :: error_string
integer(4) :: hash_tmp,hash_in,io_unit,int_tmp,ndims
integer(i4), allocatable :: dim_sizes(:)
REAL(r8) :: pt_tmp(2),scale_tmp,diff_val
real(r8), pointer :: vals_tmp(:)
logical :: success,pm_save,logical_tmp
CHARACTER(LEN=:), ALLOCATABLE :: profType
! CALL hdf5_write(self%device%fe_rep%mesh%r,filename,'mesh/R')
! CALL hdf5_write(self%device%fe_rep%mesh%lc,filename,'mesh/LC')
!---Check compatibility of mesh and FE representation
CALL hdf5_read(int_tmp,filename,'tokamaker/FE_ORDER',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read FE order.'
  RETURN
END IF
IF(int_tmp/=self%device%fe_rep%order)THEN
  error_string='FE order of equilibrium does not match current device.'
  RETURN
END IF
hash_tmp = oft_simple_hash(C_LOC(self%device%mesh%lc),INT(4*3*self%device%mesh%nc,8))
CALL hdf5_read(hash_in,filename,'tokamaker/LC_HASH',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read cell list hash.'
  RETURN
END IF
IF(hash_tmp/=hash_in)THEN
  error_string='Cell list hash for equilibrium does not match current device.'
  RETURN
END IF
hash_tmp = oft_simple_hash(C_LOC(self%device%mesh%reg),INT(4*self%device%mesh%nc,8))
CALL hdf5_read(hash_in,filename,'tokamaker/REG_HASH',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read region hash.'
  RETURN
END IF
IF(hash_tmp/=hash_in)THEN
  error_string='Region hash for equilibrium does not match current device.'
  RETURN
END IF
CALL hdf5_read(int_tmp,filename,'tokamaker/NCOILS',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read number of coils.'
  RETURN
END IF
IF(int_tmp/=self%device%ncoils)THEN
  error_string='Number of coils for equilibrium does not match current device.'
  RETURN
END IF
!---Load psi and update bounds
NULLIFY(vals_tmp)
CALL self%psi%get_local(vals_tmp)
CALL hdf5_read(vals_tmp,filename,'tokamaker/PSI',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read PSI values.'
  RETURN
END IF
CALL self%psi%restore_local(vals_tmp)
CALL hdf5_read(self%coil_currs,filename,'tokamaker/COIL_CURRENTS',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read coil currents.'
  RETURN
END IF
CALL gs_update_bounds(self)
!---Load flux functions
CALL hdf5_read(self%ffp_scale,filename,'tokamaker/FFP_SCALE',success=success)
IF(.NOT.success)THEN
  error_string="Failed to read F*F' scale."
  RETURN
END IF
CALL hdf5_read(profType,filename,'tokamaker/FFP_PROFILE/TYPE',success=success)
IF(.NOT.success)THEN
  error_string="Failed to read F*F' profile type."
  RETURN
END IF
IF(ASSOCIATED(self%I))THEN
  CALL self%I%delete()
  DEALLOCATE(self%I)
END IF
CALL gs_profile_alloc(profType,self%I)
DEALLOCATE(profType)
CALL self%I%load(filename,'tokamaker/FFP_PROFILE',success=success)
CALL self%I%update(self)
IF(.NOT.success)THEN
  error_string="Failed to load F*F' profile."
  RETURN
END IF
CALL hdf5_read(self%I%f_offset,filename,'tokamaker/F0',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read F0 values.'
  RETURN
END IF
CALL hdf5_read(self%p_scale,filename,'tokamaker/P_SCALE',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read P scale.'
  RETURN
END IF
CALL hdf5_read(profType,filename,'tokamaker/PP_PROFILE/TYPE',success=success)
IF(.NOT.success)THEN
  error_string="Failed to read P' profile type."
  RETURN
END IF
IF(ASSOCIATED(self%P))THEN
  CALL self%P%delete()
  DEALLOCATE(self%P)
END IF
CALL gs_profile_alloc(profType,self%P)
DEALLOCATE(profType)
CALL self%P%load(filename,'tokamaker/PP_PROFILE',success=success)
CALL self%P%update(self)
IF(.NOT.success)THEN
  error_string="Failed to load P' profile."
  RETURN
END IF
IF(hdf5_field_exist(filename,'tokamaker/NI_PROFILE'))THEN
  CALL hdf5_read(profType,filename,'tokamaker/NI_PROFILE/TYPE',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read non-inductive current profile type.'
    RETURN
  END IF
  IF(ASSOCIATED(self%I_NI))THEN
    CALL self%I_NI%delete()
    DEALLOCATE(self%I_NI)
  END IF
  CALL gs_profile_alloc(profType,self%I_NI)
  DEALLOCATE(profType)
  CALL self%I_NI%load(filename,'tokamaker/NI_PROFILE',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to load non-inductive current profile.'
    RETURN
  END IF
  CALL self%I_NI%update(self)
END IF
IF(hdf5_field_exist(filename,'tokamaker/ETA_PROFILE'))THEN
  CALL hdf5_read(profType,filename,'tokamaker/ETA_PROFILE/TYPE',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read ETA profile type.'
    RETURN
  END IF
  IF(ASSOCIATED(self%eta))THEN
    CALL self%eta%delete()
    DEALLOCATE(self%eta)
  END IF
  CALL gs_profile_alloc(profType,self%eta)
  DEALLOCATE(profType)
  CALL self%eta%load(filename,'tokamaker/ETA_PROFILE',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to load ETA profile.'
    RETURN
  END IF
  CALL self%eta%update(self)
END IF
! IF(hdf5_field_exist(filename,'tokamaker/P_ANI'))THEN
!   CALL hdf5_read(profType,filename,'tokamaker/P_ANI/TYPE',success=success)
!   IF(.NOT.success)GOTO 100
!   IF(ASSOCIATED(self%P_ani))THEN
!     CALL self%P_ani%delete()
!     DEALLOCATE(self%P_ani)
!   END IF
!   CALL gs_ani_alloc(profType,self%P_ani)
!   DEALLOCATE(profType)
!   CALL self%P_ani%load(filename,'tokamaker/P_ANI',success=success)
!   IF(.NOT.success)GOTO 100
!   CALL self%P_ani%update(self)
! END IF
!---Check consistency of values and warn if they differ beyond expected numerical differences
CALL hdf5_read(pt_tmp,filename,'tokamaker/PSI_BOUNDS',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read PSI bounds.'
  RETURN
END IF
diff_val=ABS(pt_tmp(2)-pt_tmp(1))*1.d-2
IF(ANY(ABS(pt_tmp-self%plasma_bounds)>diff_val))THEN
  CALL oft_warn('Plasma bounds in equilibrium file do not match recomputed values.')
END IF
CALL hdf5_read(pt_tmp,filename,'tokamaker/O_POINT',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read O-point.'
  RETURN
END IF
diff_val=self%device%mesh%hrms
IF(SQRT(SUM((pt_tmp-self%o_point)**2))>diff_val)THEN
  CALL oft_warn('O-point in equilibrium file does not match recomputed location.')
END IF
CALL hdf5_read(pt_tmp,filename,'tokamaker/LIM_POINT',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read Limiting point.'
  RETURN
END IF
diff_val=self%device%mesh%hrms
IF(SQRT(SUM((pt_tmp-self%lim_point)**2))>diff_val)THEN
  CALL oft_warn('Limiting point in equilibrium file does not match recomputed location.')
END IF
CALL hdf5_read(logical_tmp,filename,'tokamaker/DIVERTED',success=success)
IF(.NOT.success)THEN
  error_string='Failed to read Diverted status.'
  RETURN
END IF
IF(logical_tmp.NEQV.self%diverted)THEN
  CALL oft_warn('Diverted status in equilibrium file does not match recomputed value.')
END IF
!---Load targets and coil currents
IF(hdf5_field_exist(filename,'tokamaker/IP_TARGET'))THEN
  CALL hdf5_read(self%Itor_target,filename,'tokamaker/IP_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read Ip target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/IP_RATIO_TARGET'))THEN
  CALL hdf5_read(self%Ip_ratio_target,filename,'tokamaker/IP_RATIO_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read Ip ratio target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/R0_TARGET'))THEN
  CALL hdf5_read(self%R0_target,filename,'tokamaker/R0_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read R0 target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/Z0_TARGET'))THEN
  CALL hdf5_read(self%Z0_target,filename,'tokamaker/Z0_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read Z0 target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/PAX_TARGET'))THEN
  CALL hdf5_read(self%pax_target,filename,'tokamaker/PAX_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read P axis target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/ESTORE_TARGET'))THEN
  CALL hdf5_read(self%estore_target,filename,'tokamaker/ESTORE_TARGET',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read stored energy target.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/ISOFLUX_TARGETS'))THEN
  CALL hdf5_field_get_sizes(filename,'tokamaker/ISOFLUX_TARGETS',ndims,dim_sizes)
  IF(dim_sizes(1)/=5)CALL oft_abort('Invalid first dimension for isoflux targets', 'gs_load_tokamaker', __FILE__)
  self%isoflux_ntargets=dim_sizes(2)
  IF(ASSOCIATED(self%isoflux_targets))DEALLOCATE(self%isoflux_targets)
  ALLOCATE(self%isoflux_targets(5,self%isoflux_ntargets))
  DEALLOCATE(dim_sizes)
  CALL hdf5_read(self%isoflux_targets,filename,'tokamaker/ISOFLUX_TARGETS',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read isoflux targets.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/FLUX_TARGETS'))THEN
  CALL hdf5_field_get_sizes(filename,'tokamaker/FLUX_TARGETS',ndims,dim_sizes)
  IF(dim_sizes(1)/=4)CALL oft_abort('Invalid first dimension for flux targets', 'gs_load_tokamaker', __FILE__)
  self%flux_ntargets=dim_sizes(2)
  IF(ASSOCIATED(self%flux_targets))DEALLOCATE(self%flux_targets)
  ALLOCATE(self%flux_targets(4,self%flux_ntargets))
  DEALLOCATE(dim_sizes)
  CALL hdf5_read(self%flux_targets,filename,'tokamaker/FLUX_TARGETS',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read flux targets.'
    RETURN
  END IF
END IF
IF(hdf5_field_exist(filename,'tokamaker/SADDLE_TARGETS'))THEN
  CALL hdf5_field_get_sizes(filename,'tokamaker/SADDLE_TARGETS',ndims,dim_sizes)
  IF(dim_sizes(1)/=3)CALL oft_abort('Invalid first dimension for saddle targets', 'gs_load_tokamaker', __FILE__)
  self%saddle_ntargets=dim_sizes(2)
  IF(ASSOCIATED(self%saddle_targets))DEALLOCATE(self%saddle_targets)
  ALLOCATE(self%saddle_targets(3,self%saddle_ntargets))
  DEALLOCATE(dim_sizes)
  CALL hdf5_read(self%saddle_targets,filename,'tokamaker/SADDLE_TARGETS',success=success)
  IF(.NOT.success)THEN
    error_string='Failed to read saddle targets.'
    RETURN
  END IF
END IF
end subroutine gs_load_tokamaker
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine gs_save_ifile(gseq,filename,npsi,ntheta,psi_pad,lcfs_press,pack_lcfs,single_prec,error_str)
class(gs_equil), intent(inout) :: gseq !< G-S object
CHARACTER(LEN=OFT_PATH_SLEN), intent(in) :: filename !< Outpute filename
integer(4), intent(in) :: npsi !< Number of points in flux coordinate
integer(4), intent(in) :: ntheta !< Number of points in poloidal coordinate
REAL(8), intent(in) :: psi_pad !< Padding at LCFS in normalized units
REAL(8), optional, intent(in) :: lcfs_press !< LCFS pressure
LOGICAL, OPTIONAL, INTENT(in) :: pack_lcfs !< Use quadratic packing toward LCFS?
LOGICAL, OPTIONAL, INTENT(in) :: single_prec !< Save file with single precision fields?
CHARACTER(LEN=OFT_ERROR_SLEN), OPTIONAL, INTENT(out) :: error_str
type(gsinv_interp), pointer :: field
type(oft_lag_brinterp) :: psi_int
real(8) :: gop(3,3),psi_surf(1),pt_last(3)
real(8) :: raxis,zaxis,f(3),pt(3),rmax,x1,x2,xr
real(8), allocatable :: ptout(:,:)
real(8), allocatable :: rout(:,:),zout(:,:),cout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell,io_unit
LOGICAL :: do_pack,save_single
TYPE(spline_type) :: rz
type(gs_factory), pointer :: device
device=>gseq%device
!---
IF(PRESENT(error_str))error_str=""
WRITE(*,'(3A)')oft_indent,'Saving iFile: ',TRIM(filename)
CALL oft_increase_indent
do_pack=.FALSE.
save_single=.FALSE.
IF(PRESENT(pack_lcfs))do_pack=pack_lcfs
IF(PRESENT(single_prec))save_single=single_prec
!---
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
!
xr = (x2-x1)
x1 = x1 + xr*psi_pad
xr = (x2-x1)
!
psi_int%u=>gseq%psi
CALL psi_int%setup(device%fe_rep)
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  IF(device%dipole_mode.OR.device%mirror_mode)THEN
    pt=[raxis*j/REAL(100,8),0.d0,0.d0]
  ELSE
    pt=[(device%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  END IF
  CALL bmesh_findcell(device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF(device%dipole_mode.OR.device%mirror_mode)THEN
    IF(psi_surf(1)>x1)EXIT
  ELSE
    IF(psi_surf(1)<x1)EXIT
  END IF
  rmax=pt(1)
END DO
pt_last=[(.9d0*rmax+.1d0*raxis),zaxis,0.d0]
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Axis Position:'
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
  CALL oft_decrease_indent
END IF
!---Trace
call set_tracer(1)
ALLOCATE(cout(npsi,4))
ALLOCATE(rout(ntheta,npsi))
ALLOCATE(zout(ntheta,npsi))
!$omp parallel private(j,psi_surf,pt,ptout,field,rz,gop) firstprivate(pt_last)
ALLOCATE(field)
field%u=>gseq%psi
CALL field%setup(device%fe_rep)
active_tracer%neq=3
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
!$omp do schedule(dynamic,1)
do j=2,npsi
  IF(PRESENT(error_str))THEN
    IF(error_str/="")CYCLE
  END IF
  !---------------------------------------------------------------------------
  ! Trace contour
  !---------------------------------------------------------------------------
  IF(pack_lcfs)THEN
    psi_surf = xr*(1.d0-(j-1)/REAL(npsi-1,8))**2 + x1
  ELSE
    psi_surf = xr*(1.d0-(j-1)/REAL(npsi-1,8)) + x1
  END IF
  ! psi_surf(1)=x2 - psi_surf(1)
  IF(gseq%diverted.AND.ABS((psi_surf(1)-x1)/xr)<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  pt=pt_last
  !$omp critical
  CALL gs_psi2r(gseq,psi_surf(1),pt,psi_int=psi_int)
  !$omp end critical
  CALL tracinginv_fs(device%mesh,pt,ptout)
  pt_last=pt
  !---Exit if trace fails
  IF(active_tracer%status/=1)THEN
    IF(PRESENT(error_str))THEN
      !$omp critical
      WRITE(error_str,"(A,F10.4)")"Tracing failed at psi = ",1.d0-psi_surf
      !$omp end critical
      CYCLE
    ELSE
      call oft_abort('Trace did not complete.','gs_save_decon',__FILE__)
    END IF
  END IF
  !---------------------------------------------------------------------------
  ! Perform cubic spline interpolation
  !---------------------------------------------------------------------------
  !---Allocate spline
  CALL spline_alloc(rz,active_tracer%nsteps,2)
  !---Setup Spline
  rz%xs(0:active_tracer%nsteps) = ptout(1,1:active_tracer%nsteps+1)/ptout(1,active_tracer%nsteps+1)
  rz%fs(0:active_tracer%nsteps,1) = ptout(2,1:active_tracer%nsteps+1)
  rz%fs(0:active_tracer%nsteps,2) = ptout(3,1:active_tracer%nsteps+1)
  CALL spline_fit(rz,"periodic")
  !---Resample trace
  DO k=0,ntheta-1
    CALL spline_eval(rz,k/REAL(ntheta-1,8),0)
    rout(k+1,j)=rz%f(1)
    zout(k+1,j)=rz%f(2)
  END DO
  !---Destroy Spline
  CALL spline_dealloc(rz)
  !---------------------------------------------------------------------------
  ! Save DCON information
  !---------------------------------------------------------------------------
  cout(j,1)=psi_surf(1) ! Poloidal flux
  !---Toroidal flux function
  IF(gseq%mode==0)THEN
    cout(j,2)=gseq%ffp_scale*gseq%I%f(psi_surf(1))+gseq%I%f_offset
  ELSE
    cout(j,2)=SQRT(gseq%ffp_scale*gseq%I%f(psi_surf(1)) + gseq%I%f_offset**2) &
    + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  cout(j,3)=gseq%p_scale*gseq%P%f(psi_surf(1))/mu0 ! Plasma pressure
  cout(j,4)=cout(j,2)*active_tracer%v(3)/(2*pi) ! Safety Factor (q)
end do
CALL active_tracer%delete
CALL field%delete
DEALLOCATE(ptout,field)
!$omp end parallel
CALL psi_int%delete()
IF(PRESENT(error_str))THEN
  IF(error_str/="")THEN
    DEALLOCATE(cout,rout,zout)
    RETURN
  END IF
END IF
!---Information for O-point
rout(:,1)=raxis
zout(:,1)=zaxis
cout(1,1)=x2
IF(gseq%mode==0)THEN
  cout(1,2)=(gseq%ffp_scale*gseq%I%f(x2)+gseq%I%f_offset)
ELSE
  cout(1,2)=SQRT(gseq%ffp_scale*gseq%I%f(x2) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
END IF
cout(1,3)=gseq%p_scale*gseq%P%f(x2)/mu0
cout(1,4)=(cout(3,4)-cout(2,4))*(x2-cout(2,1))/(cout(3,1)-cout(2,1)) + cout(2,4)
!---Add LCFS pressure if specified
IF(PRESENT(lcfs_press))cout(:,3)=cout(:,3)+lcfs_press
!---------------------------------------------------------------------------
! Create output file
!---------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename),FORM='UNFORMATTED')
!---------------------------------------------------------------------------
! Write array lengths
!---------------------------------------------------------------------------
WRITE(io_unit)INT(npsi,4),INT(ntheta,4)
!---------------------------------------------------------------------------
! Write out flux surface quantities
!
! cout(:,1) -> psi(0:npsi)
! cout(:,2) -> f(0:npsi)
! cout(:,3) -> p(0:npsi)
! cout(:,4) -> q(0:npsi)
!---------------------------------------------------------------------------
DO j=1,4
  IF(save_single)THEN
    WRITE(io_unit)REAL(cout(:,j),4)
  ELSE
    WRITE(io_unit)cout(:,j)
  END IF
END DO
!---------------------------------------------------------------------------
! Write out inverse representation
!
! rout -> r(0:ntheta,0:npsi)
! zout -> z(0:ntheta,0:npsi)
!---------------------------------------------------------------------------
IF(save_single)THEN
  WRITE(io_unit)REAL(rout,4)
  WRITE(io_unit)REAL(zout,4)
ELSE
  WRITE(io_unit)rout
  WRITE(io_unit)zout
END IF
!---------------------------------------------------------------------------
! Close output file
!---------------------------------------------------------------------------
CLOSE(io_unit)
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,2ES11.3)')oft_indent,'Psi  = ',x1,x2
  WRITE(*,'(2A,F7.2)')oft_indent,'Qmin = ',MINVAL(cout(4,:))
  WRITE(*,'(2A,F7.2)')oft_indent,'Qmax = ',MAXVAL(cout(4,:))
  ! WRITE(*,'(2A)')oft_indent,'Done'
END IF
CALL oft_decrease_indent
!---
DEALLOCATE(cout,rout,zout)
end subroutine gs_save_ifile
!---------------------------------------------------------------------------
!> Save equilibrium to General Atomics gEQDSK file
!------------------------------------------------------------------------------
subroutine gs_save_eqdsk(gseq,filename,nr,nz,rbounds,zbounds,run_info,limiter_file,psi_pad,rcentr_in,trunc_eq,lcfs_press,cocos,error_str)
class(gs_equil), intent(inout) :: gseq !< Equilibrium to save
CHARACTER(LEN=OFT_PATH_SLEN), intent(in) :: filename !< Outpute filename
integer(4), intent(in) :: nr !< Number of radial points for flux/psi grid
integer(4), intent(in) :: nz !< Number of vertical points for flux grid
real(8), intent(in) :: rbounds(2) !< Radial extents for flux grid
real(8), intent(in) :: zbounds(2) !< Radial extents for flux grid
CHARACTER(LEN=40), intent(in) :: run_info !< Run information string [40]
CHARACTER(LEN=OFT_PATH_SLEN), intent(in) :: limiter_file !< Path to limiter file
REAL(8), intent(in) :: psi_pad !< Padding at LCFS in normalized units
REAL(8), optional, intent(in) :: rcentr_in !< Value to use for RCENTR (otherwise geometric center is used)
LOGICAL, OPTIONAL, INTENT(in) :: trunc_eq !< Truncate equilibrium at psi_pad
REAL(8), optional, intent(in) :: lcfs_press !< LCFS pressure
INTEGER, OPTIONAL, INTENT(in) :: cocos
CHARACTER(LEN=OFT_ERROR_SLEN), OPTIONAL, INTENT(out) :: error_str
!
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,xr,psi_trace
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:),rout(:),zout(:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,k,cell,io_unit,lim_max
type(gsinv_interp), pointer :: field
TYPE(spline_type) :: rz
type(gs_factory), pointer :: device
!---
INTEGER(4) :: nlim
REAL(8), ALLOCATABLE, DIMENSION(:) :: rlim,zlim
CHARACTER(LEN=48) :: eqdsk_case
INTEGER(4) :: idum
REAL(8) :: rdim,zdim,rleft,zmid,rcentr,bcentr,itor,xdum,rHFS,rLFS,zHFS,fptmp
REAL(8), ALLOCATABLE, DIMENSION(:) :: fpol,pres,ffprim,pprime,qpsi
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: psirz
LOGICAL :: do_truncate
!---
device=>gseq%device
IF(PRESENT(error_str))error_str=""
WRITE(*,'(3A)')oft_indent,'Saving gEQDSK: ',TRIM(filename)
CALL oft_increase_indent
!---
ALLOCATE(fpol(nr),pres(nr),ffprim(nr),pprime(nr),qpsi(nr))
ALLOCATE(psirz(nr,nz))
!---
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
do_truncate=.TRUE.
IF(PRESENT(trunc_eq))do_truncate=trunc_eq
xr = (x2-x1)
IF(do_truncate)THEN
  x1 = x1 + xr*psi_pad
  xr = (x2-x1)
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(device%fe_rep)
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  IF(device%dipole_mode)THEN
    pt=[raxis*j/REAL(100,8),0.d0,0.d0]
  ELSE
    pt=[(device%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  END IF
  CALL bmesh_findcell(device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_tmp)
  IF(device%dipole_mode.OR.device%mirror_mode)THEN
    IF(psi_tmp(1)>x1)EXIT
  ELSE
    IF(psi_tmp(1)<x1)EXIT
  END IF
  rmax=pt(1)
END DO
pt_last=[(.1d0*rmax+.9d0*raxis),zaxis,0.d0]
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Axis Position:'
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
  CALL oft_decrease_indent
END IF
!---Trace
call set_tracer(1)
ALLOCATE(rout(nr))
ALLOCATE(zout(nr))
!$omp parallel private(j,psi_surf,psi_trace,pt,ptout,field,fptmp) firstprivate(pt_last)
ALLOCATE(field)
field%u=>gseq%psi
CALL field%setup(device%fe_rep)
active_tracer%neq=3
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
!$omp do schedule(dynamic,1)
do j=1,nr
  IF(PRESENT(error_str))THEN
    IF(error_str/="")CYCLE
  END IF
  !------------------------------------------------------------------------------
  ! Trace contour
  !------------------------------------------------------------------------------
  psi_surf = x2 - xr*((j-1)/REAL(nr-1,8))
  psi_trace = psi_surf
  IF((.NOT.do_truncate).AND.((psi_trace-x1)/xr<psi_pad))psi_trace = x1 + xr*psi_pad
  IF(gseq%diverted.AND.(psi_trace-x1)/xr<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  IF(j>1)THEN
    pt=pt_last
    !$omp critical
    CALL gs_psi2r(gseq,psi_trace,pt,psi_int=psi_int)
    !$omp end critical
    IF(j==nr)THEN
      ALLOCATE(ptout(3,active_tracer%maxsteps+1))
      CALL tracinginv_fs(device%mesh,pt(1:2),ptout)
    ELSE
      CALL tracinginv_fs(device%mesh,pt(1:2))
    END IF
    pt_last=pt
    !---Exit if trace fails
    IF(active_tracer%status/=1)THEN
      IF(PRESENT(error_str))THEN
        !$omp critical
        WRITE(error_str,"(A,F10.4)")"Tracing failed at psi = ",1.d0-(psi_trace-x1)/xr
        !$omp end critical
        CYCLE
      ELSE
        call oft_abort('Trace did not complete.','gs_save_eqdsk',__FILE__)
      END IF
    END IF
  END IF
  IF(j==nr)THEN
    ! !---Get midplane extents
    ! rLFS=ptout(2,1)
    ! zHFS=1.d99
    ! DO k=2,active_tracer%nsteps+1
    !   IF(ABS(ptout(3,k)-zaxis)<zHFS)THEN
    !     zHFS=ABS(ptout(3,k)-zaxis)
    !     rHFS=ptout(2,k)
    !   END IF
    ! END DO
    !------------------------------------------------------------------------------
    ! Perform Cubic Spline Interpolation
    !------------------------------------------------------------------------------
    !---Allocate spline
    CALL spline_alloc(rz,active_tracer%nsteps,2)
    !---Setup Spline
    rz%xs(0:active_tracer%nsteps) = ptout(1,1:active_tracer%nsteps+1)/ptout(1,active_tracer%nsteps+1)
    rz%fs(0:active_tracer%nsteps,1) = ptout(2,1:active_tracer%nsteps+1)
    rz%fs(0:active_tracer%nsteps,2) = ptout(3,1:active_tracer%nsteps+1)
    CALL spline_fit(rz,"periodic")
    !---Resample trace
    DO k=0,nr-1
      CALL spline_eval(rz,k/REAL(nr-1,8),0)
      rout(k+1)=rz%f(1)
      zout(k+1)=rz%f(2)
    END DO
    !---Destroy Spline
    CALL spline_dealloc(rz)
    DEALLOCATE(ptout)
  END IF
  !------------------------------------------------------------------------------
  ! Compute Mercier Profiles
  !------------------------------------------------------------------------------
  !---Get flux variables
  IF(gseq%mode==0)THEN
    fptmp=gseq%ffp_scale*gseq%I%f(psi_trace)+gseq%I%f_offset
    fpol(j)=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
    ffprim(j)=gseq%I%fp(psi_surf)*((gseq%ffp_scale**2)*gseq%I%f(psi_surf)+gseq%ffp_scale*gseq%I%f_offset)
  ELSE
    fptmp=SQRT(gseq%ffp_scale*gseq%I%f(psi_trace) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
    fpol(j)=SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
    ffprim(j)=0.5d0*gseq%ffp_scale*gseq%I%fp(psi_surf)
  END IF
  pres(j)=gseq%p_scale*gseq%P%f(psi_surf)/mu0
  pprime(j)=gseq%p_scale*gseq%P%fp(psi_surf)/mu0
  !---Safety Factor (q)
  IF(j>1)qpsi(j)=fptmp*active_tracer%v(3)/(2*pi)
end do
CALL active_tracer%delete
CALL field%delete()
DEALLOCATE(field)
!$omp end parallel
IF(PRESENT(error_str))THEN
  IF(error_str/="")THEN
    DEALLOCATE(rout,zout)
    DEALLOCATE(fpol,pres,ffprim,pprime,qpsi,psirz)
    RETURN
  END IF
END IF
IF(PRESENT(lcfs_press))pres=pres+lcfs_press
!---Extrapolate q on axis
f(1) = x2
f(2) = x2 - xr*(1.d0/REAL(nr-1,8))
f(3) = x2 - xr*(2.d0/REAL(nr-1,8))
qpsi(1) = (qpsi(3)-qpsi(2))*(f(1)-f(2))/(f(3)-f(2)) + qpsi(2)
!---Extrapolate to LCFS
IF(.NOT.do_truncate)THEN
  DO i=1,nr
    pt=[rout(i),zout(i),0.d0]
    pt_last(1:2)=pt(1:2)-gseq%o_point
    pt_last=pt_last/SQRT(SUM(pt_last(1:2)**2))
    CALL gs_psi2pt(gseq,x1,pt,gseq%o_point,pt_last)
    rout(i)=pt(1)
    zout(i)=pt(2)
  END DO
END IF
!---Sample flux grid
rdim = rbounds(2)-rbounds(1)
zdim = zbounds(2)-zbounds(1)
DO i=1,nr
  cell=0
  pt(1) = (i-1)*rdim/REAL(nr-1,8) + rbounds(1)
  DO j=1,nz
    pt(2) = (j-1)*zdim/REAL(nz-1,8) + zbounds(1)
    call bmesh_findcell(device%mesh,cell,pt,f)
    call psi_int%interp(cell,f,gop,psi_tmp)
    psirz(i,j)=psi_tmp(1)
  END DO
END DO
CALL psi_int%delete()
!------------------------------------------------------------------------------
! Create output file
!------------------------------------------------------------------------------
WRITE(eqdsk_case,'(A,1X,A)')'tMaker:',run_info
rleft = rbounds(1)
zmid = (zbounds(2)+zbounds(1))/2.d0
IF(PRESENT(rcentr_in))THEN
  rcentr = rcentr_in
ELSE
  rcentr = (MAXVAL(rout)+MINVAL(rout))/2.d0 !raxis
END IF
bcentr = gseq%I%f_offset/rcentr
CALL gs_itor_nl(gseq,itor)
itor = itor/mu0
idum = 0 ! dummy variable
xdum = 0.d0 ! dummy variable
! Read or set limiting contour
IF(TRIM(limiter_file)=='')THEN
  IF(device%lim_nloops>1)THEN
    IF(nlim/=device%nlim_con+1)CALL oft_warn("Multiply-connected plasma region detected: Using largest boundary loop as limiter")
    nlim=0
    DO i=1,device%lim_nloops
      IF(device%lim_ptr(i+1)-device%lim_ptr(i)>nlim)THEN
        nlim=device%lim_ptr(i+1)-device%lim_ptr(i)
        lim_max=i
      END IF
    END DO
  ELSE
    lim_max=1
  END IF
  nlim=device%lim_ptr(lim_max+1)-device%lim_ptr(lim_max)+1
  ALLOCATE(rlim(nlim),zlim(nlim))
  DO i=device%lim_ptr(lim_max),device%lim_ptr(lim_max+1)-1
    rlim(i-device%lim_ptr(lim_max)+1)=device%mesh%r(1,device%lim_con(i))
    zlim(i-device%lim_ptr(lim_max)+1)=device%mesh%r(2,device%lim_con(i))
  END DO
  rlim(nlim)=rlim(1)
  zlim(nlim)=zlim(1)
  ! rlim=[rbounds(1),rbounds(1),rbounds(2),rbounds(2),rbounds(1)]
  ! zlim=[zbounds(1),zbounds(2),zbounds(2),zbounds(1),zbounds(1)]
ELSE
  WRITE(*,*)'  Limiter file: "',TRIM(limiter_file),'"'
  OPEN(NEWUNIT=io_unit,FILE=TRIM(limiter_file))
  READ(io_unit,*)nlim
  ALLOCATE(rlim(nlim),zlim(nlim))
  DO i=1,nlim
    READ(io_unit,*)rlim(i),zlim(i)
  END DO
  CLOSE(io_unit)
END IF
! COCOS transform
IF(PRESENT(cocos))THEN
  IF(cocos == 2)THEN
    WRITE(*,*) 'Using COCOS=2...'
    ffprim = -ffprim
    pprime = -pprime
    psirz = -psirz
    ! fpol = -fpol
    ! bcentr = -bcentr
    ! itor = -itor
    x1 = -x1
    x2 = -x2
  ELSE IF(cocos /= 7)THEN
    CALL oft_abort('Invalid COCOS version.','gs_save_eqdsk',__FILE__)
  END IF
END IF
! Write out gEQDSK file
2000 format(a48,3i4)
2020 format(5e16.9)
2022 format(2i5)
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
WRITE (io_unit,2000) eqdsk_case,0,nr,nz
WRITE (io_unit,2020) REAL([rdim,zdim,rcentr,rleft,zmid],4)
WRITE (io_unit,2020) REAL([raxis,zaxis,x2,x1,bcentr],4)
WRITE (io_unit,2020) REAL([itor,x2,xdum,raxis,xdum],4)
WRITE (io_unit,2020) REAL([zaxis,xdum,x1,xdum,xdum],4)
WRITE (io_unit,2020) (REAL(fpol(i),4),i=1,nr)
WRITE (io_unit,2020) (REAL(pres(i),4),i=1,nr)
WRITE (io_unit,2020) (REAL(ffprim(i),4),i=1,nr)
WRITE (io_unit,2020) (REAL(pprime(i),4),i=1,nr)
WRITE (io_unit,2020) ((REAL(psirz(i,j),4),i=1,nr),j=1,nz)
WRITE (io_unit,2020) (REAL(qpsi(i),4),i=1,nr)
WRITE (io_unit,2022) nr,nlim
WRITE (io_unit,2020) (REAL([rout(i),zout(i)],4),i=1,nr)
WRITE (io_unit,2020) (REAL([rlim(i),zlim(i)],4),i=1,nlim)
CLOSE (io_unit)
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,2ES11.3)')oft_indent,'Psi  = ',x1,x2
  WRITE(*,'(2A,ES11.3)')oft_indent,'Qmin = ',MINVAL(qpsi)
  WRITE(*,'(2A,ES11.3)')oft_indent,'Qmax = ',MAXVAL(qpsi)
END IF
CALL oft_decrease_indent
!---
DEALLOCATE(rout,zout,rlim,zlim)
DEALLOCATE(fpol,pres,ffprim,pprime,qpsi,psirz)
end subroutine gs_save_eqdsk
!------------------------------------------------------------------------------
!> Evaluate terms in augmented tracing ODE for computing Sauter factors (see @ref sauter_fc)
!------------------------------------------------------------------------------
subroutine sauter_apply(self,cell,f,gop,val)
class(sauter_interp), intent(inout) :: self !< Interpolation object
integer(4), intent(in) :: cell !< Cell for interpolation
real(8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(8), intent(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
real(8), intent(out) :: val(:) !< Reconstructed field at f [8]
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),d2op(6),pt(3),grad(3),tmp
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
!---
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
ELSE
  bratio = MIN(1.d0,mod_B/self%bmax)
  val(8)=val(2)*((1.d0 - SQRT(1.d0 - bratio) * (1.d0 + bratio / 2.d0)) / bratio**2)
END IF
deallocate(j)
end subroutine sauter_apply
!------------------------------------------------------------------------------
!> Compute factors required for Sauter bootstrap formula
!------------------------------------------------------------------------------
subroutine sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(4), intent(in) :: nr !< Number of flux sample points
real(8), intent(in) :: psi_q(nr) !< Location of flux sample points
real(8), intent(out) :: fc(nr) !< Trapped particle fraction \f$ f_c \f$
real(8), intent(out) :: r_avgs(nr,3) !< Flux surface averaged radial coordinates \f$<R>\f$, \f$<1/R>\f$, \f$<a>\f$
real(8), intent(out) :: modb_avgs(nr,2) !< Flux surface averaged field strength \f$<|B|>\f$, \f$<|B|^2>\f$
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,fpol,qpsi,h,h2,hf,ftu,ftl
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,cell
type(sauter_interp), target :: field
type(gs_factory), pointer :: device
!---
device=>gseq%device
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
!   x1 = x1 + (x2-x1)*2.d-2
!   x2 = x2 + (x1-x2)*2.d-2
END IF
! IF(.NOT.gseq%free)x1 = x1 + (x2-x1)*2.d-2
psi_int%u=>gseq%psi
CALL psi_int%setup(device%fe_rep)
!---Find Rmax along Zaxis
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
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Axis Position:'
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
  CALL oft_decrease_indent
END IF
!---Trace
call set_tracer(1)
! !$omp parallel private(j,psi_surf,pt,ptout,fpol,qpsi,field) firstprivate(pt_last)
field%u=>gseq%psi
field%mag_axis=gseq%o_point
CALL field%setup(device%fe_rep)
active_tracer%neq=8
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
! !$omp do schedule(dynamic,1)
do j=1,nr
  !------------------------------------------------------------------------------
  ! Trace contour
  !------------------------------------------------------------------------------
  psi_surf=psi_q(j)*(x2-x1) + x1
  IF(gseq%diverted.AND.psi_q(j)<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  !
  pt=pt_last
  ! !$omp critical
  CALL gs_psi2r(gseq,psi_surf,pt,psi_int=psi_int)
  ! !$omp end critical
  IF(gseq%mode==0)THEN
    field%f_surf=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    field%f_surf=SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  field%bmax=0.d0
  field%stage_1=.TRUE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  field%stage_1=.FALSE.
  CALL tracinginv_fs(device%mesh,pt(1:2))
  pt_last=pt
  !---Skip point if trace fails
  if(active_tracer%status/=1)THEN
    WRITE(*,*)j,pt
    CALL oft_warn("sauter_fc: Trace did not complete")
    CYCLE
  end if
  !---Get flux variables
  r_avgs(j,1)=active_tracer%v(3)/active_tracer%v(2)
  r_avgs(j,2)=active_tracer%v(4)/active_tracer%v(2)
  r_avgs(j,3)=active_tracer%v(5)/active_tracer%v(2)
  modb_avgs(j,1)=active_tracer%v(6)/active_tracer%v(2)
  modb_avgs(j,2)=active_tracer%v(7)/active_tracer%v(2)
  !---Sauter formula
  ! Equation 4
  h = modb_avgs(j,1)/field%bmax
  h2 = modb_avgs(j,2)/field%bmax**2
  ftu = 1.d0 - h2 / (h**2) * (1.d0 - SQRT(1.d0 - h) * (1.d0 + 0.5d0 * h))
  ! Equation 7
  hf = active_tracer%v(8)/active_tracer%v(2)
  ftl = 1.d0 - h2 * hf
  ! Equation 18,19
  fc(j) = 1.d0 - (0.75d0 * ftu + 0.25d0 * ftl)
end do
CALL field%delete
CALL active_tracer%delete
DEALLOCATE(ptout)
! !$omp end parallel
CALL psi_int%delete()
end subroutine sauter_fc
END MODULE oft_gs_util
