!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_gs_util.F90
!
!> Grad-Shafranov utility subroutines for TokaMaker
!!
!! @authors Chris Hansen
!! @date September 2017
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE oft_gs_util
USE oft_base
USE spline_mod
USE oft_io, ONLY: hdf5_create_file, hdf5_create_group, hdf5_write, hdf5_read
USE oft_mesh_type, ONLY: smesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_blagrange, oft_blag_geval
USE oft_blag_operators, ONLY: oft_blag_project, oft_lag_brinterp, oft_lag_bginterp, &
  oft_blag_vproject
USE tracing_2d, ONLY: active_tracer, tracinginv_fs, set_tracer
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_eq, flux_func, gs_get_cond_source, gs_get_cond_weights, &
  gs_set_cond_weights, gs_estored, gs_dflux, gs_tflux, gs_helicity, gs_itor_nl, &
  gs_psimax, gs_test_bounds, gs_b_interp, oft_indent, gs_get_qprof, gsinv_interp, &
  gs_psi2r, oft_increase_indent, oft_decrease_indent, oft_indent, gs_psi2pt
USE oft_gs_profiles
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Need docs
!---------------------------------------------------------------------------
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
CLASS(gs_eq), POINTER :: gs_fit => NULL()
CLASS(flux_func), POINTER :: ff_fit => NULL()
REAL(8), POINTER, DIMENSION(:,:) :: tmpprof => NULL()
CONTAINS
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! @param[in] filename File storing function definition
!! @param[out] f Flux function object
!------------------------------------------------------------------------------
SUBROUTINE gs_profile_load(filename,F)
CHARACTER(LEN=*), INTENT(in) :: filename
CLASS(flux_func), POINTER, INTENT(out) :: F
!---
REAL(8) :: alpha,beta,sep
REAL(8), ALLOCATABLE, DIMENSION(:) :: cofs,yvals
INTEGER(4) :: ncofs,npsi,io_unit
LOGICAL :: zero_grad
CHARACTER(LEN=10) :: profType
!---
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
READ(io_unit,*)profType
!---
SELECT CASE(TRIM(profType))
  CASE("zero")
    ALLOCATE(zero_flux_func::F)
  CASE("flat")
    CALL create_flat_f(F)
  CASE("linear")
    READ(io_unit,*)ncofs
    READ(io_unit,*)alpha
    CALL create_linear_ff(F,alpha)
    F%ncofs=ncofs
  CASE("poly")
    READ(io_unit,*)ncofs,zero_grad
    ALLOCATE(cofs(ncofs))
    READ(io_unit,*)cofs
    CALL create_poly_ff(F,ncofs,cofs,zero_grad)
    DEALLOCATE(cofs)
  CASE("linterp")
    READ(io_unit,*)ncofs,alpha
    ALLOCATE(cofs(ncofs),yvals(ncofs))
    READ(io_unit,*)cofs
    READ(io_unit,*)yvals
    CALL create_linterp_ff(F,ncofs,cofs,yvals,alpha)
    DEALLOCATE(cofs,yvals)
  CASE("idcd")
    READ(io_unit,*)ncofs
    READ(io_unit,*)sep,alpha
    CALL create_twolam_ff(F,sep,alpha)
    F%ncofs=ncofs
  CASE("step-slant")
    READ(io_unit,*)ncofs
    READ(io_unit,*)sep,alpha,beta
    CALL create_stepslant_ff(F,sep,alpha,beta)
    F%ncofs=ncofs
  ! CASE("mercier")
  !   READ(io_unit,*)npsi
  !   CALL create_mercier_ff(F,npsi)
  !   F%ncofs=0
  CASE("wesson")
    READ(io_unit,*)ncofs
    READ(io_unit,*)beta
    CALL create_wesson_ff(F,ncofs,beta)
  CASE DEFAULT
    CALL oft_abort('Invalid profile type.','gs_profile_load',__FILE__)
END SELECT
CLOSE(io_unit)
END SUBROUTINE gs_profile_load
!------------------------------------------------------------------------------
!> Save flux function object to definition file
!!
!! @param[in] filename File to store function definition
!! @param[in] f Flux function object
!------------------------------------------------------------------------------
SUBROUTINE gs_profile_save(filename,F)
CHARACTER(LEN=*), INTENT(in) :: filename
CLASS(flux_func), POINTER, INTENT(in) :: F
!---
INTEGER(4) :: io_unit
!---
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
!---
SELECT TYPE(this=>F)
  TYPE IS(zero_flux_func)
    WRITE(io_unit,*)"zero"
  TYPE IS(flat_flux_func)
    WRITE(io_unit,*)"flat"
  TYPE IS(linear_flux_func)
    WRITE(io_unit,*)"linear"
    WRITE(io_unit,*)this%ncofs
    WRITE(io_unit,*)this%alpha
  TYPE IS(poly_flux_func)
    WRITE(io_unit,*)"poly"
    IF(this%zero_grad)THEN
      WRITE(io_unit,*)this%deg+1,this%zero_grad
      WRITE(io_unit,*)1.d0,this%cofs
    ELSE
      WRITE(io_unit,*)this%deg,this%zero_grad
      WRITE(io_unit,*)this%cofs
    END IF
  TYPE IS(linterp_flux_func)
    WRITE(io_unit,*)"linterp"
    WRITE(io_unit,*)this%npsi,this%y0
    WRITE(io_unit,*)this%x
    WRITE(io_unit,*)this%yp
  TYPE IS(twolam_flux_func)
    WRITE(io_unit,*)"idcd"
    WRITE(io_unit,*)this%ncofs
    WRITE(io_unit,*)this%sep,this%alpha
  TYPE IS(stepslant_flux_func)
    WRITE(io_unit,*)"step-slant"
    WRITE(io_unit,*)this%ncofs
    WRITE(io_unit,*)this%sep,this%alpha,this%beta
  ! TYPE IS(mercier_flux_func)
  !   WRITE(io_unit,*)"mercier"
  !   WRITE(io_unit,*)this%npsi
  TYPE IS(wesson_flux_func)
    WRITE(io_unit,*)"wesson"
    WRITE(io_unit,*)this%ncofs
    WRITE(io_unit,*)this%gamma
  CLASS DEFAULT
    CALL oft_abort('Invalid profile type.','gs_profile_save',__FILE__)
END SELECT
CLOSE(io_unit)
END SUBROUTINE gs_profile_save
!---------------------------------------------------------------------------
! SUBROUTINE gs_save
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in] filename Filename for restart file
!! @param[in] mpsi_sample Number of flux (radial) sampling points (optional)
!---------------------------------------------------------------------------
SUBROUTINE gs_save(self,filename,mpsi_sample)
class(gs_eq), target, intent(inout) :: self
character(LEN=*), intent(in) :: filename
INTEGER(4), OPTIONAL, intent(in) :: mpsi_sample
integer(4) :: i,j,m,np_plot
real(8) :: x1,x2,r
real(8), allocatable, dimension(:,:) :: tmpout
integer(4), POINTER :: lctmp(:,:)
real(8), POINTER :: rtmp(:,:),cond_corr(:,:,:)
REAL(8), PARAMETER :: rst_version = 2.d0
LOGICAL :: pm_save
CLASS(oft_vector), POINTER :: a,br,bt,bz
CLASS(oft_solver), POINTER :: solver
TYPE(gs_b_interp) :: field
real(r8), POINTER :: vals_tmp(:)
m=100
IF(PRESENT(mpsi_sample))m=mpsi_sample
!---
CALL hdf5_create_file(filename)
!---Save mesh
CALL hdf5_create_group(filename,'mesh')
CALL hdf5_write(smesh%np,filename,'mesh/np')
CALL hdf5_write(smesh%nc,filename,'mesh/nc')
CALL hdf5_write(smesh%r(1:2,:),filename,'mesh/r')
CALL hdf5_write(smesh%lc,filename,'mesh/lc')
CALL hdf5_write(smesh%reg,filename,'mesh/regions')
CALL hdf5_write(oft_blagrange%order,filename,'mesh/order')
CALL smesh%tessellate(rtmp,lctmp,oft_blagrange%order)
CALL hdf5_write(rtmp,filename,'mesh/r_plot')
CALL hdf5_write(lctmp,filename,'mesh/lc_plot')
np_plot=SIZE(rtmp,DIM=2)
DEALLOCATE(rtmp,lctmp)
!---Save GS components
CALL hdf5_create_group(filename,'gs')
CALL hdf5_write(rst_version,filename,'gs/version')
CALL hdf5_write(m,filename,'gs/mpsi')
IF(self%has_plasma)THEN
  CALL hdf5_write(0.d0,filename,'gs/vac_flag')
ELSE
  CALL hdf5_write(1.d0,filename,'gs/vac_flag')
END IF
CALL hdf5_write(self%psiscale,filename,'gs/psiscale')
CALL hdf5_write(self%alam,filename,'gs/alam')
CALL hdf5_write(self%pnorm,filename,'gs/pnorm')
CALL hdf5_write(self%R0_target,filename,'gs/r0_target')
CALL hdf5_write(self%Itor_target,filename,'gs/Ip_target')
CALL hdf5_write(self%vcontrol_val,filename,'gs/vcont_val')
CALL hdf5_write(self%plasma_bounds,filename,'gs/bounds')
CALL hdf5_write(self%o_point,filename,'gs/o_point')
CALL hdf5_write(self%nx_points,filename,'/gs/nx_points')
CALL hdf5_write(self%x_points,filename,'/gs/x_points')
IF(self%isoflux_ntargets>0)THEN
  CALL hdf5_write(self%isoflux_ntargets,filename,'gs/isoflux_ntargets')
  CALL hdf5_write(self%isoflux_targets,filename,'gs/isoflux_targets')
END IF
! CALL vector_cast(psiv,self%psi)
NULLIFY(vals_tmp)
CALL self%psi%get_local(vals_tmp)
CALL hdf5_write(vals_tmp,filename,'gs/psi')
!---Save B-field
CALL self%psi%new(a)
CALL self%psi%new(br)
CALL self%psi%new(bt)
CALL self%psi%new(bz)
! CALL vector_cast(psiv,a)
field%gs=>self
CALL field%setup()
CALL create_cg_solver(solver)
solver%A=>self%mop
solver%its=-2
CALL create_diag_pre(solver%pre)
!---Project to plotting grid
pm_save=oft_env%pm; oft_env%pm=.FALSE.
ALLOCATE(tmpout(3,a%n))
CALL oft_blag_vproject(field,br,bt,bz)
CALL a%set(0.d0)
CALL solver%apply(a,br)
CALL a%get_local(vals_tmp)
tmpout(1,:)=vals_tmp
CALL a%set(0.d0)
CALL solver%apply(a,bt)
CALL a%get_local(vals_tmp)
tmpout(2,:)=vals_tmp
CALL a%set(0.d0)
CALL solver%apply(a,bz)
CALL a%get_local(vals_tmp)
tmpout(3,:)=vals_tmp
oft_env%pm=pm_save
CALL hdf5_write(tmpout,filename,'gs/B')
CALL a%delete
CALL br%delete
CALL bt%delete
CALL bz%delete
DEALLOCATE(a,br,bt,bz,tmpout,vals_tmp)
CALL field%delete()
CALL solver%pre%delete()
DEALLOCATE(solver%pre)
CALL solver%delete()
DEALLOCATE(solver)
!---Save coil/conductor info
ALLOCATE(tmpout(smesh%nc,1))
CALL gs_get_cond_source(self,tmpout(:,1))
CALL hdf5_write(tmpout(:,1),filename,'gs/cond_source')
DEALLOCATE(tmpout)
IF(self%ncoil_regs>0)THEN
  ALLOCATE(tmpout(self%ncoil_regs,1))
  tmpout=0.d0
  DO i=1,self%ncoil_regs
    DO j=1,self%ncoils
      tmpout(i,1)=tmpout(i,1) &
        + self%coil_currs(j)*self%coil_nturns(self%coil_regions(i)%id,j)
    END DO
  END DO
  CALL hdf5_write(tmpout(:,1),filename,'gs/int_coils')
  DEALLOCATE(tmpout)
END IF
IF(self%ncoils_ext>0)THEN
  ALLOCATE(tmpout(self%ncoils_ext,1))
  tmpout=0.d0
  DO i=1,self%ncoils_ext
    DO j=1,self%ncoils
      tmpout(i,1)=tmpout(i,1) &
        + self%coil_currs(j)*self%coil_nturns(smesh%nreg+i,j)
    END DO
  END DO
  CALL hdf5_write(tmpout(:,1),filename,'gs/ext_coils')
  DEALLOCATE(tmpout)
END IF
!---Get plasma bounds
x1=0.d0; x2=1.d0
IF(self%plasma_bounds(1)>-1.d98)THEN
  x1=self%plasma_bounds(1); x2=self%plasma_bounds(2)
END IF
ALLOCATE(tmpout(2,m))
!---Save F profile
CALL hdf5_create_group(filename,'gs/f')
CALL hdf5_write(self%I%f_offset,filename,'gs/f/offset')
CALL hdf5_write(self%mode,filename,'gs/f/mode')
DO i=1,m
  r=(i-1)*(x2-x1)/(m-1) + x1
  IF(i==1)r = r + (x2-x1)*1.d-10
  IF(i==m)r = r - (x2-x1)*1.d-10
  tmpout(1,i)=self%I%fp(r)
  tmpout(2,i)=self%I%f(r)
END DO
CALL hdf5_write(tmpout,filename,'gs/f/sample')
!---Save P profile
CALL hdf5_create_group(filename,'gs/p')
DO i=1,m
  r=(i-1)*(x2-x1)/(m-1) + x1
  IF(i==1)r = r + (x2-x1)*1.d-10
  IF(i==m)r = r - (x2-x1)*1.d-10
  tmpout(1,i)=self%P%fp(r)
  tmpout(2,i)=self%P%f(r)
END DO
CALL hdf5_write(tmpout,filename,'gs/p/sample')
!---
NULLIFY(cond_corr)
DO i=1,self%ncond_regs
  IF(ASSOCIATED(self%cond_regions(i)%corr_3d))THEN
    IF(.NOT.ASSOCIATED(cond_corr))THEN
      ALLOCATE(cond_corr(3,SIZE(self%cond_regions(i)%corr_3d,DIM=2),np_plot))
      cond_corr=0.d0
    END IF
    DO j=1,self%cond_regions(i)%neigs
      cond_corr=cond_corr+self%cond_regions(i)%corr_3d(:,:,:,j)*self%cond_regions(i)%weights(j)
    END DO
  END IF
END DO
IF(ASSOCIATED(cond_corr))THEN
  CALL hdf5_write(SQRT(SUM(cond_corr**2,DIM=2)/REAL(SIZE(cond_corr,DIM=2),8)),filename,'gs/cond_brms')
  DEALLOCATE(cond_corr)
END IF
!---
DEALLOCATE(tmpout)
END SUBROUTINE gs_save
!---------------------------------------------------------------------------
! SUBROUTINE gs_load
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in] filename Filename for restart file
!---------------------------------------------------------------------------
SUBROUTINE gs_load(self,filename)
class(gs_eq), target, intent(inout) :: self
character(LEN=*), intent(in) :: filename
integer(4) :: i,m
real(8) :: x1,x2,tmpval,tmp_version
real(8), allocatable, dimension(:,:) :: tmpin
REAL(8), PARAMETER :: rst_version = 1.d0
real(8), pointer :: vals_tmp(:)
!---
CALL hdf5_read(tmp_version,filename,'gs/version')
!---Load mesh
! CALL hdf5_read(smesh%np,filename,'mesh/np')
! CALL hdf5_read(smesh%nf,filename,'mesh/nc')
! CALL hdf5_read(smesh%r(1:2,:),filename,'mesh/r')
! CALL hdf5_read(smesh%lf,filename,'mesh/lc')
CALL hdf5_read(tmpval,filename,'mesh/order')
IF(INT(tmpval)/=oft_blagrange%order)CALL oft_abort("order mismatch","gs_load",__FILE__)
!---Load GS components
CALL hdf5_read(tmpval,filename,'gs/mpsi')
m=INT(tmpval)
CALL hdf5_read(self%psiscale,filename,'gs/psiscale')
CALL hdf5_read(self%alam,filename,'gs/alam')
CALL hdf5_read(self%pnorm,filename,'gs/pnorm')
CALL hdf5_read(self%R0_target,filename,'gs/r0_target')
CALL hdf5_read(self%Itor_target,filename,'gs/Ip_target')
CALL hdf5_read(self%vcontrol_val,filename,'gs/vcont_val')
CALL hdf5_read(self%plasma_bounds,filename,'gs/bounds')
! CALL vector_cast(psiv,self%psi)
NULLIFY(vals_tmp)
CALL self%psi%get_local(vals_tmp)
CALL hdf5_read(vals_tmp,filename,'gs/psi')
CALL self%psi%restore_local(vals_tmp)
IF(self%ncond_eigs>0)THEN
  ALLOCATE(tmpin(smesh%nc,1))
  CALL hdf5_read(tmpin(:,1),filename,'gs/cond_source')
  CALL cond_fit(self,tmpin)
  DEALLOCATE(tmpin)
END IF
IF(self%ncoil_regs>0)THEN
  CALL oft_abort("Not supported","gs_load",__FILE__)
  ! ALLOCATE(tmpin(self%ncoil_regs,1))
  ! CALL hdf5_read(tmpin(:,1),filename,'gs/int_coils')
  ! DO i=1,self%ncoil_regs
  !   self%coil_regions(i)%curr=tmpin(i,1)
  ! END DO
  ! DEALLOCATE(tmpin)
END IF
IF(self%ncoils_ext>0)THEN
  CALL oft_abort("Not supported","gs_load",__FILE__)
  ! ALLOCATE(tmpin(self%ncoils_ext,1))
  ! CALL hdf5_read(tmpin(:,1),filename,'gs/ext_coils')
  ! DO i=1,self%ncoils_ext
  !   self%coils_ext(i)%curr=tmpin(i,1)
  ! END DO
  ! DEALLOCATE(tmpin)
END IF
!---Get plasma bounds
x1=0.d0; x2=1.d0
IF(self%plasma_bounds(1)>-1.d98)THEN
  x1=self%plasma_bounds(1); x2=self%plasma_bounds(2)
END IF
ALLOCATE(tmpin(2,m))
!---Load F profile
CALL hdf5_read(self%I%f_offset,filename,'gs/f/offset')
CALL hdf5_read(tmpval,filename,'gs/f/mode')
IF(self%mode/=INT(tmpval))CALL oft_abort("Mode mismatch","gs_load",__FILE__)
self%mode=INT(tmpval)
CALL hdf5_read(tmpin,filename,'gs/f/sample')
CALL self%I%update(self)
CALL fit_ff(self%I,tmpin)
!---Load P profile
CALL hdf5_read(tmpin,filename,'gs/p/sample')
CALL self%P%update(self)
CALL fit_ff(self%P,tmpin)
!---
DEALLOCATE(tmpin)
END SUBROUTINE gs_load
!------------------------------------------------------------------------------
! SUBROUTINE cond_fit_error
!------------------------------------------------------------------------------
!> Needs docs
!!
!! @param[in,out] m Needs docs
!! @param[in,out] n Needs docs
!! @param[in,out] cofs Needs docs [n]
!! @param[in,out] err Needs docs [m]
!! @param[in,out] iflag Needs docs
!------------------------------------------------------------------------------
SUBROUTINE cond_fit_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
!---
CALL gs_set_cond_weights(gs_fit,cofs,.FALSE.)
CALL gs_get_cond_source(gs_fit,err)
err=err-tmpprof(:,1)
end subroutine cond_fit_error
!---------------------------------------------------------------------------
! SUBROUTINE cond_fit
!---------------------------------------------------------------------------
!> Needs docs
!!
!! @param[in,out] self Needs docs
!! @param[in,out] tmpin Needs docs
!---------------------------------------------------------------------------
subroutine cond_fit(self,tmpin)
class(gs_eq), target, intent(inout) :: self
real(8), target, intent(in) :: tmpin(:,:)
real(8), allocatable :: wttmp(:),contmp(:)
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs
integer(4), allocatable, dimension(:) :: ipvt
!---
gs_fit=>self
tmpprof=>tmpin
ALLOCATE(contmp(smesh%nc),wttmp(self%ncond_eigs))
CALL gs_get_cond_weights(gs_fit,wttmp,.TRUE.)
!---Use MINPACK to find maximum (zero gradient)
ncons=smesh%nc
ncofs=self%ncond_eigs
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
call lmdif(cond_fit_error,ncons,ncofs,wttmp,contmp, &
             ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
             nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!---
DEALLOCATE(contmp,wttmp)
end subroutine cond_fit
!------------------------------------------------------------------------------
! SUBROUTINE fit_ff_error
!------------------------------------------------------------------------------
!> Needs docs
!!
!! @param[in,out] m Needs docs
!! @param[in,out] n Needs docs
!! @param[in,out] cofs Needs docs [n]
!! @param[in,out] err Needs docs [m]
!! @param[in,out] iflag Needs docs
!------------------------------------------------------------------------------
SUBROUTINE fit_ff_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
integer(4) :: i,mtmp,ierr
real(8) :: r,x1,x2
!---
x1=0.d0; x2=1.d0
IF(ff_fit%plasma_bounds(1)>-1.d98)THEN
  x1=ff_fit%plasma_bounds(1); x2=ff_fit%plasma_bounds(2)
END IF
ierr=ff_fit%set_cofs(cofs)
!---
mtmp=m/2
DO i=1,mtmp
  r=(i-1)*(x2-x1)/(mtmp-1) + x1
  IF(i==1)r = r + (x2-x1)*1.d-10
  IF(i==mtmp)r = r - (x2-x1)*1.d-10
  err(i)=tmpprof(1,i)-ff_fit%fp(r)
  err(mtmp+i)=tmpprof(2,i)-ff_fit%f(r)
END DO
end subroutine fit_ff_error
!---------------------------------------------------------------------------
! SUBROUTINE fit_ff
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self Needs docs
!! @param[in,out] tmpin Needs docs
!---------------------------------------------------------------------------
subroutine fit_ff(self,tmpin)
class(flux_func), target, intent(inout) :: self
real(8), target, intent(in) :: tmpin(:,:)
integer(4) :: dims(2)
real(8), allocatable :: coftmp(:),contmp(:)
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs
integer(4), allocatable, dimension(:) :: ipvt
!---
ff_fit=>self
tmpprof=>tmpin
dims=SHAPE(tmpin)
ALLOCATE(contmp(2*dims(2)),coftmp(self%ncofs))
CALL self%get_cofs(coftmp)
!---Use MINPACK to find maximum (zero gradient)
ncons=2*dims(2)
ncofs=self%ncofs
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
call lmdif(fit_ff_error,ncons,ncofs,coftmp,contmp, &
             ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
             nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!---
DEALLOCATE(contmp,coftmp)
end subroutine fit_ff
!---------------------------------------------------------------------------
! SUBROUTINE gs_comp_globals
!---------------------------------------------------------------------------
!> Compute toroidal current for Grad-Shafranov equilibrium
!!
!! @param[in,out] self G-S object
!! @param[out] itor Toroidal current
!! @param[out] centroid Current centroid (optional) [2]
!---------------------------------------------------------------------------
subroutine gs_comp_globals(self,itor,centroid,vol,pvol,dflux,tflux,bp_vol)
class(gs_eq), intent(inout) :: self
real(8), intent(out) :: itor,centroid(2),vol,pvol,dflux,tflux,bp_vol
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
real(8) :: itor_loc,goptmp(3,3),v,psitmp(1),gpsitmp(3)
real(8) :: pt(3),curr_cent(2),Btor,Bpol(2)
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup
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
  do m=1,oft_blagrange%quad%np
    call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,oft_blagrange%quad%pts(:,m),goptmp,psitmp)
    IF(psitmp(1)<self%plasma_bounds(1))CYCLE
    pt=smesh%log2phys(i,oft_blagrange%quad%pts(:,m))
    !---Compute Magnetic Field
    IF(gs_test_bounds(self,pt))THEN
      IF(self%mode==0)THEN
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + (self%alam**2)*self%I%Fp(psitmp(1))*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+self%eps))
      ELSE
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + .5d0*self%alam*self%I%Fp(psitmp(1))/(pt(1)+self%eps))
      END IF
      itor = itor + itor_loc*v*oft_blagrange%quad%wts(m)
      centroid = centroid + itor_loc*pt(1:2)*v*oft_blagrange%quad%wts(m)
      pvol = pvol + (self%pnorm*self%P%F(psitmp(1)))*v*oft_blagrange%quad%wts(m)*pt(1)
      vol = vol + v*oft_blagrange%quad%wts(m)*pt(1)
      !---Compute total toroidal Field
      IF(self%mode==0)THEN
        Btor = (self%alam*(self%I%F(psitmp(1))) + self%I%f_offset)/(pt(1)+self%eps)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2))/(pt(1)+self%eps)
      END IF
      tflux = tflux + Btor*v*oft_blagrange%quad%wts(m)
      !---Compute internal inductance
      call psi_geval%interp(i,oft_blagrange%quad%pts(:,m),goptmp,gpsitmp)
      Bpol = [gpsitmp(1),gpsitmp(2)]/(pt(1)+self%eps)
      bp_vol = bp_vol + SUM(Bpol**2)*v*oft_blagrange%quad%wts(m)*pt(1)
      !---Compute differential toroidal Field
      IF(self%mode==0)THEN
        Btor = self%alam*(self%I%F(psitmp(1)))/pt(1)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2) &
        - self%I%f_offset)/pt(1)
      END IF
      dflux = dflux + Btor*v*oft_blagrange%quad%wts(m)
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
!---------------------------------------------------------------------------
!> Compute plasma loop voltage
!---------------------------------------------------------------------------
subroutine gs_calc_vloop(self,vloop)
class(gs_eq), intent(inout) :: self !< G-S object
real(8), intent(out) :: vloop !< loop voltage
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
real(8) :: itor_loc !< local toroidal current in integration
real(8) :: itor !< toroidal current
real(8) :: j_NI_loc !< local non-inductive current in integration
real(8) :: I_NI !< non-inductive F*F'
real(8) :: eta_jsq !< eta*j_NI**2 
real(8) :: goptmp(3,3) !< needs docs
real(8) :: v !< volume
real(8) :: pt(3) !< radial coordinate
real(8) :: curr_cent(2) !< needs docs
real(8) :: psitmp(1) !< magnetic flux coordinate
real(8) :: gpsitmp(3) !< needs docs
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup
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
  do m=1,oft_blagrange%quad%np
    call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,oft_blagrange%quad%pts(:,m),goptmp,psitmp)
    IF(psitmp(1)<self%plasma_bounds(1))CYCLE
    pt=smesh%log2phys(i,oft_blagrange%quad%pts(:,m))
    !---Compute toroidal current itor, and eta*j^2 eta_jsq (numerator of Vloop integral)
    IF(gs_test_bounds(self,pt))THEN
      IF(ASSOCIATED(self%I_NI))I_NI=self%I_NI%Fp(psitmp(1))
      IF(self%mode==0)THEN
        j_NI_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
          + (self%alam**2)*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+self%eps))
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
          + (self%alam**2)*self%I%Fp(psitmp(1))*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+self%eps))
      ELSE
        j_NI_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
          + (0.5d0*self%alam*self%I%Fp(psitmp(1)) - I_NI)/(pt(1)+self%eps))
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
          + .5d0*self%alam*self%I%Fp(psitmp(1))/(pt(1)+self%eps))
      END IF
      eta_jsq = eta_jsq + self%eta%fp(psitmp(1))*(j_NI_loc**2)*v*oft_blagrange%quad%wts(m)*pt(1)
      itor = itor + itor_loc*v*oft_blagrange%quad%wts(m)
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
!---------------------------------------------------------------------------
! SUBROUTINE gs_analyze
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
SUBROUTINE gs_analyze(self)
class(gs_eq), intent(inout) :: self
integer(4) :: i,io_unit
integer(4), parameter :: npsi = 50
real(8) :: Itor,centroid(2),vol,pvol,dflux,tflux,pmax,curr,bp_vol,li
real(8) :: psimax,baxis(2),prof(npsi),psi_q(npsi),dl,rbounds(2,2),zbounds(2,2)
real(8) :: beta(2),q95,tmp,psi0,psi1
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Equilibrium Statistics:'
CALL oft_increase_indent
IF(.NOT.self%has_plasma)THEN
  WRITE(*,'(2A)')oft_indent,'Vacuum only'
  RETURN
END IF
IF(self%mode==0)THEN
  WRITE(*,'(2A)')oft_indent,"Flux Spec               =   F * F'"
ELSE
  WRITE(*,'(2A)')oft_indent,"Flux Spec               =   (F^2)'"
END IF
IF(self%diverted)THEN
  WRITE(*,'(2A)')oft_indent,'Topology                =   Diverted'
ELSE
  WRITE(*,'(2A)')oft_indent,'Topology                =   Limited'
END IF
CALL gs_comp_globals(self,Itor,centroid,vol,pvol,dflux,tflux,bp_vol)
!---Get q-profile
psi0=0.02d0; psi1=0.98d0 !1.d0
! IF(self%plasma_bounds(1)>-1.d98)THEN
  ! psi0=self%plasma_bounds(1); psi1=self%plasma_bounds(2)
  ! psi0 = psi0 + (psi1-psi0)*2.d-2
  ! psi1 = psi1 + (psi0-psi1)*2.d-2
! ELSE
  ! IF(.NOT.self%free)psi0 = psi0 + (psi1-psi0)*2.d-2
! END IF
do i=1,npsi
  psi_q(i)=(psi1-psi0)*((i-1)/REAL(npsi-1,8)) + psi0
end do
CALL gs_get_qprof(self,npsi,psi_q,prof,dl,rbounds,zbounds)
IF(dl<=0.d0)CALL oft_abort('Tracing q-profile failed','gs_analyze',__FILE__)
q95=linterp(psi_q,prof,npsi,0.05d0)
!
baxis=self%o_point
psimax=1.d0
IF(self%plasma_bounds(1)>-1.d98)psimax=self%plasma_bounds(2)
pmax=self%pnorm*self%P%f(psimax)*self%psiscale*self%psiscale
beta(1) = (2.d0*pvol/vol)/(Itor/dl)**2
IF(ABS(self%I%f_offset)>0.d0)beta(2) = 2.d0*pvol/vol/(self%I%f_offset/centroid(1))**2
WRITE(*,'(2A,ES11.3)')oft_indent,'Toroidal Current [A]    = ',Itor/mu0
WRITE(*,'(2A,2F8.3)') oft_indent,'Current Centroid [m]    =',centroid
WRITE(*,'(2A,2F8.3)') oft_indent,'Magnetic Axis [m]       =',baxis
WRITE(*,'(2A,F8.3)')  oft_indent,'Elongation              =',(zbounds(2,2)-zbounds(2,1))/(rbounds(1,2)-rbounds(1,1))
IF(self%diverted)THEN
  IF(self%x_points(2,self%nx_points)<zbounds(2,1))THEN
    zbounds(:,1)=self%x_points(:,self%nx_points)
  ELSE IF(self%x_points(2,self%nx_points)>zbounds(2,2))THEN
    zbounds(:,2)=self%x_points(:,self%nx_points)
  END IF
END IF
tmp=((rbounds(1,2)+rbounds(1,1))/2.d0-(zbounds(1,2)+zbounds(1,1))/2.d0)*2.d0/(rbounds(1,2)-rbounds(1,1))
WRITE(*,'(2A,F8.3)')oft_indent,'Triangularity           =',tmp
WRITE(*,'(2A,F8.3)')oft_indent,'Plasma Volume [m^3]     =',vol*2.d0*pi
WRITE(*,'(2A,3F8.3)') oft_indent,'q_0, q_95, q_a          =',prof(1),q95,prof(npsi)
WRITE(*,'(2A,ES11.3)')oft_indent,'Peak Pressure [Pa]      = ',pmax/mu0
WRITE(*,'(2A,ES11.3)')oft_indent,'Stored Energy [J]       = ',pvol*2.d0*pi/mu0*3.d0/2.d0
WRITE(*,'(2A,F8.3)')  oft_indent,'<Beta_pol> [%]          = ',1.d2*beta(1)
IF(ABS(self%I%f_offset)>0.d0)THEN
  WRITE(*,'(2A,F8.3)')oft_indent,'<Beta_tor> [%]          = ',1.d2*beta(2)
END IF
WRITE(*,'(2A,ES11.3)')oft_indent,'Diamagnetic flux [Wb]   = ',dflux
WRITE(*,'(2A,ES11.3)')oft_indent,'Toroidal flux [Wb]      = ',tflux
WRITE(*,'(2A,ES11.3)')oft_indent,'li                      = ',(bp_vol/vol)/((Itor/dl)**2)
! IF(.NOT.self%free)THEN
!   CALL gs_helicity(self,dflux,Itor)
!   WRITE(*,'(2A,ES11.3)')oft_indent,'Magnetic Energy [J]     = ',dflux/(2.d0*mu0)
!   WRITE(*,'(2A,ES11.3)')oft_indent,'Magnetic Helicity [J-m] = ',Itor/(2.d0*mu0)
! END IF
! !---
! WRITE(*,*)
! WRITE(*,'(2A)')oft_indent,'Coil Currents [A]'
! DO i=1,self%ncoil_regs
!   curr = self%coil_regions(i)%curr & 
!     + self%coil_regions(i)%vcont_gain*self%vcontrol_val
!   WRITE(*,'(A,ES16.7)')oft_indent,curr*self%coil_regions(i)%area/mu0
! END DO
! WRITE(*,*)
CALL oft_decrease_indent
!---------------------------------------------------------------------------
! Create output file for q
!---------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE='safety_factor.dat')
WRITE(io_unit,'(A)')'# TokaMaker q-profile "Psi, q"'
DO i=1,npsi
  WRITE(io_unit,'(2ES11.3)')psi_q(i),prof(i)
END DO
CLOSE(io_unit)
END SUBROUTINE gs_analyze
!---------------------------------------------------------------------------
! SUBROUTINE gs_save_decon
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine gs_save_decon(gseq,npsi,ntheta,error_str)
class(gs_eq), intent(inout) :: gseq
integer(4), intent(in) :: npsi
integer(4), intent(in) :: ntheta
CHARACTER(LEN=OFT_ERROR_SLEN), OPTIONAL, INTENT(out) :: error_str
type(gsinv_interp), target :: field
type(oft_lag_brinterp) :: psi_int
real(8) :: gop(3,3),psi_surf(1),pt_last(3)
real(8) :: raxis,zaxis,f(3),pt(3),rmax,x1,x2,xr
real(8), allocatable :: ptout(:,:)
real(8), allocatable :: rout(:,:),zout(:,:),cout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell,io_unit
TYPE(spline_type) :: rz
!---
IF(PRESENT(error_str))error_str=""
WRITE(*,'(2A)')oft_indent,'Saving DCON file'
CALL oft_increase_indent
!---
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
xr = (x2-x1)
x1 = x1 + xr*1.d-3
x2 = x2 - xr*1.d-3
psi_int%u=>gseq%psi
CALL psi_int%setup()
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(smesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF( psi_surf(1) < x1)EXIT
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
ALLOCATE(cout(4,npsi))
ALLOCATE(rout(npsi,ntheta))
ALLOCATE(zout(npsi,ntheta))
!$omp parallel private(j,psi_surf,pt,ptout,field,rz,gop) firstprivate(pt_last)
field%u=>gseq%psi
CALL field%setup()
active_tracer%neq=3
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
!$omp do schedule(dynamic,1)
do j=1,npsi-1
  IF(PRESENT(error_str))THEN
    IF(error_str/="")CYCLE
  END IF
  !---------------------------------------------------------------------------
  ! Trace contour
  !---------------------------------------------------------------------------
  psi_surf(1)=(x2-x1)*(1.d0-j/REAL(npsi,4))**2
  psi_surf(1)=x2 - psi_surf(1)
  IF(gseq%diverted.AND.(psi_surf(1)-x1)/(x2-x1)<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  pt=pt_last
  !$omp critical
  CALL gs_psi2r(gseq,psi_surf(1),pt)
  !$omp end critical
  CALL tracinginv_fs(pt,ptout)
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
    rout(j,k+1)=rz%f(1)
    zout(j,k+1)=rz%f(2)
  END DO
  !---Destroy Spline
  CALL spline_dealloc(rz)
  !---------------------------------------------------------------------------
  ! Save DCON information
  !---------------------------------------------------------------------------
  cout(1,j)=psi_surf(1) ! Poloidal flux
  !---Toroidal flux function
  IF(gseq%mode==0)THEN
    cout(2,j)=gseq%alam*gseq%I%f(psi_surf(1))+gseq%I%f_offset
  ELSE
    cout(2,j)=SQRT(gseq%alam*gseq%I%f(psi_surf(1)) + gseq%I%f_offset**2) &
    + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  cout(3,j)=gseq%pnorm*gseq%P%f(psi_surf(1))/mu0 ! Plasma pressure
  cout(4,j)=cout(2,j)*active_tracer%v(3)/(2*pi) ! Safety Factor (q)
end do
CALL active_tracer%delete
DEALLOCATE(ptout)
!$omp end parallel
IF(PRESENT(error_str))THEN
  IF(error_str/="")THEN
    DEALLOCATE(cout,rout,zout)
    RETURN
  END IF
END IF
!---Information for O-point
rout(npsi,:)=raxis
zout(npsi,:)=zaxis
cout(1,npsi)=x2
IF(gseq%mode==0)THEN
  cout(2,npsi)=(gseq%alam*gseq%I%f(x2)+gseq%I%f_offset)
ELSE
  cout(2,npsi)=SQRT(gseq%alam*gseq%I%f(x2) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
END IF
cout(3,npsi)=gseq%pnorm*gseq%P%f(x2)/mu0
cout(4,npsi)=(cout(4,npsi-2)-cout(4,npsi-1))*(x2-cout(1,npsi-1))/(cout(1,npsi-2)-cout(1,npsi-1)) + cout(4,npsi-1)
!---------------------------------------------------------------------------
! Create output file
!---------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE='Psitri.dci',FORM='UNFORMATTED')
!---------------------------------------------------------------------------
! Write array lengths
!---------------------------------------------------------------------------
WRITE(io_unit)INT(npsi-1,4),INT(ntheta-1,4)
!---------------------------------------------------------------------------
! Write out flux surface quantities
!
! cout(1,:) -> psi(0:mpsi)
! cout(2,:) -> f(0:mpsi)
! cout(3,:) -> p(0:mpsi)
! cout(4,:) -> q(0:mpsi)
!---------------------------------------------------------------------------
DO j=1,4
  WRITE(io_unit)REAL(cout(j,:),4)
END DO
!---------------------------------------------------------------------------
! Write out inverse representation
!
! rout -> r(0:mpsi,0:mtheta)
! zout -> z(0:mpsi,0:mtheta)
!---------------------------------------------------------------------------
WRITE(io_unit)REAL(rout,4)
WRITE(io_unit)REAL(zout,4)
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
end subroutine gs_save_decon
!---------------------------------------------------------------------------
!> Save equilibrium to General Atomics gEQDSK file
!---------------------------------------------------------------------------
subroutine gs_save_eqdsk(gseq,filename,nr,nz,rbounds,zbounds,run_info,limiter_file,psi_pad,rcentr_in,trunc_eq,error_str)
class(gs_eq), intent(inout) :: gseq !< Equilibrium to save
CHARACTER(LEN=OFT_PATH_SLEN), intent(in) :: filename 
integer(4), intent(in) :: nr !< Number of radial points for flux/psi grid
integer(4), intent(in) :: nz !< Number of vertical points for flux grid
real(8), intent(in) :: rbounds(2) !< Radial extents for flux grid
real(8), intent(in) :: zbounds(2) !< Radial extents for flux grid
CHARACTER(LEN=40), intent(in) :: run_info !< Run information string [40]
CHARACTER(LEN=OFT_PATH_SLEN), intent(in) :: limiter_file !< Path to limiter file
REAL(8), intent(in) :: psi_pad !< Padding at LCFS in normalized units
REAL(8), optional, intent(in) :: rcentr_in !< Value to use for RCENTR (otherwise geometric center is used)
LOGICAL, OPTIONAL, INTENT(in) :: trunc_eq !< Truncate equilibrium at psi_pad
CHARACTER(LEN=OFT_ERROR_SLEN), OPTIONAL, INTENT(out) :: error_str
!
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,xr,psi_trace
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:),rout(:),zout(:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,k,cell,io_unit,lim_max
type(gsinv_interp), target :: field
TYPE(spline_type) :: rz
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
IF(PRESENT(error_str))error_str=""
WRITE(*,'(3A)')oft_indent,'Saving EQDSK file: ',TRIM(filename)
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
CALL psi_int%setup()
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(smesh,cell,pt,f)
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
ALLOCATE(rout(nr))
ALLOCATE(zout(nr))
!$omp parallel private(j,psi_surf,psi_trace,pt,ptout,field,fptmp) firstprivate(pt_last)
field%u=>gseq%psi
CALL field%setup()
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
  !---------------------------------------------------------------------------
  ! Trace contour
  !---------------------------------------------------------------------------
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
    CALL gs_psi2r(gseq,psi_trace,pt)
    !$omp end critical
    IF(j==nr)THEN
      ALLOCATE(ptout(3,active_tracer%maxsteps+1))
      CALL tracinginv_fs(pt(1:2),ptout)
    ELSE
      CALL tracinginv_fs(pt(1:2))
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
    !---------------------------------------------------------------------------
    ! Perform Cubic Spline Interpolation
    !---------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------
  ! Compute Mercier Profiles
  !---------------------------------------------------------------------------
  !---Get flux variables
  IF(gseq%mode==0)THEN
    fptmp=gseq%alam*gseq%I%f(psi_trace)+gseq%I%f_offset
    fpol(j)=gseq%alam*gseq%I%f(psi_surf)+gseq%I%f_offset
    ffprim(j)=(gseq%alam**2)*gseq%I%fp(psi_surf)*(gseq%I%f(psi_surf)+gseq%I%f_offset/gseq%alam)
  ELSE
    fptmp=SQRT(gseq%alam*gseq%I%f(psi_trace) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
    fpol(j)=SQRT(gseq%alam*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
    ffprim(j)=0.5d0*gseq%alam*gseq%I%fp(psi_surf)
  END IF
  pres(j)=gseq%pnorm*gseq%P%f(psi_surf)/mu0
  pprime(j)=gseq%pnorm*gseq%P%fp(psi_surf)/mu0
  !---Safety Factor (q)
  IF(j>1)qpsi(j)=fptmp*active_tracer%v(3)/(2*pi)
end do
CALL active_tracer%delete
!$omp end parallel
IF(PRESENT(error_str))THEN
  IF(error_str/="")THEN
    DEALLOCATE(rout,zout)
    DEALLOCATE(fpol,pres,ffprim,pprime,qpsi,psirz)
    RETURN
  END IF
END IF
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
    call bmesh_findcell(smesh,cell,pt,f)
    call psi_int%interp(cell,f,gop,psi_tmp)
    psirz(i,j)=psi_tmp(1)
  END DO
END DO
!---------------------------------------------------------------------------
! Create output file
!---------------------------------------------------------------------------
WRITE(eqdsk_case,'(A,X,A)')'tMaker:',run_info
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
  IF(gseq%lim_nloops>1)THEN
    IF(nlim/=gseq%nlim_con+1)CALL oft_warn("Multiply-connected plasma region detected: Using largest boundary loop as limiter")
    nlim=0
    DO i=1,gseq%lim_nloops
      IF(gseq%lim_ptr(i+1)-gseq%lim_ptr(i)>nlim)THEN
        nlim=gseq%lim_ptr(i+1)-gseq%lim_ptr(i)
        lim_max=i
      END IF
    END DO
  ELSE
    lim_max=1
  END IF
  nlim=gseq%lim_ptr(lim_max+1)-gseq%lim_ptr(lim_max)+1
  ALLOCATE(rlim(nlim),zlim(nlim))
  DO i=gseq%lim_ptr(lim_max),gseq%lim_ptr(lim_max+1)-1
    rlim(i-gseq%lim_ptr(lim_max)+1)=smesh%r(1,gseq%lim_con(i))
    zlim(i-gseq%lim_ptr(lim_max)+1)=smesh%r(2,gseq%lim_con(i))
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
! Write out gEQDSK file
2000 format(a48,3i4)
2020 format(5e16.9)
2022 format(2i5)
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
WRITE (io_unit,2000) eqdsk_case,0,nr,nz
WRITE (io_unit,2020) rdim,zdim,rcentr,rleft,zmid
WRITE (io_unit,2020) raxis,zaxis,x2,x1,bcentr
WRITE (io_unit,2020) itor,x2,xdum,raxis,xdum
WRITE (io_unit,2020) zaxis,xdum,x1,xdum,xdum
WRITE (io_unit,2020) (fpol(i),i=1,nr)
WRITE (io_unit,2020) (pres(i),i=1,nr)
WRITE (io_unit,2020) (ffprim(i),i=1,nr)
WRITE (io_unit,2020) (pprime(i),i=1,nr)
WRITE (io_unit,2020) ((psirz(i,j),i=1,nr),j=1,nz)
WRITE (io_unit,2020) (qpsi(i),i=1,nr)
WRITE (io_unit,2022) nr,nlim
WRITE (io_unit,2020) (rout(i),zout(i),i=1,nr)
WRITE (io_unit,2020) (rlim(i),zlim(i),i=1,nlim)
CLOSE (io_unit)
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,2ES11.3)')oft_indent,'Psi  = ',x1,x2
  WRITE(*,'(2A,ES11.3)')oft_indent,'Qmin = ',MINVAL(qpsi)
  WRITE(*,'(2A,ES11.3)')oft_indent,'Qmax = ',MAXVAL(qpsi)
  ! WRITE(*,'(2A)')oft_indent,'Done'
END IF
CALL oft_decrease_indent
!---
DEALLOCATE(rout,zout,rlim,zlim)
DEALLOCATE(fpol,pres,ffprim,pprime,qpsi,psirz)
end subroutine gs_save_eqdsk
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine sauter_apply(self,cell,f,gop,val)
class(sauter_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),d2op(6),pt(3),grad(3),tmp
real(8) :: s,c,Bp2,Bt2,mod_B,bratio
!---Get dofs
allocate(j(oft_blagrange%nce))
call oft_blagrange%ncdofs(cell,j)
!---Reconstruct gradient
grad=0.d0
do jc=1,oft_blagrange%nce
  call oft_blag_geval(oft_blagrange,cell,jc,f,rop,gop)
  grad=grad+self%uvals(j(jc))*rop
end do
!---Get radial position
pt=smesh%log2phys(cell,f)
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
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine sauter_fc(gseq,nr,psi_q,fc,r_avgs,modb_avgs)
class(gs_eq), intent(inout) :: gseq
integer(4), intent(in) :: nr
real(8), intent(in) :: psi_q(nr)
real(8), intent(out) :: fc(nr),r_avgs(nr,3),modb_avgs(nr,2)
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,fpol,qpsi,h,h2,hf,ftu,ftl
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,cell
type(sauter_interp), target :: field
!---
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
CALL psi_int%setup()
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(smesh,cell,pt,f)
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
CALL field%setup()
active_tracer%neq=8
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
! !$omp do schedule(dynamic,1)
do j=1,nr
  !---------------------------------------------------------------------------
  ! Trace contour
  !---------------------------------------------------------------------------
  psi_surf=psi_q(j)*(x2-x1) + x1
  IF(gseq%diverted.AND.psi_q(j)<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  !
  pt=pt_last
  ! !$omp critical
  CALL gs_psi2r(gseq,psi_surf,pt)
  ! !$omp end critical
  IF(gseq%mode==0)THEN
    field%f_surf=gseq%alam*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    field%f_surf=SQRT(gseq%alam*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
      + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  field%bmax=0.d0
  field%stage_1=.TRUE.
  CALL tracinginv_fs(pt(1:2))
  field%stage_1=.FALSE.
  CALL tracinginv_fs(pt(1:2))
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
end subroutine sauter_fc
END MODULE oft_gs_util
