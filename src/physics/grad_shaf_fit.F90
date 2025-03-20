!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_gs_fit.F90
!
!> GS Fitting implementation
!!
!! @authors Chris Hansen
!! @date March 2014
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE oft_gs_fit
use oft_base
USE oft_io, ONLY: hdf5_read, oft_file_exist
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell
!---
use oft_la_base, only: oft_vector
use oft_lu, only: lapack_matinv
!---
use oft_lag_basis, only: oft_blag_d2eval, oft_blag_geval
use oft_blag_operators, only: oft_lag_brinterp, oft_lag_bginterp
use oft_gs, only: gs_eq, gs_dflux, gs_get_cond_weights, gs_itor_nl, gs_psi2r, &
  gs_psimax, gs_set_cond_weights, gs_err_reason, gs_test_bounds, gs_get_cond_scales, &
  gs_get_qprof, gs_epsilon, gsinv_interp, oft_indent, oft_decrease_indent, oft_increase_indent
use oft_gs_profiles, only: twolam_flux_func
use tracing_2d, only: active_tracer, tracer, set_tracer
use mhd_utils, only: mu0
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
! TYPE fit_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE :: fit_constraint
  INTEGER(4) :: ncomp = 0
  REAL(8) :: wt = 1.d0
  REAL(8) :: val = 0.d0
  REAL(8), POINTER, DIMENSION(:) :: nax_corr => NULL()
  REAL(8), POINTER, DIMENSION(:,:) :: comp_r => NULL()
  REAL(8), POINTER, DIMENSION(:,:) :: comp_n => NULL()
CONTAINS
  PROCEDURE :: error => fit_dummy_error
  PROCEDURE :: eval => fit_dummy_eval
  PROCEDURE :: setup_comp => fit_dummy_setup_comp
  PROCEDURE :: get_nax => fit_dummy_nax_corr
  PROCEDURE :: is_parallel => fit_dummy_parallel
END TYPE fit_constraint
!---------------------------------------------------------------------------------
! TYPE coil_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: coil_constraint
  INTEGER(4) :: coil = 0
CONTAINS
  PROCEDURE :: error => fit_coil_error
  PROCEDURE :: eval => fit_coil_eval
END TYPE coil_constraint
!---------------------------------------------------------------------------------
! TYPE vcont_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: vcont_constraint
  INTEGER(4) :: coil = 0
CONTAINS
  PROCEDURE :: error => fit_vcont_error
  PROCEDURE :: eval => fit_vcont_eval
END TYPE vcont_constraint
!---------------------------------------------------------------------------------
! TYPE field_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: field_constraint
  INTEGER(4) :: cell = 0
  REAL(8) :: phi = 0.d0
  REAL(8) :: f(3) = 0.d0
  REAL(8) :: r(3) = [0.d0,0.d0,0.d0]
  REAL(8) :: v(3) = [0.d0,0.d0,0.d0]
CONTAINS
  PROCEDURE :: error => fit_field_error
  PROCEDURE :: eval => fit_field_eval
  PROCEDURE :: setup_comp => fit_field_setup_comp
END TYPE field_constraint
!---------------------------------------------------------------------------------
! TYPE flux_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: flux_constraint
  INTEGER(4) :: cell = 0
  REAL(8) :: f(3) = 0.d0
  REAL(8) :: r(3) = [0.d0,0.d0,0.d0]
CONTAINS
  PROCEDURE :: error => fit_flux_error
  PROCEDURE :: eval => fit_flux_eval
END TYPE flux_constraint
!---------------------------------------------------------------------------------
! TYPE saddle_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: saddle_constraint
  INTEGER(4) :: cells(2) = 0
  INTEGER(4) :: nr = 10
  INTEGER(4) :: nt = 10
  REAL(8) :: f(3,2) = 0.d0
  REAL(8) :: r(3,2) = 0.d0
  REAL(8) :: width = 0.d0
CONTAINS
  PROCEDURE :: error => fit_saddle_error
  PROCEDURE :: eval => fit_saddle_eval
  PROCEDURE :: setup_comp => fit_saddle_setup_comp
END TYPE saddle_constraint
!---------------------------------------------------------------------------------
! TYPE itor_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: itor_constraint
CONTAINS
  PROCEDURE :: error => fit_itor_error
  PROCEDURE :: eval => fit_itor_eval
END TYPE itor_constraint
!---------------------------------------------------------------------------------
! TYPE itor_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: dflux_constraint
CONTAINS
  PROCEDURE :: error => fit_dflux_error
  PROCEDURE :: eval => fit_dflux_eval
END TYPE dflux_constraint
!---------------------------------------------------------------------------------
! TYPE press_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: press_constraint
  INTEGER(4) :: cell = 0
  REAL(8) :: f(3) = 0.d0
  REAL(8) :: r(3) = [0.d0,0.d0,0.d0]
CONTAINS
  PROCEDURE :: error => fit_press_error
  PROCEDURE :: eval => fit_press_eval
END TYPE press_constraint
! !---------------------------------------------------------------------------------
! ! TYPE injlam_constraint
! !---------------------------------------------------------------------------------
! !> Needs Docs
! !---------------------------------------------------------------------------------
! TYPE, EXTENDS(fit_constraint) :: injlam_constraint
! CONTAINS
!   PROCEDURE :: error => fit_injlam_error
!   PROCEDURE :: eval => fit_injlam_eval
! END TYPE injlam_constraint
!---------------------------------------------------------------------------------
! TYPE qmin_constraint
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE, EXTENDS(fit_constraint) :: q_constraint
  INTEGER(4) :: type = 1 !< Type of constraint (1 -> min, 2 -> max, 3 -> rel_flux)
  REAL(8) :: loc = 0.d0 !< Location of constraint in normalized flux
CONTAINS
  PROCEDURE :: error => fit_q_error
  PROCEDURE :: eval => fit_q_eval
  PROCEDURE :: is_parallel => fit_q_parallel
END TYPE q_constraint
!---------------------------------------------------------------------------------
! TYPE fit_constraint_ptr
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
TYPE :: fit_constraint_ptr
  CLASS(fit_constraint), POINTER :: con => NULL()
END TYPE fit_constraint_ptr
!---
TYPE(gs_eq), POINTER, PUBLIC :: gs_active => NULL()
REAL(8), ALLOCATABLE :: cofs_best(:)
REAL(8) :: chi_best = 1.d99
REAL(8) :: alam_best = 1.d99
REAL(8) :: pnorm_best = 1.d99
REAL(8) :: vcont_best = 1.d99
CLASS(oft_vector), POINTER :: psi_best => NULL()
REAL(8), ALLOCATABLE :: cofs_scale(:)
REAL(8), ALLOCATABLE :: curr_in(:)
INTEGER(4), PRIVATE :: ncons = 0
INTEGER(4), PRIVATE :: ncofs = 0
INTEGER(4), PRIVATE :: ncond_active = 0
INTEGER(4), PRIVATE :: feval_count = 0
INTEGER(4), PRIVATE :: geval_count = 0
LOGICAL :: fit_pm = .FALSE.
LOGICAL, PRIVATE :: fit_I = .TRUE.
LOGICAL, PRIVATE :: fit_P = .TRUE.
LOGICAL, PRIVATE :: fit_Pnorm = .TRUE.
LOGICAL, PRIVATE :: fit_Alam = .FALSE.
LOGICAL, PRIVATE :: fit_R0 = .FALSE.
LOGICAL, PRIVATE :: fit_coils = .FALSE.
LOGICAL, PRIVATE :: fit_F0 = .FALSE.
LOGICAL, PRIVATE :: fit_V0 = .FALSE.
LOGICAL, PRIVATE :: fixed_centering = .FALSE.
LOGICAL, PRIVATE :: linearized_fit = .FALSE.
TYPE(fit_constraint_ptr), POINTER, DIMENSION(:) :: conlist => NULL()
!---
PUBLIC fit_constraint, field_constraint, itor_constraint, fit_pm
PUBLIC fit_constraint_ptr, fit_gs, fit_load!, create_driveinterp
CONTAINS
!---------------------------------------------------------------------------------
! SUBROUTINE fit_gs
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE fit_gs(gs,inpath,outpath,fitI,fitP,fitPnorm,fitAlam,fitR0,fitV0,fitCoils,fitF0,fixedCentering)
TYPE(gs_eq), TARGET, INTENT(inout) :: gs
CHARACTER(LEN=*), INTENT(in) :: inpath
CHARACTER(LEN=*), INTENT(in) :: outpath
LOGICAL, OPTIONAL, INTENT(in) :: fitI
LOGICAL, OPTIONAL, INTENT(in) :: fitP
LOGICAL, OPTIONAL, INTENT(in) :: fitPnorm
LOGICAL, OPTIONAL, INTENT(in) :: fitAlam
LOGICAL, OPTIONAL, INTENT(in) :: fitR0
LOGICAL, OPTIONAL, INTENT(in) :: fitV0
LOGICAL, OPTIONAL, INTENT(in) :: fitCoils
LOGICAL, OPTIONAL, INTENT(in) :: fitF0
LOGICAL, OPTIONAL, INTENT(in) :: fixedCentering
!---
real(8), allocatable :: fjac(:,:),qtf(:),error(:),cofs(:)
real(8), allocatable :: wa1(:),wa2(:),wa3(:),wa4(:)
real(8) :: ftol,xtol,gtol,epsfcn,factor
integer(4) :: maxfev,mode,nprint,info,nfev,njev,ldfjac,i,j,js,je,offset,ierr,io_unit
integer(4), allocatable :: ipvt(:)
logical :: file_exists,comp_var
NAMELIST/gs_fit_options/ftol,xtol,gtol,maxfev,epsfcn,factor,comp_var,linearized_fit
!---
IF(PRESENT(fitI))fit_I=fitI
IF(PRESENT(fitP))fit_P=fitP
IF(PRESENT(fitPnorm))fit_Pnorm=fitPnorm
IF(PRESENT(fitAlam))fit_alam=fitAlam
IF(PRESENT(fitR0))fit_R0=fitR0
IF(PRESENT(fitV0))fit_V0=fitV0
IF(PRESENT(fitCoils))fit_coils=fitCoils
IF(PRESENT(fitF0))fit_F0=fitF0
IF(PRESENT(fixedCentering))fixed_centering=fixedCentering
IF(fitPnorm.AND.fitR0)CALL oft_abort('R0 or Pnorm fitting cannot be used together', &
'fit_gs',__FILE__)
IF(fitPnorm)gs%R0_target=-1.d0
!---Load constraints
gs_active=>gs
WRITE(*,*)
WRITE(*,'(A)')'*** Loading fit constraints ***'
CALL fit_load(inpath,conlist)
!---Count coefficients
ncofs=0
IF(gs%free)THEN
  IF(fit_alam.OR.(gs_active%Itor_target>0.d0))ncofs = ncofs+1
ELSE
  ncofs=ncofs+1
  IF(fit_alam)CALL oft_abort('Lambda cannot be fit in fixed boundary mode.', &
  'fit_gs',__FILE__)
END IF
IF(gs_active%I%ncofs>0.AND.fit_I)THEN
  ncofs = ncofs+gs_active%I%ncofs
END IF
IF(fit_Pnorm.OR.fit_R0)ncofs = ncofs+1
IF(fit_V0)ncofs = ncofs + 1
IF(gs_active%P%ncofs>0.AND.fit_P)THEN
  ncofs = ncofs+gs_active%P%ncofs
END IF
ncond_active = 0
IF(gs_active%ncond_regs>0)THEN
  DO i=1,gs_active%ncond_regs
    IF(gs_active%cond_regions(i)%pair<0)THEN
      DO j=1,gs_active%cond_regions(i)%neigs
        IF(gs_active%cond_regions(i)%fixed(j))CYCLE
        ncond_active=ncond_active+1
      END DO
    END IF
  END DO
  WRITE(*,*)'Fixed',ncond_active,gs_active%ncond_eigs
  ncofs = ncofs + ncond_active !gs_active%ncond_eigs
END IF
IF(fit_coils)ncofs = ncofs + gs_active%ncoils
IF(fit_F0)ncofs = ncofs + 1
!---
ALLOCATE(cofs(ncofs),cofs_scale(ncofs))
cofs_scale=1.0
offset=0
IF(gs%free)THEN
  IF(fit_alam)THEN
    offset=1
    cofs(1)=gs_active%alam
  ELSE IF(gs_active%Itor_target>0.d0)THEN
    offset=1
    cofs(1)=gs_active%Itor_target
  END IF
ELSE
  offset=1
  IF(gs_active%Itor_target>0.d0)THEN
    cofs(1)=gs_active%Itor_target
  ELSE
    cofs(1)=gs_active%psiscale
    cofs_scale(1)=1.d0/ABS(gs_active%psiscale)
  END IF
  IF(fit_alam)CALL oft_abort('Lambda cannot be fit in fixed boundary mode.', &
  'fit_gs',__FILE__)
END IF
IF(gs_active%I%ncofs>0.AND.fit_I)THEN
  js = offset; je = offset+gs_active%I%ncofs
  CALL gs_active%I%get_cofs(cofs(js+1:je))
  offset = je
END IF
IF(fit_R0)THEN
  cofs(offset+1)=gs_active%R0_target-gs_active%rmin
  cofs_scale(offset+1)=5.d0/(gs_active%spatial_bounds(2,1)-gs_active%spatial_bounds(1,1))
  offset=offset+1
ELSE IF(fit_Pnorm)THEN
  cofs(offset+1)=gs_active%pnorm
  cofs_scale(offset+1)=1.d0/MAX(gs_active%pnorm,1.d-2)
  offset=offset+1
ELSE IF(gs_active%estore_target>0.d0)THEN
  cofs(offset+1)=gs_active%estore_target
  cofs_scale(offset+1)=1.d0/MAX(gs_active%estore_target,1.d-2)
  offset=offset+1
END IF
IF(fit_V0)THEN
  cofs(offset+1)=gs_active%V0_target
  cofs_scale(offset+1)=40.d0/(gs_active%spatial_bounds(2,2)-gs_active%spatial_bounds(1,2))
  offset=offset+1
END IF
IF(gs_active%P%ncofs>0.AND.fit_P)THEN
  js = offset; je = offset+gs_active%P%ncofs
  CALL gs_active%P%get_cofs(cofs(js+1:je))
  offset = je
END IF
IF(ncond_active>0)THEN
  js = offset; je = offset+ncond_active
  CALL gs_get_cond_weights(gs_active,cofs(js+1:je),.TRUE.)
  CALL gs_get_cond_scales(gs_active,cofs_scale(js+1:je),.TRUE.)
  offset = je
END IF
IF(fit_coils)THEN
  js = offset; je = offset+gs_active%ncoils
  ALLOCATE(curr_in(gs_active%ncoils))
  DO i=1,gs_active%ncoils
    curr_in(i)=gs_active%coil_currs(i)
    cofs(js+i)=0.d0
    cofs_scale(js+i) = 1.0/ABS(gs_active%coil_currs(i))
  END DO
  offset = je
END IF
IF(fit_F0)THEN
  cofs(offset+1) = gs_active%I%f_offset
  cofs_scale(js+1:je) = 1.0/ABS(gs_active%I%f_offset)
END IF
!---Count constraints
ncons=SIZE(conlist)
ALLOCATE(error(ncons))
error=1.d99
WRITE(*,*)
WRITE(*,'(A)')'========================================'
WRITE(*,'(2A)')oft_indent,'Starting Fit:'
CALL oft_increase_indent
WRITE(*,'(2A,I4)')oft_indent,'# of free parameters   = ',ncofs
WRITE(*,'(2A,I4)')oft_indent,'# of constraints       = ',ncons
WRITE(*,*)
CALL oft_decrease_indent
!---
allocate(fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
!---
mode = 2
factor = 1.d0
maxfev = 100
ftol = 1.d-3
xtol = 1.d-3
gtol = 1.d-3
epsfcn = 1.d-3
nprint = 0
ldfjac = ncons
comp_var = .FALSE.
!------------------------------------------------------------------------------
! Load settings
!------------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,gs_fit_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error parsing "gs_fit_options" in input file.','fit_gs',__FILE__)
!---
chi_best=1.d99
ALLOCATE(cofs_best(ncofs))
cofs_best=cofs
CALL gs_active%psi%new(psi_best)
!---Load coefficient values
INQUIRE(EXIST=file_exists,FILE=TRIM(inpath)//'_cofs')
IF(file_exists)THEN
  OPEN(NEWUNIT=io_unit,FILE=TRIM(inpath)//'_cofs')
  READ(io_unit,*)i
  IF(i==ncofs)THEN
    READ(io_unit,*)cofs
    READ(io_unit,*)gs_active%pnorm,gs_active%vcontrol_val
  END IF
  CLOSE(io_unit)
ENDIF
IF(maxfev>0)THEN
  !---Initialize
  CALL fit_error(ncons,ncofs,cofs,error,info)
  ! gs_active%Itor_target=-1.d0
  ! IF(fit_Alam)cofs(1)=gs_active%alam
  IF(gs_active%ierr<0)CALL oft_abort('Initial equilibrium solve failed to converge','fit_gs',__FILE__)
  !---
  call lmder(fit_error_grad,ncons,ncofs,cofs,error,fjac,ldfjac, &
             ftol,xtol,gtol,maxfev,cofs_scale,mode,factor,nprint,info,nfev,njev, &
             ipvt,qtf,wa1,wa2,wa3,wa4)
  WRITE(*,'(A)')'========================================'
  WRITE(*,'(2A)')oft_indent,'Fit Complete '
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'Final error          =',SQRT(SUM(error**2))
  WRITE(*,'(3A)')oft_indent,'Termination reason:  ',minpack_exit_reason(info)
  CALL oft_decrease_indent
  WRITE(*,*)
  IF(SQRT(SUM(error**2))>chi_best)THEN
    cofs=cofs_best
    gs_active%pnorm=pnorm_best
    gs_active%vcontrol_val=vcont_best
    CALL gs_active%psi%add(0.d0,1.d0,psi_best)
  END IF
END IF
!---
feval_count=-1
IF(comp_var)THEN
  CALL fit_confidence(ncons,ncofs,cofs)
  CALL fit_error(ncons,ncofs,cofs,error,info)
ELSE
  CALL fit_error(ncons,ncofs,cofs,error,i)
END IF
!---
OPEN(NEWUNIT=io_unit,FILE=TRIM(outpath))
DO i=1,ncons
  WRITE(io_unit,'(I8,4ES20.12)')i,error(i),conlist(i)%con%eval(gs_active),conlist(i)%con%val, &
  conlist(i)%con%get_nax(gs_active)
END DO
CLOSE(io_unit)
OPEN(NEWUNIT=io_unit,FILE=TRIM(outpath)//'_cofs')
WRITE(io_unit,*)ncofs
WRITE(io_unit,*)cofs
WRITE(io_unit,*)gs_active%pnorm,gs_active%vcontrol_val
CLOSE(io_unit)
!---Cleanup
CALL psi_best%delete
DEALLOCATE(cofs_best,psi_best)
DEALLOCATE(cofs_scale,fjac,qtf,wa1,wa2)
DEALLOCATE(wa3,wa4,ipvt)
CONTAINS
!------------------------------------------------------------------------------
!> Decode MINPACK exit flag info
!!
!! @param[out] info Error flag
!------------------------------------------------------------------------------
function minpack_exit_reason(info) result(exit_reason)
integer(4), intent(in) :: info
CHARACTER(LEN=OFT_ERROR_SLEN) :: exit_reason
! info = 0  improper input parameters.
!
! info = 1  both actual and predicted relative reductions
!           in the sum of squares are at most ftol.
!
! info = 2  relative error between two consecutive iterates
!           is at most xtol.
!
! info = 3  conditions for info = 1 and info = 2 both hold.
! 
! info = 4  the cosine of the angle between fvec and any
!           column of the jacobian is at most gtol in
!           absolute value.
! 
! info = 5  number of calls to fcn with iflag = 1 has
!           reached maxfev.
! 
! info = 6  ftol is too small. no further reduction in
!           the sum of squares is possible.
! 
! info = 7  xtol is too small. no further improvement in
!           the approximate solution x is possible.
! 
! info = 8  gtol is too small. fvec is orthogonal to the
!           columns of the jacobian to machine precision.
SELECT CASE(info)
  CASE(0)
    exit_reason='Improper input parameters'
  CASE(1)
    exit_reason='Actual and predicted delta_rel(err) < ftol'
  CASE(2)
    exit_reason='delta_rel(err) < xtol'
  CASE(3)
    exit_reason='Actual and predicted delta_rel(err) < ftol and delta_rel(err) < xtol'
  CASE(4)
    exit_reason='Cosine of angle between fvec and any jacobian column < gtol'
  CASE(5)
    exit_reason='nfev > maxfev'
  CASE(6)
    exit_reason='ftol is too small. No further reduction possible'
  CASE(7)
    exit_reason='xtol is too small. No further reduction possible'
  CASE(8)
    exit_reason='gtol is too small. No further reduction possible'
  CASE DEFAULT
    exit_reason='Unknown reason'
END SELECT
end function minpack_exit_reason
END SUBROUTINE fit_gs
!---------------------------------------------------------------------------------
! SUBROUTINE fit_confidence
!---------------------------------------------------------------------------------
!> Estimate confidence in the fit from linearization
!---------------------------------------------------------------------------------
SUBROUTINE fit_confidence(m,n,cofs)
INTEGER(4), INTENT(in) :: m,n
REAL(8), INTENT(in) :: cofs(n)
INTEGER(4) :: i,info
REAL(8) :: sigma
REAL(8), ALLOCATABLE :: inv_mat(:,:),cof_tmp(:),error(:),err0(:),jac_mat(:,:)
ALLOCATE(inv_mat(n,n),cof_tmp(n),error(m),jac_mat(m,n),err0(m))
cof_tmp=cofs
CALL fit_error_grad(m,n,cof_tmp,err0,jac_mat,m,1)
CALL fit_error_grad(m,n,cof_tmp,err0,jac_mat,m,2)
sigma=SUM(err0**2)/REAL(m-n,8)
inv_mat=MATMUL(TRANSPOSE(jac_mat),jac_mat)
CALL lapack_matinv(n,inv_mat,info)
inv_mat = sigma*inv_mat
WRITE(*,*)'Confidence intervals'
DO i=1,n
  WRITE(*,'(4X,I3,2ES11.3)')i,cofs(i),inv_mat(i,i)
END DO
END SUBROUTINE fit_confidence
!---------------------------------------------------------------------------------
! SUBROUTINE fit_error
!---------------------------------------------------------------------------------
!>
!---------------------------------------------------------------------------------
SUBROUTINE fit_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
REAL(8) :: jac_mat(m,n)
CALL fit_error_grad(m,n,cofs,err,jac_mat,m,1)
END SUBROUTINE fit_error
!---------------------------------------------------------------------------------
! SUBROUTINE fit_error_grad
!---------------------------------------------------------------------------------
!>
!---------------------------------------------------------------------------------
SUBROUTINE fit_error_grad(m,n,cofs,err,jac_mat,ldjac_mat,iflag)
integer(4), intent(in) :: m,n,ldjac_mat
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m),jac_mat(ldjac_mat,n)
integer(4), intent(in) :: iflag
logical :: plot_save
integer(4) :: i,j,js,je,offset,ierr
real(8) :: rel_err(m),abs_err(m),dx,dxi
real(8), allocatable :: cof_tmp(:),err_tmp(:)
real(8), save :: alam_in,pnorm_in,bounds_in(2,2),vcont_in,ip_target_in,estore_target_in
real(8), allocatable, save :: cofs_in(:)
class(oft_vector), pointer, save :: psi_center => NULL()
CHARACTER(LEN=40) :: err_reason
IF(.NOT.ASSOCIATED(psi_center))THEN
  ALLOCATE(cofs_in(n))
  cofs_in=cofs
  alam_in=gs_active%alam
  ip_target_in=gs_active%Itor_target
  pnorm_in=gs_active%pnorm
  estore_target_in=gs_active%estore_target
  ! bounds_in=gs_active%spatial_bounds
  bounds_in(:,1)=gs_active%plasma_bounds
  vcont_in=gs_active%vcontrol_val
  CALL gs_active%psi%new(psi_center)
  CALL psi_center%add(0.d0,1.d0,gs_active%psi)
END IF
!---
ALLOCATE(cof_tmp(n),err_tmp(m))
!---
IF(iflag==1)THEN
  feval_count=feval_count+1
  IF(oft_env%pm)THEN
    IF(feval_count>0)THEN
      WRITE(*,'(2A,I5)')oft_indent,'Function evaluation',feval_count
    END IF
  END IF
  CALL oft_increase_indent
  IF(ANY(isnan(cofs)))THEN
    WRITE(*,'(2A)')oft_indent,'Step failed: Bad parameters!'
    WRITE(*,*)
    CALL oft_decrease_indent
    err=1.d99
    RETURN
  END IF
  !---
  offset=0
  IF(gs_active%free)THEN
    IF(fit_alam)THEN
      offset=1
      gs_active%alam=cofs(1)
    ELSE IF(gs_active%Itor_target>0.d0)THEN
      offset=1
      gs_active%Itor_target=cofs(1)
    END IF
  ELSE
    offset=1
    IF(gs_active%Itor_target>0.d0)THEN
      gs_active%Itor_target=cofs(1)
    ELSE
      gs_active%psiscale=cofs(1)
    END IF
  END IF
  IF(gs_active%I%ncofs>0.AND.fit_I)THEN
    js = offset; je = offset+gs_active%I%ncofs
    ierr=gs_active%I%set_cofs(cofs(js+1:je))
    IF(ierr<0)THEN
      WRITE(*,'(2A)')oft_indent,'Invalid I coefficients'
      CALL oft_increase_indent
      DO i=js+1,je
        WRITE(*,'(A,ES11.3)')oft_indent,cofs(i)
      END DO
      WRITE(*,*)
      CALL oft_decrease_indent
      DO i=1,m
        err(i)=conlist(i)%con%val*conlist(i)%con%wt
      END DO
      RETURN
    END IF
    offset = je
  END IF
  IF(fit_R0)THEN
    gs_active%R0_target=cofs(offset+1)+gs_active%rmin
    offset=offset+1
  ELSE IF(fit_Pnorm)THEN
    gs_active%pnorm=cofs(offset+1)
    offset=offset+1
  ELSE IF(gs_active%estore_target>0.d0)THEN
    gs_active%estore_target=cofs(offset+1)
    offset=offset+1
  END IF
  IF(fit_V0)THEN
    gs_active%V0_target=cofs(offset+1)
    offset=offset+1
  END IF
  IF(gs_active%P%ncofs>0.AND.fit_P)THEN
    js = offset; je = offset+gs_active%P%ncofs
    ierr=gs_active%P%set_cofs(cofs(js+1:je))
    IF(ierr<0)THEN
      WRITE(*,'(2A)')oft_indent,'Invalid P coefficients'
      CALL oft_increase_indent
      DO i=js+1,je
        WRITE(*,'(A,ES11.3)')oft_indent,cofs(i)
      END DO
      WRITE(*,*)
      CALL oft_decrease_indent
      DO i=1,m
        err(i)=conlist(i)%con%val*conlist(i)%con%wt
      END DO
      RETURN
    END IF
    offset = je
  END IF
  IF(ncond_active>0)THEN
    js = offset; je = offset+ncond_active
    CALL gs_set_cond_weights(gs_active,cofs(js+1:je),.TRUE.)
    offset = je
  END IF
  IF(fit_coils)THEN
    js = offset; je = offset+gs_active%ncoils
    DO i=1,gs_active%ncoils
      gs_active%coil_currs(i)=cofs(js+i)+curr_in(i)
    END DO
    offset = je
  END IF
  IF(fit_F0)gs_active%I%f_offset = cofs(offset+1)
  ! !---Centering
  ! alam_in=gs_active%alam
  ! ip_target_in=gs_active%Itor_target
  ! pnorm_in=gs_active%pnorm
  ! estore_target_in=gs_active%estore_target
  ! bounds_in=gs_active%spatial_bounds
  ! vcont_in=gs_active%vcontrol_val
  ! CALL psi_center%add(0.d0,1.d0,gs_active%psi)
  CALL run_err(.FALSE.,err,m,ierr)
  IF(ierr/=0)THEN
    CALL reset_eq
    err_reason=gs_err_reason(ierr)
  END IF
  !---
  abs_err=0.d0
  rel_err=0.d0
  !$omp parallel do schedule(dynamic,1)
  DO i=1,m
    IF(conlist(i)%con%is_parallel())CYCLE
    abs_err(i)=(conlist(i)%con%val-conlist(i)%con%eval(gs_active))
    IF(ABS(conlist(i)%con%val)>1.d-10)rel_err(i)=abs_err(i)/conlist(i)%con%val
  END DO
  DO i=1,m
    IF(conlist(i)%con%is_parallel())THEN
      abs_err(i)=(conlist(i)%con%val-conlist(i)%con%eval(gs_active))
      IF(ABS(conlist(i)%con%val)>1.d-10)rel_err(i)=abs_err(i)/conlist(i)%con%val
    END IF
  END DO
  IF(ierr==0)THEN
    IF(SQRT(SUM(err**2))<chi_best)THEN
      chi_best=SQRT(SUM(err**2))
      cofs_best=cofs
      alam_best=gs_active%alam
      pnorm_best=gs_active%pnorm
      vcont_best=gs_active%vcontrol_val
      CALL psi_best%add(0.d0,1.d0,gs_active%psi)
    END IF
    alam_in=gs_active%alam
    ip_target_in=gs_active%Itor_target
    pnorm_in=gs_active%pnorm
    estore_target_in=gs_active%estore_target
    ! bounds_in=gs_active%spatial_bounds
    bounds_in(:,1)=gs_active%plasma_bounds
    vcont_in=gs_active%vcontrol_val
    CALL psi_center%add(0.d0,1.d0,gs_active%psi)
  END IF
  IF(oft_env%pm)THEN
    offset=0
    IF(ierr<0)THEN
      WRITE(*,'(3A)')oft_indent,'Step Failed: ',TRIM(err_reason)
    END IF
    WRITE(*,'(2A,ES11.3)')oft_indent,'Alam              =',gs_active%alam
    WRITE(*,'(2A,ES11.3)')oft_indent,'P_scale           =',gs_active%pnorm
    IF(gs_active%R0_target>0.d0)THEN
      WRITE(*,'(2A,ES11.3)')oft_indent,'R0_target         =',gs_active%R0_target
    END IF
    IF(gs_active%V0_target>-1.d98)THEN
      WRITE(*,'(2A,ES11.3)')oft_indent,'V0_target         =',gs_active%V0_target
    END IF
    IF(gs_active%free)THEN
      IF(fit_alam)THEN
        offset=offset+1
      ELSE IF(gs_active%Itor_target>0.d0)THEN
        WRITE(*,'(2A,ES11.3)')oft_indent,'Itor_target       =',gs_active%Itor_target/mu0
        offset=offset+1
      END IF
    ELSE
      IF(gs_active%Itor_target>0.d0)THEN
        WRITE(*,'(2A,ES11.3)')oft_indent,'Itor_target       =',gs_active%Itor_target/mu0
      ELSE
        WRITE(*,'(2A,ES11.3)')oft_indent,'Psi_scale         =',gs_active%psiscale
      END IF
      offset=offset+1
    END IF
    IF(gs_active%I%ncofs>0.AND.fit_I)THEN
      js = offset; je = offset+gs_active%I%ncofs
      WRITE(*,'(2A)',ADVANCE="NO")oft_indent,'F_cofs            ='
      DO i=js+1,je
        WRITE(*,'(ES11.3)',ADVANCE="NO")cofs(i)
      END DO
      WRITE(*,*)
      offset = je
    END IF
    IF(fit_pnorm.OR.fit_R0.OR.(gs_active%estore_target>0.d0))offset=offset+1
    IF(fit_V0)offset=offset+1
    IF(gs_active%P%ncofs>0.AND.fit_P)THEN
      js = offset; je = offset+gs_active%P%ncofs
      WRITE(*,'(2A)',ADVANCE="NO")oft_indent,'P_cofs            ='
      DO i=js+1,je
        WRITE(*,'(ES11.3)',ADVANCE="NO")cofs(i)
      END DO
      WRITE(*,*)
      offset=je
    END IF
    IF(ncond_active>0)THEN
      js = offset; je = offset+ncond_active
      WRITE(*,'(2A)',ADVANCE="NO")oft_indent,'Cond weights      ='
      DO i=js+1,je
        WRITE(*,'(ES11.3)',ADVANCE="NO")cofs(i)
      END DO
      WRITE(*,*)
      offset=je
    END IF
    IF(fit_coils)THEN
      js = offset; je = offset+gs_active%ncoils
      WRITE(*,'(2A)',ADVANCE="NO")oft_indent,'Coil currents [%]  ='
      DO i=js+1,je
        WRITE(*,'(100ES11.3)',ADVANCE="NO")cofs(i)/curr_in(i-js)
      END DO
      WRITE(*,*)
      offset = je
    END IF
    IF(fit_F0)THEN
      js = offset; je = offset+1
      WRITE(*,'(2A,ES11.3)')oft_indent,'f_offset          =',cofs(js+1)
      offset = je
    END IF
    IF(ierr==0)THEN
      WRITE(*,'(2A,ES11.3)')oft_indent,'Maximum Rel Error =',MAXVAL(ABS(rel_err))
      WRITE(*,'(2A,ES11.3)')oft_indent,'Maximum Abs Error =',MAXVAL(ABS(abs_err))
      WRITE(*,'(2A,ES11.3)')oft_indent,'Total Weighted Error   =',SQRT(SUM(err**2))
      WRITE(*,'(2A,ES11.3)')oft_indent,'RMS Weighted Error     =',SQRT(SUM(err**2)/REAL(m,4))
    END IF
    WRITE(*,*)
  END IF
  CALL oft_decrease_indent
  IF(ANY(isnan(err)))THEN
    WRITE(*,'(2A)')oft_indent,'Evaluation failed: NaN detected'
    err=1.d99
  END IF
ELSE
  IF(oft_env%pm)THEN
    geval_count=geval_count+1
    IF(geval_count>0)THEN
      WRITE(*,'(2A,I5)')oft_indent,'Gradient evaluation ',geval_count
    END IF
  END IF
  CALL oft_increase_indent
  plot_save=gs_active%plot_final; gs_active%plot_final=.FALSE.
  !---
  offset=0
  dxi = 1.d-2
  IF(gs_active%free)THEN
    IF(fit_alam)THEN
      CALL reset_eq
      dx = dxi/cofs_scale(offset+1)
      gs_active%alam=cofs(offset+1) + dx
      CALL run_err(.FALSE.,jac_mat(:,offset+1),m,ierr)
      jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
      gs_active%alam=cofs(offset+1)
      offset=1
    ELSE IF(gs_active%Itor_target>0.d0)THEN
      CALL reset_eq
      dx = dxi/cofs_scale(offset+1)
      gs_active%Itor_target=cofs(offset+1) + dx
      CALL run_err(.FALSE.,jac_mat(:,offset+1),m,ierr)
      jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
      gs_active%Itor_target=cofs(offset+1)
      offset=1
    END IF
  ELSE
    IF(gs_active%Itor_target>0.d0)THEN
      gs_active%Itor_target=cofs(1)
    ELSE
      CALL reset_eq
      dx = dxi/cofs_scale(offset+1)
      gs_active%psiscale=cofs(offset+1) + dx
      CALL run_err(linearized_fit,jac_mat(:,offset+1),m,ierr)
      jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
      gs_active%psiscale=cofs(offset+1)
    END IF
    offset=1
  END IF
  IF(gs_active%I%ncofs>0.AND.fit_I)THEN
    js = offset; je = offset+gs_active%I%ncofs
    cof_tmp(1:je-js)=cofs(js+1:je)
    DO j=1,je-js
      CALL reset_eq
      dx = dxi/cofs_scale(js+j)
      cof_tmp(j)=cofs(js+j) + dx
      ierr=gs_active%I%set_cofs(cof_tmp(1:je-js))
      cof_tmp(j)=cofs(js+j)
      IF(ierr<0)THEN
        WRITE(*,'(2A)')oft_indent,'Invalid I coefficients'
        CALL oft_increase_indent
        DO i=js+1,je
          WRITE(*,'(A,ES11.3)')oft_indent,cofs(i)
        END DO
        WRITE(*,*)
        CALL oft_decrease_indent
        DO i=1,m
          jac_mat(i,js+j)=conlist(i)%con%val*conlist(i)%con%wt
        END DO
      ELSE
        CALL run_err(linearized_fit,jac_mat(:,js+j),m,ierr)
      END IF
      jac_mat(:,js+j)=(jac_mat(:,js+j)-err)/dx
    END DO
    cof_tmp(1:je-js)=cofs(js+1:je)
    ierr=gs_active%I%set_cofs(cof_tmp(1:je-js))
    offset = je
  END IF
  IF(fit_R0)THEN
    CALL reset_eq
    dx = 2.0*dxi/cofs_scale(offset+1)
    gs_active%R0_target=cofs(offset+1)+gs_active%rmin + dx
    CALL run_err(.FALSE.,jac_mat(:,offset+1),m,ierr)
    jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
    gs_active%R0_target=cofs(offset+1)+gs_active%rmin
    offset=offset+1
  ELSE IF(fit_Pnorm)THEN
    CALL reset_eq
    dx = dxi/cofs_scale(offset+1)
    gs_active%pnorm=cofs(offset+1) + dx
    CALL run_err(linearized_fit,jac_mat(:,offset+1),m,ierr)
    jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
    gs_active%pnorm=cofs(offset+1)
    offset=offset+1
  ELSE IF(gs_active%estore_target>0.d0)THEN
    CALL reset_eq
    dx = dxi/cofs_scale(offset+1)
    gs_active%estore_target=cofs(offset+1) + dx
    CALL run_err(linearized_fit,jac_mat(:,offset+1),m,ierr)
    jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
    gs_active%estore_target=cofs(offset+1)
    offset=offset+1
  END IF
  IF(fit_V0)THEN
    CALL reset_eq
    dx = 50.0*dxi/cofs_scale(offset+1)
    gs_active%V0_target=cofs(offset+1) + dx
    CALL run_err(.FALSE.,jac_mat(:,offset+1),m,ierr)
    jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
    gs_active%V0_target=cofs(offset+1)
    offset=offset+1
  END IF
  IF(gs_active%P%ncofs>0.AND.fit_P)THEN
    js = offset; je = offset+gs_active%P%ncofs
    cof_tmp(1:je-js)=cofs(js+1:je)
    DO j=1,je-js
      CALL reset_eq
      dx = dxi/cofs_scale(js+j)
      cof_tmp(j)=cofs(js+j) + dx
      ierr=gs_active%P%set_cofs(cof_tmp(1:je-js))
      cof_tmp(j)=cofs(js+j)
      IF(ierr<0)THEN
        WRITE(*,'(2A)')oft_indent,'Invalid P coefficients'
        CALL oft_increase_indent
        DO i=js+1,je
          WRITE(*,'(A,ES11.3)')oft_indent,cofs(i)
        END DO
        WRITE(*,*)
        CALL oft_decrease_indent
        DO i=1,m
          jac_mat(i,js+j)=conlist(i)%con%val*conlist(i)%con%wt
        END DO
      ELSE
        CALL run_err(linearized_fit,jac_mat(:,js+j),m,ierr)
      END IF
      jac_mat(:,js+j)=(jac_mat(:,js+j)-err)/dx
    END DO
    cof_tmp(1:je-js)=cofs(js+1:je)
    ierr=gs_active%P%set_cofs(cof_tmp(1:je-js))
    offset = je
  END IF
  IF(ncond_active>0)THEN
    js = offset; je = offset+ncond_active
    cof_tmp(1:je-js)=cofs(js+1:je)
    DO j=1,je-js
      CALL reset_eq
      dx = dxi/cofs_scale(js+j)
      cof_tmp(j)=cofs(js+j) + dx
      CALL gs_set_cond_weights(gs_active,cof_tmp(1:je-js),.TRUE.)
      CALL run_err(linearized_fit,jac_mat(:,js+j),m,ierr)
      jac_mat(:,js+j)=(jac_mat(:,js+j)-err)/dx
      cof_tmp(j)=cofs(js+j)
    END DO
    cof_tmp(1:je-js)=cofs(js+1:je)
    CALL gs_set_cond_weights(gs_active,cof_tmp(1:je-js),.TRUE.)
    offset = je
  END IF
  IF(fit_coils)THEN
    js = offset; je = offset+gs_active%ncoils
    DO j=1,gs_active%ncoils
      CALL reset_eq
      dx = dxi/cofs_scale(js+j)
      gs_active%coil_currs(j)=cofs(js+j)+curr_in(j)+dx
      CALL run_err(linearized_fit,jac_mat(:,js+j),m,ierr)
      jac_mat(:,js+j)=(jac_mat(:,js+j)-err)/dx
      gs_active%coil_currs(j)=cofs(js+j)+curr_in(j)
    END DO
    offset = je
  END IF
  IF(fit_F0)THEN
    CALL reset_eq
    dx = dxi/cofs_scale(offset+1)
    gs_active%I%f_offset = cofs(offset+1) + dx
    CALL run_err(linearized_fit,jac_mat(:,offset+1),m,ierr)
    jac_mat(:,offset+1)=(jac_mat(:,offset+1)-err)/dx
    gs_active%I%f_offset = cofs(offset+1)
    offset = offset+1
  END IF
  gs_active%plot_final=plot_save
  CALL oft_decrease_indent
  IF(ANY(isnan(jac_mat)))THEN
    CALL oft_abort("Gradient failed: NaN detected", "fit_error_grad", __FILE__)
  END IF
  WRITE(*,*)
END IF
!---Centering
CALL reset_eq
!
DEALLOCATE(cof_tmp,err_tmp)
CONTAINS
SUBROUTINE reset_eq
gs_active%alam=alam_in
gs_active%Itor_target=ip_target_in
gs_active%pnorm=pnorm_in
gs_active%estore_target=estore_target_in
! gs_active%spatial_bounds=bounds_in
gs_active%plasma_bounds=bounds_in(:,1)
gs_active%vcontrol_val=vcont_in
CALL gs_active%psi%add(0.d0,1.d0,psi_center)
END SUBROUTINE reset_eq
END SUBROUTINE fit_error_grad
!
SUBROUTINE run_err(linear,err,m,ierr)
LOGICAL, INTENT(in) :: linear
REAL(8), INTENT(inout) :: err(m)
INTEGER(4), INTENT(in) :: m
INTEGER(4), INTENT(out) :: ierr
INTEGER(4) :: i
REAL(8) :: alamin
LOGICAL :: pm_save
!---
alamin=gs_active%alam
pm_save=oft_env%pm; oft_env%pm=fit_pm
IF(linear)THEN
  call gs_active%lin_solve(.TRUE.,ierr)
ELSE
  call gs_active%solve(ierr)
END IF
oft_env%pm=pm_save
IF(ierr<0)THEN
  CALL gs_active%psi%set(0.d0)
  DO i=1,m
    err(i)=conlist(i)%con%val*conlist(i)%con%wt*2.d0
  END DO
  IF(.NOT.gs_active%free)gs_active%alam=alamin
ELSE
  !$omp parallel do schedule(dynamic,1)
  DO i=1,m
    IF(conlist(i)%con%is_parallel())CYCLE
    err(i)=conlist(i)%con%error(gs_active)
  END DO
  DO i=1,m
    IF(conlist(i)%con%is_parallel())THEN
      err(i)=conlist(i)%con%error(gs_active)
    END IF
  END DO
END IF
END SUBROUTINE run_err
!---------------------------------------------------------------------------------
! SUBROUTINE fit_load
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE fit_load(filename,cons)
CHARACTER(LEN=*), INTENT(in) :: filename
TYPE(fit_constraint_ptr), POINTER, INTENT(out) :: cons(:)
!---
TYPE(coil_constraint), POINTER :: coil_con
TYPE(field_constraint), POINTER :: field_con
TYPE(flux_constraint), POINTER :: flux_con
TYPE(saddle_constraint), POINTER :: saddle_con
TYPE(itor_constraint), POINTER :: itor_con
TYPE(dflux_constraint), POINTER :: dflux_con
TYPE(press_constraint), POINTER :: press_con
! TYPE(injlam_constraint), POINTER :: injlam_con
TYPE(q_constraint), POINTER :: q_con
!---
INTEGER(4) :: i,j,k,n,m,npsi,io_unit,ncons,neddy
REAL(8) :: rtmp(2)
REAL(8), ALLOCATABLE :: nax_corr(:,:),nax_tmp(:,:)
CHARACTER(LEN=2) :: num_str
LOGICAL :: file_exists
CALL oft_increase_indent()
!---
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
READ(io_unit,*)n
!---
ncons=n
IF(fit_coils)ncons=n+gs_active%ncoils
! IF(gs_active%V0_target>-1.d98)ncons=ncons+1
ALLOCATE(cons(ncons))
j=1
!---
IF(fit_coils)THEN
  DO i=1,gs_active%ncoils
    ALLOCATE(coil_con)
    coil_con%coil=i
    coil_con%val=gs_active%coil_currs(i)/mu0
    coil_con%wt=ABS(1.d0/(.05d0*coil_con%val))
    cons(j)%con=>coil_con
    j=j+1
  END DO
END IF
! IF(gs_active%V0_target>-1.d98)THEN
!   ALLOCATE(vcont_constraint::cons(j)%con)
!   cons(j)%con%val=0.d0
!   cons(j)%con%wt=2.d-8
!   j=j+1
! END IF
!---Load
neddy=0
DO i=1,n
  READ(io_unit,*,END=300)
  READ(io_unit,*,END=300)m
  SELECT CASE(m)
    CASE(1)
      ALLOCATE(field_con)
      READ(io_unit,*,END=300)rtmp,field_con%phi
      READ(io_unit,*,END=300)field_con%v
      READ(io_unit,*,END=300)field_con%val,field_con%wt
      field_con%r(1:2)=rtmp
      neddy=neddy+1
      cons(j)%con=>field_con
    CASE(2)
      ALLOCATE(itor_con)
      READ(io_unit,*,END=300)itor_con%val,itor_con%wt
      cons(j)%con=>itor_con
    ! CASE(4)
    !   ALLOCATE(injlam_con)
    !   READ(io_unit,*,END=300)injlam_con%val,injlam_con%wt
    !   cons(j)%con=>injlam_con
    ! CASE(5)
    CASE(7)
      ALLOCATE(flux_con)
      READ(io_unit,*,END=300)rtmp
      READ(io_unit,*,END=300)flux_con%val,flux_con%wt
      flux_con%r(1:2)=rtmp
      cons(j)%con=>flux_con
    CASE(8)
      ALLOCATE(dflux_con)
      READ(io_unit,*,END=300)dflux_con%val,dflux_con%wt
      cons(j)%con=>dflux_con
    CASE(9)
      ALLOCATE(press_con)
      READ(io_unit,*,END=300)rtmp
      READ(io_unit,*,END=300)press_con%val,press_con%wt
      press_con%r(1:2)=rtmp
      cons(j)%con=>press_con
    CASE(10)
      ALLOCATE(q_con)
      READ(io_unit,*,END=300)q_con%type,q_con%loc
      READ(io_unit,*,END=300)q_con%val,q_con%wt
      cons(j)%con=>q_con
    CASE(11)
      ALLOCATE(saddle_con)
      READ(io_unit,*,END=300)rtmp
      saddle_con%r(1:2,1)=rtmp
      READ(io_unit,*,END=300)rtmp
      saddle_con%r(1:2,2)=rtmp
      READ(io_unit,*,END=300)saddle_con%width
      READ(io_unit,*,END=300)saddle_con%val,saddle_con%wt
      neddy=neddy+1
      cons(j)%con=>saddle_con
    CASE DEFAULT
      CALL oft_abort('Invalid constraint type.','fit_load',__FILE__)
  END SELECT
  j=j+1
  IF(j>ncons)EXIT
END DO
CLOSE(io_unit)
!---Load non-axisymmetric corrections
j=0
DO i=1,gs_active%ncond_regs
  j=MAX(j,gs_active%cond_regions(i)%neigs)
END DO
ALLOCATE(nax_corr(gs_active%ncond_eigs,neddy),nax_tmp(j,neddy))
nax_corr=0.d0
DO i=1,gs_active%ncond_regs
  nax_tmp=0.d0
  WRITE(num_str,'(I2.2)')i
  IF(oft_file_exist('wall_eig.rst'))THEN
    CALL hdf5_read(nax_tmp(1:gs_active%cond_regions(i)%neigs,:), 'wall_eig.rst', &
      'corr_'//num_str, success=file_exists)
    WRITE(*,'(2A,I4)')oft_indent,'Non-axisymmetric corrections found: ',i
  END IF
  DO j=1,gs_active%cond_regions(i)%neigs
    nax_corr(gs_active%cond_regions(i)%eig_map(j),:)=nax_corr(gs_active%cond_regions(i)%eig_map(j),:) &
      + nax_tmp(j,:)
  END DO
END DO
IF(ANY(ABS(nax_corr)>0.d0))THEN
  k=0
  DO i=1,ncons
    CALL cons(i)%con%setup_comp()
    IF(cons(i)%con%ncomp==0)CYCLE
    k=k+1
    IF(ANY(ABS(nax_corr(:,k))>0.d0))THEN
      IF(oft_debug_print(1))WRITE(*,'(2A,2I4)')oft_indent,'Non-ax correction: ',i,k
      ALLOCATE(cons(i)%con%nax_corr(gs_active%ncond_eigs))
      cons(i)%con%nax_corr=nax_corr(:,k)
    END IF
  END DO
END IF
DEALLOCATE(nax_corr)
DO j=1,ncons
  cons(j)%con%wt=MAX(0.d0, cons(j)%con%wt)
END DO
CALL oft_decrease_indent()
RETURN
300 CALL oft_abort('EOF reached while reading constraints', &
                   'fit_load',__FILE__)
END SUBROUTINE fit_load
!---------------------------------------------------------------------------------
! FUNCTION fit_dummy_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dummy_error(self,gs) RESULT(err)
CLASS(fit_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err
CALL oft_abort('No constraint type specified.','fit_dummy_error',__FILE__)
err=0.d0
END FUNCTION fit_dummy_error
!---------------------------------------------------------------------------------
! FUNCTION fit_dummy_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dummy_eval(self,gs) RESULT(val)
CLASS(fit_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
CALL oft_abort('No constraint type specified.','fit_dummy_eval',__FILE__)
val=0.d0
END FUNCTION fit_dummy_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_dummy_nax_corr
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dummy_nax_corr(self,gs) RESULT(val)
CLASS(fit_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
val = 0.d0
IF(ASSOCIATED(self%nax_corr))val = DOT_PRODUCT(self%nax_corr,gs%cond_weights)
END FUNCTION fit_dummy_nax_corr
!---------------------------------------------------------------------------------
! FUNCTION fit_dummy_setup_comp
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE fit_dummy_setup_comp(self)
CLASS(fit_constraint), INTENT(inout) :: self
END SUBROUTINE fit_dummy_setup_comp
!---------------------------------------------------------------------------------
! FUNCTION fit_dummy_parallel
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dummy_parallel(self) RESULT(is_parallel)
CLASS(fit_constraint), INTENT(inout) :: self
LOGICAL :: is_parallel
is_parallel=.FALSE.
END FUNCTION fit_dummy_parallel
!---------------------------------------------------------------------------------
! FUNCTION fit_coil_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_coil_error(self,gs) RESULT(err)
CLASS(coil_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err
err = (gs%coil_currs(self%coil)/mu0 - self%val)*self%wt
END FUNCTION fit_coil_error
!---------------------------------------------------------------------------------
! FUNCTION fit_coil_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_coil_eval(self,gs) RESULT(val)
CLASS(coil_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
val = gs%coil_currs(self%coil)/mu0
END FUNCTION fit_coil_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_vcont_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_vcont_error(self,gs) RESULT(err)
CLASS(vcont_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err
err = (gs%vcontrol_val - self%val)*self%wt
END FUNCTION fit_vcont_error
!---------------------------------------------------------------------------------
! FUNCTION fit_vcont_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_vcont_eval(self,gs) RESULT(val)
CLASS(vcont_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
val = gs%vcontrol_val
END FUNCTION fit_vcont_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_q_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_q_error(self,gs) RESULT(err)
CLASS(q_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err,qval
qval = self%eval(gs)
err = (qval - self%val)*self%wt
IF(self%type==-1.AND.qval>self%val)err=0.d0
END FUNCTION fit_q_error
!---------------------------------------------------------------------------------
! FUNCTION fit_q_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_q_eval(self,gs) RESULT(val)
CLASS(q_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
INTEGER(4) :: j
INTEGER(4), PARAMETER :: npsi=30
REAL(8) :: val,qmin,qmax,qprof(npsi),psi_q(npsi),dl,rbounds(2,2),zbounds(2,2),psi0,psi1
psi0=0.d0; psi1=1.d0
! IF(gs%plasma_bounds(1)>-1.d98)THEN
  ! psi0=gs%plasma_bounds(1); psi1=gs%plasma_bounds(2)
  psi0 = psi0 + (psi1-psi0)*2.d-2
  psi1 = psi1 + (psi0-psi1)*2.d-2
! ELSE
!   IF(.NOT.gs%free)psi0 = psi0 + (psi1-psi0)*2.d-2
! END IF
!
IF(ABS(self%type)==1)THEN
  do j=2,31
    psi_q(j-1)=(psi1-psi0)*((j-1)/REAL(npsi,8)) + psi0
  end do
  CALL gs_get_qprof(gs,npsi,psi_q,qprof)
  qmin = MINVAL(qprof)
  qmax = MAXVAL(qprof)
  IF(qmax<0.d0)THEN
    val = ABS(qmax)
  ELSE
    val = qmin
  END IF
ELSEIF(self%type==2)THEN
  do j=2,31
    psi_q(j-1)=(psi1-psi0)*((j-1)/REAL(npsi,8)) + psi0
  end do
  CALL gs_get_qprof(gs,npsi,psi_q,qprof)
  qmin = MINVAL(qprof)
  qmax = MAXVAL(qprof)
  IF(qmax<0.d0)THEN
    val = ABS(qmin)
  ELSE
    val = qmax
  END IF
ELSEIF(self%type==3)THEN
  psi_q(1)=self%loc
  CALL gs_get_qprof(gs,1,psi_q,qprof)
  val = qprof(1)
  ! val = ABS(linterp(qprof(2,:), qprof(1,:), 30, self%loc))
END IF
END FUNCTION fit_q_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_q_parallel
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_q_parallel(self) RESULT(is_parallel)
CLASS(q_constraint), INTENT(inout) :: self
LOGICAL :: is_parallel
is_parallel=.TRUE.
END FUNCTION fit_q_parallel
!---------------------------------------------------------------------------------
! FUNCTION fit_field_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_field_error(self,gs) RESULT(err)
CLASS(field_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err
err = self%eval(gs)
IF(ASSOCIATED(self%nax_corr))err=err+DOT_PRODUCT(self%nax_corr,gs%cond_weights)
err = (self%val - err)*self%wt
END FUNCTION fit_field_error
!---------------------------------------------------------------------------------
! FUNCTION fit_field_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_field_eval(self,gs) RESULT(val)
CLASS(field_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
TYPE(oft_lag_brinterp), TARGET :: psi_eval
TYPE(oft_lag_bginterp), TARGET :: psi_geval
REAL(8) :: goptmp(3,3),v,psi(1),gpsi(3),rmin,rdiff,btmp(3)
INTEGER(4) :: i,ip
CLASS(oft_bmesh), POINTER :: smesh
smesh=>gs%mesh
IF(self%cell==0)THEN
  call bmesh_findcell(smesh,self%cell,self%r,self%f)
  IF((minval(self%f)<-1.d-6).OR.(maxval(self%f)>1.d0+1.d-6))THEN
    ip=0
    rmin=1.d99
    DO i=1,smesh%nbp
      rdiff=SUM((self%r-smesh%r(:,smesh%lbp(i)))**2)
      IF(rdiff<rmin)THEN
        rmin=rdiff
        ip=i
      END IF
    END DO
    IF(oft_debug_print(1))WRITE(*,*)'Point projected',self%r,smesh%r(:,smesh%lbp(ip))
    self%r=smesh%r(:,smesh%lbp(ip))
    self%cell=0
    call bmesh_findcell(smesh,self%cell,self%r,self%f)
  END IF
END IF
!---
psi_eval%u=>gs%psi
CALL psi_eval%setup(gs%fe_rep)
psi_geval%u=>gs%psi
CALL psi_geval%setup(gs%fe_rep)
call smesh%jacobian(self%cell,self%f,goptmp,v)
call psi_eval%interp(self%cell,self%f,goptmp,psi)
call psi_geval%interp(self%cell,self%f,goptmp,gpsi)
btmp(1)= -gs%psiscale*gpsi(2)/self%r(1)
IF(gs%mode==0)THEN
  btmp(2)= gs%psiscale*gs%alam*(gs%I%f(psi(1))+gs%I%f_offset/gs%alam)/self%r(1)
ELSE
  btmp(2)=gs%psiscale*SQRT(gs%alam*gs%I%f(psi(1)) + gs%I%f_offset**2)/self%r(1)
END IF
btmp(3)= gs%psiscale*gpsi(1)/self%r(1)
!---
val = DOT_PRODUCT(btmp,self%v)
END FUNCTION fit_field_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_field_setup_comp
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE fit_field_setup_comp(self)
CLASS(field_constraint), INTENT(inout) :: self
self%ncomp=1
ALLOCATE(self%comp_r(3,1),self%comp_n(3,1))
self%comp_r(:,1)=(/self%r(1),self%phi,self%r(2)/)
self%comp_n(:,1)=self%v/magnitude(self%v)
END SUBROUTINE fit_field_setup_comp
!---------------------------------------------------------------------------------
! FUNCTION fit_flux_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_flux_error(self,gs) RESULT(err)
CLASS(flux_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err,psitmp
psitmp = self%eval(gs)
err = (self%val - psitmp)*self%wt
END FUNCTION fit_flux_error
!---------------------------------------------------------------------------------
! FUNCTION fit_flux_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_flux_eval(self,gs) RESULT(val)
CLASS(flux_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
TYPE(oft_lag_brinterp), TARGET :: psi_eval
REAL(8) :: goptmp(3,3),v,psi(1),rmin,rdiff
INTEGER(4) :: i,ip
CLASS(oft_bmesh), POINTER :: smesh
smesh=>gs%mesh
IF(self%cell==0)THEN
  call bmesh_findcell(smesh,self%cell,self%r,self%f)
  IF((minval(self%f)<-1.d-6).OR.(maxval(self%f)>1.d0+1.d-6))THEN
    ip=0
    rmin=1.d99
    DO i=1,smesh%nbp
      rdiff=SUM((self%r-smesh%r(:,smesh%lbp(i)))**2)
      IF(rdiff<rmin)THEN
        rmin=rdiff
        ip=i
      END IF
    END DO
    IF(oft_debug_print(1))WRITE(*,*)'Point projected',self%r,smesh%r(:,smesh%lbp(ip))
    self%r=smesh%r(:,smesh%lbp(ip))
    self%cell=0
    call bmesh_findcell(smesh,self%cell,self%r,self%f)
  END IF
END IF
!---
psi_eval%u=>gs%psi
CALL psi_eval%setup(gs%fe_rep)
call smesh%jacobian(self%cell,self%f,goptmp,v)
call psi_eval%interp(self%cell,self%f,goptmp,psi)
!---
val = psi(1)*2.d0*pi
END FUNCTION fit_flux_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_saddle_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_saddle_error(self,gs) RESULT(err)
CLASS(saddle_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err
err = self%eval(gs)
IF(ASSOCIATED(self%nax_corr))err=err+DOT_PRODUCT(self%nax_corr,gs%cond_weights)
err = (self%val - err)*self%wt
END FUNCTION fit_saddle_error
!---------------------------------------------------------------------------------
! FUNCTION fit_saddle_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_saddle_eval(self,gs) RESULT(val)
CLASS(saddle_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
TYPE(oft_lag_brinterp), TARGET :: psi_eval
REAL(8) :: goptmp(3,3),psi(1,2)
INTEGER(4) :: i
IF(self%cells(1)==0)THEN
  DO i=1,2
    call bmesh_findcell(gs%mesh,self%cells(i),self%r(:,i),self%f(:,i))
    IF((MINVAL(self%f(:,i))<-1.d-6).OR.(MAXVAL(self%f(:,i))>1.d0+1.d-6))THEN
      CALL oft_abort("Saddle coil off mesh","fit_saddle_error",__FILE__)
    END IF
  END DO
END IF
!---
psi_eval%u=>gs%psi
CALL psi_eval%setup(gs%fe_rep)
DO i=1,2
  call psi_eval%interp(self%cells(i),self%f(:,i),goptmp,psi(:,i))
END DO
!---
val = -(psi(1,1)-psi(1,2))*self%width
END FUNCTION fit_saddle_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_saddle_setup_comp
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE fit_saddle_setup_comp(self)
CLASS(saddle_constraint), INTENT(inout) :: self
INTEGER(4) :: i,j
REAL(8) :: r,z,nr,nz,rmag,theta,da,dr,dtheta
self%ncomp=self%nr*self%nt
ALLOCATE(self%comp_r(3,self%ncomp),self%comp_n(3,self%ncomp))
nr = -(self%r(2,2)-self%r(2,1))
nz = (self%r(1,2)-self%r(1,1))
rmag = SQRT(nr*nr + nz*nz)
nr = nr/rmag
nz = nz/rmag
dr = rmag/REAL(self%nr-1,8)
dtheta = self%width/REAL(self%nt-1,8)
DO i=1,self%nr
  r = (self%r(1,2) - self%r(1,1))*(i-1)/REAL(self%nr-1,8) + self%r(1,1)
  z = (self%r(2,2) - self%r(2,1))*(i-1)/REAL(self%nr-1,8) + self%r(2,1)
  DO j=1,self%nt
    theta = self%width*(j-1)/REAL(self%nt-1,8) - self%width/2.d0
    da = dr*(r*dtheta)
    self%comp_r(:,(i-1)*self%nt+j) = (/r,theta,z/)
    self%comp_n(:,(i-1)*self%nt+j) = da*(/nr,0.d0,nz/)
  END DO
END DO
END SUBROUTINE fit_saddle_setup_comp
!---------------------------------------------------------------------------------
! FUNCTION fit_itor_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_itor_error(self,gs) RESULT(err)
CLASS(itor_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err,itor
itor = self%eval(gs)
err = (self%val - itor)*self%wt
END FUNCTION fit_itor_error
!---------------------------------------------------------------------------------
! FUNCTION fit_itor_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_itor_eval(self,gs) RESULT(val)
CLASS(itor_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val,itor
CALL gs_itor_nl(gs,itor)
val=itor/mu0
END FUNCTION fit_itor_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_dflux_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dflux_error(self,gs) RESULT(err)
CLASS(dflux_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err,dflux
dflux=self%eval(gs)
err = (self%val - dflux)*self%wt
END FUNCTION fit_dflux_error
!---------------------------------------------------------------------------------
! FUNCTION fit_dflux_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_dflux_eval(self,gs) RESULT(val)
CLASS(dflux_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
val=gs_dflux(gs)
END FUNCTION fit_dflux_eval
!---------------------------------------------------------------------------------
! FUNCTION fit_press_error
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_press_error(self,gs) RESULT(err)
CLASS(press_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: err,press,psi(1),goptmp(3,4)
TYPE(oft_lag_brinterp), TARGET :: psi_eval
LOGICAL :: in_plasma
press = self%eval(gs)
!---Handle outside plasma guidance
IF(self%r(1)>0.d0)THEN
  psi_eval%u=>gs%psi
  CALL psi_eval%setup(gs%fe_rep)
  call psi_eval%interp(self%cell,self%f,goptmp,psi)
  in_plasma=.TRUE.
  IF(gs_test_bounds(gs,self%r))THEN
    IF(psi(1)<=gs%plasma_bounds(1))in_plasma=.FALSE.
  ELSE
    in_plasma=.FALSE.
  END IF
  IF(.NOT.in_plasma)THEN
    err=gs%psiscale*gs%psiscale*gs%pnorm*gs%P%f(gs%plasma_bounds(2))/mu0/(gs%plasma_bounds(2)-gs%plasma_bounds(1))
    press=err*(psi(1)-gs%plasma_bounds(1))
  END IF
  CALL psi_eval%delete()
END IF
err = (self%val - press)*self%wt
END FUNCTION fit_press_error
!---------------------------------------------------------------------------------
! FUNCTION fit_press_eval
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
FUNCTION fit_press_eval(self,gs) RESULT(val)
CLASS(press_constraint), INTENT(inout) :: self
TYPE(gs_eq), INTENT(inout) :: gs
REAL(8) :: val
TYPE(oft_lag_brinterp), TARGET :: psi_eval
REAL(8) :: goptmp(3,3),v,psi(1),rmin,rdiff
INTEGER(4) :: i,ip
CLASS(oft_bmesh), POINTER :: smesh
smesh=>gs%mesh
IF(self%r(1)<0.d0)THEN ! Magnetic axis pressure constraint
  val = gs%psiscale*gs%psiscale*gs%pnorm*gs%P%f(gs%plasma_bounds(2))/mu0
  RETURN
END IF
IF(self%cell==0)THEN
  call bmesh_findcell(smesh,self%cell,self%r,self%f)
  IF((minval(self%f)<-1.d-6).OR.(maxval(self%f)>1.d0+1.d-6))THEN
    ip=0
    rmin=1.d99
    DO i=1,smesh%nbp
      rdiff=SUM((self%r-smesh%r(:,smesh%lbp(i)))**2)
      IF(rdiff<rmin)THEN
        rmin=rdiff
        ip=i
      END IF
    END DO
    IF(oft_debug_print(1))WRITE(*,*)'Point projected',self%r,smesh%r(:,smesh%lbp(ip))
    self%r=smesh%r(:,smesh%lbp(ip))
    self%cell=0
    call bmesh_findcell(smesh,self%cell,self%r,self%f)
  END IF
END IF
!---
psi_eval%u=>gs%psi
CALL psi_eval%setup(gs%fe_rep)
call smesh%jacobian(self%cell,self%f,goptmp,v)
call psi_eval%interp(self%cell,self%f,goptmp,psi)
!---
val = gs%psiscale*gs%psiscale*gs%pnorm*gs%P%f(psi(1))/mu0
END FUNCTION fit_press_eval
! !---------------------------------------------------------------------------------
! ! FUNCTION fit_injlam_error
! !---------------------------------------------------------------------------------
! !> Needs Docs
! !---------------------------------------------------------------------------------
! FUNCTION fit_injlam_error(self,gs) RESULT(err)
! CLASS(injlam_constraint), INTENT(inout) :: self
! TYPE(gs_eq), INTENT(inout) :: gs
! REAL(8) :: err,lam
! !---
! lam=gs%I%fp(0.d0)*gs%alam
! IF(oft_debug_print(1))WRITE(*,*)'IL',lam
! err = (self%val - lam)*self%wt
! END FUNCTION fit_injlam_error
! !---------------------------------------------------------------------------------
! ! FUNCTION fit_injlam_eval
! !---------------------------------------------------------------------------------
! !> Needs Docs
! !---------------------------------------------------------------------------------
! FUNCTION fit_injlam_eval(self,gs) RESULT(val)
! CLASS(injlam_constraint), INTENT(inout) :: self
! TYPE(gs_eq), INTENT(inout) :: gs
! REAL(8) :: val
! val=gs%I%fp(0.d0)*gs%alam
! END FUNCTION fit_injlam_eval
END MODULE oft_gs_fit
