!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_h1_basis.F90
!
!> @defgroup doxy_oft_h1 Nedelec H1
!! Nedelec H1 finite element implementation for the Open FUSION Toolkit
!! @ingroup doxy_oft_fem
!
!> Base Nedelec H1 FE class and basis evaluation
!! - FE Construction
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_h1
!---------------------------------------------------------------------------
MODULE oft_h1_basis
! USE timer
USE oft_base
USE multigrid, ONLY: multigrid_mesh, multigrid_level
USE oft_la_base, ONLY: oft_matrix, oft_graph
USE fem_base, ONLY: oft_fem_type, oft_fem_ptr, oft_ml_fem_type, oft_bfem_type!, &
! oft_ml_bfem_type
USE fem_composite, ONLY: oft_fem_comp_type, oft_ml_fem_comp_type
USE oft_hcurl_basis, ONLY: ML_oft_hcurl
USE oft_h0_basis, ONLY: oft_h0_fem, oft_h0_setup_vol
IMPLICIT NONE
#include "local.h"
! !---------------------------------------------------------------------------
! !> Needs docs
! !---------------------------------------------------------------------------
! type :: h1_ops
!   class(oft_matrix), pointer :: interp => NULL() !< Interpolation graph
!   type(oft_graph), pointer :: hgrad_interp_graph => NULL() !< Interpolation graph
!   class(oft_matrix), pointer :: hgrad_interp => NULL() !< Interpolation graph
! end type h1_ops
!---Global Variables
! integer(i4) :: oft_h1_blevel = 0 !< Highest level on base meshes
! integer(i4) :: oft_h1_level = 0 !< Active FE level
! integer(i4) :: oft_h1_lin_level = 0 !< Highest linear element level
! integer(i4) :: oft_h1_minlev = 0 !<
! integer(i4) :: oft_h1_nlevels = 0 !< Number of total levels
! type(h1_ops), pointer :: oft_h1_ops !< Active operators
! type(h1_ops), pointer :: oft_h1_ops_lin !< Highest linear element operators
! type(h1_ops), pointer :: ML_oft_h1_ops(:) !< ML container for all operators
!---
! class(oft_h0_fem), pointer :: oft_hgrad !< Active FE representation
! class(oft_h0_fem), pointer :: oft_hgrad_lin !< Highest linear element representation
type(oft_ml_fem_type), TARGET :: ML_oft_hgrad !< ML container for all FE representations
!
! type(oft_fem_comp_type), pointer :: oft_h1 !< Active H1 representation
type(oft_ml_fem_comp_type), TARGET :: ML_oft_h1 !< ML container for H1 representations
contains
! !---------------------------------------------------------------------------
! !> Set the current level for Nedelec H1 finite elements
! !---------------------------------------------------------------------------
! subroutine oft_h1_set_level(level)
! integer(i4), intent(in) :: level !< Desired level
! DEBUG_STACK_PUSH
! ! if(level>oft_h1_nlevels.OR.level<=0)then
! !   write(*,*)level
! !   call oft_abort('Invalid FE level','oft_h1_set_level',__FILE__)
! ! end if
! ! if(level<mg_mesh%mgdim)then
! !   call multigrid_level(level)
! ! else
! !   call multigrid_level(mg_mesh%mgdim)
! ! end if
! CALL ML_oft_h1%set_level(level)
! ! oft_h1=>ML_oft_h1%current_level
! CALL ML_oft_hgrad%set_level(level)
! ! SELECT TYPE(this=>ML_oft_hgrad%current_level)
! !   CLASS IS(oft_h0_fem)
! !     oft_hgrad=>this
! !   CLASS DEFAULT
! !     CALL oft_abort("Error casting H0 object", "oft_h1_set_level", __FILE__)
! ! END SELECT
! ! oft_hgrad=>ML_oft_hgrad%current_level
! !---
! ! oft_h1_level=level
! ! oft_h1_ops=>ML_oft_h1_ops(level)
! ! call oft_hcurl_set_level(level)
! CALL ML_oft_hcurl%set_level(level)
! DEBUG_STACK_POP
! end subroutine oft_h1_set_level
!---------------------------------------------------------------------------
!> Construct Nedelec H1 vector FE on each mesh level
!!
!! @note Highest supported representation is quadratic
!---------------------------------------------------------------------------
subroutine oft_h1_setup(mg_mesh,order,ML_hcurl_obj,ML_h0_obj,ML_hcurl_aug_obj,minlev)
type(multigrid_mesh), target, intent(inout) :: mg_mesh
integer(i4), intent(in) :: order !< Order of representation desired
type(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_obj
type(oft_ml_fem_type), intent(inout) :: ML_h0_obj
type(oft_ml_fem_comp_type), intent(inout) :: ML_hcurl_aug_obj
integer(i4), optional, intent(in) :: minlev
integer(i4) :: i,nlevels,minlev_out
DEBUG_STACK_PUSH
minlev_out=1
IF(PRESENT(minlev))minlev_out=minlev
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'**** Creating Nedelec H1 FE space'
  WRITE(*,'(2X,A,I4)')'Order  = ',order
  WRITE(*,'(2X,A,I4)')'Minlev = ',minlev_out
END IF
ML_oft_hgrad%ml_mesh=>mg_mesh
!---Allocate multigrid operators
nlevels=mg_mesh%mgdim+(order-1)
IF(minlev_out<0)minlev_out=nlevels
! allocate(ML_oft_h1_ops(nlevels))
ML_oft_hgrad%minlev=minlev_out
ML_oft_hgrad%nlevels=nlevels
ML_hcurl_aug_obj%minlev=minlev_out
ML_hcurl_aug_obj%nlevels=nlevels

! oft_h1_minlev=minlev
! oft_h1_nlevels=nlevels
!---Set linear elements
do i=1,mg_mesh%mgdim-1
  IF(i<ML_hcurl_aug_obj%minlev)CYCLE
  ! ALLOCATE(oft_h0_fem::ML_oft_hgrad%levels(i)%fe)
  CALL multigrid_level(mg_mesh,i)
  if(mg_mesh%level==mg_mesh%nbase)THEN
    ML_hcurl_aug_obj%blevel=i
    ! oft_h1_blevel=i
  END IF
  CALL oft_h0_setup_vol(ML_oft_hgrad%levels(i)%fe,mg_mesh%mesh,2)
  ! !---
  ! oft_hgrad%mesh=>mg_mesh%mesh
  ! oft_hgrad%order=2
  ! oft_hgrad%dim=1
  ! oft_hgrad%gstruct=(/1,1,0,0/)
  ! call oft_hgrad%setup(5)
end do
call multigrid_level(mg_mesh,mg_mesh%mgdim)
!---Set high order elements
do i=1,order
  IF(mg_mesh%mgdim+i-1<ML_hcurl_aug_obj%minlev)CYCLE
  ! ALLOCATE(oft_h0_fem::ML_oft_hgrad%levels(mg_mesh%mgdim+i-1)%fe)
  ! call oft_h1_set_level(mg_mesh%mgdim+i-1)
  IF(ML_h0_obj%nlevels>=mg_mesh%mgdim+i)THEN
    ! DEALLOCATE(ML_oft_hgrad%levels(oft_h1_level)%fe)
    ! ML_oft_hgrad%levels(mg_mesh%mgdim+i-1)%fe=>ML_oft_h0%levels(mg_mesh%mgdim+i)%fe
    SELECT TYPE(this=>ML_h0_obj%levels(mg_mesh%mgdim+i)%fe)
      CLASS IS(oft_h0_fem)
        ML_oft_hgrad%levels(mg_mesh%mgdim+i-1)%fe=>this
      CLASS DEFAULT
        CALL oft_abort("Error casting H0 object", "oft_h1_setup", __FILE__)
    END SELECT
    ! oft_hgrad=>ML_oft_hgrad%levels(oft_h1_level)%fe
  ELSE
    CALL oft_h0_setup_vol(ML_oft_hgrad%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%mesh,i)
    ! !---
    ! oft_hgrad%mesh=>mg_mesh%mesh
    ! oft_hgrad%order=i+1
    ! oft_hgrad%dim=1
    ! select case(oft_hgrad%order)
    !   case(2)
    !     oft_hgrad%gstruct=(/1,1,0,0/)
    !   case(3)
    !     oft_hgrad%gstruct=(/1,2,1,0/)
    !   case(4)
    !     oft_hgrad%gstruct=(/1,3,3,1/)
    !   case(5)
    !     oft_hgrad%gstruct=(/1,4,6,4/)
    !   case default
    !     call oft_abort('Invalid polynomial degree (npmax=4)','oft_h1_setup',__FILE__)
    ! end select
    ! call oft_hgrad%setup((i+1)*2+1)
  END IF
end do
! IF(mg_mesh%mgdim>=ML_hcurl_aug_obj%minlev)THEN
!   oft_h1_lin_level=mg_mesh%mgdim
!   ! oft_h1_ops_lin=>ML_oft_h1_ops(mg_mesh%mgdim)
! ! ELSE
! !   oft_h1_lin_level=-1
! END IF
!---Setup composite structure
! ML_oft_h1%nlevels=oft_h1_nlevels
ML_hcurl_aug_obj%nfields=2
ALLOCATE(ML_hcurl_aug_obj%ml_fields(ML_hcurl_aug_obj%nfields))
ALLOCATE(ML_hcurl_aug_obj%field_tags(ML_hcurl_aug_obj%nfields))
ML_hcurl_aug_obj%ml_fields(1)%ml=>ML_hcurl_obj
ML_hcurl_aug_obj%field_tags(1)='c'
ML_hcurl_aug_obj%ml_fields(2)%ml=>ML_oft_hgrad
ML_hcurl_aug_obj%field_tags(2)='g'
call ML_hcurl_aug_obj%setup
CALL ML_hcurl_aug_obj%set_level(ML_hcurl_aug_obj%nlevels,propogate=.TRUE.)
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
end subroutine oft_h1_setup
end module oft_h1_basis
