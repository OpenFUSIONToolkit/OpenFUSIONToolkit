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
USE oft_h0_basis, ONLY: oft_h0_fem, oft_h0_setup_vol
IMPLICIT NONE
#include "local.h"
contains
!---------------------------------------------------------------------------
!> Construct Nedelec H1 vector FE on each mesh level
!!
!! @note Highest supported representation is quadratic
!---------------------------------------------------------------------------
subroutine oft_h1_setup(mg_mesh,order,ML_hcurl_obj,ML_h0_obj,ML_hcurl_aug_obj,ML_hgrad_obj,minlev)
type(multigrid_mesh), target, intent(inout) :: mg_mesh
integer(i4), intent(in) :: order !< Order of representation desired
type(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_obj
type(oft_ml_fem_type), intent(inout) :: ML_h0_obj
type(oft_ml_fem_comp_type), intent(inout) :: ML_hcurl_aug_obj
type(oft_ml_fem_type), target, intent(inout) :: ML_hgrad_obj
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
ML_hgrad_obj%ml_mesh=>mg_mesh
!---Allocate multigrid operators
nlevels=mg_mesh%mgdim+(order-1)
IF(minlev_out<0)minlev_out=nlevels
ML_hgrad_obj%minlev=minlev_out
ML_hgrad_obj%nlevels=nlevels
ML_hcurl_aug_obj%minlev=minlev_out
ML_hcurl_aug_obj%nlevels=nlevels
!---Set linear elements
do i=1,mg_mesh%mgdim-1
  IF(i<ML_hcurl_aug_obj%minlev)CYCLE
  CALL multigrid_level(mg_mesh,i)
  if(mg_mesh%level==mg_mesh%nbase)ML_hcurl_aug_obj%blevel=i
  CALL oft_h0_setup_vol(ML_hgrad_obj%levels(i)%fe,mg_mesh%mesh,2)
end do
call multigrid_level(mg_mesh,mg_mesh%mgdim)
!---Set high order elements
do i=1,order
  IF(mg_mesh%mgdim+i-1<ML_hcurl_aug_obj%minlev)CYCLE
  IF(ML_h0_obj%nlevels>=mg_mesh%mgdim+i)THEN
    SELECT TYPE(this=>ML_h0_obj%levels(mg_mesh%mgdim+i)%fe)
      CLASS IS(oft_h0_fem)
      ML_hgrad_obj%levels(mg_mesh%mgdim+i-1)%fe=>this
      CLASS DEFAULT
        CALL oft_abort("Error casting H0 object", "oft_h1_setup", __FILE__)
    END SELECT
  ELSE
    CALL oft_h0_setup_vol(ML_hgrad_obj%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%mesh,i)
  END IF
end do
!---Setup composite structure
ML_hcurl_aug_obj%nfields=2
ALLOCATE(ML_hcurl_aug_obj%ml_fields(ML_hcurl_aug_obj%nfields))
ALLOCATE(ML_hcurl_aug_obj%field_tags(ML_hcurl_aug_obj%nfields))
ML_hcurl_aug_obj%ml_fields(1)%ml=>ML_hcurl_obj
ML_hcurl_aug_obj%field_tags(1)='c'
ML_hcurl_aug_obj%ml_fields(2)%ml=>ML_hgrad_obj
ML_hcurl_aug_obj%field_tags(2)='g'
call ML_hcurl_aug_obj%setup
CALL ML_hcurl_aug_obj%set_level(ML_hcurl_aug_obj%nlevels,propogate=.TRUE.)
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
end subroutine oft_h1_setup
end module oft_h1_basis
