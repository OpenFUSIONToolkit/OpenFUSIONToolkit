!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_h1_fields.F90
!
!> Nedelec H1 FE field handling
!! - Field creation
!! - Restart file handling
!! Nedelec H1(Grad) FE field handling
!! - Field creation
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_h1
!---------------------------------------------------------------------------
MODULE oft_h1_fields
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE oft_h1_basis, ONLY: ML_oft_h1, ML_oft_hgrad
IMPLICIT NONE
#include "local.h"
contains
! !---------------------------------------------------------------------------
! ! SUBROUTINE: oft_h1_create
! !---------------------------------------------------------------------------
! !> Create a H1 field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_h1_create(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_h1%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_h1_create
! !---------------------------------------------------------------------------
! ! SUBROUTINE: oft_hgrad_create
! !---------------------------------------------------------------------------
! !> Create a H1(Grad) field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_hgrad_create(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_hgrad%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_hgrad_create
end module oft_h1_fields
