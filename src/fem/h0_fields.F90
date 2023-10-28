!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (OpenFUSIONToolkit)
!---------------------------------------------------------------------------
!> @file oft_h0_fields.F90
!
!> Nedelec H0 FE field handling
!! - Field creation
!! - Restart file handling
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_h0
!---------------------------------------------------------------------------
MODULE oft_h0_fields
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE oft_h0_basis, ONLY: oft_h0, ML_oft_h0
IMPLICIT NONE
#include "local.h"
contains
!---------------------------------------------------------------------------
! SUBROUTINE: oft_h0_create
!---------------------------------------------------------------------------
!> Create a H0 field
!!
!! @param[out] new field to create
!! @param[in] level FE level for init (optional)
!! @param[in] cache Allow caching (optional)
!! @param[in] native Force native representation (optional)
!---------------------------------------------------------------------------
subroutine oft_h0_create(new,level,cache,native)
class(oft_vector), pointer, intent(out) :: new
integer(i4), optional, intent(in) :: level
logical, optional, intent(in) :: cache
logical, optional, intent(in) :: native
DEBUG_STACK_PUSH
CALL ML_oft_h0%vec_create(new,level=level,cache=cache,native=native)
DEBUG_STACK_POP
end subroutine oft_h0_create
end module oft_h0_fields
