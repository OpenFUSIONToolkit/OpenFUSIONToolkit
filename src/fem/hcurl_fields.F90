!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_hcurl_fields.F90
!
!> Nedelec H1(Curl) FE field handling
!! - Field creation
!! - Restart file handling
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_hcurl
!---------------------------------------------------------------------------
MODULE oft_hcurl_fields
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE oft_hcurl_basis, ONLY: ML_oft_hcurl
IMPLICIT NONE
#include "local.h"
contains
! !---------------------------------------------------------------------------
! ! SUBROUTINE: oft_hcurl_create
! !---------------------------------------------------------------------------
! !> Create a H1(Curl) field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_hcurl_create(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_hcurl%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_hcurl_create
end module oft_hcurl_fields
