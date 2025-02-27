!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_lag_operators.F90
!
!> Nedelec scalar and vector Lagrange FE field handling
!! - Field creation
!! - Restart file handling
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_lag
!---------------------------------------------------------------------------
MODULE oft_lag_fields
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE oft_lag_basis, ONLY: ML_oft_lagrange, ML_oft_vlagrange, ML_oft_blagrange
IMPLICIT NONE
#include "local.h"
contains
! !---------------------------------------------------------------------------
! !> Create a Lagrange scalar field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_lag_create(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_lagrange%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_lag_create
! !---------------------------------------------------------------------------
! !> Create a boundary Lagrange scalar field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_blag_create(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_blagrange%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_blag_create
! !---------------------------------------------------------------------------
! !> Create a Lagrange vector field
! !!
! !! @param[out] new field to create
! !! @param[in] level FE level for init (optional)
! !! @param[in] cache Allow caching (optional)
! !! @param[in] native Force native representation (optional)
! !---------------------------------------------------------------------------
! subroutine oft_lag_vcreate(new,level,cache,native)
! class(oft_vector), pointer, intent(out) :: new
! integer(i4), optional, intent(in) :: level
! logical, optional, intent(in) :: cache
! logical, optional, intent(in) :: native
! DEBUG_STACK_PUSH
! CALL ML_oft_vlagrange%vec_create(new,level=level,cache=cache,native=native)
! DEBUG_STACK_POP
! end subroutine oft_lag_vcreate
end module oft_lag_fields
