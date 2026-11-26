!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Module body template for native update kernels (destination y, source x), included by
! native_kernels.F90
!
! Required macros (undefined by the including file after each instantiation):
!   KY     - Destination (and scalar) type (eg. COMPLEX(c8))
!   KX     - Source type (eg. REAL(r8)), must be convertible to KY without loss
!   KYZERO - Zero of type KY
!   KTOY   - Function-like macro converting a KX value to KY
!
! The including file provides the MODULE and END MODULE statements (so the module
! name is visible to CMake's dependency scanner). Specific procedures are private;
! only the generic names are public, so generics of the same name from each
! instantiation merge in oft_native_kernels.
!---------------------------------------------------------------------------------
USE oft_local
IMPLICIT NONE
PRIVATE
PUBLIC :: nk_copy, nk_axpby, nk_axpbypcz, nk_mult, nk_div
INTERFACE nk_copy
  MODULE PROCEDURE k_copy
END INTERFACE nk_copy
INTERFACE nk_axpby
  MODULE PROCEDURE k_axpby
END INTERFACE nk_axpby
INTERFACE nk_axpbypcz
  MODULE PROCEDURE k_axpbypcz
END INTERFACE nk_axpbypcz
INTERFACE nk_mult
  MODULE PROCEDURE k_mult
END INTERFACE nk_mult
INTERFACE nk_div
  MODULE PROCEDURE k_div
END INTERFACE nk_div
CONTAINS
!---------------------------------------------------------------------------------
!> Copy, \f$ y_i = x_i \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_copy(n,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KX, INTENT(in) :: x(n) !< Source array [n]
KY, INTENT(out) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=KTOY(x(i))
END DO
END SUBROUTINE k_copy
!---------------------------------------------------------------------------------
!> Scaled update, \f$ y = \gamma y + \alpha x \f$
!!
!! When \f$ \gamma = 0 \f$ the destination is not read, so non-finite values in `y`
!! are not propagated.
!---------------------------------------------------------------------------------
SUBROUTINE k_axpby(n,gamma,alpha,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KY, INTENT(in) :: gamma !< Scale of destination
KY, INTENT(in) :: alpha !< Scale of source
KX, INTENT(in) :: x(n) !< Source array [n]
KY, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
IF(gamma==KYZERO)THEN
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    y(i)=alpha*x(i)
  END DO
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    y(i)=gamma*y(i)+alpha*x(i)
  END DO
END IF
END SUBROUTINE k_axpby
!---------------------------------------------------------------------------------
!> Scaled update with two sources, \f$ y = \gamma y + \alpha x + \beta z \f$
!!
!! When \f$ \gamma = 0 \f$ the destination is not read, so non-finite values in `y`
!! are not propagated.
!---------------------------------------------------------------------------------
SUBROUTINE k_axpbypcz(n,gamma,alpha,x,beta,z,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KY, INTENT(in) :: gamma !< Scale of destination
KY, INTENT(in) :: alpha !< Scale of first source
KX, INTENT(in) :: x(n) !< First source array [n]
KY, INTENT(in) :: beta !< Scale of second source
KX, INTENT(in) :: z(n) !< Second source array [n]
KY, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
IF(gamma==KYZERO)THEN
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    y(i)=alpha*x(i)+beta*z(i)
  END DO
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    y(i)=gamma*y(i)+alpha*x(i)+beta*z(i)
  END DO
END IF
END SUBROUTINE k_axpbypcz
!---------------------------------------------------------------------------------
!> Elementwise product, \f$ y_i = y_i x_i \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_mult(n,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KX, INTENT(in) :: x(n) !< Source array [n]
KY, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=y(i)*x(i)
END DO
END SUBROUTINE k_mult
!---------------------------------------------------------------------------------
!> Elementwise division, \f$ y_i = y_i / x_i \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_div(n,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KX, INTENT(in) :: x(n) !< Source array [n]
KY, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=y(i)/x(i)
END DO
END SUBROUTINE k_div
