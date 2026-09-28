!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_stitching.F90
!> Regression tests for owned-entry real and complex dot products.
!!
!! Unowned entries deliberately differ from their owners: their values must not
!! affect a dot product, including after cancellation in a solver work vector.
!! @ingroup testing
!---------------------------------------------------------------------------------
PROGRAM test_stitching
USE oft_base
USE oft_stitching, ONLY: oft_seam, oft_global_dp
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
IMPLICIT NONE
! Exceed the default OpenMP vector threshold (2000).
INTEGER(i4), PARAMETER :: n=4096
REAL(r8), PARAMETER :: ghost=2.d0**40
REAL(r8) :: a(n),b(n),rank_scale,rank_sum,rank_square_sum
COMPLEX(c8) :: ca(n),cb(n),expected_complex
TYPE(oft_seam) :: seam,local_seam
INTEGER(i4) :: nfail

CALL oft_init
nfail=0
rank_scale=REAL(oft_env%rank+1,r8)
rank_sum=REAL(oft_env%nprocs*(oft_env%nprocs+1),r8)/2.d0
rank_square_sum=REAL(oft_env%nprocs*(oft_env%nprocs+1)*(2*oft_env%nprocs+1),r8)/6.d0

! Include an owned boundary entry, an interior entry, and two unowned entries.
seam%full=.TRUE.
seam%nbe=3
ALLOCATE(seam%be(n),seam%lbe(seam%nbe),seam%leo(seam%nbe))
seam%lbe=[1,2,n/2+2]
seam%leo=[.TRUE.,.FALSE.,.FALSE.]
seam%be=.FALSE.
seam%be(seam%lbe)=.TRUE.
! Leave lie unassociated: dot products do not require an interior-entry list.

a=0.d0; b=0.d0
a(1)=rank_scale; b(1)=2.d0
a(n/2+1)=3.d0*rank_scale; b(n/2+1)=4.d0
a([2,n/2+2])=ghost; b([2,n/2+2])=ghost
ca=CMPLX(0.d0,0.d0,c8); cb=CMPLX(0.d0,0.d0,c8)
ca(1)=CMPLX(1.d0,2.d0,c8)*rank_scale
cb(1)=CMPLX(3.d0,4.d0,c8)
ca(n/2+1)=CMPLX(2.d0,-1.d0,c8)*rank_scale
cb(n/2+1)=CMPLX(-1.d0,3.d0,c8)
ca(2)=CMPLX(ghost,ghost,c8); cb(2)=CMPLX(ghost,ghost,c8)
ca(n/2+2)=CMPLX(ghost,-ghost,c8); cb(n/2+2)=CMPLX(ghost,ghost,c8)
! CONJG(1+2i)*(3+4i) + CONJG(2-i)*(-1+3i) = 6+3i.
expected_complex=CMPLX(6.d0,3.d0,c8)

CALL check_real(oft_global_dp(seam,a,b,n),14.d0*rank_scale,'owned real dot')
CALL check_real(oft_global_dp(seam,a,a,n),10.d0*rank_scale**2,'owned real squared norm')
CALL check_complex(oft_global_dp(seam,ca,cb,n),expected_complex*rank_scale,'owned complex dot')
CALL check_complex(oft_global_dp(seam,ca,ca,n),CMPLX(10.d0*rank_scale**2,0.d0,c8),'owned complex squared norm')

! Some seams provide only boundary indices and ownership, without a full mask.
DEALLOCATE(seam%be)
CALL check_real(oft_global_dp(seam,a,b,n),14.d0*rank_scale,'no-mask real dot')
CALL check_complex(oft_global_dp(seam,ca,cb,n),expected_complex*rank_scale,'no-mask complex dot')

! Every rank has different owned values; no_reduce must retain its local result.
seam%full=.FALSE.
CALL check_real(oft_global_dp(seam,a,b,n,no_reduce=.TRUE.),14.d0*rank_scale,'local real dot')
CALL check_complex(oft_global_dp(seam,ca,cb,n,no_reduce=.TRUE.),expected_complex*rank_scale,'local complex dot')
CALL check_real(oft_global_dp(seam,a,b,n),14.d0*rank_sum,'global real dot')
CALL check_real(oft_global_dp(seam,a,a,n),10.d0*rank_square_sum,'global real squared norm')
CALL check_complex(oft_global_dp(seam,ca,cb,n),expected_complex*rank_sum,'global complex dot')
CALL check_complex(oft_global_dp(seam,ca,ca,n),CMPLX(10.d0*rank_square_sum,0.d0,c8),'global complex squared norm')

seam%skip=.TRUE.
CALL check_real(oft_global_dp(seam,a,b,n),0.d0,'skipped real dot')
CALL check_complex(oft_global_dp(seam,ca,cb,n),CMPLX(0.d0,0.d0,c8),'skipped complex dot')

! Local vectors without boundary entries need no allocated ownership arrays.
! In particular, leave local_seam%lie and local_seam%be unassociated.
local_seam%full=.TRUE.
a([2,n/2+2])=0.d0; b([2,n/2+2])=0.d0
ca([2,n/2+2])=CMPLX(0.d0,0.d0,c8); cb([2,n/2+2])=CMPLX(0.d0,0.d0,c8)
CALL check_real(oft_global_dp(local_seam,a,b,n),14.d0*rank_scale,'no-boundary real dot')
CALL check_complex(oft_global_dp(local_seam,ca,cb,n),expected_complex*rank_scale,'no-boundary complex dot')

DEALLOCATE(seam%lbe,seam%leo)
IF(nfail/=0)CALL oft_abort('Stitching dot product checks failed','test_stitching',__FILE__)
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Check a real result for finiteness and agreement with the owned-entry sum.
!---------------------------------------------------------------------------------
SUBROUTINE check_real(actual,expected,label)
REAL(r8), INTENT(in) :: actual !< Computed result
REAL(r8), INTENT(in) :: expected !< Independently computed result
CHARACTER(LEN=*), INTENT(in) :: label !< Description printed on failure
IF(.NOT.IEEE_IS_FINITE(actual))THEN
  nfail=nfail+1
  WRITE(*,*)TRIM(label),': nonfinite result on rank ',oft_env%rank
ELSE IF(ABS(actual-expected)>1.d-12*MAX(1.d0,ABS(expected)))THEN
  nfail=nfail+1
  WRITE(*,*)TRIM(label),': rank, actual, expected = ',oft_env%rank,actual,expected
END IF
END SUBROUTINE check_real
!---------------------------------------------------------------------------------
!> Check both components so a NaN cannot pass a complex comparison.
!---------------------------------------------------------------------------------
SUBROUTINE check_complex(actual,expected,label)
COMPLEX(c8), INTENT(in) :: actual !< Computed result
COMPLEX(c8), INTENT(in) :: expected !< Independently computed result
CHARACTER(LEN=*), INTENT(in) :: label !< Description printed on failure
CALL check_real(REAL(actual,r8),REAL(expected,r8),label//' (real)')
CALL check_real(AIMAG(actual),AIMAG(expected),label//' (imaginary)')
END SUBROUTINE check_complex
END PROGRAM test_stitching
