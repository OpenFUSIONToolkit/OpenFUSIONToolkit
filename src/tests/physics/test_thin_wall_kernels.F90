!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_thin_wall_kernels.F90
!
!> Regression test for the ThinCurr triangle potential and its gradient
!! (@ref thin_wall::tw_compute_phipot "tw_compute_phipot" and @ref thin_wall::tw_compute_phigrad "tw_compute_phigrad").
!!
!! Values are written for points off the triangle, in its plane, and on the lines of its edges beyond either end,
!! for a triangle in metre and in millimetre scale, and compared with quadrature references in
!! test_thin_wall_kernels.py.
!!
!! @authors Ian Stewart
!! @date October 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_thin_wall_kernels
USE oft_base
USE thin_wall, ONLY: tw_compute_phipot, tw_compute_phigrad
USE thin_wall_hodlr, ONLY: tw_lowrank_qrcp
IMPLICIT NONE
INTEGER(i4) :: i,j,k,io_unit,rank
REAL(r8) :: pt(3)
REAL(r8), ALLOCATABLE :: A(:,:),U(:,:),V(:,:)
REAL(r8), PARAMETER :: tri(3,3) = RESHAPE([0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0],[3,3]) !< Unit triangle
REAL(r8), PARAMETER :: scales(2) = [1.d0, 1.d-3] !< Metre and millimetre scale
!---Observation points (unit triangle): off the plane, in the plane, and collinear with each edge beyond either end
REAL(r8), PARAMETER :: pts(3,8) = RESHAPE([0.3d0,0.2d0,0.5d0, 1.5d0,0.8d0,-0.4d0, 2.d0,0.5d0,0.d0, &
  2.d0,0.d0,0.d0, -1.d0,0.d0,0.d0, 0.d0,-1.5d0,0.d0, 2.d0,-1.d0,0.d0, -0.5d0,1.5d0,0.d0],[3,8])
CALL oft_init
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_kernels.results')
  DO j=1,2
    DO i=1,8
      pt=pts(:,i)*scales(j)
      WRITE(io_unit,'(8ES24.15)')scales(j),pt,tw_compute_phipot(tri*scales(j),pt),tw_compute_phigrad(tri*scales(j),pt)
    END DO
  END DO
  CLOSE(io_unit)
  !---Low-rank compression: M, N, tolerance, rank and error of each block, then the block
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_qrcp.results')
  DO k=1,3
    CALL qrcp_block(k)
    CALL tw_lowrank_qrcp(A,1.d-8,1,U,V,rank)
    IF(rank>0)THEN
      WRITE(io_unit,'(2I6,ES24.15,I6,ES24.15)')SIZE(A,1),SIZE(A,2),1.d-8,rank,SQRT(SUM((A-MATMUL(U,V))**2))
      DEALLOCATE(U,V)
    ELSE
      WRITE(io_unit,'(2I6,ES24.15,I6,ES24.15)')SIZE(A,1),SIZE(A,2),1.d-8,rank,0.d0
    END IF
    DO i=1,SIZE(A,1)
      WRITE(io_unit,'(100ES24.15)')A(i,:)
    END DO
    DEALLOCATE(A)
  END DO
  CLOSE(io_unit)
END IF
CALL oft_finalize
CONTAINS
!> Test blocks: 1/|x-y| between separated point sets; exact rank 18 (Chebyshev products) with break-even rank 20;
!! 1/(|x-y|+0.01) between interleaved point sets (not compressible)
SUBROUTINE qrcp_block(icase)
INTEGER(i4), INTENT(in) :: icase
INTEGER(i4) :: ii,jj,kk
REAL(r8) :: x,y
SELECT CASE(icase)
CASE(1)
  ALLOCATE(A(60,50))
  DO jj=1,50
    DO ii=1,60
      A(ii,jj)=1.d0/(3.d0+REAL(ii,8)/60.d0+REAL(jj,8)/50.d0)
    END DO
  END DO
CASE(2)
  ALLOCATE(A(40,40))
  A=0.d0
  DO jj=1,40
    y=COS(3.14159265358979d0*(REAL(jj,8)-0.5d0)/40.d0)
    DO ii=1,40
      x=-0.95d0+1.9d0*REAL(ii-1,8)/39.d0
      DO kk=0,17
        A(ii,jj)=A(ii,jj)+COS(REAL(kk,8)*ACOS(x))*COS(REAL(kk,8)*ACOS(y))/REAL(kk+1,8)
      END DO
    END DO
  END DO
CASE(3)
  ALLOCATE(A(30,30))
  DO jj=1,30
    DO ii=1,30
      A(ii,jj)=1.d0/(ABS(REAL(2*ii,8)-REAL(2*jj+1,8))/60.d0+0.01d0)
    END DO
  END DO
END SELECT
END SUBROUTINE qrcp_block
END PROGRAM test_thin_wall_kernels
