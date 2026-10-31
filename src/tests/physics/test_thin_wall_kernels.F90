!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_thin_wall_kernels.F90
!
!> Regression test for ThinCurr kernels and block compression, compared with independent references in
!! test_thin_wall_kernels.py:
!! - triangle potential and its gradient (@ref thin_wall::tw_compute_phipot "tw_compute_phipot",
!!   @ref thin_wall::tw_compute_phigrad "tw_compute_phigrad") at points off the triangle, in its plane, and on the
!!   lines of its edges beyond either end, for a triangle in metre and in millimetre scale (quadrature);
!! - low-rank compression (@ref thin_wall_hodlr::tw_lowrank_qrcp "tw_lowrank_qrcp") of a smooth kernel block, a block
!!   of exact rank just below the break-even rank, a Kahan block whose pivoted QR needs more columns than the
!!   break-even rank although its SVD rank is below it, and an incompressible block (SVD);
!! - ACA+ (@ref thin_wall::tw_aca_plus "tw_aca_plus") of a smooth kernel block with absolute and relative tolerance,
!!   and of an incompressible block (dense block);
!! - straight segment pair integrals (@ref thin_wall::tw_seg_seg_int "tw_seg_seg_int") for close parallel, skew and
!!   collinear pairs (quadrature);
!! - the softened kernel of a circular coil (@ref thin_wall::tw_ring_seg_int "tw_ring_seg_int") along a polygon of
!!   the same coil (softened ring-ring mutual), and the exact kernel along straight segments passing within 1 mm or
!!   1 cm of the coil, tangent or tilted, with the closest point between or on the samples used to find it
!!   (quadrature).
!!
!! @authors Ian Stewart
!! @date October 2026
!! @ingroup testing
!------------------------------------------------------------------------------
MODULE test_thin_wall_kernels_aca
USE thin_wall, ONLY: tw_aca_source
IMPLICIT NONE
!> Dense test block given by rows and columns
TYPE, EXTENDS(tw_aca_source) :: dense_source
  REAL(8), ALLOCATABLE :: K(:,:)
CONTAINS
  PROCEDURE :: row => dense_row
  PROCEDURE :: col => dense_col
END TYPE dense_source
CONTAINS
SUBROUTINE dense_row(self,i,vec)
CLASS(dense_source), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: i
REAL(8), INTENT(out) :: vec(:)
vec=self%K(i,:)
END SUBROUTINE dense_row
SUBROUTINE dense_col(self,i,vec)
CLASS(dense_source), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: i
REAL(8), INTENT(out) :: vec(:)
vec=self%K(:,i)
END SUBROUTINE dense_col
END MODULE test_thin_wall_kernels_aca
!------------------------------------------------------------------------------
PROGRAM test_thin_wall_kernels
USE oft_base
USE thin_wall, ONLY: tw_compute_phipot, tw_compute_phigrad, tw_aca_plus, tw_seg_seg_int, tw_ring_seg_int
USE axi_green, ONLY: green
USE thin_wall_hodlr, ONLY: tw_lowrank_qrcp
USE test_thin_wall_kernels_aca, ONLY: dense_source
IMPLICIT NONE
TYPE(dense_source) :: src
INTEGER(i4) :: i,j,k,io_unit,rank
REAL(r8) :: pt(3),ua(3),pa(3),pb(3)
!---Segments near a ring (R=1, Z=0): distance from the ring, tilt from the tangent (degrees), length, and fraction of
!   the length at which the segment passes the ring
REAL(r8), PARAMETER :: rsegs(4,6) = RESHAPE([1.d-3,30.d0,1.d0,0.0625d0, 1.d-3,30.d0,0.5d0,0.0625d0, &
  1.d-3,30.d0,0.1d0,0.0625d0, 1.d-2,30.d0,1.d0,0.0625d0, 1.d-3,0.d0,1.d0,0.0625d0, 1.d-3,30.d0,1.d0,0.125d0],[4,6])
REAL(r8), ALLOCATABLE :: A(:,:),U(:,:),V(:,:)
!---Segment pairs (a0, a1, b0, b1): partly overlapping parallel, parallel at two offsets, skew, collinear and separated
REAL(r8), PARAMETER :: segs(3,4,5) = RESHAPE([0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,0.01d0,0.d0, 0.5d0,0.01d0,0.d0, &
  0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,0.01d0,0.d0, 1.d0,0.01d0,0.d0, &
  0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.d0,0.001d0,0.d0, 1.d0,0.001d0,0.d0, &
  0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 0.3d0,0.05d0,-0.2d0, 0.6d0,0.02d0,0.4d0, &
  0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.2d0,0.001d0,0.d0, 2.d0,0.001d0,0.d0],[3,4,5])
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
  DO k=1,4
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
  !---Segment pair integrals: the four points, then the integral
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_segseg.results')
  DO k=1,5
    WRITE(io_unit,'(13ES24.15)')segs(:,:,k),tw_seg_seg_int(segs(:,1,k),segs(:,2,k),segs(:,3,k),segs(:,4,k),0.d0)
  END DO
  CLOSE(io_unit)
  !---Softened ring kernel along a 400-point polygon of the same ring (R=1, Z=0, c^2=1E-4), and the softened
  !   ring-ring mutual (height difference sqrt(dz^2+c^2))
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_ringpoly.results')
  pt=0.d0
  DO k=1,400
    pt(1)=pt(1)+tw_ring_seg_int([1.d0,0.d0],[COS(2.d0*pi*(k-1)/400.d0),SIN(2.d0*pi*(k-1)/400.d0),0.d0], &
      [COS(2.d0*pi*k/400.d0),SIN(2.d0*pi*k/400.d0),0.d0],1.d-4)
  END DO
  WRITE(io_unit,'(2ES24.15)')pt(1),-8.d0*pi*pi*green(1.d0,SQRT(1.d-4),1.d0,0.d0)
  CLOSE(io_unit)
  !---Exact ring kernel along segments passing close to the ring: the two end points, then the integral
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_ringseg.results')
  DO k=1,6
    ua=[0.d0,COS(rsegs(2,k)*pi/180.d0),SIN(rsegs(2,k)*pi/180.d0)]
    pa=[1.d0+rsegs(1,k),0.d0,0.d0]-rsegs(4,k)*rsegs(3,k)*ua
    pb=pa+rsegs(3,k)*ua
    WRITE(io_unit,'(7ES24.15)')pa,pb,tw_ring_seg_int([1.d0,0.d0],pa,pb)
  END DO
  CLOSE(io_unit)
  !---ACA+: case, number of terms and relative error against the dense block
  OPEN(NEWUNIT=io_unit,FILE='thin_wall_aca.results')
  ALLOCATE(src%K(60,50))
  DO j=1,50
    DO i=1,60
      src%K(i,j)=1.d0/(3.d0+REAL(j,8)/50.d0-REAL(i,8)/60.d0)
    END DO
  END DO
  CALL tw_aca_plus(src,60,50,1.d-10,1,40,1,1,U,V,rank)
  WRITE(io_unit,'(2I6,ES24.15)')1,rank,SQRT(SUM((src%K-MATMUL(U,V))**2)/SUM(src%K**2))
  DEALLOCATE(U,V)
  CALL tw_aca_plus(src,60,50,1.d-8,1,40,1,1,U,V,rank,rel_tol=.TRUE.)
  WRITE(io_unit,'(2I6,ES24.15)')2,rank,SQRT(SUM((src%K-MATMUL(U,V))**2)/SUM(src%K**2))
  DEALLOCATE(U,V)
  DO j=1,50
    DO i=1,60
      src%K(i,j)=1.d0/(ABS(REAL(i,8)/60.d0-REAL(j,8)/50.d0)+1.d-3)
    END DO
  END DO
  CALL tw_aca_plus(src,60,50,1.d-10,1,10,1,1,U,V,rank)
  WRITE(io_unit,'(2I6,ES24.15)')3,rank,0.d0
  CLOSE(io_unit)
END IF
CALL oft_finalize
CONTAINS
!---Test blocks: 1/|x-y| between separated point sets; exact rank 18 (Chebyshev products) with break-even rank 20;
!   1/(|x-y|+0.01) between interleaved point sets (not compressible); Kahan matrix
!   diag(1,s,...,s^(n-1)) (I - c U) (U strictly upper ones, s^2+c^2=1, columns scaled
!   by (1-1E-3)^j so pivoting keeps the natural order) with n=48 and c=0.8: SVD rank 23 at the tolerance,
!   break-even rank 24, but pivoted QR needs 29 columns
SUBROUTINE qrcp_block(icase)
INTEGER(i4), INTENT(in) :: icase
INTEGER(i4) :: ii,jj,kk
REAL(r8) :: x,y,c,sk
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
CASE(4)
  ALLOCATE(A(48,48))
  c=0.8d0
  sk=SQRT(1.d0-c**2)
  A=0.d0
  DO jj=1,48
    DO ii=1,jj
      A(ii,jj)=-c
      IF(ii==jj)A(ii,jj)=1.d0
      A(ii,jj)=A(ii,jj)*sk**(ii-1)*(1.d0-1.d-3)**(jj-1)*(1.d-8/1.8d-5)
    END DO
  END DO
END SELECT
END SUBROUTINE qrcp_block
END PROGRAM test_thin_wall_kernels
