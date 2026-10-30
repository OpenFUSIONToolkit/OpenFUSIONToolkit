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
IMPLICIT NONE
INTEGER(i4) :: i,j,io_unit
REAL(r8) :: pt(3)
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
END IF
CALL oft_finalize
END PROGRAM test_thin_wall_kernels
