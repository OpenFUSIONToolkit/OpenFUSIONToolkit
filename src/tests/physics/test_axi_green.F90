!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_axi_green.F90
!
!> Regression test for the axisymmetric Green's function of a current filament
!! (@ref axi_green::green "green" and @ref axi_green::grad_green "grad_green").
!!
!! Values and gradients are written for points in the medium and far field of one filament and compared with
!! independent references in test_axi_green.py.
!!
!! @authors Ian Stewart
!! @date October 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_axi_green
USE oft_base
USE axi_green, ONLY: green, grad_green
IMPLICIT NONE
INTEGER(i4) :: i,io_unit
REAL(r8) :: fg,gg(2)
REAL(r8), PARAMETER :: rc = 1.5d0 !< Filament radius
REAL(r8), PARAMETER :: zc = 0.2d0 !< Filament height
!---Observation points (R,Z): two in the medium field, two in the far field
REAL(r8), PARAMETER :: pts(2,4) = RESHAPE([2.1d0, 0.6d0, 0.9d0, -1.4d0, 30.d0, 12.d0, 25.d0, -60.d0],[2,4])
CALL oft_init
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='axi_green.results')
  WRITE(io_unit,'(2ES24.15)')rc,zc
  DO i=1,4
    CALL grad_green(pts(1,i),pts(2,i),rc,zc,fg,gg)
    WRITE(io_unit,'(6ES24.15)')pts(:,i),green(pts(1,i),pts(2,i),rc,zc),fg,gg
  END DO
  CLOSE(io_unit)
END IF
CALL oft_finalize
END PROGRAM test_axi_green
