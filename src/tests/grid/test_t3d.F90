!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_t3d.F90
!
!> \defgroup testing Tests
!! Regression tests for the Open FUSION Toolkit
!
!> Regression tests for the CUBIT mesh interface. Test perform refinements
!! and compare the resulting mesh volume and surface area with results from
!! refinement within CUBIT.
!!
!! The current test cases are:
!! - Unit sphere
!! - Unit sphere with cuts
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_t3d
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_quadrature
USE oft_mesh_t3d, ONLY: mesh_t3d_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
IMPLICIT NONE
INTEGER(i4) :: io_unit
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_t3d_id)CALL oft_abort('Wrong mesh type, test for T3D only.','main',__FILE__)
CALL plot_file%setup("Test")
CALL mg_mesh%mesh%setup_io(plot_file,1)
IF(oft_env%head_proc)OPEN(NEWUNIT=io_unit,FILE='t3d.results')
CALL compute_volume
CALL compute_area
IF(oft_env%head_proc)CLOSE(io_unit)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Compute volume of the current mesh and output
!------------------------------------------------------------------------------
SUBROUTINE compute_volume
INTEGER(i4) :: i,m
REAL(r8) :: v,det,goptmp(3,4),volume
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
volume=0._r8
CALL mg_mesh%mesh%quad_rule(8,quad)
!---
DO i=1,mg_mesh%mesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL mg_mesh%mesh%jacobian(i,quad%pts(:,m),goptmp,v)
    det=v*quad%wts(m)
    volume=volume+det
  END DO
END DO
CALL quad%delete
volume=oft_mpi_sum(volume)
IF(oft_env%head_proc)THEN
  WRITE(*,*)'Mesh Volume =',volume
  WRITE(io_unit,*)volume
END IF
END SUBROUTINE compute_volume
!------------------------------------------------------------------------------
!> Compute surface area of the current mesh and output
!------------------------------------------------------------------------------
SUBROUTINE compute_area
INTEGER(i4) :: i,j,m
REAL(r8) :: a,det,goptmp(3,4),area
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
area=0._r8
CALL mg_mesh%smesh%quad_rule(8,quad)
DO i=1,mg_mesh%smesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL mg_mesh%smesh%jacobian(i,quad%pts(:,m),goptmp,a)
    det=a*quad%wts(m)
    area=area+det
  END DO
END DO
CALL quad%delete
area=oft_mpi_sum(area)
IF(oft_env%head_proc)THEN
  WRITE(*,*)'Mesh Area   =',area
  WRITE(io_unit,*)area
END IF
END SUBROUTINE compute_area
END PROGRAM test_t3d
