!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_gmsh.F90
!
!> Regression tests for the Gmsh mesh interface. Test perform refinements
!! and compare the resulting mesh volume and surface area to baseline case.
!!
!! The current test cases are:
!! - Unit cylinder
!!
!! @authors Chris Hansen
!! @date January 2016
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_gmsh
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: mesh, smesh
USE oft_mesh_native, ONLY: mesh_native_id
USE oft_mesh_gmsh, ONLY: mesh_gmsh_id
USE multigrid_build, ONLY: multigrid_construct
IMPLICIT NONE
INTEGER(i4) :: io_unit
TYPE(xdmf_plot_file) :: plot_file
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct
IF(ALL(mesh%cad_type/=[mesh_gmsh_id,mesh_native_id]))CALL oft_abort('Wrong mesh type.','main',__FILE__)
CALL plot_file%setup("Test")
CALL mesh%setup_io(plot_file,1)
IF(oft_env%head_proc)OPEN(NEWUNIT=io_unit,FILE='gmsh.results')
CALL compute_volume
CALL compute_area
IF(oft_env%head_proc)CLOSE(io_unit)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE compute_volume
!---------------------------------------------------------------------------
!> Compute volume of the current mesh and output
!---------------------------------------------------------------------------
SUBROUTINE compute_volume
INTEGER(i4) :: i,m
REAL(r8) :: v,det,goptmp(3,4),volume
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
volume=0._r8
CALL mesh%quad_rule(8,quad)
!---
DO i=1,mesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL mesh%jacobian(i,quad%pts(:,m),goptmp,v)
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
!---------------------------------------------------------------------------
! SUBROUTINE compute_area
!---------------------------------------------------------------------------
!> Compute surface area of the current mesh and output
!---------------------------------------------------------------------------
SUBROUTINE compute_area
INTEGER(i4) :: i,j,m
REAL(r8) :: a,det,goptmp(3,4),area
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
area=0._r8
CALL smesh%quad_rule(8,quad)
!---
DO i=1,smesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL smesh%jacobian(i,quad%pts(:,m),goptmp,a)
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
END PROGRAM test_gmsh
