!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_mapping_find.F90
!
!> Regression tests for tetmesh . Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$, with MG
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_mapping_find
USE oft_base
USE oft_quadrature
USE oft_mesh_type, ONLY: mesh_findcell
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
IMPLICIT NONE
INTEGER(i4) :: ierr
TYPE(multigrid_mesh) :: mg_mesh
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_sphere_id)CALL oft_abort('Wrong mesh type, test for SPHERE only.','main',__FILE__)
IF(oft_env%nprocs>1)CALL oft_abort('Test is for serial meshes only.','main',__FILE__)
!---Run test
CALL check_surface_points(4,1.d-6)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------
!> Validates @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell" by
!! creating points on the boundary of the mesh and comparing them to the
!! locations found by @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell".
!---------------------------------------------------------------------------
SUBROUTINE check_surface_points(order,tol)
INTEGER(i4), INTENT(in) :: order
REAL(r8), INTENT(in) :: tol
INTEGER(i4) :: i,j,cell,fail_count,io_unit
REAL(r8) :: f(4),pt_face(3),pt_cell(3)
TYPE(oft_quad_type) :: quad
CALL mg_mesh%smesh%quad_rule(order,quad)
!---
fail_count=0
!$omp parallel do private(j,cell,f,pt_face,pt_cell) reduction(+:fail_count)
DO i=1,mg_mesh%smesh%nc
  !---Get local reconstructed operators
  DO j=1,quad%np ! Loop over quadrature points
    pt_face=mg_mesh%smesh%log2phys(i,quad%pts(:,j))
    cell=0
    CALL mesh_findcell(mg_mesh%mesh,cell,pt_face,f)
    IF(cell>0)THEN
      pt_cell=mg_mesh%mesh%log2phys(cell,f)
    ELSE
      pt_cell=1.d99
    END IF
    IF(SQRT(SUM((pt_face-pt_cell)**2))>tol)fail_count=fail_count+1
  END DO
END DO
CALL quad%delete
!---
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='mapping_find.results')
  WRITE(io_unit,*)fail_count
  CLOSE(io_unit)
END IF
END SUBROUTINE check_surface_points
END PROGRAM test_mapping_find
