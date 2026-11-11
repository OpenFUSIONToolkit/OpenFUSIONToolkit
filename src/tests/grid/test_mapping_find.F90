!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_mapping_find.F90
!
!> Regression tests for @ref oft_mesh_type::mesh_findcell "mesh_findcell" on
!! tetrahedral and hexahedral meshes of a sphere or cube.
!!
!! The current test cases are:
!! - Locate points on the boundary of the mesh starting from a global search
!! - Locate points inside each cell starting from a guess cell one or two cells away
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_mapping_find
USE oft_base
USE oft_quadrature
USE oft_mesh_type, ONLY: mesh_findcell
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
IMPLICIT NONE
INTEGER(i4) :: io_unit,fail_count
TYPE(multigrid_mesh) :: mg_mesh
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF((mg_mesh%mesh%cad_type/=mesh_sphere_id).AND.(mg_mesh%mesh%cad_type/=mesh_cube_id))THEN
  CALL oft_abort('Wrong mesh type, test for SPHERE and CUBE only.','main',__FILE__)
END IF
IF(oft_env%nprocs>1)CALL oft_abort('Test is for serial meshes only.','main',__FILE__)
!---Run tests
fail_count=check_surface_points(4,1.d-6)
fail_count=fail_count+check_guess_cells(2)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='mapping_find.results')
  WRITE(io_unit,*)fail_count
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Validates @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell" by
!! creating points on the boundary of the mesh and comparing them to the
!! locations found by @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell".
!------------------------------------------------------------------------------
FUNCTION check_surface_points(order,tol) RESULT(fail_count)
INTEGER(i4), INTENT(in) :: order
REAL(r8), INTENT(in) :: tol
INTEGER(i4) :: fail_count
INTEGER(i4) :: i,j,cell
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
IF(fail_count>0)WRITE(*,*)'Surface point failures: ',fail_count
END FUNCTION check_surface_points
!------------------------------------------------------------------------------
! Locate points inside each cell (volume quadrature points) starting from guess
! cells one and two cells away, which exercises the walk from
! cell to cell in mesh_findcell. Interior points belong to a single cell, so the
! search must return that cell.
!------------------------------------------------------------------------------
FUNCTION check_guess_cells(order) RESULT(fail_count)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4) :: fail_count
INTEGER(i4) :: i,j,k,l,cell,guess1,guess2
REAL(r8) :: f(4),pt(3)
TYPE(oft_quad_type) :: quad
CALL mg_mesh%mesh%quad_rule(order,quad)
fail_count=0
!$omp parallel do private(j,k,l,cell,guess1,guess2,f,pt) reduction(+:fail_count)
DO i=1,mg_mesh%mesh%nc
  DO j=1,quad%np
    pt=mg_mesh%mesh%log2phys(i,quad%pts(:,j))
    DO k=1,mg_mesh%mesh%cell_nf
      guess1=mg_mesh%mesh%lcc(k,i)
      IF(guess1==0)CYCLE
      cell=guess1
      CALL mesh_findcell(mg_mesh%mesh,cell,pt,f)
      IF((cell/=i).OR.(.NOT.f_in_cell(f)))fail_count=fail_count+1
      DO l=1,mg_mesh%mesh%cell_nf
        guess2=mg_mesh%mesh%lcc(l,guess1)
        IF((guess2==0).OR.(guess2==i))CYCLE
        cell=guess2
        CALL mesh_findcell(mg_mesh%mesh,cell,pt,f)
        IF((cell/=i).OR.(.NOT.f_in_cell(f)))fail_count=fail_count+1
      END DO
    END DO
  END DO
END DO
CALL quad%delete
IF(fail_count>0)WRITE(*,*)'Guess cell failures: ',fail_count
END FUNCTION check_guess_cells
!------------------------------------------------------------------------------
! All 4 entries of the logical position returned by mesh_findcell lie in [0,1]
! (callers range-check all 4, including the unused 4th entry on hexes)
!------------------------------------------------------------------------------
FUNCTION f_in_cell(f) RESULT(inside)
REAL(r8), INTENT(in) :: f(4)
LOGICAL :: inside
inside=(MINVAL(f)>=-1.d-8).AND.(MAXVAL(f)<=1.d0+1.d-8)
END FUNCTION f_in_cell
END PROGRAM test_mapping_find
