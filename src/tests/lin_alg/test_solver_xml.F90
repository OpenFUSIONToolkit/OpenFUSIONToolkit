!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_solver_xml.F90
!
!> Regression tests for scalar Lagrange finite elements. Tests are performed
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
PROGRAM test_solver_xml
USE oft_base
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid_build, ONLY: multigrid_construct
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_solver_xml
!
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, &
oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, lag_inject, &
lag_interp, lag_zerob, df_lop, nu_lop, oft_lag_getlop, oft_lag_getmop, lag_getlop_pre
IMPLICIT NONE
INTEGER(i4), PARAMETER :: order=4
INTEGER(i4) :: io_unit,ierr
#if defined(HAVE_XML)
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
CALL oft_lag_setup(order)
!---Run tests
oft_env%pm=.FALSE.
CALL test_lap
!---Finalize enviroment
CALL oft_finalize
#else
WRITE(*,*)'SKIP TEST'
#endif
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_lap
!------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!------------------------------------------------------------------------------
SUBROUTINE test_lap
!---Create solver objects
CLASS(oft_solver), POINTER :: linv
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: solver_node
#endif
!---Set FE level
CALL oft_lag_set_level(oft_lagrange_nlevels)
!---Create solver fields
CALL oft_lag_create(u)
CALL oft_lag_create(v)
!---Get FE operators
CALL oft_lag_getlop(lop,'zerob')
CALL oft_lag_getmop(mop,'none')
!---Setup matrix solver
#ifdef HAVE_XML
CALL xml_get_element(oft_env%xml,"solver",solver_node,ierr)
IF(ierr==0)THEN
  CALL create_solver_xml(linv,solver_node)
ELSE
  CALL oft_abort('Could not find XML node.','test_lap',__FILE__)
END IF
#endif
linv%A=>lop
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='xml.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrices
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
!---Destory preconditioner
IF(ASSOCIATED(linv%pre))CALL linv%pre%delete
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lap
END PROGRAM test_solver_xml
