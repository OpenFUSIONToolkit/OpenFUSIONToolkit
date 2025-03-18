!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
PROGRAM test_solver_xml
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_solver_xml
!
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, &
  oft_lag_zerob, df_lop, nu_lop, oft_lag_getlop, oft_lag_getmop, lag_getlop_pre
IMPLICIT NONE
INTEGER(i4), PARAMETER :: order=4
INTEGER(i4) :: io_unit,ierr
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
#if defined(HAVE_XML)
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
!---Run tests
oft_env%pm=.FALSE.
CALL test_lap
!---Finalize enviroment
CALL oft_finalize
#else
WRITE(*,*)'SKIP TEST'
#endif
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!---------------------------------------------------------------------------------
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
CALL ML_oft_lagrange%set_level(ML_oft_lagrange%nlevels)
!---Create solver fields
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
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
CALL lag_zerob%apply(v)
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
