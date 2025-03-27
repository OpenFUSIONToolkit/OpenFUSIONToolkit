!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_h1.F90
!
!> Regression tests for scalar h1 finite elements. Tests are performed
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
PROGRAM test_h1
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_la_base, ONLY: oft_vector,oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---
USE oft_h1_basis, ONLY: oft_h1_setup
USE fem_base, ONLY: oft_ml_fem_type
USE oft_h1_operators, ONLY: h1_setup_interp, h1_mloptions, &
  oft_h1_zerob, df_lop, nu_lop, oft_h1_getlop, oft_h1_getmop, h1_getlop_pre
IMPLICIT NONE
INTEGER(i4) :: minlev
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_h1
TYPE(oft_h1_zerob), TARGET :: h1_zerob
INTEGER(i4) :: order,ierr,io_unit
LOGICAL :: mg_test
NAMELIST/test_h1_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_h1_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
minlev=2
IF(oft_env%nprocs>1)minlev=mg_mesh%nbase+1
IF(mg_mesh%mesh%type==3)minlev=mg_mesh%mgmax
CALL oft_h1_setup(mg_mesh,order,ML_oft_h1,minlev=minlev)
h1_zerob%ML_H1_rep=>ML_oft_h1
IF(mg_test)THEN
  CALL h1_setup_interp(ML_oft_h1)
  CALL h1_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(mg_test)THEN
  CALL test_lapmg
ELSE
  CALL test_lap
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_lap
!---Create solver objects
CLASS(oft_solver), POINTER :: linv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_oft_h1%set_level(ML_oft_h1%nlevels)
!---Create solver fields
CALL ML_oft_h1%vec_create(u)
CALL ML_oft_h1%vec_create(v)
!---Get FE operators
CALL oft_h1_getlop(ML_oft_h1%current_level,lop,'zerob')
CALL oft_h1_getmop(ML_oft_h1%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-3
CALL create_diag_pre(linv%pre)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL h1_zerob%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='h1.results')
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
CALL linv%pre%delete
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lap
!---------------------------------------------------------------------------------
!> Same as \ref test_h1::test_lap "test_lap" but use MG preconditioning.
!---------------------------------------------------------------------------------
SUBROUTINE test_lapmg
!---Solver object
CLASS(oft_solver), POINTER :: linv => NULL()
!---Local variables
REAL(r8) :: uu
INTEGER(i4) :: i,nlevels
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_matrix_ptr), POINTER :: ml_lop(:) => NULL()
!------------------------------------------------------------------------------
! Create ML Matrices
!------------------------------------------------------------------------------
nlevels=ML_oft_h1%nlevels-minlev+1
!---Create solver fields
CALL ML_oft_h1%vec_create(u)
CALL ML_oft_h1%vec_create(v)
!---Get FE operators
CALL oft_h1_getmop(ML_oft_h1%current_level,mop,'none')
!------------------------------------------------------------------------------
! Setup matrix solver
!------------------------------------------------------------------------------
CALL create_cg_solver(linv,force_native=.TRUE.)
linv%its=-3
linv%A=>lop
!---Setup MG preconditioner
CALL h1_getlop_pre(ML_oft_h1,linv%pre,ml_lop,'zerob',nlevels=nlevels)
lop=>ml_lop(nlevels)%M
linv%A=>lop
linv%bc=>h1_zerob
!------------------------------------------------------------------------------
! Solve system
!------------------------------------------------------------------------------
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL h1_zerob%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='h1.results')
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
CALL linv%pre%delete
DEALLOCATE(linv%pre)
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lapmg
END PROGRAM test_h1
