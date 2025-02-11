!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_lag.F90
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
PROGRAM test_lag
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: mg_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, &
  oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, lag_inject, &
  lag_interp, lag_zerob, df_lop, nu_lop, oft_lag_getlop, oft_lag_getmop, lag_getlop_pre
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
IMPLICIT NONE
INTEGER(i4) :: minlev,ierr,io_unit
TYPE(xdmf_plot_file) :: plot_file
INTEGER(i4) :: order
LOGICAL :: mg_test
NAMELIST/test_lag_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_lag_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
minlev=2
IF(mesh%type==3)minlev=mg_mesh%mgmax
CALL oft_lag_setup(order,minlev)
IF(mg_test)THEN
  CALL lag_setup_interp
  CALL lag_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(mg_test)THEN
  CALL test_lapmg
ELSE
  CALL plot_file%setup("Test")
  CALL mesh%setup_io(plot_file,order)
  CALL test_lap
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_lap
!------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!------------------------------------------------------------------------------
SUBROUTINE test_lap
!---Create solver objects
CLASS(oft_solver), POINTER :: linv => NULL()
!---Local variables
REAL(r8) :: uu
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL oft_lag_set_level(oft_lagrange_nlevels)
!---Create solver fields
CALL oft_lag_create(u)
CALL oft_lag_create(v)
!---Get FE operators
CALL oft_lag_getlop(lop,'zerob')
CALL oft_lag_getmop(mop,'none')
!---Setup matrix solver
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-3
CALL create_diag_pre(linv%pre)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
CALL u%get_local(vals)
CALL mesh%save_vertex_scalar(vals,plot_file,'T')
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
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
!------------------------------------------------------------------------------
! SUBROUTINE: test_lapmg
!------------------------------------------------------------------------------
!> Same as \ref test_lag::test_lap "test_lap" but use MG preconditioning.
!------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
nlevels=oft_lagrange_nlevels-minlev+1
!---Create solver fields
CALL oft_lag_create(u)
CALL oft_lag_create(v)
!---Get FE operators
CALL oft_lag_getmop(mop,'none')
!---------------------------------------------------------------------------
! Setup matrix solver
!---------------------------------------------------------------------------
CALL create_cg_solver(linv,force_native=.TRUE.)
linv%its=-3
!---Setup MG preconditioner
CALL lag_getlop_pre(linv%pre,ml_lop,nlevels=nlevels)
lop=>ml_lop(nlevels)%M
linv%A=>lop
linv%bc=>lag_zerob
!---------------------------------------------------------------------------
! Solve system
!---------------------------------------------------------------------------
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
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
END PROGRAM test_lag
