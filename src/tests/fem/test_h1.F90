!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_h1.F90
!
!> Regression tests for vector H1 finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ M B = M \hat{I} \f$
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
program test_h1
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
! USE oft_io
USE multigrid, ONLY: mg_mesh, multigrid_level
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_native_mlpre, create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_vcreate
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vgetmop, oft_lag_vproject
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_level, oft_hcurl_nlevels
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions, hcurl_zerob
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_nlevels, oft_h1_set_level, oft_h1_ops, &
  oft_h1_level
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1_getmop, h1_setup_interp, h1_getmop_pre, h1_mloptions, &
  oft_h1_rinterp, h1_zerob, h1grad_zerop
IMPLICIT NONE
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
CALL multigrid_construct
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---H1(Curl) subspace
CALL oft_hcurl_setup(order)
IF(mg_test)CALL hcurl_setup_interp
!---H1(Grad) subspace
CALL oft_h0_setup(order+1)
IF(mg_test)CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order)
IF(mg_test)THEN
  CALL h1_setup_interp(create_full=.TRUE.)
  CALL h1_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
IF(mg_test)THEN
  CALL test_mopmg
ELSE
  CALL test_mop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_mop
!------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H1(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL oft_h1_set_level(oft_h1_nlevels)
!---Create solver fields
CALL oft_h1_create(u)
CALL oft_h1_create(v)
!---Get FE operators
CALL h1_getmop(mop,'none')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>mop
winv%its=-3
winv%rtol=1.d-10
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(1.d0)
CALL h1grad_zerop(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='h1.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mop
!------------------------------------------------------------------------------
! SUBROUTINE: test_mopmg
!------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H1(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
TYPE(oft_native_gmres_solver) :: winv
CLASS(oft_solver), POINTER :: pre
TYPE(oft_matrix_ptr), POINTER :: mats(:) => NULL()
!---Local variables
INTEGER(i4) :: i,nlevels
REAL(r8) :: uu,r,rsq
REAL(r8), POINTER :: bvout(:,:)
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_vector), POINTER :: u,v,x,y
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_h1_rinterp) :: field
!---Set FE level
CALL oft_h1_set_level(oft_h1_nlevels)
!---Create solver fields
CALL oft_h1_create(u)
CALL oft_h1_create(v)
!---Setup matrix solver
nlevels=oft_h1_nlevels-1
CALL h1_getmop_pre(winv%pre,mats,nlevels=nlevels)
mop=>mats(nlevels)%M
winv%A=>mop
winv%nrits=20
winv%its=-3
winv%rtol=1.d-9
!---Solve
CALL u%set(1.d0)
CALL h1grad_zerop(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='h1.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mopmg
end program test_h1
