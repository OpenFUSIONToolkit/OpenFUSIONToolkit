!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_arpack.F90
!
!> Regression tests for ARPACK interface. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Compute eigenvalues of Lagrange Poisson system
!!
!! @authors Chris Hansen
!! @date August 2017
!! @ingroup testing
!--------------------------------------------------------------------------------
PROGRAM test_arpack
#if defined( HAVE_ARPACK )
USE oft_base
!--Grid
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_deriv_matrices, ONLY: create_diagmatrix
USE oft_solver_utils, ONLY: create_native_pre
USE oft_arpack, ONLY: oft_irlm_eigsolver, oft_iram_eigsolver
!---Lagrange FE space
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, oft_lag_getlop, oft_lag_zerob
IMPLICIT NONE
INTEGER(i4) :: iounit,ierr
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
INTEGER(i4) :: order=1
INTEGER(i4) :: minlev=1
NAMELIST/test_arpack_options/order,minlev
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_arpack_options,IOSTAT=ierr)
CLOSE(iounit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!------------------------------------------------------------------------------
! Run tests
!------------------------------------------------------------------------------
oft_env%pm=.FALSE.
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=minlev)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
CALL test_lop_eig()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator LAG::LOP
!------------------------------------------------------------------------------
SUBROUTINE test_lop_eig()
INTEGER(i4) :: i
REAL(r8) :: lam1,lam2
CLASS(oft_vector), POINTER :: u => NULL()
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: lop => NULL()
TYPE(oft_irlm_eigsolver) :: arsolver1
TYPE(oft_iram_eigsolver) :: arsolver2
!------------------------------------------------------------------------------
! Compute optimal smoother coefficients
!------------------------------------------------------------------------------
IF(oft_env%head_proc)OPEN(NEWUNIT=iounit,FILE='arpack.results')
DO i=ML_oft_lagrange%nlevels-order+1,ML_oft_lagrange%nlevels
  CALL ML_oft_lagrange%set_level(i)
  !---Create fields
  CALL ML_oft_lagrange%vec_create(u)
  !---Create matrices
  CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
  CALL create_diagmatrix(md,lop%D)
  !---Test Lanzcos solver
  arsolver1%A=>lop
  arsolver1%M=>md
  arsolver1%tol=1.E-5_r8
  arsolver1%bc=>lag_zerob
  arsolver1%mode=2
  CALL create_native_pre(arsolver1%Minv, "jacobi")
  arsolver1%Minv%A=>lop
  CALL arsolver1%max(u,lam1)
  !---Test Arnoldi solver
  arsolver2%A=>lop
  arsolver2%M=>md
  arsolver2%tol=1.E-5_r8
  arsolver2%bc=>lag_zerob
  arsolver2%mode=2
  arsolver2%Minv=>arsolver1%Minv
  CALL arsolver2%max(u,lam2)
  IF(oft_env%head_proc)WRITE(iounit,*)lam1,lam2
  !---Delete field and matrices
  CALL u%delete
  CALL md%delete
  CALL lop%delete
  DEALLOCATE(u,md,lop)
END DO
!---Output
IF(oft_env%head_proc)CLOSE(iounit)
END SUBROUTINE test_lop_eig
#else
WRITE(*,*)'SKIP TEST'
STOP
#endif
END PROGRAM test_arpack
