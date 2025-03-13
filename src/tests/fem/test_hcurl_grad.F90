!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_hcurl_grad.F90
!
!> Regression tests for vector H(Curl) + Grad(H^1) finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ M B = M \hat{I} \f$
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
program test_hcurl_grad
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
! USE oft_io
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_native_mlpre, create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Grad(H^1) FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_getmop, hcurl_grad_setup_interp, hcurl_grad_getmop_pre, hcurl_grad_mloptions, &
  oft_hcurl_grad_rinterp, oft_hcurl_grad_gzerop
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit,minlev
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_h1,ML_h1grad
TYPE(oft_ml_fem_type), TARGET :: ML_oft_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
LOGICAL :: mg_test
NAMELIST/test_hcurl_grad_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_hcurl_grad_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
minlev=2
IF(oft_env%nprocs>1)minlev=mg_mesh%nbase+1
IF(mg_mesh%mesh%type==3)minlev=mg_mesh%mgdim
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_oft_h1,minlev=minlev)
IF(mg_test)CALL h1_setup_interp(ML_oft_h1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl,minlev=minlev)
IF(mg_test)CALL hcurl_setup_interp(ML_oft_hcurl)
!--- Full H(Curl) space
CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,minlev=minlev)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
IF(mg_test)THEN
  CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_oft_h1,create_full=.TRUE.)
  CALL hcurl_grad_mloptions
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
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Get FE operators
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>mop
winv%its=-3
winv%rtol=1.d-10
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mop
!------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
! TYPE(oft_native_gmres_solver) :: winv
CLASS(oft_solver), POINTER :: winv => NULL()
CLASS(oft_solver), POINTER :: pre
TYPE(oft_matrix_ptr), POINTER :: mats(:) => NULL()
!---Local variables
INTEGER(i4) :: i,nlevels
REAL(r8) :: uu,r,rsq
REAL(r8), POINTER :: bvout(:,:)
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_vector), POINTER :: u,v,x,y
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_hcurl_grad_rinterp) :: field
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Setup matrix solver
nlevels=ML_hcurl_grad%nlevels-ML_hcurl_grad%minlev
CALL create_cg_solver(winv,force_native=.TRUE.)
winv%its=-3
winv%rtol=1.d-10
CALL hcurl_grad_getmop_pre(ML_hcurl_grad,winv%pre,mats,nlevels=nlevels)
mop=>mats(nlevels)%M
winv%A=>mop
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mopmg
end program test_hcurl_grad
