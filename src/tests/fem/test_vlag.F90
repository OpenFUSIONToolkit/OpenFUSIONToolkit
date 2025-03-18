!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_vlag.F90
!
!> Regression tests for vector Lagrange finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Invert the mass matrix
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_vlag
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type, oft_ml_fe_comp_vecspace
USE oft_lag_basis, ONLY: oft_lag_setup!, &
  ! ML_oft_lagrange, ML_oft_blagrange, ML_oft_vlagrange
USE oft_lag_operators, ONLY: oft_lag_vgetmop, oft_vlag_zerob, &
  df_lop, nu_lop, lag_setup_interp, lag_mloptions, oft_lag_getlop
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, oft_graph_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_la_utils, ONLY: create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit,minlev
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange,ML_oft_blagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
TYPE(oft_vlag_zerob), TARGET :: vlag_zerob
LOGICAL :: mg_test=.FALSE.
NAMELIST/test_lag_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_lag_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
minlev=2
IF(oft_env%nprocs>1)minlev=mg_mesh%nbase+1
IF(mg_mesh%mesh%type==3)minlev=mg_mesh%mgmax
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange,minlev=minlev)
vlag_zerob%ML_vlag_rep=>ML_oft_vlagrange
IF(mg_test)THEN
  CALL lag_setup_interp(ML_oft_lagrange,ML_oft_vlagrange)
  CALL lag_mloptions
END IF
!---Run tests
oft_env%pm=.TRUE.!.FALSE.
IF(mg_test)THEN
  CALL test_mopmg
ELSE
  CALL test_mop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Invert the mass matrix and output required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: minv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_oft_vlagrange%set_level(ML_oft_lagrange%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
CALL create_diag_pre(minv%pre)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
  WRITE(io_unit,*)minv%cits
  WRITE(io_unit,*)uu/u%ng
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrix
CALL mop%delete
DEALLOCATE(mop)
!---Destory preconditioner
CALL minv%pre%delete
!---Destory solver
CALL minv%delete
END SUBROUTINE test_mop
!---------------------------------------------------------------------------------
!> Invert the mass matrix and output required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
CLASS(oft_solver), POINTER :: linv => NULL()
!--- ML structures for MG-preconditioner
type(oft_matrix_ptr), pointer :: ml_int(:) => NULL()
type(oft_matrix_ptr), pointer :: ml_vlop(:) => NULL()
type(oft_matrix_ptr), pointer :: ml_lop(:) => NULL()
INTEGER(i4), ALLOCATABLE :: levels(:)
REAL(r8), ALLOCATABLE :: df(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
!---Local variables
CLASS(oft_vector), pointer :: fvec,cvec
TYPE(oft_graph_ptr) :: graphs(3,3)
TYPE(oft_matrix_ptr) :: mats(3,3)
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_ml_fe_comp_vecspace) :: ml_vecspace
REAL(r8) :: uu
INTEGER(i4) :: i,nlevels
!------------------------------------------------------------------------------
! Create ML Matrices
!------------------------------------------------------------------------------
nlevels=ML_oft_lagrange%nlevels-minlev+1
!---Get FE operators
ALLOCATE(ml_int(nlevels-1),ml_vlop(nlevels),ml_lop(nlevels))
ALLOCATE(levels(nlevels),df(nlevels),nu(nlevels))
DO i=1,nlevels
  CALL ML_oft_vlagrange%set_level(minlev+i-1,propogate=.TRUE.)
  levels(i)=minlev+(i-1)
  df(i)=df_lop(levels(i))
  nu(i)=nu_lop(levels(i))
  !---
  NULLIFY(ml_lop(i)%M)
  CALL oft_lag_getlop(ML_oft_lagrange%current_level,ml_lop(i)%M,'zerob')
!------------------------------------------------------------------------------
! Create composite matrix
!------------------------------------------------------------------------------
  ALLOCATE(graphs(1,1)%g)
  graphs(1,1)%g%nr=ML_oft_lagrange%current_level%ne
  graphs(1,1)%g%nrg=ML_oft_lagrange%current_level%global%ne
  graphs(1,1)%g%nc=ML_oft_lagrange%current_level%ne
  graphs(1,1)%g%ncg=ML_oft_lagrange%current_level%global%ne
  graphs(1,1)%g%nnz=ML_oft_lagrange%current_level%nee
  graphs(1,1)%g%kr=>ML_oft_lagrange%current_level%kee
  graphs(1,1)%g%lc=>ML_oft_lagrange%current_level%lee
  !---
  graphs(2,2)%g=>graphs(1,1)%g
  graphs(3,3)%g=>graphs(1,1)%g
  !---Get coarse and fine vectors
  CALL ML_oft_vlagrange%vec_create(cvec)
  CALL ML_oft_vlagrange%vec_create(fvec)
  NULLIFY(ml_vlop(i)%M)
  !---Construct matrix
  CALL create_matrix(ml_vlop(i)%M,graphs,fvec,cvec)
  DEALLOCATE(graphs(1,1)%g)
!------------------------------------------------------------------------------
! Combine child matrices into composite matrix
!------------------------------------------------------------------------------
  !---Specify child matrices
  mats(1,1)%m=>ml_lop(i)%M
  mats(2,2)%m=>ml_lop(i)%M
  mats(3,3)%m=>ml_lop(i)%M
  !---Combine matrices
  CALL combine_matrices(mats,3,3,ml_vlop(i)%M)
  CALL ml_vlop(i)%M%assemble(fvec)
  !---Delete temporaries
  CALL cvec%delete
  CALL fvec%delete
  DEALLOCATE(cvec,fvec)
  !---
  IF(i>1)ml_int(i-1)%M=>ML_oft_vlagrange%interp_matrices(i)%m !oft_lagrange_ops%vinterp
END DO
CALL ML_oft_vlagrange%set_level(ML_oft_lagrange%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(linv,force_native=.TRUE.)
linv%A=>ml_vlop(nlevels)%M
linv%its=-3
linv%itplot=1
ml_vecspace%ML_FE_rep=>ML_oft_vlagrange
CALL create_mlpre(linv%pre,ml_vlop(1:nlevels),levels, &
  nlevels=nlevels,ml_vecspace=ml_vecspace, &
  bc=vlag_zerob,stype=1,df=df,nu=nu)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL vlag_zerob%apply(v)
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
!---Destroy matrix
CALL mop%delete
DEALLOCATE(mop)
!---Destory preconditioner
CALL linv%pre%delete
DEALLOCATE(linv%pre)
!---Destory solver
CALL linv%delete
END SUBROUTINE test_mopmg
END PROGRAM test_vlag
