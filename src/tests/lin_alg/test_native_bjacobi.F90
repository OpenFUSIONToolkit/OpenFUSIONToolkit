!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_native_bjacobi.F90
!
!> Regression tests for native Block-Jacobi solver. Tests are performed
!! on a unit cube at different polynomial orders using SuperLU.
!!
!! The current test cases are:
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$
!!
!! @authors Chris Hansen
!! @date May 2014
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_native_bjacobi
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---LA imports
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_native_la, ONLY: oft_native_matrix
USE oft_native_solvers, ONLY: oft_bjprecond, oft_native_gmres_solver
USE oft_lu, ONLY: oft_lusolver, oft_ilusolver
!---FE imports
USE fem_base, ONLY: oft_fem_type, oft_ml_fem_type
USE fem_utils, ONLY: fem_partition
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
INTEGER(i4), PARAMETER :: order=3
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
INTEGER(i4) :: nlocal = 1
INTEGER(i4) :: sol_type = 1
LOGICAL :: use_ilu = .FALSE.
NAMELIST/test_bj_options/nlocal,sol_type,use_ilu
IF(use_petsc)THEN
  WRITE(*,*)'SKIP TEST'
  STOP
END IF
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_bj_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check available solvers
#if !defined( HAVE_SUPERLU )
IF(sol_type==1)THEN
  WRITE(*,*)'SKIP TEST: 1'
  CALL oft_finalize
END IF
#endif
#if !defined( HAVE_SUPERLU_DIST )
IF(sol_type==2)THEN
  WRITE(*,*)'SKIP TEST: 2'
  CALL oft_finalize
END IF
#endif
#if !defined( HAVE_MUMPS )
IF(sol_type==3)THEN
  WRITE(*,*)'SKIP TEST: 3'
  CALL oft_finalize
END IF
#endif
#if !defined( HAVE_UMFPACK )
IF(sol_type==4)THEN
  WRITE(*,*)'SKIP TEST: 4'
  CALL oft_finalize
END IF
#endif
#if !defined( HAVE_MKL )
IF(sol_type==5)THEN
  WRITE(*,*)'SKIP TEST: 5'
  CALL oft_finalize
END IF
#endif
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
!---Run tests
oft_env%pm=.TRUE. !.FALSE.
CALL test_lap(nlocal,sol_type)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!------------------------------------------------------------------------------
SUBROUTINE test_lap(nlocal,sol_type)
INTEGER(i4), INTENT(in) :: nlocal,sol_type
!---Create solver objects
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_bjprecond), TARGET :: linv_pre
TYPE(oft_lusolver), TARGET :: linv_pre_lu
TYPE(oft_ilusolver), TARGET :: linv_pre_ilu
!---Local variables
INTEGER(i4) :: ierr
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_native_matrix), POINTER :: lop_native
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_fem_type), POINTER :: lag_rep
!---Create solver fields
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
!---Setup matrix solver
linv%A=>lop
linv%its=-3
linv%itplot=1
linv%nrits=20
linv%pre=>linv_pre
!---Setup preconditioner
IF(nlocal>0)THEN
  linv_pre%nlocal=nlocal
ELSE
  linv_pre%nlocal=1
  SELECT TYPE(this=>ML_oft_lagrange%current_level)
    CLASS IS(oft_fem_type)
      lag_rep=>this
    CLASS DEFAULT
      CALL oft_abort("Incorrect FE type","test_lap",__FILE__)
  END SELECT
  SELECT TYPE(this=>lop)
    CLASS IS(oft_native_matrix)
      ALLOCATE(this%color(this%nr))
      CALL fem_partition(lag_rep,this%color,ABS(nlocal))
    CLASS DEFAULT
      CALL oft_abort("Incorrect matrix type","test_lap",__FILE__)
  END SELECT
END IF
IF(use_ilu)THEN
  linv_pre%pre=>linv_pre_ilu
  SELECT CASE(sol_type)
    CASE(0)
      linv_pre_ilu%package='native'
    ! CASE(1)
    !   linv_pre_ilu%package='super'
    CASE(5)
      linv_pre_ilu%package='mkl'
    CASE DEFAULT
      CALL oft_abort("Unknown ILU solver package","test_lap",__FILE__)
  END SELECT
ELSE
  linv_pre%pre=>linv_pre_lu
  SELECT CASE(sol_type)
    CASE(1)
      linv_pre_lu%package='super'
    CASE(2)
      linv_pre_lu%package='superd'
    CASE(3)
      linv_pre_lu%package='mumps'
    CASE(4)
      linv_pre_lu%package='umfpack'
    CASE(5)
      linv_pre_lu%package='mkl'
    CASE DEFAULT
      CALL oft_abort("Unknown LU solver package","test_lap",__FILE__)
  END SELECT
END IF
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob%apply(v)
CALL u%set(0.d0)
WRITE(*,*)'Solve In'
CALL linv%apply(u,v)
WRITE(*,*)'Solve Out'
uu=SQRT(u%dot(u))
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='bjacobi.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
WRITE(*,*)'Cleanup Vec'
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
WRITE(*,*)'Cleanup Mat'
!---Destroy matrices
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
WRITE(*,*)'Cleanup Pre'
!---Destory preconditioner
CALL linv_pre%delete
WRITE(*,*)'Cleanup Solver'
!---Destory solver
CALL linv%delete
WRITE(*,*)'Done'
END SUBROUTINE test_lap
END PROGRAM test_native_bjacobi
