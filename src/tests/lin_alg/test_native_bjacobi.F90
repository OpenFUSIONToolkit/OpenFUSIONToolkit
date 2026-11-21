!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_native_bjacobi.F90
!
!> Regression tests for native Block-Jacobi solver. Tests are performed
!! on a unit cube at different polynomial orders using several LU packages.
!!
!! The current test cases are:
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$
!! - Solve the same system with a direct solve of the full distributed system (`dist_lu`),
!!   then check refactorization against a fresh factorization after small and large changes
!!   in matrix values (same pattern) and, for distributed runs, after changing the sparsity pattern
!!
!! @authors Chris Hansen
!! @date May 2014
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_native_bjacobi
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---LA imports
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph_ptr
USE oft_la_utils, ONLY: create_matrix
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
LOGICAL :: dist_lu = .FALSE.
NAMELIST/test_bj_options/nlocal,sol_type,use_ilu,dist_lu
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
#if !defined( HAVE_STRUMPACK )
IF(sol_type==6)THEN
  WRITE(*,*)'SKIP TEST: 6'
  CALL oft_finalize
END IF
#elif !defined( HAVE_STRUMPACK_MPI )
IF((sol_type==6).AND.dist_lu.AND.(oft_env%nprocs>1))THEN
  WRITE(*,*)'SKIP TEST: 6 (distributed)'
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
!---------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_lap(nlocal,sol_type)
INTEGER(i4), INTENT(in) :: nlocal,sol_type
!---Create solver objects
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_bjprecond), TARGET :: linv_pre
TYPE(oft_lusolver), TARGET :: linv_pre_lu
TYPE(oft_ilusolver), TARGET :: linv_pre_ilu
!---Local variables
INTEGER(i4) :: i,j,ierr
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_native_matrix), POINTER :: lop_native
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: lop2 => NULL()
CLASS(oft_fem_type), POINTER :: lag_rep
TYPE(oft_graph_ptr) :: graphs(1,1)
!---Create solver fields
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
!---Get FE operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
SELECT TYPE(this=>ML_oft_lagrange%current_level)
  CLASS IS(oft_fem_type)
    lag_rep=>this
  CLASS DEFAULT
    CALL oft_abort("Incorrect FE type","test_lap",__FILE__)
END SELECT
SELECT TYPE(this=>lop)
  CLASS IS(oft_native_matrix)
    lop_native=>this
  CLASS DEFAULT
    CALL oft_abort("Incorrect matrix type","test_lap",__FILE__)
END SELECT
!---Setup matrix solver
linv%A=>lop
linv%its=-3
linv%itplot=1
linv%nrits=20
linv%pre=>linv_pre
IF(dist_lu)linv%pre=>linv_pre_lu
!---Setup preconditioner
IF(nlocal>0)THEN
  linv_pre%nlocal=nlocal
ELSE
  linv_pre%nlocal=1
  ALLOCATE(lop_native%color(lop_native%nr))
  CALL fem_partition(lag_rep,lop_native%color,ABS(nlocal))
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
  linv_pre_lu%distributed=dist_lu
  !---Distributed factorization must perform full setup even if not requested
  IF(dist_lu.AND.(oft_env%nprocs>1))CALL linv_pre_lu%update(.FALSE.)
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
    CASE(6)
      linv_pre_lu%package='strumpack'
    CASE DEFAULT
      CALL oft_abort("Unknown LU solver package","test_lap",__FILE__)
  END SELECT
END IF
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
uu=SQRT(u%dot(u))
!---Check refactorization against a fresh factorization of the same matrix
IF(dist_lu)THEN
  CALL u%set(1.d0)
  CALL mop%apply(u,v)
  CALL lag_zerob%apply(v)
  !---Same pattern: slightly change a symmetric subset of off-diagonal entries
  CALL scale_entries(lag_rep,lop_native,11,1.05d0)
  CALL linv_pre_lu%update(.FALSE.)
  CALL check_lu(linv_pre_lu,lop,v)
  !---Same pattern: amplify a symmetric subset of off-diagonal entries to change pivoting
  !---(STRUMPACK reuses its matching/scaling for same-pattern refactors, so large value
  !---changes require a full factorization)
  CALL scale_entries(lag_rep,lop_native,7,20.d0)
  CALL linv_pre_lu%update(sol_type==6)
  CALL check_lu(linv_pre_lu,lop,v)
  !---New pattern: drop a symmetric subset of off-diagonal entries (distributed only), which
  !---must be detected automatically as local factorizations require `update(.TRUE.)`
  IF(oft_env%nprocs>1)THEN
    CALL drop_entries(lag_rep,lop_native,graphs(1,1),lop2)
    linv_pre_lu%A=>lop2
    CALL linv_pre_lu%update(.FALSE.)
    CALL check_lu(linv_pre_lu,lop2,v)
  END IF
END IF
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='bjacobi.results')
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
IF(ASSOCIATED(lop2))THEN
  CALL lop2%delete
  DEALLOCATE(lop2,graphs(1,1)%g%kr,graphs(1,1)%g%lc)
  DEALLOCATE(graphs(1,1)%g)
END IF
!---Destory preconditioner
CALL linv_pre%delete
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lap
!---------------------------------------------------------------------------------
! Is entry (i,j) in a symmetric subset of off-diagonal entries chosen by global index?
!---------------------------------------------------------------------------------
LOGICAL FUNCTION select_entry(fe,i,j,m)
CLASS(oft_fem_type), INTENT(in) :: fe
INTEGER(i4), INTENT(in) :: i,j,m
select_entry=(i/=j).AND.(MOD(fe%global%le(i)+fe%global%le(j),INT(m,8))==0)
END FUNCTION select_entry
!---------------------------------------------------------------------------------
! Scale a symmetric subset of off-diagonal entries (see select_entry)
!---------------------------------------------------------------------------------
SUBROUTINE scale_entries(fe,A,m,factor)
CLASS(oft_fem_type), INTENT(in) :: fe
CLASS(oft_native_matrix), INTENT(inout) :: A
INTEGER(i4), INTENT(in) :: m
REAL(r8), INTENT(in) :: factor
INTEGER(i4) :: i,j
DO i=1,A%nr
  DO j=A%kr(i),A%kr(i+1)-1
    IF(select_entry(fe,i,A%lc(j),m))A%M(j)=factor*A%M(j)
  END DO
END DO
END SUBROUTINE scale_entries
!---------------------------------------------------------------------------------
! Create a copy of a matrix with a symmetric subset of off-diagonal entries removed
!---------------------------------------------------------------------------------
SUBROUTINE drop_entries(fe,A,graph,Anew)
CLASS(oft_fem_type), INTENT(inout) :: fe
CLASS(oft_native_matrix), INTENT(in) :: A
TYPE(oft_graph_ptr), INTENT(inout) :: graph ! Must persist until `Anew` is deleted
CLASS(oft_matrix), POINTER, INTENT(out) :: Anew
INTEGER(i4) :: i,j,k
CLASS(oft_vector), POINTER :: tmp_vec
ALLOCATE(graph%g)
graph%g%nr=A%nr
graph%g%nrg=fe%global%ne
graph%g%nc=A%nc
graph%g%ncg=fe%global%ne
ALLOCATE(graph%g%kr(A%nr+1),graph%g%lc(A%nnz))
k=0
DO i=1,A%nr
  graph%g%kr(i)=k+1
  DO j=A%kr(i),A%kr(i+1)-1
    IF(select_entry(fe,i,A%lc(j),5))CYCLE
    k=k+1
    graph%g%lc(k)=A%lc(j)
  END DO
END DO
graph%g%kr(A%nr+1)=k+1
graph%g%nnz=k
CALL fe%vec_create(tmp_vec)
CALL create_matrix(Anew,RESHAPE([graph],[1,1]),tmp_vec,tmp_vec)
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)
SELECT TYPE(this=>Anew)
  CLASS IS(oft_native_matrix)
    IF(this%nnz/=k)CALL oft_abort("Unexpected matrix size","drop_entries",__FILE__)
    DO i=1,A%nr
      k=this%kr(i)
      DO j=A%kr(i),A%kr(i+1)-1
        IF(select_entry(fe,i,A%lc(j),5))CYCLE
        IF(this%lc(k)/=A%lc(j))CALL oft_abort("Unexpected matrix pattern","drop_entries",__FILE__)
        this%M(k)=A%M(j)
        k=k+1
      END DO
    END DO
  CLASS DEFAULT
    CALL oft_abort("Incorrect matrix type","drop_entries",__FILE__)
END SELECT
END SUBROUTINE drop_entries
!---------------------------------------------------------------------------------
! Check the solution from a (re)factored LU solver against a fresh factorization
!---------------------------------------------------------------------------------
SUBROUTINE check_lu(lu,A,b)
CLASS(oft_lusolver), INTENT(inout) :: lu
CLASS(oft_matrix), POINTER, INTENT(in) :: A
CLASS(oft_vector), INTENT(inout) :: b
TYPE(oft_lusolver) :: lu_ref
CLASS(oft_vector), POINTER :: x,x_ref
REAL(r8) :: err,xnorm
CALL b%new(x)
CALL b%new(x_ref)
CALL lu%apply(x,b)
lu_ref%A=>A
lu_ref%package=lu%package
lu_ref%distributed=lu%distributed
CALL lu_ref%apply(x_ref,b)
xnorm=SQRT(x_ref%dot(x_ref))
CALL x%add(1.d0,-1.d0,x_ref)
err=SQRT(x%dot(x))
IF(oft_env%head_proc)WRITE(*,'(A,ES11.3)')'Refactorization relative error = ',err/xnorm
IF(err>1.d-8*xnorm)CALL oft_abort("Solution incorrect after refactorization","check_lu",__FILE__)
CALL lu_ref%delete
CALL x%delete
CALL x_ref%delete
DEALLOCATE(x,x_ref)
END SUBROUTINE check_lu
END PROGRAM test_native_bjacobi
