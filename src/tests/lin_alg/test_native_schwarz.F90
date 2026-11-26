!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_native_schwarz.F90
!
!> Regression tests for the native overlapping Schwarz preconditioner and the
!! overlapping subdomain machinery, on a third order Lagrange Laplacian on a cube.
!!
!! The current test cases are:
!! - Overlap 0 (one subdomain per process) reproduces block-Jacobi bit for bit
!! - Overlap 1 with several blocks per process reproduces block-Jacobi with
!!   `loverlap` bit for bit (restricted)
!! - Subdomain matrices with overlap across processes match the entries of the
!!   globally assembled matrix (assembled independently from all local matrices)
!! - Subdomains from a given DOF list: the DOFs of a graph-grown subdomain reproduce
!!   it exactly (numbering, matrix, gathered values), an arbitrary list (not all owned
!!   DOFs) matches the assembled matrix, global indices map back to subdomain indices,
!!   and submatrices of the subdomain matrix follow its values through updates
!! - Updating the preconditioner after the operator values change
!! - GMRES iteration counts with restricted, additive and weighted (partition of
!!   unity) additive Schwarz for 0, 1 and 2 overlap layers (written to
!!   `schwarz.results` for checking by the driver)
!!
!! All cases also run on partially periodic meshes, where a process holds several images
!! of the DOFs on periodic seams. On fully periodic meshes the Laplacian is singular
!! (constant null space), so only the subdomain matrices are checked.
!!
!! @ingroup testing
! Direct package for the subdomain solves (the exactness checks and iteration tables
! assume exact sub-solves; without any direct package the test is skipped)
#if defined(HAVE_SUPERLU)
#define SUB_LU_PACK 'super'
#elif defined(HAVE_UMFPACK)
#define SUB_LU_PACK 'umfpack'
#elif defined(HAVE_MKL)
#define SUB_LU_PACK 'mkl'
#define SUB_LU_MKL
#elif defined(HAVE_SUPERLU_DIST)
#define SUB_LU_PACK 'superd'
#else
#define SUB_LU_PACK 'none'
#define SUB_LU_MISSING
#endif
!------------------------------------------------------------------------------
! Concurrent direct solves (a module procedure: no OpenMP in internal procedures)
!------------------------------------------------------------------------------
MODULE test_native_schwarz_helpers
USE oft_base
USE oft_native_la, ONLY: oft_native_vector
USE oft_lu, ONLY: oft_lusolver
USE oft_stitching, ONLY: oft_seam
IMPLICIT NONE
CONTAINS
! Local vector of n values with a seam that skips all stitching
SUBROUTINE local_vec(vec,n,seam)
TYPE(oft_native_vector), INTENT(inout) :: vec
INTEGER(i4), INTENT(in) :: n
TYPE(oft_seam), TARGET, INTENT(inout) :: seam
ALLOCATE(vec%v(MAX(1,n)))
vec%n=n
vec%ng=INT(n,i8)
vec%nslice=n
vec%stitch_info=>seam
END SUBROUTINE local_vec
! Apply each solver to its own right-hand side, concurrently
SUBROUTINE solve_concurrent(n,lus,xs,bs)
INTEGER(i4), INTENT(in) :: n
TYPE(oft_lusolver), INTENT(inout) :: lus(n)
TYPE(oft_native_vector), INTENT(inout) :: xs(n),bs(n)
INTEGER(i4) :: i
!$omp parallel do schedule(static,1)
DO i=1,n
  CALL lus(i)%apply(xs(i),bs(i))
END DO
END SUBROUTINE solve_concurrent
END MODULE test_native_schwarz_helpers
!------------------------------------------------------------------------------
PROGRAM test_native_schwarz
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_stitching, ONLY: oft_seam
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph_ptr, OFT_ACCESS_READ
USE oft_la_utils, ONLY: create_matrix
USE oft_native_la, ONLY: oft_native_matrix, oft_native_vector, oft_native_submatrix, native_matrix_set_pattern
USE oft_native_overlap_mod, ONLY: oft_native_overlap, native_overlap_setup, native_overlap_setup_dofs, &
  native_overlap_update, native_overlap_gather, native_overlap_index, native_overlap_delete
USE oft_native_solvers, ONLY: oft_bjprecond, oft_schwarz_precond, oft_native_gmres_solver
USE oft_lu, ONLY: oft_lusolver
USE oft_solver_base, ONLY: oft_solver
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
USE test_native_schwarz_helpers, ONLY: solve_concurrent, local_vec
IMPLICIT NONE
#ifdef SUB_LU_MKL
EXTERNAL :: mkl_set_num_threads
#endif
INTEGER(i4) :: io_unit,nfail,its(11)
INTEGER(i4), PARAMETER :: order=3
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
CLASS(oft_vector), POINTER :: u,v,w,g,g0
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_native_matrix), POINTER :: A => NULL()
INTEGER(i4) :: nukey = 0
INTEGER(i8), ALLOCATABLE :: ukey(:)
REAL(r8), ALLOCATABLE :: uval(:)
IF(use_petsc)THEN
  WRITE(*,*)'SKIP TEST'
  STOP
END IF
#ifdef SUB_LU_MISSING
WRITE(*,*)'SKIP TEST: no direct solver package for the subdomain solves'
STOP
#endif
CALL oft_init
#ifdef SUB_LU_MKL
!---Separate factorizations of the same matrix by multithreaded MKL PARDISO are not
!   bitwise identical, but the checks against block-Jacobi compare solutions exactly
CALL mkl_set_num_threads(1)
#endif
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
oft_env%pm=.FALSE.
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(w)
CALL ML_oft_lagrange%vec_create(g)
CALL ML_oft_lagrange%vec_create(g0)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,g)
CALL lag_zerob%apply(g)
CALL g0%copy(g) ! Solvers overwrite the right-hand side with the residual
nfail=0
its=0
SELECT TYPE(lop)
TYPE IS(oft_native_matrix)
  A=>lop
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','main',__FILE__)
END SELECT
!---Fully periodic: singular operator, check subdomain matrices only
IF(mg_mesh%mesh%periodic%nper==3)THEN
  CALL assemble_global()
  CALL check_subdomain(1)
  CALL check_subdomain(2)
  CALL write_results()
  CALL oft_finalize
END IF
!---Overlap 0 is block-Jacobi on the owned rows
its(1)=solve(1,0,.FALSE.,.FALSE.,v)
its(2)=solve(1,0,.FALSE.,.TRUE.,w)
CALL check_exact('Overlap 0 iterations = block-Jacobi',REAL(its(2)-its(1),r8),0.d0)
CALL check_exact('Overlap 0 solution = block-Jacobi',maxdiff(v,w),0.d0)
!---Overlap 1 within a process is block-Jacobi with loverlap
its(1)=solve(4,1,.FALSE.,.FALSE.,v)
its(2)=solve(4,1,.FALSE.,.TRUE.,w)
CALL check_exact('Local overlap 1 iterations = loverlap',REAL(its(2)-its(1),r8),0.d0)
CALL check_exact('Local overlap 1 solution = loverlap',maxdiff(v,w),0.d0)
!---Subdomain matrices against the assembled global matrix
CALL assemble_global()
CALL check_subdomain(1)
CALL check_subdomain(2)
!---Subdomains from a given DOF list
CALL check_dofs(1)
CALL check_dofs(2)
!---Concurrent direct solves on submatrices sharing a graph
CALL check_shared_graph(SUB_LU_PACK)
#if defined(HAVE_MKL) && !defined(SUB_LU_MKL)
CALL check_shared_graph('mkl')
#endif
!---Iteration counts with overlap across processes
its(1)=solve(1,0,.FALSE.,.TRUE.,v)
its(2)=solve(1,1,.FALSE.,.TRUE.,v)
its(3)=solve(1,2,.FALSE.,.TRUE.,v)
its(4)=solve(1,1,.TRUE.,.TRUE.,v)
its(5)=solve(1,2,.TRUE.,.TRUE.,v)
its(9)=solve(1,1,.TRUE.,.TRUE.,v,weighted=.TRUE.)
its(10)=solve(1,2,.TRUE.,.TRUE.,v,weighted=.TRUE.)
!---Numeric update of the subdomain matrices after the operator values change
CALL check_update()
!---Operator rebuilt with another sparsity pattern and update with a new pattern
CALL check_new_pattern(1)
CALL check_new_pattern(4)
!---Sub-solver life cycle (setup, updates, delete and re-setup) through the Schwarz paths
CALL check_lifecycle(SUB_LU_PACK,1)
CALL check_lifecycle(SUB_LU_PACK,4)
#if defined(HAVE_SUPERLU_DIST)
IF(TRIM(SUB_LU_PACK)/='superd')THEN
  CALL check_lifecycle('superd',1)
  CALL check_lifecycle('superd',4)
END IF
#endif
!---Several blocks per process: no overlap, additive and restricted with overlap 1
its(6)=solve(4,0,.TRUE.,.TRUE.,v)
its(7)=solve(4,1,.TRUE.,.TRUE.,v)
its(8)=solve(4,1,.FALSE.,.TRUE.,v)
its(11)=solve(4,1,.TRUE.,.TRUE.,v,weighted=.TRUE.)
CALL write_results()
CALL oft_finalize
CONTAINS
! Write the number of failed checks and the iteration counts for the driver
SUBROUTINE write_results()
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,11I6)')'Iterations: ',its
  OPEN(NEWUNIT=io_unit,FILE='schwarz.results')
  WRITE(io_unit,'(I0)')nfail
  WRITE(io_unit,'(*(I0,1X))')its
  CLOSE(io_unit)
END IF
END SUBROUTINE write_results
! Solve the Laplace problem with GMRES preconditioned by block-Jacobi (schwarz=.FALSE.,
! with loverlap if overlap>0) or Schwarz (weighted: partition of unity additive),
! returning the iteration count
FUNCTION solve(nlocal,overlap,additive,schwarz,x,weighted) RESULT(nits)
INTEGER(i4), INTENT(in) :: nlocal,overlap
LOGICAL, INTENT(in) :: additive,schwarz
CLASS(oft_vector), INTENT(inout) :: x
LOGICAL, OPTIONAL, INTENT(in) :: weighted
INTEGER(i4) :: nits
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_bjprecond), TARGET :: bj
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
linv%A=>lop
linv%its=-3
linv%nrits=20
IF(schwarz)THEN
  sw%nlocal=nlocal
  sw%overlap=overlap
  sw%additive=additive
  IF(PRESENT(weighted))sw%weighted=weighted
  sw%pre=>sub
  linv%pre=>sw
ELSE
  bj%nlocal=nlocal
  bj%loverlap=(overlap>0)
  bj%pre=>sub
  linv%pre=>bj
END IF
CALL x%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(x,g)
nits=linv%cits
IF(schwarz)THEN
  CALL sw%delete
ELSE
  CALL bj%delete
END IF
END FUNCTION solve
! Scaling the operator and updating the preconditioner scales the solution
SUBROUTINE check_update()
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
linv%A=>lop
linv%its=-3
linv%nrits=20
linv%atol=1.d-12
sw%overlap=1
sw%pre=>sub
linv%pre=>sw
CALL v%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(v,g)
SELECT TYPE(lop)
TYPE IS(oft_native_matrix)
  lop%M=2.d0*lop%M
  lop%full_current=.FALSE.
END SELECT
CALL linv%update(.FALSE.)
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
CALL w%scale(2.d0)
CALL check_tol('Updated Schwarz solution scales with operator',maxdiff(v,w)/SQRT(v%dot(v)),1.d-8)
SELECT TYPE(lop)
TYPE IS(oft_native_matrix)
  lop%M=0.5d0*lop%M
  lop%full_current=.FALSE.
END SELECT
CALL sw%delete
END SUBROUTINE check_update
! Schwarz preconditioner (overlap 1, one subdomain per process across processes or 4 blocks
! per process) after its operator is rebuilt in place with another sparsity pattern and the
! solver is updated with a new pattern: the subdomains are rebuilt, and the solve matches a
! fresh setup on the new operator, first with more and then with fewer non-zeros than the
! previous pattern
SUBROUTINE check_new_pattern(nlocal)
INTEGER(i4), INTENT(in) :: nlocal
CLASS(oft_matrix), POINTER :: asmall,afull,opx
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
INTEGER(i4) :: k,nits,nfresh
CHARACTER(LEN=32) :: lab
NULLIFY(asmall,afull,opx)
CALL pattern_operator(asmall,.TRUE.)
CALL pattern_operator(afull,.FALSE.)
CALL pattern_operator(opx,.TRUE.)
sub%package=SUB_LU_PACK
linv%A=>opx
linv%its=-3
linv%nrits=20
sw%nlocal=nlocal
sw%overlap=1
sw%pre=>sub
linv%pre=>sw
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
DO k=1,2
  IF(k==1)THEN
    CALL replace_operator(opx,afull)
    WRITE(lab,'(A,I1,A)')' (more non-zeros, nlocal=',nlocal,')'
  ELSE
    CALL replace_operator(opx,asmall)
    WRITE(lab,'(A,I1,A)')' (fewer non-zeros, nlocal=',nlocal,')'
  END IF
  CALL linv%update(.TRUE.)
  CALL w%set(0.d0)
  CALL g%copy(g0)
  CALL linv%apply(w,g)
  nits=linv%cits
  nfresh=fresh_solve(opx,nlocal,v)
  CALL check_exact('New pattern = fresh setup (iterations)'//TRIM(lab),REAL(nits-nfresh,r8),0.d0)
  CALL check_tol('New pattern = fresh setup (solution)'//TRIM(lab),maxdiff(v,w)/SQRT(v%dot(v)),1.d-10)
END DO
CALL sw%delete
CALL opx%delete
CALL asmall%delete
CALL afull%delete
DEALLOCATE(opx,asmall,afull)
END SUBROUTINE check_new_pattern
! Solve with a freshly set up Schwarz preconditioner (overlap 1) on an operator
FUNCTION fresh_solve(op,nlocal,x) RESULT(nits)
CLASS(oft_matrix), POINTER, INTENT(in) :: op
INTEGER(i4), INTENT(in) :: nlocal
CLASS(oft_vector), INTENT(inout) :: x
INTEGER(i4) :: nits
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
linv%A=>op
linv%its=-3
linv%nrits=20
sw%nlocal=nlocal
sw%overlap=1
sw%pre=>sub
linv%pre=>sw
CALL x%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(x,g)
nits=linv%cits
CALL sw%delete
END FUNCTION fresh_solve
! Copy of the operator on its own pattern, optionally without the couplings between rows
! and columns whose global indices sum to a multiple of 5 (a symmetric rule, the same on
! every process), lumped onto the diagonal: removing a symmetric pair of entries a and
! adding |a| to both diagonal entries adds a positive semi-definite term, so the operator
! stays symmetric positive definite
SUBROUTINE pattern_operator(mat,drop)
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat
LOGICAL, INTENT(in) :: drop
TYPE(oft_graph_ptr) :: graphs(1,1)
INTEGER(i4) :: i,k,m,n
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=A%nr
graphs(1,1)%g%nrg=A%nrg
graphs(1,1)%g%nc=A%nc
graphs(1,1)%g%ncg=A%ncg
ALLOCATE(graphs(1,1)%g%kr(A%nr+1),graphs(1,1)%g%lc(MAX(1,A%nnz)))
n=0
graphs(1,1)%g%kr(1)=1
DO i=1,A%nr
  DO k=A%kr(i),A%kr(i+1)-1
    IF(dropped(i,A%lc(k),drop))CYCLE
    n=n+1
    graphs(1,1)%g%lc(n)=A%lc(k)
  END DO
  graphs(1,1)%g%kr(i+1)=n+1
END DO
graphs(1,1)%g%nnz=n
CALL create_matrix(mat,graphs,u,u,native=.TRUE.)
DEALLOCATE(graphs(1,1)%g%kr,graphs(1,1)%g%lc)
DEALLOCATE(graphs(1,1)%g)
!---Values of the operator on the kept entries, dropped entries lumped onto the diagonal
SELECT TYPE(mat)
TYPE IS(oft_native_matrix)
  DO i=1,A%nr
    DO k=A%kr(i),A%kr(i+1)-1
      IF(dropped(i,A%lc(k),drop))THEN
        DO m=mat%kr(i),mat%kr(i+1)-1
          IF(mat%lc(m)==i)mat%M(m)=mat%M(m)+ABS(A%M(k))
        END DO
        CYCLE
      END IF
      DO m=mat%kr(i),mat%kr(i+1)-1
        IF(mat%lc(m)==A%lc(k))mat%M(m)=mat%M(m)+A%M(k)
      END DO
    END DO
  END DO
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','pattern_operator',__FILE__)
END SELECT
END SUBROUTINE pattern_operator
! Coupling dropped from the operator by pattern_operator
FUNCTION dropped(ii,jj,drop) RESULT(skip)
INTEGER(i4), INTENT(in) :: ii,jj
LOGICAL, INTENT(in) :: drop
LOGICAL :: skip
INTEGER(i8) :: gi,gj
gi=A%i_map(1)%offsetg+A%i_map(1)%lge(ii)
gj=A%j_map(1)%offsetg+A%j_map(1)%lge(jj)
skip=drop.AND.(gi/=gj).AND.(MOD(gi+gj,5_i8)==0)
END FUNCTION dropped
! Rebuild an operator in place on the pattern and values of another one (as for an operator
! whose sparsity changes between solves)
SUBROUTINE replace_operator(dst,src)
CLASS(oft_matrix), POINTER, INTENT(inout) :: dst,src
SELECT TYPE(dst)
TYPE IS(oft_native_matrix)
  SELECT TYPE(src)
  TYPE IS(oft_native_matrix)
    CALL native_matrix_set_pattern(dst,src%pattern)
    DEALLOCATE(dst%M)
    ALLOCATE(dst%M(dst%nnz))
    dst%M=src%M(1:src%nnz)
    IF(ASSOCIATED(dst%Mfull))DEALLOCATE(dst%Mfull)
    dst%full_current=.FALSE.
  END SELECT
END SELECT
END SUBROUTINE replace_operator
! Schwarz preconditioner (overlap 1, one subdomain per process across processes or 4 blocks
! per process) through the life cycle of its sub-solvers with a given direct package: setup
! and factorization, numeric update, update with a new pattern, delete and re-setup, and a
! final delete. Each solve reproduces the first one. With SuperLU-DIST, sub-solvers deleted
! after the subdomain matrix they point to freed their factors with a wrong size.
SUBROUTINE check_lifecycle(package,nlocal)
CHARACTER(LEN=*), INTENT(in) :: package
INTEGER(i4), INTENT(in) :: nlocal
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
INTEGER(i4) :: k,nits(4)
REAL(r8) :: d(3)
CHARACTER(LEN=32) :: lab
WRITE(lab,'(3A,I1,A)')' (',TRIM(package),', nlocal=',nlocal,')'
sub%package=package
linv%A=>lop
linv%its=-3
linv%nrits=20
sw%nlocal=nlocal
sw%overlap=1
sw%pre=>sub
linv%pre=>sw
CALL v%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(v,g)
nits(1)=linv%cits
DO k=1,3
  SELECT CASE(k)
  CASE(1)
    CALL linv%update(.FALSE.)
  CASE(2)
    CALL linv%update(.TRUE.)
  CASE(3)
    !---Delete resets the number of blocks
    CALL sw%delete
    sw%nlocal=nlocal
  END SELECT
  CALL w%set(0.d0)
  CALL g%copy(g0)
  CALL linv%apply(w,g)
  nits(k+1)=linv%cits
  d(k)=maxdiff(v,w)/SQRT(v%dot(v))
END DO
CALL sw%delete
CALL check_exact('Sub-solver life cycle iterations'//TRIM(lab),REAL(MAXVAL(ABS(nits(2:4)-nits(1))),r8),0.d0)
CALL check_tol('Sub-solver life cycle solutions'//TRIM(lab),MAXVAL(d),1.d-10)
END SUBROUTINE check_lifecycle
! Global matrix assembled from the local (unassembled) matrices of all processes,
! as sorted unique keys (row-1)*nrg+col with summed values
SUBROUTINE assemble_global()
INTEGER(i4) :: i,j,k,nloc,ntot,ierr
INTEGER(i4), ALLOCATABLE :: counts(:),displs(:),perm(:)
INTEGER(i8), ALLOCATABLE :: gidx(:),lrow(:),lcol(:),grow(:),gcol(:),key(:)
REAL(r8), ALLOCATABLE :: lval(:),gval(:)
INTEGER(i8) :: ng
ALLOCATE(gidx(A%nr))
DO i=1,A%nr
  gidx(i)=A%i_map(1)%offsetg+A%i_map(1)%lge(i)
END DO
ng=A%nrg
nloc=A%nnz
ALLOCATE(lrow(nloc),lcol(nloc),lval(nloc))
DO i=1,A%nr
  DO j=A%kr(i),A%kr(i+1)-1
    lrow(j)=gidx(i)
    lcol(j)=gidx(A%lc(j))
    lval(j)=A%M(j)
  END DO
END DO
ALLOCATE(counts(oft_env%nprocs),displs(oft_env%nprocs))
#ifdef HAVE_MPI
CALL MPI_ALLGATHER(nloc,1,OFT_MPI_I4,counts,1,OFT_MPI_I4,oft_env%comm,ierr)
#else
counts=nloc
#endif
displs(1)=0
DO i=2,oft_env%nprocs
  displs(i)=displs(i-1)+counts(i-1)
END DO
ntot=SUM(counts)
ALLOCATE(grow(ntot),gcol(ntot),gval(ntot))
#ifdef HAVE_MPI
CALL MPI_ALLGATHERV(lrow,nloc,OFT_MPI_I8,grow,counts,displs,OFT_MPI_I8,oft_env%comm,ierr)
CALL MPI_ALLGATHERV(lcol,nloc,OFT_MPI_I8,gcol,counts,displs,OFT_MPI_I8,oft_env%comm,ierr)
CALL MPI_ALLGATHERV(lval,nloc,OFT_MPI_R8,gval,counts,displs,OFT_MPI_R8,oft_env%comm,ierr)
#else
grow=lrow
gcol=lcol
gval=lval
#endif
!---Sum duplicate entries (key = row*ng + col)
ALLOCATE(key(ntot),perm(ntot))
DO k=1,ntot
  key(k)=(grow(k)-1)*ng+gcol(k)
  perm(k)=k
END DO
CALL sort_array(key,perm,ntot)
ALLOCATE(ukey(ntot),uval(ntot))
j=0
DO k=1,ntot
  IF(j>0)THEN
    IF(key(k)==ukey(j))THEN
      uval(j)=uval(j)+gval(perm(k))
      CYCLE
    END IF
  END IF
  j=j+1
  ukey(j)=key(k)
  uval(j)=gval(perm(k))
END DO
nukey=j
END SUBROUTINE assemble_global
! Subdomain matrix entries and pattern against the assembled global matrix
SUBROUTINE check_entries(ovl,label)
TYPE(oft_native_overlap), INTENT(in) :: ovl
CHARACTER(LEN=*), INTENT(in) :: label
INTEGER(i4) :: i,k,nmiss,ir,ic
INTEGER(i4), ALLOCATABLE :: sp(:)
INTEGER(i8), ALLOCATABLE :: sg(:)
REAL(r8) :: dmax,ref
INTEGER(i8) :: ng,p8
ng=A%nrg
!---Every subdomain entry matches the global matrix
dmax=0.d0
DO i=1,ovl%n
  DO k=ovl%amat%kr(i),ovl%amat%kr(i+1)-1
    p8=search_sorted((ovl%gsub(i)-1)*ng+ovl%gsub(ovl%amat%lc(k)),ukey,INT(nukey,i8))
    ref=0.d0
    IF(p8>0)ref=uval(p8)
    dmax=MAX(dmax,ABS(ovl%amat%M(k)-ref)/MAX(1.d0,ABS(ref)))
  END DO
END DO
CALL check_tol('Subdomain entries match assembled matrix ('//label//')',oft_mpi_max(dmax),1.d-12)
!---Every non-zero global entry within the subdomain is in the subdomain matrix
ALLOCATE(sg(MAX(1,ovl%n)),sp(MAX(1,ovl%n)))
DO i=1,ovl%n
  sg(i)=ovl%gsub(i)
  sp(i)=i
END DO
IF(ovl%n>0)CALL sort_array(sg,sp,ovl%n)
nmiss=0
DO k=1,nukey
  IF(uval(k)==0.d0)CYCLE
  p8=search_sorted((ukey(k)-1)/ng+1,sg,INT(ovl%n,i8))
  IF(p8==0)CYCLE
  ir=sp(p8)
  p8=search_sorted(MOD(ukey(k)-1,ng)+1,sg,INT(ovl%n,i8))
  IF(p8==0)CYCLE
  ic=sp(p8)
  IF(.NOT.ANY(ovl%amat%lc(ovl%amat%kr(ir):ovl%amat%kr(ir+1)-1)==ic))nmiss=nmiss+1
END DO
CALL check_exact('Subdomain pattern contains all entries ('//label//')',REAL(oft_mpi_sum(nmiss),r8),0.d0)
END SUBROUTINE check_entries
! Subdomain matrix with overlap across processes against the assembled global matrix
SUBROUTINE check_subdomain(levels)
INTEGER(i4), INTENT(in) :: levels
TYPE(oft_native_overlap) :: ovl
CALL native_overlap_setup(ovl,A,u,levels)
CALL check_entries(ovl,'overlap '//CHAR(48+levels))
!---Subdomain contains the owned rows
CALL check_exact('Subdomain grows the owned rows (overlap '//CHAR(48+levels)//')', &
  REAL(oft_mpi_sum(COUNT(ovl%osub(1:ovl%nloc)==oft_env%rank)-A%nrslice),r8),0.d0)
CALL native_overlap_delete(ovl)
END SUBROUTINE check_subdomain
! Subdomains from a given DOF list: reproduce a graph-grown subdomain from its DOFs,
! check an arbitrary list, the global to subdomain index map and submatrices of the
! subdomain matrix through a numeric update
SUBROUTINE check_dofs(levels)
INTEGER(i4), INTENT(in) :: levels
TYPE(oft_native_overlap) :: ovl,ovd
TYPE(oft_native_submatrix) :: sub
TYPE(oft_seam), TARGET :: lseam
TYPE(oft_native_vector), TARGET :: uskip
INTEGER(i4) :: i,k,n,m,nbad
INTEGER(i4), ALLOCATABLE :: sp(:),so(:),idx(:),rows(:)
INTEGER(i8), ALLOCATABLE :: sg(:),gq(:)
REAL(r8), ALLOCATABLE :: x(:),xs(:),xd(:),m0(:)
CHARACTER(LEN=9) :: lab
lab='overlap '//CHAR(48+levels)
CALL native_overlap_setup(ovl,A,u,levels)
!---Sorted DOF list of the graph-grown subdomain
n=ovl%n
ALLOCATE(sg(MAX(1,n)),sp(MAX(1,n)),so(MAX(1,n)))
DO i=1,n
  sg(i)=ovl%gsub(i)
  sp(i)=i
END DO
IF(n>0)CALL sort_array(sg,sp,n)
DO i=1,n
  so(i)=ovl%osub(sp(i))
END DO
CALL native_overlap_setup_dofs(ovd,A,u,n,sg,so)
!---Same numbering, pattern and values (exact)
nbad=0
IF((ovd%n/=ovl%n).OR.(ovd%nloc/=ovl%nloc).OR.(ovd%nown/=ovl%nown))THEN
  nbad=1
ELSE IF(ovd%amat%nnz/=ovl%amat%nnz)THEN
  nbad=1
ELSE
  nbad=nbad+COUNT(ovd%gsub/=ovl%gsub)+COUNT(ovd%osub/=ovl%osub)+COUNT(ovd%lsub/=ovl%lsub)
  nbad=nbad+COUNT(ovd%amat%kr/=ovl%amat%kr)
  nbad=nbad+COUNT(ovd%amat%lc(1:ovl%amat%nnz)/=ovl%amat%lc(1:ovl%amat%nnz))
  nbad=nbad+COUNT(ovd%amat%M(1:ovl%amat%nnz)/=ovl%amat%M(1:ovl%amat%nnz))
END IF
CALL check_exact('DOF list reproduces graph subdomain ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
!---Gathered values: global index of each DOF, identical for both subdomains
ALLOCATE(x(MAX(1,A%nr)),xs(MAX(1,n)),xd(MAX(1,n)))
DO i=1,A%nr
  x(i)=REAL(A%i_map(1)%offsetg+A%i_map(1)%lge(i),r8)
END DO
xs=0.d0
xd=-1.d0
CALL native_overlap_gather(ovl,x,xs)
CALL native_overlap_gather(ovd,x,xd)
nbad=0
IF(ovd%n==n)THEN
  DO i=1,n
    IF((xs(i)/=xd(i)).OR.(xd(i)/=REAL(ovd%gsub(i),r8)))nbad=nbad+1
  END DO
END IF
CALL check_exact('DOF list gathers identical values ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
!---Global to subdomain index map (and absent indices)
ALLOCATE(gq(n+2),idx(n+2))
gq(1:n)=ovd%gsub(1:n)
gq(n+1)=0_i8
gq(n+2)=A%nrg+1
CALL native_overlap_index(ovd,n+2,gq,idx)
nbad=COUNT(idx(n+1:n+2)/=0)
DO i=1,n
  IF(idx(i)/=i)nbad=nbad+1
END DO
CALL check_exact('Subdomain index map round trip ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
!---Submatrix of all rows of the subdomain matrix, also after a numeric update
lseam%skip=.TRUE.
uskip%n=ovd%n
uskip%stitch_info=>lseam
ALLOCATE(rows(n))
rows=[(i,i=1,n)]
CALL sub%setup(ovd%amat,uskip,rows=rows)
nbad=0
IF(sub%nnz/=ovd%amat%nnz)THEN
  nbad=1
ELSE
  nbad=COUNT(sub%M(1:sub%nnz)/=ovd%amat%M(1:sub%nnz))
END IF
CALL check_exact('Submatrix of subdomain matrix ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
ALLOCATE(m0(MAX(1,ovd%amat%nnz)))
m0=ovd%amat%M
A%M=2.d0*A%M
A%full_current=.FALSE.
CALL native_overlap_update(ovd)
CALL sub%update_slice
nbad=0
IF(sub%nnz==ovd%amat%nnz)THEN
  nbad=COUNT(sub%M(1:sub%nnz)/=ovd%amat%M(1:sub%nnz))+COUNT(ovd%amat%M(1:sub%nnz)/=2.d0*m0(1:sub%nnz))
END IF
CALL check_exact('Submatrix follows subdomain update ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
A%M=0.5d0*A%M
A%full_current=.FALSE.
CALL sub%delete
CALL native_overlap_delete(ovd)
!---Arbitrary list: the graph subdomain without every third owned DOF
m=0
DO i=1,n
  IF((so(i)==oft_env%rank).AND.(MOD(sg(i),3_i8)==0))CYCLE
  m=m+1
  sg(m)=sg(i)
  so(m)=so(i)
END DO
CALL native_overlap_setup_dofs(ovd,A,u,m,sg,so)
CALL check_entries(ovd,'DOF list, '//lab)
nbad=0
IF(ovd%n/=m)THEN
  nbad=1
ELSE
  DO k=1,m
    IF(ovd%gsub(ovd%gperm(k))/=sg(k))nbad=nbad+1
  END DO
END IF
CALL check_exact('Subdomain is the given DOF list ('//lab//')',REAL(oft_mpi_sum(nbad),r8),0.d0)
CALL native_overlap_delete(ovd)
CALL native_overlap_delete(ovl)
END SUBROUTINE check_dofs
! Binary search in a sorted int8 array
FUNCTION search_sorted(item,list,n) RESULT(pos)
INTEGER(i8), INTENT(in) :: item,n
INTEGER(i8), INTENT(in) :: list(:)
INTEGER(i8) :: pos,ilow,ihigh,imid
pos=0
ilow=1
ihigh=n
DO WHILE(ilow<=ihigh)
  imid=(ilow+ihigh)/2
  IF(list(imid)==item)THEN
    pos=imid
    RETURN
  ELSE IF(list(imid)<item)THEN
    ilow=imid+1
  ELSE
    ihigh=imid-1
  END IF
END DO
END FUNCTION search_sorted
! Global maximum of |x_i - y_i| over all local entries
! Direct solves of submatrices that share a graph (the same rows, so the same cached
! extraction), factored and applied concurrently, match a solve of the same system on
! its own (MKL PARDISO modifies the graph it is given while it runs)
SUBROUTINE check_shared_graph(package)
CHARACTER(LEN=*), INTENT(in) :: package
INTEGER(i4), PARAMETER :: ncopy=4
TYPE(oft_native_submatrix), TARGET :: subs(ncopy)
TYPE(oft_lusolver), TARGET :: lus(ncopy),lref
TYPE(oft_native_vector), TARGET :: xs(ncopy),bs(ncopy),xref,bref
TYPE(oft_seam), TARGET :: lseam
INTEGER(i4) :: i,j,nr,nbad
REAL(r8) :: d
lref%package=package
IF(.NOT.lref%check_thread())RETURN
DO i=1,ncopy
  CALL subs(i)%setup(A,u)
END DO
nr=subs(1)%nr
nbad=0
DO i=2,ncopy
  IF(.NOT.(ASSOCIATED(subs(i)%kr,subs(1)%kr).AND.ASSOCIATED(subs(i)%lc,subs(1)%lc)))nbad=nbad+1
END DO
CALL check_exact('Submatrices with the same rows share a graph',REAL(oft_mpi_sum(nbad),r8),0.d0)
!---Reference: one solve on its own
lseam%skip=.TRUE.
CALL local_vec(xref,nr,lseam)
CALL local_vec(bref,nr,lseam)
xref%v=0.d0
DO j=1,nr
  bref%v(j)=1.d0+REAL(MOD(j,7),r8)
END DO
lref%A=>subs(1)
CALL lref%apply(xref,bref)
!---All copies factored and solved at the same time
DO i=1,ncopy
  CALL local_vec(xs(i),nr,lseam)
  CALL local_vec(bs(i),nr,lseam)
  xs(i)%v=0.d0
  DO j=1,nr
    bs(i)%v(j)=1.d0+REAL(MOD(j,7),r8)
  END DO
  lus(i)%package=package
  lus(i)%A=>subs(i)
END DO
CALL solve_concurrent(ncopy,lus,xs,bs)
d=0.d0
DO i=1,ncopy
  IF(nr>0)d=MAX(d,MAXVAL(ABS(xs(i)%v-xref%v))/MAX(MAXVAL(ABS(xref%v)),1.d-300))
END DO
CALL check_tol('Concurrent solves on a shared graph ('//TRIM(package)//')',oft_mpi_max(d),1.d-10)
CALL lref%delete
DEALLOCATE(xref%v,bref%v)
DO i=1,ncopy
  CALL lus(i)%delete
  CALL subs(i)%delete
  DEALLOCATE(xs(i)%v,bs(i)%v)
END DO
END SUBROUTINE check_shared_graph
FUNCTION maxdiff(x,y) RESULT(d)
CLASS(oft_vector), INTENT(inout) :: x,y
REAL(r8) :: d
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: xv,yv
NULLIFY(xv,yv)
CALL x%borrow_local(xv,OFT_ACCESS_READ)
CALL y%borrow_local(yv,OFT_ACCESS_READ)
d=0.d0
IF(x%n>0)d=MAXVAL(ABS(xv-yv))
CALL x%release_local(xv,OFT_ACCESS_READ,.FALSE.)
CALL y%release_local(yv,OFT_ACCESS_READ,.FALSE.)
d=oft_mpi_max(d)
END FUNCTION maxdiff
SUBROUTINE check_exact(name,val,ref)
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), INTENT(in) :: val,ref
IF(val/=ref)THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,'(3A,2ES25.16)')'FAILED: ',name,' ',val,ref
ELSE
  IF(oft_env%head_proc)WRITE(*,'(3A,ES25.16)')'PASSED: ',name,' ',val
END IF
END SUBROUTINE check_exact
SUBROUTINE check_tol(name,val,tol)
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), INTENT(in) :: val,tol
IF(val>tol)THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,'(3A,2ES25.16)')'FAILED: ',name,' ',val,tol
ELSE
  IF(oft_env%head_proc)WRITE(*,'(3A,ES25.16)')'PASSED: ',name,' ',val
END IF
END SUBROUTINE check_tol
END PROGRAM test_native_schwarz
