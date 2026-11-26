!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_cell_schwarz.F90
!
!> Regression tests for the cell mode of the native Schwarz preconditioner (blocks
!! grown by layers of mesh cells by a cell overlap builder attached to the matrix), on
!! a Lagrange Laplacian on a cube, with 1 and 3 blocks per process.
!!
!! The current test cases are:
!! - Exactness: with enough layers for every block to cover the whole mesh, restricted
!!   Schwarz with direct sub-solves converges in one GMRES iteration to the reference
!!   solution, and so does weighted (partition of unity) additive Schwarz; additive
!!   Schwarz converges. With one block per process the same holds in graph mode with
!!   enough graph layers for every subdomain to cover the mesh
!! - GMRES iteration counts for 0, 1 and 2 layers with point, edge and face adjacency,
!!   next to graph overlap 0, 1 and 2 (written to `cell_schwarz.results` for checking
!!   by the driver)
!! - GMRES iteration count of weighted additive Schwarz with one layer of point
!!   adjacent cells
!! - Numeric update after the operator values change matches a fresh setup
!! - Delete and re-setup (and an update with a new pattern) reproduce the first solve
!!
!! The cube may be periodic (`ref_per` in `cube_options`), so that blocks and subdomains
!! grow across periodic seams.
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
#elif defined(HAVE_SUPERLU_DIST)
#define SUB_LU_PACK 'superd'
#else
#define SUB_LU_PACK 'none'
#define SUB_LU_MISSING
#endif
!------------------------------------------------------------------------------
PROGRAM test_cell_schwarz
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_overlap_dofs, overlap_dofs_delete, OFT_ACCESS_READ, &
  OFT_CELL_ADJ_POINT, OFT_CELL_ADJ_EDGE, OFT_CELL_ADJ_FACE
USE oft_native_la, ONLY: oft_native_matrix
USE oft_native_overlap_mod, ONLY: oft_native_overlap, native_overlap_setup, native_overlap_delete
USE oft_native_solvers, ONLY: oft_schwarz_precond, oft_native_gmres_solver
USE oft_lu, ONLY: oft_lusolver
USE oft_solver_base, ONLY: oft_solver
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
USE fem_cell_overlap, ONLY: oft_fem_cell_overlap, fem_cell_overlap_create
IMPLICIT NONE
INTEGER(i4), PARAMETER :: nres=21
INTEGER(i4) :: io_unit,ierr,nfail,ip,order=3
INTEGER(i4) :: its(nres,2)
INTEGER(i8) :: psig(4,2),gsig
INTEGER(i4), PARAMETER :: nparts_list(2)=[1,3]
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(oft_fem_cell_overlap), TARGET :: cbld
CLASS(oft_vector), POINTER :: u,v,w,x,g,g0
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
NAMELIST/test_cell_schwarz_options/order
IF(use_petsc)THEN
  WRITE(*,*)'SKIP TEST'
  STOP
END IF
#ifdef SUB_LU_MISSING
WRITE(*,*)'SKIP TEST: no direct solver package for the subdomain solves'
STOP
#endif
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_cell_schwarz_options,IOSTAT=ierr)
CLOSE(io_unit)
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
oft_env%pm=.FALSE.
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(w)
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(g)
CALL ML_oft_lagrange%vec_create(g0)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,g)
CALL lag_zerob%apply(g)
CALL g0%copy(g) ! Solvers overwrite the right-hand side with the residual
IF(oft_env%head_proc)WRITE(*,'(A,I10,A,I10)')'Global DoFs: ',lop%nrg,' cells: ',mg_mesh%mesh%global%nc
nfail=0
its=0
psig=0_i8
DO ip=1,2
  CALL run_cases(nparts_list(ip),its(:,ip),psig(:,ip))
END DO
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='cell_schwarz.results')
  WRITE(io_unit,'(I0)')nfail
  DO ip=1,2
    WRITE(*,'(A,I2,A,21I5)')'nparts=',nparts_list(ip),' results: ',its(:,ip)
    WRITE(io_unit,'(*(I0,1X))')its(:,ip)
  END DO
  WRITE(io_unit,'(*(I0,1X))')psig
  CLOSE(io_unit)
END IF
CALL oft_finalize
CONTAINS
! All cases for one number of blocks per process. Results (`res`):
! 1: layers covering the mesh, 2-3: restricted/additive iterations with those layers,
! 4-12: cell mode iterations for (point, edge, face) x layers 0-2,
! 13-15: graph overlap 0-2, 16: iterations before update, 17: after update,
! 18: after delete and re-setup, 19: after update with new pattern,
! 20: weighted additive iterations with covering layers, 21: weighted additive with
! one layer of point adjacent cells
SUBROUTINE run_cases(nparts,res,sig)
INTEGER(i4), INTENT(in) :: nparts
INTEGER(i4), INTENT(out) :: res(nres)
INTEGER(i8), INTENT(out) :: sig(4) !< Hashes of the cell partition, DoF colors, graph blocks and cell ranks
INTEGER(i4) :: kfull,ia,k,nb
REAL(r8) :: rnorm,t0
TYPE(oft_overlap_dofs) :: dofs
CHARACTER(LEN=16) :: tag
WRITE(tag,'(A,I1)')' nparts=',nparts
t0=omp_get_wtime()
CALL fem_cell_overlap_create(cbld,mg_mesh,ML_oft_lagrange%current_level,nparts)
lop%cell_ovl=>cbld
nb=cbld%nblocks()
!---Partition signature (the reference iteration tables hold for one partition only)
sig(1)=part_hash(cbld%mesh%nc,cbld%cpart,cbld%mesh%global%lc)
sig(2)=part_hash(cbld%ndofs,cbld%dcolor,row_gidx())
sig(4)=part_hash(cbld%mesh%nc,[(oft_env%rank+1,k=1,cbld%mesh%nc)],cbld%mesh%global%lc)
!---Smallest number of point layers for every block to cover the whole mesh
kfull=-1
DO k=0,64
  CALL cbld%build(lop,u,k,OFT_CELL_ADJ_POINT,-1,dofs)
  IF(-oft_mpi_max(-MINVAL(dofs%blk(1:nb)%n))==lop%nrg)THEN
    kfull=k
    EXIT
  END IF
END DO
CALL overlap_dofs_delete(dofs)
IF(kfull<0)CALL oft_abort('Blocks do not cover the mesh','run_cases',__FILE__)
res(1)=kfull
CALL timing('Covering layers',tag,t0)
!---Reference solution (graph overlap, converged tightly)
res(13)=solve(nparts,1,.FALSE.,OFT_CELL_ADJ_POINT,.FALSE.,x,rnorm,rtol=1.d-12)
!---Exactness: restricted converges in one iteration to the reference solution
res(2)=solve(nparts,kfull,.TRUE.,OFT_CELL_ADJ_POINT,.FALSE.,v,rnorm)
CALL check_exact('Exact blocks, restricted: 1 iteration'//tag,REAL(res(2),r8),1.d0)
CALL check_tol('Exact blocks, restricted: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
res(3)=solve(nparts,kfull,.TRUE.,OFT_CELL_ADJ_POINT,.TRUE.,v,rnorm)
CALL check_tol('Exact blocks, additive: converged'//tag,rnorm,1.d-9)
CALL check_tol('Exact blocks, additive: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-6)
res(20)=solve(nparts,kfull,.TRUE.,OFT_CELL_ADJ_POINT,.TRUE.,v,rnorm,weighted=.TRUE.)
CALL check_exact('Exact blocks, weighted: 1 iteration'//tag,REAL(res(20),r8),1.d0)
CALL check_tol('Exact blocks, weighted: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
!---Graph mode exactness (one subdomain per process: graph growth crosses processes)
IF(nparts==1)THEN
  kfull=graph_cover()
  k=solve(nparts,kfull,.FALSE.,OFT_CELL_ADJ_POINT,.FALSE.,v,rnorm)
  CALL check_exact('Exact graph subdomains, restricted: 1 iteration'//tag,REAL(k,r8),1.d0)
  CALL check_tol('Exact graph subdomains, restricted: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
  k=solve(nparts,kfull,.FALSE.,OFT_CELL_ADJ_POINT,.TRUE.,v,rnorm)
  CALL check_tol('Exact graph subdomains, additive: converged'//tag,rnorm,1.d-9)
  CALL check_tol('Exact graph subdomains, additive: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-6)
  kfull=res(1)
END IF
CALL timing('Exactness',tag,t0)
!---Iteration counts: (point, edge, face) x layers 0-2, graph overlap 0-2
DO ia=1,3
  DO k=0,2
    res(3+(ia-1)*3+k+1)=solve(nparts,k,.TRUE.,ia,.FALSE.,v,rnorm)
    CALL check_tol('Cell mode converged'//tag,rnorm,1.d-9)
    CALL check_tol('Cell mode solution'//tag,maxdiff(v,x)/maxabs(x),1.d-5)
  END DO
END DO
res(21)=solve(nparts,1,.TRUE.,OFT_CELL_ADJ_POINT,.TRUE.,v,rnorm,weighted=.TRUE.)
CALL check_tol('Cell mode weighted converged'//tag,rnorm,1.d-9)
CALL check_tol('Cell mode weighted solution'//tag,maxdiff(v,x)/maxabs(x),1.d-5)
CALL timing('Cell mode counts',tag,t0)
gsig=0_i8
DO k=0,2
  res(13+k)=solve(nparts,k,.FALSE.,OFT_CELL_ADJ_POINT,.FALSE.,v,rnorm)
END DO
sig(3)=gsig
CALL timing('Graph mode counts',tag,t0)
!---Numeric update and delete/re-setup
CALL check_update(nparts,res(16:19))
CALL timing('Update and re-setup',tag,t0)
NULLIFY(lop%cell_ovl)
CALL cbld%delete
END SUBROUTINE run_cases
! Smallest number of graph layers for every process' subdomain to cover the whole mesh
FUNCTION graph_cover() RESULT(kg)
INTEGER(i4) :: kg,k
TYPE(oft_native_overlap) :: ovl
CLASS(oft_native_matrix), POINTER :: lnat
SELECT TYPE(lop)
CLASS IS(oft_native_matrix)
  lnat=>lop
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','graph_cover',__FILE__)
END SELECT
kg=-1
DO k=0,128
  CALL native_overlap_setup(ovl,lnat,u,k)
  IF(-oft_mpi_max(-ovl%n)==lop%nrg)kg=k
  CALL native_overlap_delete(ovl)
  IF(kg>=0)EXIT
END DO
IF(kg<0)CALL oft_abort('Graph subdomains do not cover the mesh','graph_cover',__FILE__)
IF(oft_env%head_proc)WRITE(*,'(A,I4)')'Graph layers covering the mesh: ',kg
END FUNCTION graph_cover
! Print the wall time since the last call (or `t0` set by the caller)
SUBROUTINE timing(label,tag,t0)
CHARACTER(LEN=*), INTENT(in) :: label,tag
REAL(r8), INTENT(inout) :: t0
IF(oft_env%head_proc)WRITE(*,'(3A,F10.2,A)')'Time: ',label,tag,omp_get_wtime()-t0,' s'
t0=omp_get_wtime()
END SUBROUTINE timing
! Solve the Laplace problem with GMRES preconditioned by Schwarz (cell or graph mode)
! with direct sub-solves, returning the iteration count and relative residual
! (weighted: partition of unity additive)
FUNCTION solve(nlocal,overlap,cell,adj,additive,y,rnorm,rtol,weighted) RESULT(nits)
INTEGER(i4), INTENT(in) :: nlocal,overlap,adj
LOGICAL, INTENT(in) :: cell,additive
CLASS(oft_vector), INTENT(inout) :: y
REAL(r8), INTENT(out) :: rnorm
REAL(r8), OPTIONAL, INTENT(in) :: rtol
LOGICAL, OPTIONAL, INTENT(in) :: weighted
INTEGER(i4) :: nits
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
CALL gmres_init(linv,sw,sub,nlocal,overlap,cell,adj,additive)
IF(PRESENT(rtol))linv%rtol=rtol
IF(PRESENT(weighted))sw%weighted=weighted
CALL y%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(y,g)
nits=linv%cits
rnorm=SQRT(g%dot(g)/g0%dot(g0))
!---Graph mode blocks without overlap (partition signature)
IF((.NOT.cell).AND.(overlap==0).AND.(nlocal>1))gsig=graph_hash(sw)
CALL sw%delete
CALL linv%delete
END FUNCTION solve
! Hash of a partition (part of each local item, weighted by the item's global index), summed
! over processes: identifies the partition the reference tables in the driver were made with
FUNCTION part_hash(n,part,gidx) RESULT(h)
INTEGER(i4), INTENT(in) :: n
INTEGER(i4), INTENT(in) :: part(:)
INTEGER(i8), INTENT(in) :: gidx(:)
INTEGER(i8) :: h
INTEGER(i8), PARAMETER :: p=2147483647_i8
INTEGER(i4) :: i
h=0_i8
DO i=1,n
  h=MOD(h+INT(part(i),i8)*MOD(MOD(gidx(i),p)*48271_i8,p),p)
END DO
h=MOD(oft_mpi_sum(h),p)
END FUNCTION part_hash
! Global index of each local row of the operator
FUNCTION row_gidx() RESULT(gidx)
INTEGER(i8), ALLOCATABLE :: gidx(:)
INTEGER(i4) :: i
ALLOCATE(gidx(MAX(1,lop%nr)))
DO i=1,lop%nr
  gidx(i)=lop%i_map(1)%offsetg+lop%i_map(1)%lge(i)
END DO
END FUNCTION row_gidx
! Hash of the graph mode blocks (block of each local row, 0 for rows in no block)
FUNCTION graph_hash(sw) RESULT(h)
TYPE(oft_schwarz_precond), INTENT(in) :: sw
INTEGER(i8) :: h
INTEGER(i4) :: b
INTEGER(i4), ALLOCATABLE :: rpart(:)
ALLOCATE(rpart(MAX(1,lop%nr)))
rpart=0
DO b=1,SIZE(sw%parts)
  rpart(sw%parts(b)%v)=b
END DO
h=part_hash(lop%nr,rpart,row_gidx())
END FUNCTION graph_hash
! Set up GMRES with a Schwarz preconditioner
SUBROUTINE gmres_init(linv,sw,sub,nlocal,overlap,cell,adj,additive)
TYPE(oft_native_gmres_solver), INTENT(inout) :: linv
TYPE(oft_schwarz_precond), TARGET, INTENT(inout) :: sw
CLASS(oft_solver), TARGET, INTENT(inout) :: sub
INTEGER(i4), INTENT(in) :: nlocal,overlap,adj
LOGICAL, INTENT(in) :: cell,additive
linv%A=>lop
linv%its=-3
linv%nrits=20
linv%rtol=1.d-10
sw%nlocal=nlocal
sw%overlap=overlap
sw%cell_mode=cell
sw%adjacency=adj
sw%additive=additive
sw%pre=>sub
linv%pre=>sw
END SUBROUTINE gmres_init
! Scaling the operator and updating the preconditioner matches a fresh setup; delete
! and re-setup, and an update with a new pattern, reproduce the first solve
SUBROUTINE check_update(nparts,res)
INTEGER(i4), INTENT(in) :: nparts
INTEGER(i4), INTENT(out) :: res(4)
INTEGER(i4) :: nfresh
REAL(r8) :: rnorm
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
CHARACTER(LEN=16) :: tag
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
WRITE(tag,'(A,I1)')' nparts=',nparts
CALL gmres_init(linv,sw,sub,nparts,1,.TRUE.,OFT_CELL_ADJ_POINT,.FALSE.)
CALL v%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(v,g)
res(1)=linv%cits
!---Numeric update after scaling the operator, against a fresh setup
CALL scale_op(2.d0)
CALL linv%update(.FALSE.)
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
res(2)=linv%cits
nfresh=solve(nparts,1,.TRUE.,OFT_CELL_ADJ_POINT,.FALSE.,x,rnorm)
CALL check_exact('Updated = fresh setup (iterations)'//tag,REAL(res(2),r8),REAL(nfresh,r8))
CALL check_tol('Updated = fresh setup (solution)'//tag,maxdiff(w,x)/maxabs(x),1.d-12)
CALL w%scale(2.d0)
CALL check_tol('Updated solution scales with operator'//tag,maxdiff(v,w)/maxabs(v),1.d-8)
CALL scale_op(0.5d0)
!---Delete and re-setup on the original operator
CALL sw%delete
CALL linv%update(.FALSE.)
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
res(3)=linv%cits
CALL check_exact('Re-setup iterations'//tag,REAL(res(3),r8),REAL(res(1),r8))
CALL check_tol('Re-setup solution'//tag,maxdiff(v,w)/maxabs(v),1.d-12)
!---Update with a new pattern (rebuild)
CALL linv%update(.TRUE.)
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
res(4)=linv%cits
CALL check_exact('New pattern iterations'//tag,REAL(res(4),r8),REAL(res(1),r8))
CALL check_tol('New pattern solution'//tag,maxdiff(v,w)/maxabs(v),1.d-12)
CALL sw%delete
CALL linv%delete
END SUBROUTINE check_update
! Scale the operator values
SUBROUTINE scale_op(fac)
REAL(r8), INTENT(in) :: fac
SELECT TYPE(lop)
TYPE IS(oft_native_matrix)
  lop%M=fac*lop%M
  lop%full_current=.FALSE.
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','scale_op',__FILE__)
END SELECT
END SUBROUTINE scale_op
! Global maximum of |x_i - y_i| over all local entries
FUNCTION maxdiff(a,b) RESULT(d)
CLASS(oft_vector), INTENT(inout) :: a,b
REAL(r8) :: d
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: av,bv
NULLIFY(av,bv)
CALL a%borrow_local(av,OFT_ACCESS_READ)
CALL b%borrow_local(bv,OFT_ACCESS_READ)
d=0.d0
IF(a%n>0)d=MAXVAL(ABS(av-bv))
CALL a%release_local(av,OFT_ACCESS_READ,.FALSE.)
CALL b%release_local(bv,OFT_ACCESS_READ,.FALSE.)
d=oft_mpi_max(d)
END FUNCTION maxdiff
! Global maximum of |x_i|
FUNCTION maxabs(a) RESULT(d)
CLASS(oft_vector), INTENT(inout) :: a
REAL(r8) :: d
d=a%norm(3)
END FUNCTION maxabs
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
END PROGRAM test_cell_schwarz
