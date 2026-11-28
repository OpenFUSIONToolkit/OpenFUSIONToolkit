!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_cell_schwarz_fields.F90
!
!> Regression tests for field based blocks (`nlocal=-1` and field groups `nlocal<-1`)
!! in the cell mode of the native Schwarz preconditioner, on a coupled three field
!! system on a cube (vector Lagrange representation), with 1 and 3 blocks per process.
!!
!! The operator is `C (x) L`, where `L` is the Lagrange Laplacian (Dirichlet boundary
!! rows) and `C` a symmetric positive definite 3x3 coupling matrix, so every field is
!! coupled to every other one. Field based blocks drop the coupling between groups
!! within a block, so they are not exact even when they cover the whole mesh: the
!! preconditioned operator is then `(Cg^-1 C) (x) I`, with `Cg` the block diagonal part
!! of `C` by group, and GMRES converges in at most 3 iterations.
!!
!! The current test cases are:
!! - Covering blocks: restricted (one group per field and groups 1,1,2) and weighted
!!   additive (one group per field) converge in at most 3 GMRES iterations to the
!!   reference solution
!! - GMRES iteration counts for 0, 1 and 2 layers of point adjacent cells without
!!   field groups, with one group per field and with groups 1,1,2 (written to
!!   `cell_schwarz_fields.results` for checking by the driver)
!! - Groups 1,1,1 out of 2 (group 2 empty) reproduce the blocks without field groups
!! - Without `groups`, `nlocal<-1` uses the builder's default field groups (here 1,1,2)
!! - Graph mode field based blocks (one group per field, groups 1,1,2) converge
!! - Number of blocks is the number of builder blocks times the number of groups
!! - Additive and weighted additive with one group per field converge
!! - Numeric update after the operator values change matches a fresh setup, and delete
!!   and re-setup (both keep the field groups) and an update with a new pattern
!!   reproduce the first solve
!!
!! @ingroup testing
! Direct package for the subdomain solves (without any direct package the test is skipped)
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
PROGRAM test_cell_schwarz_fields
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_overlap_dofs, overlap_dofs_delete, OFT_ACCESS_READ, &
  OFT_CELL_ADJ_POINT
USE oft_native_la, ONLY: oft_native_matrix
USE oft_native_solvers, ONLY: oft_schwarz_precond, oft_native_gmres_solver
USE oft_lu, ONLY: oft_lusolver
USE oft_solver_base, ONLY: oft_solver
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
USE fem_cell_overlap, ONLY: oft_fem_cell_overlap, fem_cell_overlap_create
IMPLICIT NONE
INTEGER(i4), PARAMETER :: nres=21,nf=3
INTEGER(i4) :: io_unit,ierr,nfail,ip,order=2
INTEGER(i4) :: its(nres,2)
INTEGER(i4), PARAMETER :: nparts_list(2)=[1,3]
!---Coupling between the fields (symmetric, diagonally dominant)
REAL(r8), PARAMETER :: cmat(nf,nf)=RESHAPE([1.d0,0.4d0,0.2d0, 0.4d0,1.d0,0.3d0, 0.2d0,0.3d0,1.d0],[nf,nf])
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_lag
TYPE(oft_ml_fem_comp_type), TARGET :: ML_vlag
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(oft_fem_cell_overlap), TARGET :: cbld
CLASS(oft_vector), POINTER :: v,w,x,g,g0
CLASS(oft_vector), POINTER :: us,gs
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: aop => NULL()
NAMELIST/test_cell_schwarz_fields_options/order
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
READ(io_unit,test_cell_schwarz_fields_options,IOSTAT=ierr)
CLOSE(io_unit)
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,order,ML_lag,ML_vlag_obj=ML_vlag,minlev=-1)
lag_zerob%ML_lag_rep=>ML_lag
oft_env%pm=.FALSE.
!---Scalar operators and right-hand side (mass of a constant, zero on the boundary)
CALL ML_lag%vec_create(us)
CALL ML_lag%vec_create(gs)
CALL oft_lag_getlop(ML_lag%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_lag%current_level,mop,'none')
CALL us%set(1.d0)
CALL mop%apply(us,gs)
CALL lag_zerob%apply(gs)
!---Coupled operator C (x) L and right-hand side (field f: f times the scalar one)
CALL ML_vlag%vec_create(v)
CALL ML_vlag%vec_create(w)
CALL ML_vlag%vec_create(x)
CALL ML_vlag%vec_create(g)
CALL ML_vlag%vec_create(g0)
CALL build_operator()
CALL set_rhs()
nfail=0
its=0
DO ip=1,2
  CALL run_cases(nparts_list(ip),its(:,ip))
END DO
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='cell_schwarz_fields.results')
  WRITE(io_unit,'(I0)')nfail
  DO ip=1,2
    WRITE(io_unit,'(*(I0,1X))')its(:,ip)
  END DO
  CLOSE(io_unit)
END IF
CALL oft_finalize
CONTAINS
! Fill the composite matrix with C (x) L from the local values of the scalar Laplacian
! (native matrices hold the local element contributions, summed over processes when
! applied, so the local values can be copied block by block)
SUBROUTINE build_operator()
INTEGER(i4) :: i,f,fg,n
INTEGER(i4) :: row(1)
CALL ML_vlag%current_level%mat_create(aop)
CALL aop%zero
SELECT TYPE(lop)
CLASS IS(oft_native_matrix)
  DO f=1,nf
    DO fg=1,nf
      DO i=1,lop%nr
        n=lop%kr(i+1)-lop%kr(i)
        IF(n==0)CYCLE
        row(1)=i
        CALL aop%add_values(row,lop%lc(lop%kr(i):lop%kr(i+1)-1), &
          RESHAPE(cmat(f,fg)*lop%M(lop%kr(i):lop%kr(i+1)-1),[1,n]),1,n,f,fg)
      END DO
    END DO
  END DO
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','build_operator',__FILE__)
END SELECT
CALL aop%assemble(v)
END SUBROUTINE build_operator
! Right-hand side: field f is f times the scalar right-hand side
SUBROUTINE set_rhs()
INTEGER(i4) :: f
REAL(r8), POINTER, DIMENSION(:) :: vals
NULLIFY(vals)
CALL gs%get_local(vals)
DO f=1,nf
  CALL g0%restore_local(REAL(f,r8)*vals,iblock=f)
END DO
DEALLOCATE(vals)
END SUBROUTINE set_rhs
! All cases for one number of blocks per process. Results (`res`):
! 1: layers covering the mesh, 2: restricted with one group per field and covering
! layers, 3: same with groups 1,1,2, 4: weighted additive with one group per field and
! covering layers, 5-7: no field groups, layers 0-2, 8-10: one group per field, layers
! 0-2, 11-13: groups 1,1,2, layers 0-2, 14: groups 1,1,1 of 2, one layer, 15: before
! update, 16: after update, 17: after delete and re-setup, 18: after update with new
! pattern (one group per field, one layer), 19: weighted additive and 20: additive with
! one group per field and one layer, 21: number of blocks with one group per field
SUBROUTINE run_cases(nparts,res)
INTEGER(i4), INTENT(in) :: nparts
INTEGER(i4), INTENT(out) :: res(nres)
INTEGER(i4) :: kfull,k,kd,nb
REAL(r8) :: rnorm
TYPE(oft_overlap_dofs) :: dofs
CHARACTER(LEN=16) :: tag
WRITE(tag,'(A,I1)')' nparts=',nparts
CALL fem_cell_overlap_create(cbld,mg_mesh,ML_vlag%current_level,nparts)
aop%cell_ovl=>cbld
cbld%field_groups=[1,1,2] ! Default field groups (used by nlocal<-1 without groups)
nb=cbld%nblocks()
!---Smallest number of point layers for every block to cover the whole mesh
kfull=-1
DO k=0,64
  CALL cbld%build(aop,v,k,OFT_CELL_ADJ_POINT,-1,dofs)
  IF(-oft_mpi_max(-MINVAL(dofs%blk(1:nb)%n))==aop%nrg)THEN
    kfull=k
    EXIT
  END IF
END DO
CALL overlap_dofs_delete(dofs)
IF(kfull<0)CALL oft_abort('Blocks do not cover the mesh','run_cases',__FILE__)
res(1)=kfull
!---Reference solution (graph overlap, converged tightly)
k=solve(1,1,.FALSE.,.FALSE.,x,rnorm,rtol=1.d-12)
!---Covering blocks: field (group) block-Jacobi, at most 3 iterations
res(2)=solve(-1,kfull,.TRUE.,.FALSE.,v,rnorm)
CALL check_tol('Covering blocks, one group per field: <=3 iterations'//tag,REAL(res(2),r8),3.d0)
CALL check_tol('Covering blocks, one group per field: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
res(3)=solve(-2,kfull,.TRUE.,.FALSE.,v,rnorm,groups=[1,1,2])
CALL check_tol('Covering blocks, groups 1,1,2: <=3 iterations'//tag,REAL(res(3),r8),3.d0)
CALL check_tol('Covering blocks, groups 1,1,2: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
res(4)=solve(-1,kfull,.TRUE.,.TRUE.,v,rnorm,weighted=.TRUE.)
CALL check_tol('Covering blocks, weighted: <=3 iterations'//tag,REAL(res(4),r8),3.d0)
CALL check_tol('Covering blocks, weighted: solution'//tag,maxdiff(v,x)/maxabs(x),1.d-8)
!---Iteration counts: no groups, one group per field, groups 1,1,2 x layers 0-2
DO k=0,2
  res(5+k)=solve(1,k,.TRUE.,.FALSE.,v,rnorm)
  CALL check_conv('No field groups'//tag,rnorm,v)
  res(8+k)=solve(-1,k,.TRUE.,.FALSE.,v,rnorm,nblocks=res(21))
  CALL check_conv('One group per field'//tag,rnorm,v)
  res(11+k)=solve(-2,k,.TRUE.,.FALSE.,v,rnorm,groups=[1,1,2])
  CALL check_conv('Groups 1,1,2'//tag,rnorm,v)
END DO
CALL check_exact('Blocks: builder blocks x fields'//tag,REAL(res(21),r8),REAL(nf*nb,r8))
!---Graph mode field based blocks (groups read per field) converge
k=solve(-1,1,.FALSE.,.FALSE.,v,rnorm)
CALL check_conv('Graph mode, one group per field'//tag,rnorm,v)
k=solve(-2,1,.FALSE.,.FALSE.,v,rnorm,groups=[1,1,2])
CALL check_conv('Graph mode, groups 1,1,2'//tag,rnorm,v)
!---Groups 1,1,1 of 2 (empty group 2): same blocks as without field groups
k=solve(1,1,.TRUE.,.FALSE.,w,rnorm)
res(14)=solve(-2,1,.TRUE.,.FALSE.,v,rnorm,groups=[1,1,1])
CALL check_exact('Empty group: iterations'//tag,REAL(res(14),r8),REAL(k,r8))
CALL check_tol('Empty group: solution'//tag,maxdiff(v,w)/maxabs(w),1.d-12)
!---No groups given: the builder's default groups (1,1,2)
k=solve(-2,1,.TRUE.,.FALSE.,w,rnorm,groups=[1,1,2])
kd=solve(-2,1,.TRUE.,.FALSE.,v,rnorm)
CALL check_exact('Default groups: iterations'//tag,REAL(kd,r8),REAL(k,r8))
CALL check_tol('Default groups: solution'//tag,maxdiff(v,w)/maxabs(w),1.d-12)
!---Additive combinations with one group per field
res(19)=solve(-1,1,.TRUE.,.TRUE.,v,rnorm,weighted=.TRUE.)
CALL check_conv('One group per field, weighted'//tag,rnorm,v)
res(20)=solve(-1,1,.TRUE.,.TRUE.,v,rnorm)
CALL check_conv('One group per field, additive'//tag,rnorm,v)
!---Numeric update and delete/re-setup
CALL check_update(nparts,res(15:18))
NULLIFY(aop%cell_ovl)
CALL cbld%delete
END SUBROUTINE run_cases
! Converged to the reference solution
SUBROUTINE check_conv(name,rnorm,y)
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), INTENT(in) :: rnorm
CLASS(oft_vector), INTENT(inout) :: y
CALL check_tol(name//': converged',rnorm,1.d-9)
CALL check_tol(name//': solution',maxdiff(y,x)/maxabs(x),1.d-5)
END SUBROUTINE check_conv
! Solve the coupled problem with GMRES preconditioned by Schwarz (cell mode with point
! adjacency, or graph mode) with direct sub-solves, returning the iteration count and
! relative residual (`groups`: field groups for nlocal<-1; `nblocks`: number of blocks
! after setup)
FUNCTION solve(nlocal,overlap,cell,additive,y,rnorm,rtol,weighted,groups,nblocks) RESULT(nits)
INTEGER(i4), INTENT(in) :: nlocal,overlap
LOGICAL, INTENT(in) :: cell,additive
CLASS(oft_vector), INTENT(inout) :: y
REAL(r8), INTENT(out) :: rnorm
REAL(r8), OPTIONAL, INTENT(in) :: rtol
LOGICAL, OPTIONAL, INTENT(in) :: weighted
INTEGER(i4), OPTIONAL, INTENT(in) :: groups(:)
INTEGER(i4), OPTIONAL, INTENT(out) :: nblocks
INTEGER(i4) :: nits
TYPE(oft_native_gmres_solver) :: linv
TYPE(oft_schwarz_precond), TARGET :: sw
TYPE(oft_lusolver), TARGET :: sub
sub%package=SUB_LU_PACK
CALL gmres_init(linv,sw,sub,nlocal,overlap,cell,additive,groups)
IF(PRESENT(rtol))linv%rtol=rtol
IF(PRESENT(weighted))sw%weighted=weighted
CALL y%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(y,g)
nits=linv%cits
rnorm=SQRT(g%dot(g)/g0%dot(g0))
IF(PRESENT(nblocks))nblocks=sw%nlocal
CALL sw%delete
CALL linv%delete
END FUNCTION solve
! Set up GMRES with a Schwarz preconditioner
SUBROUTINE gmres_init(linv,sw,sub,nlocal,overlap,cell,additive,groups)
TYPE(oft_native_gmres_solver), INTENT(inout) :: linv
TYPE(oft_schwarz_precond), TARGET, INTENT(inout) :: sw
CLASS(oft_solver), TARGET, INTENT(inout) :: sub
INTEGER(i4), INTENT(in) :: nlocal,overlap
LOGICAL, INTENT(in) :: cell,additive
INTEGER(i4), OPTIONAL, INTENT(in) :: groups(:)
linv%A=>aop
linv%its=-3
linv%nrits=20
linv%rtol=1.d-10
sw%nlocal=nlocal
IF(PRESENT(groups))sw%slice_group=groups
sw%overlap=overlap
sw%cell_mode=cell
sw%adjacency=OFT_CELL_ADJ_POINT
sw%additive=additive
sw%pre=>sub
linv%pre=>sw
END SUBROUTINE gmres_init
! Scaling the operator and updating the preconditioner matches a fresh setup; delete
! and re-setup, and an update with a new pattern, reproduce the first solve (one group
! per field, one layer)
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
CALL gmres_init(linv,sw,sub,-1,1,.TRUE.,.FALSE.)
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
nfresh=solve(-1,1,.TRUE.,.FALSE.,x,rnorm)
CALL check_exact('Updated = fresh setup (iterations)'//tag,REAL(res(2),r8),REAL(nfresh,r8))
CALL check_tol('Updated = fresh setup (solution)'//tag,maxdiff(w,x)/maxabs(x),1.d-12)
CALL w%scale(2.d0)
CALL check_tol('Updated solution scales with operator'//tag,maxdiff(v,w)/maxabs(v),1.d-8)
CALL scale_op(0.5d0)
!---Delete and re-setup on the original operator (field groups are kept)
CALL sw%delete
CALL check_exact('Delete: requested nlocal kept'//tag,REAL(sw%nlocal,r8),-1.d0)
CALL linv%update(.FALSE.)
CALL w%set(0.d0)
CALL g%copy(g0)
CALL linv%apply(w,g)
res(3)=linv%cits
CALL check_exact('Re-setup iterations'//tag,REAL(res(3),r8),REAL(res(1),r8))
CALL check_tol('Re-setup solution'//tag,maxdiff(v,w)/maxabs(v),1.d-12)
!---Update with a new pattern (rebuild with the requested field groups)
CALL linv%update(.TRUE.)
CALL check_exact('New pattern: requested nlocal restored'//tag,REAL(sw%nlocal,r8),-1.d0)
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
SELECT TYPE(aop)
CLASS IS(oft_native_matrix)
  aop%M=fac*aop%M
  aop%full_current=.FALSE.
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
END PROGRAM test_cell_schwarz_fields
