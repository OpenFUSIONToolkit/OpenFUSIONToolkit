!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_native_la.F90
!
!> Regression tests for native vector and matrix operations that are not
!! otherwise exercised by the solver tests. Tests are performed on a unit cube
!! using a second order Lagrange Laplacian made non-symmetric by row scaling.
!!
!! The current test cases are:
!! - Transpose products: adjoint identity \f$ x \cdot (A y) = (A^T x) \cdot y \f$
!! - 1-norm and inf-norm against reductions over the owned (slice) values
!! - Complex vector products and mixed real/complex dot products
!! - Borrowed (zero-copy) access with each access mode, with and without
!!   synchronization, for full vectors and blocks, against get_local/restore_local
!! - copy, maxpy and mdot for real and complex vectors against add and dot
!! - Boundary condition operators that skip synchronization against the copying path
!! - Summation stitching gives bitwise identical shared values (eg. periodic images)
!! - Complex stitching (all methods) against stitching real and imaginary parts
!! - Accumulate access for blocks and full vectors, real and complex, against
!!   get_local/restore_local(add=.TRUE.)
!! - Shared matrix non-zero patterns: sharing by FE operators, identical results for
!!   independent and pattern-sharing matrices (real and complex), linkage built once,
!!   and patterns freed with their last reference
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_native_la
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix, oft_vector_ptr, oft_cvector_ptr, &
  OFT_ACCESS_READ, OFT_ACCESS_WRITE, OFT_ACCESS_READWRITE, OFT_ACCESS_ACCUMULATE
USE oft_native_la, ONLY: oft_native_matrix, oft_native_vector, oft_native_cvector, &
  oft_native_cmatrix, oft_native_pattern, native_matrix_share, native_matrix_setup_full, &
  native_pattern_live
USE oft_la_base, ONLY: oft_cmatrix, oft_graph_ptr
USE oft_la_utils, ONLY: create_matrix
USE fem_base, ONLY: oft_ml_fem_type, oft_afem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_getlop, oft_lag_getmop, oft_lag_zerob, oft_vlag_zerob
IMPLICIT NONE
INTEGER(i4) :: io_unit,nfail
REAL(r8), PARAMETER :: tol=1.d-12
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange,minlev=-1)
!---Run tests
nfail=0
CALL test_ops()
CALL test_access()
CALL test_multi()
CALL test_bc()
CALL test_stitch()
CALL test_accumulate()
CALL test_pattern()
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='native_la.results')
  WRITE(io_unit,*)nfail
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
! Check relative agreement of two values and record failures
SUBROUTINE check(name,val,ref)
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), INTENT(in) :: val,ref
IF(ABS(val-ref)>tol*MAX(ABS(ref),1.d0))THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,'(3A,2ES25.16)')'FAILED: ',name,' ',val,ref
ELSE
  IF(oft_env%head_proc)WRITE(*,'(3A,ES25.16)')'PASSED: ',name,' ',val
END IF
END SUBROUTINE check
! Check that a value is exactly equal to a reference (bitwise reproducible operations)
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
! Global maximum of |x_i - y_i| over all local entries (including shared entries)
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
! Global maximum of |x_i - y_i| over all local entries (complex)
FUNCTION cmaxdiff(x,y) RESULT(d)
CLASS(oft_cvector), INTENT(inout) :: x,y
REAL(r8) :: d
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: xv,yv
NULLIFY(xv,yv)
CALL x%borrow_local(xv,OFT_ACCESS_READ)
CALL y%borrow_local(yv,OFT_ACCESS_READ)
d=0.d0
IF(x%n>0)d=MAXVAL(ABS(xv-yv))
CALL x%release_local(xv,OFT_ACCESS_READ,.FALSE.)
CALL y%release_local(yv,OFT_ACCESS_READ,.FALSE.)
d=oft_mpi_max(d)
END FUNCTION cmaxdiff
! Global maximum difference between shared values and their owner's value
FUNCTION inconsistency(x) RESULT(d)
CLASS(oft_vector), INTENT(inout) :: x
REAL(r8) :: d
CLASS(oft_vector), POINTER :: y
CALL x%new(y)
CALL y%copy(x)
CALL y%stitch(0)
d=maxdiff(x,y)
CALL y%delete
DEALLOCATE(y)
END FUNCTION inconsistency
! Reproducible value for a global index
PURE FUNCTION gval(ig,a,b) RESULT(v)
INTEGER(i8), INTENT(in) :: ig
REAL(r8), INTENT(in) :: a,b
REAL(r8) :: v
v=SIN(a*REAL(ig,r8))+b
END FUNCTION gval
! Reproducible complex values from the global index of each local entry, set without
! stitching (complex vector seams are not supported in parallel)
SUBROUTINE cfill(x,a,b)
CLASS(oft_cvector), INTENT(inout) :: x
REAL(r8), INTENT(in) :: a,b
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: vals
INTEGER(i4) :: i
INTEGER(i8) :: ig
NULLIFY(vals)
CALL x%borrow_local(vals,OFT_ACCESS_WRITE)
DO i=1,x%n
  ig=x%map(1)%offsetg+x%map(1)%lge(i)
  vals(i)=CMPLX(gval(ig,a,b),gval(ig,b+0.3d0,-a),c8)
END DO
CALL x%release_local(vals,OFT_ACCESS_WRITE,.FALSE.)
END SUBROUTINE cfill
! Reproducible values from the global index of each local entry
SUBROUTINE fill(x,a,b)
CLASS(oft_vector), INTENT(inout) :: x
REAL(r8), INTENT(in) :: a,b
REAL(r8), POINTER, DIMENSION(:) :: vals
INTEGER(i4) :: i
NULLIFY(vals)
CALL x%get_local(vals)
DO i=1,x%n
  vals(i)=SIN(a*REAL(x%map(1)%offsetg+x%map(1)%lge(i),r8))+b
END DO
CALL x%restore_local(vals)
DEALLOCATE(vals)
END SUBROUTINE fill
! Exercise native matrix and vector operations
SUBROUTINE test_ops()
INTEGER(i4) :: i,j
REAL(r8) :: v1,v2
REAL(r8), POINTER, DIMENSION(:) :: vals,xr,xi
COMPLEX(c8) :: c1,c2
COMPLEX(c8), POINTER, DIMENSION(:) :: cvals
CLASS(oft_vector), POINTER :: x,y,z,w
CLASS(oft_cvector), POINTER :: cx,cz
CLASS(oft_matrix), POINTER :: lop => NULL()
!---Create vectors and operator
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL ML_oft_lagrange%vec_create(z)
CALL ML_oft_lagrange%vec_create(w)
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'none')
!---Scale rows by a function of the global row index to make the operator non-symmetric
SELECT TYPE(lop)
CLASS IS(oft_native_matrix)
  DO i=1,lop%nr
    v1=1.d0+0.5d0*SIN(0.7d0*REAL(x%map(1)%offsetg+x%map(1)%lge(i),r8))
    DO j=lop%kr(i),lop%kr(i+1)-1
      lop%M(j)=lop%M(j)*v1
    END DO
  END DO
CLASS DEFAULT
  CALL oft_abort('Native matrix required','test_ops',__FILE__)
END SELECT
CALL fill(x,0.37d0,0.1d0)
CALL fill(y,0.11d0,-0.2d0)
!---Operator must be non-symmetric for the transpose test to be meaningful
CALL lop%apply(x,z)
CALL lop%applyt(x,w)
CALL w%add(1.d0,-1.d0,z)
v1=w%norm(2)
IF(v1<1.d-3*z%norm(2))THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,*)'FAILED: Operator is symmetric'
END IF
!---Adjoint identity
CALL lop%apply(y,z)
v1=x%dot(z)
CALL lop%applyt(x,w)
v2=w%dot(y)
CALL check('x.(A y) = (A^T x).y',v2,v1)
!---Repeated transpose product (cached index) gives the same result
CALL lop%applyt(x,z)
CALL z%add(1.d0,-1.d0,w)
CALL check('Repeated A^T x',z%norm(3),0.d0)
!---Norms against reductions over owned values
NULLIFY(vals)
CALL w%get_slice(vals)
v1=0.d0
v2=0.d0
DO i=1,w%nslice
  v1=v1+ABS(vals(i))
  v2=MAX(v2,ABS(vals(i)))
END DO
DEALLOCATE(vals)
CALL check('1-norm',w%norm(1),oft_mpi_sum(v1))
CALL check('Inf-norm',w%norm(3),oft_mpi_max(v2))
!---Complex vectors
CALL x%new(cx)
CALL x%new(cz)
NULLIFY(xr,xi,cvals)
CALL x%get_local(xr)
CALL y%get_local(xi)
ALLOCATE(cvals(x%n))
cvals=CMPLX(xr,xi,c8)
CALL cx%restore_local(cvals)
!---Transpose product matches real and imaginary parts applied separately
CALL lop%applyt(cx,cz)
CALL cz%get_local(cvals)
CALL lop%applyt(x,z)
CALL lop%applyt(y,w)
CALL z%get_local(xr)
CALL w%get_local(xi)
CALL check('Re(A^T cx)',MAXVAL(ABS(REAL(cvals,r8)-xr)),0.d0)
CALL check('Im(A^T cx)',MAXVAL(ABS(AIMAG(cvals)-xi)),0.d0)
!---Mixed real/complex dot products
c1=x%dot(cz)
c2=cz%dot(x)
CALL check('Re(x.cz)',REAL(c1,r8),x%dot(z))
CALL check('Im(x.cz)',AIMAG(c1),x%dot(w))
CALL check('Re(cz.x)',REAL(c2,r8),x%dot(z))
CALL check('Im(cz.x)',AIMAG(c2),-x%dot(w))
!---Complex norms (owned values only)
DEALLOCATE(cvals)
CALL cz%get_slice(cvals)
v1=0.d0
v2=0.d0
DO i=1,cz%nslice
  v1=v1+ABS(cvals(i))
  v2=MAX(v2,ABS(cvals(i)))
END DO
CALL check('Complex 1-norm',cz%norm(1),oft_mpi_sum(v1))
CALL check('Complex inf-norm',cz%norm(3),oft_mpi_max(v2))
DEALLOCATE(xr,xi,cvals)
CALL cx%delete
CALL cz%delete
DEALLOCATE(cx,cz)
!---Cleanup
CALL x%delete
CALL y%delete
CALL z%delete
CALL w%delete
DEALLOCATE(x,y,z,w)
CALL lop%delete
DEALLOCATE(lop)
END SUBROUTINE test_ops
! Borrowed access to local values with each access mode
SUBROUTINE test_access()
INTEGER(i4) :: i,iblk
REAL(r8) :: v1
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: bvals
CLASS(oft_vector), POINTER :: x,y,z
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL ML_oft_lagrange%vec_create(z)
CALL fill(x,0.37d0,0.1d0)
!---Read access sees the same values as get_local
NULLIFY(vals,bvals)
CALL x%get_local(vals)
CALL x%borrow_local(bvals,OFT_ACCESS_READ)
v1=0.d0
IF(x%n>0)v1=MAXVAL(ABS(bvals-vals))
CALL x%release_local(bvals,OFT_ACCESS_READ,.FALSE.)
CALL check_exact('Borrow read = get_local',oft_mpi_max(v1),0.d0)
IF(ASSOCIATED(bvals))THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,*)'FAILED: Borrowed array not nullified on release'
END IF
!---Write access (every entry assigned) matches restore_local of the same values
CALL y%borrow_local(bvals,OFT_ACCESS_WRITE)
DO i=1,y%n
  bvals(i)=vals(i)
END DO
CALL y%release_local(bvals,OFT_ACCESS_WRITE,.FALSE.)
CALL check_exact('Borrow write = restore_local',maxdiff(x,y),0.d0)
!---Read-write access modifies the vector in place
CALL y%borrow_local(bvals,OFT_ACCESS_READWRITE)
DO i=1,y%n
  bvals(i)=2.d0*bvals(i)
END DO
CALL y%release_local(bvals,OFT_ACCESS_READWRITE,.FALSE.)
CALL z%copy(x)
CALL z%scale(2.d0)
CALL check_exact('Borrow read-write in place',maxdiff(y,z),0.d0)
!---Process dependent values: sync=.TRUE. takes the owner's values, as restore_local
DO i=1,x%n
  vals(i)=vals(i)+REAL(oft_env%rank,r8)
END DO
CALL z%restore_local(vals)
CALL y%borrow_local(bvals,OFT_ACCESS_WRITE)
DO i=1,y%n
  bvals(i)=vals(i)
END DO
CALL y%release_local(bvals,OFT_ACCESS_WRITE,.TRUE.)
CALL check_exact('Borrow write with sync = restore_local',maxdiff(y,z),0.d0)
CALL check_exact('Shared values consistent after sync',inconsistency(y),0.d0)
!---Without sync a modification made identically on all processes keeps shared values consistent
CALL z%copy(y)
CALL y%borrow_local(bvals,OFT_ACCESS_READWRITE)
DO i=1,y%n
  bvals(i)=bvals(i)-0.5d0*REAL(MOD(x%map(1)%offsetg+x%map(1)%lge(i),7_i8),r8)
END DO
CALL y%release_local(bvals,OFT_ACCESS_READWRITE,.FALSE.)
CALL z%borrow_local(bvals,OFT_ACCESS_READWRITE)
DO i=1,z%n
  bvals(i)=bvals(i)-0.5d0*REAL(MOD(x%map(1)%offsetg+x%map(1)%lge(i),7_i8),r8)
END DO
CALL z%release_local(bvals,OFT_ACCESS_READWRITE,.TRUE.)
CALL check_exact('Borrow read-write without sync = with sync',maxdiff(y,z),0.d0)
CALL check_exact('Shared values consistent without sync',inconsistency(y),0.d0)
DEALLOCATE(vals)
CALL x%delete
CALL y%delete
CALL z%delete
DEALLOCATE(x,y,z)
!---Block access on a 3-block (vector Lagrange) field
CALL ML_oft_vlagrange%vec_create(x)
CALL ML_oft_vlagrange%vec_create(y)
CALL x%set(0.d0,random=.TRUE.)
CALL y%copy(x)
DO iblk=1,3
  NULLIFY(vals)
  CALL x%get_local(vals,iblk)
  CALL x%borrow_local(bvals,OFT_ACCESS_READ,iblk)
  v1=0.d0
  IF(SIZE(vals)>0)v1=MAXVAL(ABS(bvals-vals))
  IF(SIZE(bvals)/=x%map(iblk)%n)v1=1.d0
  CALL x%release_local(bvals,OFT_ACCESS_READ,.FALSE.,iblk)
  CALL check_exact('Block borrow read = get_local',oft_mpi_max(v1),0.d0)
  !---Scale one block in place, compare with block restore_local
  CALL y%borrow_local(bvals,OFT_ACCESS_READWRITE,iblk)
  DO i=1,SIZE(bvals)
    bvals(i)=REAL(iblk+1,r8)*bvals(i)
  END DO
  CALL y%release_local(bvals,OFT_ACCESS_READWRITE,.FALSE.,iblk)
  vals=REAL(iblk+1,r8)*vals
  CALL x%restore_local(vals,iblk)
  DEALLOCATE(vals)
END DO
CALL check_exact('Block borrow read-write = block restore_local',maxdiff(x,y),0.d0)
CALL x%delete
CALL y%delete
DEALLOCATE(x,y)
END SUBROUTINE test_access
! copy, maxpy and mdot against add and dot (real and complex)
SUBROUTINE test_multi()
INTEGER(i4), PARAMETER :: nv=4
INTEGER(i4) :: i,k
REAL(r8) :: alpha(nv),dots(nv),v1
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: bvals
COMPLEX(c8) :: calpha(nv),cdots(nv)
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: cvals
TYPE(oft_vector_ptr) :: vs(nv),ws(nv)
TYPE(oft_cvector_ptr) :: cvs(nv)
CLASS(oft_vector), POINTER :: x,y
CLASS(oft_cvector), POINTER :: cx,cy
DO k=1,nv
  CALL ML_oft_lagrange%vec_create(vs(k)%f)
  CALL fill(vs(k)%f,0.13d0*REAL(k,r8),0.1d0*REAL(k-2,r8))
  alpha(k)=0.3d0*REAL(k,r8)-0.7d0
END DO
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL fill(x,0.29d0,-0.05d0)
!---copy does not read the destination
NULLIFY(bvals)
CALL y%borrow_local(bvals,OFT_ACCESS_WRITE)
bvals=ieee_value(1.d0,ieee_quiet_nan)
CALL y%release_local(bvals,OFT_ACCESS_WRITE,.FALSE.)
CALL y%copy(x)
CALL check_exact('Copy',maxdiff(x,y),0.d0)
CALL y%borrow_local(bvals,OFT_ACCESS_WRITE)
bvals=ieee_value(1.d0,ieee_quiet_nan)
CALL y%release_local(bvals,OFT_ACCESS_WRITE,.FALSE.)
CALL y%add(0.d0,1.d0,x)
CALL check_exact('Add with gamma=0',maxdiff(x,y),0.d0)
!---maxpy is bitwise identical to successive add calls
CALL y%copy(x)
CALL y%maxpy(alpha,vs,nv)
DO k=1,nv
  CALL x%add(1.d0,alpha(k),vs(k)%f)
END DO
CALL check_exact('maxpy = successive add',maxdiff(x,y),0.d0)
!---mdot against individual dot products (self included in the list)
CALL fill(x,0.29d0,-0.05d0)
ws(1:nv-1)=vs(1:nv-1)
ws(nv)%f=>x
dots=x%mdot(ws,nv)
DO k=1,nv
  CALL check('mdot',dots(k),x%dot(ws(k)%f))
END DO
!---Complex vectors, filled without stitching
CALL x%new(cx)
CALL x%new(cy)
DO k=1,nv
  CALL x%new(cvs(k)%f)
  CALL cfill(cvs(k)%f,0.17d0*REAL(k,r8),0.2d0*REAL(k-2,r8))
  calpha(k)=CMPLX(alpha(k),0.1d0*REAL(k,r8),c8)
END DO
CALL cfill(cx,0.31d0,0.07d0)
NULLIFY(cvals)
CALL cy%borrow_local(cvals,OFT_ACCESS_WRITE)
cvals=CMPLX(ieee_value(1.d0,ieee_quiet_nan),0.d0,c8)
CALL cy%release_local(cvals,OFT_ACCESS_WRITE,.FALSE.)
CALL cy%copy(cx)
CALL check_exact('Complex copy',cmaxdiff(cx,cy),0.d0)
CALL cy%maxpy(calpha,cvs,nv)
DO k=1,nv
  CALL cx%add((1.d0,0.d0),calpha(k),cvs(k)%f)
END DO
CALL check_exact('Complex maxpy = successive add',cmaxdiff(cx,cy),0.d0)
!---Complex mdot (conjugates self) against individual dot products
cdots=cx%mdot(cvs,nv)
DO k=1,nv
  CALL check('Re(complex mdot)',REAL(cdots(k),r8),REAL(cx%dot(cvs(k)%f),r8))
  CALL check('Im(complex mdot)',AIMAG(cdots(k)),AIMAG(cx%dot(cvs(k)%f)))
END DO
!---Mixed real/complex mdot
cdots=cx%mdot(vs,nv)
DO k=1,nv
  CALL check('Re(complex.real mdot)',REAL(cdots(k),r8),REAL(cx%dot(vs(k)%f),r8))
  CALL check('Im(complex.real mdot)',AIMAG(cdots(k)),AIMAG(cx%dot(vs(k)%f)))
END DO
cdots=x%mdot(cvs,nv)
DO k=1,nv
  CALL check('Re(real.complex mdot)',REAL(cdots(k),r8),REAL(x%dot(cvs(k)%f),r8))
  CALL check('Im(real.complex mdot)',AIMAG(cdots(k)),AIMAG(x%dot(cvs(k)%f)))
END DO
!---Complex copy from a real vector
CALL cy%copy(x)
CALL cy%borrow_local(cvals,OFT_ACCESS_READ)
CALL x%borrow_local(bvals,OFT_ACCESS_READ)
v1=0.d0
DO i=1,x%n
  v1=MAX(v1,ABS(cvals(i)-bvals(i)))
END DO
CALL x%release_local(bvals,OFT_ACCESS_READ,.FALSE.)
CALL cy%release_local(cvals,OFT_ACCESS_READ,.FALSE.)
CALL check_exact('Complex copy from real',oft_mpi_max(v1),0.d0)
!---Cleanup
DO k=1,nv
  CALL cvs(k)%f%delete
  DEALLOCATE(cvs(k)%f)
END DO
DO k=1,nv
  CALL vs(k)%f%delete
  DEALLOCATE(vs(k)%f)
END DO
CALL cx%delete
CALL cy%delete
CALL x%delete
CALL y%delete
DEALLOCATE(cx,cy,x,y)
END SUBROUTINE test_multi
! Boundary condition operators that release without synchronization give the same result
! as zeroing a copy and restoring it (which synchronizes)
SUBROUTINE test_bc()
INTEGER(i4) :: i,j,iblk
REAL(r8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: x,y
CLASS(oft_afem_type), POINTER :: lag_rep
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(oft_vlag_zerob), TARGET :: vlag_zerob
lag_zerob%ML_lag_rep=>ML_oft_lagrange
vlag_zerob%ML_vlag_rep=>ML_oft_vlagrange
!---Scalar field
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL fill(x,0.37d0,0.1d0)
CALL y%copy(x)
CALL lag_zerob%apply(x)
NULLIFY(vals)
CALL y%get_local(vals)
DO i=1,ML_oft_lagrange%current_level%nbe
  j=ML_oft_lagrange%current_level%lbe(i)
  IF(ML_oft_lagrange%current_level%global%gbe(j))vals(j)=0.d0
END DO
CALL y%restore_local(vals)
DEALLOCATE(vals)
CALL check_exact('Lagrange zerob',maxdiff(x,y),0.d0)
CALL check_exact('Lagrange zerob consistent',inconsistency(x),0.d0)
CALL x%delete
CALL y%delete
DEALLOCATE(x,y)
!---Vector field
CALL ML_oft_vlagrange%vec_create(x)
CALL ML_oft_vlagrange%vec_create(y)
CALL x%set(1.d0,random=.TRUE.)
CALL y%copy(x)
CALL vlag_zerob%apply(x)
lag_rep=>ML_oft_vlagrange%current_level%fields(1)%fe
DO iblk=1,3
  NULLIFY(vals)
  CALL y%get_local(vals,iblk)
  DO i=1,lag_rep%nbe
    j=lag_rep%lbe(i)
    IF(lag_rep%global%gbe(j))vals(j)=0.d0
  END DO
  CALL y%restore_local(vals,iblk)
  DEALLOCATE(vals)
END DO
CALL check_exact('Vector Lagrange zerob',maxdiff(x,y),0.d0)
CALL check_exact('Vector Lagrange zerob consistent',inconsistency(x),0.d0)
CALL x%delete
CALL y%delete
DEALLOCATE(x,y)
END SUBROUTINE test_bc
! Summation stitching of values that differ between copies of a shared entry gives every copy
! the same value, bit for bit, and conserves the global sum. Run on periodic meshes this covers
! several images of an entry on one process.
SUBROUTINE test_stitch()
INTEGER(i4) :: i,m
REAL(r8) :: s_raw,s_sum,d
CLASS(oft_vector), POINTER :: x,y,ones
CLASS(oft_cvector), POINTER :: cx,cy
!---Real
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(ones)
CALL ones%set(1.d0)
s_raw=0.d0
SELECT TYPE(x)
TYPE IS(oft_native_vector)
  DO i=1,x%n
    x%v(i)=gval(INT(i+1000*oft_env%rank,i8),0.71d0,0.1d0)/3.d0
    s_raw=s_raw+x%v(i)
  END DO
CLASS DEFAULT
  CALL oft_abort('Native vector expected','test_stitch',__FILE__)
END SELECT
s_raw=oft_mpi_sum(s_raw)
CALL x%stitch(1)
CALL check_exact('Summation stitch consistent',inconsistency(x),0.d0)
s_sum=x%dot(ones)
CALL check('Summation stitch conserves sum',s_sum,s_raw)
!---Complex: summation is consistent across copies and each stitching method matches
!---stitching the real and imaginary parts separately
CALL x%new(cx)
CALL x%new(cy)
CALL x%new(y)
DO m=0,2
  CALL cfill_raw(cx)
  SELECT TYPE(cx)
  TYPE IS(oft_native_cvector)
    SELECT TYPE(x)
    TYPE IS(oft_native_vector)
      SELECT TYPE(y)
      TYPE IS(oft_native_vector)
        x%v=REAL(cx%v,r8)
        y%v=AIMAG(cx%v)
        CALL cx%stitch(m)
        CALL x%stitch(m)
        CALL y%stitch(m)
        IF(m==1)THEN
          ! Elements with three or more contributions are summed in a different order
          d=MAXVAL(ABS(cx%v-CMPLX(x%v,y%v,c8)))
          CALL check('Complex stitch(1) matches real/imaginary parts',oft_mpi_max(d),0.d0)
          CALL cy%add((0.d0,0.d0),(1.d0,0.d0),cx)
          CALL cy%stitch(0)
          CALL check_exact('Complex summation stitch consistent',cmaxdiff(cx,cy),0.d0)
        ELSE
          d=MAXVAL(ABS(cx%v-CMPLX(x%v,y%v,c8)))
          CALL check_exact('Complex stitch('//CHAR(48+m)//') matches real/imaginary parts',oft_mpi_max(d),0.d0)
        END IF
      END SELECT
    END SELECT
  END SELECT
END DO
!---Complex summation of a uniform field counts every copy once
CALL cx%set((1.d0,-2.d0))
CALL cx%stitch(1)
CALL ones%set(1.d0)
CALL ones%stitch(1)
SELECT TYPE(cx)
TYPE IS(oft_native_cvector)
  SELECT TYPE(ones)
  TYPE IS(oft_native_vector)
    CALL check_exact('Complex summation multiplicity',oft_mpi_max(MAXVAL(ABS(cx%v-CMPLX(1.d0,-2.d0,c8)*ones%v))),0.d0)
  END SELECT
END SELECT
CALL cx%delete
CALL cy%delete
CALL y%delete
DEALLOCATE(cx,cy,y)
CALL x%delete
CALL ones%delete
DEALLOCATE(x,ones)
END SUBROUTINE test_stitch
! Distinct values on every copy of an entry (local index and rank dependent)
SUBROUTINE cfill_raw(c)
CLASS(oft_cvector), INTENT(inout) :: c
INTEGER(i4) :: i
SELECT TYPE(c)
TYPE IS(oft_native_cvector)
  DO i=1,c%n
    c%v(i)=CMPLX(gval(INT(i+1000*oft_env%rank,i8),0.71d0,0.1d0), &
      gval(INT(i+1000*oft_env%rank,i8),0.53d0,-0.2d0),c8)/3.d0
  END DO
CLASS DEFAULT
  CALL oft_abort('Native complex vector expected','cfill_raw',__FILE__)
END SELECT
END SUBROUTINE cfill_raw
! Shared non-zero patterns of native CRS matrices
! Accumulate access (zeroed storage, summed across copies and added on completion) against the
! get_local/restore_local(add=.TRUE.) pattern it replaces, for a blocked vector and a full vector
SUBROUTINE test_accumulate()
INTEGER(i4) :: i,k,iblk
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: avals
COMPLEX(c8), POINTER, DIMENSION(:) :: cvals
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: cavals
CLASS(oft_vector), POINTER :: x,y,z
CLASS(oft_cvector), POINTER :: cy,cz
!---Blocked (vector Lagrange) field, starting from non-zero values
CALL ML_oft_vlagrange%vec_create(x)
CALL x%new(y)
CALL x%new(z)
CALL x%set(1.d0,random=.TRUE.)
DO k=1,2 ! Second pass reuses the accumulation storage
  CALL y%copy(x)
  CALL z%copy(x)
  NULLIFY(vals,avals)
  DO iblk=1,3
    CALL z%get_local(vals,iblk)
    DO i=1,SIZE(vals)
      vals(i)=contrib(i,iblk+k)
    END DO
    CALL z%restore_local(vals,iblk,add=.TRUE.,wait=(iblk/=3))
    DEALLOCATE(vals)
  END DO
  DO iblk=3,1,-1 ! Any order; the last release completes the accumulation
    CALL y%borrow_local(avals,OFT_ACCESS_ACCUMULATE,iblk)
    IF(MAXVAL(ABS(avals))/=0.d0)THEN
      nfail=nfail+1
      WRITE(*,*)'FAILED: Accumulation storage not zeroed'
    END IF
    DO i=1,SIZE(avals)
      avals(i)=contrib(i,iblk+k)
    END DO
    CALL y%release_local(avals,OFT_ACCESS_ACCUMULATE,(iblk==1),iblk)
  END DO
  CALL check_exact('Blocked accumulate = restore_local(add)',maxdiff(y,z),0.d0)
END DO
CALL check_exact('Accumulated values consistent',inconsistency(y),0.d0)
!---Only borrowed blocks change
CALL y%copy(x)
CALL y%borrow_local(avals,OFT_ACCESS_ACCUMULATE,2)
avals=1.d0
CALL y%release_local(avals,OFT_ACCESS_ACCUMULATE,.TRUE.,2)
CALL z%copy(y)
CALL z%add(1.d0,-1.d0,x)
NULLIFY(vals)
CALL z%get_local(vals,1)
CALL check_exact('Accumulate leaves other blocks',oft_mpi_max(MAXVAL(ABS(vals))),0.d0)
DEALLOCATE(vals)
CALL x%delete
CALL y%delete
CALL z%delete
DEALLOCATE(x,y,z)
!---Full (scalar) vector, real and complex
CALL ML_oft_lagrange%vec_create(y)
CALL y%new(z)
CALL fill(y,0.21d0,0.3d0)
CALL z%copy(y)
NULLIFY(vals,avals)
ALLOCATE(vals(z%n))
DO i=1,z%n
  vals(i)=contrib(i,7)
END DO
CALL z%restore_local(vals,add=.TRUE.)
DEALLOCATE(vals)
CALL y%borrow_local(avals,OFT_ACCESS_ACCUMULATE)
DO i=1,y%n
  avals(i)=contrib(i,7)
END DO
CALL y%release_local(avals,OFT_ACCESS_ACCUMULATE,.TRUE.)
CALL check_exact('Full accumulate = restore_local(add)',maxdiff(y,z),0.d0)
CALL y%new(cy)
CALL y%new(cz)
CALL cfill(cy,0.17d0,0.2d0)
CALL cz%copy(cy)
NULLIFY(cvals,cavals)
ALLOCATE(cvals(cz%n))
DO i=1,cz%n
  cvals(i)=CMPLX(contrib(i,8),contrib(i,9),c8)
END DO
CALL cz%restore_local(cvals,add=.TRUE.)
DEALLOCATE(cvals)
CALL cy%borrow_local(cavals,OFT_ACCESS_ACCUMULATE)
DO i=1,cy%n
  cavals(i)=CMPLX(contrib(i,8),contrib(i,9),c8)
END DO
CALL cy%release_local(cavals,OFT_ACCESS_ACCUMULATE,.TRUE.)
CALL check_exact('Complex accumulate = restore_local(add)',cmaxdiff(cy,cz),0.d0)
CALL y%delete
CALL z%delete
CALL cy%delete
CALL cz%delete
DEALLOCATE(y,z,cy,cz)
END SUBROUTINE test_accumulate
! Contribution that differs between copies of an entry (local index and rank dependent)
PURE FUNCTION contrib(i,k) RESULT(v)
INTEGER(i4), INTENT(in) :: i,k
REAL(r8) :: v
v=COS(0.37d0*REAL(i+1000*oft_env%rank,r8)+REAL(k,r8))/7.d0
END FUNCTION contrib
SUBROUTINE test_pattern()
INTEGER(i4) :: i,nlive0,nlive1,refs0
REAL(r8) :: v1
REAL(r8), POINTER, DIMENSION(:) :: xr,xi
COMPLEX(c8), POINTER, DIMENSION(:) :: cvals
CLASS(oft_vector), POINTER :: x,y,z,w
CLASS(oft_cvector), POINTER :: cx,cz
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: indep => NULL()
CLASS(oft_matrix), POINTER :: shared => NULL()
CLASS(oft_matrix), POINTER :: shared2 => NULL()
CLASS(oft_cmatrix), POINTER :: cshared => NULL()
TYPE(oft_graph_ptr) :: graphs(1,1)
TYPE(oft_native_matrix), POINTER :: A,B,B2,L,Mm
TYPE(oft_native_cmatrix), POINTER :: C
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL ML_oft_lagrange%vec_create(z)
CALL ML_oft_lagrange%vec_create(w)
CALL fill(x,0.37d0,0.1d0)
CALL fill(y,0.11d0,-0.2d0)
nlive0=native_pattern_live
!---Operators created by the same FE representation share its cached pattern
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'none')
L=>native_cast(lop)
refs0=L%pattern%refcount
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
Mm=>native_cast(mop)
CALL check_exact('FE operators share pattern',MERGE(1.d0,0.d0,ASSOCIATED(L%pattern,Mm%pattern)),1.d0)
CALL check_exact('Shared pattern reference count',REAL(Mm%pattern%refcount-refs0,r8),1.d0)
CALL check_exact('No pattern created for second operator',REAL(native_pattern_live-nlive0,r8),0.d0)
!---Independent matrix with the same structure (no pattern cache) and the operator values
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=ML_oft_lagrange%current_level%ne
graphs(1,1)%g%nrg=ML_oft_lagrange%current_level%global%ne
graphs(1,1)%g%nc=ML_oft_lagrange%current_level%ne
graphs(1,1)%g%ncg=ML_oft_lagrange%current_level%global%ne
graphs(1,1)%g%nnz=ML_oft_lagrange%current_level%nee
graphs(1,1)%g%kr=>ML_oft_lagrange%current_level%kee
graphs(1,1)%g%lc=>ML_oft_lagrange%current_level%lee
CALL create_matrix(indep,graphs,x,x)
DEALLOCATE(graphs(1,1)%g)
A=>native_cast(indep)
CALL check_exact('Independent matrix has new pattern',REAL(native_pattern_live-nlive0,r8),1.d0)
CALL check_exact('Independent pattern not shared',MERGE(1.d0,0.d0,ASSOCIATED(A%pattern,L%pattern)),0.d0)
A%M=L%M
!---Matrix sharing the pattern of the independent matrix gives identical results
CALL native_matrix_share(indep,shared)
B=>native_cast(shared)
CALL check_exact('Share uses source pattern',MERGE(1.d0,0.d0,ASSOCIATED(A%pattern,B%pattern)),1.d0)
CALL check_exact('Share reference count',REAL(A%pattern%refcount,r8),2.d0)
CALL check_exact('Share zero-initialized',oft_mpi_max(MAXVAL(ABS(B%M))),0.d0)
B%M=A%M
CALL indep%apply(x,z)
CALL shared%apply(x,w)
CALL check_exact('Shared A x = independent A x',maxdiff(z,w),0.d0)
CALL indep%applyt(y,z)
CALL shared%applyt(y,w)
CALL check_exact('Shared A^T y = independent A^T y',maxdiff(z,w),0.d0)
CALL lop%apply(x,w)
CALL indep%apply(x,z)
CALL check_exact('FE operator A x = independent A x',maxdiff(z,w),0.d0)
!---Deleting one of several matrices sharing a pattern keeps it for the others
nlive1=native_pattern_live
CALL indep%delete
DEALLOCATE(indep)
CALL check_exact('Pattern kept after deleting one user',REAL(native_pattern_live-nlive1,r8),0.d0)
CALL check_exact('Reference count after delete',REAL(B%pattern%refcount,r8),1.d0)
CALL shared%apply(x,w)
CALL check_exact('A x unchanged after deleting source',maxdiff(z,w),0.d0)
!---Linkage for full local matrices is built once per pattern
CALL native_matrix_share(shared,shared2)
B2=>native_cast(shared2)
B2%M=B%M
CALL native_matrix_setup_full(B,x)
CALL native_matrix_setup_full(B2,x)
CALL check_exact('Linkage shared',MERGE(1.d0,0.d0,ASSOCIATED(B%linkage,B2%linkage)),1.d0)
CALL check_exact('Linkage owned by pattern',MERGE(1.d0,0.d0,ASSOCIATED(B2%linkage,B2%pattern%linkage)),1.d0)
CALL check_exact('Redundant rows shared',REAL(B2%nred-B%nred,r8),0.d0)
CALL B%update_slice
CALL B2%update_slice
v1=0.d0
IF(B%nnz>0)v1=MAXVAL(ABS(B%Mfull-B2%Mfull))
CALL check_exact('Full local matrices equal',oft_mpi_max(v1),0.d0)
CALL native_matrix_setup_full(L,x)
CALL native_matrix_setup_full(Mm,x)
CALL check_exact('FE operator linkage shared',MERGE(1.d0,0.d0,ASSOCIATED(L%linkage,Mm%linkage)),1.d0)
!---Complex matrix sharing a real pattern
CALL native_matrix_share(shared,cshared)
SELECT TYPE(cshared)
CLASS IS(oft_native_cmatrix)
  C=>cshared
CLASS DEFAULT
  CALL oft_abort('Native complex matrix expected','test_pattern',__FILE__)
END SELECT
CALL check_exact('Complex share uses real pattern',MERGE(1.d0,0.d0,ASSOCIATED(C%pattern,B%pattern)),1.d0)
C%M=CMPLX(B%M,0.d0,c8)
CALL x%new(cx)
CALL x%new(cz)
NULLIFY(xr,xi,cvals)
CALL x%get_local(xr)
CALL y%get_local(xi)
ALLOCATE(cvals(x%n))
cvals=CMPLX(xr,xi,c8)
CALL cx%restore_local(cvals)
CALL cshared%apply(cx,cz)
CALL cz%get_local(cvals)
CALL shared%apply(x,z)
CALL shared%apply(y,w)
CALL z%get_local(xr)
CALL w%get_local(xi)
! Products are bitwise equal before stitching; entries with three or more contributions are
! summed in ascending order of the complex value, not of each part, so only round-off agrees
CALL check('Re(complex shared A cx)',oft_mpi_max(MAXVAL(ABS(REAL(cvals,r8)-xr))),0.d0)
CALL check('Im(complex shared A cx)',oft_mpi_max(MAXVAL(ABS(AIMAG(cvals)-xi))),0.d0)
DEALLOCATE(xr,xi,cvals)
CALL cx%delete
CALL cz%delete
DEALLOCATE(cx,cz)
!---Pattern is freed with its last reference
CALL check_exact('Reference count with three users',REAL(B%pattern%refcount,r8),3.d0)
CALL shared%delete
DEALLOCATE(shared)
CALL cshared%delete
DEALLOCATE(cshared)
CALL check_exact('Pattern kept for last user',REAL(native_pattern_live-nlive1,r8),0.d0)
CALL shared2%delete
DEALLOCATE(shared2)
CALL check_exact('Pattern freed with last user',REAL(native_pattern_live-nlive1,r8),-1.d0)
!---FE pattern stays cached after its operators are deleted
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
CALL check_exact('FE pattern kept by cache',REAL(native_pattern_live-nlive0,r8),0.d0)
CALL x%delete
CALL y%delete
CALL z%delete
CALL w%delete
DEALLOCATE(x,y,z,w)
END SUBROUTINE test_pattern
! Cast a matrix to a native CRS matrix pointer (aborts otherwise)
FUNCTION native_cast(mat) RESULT(nmat)
CLASS(oft_matrix), POINTER, INTENT(in) :: mat
TYPE(oft_native_matrix), POINTER :: nmat
SELECT TYPE(mat)
TYPE IS(oft_native_matrix)
  nmat=>mat
CLASS DEFAULT
  CALL oft_abort('Native matrix expected','native_cast',__FILE__)
END SELECT
END FUNCTION native_cast
END PROGRAM test_native_la
