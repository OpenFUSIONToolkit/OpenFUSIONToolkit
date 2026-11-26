!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_replicated_seam.F90
!
!> Regression tests for replicated vectors (every processor stores all entries, each entry
!! owned by one processor).
!!
!! The current test cases are:
!! - Synchronization by ownership (stitch method 0)
!! - Summation of contributions (stitch method 1)
!! - Dot products for real and complex vectors (single and multiple)
!! - Slice extraction/restoration
!! - Eigenvalues of a diagonal system using (P)ARPACK
!!
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
PROGRAM test_replicated_seam
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_vector_ptr, oft_cvector_ptr
USE oft_native_la, ONLY: oft_native_vector, oft_native_matrix, native_vector_setup_replicated
USE oft_stitching, ONLY: oft_stitch_check
#ifdef HAVE_ARPACK
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
IMPLICIT NONE
INTEGER(i4) :: i,n_owned,iounit,ierr,nfail
INTEGER(i4), ALLOCATABLE :: owner(:)
REAL(r8) :: dval,dref
REAL(r8), POINTER, DIMENSION(:) :: vals,slice
REAL(r8), ALLOCATABLE :: dots(:)
COMPLEX(c8) :: cdval
COMPLEX(c8), POINTER, DIMENSION(:) :: cvals
COMPLEX(c8), ALLOCATABLE :: cdots(:)
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_cvector), POINTER :: cu,cv
TYPE(oft_vector_ptr) :: vlist(2)
TYPE(oft_cvector_ptr) :: cvlist(2)
INTEGER(i4) :: n = 1000
LOGICAL :: one_owner = .FALSE. ! All entries owned by the last processor (others own none)
LOGICAL :: skip_eigs = .FALSE. ! Skip the eigenvalue test
NAMELIST/test_replicated_options/n,one_owner,skip_eigs
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_replicated_options,IOSTAT=ierr)
CLOSE(iounit)
oft_env%pm=.FALSE.
nfail=0
!---Interleaved (non-contiguous) ownership, identical on all processors
ALLOCATE(owner(n))
DO i=1,n
  owner(i)=MOD((i-1)/7+3*(i-1),oft_env%nprocs)
  IF(one_owner)owner(i)=oft_env%nprocs-1
END DO
n_owned=COUNT(owner==oft_env%rank)
!---Create replicated vector
ALLOCATE(oft_native_vector::u)
SELECT TYPE(this=>u)
CLASS IS(oft_native_vector)
  CALL native_vector_setup_replicated(this,owner)
END SELECT
CALL oft_stitch_check(u%stitch_info)
CALL u%new(v)
IF(u%nslice/=n_owned.OR.u%ng/=n.OR.u%n/=n)CALL fail('Incorrect vector sizes')
NULLIFY(vals,slice,cvals)
ALLOCATE(vals(n))
!------------------------------------------------------------------------------
! Stitch method 0: owned values propagated to all processors
!------------------------------------------------------------------------------
DO i=1,n
  vals(i)=ref_val(i)
  IF(owner(i)/=oft_env%rank)vals(i)=-1.d3*(oft_env%rank+1) ! Garbage on non-owned entries
END DO
CALL u%restore_local(vals) ! Calls stitch(0)
CALL u%get_local(vals)
IF(ANY(vals/=[(ref_val(i),i=1,n)]))CALL fail('Stitch method 0')
!------------------------------------------------------------------------------
! Stitch method 1: sum of contributions from all processors
!------------------------------------------------------------------------------
DO i=1,n
  vals(i)=(oft_env%rank+1)*ref_val(i)
END DO
CALL v%set(0.d0)
CALL v%restore_local(vals,add=.TRUE.) ! Calls stitch(1)
CALL v%get_local(vals)
dref=REAL(oft_env%nprocs*(oft_env%nprocs+1)/2,8)
DO i=1,n
  IF(ABS(vals(i)-dref*ref_val(i))>1.d-12*ABS(dref*ref_val(i)))THEN
    CALL fail('Stitch method 1')
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Dot products
!------------------------------------------------------------------------------
dref=0.d0
DO i=1,n
  dref=dref+ref_val(i)**2
END DO
dval=u%dot(u)
IF(ABS(dval-dref)>1.d-12*dref)CALL fail('Real dot product')
CALL v%set(1.d0)
vlist(1)%f=>u
vlist(2)%f=>v
ALLOCATE(dots(2))
dots=u%mdot(vlist,2)
IF(ABS(dots(1)-dref)>1.d-12*dref)CALL fail('Real mdot (1)')
IF(ABS(dots(2)-SUM([(ref_val(i),i=1,n)]))>1.d-12*ABS(SUM([(ref_val(i),i=1,n)])))CALL fail('Real mdot (2)')
IF(ABS(v%sum()-REAL(n,8))>1.d-12*n)CALL fail('Real sum')
IF(ABS(v%norm(2)-SQRT(REAL(n,8)))>1.d-12*SQRT(REAL(n,8)))CALL fail('Real 2-norm')
!---Complex vectors
CALL u%new(cu)
CALL u%new(cv)
ALLOCATE(cvals(n))
DO i=1,n
  cvals(i)=CMPLX(ref_val(i),2.d0*ref_val(i),c8)
END DO
CALL cu%restore_local(cvals)
CALL cv%set((0.d0,1.d0))
cdval=cu%dot(cu)
IF(ABS(cdval-5.d0*dref)>5.d-12*dref)CALL fail('Complex dot product')
cvlist(1)%f=>cu
cvlist(2)%f=>cv
ALLOCATE(cdots(2))
cdots=cu%mdot(cvlist,2)
IF(ABS(cdots(1)-5.d0*dref)>5.d-12*dref)CALL fail('Complex mdot (1)')
cdval=SUM(CONJG(cvals))*(0.d0,1.d0)
IF(ABS(cdots(2)-cdval)>1.d-12*ABS(cdval))CALL fail('Complex mdot (2)')
!---Complex stitch method 1
DO i=1,n
  cvals(i)=(oft_env%rank+1)*CMPLX(ref_val(i),-ref_val(i),c8)
END DO
CALL cv%set((0.d0,0.d0))
CALL cv%restore_local(cvals,add=.TRUE.) ! Calls stitch(1)
CALL cv%get_local(cvals)
dref=REAL(oft_env%nprocs*(oft_env%nprocs+1)/2,8)
DO i=1,n
  IF(ABS(cvals(i)-dref*CMPLX(ref_val(i),-ref_val(i),c8))>1.d-12*ABS(dref*ref_val(i)))THEN
    CALL fail('Complex stitch method 1')
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Slices: owned entries only, full vector restored on all processors
!------------------------------------------------------------------------------
CALL u%get_slice(slice)
IF(SIZE(slice)/=n_owned)CALL fail('Slice size')
CALL u%set(0.d0)
CALL u%restore_slice(slice)
CALL u%get_local(vals)
IF(ANY(vals/=[(ref_val(i),i=1,n)]))CALL fail('Slice round trip')
!------------------------------------------------------------------------------
! Eigenvalues of a diagonal system
!------------------------------------------------------------------------------
#ifdef HAVE_ARPACK
IF((n>=20).AND.(.NOT.skip_eigs))CALL test_eigs() ! Ritz space must fit in the system
#endif
!---Report
nfail=oft_mpi_max(nfail)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=iounit,FILE='replicated.results')
  IF(nfail==0)THEN
    WRITE(iounit,'(A)')'PASS'
  ELSE
    WRITE(iounit,'(A,I4)')'FAIL',nfail
  END IF
  CLOSE(iounit)
END IF
CALL oft_finalize
CONTAINS
!> Reference value for each entry
FUNCTION ref_val(i) RESULT(val)
INTEGER(i4), INTENT(in) :: i
REAL(r8) :: val
val=SIN(REAL(i,8))+1.5d0
END FUNCTION ref_val
!> Record a failure
SUBROUTINE fail(msg)
CHARACTER(LEN=*), INTENT(in) :: msg
WRITE(*,'(A,I4,2A)')'FAILED (rank',oft_env%rank,'): ',msg
nfail=nfail+1
END SUBROUTINE fail
#ifdef HAVE_ARPACK
!> Generalized eigenvalues of A x = lam M x, with A diagonal and M=2*I (replicated, local matrices)
SUBROUTINE test_eigs()
INTEGER(i4) :: j
REAL(r8) :: lam
TYPE(oft_native_matrix), TARGET :: amat,mmat
CLASS(oft_solver), POINTER :: minv
TYPE(oft_irlm_eigsolver) :: arsolver
!---Four well separated eigenvalues (1000, 2000, 3000, 4000) above a cluster in (0,1]
CALL setup_diag(amat,[(REAL(i,8)/REAL(n,8),i=1,n-4),1.d3,2.d3,3.d3,4.d3])
CALL setup_diag(mmat,[(2.d0,i=1,n)])
CALL create_cg_solver(minv)
minv%A=>mmat
minv%its=-2
minv%atol=1.d-14
minv%rtol=1.d-14
CALL create_diag_pre(minv%pre)
arsolver%A=>amat
arsolver%M=>mmat
arsolver%Minv=>minv
arsolver%mode=2
arsolver%nev=4
arsolver%its=1000
arsolver%tol=1.d-10
CALL u%set(1.d0)
CALL arsolver%apply(u,lam)
IF(arsolver%info<0)THEN
  CALL fail('ARPACK error')
  RETURN
END IF
DO j=1,4
  dref=REAL(j,8)*5.d2
  IF(MINVAL(ABS(arsolver%eig_val(1,1:4)-dref))>1.d-8*dref)THEN
    CALL fail('Eigenvalue')
    WRITE(*,*)'  Computed: ',arsolver%eig_val(1,1:4)
  END IF
END DO
!---Eigenvectors are returned as owned slices
IF(SIZE(arsolver%eig_vec,DIM=1)/=n_owned)CALL fail('Eigenvector slice size')
CALL arsolver%delete()
CALL minv%pre%delete()
CALL minv%delete()
END SUBROUTINE test_eigs
!> Setup a replicated (local) diagonal matrix
SUBROUTINE setup_diag(mat,d)
TYPE(oft_native_matrix), INTENT(inout) :: mat
REAL(r8), INTENT(in) :: d(:)
INTEGER(i4) :: j
mat%nr=n; mat%nrg=n
mat%nc=n; mat%ncg=n
mat%nnz=n
ALLOCATE(mat%kr(n+1),mat%lc(n),mat%M(n))
DO j=1,n
  mat%kr(j)=j
  mat%lc(j)=j
  mat%M(j)=d(j)
END DO
mat%kr(n+1)=n+1
mat%force_local=.TRUE. ! Full rows are stored on every processor
CALL mat%assemble(u)
END SUBROUTINE setup_diag
#endif
END PROGRAM test_replicated_seam
