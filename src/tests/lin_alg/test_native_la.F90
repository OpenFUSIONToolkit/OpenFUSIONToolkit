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
!! - Complex vector products and mixed real/complex dot products (serial only)
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_native_la
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix
USE oft_native_la, ONLY: oft_native_matrix
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_getlop
IMPLICIT NONE
INTEGER(i4) :: io_unit,nfail
REAL(r8), PARAMETER :: tol=1.d-12
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,2,ML_oft_lagrange,minlev=-1)
!---Run tests
nfail=0
CALL test_ops()
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
!---Complex vectors (complex stitching is only supported for local vectors)
IF(oft_env%nprocs==1)THEN
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
  !---Complex norms
  v1=0.d0
  v2=0.d0
  DO i=1,cz%n
    v1=v1+ABS(cvals(i))
    v2=MAX(v2,ABS(cvals(i)))
  END DO
  CALL check('Complex 1-norm',cz%norm(1),v1)
  CALL check('Complex inf-norm',cz%norm(3),v2)
  DEALLOCATE(xr,xi,cvals)
  CALL cx%delete
  CALL cz%delete
  DEALLOCATE(cx,cz)
END IF
!---Cleanup
CALL x%delete
CALL y%delete
CALL z%delete
CALL w%delete
DEALLOCATE(x,y,z,w)
CALL lop%delete
DEALLOCATE(lop)
END SUBROUTINE test_ops
END PROGRAM test_native_la
