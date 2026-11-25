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
!! - Zeroing rows and columns of block matrices (real and complex), with and without
!!   keeping diagonal entries, against the expected entries of the CSR structure
!! - Synchronization by ownership after ownership is set or changed on a seam whose
!!   stitching plan was already built
!! - maxpy and mdot with non-native vectors in the list (copied) against native vectors
!!
!! @ingroup testing
!------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! Minimal non-native vector types for testing the fallbacks of native operations
! (each procedure forwards to a wrapped native vector)
!---------------------------------------------------------------------------------
MODULE test_wrap_vectors
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_vector_ptr, oft_cvector_ptr
IMPLICIT NONE
PRIVATE
!> Real vector that is not a native vector
TYPE, PUBLIC, EXTENDS(oft_vector) :: wrap_vector
  CLASS(oft_vector), POINTER :: inner => NULL() !< Wrapped vector
CONTAINS
  PROCEDURE :: new_real => wv_new_real
  PROCEDURE :: new_complex => wv_new_complex
  PROCEDURE :: set => wv_set
  PROCEDURE :: get_local => wv_get_local
  PROCEDURE :: restore_local => wv_restore_local
  PROCEDURE :: get_slice => wv_get_slice
  PROCEDURE :: restore_slice => wv_restore_slice
  PROCEDURE :: add => wv_add
  PROCEDURE :: mult => wv_mult
  PROCEDURE :: dot_real => wv_dot_real
  PROCEDURE :: dot_complex => wv_dot_complex
  PROCEDURE :: mdot_real => wv_mdot_real
  PROCEDURE :: mdot_complex => wv_mdot_complex
  PROCEDURE :: scale => wv_scale
  PROCEDURE :: sum => wv_sum
  PROCEDURE :: norm => wv_norm
  PROCEDURE :: stitch => wv_stitch
END TYPE wrap_vector
!> Complex vector that is not a native vector
TYPE, PUBLIC, EXTENDS(oft_cvector) :: wrap_cvector
  CLASS(oft_cvector), POINTER :: inner => NULL() !< Wrapped vector
CONTAINS
  PROCEDURE :: new_real => wc_new_real
  PROCEDURE :: new_complex => wc_new_complex
  PROCEDURE :: set => wc_set
  PROCEDURE :: get_local => wc_get_local
  PROCEDURE :: restore_local => wc_restore_local
  PROCEDURE :: get_slice => wc_get_slice
  PROCEDURE :: restore_slice => wc_restore_slice
  PROCEDURE :: add_real => wc_add_real
  PROCEDURE :: add_complex => wc_add_complex
  PROCEDURE :: mult_real => wc_mult_real
  PROCEDURE :: mult_complex => wc_mult_complex
  PROCEDURE :: dot_real => wc_dot_real
  PROCEDURE :: dot_complex => wc_dot_complex
  PROCEDURE :: mdot_real => wc_mdot_real
  PROCEDURE :: mdot_complex => wc_mdot_complex
  PROCEDURE :: scale => wc_scale
  PROCEDURE :: sum => wc_sum
  PROCEDURE :: norm => wc_norm
  PROCEDURE :: stitch => wc_stitch
END TYPE wrap_cvector
PUBLIC wrap_vector_create, wrap_cvector_create
CONTAINS
! Wrap a real vector (sizes and seam are taken from it)
SUBROUTINE wrap_vector_create(inner,new)
CLASS(oft_vector), POINTER, INTENT(in) :: inner
CLASS(oft_vector), POINTER, INTENT(out) :: new
ALLOCATE(wrap_vector::new)
SELECT TYPE(new)
TYPE IS(wrap_vector)
  new%inner=>inner
  new%n=inner%n
  new%nslice=inner%nslice
  new%nblocks=inner%nblocks
  new%ng=inner%ng
  new%map=>inner%map
  new%stitch_info=>inner%stitch_info
END SELECT
END SUBROUTINE wrap_vector_create
! Wrap a complex vector (sizes and seam are taken from it)
SUBROUTINE wrap_cvector_create(inner,new)
CLASS(oft_cvector), POINTER, INTENT(in) :: inner
CLASS(oft_cvector), POINTER, INTENT(out) :: new
ALLOCATE(wrap_cvector::new)
SELECT TYPE(new)
TYPE IS(wrap_cvector)
  new%inner=>inner
  new%n=inner%n
  new%nslice=inner%nslice
  new%nblocks=inner%nblocks
  new%ng=inner%ng
  new%map=>inner%map
  new%stitch_info=>inner%stitch_info
END SELECT
END SUBROUTINE wrap_cvector_create
!---Real wrapper
SUBROUTINE wv_new_real(self,new)
CLASS(wrap_vector), INTENT(in) :: self
CLASS(oft_vector), POINTER, INTENT(out) :: new
CLASS(oft_vector), POINTER :: tmp
CALL self%inner%new(tmp)
CALL wrap_vector_create(tmp,new)
END SUBROUTINE wv_new_real
SUBROUTINE wv_new_complex(self,new)
CLASS(wrap_vector), INTENT(in) :: self
CLASS(oft_cvector), POINTER, INTENT(out) :: new
CALL self%inner%new(new)
END SUBROUTINE wv_new_complex
SUBROUTINE wv_set(self,alpha,iblock,random)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), INTENT(in) :: alpha
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: random
CALL self%inner%set(alpha,iblock,random)
END SUBROUTINE wv_set
SUBROUTINE wv_get_local(self,array,iblock)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), POINTER, INTENT(inout) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
CALL self%inner%get_local(array,iblock)
END SUBROUTINE wv_get_local
SUBROUTINE wv_restore_local(self,array,iblock,add,wait)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), INTENT(in) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: add,wait
CALL self%inner%restore_local(array,iblock,add,wait)
END SUBROUTINE wv_restore_local
SUBROUTINE wv_get_slice(self,array,iblock)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), POINTER, INTENT(inout) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
CALL self%inner%get_slice(array,iblock)
END SUBROUTINE wv_get_slice
SUBROUTINE wv_restore_slice(self,array,iblock,wait)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), INTENT(in) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: wait
CALL self%inner%restore_slice(array,iblock,wait)
END SUBROUTINE wv_restore_slice
SUBROUTINE wv_add(self,gamma,alpha,a,beta,b)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), INTENT(in) :: gamma,alpha
CLASS(oft_vector), TARGET, INTENT(inout) :: a
REAL(r8), OPTIONAL, INTENT(in) :: beta
CLASS(oft_vector), TARGET, OPTIONAL, INTENT(inout) :: b
CALL self%inner%add(gamma,alpha,a,beta,b)
END SUBROUTINE wv_add
SUBROUTINE wv_mult(self,a,div_flag)
CLASS(wrap_vector), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a
LOGICAL, OPTIONAL, INTENT(in) :: div_flag
CALL self%inner%mult(a,div_flag)
END SUBROUTINE wv_mult
FUNCTION wv_dot_real(self,a) RESULT(dot)
CLASS(wrap_vector), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a
REAL(r8) :: dot
dot=self%inner%dot(a)
END FUNCTION wv_dot_real
FUNCTION wv_dot_complex(self,a) RESULT(dot)
CLASS(wrap_vector), INTENT(inout) :: self
CLASS(oft_cvector), TARGET, INTENT(inout) :: a
COMPLEX(c8) :: dot
dot=self%inner%dot(a)
END FUNCTION wv_dot_complex
FUNCTION wv_mdot_real(self,a,n) RESULT(dots)
CLASS(wrap_vector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: n
TYPE(oft_vector_ptr), INTENT(inout) :: a(n)
REAL(r8) :: dots(n)
dots=self%inner%mdot(a,n)
END FUNCTION wv_mdot_real
FUNCTION wv_mdot_complex(self,a,n) RESULT(dots)
CLASS(wrap_vector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: n
TYPE(oft_cvector_ptr), INTENT(inout) :: a(n)
COMPLEX(c8) :: dots(n)
dots=self%inner%mdot(a,n)
END FUNCTION wv_mdot_complex
SUBROUTINE wv_scale(self,alpha)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8), INTENT(in) :: alpha
CALL self%inner%scale(alpha)
END SUBROUTINE wv_scale
FUNCTION wv_sum(self) RESULT(sum)
CLASS(wrap_vector), INTENT(inout) :: self
REAL(r8) :: sum
sum=self%inner%sum()
END FUNCTION wv_sum
FUNCTION wv_norm(self,itype) RESULT(norm)
CLASS(wrap_vector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: itype
REAL(r8) :: norm
norm=self%inner%norm(itype)
END FUNCTION wv_norm
SUBROUTINE wv_stitch(self,up_method)
CLASS(wrap_vector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: up_method
CALL self%inner%stitch(up_method)
END SUBROUTINE wv_stitch
!---Complex wrapper
SUBROUTINE wc_new_real(self,new)
CLASS(wrap_cvector), INTENT(in) :: self
CLASS(oft_vector), POINTER, INTENT(out) :: new
CALL self%inner%new(new)
END SUBROUTINE wc_new_real
SUBROUTINE wc_new_complex(self,new)
CLASS(wrap_cvector), INTENT(in) :: self
CLASS(oft_cvector), POINTER, INTENT(out) :: new
CLASS(oft_cvector), POINTER :: tmp
CALL self%inner%new(tmp)
CALL wrap_cvector_create(tmp,new)
END SUBROUTINE wc_new_complex
SUBROUTINE wc_set(self,alpha,iblock,random)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: alpha
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: random
CALL self%inner%set(alpha,iblock,random)
END SUBROUTINE wc_set
SUBROUTINE wc_get_local(self,array,iblock)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), POINTER, INTENT(inout) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
CALL self%inner%get_local(array,iblock)
END SUBROUTINE wc_get_local
SUBROUTINE wc_restore_local(self,array,iblock,add,wait)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: add,wait
CALL self%inner%restore_local(array,iblock,add,wait)
END SUBROUTINE wc_restore_local
SUBROUTINE wc_get_slice(self,array,iblock)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), POINTER, INTENT(inout) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
CALL self%inner%get_slice(array,iblock)
END SUBROUTINE wc_get_slice
SUBROUTINE wc_restore_slice(self,array,iblock,wait)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: array(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
LOGICAL, OPTIONAL, INTENT(in) :: wait
CALL self%inner%restore_slice(array,iblock,wait)
END SUBROUTINE wc_restore_slice
SUBROUTINE wc_add_real(self,gamma,alpha,a,beta,b)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: gamma,alpha
CLASS(oft_vector), TARGET, INTENT(inout) :: a
COMPLEX(c8), OPTIONAL, INTENT(in) :: beta
CLASS(oft_vector), TARGET, OPTIONAL, INTENT(inout) :: b
CALL self%inner%add(gamma,alpha,a,beta,b)
END SUBROUTINE wc_add_real
SUBROUTINE wc_add_complex(self,gamma,alpha,a,beta,b)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: gamma,alpha
CLASS(oft_cvector), TARGET, INTENT(inout) :: a
COMPLEX(c8), OPTIONAL, INTENT(in) :: beta
CLASS(oft_cvector), TARGET, OPTIONAL, INTENT(inout) :: b
CALL self%inner%add(gamma,alpha,a,beta,b)
END SUBROUTINE wc_add_complex
SUBROUTINE wc_mult_real(self,a,div_flag)
CLASS(wrap_cvector), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a
LOGICAL, OPTIONAL, INTENT(in) :: div_flag
CALL self%inner%mult(a,div_flag)
END SUBROUTINE wc_mult_real
SUBROUTINE wc_mult_complex(self,a,div_flag)
CLASS(wrap_cvector), INTENT(inout) :: self
CLASS(oft_cvector), TARGET, INTENT(inout) :: a
LOGICAL, OPTIONAL, INTENT(in) :: div_flag
CALL self%inner%mult(a,div_flag)
END SUBROUTINE wc_mult_complex
FUNCTION wc_dot_real(self,a) RESULT(dot)
CLASS(wrap_cvector), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a
COMPLEX(c8) :: dot
dot=self%inner%dot(a)
END FUNCTION wc_dot_real
FUNCTION wc_dot_complex(self,a) RESULT(dot)
CLASS(wrap_cvector), INTENT(inout) :: self
CLASS(oft_cvector), TARGET, INTENT(inout) :: a
COMPLEX(c8) :: dot
dot=self%inner%dot(a)
END FUNCTION wc_dot_complex
FUNCTION wc_mdot_real(self,a,n) RESULT(dots)
CLASS(wrap_cvector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: n
TYPE(oft_vector_ptr), INTENT(inout) :: a(n)
COMPLEX(c8) :: dots(n)
dots=self%inner%mdot(a,n)
END FUNCTION wc_mdot_real
FUNCTION wc_mdot_complex(self,a,n) RESULT(dots)
CLASS(wrap_cvector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: n
TYPE(oft_cvector_ptr), INTENT(inout) :: a(n)
COMPLEX(c8) :: dots(n)
dots=self%inner%mdot(a,n)
END FUNCTION wc_mdot_complex
SUBROUTINE wc_scale(self,alpha)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8), INTENT(in) :: alpha
CALL self%inner%scale(alpha)
END SUBROUTINE wc_scale
FUNCTION wc_sum(self) RESULT(sum)
CLASS(wrap_cvector), INTENT(inout) :: self
COMPLEX(c8) :: sum
sum=self%inner%sum()
END FUNCTION wc_sum
FUNCTION wc_norm(self,itype) RESULT(norm)
CLASS(wrap_cvector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: itype
REAL(r8) :: norm
norm=self%inner%norm(itype)
END FUNCTION wc_norm
SUBROUTINE wc_stitch(self,up_method)
CLASS(wrap_cvector), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: up_method
CALL self%inner%stitch(up_method)
END SUBROUTINE wc_stitch
END MODULE test_wrap_vectors
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
USE oft_stitching, ONLY: oft_seam, oft_global_stitch, seam_free_plan
USE fem_base, ONLY: oft_ml_fem_type, oft_afem_type
USE test_wrap_vectors, ONLY: wrap_vector_create, wrap_cvector_create
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
CALL test_nonnative()
CALL test_bc()
CALL test_stitch()
CALL test_stitch_owner()
CALL test_accumulate()
CALL test_pattern()
CALL test_zero()
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
! Check that a difference is within a few ulp of a scale (results whose rounding depends on
! how the compiler contracts or vectorizes the arithmetic, eg. FMA)
SUBROUTINE check_ulp(name,diff,scale)
CHARACTER(LEN=*), INTENT(in) :: name
REAL(r8), INTENT(in) :: diff,scale
IF(.NOT.(diff<=8.d0*EPSILON(1.d0)*scale))THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,'(3A,2ES25.16)')'FAILED: ',name,' ',diff,scale
ELSE
  IF(oft_env%head_proc)WRITE(*,'(3A,ES25.16)')'PASSED: ',name,' ',diff
END IF
END SUBROUTINE check_ulp
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
! Maximum absolute value of a (possibly empty) array
PURE FUNCTION maxabs(d) RESULT(m)
REAL(r8), INTENT(in) :: d(:)
REAL(r8) :: m
m=0.d0
IF(SIZE(d)>0)m=MAXVAL(ABS(d))
END FUNCTION maxabs
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
!---maxpy matches successive add calls to within rounding (same order of operations, but
!---compilers may contract the multiply-adds differently)
CALL y%copy(x)
CALL y%maxpy(alpha,vs,nv)
DO k=1,nv
  CALL x%add(1.d0,alpha(k),vs(k)%f)
END DO
CALL check_ulp('maxpy = successive add',maxdiff(x,y),x%norm(3))
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
CALL check_ulp('Complex maxpy = successive add',cmaxdiff(cx,cy),cx%norm(3))
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
        ELSE IF(m==2)THEN
          ! Sums of remote contributions, which the compiler may vectorize differently for
          ! real and complex values
          d=MAXVAL(ABS(cx%v-CMPLX(x%v,y%v,c8)))
          CALL check_ulp('Complex stitch(2) matches real/imaginary parts',oft_mpi_max(d),cx%norm(3))
        ELSE
          d=MAXVAL(ABS(cx%v-CMPLX(x%v,y%v,c8)))
          CALL check_exact('Complex stitch(0) matches real/imaginary parts',oft_mpi_max(d),0.d0)
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
! maxpy and mdot with non-native vectors in the list, which are copied, give the same
! results as with the native vectors they wrap (same kernels on the same values)
SUBROUTINE test_nonnative()
INTEGER(i4), PARAMETER :: nv=4
INTEGER(i4) :: k
REAL(r8) :: alpha(nv),dref(nv),dmix(nv)
COMPLEX(c8) :: calpha(nv),cdref(nv),cdmix(nv)
TYPE(oft_vector_ptr) :: vs(nv),ws(nv)
TYPE(oft_cvector_ptr) :: cvs(nv),cws(nv)
CLASS(oft_vector), POINTER :: x,y
CLASS(oft_cvector), POINTER :: cx,cy
!---Real vectors; every other one wrapped in a non-native vector type
DO k=1,nv
  CALL ML_oft_lagrange%vec_create(vs(k)%f)
  CALL fill(vs(k)%f,0.13d0*REAL(k,r8),0.1d0*REAL(k-2,r8))
  alpha(k)=0.3d0*REAL(k,r8)-0.7d0
  ws(k)%f=>vs(k)%f
END DO
DO k=2,nv,2
  CALL wrap_vector_create(vs(k)%f,ws(k)%f)
END DO
CALL ML_oft_lagrange%vec_create(x)
CALL ML_oft_lagrange%vec_create(y)
CALL fill(x,0.29d0,-0.05d0)
CALL y%copy(x)
CALL x%maxpy(alpha,vs,nv)
CALL y%maxpy(alpha,ws,nv)
CALL check_exact('maxpy with non-native vectors',maxdiff(x,y),0.d0)
dref=x%mdot(vs,nv)
dmix=x%mdot(ws,nv)
CALL check_exact('mdot with non-native vectors',MAXVAL(ABS(dmix-dref)),0.d0)
!---Complex vectors; every other one wrapped
DO k=1,nv
  CALL x%new(cvs(k)%f)
  CALL cfill(cvs(k)%f,0.17d0*REAL(k,r8),0.2d0*REAL(k-2,r8))
  calpha(k)=CMPLX(alpha(k),0.1d0*REAL(k,r8),c8)
  cws(k)%f=>cvs(k)%f
END DO
DO k=2,nv,2
  CALL wrap_cvector_create(cvs(k)%f,cws(k)%f)
END DO
CALL x%new(cx)
CALL x%new(cy)
CALL cfill(cx,0.31d0,0.07d0)
CALL cy%copy(cx)
CALL cx%maxpy(calpha,cvs,nv)
CALL cy%maxpy(calpha,cws,nv)
CALL check_exact('Complex maxpy with non-native vectors',cmaxdiff(cx,cy),0.d0)
cdref=cx%mdot(cvs,nv)
cdmix=cx%mdot(cws,nv)
CALL check_exact('Complex mdot with non-native vectors',MAXVAL(ABS(cdmix-cdref)),0.d0)
!---Mixed real/complex mdot
cdref=x%mdot(cvs,nv)
cdmix=x%mdot(cws,nv)
CALL check_exact('Real-complex mdot with non-native vectors',MAXVAL(ABS(cdmix-cdref)),0.d0)
cdref=cx%mdot(vs,nv)
cdmix=cx%mdot(ws,nv)
CALL check_exact('Complex-real mdot with non-native vectors',MAXVAL(ABS(cdmix-cdref)),0.d0)
!---Cleanup (wrappers only reference the native vectors)
DO k=2,nv,2
  DEALLOCATE(ws(k)%f,cws(k)%f)
END DO
DO k=1,nv
  CALL vs(k)%f%delete
  CALL cvs(k)%f%delete
  DEALLOCATE(vs(k)%f,cvs(k)%f)
END DO
CALL x%delete
CALL y%delete
CALL cx%delete
CALL cy%delete
DEALLOCATE(x,y,cx,cy)
END SUBROUTINE test_nonnative
! Synchronization by ownership uses the ownership current at the time of the stitch, also
! when the stitching plan was built before ownership was set or changed
SUBROUTINE test_stitch_owner()
INTEGER(i4) :: i
REAL(r8), ALLOCATABLE, DIMENSION(:) :: raw,ref,a
LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: leo_orig,leo_tmp
TYPE(oft_seam) :: sc
CLASS(oft_vector), POINTER :: x
CALL ML_oft_lagrange%vec_create(x)
SELECT TYPE(x)
TYPE IS(oft_native_vector)
  ALLOCATE(raw(x%n),ref(x%n),a(x%n))
  DO i=1,x%n
    raw(i)=gval(INT(i+1000*oft_env%rank,i8),0.37d0,0.2d0)
  END DO
  !---Reference: synchronization with the vector's own seam
  ref=raw
  CALL oft_global_stitch(x%stitch_info,ref,0)
  !---Copy of the seam (shares linkage) whose ownership is changed below
  sc=x%stitch_info
  CALL seam_free_plan(sc)
  leo_orig=>x%stitch_info%leo
  !---Plan built by summation before ownership is set
  NULLIFY(sc%leo)
  a=raw
  CALL oft_global_stitch(sc,a,1)
  sc%leo=>leo_orig
  a=raw
  CALL oft_global_stitch(sc,a,0)
  CALL check_exact('Stitch(0) after ownership set',oft_mpi_max(maxabs(a-ref)),0.d0)
  !---Ownership changed in place after the plan was built
  IF(ASSOCIATED(leo_orig))THEN
    ALLOCATE(leo_tmp(SIZE(leo_orig)))
    leo_tmp=.TRUE.
    sc%leo=>leo_tmp
    a=raw
    CALL oft_global_stitch(sc,a,0)
    leo_tmp=leo_orig
    a=raw
    CALL oft_global_stitch(sc,a,0)
    CALL check_exact('Stitch(0) after ownership changed in place',oft_mpi_max(maxabs(a-ref)),0.d0)
    DEALLOCATE(leo_tmp)
  END IF
  CALL seam_free_plan(sc)
  DEALLOCATE(raw,ref,a)
CLASS DEFAULT
  CALL oft_abort('Native vector expected','test_stitch_owner',__FILE__)
END SELECT
CALL x%delete
DEALLOCATE(x)
END SUBROUTINE test_stitch_owner
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
! Zero rows and columns of block CSR matrices and compare every entry with the expected
! result: zeroed entries in the selected rows/columns of the selected block, original values
! elsewhere and, with keep_diag, on the diagonal
SUBROUTINE test_zero()
INTEGER(i4) :: i,j,k,ne,nsel,roff,coff
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: isel
REAL(r8), ALLOCATABLE, DIMENSION(:) :: mref
CLASS(oft_vector), POINTER :: xs,xv
CLASS(oft_matrix), POINTER :: sq => NULL()
CLASS(oft_matrix), POINTER :: rect => NULL()
CLASS(oft_cmatrix), POINTER :: csq => NULL()
TYPE(oft_graph_ptr) :: graphs(3,3)
TYPE(oft_native_matrix), POINTER :: A,R
TYPE(oft_native_cmatrix), POINTER :: C
CALL ML_oft_lagrange%vec_create(xs)
CALL ML_oft_vlagrange%vec_create(xv)
ne=ML_oft_lagrange%current_level%ne
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=ne
graphs(1,1)%g%nrg=ML_oft_lagrange%current_level%global%ne
graphs(1,1)%g%nc=ne
graphs(1,1)%g%ncg=ML_oft_lagrange%current_level%global%ne
graphs(1,1)%g%nnz=ML_oft_lagrange%current_level%nee
graphs(1,1)%g%kr=>ML_oft_lagrange%current_level%kee
graphs(1,1)%g%lc=>ML_oft_lagrange%current_level%lee
DO i=1,3
  DO j=1,3
    graphs(i,j)%g=>graphs(1,1)%g
  END DO
END DO
!---Square 3x3 block matrix (vector Lagrange) and rectangular 1x3 block matrix
CALL create_matrix(sq,graphs,xv,xv)
CALL create_matrix(rect,graphs(1:1,:),xs,xv)
A=>native_cast(sq)
R=>native_cast(rect)
!---Rows and columns selected in each block (block-local indices)
nsel=MIN(5,ne)
ALLOCATE(isel(nsel))
DO i=1,nsel
  isel(i)=1+((i-1)*ne)/nsel
END DO
ALLOCATE(mref(A%nnz))
DO i=1,A%nnz
  mref(i)=1.d0+REAL(i,r8)
END DO
!---Rows of block 2, keeping diagonal entries
roff=A%i_map(2)%offset
A%M=mref
CALL sq%zero_rows(nsel,isel,2,keep_diag=.TRUE.)
CALL check_exact('zero_rows block, keep_diag',REAL(oft_mpi_sum(zero_errs(A%nr,A%kr,A%lc,mref,A%M,roff,nsel,isel,.TRUE.,.TRUE.)),r8),0.d0)
!---Rows of block 2, including diagonal entries
A%M=mref
CALL sq%zero_rows(nsel,isel,2)
CALL check_exact('zero_rows block',REAL(oft_mpi_sum(zero_errs(A%nr,A%kr,A%lc,mref,A%M,roff,nsel,isel,.TRUE.,.FALSE.)),r8),0.d0)
!---Columns of block 3, keeping diagonal entries
coff=A%j_map(3)%offset
A%M=mref
CALL A%zero_cols(nsel,isel,3,keep_diag=.TRUE.)
CALL check_exact('zero_cols block, keep_diag',REAL(oft_mpi_sum(zero_errs(A%nr,A%kr,A%lc,mref,A%M,coff,nsel,isel,.FALSE.,.TRUE.)),r8),0.d0)
!---Columns of block 2 of the rectangular matrix (column offset differs from any row offset)
DEALLOCATE(mref)
ALLOCATE(mref(R%nnz))
DO i=1,R%nnz
  mref(i)=1.d0+REAL(i,r8)
END DO
coff=R%j_map(2)%offset
R%M=mref
CALL R%zero_cols(nsel,isel,2)
CALL check_exact('zero_cols rectangular block',REAL(oft_mpi_sum(zero_errs(R%nr,R%kr,R%lc,mref,R%M,coff,nsel,isel,.FALSE.,.FALSE.)),r8),0.d0)
!---Complex matrix with the same pattern
DEALLOCATE(mref)
ALLOCATE(mref(A%nnz))
DO i=1,A%nnz
  mref(i)=1.d0+REAL(i,r8)
END DO
CALL native_matrix_share(sq,csq)
SELECT TYPE(csq)
CLASS IS(oft_native_cmatrix)
  C=>csq
CLASS DEFAULT
  CALL oft_abort('Native complex matrix expected','test_zero',__FILE__)
END SELECT
roff=C%i_map(2)%offset
C%M=CMPLX(mref,-mref,c8)
CALL csq%zero_rows(nsel,isel,2,keep_diag=.TRUE.)
k=zero_errs(C%nr,C%kr,C%lc,mref,REAL(C%M,r8),roff,nsel,isel,.TRUE.,.TRUE.) &
  +zero_errs(C%nr,C%kr,C%lc,-mref,AIMAG(C%M),roff,nsel,isel,.TRUE.,.TRUE.)
CALL check_exact('Complex zero_rows block, keep_diag',REAL(oft_mpi_sum(k),r8),0.d0)
coff=C%j_map(3)%offset
C%M=CMPLX(mref,-mref,c8)
CALL C%zero_cols(nsel,isel,3,keep_diag=.TRUE.)
k=zero_errs(C%nr,C%kr,C%lc,mref,REAL(C%M,r8),coff,nsel,isel,.FALSE.,.TRUE.) &
  +zero_errs(C%nr,C%kr,C%lc,-mref,AIMAG(C%M),coff,nsel,isel,.FALSE.,.TRUE.)
CALL check_exact('Complex zero_cols block, keep_diag',REAL(oft_mpi_sum(k),r8),0.d0)
DEALLOCATE(isel,mref)
CALL csq%delete
CALL sq%delete
CALL rect%delete
DEALLOCATE(csq,sq,rect)
DEALLOCATE(graphs(1,1)%g)
CALL xs%delete
CALL xv%delete
DEALLOCATE(xs,xv)
END SUBROUTINE test_zero
! Number of CSR entries that differ from the expected result of zeroing rows (`rows`) or
! columns `off+isel(:)`, with the diagonal entries kept if `keep`
FUNCTION zero_errs(nr,kr,lc,mref,mval,off,nsel,isel,rows,keep) RESULT(nerr)
INTEGER(i4), INTENT(in) :: nr,off,nsel
INTEGER(i4), INTENT(in) :: kr(:),lc(:),isel(:)
REAL(r8), INTENT(in) :: mref(:),mval(:)
LOGICAL, INTENT(in) :: rows,keep
INTEGER(i4) :: nerr
INTEGER(i4) :: i,j
LOGICAL :: zeroed
nerr=0
DO i=1,nr
  DO j=kr(i),kr(i+1)-1
    IF(rows)THEN
      zeroed=ANY(off+isel(1:nsel)==i)
    ELSE
      zeroed=ANY(off+isel(1:nsel)==lc(j))
    END IF
    IF(keep.AND.(lc(j)==i))zeroed=.FALSE.
    IF(zeroed)THEN
      IF(mval(j)/=0.d0)nerr=nerr+1
    ELSE
      IF(mval(j)/=mref(j))nerr=nerr+1
    END IF
  END DO
END DO
END FUNCTION zero_errs
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
