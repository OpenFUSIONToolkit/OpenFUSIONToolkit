!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_native_la.F90
!
!> @defgroup doxy_oft_native_la Native backend
!! Native vector and matrix backend
!! @ingroup doxy_oft_lin_alg
!
!> Abstract field interfaces and native vector implementations
!!
!! Abstract interface definitions
!! - Abstract field class
!!
!! Native vector implementations
!! - Vector class
!! - Composite field class
!!
!! @sa oft_petsc_vectors
!!
!! @authors Chris Hansen
!! @date Feburary 2012
!! @ingroup doxy_oft_native_la
!---------------------------------------------------------------------------------
MODULE oft_native_la
USE, INTRINSIC :: iso_c_binding, only: c_int
USE oft_local
USE oft_base
USE oft_sort, ONLY: search_array, sort_array
USE oft_stitching, ONLY: oft_seam, oft_global_stitch, oft_global_dp, &
  oft_global_reduction, oft_stitch_check
USE oft_io
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_cvector_ptr, &
  oft_matrix, oft_cmatrix, oft_matrix_map, oft_graph
IMPLICIT NONE
#include "local.h"
private
!---------------------------------------------------------------------------------
!> Native vector class
!!
!! Used for implementing global vector operations in OFT
!---------------------------------------------------------------------------------
type, public, extends(oft_vector) :: oft_native_vector
  real(r8), pointer, contiguous, dimension(:) :: v => NULL() !< Vector values
  real(r8), pointer, contiguous, dimension(:) :: local_tmp => NULL() !< Local value cache
contains
  procedure :: new_real => vec_new_vec
  procedure :: new_complex => vec_new_cvec
  !> Set all elements to a scalar
  procedure :: set => vec_set
  !> Get local portion of field
  procedure :: get_local => vec_get_local
  !> Restore local portion of field
  procedure :: restore_local => vec_restore_local
  !> Get owned portion of field
  procedure :: get_slice => vec_get_slice
  !> Restore owned portion of field
  procedure :: restore_slice => vec_restore_slice
  !> Add vectors
  procedure :: add => vec_add
  !> Multiply fields element by element
  procedure :: mult => vec_mult
  procedure :: dot_real => vec_dot_vec
  procedure :: dot_complex => vec_dot_cvec
  procedure :: mdot_real => vec_mdot_vec
  procedure :: mdot_complex => vec_mdot_cvec
  !> Scale vector by a scalar
  procedure :: scale => vec_scale
  !> Sum reduction over vector
  procedure :: sum => vec_sum
  !> Norm of field
  procedure :: norm => vec_norm
  !> Perform global stitching
  procedure :: stitch => vec_stitch
  !> Delete vector
  procedure :: delete => vec_delete
end type oft_native_vector
!---------------------------------------------------------------------------------
!> Native complex vector class
!!
!! Used for implementing global vector operations in OFT
!---------------------------------------------------------------------------------
type, public, extends(oft_cvector) :: oft_native_cvector
  complex(c8), pointer, contiguous, dimension(:) :: v => NULL() !< Vector values
  complex(c8), pointer, contiguous, dimension(:) :: local_tmp => NULL() !< Local value cache
contains
  procedure :: new_real => cvec_new_vec
  procedure :: new_complex => cvec_new_cvec
  !> Set all elements to a scalar
  procedure :: set => cvec_set
  !> Get local portion of field
  procedure :: get_local => cvec_get_local
  !> Restore local portion of field
  procedure :: restore_local => cvec_restore_local
  !> Get owned portion of field
  procedure :: get_slice => cvec_get_slice
  !> Restore owned portion of field
  procedure :: restore_slice => cvec_restore_slice
  procedure :: add_real => cvec_add_vec
  procedure :: add_complex => cvec_add_cvec
  procedure :: mult_real => cvec_mult_vec
  procedure :: mult_complex => cvec_mult_cvec
  procedure :: dot_real => cvec_dot_vec
  procedure :: dot_complex => cvec_dot_cvec
  procedure :: mdot_real => cvec_mdot_vec
  procedure :: mdot_complex => cvec_mdot_cvec
  !> Scale vector by a scalar
  procedure :: scale => cvec_scale
  !> Sum reduction over vector
  procedure :: sum => cvec_sum
  !> Norm of field
  procedure :: norm => cvec_norm
  !> Perform global stitching
  procedure :: stitch => cvec_stitch
  !> Delete vector
  procedure :: delete => cvec_delete
end type oft_native_cvector
!---------------------------------------------------------------------------------
!> Native CRS matrix class
!---------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_native_matrix
  logical :: full_current = .FALSE. !< Is full matrix current?
  integer(i4) :: nnz = 0 !< Number of non-zero entries
  integer(i4) :: nred = 0 !< Number of redundant rows
  integer(i4), pointer, contiguous, dimension(:) :: color => NULL() !< Element coloring
  integer(i4), pointer, contiguous, dimension(:) :: redind => NULL() !< Redundant rows
  integer(i4), pointer, contiguous, dimension(:) :: kr => NULL() !< Row pointer to column list
  integer(i4), pointer, contiguous, dimension(:) :: lc => NULL() !< Column list
  integer(2), pointer, contiguous, dimension(:) :: lc_small => NULL() !< 2-byte copy of column list
  real(r8), pointer, contiguous, dimension(:) :: M => NULL() !< Matrix entries
  real(r8), pointer, contiguous, dimension(:) :: Md => NULL() !< Diagonal matrix entries
  real(r8), pointer, contiguous, dimension(:) :: Mfull => NULL() !< Full local matrix entries
  type(oft_matrix_map), pointer, dimension(:,:) :: map => NULL() !< Global context
  type(oft_seam), pointer :: linkage => NULL() !< Global linkage information
contains
  procedure :: apply_real => mat_apply_vec
  procedure :: apply_complex => mat_apply_cvec
  procedure :: applyt_real => mat_applyt_vec
  procedure :: applyt_complex => mat_applyt_cvec
  !> Set values of the matrix
  procedure :: set_values => mat_set_values
  !> Add values to the matrix
  procedure :: add_values => mat_add_values
  !> Add values atomically to the matrix
  procedure :: atomic_add_values => mat_atomic_add_values
  !> Complete matrix assembly
  procedure :: assemble => mat_assemble
  !> Zero all elements
  procedure :: zero => mat_zero
  !> Zero all elements in a given row
  procedure :: zero_rows => mat_zero_rows
  !> Zero all elements in a given column
  procedure :: zero_cols => mat_zero_cols
  !> Update full matrix with adjacent processor information
  procedure :: update_slice => matrix_update_slice
  !> Delete matrix
  procedure :: delete => mat_delete
end type oft_native_matrix
!---------------------------------------------------------------------------------
!> Native CRS complex matrix class
!---------------------------------------------------------------------------------
type, public, extends(oft_cmatrix) :: oft_native_cmatrix
  logical :: full_current = .FALSE. !< Is full matrix current?
  integer(i4) :: nnz = 0 !< Number of non-zero entries
  integer(i4) :: nred = 0 !< Number of redundant rows
  integer(i4), pointer, contiguous, dimension(:) :: color => NULL() !< Element coloring
  integer(i4), pointer, contiguous, dimension(:) :: redind => NULL() !< Redundant rows
  integer(i4), pointer, contiguous, dimension(:) :: kr => NULL() !< Row pointer to column list
  integer(i4), pointer, contiguous, dimension(:) :: lc => NULL() !< Column list
  integer(2), pointer, contiguous, dimension(:) :: lc_small => NULL() !< 2-byte copy of column list
  complex(c8), pointer, contiguous, dimension(:) :: M => NULL() !< Matrix entries
  complex(c8), pointer, contiguous, dimension(:) :: Md => NULL() !< Diagonal matrix entries
  complex(c8), pointer, contiguous, dimension(:) :: Mfull => NULL() !< Full local matrix entries
  type(oft_matrix_map), pointer, dimension(:,:) :: map => NULL() !< Global context
  type(oft_seam), pointer :: linkage => NULL() !< Global linkage information
contains
  procedure :: apply_real => cmat_apply_vec
  procedure :: apply_complex => cmat_apply_cvec
  procedure :: applyt_real => cmat_applyt_vec
  procedure :: applyt_complex => cmat_applyt_cvec
  !> Set values of the matrix
  procedure :: set_values => cmat_set_values
  !> Add values to the matrix
  procedure :: add_values => cmat_add_values
  !> Add values atomically to the matrix
  procedure :: atomic_add_values => cmat_atomic_add_values
  !> Complete matrix assembly
  procedure :: assemble => cmat_assemble
  !> Zero all elements
  procedure :: zero => cmat_zero
  !> Zero all elements in a given row
  procedure :: zero_rows => cmat_zero_rows
  !> Zero all elements in a given column
  procedure :: zero_cols => cmat_zero_cols
  !> Delete matrix
  procedure :: delete => cmat_delete
end type oft_native_cmatrix
!---------------------------------------------------------------------------------
!> Class for extracting a submatrix of @ref oft_native_matrix (eg. for block solvers)
!---------------------------------------------------------------------------------
type, public, extends(oft_native_matrix) :: oft_native_submatrix
  integer(i4) :: slice = -1 !< Slice of parent matrix
  integer(i4), pointer, dimension(:) :: lcmap => NULL() !< Column list
  integer(i4), pointer, dimension(:) :: part => NULL() !< Column list
  class(oft_native_matrix), pointer :: parent => NULL() !< Parent matrix
contains
  !> Build submatrix
  procedure :: setup => submatrix_setup
  !> Update matrix with values from parent matrix
  procedure :: update_slice => submatrix_update_slice
  !> Delete matrix
  procedure :: delete => submatrix_delete
end type oft_native_submatrix
!---------------------------------------------------------------------------------
!> Native dense matrix implementation
!---------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_native_dense_matrix
  real(r8), pointer, CONTIGUOUS, DIMENSION(:,:) :: M => NULL() !< Matrix values
contains
  procedure :: apply_real => dense_mat_apply_vec
  procedure :: apply_complex => dense_mat_apply_cvec
  procedure :: applyt_real => dense_mat_applyt_vec
  procedure :: applyt_complex => dense_mat_applyt_cvec
  !> Set values of the matrix
  procedure :: set_values => dense_mat_set_values
  !> Add values to the matrix
  procedure :: add_values => dense_mat_add_values
  !> Add values atomically to the matrix
  procedure :: atomic_add_values => dense_mat_atomic_add_values
  !> Complete matrix assembly
  procedure :: assemble => dense_mat_assemble
  !> Zero all elements
  procedure :: zero => dense_mat_zero
  !> Zero all elements in a given row
  procedure :: zero_rows => dense_mat_zero_rows
end type oft_native_dense_matrix
!---------------------------------------------------------------------------------
!> Native dense complex matrix implementation
!---------------------------------------------------------------------------------
type, public, extends(oft_cmatrix) :: oft_native_dense_cmatrix
  complex(c8), pointer, CONTIGUOUS, DIMENSION(:,:) :: M => NULL() !< Matrix values
contains
  procedure :: apply_real => dense_cmat_apply_vec
  procedure :: apply_complex => dense_cmat_apply_cvec
  procedure :: applyt_real => dense_cmat_applyt_vec
  procedure :: applyt_complex => dense_cmat_applyt_cvec
  !> Set values of the matrix
  procedure :: set_values => dense_cmat_set_values
  !> Add values to the matrix
  procedure :: add_values => dense_cmat_add_values
  !> Add values atomically to the matrix
  procedure :: atomic_add_values => dense_cmat_atomic_add_values
  !> Complete matrix assembly
  procedure :: assemble => dense_cmat_assemble
  !> Zero all elements
  procedure :: zero => dense_cmat_zero
  !> Zero all elements in a given row
  procedure :: zero_rows => dense_cmat_zero_rows
end type oft_native_dense_cmatrix
INTERFACE
!------------------------------------------------------------------------------
!> Parition CRS graph using METIS
!------------------------------------------------------------------------------
  SUBROUTINE oft_metis_partgraph(nr,nnz,kr,lc,npart,part,type,info) BIND(C,NAME="oft_metis_partGraph")
  IMPORT c_int
  INTEGER(c_int), INTENT(in) :: nr !< Number of rows
  INTEGER(c_int), INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), DIMENSION(nr+1), INTENT(inout) :: kr !< Row pointer into lc [nr+1]
  INTEGER(c_int), DIMENSION(nnz), INTENT(inout) :: lc !< Column indices in CRS format [nnz]
  INTEGER(c_int), INTENT(in) :: npart !< Number of partitions
  INTEGER(c_int), DIMENSION(npart), INTENT(inout) :: part !< Partition for each row [n]
  INTEGER(c_int), INTENT(in) :: type !< Type of partitioning to use (1->ML_Recursive,2->ML_KWay)
  INTEGER(c_int), INTENT(inout) :: info !< Parition return status
  END SUBROUTINE oft_metis_partgraph
END INTERFACE
INTEGER(i4), PARAMETER :: lc_offset = HUGE(INT(1,2))
!---Declare public entities
public native_vector_cast, native_cvector_cast, native_matrix_cast, native_cmatrix_cast
public partition_graph, native_matrix_setup_full
public native_vector_slice_push, native_vector_slice_pop
contains
!---------------------------------------------------------------------------------
!> Cast an abstract vector object to @ref oft_native_vector
!!
!! The source vector must be @ref oft_native_vector or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!---------------------------------------------------------------------------------
FUNCTION native_vector_cast(self,source) result(success)
class(oft_native_vector), pointer, intent(out) :: self !< Cast pointer
class(oft_vector), target, intent(in) :: source !< Abstract vector object
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_vector)
    self=>source
    success=.TRUE.
  class default
    NULLIFY(self)
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION native_vector_cast
!---------------------------------------------------------------------------------
!> Create a new vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine vec_new_vec(self,new)
class(oft_native_vector), intent(in) :: self !< Vector object
class(oft_vector), pointer, intent(out) :: new !< New vector
DEBUG_STACK_PUSH
allocate(oft_native_vector::new)
select type(this=>new)
  class is(oft_native_vector)
    this%n=self%n
    this%nslice=self%nslice
    this%ng=self%ng
    this%nblocks=self%nblocks
    this%map=>self%map
    this%stitch_info=>self%stitch_info
    allocate(this%v(this%n))
    call this%set(0.d0)
  class default
    CALL oft_abort('Failure to allocate vector.','vec_new_vec',__FILE__)
end select
DEBUG_STACK_POP
end subroutine vec_new_vec
!---------------------------------------------------------------------------------
!> Create a new complex vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine vec_new_cvec(self,new)
class(oft_native_vector), intent(in) :: self !< Vector object
class(oft_cvector), pointer, intent(out) :: new !< New vector
DEBUG_STACK_PUSH
allocate(oft_native_cvector::new)
select type(this=>new)
  class is(oft_native_cvector)
    this%n=self%n
    this%nslice=self%nslice
    this%ng=self%ng
    this%nblocks=self%nblocks
    this%map=>self%map
    this%stitch_info=>self%stitch_info
    allocate(this%v(this%n))
    call this%set((0.d0,0.d0))
  class default
    CALL oft_abort('Failure to allocate vector.','vec_new_cvec',__FILE__)
end select
DEBUG_STACK_POP
end subroutine vec_new_cvec
!---------------------------------------------------------------------------------
!> Set all elements to a scalar
!---------------------------------------------------------------------------------
subroutine vec_set(self,alpha,iblock,random)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: alpha !< Updated vector value
integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
logical :: random_flag
integer(i4) :: i
DEBUG_STACK_PUSH
random_flag=.FALSE.
IF(PRESENT(random))random_flag=random
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_set',__FILE__)
  !---
  IF(random_flag)THEN
    CALL oft_random_number(self%v(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n),self%map(iblock)%n)
  ELSE
    do i=1,self%map(iblock)%n
      self%v(self%map(iblock)%offset+i) = alpha
    end do
  END IF
  CALL self%stitch(0)
ELSE
  IF(random_flag)THEN
    CALL oft_random_number(self%v,self%n)
    CALL self%stitch(0)
  ELSE
    !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
    do i=1,self%n
      self%v(i) = alpha
    end do
  END IF
END IF
DEBUG_STACK_POP
end subroutine vec_set
!---------------------------------------------------------------------------------
!> Get values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine vec_get_slice(self,array,iblock)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), pointer, intent(inout) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,j,nslice
DEBUG_STACK_PUSH
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_get_slice',__FILE__)
  !---
  nslice=self%map(iblock)%nslice
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(nslice))
  !$omp simd
  do j=1,self%map(iblock)%nslice
    array(j) = self%v(self%map(iblock)%offset+self%map(iblock)%slice(j))
  end do
ELSE
  !---
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%nslice))
  IF(self%n==self%nslice)THEN
    array=self%v
  ELSE
    !$omp parallel do private(j) if(self%n>OFT_OMP_VTHRESH)
    do i=1,self%nblocks
      !$omp simd
      do j=1,self%map(i)%nslice
        array(self%map(i)%soffset+j) = self%v(self%map(i)%offset+self%map(i)%slice(j))
      end do
    end do
  END IF
END IF
DEBUG_STACK_POP
end subroutine vec_get_slice
!---------------------------------------------------------------------------------
!> Set/add values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine vec_restore_slice(self,array,iblock,wait)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: wait !< Wait to perform global sync?
integer(i4) :: i,j
logical :: do_wait
DEBUG_STACK_PUSH
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_slice',__FILE__)
  !$omp simd
  do j=1,self%map(iblock)%nslice
    self%v(self%map(iblock)%offset+self%map(iblock)%slice(j)) = array(j)
  end do
ELSE
  IF(self%n==self%nslice)THEN
    self%v=array
  ELSE
    !$omp parallel do private(j) if(self%n>OFT_OMP_VTHRESH)
    do i=1,self%nblocks
      !$omp simd
      do j=1,self%map(i)%nslice
        self%v(self%map(i)%offset+self%map(i)%slice(j)) = array(self%map(i)%soffset+j)
      end do
    end do
  END IF
END IF
IF(.NOT.do_wait)call self%stitch(0)
DEBUG_STACK_POP
end subroutine vec_restore_slice
!---------------------------------------------------------------------------------
!> Get local values from vector
!---------------------------------------------------------------------------------
subroutine vec_get_local(self,array,iblock)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), pointer, intent(inout) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,j
DEBUG_STACK_PUSH
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_get_local',__FILE__)
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%n))
  !$omp do simd
  do j=1,self%map(iblock)%n
    array(j) = self%v(self%map(iblock)%offset+j)
  end do
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%n))
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  DO i=1,self%n
    array(i)=self%v(i)
  END DO
  !$omp end parallel
END IF
DEBUG_STACK_POP
end subroutine vec_get_local
!---------------------------------------------------------------------------------
!> Set/add local values to vector
!---------------------------------------------------------------------------------
subroutine vec_restore_local(self,array,iblock,add,wait)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: add !< Add values instead of replace
logical, optional, intent(in) :: wait !< Wait to perform global sync
logical :: do_add,do_wait
integer(i4) :: j
DEBUG_STACK_PUSH
do_add=.FALSE.
IF(PRESENT(add))do_add=add
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_local',__FILE__)
  IF(do_add)THEN
    IF(.NOT.ASSOCIATED(self%local_tmp))THEN
      ALLOCATE(self%local_tmp(self%n))
      self%local_tmp=0.d0
    END IF
    !$omp simd
    do j=1,self%map(iblock)%n
      self%local_tmp(self%map(iblock)%offset+j) = array(j)
    end do
    IF(.NOT.do_wait)THEN
      CALL full_restore_local(self%local_tmp)
      DEALLOCATE(self%local_tmp)
    END IF
  ELSE
    !$omp simd
    do j=1,self%map(iblock)%n
      self%v(self%map(iblock)%offset+j) = array(j)
    end do
    IF(.NOT.do_wait)CALL self%stitch(0)
  END IF
ELSE
  IF(ASSOCIATED(self%local_tmp))CALL oft_abort("Existing partial restore found when restoring full vector", &
    "vec_restore_local",__FILE__)
  IF(do_wait)CALL oft_abort("Wait not supported when restoring full vector","vec_restore_local",__FILE__)
  ALLOCATE(self%local_tmp(self%n))
  self%local_tmp=array
  CALL full_restore_local(self%local_tmp)
  DEALLOCATE(self%local_tmp)
END IF
DEBUG_STACK_POP
contains
!
subroutine full_restore_local(array_loc)
real(r8), intent(inout) :: array_loc(:)
integer(i4) :: i
IF(do_add)THEN
  call oft_global_stitch(self%stitch_info,array_loc,1)
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  DO i=1,self%n
    self%v(i)=self%v(i)+array_loc(i)
  END DO
  !$omp end parallel
ELSE
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  DO i=1,self%n
    self%v(i)=array_loc(i)
  END DO
  !$omp end parallel
  IF(.NOT.do_wait)CALL self%stitch(0)
END IF
end subroutine full_restore_local
end subroutine vec_restore_local
!---------------------------------------------------------------------------------
!> Add vectors
!!
!! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
!---------------------------------------------------------------------------------
subroutine vec_add(self,gamma,alpha,a,beta,b)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: gamma !< Scale of source vector
real(r8), intent(in) :: alpha !< Scale of first vector
class(oft_vector), target, intent(inout) :: a !< First vector to add
real(r8), optional, intent(in) :: beta !< Scale of second vector (optional)
class(oft_vector), target, optional, intent(inout) :: b !< Second vector to add (optional)
logical :: dealloc_flag(2)
integer(i4) :: i
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
dealloc_flag=.FALSE.
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','vec_add_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag(1)=.TRUE.
END SELECT
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','vec_add_vec',__FILE__)
  SELECT TYPE(b)
  CLASS IS(oft_native_vector)
    btmp=>b%v
  CLASS DEFAULT
    NULLIFY(btmp)
    CALL b%get_local(btmp)
    dealloc_flag(2)=.TRUE.
  END SELECT
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i) + beta*btmp(i)
  end do
  !$omp end parallel
else ! Add one vector to the source
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i)
  end do
  !$omp end parallel
end if
IF(dealloc_flag(1))DEALLOCATE(atmp)
IF(dealloc_flag(2))DEALLOCATE(btmp)
DEBUG_STACK_POP
end subroutine vec_add
!---------------------------------------------------------------------------------
!> Elementwise multiplication with another vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine vec_mult(self,a,div_flag)
class(oft_native_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
logical :: div,dealloc_flag
integer(i4) :: i
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
!---
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','vec_mult_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flag=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag=.TRUE.
END SELECT
!---
if(div)then
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,self%n
    self%v(i) = self%v(i)/atmp(i)
  end do
  !$omp end parallel
else
  !$omp parallel if(self%n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,self%n
    self%v(i) = self%v(i)*atmp(i)
  end do
  !$omp end parallel
end if
IF(dealloc_flag)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine vec_mult
!---------------------------------------------------------------------------------
!> Scale vector by a scalar
!---------------------------------------------------------------------------------
subroutine vec_scale(self,alpha)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: alpha !< Scale factor
integer(i4) :: i
DEBUG_STACK_PUSH
!$omp parallel if(self%n>OFT_OMP_VTHRESH)
!$omp do simd
do i=1,self%n
  self%v(i) = alpha*self%v(i)
end do
!$omp end parallel
DEBUG_STACK_POP
end subroutine vec_scale
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function vec_dot_vec(self,a) result(dot)
class(oft_native_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
real(r8) :: dot !< \f$ \sum_i self_i a_i \f$
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_dot_vec',__FILE__)
    dot=oft_global_dp(self%stitch_info,self%v,this%v,this%n)
  class default
    CALL oft_abort('"a" is not a native vector.','vec_dot_vec',__FILE__)
end select
DEBUG_STACK_POP
end function vec_dot_vec
!---------------------------------------------------------------------------------
!> Dot product with a complex vector
!---------------------------------------------------------------------------------
function vec_dot_cvec(self,a) result(dot)
class(oft_native_vector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
complex(c8) :: dot !< \f$ \sum_i self_i a_i \f$
complex(c8), allocatable, dimension(:) :: vtmp
DEBUG_STACK_PUSH
allocate(vtmp(self%n))
vtmp=(1.d0,0.d0)*self%v
select type(this=>a)
  class is(oft_native_cvector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_dot_cvec',__FILE__)
    dot=oft_global_dp(self%stitch_info,vtmp,this%v,this%n)
  class default
    CALL oft_abort('"a" is not a native vector.','vec_dot_cvec',__FILE__)
end select
deallocate(vtmp)
DEBUG_STACK_POP
end function vec_dot_cvec
!---------------------------------------------------------------------------------
!> Dot product with an array of vectors
!---------------------------------------------------------------------------------
function vec_mdot_vec(self,a,n) result(dots)
class(oft_native_vector), intent(inout) :: self !< Vector object
type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
real(r8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
integer(i4) :: i
DEBUG_STACK_PUSH
DO i=1,n
  select type(this=>a(i)%f)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_mdot_vec',__FILE__)
    dots(i)=oft_global_dp(self%stitch_info,self%v,this%v,this%n,no_reduce=.TRUE.)
  class default
    CALL oft_abort('"a" is not a native vector.','vec_mdot_vec',__FILE__)
  end select
END DO
IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
DEBUG_STACK_POP
end function vec_mdot_vec
!---------------------------------------------------------------------------------
!> Dot product with an array of complex vectors
!---------------------------------------------------------------------------------
function vec_mdot_cvec(self,a,n) result(dots)
class(oft_native_vector), intent(inout) :: self !< Vector object
type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
integer(i4) :: i
complex(c8), allocatable, dimension(:) :: vtmp
DEBUG_STACK_PUSH
allocate(vtmp(self%n))
vtmp=(1.d0,0.d0)*self%v
DO i=1,n
  select type(this=>a(i)%f)
  class is(oft_native_cvector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_mdot_cvec',__FILE__)
    dots(i)=oft_global_dp(self%stitch_info,vtmp,this%v,this%n,no_reduce=.TRUE.)
  class default
    CALL oft_abort('"a" is not a native complex vector.','vec_mdot_cvec',__FILE__)
  end select
END DO
deallocate(vtmp)
IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
DEBUG_STACK_POP
end function vec_mdot_cvec
!---------------------------------------------------------------------------------
!> Sum reduction over vector
!---------------------------------------------------------------------------------
function vec_sum(self) result(sum)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8) :: sum !< Sum of vector elements
DEBUG_STACK_PUSH
sum=oft_global_reduction(self%stitch_info,self%v,self%n)
DEBUG_STACK_POP
end function vec_sum
!---------------------------------------------------------------------------------
!> Compute norm of vector
!---------------------------------------------------------------------------------
function vec_norm(self,itype) result(norm)
class(oft_native_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
real(r8) :: norm !< Specified norm of vector
integer(i4) :: i
real(r8) :: nloc
real(r8), allocatable, dimension(:) :: vtmp
DEBUG_STACK_PUSH
SELECT CASE(itype)
  CASE(1)
    ALLOCATE(vtmp(self%n))
    vtmp=ABS(self%v)
    norm=oft_global_reduction(self%stitch_info,vtmp,self%n)
    DEALLOCATE(vtmp)
  CASE(2)
    norm=SQRT(self%dot(self))
  CASE(3)
    nloc=0.d0
    DO i=1,self%n
      nloc=MAX(nloc,ABS(self%v(i)))
    END DO
    norm=oft_mpi_max(nloc)
  CASE DEFAULT
    CALL oft_abort("Invalid norm type","vec_norm",__FILE__)
END SELECT
DEBUG_STACK_POP
end function vec_norm
!---------------------------------------------------------------------------------
!> Perform global stitching
!---------------------------------------------------------------------------------
subroutine vec_stitch(self,up_method)
class(oft_native_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: up_method !< Type of stitching to perform
DEBUG_STACK_PUSH
call oft_global_stitch(self%stitch_info,self%v,up_method)
DEBUG_STACK_POP
end subroutine vec_stitch
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine vec_delete(self)
class(oft_native_vector), intent(inout) :: self !< Vector object
DEBUG_STACK_PUSH
if(associated(self%v))deallocate(self%v)
if(associated(self%local_tmp))deallocate(self%local_tmp)
self%n=-1
self%nslice=-1
self%ng=-1
self%nblocks=1
NULLIFY(self%map,self%stitch_info)
DEBUG_STACK_POP
end subroutine vec_delete
!------------------------------------------------------------------------------
!> Insert vector data into a restart structure for output
!!
!! Vector data and associated offsets are copied into the restart structure
!! for use with \ref oft_io::hdf5_rst_write "hdf5_rst_write".
!------------------------------------------------------------------------------
subroutine native_vector_slice_push(self,ig,rst_info,alloc_only)
class(oft_native_vector), pointer, intent(in) :: self !< Source vector
integer(i8), intent(in) :: ig(:) !< Global indices for `self`
type(hdf5_rst), intent(out) :: rst_info !< Restart structure for output
logical, optional, intent(in) :: alloc_only !< Allocate and setup data structures only (optional)
integer(i4) :: nemax,ierr,i,j,k,l,n
integer(i8) :: ii
integer(i8), allocatable, dimension(:) :: t1i,t2i,ind
real(r8) :: fcount
real(r8), allocatable, dimension(:) :: t1v,t2v
#ifdef OFT_MPI_F08
type(mpi_request) :: requests(4)
#else
integer(i4) :: requests(4)
#endif
DEBUG_STACK_PUSH
if(self%stitch_info%full)then
  rst_info%offset=0
  rst_info%count=INT(self%ng,4)
  rst_info%dim=self%ng
  rst_info%full=.TRUE.
  allocate(rst_info%data(rst_info%count))
  rst_info%data=0.d0
  IF(PRESENT(alloc_only))THEN
    IF(alloc_only)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
  do i=1,self%n
    rst_info%data(ABS(ig(i)))=self%v(i)
  end do
else
#ifdef HAVE_MPI
  fcount=REAL(self%ng,8)/oft_env%nprocs
  rst_info%dim=self%ng
  rst_info%offset=CEILING(fcount*oft_env%rank)
  rst_info%count=INT(CEILING(fcount*(oft_env%rank+1))-rst_info%offset,4)
  IF(oft_env%rank==oft_env%nprocs-1)rst_info%count=INT(self%ng-rst_info%offset,4)
  rst_info%full=.FALSE.
  allocate(rst_info%data(rst_info%count))
  rst_info%data=0.d0
  IF(PRESENT(alloc_only))THEN
    IF(alloc_only)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
  !---Create temporary field storage
  nemax=oft_mpi_max(self%nslice)
  allocate(t1i(nemax+1),t1v(nemax))
  allocate(t2i(nemax+1),t2v(nemax))
  t1i=0; t1v=0.d0
  t2i=0; t2v=0.d0
  !---Extract values owned by current processor
  j=0
  DO i=1,self%n
    IF(.NOT.self%stitch_info%be(i))THEN
      ii=ABS(ig(i))
      j=j+1
      t1i(j) = ii
      t1v(j) = self%v(i)
    END IF
  END DO
  DO i=1,self%stitch_info%nbe
    IF(self%stitch_info%leo(i))THEN
      l=self%stitch_info%lbe(i)
      ii=ABS(ig(l))
      j=j+1
      t1i(j) = ii
      t1v(j) = self%v(l)
    END IF
  END DO
  !---Order entries by global index
  n=j
  ALLOCATE(ind(n))
  ind=[(i,i=1,n)]
  t2i=t1i
  t2v=t1v
  CALL sort_array(t1i,ind,INT(n,8))
  t1v=0.d0
  DO i=1,n
    t1v(i)=t2v(ind(i))
  END DO
  t1i(nemax+1)=n ! Set extent
  DEALLOCATE(ind)
  !---Save elements to output array (from current processor)
  do i=1,n
    IF(t1i(i)<=rst_info%offset)CYCLE
    IF(t1i(i)>rst_info%offset+rst_info%count)EXIT
    rst_info%data(t1i(i)-rst_info%offset)=t1v(i)
  enddo
  !---Set adjacent processors
  j=oft_env%rank + 1
  if(j > oft_env%nprocs-1)j = 0
  k=oft_env%rank - 1
  if(k < 0)k = oft_env%nprocs-1
  !---Setup initial sends
  CALL MPI_ISEND(t1i, nemax+1, OFT_MPI_I8, j, 1, oft_env%comm, requests(1), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
  CALL MPI_ISEND(t1v, nemax, OFT_MPI_R8, j, 2, oft_env%comm, requests(2), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
  !---Setup initial recvs
  CALL MPI_IRECV(t2i, nemax+1, OFT_MPI_I8, k, 1, oft_env%comm, requests(3), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
  CALL MPI_IRECV(t2v, nemax, OFT_MPI_R8, k, 2, oft_env%comm, requests(4), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
  !---Ring communication
  do l=1,oft_env%nprocs-1
    !---Check for transfer completion
    CALL oft_mpi_waitall(4, requests, ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','vector_slice_push',__FILE__)
    t1i=t2i
    t1v=t2v
    !---Setup next transfer
    IF(l<oft_env%nprocs-1)THEN
      !---Sends
      CALL MPI_ISEND(t1i, nemax+1, OFT_MPI_I8, j, 1, oft_env%comm, requests(1), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
      CALL MPI_ISEND(t1v, nemax, OFT_MPI_R8, j, 2, oft_env%comm, requests(2), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
      !---Recvs
      CALL MPI_IRECV(t2i, nemax+1, OFT_MPI_I8, k, 1, oft_env%comm, requests(3), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
      CALL MPI_IRECV(t2v, nemax, OFT_MPI_R8, k, 2, oft_env%comm, requests(4), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
    END IF
    !---Save elements to output array (from next processor)
    do i=1,INT(t1i(nemax+1),4)
      IF(t1i(i)<=rst_info%offset)CYCLE
      IF(t1i(i)>rst_info%offset+rst_info%count)EXIT
      rst_info%data(t1i(i)-rst_info%offset)=t1v(i)
    enddo
  enddo
  deallocate(t1i,t1v,t2i,t2v)
#else
CALL oft_abort("Distributed vectors require MPI","native_vector_slice_push",__FILE__)
#endif
end if
DEBUG_STACK_POP
end subroutine native_vector_slice_push
!------------------------------------------------------------------------------
!> Insert data from a restart structure into a vector
!!
!! Vector data is copied from the restart structure following a call to
!! \ref oft_io::hdf5_rst_read "hdf5_rst_read".
!!
!! @note The restart structure should be setup first using a call to
!! \ref oft_native_vectors::vector_slice_push "vector_slice_push"
!------------------------------------------------------------------------------
subroutine native_vector_slice_pop(self,ig,rst_info)
class(oft_native_vector), pointer, intent(inout) :: self !< Destination vector
integer(i8), intent(in) :: ig(:) !< Global indices for `self`
type(hdf5_rst), intent(inout) :: rst_info !< Restart structure for output
integer(i4) :: nemax,ierr,i,j,k,l
integer(i8) :: ii,t1i(2),t2i(2)
real(r8), allocatable, dimension(:) :: t1v,t2v
#ifdef OFT_MPI_F08
type(mpi_request) :: requests(4)
#else
integer(i4) :: requests(4)
#endif
DEBUG_STACK_PUSH
if(self%stitch_info%full)then
  do i=1,self%n
    self%v(i)=rst_info%data(ABS(ig(i)))
  end do
else
#ifdef HAVE_MPI
  !---Create temporary field storage
  nemax=oft_mpi_max(rst_info%count)
  ALLOCATE(t1v(nemax),t2v(nemax))
  t1v=0.d0; t2v=0.d0
  t1i=(/rst_info%offset,rst_info%offset+rst_info%count/)
  do i=1,rst_info%count
    t1v(i)=rst_info%data(i)
  end do
  !---Fetch elements from local chunk
  self%v=0.d0
  DO i=1,self%n
    ii=ABS(ig(i))
    IF(ii>t1i(1).AND.ii<=t1i(2))self%v(i)=t1v(ii-t1i(1))
  END DO
  !---Set adjacent processors
  j=oft_env%rank + 1
  if(j > oft_env%nprocs-1)j = 0
  k=oft_env%rank - 1
  if(k < 0)k = oft_env%nprocs-1
  !---Setup initial sends
  CALL MPI_ISEND(t1i, 2, OFT_MPI_I8, j, 1, oft_env%comm, requests(1), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
  CALL MPI_ISEND(t1v, nemax, OFT_MPI_R8, j, 2, oft_env%comm, requests(2), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
  !---Setup initial recvs
  CALL MPI_IRECV(t2i, 2, OFT_MPI_I8, k, 1, oft_env%comm, requests(3), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
  CALL MPI_IRECV(t2v, nemax, OFT_MPI_R8, k, 2, oft_env%comm, requests(4), ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
  !---Ring communication
  do l=1,oft_env%nprocs-1
    !---Check for transfer completion
    CALL oft_mpi_waitall(4, requests, ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','vector_slice_pop',__FILE__)
    t1i=t2i
    t1v=t2v
    !---Setup next transfer
    IF(l<oft_env%nprocs-1)THEN
      !---Sends
      CALL MPI_ISEND(t1i, 2, OFT_MPI_I8, j, 1, oft_env%comm, requests(1), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
      CALL MPI_ISEND(t1v, nemax, OFT_MPI_R8, j, 2, oft_env%comm, requests(2), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','vector_slice_push',__FILE__)
      !---Recvs
      CALL MPI_IRECV(t2i, 2, OFT_MPI_I8, k, 1, oft_env%comm, requests(3), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
      CALL MPI_IRECV(t2v, nemax, OFT_MPI_R8, k, 2, oft_env%comm, requests(4), ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','vector_slice_push',__FILE__)
    END IF
    !---Fetch elements from next chunk
    do i=1,self%n
      ii=ABS(ig(i))
      IF(ii>t1i(1).AND.ii<=t1i(2))self%v(i)=t1v(ii-t1i(1))
    end do
  enddo
  DEALLOCATE(t1v,t2v)
#else
CALL oft_abort("Distributed vectors require MPI","native_vector_slice_pop",__FILE__)
#endif
end if
DEBUG_STACK_POP
end subroutine native_vector_slice_pop
!---------------------------------------------------------------------------------
!> Cast an abstract vector object to @ref oft_native_cvector
!!
!! The source vector must be @ref oft_native_vector or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!---------------------------------------------------------------------------------
FUNCTION native_cvector_cast(self,source) result(success)
class(oft_native_cvector), pointer, intent(out) :: self !< Cast pointer
class(oft_cvector), target, intent(in) :: source !< Abstract vector object
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cvector)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION native_cvector_cast
!---------------------------------------------------------------------------------
!> Create a new vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine cvec_new_vec(self,new)
class(oft_native_cvector), intent(in) :: self !< Vector object
class(oft_vector), pointer, intent(out) :: new !< New vector
DEBUG_STACK_PUSH
allocate(oft_native_vector::new)
select type(this=>new)
  class is(oft_native_vector)
    this%n=self%n
    this%nslice=self%nslice
    this%ng=self%ng
    this%nblocks=self%nblocks
    this%map=>self%map
    this%stitch_info=>self%stitch_info
    allocate(this%v(this%n))
    call this%set(0.d0)
  class default
    CALL oft_abort('Failure to allocate native real vector.','cvec_new_vec',__FILE__)
end select
DEBUG_STACK_POP
end subroutine cvec_new_vec
!---------------------------------------------------------------------------------
!> Create a new vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine cvec_new_cvec(self,new)
class(oft_native_cvector), intent(in) :: self !< Vector object
class(oft_cvector), pointer, intent(out) :: new !< New vector
DEBUG_STACK_PUSH
allocate(oft_native_cvector::new)
select type(this=>new)
  class is(oft_native_cvector)
    this%n=self%n
    this%nslice=self%nslice
    this%ng=self%ng
    this%nblocks=self%nblocks
    this%map=>self%map
    this%stitch_info=>self%stitch_info
    allocate(this%v(this%n))
    call this%set((0.d0,0.d0))
  class default
    CALL oft_abort('Failure to allocate vector.','vec_new_cvec',__FILE__)
end select
DEBUG_STACK_POP
end subroutine cvec_new_cvec
!---------------------------------------------------------------------------------
!> Set all elements to a scalar
!---------------------------------------------------------------------------------
subroutine cvec_set(self,alpha,iblock,random)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: alpha !< Updated vector value
integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
logical :: random_flag
integer(i4) :: i
DEBUG_STACK_PUSH
random_flag=.FALSE.
IF(PRESENT(random))random_flag=random
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_set',__FILE__)
  !---
  IF(random_flag)THEN
    CALL oft_random_number(self%v(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n),self%map(iblock)%n)
  ELSE
    do i=1,self%map(iblock)%n
      self%v(self%map(iblock)%offset+i) = alpha
    end do
  END IF
  CALL self%stitch(0)
ELSE
  IF(random_flag)THEN
    CALL oft_random_number(self%v,self%n)
    CALL self%stitch(0)
  ELSE
    !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
    do i=1,self%n
      self%v(i) = alpha
    end do
  END IF
END IF
DEBUG_STACK_POP
end subroutine cvec_set
!---------------------------------------------------------------------------------
!> Get values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine cvec_get_slice(self,array,iblock)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), pointer, intent(inout) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,j,nslice
DEBUG_STACK_PUSH
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_get_slice',__FILE__)
  !---
  nslice=self%map(iblock)%nslice
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(nslice))
  do j=1,self%map(iblock)%nslice
    array(j) = self%v(self%map(iblock)%offset+self%map(iblock)%slice(j))
  end do
ELSE
  !---
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%nslice))
  IF(self%n==self%nslice)THEN
    array=self%v
  ELSE
    do i=1,self%nblocks
      do j=1,self%map(i)%nslice
        array(self%map(i)%soffset+j) = self%v(self%map(i)%offset+self%map(i)%slice(j))
      end do
    end do
  END IF
END IF
DEBUG_STACK_POP
end subroutine cvec_get_slice
!---------------------------------------------------------------------------------
!> Set/add values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine cvec_restore_slice(self,array,iblock,wait)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: wait !< Wait to perform global sync?
integer(i4) :: i,j
logical :: do_wait
DEBUG_STACK_PUSH
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_restore_slice',__FILE__)
  do j=1,self%map(iblock)%nslice
    self%v(self%map(iblock)%offset+self%map(iblock)%slice(j)) = array(j)
  end do
ELSE
  IF(self%n==self%nslice)THEN
    self%v=array
  ELSE
    do i=1,self%nblocks
      do j=1,self%map(i)%nslice
        self%v(self%map(i)%offset+self%map(i)%slice(j)) = array(self%map(i)%soffset+j)
      end do
    end do
  END IF
END IF
IF(.NOT.do_wait)call self%stitch(0)
DEBUG_STACK_POP
end subroutine cvec_restore_slice
!---------------------------------------------------------------------------------
!> Get local values from vector
!---------------------------------------------------------------------------------
subroutine cvec_get_local(self,array,iblock)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), pointer, intent(inout) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_get_local',__FILE__)
  !---
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%n))
  do j=1,self%map(iblock)%n
    array(j) = self%v(self%map(iblock)%offset+j)
  end do
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%n))
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  DO i=1,self%n
    array(i)=self%v(i)
  END DO
END IF
DEBUG_STACK_POP
end subroutine cvec_get_local
!---------------------------------------------------------------------------------
!> Set/add local values to vector
!---------------------------------------------------------------------------------
subroutine cvec_restore_local(self,array,iblock,add,wait)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: add !< Add values instead of replace
logical, optional, intent(in) :: wait !< Wait to perform global sync?
integer(i4) :: j
logical :: do_add,do_wait
DEBUG_STACK_PUSH
do_add=.FALSE.
IF(PRESENT(add))do_add=add
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_restore_local',__FILE__)
  !---
  IF(do_add)THEN
    IF(.NOT.ASSOCIATED(self%local_tmp))THEN
      ALLOCATE(self%local_tmp(self%n))
      self%local_tmp=(0.d0,0.d0)
    END IF
    !$omp simd
    do j=1,self%map(iblock)%n
      self%local_tmp(self%map(iblock)%offset+j) = array(j)
    end do
    IF(.NOT.do_wait)THEN
      CALL full_restore_local(self%local_tmp)
      DEALLOCATE(self%local_tmp)
    END IF
  ELSE
    !$omp simd
    do j=1,self%map(iblock)%n
      self%v(self%map(iblock)%offset+j) = array(j)
    end do
    IF(.NOT.do_wait)CALL self%stitch(0)
  END IF
ELSE
  IF(ASSOCIATED(self%local_tmp))CALL oft_abort("Existing partial restore found when restoring full vector", &
    "cvec_restore_local",__FILE__)
  IF(do_wait)CALL oft_abort("Wait not supported when restoring full vector", &
    "cvec_restore_local",__FILE__)
  ALLOCATE(self%local_tmp(self%n))
  self%local_tmp=array
  CALL full_restore_local(self%local_tmp)
  DEALLOCATE(self%local_tmp)
END IF
DEBUG_STACK_POP
contains
!
subroutine full_restore_local(array_loc)
complex(c8), intent(inout) :: array_loc(:)
integer(i4) :: i
IF(do_add)THEN
  call oft_global_stitch(self%stitch_info,array_loc,1)
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  DO i=1,self%n
    self%v(i)=self%v(i)+array_loc(i)
  END DO
ELSE
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  DO i=1,self%n
    self%v(i)=array_loc(i)
  END DO
  IF(.NOT.do_wait)CALL self%stitch(0)
END IF
end subroutine full_restore_local
end subroutine cvec_restore_local
!---------------------------------------------------------------------------------
!> Add vectors
!!
!! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
!---------------------------------------------------------------------------------
subroutine cvec_add_vec(self,gamma,alpha,a,beta,b)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: gamma !< Scale of source vector
complex(c8), intent(in) :: alpha !< Scale of first vector
class(oft_vector), target, intent(inout) :: a !< First vector to add
complex(c8), optional, intent(in) :: beta !< Scale of second vector (optional)
class(oft_vector), target, optional, intent(inout) :: b !< Second vector to add (optional)
logical :: dealloc_flag(2)
integer(i4) :: i
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
dealloc_flag=.FALSE.
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','cvec_add_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag(1)=.TRUE.
END SELECT
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','cvec_add_vec',__FILE__)
  SELECT TYPE(b)
  CLASS IS(oft_native_vector)
    btmp=>b%v
  CLASS DEFAULT
    NULLIFY(btmp)
    CALL b%get_local(btmp)
    dealloc_flag(2)=.TRUE.
  END SELECT
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i) + beta*btmp(i)
  end do
else ! Add one vector to the source
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i)
  end do
end if
IF(dealloc_flag(1))DEALLOCATE(atmp)
IF(dealloc_flag(2))DEALLOCATE(btmp)
DEBUG_STACK_POP
end subroutine cvec_add_vec
!---------------------------------------------------------------------------------
!> Add vectors
!!
!! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
!---------------------------------------------------------------------------------
subroutine cvec_add_cvec(self,gamma,alpha,a,beta,b)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: gamma !< Scale of source vector
complex(c8), intent(in) :: alpha !< Scale of first vector
class(oft_cvector), target, intent(inout) :: a !< First vector to add
complex(c8), optional, intent(in) :: beta !< Scale of second vector (optional)
class(oft_cvector), target, optional, intent(inout) :: b !< Second vector to add (optional)
logical :: dealloc_flag(2)
integer(i4) :: i
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
dealloc_flag=.FALSE.
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','cvec_add_cvec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag(1)=.TRUE.
END SELECT
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','cvec_add_cvec',__FILE__)
  SELECT TYPE(b)
  CLASS IS(oft_native_cvector)
    btmp=>b%v
  CLASS DEFAULT
    NULLIFY(btmp)
    CALL b%get_local(btmp)
    dealloc_flag(2)=.TRUE.
  END SELECT
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i) + beta*btmp(i)
  end do
else ! Add one vector to the source
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = gamma*self%v(i) + alpha*atmp(i)
  end do
end if
IF(dealloc_flag(1))DEALLOCATE(atmp)
IF(dealloc_flag(2))DEALLOCATE(btmp)
DEBUG_STACK_POP
end subroutine cvec_add_cvec
!---------------------------------------------------------------------------------
!> Elementwise multiplication with another vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine cvec_mult_vec(self,a,div_flag)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
logical :: div,dealloc_flag
integer(i4) :: i
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
!---
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','cvec_mult_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flag=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag=.TRUE.
END SELECT
!---
if(div)then
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = self%v(i)/atmp(i)
  end do
else
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = self%v(i)*atmp(i)
  end do
end if
IF(dealloc_flag)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_mult_vec
!---------------------------------------------------------------------------------
!> Elementwise multiplication with another vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine cvec_mult_cvec(self,a,div_flag)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
logical :: div,dealloc_flag
integer(i4) :: i
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
!---
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','cvec_mult_cvec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
  dealloc_flag=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flag=.TRUE.
END SELECT
!---
if(div)then
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = self%v(i)/atmp(i)
  end do
else
  !$omp parallel do if(self%n>OFT_OMP_VTHRESH)
  do i=1,self%n
    self%v(i) = self%v(i)*atmp(i)
  end do
end if
IF(dealloc_flag)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_mult_cvec
!---------------------------------------------------------------------------------
!> Scale vector by a scalar
!---------------------------------------------------------------------------------
subroutine cvec_scale(self,alpha)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: alpha !< Scale factor
integer(i4) :: i
DEBUG_STACK_PUSH
!$omp parallel do if(self%n>OFT_OMP_VTHRESH)
do i=1,self%n
  self%v(i) = alpha*self%v(i)
end do
DEBUG_STACK_POP
end subroutine cvec_scale
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function cvec_dot_vec(self,a) result(dot)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
complex(c8) :: dot !< \f$ \sum_i self_i a_i \f$
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
ALLOCATE(atmp(self%n))
select type(this=>a)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','cvec_dot_vec',__FILE__)
    atmp=(1.d0,0.d0)*this%v
    dot=oft_global_dp(self%stitch_info,self%v,atmp,this%n)
  class default
    CALL oft_abort('"a" is not a native vector.','cvec_dot_vec',__FILE__)
end select
DEALLOCATE(atmp)
DEBUG_STACK_POP
end function cvec_dot_vec
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function cvec_dot_cvec(self,a) result(dot)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
complex(c8) :: dot !< \f$ \sum_i self_i a_i \f$
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_cvector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_dot_cvec',__FILE__)
    dot=oft_global_dp(self%stitch_info,self%v,this%v,this%n)
  class default
    CALL oft_abort('"a" is not a native vector.','vec_dot_cvec',__FILE__)
end select
DEBUG_STACK_POP
end function cvec_dot_cvec
!---------------------------------------------------------------------------------
!> Dot product with an array of vectors
!---------------------------------------------------------------------------------
function cvec_mdot_vec(self,a,n) result(dots)
class(oft_native_cvector), intent(inout) :: self !< Vector object
type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
integer(i4) :: i
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
ALLOCATE(atmp(self%n))
DO i=1,n
  select type(this=>a(i)%f)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','cvec_mdot_vec',__FILE__)
    atmp=(1.d0,0.d0)*this%v
    dots(i)=oft_global_dp(self%stitch_info,self%v,atmp,this%n,no_reduce=.TRUE.)
  class default
    CALL oft_abort('"a" is not a native vector.','cvec_mdot_vec',__FILE__)
  end select
END DO
DEALLOCATE(atmp)
dots=oft_mpi_sum(dots,n)
DEBUG_STACK_POP
end function cvec_mdot_vec
!---------------------------------------------------------------------------------
!> Dot product with an array of vectors
!---------------------------------------------------------------------------------
function cvec_mdot_cvec(self,a,n) result(dots)
class(oft_native_cvector), intent(inout) :: self !< Vector object
type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
integer(i4) :: i
DEBUG_STACK_PUSH
DO i=1,n
  select type(this=>a(i)%f)
  class is(oft_native_cvector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','cvec_mdot_cvec',__FILE__)
    dots(i)=oft_global_dp(self%stitch_info,self%v,this%v,this%n,no_reduce=.TRUE.)
  class default
    CALL oft_abort('"a" is not a native complex vector.','cvec_mdot_cvec',__FILE__)
  end select
END DO
dots=oft_mpi_sum(dots,n)
DEBUG_STACK_POP
end function cvec_mdot_cvec
!---------------------------------------------------------------------------------
!> Sum reduction over vector
!---------------------------------------------------------------------------------
function cvec_sum(self) result(sum)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8) :: sum !< Sum of vector elements
DEBUG_STACK_PUSH
sum=oft_global_reduction(self%stitch_info,self%v,self%n)
DEBUG_STACK_POP
end function cvec_sum
!---------------------------------------------------------------------------------
!> Compute norm of vector
!---------------------------------------------------------------------------------
function cvec_norm(self,itype) result(norm)
class(oft_native_cvector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
REAL(r8) :: norm !< Specified norm of vector
integer(i4) :: i
REAL(r8) :: nloc
REAL(r8), allocatable, dimension(:) :: vtmp
DEBUG_STACK_PUSH
SELECT CASE(itype)
  CASE(1)
    ALLOCATE(vtmp(self%n))
    DO i=1,self%n
      vtmp(i)=ABS(self%v(i))
    END DO
    norm=oft_global_reduction(self%stitch_info,vtmp,self%n)
    DEALLOCATE(vtmp)
  CASE(2)
    norm=SQRT(REAL(self%dot(self),r8))
  CASE(3)
    nloc=0.d0
    DO i=1,self%n
      nloc=MAX(nloc,ABS(self%v(i)))
    END DO
    norm=oft_mpi_max(nloc)
  CASE DEFAULT
    CALL oft_abort("Invalid norm type","cvec_norm",__FILE__)
END SELECT
DEBUG_STACK_POP
end function cvec_norm
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine cvec_stitch(self,up_method)
class(oft_native_cvector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: up_method
DEBUG_STACK_PUSH
call oft_global_stitch(self%stitch_info,self%v,up_method)
DEBUG_STACK_POP
end subroutine cvec_stitch
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine cvec_delete(self)
class(oft_native_cvector), intent(inout) :: self !< Vector object
DEBUG_STACK_PUSH
if(associated(self%v))deallocate(self%v)
if(associated(self%local_tmp))deallocate(self%local_tmp)
self%n=-1
self%nslice=-1
self%ng=-1
self%nblocks=1
NULLIFY(self%map,self%stitch_info)
DEBUG_STACK_POP
end subroutine cvec_delete
!---------------------------------------------------------------------------------
!> Cast an abstract matrix object to @ref oft_native_matrix
!!
!! The source matrix must be @ref oft_native_matrix or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!---------------------------------------------------------------------------------
FUNCTION native_matrix_cast(self,source) result(success)
class(oft_native_matrix), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_matrix), target, intent(in) :: source !< Abstract vector object
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_matrix)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION native_matrix_cast
!---------------------------------------------------------------------------------
!> Compute matrix-vector product
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine mat_apply_vec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j
real(r8) :: tmp
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','mat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','mat_apply_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_vector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=0.d0
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=0.d0
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(INT(self%lc_small(j),4)+lc_offset)
    end do
    btmp(i)=tmp
  end do
ELSE
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
   do i=1,self%nr
    tmp=0.d0
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(self%lc(j))
    end do
    btmp(i)=tmp
  end do
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_apply_vec
!---------------------------------------------------------------------------------
!> Compute matrix-vector product (complex)
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine mat_apply_cvec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source vector
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j
complex(c8) :: tmp
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','mat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','mat_apply_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(INT(self%lc_small(j),4)+lc_offset)
    end do
    btmp(i)=tmp
  end do
ELSE
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(self%lc(j))
    end do
    btmp(i)=tmp
  end do
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine mat_applyt_vec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j,k
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','mat_applyt_vec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','mat_applyt_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_vector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=0.d0
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=INT(self%lc_small(j),4)+lc_offset
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
ELSE
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=self%lc(j)
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose (complex vector)
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine mat_applyt_cvec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source vector
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j,k
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','mat_applyt_cvec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','mat_applyt_cvec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=INT(self%lc_small(j),4)+lc_offset
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
ELSE
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=self%lc(j)
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4) :: i,j,jp,jn,jk
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ANY(i_inds>self%nr))call oft_abort('Index exceeds graph','mat_set_values',__FILE__)
if(ANY(j_inds>self%nc))call oft_abort('Index exceeds graph','mat_set_values',__FILE__)
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  do i=1,n
    jp=self%map(iblock,jblock)%kr(i_inds(i))
    jn=self%map(iblock,jblock)%kr(i_inds(i)+1)-jp
    do j=1,m
      jk=search_array(j_inds(j),self%map(iblock,jblock)%lc(jp:jp+jn-1),jn)+jp-1
      self%M(self%map(iblock,jblock)%lc_map(jk))=b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
    'mat_set_values',__FILE__)
  do i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    do j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      self%M(jk)=b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine mat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,jp,jn,jk,jktmp(n,m)
logical :: fail_test,use_cache
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
self%full_current=.FALSE.
fail_test=.FALSE.
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  use_cache=.FALSE.
  IF(PRESENT(loc_cache))THEN
    IF(loc_cache(1,1)>0)THEN
      use_cache=.TRUE.
      jktmp=loc_cache
    END IF
  END IF
  IF(.NOT.use_cache)THEN
    DO i=1,n
      jp=self%map(iblock,jblock)%ext(1,i_inds(i))
      jn=self%map(iblock,jblock)%ext(2,i_inds(i))
      DO j=1,m
        jktmp(i,j)=search_array(self%j_map(jblock)%offset+j_inds(j),self%lc(jp:jn),jn-jp+1)
        IF(jktmp(i,j)==0)fail_test=.TRUE.
      END DO
    END DO
    IF(PRESENT(loc_cache))loc_cache=jktmp
  END IF
  DO i=1,n
    jp=self%map(iblock,jblock)%ext(1,i_inds(i))
    DO j=1,m
      jk=jktmp(i,j)+jp-1
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
    'matrix_add_values',__FILE__)
  DO i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    DO j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      IF(jk==jp-1)fail_test=.TRUE.
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
END IF
IF(fail_test)CALL oft_abort('Entry not found!','mat_add_values',__FILE__)
DEBUG_STACK_POP
end subroutine mat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine mat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,jp,jn,jk,jktmp(n,m)
logical :: fail_test,use_cache
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
self%full_current=.FALSE.
fail_test=.FALSE.
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  use_cache=.FALSE.
  IF(PRESENT(loc_cache))THEN
    IF(loc_cache(1,1)>0)THEN
      use_cache=.TRUE.
      jktmp=loc_cache
    END IF
  END IF
  IF(.NOT.use_cache)THEN
    DO i=1,n
      jp=self%map(iblock,jblock)%ext(1,i_inds(i))
      jn=self%map(iblock,jblock)%ext(2,i_inds(i))
      DO j=1,m
        jktmp(i,j)=search_array(self%j_map(jblock)%offset+j_inds(j),self%lc(jp:jn),jn-jp+1)
        IF(jktmp(i,j)==0)fail_test=.TRUE.
      END DO
    END DO
    IF(PRESENT(loc_cache))loc_cache=jktmp
  END IF
  DO i=1,n
    jp=self%map(iblock,jblock)%ext(1,i_inds(i))
    DO j=1,m
      jk=jktmp(i,j)+jp-1
      !$omp atomic
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
    'matrix_add_values',__FILE__)
  DO i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    DO j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      IF(jk==jp-1)fail_test=.TRUE.
      !$omp atomic
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
END IF
IF(fail_test)CALL oft_abort('Entry not found!','mat_add_values',__FILE__)
DEBUG_STACK_POP
end subroutine mat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine mat_assemble(self,diag)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
integer(i4) :: i,jp,jn,jr
DEBUG_STACK_PUSH
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_vector)
    !---Get diagonal matrix values
    !$omp parallel do private(jp,jn,jr)
    do i=1,self%nr ! loop over dof
      jp=self%kr(i); jn=self%kr(i+1)-jp
      jr=search_array(i,self%lc(jp:jp+jn-1),jn)+jp-1
      this%v(i)=self%M(jr)
    enddo
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a vector object.','mat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add(0.d0,1.d0,self%D)
end if
!---Common assembly tasks
IF(self%nc<2*lc_offset)THEN
  IF(.NOT.ASSOCIATED(self%lc_small))THEN
    ALLOCATE(self%lc_small(self%nnz))
    !$omp simd
    DO i=1,self%nnz
      self%lc_small(i)=INT(self%lc(i)-lc_offset,2)
    END DO
  END IF
END IF
DEBUG_STACK_POP
end subroutine mat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine mat_zero(self)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
DEBUG_STACK_PUSH
self%full_current=.FALSE.
self%M=0.d0
DEBUG_STACK_POP
end subroutine mat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine mat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
integer(i4) :: i,j,ir
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(iblock))THEN
  !$omp parallel do private(ir,j)
  do i=1,nrows
    ir=self%i_map(iblock)%offset+irows(i)
    do j=self%kr(ir),self%kr(ir+1)-1
      IF(zero_diag.OR.(self%lc(j)/=irows(i)))self%M(j)=0.d0
    end do
  end do
ELSE
  !$omp parallel do private(j)
  DO i=1,nrows
    DO j=self%kr(irows(i)),self%kr(irows(i)+1)-1
      IF(zero_diag.OR.(self%lc(j)/=irows(i)))THEN
        self%M(j)=0.d0
      END IF
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine mat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine mat_zero_cols(self,ncols,icols,jblock,keep_diag)
class(oft_native_matrix), intent(inout) :: self
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
!---
integer(i4) :: i,j
logical :: zero_diag
!---
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(jblock))THEN
  !$omp parallel do private(j)
  DO i=1,self%nr
    do j=self%kr(i),self%kr(i+1)-1
      IF(ANY(self%lc(j)==icols+self%i_map(jblock)%offset))THEN
        IF(zero_diag.OR.(self%lc(j)/=i))self%M(j)=0.d0
      END IF
    end do
  end do
ELSE
  !$omp parallel do private(j)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      IF(ANY(self%lc(j)==icols))THEN
        IF(zero_diag.OR.(self%lc(j)/=i))self%M(j)=0.d0
      END IF
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine mat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine mat_delete(self)
class(oft_native_matrix), intent(inout) :: self
integer(i4) :: i,j
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ASSOCIATED(self%M))deallocate(self%M)
if(ASSOCIATED(self%Md))deallocate(self%Md)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
IF(ASSOCIATED(self%kr))DEALLOCATE(self%kr)
IF(ASSOCIATED(self%lc))DEALLOCATE(self%lc)
IF(ASSOCIATED(self%lc_small))DEALLOCATE(self%lc_small)
!---
if(ASSOCIATED(self%map))THEN
  DO i=1,self%ni
    DO j=1,self%nj
      IF(ASSOCIATED(self%map(i,j)%lc_map))DEALLOCATE(self%map(i,j)%lc_map)
      IF(ASSOCIATED(self%map(i,j)%ext))DEALLOCATE(self%map(i,j)%ext)
    END DO
  END DO
  DEALLOCATE(self%map,self%i_map,self%j_map)
END IF
IF(ASSOCIATED(self%linkage))THEN
  IF(ASSOCIATED(self%linkage%send))THEN
    DO i=0,self%linkage%nproc_con
      IF(ASSOCIATED(self%linkage%send(i)%v))DEALLOCATE(self%linkage%send(i)%v,self%linkage%recv(i)%v)
      IF(ASSOCIATED(self%linkage%csend(i)%v))DEALLOCATE(self%linkage%csend(i)%v,self%linkage%crecv(i)%v)
    END DO
    DEALLOCATE(self%linkage%send,self%linkage%recv)
    DEALLOCATE(self%linkage%csend,self%linkage%crecv)
  END IF
  IF(ASSOCIATED(self%linkage%kle))DEALLOCATE(self%linkage%kle)
  IF(ASSOCIATED(self%linkage%lbe))DEALLOCATE(self%linkage%lbe)
  IF(ASSOCIATED(self%linkage%lle))DEALLOCATE(self%linkage%lle)
  DEALLOCATE(self%linkage)
END IF
IF(ASSOCIATED(self%Mfull))DEALLOCATE(self%Mfull)
IF(ASSOCIATED(self%redind))DEALLOCATE(self%redind)
IF(ASSOCIATED(self%color))DEALLOCATE(self%color)
self%nred=0
DEBUG_STACK_POP
end subroutine mat_delete
!---------------------------------------------------------------------------------
!> Build full local submatrix for matrix
!! - Create stitching structures
!---------------------------------------------------------------------------------
subroutine native_matrix_setup_full(self,u)
class(oft_native_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: u !< Constructed vector for parent row space
class(oft_native_vector), pointer :: uv
type :: eout
  integer(i8), pointer, dimension(:,:) :: le => NULL()
end type eout
type :: emap
  integer(i4), pointer, dimension(:) :: le => NULL()
  integer(i4), pointer, dimension(:,:) :: linktmp => NULL()
end type emap
type(eout), allocatable, dimension(:) :: lerecv,lesend
type(emap), allocatable, dimension(:) :: lrmaps
integer(i4) :: i,j,k,m,mm,ierr,nbetmp,nbemax,jp,jn,jk,neout,netmp,neel
integer(i4), allocatable, dimension(:) :: isort,lsort,ind_sort
integer(i4), allocatable, dimension(:) :: lle_tmp,ncon,nrecv
integer(i8) :: last
integer(i8), allocatable, dimension(:) :: lge,sort_list
integer(i8), pointer, dimension(:,:) :: letmp,leout
logical :: periodic
logical, allocatable, dimension(:) :: glob_inc
logical, allocatable, dimension(:,:) :: skips
IF(ASSOCIATED(self%linkage))RETURN
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(2))WRITE(*,*)'Creating CRS stitching structures'
SELECT TYPE(u)
CLASS IS(oft_native_vector)
  uv=>u
CLASS DEFAULT
  CALL oft_abort('"u" is not a vector object.','native_matrix_setup_full',__FILE__)
END SELECT
IF(uv%stitch_info%skip)THEN
  ALLOCATE(self%linkage)
  self%linkage%full=.TRUE.
  self%linkage%skip=.TRUE.
  DEBUG_STACK_POP
  RETURN
END IF
!---Get boundary entry count
periodic=.FALSE.
ALLOCATE(lge(self%nr))
DO i=1,self%ni
  periodic=(periodic.OR.self%i_map(i)%per)
  !$omp parallel do if(self%i_map(i)%n>OFT_OMP_VTHRESH)
  DO j=1,self%i_map(i)%n
    lge(self%i_map(i)%offset+j)=self%i_map(i)%offsetg+self%i_map(i)%lge(j)
  END DO
END DO
IF(MINVAL(lge)<=0)CALL oft_abort('BAD LGE value','native_matrix_setup_full',__FILE__)
!---
ALLOCATE(ncon(0:uv%stitch_info%nproc_con),nrecv(uv%stitch_info%nproc_con))
ncon=0
ALLOCATE(self%linkage)
self%linkage%nproc_con=uv%stitch_info%nproc_con
self%linkage%proc_split=uv%stitch_info%proc_split
self%linkage%proc_con=>uv%stitch_info%proc_con
self%linkage%send_reqs=>uv%stitch_info%send_reqs
self%linkage%recv_reqs=>uv%stitch_info%recv_reqs
ALLOCATE(self%linkage%kle(0:self%linkage%nproc_con+1)) ! Allocate point linkage arrays
self%linkage%kle=0
nbetmp=0
ALLOCATE(glob_inc(uv%n))
glob_inc=.FALSE.
IF(periodic)THEN
  DO i=uv%stitch_info%kle(0),uv%stitch_info%kle(1)-1
    k = uv%stitch_info%lbe(uv%stitch_info%lle(2,i))
    glob_inc(k)=.TRUE.
  END DO
END IF
!$omp parallel do private(j,k) reduction(+:nbetmp) schedule(guided)
DO i=1,uv%stitch_info%nbe
  k=uv%stitch_info%lbe(i)
  IF(glob_inc(k))THEN
    nbetmp=nbetmp+(self%kr(k+1)-self%kr(k))
  ELSE
    DO j=self%kr(k),self%kr(k+1)-1
      IF(uv%stitch_info%be(self%lc(j)))nbetmp=nbetmp+1
    END DO
  END IF
END DO
IF(periodic)THEN
  !$omp parallel do private(j) reduction(+:nbetmp) schedule(guided)
  DO i=1,uv%n
    IF(uv%stitch_info%be(i))CYCLE
    DO j=self%kr(i),self%kr(i+1)-1
      IF(glob_inc(self%lc(j)))nbetmp=nbetmp+1
    END DO
  END DO
END IF
self%linkage%nbe=nbetmp
self%linkage%full=uv%stitch_info%full
nbemax=nbetmp
IF(.NOT.uv%stitch_info%full)nbemax=oft_mpi_max(nbetmp)
IF(oft_debug_print(3))WRITE(*,'(2A,I8)')oft_indent,'Max # of boundary elements',nbemax
self%linkage%nbemax=nbemax
!---Allocate temporary Send/Recv arrays
ALLOCATE(lerecv(0:self%linkage%nproc_con),lesend(0:self%linkage%nproc_con),lrmaps(0:self%linkage%nproc_con))
ALLOCATE(self%linkage%lbe(self%linkage%nbe))
ALLOCATE(lesend(0)%le(2,nbemax))
ALLOCATE(lrmaps(0)%le(nbemax))
ALLOCATE(lrmaps(0)%linktmp(2,2*nbemax))
lrmaps(0)%linktmp=0
!---
allocate(lle_tmp(uv%stitch_info%nle))
lle_tmp=uv%stitch_info%lle(1,:)
DO j=0,self%linkage%nproc_con
  jp=uv%stitch_info%kle(j)
  jn=uv%stitch_info%kle(j+1)-jp
  IF(jn>0)call sort_array(lle_tmp(jp:jp+jn-1),jn)
END DO
!---Point dummy output array to Send array
ALLOCATE(skips(0:self%linkage%nproc_con,self%nr))
skips=.TRUE.
DO i=1,uv%stitch_info%nbe
  DO m=0,self%linkage%nproc_con
    jp=uv%stitch_info%kle(m)
    jn=uv%stitch_info%kle(m+1)-jp
    IF(jn>0)THEN
      jk=search_array(i,lle_tmp(jp:jp+jn-1),jn)
      IF(jk/=0)skips(m,uv%stitch_info%lbe(i))=.FALSE.
    END IF
  END DO
END DO
DEALLOCATE(lle_tmp)
!---
ALLOCATE(sort_list(self%linkage%nbe),ind_sort(self%linkage%nbe),letmp(2,self%linkage%nbe))
nbetmp=0
ind_sort=(/(i,i=1,self%linkage%nbe)/)
!$omp parallel do private(j,k,m) schedule(guided)
DO i=1,uv%stitch_info%nbe
  k=uv%stitch_info%lbe(i)
  IF(glob_inc(k))THEN
    !$omp atomic capture
    m=nbetmp
    nbetmp=nbetmp+self%kr(k+1)-self%kr(k)
    !$omp end atomic
    DO j=self%kr(k),self%kr(k+1)-1
      m=m+1
      self%linkage%lbe(m) = j
      letmp(:,m) = (/k,self%lc(j)/)
      sort_list(m) = lge(k)
    END DO
  ELSE
    DO j=self%kr(k),self%kr(k+1)-1
      IF(uv%stitch_info%be(self%lc(j)))THEN
        !$omp atomic capture
        nbetmp=nbetmp+1
        m=nbetmp
        !$omp end atomic
        self%linkage%lbe(m) = j
        letmp(:,m) = (/k,self%lc(j)/)
        sort_list(m) = lge(k)
      END IF
    END DO
  END IF
END DO
IF(periodic)THEN
  !$omp parallel do private(j,m) schedule(guided)
  DO i=1,uv%n
    IF(uv%stitch_info%be(i))CYCLE
    DO j=self%kr(i),self%kr(i+1)-1
      IF(glob_inc(self%lc(j)))THEN
        !$omp atomic capture
        nbetmp=nbetmp+1
        m=nbetmp
        !$omp end atomic
        self%linkage%lbe(m) = j
        letmp(:,m) = (/i,self%lc(j)/)
        sort_list(m) = lge(i)
      END IF
    END DO
  END DO
END IF
DEALLOCATE(glob_inc)
CALL sort_array(sort_list,ind_sort,self%linkage%nbe)
sort_list(1)=lge(letmp(2,ind_sort(1)))
m=1
DO i=2,self%linkage%nbe
  sort_list(i)=lge(letmp(2,ind_sort(i)))
  IF(lge(letmp(1,ind_sort(m)))==lge(letmp(1,ind_sort(i))))CYCLE
  CALL sort_array(sort_list(m:i-1),ind_sort(m:i-1),i-m)
  m=i
END DO
CALL sort_array(sort_list(m:i-1),ind_sort(m:i-1),i-m)
!---
DO i=1,self%linkage%nbe
  k=ind_sort(i)
  lesend(0)%le(:,i) = lge(letmp(:,k))
  lrmaps(0)%le(i) = k
  DO m=1,self%linkage%nproc_con
    IF(skips(m,letmp(1,k)).OR.skips(m,letmp(2,k)))CYCLE
    ncon(m)=ncon(m)+1
  END DO
END DO
#ifdef HAVE_MPI
IF(.NOT.uv%stitch_info%full)THEN
  DO j=1,self%linkage%nproc_con
    CALL MPI_ISEND(ncon(j),1,OFT_MPI_I4,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','native_matrix_setup_full',__FILE__)
    CALL MPI_IRECV(nrecv(j),1,OFT_MPI_I4,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','native_matrix_setup_full',__FILE__)
  END DO
  CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%recv_reqs,ierr)
  DO i=1,self%linkage%nproc_con
    ALLOCATE(lesend(i)%le(2,ncon(i)))
    lesend(i)%le=0
    ALLOCATE(lrmaps(i)%le(ncon(i)))
    ALLOCATE(lrmaps(i)%linktmp(2,2*ncon(i)))
    lrmaps(i)%linktmp=0
    ALLOCATE(lerecv(i)%le(2,nrecv(i)))
  END DO
  CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%send_reqs,ierr)
  ncon=0
  DO i=1,self%linkage%nbe
    k=ind_sort(i)
    DO m=1,self%linkage%nproc_con
      IF(skips(m,letmp(1,k)).OR.skips(m,letmp(2,k)))CYCLE
      ncon(m)=ncon(m)+1
      lesend(m)%le(:,ncon(m)) = lge(letmp(:,k))
      lrmaps(m)%le(ncon(m)) = k
    END DO
  END DO
END IF
#endif
DEALLOCATE(letmp,skips,lge)
!---
#ifdef HAVE_MPI
IF(.NOT.uv%stitch_info%full)THEN
  CALL oft_mpi_barrier(ierr) ! Wait for all processes
  !---Create Send and Recv calls
  do j=1,self%linkage%nproc_con
    CALL MPI_ISEND(lesend(j)%le,2*ncon(j),OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','native_matrix_setup_full',__FILE__)
    CALL MPI_IRECV(lerecv(j)%le,2*nrecv(j),OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','native_matrix_setup_full',__FILE__)
  end do
  !---Loop over each connected processor
  do while(.TRUE.)
    IF(oft_mpi_check_reqs(self%linkage%nproc_con,self%linkage%recv_reqs))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(self%linkage%nproc_con,self%linkage%recv_reqs,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','native_matrix_setup_full',__FILE__)
    letmp=>lerecv(j)%le ! Point dummy input array to current Recv array
    leout=>lesend(j)%le ! Point dummy out array to current Send array
    !---Get extents of transfer arrays
    neout=ncon(j)
    netmp=nrecv(j)
    !---Determine location of boundary points on other processors
    ncon(j)=0
    neel=0
    m=1
    jp=1
    last=0
    do mm=1,netmp ! Loop over input points
      IF(letmp(1,mm)==last.AND.periodic)THEN
        m=jp
      ELSE
        last=letmp(1,mm)
      END IF
      IF(m>neout)EXIT
      do while(leout(1,m)<=letmp(1,mm))
        IF(leout(1,m)<letmp(1,mm))jp=m+1
        IF(leout(1,m)==letmp(1,mm))THEN
          IF(leout(2,m)==letmp(2,mm))THEN ! Found match
            ncon(j)=ncon(j)+1
            lrmaps(j)%linktmp(:,ncon(j))=(/mm,m/)
            neel=neel+1
          ELSE IF(leout(2,m)>letmp(2,mm))THEN
            EXIT
          END IF
        END IF
        m=m+1
        IF(m>neout)EXIT
      end do
    end do
    self%linkage%kle(j)=neel
  end do
END IF
#endif
!---Local connectivity
!---Determine location of boundary points on other processors
ncon(0)=0
neel=0
leout=>lesend(0)%le
jp=uv%stitch_info%kle(0)
jn=uv%stitch_info%kle(1)-jp
IF(jn>0)THEN
  DO m=1,self%linkage%nbe
    IF(leout(1,m)==0)EXIT
  END DO
  neout=m-1
  m=1
  jp=1
  last=0
  do mm=1,neout ! Loop over input points
    IF(leout(1,mm)==last.AND.periodic)THEN
      m=jp
    ELSE
      last=leout(1,mm)
    END IF
    IF(m>neout)EXIT
    do while(leout(1,m)<=leout(1,mm))
      IF(leout(1,m)<leout(1,mm))jp=m+1
      IF((leout(1,m)==leout(1,mm)).AND.(m/=mm))THEN
        IF(leout(2,m)==leout(2,mm))THEN ! Found match
          ncon(0)=ncon(0)+1
          lrmaps(0)%linktmp(:,ncon(0))=(/mm,m/)
          neel=neel+1
        ELSE IF(leout(2,m)>leout(2,mm))THEN
          EXIT
        END IF
      END IF
      m=m+1
      IF(m>neout)EXIT
    end do
  end do
END IF
self%linkage%kle(0)=neel
!---Condense linkage to sparse rep
self%linkage%nle=SUM(self%linkage%kle)
self%linkage%kle(self%linkage%nproc_con+1)=self%linkage%nle+1
do i=self%linkage%nproc_con,0,-1 ! cumulative unique point linkage count
  self%linkage%kle(i)=self%linkage%kle(i+1)-self%linkage%kle(i)
end do
!---
if(self%linkage%kle(0)/=1)call oft_abort('Bad element linkage count','native_matrix_setup_full', &
__FILE__)
!---
allocate(self%linkage%lle(2,self%linkage%nle))
allocate(self%linkage%send(0:self%linkage%nproc_con),self%linkage%recv(0:self%linkage%nproc_con))
allocate(self%linkage%csend(0:self%linkage%nproc_con),self%linkage%crecv(0:self%linkage%nproc_con))
!---
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do schedule(dynamic,1)
do i=0,self%linkage%nproc_con
  !---
  DO j=1,ncon(i)
    lsort(j)=lrmaps(i)%linktmp(2,j)
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%linkage%proc_con(i)<oft_env%rank)lsort(j)=lrmaps(i)%linktmp(1,j)
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  IF(i==0)THEN
    DO m=0,ncon(i)-1
      self%linkage%lle(:,m+self%linkage%kle(i)) = &
      (/lrmaps(i)%le(lrmaps(i)%linktmp(2,isort(m+1))),lrmaps(i)%le(lrmaps(i)%linktmp(1,isort(m+1)))/)
    END DO
  ELSE
    DO m=0,ncon(i)-1
      self%linkage%lle(:,m+self%linkage%kle(i)) = &
      (/lrmaps(i)%le(lrmaps(i)%linktmp(2,isort(m+1))),-1/)
    END DO
  END IF
  !---Allocate permanent stitching arrays
  self%linkage%send(i)%n=ncon(i); self%linkage%recv(i)%n=ncon(i)
  allocate(self%linkage%send(i)%v(self%linkage%send(i)%n))
  allocate(self%linkage%recv(i)%v(self%linkage%recv(i)%n))
  self%linkage%csend(i)%n=ncon(i); self%linkage%crecv(i)%n=ncon(i)
  allocate(self%linkage%csend(i)%v(self%linkage%csend(i)%n))
  allocate(self%linkage%crecv(i)%v(self%linkage%crecv(i)%n))
end do
DEALLOCATE(lsort,isort)
!$omp end parallel
DEALLOCATE(ncon,nrecv)
DEALLOCATE(lesend(0)%le)
DEALLOCATE(lrmaps(0)%le)
DEALLOCATE(lrmaps(0)%linktmp)
IF(.NOT.uv%stitch_info%full)THEN
  CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','native_matrix_setup_full',__FILE__)
  CALL oft_mpi_barrier(ierr) ! Wait for all processes
  DO i=1,self%linkage%nproc_con
    DEALLOCATE(lesend(i)%le)
    DEALLOCATE(lrmaps(i)%le)
    DEALLOCATE(lrmaps(i)%linktmp)
    DEALLOCATE(lerecv(i)%le)
  END DO
  DEALLOCATE(lesend,lrmaps,lerecv)
END IF
!---
CALL oft_stitch_check(self%linkage)
!---Flag redundant rows (from periodicity)
IF(periodic)THEN
  IF(oft_debug_print(3))WRITE(*,*)'  Flagging redundant rows'
  ALLOCATE(glob_inc(uv%n))
  glob_inc=.FALSE.
  DO i=uv%stitch_info%kle(0),uv%stitch_info%kle(1)-1
    k=uv%stitch_info%lbe(uv%stitch_info%lle(1,i))
    j=uv%stitch_info%lbe(uv%stitch_info%lle(2,i))
    IF(k>j)glob_inc(k)=.TRUE.
  END DO
  self%nred=COUNT(glob_inc)
  ALLOCATE(self%redind(self%nred))
  self%nred=0
  DO i=1,uv%stitch_info%nbe
    IF(glob_inc(uv%stitch_info%lbe(i)))THEN
      self%nred=self%nred+1
      self%redind(self%nred)=uv%stitch_info%lbe(i)
    END IF
  END DO
  DEALLOCATE(glob_inc)
END IF
! IF(oft_debug_print(3))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine native_matrix_setup_full
!---------------------------------------------------------------------------------
!> Update full local matrix with values from other processors
!---------------------------------------------------------------------------------
subroutine matrix_update_slice(self)
class(oft_native_matrix), intent(inout) :: self
integer(i4) :: i,j,k
IF(self%full_current)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(3))WRITE(*,*)'  Get local slice'
IF(.NOT.ASSOCIATED(self%Mfull))ALLOCATE(self%Mfull(self%nnz))
!---Get values
DO i=1,self%nnz
  self%Mfull(i)=self%M(i)
END DO
CALL oft_global_stitch(self%linkage,self%Mfull,1)
!---Zero out redundant rows
IF(ASSOCIATED(self%redind))THEN
  !$omp parallel do private(j,k)
  DO i=1,self%nred
    k=self%redind(i)
    DO j=self%kr(k),self%kr(k+1)-1
      self%Mfull(j)=0.d0
      IF(self%lc(j)==k)self%Mfull(j)=1.d0
    END DO
  END DO
END IF
self%full_current=.TRUE.
! IF(oft_debug_print(3))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine matrix_update_slice
!---------------------------------------------------------------------------------
!> Setup local submatrix from parent matrix
!! - Create stitching structures for parent matrix if necessary
!---------------------------------------------------------------------------------
subroutine submatrix_setup(self,mat,u,slice,part)
class(oft_native_submatrix), intent(inout) :: self
class(oft_matrix), intent(inout) :: mat !< Parent matrix
class(oft_vector), intent(inout) :: u !< Constructed vector for parent row space
integer(i4), optional, intent(in) :: slice !< Index of field slice to use (optional)
integer(i4), optional, target, intent(in) :: part(:) !< Array of elements to define submatrix [:] (optional)
class(oft_native_matrix), pointer :: parent
DEBUG_STACK_PUSH
!---
IF(PRESENT(part).AND.PRESENT(slice))CALL oft_abort('"slice" and "part" cannot both be specified', &
  'submatrix_setup',__FILE__)
IF(PRESENT(part))self%part=>part
IF(PRESENT(slice))self%slice=slice
IF(.NOT.native_matrix_cast(parent,mat))CALL oft_abort('Native matrix required', &
  'submatrix_setup',__FILE__)
self%parent=>parent
IF(self%slice<0.AND.ASSOCIATED(self%part))THEN
  self%nr=SIZE(self%part)
ELSE IF(self%slice<0.AND.(.NOT.ASSOCIATED(self%part)))THEN
  self%nr=self%parent%nrslice
ELSE IF(self%slice>0)THEN
  self%nr=parent%i_map(self%slice)%nslice
END IF
self%nrg=self%nr
self%nc=self%nr
self%ncg=self%nr
self%ni=1
self%nj=1
self%nrslice=self%nr
self%ncslice=self%nr
ALLOCATE(self%i_map(1),self%j_map(1),self%map(1,1))
self%i_map(1)%n=self%nr
self%j_map(1)%n=self%nr
!---
CALL native_matrix_setup_full(self%parent,u)
CALL self%update_slice
! IF(oft_debug_print(3))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine submatrix_setup
!---------------------------------------------------------------------------------
!> Update matrix with values from parent matrix
!---------------------------------------------------------------------------------
subroutine submatrix_update_slice(self)
class(oft_native_submatrix), intent(inout) :: self
class(oft_native_matrix), pointer :: parent
integer(i4) :: i,j,k,jp,jn
integer(i4), allocatable, dimension(:) :: slice_map
DEBUG_STACK_PUSH
IF(oft_debug_print(3))WRITE(*,*)'  Get local slice'
parent=>self%parent
!---Get graph
IF(.NOT.ASSOCIATED(self%lcmap))THEN
  !---
  ALLOCATE(slice_map(parent%nr))
  slice_map=0
  IF(self%slice>0)THEN
    i=self%slice
    !$omp parallel do
    DO j=1,parent%i_map(i)%nslice
      slice_map(parent%i_map(i)%offset+parent%i_map(i)%slice(j))=j
    END DO
  ELSE
    IF(ASSOCIATED(self%part))THEN
      !$omp parallel do
      DO i=1,self%nr
        slice_map(self%part(i))=i
      END DO
    ELSE
      DO i=1,parent%ni
        !$omp parallel do
        DO j=1,parent%i_map(i)%nslice
          slice_map(parent%i_map(i)%offset+parent%i_map(i)%slice(j))=parent%i_map(i)%soffset+j
        END DO
      END DO
    END IF
  END IF
  !---
  ALLOCATE(self%kr(self%nr+1))
  self%kr=0
  !$omp parallel do private(j) schedule(static,20)
  do i=1,parent%nr
    IF(slice_map(i)==0)CYCLE
    do j=parent%kr(i),parent%kr(i+1)-1
      IF(slice_map(parent%lc(j))/=0)self%kr(slice_map(i))=self%kr(slice_map(i))+1
    end do
  end do
  !---
  self%nnz=sum(self%kr)
  self%kr(self%nr+1)=self%nnz+1
  do i=self%nr,1,-1 ! cumulative point to point count
    self%kr(i)=self%kr(i+1)-self%kr(i)
  end do
  if(self%kr(1)/=1)call oft_abort('Bad element to element count', &
    'submatrix_update_slice',__FILE__)
  !---
  ALLOCATE(self%lcmap(self%nnz))
  ALLOCATE(self%lc(self%nnz))
  do i=1,parent%nr
    IF(slice_map(i)==0)CYCLE
    k=self%kr(slice_map(i))
    do j=parent%kr(i),parent%kr(i+1)-1
      IF(slice_map(parent%lc(j))/=0)THEN
        self%lc(k)=slice_map(parent%lc(j))
        self%lcmap(k)=j
        k=k+1
      END IF
    end do
    jp=self%kr(slice_map(i))
    jn=self%kr(slice_map(i)+1)-1
    CALL sort_array(self%lc(jp:jn),self%lcmap(jp:jn),jn-jp+1)
  end do
  DEALLOCATE(slice_map)
  ALLOCATE(self%M(self%nnz))
END IF
!---Get values
CALL matrix_update_slice(parent)
!$omp parallel do
DO i=1,self%nnz
  self%M(i)=parent%Mfull(self%lcmap(i))
END DO
self%full_current=.FALSE.
! IF(oft_debug_print(3))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine submatrix_update_slice
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine submatrix_delete(self)
class(oft_native_submatrix), intent(inout) :: self
DEBUG_STACK_PUSH
CALL mat_delete(self)
IF(ASSOCIATED(self%lcmap))DEALLOCATE(self%lcmap)
IF(ASSOCIATED(self%part))NULLIFY(self%part)
NULLIFY(self%parent)
DEBUG_STACK_POP
end subroutine submatrix_delete
!---------------------------------------------------------------------------------
!> Perform graph partitioning (METIS)
!---------------------------------------------------------------------------------
subroutine partition_graph(graph,n,part)
type(oft_graph), intent(inout) :: graph !< Graph to partition
integer(i4), intent(in) :: n !< Desired number of partitions
integer(i4), intent(inout) :: part(:) !< Array of partition ID for each element [graph%nr]
integer(i4) :: ptype,ierr
DEBUG_STACK_PUSH
ptype=2 ! Use ML K-Way
CALL oft_metis_partGraph(graph%nr,graph%nnz,graph%kr,graph%lc,n,part,ptype,ierr)
IF(ierr<0)CALL oft_abort('Graph paritioning failed','partition_graph',__FILE__)
DEBUG_STACK_POP
end subroutine partition_graph
!---------------------------------------------------------------------------------
!> Cast a matrix object to a oft_native_matrix
!!
!! The source matrix must be oft_native_matrix or a child class, otherwise an error will be thrown
!---------------------------------------------------------------------------------
FUNCTION native_cmatrix_cast(self,source) result(ierr)
class(oft_native_cmatrix), pointer, intent(out) :: self !< Pointer to cast crsmatrix
class(oft_cmatrix), target, intent(in) :: source !< Source matrix to cast
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cmatrix)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION native_cmatrix_cast
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine cmat_apply_vec(self,a,b)
class(oft_native_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
complex(c8) :: tmp
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: btmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','cmat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','cmat_apply_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(INT(self%lc_small(j),4)+lc_offset)
    end do
    btmp(i)=tmp
  end do
ELSE
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(self%lc(j))
    end do
    btmp(i)=tmp
  end do
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine cmat_apply_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine cmat_apply_cvec(self,a,b)
class(oft_native_cmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j
complex(c8) :: tmp
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','cmat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','cmat_apply_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(INT(self%lc_small(j),4)+lc_offset)
    end do
    btmp(i)=tmp
  end do
ELSE
  !$omp parallel do private(j,tmp) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  do i=1,self%nr
    tmp=(0.d0,0.d0)
    do j=self%kr(i),self%kr(i+1)-1
      tmp=tmp+self%M(j)*atmp(self%lc(j))
    end do
    btmp(i)=tmp
  end do
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine cmat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine cmat_applyt_vec(self,a,b)
class(oft_native_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j,k
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','cmat_applyt_vec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','cmat_applyt_vec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=INT(self%lc_small(j),4)+lc_offset
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
ELSE
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=self%lc(j)
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine cmat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine cmat_applyt_cvec(self,a,b)
class(oft_native_cmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
LOGICAL :: dealloc_flags(2)
integer(i4) :: i,j,k
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','cmat_applyt_cvec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','cmat_applyt_cvec',__FILE__)
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  atmp=>a%v
  dealloc_flags(1)=.FALSE.
CLASS DEFAULT
  NULLIFY(atmp)
  CALL a%get_local(atmp)
  dealloc_flags(1)=.TRUE.
END SELECT
SELECT TYPE(b)
CLASS IS(oft_native_cvector)
  btmp=>b%v
  dealloc_flags(2)=.FALSE.
CLASS DEFAULT
  NULLIFY(btmp)
  CALL b%get_local(btmp)
  dealloc_flags(2)=.TRUE.
END SELECT
!---Apply operator
btmp=(0.d0,0.d0)
IF(ASSOCIATED(self%lc_small))THEN
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=INT(self%lc_small(j),4)+lc_offset
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
ELSE
  !$omp parallel do private(j,k) schedule(static,20) if(self%nnz>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      k=self%lc(j)
      !$omp atomic
      btmp(k) = btmp(k) + self%M(j)*atmp(i)
    END DO
  END DO
END IF
IF(dealloc_flags(1))DEALLOCATE(atmp)
IF(dealloc_flags(2))THEN
  CALL b%restore_local(btmp,wait=.TRUE.)
  DEALLOCATE(btmp)
END IF
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine cmat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4) :: i,j,jp,jn,jk
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ANY(i_inds>self%nr))call oft_abort('Index exceeds graph','cmat_set_values',__FILE__)
if(ANY(j_inds>self%nc))call oft_abort('Index exceeds graph','cmat_set_values',__FILE__)
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  do i=1,n
    jp=self%map(iblock,jblock)%kr(i_inds(i))
    jn=self%map(iblock,jblock)%kr(i_inds(i)+1)-jp
    do j=1,m
      jk=search_array(j_inds(j),self%map(iblock,jblock)%lc(jp:jp+jn-1),jn)+jp-1
      self%M(self%map(iblock,jblock)%lc_map(jk))=b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'cmat_set_values',__FILE__)
  do i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    do j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      self%M(jk)=b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine cmat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,jp,jn,jk,jktmp(n,m)
logical :: fail_test,use_cache
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','cmat_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','cmat_add_values',__FILE__)
self%full_current=.FALSE.
fail_test=.FALSE.
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  use_cache=.FALSE.
  IF(PRESENT(loc_cache))THEN
    IF(loc_cache(1,1)>0)THEN
      use_cache=.TRUE.
      jktmp=loc_cache
    END IF
  END IF
  IF(.NOT.use_cache)THEN
    DO i=1,n
      jp=self%map(iblock,jblock)%ext(1,i_inds(i))
      jn=self%map(iblock,jblock)%ext(2,i_inds(i))
      DO j=1,m
        jktmp(i,j)=search_array(self%j_map(jblock)%offset+j_inds(j),self%lc(jp:jn),jn-jp+1)
        IF(jktmp(i,j)==0)fail_test=.TRUE.
      END DO
    END DO
    IF(PRESENT(loc_cache))loc_cache=jktmp
  END IF
  DO i=1,n
    jp=self%map(iblock,jblock)%ext(1,i_inds(i))
    DO j=1,m
      jk=jktmp(i,j)+jp-1
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'cmat_add_values',__FILE__)
  DO i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    DO j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      IF(jk==jp-1)fail_test=.TRUE.
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
END IF
IF(fail_test)CALL oft_abort('Entry not found!','cmat_add_values',__FILE__)
DEBUG_STACK_POP
end subroutine cmat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine cmat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,jp,jn,jk,jktmp(n,m)
logical :: fail_test,use_cache
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','cmat_atomic_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','cmat_atomic_add_values',__FILE__)
self%full_current=.FALSE.
fail_test=.FALSE.
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  use_cache=.FALSE.
  IF(PRESENT(loc_cache))THEN
    IF(loc_cache(1,1)>0)THEN
      use_cache=.TRUE.
      jktmp=loc_cache
    END IF
  END IF
  IF(.NOT.use_cache)THEN
    DO i=1,n
      jp=self%map(iblock,jblock)%ext(1,i_inds(i))
      jn=self%map(iblock,jblock)%ext(2,i_inds(i))
      DO j=1,m
        jktmp(i,j)=search_array(self%j_map(jblock)%offset+j_inds(j),self%lc(jp:jn),jn-jp+1)
        IF(jktmp(i,j)==0)fail_test=.TRUE.
      END DO
    END DO
    IF(PRESENT(loc_cache))loc_cache=jktmp
  END IF
  DO i=1,n
    jp=self%map(iblock,jblock)%ext(1,i_inds(i))
    DO j=1,m
      jk=jktmp(i,j)+jp-1
      !$omp atomic
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'cmat_atomic_add_values',__FILE__)
  DO i=1,n
    jp=self%kr(i_inds(i))
    jn=self%kr(i_inds(i)+1)-jp
    DO j=1,m
      jk=search_array(j_inds(j),self%lc(jp:jp+jn-1),jn)+jp-1
      IF(jk==jp-1)fail_test=.TRUE.
      !$omp atomic
      self%M(jk) = self%M(jk) + b(i,j)
    END DO
  END DO
END IF
IF(fail_test)CALL oft_abort('Entry not found!','cmat_atomic_add_values',__FILE__)
DEBUG_STACK_POP
end subroutine cmat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine cmat_assemble(self,diag)
class(oft_native_cmatrix), intent(inout) :: self
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
integer(i4) :: i,jp,jn,jr
DEBUG_STACK_PUSH
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_cvector)
    !---Get diagonal matrix values
    !$omp parallel do private(jp,jn,jr)
    do i=1,self%nr ! loop over dof
      jp=self%kr(i); jn=self%kr(i+1)-jp
      jr=search_array(i,self%lc(jp:jp+jn-1),jn)+jp-1
      this%v(i)=self%M(jr)
    enddo
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a vector object.','cmat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add((0.d0,0.d0),(1.d0,0.d0),self%D)
end if
!---Common assembly tasks
IF(self%nc<2*lc_offset)THEN
  IF(.NOT.ASSOCIATED(self%lc_small))THEN
    ALLOCATE(self%lc_small(self%nnz))
    DO i=1,self%nnz
      self%lc_small(i)=INT(self%lc(i)-lc_offset,2)
    END DO
  END IF
END IF
DEBUG_STACK_POP
end subroutine cmat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine cmat_zero(self)
class(oft_native_cmatrix), intent(inout) :: self
DEBUG_STACK_PUSH
self%full_current=.FALSE.
self%M=(0.d0,0.d0)
DEBUG_STACK_POP
end subroutine cmat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine cmat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
integer(i4) :: i,j,ir
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(iblock))THEN
  !$omp parallel do private(ir,j)
  do i=1,nrows
    ir=self%i_map(iblock)%offset+irows(i)
    do j=self%kr(ir),self%kr(ir+1)-1
      IF(zero_diag.OR.(self%lc(j)/=irows(i)))self%M(j)=(0.d0,0.d0)
    end do
  end do
ELSE
  !$omp parallel do private(j)
  DO i=1,nrows
    DO j=self%kr(irows(i)),self%kr(irows(i)+1)-1
      IF(zero_diag.OR.(self%lc(j)/=irows(i)))THEN
        self%M(j)=(0.d0,0.d0)
      END IF
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine cmat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine cmat_zero_cols(self,ncols,icols,jblock,keep_diag)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
!---
integer(i4) :: i,j
logical :: zero_diag
!---
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(jblock))THEN
  !$omp parallel do private(j)
  DO i=1,self%nr
    do j=self%kr(i),self%kr(i+1)-1
      IF(ANY(self%lc(j)==icols+self%i_map(jblock)%offset))THEN
        IF(zero_diag.OR.(self%lc(j)/=i))self%M(j)=(0.d0,0.d0)
      END IF
    end do
  end do
ELSE
  !$omp parallel do private(j)
  DO i=1,self%nr
    DO j=self%kr(i),self%kr(i+1)-1
      IF(ANY(self%lc(j)==icols))THEN
        IF(zero_diag.OR.(self%lc(j)/=i))self%M(j)=(0.d0,0.d0)
      END IF
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine cmat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine cmat_delete(self)
class(oft_native_cmatrix), intent(inout) :: self
integer(i4) :: i,j
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ASSOCIATED(self%M))deallocate(self%M)
if(ASSOCIATED(self%Md))deallocate(self%Md)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
IF(ASSOCIATED(self%kr))DEALLOCATE(self%kr)
IF(ASSOCIATED(self%lc))DEALLOCATE(self%lc)
IF(ASSOCIATED(self%lc_small))DEALLOCATE(self%lc_small)
!---
if(ASSOCIATED(self%map))THEN
  DO i=1,self%ni
    DO j=1,self%nj
      IF(ASSOCIATED(self%map(i,j)%lc_map))DEALLOCATE(self%map(i,j)%lc_map)
      IF(ASSOCIATED(self%map(i,j)%ext))DEALLOCATE(self%map(i,j)%ext)
    END DO
  END DO
  DEALLOCATE(self%map,self%i_map,self%j_map)
END IF
IF(ASSOCIATED(self%linkage))THEN
  IF(ASSOCIATED(self%linkage%send))THEN
    DO i=0,self%linkage%nproc_con
      IF(ASSOCIATED(self%linkage%send(i)%v))DEALLOCATE(self%linkage%send(i)%v,self%linkage%recv(i)%v)
      IF(ASSOCIATED(self%linkage%csend(i)%v))DEALLOCATE(self%linkage%csend(i)%v,self%linkage%crecv(i)%v)
    END DO
    DEALLOCATE(self%linkage%send,self%linkage%recv)
    DEALLOCATE(self%linkage%csend,self%linkage%crecv)
  END IF
  IF(ASSOCIATED(self%linkage%kle))DEALLOCATE(self%linkage%kle)
  IF(ASSOCIATED(self%linkage%lbe))DEALLOCATE(self%linkage%lbe)
  IF(ASSOCIATED(self%linkage%lle))DEALLOCATE(self%linkage%lle)
  DEALLOCATE(self%linkage)
END IF
IF(ASSOCIATED(self%Mfull))DEALLOCATE(self%Mfull)
IF(ASSOCIATED(self%redind))DEALLOCATE(self%redind)
IF(ASSOCIATED(self%color))DEALLOCATE(self%color)
self%nred=0
DEBUG_STACK_POP
end subroutine cmat_delete
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_apply_vec(self,a,b)
class(oft_native_dense_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,i1,i2
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
#ifdef BLAS_THREADED
CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
#else
bvals=0.d0
!$omp parallel private(i1,i2,i,j)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nr,i1,i2)
DO j=1,self%nc
  !$omp simd
  DO i=i1,i2
    bvals(i)=bvals(i)+self%M(i,j)*avals(j)
  END DO
END DO
!$omp end parallel
#endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_mat_apply_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_apply_cvec(self,a,b)
class(oft_native_dense_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,i1,i2
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(i1,i2,i,j)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nr,i1,i2)
DO j=1,self%nc
  !$omp simd
  DO i=i1,i2
    bvals(i)=bvals(i)+self%M(i,j)*avals(j)
  END DO
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_mat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_applyt_vec(self,a,b)
class(oft_native_dense_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,j1,j2
real(r8) :: tmp
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
#ifdef BLAS_THREADED
CALL dgemv('T',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
#else
bvals=0.d0
!$omp parallel private(j1,j2,i,j,tmp)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nc,j1,j2)
DO j=j1,j2
  tmp=0.d0
  !$omp simd reduction(+:tmp)
  DO i=1,self%nr
    tmp=tmp+self%M(i,j)*avals(i)
  END DO
  bvals(j)=tmp
END DO
!$omp end parallel
#endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_mat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_applyt_cvec(self,a,b)
class(oft_native_dense_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,j1,j2
COMPLEX(c8) :: tmp
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('T',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(j1,j2,i,j,tmp)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nc,j1,j2)
DO j=j1,j2
  tmp=0.d0
  !$omp simd reduction(+:tmp)
  DO i=1,self%nr
    tmp=tmp+self%M(i,j)*avals(i)
  END DO
  bvals(j)=tmp
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_mat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine dense_mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_native_dense_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      self%M(ii,jj)=b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'mat_set_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  do i=1,n
    do j=1,m
      self%M(i_inds(i),j_inds(j))=b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_mat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine dense_mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_dense_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      self%M(ii,jj)=self%M(ii,jj)+b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'mat_set_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_set_values',__FILE__)
  do i=1,n
    do j=1,m
      self%M(i_inds(i),j_inds(j))=self%M(i_inds(i),j_inds(j))+b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_mat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine dense_mat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_dense_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph', &
    'dense_mat_atomic_add_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph', &
    'dense_mat_atomic_add_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      !$omp atomic
      self%M(ii,jj)=self%M(ii,jj)+b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'dense_mat_atomic_add_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_mat_atomic_add_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_mat_atomic_add_values',__FILE__)
  do i=1,n
    do j=1,m
      !$omp atomic
      self%M(i_inds(i),j_inds(j))=self%M(i_inds(i),j_inds(j))+b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_mat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine dense_mat_assemble(self,diag)
class(oft_native_dense_matrix), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
integer(i4) :: i
DEBUG_STACK_PUSH
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_vector)
    !---Get diagonal matrix values
    !$omp parallel do
    do i=1,self%nr ! loop over dof
      this%v(i)=self%M(i,i)
    enddo
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a native vector.','dense_mat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add(0.d0,1.d0,self%D)
end if
DEBUG_STACK_POP
end subroutine dense_mat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine dense_mat_zero(self)
class(oft_native_dense_matrix), intent(inout) :: self
DEBUG_STACK_PUSH
self%M=0.d0
DEBUG_STACK_POP
end subroutine dense_mat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine dense_mat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_native_dense_matrix), intent(inout) :: self
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
integer(i4) :: i,j,ir
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(iblock))THEN
  !$omp parallel do private(ir,j)
  do i=1,nrows
    ir=self%i_map(iblock)%offset+irows(i)
    do j=1,self%nc
      IF(zero_diag.OR.(ir/=j))self%M(ir,j)=0.d0
    end do
  end do
ELSE
  !$omp parallel do private(ir,j)
  DO i=1,nrows
    ir=irows(i)
    DO j=1,self%nc
      IF(zero_diag.OR.(ir/=j))self%M(ir,j)=0.d0
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine dense_mat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine dense_mat_zero_cols(self,ncols,icols,jblock,keep_diag)
class(oft_native_dense_matrix), intent(inout) :: self
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
!---
integer(i4) :: i,j,ic
logical :: zero_diag
!---
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(jblock))THEN
  !$omp parallel do private(ic,j)
  do i=1,ncols
    ic=self%j_map(jblock)%offset+icols(i)
    do j=1,self%nr
      IF(zero_diag.OR.(ic/=j))self%M(j,ic)=0.d0
    end do
  end do
ELSE
  !$omp parallel do private(ic,j)
  DO i=1,ncols
    ic=icols(i)
    DO j=1,self%nr
      IF(zero_diag.OR.(ic/=j))self%M(j,ic)=0.d0
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine dense_mat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine dense_mat_delete(self)
class(oft_native_dense_matrix), intent(inout) :: self
DEBUG_STACK_PUSH
if(ASSOCIATED(self%M))deallocate(self%M)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
DEBUG_STACK_POP
end subroutine dense_mat_delete
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_apply_vec(self,a,b)
class(oft_native_dense_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,i1,i2
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(i1,i2,i,j)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nr,i1,i2)
DO j=1,self%nc
  !$omp simd
  DO i=i1,i2
    bvals(i)=bvals(i)+self%M(i,j)*avals(j)
  END DO
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_cmat_apply_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_apply_cvec(self,a,b)
class(oft_native_dense_cmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,i1,i2
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(i1,i2,i,j)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nr,i1,i2)
DO j=1,self%nc
  !$omp simd
  DO i=i1,i2
    bvals(i)=bvals(i)+self%M(i,j)*avals(j)
  END DO
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_cmat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_applyt_vec(self,a,b)
class(oft_native_dense_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !<  Result of matrix product
integer(i4) :: i,j,j1,j2
COMPLEX(c8) :: tmp
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(j1,j2,i,j,tmp)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nc,j1,j2)
DO j=j1,j2
  tmp=0.d0
  !$omp simd reduction(+:tmp)
  DO i=1,self%nr
    tmp=tmp+self%M(i,j)*avals(i)
  END DO
  bvals(j)=tmp
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_cmat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_applyt_cvec(self,a,b)
class(oft_native_dense_cmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
integer(i4) :: i,j,j1,j2
COMPLEX(c8) :: tmp
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
!---
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! #ifdef BLAS_THREADED
! CALL dgemv('N',self%nr,self%nc,1.d0,self%M,self%nr,avals,1,0.d0,bvals,1)
! #else
bvals=(0.d0,0.d0)
!$omp parallel private(j1,j2,i,j,tmp)
CALL oft_thread_slice(oft_tid,oft_env%nthreads,self%nc,j1,j2)
DO j=j1,j2
  tmp=0.d0
  !$omp simd reduction(+:tmp)
  DO i=1,self%nr
    tmp=tmp+self%M(i,j)*avals(i)
  END DO
  bvals(j)=tmp
END DO
!$omp end parallel
! #endif
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine dense_cmat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine dense_cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_native_dense_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      self%M(ii,jj)=b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'dense_cmat_set_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
  do i=1,n
    do j=1,m
      self%M(i_inds(i),j_inds(j))=b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine dense_cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_dense_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      self%M(ii,jj)=self%M(ii,jj)+b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'dense_cmat_add_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
  do i=1,n
    do j=1,m
      self%M(i_inds(i),j_inds(j))=self%M(i_inds(i),j_inds(j))+b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine dense_cmat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_dense_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: i,j,ii,jj
DEBUG_STACK_PUSH
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  IF(ANY(self%i_map(iblock)%offset+i_inds>self%nr))CALL oft_abort('Index exceeds graph', &
    'dense_cmat_atomic_add_values',__FILE__)
  IF(ANY(self%j_map(jblock)%offset+j_inds>self%nc))CALL oft_abort('Index exceeds graph', &
    'dense_cmat_atomic_add_values',__FILE__)
  do i=1,n
    ii=self%i_map(iblock)%offset+i_inds(i)
    do j=1,m
      jj=self%j_map(jblock)%offset+j_inds(j)
      !$omp atomic
      self%M(ii,jj)=self%M(ii,jj)+b(i,j)
    end do
  end do
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
                                                        'dense_cmat_atomic_add_values',__FILE__)
  IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_atomic_add_values',__FILE__)
  IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_atomic_add_values',__FILE__)
  do i=1,n
    do j=1,m
      !$omp atomic
      self%M(i_inds(i),j_inds(j))=self%M(i_inds(i),j_inds(j))+b(i,j)
    end do
  end do
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine dense_cmat_assemble(self,diag)
class(oft_native_dense_cmatrix), intent(inout) :: self
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
integer(i4) :: i
DEBUG_STACK_PUSH
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_cvector)
    !---Get diagonal matrix values
    !$omp parallel do
    do i=1,self%nr ! loop over dof
      this%v(i)=self%M(i,i)
    enddo
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a native vector.','dense_cmat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add((0.d0,0.d0),(1.d0,0.d0),self%D)
end if
DEBUG_STACK_POP
end subroutine dense_cmat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine dense_cmat_zero(self)
class(oft_native_dense_cmatrix), intent(inout) :: self
DEBUG_STACK_PUSH
self%M=(0.d0,0.d0)
DEBUG_STACK_POP
end subroutine dense_cmat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine dense_cmat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_native_dense_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
integer(i4) :: i,j,ir
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(iblock))THEN
  !$omp parallel do private(ir,j)
  do i=1,nrows
    ir=self%i_map(iblock)%offset+irows(i)
    do j=1,self%nc
      IF(zero_diag.OR.(ir/=j))self%M(ir,j)=(0.d0,0.d0)
    end do
  end do
ELSE
  !$omp parallel do private(ir,j)
  DO i=1,nrows
    ir=irows(i)
    DO j=1,self%nc
      IF(zero_diag.OR.(ir/=j))self%M(ir,j)=(0.d0,0.d0)
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine dense_cmat_zero_cols(self,ncols,icols,jblock,keep_diag)
class(oft_native_dense_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
!---
integer(i4) :: i,j,ic
logical :: zero_diag
!---
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
!---
IF(PRESENT(jblock))THEN
  !$omp parallel do private(ic,j)
  do i=1,ncols
    ic=self%j_map(jblock)%offset+icols(i)
    do j=1,self%nr
      IF(zero_diag.OR.(ic/=j))self%M(j,ic)=(0.d0,0.d0)
    end do
  end do
ELSE
  !$omp parallel do private(ic,j)
  DO i=1,ncols
    ic=icols(i)
    DO j=1,self%nr
      IF(zero_diag.OR.(ic/=j))self%M(j,ic)=(0.d0,0.d0)
    END DO
  END DO
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine dense_cmat_delete(self)
class(oft_native_dense_cmatrix), intent(inout) :: self
DEBUG_STACK_PUSH
if(ASSOCIATED(self%M))deallocate(self%M)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
DEBUG_STACK_POP
end subroutine dense_cmat_delete
end module oft_native_la
