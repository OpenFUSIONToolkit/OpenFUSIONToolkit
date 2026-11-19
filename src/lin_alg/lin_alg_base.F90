!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_la_base.F90
!
!> @defgroup doxy_oft_lin_alg Linera Algebra
!! Linear Algebra framework for the Open FUSION Toolkit
!
!> Abstract vector and matrix interfaces
!!
!! @sa oft_native_la
!! @sa oft_petsc_la
!!
!! @authors Chris Hansen
!! @date Feburary 2012
!! @ingroup doxy_oft_lin_alg
!---------------------------------------------------------------------------------
MODULE oft_la_base
USE oft_local
USE oft_base
USE oft_stitching, ONLY: oft_seam
USE oft_scatter, ONLY: oft_scatter_plan
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> @name Vector access modes
!!
!! Access modes for borrowing the local values of a vector without a copy
!! (`borrow_local`/`release_local` on @ref oft_vector and @ref oft_cvector). The same
!! mode must be passed to the borrow and the matching release.
!!
!! The modes are bit flags (`OFT_ACCESS_READWRITE = IOR(OFT_ACCESS_READ,OFT_ACCESS_WRITE)`).
!! Besides documenting intent they say which copy of the data becomes stale for a
!! backend that keeps values in more than one memory space (eg. host and device):
!! - `OFT_ACCESS_READ`: current values are needed on borrow; nothing becomes stale on release.
!! - `OFT_ACCESS_WRITE`: values on borrow are undefined (no transfer is needed before the
!!   borrow) and the caller must assign every entry of the borrowed range; after release
!!   the borrowed copy is authoritative and any other copy is stale.
!! - `OFT_ACCESS_READWRITE`: current values are needed on borrow; after release the
!!   borrowed copy is authoritative and any other copy is stale.
!! - `OFT_ACCESS_ACCUMULATE`: the borrowed array is zeroed accumulation storage, not the
!!   vector's values (no transfer is needed before the borrow). Contributions written
!!   into it are summed across processes (`stitch(1)`) and added to the vector when the
!!   accumulation is completed by a release with `sync=.TRUE.`, for operators that
!!   assemble a result element by element (eg. matrix-free applies). Several blocks may
!!   be borrowed in one accumulation; all but the last are released with
!!   `sync=.FALSE.` and the summation is done once, on the last release.
!! @{
INTEGER(i4), PARAMETER, PUBLIC :: OFT_ACCESS_READ = 1 !< Read-only access
INTEGER(i4), PARAMETER, PUBLIC :: OFT_ACCESS_WRITE = 2 !< Write-only access (all entries must be set)
INTEGER(i4), PARAMETER, PUBLIC :: OFT_ACCESS_READWRITE = 3 !< Read-write access
INTEGER(i4), PARAMETER, PUBLIC :: OFT_ACCESS_ACCUMULATE = 4 !< Accumulate contributions (added on completion)
!! @}
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
type, public :: oft_local_mat
  real(r8), pointer, contiguous, dimension(:,:) :: m => NULL() !< Needs docs
  integer(i4), pointer, contiguous, dimension(:,:) :: ind => NULL() !< Needs docs
end type oft_local_mat
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
type, public :: oft_local_cmat
  real(r8), pointer, contiguous, dimension(:,:) :: mreal => NULL() !< Needs docs
  real(r8), pointer, contiguous, dimension(:,:) :: mcomp => NULL() !< Needs docs
  integer(i4), pointer, contiguous, dimension(:,:) :: ind => NULL() !< Needs docs
end type oft_local_cmat
!---------------------------------------------------------------------------------
!> Block mapping structure
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_map
  LOGICAL :: local=.FALSE. !< Map corresponds to local vector?
  LOGICAL :: per=.FALSE. !< Map has periodic elements?
  INTEGER(i4) :: n = 0 !< Local size of vector
  INTEGER(i4) :: nslice = -1 !< Number of local elements
  INTEGER(i4) :: offset = 0 !< Offset within full vector
  INTEGER(i4) :: soffset = 0 !< Offset of local slice within full vector
  INTEGER(i8) :: ng = 0 !< Local size of vector
  INTEGER(i8) :: offsetg = 0 !< Global offset within full vector
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: gbe => NULL() !< Global boundary element flag
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:)  :: slice => NULL() !< List of boundary elements in full vector
  INTEGER(i8), POINTER, CONTIGUOUS, DIMENSION(:)  :: lge => NULL() !< List of global indices
  !> Local scatter plan gathering the owned entries (`slice`) of this block (built on
  !! first use by native vectors, see `get_slice`)
  TYPE(oft_scatter_plan) :: slice_plan
END TYPE oft_map
!---------------------------------------------------------------------------------
!> Block mapping pointer
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: map_list
  TYPE(oft_map), POINTER :: m => NULL() !< Needs docs
END TYPE map_list
!---------------------------------------------------------------------------------
!> Abstract vector class
!!
!! Basic class for OFT vectors
!---------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_vector
  INTEGER(i4) :: n = -1 !< Local size of vector
  INTEGER(i4) :: nslice = -1 !< Number of locally-owned values
  INTEGER(i4) :: nblocks = 1 !< Number of blocks
  INTEGER(i8) :: ng = -1 !< Gobal size of vector
  TYPE(oft_map), POINTER, DIMENSION(:) :: map => NULL() !< Block mapping
  TYPE(oft_seam), POINTER :: stitch_info => NULL() !< Seam information
CONTAINS
  procedure(vec_new_vec), deferred :: new_real
  procedure(vec_new_cvec), deferred :: new_complex
  !> Create a new vector as a bare copy
  generic :: new => new_real, new_complex
  !> Set all elements to a scalar
  procedure(vec_set), deferred :: set
  !> Get local portion of vector
  procedure(vec_get_local), deferred :: get_local
  !> Restore local portion of vector
  procedure(vec_restore_local), deferred :: restore_local
  !> Get locally-owned portion of vector
  procedure(vec_get_slice), deferred :: get_slice
  !> Restore locally-owned portion of vector
  procedure(vec_restore_slice), deferred :: restore_slice
  !> Add vectors
  procedure(vec_add), deferred :: add
  !> Multiply vectors element by element
  procedure(vec_mult), deferred :: mult
  procedure(vec_dot_vec), deferred :: dot_real
  procedure(vec_dot_cvec), deferred :: dot_complex
  !> Dot product with a second vector
  generic :: dot => dot_real, dot_complex
  procedure(vec_mdot_vec), deferred :: mdot_real
  procedure(vec_mdot_cvec), deferred :: mdot_complex
  !> Dot product with a set of vectors
  generic :: mdot => mdot_real, mdot_complex
  !> Scale vector by a scalar
  procedure(vec_scale), deferred :: scale
  !> Sum reduction over vector
  procedure(vec_sum), deferred :: sum
  !> Norm of vector
  procedure(vec_norm), deferred :: norm
  !> Perform global sitching
  procedure(vec_stitch), deferred :: stitch
  !> Borrow local values (no copy where supported), see @ref vector_borrow_local
  procedure :: borrow_local => vector_borrow_local
  !> Release local values obtained with `borrow_local`, see @ref vector_release_local
  procedure :: release_local => vector_release_local
  !> Copy values from a second vector
  procedure :: copy => vector_copy
  !> Add a linear combination of vectors
  procedure :: maxpy => vector_maxpy
  !> Delete vector
  procedure :: delete => vector_delete
END TYPE oft_vector
!---------------------------------------------------------------------------------
!> Vector pointer
!---------------------------------------------------------------------------------
type, public :: oft_vector_ptr
  class(oft_vector), pointer :: f => NULL() !< Needs docs
end type oft_vector_ptr
!---------------------------------------------------------------------------------
!> Abstract complex vector class
!!
!! Basic class for complex OFT vector
!---------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_cvector
  INTEGER(i4) :: n = -1 !< Local size of vector
  INTEGER(i4) :: nslice = -1 !< Number of locally-owned values
  INTEGER(i4) :: nblocks = 1 !< Number of blocks
  INTEGER(i8) :: ng = -1 !< Gobal size of vector
  TYPE(oft_map), POINTER, DIMENSION(:) :: map => NULL() !< Block mapping
  TYPE(oft_seam), POINTER :: stitch_info => NULL() !< Seam information
CONTAINS
  procedure(cvec_new_vec), deferred :: new_real
  procedure(cvec_new_cvec), deferred :: new_complex
  !> Create a new vector as a bare copy
  generic :: new => new_real, new_complex
  !> Set all elements to a scalar
  procedure(cvec_set), deferred :: set
  !> Get local portion of vector
  procedure(cvec_get_local), deferred :: get_local
  !> Restore local portion of vector
  procedure(cvec_restore_local), deferred :: restore_local
  !> Get locally-owned portion of vector
  procedure(cvec_get_slice), deferred :: get_slice
  !> Restore locally-owned portion of vector
  procedure(cvec_restore_slice), deferred :: restore_slice
  procedure(cvec_add_vec), deferred :: add_real
  procedure(cvec_add_cvec), deferred :: add_complex
  !> Add vectors
  generic :: add => add_real, add_complex
  procedure(cvec_mult_vec), deferred :: mult_real
  procedure(cvec_mult_cvec), deferred :: mult_complex
  !> Multiply vectors element by element
  generic :: mult => mult_real, mult_complex
  procedure(cvec_dot_vec), deferred :: dot_real
  procedure(cvec_dot_cvec), deferred :: dot_complex
  !> Dot product with a second vector
  generic :: dot => dot_real, dot_complex
  procedure(cvec_mdot_vec), deferred :: mdot_real
  procedure(cvec_mdot_cvec), deferred :: mdot_complex
  !> Dot product with a set of vectors
  generic :: mdot => mdot_real, mdot_complex
  !> Scale vector by a scalar
  procedure(cvec_scale), deferred :: scale
  !> Sum reduction over vector
  procedure(cvec_sum), deferred :: sum
  !> Norm of vector
  procedure(cvec_norm), deferred :: norm
  !> Perform global sitching
  procedure(cvec_stitch), deferred :: stitch
  !> Borrow local values (no copy where supported), see @ref cvector_borrow_local
  procedure :: borrow_local => cvector_borrow_local
  !> Release local values obtained with `borrow_local`, see @ref cvector_release_local
  procedure :: release_local => cvector_release_local
  procedure :: copy_real => cvector_copy_vec
  procedure :: copy_complex => cvector_copy_cvec
  !> Copy values from a second vector
  generic :: copy => copy_real, copy_complex
  !> Add a linear combination of vectors
  procedure :: maxpy => cvector_maxpy
  !> Delete vector
  procedure :: delete => cvector_delete
END TYPE oft_cvector
!---------------------------------------------------------------------------------
!> Vector pointer (complex)
!---------------------------------------------------------------------------------
type, public :: oft_cvector_ptr
  class(oft_cvector), pointer :: f => NULL() !< Needs docs
end type oft_cvector_ptr
!--- Real interfaces
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------------
  !> Create a new vector as a bare copy of an existing vector
  !---------------------------------------------------------------------------------
  subroutine vec_new_vec(self,new)
  import oft_vector
  class(oft_vector), intent(in) :: self !< Vector object
  class(oft_vector), pointer, intent(out) :: new !< New vector
  end subroutine vec_new_vec
  !---------------------------------------------------------------------------------
  !> Create a new vector as a bare copy of an existing complex vector
  !---------------------------------------------------------------------------------
  subroutine vec_new_cvec(self,new)
  import oft_vector, oft_cvector
  class(oft_vector), intent(in) :: self !< Vector object
  class(oft_cvector), pointer, intent(out) :: new !< New vector
  end subroutine vec_new_cvec
  !---------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !---------------------------------------------------------------------------------
  subroutine vec_set(self,alpha,iblock,random)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), intent(in) :: alpha !< Updated vector value
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
  logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
  end subroutine vec_set
  !---------------------------------------------------------------------------------
  !> Get local values from vector
  !---------------------------------------------------------------------------------
  subroutine vec_get_local(self,array,iblock)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), pointer, intent(inout) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
  end subroutine vec_get_local
  !---------------------------------------------------------------------------------
  !> Set/add local values to vector
  !---------------------------------------------------------------------------------
  subroutine vec_restore_local(self,array,iblock,add,wait)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), intent(in) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
  logical, optional, intent(in) :: add !< Add values instead of replace
  logical, optional, intent(in) :: wait !< Wait to perform global sync
  end subroutine vec_restore_local
  !---------------------------------------------------------------------------------
  !> Get values for locally-owned portion of vector (slice)
  !---------------------------------------------------------------------------------
  subroutine vec_get_slice(self,array,iblock)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), pointer, intent(inout) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
  end subroutine vec_get_slice
  !---------------------------------------------------------------------------------
  !> Set/add values for locally-owned portion of vector (slice)
  !---------------------------------------------------------------------------------
  subroutine vec_restore_slice(self,array,iblock,wait)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), intent(in) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
  logical, optional, intent(in) :: wait !< Wait to perform global sync?
  end subroutine vec_restore_slice
  !---------------------------------------------------------------------------------
  !> Add vectors
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !---------------------------------------------------------------------------------
  subroutine vec_add(self,gamma,alpha,a,beta,b)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), intent(in) :: gamma !< Scale of source vector
  real(r8), intent(in) :: alpha !< Scale of first vector
  class(oft_vector), target, intent(inout) :: a !< First vector to add
  real(r8), optional, intent(in) :: beta !< Scale of second vector (optional)
  class(oft_vector), target, optional, intent(inout) :: b !< Second vector to add (optional)
  end subroutine vec_add
  !---------------------------------------------------------------------------------
  !> Multiply vectors element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !---------------------------------------------------------------------------------
  subroutine vec_mult(self,a,div_flag)
  import oft_vector
  class(oft_vector), intent(inout) :: self !< Vector object
  class(oft_vector), target, intent(inout) :: a !< vector for multiplication
  logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
  end subroutine vec_mult
  !---------------------------------------------------------------------------------
  !> Scale vector by a scalar
  !---------------------------------------------------------------------------------
  subroutine vec_scale(self,alpha)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8), intent(in) :: alpha !< Scale factor
  end subroutine vec_scale
  !---------------------------------------------------------------------------------
  !> Dot product with a second vector
  !---------------------------------------------------------------------------------
  function vec_dot_vec(self,a) result(dot)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
  real(r8) :: dot !< dot(self,a)
  end function vec_dot_vec
  !---------------------------------------------------------------------------------
  !> Dot product with a second vector (complex)
  !---------------------------------------------------------------------------------
  function vec_dot_cvec(self,a) result(dot)
  import oft_vector, oft_cvector, c8
  class(oft_vector), intent(inout) :: self !< Vector object
  class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
  complex(c8) :: dot !< dot(self,a)
  end function vec_dot_cvec
  !---------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !---------------------------------------------------------------------------------
  function vec_mdot_vec(self,a,n) result(dots)
  import oft_vector, oft_vector_ptr, r8, i4
  class(oft_vector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: n !< Length of vector array
  type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  real(r8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
  end function vec_mdot_vec
  !---------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors (complex)
  !---------------------------------------------------------------------------------
  function vec_mdot_cvec(self,a,n) result(dots)
  import oft_vector, oft_cvector_ptr, c8, i4
  class(oft_vector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: n !< Length of vector array
  type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
  end function vec_mdot_cvec
  !---------------------------------------------------------------------------------
  !> Sum reduction over a vector
  !---------------------------------------------------------------------------------
  function vec_sum(self) result(sum)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  real(r8) :: sum !< Sum of vector elements
  end function vec_sum
  !---------------------------------------------------------------------------------
  !> Norm of a vector
  !---------------------------------------------------------------------------------
  function vec_norm(self,itype) result(norm)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
  real(r8) :: norm !< Specified norm of vector
  end function vec_norm
  !---------------------------------------------------------------------------------
  !> Perform global stitching
  !---------------------------------------------------------------------------------
  subroutine vec_stitch(self,up_method)
  import oft_vector, i4
  class(oft_vector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: up_method !< Type of stitching to perform
  end subroutine vec_stitch
END INTERFACE
!--- Complex interfaces
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------------
  !> Create a new vector as a bare copy of an existing real vector
  !---------------------------------------------------------------------------------
  subroutine cvec_new_vec(self,new)
  import oft_cvector, oft_vector
  class(oft_cvector), intent(in) :: self !< Vector object
  class(oft_vector), pointer, intent(out) :: new !< New vector
  end subroutine cvec_new_vec
  !---------------------------------------------------------------------------------
  !> Create a new vector as a bare copy of an existing vector
  !---------------------------------------------------------------------------------
  subroutine cvec_new_cvec(self,new)
  import oft_cvector
  class(oft_cvector), intent(in) :: self !< Vector object
  class(oft_cvector), pointer, intent(out) :: new !< New vector
  end subroutine cvec_new_cvec
  !---------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !---------------------------------------------------------------------------------
  subroutine cvec_set(self,alpha,iblock,random)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: alpha !< Updated vector value
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
  logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
  end subroutine cvec_set
  !---------------------------------------------------------------------------------
  !> Get local values from vector
  !---------------------------------------------------------------------------------
  subroutine cvec_get_local(self,array,iblock)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), pointer, intent(inout) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
  end subroutine cvec_get_local
  !---------------------------------------------------------------------------------
  !> Set/add local values to vector
  !---------------------------------------------------------------------------------
  subroutine cvec_restore_local(self,array,iblock,add,wait)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
  logical, optional, intent(in) :: add !< Add values instead of replace
  logical, optional, intent(in) :: wait !< Wait to perform global sync
  end subroutine cvec_restore_local
  !---------------------------------------------------------------------------------
  !> Get values for locally-owned portion of vector (slice)
  !---------------------------------------------------------------------------------
  subroutine cvec_get_slice(self,array,iblock)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), pointer, intent(inout) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
  end subroutine cvec_get_slice
  !---------------------------------------------------------------------------------
  !> Set/add values for locally-owned portion of vector (slice)
  !---------------------------------------------------------------------------------
  subroutine cvec_restore_slice(self,array,iblock,wait)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
  logical, optional, intent(in) :: wait !< Wait to perform global sync
  end subroutine cvec_restore_slice
  !---------------------------------------------------------------------------------
  !> Add vectors
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !---------------------------------------------------------------------------------
  subroutine cvec_add_vec(self,gamma,alpha,a,beta,b)
  import oft_cvector, oft_vector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: gamma !< Scale of source vector
  complex(c8), intent(in) :: alpha !< Scale of first vector
  class(oft_vector), target, intent(inout) :: a !< First vector to add
  complex(c8), optional, intent(in) :: beta !< Scale of second vector (optional)
  class(oft_vector), target, optional, intent(inout) :: b !< Second vector to add (optional)
  end subroutine cvec_add_vec
  !---------------------------------------------------------------------------------
  !> Add vectors
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !---------------------------------------------------------------------------------
  subroutine cvec_add_cvec(self,gamma,alpha,a,beta,b)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: gamma !< Scale of source vector
  complex(c8), intent(in) :: alpha !< Scale of first vector
  class(oft_cvector), target, intent(inout) :: a !< First vector to add
  complex(c8), optional, intent(in) :: beta !< Scale of second vector (optional)
  class(oft_cvector), target, optional, intent(inout) :: b !< Second vector to add (optional)
  end subroutine cvec_add_cvec
  !---------------------------------------------------------------------------------
  !> Multiply vectors element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !---------------------------------------------------------------------------------
  subroutine cvec_mult_vec(self,a,div_flag)
  import oft_cvector, oft_vector
  class(oft_cvector), intent(inout) :: self !< Vector object
  class(oft_vector), target, intent(inout) :: a !< vector for multiplication
  logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
  end subroutine cvec_mult_vec
  !---------------------------------------------------------------------------------
  !> Multiply vectors element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !---------------------------------------------------------------------------------
  subroutine cvec_mult_cvec(self,a,div_flag)
  import oft_cvector
  class(oft_cvector), intent(inout) :: self !< Vector object
  class(oft_cvector), target, intent(inout) :: a !< vector for multiplication
  logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
  end subroutine cvec_mult_cvec
  !---------------------------------------------------------------------------------
  !> Scale vector by a scalar
  !---------------------------------------------------------------------------------
  subroutine cvec_scale(self,alpha)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8), intent(in) :: alpha !< Scale factor
  end subroutine cvec_scale
  !---------------------------------------------------------------------------------
  !> Dot product with a second vector (real)
  !---------------------------------------------------------------------------------
  function cvec_dot_vec(self,a) result(dot)
  import oft_cvector, oft_vector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
  complex(c8) :: dot
  end function cvec_dot_vec
  !---------------------------------------------------------------------------------
  !> Dot product with a second vector
  !---------------------------------------------------------------------------------
  function cvec_dot_cvec(self,a) result(dot)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
  complex(c8) :: dot
  end function cvec_dot_cvec
  !---------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors (real)
  !---------------------------------------------------------------------------------
  function cvec_mdot_vec(self,a,n) result(dots)
  import oft_cvector, oft_vector_ptr, c8, i4
  class(oft_cvector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: n !< Length of vector array
  type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
  end function cvec_mdot_vec
  !---------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !---------------------------------------------------------------------------------
  function cvec_mdot_cvec(self,a,n) result(dots)
  import oft_cvector, oft_cvector_ptr, c8, i4
  class(oft_cvector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: n !< Length of vector array
  type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
  end function cvec_mdot_cvec
  !---------------------------------------------------------------------------------
  !> Sum reduction over a vector
  !---------------------------------------------------------------------------------
  function cvec_sum(self) result(sum)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self !< Vector object
  complex(c8) :: sum !< Sum of vector elements
  end function cvec_sum
  !---------------------------------------------------------------------------------
  !> Norm of a vector
  !---------------------------------------------------------------------------------
  function cvec_norm(self,itype) result(norm)
  import oft_cvector, i4, r8
  class(oft_cvector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
  REAL(r8) :: norm !< Specified norm of vector
  end function cvec_norm
  !---------------------------------------------------------------------------------
  !> Perform global stitching
  !---------------------------------------------------------------------------------
  subroutine cvec_stitch(self,up_method)
  import oft_cvector, i4
  class(oft_cvector), intent(inout) :: self !< Vector object
  integer(i4), intent(in) :: up_method !< Type of stitching to perform
  end subroutine cvec_stitch
END INTERFACE
!---------------------------------------------------------------------------------
!> Sub-matrix mapping
!---------------------------------------------------------------------------------
type, public :: oft_matrix_map
  integer(i4) :: nnz = 0 !< Number of non-zeros in block
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: kr => NULL() !< Row pointer to column list
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: lc => NULL() !< Column list
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: lc_map => NULL() !< List of entries in full matrix
  integer(i4), pointer, CONTIGUOUS, dimension(:,:) :: ext => NULL() !< Row extents within full matrix
end type oft_matrix_map
!---------------------------------------------------------------------------------
!> CSR graph representation
!---------------------------------------------------------------------------------
type, public :: oft_graph
  INTEGER(i4) :: nr = 0 !< Local number of rows
  INTEGER(i4) :: nc = 0 !< Local number of columns
  INTEGER(i8) :: nrg = 0 !< Gobal number of rows
  INTEGER(i8) :: ncg = 0 !< Gobal number of columns
  INTEGER(i4) :: nnz = 0 !< Number of local non-zeros
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: kr => NULL() !< Row pointer to column list
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: lc => NULL() !< Column list
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: M => NULL() !< Matrix values (optional)
end type oft_graph
!---------------------------------------------------------------------------------
!> Graph pointer
!---------------------------------------------------------------------------------
type, public :: oft_graph_ptr
  type(oft_graph), pointer :: g => NULL() !< Graph
end type oft_graph_ptr
!---------------------------------------------------------------------------------
!> Abstract matrix class
!!
!! Basic class for OFT matices (ex. fem operators)
!---------------------------------------------------------------------------------
type, abstract, public :: oft_matrix
  logical :: force_local = .FALSE. !< Do not stitch resulting vector? (Native ONLY)
  integer(i4) :: nr !< Local number of rows
  integer(i4) :: nc !< Local number of columns
  integer(i8) :: nrg !< Gobal number of rows
  integer(i8) :: ncg !< Gobal number of columns
  integer(i4) :: ni = 0 !< Number of row blocks
  integer(i4) :: nj = 0 !< Number of column blocks
  integer(i4) :: nrslice = 0 !< Number of owned rows
  integer(i4) :: ncslice = 0 !< Number of owned columns
  integer(i4), pointer, contiguous, dimension(:) :: color => NULL() !< Element coloring
  type(oft_map), pointer, dimension(:) :: i_map => NULL() !< Row block mapping
  type(oft_map), pointer, dimension(:) :: j_map => NULL() !< Column block mapping
  class(oft_vector), pointer :: D => NULL() !< Diagonal entries for scaling
contains
  !> Compute matrix-vector product
  generic :: apply => apply_real, apply_complex
  procedure(mat_apply_vec), deferred :: apply_real
  procedure :: apply_complex => matrix_apply_cvec
  !> Compute matrix-vector product for matrix transpose
  generic :: applyt => applyt_real, applyt_complex
  procedure(mat_apply_vec), deferred :: applyt_real
  procedure :: applyt_complex => matrix_applyt_cvec
  !> Set values of the matrix
  procedure(mat_set_values), deferred :: set_values
  !> Add values to the matrix
  procedure(mat_add_values), deferred :: add_values
  !> Add values atomically to the matrix
  procedure(mat_add_values), deferred :: atomic_add_values
  !> Complete matrix assembly
  procedure(mat_assemble), deferred :: assemble
  !> Zero all elements
  procedure(mat_zero), deferred :: zero
  !> Zero all elements in a given row
  procedure(mat_zero_rows), deferred :: zero_rows
  !> Save matrix elements in CSR form
  procedure :: save => matrix_save
  !> Delete matrix
  procedure :: delete => matrix_delete
end type oft_matrix
!---------------------------------------------------------------------------------
!> Matrix pointer
!---------------------------------------------------------------------------------
type, public :: oft_matrix_ptr
  class(oft_matrix), pointer :: m => NULL() !< Matrix
end type oft_matrix_ptr
!---------------------------------------------------------------------------------
!> Abstract complex matrix class
!!
!! Basic class for OFT matices (ex. fem operators)
!---------------------------------------------------------------------------------
type, abstract, public :: oft_cmatrix
  logical :: force_local = .FALSE. !< Do not stitch resulting vector (Native ONLY)
  integer(i4) :: nr !< Local number of rows
  integer(i4) :: nc !< Local number of columns
  integer(i8) :: nrg !< Gobal number of rows
  integer(i8) :: ncg !< Gobal number of columns
  integer(i4) :: ni = 0 !< Number of row blocks
  integer(i4) :: nj = 0 !< Number of column blocks
  integer(i4) :: nrslice = 0 !< Number of owned rows
  integer(i4) :: ncslice = 0 !< Number of owned columns
  type(oft_map), pointer, dimension(:) :: i_map => NULL() !< Row block mapping
  type(oft_map), pointer, dimension(:) :: j_map => NULL() !< Column block mapping
  class(oft_cvector), pointer :: D => NULL() !< Diagonal entries for scaling
contains
  !> Compute matrix-vector product
  generic :: apply => apply_real, apply_complex
  procedure :: apply_real => cmatrix_apply_vec
  procedure(cmat_apply_cvec), deferred :: apply_complex
  !> Compute matrix-vector product for matrix transpose
  generic :: applyt => applyt_real, applyt_complex
  procedure :: applyt_real => cmatrix_applyt_vec
  procedure(cmat_apply_cvec), deferred :: applyt_complex
  !> Set values of the matrix
  procedure(cmat_set_values), deferred :: set_values
  !> Add values to the matrix
  procedure(cmat_add_values), deferred :: add_values
  !> Add values atomically to the matrix
  procedure(cmat_add_values), deferred :: atomic_add_values
  !> Complete matrix assembly
  procedure(cmat_assemble), deferred :: assemble
  !> Zero all elements
  procedure(cmat_zero), deferred :: zero
  !> Zero all elements in a given row
  procedure(cmat_zero_rows), deferred :: zero_rows
  !> Delete matrix
  procedure :: delete => cmatrix_delete
end type oft_cmatrix
!---------------------------------------------------------------------------------
!> Matrix pointer
!---------------------------------------------------------------------------------
type, public :: oft_cmatrix_ptr
  class(oft_cmatrix), pointer :: m => NULL() !< Matrix
end type oft_cmatrix_ptr
!--- Real interfaces
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------------
  !> Compute matrix-vector product
  !!
  !! b = self * a
  !---------------------------------------------------------------------------------
  subroutine mat_apply_vec(self,a,b)
  import oft_vector, oft_matrix
  class(oft_matrix), intent(inout) :: self !< Matrix object
  class(oft_vector), target, intent(inout) :: a !< Source vector
  class(oft_vector), intent(inout) :: b !< Result of matrix product
  end subroutine mat_apply_vec
  !---------------------------------------------------------------------------------
  !> Compute matrix-vector product (complex vector)
  !!
  !! b = self * a
  !---------------------------------------------------------------------------------
  subroutine mat_apply_cvec(self,a,b)
  import oft_cvector, oft_matrix
  class(oft_matrix), intent(inout) :: self !< Matrix object
  class(oft_cvector), target, intent(inout) :: a !< Source vector
  class(oft_cvector), intent(inout) :: b !< Result of matrix product
  end subroutine mat_apply_cvec
  !---------------------------------------------------------------------------------
  !> Set values of a matrix
  !---------------------------------------------------------------------------------
  subroutine mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
  import oft_matrix, i4, r8
  class(oft_matrix), intent(inout) :: self !< Matrix object
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
  real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  end subroutine mat_set_values
  !---------------------------------------------------------------------------------
  !> Add values to a matrix
  !---------------------------------------------------------------------------------
  subroutine mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
  import oft_matrix, i4, r8
  class(oft_matrix), intent(inout) :: self !< Matrix object
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
  real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
  end subroutine mat_add_values
  !---------------------------------------------------------------------------------
  !> Finish assembly of matrix and optionally extract real diagonal
  !---------------------------------------------------------------------------------
  subroutine mat_assemble(self,diag)
  import oft_matrix, oft_vector
  class(oft_matrix), intent(inout) :: self !< Matrix object
  class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
  end subroutine mat_assemble
  !---------------------------------------------------------------------------------
  !> Zero all entries in matrix
  !---------------------------------------------------------------------------------
  subroutine mat_zero(self)
  import oft_matrix
  class(oft_matrix), intent(inout) :: self !< Matrix object
  end subroutine mat_zero
  !---------------------------------------------------------------------------------
  !> Zero all entries in the specified rows
  !---------------------------------------------------------------------------------
  subroutine mat_zero_rows(self,nrows,irows,iblock,keep_diag)
  import oft_matrix, i4
  class(oft_matrix), intent(inout) :: self !< Matrix object
  integer(i4), intent(in) :: nrows !< Number of rows to zero
  integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
  end subroutine mat_zero_rows
END INTERFACE
!--- Complex interfaces
ABSTRACT INTERFACE
  ! !---------------------------------------------------------------------------------
  ! !> Compute matrix-vector product (real vector)
  ! !!
  ! !! b = self * a
  ! !---------------------------------------------------------------------------------
  ! subroutine cmat_apply_vec(self,a,b)
  ! import oft_cvector, oft_vector, oft_cmatrix
  ! class(oft_cmatrix), intent(inout) :: self !< Matrix object
  ! class(oft_vector), target, intent(inout) :: a !< Source vector
  ! class(oft_cvector), intent(inout) :: b !< Result of matrix product
  ! end subroutine cmat_apply_vec
  !---------------------------------------------------------------------------------
  !> Compute matrix-vector product
  !!
  !! b = self * a
  !---------------------------------------------------------------------------------
  subroutine cmat_apply_cvec(self,a,b)
  import oft_cvector, oft_cmatrix
  class(oft_cmatrix), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a !< Source vector
  class(oft_cvector), intent(inout) :: b !< Result of matrix product
  end subroutine cmat_apply_cvec
  !---------------------------------------------------------------------------------
  !> Set values of a matrix
  !---------------------------------------------------------------------------------
  subroutine cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
  import oft_cmatrix, i4, c8
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
  complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  end subroutine cmat_set_values
  !---------------------------------------------------------------------------------
  !> Add values to a matrix
  !---------------------------------------------------------------------------------
  subroutine cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
  import oft_cmatrix, i4, c8
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
  complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
  end subroutine cmat_add_values
  !---------------------------------------------------------------------------------
  !> Finish assembly of matrix and optionally extract diagonals
  !---------------------------------------------------------------------------------
  subroutine cmat_assemble(self,diag)
  import oft_cmatrix, oft_cvector
  class(oft_cmatrix), intent(inout) :: self
  class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
  end subroutine cmat_assemble
  !---------------------------------------------------------------------------------
  !> Zero all entries in matrix
  !---------------------------------------------------------------------------------
  subroutine cmat_zero(self)
  import oft_cmatrix
  class(oft_cmatrix), intent(inout) :: self
  end subroutine cmat_zero
  !---------------------------------------------------------------------------------
  !> Zero all entries in the specified rows
  !---------------------------------------------------------------------------------
  subroutine cmat_zero_rows(self,nrows,irows,iblock,keep_diag)
  import oft_cmatrix, i4
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: nrows !< Number of rows to zero
  integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
  end subroutine cmat_zero_rows
END INTERFACE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_ml_vecspace
CONTAINS
  !> Create vector on specified level
  PROCEDURE(oft_veccreate_proto), DEFERRED :: vec_create
  !> Inject vector one level down (restriction)
  PROCEDURE(oft_inject_proto), DEFERRED :: inject
  !> Interpolate vector one level up (prolongation)
  PROCEDURE(oft_interp_proto), DEFERRED :: interp
END TYPE oft_ml_vecspace
!--- Interfaces
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Inject vector one level down (restriction)
  !------------------------------------------------------------------------------
  SUBROUTINE oft_inject_proto(self,afine,acors)
    IMPORT oft_ml_vecspace, oft_vector
    CLASS(oft_ml_vecspace), INTENT(inout) :: self !< ML vector space object
    CLASS(oft_vector), INTENT(inout) :: afine !< Source vector
    CLASS(oft_vector), INTENT(inout) :: acors !< Injected vector
  END SUBROUTINE oft_inject_proto
  !------------------------------------------------------------------------------
  !> Interpolate vector one level up (prolongation)
  !------------------------------------------------------------------------------
  SUBROUTINE oft_interp_proto(self,acors,afine)
    IMPORT oft_ml_vecspace, oft_vector
    CLASS(oft_ml_vecspace), INTENT(inout) :: self !< ML vector space object
    CLASS(oft_vector), INTENT(inout) :: acors !< Source vector
    CLASS(oft_vector), INTENT(inout) :: afine !< Interpolated vector
  END SUBROUTINE oft_interp_proto
  !------------------------------------------------------------------------------
  !> Create vector on specified level
  !------------------------------------------------------------------------------
  SUBROUTINE oft_veccreate_proto(self,new,level,cache,native)
    IMPORT oft_ml_vecspace, oft_vector, i4
    CLASS(oft_ml_vecspace), INTENT(inout) :: self !< ML vector space object
    CLASS(oft_vector), POINTER, INTENT(out) :: new !< Vector to create
    INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Vectorspace level for init (optional)
    LOGICAL, OPTIONAL, INTENT(in) :: cache !< Allow caching (optional)
    LOGICAL, OPTIONAL, INTENT(in) :: native !< Force native representation (optional)
  END SUBROUTINE oft_veccreate_proto
END INTERFACE
!---Declare public entities
public vector_extrapolate
contains
!---------------------------------------------------------------------------------
!> Borrow the local values (or a block of them) of a vector
!!
!! On return `array` refers to the local values of the vector, including values
!! shared with other processes (as for `get_local`), or of block `iblock` only.
!! Backends that store local values contiguously (native vectors) return a pointer
!! into the vector's own storage without copying; this default implementation,
!! used by other backends, returns a copy obtained with `get_local` (or a zeroed
!! array for @ref OFT_ACCESS_ACCUMULATE, without reading the vector). Callers must
!! work correctly in both cases.
!!
!! Contract:
!! - `array` is owned by the vector: do not allocate, deallocate or reassociate it,
!!   and always give it back with `release_local` using the same `mode` and `iblock`.
!! - Between borrow and release do not call any other method of the vector (or of a
!!   vector created as an alias of it) that reads or modifies its values, including
!!   `stitch`, and do not delete it. Borrowing several different vectors at once is
!!   fine, as is borrowing the same vector more than once for reading only.
!! - Writing through `array` is only allowed if `mode` includes @ref OFT_ACCESS_WRITE.
!!   With @ref OFT_ACCESS_WRITE alone the initial contents are undefined and every
!!   entry must be assigned before release.
!! - `array` is invalid after release (it is nullified).
!! - With @ref OFT_ACCESS_ACCUMULATE `array` starts as zero and holds contributions to
!!   be added to the vector (see the access modes). Between the first accumulate borrow
!!   and the completing release (`sync=.TRUE.`) no other access to the vector is
!!   allowed, and blocks not borrowed in the accumulation are left unchanged.
!---------------------------------------------------------------------------------
subroutine vector_borrow_local(self,array,mode,iblock)
class(oft_vector), intent(inout) :: self !< Vector object
real(r8), pointer, contiguous, intent(inout) :: array(:) !< Local values (output)
integer(i4), intent(in) :: mode !< Access mode (@ref OFT_ACCESS_READ, @ref OFT_ACCESS_WRITE, @ref OFT_ACCESS_READWRITE or @ref OFT_ACCESS_ACCUMULATE)
integer(i4), optional, intent(in) :: iblock !< Sub-block to borrow (optional)
real(r8), pointer, dimension(:) :: vtmp
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'vector_borrow_local',__FILE__)
NULLIFY(vtmp)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Zeroed contributions; values are not read, as get_local could overwrite contributions
  !---already held by the backend for the blocks released earlier in this accumulation
  IF(PRESENT(iblock))THEN
    IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vector_borrow_local',__FILE__)
    ALLOCATE(vtmp(self%map(iblock)%n))
  ELSE
    ALLOCATE(vtmp(self%n))
  END IF
  vtmp=0.d0
ELSE
  CALL self%get_local(vtmp,iblock)
END IF
array=>vtmp
DEBUG_STACK_POP
end subroutine vector_borrow_local
!---------------------------------------------------------------------------------
!> Release local values obtained with @ref vector_borrow_local
!!
!! If `mode` includes @ref OFT_ACCESS_WRITE the values in `array` become the local
!! values of the vector. `sync` then selects whether values shared between processes
!! are synchronized from their owning process (`stitch(0)`, as done by `restore_local`):
!! - `sync=.TRUE.` is always safe and must be used when processes sharing a value
!!   may have written different values to it.
!! - `sync=.FALSE.` skips the MPI exchange. Use it only if shared values are
!!   consistent across processes when the array is released, ie. the vector was
!!   consistent on borrow and every process sharing a value made the same
!!   modification (eg. zeroing based on a global boundary flag). Native vectors are
!!   consistent after any stitching, including summation, and check this condition
!!   in debug builds (`OFT_DEBUG`). To synchronize several blocks once, release
!!   them all with `sync=.TRUE.` or release the whole vector instead.
!! Backends that keep a single copy of each shared value (eg. PETSc) only retain
!! the owned values, which is the same result as `sync=.TRUE.`.
!!
!! For @ref OFT_ACCESS_READ `sync` must be `.FALSE.`. For @ref OFT_ACCESS_ACCUMULATE,
!! `sync=.TRUE.` completes the accumulation (contributions are summed across processes
!! with `stitch(1)` and added to the vector) and `sync=.FALSE.` defers this to the
!! release of a later block of the same accumulation. `array` is nullified on return.
!---------------------------------------------------------------------------------
subroutine vector_release_local(self,array,mode,sync,iblock)
class(oft_vector), intent(inout) :: self !< Vector object
real(r8), pointer, contiguous, intent(inout) :: array(:) !< Borrowed local values (nullified on return)
integer(i4), intent(in) :: mode !< Access mode used in the matching borrow
logical, intent(in) :: sync !< Synchronize shared values from their owner after writing?
integer(i4), optional, intent(in) :: iblock !< Sub-block used in the matching borrow (optional)
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'vector_release_local',__FILE__)
IF(sync.AND.(mode==OFT_ACCESS_READ))CALL oft_abort('Synchronization requested for read-only access', &
  'vector_release_local',__FILE__)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Summed into a full-size cache by restore_local, stitched and added on the last block
  IF(PRESENT(iblock))THEN
    CALL self%restore_local(array,iblock,add=.TRUE.,wait=(.NOT.sync))
  ELSE
    IF(.NOT.sync)CALL oft_abort('Accumulation of the full vector must be completed on release', &
      'vector_release_local',__FILE__)
    CALL self%restore_local(array,add=.TRUE.)
  END IF
ELSE IF(IAND(mode,OFT_ACCESS_WRITE)/=0)THEN
  IF(PRESENT(iblock))THEN
    CALL self%restore_local(array,iblock,wait=(.NOT.sync))
  ELSE
    CALL self%restore_local(array)
  END IF
END IF
DEALLOCATE(array)
DEBUG_STACK_POP
end subroutine vector_release_local
!---------------------------------------------------------------------------------
!> Copy values from a second vector
!!
!! self = a
!!
!! @note This default implementation uses `add`, so non-finite values in `self` may
!! propagate; backends should override it with a direct copy.
!---------------------------------------------------------------------------------
subroutine vector_copy(self,a)
class(oft_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Source vector
DEBUG_STACK_PUSH
CALL self%add(0.d0,1.d0,a)
DEBUG_STACK_POP
end subroutine vector_copy
!---------------------------------------------------------------------------------
!> Add a linear combination of vectors
!!
!! self = self + \f$ \sum_k \alpha_k a_k \f$
!!
!! Each entry is updated by adding the terms in order of increasing `k`, so the result
!! is the same as `n` successive calls to `add`. This default implementation does
!! exactly that.
!---------------------------------------------------------------------------------
subroutine vector_maxpy(self,alpha,a,n)
class(oft_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: n !< Number of vectors
real(r8), intent(in) :: alpha(n) !< Scale factors [n]
type(oft_vector_ptr), intent(inout) :: a(n) !< Vectors to add [n]
integer(i4) :: k
DEBUG_STACK_PUSH
DO k=1,n
  CALL self%add(1.d0,alpha(k),a(k)%f)
END DO
DEBUG_STACK_POP
end subroutine vector_maxpy
!---------------------------------------------------------------------------------
!> Borrow the local values (or a block of them) of a complex vector
!!
!! See @ref vector_borrow_local for the contract.
!---------------------------------------------------------------------------------
subroutine cvector_borrow_local(self,array,mode,iblock)
class(oft_cvector), intent(inout) :: self !< Vector object
complex(c8), pointer, contiguous, intent(inout) :: array(:) !< Local values (output)
integer(i4), intent(in) :: mode !< Access mode (@ref OFT_ACCESS_READ, @ref OFT_ACCESS_WRITE, @ref OFT_ACCESS_READWRITE or @ref OFT_ACCESS_ACCUMULATE)
integer(i4), optional, intent(in) :: iblock !< Sub-block to borrow (optional)
complex(c8), pointer, dimension(:) :: vtmp
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'cvector_borrow_local',__FILE__)
NULLIFY(vtmp)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Zeroed contributions; values are not read, as get_local could overwrite contributions
  !---already held by the backend for the blocks released earlier in this accumulation
  IF(PRESENT(iblock))THEN
    IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvector_borrow_local',__FILE__)
    ALLOCATE(vtmp(self%map(iblock)%n))
  ELSE
    ALLOCATE(vtmp(self%n))
  END IF
  vtmp=(0.d0,0.d0)
ELSE
  CALL self%get_local(vtmp,iblock)
END IF
array=>vtmp
DEBUG_STACK_POP
end subroutine cvector_borrow_local
!---------------------------------------------------------------------------------
!> Release local values obtained with @ref cvector_borrow_local
!!
!! See @ref vector_release_local for the contract.
!---------------------------------------------------------------------------------
subroutine cvector_release_local(self,array,mode,sync,iblock)
class(oft_cvector), intent(inout) :: self !< Vector object
complex(c8), pointer, contiguous, intent(inout) :: array(:) !< Borrowed local values (nullified on return)
integer(i4), intent(in) :: mode !< Access mode used in the matching borrow
logical, intent(in) :: sync !< Synchronize shared values from their owner after writing?
integer(i4), optional, intent(in) :: iblock !< Sub-block used in the matching borrow (optional)
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'cvector_release_local',__FILE__)
IF(sync.AND.(mode==OFT_ACCESS_READ))CALL oft_abort('Synchronization requested for read-only access', &
  'cvector_release_local',__FILE__)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Summed into a full-size cache by restore_local, stitched and added on the last block
  IF(PRESENT(iblock))THEN
    CALL self%restore_local(array,iblock,add=.TRUE.,wait=(.NOT.sync))
  ELSE
    IF(.NOT.sync)CALL oft_abort('Accumulation of the full vector must be completed on release', &
      'cvector_release_local',__FILE__)
    CALL self%restore_local(array,add=.TRUE.)
  END IF
ELSE IF(IAND(mode,OFT_ACCESS_WRITE)/=0)THEN
  IF(PRESENT(iblock))THEN
    CALL self%restore_local(array,iblock,wait=(.NOT.sync))
  ELSE
    CALL self%restore_local(array)
  END IF
END IF
DEALLOCATE(array)
DEBUG_STACK_POP
end subroutine cvector_release_local
!---------------------------------------------------------------------------------
!> Copy values from a real vector
!!
!! self = a
!!
!! @note This default implementation uses `add`, see @ref vector_copy.
!---------------------------------------------------------------------------------
subroutine cvector_copy_vec(self,a)
class(oft_cvector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Source vector
DEBUG_STACK_PUSH
CALL self%add((0.d0,0.d0),(1.d0,0.d0),a)
DEBUG_STACK_POP
end subroutine cvector_copy_vec
!---------------------------------------------------------------------------------
!> Copy values from a second complex vector
!!
!! self = a
!!
!! @note This default implementation uses `add`, see @ref vector_copy.
!---------------------------------------------------------------------------------
subroutine cvector_copy_cvec(self,a)
class(oft_cvector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< Source vector
DEBUG_STACK_PUSH
CALL self%add((0.d0,0.d0),(1.d0,0.d0),a)
DEBUG_STACK_POP
end subroutine cvector_copy_cvec
!---------------------------------------------------------------------------------
!> Add a linear combination of complex vectors
!!
!! self = self + \f$ \sum_k \alpha_k a_k \f$
!!
!! See @ref vector_maxpy.
!---------------------------------------------------------------------------------
subroutine cvector_maxpy(self,alpha,a,n)
class(oft_cvector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: n !< Number of vectors
complex(c8), intent(in) :: alpha(n) !< Scale factors [n]
type(oft_cvector_ptr), intent(inout) :: a(n) !< Vectors to add [n]
integer(i4) :: k
DEBUG_STACK_PUSH
DO k=1,n
  CALL self%add((1.d0,0.d0),alpha(k),a(k)%f)
END DO
DEBUG_STACK_POP
end subroutine cvector_maxpy
!---------------------------------------------------------------------------------
!> Finalize vector
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized vectors
!---------------------------------------------------------------------------------
subroutine vector_delete(self)
class(oft_vector), intent(inout) :: self
call oft_warn('Finalizing general real vector, this may indicate an error.')
end subroutine vector_delete
!---------------------------------------------------------------------------------
!> Finalize vector
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized vectors
!---------------------------------------------------------------------------------
subroutine cvector_delete(self)
class(oft_cvector), intent(inout) :: self
call oft_warn('Finalizing general complex vector, this may indicate an error.')
end subroutine cvector_delete
!---------------------------------------------------------------------------------
!> Compute matrix vector product (complex)
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine matrix_apply_cvec(self,a,b)
class(oft_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: atmp,btmp
complex(c8), pointer :: avals(:)
CALL a%new(atmp)
CALL b%new(btmp)
CALL b%set((0.d0,0.d0))
nullify(avals)
CALL a%get_slice(avals)
!---Apply real part
CALL atmp%restore_slice(REAL(avals))
CALL self%apply(atmp,btmp)
CALL b%copy(btmp)
!---Apply imaginary part
CALL atmp%restore_slice(AIMAG(avals))
CALL self%apply(atmp,btmp)
CALL b%add((1.d0,0.d0),(0.d0,1.d0),btmp)
!---Cleanup
CALL atmp%delete()
CALL btmp%delete()
DEALLOCATE(atmp,btmp,avals)
end subroutine matrix_apply_cvec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose (complex vector)
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine matrix_applyt_cvec(self,a,b)
class(oft_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: atmp,btmp
complex(c8), pointer :: avals(:)
CALL a%new(atmp)
CALL b%new(btmp)
CALL b%set((0.d0,0.d0))
nullify(avals)
CALL a%get_slice(avals)
!---Apply real part
CALL atmp%restore_slice(REAL(avals,r8))
CALL self%applyt(atmp,btmp)
CALL b%copy(btmp)
!---Apply imaginary part
CALL atmp%restore_slice(AIMAG(avals))
CALL self%applyt(atmp,btmp)
CALL b%add((1.d0,0.d0),(0.d0,1.d0),btmp)
!---Cleanup
CALL atmp%delete()
CALL btmp%delete()
DEALLOCATE(atmp,btmp,avals)
end subroutine matrix_applyt_cvec
!---------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!---------------------------------------------------------------------------------
subroutine matrix_delete(self)
class(oft_matrix), intent(inout) :: self
call oft_warn('Finalizing general real matrix, this may indicate an error.')
end subroutine matrix_delete
!---------------------------------------------------------------------------------
!> Save matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!---------------------------------------------------------------------------------
subroutine matrix_save(self, filename, path, bc_flags, nfields)
class(oft_matrix), intent(inout) :: self
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
logical, intent(in), optional, dimension(:) :: bc_flags
integer, intent(in), optional :: nfields
call oft_abort('No matrix save implemented for this matrix type', 'matrix_save', __FILE__)
end subroutine matrix_save
!---------------------------------------------------------------------------------
!> Compute matrix vector product
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine cmatrix_apply_vec(self,a,b)
class(oft_cmatrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: atmp
CALL a%new(atmp)
CALL atmp%copy(a)
CALL self%apply_complex(atmp,b)
CALL atmp%delete()
DEALLOCATE(atmp)
end subroutine cmatrix_apply_vec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine cmatrix_applyt_vec(self,a,b)
class(oft_cmatrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: atmp
CALL a%new(atmp)
CALL atmp%copy(a)
CALL self%applyt_complex(atmp,b)
CALL atmp%delete()
DEALLOCATE(atmp)
end subroutine cmatrix_applyt_vec
!---------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!---------------------------------------------------------------------------------
subroutine cmatrix_delete(self)
class(oft_cmatrix), intent(inout) :: self
call oft_warn('Finalizing general complex matrix, this may indicate an error.')
end subroutine cmatrix_delete
!------------------------------------------------------------------------------
!> Extrapolate a vector using a polynomial fit to previous vectors
!------------------------------------------------------------------------------
subroutine vector_extrapolate(x,vectors,n,xe,output)
real(r8), intent(in) :: x(:) !< Array of previous positions [n]
type(oft_vector_ptr), intent(inout) :: vectors(:) !< Array of previous vectors [n]
integer(i4), intent(in) :: n !< Number of vectors to use for interpolation
real(r8), intent(in) :: xe !< Position to extrapolate to
class(oft_vector), intent(inout) :: output !< Extrapolated vector
!---
integer(i4) :: i,j
real(r8) :: y
DEBUG_STACK_PUSH
!---
CALL output%set(0.d0)
DO i=1,n
  y=1.d0
  DO j=1,n
    IF(i==j)CYCLE
    y=y*(xe-x(j))/(x(i)-x(j))
  END DO
  CALL output%add(1.d0,y,vectors(i)%f)
END DO
DEBUG_STACK_POP
end subroutine vector_extrapolate
end module oft_la_base
