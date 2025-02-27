!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_la_base.F90
!
!> @defgroup doxy_oft_lin_alg Linera Algebra
!! Linear Algebra framework for the Open FUSION Toolkit
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
!! @ingroup doxy_oft_lin_alg
!------------------------------------------------------------------------------
MODULE oft_la_base
USE oft_local
USE oft_base
USE oft_stitching, ONLY: oft_seam
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, public :: oft_local_mat
  real(r8), pointer, contiguous, dimension(:,:) :: m => NULL()
  integer(i4), pointer, contiguous, dimension(:,:) :: ind => NULL()
end type oft_local_mat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, public :: oft_local_cmat
  real(r8), pointer, contiguous, dimension(:,:) :: mreal => NULL()
  real(r8), pointer, contiguous, dimension(:,:) :: mcomp => NULL()
  integer(i4), pointer, contiguous, dimension(:,:) :: ind => NULL()
end type oft_local_cmat
!------------------------------------------------------------------------------
!> Block mapping structure
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_map
  LOGICAL :: local=.FALSE. !< Map corresponds to local field
  LOGICAL :: per=.FALSE. !< Flag indicating map has periodic elements
  INTEGER(i4) :: n = 0 !< Local size of field
  INTEGER(i4) :: nslice = -1 !< Number of local elements
  INTEGER(i4) :: offset = 0 !< Offset within full field
  INTEGER(i4) :: soffset = 0 !< Offset of local slice within full field
  INTEGER(i8) :: ng = 0 !< Local size of field
  INTEGER(i8) :: offsetg = 0 !< Global offset within full field
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: gbe => NULL() !< Global boundary element flag
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:)  :: slice => NULL() !< List of boundary elements in full vector
  INTEGER(i8), POINTER, CONTIGUOUS, DIMENSION(:)  :: lge => NULL() !< List of global indices
END TYPE oft_map
!------------------------------------------------------------------------------
!> Block mapping pointer
!------------------------------------------------------------------------------
TYPE, PUBLIC :: map_list
  TYPE(oft_map), POINTER :: m => NULL()
END TYPE map_list
!------------------------------------------------------------------------------
!> Abstract field class
!!
!! Basic class for OFT fields (ex. fem vectors)
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_vector
  INTEGER(i4) :: n = -1 !< Local size of field
  INTEGER(i4) :: nslice = -1 !< Local size of field
  INTEGER(i4) :: nblocks = 1 !< Number of blocks
  INTEGER(i8) :: ng = -1 !< Gobal size of field
  TYPE(oft_map), POINTER, DIMENSION(:) :: map => NULL() !< Block mapping
  TYPE(oft_seam), POINTER :: stitch_info => NULL() !< Seam information
CONTAINS
  procedure(vec_new_vec), deferred :: new_real
  procedure(vec_new_cvec), deferred :: new_complex
  !> Create a new field as a bare copy
  generic :: new => new_real, new_complex
  !> Set all elements to a scalar
  procedure(vec_set), deferred :: set
  !> Get local portion of field
  procedure(vec_get_local), deferred :: get_local
  !> Restore local portion of field
  procedure(vec_restore_local), deferred :: restore_local
  !> Get owned portion of field
  procedure(vec_get_slice), deferred :: get_slice
  !> Restore owned portion of field
  procedure(vec_restore_slice), deferred :: restore_slice
  !> Add fields
  procedure(vec_add), deferred :: add
  !> Multiply fields element by element
  procedure(vec_mult), deferred :: mult
  procedure(vec_dot_vec), deferred :: dot_real
  procedure(vec_dot_cvec), deferred :: dot_complex
  !> Dot product with a second field
  generic :: dot => dot_real, dot_complex
  procedure(vec_mdot_vec), deferred :: mdot_real
  procedure(vec_mdot_cvec), deferred :: mdot_complex
  !> Dot product with a set of fields
  generic :: mdot => mdot_real, mdot_complex
  !> Scale field by a scalar
  procedure(vec_scale), deferred :: scale
  !> Sum reduction over field
  procedure(vec_sum), deferred :: sum
  !> Norm of field
  procedure(vec_norm), deferred :: norm
  !> Perform global sitching
  procedure(vec_stitch), deferred :: stitch
  !> Delete field
  procedure :: delete => vector_delete
END TYPE oft_vector
!------------------------------------------------------------------------------
!> Field pointer
!------------------------------------------------------------------------------
type, public :: oft_vector_ptr
  class(oft_vector), pointer :: f => NULL()
end type oft_vector_ptr
!------------------------------------------------------------------------------
!> Abstract complex vector class
!!
!! Basic class for OFT fields (ex. fem vectors)
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_cvector
  INTEGER(i4) :: n = -1 !< Local size of field
  INTEGER(i4) :: nslice = -1 !< Local size of field
  INTEGER(i4) :: nblocks = 1 !< Number of blocks
  INTEGER(i8) :: ng = -1 !< Gobal size of field
  TYPE(oft_map), POINTER, DIMENSION(:) :: map => NULL() !< Block mapping
  TYPE(oft_seam), POINTER :: stitch_info => NULL() !< Seam information
CONTAINS
  procedure(cvec_new_vec), deferred :: new_real
  procedure(cvec_new_cvec), deferred :: new_complex
  !> Create a new field as a bare copy
  generic :: new => new_real, new_complex
  !> Set all elements to a scalar
  procedure(cvec_set), deferred :: set
  !> Get local portion of field
  procedure(cvec_get_local), deferred :: get_local
  !> Restore local portion of field
  procedure(cvec_restore_local), deferred :: restore_local
  !> Get owned portion of field
  procedure(cvec_get_slice), deferred :: get_slice
  !> Restore owned portion of field
  procedure(cvec_restore_slice), deferred :: restore_slice
  procedure(cvec_add_vec), deferred :: add_real
  procedure(cvec_add_cvec), deferred :: add_complex
  !> Add fields
  generic :: add => add_real, add_complex
  procedure(cvec_mult_vec), deferred :: mult_real
  procedure(cvec_mult_cvec), deferred :: mult_complex
  !> Multiply fields element by element
  generic :: mult => mult_real, mult_complex
  procedure(cvec_dot_vec), deferred :: dot_real
  procedure(cvec_dot_cvec), deferred :: dot_complex
  !> Dot product with a second field
  generic :: dot => dot_real, dot_complex
  procedure(cvec_mdot_vec), deferred :: mdot_real
  procedure(cvec_mdot_cvec), deferred :: mdot_complex
  !> Dot product with a set of fields
  generic :: mdot => mdot_real, mdot_complex
  !> Scale field by a scalar
  procedure(cvec_scale), deferred :: scale
  !> Sum reduction over field
  procedure(cvec_sum), deferred :: sum
  !> Norm of field
  procedure(cvec_norm), deferred :: norm
  !> Perform global sitching
  procedure(cvec_stitch), deferred :: stitch
  !> Delete field
  procedure :: delete => cvector_delete
END TYPE oft_cvector
!------------------------------------------------------------------------------
!> Field pointer
!------------------------------------------------------------------------------
type, public :: oft_cvector_ptr
  class(oft_cvector), pointer :: f => NULL()
end type oft_cvector_ptr
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Create a new field as a bare copy
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_new_vec(self,new)
  import oft_vector
  class(oft_vector), intent(in) :: self
  class(oft_vector), pointer, intent(out) :: new !< New field
  end subroutine vec_new_vec
  !------------------------------------------------------------------------------
  !> Create a new field as a bare copy
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_new_cvec(self,new)
  import oft_vector, oft_cvector
  class(oft_vector), intent(in) :: self
  class(oft_cvector), pointer, intent(out) :: new !< New field
  end subroutine vec_new_cvec
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_set(self,alpha,iblock,random)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  real(r8), intent(in) :: alpha !< Updated field value
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
  logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
  end subroutine vec_set
  !------------------------------------------------------------------------------
  !> Needs docs
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_get_local(self,array,iblock)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  real(r8), pointer, intent(inout) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block
  end subroutine vec_get_local
  !------------------------------------------------------------------------------
  !> Needs docs
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_restore_local(self,array,iblock,add,wait)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  real(r8), intent(in) :: array(:) !< Local values
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block
  logical, optional, intent(in) :: add !< Add values instead of replace
  logical, optional, intent(in) :: wait !< Wait to perform global sync
  end subroutine vec_restore_local
  !------------------------------------------------------------------------------
  !> Needs docs
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_get_slice(self,array,iblock)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  real(r8), pointer, intent(inout) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block
  end subroutine vec_get_slice
  !------------------------------------------------------------------------------
  !> Needs docs
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_restore_slice(self,array,iblock,wait)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  real(r8), intent(in) :: array(:) !< Slice values
  integer(i4), optional, intent(in) :: iblock !< Vector sub-block
  logical, optional, intent(in) :: wait !< Wait to perform global sync
  end subroutine vec_restore_slice
  !------------------------------------------------------------------------------
  !> Add fields
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_add(self,gamma,alpha,a,beta,b)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self
  real(r8), intent(in) :: gamma !< Scale of source field
  real(r8), intent(in) :: alpha !< Scale of first field
  class(oft_vector), target, intent(inout) :: a !< First field to add
  real(r8), optional, intent(in) :: beta !< Scale of second field (optional)
  class(oft_vector), target, optional, intent(inout) :: b !< Second field to add (optional)
  end subroutine vec_add
  !------------------------------------------------------------------------------
  !> Multiply fields element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_mult(self,a,div_flag)
  import oft_vector
  class(oft_vector), intent(inout) :: self
  class(oft_vector), target, intent(inout) :: a !< Field for multiplication
  logical, optional, intent(in) :: div_flag
  end subroutine vec_mult
  !------------------------------------------------------------------------------
  !> Scale field by a scalar
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_scale(self,alpha)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self
  real(r8), intent(in) :: alpha !< Scale factor
  end subroutine vec_scale
  !------------------------------------------------------------------------------
  !> Dot product with a second field (real)
  !!
  !! @result dot(self,a)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  function vec_dot_vec(self,a) result(dot)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self
  class(oft_vector), target, intent(inout) :: a !< Second field for dot product
  real(r8) :: dot
  end function vec_dot_vec
  !------------------------------------------------------------------------------
  !> Dot product with a second field (complex)
  !!
  !! @result dot(self,a)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  function vec_dot_cvec(self,a) result(dot)
  import oft_vector, oft_cvector, c8
  class(oft_vector), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a !< Second field for dot product
  complex(c8) :: dot
  end function vec_dot_cvec
  !------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !!
  !! @result \f$ \sum_i self_i a(j)_i \f$
  !------------------------------------------------------------------------------
  function vec_mdot_vec(self,a,n) result(dots)
  import oft_vector, oft_vector_ptr, r8, i4
  class(oft_vector), intent(inout) :: self
  type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  integer(i4), intent(in) :: n !< Length of vector array
  real(r8) :: dots(n)
  end function vec_mdot_vec
  !------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !!
  !! @result \f$ \sum_i self_i a(j)_i \f$
  !------------------------------------------------------------------------------
  function vec_mdot_cvec(self,a,n) result(dots)
  import oft_vector, oft_cvector_ptr, c8, i4
  class(oft_vector), intent(inout) :: self
  type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
  integer(i4), intent(in) :: n !< Length of vector array
  complex(c8) :: dots(n)
  end function vec_mdot_cvec
  !------------------------------------------------------------------------------
  !> Sum reduction over a field
  !!
  !! @result sum(self)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  function vec_sum(self) result(sum)
  import oft_vector, r8
  class(oft_vector), intent(inout) :: self
  real(r8) :: sum
  end function vec_sum
  !------------------------------------------------------------------------------
  !> Norm of a field
  !!
  !! @result norm(self)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  function vec_norm(self,itype) result(norm)
  import oft_vector, i4, r8
  class(oft_vector), intent(inout) :: self
  integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
  real(r8) :: norm
  end function vec_norm
  !------------------------------------------------------------------------------
  !> Perform global stitching
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields
  !------------------------------------------------------------------------------
  subroutine vec_stitch(self,up_method)
  import oft_vector, i4
  class(oft_vector), intent(inout) :: self
  integer(i4), intent(in) :: up_method !< Type of stitching to perform
  end subroutine vec_stitch
END INTERFACE
!
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Create a new field as a bare copy
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !!
  !! @param[out] new New field
  !------------------------------------------------------------------------------
  subroutine cvec_new_vec(self,new)
  import oft_cvector, oft_vector
  class(oft_cvector), intent(in) :: self
  class(oft_vector), pointer, intent(out) :: new
  end subroutine cvec_new_vec
  !------------------------------------------------------------------------------
  !> Create a new field as a bare copy
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !!
  !! @param[out] new New field
  !------------------------------------------------------------------------------
  subroutine cvec_new_cvec(self,new)
  import oft_cvector
  class(oft_cvector), intent(in) :: self
  class(oft_cvector), pointer, intent(out) :: new
  end subroutine cvec_new_cvec
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @param[in] alpha Updated field value
  !! @param[in] random Set to random number, if true alpha is ignored (optional)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_set(self,alpha,iblock,random)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: alpha
  integer(i4), optional, intent(in) :: iblock
  logical, optional, intent(in) :: random
  end subroutine cvec_set
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @param[in] alpha Updated field value
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_get_local(self,array,iblock)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), pointer, intent(inout) :: array(:)
  integer(i4), optional, intent(in) :: iblock
  end subroutine cvec_get_local
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @param[in] alpha Updated field value
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_restore_local(self,array,iblock,add,wait)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: array(:)
  integer(i4), optional, intent(in) :: iblock
  logical, optional, intent(in) :: add
  logical, optional, intent(in) :: wait
  end subroutine cvec_restore_local
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @param[in] alpha Updated field value
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_get_slice(self,array,iblock)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), pointer, intent(inout) :: array(:)
  integer(i4), optional, intent(in) :: iblock
  end subroutine cvec_get_slice
  !------------------------------------------------------------------------------
  !> Set all elements to a scalar
  !!
  !! @param[in] alpha Updated field value
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_restore_slice(self,array,iblock,wait)
  import oft_cvector, i4, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: array(:)
  integer(i4), optional, intent(in) :: iblock
  logical, optional, intent(in) :: wait
  end subroutine cvec_restore_slice
  !------------------------------------------------------------------------------
  !> Add fields
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !! @param[in] gamma Scale of source field
  !! @param[in] alpha Scale of first field
  !! @param[in] a First field to add
  !! @param[in] beta Scale of second field (optional)
  !! @param[in] b Second field to add (optional)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_add_vec(self,gamma,alpha,a,beta,b)
  import oft_cvector, oft_vector, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: gamma
  complex(c8), intent(in) :: alpha
  class(oft_vector), target, intent(inout) :: a
  complex(c8), optional, intent(in) :: beta
  class(oft_vector), target, optional, intent(inout) :: b
  end subroutine cvec_add_vec
  !------------------------------------------------------------------------------
  !> Add fields
  !!
  !! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
  !! @param[in] gamma Scale of source field
  !! @param[in] alpha Scale of first field
  !! @param[in] a First field to add
  !! @param[in] beta Scale of second field (optional)
  !! @param[in] b Second field to add (optional)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_add_cvec(self,gamma,alpha,a,beta,b)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: gamma
  complex(c8), intent(in) :: alpha
  class(oft_cvector), target, intent(inout) :: a
  complex(c8), optional, intent(in) :: beta
  class(oft_cvector), target, optional, intent(inout) :: b
  end subroutine cvec_add_cvec
  !------------------------------------------------------------------------------
  !> Multiply fields element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !! @param[in] a First field to add
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_mult_vec(self,a,div_flag)
  import oft_cvector, oft_vector
  class(oft_cvector), intent(inout) :: self
  class(oft_vector), target, intent(inout) :: a
  logical, optional, intent(in) :: div_flag
  end subroutine cvec_mult_vec
  !------------------------------------------------------------------------------
  !> Multiply fields element by element
  !!
  !! \f$ self_i = self_i * a_i \f$
  !! @param[in] a First field to add
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_mult_cvec(self,a,div_flag)
  import oft_cvector
  class(oft_cvector), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a
  logical, optional, intent(in) :: div_flag
  end subroutine cvec_mult_cvec
  !------------------------------------------------------------------------------
  !> Scale field by a scalar
  !!
  !! @param[in] alpha Factor to scale field
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_scale(self,alpha)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8), intent(in) :: alpha
  end subroutine cvec_scale
  !------------------------------------------------------------------------------
  !> Dot product with a second field
  !!
  !! @param[in] a Second field for dot product
  !! @result dot(self,a)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  function cvec_dot_vec(self,a) result(dot)
  import oft_cvector, oft_vector, c8
  class(oft_cvector), intent(inout) :: self
  class(oft_vector), target, intent(inout) :: a
  complex(c8) :: dot
  end function cvec_dot_vec
  !------------------------------------------------------------------------------
  !> Dot product with a second field
  !!
  !! @param[in] a Second field for dot product
  !! @result dot(self,a)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  function cvec_dot_cvec(self,a) result(dot)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a
  complex(c8) :: dot
  end function cvec_dot_cvec
  !------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !!
  !! @param[in] a Array of vectors for dot product [n]
  !! @param[in] n Length of vector array
  !! @result \f$ \sum_i self_i a(j)_i \f$
  !------------------------------------------------------------------------------
  function cvec_mdot_vec(self,a,n) result(dots)
  import oft_cvector, oft_vector_ptr, c8, i4
  class(oft_cvector), intent(inout) :: self
  type(oft_vector_ptr), intent(inout) :: a(n)
  integer(i4), intent(in) :: n
  complex(c8) :: dots(n)
  end function cvec_mdot_vec
  !------------------------------------------------------------------------------
  !> Dot product with an arrays of vectors
  !!
  !! @param[in] a Array of vectors for dot product [n]
  !! @param[in] n Length of vector array
  !! @result \f$ \sum_i self_i a(j)_i \f$
  !------------------------------------------------------------------------------
  function cvec_mdot_cvec(self,a,n) result(dots)
  import oft_cvector, oft_cvector_ptr, c8, i4
  class(oft_cvector), intent(inout) :: self
  type(oft_cvector_ptr), intent(inout) :: a(n)
  integer(i4), intent(in) :: n
  complex(c8) :: dots(n)
  end function cvec_mdot_cvec
  !------------------------------------------------------------------------------
  !> Sum reduction over a field
  !!
  !! @result sum(self)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  function cvec_sum(self) result(sum)
  import oft_cvector, c8
  class(oft_cvector), intent(inout) :: self
  complex(c8) :: sum
  end function cvec_sum
  !------------------------------------------------------------------------------
  !> Norm of a field
  !!
  !! @param[in] itype Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
  !! @result norm(self)
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  function cvec_norm(self,itype) result(norm)
  import oft_cvector, i4, r8
  class(oft_cvector), intent(inout) :: self
  integer(i4), intent(in) :: itype
  REAL(r8) :: norm
  end function cvec_norm
  !------------------------------------------------------------------------------
  !> Perform global stitching
  !!
  !! @param[in] up_method Type of stitching to perform
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized fields.
  !------------------------------------------------------------------------------
  subroutine cvec_stitch(self,up_method)
  import oft_cvector, i4
  class(oft_cvector), intent(inout) :: self
  integer(i4), intent(in) :: up_method
  end subroutine cvec_stitch
END INTERFACE
!------------------------------------------------------------------------------
!> Sub-matrix mapping
!------------------------------------------------------------------------------
type, public :: oft_matrix_map
  integer(i4) :: nnz = 0 !< Number of non-zeros in block
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: kr => NULL() !< Row pointer to column list
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: lc => NULL() !< Column list
  integer(i4), pointer, CONTIGUOUS, dimension(:) :: lc_map => NULL() !< List of entries in full matrix
  integer(i4), pointer, CONTIGUOUS, dimension(:,:) :: ext => NULL() !< Row extents within full matrix
end type oft_matrix_map
!------------------------------------------------------------------------------
!> CSR graph representation
!------------------------------------------------------------------------------
type, public :: oft_graph
  INTEGER(i4) :: nr = 0 !< Local number of rows
  INTEGER(i4) :: nc = 0 !< Local number of columns
  INTEGER(i8) :: nrg = 0 !< Gobal number of rows
  INTEGER(i8) :: ncg = 0 !< Gobal number of columns
  INTEGER(i4) :: nnz = 0 !< Number of local non-zeros
  INTEGER(i4), CONTIGUOUS, POINTER, dimension(:) :: kr => NULL() !< Row pointer to column list
  INTEGER(i4), CONTIGUOUS, POINTER, dimension(:) :: lc => NULL() !< Column list
end type oft_graph
!------------------------------------------------------------------------------
!> Graph pointer
!------------------------------------------------------------------------------
type, public :: oft_graph_ptr
  type(oft_graph), pointer :: g => NULL() !< Graph
end type oft_graph_ptr
!------------------------------------------------------------------------------
!> Abstract matrix class
!!
!! Basic class for OFT matices (ex. fem operators)
!------------------------------------------------------------------------------
type, abstract, public :: oft_matrix
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
  class(oft_vector), pointer :: D => NULL() !< Diagonal entries for scaling
contains
  procedure(mat_apply_vec), deferred :: apply_real
  procedure :: apply_complex => matrix_apply_cvec
  !> Apply the matrix
  generic :: apply => apply_real, apply_complex
  procedure(mat_apply_vec), deferred :: applyt_real
  procedure :: applyt_complex => matrix_applyt_cvec
  !> Apply the matrix transpose
  generic :: applyt => applyt_real, applyt_complex
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
  !> Delete matrix
  procedure :: delete => matrix_delete
end type oft_matrix
!------------------------------------------------------------------------------
!> Matrix pointer
!------------------------------------------------------------------------------
type, public :: oft_matrix_ptr
  class(oft_matrix), pointer :: m => NULL() !< Matrix
end type oft_matrix_ptr
!------------------------------------------------------------------------------
!> Abstract complex matrix class
!!
!! Basic class for OFT matices (ex. fem operators)
!------------------------------------------------------------------------------
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
  procedure :: apply_real => cmatrix_apply_vec
  procedure(cmat_apply_cvec), deferred :: apply_complex
  !> Apply the matrix
  generic :: apply => apply_real, apply_complex
  procedure :: applyt_real => cmatrix_applyt_vec
  procedure(cmat_apply_cvec), deferred :: applyt_complex
  !> Apply the matrix transpose
  generic :: applyt => applyt_real, applyt_complex
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
!------------------------------------------------------------------------------
!> Matrix pointer
!------------------------------------------------------------------------------
type, public :: oft_cmatrix_ptr
  class(oft_cmatrix), pointer :: m => NULL() !< Matrix
end type oft_cmatrix_ptr
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Apply the matrix to a field (real)
  !!
  !! b = self * a
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_apply_vec(self,a,b)
  import oft_vector, oft_matrix
  class(oft_matrix), intent(inout) :: self
  class(oft_vector), target, intent(inout) :: a !< Source vector
  class(oft_vector), intent(inout) :: b !< Result of matrix product
  end subroutine mat_apply_vec
  !------------------------------------------------------------------------------
  !> Apply the matrix to a field (complex)
  !!
  !! b = self * a
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_apply_cvec(self,a,b)
  import oft_cvector, oft_matrix
  class(oft_matrix), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a !< Source vector
  class(oft_cvector), intent(inout) :: b !< Result of matrix product
  end subroutine mat_apply_cvec
  !------------------------------------------------------------------------------
  !> Set values of a matrix
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
  import oft_matrix, i4, r8
  class(oft_matrix), intent(inout) :: self
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
  real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  end subroutine mat_set_values
  !------------------------------------------------------------------------------
  !> Add values to a matrix
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
  import oft_matrix, i4, r8
  class(oft_matrix), intent(inout) :: self
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
  real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
  end subroutine mat_add_values
  !------------------------------------------------------------------------------
  !> Finish assembly of matrix and optionally extract real diagonal
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_assemble(self,diag)
  import oft_matrix, oft_vector
  class(oft_matrix), intent(inout) :: self
  class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
  end subroutine mat_assemble
  !------------------------------------------------------------------------------
  !> Zero all entries in matrix
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_zero(self)
  import oft_matrix
  class(oft_matrix), intent(inout) :: self
  end subroutine mat_zero
  !------------------------------------------------------------------------------
  !> Zero all entries in the specified rows
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine mat_zero_rows(self,nrows,irows,iblock,keep_diag)
  import oft_matrix, i4
  class(oft_matrix), intent(inout) :: self
  integer(i4), intent(in) :: nrows !< Number of rows to zero
  integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
  end subroutine mat_zero_rows
END INTERFACE
!
ABSTRACT INTERFACE
  ! !------------------------------------------------------------------------------
  ! !> Apply the matrix to a field and optionally add it to an existing field.
  ! !!
  ! !! b = self * a
  ! !!
  ! !! @note This subroutine is a dummy routine used to specify the interface
  ! !! of the member function and catch errors in uninitialized matrices.
  ! !!
  ! !! @param[in] a Source field
  ! !! @param[in,out] b Result of matrix product
  ! !! @param[in] beta Factor for matrix addition
  ! !------------------------------------------------------------------------------
  ! subroutine cmat_apply_vec(self,a,b)
  ! import oft_cvector, oft_vector, oft_cmatrix
  ! class(oft_cmatrix), intent(inout) :: self
  ! class(oft_vector), target, intent(inout) :: a
  ! class(oft_cvector), intent(inout) :: b
  ! end subroutine cmat_apply_vec
  !------------------------------------------------------------------------------
  !> Apply the matrix to a field
  !!
  !! b = self * a
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_apply_cvec(self,a,b)
  import oft_cvector, oft_cmatrix
  class(oft_cmatrix), intent(inout) :: self
  class(oft_cvector), target, intent(inout) :: a !< Source vector
  class(oft_cvector), intent(inout) :: b !< Result of matrix product
  end subroutine cmat_apply_cvec
  !------------------------------------------------------------------------------
  !> Set values of a matrix
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
  import oft_cmatrix, i4, c8
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
  complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  end subroutine cmat_set_values
  !------------------------------------------------------------------------------
  !> Add values to a matrix.
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
  import oft_cmatrix, i4, c8
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
  integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
  complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
  integer(i4), intent(in) :: n !< Number of rows in local matrix
  integer(i4), intent(in) :: m !< Number of columns in local matrix
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  integer(i4), optional, intent(in) :: jblock !< Column block (optional)
  integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
  end subroutine cmat_add_values
  !------------------------------------------------------------------------------
  !> Finish assembly of matrix and optionally extract diagonals
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_assemble(self,diag)
  import oft_cmatrix, oft_cvector
  class(oft_cmatrix), intent(inout) :: self
  class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
  end subroutine cmat_assemble
  !------------------------------------------------------------------------------
  !> Zero all entries in matrix
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_zero(self)
  import oft_cmatrix
  class(oft_cmatrix), intent(inout) :: self
  end subroutine cmat_zero
  !------------------------------------------------------------------------------
  !> Zero all entries in the specified rows
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized matrices
  !------------------------------------------------------------------------------
  subroutine cmat_zero_rows(self,nrows,irows,iblock,keep_diag)
  import oft_cmatrix, i4
  class(oft_cmatrix), intent(inout) :: self
  integer(i4), intent(in) :: nrows !< Number of rows to zero
  integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
  integer(i4), optional, intent(in) :: iblock !< Row block (optional)
  logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
  end subroutine cmat_zero_rows
END INTERFACE
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_ml_vecspace
CONTAINS
  !> Needs docs
  PROCEDURE(oft_veccreate_proto), DEFERRED :: vec_create
  !> Needs docs
  PROCEDURE(oft_interp_proto), DEFERRED :: interp
  !> Needs docs
  PROCEDURE(oft_inject_proto), DEFERRED :: inject
END TYPE oft_ml_vecspace
!---Interfaces
ABSTRACT INTERFACE
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE oft_inject_proto(self,afine,acors)
    IMPORT oft_ml_vecspace, oft_vector
    CLASS(oft_ml_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: afine
    CLASS(oft_vector), INTENT(inout) :: acors
  END SUBROUTINE oft_inject_proto
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE oft_interp_proto(self,acors,afine)
    IMPORT oft_ml_vecspace, oft_vector
    CLASS(oft_ml_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: acors
    CLASS(oft_vector), INTENT(inout) :: afine
  END SUBROUTINE oft_interp_proto
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE oft_veccreate_proto(self,new,level,cache,native)
    IMPORT oft_ml_vecspace, oft_vector, i4
    CLASS(oft_ml_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), POINTER, INTENT(out) :: new
    INTEGER(i4), OPTIONAL, INTENT(in) :: level
    LOGICAL, OPTIONAL, INTENT(in) :: cache
    LOGICAL, OPTIONAL, INTENT(in) :: native
  END SUBROUTINE oft_veccreate_proto
END INTERFACE
!---Declare public entities
public vector_extrapolate
contains
!------------------------------------------------------------------------------
!> Finalize field
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized fields.
!------------------------------------------------------------------------------
subroutine vector_delete(self)
class(oft_vector), intent(inout) :: self
call oft_warn('Finalizing general real vector, this may indicate an error.')
end subroutine vector_delete
!------------------------------------------------------------------------------
!> Finalize field
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized fields.
!------------------------------------------------------------------------------
subroutine cvector_delete(self)
class(oft_cvector), intent(inout) :: self
call oft_warn('Finalizing general complex vector, this may indicate an error.')
end subroutine cvector_delete
!------------------------------------------------------------------------------
!> Apply the matrix to a field and optionally add it to an existing field.
!!
!! b = self^T * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[in,out] b Result of matrix product
!! @param[in] beta Factor for matrix addition
!------------------------------------------------------------------------------
subroutine matrix_apply_cvec(self,a,b)
class(oft_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
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
CALL b%add((0.d0,0.d0),(1.d0,0.d0),btmp)
!---Apply imaginary part
CALL atmp%restore_slice(AIMAG(avals))
CALL self%apply(atmp,btmp)
CALL b%add((1.d0,0.d0),(0.d0,1.d0),btmp)
!---Cleanup
CALL atmp%delete()
CALL btmp%delete()
DEALLOCATE(atmp,btmp,avals)
end subroutine matrix_apply_cvec
!------------------------------------------------------------------------------
!> Apply the matrix to a field and optionally add it to an existing field.
!!
!! b = self^T * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[in,out] b Result of matrix product
!! @param[in] beta Factor for matrix addition
!------------------------------------------------------------------------------
subroutine matrix_applyt_cvec(self,a,b)
class(oft_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
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
CALL b%add((0.d0,0.d0),(1.d0,0.d0),btmp)
!---Apply imaginary part
CALL atmp%restore_slice(AIMAG(avals))
CALL self%applyt(atmp,btmp)
CALL b%add((1.d0,0.d0),(0.d0,1.d0),btmp)
!---Cleanup
CALL atmp%delete()
CALL btmp%delete()
DEALLOCATE(atmp,btmp,avals)
end subroutine matrix_applyt_cvec
!------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine matrix_delete(self)
class(oft_matrix), intent(inout) :: self
call oft_warn('Finalizing general real matrix, this may indicate an error.')
end subroutine matrix_delete
!------------------------------------------------------------------------------
!> Apply the matrix to a field and optionally add it to an existing field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[in,out] b Result of matrix product
!! @param[in] beta Factor for matrix addition
!------------------------------------------------------------------------------
subroutine cmatrix_apply_vec(self,a,b)
class(oft_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
class(oft_vector), pointer :: atmp
CALL a%new(atmp)
CALL atmp%add(0.d0,1.d0,a)
CALL self%apply(atmp,b)
CALL atmp%delete()
DEALLOCATE(atmp)
end subroutine cmatrix_apply_vec
!------------------------------------------------------------------------------
!> Apply the matrix to a field and optionally add it to an existing field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[in,out] b Result of matrix product
!! @param[in] beta Factor for matrix addition
!------------------------------------------------------------------------------
subroutine cmatrix_applyt_vec(self,a,b)
class(oft_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
class(oft_vector), pointer :: atmp
CALL a%new(atmp)
CALL atmp%add(0.d0,1.d0,a)
CALL self%applyt(atmp,b)
CALL atmp%delete()
DEALLOCATE(atmp)
end subroutine cmatrix_applyt_vec
!------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine cmatrix_delete(self)
class(oft_cmatrix), intent(inout) :: self
call oft_warn('Finalizing general complex matrix, this may indicate an error.')
end subroutine cmatrix_delete
!---------------------------------------------------------------------------
! SUBROUTINE: vector_extrapolate
!---------------------------------------------------------------------------
!> Extrapolate a field using a polynomial fit to previous fields
!!
!! @param[in] x Array of previous positions [n]
!! @param[in,out] fields Array of previous fields [n]
!! @param[in] n Number of fields to use for interpolation
!! @param[in] xe Position to extrapolate to
!! @param[in,out] output Extrapolated field
!---------------------------------------------------------------------------
subroutine vector_extrapolate(x,fields,n,xe,output)
real(r8), intent(in) :: x(:)
type(oft_vector_ptr), intent(inout) :: fields(:)
integer(i4), intent(in) :: n
real(r8), intent(in) :: xe
class(oft_vector), intent(inout) :: output
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
  CALL output%add(1.d0,y,fields(i)%f)
END DO
DEBUG_STACK_POP
end subroutine vector_extrapolate
end module oft_la_base
