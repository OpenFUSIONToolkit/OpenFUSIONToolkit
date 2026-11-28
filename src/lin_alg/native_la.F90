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
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int8_t
USE oft_local
USE oft_base
USE oft_sort, ONLY: search_array, sort_array
USE oft_stitching, ONLY: oft_seam, oft_global_stitch, oft_stitch_check, seam_free_plan, &
  seam_build_own
USE oft_io
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_cvector_ptr, &
  oft_matrix, oft_cmatrix, oft_map, oft_matrix_map, oft_graph, OFT_ACCESS_READ, OFT_ACCESS_WRITE, &
  OFT_ACCESS_READWRITE, OFT_ACCESS_ACCUMULATE
USE oft_native_kernels
USE oft_scatter, ONLY: oft_scatter_plan, scatter_plan_setup_local, scatter_plan_delete, scatter_local, &
  OFT_SCATTER_INSERT, OFT_SCATTER_ADD
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
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: acc => NULL() !< Accumulation storage (`OFT_ACCESS_ACCUMULATE`, persistent)
  LOGICAL :: acc_open = .FALSE. !< Accumulation in progress?
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
  !> Borrow local values without a copy
  PROCEDURE :: borrow_local => vec_borrow_local
  !> Release borrowed local values
  PROCEDURE :: release_local => vec_release_local
  !> Copy values from a second vector
  PROCEDURE :: copy => vec_copy
  !> Add a linear combination of vectors
  PROCEDURE :: maxpy => vec_maxpy
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
  COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: acc => NULL() !< Accumulation storage (`OFT_ACCESS_ACCUMULATE`, persistent)
  LOGICAL :: acc_open = .FALSE. !< Accumulation in progress?
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
  !> Borrow local values without a copy
  PROCEDURE :: borrow_local => cvec_borrow_local
  !> Release borrowed local values
  PROCEDURE :: release_local => cvec_release_local
  PROCEDURE :: copy_real => cvec_copy_vec
  PROCEDURE :: copy_complex => cvec_copy_cvec
  !> Add a linear combination of vectors
  PROCEDURE :: maxpy => cvec_maxpy
  !> Delete vector
  procedure :: delete => cvec_delete
end type oft_native_cvector
!---------------------------------------------------------------------------------
!> Transposed (CSC) index for the non-zero pattern of a native CRS matrix
!!
!! Built on first use by the transpose products and released with the pattern. Entries
!! within each column are ordered by increasing row index.
!---------------------------------------------------------------------------------
TYPE :: native_csc_index
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kc => NULL() !< Column pointer to row list [nc+1]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lr => NULL() !< Row index of each entry [nnz]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lm => NULL() !< Location of each entry in CRS value array [nnz]
END TYPE native_csc_index
!---------------------------------------------------------------------------------
!> Location of the diagonal entry in each row of a native CRS matrix
!!
!! Built on first use by matrix assembly and released with the pattern.
!---------------------------------------------------------------------------------
TYPE :: native_diag_index
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: loc => NULL() !< Location of diagonal in CRS value array (0 if absent) [nr]
END TYPE native_diag_index
!---------------------------------------------------------------------------------
!> Non-zero pattern of native CRS matrices, shared between matrices by reference counting
!!
!! Holds everything that depends only on the location of the non-zeros: the CRS
!! structure (`kr`, `lc`), the block maps, the 2-byte column list, the transposed
!! (CSC) index, the diagonal locations and the seam and redundant-row data used to
!! assemble full local matrices (@ref native_matrix_setup_full). Matrices with the
!! same structure (eg. several operators on one FE space) can share one pattern, each
!! holding only its own values.
!!
!! Ownership rules:
!! - A pattern is created by @ref native_pattern_create, which takes ownership of the
!!   CRS arrays and block maps passed to it, or by @ref native_matrix_adopt_pattern
!!   from a matrix whose `kr`/`lc` were set directly.
!! - Every matrix using a pattern holds one reference (@ref native_matrix_set_pattern,
!!   @ref native_matrix_share). Other holders (eg. a pattern cache of an FE
!!   representation) take a reference with @ref native_pattern_retain.
!! - @ref native_pattern_release drops one reference; the pattern and everything it
!!   owns are freed when the last reference is released. Deleting a matrix releases
!!   its reference.
!! - Reference count updates are atomic, so references may be taken and released from
!!   different threads (eg. block matrices created in a parallel loop). Only the
!!   thread releasing the last reference frees the pattern.
!! - A pattern must not be modified once created; derived data (2-byte column list,
!!   CSC index, diagonal locations, linkage) is built on first use and is valid for
!!   the lifetime of the pattern.
!!
!! `kr`, `lc`, `map`, `linkage`, `redind`, `nred` and `nnz` of a matrix are aliases of
!! the pattern components and must be treated as read-only.
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_native_pattern
  INTEGER(i4) :: refcount = 0 !< Number of references held (matrices and other holders)
  INTEGER(i4) :: nr = 0 !< Number of rows
  INTEGER(i4) :: nc = 0 !< Number of columns
  INTEGER(i4) :: nnz = 0 !< Number of non-zero entries
  INTEGER(i4) :: ni = 0 !< Number of row blocks in `map`
  INTEGER(i4) :: nj = 0 !< Number of column blocks in `map`
  INTEGER(i4) :: nred = 0 !< Number of redundant rows (set by @ref native_matrix_setup_full)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr => NULL() !< Row pointer to column list [nr+1]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lc => NULL() !< Column list [nnz]
  INTEGER(2), POINTER, CONTIGUOUS, DIMENSION(:) :: lc_small => NULL() !< 2-byte copy of column list (built by assembly)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: redind => NULL() !< Redundant rows [nred]
  TYPE(oft_matrix_map), POINTER, DIMENSION(:,:) :: map => NULL() !< Block maps [ni,nj] (optional)
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Seam for full local matrices (built by @ref native_matrix_setup_full)
  TYPE(native_csc_index), PRIVATE :: csc !< Transposed index (built by transpose products)
  TYPE(native_diag_index), PRIVATE :: diag_index !< Diagonal locations (built by assembly)
  TYPE(native_extraction), POINTER :: extracts => NULL() !< Cached submatrix extractions (see @ref native_extraction_get)
END TYPE oft_native_pattern
!---------------------------------------------------------------------------------
!> Extraction of a submatrix (row set by column set) from a native CRS pattern
!!
!! Holds the symbolic part of the extraction: the non-zero pattern of the submatrix
!! (an @ref oft_native_pattern, shared by every submatrix using this extraction) and
!! the location in the parent's value array of each submatrix entry, as a local
!! scatter plan. Extractions are cached in a list held by the parent pattern (see
!! @ref native_extraction_get), so submatrices of matrices sharing a pattern, or
!! rebuilt for new values, reuse the symbolic work; they are freed with the parent
!! pattern.
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_extraction
  INTEGER(i4) :: nr = 0 !< Number of rows of the submatrix
  INTEGER(i4) :: nc = 0 !< Number of columns of the submatrix
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: rows => NULL() !< Parent row of each submatrix row [nr]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: cols => NULL() !< Parent column of each submatrix column [nc]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lcmap => NULL() !< Parent value location of each submatrix entry [nnz]
  TYPE(oft_native_pattern), POINTER :: pattern => NULL() !< Pattern of the submatrix (one reference held)
  TYPE(oft_scatter_plan) :: vplan !< Value gather plan (parent values to submatrix values)
  TYPE(native_extraction), POINTER :: next => NULL() !< Next extraction in the parent's cache
END TYPE native_extraction
!---------------------------------------------------------------------------------
!> Native CRS matrix class
!---------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_native_matrix
  logical :: full_current = .FALSE. !< Is full matrix current?
  INTEGER(i4) :: nnz = 0 !< Number of non-zero entries (alias of `pattern%nnz`)
  INTEGER(i4) :: nred = 0 !< Number of redundant rows (alias of `pattern%nred`)
  TYPE(oft_native_pattern), POINTER :: pattern => NULL() !< Non-zero pattern (possibly shared)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: redind => NULL() !< Redundant rows (alias of `pattern%redind`)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr => NULL() !< Row pointer to column list (alias of `pattern%kr`)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lc => NULL() !< Column list (alias of `pattern%lc`)
  real(r8), pointer, contiguous, dimension(:) :: M => NULL() !< Matrix entries
  real(r8), pointer, contiguous, dimension(:) :: Md => NULL() !< Diagonal matrix entries
  real(r8), pointer, contiguous, dimension(:) :: Mfull => NULL() !< Full local matrix entries
  TYPE(oft_matrix_map), POINTER, DIMENSION(:,:) :: map => NULL() !< Block maps (alias of `pattern%map`)
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Seam for full local matrix (alias of `pattern%linkage`)
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
  !> Save matrix
  procedure :: save => mat_save
end type oft_native_matrix
!---------------------------------------------------------------------------------
!> Native CRS complex matrix class
!---------------------------------------------------------------------------------
type, public, extends(oft_cmatrix) :: oft_native_cmatrix
  logical :: full_current = .FALSE. !< Is full matrix current?
  INTEGER(i4) :: nnz = 0 !< Number of non-zero entries (alias of `pattern%nnz`)
  INTEGER(i4) :: nred = 0 !< Number of redundant rows (alias of `pattern%nred`)
  integer(i4), pointer, contiguous, dimension(:) :: color => NULL() !< Element coloring
  TYPE(oft_native_pattern), POINTER :: pattern => NULL() !< Non-zero pattern (possibly shared)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: redind => NULL() !< Redundant rows (alias of `pattern%redind`)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr => NULL() !< Row pointer to column list (alias of `pattern%kr`)
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lc => NULL() !< Column list (alias of `pattern%lc`)
  complex(c8), pointer, contiguous, dimension(:) :: M => NULL() !< Matrix entries
  complex(c8), pointer, contiguous, dimension(:) :: Md => NULL() !< Diagonal matrix entries
  complex(c8), pointer, contiguous, dimension(:) :: Mfull => NULL() !< Full local matrix entries
  TYPE(oft_matrix_map), POINTER, DIMENSION(:,:) :: map => NULL() !< Block maps (alias of `pattern%map`)
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Seam for full local matrix (alias of `pattern%linkage`)
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
  INTEGER(i4), POINTER, DIMENSION(:) :: lcmap => NULL() !< Location in the parent's full values of each entry [nnz] (alias of `ext%lcmap`)
  INTEGER(i4), POINTER, DIMENSION(:) :: part => NULL() !< Parent rows (and columns) of the submatrix
  class(oft_native_matrix), pointer :: parent => NULL() !< Parent matrix
  TYPE(oft_native_pattern), POINTER :: parent_pattern => NULL() !< Pattern of the parent (reference held, owns `ext`)
  TYPE(native_extraction), POINTER :: ext => NULL() !< Extraction (symbolic map) used
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
  !> Delete matrix
  procedure :: delete => dense_mat_delete
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
  !> Delete matrix
  procedure :: delete => dense_cmat_delete
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
!> Empty ownership weights for seams without shared entries (see @ref native_seam_owned)
INTEGER(c_int8_t), TARGET :: native_no_own(1) = 1_c_int8_t
!---------------------------------------------------------------------------------
!> Get a contiguous pointer to the local values of a vector (in place for native vectors)
!---------------------------------------------------------------------------------
INTERFACE native_get_vals
  MODULE PROCEDURE native_get_rvals
  MODULE PROCEDURE native_get_cvals
END INTERFACE native_get_vals
!---------------------------------------------------------------------------------
!> Finish writing values obtained with @ref native_get_vals
!---------------------------------------------------------------------------------
INTERFACE native_put_vals
  MODULE PROCEDURE native_put_rvals
  MODULE PROCEDURE native_put_cvals
END INTERFACE native_put_vals
!> Number of native CRS patterns currently allocated (for diagnostics and testing)
INTEGER(i4), PUBLIC, PROTECTED :: native_pattern_live = 0
!---------------------------------------------------------------------------------
!> Use a non-zero pattern for a native CRS matrix (real or complex)
!---------------------------------------------------------------------------------
INTERFACE native_matrix_set_pattern
  MODULE PROCEDURE native_matrix_set_pattern_real
  MODULE PROCEDURE native_matrix_set_pattern_comp
END INTERFACE native_matrix_set_pattern
!---------------------------------------------------------------------------------
!> Create a pattern from the CRS structure set directly on a native matrix (real or complex)
!---------------------------------------------------------------------------------
INTERFACE native_matrix_adopt_pattern
  MODULE PROCEDURE native_matrix_adopt_pattern_real
  MODULE PROCEDURE native_matrix_adopt_pattern_comp
END INTERFACE native_matrix_adopt_pattern
!---------------------------------------------------------------------------------
!> Create a native CRS matrix (real or complex) sharing the pattern of an existing one
!---------------------------------------------------------------------------------
INTERFACE native_matrix_share
  MODULE PROCEDURE native_matrix_share_real
  MODULE PROCEDURE native_matrix_share_comp
END INTERFACE native_matrix_share
!---Declare public entities
public native_vector_cast, native_cvector_cast, native_matrix_cast, native_cmatrix_cast
public partition_graph, native_matrix_setup_full
public native_vector_slice_push, native_vector_slice_pop
PUBLIC native_pattern_create, native_pattern_retain, native_pattern_release
PUBLIC native_matrix_set_pattern, native_matrix_adopt_pattern, native_matrix_share
PUBLIC native_extraction_get
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
!> Get the ownership weights of a seam for reductions over owned entries
!!
!! Returns `.FALSE.` if reductions are skipped for this seam (`seam%skip`). When the
!! seam has shared entries, `w` holds the 0/1 ownership weight of each of the `nw=n`
!! local entries (built on first use with `seam_build_own`, which must happen outside
!! of OpenMP parallel regions); otherwise `nw=0` and every entry is owned. The result
!! can always be passed to the `nk_owned_*` kernels.
!---------------------------------------------------------------------------------
FUNCTION native_seam_owned(seam,n,nw,w) RESULT(active)
TYPE(oft_seam), INTENT(inout) :: seam !< Seam structure
INTEGER(i4), INTENT(in) :: n !< Number of local entries
INTEGER(i4), INTENT(out) :: nw !< Number of weights (`n` or 0)
INTEGER(c_int8_t), POINTER, CONTIGUOUS, INTENT(out) :: w(:) !< Ownership weight of each entry [nw]
LOGICAL :: active !< Perform reduction?
active=(.NOT.seam%skip)
nw=0
w=>native_no_own
IF(active.AND.(seam%nbe>0))THEN
  IF(seam%nown/=n)CALL seam_build_own(seam,n)
  nw=n
  w=>seam%own
END IF
END FUNCTION native_seam_owned
!---------------------------------------------------------------------------------
!> Get a contiguous pointer to the local values of a real vector
!!
!! Native vectors are accessed in place, other vectors are copied with `get_local`.
!---------------------------------------------------------------------------------
FUNCTION native_get_rvals(a,vals) RESULT(copied)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Vector
REAL(r8), POINTER, CONTIGUOUS, INTENT(out) :: vals(:) !< Local values
LOGICAL :: copied !< Values were copied (deallocate when done)?
SELECT TYPE(a)
CLASS IS(oft_native_vector)
  vals=>a%v
  copied=.FALSE.
CLASS DEFAULT
  NULLIFY(vals)
  CALL a%get_local(vals)
  copied=.TRUE.
END SELECT
END FUNCTION native_get_rvals
!---------------------------------------------------------------------------------
!> Get a contiguous pointer to the local values of a complex vector
!!
!! Native vectors are accessed in place, other vectors are copied with `get_local`.
!---------------------------------------------------------------------------------
FUNCTION native_get_cvals(a,vals) RESULT(copied)
CLASS(oft_cvector), TARGET, INTENT(inout) :: a !< Vector
COMPLEX(c8), POINTER, CONTIGUOUS, INTENT(out) :: vals(:) !< Local values
LOGICAL :: copied !< Values were copied (deallocate when done)?
SELECT TYPE(a)
CLASS IS(oft_native_cvector)
  vals=>a%v
  copied=.FALSE.
CLASS DEFAULT
  NULLIFY(vals)
  CALL a%get_local(vals)
  copied=.TRUE.
END SELECT
END FUNCTION native_get_cvals
!---------------------------------------------------------------------------------
!> Finish writing values obtained with @ref native_get_rvals
!!
!! Copied values are restored (without synchronization) and deallocated.
!---------------------------------------------------------------------------------
SUBROUTINE native_put_rvals(a,vals,copied)
CLASS(oft_vector), INTENT(inout) :: a !< Vector
REAL(r8), POINTER, CONTIGUOUS, INTENT(inout) :: vals(:) !< Local values
LOGICAL, INTENT(in) :: copied !< Values were copied?
IF(copied)THEN
  CALL a%restore_local(vals,wait=.TRUE.)
  DEALLOCATE(vals)
END IF
END SUBROUTINE native_put_rvals
!---------------------------------------------------------------------------------
!> Finish writing values obtained with @ref native_get_cvals
!!
!! Copied values are restored (without synchronization) and deallocated.
!---------------------------------------------------------------------------------
SUBROUTINE native_put_cvals(a,vals,copied)
CLASS(oft_cvector), INTENT(inout) :: a !< Vector
COMPLEX(c8), POINTER, CONTIGUOUS, INTENT(inout) :: vals(:) !< Local values
LOGICAL, INTENT(in) :: copied !< Values were copied?
IF(copied)THEN
  CALL a%restore_local(vals,wait=.TRUE.)
  DEALLOCATE(vals)
END IF
END SUBROUTINE native_put_cvals
!---------------------------------------------------------------------------------
!> Locate a block of entries in a native CRS non-zero pattern
!!
!! `loc(i,j)` is set to the location of entry (`i_inds(i)`,`j_inds(j)`) in the CRS value
!! array, or 0 if the entry is not in the pattern. When `iblock` and `jblock` are given
!! the indices are local to that block and the row extents of the block map are
!! searched. In that case `loc_cache` may hold the positions within each row extent
!! from a previous call with the same indices (`loc_cache(1,1)>0`), which are then
!! used without searching, or is filled for later calls.
!---------------------------------------------------------------------------------
SUBROUTINE native_locate(kr,lc,map,j_map,n,m,i_inds,j_inds,loc,iblock,jblock,loc_cache)
INTEGER(i4), INTENT(in) :: kr(:) !< CRS row pointer
INTEGER(i4), INTENT(in) :: lc(:) !< CRS column indices
TYPE(oft_matrix_map), POINTER, INTENT(in) :: map(:,:) !< Block maps (used with `iblock`/`jblock`)
TYPE(oft_map), POINTER, INTENT(in) :: j_map(:) !< Column block mapping (used with `jblock`)
INTEGER(i4), INTENT(in) :: n !< Number of rows
INTEGER(i4), INTENT(in) :: m !< Number of columns
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices [m]
INTEGER(i4), INTENT(out) :: loc(n,m) !< Locations in CRS value array (0 if not found) [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of positions within row extents (optional)
INTEGER(i4) :: i,j,jp,jn,jk,joff
LOGICAL :: use_cache
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  use_cache=.FALSE.
  IF(PRESENT(loc_cache))use_cache=(loc_cache(1,1)>0)
  joff=j_map(jblock)%offset
  DO i=1,n
    jp=map(iblock,jblock)%ext(1,i_inds(i))
    jn=map(iblock,jblock)%ext(2,i_inds(i))
    DO j=1,m
      IF(use_cache)THEN
        jk=loc_cache(i,j)
      ELSE
        jk=search_array(joff+j_inds(j),lc(jp:jn),jn-jp+1)
        IF(PRESENT(loc_cache))loc_cache(i,j)=jk
      END IF
      loc(i,j)=0
      IF(jk>0)loc(i,j)=jk+jp-1
    END DO
  END DO
ELSE
  IF(PRESENT(iblock).OR.PRESENT(jblock))CALL oft_abort('Only one block index was supplied', &
    'native_locate',__FILE__)
  DO i=1,n
    jp=kr(i_inds(i))
    jn=kr(i_inds(i)+1)-jp
    DO j=1,m
      jk=search_array(j_inds(j),lc(jp:jp+jn-1),jn)
      loc(i,j)=0
      IF(jk>0)loc(i,j)=jk+jp-1
    END DO
  END DO
END IF
END SUBROUTINE native_locate
!---------------------------------------------------------------------------------
!> Build 2-byte copy of a CRS column list, if the number of columns allows it
!!
!! Columns are stored as `lc-lc_offset` (see @ref oft_native_kernels_base::lc_offset
!! "lc_offset"). Nothing is done if the copy already exists or if there are
!! `2*lc_offset` or more columns.
!---------------------------------------------------------------------------------
SUBROUTINE native_lc_small_build(lc_small,lc,nnz,nc)
INTEGER(2), POINTER, CONTIGUOUS, INTENT(inout) :: lc_small(:) !< 2-byte column list [nnz]
INTEGER(i4), INTENT(in) :: lc(:) !< CRS column list [nnz]
INTEGER(i4), INTENT(in) :: nnz !< Number of non-zeros
INTEGER(i4), INTENT(in) :: nc !< Number of columns
INTEGER(i4) :: i
IF(nc>=2*lc_offset)RETURN
IF(ASSOCIATED(lc_small))RETURN
ALLOCATE(lc_small(nnz))
!$omp simd
DO i=1,nnz
  lc_small(i)=INT(lc(i)-lc_offset,2)
END DO
END SUBROUTINE native_lc_small_build
!---------------------------------------------------------------------------------
!> Build the owned-entry (slice) scatter plan of a vector block, if not already built
!!
!! The plan gathers `slice(j)` of the block into entry `j` of the slice. Maps are shared
!! by all vectors of a space, so the plan is built once, on first use. The build is safe
!! when called concurrently from several threads (double-checked initialization): the
!! `slice_plan_ready` flag is read with a sequentially consistent atomic read, the plan
!! is built inside a named critical section after checking the flag again, and the flag
!! is set with a sequentially consistent atomic write only after the plan is complete,
!! so a thread that sees the flag set also sees the complete plan.
!---------------------------------------------------------------------------------
SUBROUTINE native_slice_plan(mp)
TYPE(oft_map), INTENT(inout) :: mp !< Block mapping
LOGICAL :: ready
!$omp atomic read seq_cst
ready=mp%slice_plan_ready
IF(ready)RETURN
!$omp critical (native_slice_plan_build)
!$omp atomic read seq_cst
ready=mp%slice_plan_ready
IF(.NOT.ready)THEN
  CALL scatter_plan_setup_local(mp%slice_plan,mp%nslice,mp%slice)
  !$omp atomic write seq_cst
  mp%slice_plan_ready=.TRUE.
END IF
!$omp end critical (native_slice_plan_build)
END SUBROUTINE native_slice_plan
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
INTEGER(i4) :: i1,i2
DEBUG_STACK_PUSH
random_flag=.FALSE.
IF(PRESENT(random))random_flag=random
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_set',__FILE__)
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  IF(random_flag)THEN
    CALL oft_random_number(self%v(i1:i2),self%map(iblock)%n)
  ELSE
    CALL nk_set(self%map(iblock)%n,alpha,self%v(i1:i2))
  END IF
  CALL self%stitch(0)
ELSE
  IF(random_flag)THEN
    CALL oft_random_number(self%v,self%n)
    CALL self%stitch(0)
  ELSE
    CALL nk_set(self%n,alpha,self%v)
  END IF
END IF
DEBUG_STACK_POP
end subroutine vec_set
!---------------------------------------------------------------------------------
!> Get values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
SUBROUTINE vec_get_slice(self,array,iblock)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
REAL(r8), POINTER, INTENT(inout) :: array(:) !< Slice values
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block to retrieve
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_get_slice',__FILE__)
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%nslice))
  ASSOCIATE(mp=>self%map(iblock))
    IF(mp%nslice>0)THEN
      CALL native_slice_plan(mp)
      CALL scatter_local(mp%slice_plan,self%v(mp%offset+1:mp%offset+mp%n),array(1:mp%nslice),OFT_SCATTER_INSERT)
    END IF
  END ASSOCIATE
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%nslice))
  IF(self%n==self%nslice)THEN
    CALL nk_copy(self%n,self%v,array)
  ELSE
    DO i=1,self%nblocks
      ASSOCIATE(mp=>self%map(i))
        IF(mp%nslice>0)THEN
          CALL native_slice_plan(mp)
          CALL scatter_local(mp%slice_plan,self%v(mp%offset+1:mp%offset+mp%n), &
            array(mp%soffset+1:mp%soffset+mp%nslice),OFT_SCATTER_INSERT)
        END IF
      END ASSOCIATE
    END DO
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE vec_get_slice
!---------------------------------------------------------------------------------
!> Set/add values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine vec_restore_slice(self,array,iblock,wait)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: wait !< Wait to perform global sync?
INTEGER(i4) :: i
logical :: do_wait
DEBUG_STACK_PUSH
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_slice',__FILE__)
  ASSOCIATE(mp=>self%map(iblock))
    IF(mp%nslice>0)THEN
      CALL native_slice_plan(mp)
      CALL scatter_local(mp%slice_plan,array(1:mp%nslice),self%v(mp%offset+1:mp%offset+mp%n), &
        OFT_SCATTER_INSERT,reverse=.TRUE.)
    END IF
  END ASSOCIATE
ELSE
  IF(self%n==self%nslice)THEN
    CALL nk_copy(self%n,array,self%v)
  ELSE
    DO i=1,self%nblocks
      ASSOCIATE(mp=>self%map(i))
        IF(mp%nslice>0)THEN
          CALL native_slice_plan(mp)
          CALL scatter_local(mp%slice_plan,array(mp%soffset+1:mp%soffset+mp%nslice), &
            self%v(mp%offset+1:mp%offset+mp%n),OFT_SCATTER_INSERT,reverse=.TRUE.)
        END IF
      END ASSOCIATE
    END DO
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
INTEGER(i4) :: i1,i2
DEBUG_STACK_PUSH
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_get_local',__FILE__)
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%n))
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  CALL nk_copy(self%map(iblock)%n,self%v(i1:i2),array)
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%n))
  CALL nk_copy(self%n,self%v,array)
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
INTEGER(i4) :: i1,i2
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
do_add=.FALSE.
IF(PRESENT(add))do_add=add
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_local',__FILE__)
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  IF(do_add)THEN
    IF(.NOT.ASSOCIATED(self%local_tmp))THEN
      ALLOCATE(self%local_tmp(self%n))
      CALL nk_set(self%n,0.d0,self%local_tmp)
    END IF
    CALL nk_copy(self%map(iblock)%n,array,self%local_tmp(i1:i2))
    IF(.NOT.do_wait)THEN
      CALL full_restore_local(self%local_tmp)
      DEALLOCATE(self%local_tmp)
    END IF
  ELSE
    CALL nk_copy(self%map(iblock)%n,array,self%v(i1:i2))
    IF(.NOT.do_wait)CALL self%stitch(0)
  END IF
ELSE
  IF(ASSOCIATED(self%local_tmp))CALL oft_abort("Existing partial restore found when restoring full vector", &
    "vec_restore_local",__FILE__)
  IF(do_wait)CALL oft_abort("Wait not supported when restoring full vector","vec_restore_local",__FILE__)
  IF(do_add)THEN
    !---Contributions must be summed before adding, which needs a modifiable copy
    ALLOCATE(atmp(self%n))
    CALL nk_copy(self%n,array,atmp)
    CALL full_restore_local(atmp)
    DEALLOCATE(atmp)
  ELSE
    CALL nk_copy(self%n,array,self%v)
    CALL self%stitch(0)
  END IF
END IF
DEBUG_STACK_POP
contains
!
subroutine full_restore_local(array_loc)
real(r8), intent(inout) :: array_loc(:)
CALL oft_global_stitch(self%stitch_info,array_loc,1)
CALL nk_axpby(self%n,1.d0,1.d0,array_loc,self%v)
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
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','vec_add_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','vec_add_vec',__FILE__)
  copied(2)=native_get_vals(b,btmp)
  CALL nk_axpbypcz(self%n,gamma,alpha,atmp,beta,btmp,self%v)
  IF(copied(2))DEALLOCATE(btmp)
else ! Add one vector to the source
  CALL nk_axpby(self%n,gamma,alpha,atmp,self%v)
end if
IF(copied(1))DEALLOCATE(atmp)
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
LOGICAL :: div,copied
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','vec_mult_vec',__FILE__)
copied=native_get_vals(a,atmp)
if(div)then
  CALL nk_div(self%n,atmp,self%v)
else
  CALL nk_mult(self%n,atmp,self%v)
end if
IF(copied)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine vec_mult
!---------------------------------------------------------------------------------
!> Scale vector by a scalar
!---------------------------------------------------------------------------------
subroutine vec_scale(self,alpha)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: alpha !< Scale factor
DEBUG_STACK_PUSH
CALL nk_scale(self%n,alpha,self%v)
DEBUG_STACK_POP
end subroutine vec_scale
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function vec_dot_vec(self,a) result(dot)
class(oft_native_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
real(r8) :: dot !< \f$ \sum_i self_i a_i \f$
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_dot_vec',__FILE__)
    dot=0.d0
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))dot=nk_owned_dot(self%n,self%v,this%v,nw,w)
    IF(.NOT.self%stitch_info%full)dot=oft_mpi_sum(dot)
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
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_cvector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','vec_dot_cvec',__FILE__)
    dot=(0.d0,0.d0)
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))dot=nk_owned_dot(self%n,self%v,this%v,nw,w)
    IF(.NOT.self%stitch_info%full)dot=oft_mpi_sum(dot)
  class default
    CALL oft_abort('"a" is not a native vector.','vec_dot_cvec',__FILE__)
end select
DEBUG_STACK_POP
end function vec_dot_cvec
!---------------------------------------------------------------------------------
!> Dot products with an array of vectors
!!
!! All products are computed in a single pass over row blocks (see
!! @ref oft_native_kernels_prod_rr::k_owned_mdot "nk_owned_mdot"), so each block of `self` is
!! loaded once for all `n` products. Values shared with other processes are counted
!! only on their owner and the results are combined with a single MPI reduction,
!! which is skipped when the local mesh is complete (`stitch_info%full`). Non-native
!! vectors in `a` are copied with `get_local` first.
!---------------------------------------------------------------------------------
FUNCTION vec_mdot_vec(self,a,n) RESULT(dots)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Length of vector array
TYPE(oft_vector_ptr), INTENT(inout) :: a(n) !< Array of vectors for dot product [n]
REAL(r8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
INTEGER(i4) :: k,nw
LOGICAL :: copied(n)
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
TYPE(native_rvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','vec_mdot_vec',__FILE__)
  copied(k)=native_get_rvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
dots=0.d0
IF(n>0)THEN
  IF(native_seam_owned(self%stitch_info,self%n,nw,w))CALL nk_owned_mdot(self%n,n,self%v,vals,nw,w,dots)
  IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
END IF
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END FUNCTION vec_mdot_vec
!---------------------------------------------------------------------------------
!> Dot products with an array of complex vectors
!!
!! Single pass over row blocks with a single MPI reduction, see @ref vec_mdot_vec.
!---------------------------------------------------------------------------------
FUNCTION vec_mdot_cvec(self,a,n) RESULT(dots)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Length of vector array
TYPE(oft_cvector_ptr), INTENT(inout) :: a(n) !< Array of vectors for dot product [n]
COMPLEX(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
INTEGER(i4) :: k,nw
LOGICAL :: copied(n)
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
TYPE(native_cvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','vec_mdot_cvec',__FILE__)
  copied(k)=native_get_cvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
dots=(0.d0,0.d0)
IF(n>0)THEN
  IF(native_seam_owned(self%stitch_info,self%n,nw,w))CALL nk_owned_mdot(self%n,n,self%v,vals,nw,w,dots)
  IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
END IF
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END FUNCTION vec_mdot_cvec
!---------------------------------------------------------------------------------
!> Sum reduction over vector
!---------------------------------------------------------------------------------
function vec_sum(self) result(sum)
class(oft_native_vector), intent(inout) :: self !< Vector object
real(r8) :: sum !< Sum of vector elements
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
sum=0.d0
IF(native_seam_owned(self%stitch_info,self%n,nw,w))sum=nk_owned_sum(self%n,self%v,nw,w)
IF(.NOT.self%stitch_info%full)sum=oft_mpi_sum(sum)
DEBUG_STACK_POP
end function vec_sum
!---------------------------------------------------------------------------------
!> Compute norm of vector
!---------------------------------------------------------------------------------
function vec_norm(self,itype) result(norm)
class(oft_native_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
real(r8) :: norm !< Specified norm of vector
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
SELECT CASE(itype)
  CASE(1)
    norm=0.d0
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))norm=nk_owned_norm1(self%n,self%v,nw,w)
    IF(.NOT.self%stitch_info%full)norm=oft_mpi_sum(norm)
  CASE(2)
    norm=SQRT(clamp_nonneg(self%dot(self))) ! Clamp negative round-off (safety net, NaN preserved)
  CASE(3)
    norm=oft_mpi_max(nk_norm_inf(self%n,self%v))
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
#ifdef OFT_DEBUG
!---------------------------------------------------------------------------------
!> Check that values shared between processes agree with their owner's value
!!
!! Debug builds only (`OFT_DEBUG`). Used by `release_local` with `sync=.FALSE.`, which
!! requires shared values to be consistent across processes at release. The values
!! are compared with the result of `stitch(0)` on a copy, and all processes abort if
!! any shared value differs from its owner's value by more than a relative 1.E-12.
!! This is a collective operation unless the seam is local to the process (`skip`, or
!! `full` when only periodic images on this process are compared), so it can be used
!! on local vectors from worker threads (eg. block sub-solvers) without MPI calls.
!---------------------------------------------------------------------------------
SUBROUTINE native_check_shared(seam,vals,caller)
TYPE(oft_seam), INTENT(inout) :: seam !< Seam structure of the vector
REAL(r8), INTENT(in) :: vals(:) !< Local values of the vector
CHARACTER(LEN=*), INTENT(in) :: caller !< Name of calling routine (for error message)
INTEGER(i4) :: i,j,nbad
REAL(r8) :: dmax
REAL(r8), ALLOCATABLE, DIMENSION(:) :: tmp
IF(seam%skip)RETURN
ALLOCATE(tmp(SIZE(vals)))
tmp=vals
CALL oft_global_stitch(seam,tmp,0)
nbad=0
dmax=0.d0
DO i=1,seam%nbe
  j=seam%lbe(i)
  !---A NaN on only one side compares false above, so check NaN mismatches separately
  IF((ABS(tmp(j)-vals(j))>1.d-12*MAX(ABS(tmp(j)),ABS(vals(j)))) &
    .OR.(IEEE_IS_NAN(tmp(j)).NEQV.IEEE_IS_NAN(vals(j))))THEN
    nbad=nbad+1
    dmax=MAX(dmax,ABS(tmp(j)-vals(j)))
  END IF
END DO
DEALLOCATE(tmp)
IF(nbad>0)WRITE(*,'(A,I8,A,I8,A,ES12.4)')'Rank ',oft_env%rank,': ',nbad, &
  ' shared values differ from owner, max abs difference ',dmax
IF(.NOT.seam%full)nbad=oft_mpi_sum(nbad)
IF(nbad>0)CALL oft_abort('Shared values inconsistent after release_local with sync=.FALSE. ('//caller//')', &
  'native_check_shared',__FILE__)
END SUBROUTINE native_check_shared
#endif
!---------------------------------------------------------------------------------
!> Borrow local values without a copy
!!
!! Returns a pointer to the vector's own storage (or to the range of block `iblock`).
!! See @ref oft_la_base::vector_borrow_local "vector_borrow_local" for the contract.
!---------------------------------------------------------------------------------
SUBROUTINE vec_borrow_local(self,array,mode,iblock)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
REAL(r8), POINTER, CONTIGUOUS, INTENT(inout) :: array(:) !< Local values (output)
INTEGER(i4), INTENT(in) :: mode !< Access mode (@ref oft_la_base::oft_access_read "OFT_ACCESS_READ", etc.)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block to borrow (optional)
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'vec_borrow_local',__FILE__)
IF(PRESENT(iblock))THEN
  IF(iblock<1.OR.iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_borrow_local',__FILE__)
END IF
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Zeroed accumulation storage, kept by the vector between accumulations
  IF(.NOT.self%acc_open)THEN
    IF(.NOT.ASSOCIATED(self%acc))ALLOCATE(self%acc(self%n))
    CALL nk_set(self%n,0.d0,self%acc)
    self%acc_open=.TRUE.
  END IF
  IF(PRESENT(iblock))THEN
    array=>self%acc(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n)
  ELSE
    array=>self%acc
  END IF
ELSE
  IF(PRESENT(iblock))THEN
    array=>self%v(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n)
  ELSE
    array=>self%v
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE vec_borrow_local
!---------------------------------------------------------------------------------
!> Release local values obtained with `borrow_local`
!!
!! Values were modified in place, so only the optional synchronization (`stitch(0)`)
!! is performed. See @ref oft_la_base::vector_release_local "vector_release_local"
!! for the contract.
!---------------------------------------------------------------------------------
SUBROUTINE vec_release_local(self,array,mode,sync,iblock)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
REAL(r8), POINTER, CONTIGUOUS, INTENT(inout) :: array(:) !< Borrowed local values (nullified on return)
INTEGER(i4), INTENT(in) :: mode !< Access mode used in the matching borrow
LOGICAL, INTENT(in) :: sync !< Synchronize shared values from their owner after writing?
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block used in the matching borrow (optional)
LOGICAL :: match
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: vals
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'vec_release_local',__FILE__)
IF(sync.AND.(mode==OFT_ACCESS_READ))CALL oft_abort('Synchronization requested for read-only access', &
  'vec_release_local',__FILE__)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  IF(.NOT.self%acc_open)CALL oft_abort('No accumulation in progress','vec_release_local',__FILE__)
  vals=>self%acc
ELSE
  vals=>self%v
END IF
!---Check that the array was borrowed from this vector (zero-size targets never match)
IF(PRESENT(iblock))THEN
  IF(iblock<1.OR.iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_release_local',__FILE__)
  match=(self%map(iblock)%n==0)
  IF(.NOT.match)match=ASSOCIATED(array,vals(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n))
ELSE
  match=(self%n==0)
  IF(.NOT.match)match=ASSOCIATED(array,vals)
END IF
IF(.NOT.match)CALL oft_abort('Array was not borrowed from this vector','vec_release_local',__FILE__)
NULLIFY(array)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Complete accumulation: sum contributions across copies once and add
  IF(sync)THEN
    CALL oft_global_stitch(self%stitch_info,self%acc,1)
    CALL nk_axpby(self%n,1.d0,1.d0,self%acc,self%v)
    self%acc_open=.FALSE.
  END IF
  DEBUG_STACK_POP
  RETURN
END IF
IF(sync.AND.(IAND(mode,OFT_ACCESS_WRITE)/=0))CALL self%stitch(0)
#ifdef OFT_DEBUG
IF((.NOT.sync).AND.(IAND(mode,OFT_ACCESS_WRITE)/=0))CALL native_check_shared(self%stitch_info,self%v,'vec_release_local')
#endif
DEBUG_STACK_POP
END SUBROUTINE vec_release_local
!---------------------------------------------------------------------------------
!> Copy values from a second vector
!!
!! self = a
!---------------------------------------------------------------------------------
SUBROUTINE vec_copy(self,a)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Source vector
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('Source does not match destination size.','vec_copy',__FILE__)
NULLIFY(atmp)
CALL a%borrow_local(atmp,OFT_ACCESS_READ)
CALL nk_copy(self%n,atmp,self%v)
CALL a%release_local(atmp,OFT_ACCESS_READ,.FALSE.)
DEBUG_STACK_POP
END SUBROUTINE vec_copy
!---------------------------------------------------------------------------------
!> Add a linear combination of vectors
!!
!! self = self + \f$ \sum_k \alpha_k a_k \f$
!!
!! All vectors are processed in a single pass over row blocks (see
!! @ref oft_native_kernels_vec_r8::k_maxpy "nk_maxpy"), so each block of `self` is loaded
!! once for all `n` updates. Terms are added in order of increasing `k`, as with `n`
!! successive calls to `add`; results agree to rounding, but need not be bitwise equal
!! (the compiler may contract the multiply-adds differently).
!! Non-native vectors in `a` are copied with `get_local` first.
!---------------------------------------------------------------------------------
SUBROUTINE vec_maxpy(self,alpha,a,n)
CLASS(oft_native_vector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Number of vectors
REAL(r8), INTENT(in) :: alpha(n) !< Scale factors [n]
TYPE(oft_vector_ptr), INTENT(inout) :: a(n) !< Vectors to add [n]
INTEGER(i4) :: k
LOGICAL :: copied(n)
TYPE(native_rvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','vec_maxpy',__FILE__)
  copied(k)=native_get_rvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
CALL nk_maxpy(self%n,n,alpha,vals,self%v)
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END SUBROUTINE vec_maxpy
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine vec_delete(self)
class(oft_native_vector), intent(inout) :: self !< Vector object
DEBUG_STACK_PUSH
if(associated(self%v))deallocate(self%v)
if(associated(self%local_tmp))deallocate(self%local_tmp)
IF(ASSOCIATED(self%acc))DEALLOCATE(self%acc)
self%acc_open=.FALSE.
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
INTEGER(i4) :: i1,i2
DEBUG_STACK_PUSH
random_flag=.FALSE.
IF(PRESENT(random))random_flag=random
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_set',__FILE__)
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  IF(random_flag)THEN
    CALL oft_random_number(self%v(i1:i2),self%map(iblock)%n)
  ELSE
    CALL nk_set(self%map(iblock)%n,alpha,self%v(i1:i2))
  END IF
  CALL self%stitch(0)
ELSE
  IF(random_flag)THEN
    CALL oft_random_number(self%v,self%n)
    CALL self%stitch(0)
  ELSE
    CALL nk_set(self%n,alpha,self%v)
  END IF
END IF
DEBUG_STACK_POP
end subroutine cvec_set
!---------------------------------------------------------------------------------
!> Get values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
SUBROUTINE cvec_get_slice(self,array,iblock)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
COMPLEX(c8), POINTER, INTENT(inout) :: array(:) !< Slice values
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block to retrieve
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_get_slice',__FILE__)
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%nslice))
  ASSOCIATE(mp=>self%map(iblock))
    IF(mp%nslice>0)THEN
      CALL native_slice_plan(mp)
      CALL scatter_local(mp%slice_plan,self%v(mp%offset+1:mp%offset+mp%n),array(1:mp%nslice),OFT_SCATTER_INSERT)
    END IF
  END ASSOCIATE
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%nslice))
  IF(self%n==self%nslice)THEN
    CALL nk_copy(self%n,self%v,array)
  ELSE
    DO i=1,self%nblocks
      ASSOCIATE(mp=>self%map(i))
        IF(mp%nslice>0)THEN
          CALL native_slice_plan(mp)
          CALL scatter_local(mp%slice_plan,self%v(mp%offset+1:mp%offset+mp%n), &
            array(mp%soffset+1:mp%soffset+mp%nslice),OFT_SCATTER_INSERT)
        END IF
      END ASSOCIATE
    END DO
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE cvec_get_slice
!---------------------------------------------------------------------------------
!> Set/add values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine cvec_restore_slice(self,array,iblock,wait)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: wait !< Wait to perform global sync?
INTEGER(i4) :: i
logical :: do_wait
DEBUG_STACK_PUSH
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_restore_slice',__FILE__)
  ASSOCIATE(mp=>self%map(iblock))
    IF(mp%nslice>0)THEN
      CALL native_slice_plan(mp)
      CALL scatter_local(mp%slice_plan,array(1:mp%nslice),self%v(mp%offset+1:mp%offset+mp%n), &
        OFT_SCATTER_INSERT,reverse=.TRUE.)
    END IF
  END ASSOCIATE
ELSE
  IF(self%n==self%nslice)THEN
    CALL nk_copy(self%n,array,self%v)
  ELSE
    DO i=1,self%nblocks
      ASSOCIATE(mp=>self%map(i))
        IF(mp%nslice>0)THEN
          CALL native_slice_plan(mp)
          CALL scatter_local(mp%slice_plan,array(mp%soffset+1:mp%soffset+mp%nslice), &
            self%v(mp%offset+1:mp%offset+mp%n),OFT_SCATTER_INSERT,reverse=.TRUE.)
        END IF
      END ASSOCIATE
    END DO
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
INTEGER(i4) :: i1,i2
DEBUG_STACK_PUSH
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_get_local',__FILE__)
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%n))
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  CALL nk_copy(self%map(iblock)%n,self%v(i1:i2),array)
ELSE
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%n))
  CALL nk_copy(self%n,self%v,array)
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
LOGICAL, OPTIONAL, INTENT(in) :: wait !< Wait to perform global sync
logical :: do_add,do_wait
INTEGER(i4) :: i1,i2
COMPLEX(c8), ALLOCATABLE, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
do_add=.FALSE.
IF(PRESENT(add))do_add=add
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_restore_local',__FILE__)
  i1=self%map(iblock)%offset+1
  i2=self%map(iblock)%offset+self%map(iblock)%n
  IF(do_add)THEN
    IF(.NOT.ASSOCIATED(self%local_tmp))THEN
      ALLOCATE(self%local_tmp(self%n))
      CALL nk_set(self%n,(0.d0,0.d0),self%local_tmp)
    END IF
    CALL nk_copy(self%map(iblock)%n,array,self%local_tmp(i1:i2))
    IF(.NOT.do_wait)THEN
      CALL full_restore_local(self%local_tmp)
      DEALLOCATE(self%local_tmp)
    END IF
  ELSE
    CALL nk_copy(self%map(iblock)%n,array,self%v(i1:i2))
    IF(.NOT.do_wait)CALL self%stitch(0)
  END IF
ELSE
  IF(ASSOCIATED(self%local_tmp))CALL oft_abort("Existing partial restore found when restoring full vector", &
    "cvec_restore_local",__FILE__)
  IF(do_wait)CALL oft_abort("Wait not supported when restoring full vector","cvec_restore_local",__FILE__)
  IF(do_add)THEN
    !---Contributions must be summed before adding, which needs a modifiable copy
    ALLOCATE(atmp(self%n))
    CALL nk_copy(self%n,array,atmp)
    CALL full_restore_local(atmp)
    DEALLOCATE(atmp)
  ELSE
    CALL nk_copy(self%n,array,self%v)
    CALL self%stitch(0)
  END IF
END IF
DEBUG_STACK_POP
contains
!
subroutine full_restore_local(array_loc)
complex(c8), intent(inout) :: array_loc(:)
CALL oft_global_stitch(self%stitch_info,array_loc,1)
CALL nk_axpby(self%n,(1.d0,0.d0),(1.d0,0.d0),array_loc,self%v)
end subroutine full_restore_local
end subroutine cvec_restore_local
!---------------------------------------------------------------------------------
!> Add real vectors
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
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','cvec_add_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','cvec_add_vec',__FILE__)
  copied(2)=native_get_vals(b,btmp)
  CALL nk_axpbypcz(self%n,gamma,alpha,atmp,beta,btmp,self%v)
  IF(copied(2))DEALLOCATE(btmp)
else ! Add one vector to the source
  CALL nk_axpby(self%n,gamma,alpha,atmp,self%v)
end if
IF(copied(1))DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_add_vec
!---------------------------------------------------------------------------------
!> Add complex vectors
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
LOGICAL :: copied(2)
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('First source does not match destination size.','cvec_add_cvec',__FILE__)
copied(1)=native_get_vals(a,atmp)
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(self%n/=b%n)CALL oft_abort('Second source does not match destination size.','cvec_add_cvec',__FILE__)
  copied(2)=native_get_vals(b,btmp)
  CALL nk_axpbypcz(self%n,gamma,alpha,atmp,beta,btmp,self%v)
  IF(copied(2))DEALLOCATE(btmp)
else ! Add one vector to the source
  CALL nk_axpby(self%n,gamma,alpha,atmp,self%v)
end if
IF(copied(1))DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_add_cvec
!---------------------------------------------------------------------------------
!> Elementwise multiplication with a real vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine cvec_mult_vec(self,a,div_flag)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
LOGICAL :: div,copied
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','cvec_mult_vec',__FILE__)
copied=native_get_vals(a,atmp)
if(div)then
  CALL nk_div(self%n,atmp,self%v)
else
  CALL nk_mult(self%n,atmp,self%v)
end if
IF(copied)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_mult_vec
!---------------------------------------------------------------------------------
!> Elementwise multiplication with a complex vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine cvec_mult_cvec(self,a,div_flag)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
LOGICAL :: div,copied
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
div=.FALSE.
if(present(div_flag))div=div_flag
IF(self%n/=a%n)CALL oft_abort('Vector sizes do not match.','cvec_mult_cvec',__FILE__)
copied=native_get_vals(a,atmp)
if(div)then
  CALL nk_div(self%n,atmp,self%v)
else
  CALL nk_mult(self%n,atmp,self%v)
end if
IF(copied)DEALLOCATE(atmp)
DEBUG_STACK_POP
end subroutine cvec_mult_cvec
!---------------------------------------------------------------------------------
!> Scale vector by a scalar
!---------------------------------------------------------------------------------
subroutine cvec_scale(self,alpha)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8), intent(in) :: alpha !< Scale factor
DEBUG_STACK_PUSH
CALL nk_scale(self%n,alpha,self%v)
DEBUG_STACK_POP
end subroutine cvec_scale
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function cvec_dot_vec(self,a) result(dot)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
COMPLEX(c8) :: dot !< \f$ \sum_i \overline{self_i} a_i \f$
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_vector)
    if(self%n/=this%n)call oft_abort('Vector lengths do not match.','cvec_dot_vec',__FILE__)
    dot=(0.d0,0.d0)
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))dot=nk_owned_dot(self%n,self%v,this%v,nw,w)
    IF(.NOT.self%stitch_info%full)dot=oft_mpi_sum(dot)
  class default
    CALL oft_abort('"a" is not a native vector.','cvec_dot_vec',__FILE__)
end select
DEBUG_STACK_POP
end function cvec_dot_vec
!---------------------------------------------------------------------------------
!> Dot product with a complex vector
!---------------------------------------------------------------------------------
function cvec_dot_cvec(self,a) result(dot)
class(oft_native_cvector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
COMPLEX(c8) :: dot !< \f$ \sum_i \overline{self_i} a_i \f$
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
select type(this=>a)
  class is(oft_native_cvector)
    IF(self%n/=this%n)CALL oft_abort('Vector lengths do not match.','cvec_dot_cvec',__FILE__)
    dot=(0.d0,0.d0)
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))dot=nk_owned_dot(self%n,self%v,this%v,nw,w)
    IF(.NOT.self%stitch_info%full)dot=oft_mpi_sum(dot)
  class default
    CALL oft_abort('"a" is not a native vector.','cvec_dot_cvec',__FILE__)
end select
DEBUG_STACK_POP
end function cvec_dot_cvec
!---------------------------------------------------------------------------------
!> Dot products with an array of real vectors
!!
!! Single pass over row blocks with a single MPI reduction, see @ref vec_mdot_vec.
!---------------------------------------------------------------------------------
FUNCTION cvec_mdot_vec(self,a,n) RESULT(dots)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Length of vector array
TYPE(oft_vector_ptr), INTENT(inout) :: a(n) !< Array of vectors for dot product [n]
COMPLEX(c8) :: dots(n) !< \f$ \sum_i \overline{self_i} a(j)_i \f$
INTEGER(i4) :: k,nw
LOGICAL :: copied(n)
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
TYPE(native_rvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','cvec_mdot_vec',__FILE__)
  copied(k)=native_get_rvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
dots=(0.d0,0.d0)
IF(n>0)THEN
  IF(native_seam_owned(self%stitch_info,self%n,nw,w))CALL nk_owned_mdot(self%n,n,self%v,vals,nw,w,dots)
  IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
END IF
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END FUNCTION cvec_mdot_vec
!---------------------------------------------------------------------------------
!> Dot products with an array of complex vectors
!!
!! Single pass over row blocks with a single MPI reduction, see @ref vec_mdot_vec.
!---------------------------------------------------------------------------------
FUNCTION cvec_mdot_cvec(self,a,n) RESULT(dots)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Length of vector array
TYPE(oft_cvector_ptr), INTENT(inout) :: a(n) !< Array of vectors for dot product [n]
COMPLEX(c8) :: dots(n) !< \f$ \sum_i \overline{self_i} a(j)_i \f$
INTEGER(i4) :: k,nw
LOGICAL :: copied(n)
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
TYPE(native_cvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','cvec_mdot_cvec',__FILE__)
  copied(k)=native_get_cvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
dots=(0.d0,0.d0)
IF(n>0)THEN
  IF(native_seam_owned(self%stitch_info,self%n,nw,w))CALL nk_owned_mdot(self%n,n,self%v,vals,nw,w,dots)
  IF(.NOT.self%stitch_info%full)dots=oft_mpi_sum(dots,n)
END IF
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END FUNCTION cvec_mdot_cvec
!---------------------------------------------------------------------------------
!> Sum reduction over vector
!---------------------------------------------------------------------------------
function cvec_sum(self) result(sum)
class(oft_native_cvector), intent(inout) :: self !< Vector object
complex(c8) :: sum !< Sum of vector elements
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
sum=(0.d0,0.d0)
IF(native_seam_owned(self%stitch_info,self%n,nw,w))sum=nk_owned_sum(self%n,self%v,nw,w)
IF(.NOT.self%stitch_info%full)sum=oft_mpi_sum(sum)
DEBUG_STACK_POP
end function cvec_sum
!---------------------------------------------------------------------------------
!> Compute norm of vector
!---------------------------------------------------------------------------------
function cvec_norm(self,itype) result(norm)
class(oft_native_cvector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
REAL(r8) :: norm !< Specified norm of vector
INTEGER(i4) :: nw
INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: w
DEBUG_STACK_PUSH
SELECT CASE(itype)
  CASE(1)
    norm=0.d0
    IF(native_seam_owned(self%stitch_info,self%n,nw,w))norm=nk_owned_norm1(self%n,self%v,nw,w)
    IF(.NOT.self%stitch_info%full)norm=oft_mpi_sum(norm)
  CASE(2)
    norm=SQRT(clamp_nonneg(REAL(self%dot(self),r8))) ! Clamp negative round-off (safety net, NaN preserved)
  CASE(3)
    norm=oft_mpi_max(nk_norm_inf(self%n,self%v))
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
!> Borrow local values without a copy
!!
!! Returns a pointer to the vector's own storage (or to the range of block `iblock`).
!! See @ref oft_la_base::vector_borrow_local "vector_borrow_local" for the contract.
!---------------------------------------------------------------------------------
SUBROUTINE cvec_borrow_local(self,array,mode,iblock)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
COMPLEX(c8), POINTER, CONTIGUOUS, INTENT(inout) :: array(:) !< Local values (output)
INTEGER(i4), INTENT(in) :: mode !< Access mode (@ref oft_la_base::oft_access_read "OFT_ACCESS_READ", etc.)
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block to borrow (optional)
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'cvec_borrow_local',__FILE__)
IF(PRESENT(iblock))THEN
  IF(iblock<1.OR.iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_borrow_local',__FILE__)
END IF
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Zeroed accumulation storage, kept by the vector between accumulations
  IF(.NOT.self%acc_open)THEN
    IF(.NOT.ASSOCIATED(self%acc))ALLOCATE(self%acc(self%n))
    CALL nk_set(self%n,(0.d0,0.d0),self%acc)
    self%acc_open=.TRUE.
  END IF
  IF(PRESENT(iblock))THEN
    array=>self%acc(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n)
  ELSE
    array=>self%acc
  END IF
ELSE
  IF(PRESENT(iblock))THEN
    array=>self%v(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n)
  ELSE
    array=>self%v
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE cvec_borrow_local
!---------------------------------------------------------------------------------
!> Release local values obtained with `borrow_local`
!!
!! Values were modified in place, so only the optional synchronization (`stitch(0)`)
!! is performed. See @ref oft_la_base::vector_release_local "vector_release_local"
!! for the contract.
!---------------------------------------------------------------------------------
SUBROUTINE cvec_release_local(self,array,mode,sync,iblock)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
COMPLEX(c8), POINTER, CONTIGUOUS, INTENT(inout) :: array(:) !< Borrowed local values (nullified on return)
INTEGER(i4), INTENT(in) :: mode !< Access mode used in the matching borrow
LOGICAL, INTENT(in) :: sync !< Synchronize shared values from their owner after writing?
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Sub-block used in the matching borrow (optional)
LOGICAL :: match
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: vals
DEBUG_STACK_PUSH
IF(mode<OFT_ACCESS_READ.OR.mode>OFT_ACCESS_ACCUMULATE)CALL oft_abort('Invalid access mode', &
  'cvec_release_local',__FILE__)
IF(sync.AND.(mode==OFT_ACCESS_READ))CALL oft_abort('Synchronization requested for read-only access', &
  'cvec_release_local',__FILE__)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  IF(.NOT.self%acc_open)CALL oft_abort('No accumulation in progress','cvec_release_local',__FILE__)
  vals=>self%acc
ELSE
  vals=>self%v
END IF
!---Check that the array was borrowed from this vector (zero-size targets never match)
IF(PRESENT(iblock))THEN
  IF(iblock<1.OR.iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','cvec_release_local',__FILE__)
  match=(self%map(iblock)%n==0)
  IF(.NOT.match)match=ASSOCIATED(array,vals(self%map(iblock)%offset+1:self%map(iblock)%offset+self%map(iblock)%n))
ELSE
  match=(self%n==0)
  IF(.NOT.match)match=ASSOCIATED(array,vals)
END IF
IF(.NOT.match)CALL oft_abort('Array was not borrowed from this vector','cvec_release_local',__FILE__)
NULLIFY(array)
IF(mode==OFT_ACCESS_ACCUMULATE)THEN
  !---Complete accumulation: sum contributions across copies once and add
  IF(sync)THEN
    CALL oft_global_stitch(self%stitch_info,self%acc,1)
    CALL nk_axpby(self%n,(1.d0,0.d0),(1.d0,0.d0),self%acc,self%v)
    self%acc_open=.FALSE.
  END IF
  DEBUG_STACK_POP
  RETURN
END IF
IF(sync.AND.(IAND(mode,OFT_ACCESS_WRITE)/=0))CALL self%stitch(0)
#ifdef OFT_DEBUG
!---Real and imaginary parts are checked separately
IF((.NOT.sync).AND.(IAND(mode,OFT_ACCESS_WRITE)/=0))THEN
  CALL native_check_shared(self%stitch_info,REAL(self%v,r8),'cvec_release_local (real part)')
  CALL native_check_shared(self%stitch_info,AIMAG(self%v),'cvec_release_local (imaginary part)')
END IF
#endif
DEBUG_STACK_POP
END SUBROUTINE cvec_release_local
!---------------------------------------------------------------------------------
!> Copy values from a real vector
!!
!! self = a
!---------------------------------------------------------------------------------
SUBROUTINE cvec_copy_vec(self,a)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Source vector
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('Source does not match destination size.','cvec_copy_vec',__FILE__)
NULLIFY(atmp)
CALL a%borrow_local(atmp,OFT_ACCESS_READ)
CALL nk_copy(self%n,atmp,self%v)
CALL a%release_local(atmp,OFT_ACCESS_READ,.FALSE.)
DEBUG_STACK_POP
END SUBROUTINE cvec_copy_vec
!---------------------------------------------------------------------------------
!> Copy values from a second complex vector
!!
!! self = a
!---------------------------------------------------------------------------------
SUBROUTINE cvec_copy_cvec(self,a)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
CLASS(oft_cvector), TARGET, INTENT(inout) :: a !< Source vector
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
DEBUG_STACK_PUSH
IF(self%n/=a%n)CALL oft_abort('Source does not match destination size.','cvec_copy_cvec',__FILE__)
NULLIFY(atmp)
CALL a%borrow_local(atmp,OFT_ACCESS_READ)
CALL nk_copy(self%n,atmp,self%v)
CALL a%release_local(atmp,OFT_ACCESS_READ,.FALSE.)
DEBUG_STACK_POP
END SUBROUTINE cvec_copy_cvec
!---------------------------------------------------------------------------------
!> Add a linear combination of complex vectors
!!
!! self = self + \f$ \sum_k \alpha_k a_k \f$
!!
!! Single pass over row blocks, see @ref vec_maxpy.
!---------------------------------------------------------------------------------
SUBROUTINE cvec_maxpy(self,alpha,a,n)
CLASS(oft_native_cvector), INTENT(inout) :: self !< Vector object
INTEGER(i4), INTENT(in) :: n !< Number of vectors
COMPLEX(c8), INTENT(in) :: alpha(n) !< Scale factors [n]
TYPE(oft_cvector_ptr), INTENT(inout) :: a(n) !< Vectors to add [n]
INTEGER(i4) :: k
LOGICAL :: copied(n)
TYPE(native_cvals_ptr) :: vals(n)
DEBUG_STACK_PUSH
DO k=1,n
  IF(self%n/=a(k)%f%n)CALL oft_abort('Vector lengths do not match.','cvec_maxpy',__FILE__)
  copied(k)=native_get_cvals(a(k)%f,vals(k)%v) ! Non-native vectors are copied
END DO
CALL nk_maxpy(self%n,n,alpha,vals,self%v)
DO k=1,n
  IF(copied(k))DEALLOCATE(vals(k)%v)
END DO
DEBUG_STACK_POP
END SUBROUTINE cvec_maxpy
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine cvec_delete(self)
class(oft_native_cvector), intent(inout) :: self !< Vector object
DEBUG_STACK_PUSH
if(associated(self%v))deallocate(self%v)
if(associated(self%local_tmp))deallocate(self%local_tmp)
IF(ASSOCIATED(self%acc))DEALLOCATE(self%acc)
self%acc_open=.FALSE.
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
!> Create a native CRS non-zero pattern
!!
!! The new pattern takes ownership of `kr`, `lc` and (if present) `map`, whose
!! pointers are nullified on return. Its reference count is zero; it is retained by
!! @ref native_matrix_set_pattern or @ref native_pattern_retain.
!---------------------------------------------------------------------------------
SUBROUTINE native_pattern_create(pattern,nr,nc,kr,lc,map)
TYPE(oft_native_pattern), POINTER, INTENT(out) :: pattern !< New pattern
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: nc !< Number of columns
INTEGER(i4), POINTER, CONTIGUOUS, INTENT(inout) :: kr(:) !< Row pointer [nr+1] (ownership is transferred)
INTEGER(i4), POINTER, CONTIGUOUS, INTENT(inout) :: lc(:) !< Column list [nnz] (ownership is transferred)
TYPE(oft_matrix_map), POINTER, OPTIONAL, INTENT(inout) :: map(:,:) !< Block maps (ownership is transferred)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(kr))CALL oft_abort('Row pointer not allocated','native_pattern_create',__FILE__)
IF(SIZE(kr)/=nr+1)CALL oft_abort('Row pointer size does not match number of rows', &
  'native_pattern_create',__FILE__)
ALLOCATE(pattern)
pattern%nr=nr
pattern%nc=nc
pattern%nnz=kr(nr+1)-1
pattern%kr=>kr
pattern%lc=>lc
NULLIFY(kr,lc)
IF(PRESENT(map))THEN
  IF(ASSOCIATED(map))THEN
    pattern%map=>map
    pattern%ni=SIZE(map,DIM=1)
    pattern%nj=SIZE(map,DIM=2)
    NULLIFY(map)
  END IF
END IF
!$omp atomic
native_pattern_live=native_pattern_live+1
DEBUG_STACK_POP
END SUBROUTINE native_pattern_create
!---------------------------------------------------------------------------------
!> Take a reference to a native CRS non-zero pattern
!---------------------------------------------------------------------------------
SUBROUTINE native_pattern_retain(pattern)
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Pattern
IF(.NOT.ASSOCIATED(pattern))CALL oft_abort('Pattern not associated','native_pattern_retain',__FILE__)
!$omp atomic update
pattern%refcount=pattern%refcount+1
END SUBROUTINE native_pattern_retain
!---------------------------------------------------------------------------------
!> Release a reference to a native CRS non-zero pattern
!!
!! The pattern and everything it owns are freed when the last reference is
!! released. `pattern` is nullified on return in either case.
!---------------------------------------------------------------------------------
RECURSIVE SUBROUTINE native_pattern_release(pattern)
TYPE(oft_native_pattern), POINTER, INTENT(inout) :: pattern !< Pattern
INTEGER(i4) :: i,j,nleft
TYPE(native_extraction), POINTER :: ext
IF(.NOT.ASSOCIATED(pattern))RETURN
DEBUG_STACK_PUSH
!$omp atomic capture
pattern%refcount=pattern%refcount-1
nleft=pattern%refcount
!$omp end atomic
IF(nleft<=0)THEN
  IF(ASSOCIATED(pattern%kr))DEALLOCATE(pattern%kr)
  IF(ASSOCIATED(pattern%lc))DEALLOCATE(pattern%lc)
  IF(ASSOCIATED(pattern%lc_small))DEALLOCATE(pattern%lc_small)
  IF(ASSOCIATED(pattern%redind))DEALLOCATE(pattern%redind)
  DO WHILE(ASSOCIATED(pattern%extracts))
    ext=>pattern%extracts
    pattern%extracts=>ext%next
    CALL native_extraction_free(ext)
  END DO
  CALL native_csc_delete(pattern%csc)
  CALL native_diag_delete(pattern%diag_index)
  IF(ASSOCIATED(pattern%map))THEN
    DO i=1,SIZE(pattern%map,DIM=1)
      DO j=1,SIZE(pattern%map,DIM=2)
        IF(ASSOCIATED(pattern%map(i,j)%lc_map))DEALLOCATE(pattern%map(i,j)%lc_map)
        IF(ASSOCIATED(pattern%map(i,j)%ext))DEALLOCATE(pattern%map(i,j)%ext)
      END DO
    END DO
    DEALLOCATE(pattern%map)
  END IF
  IF(ASSOCIATED(pattern%linkage))THEN
    CALL seam_free_plan(pattern%linkage)
    IF(ASSOCIATED(pattern%linkage%kle))DEALLOCATE(pattern%linkage%kle)
    IF(ASSOCIATED(pattern%linkage%lbe))DEALLOCATE(pattern%linkage%lbe)
    IF(ASSOCIATED(pattern%linkage%lle))DEALLOCATE(pattern%linkage%lle)
    DEALLOCATE(pattern%linkage)
  END IF
  DEALLOCATE(pattern)
  !$omp atomic
  native_pattern_live=native_pattern_live-1
END IF
NULLIFY(pattern)
DEBUG_STACK_POP
END SUBROUTINE native_pattern_release
!---------------------------------------------------------------------------------
!> Point the pattern aliases of a matrix at the components of its pattern
!!
!! Shared by the real and complex CRS matrix classes; the arguments are the alias
!! components of the matrix.
!---------------------------------------------------------------------------------
SUBROUTINE native_pattern_alias(pattern,kr,lc,map,linkage,redind,nred,nnz)
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Pattern (may be null)
INTEGER(i4), POINTER, CONTIGUOUS, INTENT(inout) :: kr(:) !< Row pointer alias
INTEGER(i4), POINTER, CONTIGUOUS, INTENT(inout) :: lc(:) !< Column list alias
TYPE(oft_matrix_map), POINTER, INTENT(inout) :: map(:,:) !< Block map alias
TYPE(oft_seam), POINTER, INTENT(inout) :: linkage !< Linkage alias
INTEGER(i4), POINTER, CONTIGUOUS, INTENT(inout) :: redind(:) !< Redundant row alias
INTEGER(i4), INTENT(out) :: nred !< Number of redundant rows
INTEGER(i4), INTENT(out) :: nnz !< Number of non-zeros
IF(ASSOCIATED(pattern))THEN
  kr=>pattern%kr
  lc=>pattern%lc
  map=>pattern%map
  linkage=>pattern%linkage
  redind=>pattern%redind
  nred=pattern%nred
  nnz=pattern%nnz
ELSE
  NULLIFY(kr,lc,map,linkage,redind)
  nred=0
  nnz=0
END IF
END SUBROUTINE native_pattern_alias
!---------------------------------------------------------------------------------
!> Use a non-zero pattern for a real CRS matrix
!!
!! Any previous pattern is released, a reference to `pattern` is taken and the
!! pattern aliases (`kr`, `lc`, `map`, etc.) are updated. Matrix values are not
!! allocated or changed.
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_set_pattern_real(self,pattern)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Pattern to use
TYPE(oft_native_pattern), POINTER :: old
IF(pattern%nr/=self%nr)CALL oft_abort('Pattern rows do not match matrix','native_matrix_set_pattern',__FILE__)
IF(pattern%nc/=self%nc)CALL oft_abort('Pattern columns do not match matrix','native_matrix_set_pattern',__FILE__)
CALL native_pattern_retain(pattern)
old=>self%pattern
self%pattern=>pattern
CALL native_pattern_release(old)
self%full_current=.FALSE.
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
END SUBROUTINE native_matrix_set_pattern_real
!---------------------------------------------------------------------------------
!> Use a non-zero pattern for a complex CRS matrix, see @ref native_matrix_set_pattern_real
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_set_pattern_comp(self,pattern)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Pattern to use
TYPE(oft_native_pattern), POINTER :: old
IF(pattern%nr/=self%nr)CALL oft_abort('Pattern rows do not match matrix','native_matrix_set_pattern',__FILE__)
IF(pattern%nc/=self%nc)CALL oft_abort('Pattern columns do not match matrix','native_matrix_set_pattern',__FILE__)
CALL native_pattern_retain(pattern)
old=>self%pattern
self%pattern=>pattern
CALL native_pattern_release(old)
self%full_current=.FALSE.
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
END SUBROUTINE native_matrix_set_pattern_comp
!---------------------------------------------------------------------------------
!> Create a pattern from the `kr`, `lc` and `map` set directly on a real CRS matrix
!!
!! For code that builds the CRS structure of a matrix itself: after setting `nr`,
!! `nc`, `kr` and `lc` (and optionally `map`) on a matrix without a pattern, this
!! moves them into a new pattern used by the matrix. Native matrix methods call this
!! automatically for matrices without a pattern.
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_adopt_pattern_real(self)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix
TYPE(oft_native_pattern), POINTER :: pattern
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr,lc
TYPE(oft_matrix_map), POINTER, DIMENSION(:,:) :: map
IF(ASSOCIATED(self%pattern))CALL oft_abort('Matrix already has a pattern','native_matrix_adopt_pattern',__FILE__)
kr=>self%kr
lc=>self%lc
map=>self%map
CALL native_pattern_create(pattern,self%nr,self%nc,kr,lc,map)
CALL native_matrix_set_pattern_real(self,pattern)
END SUBROUTINE native_matrix_adopt_pattern_real
!---------------------------------------------------------------------------------
!> Create a pattern from the `kr`, `lc` and `map` set directly on a complex CRS matrix,
!! see @ref native_matrix_adopt_pattern_real
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_adopt_pattern_comp(self)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix
TYPE(oft_native_pattern), POINTER :: pattern
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr,lc
TYPE(oft_matrix_map), POINTER, DIMENSION(:,:) :: map
IF(ASSOCIATED(self%pattern))CALL oft_abort('Matrix already has a pattern','native_matrix_adopt_pattern',__FILE__)
kr=>self%kr
lc=>self%lc
map=>self%map
CALL native_pattern_create(pattern,self%nr,self%nc,kr,lc,map)
CALL native_matrix_set_pattern_comp(self,pattern)
END SUBROUTINE native_matrix_adopt_pattern_comp
!---------------------------------------------------------------------------------
!> Get the dimensions, block mappings and pattern of a native CRS matrix (real or complex)
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_layout(source,pattern,nr,nc,nrg,ncg,nrslice,ncslice,ni,nj,i_map,j_map,force_local)
CLASS(*), TARGET, INTENT(inout) :: source !< Native real or complex CRS matrix
TYPE(oft_native_pattern), POINTER, INTENT(out) :: pattern !< Pattern of `source`
INTEGER(i4), INTENT(out) :: nr !< Local number of rows
INTEGER(i4), INTENT(out) :: nc !< Local number of columns
INTEGER(i8), INTENT(out) :: nrg !< Global number of rows
INTEGER(i8), INTENT(out) :: ncg !< Global number of columns
INTEGER(i4), INTENT(out) :: nrslice !< Number of owned rows
INTEGER(i4), INTENT(out) :: ncslice !< Number of owned columns
INTEGER(i4), INTENT(out) :: ni !< Number of row blocks
INTEGER(i4), INTENT(out) :: nj !< Number of column blocks
TYPE(oft_map), POINTER, INTENT(out) :: i_map(:) !< Row block mapping
TYPE(oft_map), POINTER, INTENT(out) :: j_map(:) !< Column block mapping
LOGICAL, INTENT(out) :: force_local !< Local-only flag
SELECT TYPE(source)
CLASS IS(oft_native_matrix)
  IF(.NOT.ASSOCIATED(source%pattern))CALL native_matrix_adopt_pattern_real(source)
  pattern=>source%pattern
  nr=source%nr; nc=source%nc; nrg=source%nrg; ncg=source%ncg
  nrslice=source%nrslice; ncslice=source%ncslice; ni=source%ni; nj=source%nj
  i_map=>source%i_map; j_map=>source%j_map; force_local=source%force_local
CLASS IS(oft_native_cmatrix)
  IF(.NOT.ASSOCIATED(source%pattern))CALL native_matrix_adopt_pattern_comp(source)
  pattern=>source%pattern
  nr=source%nr; nc=source%nc; nrg=source%nrg; ncg=source%ncg
  nrslice=source%nrslice; ncslice=source%ncslice; ni=source%ni; nj=source%nj
  i_map=>source%i_map; j_map=>source%j_map; force_local=source%force_local
CLASS DEFAULT
  CALL oft_abort('Source is not a native CRS matrix','native_matrix_layout',__FILE__)
END SELECT
END SUBROUTINE native_matrix_layout
!---------------------------------------------------------------------------------
!> Create a real CRS matrix sharing the non-zero pattern of an existing matrix
!!
!! `source` may be a real or complex native CRS matrix. The new matrix has the same
!! dimensions and block mappings, references the same @ref oft_native_pattern (so
!! the CRS structure, block maps, cached indices and linkage are not duplicated or
!! rebuilt) and has its own values, initialized to zero.
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_share_real(source,new)
CLASS(*), TARGET, INTENT(inout) :: source !< Matrix whose pattern is shared
CLASS(oft_matrix), POINTER, INTENT(out) :: new !< New matrix
TYPE(oft_native_pattern), POINTER :: pattern
TYPE(oft_map), POINTER :: i_map(:),j_map(:)
DEBUG_STACK_PUSH
ALLOCATE(oft_native_matrix::new)
SELECT TYPE(this=>new)
CLASS IS(oft_native_matrix)
  CALL native_matrix_layout(source,pattern,this%nr,this%nc,this%nrg,this%ncg,this%nrslice,this%ncslice, &
    this%ni,this%nj,i_map,j_map,this%force_local)
  IF(ASSOCIATED(pattern%map))THEN
    ALLOCATE(this%i_map(this%ni),this%j_map(this%nj))
    this%i_map=i_map
    this%j_map=j_map
  END IF
  CALL native_matrix_set_pattern_real(this,pattern)
  ALLOCATE(this%M(this%nnz))
  CALL nk_set(this%nnz,0.d0,this%M)
END SELECT
DEBUG_STACK_POP
END SUBROUTINE native_matrix_share_real
!---------------------------------------------------------------------------------
!> Create a complex CRS matrix sharing the non-zero pattern of an existing matrix
!!
!! See @ref native_matrix_share_real.
!---------------------------------------------------------------------------------
SUBROUTINE native_matrix_share_comp(source,new)
CLASS(*), TARGET, INTENT(inout) :: source !< Matrix whose pattern is shared
CLASS(oft_cmatrix), POINTER, INTENT(out) :: new !< New matrix
TYPE(oft_native_pattern), POINTER :: pattern
TYPE(oft_map), POINTER :: i_map(:),j_map(:)
DEBUG_STACK_PUSH
ALLOCATE(oft_native_cmatrix::new)
SELECT TYPE(this=>new)
CLASS IS(oft_native_cmatrix)
  CALL native_matrix_layout(source,pattern,this%nr,this%nc,this%nrg,this%ncg,this%nrslice,this%ncslice, &
    this%ni,this%nj,i_map,j_map,this%force_local)
  IF(ASSOCIATED(pattern%map))THEN
    ALLOCATE(this%i_map(this%ni),this%j_map(this%nj))
    this%i_map=i_map
    this%j_map=j_map
  END IF
  CALL native_matrix_set_pattern_comp(this,pattern)
  ALLOCATE(this%M(this%nnz))
  CALL nk_set(this%nnz,(0.d0,0.d0),this%M)
END SELECT
DEBUG_STACK_POP
END SUBROUTINE native_matrix_share_comp
!---------------------------------------------------------------------------------
!> Build transposed (CSC) index for a CRS non-zero pattern, if not already built
!!
!! The index is owned by the pattern, which does not change once created. Callers that may run concurrently on the same matrix must serialize calls
!! (eg. using a named critical section).
!---------------------------------------------------------------------------------
SUBROUTINE native_csc_build(csc,kr,lc,nr,nc)
TYPE(native_csc_index), INTENT(inout) :: csc !< Transposed index
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: nc !< Number of columns
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CRS row pointer [nr+1]
INTEGER(i4), INTENT(in) :: lc(:) !< CRS column list [nnz]
INTEGER(i4) :: i,j,k,nnz
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: knext
nnz=kr(nr+1)-1
IF(ASSOCIATED(csc%kc))RETURN
ALLOCATE(csc%kc(nc+1),csc%lr(nnz),csc%lm(nnz))
!---Count entries in each column
csc%kc=0
DO j=1,nnz
  csc%kc(lc(j)+1)=csc%kc(lc(j)+1)+1
END DO
csc%kc(1)=1
DO k=1,nc
  csc%kc(k+1)=csc%kc(k+1)+csc%kc(k)
END DO
!---Fill entries, traversing rows in order so rows are increasing within each column
ALLOCATE(knext(nc))
knext=csc%kc(1:nc)
DO i=1,nr
  DO j=kr(i),kr(i+1)-1
    k=lc(j)
    csc%lr(knext(k))=i
    csc%lm(knext(k))=j
    knext(k)=knext(k)+1
  END DO
END DO
DEALLOCATE(knext)
END SUBROUTINE native_csc_build
!---------------------------------------------------------------------------------
!> Release transposed (CSC) index
!---------------------------------------------------------------------------------
SUBROUTINE native_csc_delete(csc)
TYPE(native_csc_index), INTENT(inout) :: csc !< Transposed index
IF(ASSOCIATED(csc%kc))DEALLOCATE(csc%kc)
IF(ASSOCIATED(csc%lr))DEALLOCATE(csc%lr)
IF(ASSOCIATED(csc%lm))DEALLOCATE(csc%lm)
END SUBROUTINE native_csc_delete
!---------------------------------------------------------------------------------
!> Build index of diagonal entry locations for a CRS non-zero pattern, if not already built
!!
!! Rows without a diagonal entry are marked with location 0. The index is owned by the
!! pattern, which does not change once created. Callers that may run concurrently on
!! the same pattern must serialize calls.
!---------------------------------------------------------------------------------
SUBROUTINE native_diag_build(dind,kr,lc,nr)
TYPE(native_diag_index), INTENT(inout) :: dind !< Diagonal index
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CRS row pointer [nr+1]
INTEGER(i4), INTENT(in) :: lc(:) !< CRS column list [nnz]
INTEGER(i4) :: i,jp,jn,jr
IF(ASSOCIATED(dind%loc))RETURN
ALLOCATE(dind%loc(nr))
!$omp parallel do private(jp,jn,jr) if(nr>OFT_OMP_VTHRESH)
DO i=1,nr
  jp=kr(i)
  jn=kr(i+1)-jp
  jr=search_array(i,lc(jp:jp+jn-1),jn)
  IF(jr>0)jr=jr+jp-1
  dind%loc(i)=jr
END DO
END SUBROUTINE native_diag_build
!---------------------------------------------------------------------------------
!> Release index of diagonal entry locations
!---------------------------------------------------------------------------------
SUBROUTINE native_diag_delete(dind)
TYPE(native_diag_index), INTENT(inout) :: dind !< Diagonal index
IF(ASSOCIATED(dind%loc))DEALLOCATE(dind%loc)
END SUBROUTINE native_diag_delete
!---------------------------------------------------------------------------------
!> Compute matrix-vector product
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine mat_apply_vec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','mat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','mat_apply_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
IF(ASSOCIATED(self%pattern%lc_small))THEN
  CALL nk_csr_matvec(self%nr,self%kr,self%pattern%lc_small,self%M,atmp,btmp)
ELSE
  CALL nk_csr_matvec(self%nr,self%kr,self%lc,self%M,atmp,btmp)
END IF
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_apply_vec
!---------------------------------------------------------------------------------
!> Compute matrix-vector product (complex)
!!
!! b = self * a
!---------------------------------------------------------------------------------
SUBROUTINE mat_apply_cvec(self,a,b)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_cvector), TARGET, INTENT(inout) :: a !< Source vector
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(b%n/=self%nr)CALL oft_abort('Row mismatch','mat_apply_cvec',__FILE__)
IF(a%n/=self%nc)CALL oft_abort('Col mismatch','mat_apply_cvec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
IF(ASSOCIATED(self%pattern%lc_small))THEN
  CALL nk_csr_matvec(self%nr,self%kr,self%pattern%lc_small,self%M,atmp,btmp)
ELSE
  CALL nk_csr_matvec(self%nr,self%kr,self%lc,self%M,atmp,btmp)
END IF
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)CALL b%stitch(1)
DEBUG_STACK_POP
END SUBROUTINE mat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine mat_applyt_vec(self,a,b)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','mat_applyt_vec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','mat_applyt_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!$omp critical (oft_native_csc_build)
CALL native_csc_build(self%pattern%csc,self%kr,self%lc,self%nr,self%nc)
!$omp end critical (oft_native_csc_build)
CALL nk_csc_matvec(self%nc,self%pattern%csc%kc,self%pattern%csc%lr,self%pattern%csc%lm,self%M,atmp,btmp)
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
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
LOGICAL :: copied(2)
complex(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Col mismatch','mat_applyt_cvec',__FILE__)
if(a%n/=self%nr)call oft_abort('Row mismatch','mat_applyt_cvec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!$omp critical (oft_native_csc_build)
CALL native_csc_build(self%pattern%csc,self%kr,self%lc,self%nr,self%nc)
!$omp end critical (oft_native_csc_build)
CALL nk_csc_matvec(self%nc,self%pattern%csc%kc,self%pattern%csc%lr,self%pattern%csc%lm,self%M,atmp,btmp)
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)call b%stitch(1)
DEBUG_STACK_POP
end subroutine mat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ANY(i_inds>self%nr))call oft_abort('Index exceeds graph','mat_set_values',__FILE__)
if(ANY(j_inds>self%nc))call oft_abort('Index exceeds graph','mat_set_values',__FILE__)
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','mat_set_values',__FILE__)
CALL nk_set_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
end subroutine mat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_add_values',__FILE__)
self%full_current=.FALSE.
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock,loc_cache)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','mat_add_values',__FILE__)
CALL nk_add_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
end subroutine mat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix atomically (for concurrent assembly)
!---------------------------------------------------------------------------------
subroutine mat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','mat_atomic_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','mat_atomic_add_values',__FILE__)
self%full_current=.FALSE.
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock,loc_cache)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','mat_atomic_add_values',__FILE__)
CALL nk_atomic_add_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
end subroutine mat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine mat_assemble(self,diag)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_vector)
    !---Get diagonal matrix values (rows without a diagonal entry get zero)
    !$omp critical (oft_native_diag_build)
    CALL native_diag_build(self%pattern%diag_index,self%kr,self%lc,self%nr)
    !$omp end critical (oft_native_diag_build)
    CALL nk_csr_get_diag(self%nr,self%pattern%diag_index%loc,self%M,this%v)
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a vector object.','mat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add(0.d0,1.d0,self%D)
end if
!---Common assembly tasks
!$omp critical (oft_native_lc_small_build)
CALL native_lc_small_build(self%pattern%lc_small,self%lc,self%nnz,self%nc)
!$omp end critical (oft_native_lc_small_build)
DEBUG_STACK_POP
end subroutine mat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine mat_zero(self)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
DEBUG_STACK_PUSH
self%full_current=.FALSE.
CALL nk_set(self%nnz,0.d0,self%M)
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
INTEGER(i4) :: ioff
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
ioff=0
IF(PRESENT(iblock))ioff=self%i_map(iblock)%offset
CALL nk_csr_zero_rows(nrows,irows,ioff,self%kr,self%lc,zero_diag,self%M)
DEBUG_STACK_POP
end subroutine mat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine mat_zero_cols(self,ncols,icols,jblock,keep_diag)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: coff
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
coff=0
IF(PRESENT(jblock))coff=self%j_map(jblock)%offset
CALL nk_csr_zero_cols(self%nr,ncols,icols,coff,self%kr,self%lc,zero_diag,self%M)
DEBUG_STACK_POP
end subroutine mat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!!
!! Values are freed and the reference to the non-zero pattern is released, so
!! the pattern is only freed if no other matrix (or holder) uses it.
!---------------------------------------------------------------------------------
subroutine mat_delete(self)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix object
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ASSOCIATED(self%M))deallocate(self%M)
if(ASSOCIATED(self%Md))deallocate(self%Md)
IF(ASSOCIATED(self%Mfull))DEALLOCATE(self%Mfull)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
!---Block mappings belong to the matrix, block maps to the pattern
IF(ASSOCIATED(self%map))DEALLOCATE(self%i_map,self%j_map)
!---Matrices built without a pattern own their structure
IF((.NOT.ASSOCIATED(self%pattern)).AND.ASSOCIATED(self%kr))CALL native_matrix_adopt_pattern(self)
CALL native_pattern_release(self%pattern)
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
IF(ASSOCIATED(self%color))DEALLOCATE(self%color)
!---Cell overlap builder is not owned by the matrix
NULLIFY(self%cell_ovl)
DEBUG_STACK_POP
end subroutine mat_delete
!---------------------------------------------------------------------------------
!> Save matrix in CSR format
!---------------------------------------------------------------------------------
subroutine mat_save(self,filename, path, bc_flags, nfields)
class(oft_native_matrix), intent(inout) :: self !< Matrix object
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:) :: bc_flags
INTEGER, INTENT(IN), OPTIONAL :: nfields
integer(i4), dimension(:), allocatable :: bit_flags
integer(i4) :: i, j
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: lg !< global index of each element
DEBUG_STACK_PUSH
IF(oft_env%nprocs>1)CALL oft_abort('Not supported for nprocs>1','mat_save',__FILE__)
ALLOCATE(lg(self%nr))
IF(PRESENT(bc_flags))THEN
  ALLOCATE(bit_flags(size(bc_flags)))
  bit_flags=0
  bit_flags=MERGE(1,0,bc_flags)
END IF
DO i=1,self%ni
  DO j=1, self%i_map(i)%n
    lg(self%i_map(i)%offset+j)=self%i_map(i)%offsetg+self%i_map(i)%lge(j)
  END DO
END DO
!---
CALL hdf5_create_file(filename)
!---Save matrix
CALL hdf5_create_group(filename, path)
CALL hdf5_write(self%nr,filename,TRIM(path)//'/nr')
CALL hdf5_write(self%nc,filename,TRIM(path)//'/nc')
CALL hdf5_write(self%kr,filename,TRIM(path)//'/kr')
CALL hdf5_write(self%lc,filename,TRIM(path)//'/lc')
CALL hdf5_write(self%nrg,filename,TRIM(path)//'/nrg')
CALL hdf5_write(self%ncg,filename,TRIM(path)//'/ncg')
CALL hdf5_write(self%M,filename,TRIM(path)//'/M')
CALL hdf5_write(lg,filename,TRIM(path)//'/lg')
If(PRESENT(bc_flags))CALL hdf5_write(bit_flags,filename,TRIM(path)//'/bc_flags')
IF(PRESENT(nfields))CALL hdf5_write(nfields,filename,TRIM(path)//'/nfields')
DEALLOCATE(lg)
DEBUG_STACK_POP
end subroutine mat_save
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
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lle_tmp,ncon,nrecv,bemap
integer(i8) :: last
integer(i8), allocatable, dimension(:) :: lge,sort_list
integer(i8), pointer, dimension(:,:) :: letmp,leout
logical :: periodic
logical, allocatable, dimension(:) :: glob_inc
logical, allocatable, dimension(:,:) :: skips
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!---Linkage is built once per pattern and shared by all matrices using it
IF(ASSOCIATED(self%pattern%linkage))THEN
  CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
  RETURN
END IF
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
  ALLOCATE(self%pattern%linkage)
  self%pattern%linkage%full=.TRUE.
  self%pattern%linkage%skip=.TRUE.
  CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
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
!---Build into the pattern; the matrix aliases point at the pattern components
ALLOCATE(self%pattern%linkage)
self%linkage=>self%pattern%linkage
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
!---Flag connections for each boundary row (skips(:,0) is used for non-boundary rows)
ALLOCATE(bemap(self%nr))
bemap=0
DO i=1,uv%stitch_info%nbe
  bemap(uv%stitch_info%lbe(i))=i
END DO
ALLOCATE(skips(0:self%linkage%nproc_con,0:uv%stitch_info%nbe))
skips=.TRUE.
DO i=1,uv%stitch_info%nbe
  DO m=0,self%linkage%nproc_con
    jp=uv%stitch_info%kle(m)
    jn=uv%stitch_info%kle(m+1)-jp
    IF(jn>0)THEN
      jk=search_array(i,lle_tmp(jp:jp+jn-1),jn)
      IF(jk/=0)skips(m,i)=.FALSE.
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
    j = self%kr(k+1)-self%kr(k)
    !$omp atomic capture
    m=nbetmp
    nbetmp=nbetmp+j
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
    IF(skips(m,bemap(letmp(1,k))).OR.skips(m,bemap(letmp(2,k))))CYCLE
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
      IF(skips(m,bemap(letmp(1,k))).OR.skips(m,bemap(letmp(2,k))))CYCLE
      ncon(m)=ncon(m)+1
      lesend(m)%le(:,ncon(m)) = lge(letmp(:,k))
      lrmaps(m)%le(ncon(m)) = k
    END DO
  END DO
END IF
#endif
DEALLOCATE(letmp,skips,bemap,lge)
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
    !---Count matches first, then store them (with periodicity, an entry may match many remote copies)
    ncon(j)=0
    DO k=1,2
      IF(k==2)THEN
        DEALLOCATE(lrmaps(j)%linktmp)
        ALLOCATE(lrmaps(j)%linktmp(2,MAX(1,ncon(j))))
        ncon(j)=0
      END IF
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
              IF(k==2)lrmaps(j)%linktmp(:,ncon(j))=(/mm,m/)
              neel=neel+1
            ELSE IF(leout(2,m)>letmp(2,mm))THEN
              EXIT
            END IF
          END IF
          m=m+1
          IF(m>neout)EXIT
        end do
      end do
    END DO
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
  !---Count matches first, then store them (with periodicity, an entry may match many local copies)
  DO k=1,2
    IF(k==2)THEN
      DEALLOCATE(lrmaps(0)%linktmp)
      ALLOCATE(lrmaps(0)%linktmp(2,MAX(1,ncon(0))))
      ncon(0)=0
      neel=0
    END IF
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
            IF(k==2)lrmaps(0)%linktmp(:,ncon(0))=(/mm,m/)
            neel=neel+1
          ELSE IF(leout(2,m)>leout(2,mm))THEN
            EXIT
          END IF
        END IF
        m=m+1
        IF(m>neout)EXIT
      end do
    end do
  END DO
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
  self%pattern%nred=COUNT(glob_inc)
  ALLOCATE(self%pattern%redind(self%pattern%nred))
  self%pattern%nred=0
  DO i=1,uv%stitch_info%nbe
    IF(glob_inc(uv%stitch_info%lbe(i)))THEN
      self%pattern%nred=self%pattern%nred+1
      self%pattern%redind(self%pattern%nred)=uv%stitch_info%lbe(i)
    END IF
  END DO
  DEALLOCATE(glob_inc)
END IF
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
! IF(oft_debug_print(3))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine native_matrix_setup_full
!---------------------------------------------------------------------------------
!> Update full local matrix with values from other processors
!---------------------------------------------------------------------------------
SUBROUTINE matrix_update_slice(self)
CLASS(oft_native_matrix), INTENT(inout) :: self !< Matrix object
IF(self%full_current)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(3))WRITE(*,*)'  Get local slice'
IF(.NOT.ASSOCIATED(self%Mfull))ALLOCATE(self%Mfull(self%nnz))
IF(.NOT.ASSOCIATED(self%pattern))CALL oft_abort('Linkage not built (see native_matrix_setup_full)', &
  'matrix_update_slice',__FILE__)
IF(.NOT.ASSOCIATED(self%pattern%linkage))CALL oft_abort('Linkage not built (see native_matrix_setup_full)', &
  'matrix_update_slice',__FILE__)
!---Linkage may have been built through another matrix sharing the pattern
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
!---Get values
CALL nk_copy(self%nnz,self%M,self%Mfull)
CALL oft_global_stitch(self%pattern%linkage,self%Mfull,1)
!---Replace redundant rows with identity rows
IF(self%pattern%nred>0)CALL nk_csr_identity_rows(self%pattern%nred,self%pattern%redind,self%kr,self%lc,self%Mfull)
self%full_current=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE matrix_update_slice
!---------------------------------------------------------------------------------
!> Setup local submatrix from parent matrix
!!
!! The submatrix holds the entries of the parent's full local matrix (see
!! @ref native_matrix_setup_full) in a set of rows and a set of columns, numbered in
!! the order given:
!! - `rows` (and optionally `cols`, default equal to `rows`): parent local indices;
!! - `part`: rows and columns `part` (kept as a pointer);
!! - `slice`: owned rows and columns of one block, in slice order;
!! - none of these: all owned rows and columns, in slice order.
!!
!! The symbolic extraction is cached by the parent's pattern (@ref native_extraction_get).
!---------------------------------------------------------------------------------
SUBROUTINE submatrix_setup(self,mat,u,slice,part,rows,cols)
class(oft_native_submatrix), intent(inout) :: self
class(oft_matrix), intent(inout) :: mat !< Parent matrix
class(oft_vector), intent(inout) :: u !< Constructed vector for parent row space
integer(i4), optional, intent(in) :: slice !< Index of field slice to use (optional)
integer(i4), optional, target, intent(in) :: part(:) !< Array of elements to define submatrix [:] (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: rows(:) !< Parent rows of the submatrix (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: cols(:) !< Parent columns of the submatrix (optional, requires `rows`)
class(oft_native_matrix), pointer :: parent
INTEGER(i4) :: i,j,nr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: rlist
DEBUG_STACK_PUSH
!---
IF(COUNT([PRESENT(part),PRESENT(slice),PRESENT(rows)])>1)CALL oft_abort('Only one of "slice", "part" and "rows" '// &
  'may be specified','submatrix_setup',__FILE__)
IF(PRESENT(cols).AND.(.NOT.PRESENT(rows)))CALL oft_abort('"cols" requires "rows"','submatrix_setup',__FILE__)
IF(PRESENT(part))self%part=>part
IF(PRESENT(slice))self%slice=slice
IF(.NOT.native_matrix_cast(parent,mat))CALL oft_abort('Native matrix required', &
  'submatrix_setup',__FILE__)
self%parent=>parent
!---Parent rows of the submatrix
IF(PRESENT(rows))THEN
  nr=SIZE(rows)
  ALLOCATE(rlist(nr))
  rlist=rows
ELSE IF(ASSOCIATED(self%part))THEN
  nr=SIZE(self%part)
  ALLOCATE(rlist(nr))
  rlist=self%part
ELSE IF(self%slice>0)THEN
  ASSOCIATE(mp=>parent%i_map(self%slice))
    nr=mp%nslice
    ALLOCATE(rlist(nr))
    DO j=1,nr
      rlist(j)=mp%offset+mp%slice(j)
    END DO
  END ASSOCIATE
ELSE
  nr=parent%nrslice
  ALLOCATE(rlist(nr))
  DO i=1,parent%ni
    ASSOCIATE(mp=>parent%i_map(i))
      DO j=1,mp%nslice
        rlist(mp%soffset+j)=mp%offset+mp%slice(j)
      END DO
    END ASSOCIATE
  END DO
END IF
!---Full local matrix of the parent and symbolic extraction (cached by the parent pattern)
CALL native_matrix_setup_full(self%parent,u)
IF(PRESENT(cols))THEN
  self%ext=>native_extraction_get(parent%pattern,nr,rlist,SIZE(cols),cols)
ELSE
  self%ext=>native_extraction_get(parent%pattern,nr,rlist,nr,rlist)
END IF
DEALLOCATE(rlist)
CALL native_pattern_retain(parent%pattern)
self%parent_pattern=>parent%pattern
self%nr=self%ext%nr
self%nc=self%ext%nc
self%nrg=self%nr
self%ncg=self%nc
self%ni=1
self%nj=1
self%nrslice=self%nr
self%ncslice=self%nc
ALLOCATE(self%i_map(1),self%j_map(1))
self%i_map(1)%n=self%nr
self%j_map(1)%n=self%nc
CALL native_matrix_set_pattern(self,self%ext%pattern)
self%lcmap=>self%ext%lcmap
ALLOCATE(self%M(self%nnz))
CALL self%update_slice
DEBUG_STACK_POP
end subroutine submatrix_setup
!---------------------------------------------------------------------------------
!> Update matrix with values from parent matrix
!---------------------------------------------------------------------------------
subroutine submatrix_update_slice(self)
class(oft_native_submatrix), intent(inout) :: self
DEBUG_STACK_PUSH
IF(oft_debug_print(3))WRITE(*,*)'  Get local slice'
IF(.NOT.ASSOCIATED(self%ext))CALL oft_abort('Submatrix not set up','submatrix_update_slice',__FILE__)
CALL matrix_update_slice(self%parent)
CALL scatter_local(self%ext%vplan,self%parent%Mfull,self%M,OFT_SCATTER_INSERT)
self%full_current=.FALSE.
DEBUG_STACK_POP
end subroutine submatrix_update_slice
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine submatrix_delete(self)
class(oft_native_submatrix), intent(inout) :: self
DEBUG_STACK_PUSH
CALL mat_delete(self)
NULLIFY(self%lcmap,self%ext)
IF(ASSOCIATED(self%part))NULLIFY(self%part)
CALL native_pattern_release(self%parent_pattern)
NULLIFY(self%parent)
DEBUG_STACK_POP
end subroutine submatrix_delete
!---------------------------------------------------------------------------------
!> Get the extraction of a submatrix (rows by columns) from a pattern, building it if needed
!!
!! Extractions are cached in a list held by `pattern` and matched by their row and
!! column lists, so the symbolic work is done once per pattern and set of rows and
!! columns. Submatrix rows and columns are numbered in the order of `rows` and `cols`
!! (parent local indices, each list without repeats); entries of each submatrix row
!! are sorted by column. The build is threaded over rows and may run concurrently
!! for different extractions; cache access is serialized.
!---------------------------------------------------------------------------------
FUNCTION native_extraction_get(pattern,nr,rows,nc,cols) RESULT(ext)
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Parent pattern
INTEGER(i4), INTENT(in) :: nr !< Number of submatrix rows
INTEGER(i4), INTENT(in) :: rows(nr) !< Parent row of each submatrix row [nr]
INTEGER(i4), INTENT(in) :: nc !< Number of submatrix columns
INTEGER(i4), INTENT(in) :: cols(nc) !< Parent column of each submatrix column [nc]
TYPE(native_extraction), POINTER :: ext !< Extraction
TYPE(native_extraction), POINTER :: new
!---Look for a cached extraction
!$omp critical (native_extraction_cache)
ext=>native_extraction_find(pattern,nr,rows,nc,cols)
!$omp end critical (native_extraction_cache)
IF(ASSOCIATED(ext))RETURN
!---Build outside the lock and insert, unless another thread inserted the same extraction
CALL native_extraction_build(pattern,nr,rows,nc,cols,new)
!$omp critical (native_extraction_cache)
ext=>native_extraction_find(pattern,nr,rows,nc,cols)
IF(.NOT.ASSOCIATED(ext))THEN
  new%next=>pattern%extracts
  pattern%extracts=>new
  ext=>new
  NULLIFY(new)
END IF
!$omp end critical (native_extraction_cache)
IF(ASSOCIATED(new))CALL native_extraction_free(new)
END FUNCTION native_extraction_get
!---------------------------------------------------------------------------------
!> Find a cached extraction with the given rows and columns (NULL if none)
!---------------------------------------------------------------------------------
FUNCTION native_extraction_find(pattern,nr,rows,nc,cols) RESULT(ext)
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Parent pattern
INTEGER(i4), INTENT(in) :: nr !< Number of submatrix rows
INTEGER(i4), INTENT(in) :: rows(nr) !< Parent row of each submatrix row [nr]
INTEGER(i4), INTENT(in) :: nc !< Number of submatrix columns
INTEGER(i4), INTENT(in) :: cols(nc) !< Parent column of each submatrix column [nc]
TYPE(native_extraction), POINTER :: ext !< Matching extraction (NULL if none)
ext=>pattern%extracts
DO WHILE(ASSOCIATED(ext))
  IF((ext%nr==nr).AND.(ext%nc==nc))THEN
    IF(ALL(ext%rows==rows).AND.ALL(ext%cols==cols))RETURN
  END IF
  ext=>ext%next
END DO
END FUNCTION native_extraction_find
!---------------------------------------------------------------------------------
!> Build the symbolic extraction of a submatrix (threaded over submatrix rows)
!---------------------------------------------------------------------------------
SUBROUTINE native_extraction_build(pattern,nr,rows,nc,cols,ext)
TYPE(oft_native_pattern), POINTER, INTENT(in) :: pattern !< Parent pattern
INTEGER(i4), INTENT(in) :: nr !< Number of submatrix rows
INTEGER(i4), INTENT(in) :: rows(nr) !< Parent row of each submatrix row [nr]
INTEGER(i4), INTENT(in) :: nc !< Number of submatrix columns
INTEGER(i4), INTENT(in) :: cols(nc) !< Parent column of each submatrix column [nc]
TYPE(native_extraction), POINTER, INTENT(out) :: ext !< New extraction (not cached)
INTEGER(i4) :: i,j,k,r,nnz
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cmap
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr,lc
ALLOCATE(ext)
ext%nr=nr
ext%nc=nc
ALLOCATE(ext%rows(nr),ext%cols(nc))
ext%rows=rows
ext%cols=cols
!---Submatrix column of each parent column
ALLOCATE(cmap(pattern%nc))
cmap=0
DO k=1,nc
  IF(cmap(cols(k))/=0)CALL oft_abort('Repeated column','native_extraction_build',__FILE__)
  cmap(cols(k))=k
END DO
!---Count entries of each row
ALLOCATE(kr(nr+1))
kr(1)=1
!$omp parallel do private(i,j,k) if(nr>OFT_OMP_VTHRESH/10)
DO r=1,nr
  i=rows(r)
  k=0
  DO j=pattern%kr(i),pattern%kr(i+1)-1
    IF(cmap(pattern%lc(j))/=0)k=k+1
  END DO
  kr(r+1)=k
END DO
DO r=1,nr
  kr(r+1)=kr(r+1)+kr(r)
END DO
nnz=kr(nr+1)-1
!---Fill and sort each row by submatrix column
ALLOCATE(lc(nnz),ext%lcmap(nnz))
!$omp parallel do private(i,j,k) schedule(guided) if(nr>OFT_OMP_VTHRESH/10)
DO r=1,nr
  i=rows(r)
  k=kr(r)
  DO j=pattern%kr(i),pattern%kr(i+1)-1
    IF(cmap(pattern%lc(j))/=0)THEN
      lc(k)=cmap(pattern%lc(j))
      ext%lcmap(k)=j
      k=k+1
    END IF
  END DO
  IF(kr(r+1)-kr(r)>1)CALL sort_array(lc(kr(r):kr(r+1)-1),ext%lcmap(kr(r):kr(r+1)-1),kr(r+1)-kr(r))
END DO
DEALLOCATE(cmap)
CALL native_pattern_create(ext%pattern,nr,nc,kr,lc)
CALL native_pattern_retain(ext%pattern)
CALL scatter_plan_setup_local(ext%vplan,nnz,ext%lcmap)
END SUBROUTINE native_extraction_build
!---------------------------------------------------------------------------------
!> Free an extraction (not the rest of its list)
!---------------------------------------------------------------------------------
RECURSIVE SUBROUTINE native_extraction_free(ext)
TYPE(native_extraction), POINTER, INTENT(inout) :: ext !< Extraction
IF(.NOT.ASSOCIATED(ext))RETURN
IF(ASSOCIATED(ext%rows))DEALLOCATE(ext%rows)
IF(ASSOCIATED(ext%cols))DEALLOCATE(ext%cols)
IF(ASSOCIATED(ext%lcmap))DEALLOCATE(ext%lcmap)
CALL scatter_plan_delete(ext%vplan)
CALL native_pattern_release(ext%pattern)
DEALLOCATE(ext)
END SUBROUTINE native_extraction_free
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
!> Compute matrix-vector product with a real vector
!!
!! b = self * a
!---------------------------------------------------------------------------------
SUBROUTINE cmat_apply_vec(self,a,b)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Source vector
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: btmp
DEBUG_STACK_PUSH
IF(b%n/=self%nr)CALL oft_abort('Row mismatch','cmat_apply_vec',__FILE__)
IF(a%n/=self%nc)CALL oft_abort('Col mismatch','cmat_apply_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
IF(ASSOCIATED(self%pattern%lc_small))THEN
  CALL nk_csr_matvec(self%nr,self%kr,self%pattern%lc_small,self%M,atmp,btmp)
ELSE
  CALL nk_csr_matvec(self%nr,self%kr,self%lc,self%M,atmp,btmp)
END IF
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)CALL b%stitch(1)
DEBUG_STACK_POP
END SUBROUTINE cmat_apply_vec
!---------------------------------------------------------------------------------
!> Compute matrix-vector product with a complex vector
!!
!! b = self * a
!---------------------------------------------------------------------------------
SUBROUTINE cmat_apply_cvec(self,a,b)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_cvector), TARGET, INTENT(inout) :: a !< Source vector
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(b%n/=self%nr)CALL oft_abort('Row mismatch','cmat_apply_cvec',__FILE__)
IF(a%n/=self%nc)CALL oft_abort('Col mismatch','cmat_apply_cvec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
IF(ASSOCIATED(self%pattern%lc_small))THEN
  CALL nk_csr_matvec(self%nr,self%kr,self%pattern%lc_small,self%M,atmp,btmp)
ELSE
  CALL nk_csr_matvec(self%nr,self%kr,self%lc,self%M,atmp,btmp)
END IF
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)CALL b%stitch(1)
DEBUG_STACK_POP
END SUBROUTINE cmat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose with a real vector
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
SUBROUTINE cmat_applyt_vec(self,a,b)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Source vector
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: btmp
DEBUG_STACK_PUSH
IF(b%n/=self%nc)CALL oft_abort('Col mismatch','cmat_applyt_vec',__FILE__)
IF(a%n/=self%nr)CALL oft_abort('Row mismatch','cmat_applyt_vec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!$omp critical (oft_native_csc_build)
CALL native_csc_build(self%pattern%csc,self%kr,self%lc,self%nr,self%nc)
!$omp end critical (oft_native_csc_build)
CALL nk_csc_matvec(self%nc,self%pattern%csc%kc,self%pattern%csc%lr,self%pattern%csc%lm,self%M,atmp,btmp)
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)CALL b%stitch(1)
DEBUG_STACK_POP
END SUBROUTINE cmat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose with a complex vector
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
SUBROUTINE cmat_applyt_cvec(self,a,b)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_cvector), TARGET, INTENT(inout) :: a !< Source vector
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
LOGICAL :: copied(2)
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: atmp,btmp
DEBUG_STACK_PUSH
IF(b%n/=self%nc)CALL oft_abort('Col mismatch','cmat_applyt_cvec',__FILE__)
IF(a%n/=self%nr)CALL oft_abort('Row mismatch','cmat_applyt_cvec',__FILE__)
copied(1)=native_get_vals(a,atmp)
copied(2)=native_get_vals(b,btmp)
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!$omp critical (oft_native_csc_build)
CALL native_csc_build(self%pattern%csc,self%kr,self%lc,self%nr,self%nc)
!$omp end critical (oft_native_csc_build)
CALL nk_csc_matvec(self%nc,self%pattern%csc%kc,self%pattern%csc%lr,self%pattern%csc%lm,self%M,atmp,btmp)
IF(copied(1))DEALLOCATE(atmp)
CALL native_put_vals(b,btmp,copied(2))
IF(.NOT.self%force_local)CALL b%stitch(1)
DEBUG_STACK_POP
END SUBROUTINE cmat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ANY(i_inds>self%nr))call oft_abort('Index exceeds graph','cmat_set_values',__FILE__)
if(ANY(j_inds>self%nc))call oft_abort('Index exceeds graph','cmat_set_values',__FILE__)
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','cmat_set_values',__FILE__)
CALL nk_set_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
end subroutine cmat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
SUBROUTINE cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries to add [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries to add [m]
COMPLEX(c8), INTENT(in) :: b(n,m) !< Values to set [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of entry locations
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','cmat_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','cmat_add_values',__FILE__)
self%full_current=.FALSE.
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock,loc_cache)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','cmat_add_values',__FILE__)
CALL nk_add_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
END SUBROUTINE cmat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix atomically (for concurrent assembly)
!---------------------------------------------------------------------------------
subroutine cmat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
COMPLEX(c8), INTENT(in) :: b(n,m) !< Values to set [n,m]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
INTEGER(i4) :: loc(n,m)
DEBUG_STACK_PUSH
IF(ANY(i_inds>self%nr))CALL oft_abort('Index exceeds graph','cmat_atomic_add_values',__FILE__)
IF(ANY(j_inds>self%nc))CALL oft_abort('Index exceeds graph','cmat_atomic_add_values',__FILE__)
self%full_current=.FALSE.
CALL native_locate(self%kr,self%lc,self%map,self%j_map,n,m,i_inds,j_inds,loc,iblock,jblock,loc_cache)
IF(ANY(loc==0))CALL oft_abort('Entry not found!','cmat_atomic_add_values',__FILE__)
CALL nk_atomic_add_at(n,m,loc,b,self%M)
DEBUG_STACK_POP
end subroutine cmat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine cmat_assemble(self,diag)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%pattern))CALL native_matrix_adopt_pattern(self)
!---Setup diagonal scaling
if(present(diag))then
  IF(associated(self%D))THEN
    call self%D%delete
    DEALLOCATE(self%D)
  END IF
  call diag%new(self%D)
  SELECT TYPE(this=>self%D)
  CLASS IS(oft_native_cvector)
    !---Get diagonal matrix values (rows without a diagonal entry get zero)
    !$omp critical (oft_native_diag_build)
    CALL native_diag_build(self%pattern%diag_index,self%kr,self%lc,self%nr)
    !$omp end critical (oft_native_diag_build)
    CALL nk_csr_get_diag(self%nr,self%pattern%diag_index%loc,self%M,this%v)
  CLASS DEFAULT
    CALL oft_abort('"diag" is not a vector object.','cmat_assemble',__FILE__)
  END SELECT
  call self%D%stitch(1)
  call diag%add((0.d0,0.d0),(1.d0,0.d0),self%D)
end if
!---Common assembly tasks
!$omp critical (oft_native_lc_small_build)
CALL native_lc_small_build(self%pattern%lc_small,self%lc,self%nnz,self%nc)
!$omp end critical (oft_native_lc_small_build)
DEBUG_STACK_POP
end subroutine cmat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine cmat_zero(self)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
DEBUG_STACK_PUSH
self%full_current=.FALSE.
CALL nk_set(self%nnz,(0.d0,0.d0),self%M)
DEBUG_STACK_POP
end subroutine cmat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine cmat_zero_rows(self,nrows,irows,iblock,keep_diag)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: ioff
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
ioff=0
IF(PRESENT(iblock))ioff=self%i_map(iblock)%offset
CALL nk_csr_zero_rows(nrows,irows,ioff,self%kr,self%lc,zero_diag,self%M)
DEBUG_STACK_POP
end subroutine cmat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine cmat_zero_cols(self,ncols,icols,jblock,keep_diag)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: coff
DEBUG_STACK_PUSH
self%full_current=.FALSE.
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
coff=0
IF(PRESENT(jblock))coff=self%j_map(jblock)%offset
CALL nk_csr_zero_cols(self%nr,ncols,icols,coff,self%kr,self%lc,zero_diag,self%M)
DEBUG_STACK_POP
end subroutine cmat_zero_cols
!---------------------------------------------------------------------------------
!> Delete matrix
!!
!! Values are freed and the reference to the non-zero pattern is released, so
!! the pattern is only freed if no other matrix (or holder) uses it.
!---------------------------------------------------------------------------------
subroutine cmat_delete(self)
CLASS(oft_native_cmatrix), INTENT(inout) :: self !< Matrix object
DEBUG_STACK_PUSH
self%full_current=.FALSE.
if(ASSOCIATED(self%M))deallocate(self%M)
if(ASSOCIATED(self%Md))deallocate(self%Md)
IF(ASSOCIATED(self%Mfull))DEALLOCATE(self%Mfull)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  DEALLOCATE(self%D)
END IF
!---Block mappings belong to the matrix, block maps to the pattern
IF(ASSOCIATED(self%map))DEALLOCATE(self%i_map,self%j_map)
!---Matrices built without a pattern own their structure
IF((.NOT.ASSOCIATED(self%pattern)).AND.ASSOCIATED(self%kr))CALL native_matrix_adopt_pattern(self)
CALL native_pattern_release(self%pattern)
CALL native_pattern_alias(self%pattern,self%kr,self%lc,self%map,self%linkage,self%redind,self%nred,self%nnz)
IF(ASSOCIATED(self%color))DEALLOCATE(self%color)
DEBUG_STACK_POP
end subroutine cmat_delete
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_apply_vec(self,a,b)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvec(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_mat_apply_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_mat_apply_cvec(self,a,b)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvec(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_mat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply the transposed matrix to a field.
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine dense_mat_applyt_vec(self,a,b)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvect(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_mat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply the transposed matrix to a field.
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine dense_mat_applyt_cvec(self,a,b)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvect(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_mat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
SUBROUTINE dense_mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
REAL(r8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_mat_set_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_mat_set_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_mat_set_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,0)
DEBUG_STACK_POP
END SUBROUTINE dense_mat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
SUBROUTINE dense_mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
REAL(r8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of entry locations (unused)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_mat_add_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_mat_add_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_mat_add_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,1)
DEBUG_STACK_POP
END SUBROUTINE dense_mat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix atomically (for concurrent assembly)
!---------------------------------------------------------------------------------
SUBROUTINE dense_mat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
REAL(r8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of entry locations (unused)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_mat_atomic_add_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_mat_atomic_add_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_mat_atomic_add_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,2)
DEBUG_STACK_POP
END SUBROUTINE dense_mat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine dense_mat_assemble(self,diag)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
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
    CALL nk_dense_get_diag(self%nr,self%nc,self%M,this%v)
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
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: offset
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
offset=0
IF(PRESENT(iblock))offset=self%i_map(iblock)%offset
CALL nk_dense_zero_rows(self%nr,self%nc,self%M,nrows,irows,offset,zero_diag)
DEBUG_STACK_POP
end subroutine dense_mat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine dense_mat_zero_cols(self,ncols,icols,jblock,keep_diag)
CLASS(oft_native_dense_matrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: offset
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
offset=0
IF(PRESENT(jblock))offset=self%j_map(jblock)%offset
CALL nk_dense_zero_cols(self%nr,self%nc,self%M,ncols,icols,offset,zero_diag)
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
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvec(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_cmat_apply_vec
!---------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_apply_cvec(self,a,b)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvec(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_cmat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply the transposed matrix to a field.
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_applyt_vec(self,a,b)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
CLASS(oft_cvector), INTENT(inout) :: b !< Result of matrix product
real(r8), pointer, CONTIGUOUS, dimension(:) :: avals
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvect(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_cmat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply the transposed matrix to a field.
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine dense_cmat_applyt_cvec(self,a,b)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
COMPLEX(c8), pointer, CONTIGUOUS, dimension(:) :: avals,bvals
DEBUG_STACK_PUSH
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
CALL nk_dense_matvect(self%nr,self%nc,self%M,avals,bvals)
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine dense_cmat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
SUBROUTINE dense_cmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
COMPLEX(c8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_cmat_set_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_set_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,0)
DEBUG_STACK_POP
END SUBROUTINE dense_cmat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
SUBROUTINE dense_cmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
COMPLEX(c8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of entry locations (unused)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_cmat_add_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_add_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,1)
DEBUG_STACK_POP
END SUBROUTINE dense_cmat_add_values
!---------------------------------------------------------------------------------
!> Add values to a matrix atomically (for concurrent assembly)
!---------------------------------------------------------------------------------
SUBROUTINE dense_cmat_atomic_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
INTEGER(i4), INTENT(in) :: n !< Number of rows in local matrix
INTEGER(i4), INTENT(in) :: m !< Number of columns in local matrix
INTEGER(i4), INTENT(in) :: i_inds(n) !< Row indices of entries [n]
INTEGER(i4), INTENT(in) :: j_inds(m) !< Column indices of entries [m]
COMPLEX(c8), INTENT(in) :: b(n,m) !< Values [n,m]
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Row block (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: jblock !< Column block (optional)
INTEGER(i4), OPTIONAL, INTENT(inout) :: loc_cache(n,m) !< Cache of entry locations (unused)
INTEGER(i4) :: ioff,joff
DEBUG_STACK_PUSH
ioff=0
joff=0
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  ioff=self%i_map(iblock)%offset
  joff=self%j_map(jblock)%offset
ELSE IF(PRESENT(iblock).OR.PRESENT(jblock))THEN
  CALL oft_abort('Only one block index was supplied','dense_cmat_atomic_add_values',__FILE__)
END IF
IF(ANY(ioff+i_inds>self%nr))CALL oft_abort('Index exceeds graph','dense_cmat_atomic_add_values',__FILE__)
IF(ANY(joff+j_inds>self%nc))CALL oft_abort('Index exceeds graph','dense_cmat_atomic_add_values',__FILE__)
CALL nk_dense_insert(self%nr,self%nc,self%M,n,m,i_inds,j_inds,ioff,joff,b,2)
DEBUG_STACK_POP
END SUBROUTINE dense_cmat_atomic_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine dense_cmat_assemble(self,diag)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
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
    CALL nk_dense_get_diag(self%nr,self%nc,self%M,this%v)
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
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: offset
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
offset=0
IF(PRESENT(iblock))offset=self%i_map(iblock)%offset
CALL nk_dense_zero_rows(self%nr,self%nc,self%M,nrows,irows,offset,zero_diag)
DEBUG_STACK_POP
end subroutine dense_cmat_zero_rows
!---------------------------------------------------------------------------------
!> Zero all entries in the specified columns
!---------------------------------------------------------------------------------
subroutine dense_cmat_zero_cols(self,ncols,icols,jblock,keep_diag)
CLASS(oft_native_dense_cmatrix), INTENT(inout) :: self !< Matrix object
integer(i4), intent(in) :: ncols !< Number of columns to zero
integer(i4), intent(in) :: icols(ncols) !< Indices of columns to zero [ncols]
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
logical :: zero_diag
INTEGER(i4) :: offset
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=(.NOT.keep_diag)
offset=0
IF(PRESENT(jblock))offset=self%j_map(jblock)%offset
CALL nk_dense_zero_cols(self%nr,self%nc,self%M,ncols,icols,offset,zero_diag)
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
