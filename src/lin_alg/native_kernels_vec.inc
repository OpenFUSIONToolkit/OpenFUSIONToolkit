!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Module body template for single-type native kernels, included by native_kernels.F90
!
! Required macros (undefined by the including file after each instantiation):
!   KT     - Value type (eg. REAL(r8))
!   KZERO  - Zero of type KT
!   KONE   - One of type KT
!   KPTR   - Pointer-to-values type for KT (native_rvals_ptr or native_cvals_ptr)
!
! The including file provides the MODULE and END MODULE statements (so the module
! name is visible to CMake's dependency scanner). Specific procedures are private;
! only the generic names are public, so generics of the same name from each
! instantiation merge in oft_native_kernels.
!---------------------------------------------------------------------------------
USE omp_lib
USE oft_local
USE oft_base, ONLY: oft_thread_slice
USE oft_native_kernels_base
IMPLICIT NONE
PRIVATE
PUBLIC :: nk_set, nk_scale, nk_maxpy, nk_norm_inf, nk_owned_sum, nk_owned_norm1
PUBLIC :: nk_gather, nk_scatter, nk_csr_zero_rows, nk_csr_zero_cols, nk_csr_identity_rows
PUBLIC :: nk_csr_get_diag, nk_set_at, nk_add_at, nk_atomic_add_at
PUBLIC :: nk_dense_insert, nk_dense_zero_rows, nk_dense_zero_cols, nk_dense_get_diag
INTERFACE nk_set
  MODULE PROCEDURE k_set
END INTERFACE nk_set
INTERFACE nk_scale
  MODULE PROCEDURE k_scale
END INTERFACE nk_scale
INTERFACE nk_maxpy
  MODULE PROCEDURE k_maxpy
END INTERFACE nk_maxpy
INTERFACE nk_norm_inf
  MODULE PROCEDURE k_norm_inf
END INTERFACE nk_norm_inf
INTERFACE nk_owned_sum
  MODULE PROCEDURE k_owned_sum
END INTERFACE nk_owned_sum
INTERFACE nk_owned_norm1
  MODULE PROCEDURE k_owned_norm1
END INTERFACE nk_owned_norm1
INTERFACE nk_gather
  MODULE PROCEDURE k_gather
END INTERFACE nk_gather
INTERFACE nk_scatter
  MODULE PROCEDURE k_scatter
END INTERFACE nk_scatter
INTERFACE nk_csr_zero_rows
  MODULE PROCEDURE k_csr_zero_rows
END INTERFACE nk_csr_zero_rows
INTERFACE nk_csr_zero_cols
  MODULE PROCEDURE k_csr_zero_cols
END INTERFACE nk_csr_zero_cols
INTERFACE nk_csr_identity_rows
  MODULE PROCEDURE k_csr_identity_rows
END INTERFACE nk_csr_identity_rows
INTERFACE nk_csr_get_diag
  MODULE PROCEDURE k_csr_get_diag
END INTERFACE nk_csr_get_diag
INTERFACE nk_set_at
  MODULE PROCEDURE k_set_at
END INTERFACE nk_set_at
INTERFACE nk_add_at
  MODULE PROCEDURE k_add_at
END INTERFACE nk_add_at
INTERFACE nk_atomic_add_at
  MODULE PROCEDURE k_atomic_add_at
END INTERFACE nk_atomic_add_at
INTERFACE nk_dense_insert
  MODULE PROCEDURE k_dense_insert
END INTERFACE nk_dense_insert
INTERFACE nk_dense_zero_rows
  MODULE PROCEDURE k_dense_zero_rows
END INTERFACE nk_dense_zero_rows
INTERFACE nk_dense_zero_cols
  MODULE PROCEDURE k_dense_zero_cols
END INTERFACE nk_dense_zero_cols
INTERFACE nk_dense_get_diag
  MODULE PROCEDURE k_dense_get_diag
END INTERFACE nk_dense_get_diag
CONTAINS
!---------------------------------------------------------------------------------
!> Set all entries to a scalar, \f$ y_i = \alpha \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_set(n,alpha,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KT, INTENT(in) :: alpha !< Value
KT, INTENT(out) :: y(n) !< Array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=alpha
END DO
END SUBROUTINE k_set
!---------------------------------------------------------------------------------
!> Scale all entries, \f$ y_i = \alpha y_i \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_scale(n,alpha,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KT, INTENT(in) :: alpha !< Scale factor
KT, INTENT(inout) :: y(n) !< Array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=alpha*y(i)
END DO
END SUBROUTINE k_scale
!---------------------------------------------------------------------------------
!> Add a linear combination of arrays, \f$ y = y + \sum_k \alpha_k x_k \f$
!!
!! All arrays are processed in a single pass over row blocks of
!! @ref oft_native_kernels_base::native_vec_chunk "native_vec_chunk" entries, so each
!! block of `y` is loaded once for all `m` updates. Each thread processes one
!! contiguous range of rows. Terms are added in order of increasing `k`, giving the
!! same result as `m` successive `axpy` operations.
!---------------------------------------------------------------------------------
SUBROUTINE k_maxpy(n,m,alpha,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: m !< Number of arrays to add
KT, INTENT(in) :: alpha(m) !< Scale factors [m]
TYPE(KPTR), INTENT(in) :: x(m) !< Arrays to add, each [n] [m]
KT, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i,k,ib,ie,i1,i2
KT, POINTER, CONTIGUOUS, DIMENSION(:) :: xv
IF(m<1)RETURN
!$omp parallel private(i1,i2,ib,ie,k,i,xv) if(n>OFT_OMP_VTHRESH)
CALL oft_thread_slice(omp_get_thread_num(),omp_get_num_threads(),n,i1,i2)
DO ib=i1,i2,native_vec_chunk
  ie=MIN(i2,ib+native_vec_chunk-1)
  DO k=1,m
    xv=>x(k)%v
    !$omp simd
    DO i=ib,ie
      y(i)=y(i)+alpha(k)*xv(i)
    END DO
  END DO
END DO
!$omp end parallel
END SUBROUTINE k_maxpy
!---------------------------------------------------------------------------------
!> Local infinity norm, \f$ \max_i |x_i| \f$ (no reduction across processes)
!---------------------------------------------------------------------------------
FUNCTION k_norm_inf(n,x) RESULT(norm)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KT, INTENT(in) :: x(n) !< Array [n]
REAL(r8) :: norm !< Maximum absolute value (0 if `n==0`)
INTEGER(i4) :: i
norm=0.d0
!$omp parallel do simd reduction(max:norm) if(n>OFT_OMP_VTHRESH)
DO i=1,n
  norm=MAX(norm,ABS(x(i)))
END DO
END FUNCTION k_norm_inf
!---------------------------------------------------------------------------------
!> Local sum over owned entries, \f$ \sum_{i \in owned} x_i \f$
!!
!! Shared entries `lbe(i)` with `leo(i)==.FALSE.` are owned by another process and
!! are excluded (no reduction across processes).
!---------------------------------------------------------------------------------
FUNCTION k_owned_sum(n,x,nbe,lbe,leo) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KT, INTENT(in) :: x(n) !< Array [n]
INTEGER(i4), INTENT(in) :: nbe !< Number of shared entries
INTEGER(i4), INTENT(in) :: lbe(nbe) !< Location of shared entries [nbe]
LOGICAL, INTENT(in) :: leo(nbe) !< Ownership flag of shared entries [nbe]
KT :: c !< Local sum
INTEGER(i4) :: i
c=KZERO
!$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
!$omp do simd
DO i=1,n
  c=c+x(i)
END DO
!$omp do simd
DO i=1,nbe
  IF(.NOT.leo(i))c=c-x(lbe(i))
END DO
!$omp end parallel
END FUNCTION k_owned_sum
!---------------------------------------------------------------------------------
!> Local 1-norm over owned entries, \f$ \sum_{i \in owned} |x_i| \f$
!!
!! See @ref k_owned_sum for the treatment of shared entries.
!---------------------------------------------------------------------------------
FUNCTION k_owned_norm1(n,x,nbe,lbe,leo) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KT, INTENT(in) :: x(n) !< Array [n]
INTEGER(i4), INTENT(in) :: nbe !< Number of shared entries
INTEGER(i4), INTENT(in) :: lbe(nbe) !< Location of shared entries [nbe]
LOGICAL, INTENT(in) :: leo(nbe) !< Ownership flag of shared entries [nbe]
REAL(r8) :: c !< Local 1-norm
INTEGER(i4) :: i
c=0.d0
!$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
!$omp do simd
DO i=1,n
  c=c+ABS(x(i))
END DO
!$omp do simd
DO i=1,nbe
  IF(.NOT.leo(i))c=c-ABS(x(lbe(i)))
END DO
!$omp end parallel
END FUNCTION k_owned_norm1
!---------------------------------------------------------------------------------
!> Gather entries by index, \f$ y_i = x_{idx(i)} \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_gather(n,idx,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries to gather
INTEGER(i4), INTENT(in) :: idx(n) !< Source index of each entry [n]
KT, INTENT(in) :: x(:) !< Source array
KT, INTENT(inout) :: y(n) !< Destination array [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(i)=x(idx(i))
END DO
END SUBROUTINE k_gather
!---------------------------------------------------------------------------------
!> Scatter entries by index, \f$ y_{idx(i)} = x_i \f$
!!
!! Indices must be unique.
!---------------------------------------------------------------------------------
SUBROUTINE k_scatter(n,idx,x,y)
INTEGER(i4), INTENT(in) :: n !< Number of entries to scatter
INTEGER(i4), INTENT(in) :: idx(n) !< Destination index of each entry [n]
KT, INTENT(in) :: x(n) !< Source array [n]
KT, INTENT(inout) :: y(:) !< Destination array
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  y(idx(i))=x(i)
END DO
END SUBROUTINE k_scatter
!---------------------------------------------------------------------------------
!> Zero entries in a set of CSR rows
!!
!! Row `ioff+irows(i)` is zeroed for each `i`. Unless `zero_diag` is set, entries in
!! column `irows(i)` are kept.
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_zero_rows(nrows,irows,ioff,kr,lc,zero_diag,vals)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows to zero
INTEGER(i4), INTENT(in) :: irows(nrows) !< Row indices [nrows]
INTEGER(i4), INTENT(in) :: ioff !< Offset added to row indices
INTEGER(i4), INTENT(in) :: kr(:) !< CSR row pointer
INTEGER(i4), INTENT(in) :: lc(:) !< CSR column indices
LOGICAL, INTENT(in) :: zero_diag !< Zero entries in column `irows(i)` too?
KT, INTENT(inout) :: vals(:) !< CSR values
INTEGER(i4) :: i,j,ir
!$omp parallel do private(ir,j)
DO i=1,nrows
  ir=ioff+irows(i)
  DO j=kr(ir),kr(ir+1)-1
    IF(zero_diag.OR.(lc(j)/=irows(i)))vals(j)=KZERO
  END DO
END DO
END SUBROUTINE k_csr_zero_rows
!---------------------------------------------------------------------------------
!> Zero entries in a set of CSR columns
!!
!! Entries in columns `coff+icols(:)` are zeroed. Unless `zero_diag` is set,
!! diagonal entries are kept.
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_zero_cols(nr,ncols,icols,coff,kr,lc,zero_diag,vals)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: ncols !< Number of columns to zero
INTEGER(i4), INTENT(in) :: icols(ncols) !< Column indices [ncols]
INTEGER(i4), INTENT(in) :: coff !< Offset added to column indices
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(in) :: lc(:) !< CSR column indices
LOGICAL, INTENT(in) :: zero_diag !< Zero diagonal entries too?
KT, INTENT(inout) :: vals(:) !< CSR values
INTEGER(i4) :: i,j
!$omp parallel do private(j)
DO i=1,nr
  DO j=kr(i),kr(i+1)-1
    IF(ANY(lc(j)==icols+coff))THEN
      IF(zero_diag.OR.(lc(j)/=i))vals(j)=KZERO
    END IF
  END DO
END DO
END SUBROUTINE k_csr_zero_cols
!---------------------------------------------------------------------------------
!> Replace a set of CSR rows with rows of the identity matrix
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_identity_rows(nrows,irows,kr,lc,vals)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows
INTEGER(i4), INTENT(in) :: irows(nrows) !< Row indices [nrows]
INTEGER(i4), INTENT(in) :: kr(:) !< CSR row pointer
INTEGER(i4), INTENT(in) :: lc(:) !< CSR column indices
KT, INTENT(inout) :: vals(:) !< CSR values
INTEGER(i4) :: i,j,k
!$omp parallel do private(j,k)
DO i=1,nrows
  k=irows(i)
  DO j=kr(k),kr(k+1)-1
    vals(j)=KZERO
    IF(lc(j)==k)vals(j)=KONE
  END DO
END DO
END SUBROUTINE k_csr_identity_rows
!---------------------------------------------------------------------------------
!> Extract diagonal entries, \f$ d_i = M_{dloc(i)} \f$ (zero where `dloc(i)==0`)
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_get_diag(nr,dloc,vals,d)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: dloc(nr) !< Location of diagonal entry in each row (0 if absent) [nr]
KT, INTENT(in) :: vals(:) !< CSR values
KT, INTENT(out) :: d(nr) !< Diagonal entries [nr]
INTEGER(i4) :: i
!$omp parallel do if(nr>OFT_OMP_VTHRESH)
DO i=1,nr
  IF(dloc(i)>0)THEN
    d(i)=vals(dloc(i))
  ELSE
    d(i)=KZERO
  END IF
END DO
END SUBROUTINE k_csr_get_diag
!---------------------------------------------------------------------------------
!> Set entries at given locations, \f$ M_{loc(i,j)} = b_{i,j} \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_set_at(n,m,loc,b,vals)
INTEGER(i4), INTENT(in) :: n !< Number of rows in block
INTEGER(i4), INTENT(in) :: m !< Number of columns in block
INTEGER(i4), INTENT(in) :: loc(n,m) !< Locations [n,m]
KT, INTENT(in) :: b(n,m) !< Values [n,m]
KT, INTENT(inout) :: vals(:) !< Destination array
INTEGER(i4) :: i,j
DO j=1,m
  DO i=1,n
    vals(loc(i,j))=b(i,j)
  END DO
END DO
END SUBROUTINE k_set_at
!---------------------------------------------------------------------------------
!> Add to entries at given locations, \f$ M_{loc(i,j)} = M_{loc(i,j)} + b_{i,j} \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_add_at(n,m,loc,b,vals)
INTEGER(i4), INTENT(in) :: n !< Number of rows in block
INTEGER(i4), INTENT(in) :: m !< Number of columns in block
INTEGER(i4), INTENT(in) :: loc(n,m) !< Locations [n,m]
KT, INTENT(in) :: b(n,m) !< Values [n,m]
KT, INTENT(inout) :: vals(:) !< Destination array
INTEGER(i4) :: i,j
DO j=1,m
  DO i=1,n
    vals(loc(i,j))=vals(loc(i,j))+b(i,j)
  END DO
END DO
END SUBROUTINE k_add_at
!---------------------------------------------------------------------------------
!> Atomically add to entries at given locations (for concurrent assembly)
!---------------------------------------------------------------------------------
SUBROUTINE k_atomic_add_at(n,m,loc,b,vals)
INTEGER(i4), INTENT(in) :: n !< Number of rows in block
INTEGER(i4), INTENT(in) :: m !< Number of columns in block
INTEGER(i4), INTENT(in) :: loc(n,m) !< Locations [n,m]
KT, INTENT(in) :: b(n,m) !< Values [n,m]
KT, INTENT(inout) :: vals(:) !< Destination array
INTEGER(i4) :: i,j
DO j=1,m
  DO i=1,n
    !$omp atomic
    vals(loc(i,j))=vals(loc(i,j))+b(i,j)
  END DO
END DO
END SUBROUTINE k_atomic_add_at
!---------------------------------------------------------------------------------
!> Insert a block of values into a dense (column-major) matrix
!!
!! Entry (`ioff+irows(i)`,`joff+jcols(j)`) is set to (`mode=0`), incremented by
!! (`mode=1`) or atomically incremented by (`mode=2`) `b(i,j)`.
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_insert(nr,nc,vals,n,m,irows,jcols,ioff,joff,b,mode)
INTEGER(i4), INTENT(in) :: nr !< Number of matrix rows
INTEGER(i4), INTENT(in) :: nc !< Number of matrix columns
KT, INTENT(inout) :: vals(nr,nc) !< Matrix values [nr,nc]
INTEGER(i4), INTENT(in) :: n !< Number of rows in block
INTEGER(i4), INTENT(in) :: m !< Number of columns in block
INTEGER(i4), INTENT(in) :: irows(n) !< Row indices [n]
INTEGER(i4), INTENT(in) :: jcols(m) !< Column indices [m]
INTEGER(i4), INTENT(in) :: ioff !< Offset added to row indices
INTEGER(i4), INTENT(in) :: joff !< Offset added to column indices
KT, INTENT(in) :: b(n,m) !< Block values [n,m]
INTEGER(i4), INTENT(in) :: mode !< Set (0), add (1) or atomically add (2)
INTEGER(i4) :: i,j,ii,jj
DO j=1,m
  jj=joff+jcols(j)
  DO i=1,n
    ii=ioff+irows(i)
    SELECT CASE(mode)
    CASE(0)
      vals(ii,jj)=b(i,j)
    CASE(1)
      vals(ii,jj)=vals(ii,jj)+b(i,j)
    CASE DEFAULT
      !$omp atomic
      vals(ii,jj)=vals(ii,jj)+b(i,j)
    END SELECT
  END DO
END DO
END SUBROUTINE k_dense_insert
!---------------------------------------------------------------------------------
!> Zero rows `ioff+irows(:)` of a dense (column-major) matrix
!!
!! Unless `zero_diag` is set, diagonal entries are kept.
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_zero_rows(nr,nc,vals,nrows,irows,ioff,zero_diag)
INTEGER(i4), INTENT(in) :: nr !< Number of matrix rows
INTEGER(i4), INTENT(in) :: nc !< Number of matrix columns
KT, INTENT(inout) :: vals(nr,nc) !< Matrix values [nr,nc]
INTEGER(i4), INTENT(in) :: nrows !< Number of rows to zero
INTEGER(i4), INTENT(in) :: irows(nrows) !< Row indices [nrows]
INTEGER(i4), INTENT(in) :: ioff !< Offset added to row indices
LOGICAL, INTENT(in) :: zero_diag !< Zero diagonal entries too?
INTEGER(i4) :: i,j,ir
!$omp parallel do private(ir,j)
DO i=1,nrows
  ir=ioff+irows(i)
  DO j=1,nc
    IF(zero_diag.OR.(ir/=j))vals(ir,j)=KZERO
  END DO
END DO
END SUBROUTINE k_dense_zero_rows
!---------------------------------------------------------------------------------
!> Zero columns `joff+icols(:)` of a dense (column-major) matrix
!!
!! Unless `zero_diag` is set, diagonal entries are kept.
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_zero_cols(nr,nc,vals,ncols,icols,joff,zero_diag)
INTEGER(i4), INTENT(in) :: nr !< Number of matrix rows
INTEGER(i4), INTENT(in) :: nc !< Number of matrix columns
KT, INTENT(inout) :: vals(nr,nc) !< Matrix values [nr,nc]
INTEGER(i4), INTENT(in) :: ncols !< Number of columns to zero
INTEGER(i4), INTENT(in) :: icols(ncols) !< Column indices [ncols]
INTEGER(i4), INTENT(in) :: joff !< Offset added to column indices
LOGICAL, INTENT(in) :: zero_diag !< Zero diagonal entries too?
INTEGER(i4) :: i,j,ic
!$omp parallel do private(ic,j)
DO i=1,ncols
  ic=joff+icols(i)
  DO j=1,nr
    IF(zero_diag.OR.(ic/=j))vals(j,ic)=KZERO
  END DO
END DO
END SUBROUTINE k_dense_zero_cols
!---------------------------------------------------------------------------------
!> Extract the diagonal of a dense (column-major) matrix, \f$ d_i = M_{i,i} \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_get_diag(nr,nc,vals,d)
INTEGER(i4), INTENT(in) :: nr !< Number of matrix rows
INTEGER(i4), INTENT(in) :: nc !< Number of matrix columns
KT, INTENT(in) :: vals(nr,nc) !< Matrix values [nr,nc]
KT, INTENT(out) :: d(nr) !< Diagonal entries [nr]
INTEGER(i4) :: i
!$omp parallel do
DO i=1,nr
  d(i)=vals(i,i)
END DO
END SUBROUTINE k_dense_get_diag
