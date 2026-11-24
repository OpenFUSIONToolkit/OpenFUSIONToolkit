!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Module body template for native product kernels (left operand a or matrix M, right operand x,
! result y), included by native_kernels.F90
!
! Required macros (undefined by the including file after each instantiation):
!   KA     - Type of left operand or matrix values (eg. REAL(r8))
!   KB     - Type of right operand (eg. COMPLEX(c8))
!   KR     - Result type (COMPLEX(c8) if either KA or KB is complex)
!   KRZERO - Zero of type KR
!   KBPTR  - Pointer-to-values type for KB (native_rvals_ptr or native_cvals_ptr)
!   KCONJ  - Function-like macro conjugating a KA value in dot products (identity if real)
!   KREAL  - Defined only for the all-real instantiation
!
! The including file provides the MODULE and END MODULE statements (so the module
! name is visible to CMake's dependency scanner). Specific procedures are private;
! only the generic names are public, so generics of the same name from each
! instantiation merge in oft_native_kernels.
!---------------------------------------------------------------------------------
USE, INTRINSIC :: iso_c_binding, ONLY: c_int8_t
USE omp_lib
USE oft_local
USE oft_base, ONLY: oft_thread_slice
USE oft_native_kernels_base
IMPLICIT NONE
PRIVATE
PUBLIC :: nk_owned_dot, nk_owned_mdot, nk_csr_matvec, nk_csc_matvec
PUBLIC :: nk_dense_matvec, nk_dense_matvect
INTERFACE nk_owned_dot
  MODULE PROCEDURE k_owned_dot
END INTERFACE nk_owned_dot
INTERFACE nk_owned_mdot
  MODULE PROCEDURE k_owned_mdot
END INTERFACE nk_owned_mdot
INTERFACE nk_csr_matvec
  MODULE PROCEDURE k_csr_matvec
  MODULE PROCEDURE k_csr_matvec_small
END INTERFACE nk_csr_matvec
INTERFACE nk_csc_matvec
  MODULE PROCEDURE k_csc_matvec
END INTERFACE nk_csc_matvec
INTERFACE nk_dense_matvec
  MODULE PROCEDURE k_dense_matvec
END INTERFACE nk_dense_matvec
INTERFACE nk_dense_matvect
  MODULE PROCEDURE k_dense_matvect
END INTERFACE nk_dense_matvect
CONTAINS
!---------------------------------------------------------------------------------
!> Local dot product over owned entries, \f$ \sum_{i \in owned} \overline{a_i} x_i \f$
!!
!! The left operand is conjugated when it is complex. With `nw=n`, entry `i` is
!! weighted by the ownership weight `w(i)` (1 if owned by this process, 0 if owned by
!! another), so values shared between processes are counted once overall; with
!! `nw=0` all entries are owned. Only owned entries are summed, rather than summing
!! all entries and subtracting the shared ones, which avoids catastrophic
!! cancellation when shared entries dominate. No reduction across processes.
!---------------------------------------------------------------------------------
FUNCTION k_owned_dot(n,a,x,nw,w) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Number of entries
KA, INTENT(in) :: a(n) !< Left operand [n]
KB, INTENT(in) :: x(n) !< Right operand [n]
INTEGER(i4), INTENT(in) :: nw !< Number of weights (`n`, or 0 if all entries are owned)
INTEGER(c_int8_t), INTENT(in) :: w(nw) !< Ownership weight (1/0) of each entry [nw]
KR :: c !< Local dot product
INTEGER(i4) :: i
c=KRZERO
IF(nw>0)THEN
  !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    c=c+KCONJ(a(i))*x(i)*REAL(w(i),r8)
  END DO
ELSE
  !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    c=c+KCONJ(a(i))*x(i)
  END DO
END IF
END FUNCTION k_owned_dot
!---------------------------------------------------------------------------------
!> Local dot products with several arrays over owned entries
!!
!! \f$ c_k = \sum_{i \in owned} \overline{a_i} x_{k,i} \f$, computed in a single pass
!! over row blocks of @ref oft_native_kernels_base::native_vec_chunk "native_vec_chunk"
!! entries, so each block of `a` is loaded once for all `m` products. Each thread
!! processes one contiguous range of rows. See @ref k_owned_dot for the ownership
!! weights.
!---------------------------------------------------------------------------------
SUBROUTINE k_owned_mdot(n,m,a,x,nw,w,c)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: m !< Number of right operands
KA, INTENT(in) :: a(n) !< Left operand [n]
TYPE(KBPTR), INTENT(in) :: x(m) !< Right operands, each [n] [m]
INTEGER(i4), INTENT(in) :: nw !< Number of weights (`n`, or 0 if all entries are owned)
INTEGER(c_int8_t), INTENT(in) :: w(nw) !< Ownership weight (1/0) of each entry [nw]
KR, INTENT(out) :: c(m) !< Local dot products [m]
INTEGER(i4) :: i,k,ib,ie,i1,i2
KB, POINTER, CONTIGUOUS, DIMENSION(:) :: xv
KR :: tmp
c=KRZERO
IF(m<1)RETURN
!$omp parallel private(i1,i2,ib,ie,k,i,tmp,xv) reduction(+:c) if(n>OFT_OMP_VTHRESH)
CALL oft_thread_slice(omp_get_thread_num(),omp_get_num_threads(),n,i1,i2)
DO ib=i1,i2,native_vec_chunk
  ie=MIN(i2,ib+native_vec_chunk-1)
  DO k=1,m
    xv=>x(k)%v
    tmp=KRZERO
    IF(nw>0)THEN
      !$omp simd reduction(+:tmp)
      DO i=ib,ie
        tmp=tmp+KCONJ(a(i))*xv(i)*REAL(w(i),r8)
      END DO
    ELSE
      !$omp simd reduction(+:tmp)
      DO i=ib,ie
        tmp=tmp+KCONJ(a(i))*xv(i)
      END DO
    END IF
    c(k)=c(k)+tmp
  END DO
END DO
!$omp end parallel
END SUBROUTINE k_owned_mdot
!---------------------------------------------------------------------------------
!> CSR matrix-vector product, \f$ y = M x \f$
!!
!! Every entry of `y` is overwritten. Rows are split into one contiguous range per
!! thread with similar numbers of non-zeros (@ref oft_native_kernels_base::native_thread_rows
!! "native_thread_rows").
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_matvec(nr,kr,lc,vals,x,y)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(in) :: lc(:) !< CSR column indices [nnz]
KA, INTENT(in) :: vals(:) !< CSR values [nnz]
KB, INTENT(in) :: x(:) !< Source vector [nc]
KR, INTENT(out) :: y(nr) !< Result [nr]
INTEGER(i4) :: i,j,ifirst,ilast
KR :: tmp
!$omp parallel private(i,j,tmp,ifirst,ilast) if(kr(nr+1)-1>OFT_OMP_VTHRESH)
CALL native_thread_rows(kr,nr,ifirst,ilast)
DO i=ifirst,ilast
  tmp=KRZERO
  DO j=kr(i),kr(i+1)-1
    tmp=tmp+vals(j)*x(lc(j))
  END DO
  y(i)=tmp
END DO
!$omp end parallel
END SUBROUTINE k_csr_matvec
!---------------------------------------------------------------------------------
!> CSR matrix-vector product with 2-byte column indices, \f$ y = M x \f$
!!
!! Column `j` is stored as `j-lc_offset` (see @ref oft_native_kernels_base::lc_offset
!! "lc_offset"), which halves the index traffic for matrices with fewer than
!! `2*lc_offset` columns. Otherwise identical to the 4-byte version.
!---------------------------------------------------------------------------------
SUBROUTINE k_csr_matvec_small(nr,kr,lc,vals,x,y)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(2), INTENT(in) :: lc(:) !< Offset 2-byte CSR column indices [nnz]
KA, INTENT(in) :: vals(:) !< CSR values [nnz]
KB, INTENT(in) :: x(:) !< Source vector [nc]
KR, INTENT(out) :: y(nr) !< Result [nr]
INTEGER(i4) :: i,j,ifirst,ilast
KR :: tmp
!$omp parallel private(i,j,tmp,ifirst,ilast) if(kr(nr+1)-1>OFT_OMP_VTHRESH)
CALL native_thread_rows(kr,nr,ifirst,ilast)
DO i=ifirst,ilast
  tmp=KRZERO
  DO j=kr(i),kr(i+1)-1
    tmp=tmp+vals(j)*x(INT(lc(j),4)+lc_offset)
  END DO
  y(i)=tmp
END DO
!$omp end parallel
END SUBROUTINE k_csr_matvec_small
!---------------------------------------------------------------------------------
!> Transposed CSR matrix-vector product using a CSC index, \f$ y = M^T x \f$
!!
!! Column `k` of the CSR matrix holds entries `vals(lm(j))` in rows `lr(j)` for
!! `j=kc(k),...,kc(k+1)-1`. Every entry of `y` is overwritten, and columns are split
!! between threads as rows are in @ref k_csr_matvec.
!---------------------------------------------------------------------------------
SUBROUTINE k_csc_matvec(nc,kc,lr,lm,vals,x,y)
INTEGER(i4), INTENT(in) :: nc !< Number of columns
INTEGER(i4), INTENT(in) :: kc(nc+1) !< CSC column pointer [nc+1]
INTEGER(i4), INTENT(in) :: lr(:) !< Row index of each entry [nnz]
INTEGER(i4), INTENT(in) :: lm(:) !< Location of each entry in CSR values [nnz]
KA, INTENT(in) :: vals(:) !< CSR values [nnz]
KB, INTENT(in) :: x(:) !< Source vector [nr]
KR, INTENT(out) :: y(nc) !< Result [nc]
INTEGER(i4) :: j,k,ifirst,ilast
KR :: tmp
!$omp parallel private(j,k,tmp,ifirst,ilast) if(kc(nc+1)-1>OFT_OMP_VTHRESH)
CALL native_thread_rows(kc,nc,ifirst,ilast)
DO k=ifirst,ilast
  tmp=KRZERO
  DO j=kc(k),kc(k+1)-1
    tmp=tmp+vals(lm(j))*x(lr(j))
  END DO
  y(k)=tmp
END DO
!$omp end parallel
END SUBROUTINE k_csc_matvec
!---------------------------------------------------------------------------------
!> Dense (column-major) matrix-vector product, \f$ y = M x \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_matvec(nr,nc,vals,x,y)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: nc !< Number of columns
KA, INTENT(in) :: vals(nr,nc) !< Matrix values [nr,nc]
KB, INTENT(in) :: x(nc) !< Source vector [nc]
KR, INTENT(out) :: y(nr) !< Result [nr]
INTEGER(i4) :: i,j,i1,i2
#if defined(BLAS_THREADED) && defined(KREAL)
CALL dgemv('N',nr,nc,1.d0,vals,nr,x,1,0.d0,y,1)
#else
!$omp parallel private(i1,i2,i,j)
CALL oft_thread_slice(omp_get_thread_num(),omp_get_num_threads(),nr,i1,i2)
DO i=i1,i2
  y(i)=KRZERO
END DO
DO j=1,nc
  !$omp simd
  DO i=i1,i2
    y(i)=y(i)+vals(i,j)*x(j)
  END DO
END DO
!$omp end parallel
#endif
END SUBROUTINE k_dense_matvec
!---------------------------------------------------------------------------------
!> Transposed dense (column-major) matrix-vector product, \f$ y = M^T x \f$
!---------------------------------------------------------------------------------
SUBROUTINE k_dense_matvect(nr,nc,vals,x,y)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: nc !< Number of columns
KA, INTENT(in) :: vals(nr,nc) !< Matrix values [nr,nc]
KB, INTENT(in) :: x(nr) !< Source vector [nr]
KR, INTENT(out) :: y(nc) !< Result [nc]
INTEGER(i4) :: i,j,j1,j2
KR :: tmp
#if defined(BLAS_THREADED) && defined(KREAL)
CALL dgemv('T',nr,nc,1.d0,vals,nr,x,1,0.d0,y,1)
#else
!$omp parallel private(j1,j2,i,j,tmp)
CALL oft_thread_slice(omp_get_thread_num(),omp_get_num_threads(),nc,j1,j2)
DO j=j1,j2
  tmp=KRZERO
  !$omp simd reduction(+:tmp)
  DO i=1,nr
    tmp=tmp+vals(i,j)*x(i)
  END DO
  y(j)=tmp
END DO
!$omp end parallel
#endif
END SUBROUTINE k_dense_matvect
