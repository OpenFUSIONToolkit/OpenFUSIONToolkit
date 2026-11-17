!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file native_kernels.F90
!
!> Plain-array compute kernels for the native linear algebra backend
!!
!! All inner loops over vector values and matrix entries used by the native vector,
!! matrix and solver classes live here, so each loop exists once and threading,
!! vectorization and (later) offload annotations are applied in one place. Kernels
!! take explicit sizes and contiguous arrays and know nothing about the vector or
!! matrix classes. They contain no MPI communication; reductions return the local
!! (owned) contribution only.
!!
!! Real, complex and mixed real/complex versions are generated from a single source
!! by including a template file once per type combination, each time inside its own
!! module with different type macros:
!! - `native_kernels_vec.inc`: one value type (`nk_set`, `nk_scale`, `nk_maxpy`,
!!   `nk_norm_inf`, `nk_owned_sum`, `nk_owned_norm1`, `nk_gather`, `nk_scatter`,
!!   CSR row/column/diagonal operations and value insertion) for real and complex.
!! - `native_kernels_upd.inc`: destination and source types (`nk_copy`, `nk_axpby`,
!!   `nk_axpbypcz`, `nk_mult`, `nk_div`) for real/real, complex/real and
!!   complex/complex.
!! - `native_kernels_prod.inc`: left (or matrix), right and result types
!!   (`nk_owned_dot`, `nk_owned_mdot`, `nk_csr_matvec`, `nk_csc_matvec`,
!!   `nk_dense_matvec`, `nk_dense_matvect`) for all four real/complex combinations.
!!
!! Specific procedures in each instantiation are private and share their names; only
!! generic names are public, and the generics of each instantiation merge in
!! @ref oft_native_kernels, which is the module to use. Kernels that only exist for
!! real values (Jacobi smoother sweeps) and the type-independent helpers are in
!! @ref oft_native_kernels_base.
!!
!! Threading follows the rest of the native backend: each kernel opens its own
!! OpenMP parallel region, conditional on the problem size exceeding
!! `OFT_OMP_VTHRESH`, so kernels must be called from serial code.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_native_la
!---------------------------------------------------------------------------------
#include "local.h"
!---------------------------------------------------------------------------------
!> Type-independent helpers and real-only kernels for the native backend
!---------------------------------------------------------------------------------
MODULE oft_native_kernels_base
USE omp_lib
USE oft_local
IMPLICIT NONE
PRIVATE
!> Offset used to store column indices in 2 bytes (`lc_small = lc - lc_offset`)
INTEGER(i4), PUBLIC, PARAMETER :: lc_offset = HUGE(INT(1,2))
!> Number of rows processed together by multi-vector kernels (`maxpy`, `mdot`), chosen
!! so that one row block of every vector involved stays in cache
INTEGER(i4), PUBLIC, PARAMETER :: native_vec_chunk = 512
!---------------------------------------------------------------------------------
!> Contiguous pointer to real values (for multi-vector kernels)
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_rvals_ptr
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: v => NULL() !< Values
END TYPE native_rvals_ptr
!---------------------------------------------------------------------------------
!> Contiguous pointer to complex values (for multi-vector kernels)
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: native_cvals_ptr
  COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: v => NULL() !< Values
END TYPE native_cvals_ptr
PUBLIC native_row_split, native_thread_rows, nk_jacobi_start, nk_jacobi_sweep, nk_jacobi_restart
CONTAINS
!---------------------------------------------------------------------------------
!> Get the first row of one part of a non-zero balanced, contiguous row partition
!!
!! Rows are split into `nparts` contiguous ranges that each contain approximately
!! the same number of non-zeros. Part `ipart` covers rows
!! `native_row_split(kr,nr,ipart,nparts)` to `native_row_split(kr,nr,ipart+1,nparts)-1`,
!! so the parts are disjoint and together cover all rows.
!---------------------------------------------------------------------------------
PURE FUNCTION native_row_split(kr,nr,ipart,nparts) RESULT(irow)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(in) :: ipart !< Part index (0-based, 0 <= ipart <= nparts)
INTEGER(i4), INTENT(in) :: nparts !< Number of parts
INTEGER(i4) :: irow !< First row in part (`nr+1` when `ipart>=nparts`)
INTEGER(i4) :: ilow,ihigh,imid
INTEGER(i8) :: nnz_target
IF(ipart<=0)THEN
  irow=1
  RETURN
ELSE IF(ipart>=nparts)THEN
  irow=nr+1
  RETURN
END IF
nnz_target=INT(kr(1),i8)+(INT(kr(nr+1)-kr(1),i8)*INT(ipart,i8))/INT(nparts,i8)
!---First row whose entries start at or after the target
ilow=1
ihigh=nr+1
DO WHILE(ilow<ihigh)
  imid=(ilow+ihigh)/2
  IF(INT(kr(imid),i8)<nnz_target)THEN
    ilow=imid+1
  ELSE
    ihigh=imid
  END IF
END DO
irow=ilow
END FUNCTION native_row_split
!---------------------------------------------------------------------------------
!> Get the range of rows to be processed by the calling OpenMP thread
!!
!! Uses @ref native_row_split with one part per thread in the current team, so each
!! thread processes a contiguous block of rows with a similar number of non-zeros.
!! Outside of an active parallel region all rows are returned. The same splitting
!! applies to columns when a CSC column pointer is passed instead.
!---------------------------------------------------------------------------------
SUBROUTINE native_thread_rows(kr,nr,ifirst,ilast)
INTEGER(i4), INTENT(in) :: nr !< Number of rows
INTEGER(i4), INTENT(in) :: kr(nr+1) !< CSR row pointer [nr+1]
INTEGER(i4), INTENT(out) :: ifirst !< First row for this thread
INTEGER(i4), INTENT(out) :: ilast !< Last row for this thread (`ilast<ifirst` if none)
INTEGER(i4) :: ithread,nthreads
ithread=omp_get_thread_num()
nthreads=omp_get_num_threads()
ifirst=native_row_split(kr,nr,ithread,nthreads)
ilast=native_row_split(kr,nr,ithread+1,nthreads)-1
END SUBROUTINE native_thread_rows
!---------------------------------------------------------------------------------
!> Start a damped Jacobi iteration from a zero guess
!!
!! \f$ u = 0 \f$, \f$ p = \omega D^{-1} g \f$
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_start(n,df,d,g,u,p)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: df !< Damping factor \f$ \omega \f$
REAL(r8), INTENT(in) :: d(n) !< Diagonal \f$ D \f$ [n]
REAL(r8), INTENT(in) :: g(n) !< Residual [n]
REAL(r8), INTENT(out) :: u(n) !< Solution [n]
REAL(r8), INTENT(out) :: p(n) !< Correction [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
DO i=1,n
  u(i)=0.d0
  p(i)=df*g(i)/d(i)
END DO
END SUBROUTINE nk_jacobi_start
!---------------------------------------------------------------------------------
!> One fused damped Jacobi sweep, given \f$ q = A p \f$ for the current correction
!!
!! \f$ u = u + p \f$ (only if `u` is present), \f$ g = g - q \f$,
!! \f$ p = \omega D^{-1} g \f$, and optionally saves the updated `u` and `g`.
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_sweep(n,df,d,q,g,p,u,usave,gsave)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: df !< Damping factor \f$ \omega \f$
REAL(r8), INTENT(in) :: d(n) !< Diagonal \f$ D \f$ [n]
REAL(r8), INTENT(in) :: q(n) !< Product of matrix with correction [n]
REAL(r8), INTENT(inout) :: g(n) !< Residual [n]
REAL(r8), INTENT(inout) :: p(n) !< Correction [n]
REAL(r8), OPTIONAL, INTENT(inout) :: u(n) !< Solution [n] (optional)
REAL(r8), OPTIONAL, INTENT(out) :: usave(n) !< Copy of updated solution [n] (optional, requires `u`)
REAL(r8), OPTIONAL, INTENT(out) :: gsave(n) !< Copy of updated residual [n] (optional, with `usave`)
INTEGER(i4) :: i
IF(PRESENT(u))THEN
  IF(PRESENT(usave))THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
    DO i=1,n
      u(i)=u(i)+p(i)
      g(i)=g(i)-q(i)
      p(i)=df*g(i)/d(i)
      usave(i)=u(i)
      gsave(i)=g(i)
    END DO
  ELSE
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
    DO i=1,n
      u(i)=u(i)+p(i)
      g(i)=g(i)-q(i)
      p(i)=df*g(i)/d(i)
    END DO
  END IF
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH/2)
  DO i=1,n
    g(i)=g(i)-q(i)
    p(i)=df*g(i)/d(i)
  END DO
END IF
END SUBROUTINE nk_jacobi_sweep
!---------------------------------------------------------------------------------
!> Restart a Jacobi iteration from a saved state
!!
!! \f$ p = u - u_{save} \f$, \f$ g = g_{save} \f$
!---------------------------------------------------------------------------------
SUBROUTINE nk_jacobi_restart(n,u,usave,gsave,g,p)
INTEGER(i4), INTENT(in) :: n !< Number of entries
REAL(r8), INTENT(in) :: u(n) !< Solution [n]
REAL(r8), INTENT(in) :: usave(n) !< Saved solution [n]
REAL(r8), INTENT(in) :: gsave(n) !< Saved residual [n]
REAL(r8), INTENT(out) :: g(n) !< Residual [n]
REAL(r8), INTENT(out) :: p(n) !< Correction [n]
INTEGER(i4) :: i
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  p(i)=u(i)-usave(i)
  g(i)=gsave(i)
END DO
END SUBROUTINE nk_jacobi_restart
END MODULE oft_native_kernels_base
!---------------------------------------------------------------------------------
! Single-type kernels
!---------------------------------------------------------------------------------
#define KT REAL(r8)
#define KZERO 0.d0
#define KONE 1.d0
#define KPTR native_rvals_ptr
MODULE oft_native_kernels_vec_r8
#include "native_kernels_vec.inc"
END MODULE oft_native_kernels_vec_r8
#undef KT
#undef KZERO
#undef KONE
#undef KPTR
#define KT COMPLEX(c8)
#define KZERO (0.d0,0.d0)
#define KONE (1.d0,0.d0)
#define KPTR native_cvals_ptr
MODULE oft_native_kernels_vec_c8
#include "native_kernels_vec.inc"
END MODULE oft_native_kernels_vec_c8
#undef KT
#undef KZERO
#undef KONE
#undef KPTR
!---------------------------------------------------------------------------------
! Update kernels (destination/source)
!---------------------------------------------------------------------------------
#define KY REAL(r8)
#define KX REAL(r8)
#define KYZERO 0.d0
#define KTOY(a) a
MODULE oft_native_kernels_upd_rr
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_rr
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
#define KY COMPLEX(c8)
#define KX REAL(r8)
#define KYZERO (0.d0,0.d0)
#define KTOY(a) CMPLX(a,0.d0,c8)
MODULE oft_native_kernels_upd_cr
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_cr
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
#define KY COMPLEX(c8)
#define KX COMPLEX(c8)
#define KYZERO (0.d0,0.d0)
#define KTOY(a) a
MODULE oft_native_kernels_upd_cc
#include "native_kernels_upd.inc"
END MODULE oft_native_kernels_upd_cc
#undef KY
#undef KX
#undef KYZERO
#undef KTOY
!---------------------------------------------------------------------------------
! Product kernels (left or matrix/right/result)
!---------------------------------------------------------------------------------
#define KA REAL(r8)
#define KB REAL(r8)
#define KR REAL(r8)
#define KRZERO 0.d0
#define KBPTR native_rvals_ptr
#define KCONJ(a) a
#define KREAL
MODULE oft_native_kernels_prod_rr
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_rr
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#undef KREAL
#define KA REAL(r8)
#define KB COMPLEX(c8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_cvals_ptr
#define KCONJ(a) a
MODULE oft_native_kernels_prod_rc
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_rc
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#define KA COMPLEX(c8)
#define KB REAL(r8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_rvals_ptr
#define KCONJ(a) CONJG(a)
MODULE oft_native_kernels_prod_cr
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_cr
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
#define KA COMPLEX(c8)
#define KB COMPLEX(c8)
#define KR COMPLEX(c8)
#define KRZERO (0.d0,0.d0)
#define KBPTR native_cvals_ptr
#define KCONJ(a) CONJG(a)
MODULE oft_native_kernels_prod_cc
#include "native_kernels_prod.inc"
END MODULE oft_native_kernels_prod_cc
#undef KA
#undef KB
#undef KR
#undef KRZERO
#undef KBPTR
#undef KCONJ
!---------------------------------------------------------------------------------
!> Native backend kernels for all supported value types
!!
!! Collects the generic kernel interfaces of every instantiation (see the file
!! description of native_kernels.F90) together with @ref oft_native_kernels_base.
!! Generic kernels:
!! - Single type: `nk_set`, `nk_scale`, `nk_maxpy`, `nk_norm_inf`, `nk_owned_sum`,
!!   `nk_owned_norm1`, `nk_gather`, `nk_scatter`, `nk_csr_zero_rows`,
!!   `nk_csr_zero_cols`, `nk_csr_identity_rows`, `nk_csr_get_diag`, `nk_set_at`,
!!   `nk_add_at`, `nk_atomic_add_at` (real and complex)
!! - Destination/source: `nk_copy`, `nk_axpby`, `nk_axpbypcz`, `nk_mult`, `nk_div`
!!   (real/real, complex/real, complex/complex)
!! - Products: `nk_owned_dot`, `nk_owned_mdot`, `nk_csr_matvec` (4- and 2-byte
!!   column indices), `nk_csc_matvec`, `nk_dense_matvec`, `nk_dense_matvect`
!!   (all real/complex combinations; the left operand of dot products is conjugated)
!! - Real only: `nk_jacobi_start`, `nk_jacobi_sweep`, `nk_jacobi_restart`
!---------------------------------------------------------------------------------
MODULE oft_native_kernels
USE oft_native_kernels_base
USE oft_native_kernels_vec_r8
USE oft_native_kernels_vec_c8
USE oft_native_kernels_upd_rr
USE oft_native_kernels_upd_cr
USE oft_native_kernels_upd_cc
USE oft_native_kernels_prod_rr
USE oft_native_kernels_prod_rc
USE oft_native_kernels_prod_cr
USE oft_native_kernels_prod_cc
IMPLICIT NONE
PUBLIC
END MODULE oft_native_kernels
