/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* UMFPACK (SuitSparse) LU solver interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_UMFPACK
#include "umfpack.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Integer interface: umfpack_di_* (32-bit, default) or umfpack_dl_* (64-bit).
   The 32-bit interface indexes UMFPACK's internal workspace with int, which
   limits the factors to roughly 16 GB (2^31 8-byte units); larger problems fail
   with UMFPACK_ERROR_out_of_memory. With OFT_UMFPACK_LONG the CSR graph passed
   from Fortran (always 32-bit) is copied to 64-bit arrays at each full
   factorization. */
#ifdef OFT_UMFPACK_LONG
typedef SuiteSparse_long umf_int;
#define UMF_FUN(name) umfpack_dl_##name
#else
typedef int umf_int;
#define UMF_FUN(name) umfpack_di_##name
#endif

typedef struct {
  void *Symbolic;
  void *Numeric;
  umf_int *Wi;
  double *W;
  umf_int *Ap; /* 64-bit copy of rowptr (OFT_UMFPACK_LONG only) */
  umf_int *Ai; /* 64-bit copy of colind (OFT_UMFPACK_LONG only) */
  double Control [UMFPACK_CONTROL];
} factors_t;

void
oft_umfpack_dgssv_c(int iopt, int n, int nnz, int nrhs,
                    double *values, int *colind, int *rowptr,
                    double *b, int ldb,
                    factors_t **f_factors,
                    int col_perm, bool iter_refine, int *info)
/*
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs full LU decomposition
 *      = 2, performs triangular solve
 *      = 3, performs numerical factorization only
 *      = 4, free all the storage in the end
 *
 * f_factors (input/output) factors_t**
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
{
  int i,j;
  int status = 0;
  double Info [UMFPACK_INFO];
  factors_t *LUfactors;
  umf_int *Ap, *Ai;
  bool refactor;

    if ( iopt == 1 || iopt == 3 ) { /* Perform full factorization */

      // Create or unpack storage structure
      if ( !*f_factors ) {
        if ( iopt == 3) {
          printf("UMFPACK: Refactor requested without symbolic factors\n");
          *info = -1;
          return;
        }
        LUfactors = (factors_t*) malloc(sizeof(factors_t));
        if ( !LUfactors ) {
          printf("UMFPACK: Unable to allocate factor structure\n");
          *info = UMFPACK_ERROR_out_of_memory;
          return;
        }
        LUfactors->Symbolic = NULL;
        LUfactors->Numeric = NULL;
        LUfactors->Ap = NULL;
        LUfactors->Ai = NULL;
        LUfactors->Wi = (umf_int *) malloc(n * sizeof(umf_int));
        /* Set the default input options. */
        UMF_FUN(defaults)( LUfactors->Control);
        if (iter_refine) {
          LUfactors->W  = (double *) malloc(5*n * sizeof(double));
        } else {
          LUfactors->Control[UMFPACK_IRSTEP] = 0;
          LUfactors->W  = (double *) malloc(n * sizeof(double));
        }
        if ( !LUfactors->Wi || !LUfactors->W ) {
          printf("UMFPACK: Unable to allocate work arrays\n");
          free(LUfactors->Wi);
          free(LUfactors->W);
          free(LUfactors);
          *info = UMFPACK_ERROR_out_of_memory;
          return;
        }
        *f_factors = LUfactors;
        refactor = false;
      } else {
        LUfactors = *f_factors;
        refactor = true;
      }

#ifdef OFT_UMFPACK_LONG
      // Copy (possibly updated) graph to 64-bit arrays
      if ( iopt == 1 ) {
        free(LUfactors->Ap);
        free(LUfactors->Ai);
        LUfactors->Ap = (umf_int *) malloc((n+1) * sizeof(umf_int));
        LUfactors->Ai = (umf_int *) malloc(nnz * sizeof(umf_int));
        if ( !LUfactors->Ap || !LUfactors->Ai ) {
          printf("UMFPACK: Unable to allocate 64-bit index arrays\n");
          free(LUfactors->Ap);
          free(LUfactors->Ai);
          LUfactors->Ap = NULL;
          LUfactors->Ai = NULL;
          *info = UMFPACK_ERROR_out_of_memory;
          return;
        }
        for (i = 0; i <= n; i++) LUfactors->Ap[i] = rowptr[i];
        for (i = 0; i < nnz; i++) LUfactors->Ai[i] = colind[i];
      }
      Ap = LUfactors->Ap;
      Ai = LUfactors->Ai;
#else
      Ap = rowptr;
      Ai = colind;
#endif

      // Perform factorization
      if ( iopt == 1 ) {
        if (refactor) {
          UMF_FUN(free_symbolic)( &LUfactors->Symbolic);
        }
        status = UMF_FUN(symbolic)( n, n, Ap, Ai, values, &LUfactors->Symbolic, LUfactors->Control, Info);
        if ( status != 0 ){
          printf("UMFPACK: Symbolic factorization failed %d\n", status);
          *info = status;
          return;
        }
      }
      if (refactor) {
        UMF_FUN(free_numeric)( &LUfactors->Numeric);
      }
      status = UMF_FUN(numeric)( Ap, Ai, values, LUfactors->Symbolic, &LUfactors->Numeric, LUfactors->Control, Info);
      *info = status;
      if ( status != UMFPACK_OK) {
        if ( status == UMFPACK_ERROR_out_of_memory) {
          printf("UMFPACK: Out of Memory in factorization\n");
        } else if ( status == UMFPACK_WARNING_singular_matrix ) {
          printf("UMFPACK: Singular matrix detected in factorization\n");
        } else if ( status == UMFPACK_WARNING_determinant_underflow ) {
          printf("UMFPACK: Determinant underflow detected in factorization\n");
        } else if ( status == UMFPACK_WARNING_determinant_overflow ) {
          printf("UMFPACK: Determinant overflow detected in factorization\n");
        }
      }
      // UMF_FUN(report_info)(LUfactors->Control, Info);

    } else if ( iopt == 2 ) { /* Perform solve */

      LUfactors = *f_factors;
#ifdef OFT_UMFPACK_LONG
      Ap = LUfactors->Ap;
      Ai = LUfactors->Ai;
#else
      Ap = rowptr;
      Ai = colind;
#endif
      double *x;
      x  = (double *) malloc(n * sizeof(double));
      for (j = 0; j < nrhs; j++)  {
        for (i = 0; i < n; i++) x[i] = b[i + j*n];
        // Reversed order of x -> b
        status = UMF_FUN(wsolve)( UMFPACK_At, Ap, Ai, values, b+j*n, x, LUfactors->Numeric, LUfactors->Control, Info, LUfactors->Wi, LUfactors->W);
        if ( status != 0 ) break;
      }
      free(x);
      *info = status;

    } else if ( iopt == 4 ) { /* free storage */

      LUfactors = *f_factors;
      UMF_FUN(free_symbolic)( &LUfactors->Symbolic);
      UMF_FUN(free_numeric)( &LUfactors->Numeric);
      free(LUfactors->Ap);
      free(LUfactors->Ai);
      free(LUfactors->Wi);
      free(LUfactors->W);
      free(LUfactors);
      *f_factors = NULL;

    }
}
#endif
