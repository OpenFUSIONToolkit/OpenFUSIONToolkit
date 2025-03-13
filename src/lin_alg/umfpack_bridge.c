#ifdef HAVE_UMFPACK
#include "umfpack.h"
#include <stdbool.h>

/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */

typedef struct {
  void *Symbolic;
  void *Numeric;
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
 *      = 3, free all the storage in the end
 *      = 4, performs numerical factorization only
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
  bool refactor;

    if ( iopt == 1 || iopt == 3 ) { /* Perform full factorization */

      // Create or unpack storage structure
      if ( !*f_factors ) {
        LUfactors = (factors_t*) malloc(sizeof(factors_t));
        /* Set the default input options. */
        umfpack_di_defaults( LUfactors->Control);
        if (!iter_refine) LUfactors->Control[UMFPACK_IRSTEP] = 0;
        *f_factors = LUfactors;
        refactor = false;
      } else {
        LUfactors = *f_factors;
        refactor = true;
      }

      // Perform factorization
      if ( iopt == 1 ) {
        if (refactor) {
          umfpack_di_free_symbolic( &LUfactors->Symbolic);
        }
        status = umfpack_di_symbolic( n, n, rowptr, colind, values, &LUfactors->Symbolic, LUfactors->Control, Info);
        if ( status != 0 ){
          *info = status;
          return;
        }
      }
      if (refactor) {
        umfpack_di_free_numeric( &LUfactors->Numeric);
      }
      status = umfpack_di_numeric( rowptr, colind, values, LUfactors->Symbolic, &LUfactors->Numeric, LUfactors->Control, Info);
      *info = status;

    } else if ( iopt == 2 ) { /* Perform solve */

      LUfactors = *f_factors;
      double *x;
      x  = (double *) malloc(n * sizeof(double));
      for (j = 0; j < nrhs; j++)  {
        for (i = 0; i < n; i++) x[i] = b[i + j*n];
        // Reversed order of x -> b
        status = umfpack_di_solve( UMFPACK_At, rowptr, colind, values, b+j*n, x, LUfactors->Numeric, LUfactors->Control, Info);
        if ( status != 0 ) break;
      }
      free(x);
      *info = status;

    } else if ( iopt == 4 ) { /* free storage */

      LUfactors = *f_factors;
      umfpack_di_free_symbolic( &LUfactors->Symbolic);
      umfpack_di_free_numeric( &LUfactors->Numeric);
      free(LUfactors);
      *f_factors = NULL;

    }
}
#endif
