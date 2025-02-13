#ifdef HAVE_MKL
/*
 *
 */

#include <stdlib.h>
#include "mkl_spblas.h"
#include "mkl_rci.h"

typedef struct {
	double *ilu_vals;
	sparse_matrix_t A_iLU;
} factors_t;

void
oft_mkl_ilu_c(int iopt, int n, int nnz, double *values, int *colind,
			  int *rowptr, double *b, factors_t **f_factors, int *info)

{
/*
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs iLU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) ptr*
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
factors_t *LUfactors;
if ( iopt == 1 ) { /* iLU0 decomposition */
	// Create or unpack storage structure
	if ( !*f_factors ) {
		LUfactors = (factors_t*) malloc(sizeof(factors_t));
		LUfactors->ilu_vals = (double *) malloc(nnz * sizeof(double));
		LUfactors->A_iLU = NULL;
		*f_factors = LUfactors;
	} else {
		LUfactors = *f_factors;
	}
	// Setup default settings
	int ipar[128] = {0};
	double dpar[128] = {0.0};
	ipar[30] = 1;
	dpar[30] = 1.E-20;
	dpar[31] = 1.E-16;
	// Factor matrix
	dcsrilu0(&n, values, rowptr, colind, LUfactors->ilu_vals, ipar, dpar, info);
	// Destroy old matrix if necessary
	if (LUfactors->A_iLU)
		mkl_sparse_destroy(LUfactors->A_iLU);
	// Construct factored matrix obect for backsolves
	mkl_sparse_d_create_csr(&LUfactors->A_iLU, SPARSE_INDEX_BASE_ONE, n, n, rowptr, rowptr+1, colind, LUfactors->ilu_vals);
	// Optimize matrix representation
	struct matrix_descr descr_A;
	descr_A.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr_A.mode = SPARSE_FILL_MODE_LOWER;
	descr_A.diag = SPARSE_DIAG_UNIT;
	mkl_sparse_set_mv_hint(LUfactors->A_iLU, SPARSE_OPERATION_NON_TRANSPOSE, descr_A, 100);
	descr_A.mode = SPARSE_FILL_MODE_UPPER;
	descr_A.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_set_mv_hint(LUfactors->A_iLU, SPARSE_OPERATION_NON_TRANSPOSE, descr_A, 100);
	mkl_sparse_optimize(LUfactors->A_iLU);
} else if ( iopt == 2 ) { /* Perform triangular solve */
	LUfactors = *f_factors;
	// Create temporary working vector
	int i;
	double *x;
	x = (double *) malloc(n * sizeof(double));
	for (i = 0; i < n; i++) x[i] = b[i];
	// Apply triangular solve
	struct matrix_descr descr_A;
	descr_A.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr_A.mode = SPARSE_FILL_MODE_LOWER;
	descr_A.diag = SPARSE_DIAG_UNIT;
	mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, LUfactors->A_iLU, descr_A, b, x);
	descr_A.mode = SPARSE_FILL_MODE_UPPER;
	descr_A.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, LUfactors->A_iLU, descr_A, x, b);
	// Free temporary storage
	free(x);
} else if ( iopt == 3 ) { /* Free storage */
	LUfactors = *f_factors;
	if (LUfactors->A_iLU)
		mkl_sparse_destroy(LUfactors->A_iLU);
	free(LUfactors->ilu_vals);
	free(LUfactors);
	*f_factors = NULL;
}

}
#endif