#ifdef HAVE_SUPERLU_DIST
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified for use with Open FUSION Toolkit by Chris Hansen (December 2014)
 */
#include "superlu_ddefs.h"

void
oft_superlu_dist_slugrid_c(int *iopt, MPI_Comm *slu_comm, int *nprow, int *npcol,
				   gridinfo_t **grid_handle, int *info)
/*
 * This routine provides a fortran call for initializing and
 * freeing the SuperLU_DIST processor grid.  The pointer for the grid
 * structure is returned in grid_handle.
 *
 * The input option, iopt, controls the functionality:
 *   iopt=1:  allocate and define a new process grid
 *   iopt=2:  free an existing process grid
 *
 * slu_comm is the base communication handle
 * nprow is the number of processors per process grid row
 * npcol is the number of processors per process grid column
 */

{
	gridinfo_t *grid;

	if ( *iopt == 1 ) {
		/* Allocate the grid structure. */
		grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));

		/* Initialize the process grid. */
		superlu_gridinit(*slu_comm, *nprow, *npcol, grid);

		/* Set the handle passed from fortran, so that the
		* process grid can be reused. */
		*grid_handle = grid;
		*info = 0;
	} else if ( *iopt == 2 ) {
		/* Locate and free the process grid. */
		grid = *grid_handle;
		superlu_gridexit(grid);
		SUPERLU_FREE(grid);
		*info = 0;
	} else {
		fprintf(stderr, "Invalid iopt=%d passed to oft_superlu_dist_slugrid_c()\n", *iopt);
		*info = -1;
	}
}

typedef struct {
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
	dScalePermstruct_t *ScalePermstruct;
	dLUstruct_t *LUstruct;
#else
	ScalePermstruct_t *ScalePermstruct;
	LUstruct_t *LUstruct;
#endif
	SuperMatrix A;
} factors_dist_t;

void
oft_superlu_dist_dgssv_c(int *iopt, int_t *n, int_t *nnz, int *nrhs,
				double *values, int_t *colind, int_t *rowptr,
				double *b, int *ldb, gridinfo_t **grid_handle,
				factors_dist_t **f_factors, int *perm_spec, int *info)

{
/*
 * Purpose
 * =======
 *
 * This is a Fortran wrapper to use pdgssvx_ABglobal().
 *
 * Arguments
 * =========
 *
 * iopt (input) int
 *	  Specifies the operation to be performed:
 *	  = 1, performs LU decomposition for the first time
 *	  = 2, performs a subsequent LU decomposition for a new matrix
 *		   with the same sparsity pattern
 *	  = 3, performs triangular solve
 *	  = 4, frees all the storage in the end
 *
 * n	(input) int, order of the matrix A
 *
 * nnz  (input) int, number of nonzeros in matrix A
 *
 * nrhs (input) int, number of right-hand sides in the system AX = B
 *
 * values/colind/rowptr (input) column compressed data structure for A
 *
 * b	(input/output) double
 *	  On input, the right-hand side matrix of dimension (ldb, nrhs)
 *	  On output, the solution matrix
 *
 * ldb  (input) int, leading dimension of the matrix B
 *
 * grid_handle (input) int array of size 8, holds a pointer to the process
 *	  grid structure, which is created and freed separately.
 *
 * berr  (output) double, the backward error of each right-hand side
 *
 * factors (input/output) int array of size 8
 *	  If iopt == 1, it is an output and contains the pointer pointing to
 *					the structure of the factored matrices.
 *	  Otherwise, it it an input.
 *
 * info (output) int
 *
 */
	superlu_dist_options_t options;
	SuperLUStat_t stat;
	SuperMatrix A;
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
	dScalePermstruct_t *ScalePermstruct;
	dLUstruct_t *LUstruct;
#else
	ScalePermstruct_t *ScalePermstruct;
	LUstruct_t *LUstruct;
#endif
	gridinfo_t *grid;
	factors_dist_t *LUfactors;
	double *berr;

	/* Locate the process grid. */
	grid = *grid_handle;
	if ( *iopt == 1 || *iopt == 3 ) { /* LU decomposition */
		if ( !*f_factors ) {
            /*
            * Get column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A' * A
            *   permc_spec = 2: minimum degree on structure of A' + A
            *   permc_spec = 3: approximate minimum degree column ordering
            */
            // if ( *perm_spec >= 0 ) options.ColPerm = *perm_spec;
            /* */
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
			ScalePermstruct = (dScalePermstruct_t *) SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
			dScalePermstructInit(*n, *n, ScalePermstruct);
			LUstruct = (dLUstruct_t *) SUPERLU_MALLOC(sizeof(dLUstruct_t));
			dLUstructInit(*n, LUstruct);
#else
			ScalePermstruct = (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
			ScalePermstructInit(*n, *n, ScalePermstruct);
			LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
			LUstructInit(*n, LUstruct);
#endif	
			dCreate_CompCol_Matrix_dist(&A, *n, *n, *nnz, values, colind, rowptr,
						SLU_NC, SLU_D, SLU_GE);
			// Setup saved info
			LUfactors = (factors_dist_t*) SUPERLU_MALLOC(sizeof(factors_dist_t));
			*f_factors = LUfactors;
        } else {
			// Get saved info
            LUfactors = *f_factors;
			// Update matrix
			A = LUfactors->A;
			Destroy_SuperMatrix_Store_dist(&A);
			dCreate_CompCol_Matrix_dist(&A, *n, *n, *nnz, values, colind, rowptr,
						SLU_NC, SLU_D, SLU_GE);
			// Destroy and recreate
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
			dDestroy_LU(*n, grid, LUfactors->LUstruct);
			dLUstructFree(LUfactors->LUstruct);
			SUPERLU_FREE(LUfactors->LUstruct);
			LUstruct = (dLUstruct_t *) SUPERLU_MALLOC(sizeof(dLUstruct_t));
			dLUstructInit(*n, LUstruct);
#else
			Destroy_LU(*n, grid, LUfactors->LUstruct);
			LUstructFree(LUfactors->LUstruct);
			SUPERLU_FREE(LUfactors->LUstruct);
			LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
			LUstructInit(*n, LUstruct);
#endif
			if ( *iopt == 1 ) {
				// Destroy and recreate
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
				dScalePermstructFree(LUfactors->ScalePermstruct);
				SUPERLU_FREE(LUfactors->ScalePermstruct);
				ScalePermstruct = (dScalePermstruct_t *) SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
				dScalePermstructInit(*n, *n, ScalePermstruct);
#else
				ScalePermstructFree(LUfactors->ScalePermstruct);
				SUPERLU_FREE(LUfactors->ScalePermstruct);
				ScalePermstruct = (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
				ScalePermstructInit(*n, *n, ScalePermstruct);
#endif			
			} else {
				ScalePermstruct = LUfactors->ScalePermstruct;
			}
        }

		// Set options
		set_default_options_dist(&options);
		options.Equil = NO;
#if (SUPERLU_DIST_MAJOR_VERSION > 5) || ((SUPERLU_DIST_MAJOR_VERSION == 5) && (SUPERLU_DIST_MINOR_VERSION > 3))
                options.RowPerm = LargeDiag_MC64;
#else
                options.RowPerm = LargeDiag;
#endif
		options.ColPerm = MMD_AT_PLUS_A;
		options.ReplaceTinyPivot = NO;
		options.IterRefine = NO;
		options.Trans = NOTRANS;
		options.SolveInitialized = NO;
		options.RefineInitialized = NO;
		options.PrintStat = NO;
		//
		options.Fact = DOFACT;
		if ( *iopt == 3 ) options.Fact = SamePattern;

		/* Initialize the statistics variables. */
		PStatInit(&stat);

		/* Call global routine with nrhs=0 to perform the factorization. */
		pdgssvx_ABglobal(&options, &A, ScalePermstruct, NULL, *ldb, 0,
						grid, LUstruct, berr, &stat, info);

		LUfactors->ScalePermstruct = ScalePermstruct;
		LUfactors->LUstruct = LUstruct;
		LUfactors->A = A;

		/* Free un-wanted storage */
		PStatFree(&stat);

	} else if ( *iopt == 2 ) { /* Triangular solve */

		/* Extract the LU factors in the factors handle */
		LUfactors = *f_factors;
		ScalePermstruct = LUfactors->ScalePermstruct;
		LUstruct = LUfactors->LUstruct;
		A = LUfactors->A;

		PStatInit(&stat);

		// Set options
		set_default_options_dist(&options);
		options.Equil = NO;
#if (SUPERLU_DIST_MAJOR_VERSION > 5) || ((SUPERLU_DIST_MAJOR_VERSION == 5) && (SUPERLU_DIST_MINOR_VERSION > 3))
		options.RowPerm = LargeDiag_MC64;
#else
		options.RowPerm = LargeDiag;
#endif
		options.ColPerm = MMD_AT_PLUS_A;
		options.ReplaceTinyPivot = NO;
		options.IterRefine = NO;
		options.Trans = NOTRANS;
		options.SolveInitialized = NO;
		options.RefineInitialized = NO;
		options.PrintStat = NO;
		//
		options.Fact = FACTORED;

		berr = (double *) malloc(*n * sizeof(double));

		/* Solve the system A*X=B, overwriting B with X. */
		pdgssvx_ABglobal(&options, &A, ScalePermstruct, b, *ldb, *nrhs,
						grid, LUstruct, berr, &stat, info);

		/* Free un-wanted storage */
		PStatFree(&stat);
		SUPERLU_FREE(berr);

	} else if ( *iopt == 4 ) { /* Free storage */

		/* Free the LU factors in the factors handle */
		LUfactors = *f_factors;
		A = LUfactors->A;
		Destroy_SuperMatrix_Store_dist(&A);
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
		dDestroy_LU(*n, grid, LUfactors->LUstruct);
		dLUstructFree(LUfactors->LUstruct);
		dScalePermstructFree(LUfactors->ScalePermstruct);
#else
		Destroy_LU(*n, grid, LUfactors->LUstruct);
		LUstructFree(LUfactors->LUstruct);
		ScalePermstructFree(LUfactors->ScalePermstruct);
#endif
		SUPERLU_FREE(LUfactors->ScalePermstruct);
		SUPERLU_FREE(LUfactors->LUstruct);
		SUPERLU_FREE(LUfactors);

	} else {
		fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dist_dgsisx_c()\n",*iopt);
		*info = -1;
	}
}
#endif
