/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* SuperLU-DIST LU solver interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_SUPERLU_DIST
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified for use with Open FUSION Toolkit by Chris Hansen (December 2014)
 */
#include <stdbool.h>
#include "superlu_ddefs.h"

void
oft_superlu_dist_slugrid_c(int iopt, int *slu_comm, int nprow, int npcol,
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

	if ( iopt == 1 ) {
		/* Allocate the grid structure. */
		grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));

		/* Initialize the process grid. */
		MPI_Comm cComm = MPI_Comm_f2c(*slu_comm);
		superlu_gridinit(cComm, nprow, npcol, grid);

		/* Set the handle passed from fortran, so that the
		* process grid can be reused. */
		*grid_handle = grid;
		*info = 0;
	} else if ( iopt == 2 ) {

		if ( !*grid_handle ) {
			*info = -2;
			return;
		}

		/* Locate and free the process grid. */
		grid = *grid_handle;
		superlu_gridexit(grid);
		SUPERLU_FREE(grid);
		*grid_handle = NULL;
		*info = 0;
	} else {
		printf("Invalid iopt=%d passed to oft_superlu_dist_slugrid_c()\n", iopt);
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
oft_superlu_dist_dgssv_c(int iopt, int_t n, int_t nnz, int nrhs,
				double *values, int_t *colind, int_t *rowptr,
				double *b, int ldb, gridinfo_t *grid,
				factors_dist_t **f_factors, int perm_spec, bool iter_refine, int *info)
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
{
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
	factors_dist_t *LUfactors;

	/* Locate the process grid. */
	if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
		if ( !*f_factors ) {
            /*
            * Get column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A' * A
            *   permc_spec = 2: minimum degree on structure of A' + A
            *   permc_spec = 3: approximate minimum degree column ordering
            */
            // if ( *perm_spec >= 0 ) options.ColPerm = *perm_spec;
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
			ScalePermstruct = (dScalePermstruct_t *) SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
			dScalePermstructInit(n, n, ScalePermstruct);
			LUstruct = (dLUstruct_t *) SUPERLU_MALLOC(sizeof(dLUstruct_t));
			dLUstructInit(n, LUstruct);
#else
			ScalePermstruct = (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
			ScalePermstructInit(n, n, ScalePermstruct);
			LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
			LUstructInit(n, LUstruct);
#endif	
			dCreate_CompCol_Matrix_dist(&A, n, n, nnz, values, colind, rowptr,
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
			dCreate_CompCol_Matrix_dist(&A, n, n, nnz, values, colind, rowptr,
						SLU_NC, SLU_D, SLU_GE);
			// Destroy and recreate
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
			dDestroy_LU(n, grid, LUfactors->LUstruct);
			dLUstructFree(LUfactors->LUstruct);
			SUPERLU_FREE(LUfactors->LUstruct);
			LUstruct = (dLUstruct_t *) SUPERLU_MALLOC(sizeof(dLUstruct_t));
			dLUstructInit(n, LUstruct);
#else
			Destroy_LU(n, grid, LUfactors->LUstruct);
			LUstructFree(LUfactors->LUstruct);
			SUPERLU_FREE(LUfactors->LUstruct);
			LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
			LUstructInit(n, LUstruct);
#endif
			if ( iopt == 1 ) {
				// Destroy and recreate
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
				dScalePermstructFree(LUfactors->ScalePermstruct);
				SUPERLU_FREE(LUfactors->ScalePermstruct);
				ScalePermstruct = (dScalePermstruct_t *) SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
				dScalePermstructInit(n, n, ScalePermstruct);
#else
				ScalePermstructFree(LUfactors->ScalePermstruct);
				SUPERLU_FREE(LUfactors->ScalePermstruct);
				ScalePermstruct = (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
				ScalePermstructInit(n, n, ScalePermstruct);
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
		if (!iter_refine) options.IterRefine = NO;
		options.Trans = NOTRANS;
		options.SolveInitialized = NO;
		options.RefineInitialized = NO;
		options.PrintStat = NO;
		//
		options.Fact = DOFACT;
		if ( iopt == 3 ) options.Fact = SamePattern;

		/* Initialize the statistics variables. */
		PStatInit(&stat);

		/* Call global routine with nrhs=0 to perform the factorization. */
		double *berr = NULL;
		pdgssvx_ABglobal(&options, &A, ScalePermstruct, NULL, ldb, 0,
						grid, LUstruct, berr, &stat, info);

		LUfactors->ScalePermstruct = ScalePermstruct;
		LUfactors->LUstruct = LUstruct;
		LUfactors->A = A;

		/* Free un-wanted storage */
		PStatFree(&stat);

	} else if ( iopt == 2 ) { /* Triangular solve */

		if ( !*f_factors ) {
			*info = -2;
			return;
		}

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
		if (!iter_refine) options.IterRefine = NO;
		options.Trans = NOTRANS;
		options.SolveInitialized = NO;
		options.RefineInitialized = NO;
		options.PrintStat = NO;
		//
		options.Fact = FACTORED;

		double *berr = (double *) malloc(nrhs * sizeof(double));

		/* Solve the system A*X=B, overwriting B with X. */
		pdgssvx_ABglobal(&options, &A, ScalePermstruct, b, ldb, nrhs,
						grid, LUstruct, berr, &stat, info);

		/* Free un-wanted storage */
		PStatFree(&stat);
		SUPERLU_FREE(berr);

	} else if ( iopt == 4 ) { /* Free storage */

		if ( !*f_factors ) {
			*info = -2;
			return;
		}

		/* Free the LU factors in the factors handle */
		LUfactors = *f_factors;
		A = LUfactors->A;
		Destroy_SuperMatrix_Store_dist(&A);
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
		dDestroy_LU(n, grid, LUfactors->LUstruct);
		dLUstructFree(LUfactors->LUstruct);
		dScalePermstructFree(LUfactors->ScalePermstruct);
#else
		Destroy_LU(n, grid, LUfactors->LUstruct);
		LUstructFree(LUfactors->LUstruct);
		ScalePermstructFree(LUfactors->ScalePermstruct);
#endif
		SUPERLU_FREE(LUfactors->ScalePermstruct);
		SUPERLU_FREE(LUfactors->LUstruct);
		SUPERLU_FREE(LUfactors);
		*f_factors = NULL;

	} else {
		fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dist_dgssv_c()\n",iopt);
		*info = -1;
	}
}

/* Type and routine names for distributed (NR_loc) solver structures */
#if (SUPERLU_DIST_MAJOR_VERSION > 6) || ((SUPERLU_DIST_MAJOR_VERSION == 6) && (SUPERLU_DIST_MINOR_VERSION > 2))
#define OFT_SLUD_ScalePermstruct_t dScalePermstruct_t
#define OFT_SLUD_LUstruct_t dLUstruct_t
#define OFT_SLUD_SOLVEstruct_t dSOLVEstruct_t
#define OFT_SLUD_ScalePermstructInit dScalePermstructInit
#define OFT_SLUD_ScalePermstructFree dScalePermstructFree
#define OFT_SLUD_LUstructInit dLUstructInit
#define OFT_SLUD_LUstructFree dLUstructFree
#define OFT_SLUD_Destroy_LU dDestroy_LU
#else
#define OFT_SLUD_ScalePermstruct_t ScalePermstruct_t
#define OFT_SLUD_LUstruct_t LUstruct_t
#define OFT_SLUD_SOLVEstruct_t SOLVEstruct_t
#define OFT_SLUD_ScalePermstructInit ScalePermstructInit
#define OFT_SLUD_ScalePermstructFree ScalePermstructFree
#define OFT_SLUD_LUstructInit LUstructInit
#define OFT_SLUD_LUstructFree LUstructFree
#define OFT_SLUD_Destroy_LU Destroy_LU
#endif

typedef struct {
	int_t n;
	superlu_dist_options_t options;
	OFT_SLUD_ScalePermstruct_t ScalePermstruct;
	OFT_SLUD_LUstruct_t LUstruct;
	OFT_SLUD_SOLVEstruct_t SOLVEstruct;
	SuperMatrix A;
	bool A_set;
} factors_dist_loc_t;

static void oft_superlu_dist_loc_free(factors_dist_loc_t *LUfactors, gridinfo_t *grid)
{
	if ( LUfactors->A_set ) Destroy_CompRowLoc_Matrix_dist(&LUfactors->A);
	OFT_SLUD_Destroy_LU(LUfactors->n, grid, &LUfactors->LUstruct);
	OFT_SLUD_ScalePermstructFree(&LUfactors->ScalePermstruct);
	OFT_SLUD_LUstructFree(&LUfactors->LUstruct);
	if ( LUfactors->options.SolveInitialized == YES ) {
		dSolveFinalize(&LUfactors->options, &LUfactors->SOLVEstruct);
	}
	SUPERLU_FREE(LUfactors);
}

void
oft_superlu_dist_pdgssvx_c(int iopt, int_t n, int_t nnz_loc, int_t m_loc, int_t fst_row,
				int nrhs, double *values, int_t *colind, int_t *rowptr,
				double *b, int ldb, gridinfo_t *grid,
				factors_dist_loc_t **f_factors, bool iter_refine, int *info)
/*
 * Purpose
 * =======
 *
 * Fortran wrapper to use pdgssvx() with a matrix distributed by rows
 * across all processes in the grid (SLU_NR_loc format).
 *
 * Arguments
 * =========
 *
 * iopt (input) int
 *	  Specifies the operation to be performed:
 *	  = 1, performs LU decomposition (new sparsity pattern)
 *	  = 2, performs triangular solve
 *	  = 3, performs LU decomposition for a new matrix with the same sparsity pattern
 *	       (orderings are recomputed, see below)
 *	  = 4, frees all the storage in the end
 *
 * n	(input) int, global order of the matrix A
 *
 * nnz_loc (input) int, number of nonzeros in the local rows of A
 *
 * m_loc (input) int, number of rows of A stored on this process
 *
 * fst_row (input) int, global index (0-based) of the first local row; the
 *	  local rows of each process must be contiguous and ordered by rank
 *
 * nrhs (input) int, number of right-hand sides in the system AX = B
 *
 * values/colind/rowptr (input) compressed row structure for the local rows of A
 *	  using 1-based (Fortran) indexing and global column indices. These arrays
 *	  are copied, as pdgssvx() overwrites the input matrix.
 *
 * b	(input/output) double
 *	  On input, the local rows of the right-hand side matrix (ldb, nrhs)
 *	  On output, the local rows of the solution matrix
 *
 * ldb  (input) int, leading dimension of the matrix B (>= m_loc)
 *
 * grid (input) SuperLU-DIST process grid
 *
 * f_factors (input/output) Pointer to LU factors and supporting structures
 *
 * iter_refine (input) bool, perform iterative refinement?
 *
 * info (output) int, solver return status
 *
 */
{
	factors_dist_loc_t *LUfactors;
	SuperLUStat_t stat;
	double *berr;
	double *a_vals;
	int_t *a_colind, *a_rowptr;
	int_t i;

	*info = 0;
	if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
		LUfactors = *f_factors;
		if ( LUfactors && iopt == 1 ) { /* New pattern: discard everything */
			oft_superlu_dist_loc_free(LUfactors, grid);
			LUfactors = NULL;
			*f_factors = NULL;
		}
		if ( !LUfactors ) {
			LUfactors = (factors_dist_loc_t *) SUPERLU_MALLOC(sizeof(factors_dist_loc_t));
			LUfactors->n = n;
			LUfactors->A_set = false;
			set_default_options_dist(&LUfactors->options);
			LUfactors->options.Equil = NO;
#if (SUPERLU_DIST_MAJOR_VERSION > 5) || ((SUPERLU_DIST_MAJOR_VERSION == 5) && (SUPERLU_DIST_MINOR_VERSION > 3))
			LUfactors->options.RowPerm = LargeDiag_MC64;
#else
			LUfactors->options.RowPerm = LargeDiag;
#endif
			LUfactors->options.ReplaceTinyPivot = NO;
			if (!iter_refine) LUfactors->options.IterRefine = NOREFINE;
			LUfactors->options.Trans = NOTRANS;
			LUfactors->options.SolveInitialized = NO;
			LUfactors->options.RefineInitialized = NO;
			LUfactors->options.PrintStat = NO;
			LUfactors->options.Fact = DOFACT;
			OFT_SLUD_ScalePermstructInit(n, n, &LUfactors->ScalePermstruct);
			OFT_SLUD_LUstructInit(n, &LUfactors->LUstruct);
			*f_factors = LUfactors;
		} else { /* Same pattern: reuse solver structures, but not orderings */
			OFT_SLUD_Destroy_LU(LUfactors->n, grid, &LUfactors->LUstruct);
			/* Fact = SamePattern (reuse of perm_c and etree) gives incorrect factors
			   when MC64 selects a different row permutation for the new values
			   (observed with SuperLU-DIST 8.2.1 and 9.2.1), so recompute all
			   orderings. pdgssvx() re-initializes the solve structures. */
			LUfactors->options.Fact = DOFACT;
		}

		/* Create local matrix from copies of the input arrays */
		if ( LUfactors->A_set ) Destroy_CompRowLoc_Matrix_dist(&LUfactors->A);
		a_vals = doubleMalloc_dist(nnz_loc);
		a_colind = intMalloc_dist(nnz_loc);
		a_rowptr = intMalloc_dist(m_loc+1);
		for ( i = 0; i < nnz_loc; i++ ) {
			a_vals[i] = values[i];
			a_colind[i] = colind[i] - 1;
		}
		for ( i = 0; i <= m_loc; i++ ) a_rowptr[i] = rowptr[i] - 1;
		dCreate_CompRowLoc_Matrix_dist(&LUfactors->A, n, n, nnz_loc, m_loc, fst_row,
						a_vals, a_colind, a_rowptr, SLU_NR_loc, SLU_D, SLU_GE);
		LUfactors->A_set = true;

		/* Call global routine with nrhs=0 to perform the factorization. */
		PStatInit(&stat);
		berr = NULL;
		pdgssvx(&LUfactors->options, &LUfactors->A, &LUfactors->ScalePermstruct, NULL, ldb, 0,
				grid, &LUfactors->LUstruct, &LUfactors->SOLVEstruct, berr, &stat, info);
		PStatFree(&stat);

	} else if ( iopt == 2 ) { /* Triangular solve */

		if ( !*f_factors ) {
			*info = -2;
			return;
		}
		LUfactors = *f_factors;
		LUfactors->options.Fact = FACTORED;

		/* Solve the system A*X=B, overwriting B with X. */
		PStatInit(&stat);
		berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double));
		pdgssvx(&LUfactors->options, &LUfactors->A, &LUfactors->ScalePermstruct, b, ldb, nrhs,
				grid, &LUfactors->LUstruct, &LUfactors->SOLVEstruct, berr, &stat, info);
		SUPERLU_FREE(berr);
		PStatFree(&stat);

	} else if ( iopt == 4 ) { /* Free storage */

		if ( !*f_factors ) {
			*info = -2;
			return;
		}
		oft_superlu_dist_loc_free(*f_factors, grid);
		*f_factors = NULL;

	} else {
		fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dist_pdgssvx_c()\n",iopt);
		*info = -1;
	}
}
#endif
