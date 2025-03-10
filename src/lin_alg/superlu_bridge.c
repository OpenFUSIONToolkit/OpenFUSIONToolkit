#ifdef HAVE_SUPERLU
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified for use with Open FUSION Toolkit by Chris Hansen (December 2014)
 */
#include <stdbool.h>
#include "slu_ddefs.h"

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
    int *etree;
    superlu_options_t options;
#if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
#endif
} factors_t;

void
oft_superlu_dgssv_c(int iopt, int n, int nnz, int nrhs,
                    double *values, int *colind, int *rowptr,
                    double *b, int ldb,
                    factors_t **f_factors,
                    int perm_spec, bool iter_refine, int *info)
/*
* This routine can be called from Fortran.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs LU decomposition for the first time
*      = 2, performs triangular solve
*      = 3, refactor matrix with the same non-zero pattern
*      = 4, free all the storage in the end
*
* f_factors (input/output) ptr*
*      If iopt == 1, it is an output and contains the pointer pointing to
*                    the structure of the factored matrices.
*      Otherwise, it it an input.
*
*/
{
    SuperMatrix A, AC, B, X;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    char equed;
    double *R;
    double *C;
    SCformat *Lstore;
    NCformat *Ustore;
    int i, panel_size, permc_spec, relax;
    trans_t trans;
    mem_usage_t mem_usage;
    superlu_options_t options;
    #if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
    #endif
    SuperLUStat_t stat;
    factors_t *LUfactors;

    double *recip_pivot_growth;
    double *rcond;
    double *ferr;
    double *berr;
    double *x;

    trans = TRANS; // Row storage format is passed so we solve the transposed system

    if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
        if ( !*f_factors ) {
            /* Set the default input options. */
            set_default_options(&options);
            /*
            * Get column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A'*A
            *   permc_spec = 2: minimum degree on structure of A'+A
            *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
            */
            if ( perm_spec >= 0 ) options.ColPerm = perm_spec;
            if (!iter_refine) options.IterRefine = NO;
            /* */
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
            if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
        } else {
            LUfactors = *f_factors;
            Destroy_SuperNode_Matrix(LUfactors->L);
            Destroy_CompCol_Matrix(LUfactors->U);
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            perm_c = LUfactors->perm_c;
            perm_r = LUfactors->perm_r;
            etree = LUfactors->etree;
            options = LUfactors->options;
#if SUPERLU_VER_MAJOR >= 5
            Glu = LUfactors->Glu;
#endif
        }
        //
        options.Fact = DOFACT;
        if ( iopt == 3 ) options.Fact = SamePattern;
        /* Initialize the statistics variables. */
        StatInit(&stat);
        //
        dCreate_CompCol_Matrix(&A, n, n, nnz, values, colind, rowptr,
                               SLU_NC, SLU_D, SLU_GE);
        if ( iopt == 3 ) {
            //
            sp_preorder(&options, &A, perm_c, etree, &AC);
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
            //
#if SUPERLU_VER_MAJOR >= 5
            dgstrf(&options, &AC, relax, panel_size, etree,
                   NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgstrf(&options, &AC, relax, panel_size, etree,
                   NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        } else {
            permc_spec = options.ColPerm;
            get_perm_c(permc_spec, &A, perm_c);
            //
            sp_preorder(&options, &A, perm_c, etree, &AC);
            //
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
            //
#if SUPERLU_VER_MAJOR >= 5
            dgstrf(&options, &AC, relax, panel_size, etree,
            	   NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgstrf(&options, &AC, relax, panel_size, etree,
            	   NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        }
        //
//        if ( *info == 0 ) {
//            Lstore = (SCformat *) L->Store;
//            Ustore = (NCformat *) U->Store;
//            //printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
//            //printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
//            //printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
//            //dQuerySpace(L, U, &mem_usage);
//            //printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
//            //       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
//        } else {
//            printf("dgstrf() error returns INFO= %d\n", *info);
//            if ( *info <= n ) { /* factorization completes */
//                dQuerySpace(L, U, &mem_usage);
//                printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
//                       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
//            }
//        }

        if ( !*f_factors ) {
            /* Save the LU factors in the factors handle */
        	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
            LUfactors->perm_c = perm_c;
            LUfactors->perm_r = perm_r;
            LUfactors->etree = etree;
            LUfactors->options = options;
#if SUPERLU_VER_MAJOR >= 5
            LUfactors->Glu = Glu;
#endif
            *f_factors = LUfactors;
        }
        LUfactors->L = L;
        LUfactors->U = U;

        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        StatFree(&stat);

    } else if ( iopt == 2 ) { /* Triangular solve */
        /* Initialize the statistics variables. */
        StatInit(&stat);

        /* Extract the LU factors in the factors handle */
        LUfactors = *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;

        dCreate_Dense_Matrix(&B, n, nrhs, b, ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

        Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( iopt == 4 ) { /* Free storage */
        /* Free the LU factors in the factors handle */
        LUfactors = *f_factors;
        SUPERLU_FREE (LUfactors->etree);
        SUPERLU_FREE (LUfactors->perm_r);
        SUPERLU_FREE (LUfactors->perm_c);
        Destroy_SuperNode_Matrix(LUfactors->L);
        Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
        SUPERLU_FREE (LUfactors);
        f_factors = NULL;
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dgssv_c()\n",iopt);
        *info = -1;
    }
}

void
oft_superlu_dgsisx_c(int iopt, int n, int nnz, int nrhs,
                    double *values, int *colind, int *rowptr,
                    double *b, int ldb,
                    factors_t **f_factors,
                    int perm_spec, double fill_tol, int *info)
/*
* This routine can be called from Fortran.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs LU decomposition for the first time
*      = 2, performs triangular solve
*      = 3, refactor matrix with the same non-zero pattern
*      = 4, free all the storage in the end
*
* f_factors (input/output) ptr*
*      If iopt == 1, it is an output and contains the pointer pointing to
*                    the structure of the factored matrices.
*      Otherwise, it it an input.
*
*/
{
    SuperMatrix A, AC, B, X;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    char equed;
    double *R;
    double *C;
    SCformat *Lstore;
    NCformat *Ustore;
    int i, panel_size, permc_spec, relax;
    trans_t trans;
    mem_usage_t mem_usage;
    superlu_options_t options;
    #if SUPERLU_VER_MAJOR >= 5
    GlobalLU_t Glu;
    #endif
    SuperLUStat_t stat;
    factors_t *LUfactors;

    double *recip_pivot_growth;
    double *rcond;
    double *ferr;
    double *berr;
    double *x;

    trans = TRANS; // Row storage format is passed so we solve the transposed system

    if ( iopt == 1 || iopt == 3 ) { /* iLU decomposition */
        if ( !*f_factors ) {
            /* Set the default input options. */
            ilu_set_default_options(&options);
            // options.DiagPivotThresh = 1.0;
            // options.RowPerm = NO;
            options.ILU_MILU = SMILU_2;
            options.ILU_DropTol = fill_tol; //1.e-7;
            // options.ILU_FillTol = 1.e-3;
            /*
            * Get column permutation vector perm_c[], according to permc_spec:
            *   permc_spec = 0: natural ordering
            *   permc_spec = 1: minimum degree on structure of A'*A
            *   permc_spec = 2: minimum degree on structure of A'+A
            *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
            */
            if ( perm_spec >= 0 ) options.ColPerm = perm_spec;
            /* */
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
            if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
        } else {
            LUfactors = *f_factors;
            Destroy_SuperNode_Matrix(LUfactors->L);
            Destroy_CompCol_Matrix(LUfactors->U);
            L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
            perm_c = LUfactors->perm_c;
            perm_r = LUfactors->perm_r;
            etree = LUfactors->etree;
            options = LUfactors->options;
#if SUPERLU_VER_MAJOR >= 5
            Glu = LUfactors->Glu;
#endif
        }
        //
        options.Fact = DOFACT;
        if ( iopt == 3 ) options.Fact = SamePattern;
        /* Initialize the statistics variables. */
        StatInit(&stat);
        //
        dCreate_CompCol_Matrix(&A, n, n, nnz, values, colind, rowptr,
                                SLU_NC, SLU_D, SLU_GE);
        if ( iopt == 3 ) {
            sp_preorder(&options, &A, perm_c, etree, &AC);
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
#if SUPERLU_VER_MAJOR >= 5
            dgsitrf(&options, &AC, relax, panel_size, etree,
                    NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgsitrf(&options, &AC, relax, panel_size, etree,
                    NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        } else {
            permc_spec = options.ColPerm;
            get_perm_c(permc_spec, &A, perm_c);
            sp_preorder(&options, &A, perm_c, etree, &AC);
            panel_size = sp_ienv(1);
            relax = sp_ienv(2);
#if SUPERLU_VER_MAJOR >= 5
            dgsitrf(&options, &AC, relax, panel_size, etree,
                    NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);
#else
            dgsitrf(&options, &AC, relax, panel_size, etree,
                    NULL, 0, perm_c, perm_r, L, U, &stat, info);
#endif
        }
        //
        // if ( *info == 0 ) {
        //     Lstore = (SCformat *) L->Store;
        //     Ustore = (NCformat *) U->Store;
        //     printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
        //     printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
        //     printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
        //     dQuerySpace(L, U, &mem_usage);
        //     printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
        //             mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        // } else {
        //     printf("dgsitrf() error returns INFO= %d\n", *info);
        //     if ( *info <= n ) { /* factorization completes */
        //         dQuerySpace(L, U, &mem_usage);
        //         printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
        //                 mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        //     }
        // }

        if ( !*f_factors ) {
            /* Save the LU factors in the factors handle */
            LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
            LUfactors->perm_c = perm_c;
            LUfactors->perm_r = perm_r;
            LUfactors->etree = etree;
            LUfactors->options = options;
#if SUPERLU_VER_MAJOR >= 5
            LUfactors->Glu = Glu;
#endif
            *f_factors = LUfactors;
        }
        LUfactors->L = L;
        LUfactors->U = U;

        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        StatFree(&stat);

    } else if ( iopt == 2 ) { /* Triangular solve */
        /* Initialize the statistics variables. */
        StatInit(&stat);

        /* Extract the LU factors in the factors handle */
        LUfactors = *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;

        dCreate_Dense_Matrix(&B, n, nrhs, b, ldb, SLU_DN, SLU_D, SLU_GE);

        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs(trans, L, U, perm_c, perm_r, &B, &stat, info);

        Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( iopt == 4 ) { /* Free storage */
        /* Free the LU factors in the factors handle */
        LUfactors = *f_factors;
        SUPERLU_FREE (LUfactors->etree);
        SUPERLU_FREE (LUfactors->perm_r);
        SUPERLU_FREE (LUfactors->perm_c);
        Destroy_SuperNode_Matrix(LUfactors->L);
        Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
        SUPERLU_FREE (LUfactors);
        f_factors = NULL;
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to oft_superlu_dgsisx_c()\n",iopt);
        *info = -1;
    }
}
#endif
