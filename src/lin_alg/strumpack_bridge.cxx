/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* STRUMPACK sparse (multifrontal) LU solver interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_STRUMPACK
#include <cstdio>
#include <exception>
#include <memory>
#include <vector>
#include "StrumpackSparseSolver.hpp"

using namespace strumpack;

typedef struct {
    std::unique_ptr<SparseSolver<double,int>> sp; // Sparse (multifrontal) factorization
} factors_t;

extern "C" void
oft_strumpack_sparse_c(int iopt, int n, int nnz, int nrhs,
                       double *values, int *colind, int *rowptr,
                       double *b, int ldb,
                       factors_t **f_factors,
                       int comp_type, double rel_tol, double abs_tol,
                       int leaf_size, int max_rank, int min_sep_size,
                       bool iter_refine, int verbose, int *info)
/*
* This routine can be called from Fortran.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs reordering and LU decomposition for the first time
*      = 2, performs solve using existing factorization
*      = 3, refactor matrix with the same non-zero pattern
*      = 4, free all the storage in the end
*
* values, colind, rowptr (input)
*      Matrix in CSR format with 0-based indices (copied by STRUMPACK)
*
* comp_type (input) int
*      Compression type for frontal matrices: 0 -> none, 1 -> HSS,
*      2 -> BLR, -1 -> STRUMPACK default (none)
*
* rel_tol, abs_tol, leaf_size, max_rank, min_sep_size (input)
*      Compression options, values <= 0 select the STRUMPACK default
*
* iter_refine (input) bool
*      Use STRUMPACK's default outer solver (iterative refinement, or
*      GMRES when compression is used) instead of a single direct solve
*
* f_factors (input/output) ptr*
*      If iopt == 1, it is an output and contains the pointer pointing to
*                    the structure of the factored matrices.
*      Otherwise, it it an input.
*
*/
{
    ReturnCode ierr;
    *info = 0;
    if ( iopt == 1 ) { /* Reordering and LU decomposition */
        if ( !*f_factors ) *f_factors = new factors_t;
        try {
            (*f_factors)->sp.reset(new SparseSolver<double,int>(verbose != 0));
            SPOptions<double>& opts = (*f_factors)->sp->options();
            if ( comp_type == 0 ) {
                opts.set_compression(CompressionType::NONE);
            } else if ( comp_type == 1 ) {
                opts.set_compression(CompressionType::HSS);
            } else if ( comp_type == 2 ) {
                opts.set_compression(CompressionType::BLR);
            } else if ( comp_type != -1 ) {
                fprintf(stderr,"Invalid comp_type=%d passed to oft_strumpack_sparse_c()\n",comp_type);
                (*f_factors)->sp.reset();
                *info = -1;
                return;
            }
            if ( rel_tol > 0.0 ) opts.set_compression_rel_tol(rel_tol);
            if ( abs_tol > 0.0 ) opts.set_compression_abs_tol(abs_tol);
            if ( leaf_size > 0 ) opts.set_compression_leaf_size(leaf_size);
            if ( max_rank > 0 ) {
                opts.HSS_options().set_max_rank(max_rank);
                opts.BLR_options().set_max_rank(max_rank);
            }
            if ( min_sep_size > 0 ) opts.set_compression_min_sep_size(min_sep_size);
            if ( !iter_refine ) opts.set_Krylov_solver(KrylovSolver::DIRECT);
            (*f_factors)->sp->set_csr_matrix(n, rowptr, colind, values, false);
            ierr = (*f_factors)->sp->factor();
            if ( ierr != ReturnCode::SUCCESS ) {
                fprintf(stderr,"STRUMPACK factorization failed: %s\n",
                        ierr == ReturnCode::REORDERING_ERROR ? "reordering error" : "zero pivot");
                *info = 1;
            }
        } catch (std::exception& e) {
            fprintf(stderr,"STRUMPACK factorization failed: %s\n",e.what());
            (*f_factors)->sp.reset();
            *info = 1;
        }
    } else if ( iopt == 3 ) { /* Refactor with same non-zero pattern */
        if ( !*f_factors || !(*f_factors)->sp ) {
            fprintf(stderr,"Refactor requested before factorization in oft_strumpack_sparse_c()\n");
            *info = -1;
            return;
        }
        try {
            (*f_factors)->sp->update_matrix_values(n, rowptr, colind, values, false);
            ierr = (*f_factors)->sp->factor();
            if ( ierr != ReturnCode::SUCCESS ) {
                fprintf(stderr,"STRUMPACK refactorization failed: %s\n",
                        ierr == ReturnCode::MATRIX_NOT_SET ? "matrix not set" : "zero pivot");
                *info = 1;
            }
        } catch (std::exception& e) {
            fprintf(stderr,"STRUMPACK refactorization failed: %s\n",e.what());
            *info = 1;
        }
    } else if ( iopt == 2 ) { /* Solve */
        if ( !*f_factors || !(*f_factors)->sp ) {
            fprintf(stderr,"Solve requested before factorization in oft_strumpack_sparse_c()\n");
            *info = -1;
            return;
        }
        try {
            /* Solve the system A*X=B, overwriting B with X (RHS copy required) */
            std::vector<double> rhs(b, b + (size_t)ldb*nrhs);
            ierr = (*f_factors)->sp->solve(nrhs, rhs.data(), ldb, b, ldb);
            if ( ierr != ReturnCode::SUCCESS ) {
                fprintf(stderr,"STRUMPACK solve failed: %s\n",
                        ierr == ReturnCode::NO_CONVERGENCE ? "no convergence" : "matrix not set");
                *info = 1;
            }
        } catch (std::exception& e) {
            fprintf(stderr,"STRUMPACK solve failed: %s\n",e.what());
            *info = 1;
        }
    } else if ( iopt == 4 ) { /* Free storage */
        delete *f_factors;
        *f_factors = NULL;
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to oft_strumpack_sparse_c()\n",iopt);
        *info = -1;
    }
}
#endif
