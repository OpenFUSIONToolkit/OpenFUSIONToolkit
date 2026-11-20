/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* STRUMPACK sparse (multifrontal) LU solver interface for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#ifdef HAVE_STRUMPACK
#include <exception>
#include <iostream>
#include <memory>
#include <vector>
#include "StrumpackSparseSolver.hpp"
#ifdef HAVE_STRUMPACK_MPI
#include "StrumpackSparseSolverMPIDist.hpp"
#endif

using namespace strumpack;

/* Solver options used for a factorization */
struct options_t {
    int comp_type;
    int leaf_size;
    int max_rank;
    int min_sep_size;
    int verbose;
    double rel_tol;
    double abs_tol;
    bool operator==(const options_t& o) const {
        return comp_type == o.comp_type && leaf_size == o.leaf_size && max_rank == o.max_rank
            && min_sep_size == o.min_sep_size && verbose == o.verbose
            && rel_tol == o.rel_tol && abs_tol == o.abs_tol;
    }
};

/* Factorization and the matrix size/options it was computed with */
struct factors_t {
    std::unique_ptr<SparseSolver<double,int>> sp;
    int n = 0;
    int nnz = 0;
    options_t opts = {};
};

/* Convert compression type, returns false for invalid options */
static bool get_compression(const options_t& o, CompressionType& comp) {
    if ( o.comp_type == 0 || o.comp_type == -1 ) {
        comp = CompressionType::NONE;
    } else if ( o.comp_type == 1 ) {
        comp = CompressionType::HSS;
    } else if ( o.comp_type == 2 ) {
        comp = CompressionType::BLR;
    } else {
        std::cerr << "Invalid comp_type=" << o.comp_type << " passed to STRUMPACK interface" << std::endl;
        return false;
    }
    return true;
}

/* Apply options to a STRUMPACK solver */
static void apply_options(SPOptions<double>& opts, const options_t& o, CompressionType comp) {
    opts.set_compression(comp);
    if ( o.rel_tol > 0.0 ) opts.set_compression_rel_tol(o.rel_tol);
    if ( o.abs_tol > 0.0 ) opts.set_compression_abs_tol(o.abs_tol);
    if ( o.leaf_size > 0 ) opts.set_compression_leaf_size(o.leaf_size);
    if ( o.max_rank > 0 ) {
        opts.HSS_options().set_max_rank(o.max_rank);
        opts.BLR_options().set_max_rank(o.max_rank);
    }
    if ( o.min_sep_size > 0 ) opts.set_compression_min_sep_size(o.min_sep_size);
    /* Single application of the (approximate) factorization, no outer iterative solver */
    opts.set_Krylov_solver(KrylovSolver::DIRECT);
}

/* Create a new solver in f and apply options, returns false for invalid options */
static bool setup_solver(factors_t *f, const options_t& o) {
    CompressionType comp;
    if ( !get_compression(o, comp) ) return false;
    f->sp.reset(new SparseSolver<double,int>(o.verbose != 0));
    apply_options(f->sp->options(), o, comp);
    f->opts = o;
    return true;
}

extern "C" void
oft_strumpack_sparse_c(int iopt, int n, int nnz, int nrhs,
                       double *values, int *colind, int *rowptr,
                       double *b, int ldb,
                       factors_t **f_factors,
                       int comp_type, double rel_tol, double abs_tol,
                       int leaf_size, int max_rank, int min_sep_size,
                       int verbose, int *info)
/*
* This routine can be called from Fortran.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs reordering and LU decomposition for the first time
*      = 2, performs solve using existing factorization
*      = 3, refactor matrix with the same non-zero pattern (reuses the existing
*           reordering and MC64 matching/scaling, a full factorization is
*           performed instead if the solver options have changed)
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
* f_factors (input/output) ptr*
*      If iopt == 1, it is an output and contains the pointer pointing to
*                    the structure of the factored matrices.
*      Otherwise, it it an input.
*
*/
{
    const options_t o = {comp_type, leaf_size, max_rank, min_sep_size, verbose, rel_tol, abs_tol};
    *info = 0;
    try {
        if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
            if ( !*f_factors ) *f_factors = new factors_t;
            factors_t *f = *f_factors;
            const bool full = ( iopt == 1 ) || !f->sp || !( f->opts == o );
            if ( full ) {
                if ( !setup_solver(f, o) ) {
                    f->sp.reset();
                    *info = -1;
                    return;
                }
                f->sp->set_csr_matrix(n, rowptr, colind, values, false);
                f->n = n;
                f->nnz = nnz;
            } else {
                /* STRUMPACK does not report a changed non-zero pattern, so check size here */
                if ( n != f->n || nnz != f->nnz ) {
                    std::cerr << "Matrix size changed on refactor in oft_strumpack_sparse_c()" << std::endl;
                    *info = -1;
                    return;
                }
                f->sp->update_matrix_values(n, rowptr, colind, values, false);
            }
            ReturnCode ierr = f->sp->factor();
            if ( ierr != ReturnCode::SUCCESS ) {
                std::cerr << "STRUMPACK factorization failed: " << ierr << std::endl;
                *info = 1;
            }
        } else if ( iopt == 2 ) { /* Solve */
            factors_t *f = *f_factors;
            if ( !f || !f->sp ) {
                std::cerr << "Solve requested before factorization in oft_strumpack_sparse_c()" << std::endl;
                *info = -1;
                return;
            }
            if ( n != f->n ) {
                std::cerr << "Matrix size changed since factorization in oft_strumpack_sparse_c()" << std::endl;
                *info = -1;
                return;
            }
            /* Solve the system A*X=B, overwriting B with X (RHS copy required) */
            std::vector<double> rhs(b, b + (size_t)ldb*nrhs);
            ReturnCode ierr = f->sp->solve(nrhs, rhs.data(), ldb, b, ldb);
            if ( ierr != ReturnCode::SUCCESS ) {
                std::cerr << "STRUMPACK solve failed: " << ierr << std::endl;
                *info = 1;
            }
        } else if ( iopt == 4 ) { /* Free storage */
            delete *f_factors;
            *f_factors = nullptr;
        } else {
            std::cerr << "Invalid iopt=" << iopt << " passed to oft_strumpack_sparse_c()" << std::endl;
            *info = -1;
        }
    } catch (std::exception& e) {
        std::cerr << "STRUMPACK operation (iopt=" << iopt << ") failed: " << e.what() << std::endl;
        if ( *f_factors && iopt != 2 ) (*f_factors)->sp.reset();
        *info = 1;
    } catch (...) {
        std::cerr << "STRUMPACK operation (iopt=" << iopt << ") failed: unknown exception" << std::endl;
        if ( *f_factors && iopt != 2 ) (*f_factors)->sp.reset();
        *info = 1;
    }
}

#ifdef HAVE_STRUMPACK_MPI
/* Distributed factorization, the stored (0-based) local graph and the options it was computed with */
struct factors_dist_t {
    std::unique_ptr<SparseSolverMPIDist<double,int>> sp;
    MPI_Comm comm = MPI_COMM_NULL;
    int n = 0;
    std::vector<int> rowptr;
    std::vector<int> colind;
    std::vector<int> dist;
    options_t opts = {};
};

extern "C" void
oft_strumpack_dist_c(int iopt, int n, int nnz_loc, int m_loc,
                     double *values, int *colind, int *rowptr, int *dist,
                     double *b, int *f_comm,
                     factors_dist_t **f_factors,
                     int comp_type, double rel_tol, double abs_tol,
                     int leaf_size, int max_rank, int min_sep_size,
                     int verbose, int *info)
/*
* Distributed (row-block) factorization of the full system using STRUMPACK's
* SparseSolverMPIDist. This routine can be called from Fortran and must be
* called collectively by all processes in the communicator.
*
* iopt (input) int
*      Specifies the operation:
*      = 1, performs reordering and LU decomposition for the first time
*      = 2, performs solve using existing factorization
*      = 3, refactor matrix with the same non-zero pattern (reuses the existing
*           reordering and MC64 matching/scaling, a full factorization is
*           performed instead if the solver options have changed)
*      = 4, free all the storage in the end
*
* n (input) int, global number of rows
*
* nnz_loc, m_loc (input) int
*      Number of non-zeros and rows owned by this process
*
* values, colind, rowptr (input)
*      Locally-owned rows in CSR format with 1-based indices. Row pointers are
*      local, column indices are global.
*
* dist (input) int array of size nprocs+1
*      Row distribution, process p owns rows [dist[p],dist[p+1]) (0-based)
*
* b (input/output) double array of size m_loc
*      On input the locally-owned rows of the RHS, on output the solution (iopt = 2)
*
* f_comm (input) Fortran MPI communicator handle (used when iopt = 1)
*
* comp_type, rel_tol, abs_tol, leaf_size, max_rank, min_sep_size (input)
*      See oft_strumpack_sparse_c
*/
{
    const options_t o = {comp_type, leaf_size, max_rank, min_sep_size, verbose, rel_tol, abs_tol};
    *info = 0;
    try {
        if ( iopt == 1 || iopt == 3 ) { /* LU decomposition */
            if ( !*f_factors ) *f_factors = new factors_dist_t;
            factors_dist_t *f = *f_factors;
            if ( f->comm == MPI_COMM_NULL ) f->comm = MPI_Comm_f2c(*f_comm);
            int nprocs;
            MPI_Comm_size(f->comm, &nprocs);
            const bool full = ( iopt == 1 ) || !f->sp || !( f->opts == o );
            if ( iopt == 1 ) {
                /* Store 0-based copy of local graph and distribution */
                f->rowptr.resize(m_loc+1);
                for (int i = 0; i <= m_loc; ++i) f->rowptr[i] = rowptr[i] - 1;
                f->colind.resize(nnz_loc);
                for (int i = 0; i < nnz_loc; ++i) f->colind[i] = colind[i] - 1;
                f->dist.assign(dist, dist + nprocs + 1);
                f->n = n;
            } else if ( n != f->n || nnz_loc != (int)f->colind.size() || m_loc+1 != (int)f->rowptr.size() ) {
                std::cerr << "Matrix size changed on refactor in oft_strumpack_dist_c()" << std::endl;
                *info = -1;
                return;
            }
            if ( full ) {
                CompressionType comp;
                if ( !get_compression(o, comp) ) {
                    f->sp.reset();
                    *info = -1;
                    return;
                }
                f->sp.reset(new SparseSolverMPIDist<double,int>(f->comm, o.verbose != 0));
                apply_options(f->sp->options(), o, comp);
                f->opts = o;
                f->sp->set_distributed_csr_matrix(m_loc, f->rowptr.data(), f->colind.data(), values,
                                                  f->dist.data(), false);
            } else {
                f->sp->update_matrix_values(m_loc, f->rowptr.data(), f->colind.data(), values,
                                            f->dist.data(), false);
            }
            ReturnCode ierr = f->sp->factor();
            if ( ierr != ReturnCode::SUCCESS ) {
                std::cerr << "STRUMPACK factorization failed: " << ierr << std::endl;
                *info = 1;
            }
        } else if ( iopt == 2 ) { /* Solve */
            factors_dist_t *f = *f_factors;
            if ( !f || !f->sp ) {
                std::cerr << "Solve requested before factorization in oft_strumpack_dist_c()" << std::endl;
                *info = -1;
                return;
            }
            /* Solve the system A*X=B, overwriting B with X (RHS copy required) */
            std::vector<double> rhs(b, b + m_loc);
            ReturnCode ierr = f->sp->solve(rhs.data(), b);
            if ( ierr != ReturnCode::SUCCESS ) {
                std::cerr << "STRUMPACK solve failed: " << ierr << std::endl;
                *info = 1;
            }
        } else if ( iopt == 4 ) { /* Free storage */
            delete *f_factors;
            *f_factors = nullptr;
        } else {
            std::cerr << "Invalid iopt=" << iopt << " passed to oft_strumpack_dist_c()" << std::endl;
            *info = -1;
        }
    } catch (std::exception& e) {
        std::cerr << "STRUMPACK operation (iopt=" << iopt << ") failed: " << e.what() << std::endl;
        if ( *f_factors && iopt != 2 ) (*f_factors)->sp.reset();
        *info = 1;
    } catch (...) {
        std::cerr << "STRUMPACK operation (iopt=" << iopt << ") failed: unknown exception" << std::endl;
        if ( *f_factors && iopt != 2 ) (*f_factors)->sp.reset();
        *info = 1;
    }
}
#endif
#endif
