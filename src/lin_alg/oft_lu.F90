!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_lu.F90
!
!> Native module for direct solves.
!! Available interfaces:
!! - LAPACK (dgetrf, dgetrs)
!! - SuperLU
!! - MUMPS
!! - PARDISO (MKL interface)
!! - STRUMPACK
!!
!! @note Local solves are supported by all packages. Solves of the full distributed
!! system are supported by SuperLU-DIST, MUMPS, and STRUMPACK (if built with MPI support)
!! (see @ref oft_lusolver::distributed).
!!
!! @authors Chris Hansen
!! @date May 2014
!! @ingroup doxy_oft_lin_alg
!------------------------------------------------------------------------------
MODULE oft_lu
USE, INTRINSIC :: iso_c_binding, ONLY: c_bool, c_int, c_double, c_ptr, c_null_ptr, c_f_pointer, c_loc, &
  c_associated
USE oft_local
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_graph
USE oft_solver_base, ONLY: oft_solver
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
IMPLICIT NONE
#include "local.h"
#ifdef HAVE_MUMPS
#include "dmumps_struc.h"
#endif
#ifdef HAVE_MKL_PARDISO
#define HAVE_MKL
#endif
#ifdef HAVE_MKL
#define DEF_LU_PACK "mkl"
#define DEF_ILU_PACK "mkl"
#else
#define DEF_ILU_PACK "native"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_UMFPACK)
#define DEF_LU_PACK "umfpack"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_MUMPS)
#define DEF_LU_PACK "mumps"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU)
#define DEF_LU_PACK "super"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU_DIST)
#define DEF_LU_PACK "superd"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_STRUMPACK)
#define DEF_LU_PACK "strumpack"
#endif
#if !defined(DEF_LU_PACK)
#define DEF_LU_PACK "none"
#endif
PRIVATE
PUBLIC lapack_matinv, lapack_cholesky
!------------------------------------------------------------------------------
!> Native ILU(0) preconditioner information object
!------------------------------------------------------------------------------
TYPE :: native_ilu_struc
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: ju => NULL() !< Diagonal entry indices
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: jlu => NULL() !< Column values in iLU factors
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: alu => NULL() !< Values of iLU factors
END TYPE native_ilu_struc
!------------------------------------------------------------------------------
!> SuperLU solver information
!------------------------------------------------------------------------------
TYPE :: superlu_struc
  INTEGER(i4) :: col_perm = -1 !< Column permutation to use (-1 for default)
#ifdef OFT_MPI_F08
  TYPE(mpi_comm) :: comm = MPI_COMM_NULL !< Private MPI communicator for SuperLU-DIST
#else
  INTEGER(i4) :: comm = MPI_COMM_NULL !< Private MPI communicator for SuperLU-DIST
#endif
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< SuperLU factors pointer
  TYPE(C_PTR) :: grid_handle = C_NULL_PTR !< SuperLU-DIST processor grid pointer
  INTEGER(i4), POINTER, DIMENSION(:) :: kr => NULL() !< Row pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL() !< Col pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: csc_map => NULL() !< Col pointer (0-based indices)
  REAL(r8), POINTER, DIMENSION(:) :: csc_vals => NULL() !< Col pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc_orig => NULL() !< Unpermuted copy of `lc` (SuperLU-DIST only)
END TYPE superlu_struc
!------------------------------------------------------------------------------
!> STRUMPACK solver information
!!
!! Numerical options left at their default values (`-1`) use the STRUMPACK defaults.
!!
!! @note When compression is enabled it is applied to the frontal matrices, so the
!! factorization (and each solve) is approximate. This is intended for use as a
!! preconditioner (eg. block-Jacobi sub-solver). STRUMPACK's iterative refinement and
!! Krylov outer solvers are not supported, so `iter_refine` must be `.FALSE.`.
!!
!! @note A refactorization with the same non-zero pattern (`update(.FALSE.)`) reuses the
!! fill-reducing reordering and MC64 matching/scaling computed at the last full
!! factorization, for both local and distributed (`distributed`) solves. STRUMPACK only
!! pivots within each front, so after large changes in matrix values the solution can be
!! inaccurate unless a full factorization is requested with `update(.TRUE.)`, even though
!! the non-zero pattern has not changed. A full factorization is also performed
!! automatically if the non-zero pattern or any of the options below has changed since
!! the last factorization.
!------------------------------------------------------------------------------
TYPE :: strumpack_struc
  INTEGER(i4) :: comp_type = -1 !< Compression type (0 -> none, 1 -> HSS, 2 -> BLR, -1 -> default [none])
  INTEGER(i4) :: leaf_size = -1 !< Size of leaf/tile blocks for compression
  INTEGER(i4) :: max_rank = -1 !< Maximum rank of compressed blocks
  INTEGER(i4) :: min_sep_size = -1 !< Minimum separator size for compression
  LOGICAL :: verbose = .FALSE. !< Print STRUMPACK compression/factorization information
  REAL(r8) :: comp_rel_tol = -1.d0 !< Relative compression tolerance
  REAL(r8) :: comp_abs_tol = -1.d0 !< Absolute compression tolerance
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< STRUMPACK factorization pointer
  INTEGER(i4), POINTER, DIMENSION(:) :: kr => NULL() !< Row pointer (0-based indices)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL() !< Column indices (0-based indices)
END TYPE strumpack_struc
!------------------------------------------------------------------------------
!> Pardiso solver information
!------------------------------------------------------------------------------
TYPE :: pardiso_struc
  INTEGER(i4) :: mtype = 11 !< Matrix type (1 -> real/symmetric, 11 -> real/nonsymmetric)
  INTEGER(i4) :: msglvl = 0 !< Controls message output (0 -> no output, 1 -> print some stats)
  INTEGER(i8) :: pt(64) = 0 !< Handle to internal data structure
  INTEGER(i4) :: iparm(64) = 0 !< This array is used to pass various parameters to Intel MKL PARDISO
  INTEGER(i4), POINTER, DIMENSION(:) :: perm => NULL() !< Permutation vector
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< ILU0 factors pointer
END TYPE pardiso_struc
!------------------------------------------------------------------------------
!> Row-distributed matrix information for parallel direct solves
!!
!! Each process holds the full rows corresponding to its locally-owned vector
!! entries, which are numbered contiguously in rank order.
!------------------------------------------------------------------------------
TYPE :: lu_dist_struc
  LOGICAL :: active = .FALSE. !< Is a distributed solve being performed?
#ifdef OFT_MPI_F08
  TYPE(mpi_comm) :: comm = MPI_COMM_NULL !< Private MPI communicator for solver package
#else
  INTEGER(i4) :: comm = MPI_COMM_NULL !< Private MPI communicator for solver package
#endif
  INTEGER(i4) :: ng = 0 !< Global number of rows
  INTEGER(i4) :: nloc = 0 !< Number of locally-owned rows
  INTEGER(i4) :: offset = 0 !< Number of rows owned by lower ranks
  INTEGER(i4) :: nnz = 0 !< Number of non-zeros in locally-owned rows
  INTEGER(i4) :: nproc_con = 0 !< Number of neighboring processes
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: proc_con => NULL() !< Ranks of neighboring processes [nproc_con]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: counts => NULL() !< Number of rows owned by each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: displs => NULL() !< Number of rows owned by lower ranks for each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: rows => NULL() !< Global indices (1-based) of locally-owned rows [nloc]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: kr => NULL() !< Row pointer (1-based) [nloc+1]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: lc => NULL() !< Global column indices (1-based) [nnz]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: irn => NULL() !< Global row indices (1-based) [nnz]
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: vals => NULL() !< Matrix values [nnz]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: loc_src => NULL() !< Local matrix entries in locally-owned rows
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: loc_dst => NULL() !< Location of `loc_src` entries in `vals`
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: send_map => NULL() !< Local matrix entries to send to each neighbor
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: recv_map => NULL() !< Location in `vals` of entries received from each neighbor
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: src_kr => NULL() !< Row pointer of native matrix used in setup [nr+1]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: src_lc => NULL() !< Column indices of native matrix used in setup [nnz]
  !---Distributed solution returned by the solver package ("mumps" only)
  LOGICAL :: sol_map_valid = .FALSE. !< Is the solution redistribution map current?
  INTEGER(i4) :: nsol = 0 !< Number of solution entries returned on this process
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_loc => NULL() !< Solution entries returned on this process [nsol]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: isol_loc => NULL() !< Global indices (1-based) of `sol_loc` entries [nsol]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_perm => NULL() !< `sol_loc` entries ordered by destination process [nsol]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_scounts => NULL() !< Solution entries sent to each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_sdispls => NULL() !< Offsets of entries sent to each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_rcounts => NULL() !< Solution entries received from each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_rdispls => NULL() !< Offsets of entries received from each process [nprocs]
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: sol_ridx => NULL() !< Locally-owned row of each received entry [nloc]
END TYPE lu_dist_struc
!------------------------------------------------------------------------------
!> LU solver class
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_lusolver
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  LOGICAL :: update_graph = .TRUE. !< Perform full factorization including symbolic steps
  LOGICAL(kind=c_bool) :: iter_refine = .FALSE. !< Perform iterative refinement
  !> Factor the full distributed system across all processes ("superd", "mumps", and "strumpack" only).
  !! Ignored when the system is not distributed (eg. serial runs or local blocks). Changes to the
  !! matrix's sparsity pattern are detected automatically and trigger a full factorization.
  LOGICAL :: distributed = .FALSE.
  INTEGER(i4) :: nrhs = 1 !< Number of right hand sides
  INTEGER(i4), POINTER, DIMENSION(:) :: ipiv => NULL() !< Row pivot array ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: atmp => NULL() !< Local dense matrix ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: sec_rhs => NULL() !< Storage for additional RHS data
  CHARACTER(LEN=9) :: package = DEF_LU_PACK !< Factorization package
#if defined( HAVE_SUPERLU ) || defined( HAVE_SUPERLU_DIST ) || defined( HAVE_UMFPACK )
  TYPE(superlu_struc) :: superlu_struct !< SuperLU solver information
#endif
#ifdef HAVE_STRUMPACK
  TYPE(strumpack_struc) :: strumpack_struct !< STRUMPACK solver information
#endif
#ifdef HAVE_MKL
  TYPE(pardiso_struc) :: pardiso_struct !< MKL/PARDISO solver information
#endif
#ifdef HAVE_MUMPS
  TYPE(DMUMPS_STRUC) :: mumps_struct !< MUMPS solver information
#endif
  TYPE(lu_dist_struc) :: dist !< Distributed matrix information
CONTAINS
  !> Solve linear system
  PROCEDURE :: apply => lusolver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => lusolver_setup_xml
  !> Update solver with new settings/operators and refactor
  PROCEDURE :: update => lusolver_update
  !> Check thread safety
  PROCEDURE :: check_thread => lusolver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => lusolver_delete
END TYPE oft_lusolver
!------------------------------------------------------------------------------
!> ILU solver class
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_ilusolver
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  LOGICAL :: update_graph = .TRUE. !< Perform full factorization including symbolic steps
  INTEGER(i4) :: nrhs = 1 !< Number of right hand sides
  REAL(r8), POINTER, DIMENSION(:,:) :: sec_rhs => NULL() !< Storage for additional RHS data
  CHARACTER(LEN=7) :: package = DEF_ILU_PACK !< Factorization package
  TYPE(native_ilu_struc) :: native_struct !< Data for native iLU solver
#if defined( HAVE_SUPERLU )
  TYPE(superlu_struc) :: superlu_struct !< SuperLU solver information
#endif
#ifdef HAVE_MKL
  TYPE(pardiso_struc) :: pardiso_struct !< MKL/PARDISO solver information
#endif
CONTAINS
  !> Solve linear system
  PROCEDURE :: apply => ilusolver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => ilusolver_setup_xml
  !> Update solver with new settings/operators and refactor
  PROCEDURE :: update => ilusolver_update
  !> Check thread safety
  PROCEDURE :: check_thread => ilusolver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => ilusolver_delete
END TYPE oft_ilusolver
INTERFACE
#ifdef HAVE_SUPERLU
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_superlu_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dgssv
!------------------------------------------------------------------------------
!> Interface to dgsisx (iLUT) from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dgsisx(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,fill_tol,info) BIND(C,NAME="oft_superlu_dgsisx_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  REAL(c_double), VALUE, INTENT(in) :: fill_tol !< Fill tolerance for tresholded iLU
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dgsisx
#endif
#ifdef HAVE_SUPERLU_DIST
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU-DIST (local only)
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dist_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,grid_handle,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_superlu_dist_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), VALUE, INTENT(in) :: grid_handle !< Pointer to SuperLU-DIST communication object
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU-DIST internal data storage
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dist_dgssv
!------------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dist_slugrid(iopt,slu_comm,nprow,npcol,grid_handle,info) BIND(C,NAME="oft_superlu_dist_slugrid_c")
  IMPORT c_int, c_double, c_ptr
#ifdef OFT_MPI_F08
  IMPORT mpi_comm
#endif
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> Setup grid, 2 -> Clear/Release grid)
#ifdef OFT_MPI_F08
  TYPE(mpi_comm), INTENT(in) :: slu_comm !< MPI communicator for processor grid (must be MPI_COMM_SELF)
#else
  INTEGER(c_int), INTENT(in) :: slu_comm
#endif
  INTEGER(c_int), VALUE, INTENT(in) :: nprow !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: npcol !< Number of non-zeros
  TYPE(c_ptr), INTENT(inout) :: grid_handle !< Pointer to SuperLU-DIST processor grid
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dist_slugrid
!------------------------------------------------------------------------------
!> Interface to pdgssvx from SuperLU-DIST (row-distributed matrix)
!------------------------------------------------------------------------------
  SUBROUTINE oft_superlu_dist_pdgssvx(iopt,n,nnz_loc,m_loc,fst_row,nrhs,values,colind, &
    rowptr,b,ldb,grid_handle,f_factors,iter_refine,info) BIND(C,NAME="oft_superlu_dist_pdgssvx_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition with
  !! a new non-zero pattern, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Global number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz_loc !< Number of non-zeros in local rows
  INTEGER(c_int), VALUE, INTENT(in) :: m_loc !< Number of local rows
  INTEGER(c_int), VALUE, INTENT(in) :: fst_row !< Global index of first local row (0-based)
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve
  REAL(c_double), DIMENSION(*), INTENT(in) :: values !< Local matrix values in CRS format [nnz_loc]
  INTEGER(c_int), DIMENSION(*), INTENT(in) :: colind !< Global column indices (1-based) in CRS format [nnz_loc]
  INTEGER(c_int), DIMENSION(*), INTENT(in) :: rowptr !< Row pointer (1-based) into colind [m_loc+1]
  REAL(c_double), DIMENSION(*), INTENT(inout) :: b !< Local rows of RHS -> overwritten with solution [ldb,nrhs]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Leading dimension of `b` (>= m_loc)
  TYPE(c_ptr), VALUE, INTENT(in) :: grid_handle !< Pointer to SuperLU-DIST processor grid
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU-DIST internal data storage
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), INTENT(out) :: info !< Solver return status
  END SUBROUTINE oft_superlu_dist_pdgssvx
#endif
#ifdef HAVE_STRUMPACK
!------------------------------------------------------------------------------
!> Interface to sparse (multifrontal) LU factorization from STRUMPACK
!------------------------------------------------------------------------------
  SUBROUTINE oft_strumpack_sparse(iopt,n,nnz,nrhs,values,colind,rowptr,b,ldb, &
    f_factors,comp_type,rel_tol,abs_tol,leaf_size,max_rank,min_sep_size, &
    verbose,info) BIND(C,NAME="oft_strumpack_sparse_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs solve, 3 -> refactor matrix with the
  !! same non-zero pattern [full factorization if options changed], 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format (0-based) [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind (0-based) [n+1]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  REAL(c_double), DIMENSION(ldb,nrhs), INTENT(inout) :: b !< Right hand side -> overwritten with solution [ldb,nrhs]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to STRUMPACK internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: comp_type !< Compression type (0 -> none, 1 -> HSS, 2 -> BLR, -1 -> default)
  REAL(c_double), VALUE, INTENT(in) :: rel_tol !< Relative compression tolerance (<=0 -> default)
  REAL(c_double), VALUE, INTENT(in) :: abs_tol !< Absolute compression tolerance (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: leaf_size !< Size of leaf/tile blocks (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: max_rank !< Maximum rank of compressed blocks (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: min_sep_size !< Minimum separator size for compression (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: verbose !< Print STRUMPACK information (0 -> no, 1 -> yes)
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (0 -> success)
  END SUBROUTINE oft_strumpack_sparse
#ifdef HAVE_STRUMPACK_MPI
!------------------------------------------------------------------------------
!> Interface to distributed sparse (multifrontal) LU factorization from STRUMPACK
!!
!! Must be called collectively by all processes, with locally-owned rows in CRS format
!------------------------------------------------------------------------------
  SUBROUTINE oft_strumpack_dist(iopt,n,nnz_loc,m_loc,values,colind,rowptr,dist,b, &
    comm,f_factors,comp_type,rel_tol,abs_tol,leaf_size,max_rank,min_sep_size, &
    verbose,info) BIND(C,NAME="oft_strumpack_dist_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs solve, 3 -> refactor matrix with the
  !! same non-zero pattern [full factorization if options changed], 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Global number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz_loc !< Number of non-zeros in locally-owned rows
  INTEGER(c_int), VALUE, INTENT(in) :: m_loc !< Number of locally-owned rows
  REAL(c_double), DIMENSION(*), INTENT(in) :: values !< Matrix values in CRS format [nnz_loc]
  INTEGER(c_int), DIMENSION(*), INTENT(in) :: colind !< Global column indices in CRS format (1-based) [nnz_loc]
  INTEGER(c_int), DIMENSION(*), INTENT(in) :: rowptr !< Local row pointer into colind (1-based) [m_loc+1]
  INTEGER(c_int), DIMENSION(*), INTENT(in) :: dist !< First global row (0-based) on each process [nprocs+1]
  REAL(c_double), DIMENSION(*), INTENT(inout) :: b !< Right hand side -> overwritten with solution [m_loc]
  INTEGER(c_int), INTENT(in) :: comm !< Fortran MPI communicator handle
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to STRUMPACK internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: comp_type !< Compression type (0 -> none, 1 -> HSS, 2 -> BLR, -1 -> default)
  REAL(c_double), VALUE, INTENT(in) :: rel_tol !< Relative compression tolerance (<=0 -> default)
  REAL(c_double), VALUE, INTENT(in) :: abs_tol !< Absolute compression tolerance (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: leaf_size !< Size of leaf/tile blocks (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: max_rank !< Maximum rank of compressed blocks (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: min_sep_size !< Minimum separator size for compression (<=0 -> default)
  INTEGER(c_int), VALUE, INTENT(in) :: verbose !< Print STRUMPACK information (0 -> no, 1 -> yes)
  INTEGER(c_int), INTENT(out) :: info !< Solver return status (0 -> success)
  END SUBROUTINE oft_strumpack_dist
#endif
#endif
!------------------------------------------------------------------------------
!> Interface to dgssv from UMFPACK
!------------------------------------------------------------------------------
#ifdef HAVE_UMFPACK
  SUBROUTINE oft_umfpack_dgssv(iopt,n,nnz,nrhs,values,colind, &
    rowptr,b,ldb,f_factors,col_perm,iter_refine,info) BIND(C,NAME="oft_umfpack_dgssv_c")
  IMPORT c_bool, c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !<  Operation control (1 -> performs LU decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> refactor matrix with the
  !! same non-zero pattern, 4 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  INTEGER(c_int), VALUE, INTENT(in) :: nrhs !< Number of RHS to solve (nrhs=1 only)
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  INTEGER(c_int), VALUE, INTENT(in) :: ldb !< Lowest dimension of vector `b` [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to UMFPACK internal data storage
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_umfpack_dgssv
#endif
#ifdef HAVE_MKL
!------------------------------------------------------------------------------
!> Interface to MKL iLU0 solver
!------------------------------------------------------------------------------
  SUBROUTINE oft_mkl_ilu(iopt,n,nnz,values,colind, &
    rowptr,b,f_factors,info) BIND(C,NAME="oft_mkl_ilu_c")
  IMPORT c_int, c_double, c_ptr
  INTEGER(c_int), VALUE, INTENT(in) :: iopt !< Operation control (1 -> performs iLU0 decomposition for the
  !! first time, 2 -> performs triangular solve, 3 -> free all the storage in the end)
  INTEGER(c_int), VALUE, INTENT(in) :: n !< Number of rows
  INTEGER(c_int), VALUE, INTENT(in) :: nnz !< Number of non-zeros
  REAL(c_double), DIMENSION(nnz), INTENT(in) :: values !< Matrix values in CRS format [nnz]
  INTEGER(c_int), DIMENSION(nnz), INTENT(in) :: colind !< Column indices in CRS format [nnz]
  INTEGER(c_int), DIMENSION(n+1), INTENT(in) :: rowptr !< Row pointer into colind [n+1]
  REAL(c_double), DIMENSION(n), INTENT(inout) :: b !< Right hand side -> overwritten with solution [n]
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to MKL internal data storage
  INTEGER(c_int), INTENT(inout) :: info !< Pointer to MKL internal data storage
  END SUBROUTINE oft_mkl_ilu
#endif
END INTERFACE
INTERFACE lapack_matinv
  MODULE PROCEDURE lapack_matinv_real
  MODULE PROCEDURE lapack_matinv_complex
  MODULE PROCEDURE lapack_matinv_complex4
END INTERFACE lapack_matinv
INTERFACE lapack_cholesky
  MODULE PROCEDURE lapack_cholesky_real
  MODULE PROCEDURE lapack_cholesky_real4
  ! MODULE PROCEDURE lapack_cholesky_complex
END INTERFACE lapack_cholesky
CONTAINS
!------------------------------------------------------------------------------
!> Convert csr graph to csc storage
!------------------------------------------------------------------------------
SUBROUTINE convert_csr_to_csc(csr_graph,csc_graph,map)
TYPE(oft_graph), INTENT(inout) :: csr_graph !< Original graph in csr format
TYPE(oft_graph), INTENT(inout) :: csc_graph !< Resulting graph in csc format
INTEGER(i4), INTENT(inout) :: map(:) !< Mapping from csr entries to csc entries
INTEGER(i4) :: i,j,k
INTEGER(i4), POINTER :: kr_tmp(:)
DEBUG_STACK_PUSH
!---Convert to CSC format
csc_graph%nr=csr_graph%nr; csc_graph%nc=csr_graph%nc
csc_graph%nrg=-1; csc_graph%ncg=-1
csc_graph%nnz=csr_graph%nnz
!---
ALLOCATE(csc_graph%kr(csc_graph%nr+1))
csc_graph%kr=0
DO i=1,csr_graph%nr
  DO j=csr_graph%kr(i),csr_graph%kr(i+1)-1
    csc_graph%kr(csr_graph%lc(j))=csc_graph%kr(csr_graph%lc(j))+1
  END DO
END DO
!---
csc_graph%kr(csc_graph%nr+1)=csc_graph%nnz+1
DO i=csc_graph%nr,1,-1 ! cumulative point to point count
  csc_graph%kr(i)=csc_graph%kr(i+1)-csc_graph%kr(i)
END DO
IF(csc_graph%kr(1)/=1)CALL oft_abort('Bad element to element count', &
'convert_csr_to_csc',__FILE__)
!---
ALLOCATE(kr_tmp(csc_graph%nr))
ALLOCATE(csc_graph%lc(csc_graph%nnz))
kr_tmp=0
DO i=1,csr_graph%nr
  DO j=csr_graph%kr(i),csr_graph%kr(i+1)-1
    k=csc_graph%kr(csr_graph%lc(j))+kr_tmp(csr_graph%lc(j))
    csc_graph%lc(k)=i
    map(k)=j!csr_map(j)
    kr_tmp(csr_graph%lc(j))=kr_tmp(csr_graph%lc(j))+1
  END DO
END DO
DEALLOCATE(kr_tmp)
DEBUG_STACK_POP
END SUBROUTINE convert_csr_to_csc
!------------------------------------------------------------------------------
!> Solve a linear system using a direct solve
!!
!! @note Do not call with `superd`/`mumps` inside an OpenMP parallel region, even an inactive one (can deadlock on macOS)
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE lusolver_apply(self,u,g)
CLASS(oft_lusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: mode,nrhs,ldb,ierr,i,j,k,info
INTEGER(i4), POINTER :: csr_map(:),kr_tmp(:)
LOGICAL :: graph_changed
REAL(r8), POINTER, DIMENSION(:) :: mat_vals,csc_vals,vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vals,b
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: mytimer
TYPE(oft_graph) :: csr_graph
TYPE(oft_graph) :: csc_graph
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting LU solver: ',self%package,self%refactor
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'lusolver_apply',__FILE__)
!---Check call threading for MUMPS
IF(self%package(1:5)=='mumps')THEN
  IF(omp_get_num_threads()>1)CALL oft_abort('MUMPS used by multiple threads', &
    'lusolver_apply',__FILE__)
END IF
!---Check for distributed solve
IF(.NOT.self%initialized)THEN
  self%dist%active=.FALSE.
  IF(self%distributed.AND.(oft_env%nprocs>1).AND.ASSOCIATED(g%stitch_info))THEN
    self%dist%active=.NOT.(g%stitch_info%skip.OR.g%stitch_info%full)
  END IF
END IF
IF(self%dist%active)THEN
  CALL lusolver_apply_dist(self,A_native,u,g)
  DEBUG_STACK_POP
  RETURN
END IF
IF(ASSOCIATED(A_native%Mfull))THEN
  mat_vals=>A_native%Mfull
ELSE
  mat_vals=>A_native%M
END IF
!---Initialize solvers
IF(.NOT.self%initialized)THEN
  !---
  SELECT CASE(TRIM(self%package))
    CASE DEFAULT
      CALL oft_abort('Unknown LU package','lusolver_apply', __FILE__)
    CASE("lapack")
      ! Do nothing
    CASE("super")
#ifdef HAVE_SUPERLU
      ALLOCATE(self%superlu_struct%kr(A_native%nr+1))
      ALLOCATE(self%superlu_struct%lc(A_native%nnz))
      self%superlu_struct%kr=A_native%kr-1
      self%superlu_struct%lc=A_native%lc-1
#else
      CALL oft_abort('OFT not compiled with SUPERLU','lusolver_apply',__FILE__)
#endif
    CASE("superd")
#ifdef HAVE_SUPERLU_DIST
      ! CSC graph is built at the first factorization
      mode=1
      nrhs=1
      CALL MPI_Comm_dup(MPI_COMM_SELF,self%superlu_struct%comm,ierr)
      CALL oft_superlu_dist_slugrid(mode,self%superlu_struct%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
#else
      CALL oft_abort('OFT not compiled with SUPERLU-DIST','lusolver_apply',__FILE__)
#endif
    CASE("umfpack")
#ifdef HAVE_UMFPACK
      ALLOCATE(self%superlu_struct%kr(A_native%nr+1))
      ALLOCATE(self%superlu_struct%lc(A_native%nnz))
      self%superlu_struct%kr=A_native%kr-1
      self%superlu_struct%lc=A_native%lc-1
#else
      CALL oft_abort('OFT not compiled with UMFPACK','lusolver_apply',__FILE__)
#endif
    CASE("strumpack")
#if !defined(HAVE_STRUMPACK)
      CALL oft_abort('OFT not compiled with STRUMPACK','lusolver_apply',__FILE__)
#endif
    CASE("mkl")
#ifdef HAVE_MKL
      self%pardiso_struct%pt=0
      self%pardiso_struct%iparm=0
      call pardisoinit(self%pardiso_struct%pt,self%pardiso_struct%mtype,self%pardiso_struct%iparm)
      ALLOCATE(self%pardiso_struct%perm(A_native%nr))
      self%pardiso_struct%iparm(1)=1
      self%pardiso_struct%iparm(27)=1
#else
      CALL oft_abort('OFT not compiled with MKL-PARDISO','lusolver_apply',__FILE__)
#endif
    CASE("mumps")
#ifdef HAVE_MUMPS
      self%mumps_struct%job = -1
#ifdef OFT_MPI_F08
      self%mumps_struct%comm = MPI_COMM_SELF%MPI_VAL
#else
      self%mumps_struct%comm = MPI_COMM_SELF
#endif
      self%mumps_struct%sym = 0
      self%mumps_struct%par = 1
      CALL dmumps(self%mumps_struct)
      self%mumps_struct%icntl((/5,18/))=0 ! Set centralized assembly
      self%mumps_struct%icntl((/20,21/))=0 ! Set centralized dense RHS
      self%mumps_struct%icntl(3:4)=-1 ! Surpress information messages
      self%mumps_struct%n=A_native%nr
      self%mumps_struct%nz=A_native%nnz
      self%mumps_struct%jcn=>A_native%lc
      self%mumps_struct%A=>mat_vals
      ALLOCATE(self%mumps_struct%irn(A_native%nnz))
      DO i=1,A_native%nr
        self%mumps_struct%irn(A_native%kr(i):A_native%kr(i+1)-1)=i
      END DO
#else
      CALL oft_abort('OFT not compiled with MUMPS','lusolver_apply',__FILE__)
#endif
  END SELECT
  self%initialized=.TRUE.
END IF
!---Call factorization
ALLOCATE(vals(u%n,self%nrhs))
nrhs=1
vtmp=>vals(:,1)
CALL g%get_local(vtmp)
IF(self%nrhs>1)THEN
  nrhs=self%nrhs
  DO i=2,self%nrhs
    vals(:,i) = self%sec_rhs(:,i-1)
  END DO
END IF
ldb=A_native%nr
! NULLIFY(vals)
! CALL g%get_local(vals)
SELECT CASE(TRIM(self%package))
  CASE("super")
#ifdef HAVE_SUPERLU
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)mode=1
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
      !$omp critical (superlu_solve)
#endif
      CALL oft_superlu_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (SuperLU)','lusolver_apply',__FILE__)
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
      !$omp end critical (superlu_solve)
#endif
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
    !$omp critical (superlu_solve)
#endif
    CALL oft_superlu_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (SuperLU)','lusolver_apply',__FILE__)
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
    !$omp end critical (superlu_solve)
#endif
#else
    CALL oft_abort('OFT not compiled with SuperLU','lusolver_apply',__FILE__)
#endif
CASE("superd")
#ifdef HAVE_SUPERLU_DIST
    IF(self%refactor)THEN
      IF(self%update_graph.OR.(.NOT.ASSOCIATED(self%superlu_struct%lc_orig)))THEN
        !---Build CSC graph and perform full factorization
        IF(ASSOCIATED(self%superlu_struct%lc_orig))THEN
          DEALLOCATE(self%superlu_struct%csc_map,self%superlu_struct%csc_vals)
          DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc,self%superlu_struct%lc_orig)
        END IF
        csr_graph%nr=A_native%nr; csr_graph%nc=A_native%nc; csr_graph%nnz=A_native%nnz
        csr_graph%kr=>A_native%kr; csr_graph%lc=>A_native%lc
        ALLOCATE(self%superlu_struct%csc_map(A_native%nnz))
        ALLOCATE(self%superlu_struct%csc_vals(A_native%nnz))
        CALL convert_csr_to_csc(csr_graph,csc_graph,self%superlu_struct%csc_map)
        csc_graph%kr=csc_graph%kr-1; csc_graph%lc=csc_graph%lc-1
        self%superlu_struct%kr=>csc_graph%kr; self%superlu_struct%lc=>csc_graph%lc
        ALLOCATE(self%superlu_struct%lc_orig(A_native%nnz))
        self%superlu_struct%lc_orig=self%superlu_struct%lc
        mode=1
      ELSE
        !---Restore row indices, which the previous factorization (pdgssvx_ABglobal)
        !---overwrote with their row-permuted values, and refactor with the same pattern
        self%superlu_struct%lc=self%superlu_struct%lc_orig
        mode=3
      END IF
      self%superlu_struct%csc_vals=mat_vals(self%superlu_struct%csc_map)
      !!$omp critical (superlu_solve)
      CALL oft_superlu_dist_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        self%superlu_struct%csc_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%grid_handle,self%superlu_struct%f_factors, &
        self%superlu_struct%col_perm,self%iter_refine,ierr)
      !!$omp end critical (superlu_solve)
      IF(ierr/=0)CALL oft_abort('Factorization failed (SuperLU-DIST)','lusolver_apply',__FILE__)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    !!$omp critical (superlu_solve)
    CALL oft_superlu_dist_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      self%superlu_struct%csc_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%grid_handle,self%superlu_struct%f_factors, &
      self%superlu_struct%col_perm,self%iter_refine,ierr)
    !!$omp end critical (superlu_solve)
    IF(ierr/=0)CALL oft_abort('Solve failed (SuperLU-DIST)','lusolver_apply',__FILE__)
#else
    CALL oft_abort('OFT not compiled with SuperLU-DIST','lusolver_apply',__FILE__)
#endif
  CASE("umfpack")
#ifdef HAVE_UMFPACK
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)mode=1
      !!$omp critical (umfpack_solve)
      CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (UMFPACK)','lusolver_apply',__FILE__)
      !!$omp end critical (umfpack_solve)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    !!$omp critical (umfpack_solve)
    CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (UMFPACK)','lusolver_apply',__FILE__)
    !!$omp end critical (umfpack_solve)
#else
    CALL oft_abort('OFT not compiled with UMFPACK','lusolver_apply',__FILE__)
#endif
  CASE("strumpack")
#ifdef HAVE_STRUMPACK
    IF(self%iter_refine)CALL oft_abort('Iterative refinement not supported (STRUMPACK)','lusolver_apply',__FILE__)
    IF(self%refactor)THEN
      !---Full factorization if requested or the non-zero pattern has changed
      CALL strumpack_set_graph(self%strumpack_struct,A_native,graph_changed)
      mode=3
      IF(self%update_graph.OR.graph_changed)mode=1
      CALL strumpack_sparse_call(self%strumpack_struct,mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,vals,ldb,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (STRUMPACK)','lusolver_apply',__FILE__)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    CALL strumpack_sparse_call(self%strumpack_struct,mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,vals,ldb,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (STRUMPACK)','lusolver_apply',__FILE__)
#else
    CALL oft_abort('OFT not compiled with STRUMPACK','lusolver_apply',__FILE__)
#endif
  CASE("mkl")
#ifdef HAVE_MKL
    ALLOCATE(b(A_native%nr,self%nrhs))
    b=vals
    IF(self%refactor)THEN
      !---Update matrix
      mode=22
      IF(self%update_graph)mode=12
      CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,A_native%nr, &
        mat_vals,A_native%kr,A_native%lc,self%pardiso_struct%perm, &
        nrhs,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,b,vals,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (MKL-PARDISO)','lusolver_apply',__FILE__)
      self%pardiso_struct%iparm(1)=0
      self%pardiso_struct%iparm(27)=0
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=33
    IF(self%iter_refine)THEN
      self%pardiso_struct%iparm(8)=0
    ELSE
      self%pardiso_struct%iparm(8)=1
    END IF
    CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,A_native%nr, &
      mat_vals,A_native%kr,A_native%lc,self%pardiso_struct%perm, &
      nrhs,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,b,vals,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (MKL-PARDISO)','lusolver_apply',__FILE__)
    DEALLOCATE(b)
#else
    CALL oft_abort('OFT not compiled with MKL-PARDISO','lusolver_apply',__FILE__)
#endif
  CASE("mumps")
#ifdef HAVE_MUMPS
    IF(self%refactor)THEN
      ! (1 for analysis, 2 for factorization, 4 for both)
      self%mumps_struct%job = 2
      IF(self%update_graph)self%mumps_struct%job = 4
      CALL dmumps(self%mumps_struct)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    CALL c_f_pointer(C_LOC(vals),self%mumps_struct%rhs,[nrhs*g%n])
    ! self%mumps_struct%rhs => vals(1:g%n,1:nrhs)
    self%mumps_struct%lrhs = g%n
    self%mumps_struct%nrhs = nrhs
    !---
    self%mumps_struct%job = 3 ! Solve phase
    CALL dmumps(self%mumps_struct)
#else
    CALL oft_abort('OFT not compiled with MUMPS','lusolver_apply',__FILE__)
#endif
  CASE("lapack")
    IF(.NOT.ASSOCIATED(self%atmp))THEN
      ALLOCATE(self%ipiv(A_native%nr))
      ALLOCATE(self%atmp(A_native%nr,A_native%nc))
    END IF
    IF(self%refactor)THEN
      !---Update matrix
      self%atmp=0.d0
      DO i=1,A_native%nr
        DO j=A_native%kr(i),A_native%kr(i+1)-1
          self%atmp(i,A_native%lc(j))=mat_vals(j)
        END DO
      END DO
      CALL dgetrf(A_native%nr,A_native%nr,self%atmp,A_native%nr,self%ipiv,info)
      IF(info/=0)CALL oft_abort('Factorization failed (LAPACK)','lusolver_apply',__FILE__)
      self%refactor=.FALSE.
    END IF
    CALL dgetrs('N',A_native%nr,nrhs,self%atmp,A_native%nr,self%ipiv,vals,ldb,info)
    IF(info/=0)CALL oft_abort('Solve failed (LAPACK)','lusolver_apply',__FILE__)
  CASE DEFAULT
    CALL oft_abort('Unknown factorization package','lusolver_apply',__FILE__)
END SELECT
CALL u%restore_local(vals(:,1))
IF(self%nrhs>1)THEN
  DO i=2,self%nrhs
    self%sec_rhs(:,i-1) = vals(:,i)
  END DO
END IF
DEALLOCATE(vals)
DEBUG_STACK_POP
END SUBROUTINE lusolver_apply
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine lusolver_update(self,new_pattern)
class(oft_lusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))self%update_graph=new_pattern
end subroutine lusolver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition
!------------------------------------------------------------------------------
subroutine lusolver_setup_xml(self,solver_node,level)
CLASS(oft_lusolver), INTENT(inout) :: self
TYPE(xml_node), INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
!---
CHARACTER(LEN=:), ALLOCATABLE :: factor_package
CHARACTER(LEN=3) :: fac_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---
CALL xml_get_element(solver_node,"package",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,factor_package,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `package` node","lusolver_setup_xml",__FILE__)
  IF(ALLOCATED(factor_package))THEN
    IF(LEN(factor_package)>LEN(self%package))CALL oft_abort('Factorization package name too long','lusolver_setup_xml',__FILE__)
    self%package=factor_package
  END IF
END IF
CALL xml_get_element(solver_node,"distributed",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%distributed,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `distributed` node","lusolver_setup_xml",__FILE__)
END IF
#ifdef HAVE_STRUMPACK
!---STRUMPACK options
CALL xml_get_element(solver_node,"strumpack",current_node,ierr)
IF(ierr==0)CALL strumpack_setup_xml(self%strumpack_struct,current_node)
#endif
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'LU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:  ',self%package
  WRITE(*,'(2X,A,L1)')'- Distributed:  ',self%distributed
#ifdef HAVE_STRUMPACK
  IF(TRIM(self%package)=='strumpack')THEN
    WRITE(*,'(2X,A,I4)')     '- Comp type:   ',self%strumpack_struct%comp_type
    WRITE(*,'(2X,A,ES11.3)') '- Comp rtol:   ',self%strumpack_struct%comp_rel_tol
    WRITE(*,'(2X,A,ES11.3)') '- Comp atol:   ',self%strumpack_struct%comp_abs_tol
    WRITE(*,'(2X,A,I6)')     '- Leaf size:   ',self%strumpack_struct%leaf_size
    WRITE(*,'(2X,A,I6)')     '- Max rank:    ',self%strumpack_struct%max_rank
    WRITE(*,'(2X,A,I6)')     '- Min sep size:',self%strumpack_struct%min_sep_size
  END IF
#endif
END IF
DEBUG_STACK_POP
end subroutine lusolver_setup_xml
#ifdef HAVE_STRUMPACK
!------------------------------------------------------------------------------
!> Read STRUMPACK solver options from XML definition
!!
!! Supported elements: `comp_type` {none, hss, blr}, `comp_rel_tol`, `comp_abs_tol`,
!! `leaf_size`, `max_rank`, `min_sep_size`, `verbose`
!------------------------------------------------------------------------------
subroutine strumpack_setup_xml(self,strumpack_node)
TYPE(strumpack_struc), INTENT(inout) :: self
TYPE(xml_node), INTENT(in) :: strumpack_node !< XML node containing STRUMPACK options
!---
TYPE(xml_node) :: current_node
CHARACTER(LEN=:), ALLOCATABLE :: comp_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
CALL xml_get_element(strumpack_node,"comp_type",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,comp_type,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `comp_type` node","strumpack_setup_xml",__FILE__)
  SELECT CASE(TRIM(comp_type))
    CASE("none")
      self%comp_type=0
    CASE("hss")
      self%comp_type=1
    CASE("blr")
      self%comp_type=2
    CASE DEFAULT
      CALL oft_abort('Unknown compression type, must be one of {none, hss, blr}','strumpack_setup_xml',__FILE__)
  END SELECT
END IF
CALL xml_get_element(strumpack_node,"comp_rel_tol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%comp_rel_tol,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `comp_rel_tol` node","strumpack_setup_xml",__FILE__)
END IF
CALL xml_get_element(strumpack_node,"comp_abs_tol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%comp_abs_tol,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `comp_abs_tol` node","strumpack_setup_xml",__FILE__)
END IF
CALL xml_get_element(strumpack_node,"leaf_size",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%leaf_size,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `leaf_size` node","strumpack_setup_xml",__FILE__)
END IF
CALL xml_get_element(strumpack_node,"max_rank",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%max_rank,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `max_rank` node","strumpack_setup_xml",__FILE__)
END IF
CALL xml_get_element(strumpack_node,"min_sep_size",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%min_sep_size,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `min_sep_size` node","strumpack_setup_xml",__FILE__)
END IF
CALL xml_get_element(strumpack_node,"verbose",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%verbose,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `verbose` node","strumpack_setup_xml",__FILE__)
END IF
DEBUG_STACK_POP
end subroutine strumpack_setup_xml
!------------------------------------------------------------------------------
!> Store a 0-based copy of the matrix graph for STRUMPACK
!------------------------------------------------------------------------------
SUBROUTINE strumpack_set_graph(self,A,changed)
TYPE(strumpack_struc), INTENT(inout) :: self
CLASS(oft_native_matrix), INTENT(in) :: A !< Matrix
LOGICAL, INTENT(out) :: changed !< Graph differs from the stored copy?
DEBUG_STACK_PUSH
changed=.TRUE.
IF(ASSOCIATED(self%kr))THEN
  IF((SIZE(self%kr)==A%nr+1).AND.(SIZE(self%lc)==A%nnz))THEN
    changed=ANY(self%kr/=A%kr(1:A%nr+1)-1)
    IF(.NOT.changed)changed=ANY(self%lc/=A%lc(1:A%nnz)-1)
  ELSE
    DEALLOCATE(self%kr,self%lc)
  END IF
END IF
IF(changed)THEN
  IF(.NOT.ASSOCIATED(self%kr))ALLOCATE(self%kr(A%nr+1),self%lc(A%nnz))
  self%kr=A%kr(1:A%nr+1)-1
  self%lc=A%lc(1:A%nnz)-1
END IF
DEBUG_STACK_POP
END SUBROUTINE strumpack_set_graph
!------------------------------------------------------------------------------
!> Call STRUMPACK sparse solver interface using stored graph and options
!------------------------------------------------------------------------------
SUBROUTINE strumpack_sparse_call(self,mode,n,nnz,nrhs,values,b,ldb,info)
TYPE(strumpack_struc), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: mode !< Operation (see @ref oft_strumpack_sparse)
INTEGER(i4), INTENT(in) :: n !< Number of rows
INTEGER(i4), INTENT(in) :: nnz !< Number of non-zeros
INTEGER(i4), INTENT(in) :: nrhs !< Number of RHS
REAL(r8), INTENT(in) :: values(nnz) !< Matrix values in CRS format [nnz]
INTEGER(i4), INTENT(in) :: ldb !< Lowest dimension of `b`
REAL(r8), INTENT(inout) :: b(ldb,nrhs) !< Right hand side -> overwritten with solution [ldb,nrhs]
INTEGER(i4), INTENT(out) :: info !< Solver return status (0 -> success)
INTEGER(i4) :: idummy(1)
DEBUG_STACK_PUSH
IF(mode==4)THEN
  idummy=0
  CALL oft_strumpack_sparse(mode,0,0,0,values,idummy,idummy,b,1,self%f_factors, &
    self%comp_type,self%comp_rel_tol,self%comp_abs_tol,self%leaf_size,self%max_rank, &
    self%min_sep_size,MERGE(1,0,self%verbose),info)
ELSE
  CALL oft_strumpack_sparse(mode,n,nnz,nrhs,values,self%lc,self%kr,b,ldb,self%f_factors, &
    self%comp_type,self%comp_rel_tol,self%comp_abs_tol,self%leaf_size,self%max_rank, &
    self%min_sep_size,MERGE(1,0,self%verbose),info)
END IF
DEBUG_STACK_POP
END SUBROUTINE strumpack_sparse_call
#endif
!------------------------------------------------------------------------------
!> Check for thread safety
!------------------------------------------------------------------------------
recursive function lusolver_check_thread(self) result(thread_safe)
class(oft_lusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
IF(self%package(1:5)=='mumps')thread_safe=.FALSE.
IF(self%package(1:6)=='superd')thread_safe=.FALSE.
end function lusolver_check_thread
!------------------------------------------------------------------------------
!> Destroy direct solver and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine lusolver_delete(self,propogate)
class(oft_lusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: mode,nrhs,ierr,ldb,nr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ivals
REAL(r8), ALLOCATABLE, DIMENSION(:) :: rvals
self%refactor=.TRUE.
self%update_graph=.TRUE.
IF(.NOT.self%initialized)THEN
  NULLIFY(self%A)
  RETURN
END IF
IF(self%dist%active)THEN
  CALL lusolver_delete_dist(self)
  IF(ASSOCIATED(self%sec_rhs))DEALLOCATE(self%sec_rhs)
  NULLIFY(self%A)
  self%initialized=.FALSE.
  RETURN
END IF
nrhs=0
ldb=0
ALLOCATE(ivals(1),rvals(1))
SELECT CASE(TRIM(self%package))
#ifdef HAVE_SUPERLU
  CASE("super")
    mode=4
    CALL oft_superlu_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
    self%superlu_struct%f_factors=C_NULL_PTR
#endif
#ifdef HAVE_SUPERLU_DIST
  CASE("superd")
    mode=4
    CALL oft_superlu_dist_dgssv(mode,self%A%nr,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%grid_handle,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    self%superlu_struct%f_factors=C_NULL_PTR
    mode=2
    nrhs=1
    CALL oft_superlu_dist_slugrid(mode,self%superlu_struct%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
    self%superlu_struct%grid_handle=C_NULL_PTR
    CALL MPI_COMM_FREE(self%superlu_struct%comm, ierr)
    IF(ASSOCIATED(self%superlu_struct%lc_orig))THEN
      DEALLOCATE(self%superlu_struct%csc_vals,self%superlu_struct%csc_map)
      DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc,self%superlu_struct%lc_orig)
    END IF
    self%superlu_struct%comm=MPI_COMM_NULL
#endif
#ifdef HAVE_UMFPACK
  CASE("umfpack")
    mode=4
    CALL oft_umfpack_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
    self%superlu_struct%f_factors=C_NULL_PTR
#endif
#ifdef HAVE_STRUMPACK
  CASE("strumpack")
    mode=4
    CALL strumpack_sparse_call(self%strumpack_struct,mode,nrhs,nrhs,nrhs,rvals,rvals,1,ierr)
    IF(ASSOCIATED(self%strumpack_struct%kr))DEALLOCATE(self%strumpack_struct%kr,self%strumpack_struct%lc)
    self%strumpack_struct%f_factors=C_NULL_PTR
#endif
#ifdef HAVE_MUMPS
  CASE("mumps")
    self%mumps_struct%job = -2
    CALL dmumps(self%mumps_struct)
#endif
#ifdef HAVE_MKL
  CASE("mkl")
    mode=-1
    CALL pardiso(self%pardiso_struct%pt,1,1,self%pardiso_struct%mtype,mode,nr, &
      rvals,ivals,ivals,self%pardiso_struct%perm, &
      1,self%pardiso_struct%iparm,self%pardiso_struct%msglvl,rvals,rvals,ierr)
#endif
  CASE("lapack")
    DEALLOCATE(self%ipiv,self%atmp)
END SELECT
DEALLOCATE(ivals,rvals)
IF(ASSOCIATED(self%sec_rhs))DEALLOCATE(self%sec_rhs)
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine lusolver_delete
!------------------------------------------------------------------------------
!> Solve the full distributed linear system using a parallel direct solver
!!
!! The system is factored across all processes using SuperLU-DIST, MUMPS, or
!! STRUMPACK with rows distributed by ownership (see @ref lu_dist_setup).
!------------------------------------------------------------------------------
SUBROUTINE lusolver_apply_dist(self,A_native,u,g)
CLASS(oft_lusolver), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), INTENT(inout) :: A_native !< Native matrix to factor
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
#ifdef HAVE_MPI
INTEGER(i4) :: mode,nrhs,ldb,ierr,nprow,npcol
LOGICAL :: new_graph
REAL(r8) :: rdummy(1)
REAL(r8), POINTER, DIMENSION(:) :: bvals
DEBUG_STACK_PUSH
IF(self%nrhs>1)CALL oft_abort('Multiple RHS not supported for distributed solves', &
  'lusolver_apply_dist',__FILE__)
!---Initialize solvers
IF(.NOT.self%initialized)THEN
  CALL MPI_Comm_dup(oft_env%COMM,self%dist%comm,ierr)
  SELECT CASE(TRIM(self%package))
    CASE("superd")
#ifdef HAVE_SUPERLU_DIST
      !---Use the most square processor grid available
      nprow=INT(SQRT(REAL(oft_env%nprocs,8)),4)
      DO WHILE(MOD(oft_env%nprocs,nprow)/=0)
        nprow=nprow-1
      END DO
      npcol=oft_env%nprocs/nprow
      mode=1
      CALL oft_superlu_dist_slugrid(mode,self%dist%comm,nprow,npcol,self%superlu_struct%grid_handle,ierr)
#else
      CALL oft_abort('OFT not compiled with SUPERLU-DIST','lusolver_apply_dist',__FILE__)
#endif
    CASE("mumps")
#ifdef HAVE_MUMPS
      self%mumps_struct%job = -1
#ifdef OFT_MPI_F08
      self%mumps_struct%comm = self%dist%comm%MPI_VAL
#else
      self%mumps_struct%comm = self%dist%comm
#endif
      self%mumps_struct%sym = 0
      self%mumps_struct%par = 1
      CALL dmumps(self%mumps_struct)
      self%mumps_struct%icntl(3:4)=-1 ! Surpress information messages
      self%mumps_struct%icntl(5)=0 ! Assembled matrix
      self%mumps_struct%icntl(18)=3 ! Distributed matrix entries
      self%mumps_struct%icntl(20)=10 ! Distributed dense RHS
      self%mumps_struct%icntl(21)=1 ! Distributed solution
#else
      CALL oft_abort('OFT not compiled with MUMPS','lusolver_apply_dist',__FILE__)
#endif
    CASE("strumpack")
#ifdef HAVE_STRUMPACK_MPI
      IF(self%iter_refine)CALL oft_abort('Iterative refinement not supported (STRUMPACK)', &
        'lusolver_apply_dist',__FILE__)
#else
      CALL oft_abort('OFT not compiled with MPI support for STRUMPACK','lusolver_apply_dist',__FILE__)
#endif
    CASE DEFAULT
      CALL oft_abort('Distributed solves require "superd", "mumps", or "strumpack" package', &
        'lusolver_apply_dist',__FILE__)
  END SELECT
  self%initialized=.TRUE.
END IF
!---Factor matrix
IF(self%refactor)THEN
  !---Full setup if requested, not yet performed, or the non-zero pattern has changed
  new_graph=self%update_graph
  IF(.NOT.new_graph)new_graph=lu_dist_graph_changed(self%dist,A_native)
  IF(new_graph)CALL lu_dist_setup(self%dist,A_native,g)
  CALL lu_dist_update(self%dist,A_native)
  SELECT CASE(TRIM(self%package))
    CASE("superd")
#ifdef HAVE_SUPERLU_DIST
      mode=3
      IF(new_graph)mode=1
      nrhs=0
      ldb=MAX(1,self%dist%nloc)
      CALL oft_superlu_dist_pdgssvx(mode,self%dist%ng,self%dist%nnz,self%dist%nloc,self%dist%offset, &
        nrhs,self%dist%vals,self%dist%lc,self%dist%kr,rdummy,ldb,self%superlu_struct%grid_handle, &
        self%superlu_struct%f_factors,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (SuperLU-DIST)','lusolver_apply_dist',__FILE__)
#endif
    CASE("mumps")
#ifdef HAVE_MUMPS
      IF(new_graph)THEN
        self%mumps_struct%n=self%dist%ng
        self%mumps_struct%nnz_loc=INT(self%dist%nnz,8)
        self%mumps_struct%irn_loc=>self%dist%irn
        self%mumps_struct%jcn_loc=>self%dist%lc
        self%mumps_struct%job=4 ! Analysis and factorization
      ELSE
        self%mumps_struct%job=2 ! Factorization only
      END IF
      self%mumps_struct%a_loc=>self%dist%vals
      CALL dmumps(self%mumps_struct)
      IF(self%mumps_struct%infog(1)<0)CALL oft_abort('Factorization failed (MUMPS)', &
        'lusolver_apply_dist',__FILE__)
      !---Allocate storage for distributed solution, whose layout is set by the factorization
      IF(ASSOCIATED(self%dist%sol_loc))DEALLOCATE(self%dist%sol_loc,self%dist%isol_loc)
      self%dist%nsol=self%mumps_struct%info(23)
      ALLOCATE(self%dist%sol_loc(MAX(1,self%dist%nsol)),self%dist%isol_loc(MAX(1,self%dist%nsol)))
      self%dist%sol_map_valid=.FALSE.
#endif
    CASE("strumpack")
#ifdef HAVE_STRUMPACK_MPI
      mode=3
      IF(new_graph)mode=1
      CALL strumpack_dist_call(self,mode,rdummy,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed (STRUMPACK)','lusolver_apply_dist',__FILE__)
#endif
  END SELECT
  self%refactor=.FALSE.
  self%update_graph=.FALSE.
END IF
!---Solve using locally-owned rows
NULLIFY(bvals)
CALL g%get_slice(bvals)
SELECT CASE(TRIM(self%package))
  CASE("superd")
#ifdef HAVE_SUPERLU_DIST
    mode=2
    nrhs=1
    ldb=MAX(1,self%dist%nloc)
    CALL oft_superlu_dist_pdgssvx(mode,self%dist%ng,self%dist%nnz,self%dist%nloc,self%dist%offset, &
      nrhs,self%dist%vals,self%dist%lc,self%dist%kr,bvals,ldb,self%superlu_struct%grid_handle, &
      self%superlu_struct%f_factors,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (SuperLU-DIST)','lusolver_apply_dist',__FILE__)
#endif
  CASE("mumps")
#ifdef HAVE_MUMPS
    !---Solution is distributed by MUMPS and then sent to the owning processes
    self%mumps_struct%nrhs=1
    self%mumps_struct%nloc_rhs=self%dist%nloc
    self%mumps_struct%lrhs_loc=MAX(1,self%dist%nloc)
    self%mumps_struct%irhs_loc=>self%dist%rows
    self%mumps_struct%rhs_loc=>bvals
    self%mumps_struct%lsol_loc=MAX(1,self%dist%nsol)
    self%mumps_struct%sol_loc=>self%dist%sol_loc
    self%mumps_struct%isol_loc=>self%dist%isol_loc
    self%mumps_struct%job=3 ! Solve phase
    CALL dmumps(self%mumps_struct)
    IF(self%mumps_struct%infog(1)<0)CALL oft_abort('Solve failed (MUMPS)','lusolver_apply_dist',__FILE__)
    NULLIFY(self%mumps_struct%rhs_loc,self%mumps_struct%irhs_loc)
    NULLIFY(self%mumps_struct%sol_loc,self%mumps_struct%isol_loc)
    CALL lu_dist_sol_scatter(self%dist,bvals)
#endif
  CASE("strumpack")
#ifdef HAVE_STRUMPACK_MPI
    mode=2
    CALL strumpack_dist_call(self,mode,bvals,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed (STRUMPACK)','lusolver_apply_dist',__FILE__)
#endif
END SELECT
CALL u%restore_slice(bvals)
DEALLOCATE(bvals)
DEBUG_STACK_POP
#else
CALL oft_abort('Distributed solves require MPI','lusolver_apply_dist',__FILE__)
#endif
END SUBROUTINE lusolver_apply_dist
!------------------------------------------------------------------------------
!> Destroy distributed direct solver and deallocate all internal storage
!------------------------------------------------------------------------------
SUBROUTINE lusolver_delete_dist(self)
CLASS(oft_lusolver), INTENT(inout) :: self !< Solver object
#ifdef HAVE_MPI
INTEGER(i4) :: mode,nrhs,ldb,ierr,ivals(1)
REAL(r8) :: rvals(1)
DEBUG_STACK_PUSH
SELECT CASE(TRIM(self%package))
#ifdef HAVE_SUPERLU_DIST
  CASE("superd")
    nrhs=0
    ldb=1
    IF(C_ASSOCIATED(self%superlu_struct%f_factors))THEN
      mode=4
      CALL oft_superlu_dist_pdgssvx(mode,nrhs,nrhs,nrhs,nrhs,nrhs,rvals,ivals,ivals,rvals,ldb, &
        self%superlu_struct%grid_handle,self%superlu_struct%f_factors,self%iter_refine,ierr)
      self%superlu_struct%f_factors=C_NULL_PTR
    END IF
    IF(C_ASSOCIATED(self%superlu_struct%grid_handle))THEN
      mode=2
      nrhs=1
      CALL oft_superlu_dist_slugrid(mode,self%dist%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
      self%superlu_struct%grid_handle=C_NULL_PTR
    END IF
#endif
#ifdef HAVE_MUMPS
  CASE("mumps")
    self%mumps_struct%job = -2
    CALL dmumps(self%mumps_struct)
    NULLIFY(self%mumps_struct%irn_loc,self%mumps_struct%jcn_loc,self%mumps_struct%a_loc)
#endif
#ifdef HAVE_STRUMPACK_MPI
  CASE("strumpack")
    IF(C_ASSOCIATED(self%strumpack_struct%f_factors))THEN
      mode=4
      CALL strumpack_dist_call(self,mode,rvals,ierr)
    END IF
#endif
END SELECT
CALL MPI_Comm_free(self%dist%comm,ierr)
self%dist%comm=MPI_COMM_NULL
CALL lu_dist_delete(self%dist)
self%dist%active=.FALSE.
DEBUG_STACK_POP
#endif
END SUBROUTINE lusolver_delete_dist
#ifdef HAVE_STRUMPACK_MPI
!------------------------------------------------------------------------------
!> Call distributed STRUMPACK interface using stored row distribution and options
!------------------------------------------------------------------------------
SUBROUTINE strumpack_dist_call(self,mode,b,info)
CLASS(oft_lusolver), INTENT(inout) :: self !< Solver object
INTEGER(i4), INTENT(in) :: mode !< Operation (see @ref oft_strumpack_dist)
REAL(r8), INTENT(inout) :: b(*) !< Locally-owned rows of RHS -> overwritten with solution (`mode=2`)
INTEGER(i4), INTENT(out) :: info !< Solver return status (0 -> success)
INTEGER(i4) :: fcomm,idummy(1)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: dist_rows
DEBUG_STACK_PUSH
#ifdef OFT_MPI_F08
fcomm=self%dist%comm%MPI_VAL
#else
fcomm=self%dist%comm
#endif
idummy=0
SELECT CASE(mode)
  CASE(1,3)
    !---Row distribution (first row on each process and total)
    ALLOCATE(dist_rows(oft_env%nprocs+1))
    dist_rows(1:oft_env%nprocs)=self%dist%displs
    dist_rows(oft_env%nprocs+1)=self%dist%ng
    CALL oft_strumpack_dist(mode,self%dist%ng,self%dist%nnz,self%dist%nloc,self%dist%vals, &
      self%dist%lc,self%dist%kr,dist_rows,b,fcomm, &
      self%strumpack_struct%f_factors,self%strumpack_struct%comp_type,self%strumpack_struct%comp_rel_tol, &
      self%strumpack_struct%comp_abs_tol,self%strumpack_struct%leaf_size,self%strumpack_struct%max_rank, &
      self%strumpack_struct%min_sep_size,MERGE(1,0,self%strumpack_struct%verbose),info)
    DEALLOCATE(dist_rows)
  CASE(2)
    CALL oft_strumpack_dist(mode,self%dist%ng,self%dist%nnz,self%dist%nloc,self%dist%vals, &
      self%dist%lc,self%dist%kr,idummy,b,fcomm, &
      self%strumpack_struct%f_factors,self%strumpack_struct%comp_type,self%strumpack_struct%comp_rel_tol, &
      self%strumpack_struct%comp_abs_tol,self%strumpack_struct%leaf_size,self%strumpack_struct%max_rank, &
      self%strumpack_struct%min_sep_size,MERGE(1,0,self%strumpack_struct%verbose),info)
  CASE DEFAULT
    CALL oft_strumpack_dist(mode,0,0,0,b,idummy,idummy,idummy,b,fcomm, &
      self%strumpack_struct%f_factors,self%strumpack_struct%comp_type,self%strumpack_struct%comp_rel_tol, &
      self%strumpack_struct%comp_abs_tol,self%strumpack_struct%leaf_size,self%strumpack_struct%max_rank, &
      self%strumpack_struct%min_sep_size,MERGE(1,0,self%strumpack_struct%verbose),info)
END SELECT
DEBUG_STACK_POP
END SUBROUTINE strumpack_dist_call
#endif
!------------------------------------------------------------------------------
!> Setup row-distributed representation of a native matrix
!!
!! Each process takes ownership of the rows corresponding to its locally-owned
!! vector entries, which are numbered contiguously in rank order. Contributions
!! to these rows from other processes (partial sums and couplings to entries not
!! present locally) are gathered from neighboring processes. The resulting
!! communication pattern is saved for value updates in @ref lu_dist_update.
!!
!! @note All communication uses the private communicator `self%comm`, which must be set
!------------------------------------------------------------------------------
SUBROUTINE lu_dist_setup(self,A,g)
TYPE(lu_dist_struc), INTENT(inout) :: self !< Distributed matrix structure
CLASS(oft_native_matrix), INTENT(in) :: A !< Native matrix
CLASS(oft_vector), INTENT(inout) :: g !< Vector defining row ownership and linkage
#ifdef HAVE_MPI
INTEGER(i4), PARAMETER :: tag_count = 7301, tag_pattern = 7302
INTEGER(i4) :: i,j,k,m,jp,jn,last_col,nlocal_cand,ncand,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: gind,owner,con_ind,nsend,nrecv,cand_row,cand_col, &
  cand_id,cand_dst,row_ptr,row_fill,sort_col
INTEGER(i8) :: ng8
REAL(r8), POINTER, DIMENSION(:) :: vtmp
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_1d_int), ALLOCATABLE, DIMENSION(:) :: send_ij,recv_ij
#ifdef OFT_MPI_F08
TYPE(mpi_request), ALLOCATABLE, DIMENSION(:) :: reqs
#else
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: reqs
#endif
DEBUG_STACK_PUSH
CALL lu_dist_delete(self)
IF(A%nr/=g%n)CALL oft_abort('Matrix and vector sizes do not match','lu_dist_setup',__FILE__)
!---Number locally-owned rows contiguously in rank order
self%nloc=g%nslice
ALLOCATE(self%counts(oft_env%nprocs),self%displs(oft_env%nprocs))
CALL MPI_Allgather(self%nloc,1,OFT_MPI_I4,self%counts,1,OFT_MPI_I4,self%comm,ierr)
ng8=SUM(INT(self%counts,8))
IF(ng8>INT(HUGE(self%ng),8))CALL oft_abort('Global system too large for 32-bit indexing', &
  'lu_dist_setup',__FILE__)
self%ng=INT(ng8,4)
self%displs(1)=0
DO i=2,oft_env%nprocs
  self%displs(i)=self%displs(i-1)+self%counts(i-1)
END DO
self%offset=self%displs(oft_env%rank+1)
ALLOCATE(self%rows(self%nloc))
DO i=1,self%nloc
  self%rows(i)=self%offset+i
END DO
!---Propagate global indices from owners to all local copies
ALLOCATE(vtmp(self%nloc))
vtmp=REAL(self%rows,8)
CALL g%new(tmp_vec)
CALL tmp_vec%restore_slice(vtmp)
DEALLOCATE(vtmp)
NULLIFY(vtmp)
CALL tmp_vec%get_local(vtmp)
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)
ALLOCATE(gind(A%nr))
gind=NINT(vtmp,4)
DEALLOCATE(vtmp)
IF(ANY(gind<1).OR.ANY(gind>self%ng))CALL oft_abort('Invalid global row index','lu_dist_setup',__FILE__)
!---Determine owner of each local row
ALLOCATE(owner(A%nr))
DO i=1,A%nr
  owner(i)=lu_dist_owner(self%displs,gind(i))
END DO
!---Map neighboring ranks to connection index
self%nproc_con=g%stitch_info%nproc_con
ALLOCATE(self%proc_con(self%nproc_con))
ALLOCATE(con_ind(0:oft_env%nprocs-1))
con_ind=0
DO j=1,self%nproc_con
  self%proc_con(j)=g%stitch_info%proc_con(j)
  con_ind(self%proc_con(j))=j
END DO
!---Count entries in locally-owned rows and rows owned by each neighbor
ALLOCATE(nsend(self%nproc_con),nrecv(self%nproc_con))
nsend=0
nrecv=0
nlocal_cand=0
DO i=1,A%nr
  IF(owner(i)==oft_env%rank)THEN
    nlocal_cand=nlocal_cand+A%kr(i+1)-A%kr(i)
  ELSE
    j=con_ind(owner(i))
    IF(j==0)CALL oft_abort('Row owner is not a neighboring process','lu_dist_setup',__FILE__)
    nsend(j)=nsend(j)+A%kr(i+1)-A%kr(i)
  END IF
END DO
ALLOCATE(reqs(2*self%nproc_con))
DO j=1,self%nproc_con
  CALL MPI_Irecv(nrecv(j),1,OFT_MPI_I4,self%proc_con(j),tag_count,self%comm,reqs(j),ierr)
  CALL MPI_Isend(nsend(j),1,OFT_MPI_I4,self%proc_con(j),tag_count,self%comm, &
    reqs(self%nproc_con+j),ierr)
END DO
CALL MPI_Waitall(2*self%nproc_con,reqs,MPI_STATUSES_IGNORE,ierr)
!---Build candidate entries for locally-owned rows and send lists for other rows
ncand=nlocal_cand+SUM(nrecv)
ALLOCATE(cand_row(ncand),cand_col(ncand),self%loc_src(nlocal_cand))
ALLOCATE(self%send_map(self%nproc_con),send_ij(self%nproc_con))
DO j=1,self%nproc_con
  self%send_map(j)%n=nsend(j)
  ALLOCATE(self%send_map(j)%v(nsend(j)),send_ij(j)%v(2*nsend(j)))
END DO
nsend=0
k=0
DO i=1,A%nr
  IF(owner(i)==oft_env%rank)THEN
    DO m=A%kr(i),A%kr(i+1)-1
      k=k+1
      self%loc_src(k)=m
      cand_row(k)=gind(i)-self%offset
      cand_col(k)=gind(A%lc(m))
    END DO
  ELSE
    j=con_ind(owner(i))
    DO m=A%kr(i),A%kr(i+1)-1
      nsend(j)=nsend(j)+1
      self%send_map(j)%v(nsend(j))=m
      send_ij(j)%v(2*nsend(j)-1)=gind(i)
      send_ij(j)%v(2*nsend(j))=gind(A%lc(m))
    END DO
  END IF
END DO
DEALLOCATE(gind,owner,con_ind)
!---Exchange (row,column) pattern with neighbors
ALLOCATE(self%recv_map(self%nproc_con),recv_ij(self%nproc_con))
DO j=1,self%nproc_con
  self%recv_map(j)%n=nrecv(j)
  ALLOCATE(self%recv_map(j)%v(nrecv(j)),recv_ij(j)%v(2*nrecv(j)))
END DO
reqs=MPI_REQUEST_NULL
DO j=1,self%nproc_con
  IF(nrecv(j)>0)CALL MPI_Irecv(recv_ij(j)%v,2*nrecv(j),OFT_MPI_I4,self%proc_con(j),tag_pattern, &
    self%comm,reqs(j),ierr)
  IF(nsend(j)>0)CALL MPI_Isend(send_ij(j)%v,2*nsend(j),OFT_MPI_I4,self%proc_con(j),tag_pattern, &
    self%comm,reqs(self%nproc_con+j),ierr)
END DO
CALL MPI_Waitall(2*self%nproc_con,reqs,MPI_STATUSES_IGNORE,ierr)
DEALLOCATE(reqs)
k=nlocal_cand
DO j=1,self%nproc_con
  DO m=1,nrecv(j)
    k=k+1
    cand_row(k)=recv_ij(j)%v(2*m-1)-self%offset
    cand_col(k)=recv_ij(j)%v(2*m)
  END DO
  DEALLOCATE(send_ij(j)%v,recv_ij(j)%v)
END DO
DEALLOCATE(send_ij,recv_ij)
IF(ANY(cand_row<1).OR.ANY(cand_row>self%nloc))CALL oft_abort('Received entry for row not owned locally', &
  'lu_dist_setup',__FILE__)
!---Bucket candidates by row
ALLOCATE(row_ptr(self%nloc+1),row_fill(self%nloc),sort_col(ncand),cand_id(ncand),cand_dst(ncand))
row_ptr=0
DO k=1,ncand
  row_ptr(cand_row(k)+1)=row_ptr(cand_row(k)+1)+1
END DO
row_ptr(1)=1
DO i=1,self%nloc
  row_ptr(i+1)=row_ptr(i+1)+row_ptr(i)
END DO
row_fill=row_ptr(1:self%nloc)
DO k=1,ncand
  i=cand_row(k)
  sort_col(row_fill(i))=cand_col(k)
  cand_id(row_fill(i))=k
  row_fill(i)=row_fill(i)+1
END DO
DEALLOCATE(cand_row,cand_col,row_fill)
!---Sort columns in each row and merge duplicate entries
ALLOCATE(self%kr(self%nloc+1))
self%nnz=0
DO i=1,self%nloc
  jp=row_ptr(i)
  jn=row_ptr(i+1)-jp
  IF(jn>1)CALL sort_array(sort_col(jp:jp+jn-1),cand_id(jp:jp+jn-1),jn)
  self%kr(i)=self%nnz+1
  last_col=-1
  DO m=jp,jp+jn-1
    IF(sort_col(m)/=last_col)THEN
      self%nnz=self%nnz+1
      sort_col(self%nnz)=sort_col(m) ! Compact in place (nnz <= m)
      last_col=sort_col(m)
    END IF
    cand_dst(cand_id(m))=self%nnz
  END DO
END DO
self%kr(self%nloc+1)=self%nnz+1
ALLOCATE(self%lc(self%nnz),self%irn(self%nnz),self%vals(self%nnz))
self%lc=sort_col(1:self%nnz)
DO i=1,self%nloc
  self%irn(self%kr(i):self%kr(i+1)-1)=self%rows(i)
END DO
!---Save destinations for local and received values
ALLOCATE(self%loc_dst(nlocal_cand))
self%loc_dst=cand_dst(1:nlocal_cand)
k=nlocal_cand
DO j=1,self%nproc_con
  self%recv_map(j)%v=cand_dst(k+1:k+nrecv(j))
  k=k+nrecv(j)
END DO
DEALLOCATE(row_ptr,sort_col,cand_id,cand_dst,nsend,nrecv)
!---Save native matrix graph to detect changes in non-zero pattern
ALLOCATE(self%src_kr(A%nr+1),self%src_lc(A%nnz))
self%src_kr=A%kr(1:A%nr+1)
self%src_lc=A%lc(1:A%nnz)
DEBUG_STACK_POP
#endif
END SUBROUTINE lu_dist_setup
!------------------------------------------------------------------------------
!> Update values of row-distributed matrix from a native matrix
!!
!! Uses communication pattern constructed by @ref lu_dist_setup
!------------------------------------------------------------------------------
SUBROUTINE lu_dist_update(self,A)
TYPE(lu_dist_struc), INTENT(inout) :: self !< Distributed matrix structure
CLASS(oft_native_matrix), INTENT(in) :: A !< Native matrix
#ifdef HAVE_MPI
INTEGER(i4), PARAMETER :: tag_vals = 7303
INTEGER(i4) :: j,k,ierr
TYPE(oft_1d_real), ALLOCATABLE, DIMENSION(:) :: sbuf,rbuf
#ifdef OFT_MPI_F08
TYPE(mpi_request), ALLOCATABLE, DIMENSION(:) :: reqs
#else
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: reqs
#endif
DEBUG_STACK_PUSH
ALLOCATE(sbuf(self%nproc_con),rbuf(self%nproc_con),reqs(2*self%nproc_con))
reqs=MPI_REQUEST_NULL
DO j=1,self%nproc_con
  rbuf(j)%n=self%recv_map(j)%n
  ALLOCATE(rbuf(j)%v(rbuf(j)%n))
  IF(rbuf(j)%n>0)CALL MPI_Irecv(rbuf(j)%v,rbuf(j)%n,OFT_MPI_R8,self%proc_con(j),tag_vals, &
    self%comm,reqs(j),ierr)
END DO
DO j=1,self%nproc_con
  sbuf(j)%n=self%send_map(j)%n
  ALLOCATE(sbuf(j)%v(sbuf(j)%n))
  DO k=1,sbuf(j)%n
    sbuf(j)%v(k)=A%M(self%send_map(j)%v(k))
  END DO
  IF(sbuf(j)%n>0)CALL MPI_Isend(sbuf(j)%v,sbuf(j)%n,OFT_MPI_R8,self%proc_con(j),tag_vals, &
    self%comm,reqs(self%nproc_con+j),ierr)
END DO
!---Add local contributions (unstitched partial sums)
self%vals=0.d0
DO k=1,SIZE(self%loc_src)
  self%vals(self%loc_dst(k))=self%vals(self%loc_dst(k))+A%M(self%loc_src(k))
END DO
!---Add contributions from neighbors
CALL MPI_Waitall(2*self%nproc_con,reqs,MPI_STATUSES_IGNORE,ierr)
DO j=1,self%nproc_con
  DO k=1,rbuf(j)%n
    self%vals(self%recv_map(j)%v(k))=self%vals(self%recv_map(j)%v(k))+rbuf(j)%v(k)
  END DO
  DEALLOCATE(rbuf(j)%v,sbuf(j)%v)
END DO
DEALLOCATE(sbuf,rbuf,reqs)
DEBUG_STACK_POP
#endif
END SUBROUTINE lu_dist_update
!------------------------------------------------------------------------------
!> Get the process that owns a given global row
!------------------------------------------------------------------------------
PURE FUNCTION lu_dist_owner(displs,grow) RESULT(owner)
INTEGER(i4), INTENT(in) :: displs(:) !< Number of rows owned by lower ranks for each process [nprocs]
INTEGER(i4), INTENT(in) :: grow !< Global row index (1-based)
INTEGER(i4) :: owner !< Rank of owning process
INTEGER(i4) :: lo,hi,mid
lo=1
hi=SIZE(displs)
DO WHILE(lo<hi)
  mid=(lo+hi+1)/2
  IF(displs(mid)<grow)THEN
    lo=mid
  ELSE
    hi=mid-1
  END IF
END DO
owner=lo-1
END FUNCTION lu_dist_owner
!------------------------------------------------------------------------------
!> Check if the non-zero pattern of a native matrix differs from the one used in @ref lu_dist_setup
!!
!! @note Collective over `self%comm`; the result is true on all processes if the pattern
!! has changed on any process or setup has not been performed
!------------------------------------------------------------------------------
FUNCTION lu_dist_graph_changed(self,A) RESULT(changed)
TYPE(lu_dist_struc), INTENT(in) :: self !< Distributed matrix structure
CLASS(oft_native_matrix), INTENT(in) :: A !< Native matrix
LOGICAL :: changed !< Has the pattern changed?
#ifdef HAVE_MPI
LOGICAL :: changed_loc
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
changed_loc=.TRUE.
IF(ASSOCIATED(self%src_kr))THEN
  IF((SIZE(self%src_kr)==A%nr+1).AND.(SIZE(self%src_lc)==A%nnz))THEN
    changed_loc=ANY(self%src_kr/=A%kr(1:A%nr+1))
    IF(.NOT.changed_loc)changed_loc=ANY(self%src_lc/=A%lc(1:A%nnz))
  END IF
END IF
CALL MPI_Allreduce(changed_loc,changed,1,OFT_MPI_LOGICAL,MPI_LOR,self%comm,ierr)
DEBUG_STACK_POP
#else
changed=.TRUE.
#endif
END FUNCTION lu_dist_graph_changed
!------------------------------------------------------------------------------
!> Send entries of a distributed solution returned by the solver package to the owning processes
!!
!! The redistribution map is built on the first call after each factorization, once the global
!! indices of the returned entries (`isol_loc`) have been set by the solver package.
!------------------------------------------------------------------------------
SUBROUTINE lu_dist_sol_scatter(self,b)
TYPE(lu_dist_struc), INTENT(inout) :: self !< Distributed matrix structure
REAL(r8), INTENT(inout) :: b(:) !< Locally-owned rows of solution [nloc]
#ifdef HAVE_MPI
INTEGER(i4) :: i,k,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: dest,fill,isbuf,irbuf
REAL(r8), ALLOCATABLE, DIMENSION(:) :: sbuf,rbuf
DEBUG_STACK_PUSH
IF(.NOT.self%sol_map_valid)THEN
  CALL lu_dist_sol_map_delete(self)
  ALLOCATE(self%sol_perm(self%nsol),self%sol_ridx(self%nloc))
  ALLOCATE(self%sol_scounts(oft_env%nprocs),self%sol_sdispls(oft_env%nprocs))
  ALLOCATE(self%sol_rcounts(oft_env%nprocs),self%sol_rdispls(oft_env%nprocs))
  !---Order returned entries by owning process
  ALLOCATE(dest(self%nsol),fill(oft_env%nprocs))
  self%sol_scounts=0
  DO k=1,self%nsol
    dest(k)=lu_dist_owner(self%displs,self%isol_loc(k))+1
    self%sol_scounts(dest(k))=self%sol_scounts(dest(k))+1
  END DO
  self%sol_sdispls(1)=0
  DO i=2,oft_env%nprocs
    self%sol_sdispls(i)=self%sol_sdispls(i-1)+self%sol_scounts(i-1)
  END DO
  fill=self%sol_sdispls
  DO k=1,self%nsol
    fill(dest(k))=fill(dest(k))+1
    self%sol_perm(fill(dest(k)))=k
  END DO
  DEALLOCATE(dest,fill)
  !---Exchange global indices to map received entries to locally-owned rows
  CALL MPI_Alltoall(self%sol_scounts,1,OFT_MPI_I4,self%sol_rcounts,1,OFT_MPI_I4,self%comm,ierr)
  IF(SUM(self%sol_rcounts)/=self%nloc)CALL oft_abort('Incomplete distributed solution', &
    'lu_dist_sol_scatter',__FILE__)
  self%sol_rdispls(1)=0
  DO i=2,oft_env%nprocs
    self%sol_rdispls(i)=self%sol_rdispls(i-1)+self%sol_rcounts(i-1)
  END DO
  ALLOCATE(isbuf(MAX(1,self%nsol)),irbuf(MAX(1,self%nloc)))
  isbuf(1:self%nsol)=self%isol_loc(self%sol_perm)
  CALL MPI_Alltoallv(isbuf,self%sol_scounts,self%sol_sdispls,OFT_MPI_I4, &
    irbuf,self%sol_rcounts,self%sol_rdispls,OFT_MPI_I4,self%comm,ierr)
  self%sol_ridx=irbuf(1:self%nloc)-self%offset
  DEALLOCATE(isbuf,irbuf)
  IF(ANY(self%sol_ridx<1).OR.ANY(self%sol_ridx>self%nloc))CALL oft_abort( &
    'Received solution entry for row not owned locally','lu_dist_sol_scatter',__FILE__)
  self%sol_map_valid=.TRUE.
END IF
!---Send solution entries to owning processes
ALLOCATE(sbuf(MAX(1,self%nsol)),rbuf(MAX(1,self%nloc)))
sbuf(1:self%nsol)=self%sol_loc(self%sol_perm)
CALL MPI_Alltoallv(sbuf,self%sol_scounts,self%sol_sdispls,OFT_MPI_R8, &
  rbuf,self%sol_rcounts,self%sol_rdispls,OFT_MPI_R8,self%comm,ierr)
b(self%sol_ridx)=rbuf(1:self%nloc)
DEALLOCATE(sbuf,rbuf)
DEBUG_STACK_POP
#endif
END SUBROUTINE lu_dist_sol_scatter
!------------------------------------------------------------------------------
!> Deallocate solution redistribution map built by @ref lu_dist_sol_scatter
!------------------------------------------------------------------------------
SUBROUTINE lu_dist_sol_map_delete(self)
TYPE(lu_dist_struc), INTENT(inout) :: self !< Distributed matrix structure
IF(ASSOCIATED(self%sol_perm))DEALLOCATE(self%sol_perm)
IF(ASSOCIATED(self%sol_scounts))DEALLOCATE(self%sol_scounts)
IF(ASSOCIATED(self%sol_sdispls))DEALLOCATE(self%sol_sdispls)
IF(ASSOCIATED(self%sol_rcounts))DEALLOCATE(self%sol_rcounts)
IF(ASSOCIATED(self%sol_rdispls))DEALLOCATE(self%sol_rdispls)
IF(ASSOCIATED(self%sol_ridx))DEALLOCATE(self%sol_ridx)
self%sol_map_valid=.FALSE.
END SUBROUTINE lu_dist_sol_map_delete
!------------------------------------------------------------------------------
!> Deallocate internal storage for row-distributed matrix
!!
!! @note The private communicator is not freed
!------------------------------------------------------------------------------
SUBROUTINE lu_dist_delete(self)
TYPE(lu_dist_struc), INTENT(inout) :: self !< Distributed matrix structure
INTEGER(i4) :: j
IF(ASSOCIATED(self%proc_con))DEALLOCATE(self%proc_con)
IF(ASSOCIATED(self%counts))DEALLOCATE(self%counts)
IF(ASSOCIATED(self%displs))DEALLOCATE(self%displs)
IF(ASSOCIATED(self%rows))DEALLOCATE(self%rows)
IF(ASSOCIATED(self%kr))DEALLOCATE(self%kr)
IF(ASSOCIATED(self%lc))DEALLOCATE(self%lc)
IF(ASSOCIATED(self%irn))DEALLOCATE(self%irn)
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
IF(ASSOCIATED(self%loc_src))DEALLOCATE(self%loc_src)
IF(ASSOCIATED(self%loc_dst))DEALLOCATE(self%loc_dst)
IF(ASSOCIATED(self%send_map))THEN
  DO j=1,SIZE(self%send_map)
    IF(ASSOCIATED(self%send_map(j)%v))DEALLOCATE(self%send_map(j)%v)
  END DO
  DEALLOCATE(self%send_map)
END IF
IF(ASSOCIATED(self%recv_map))THEN
  DO j=1,SIZE(self%recv_map)
    IF(ASSOCIATED(self%recv_map(j)%v))DEALLOCATE(self%recv_map(j)%v)
  END DO
  DEALLOCATE(self%recv_map)
END IF
IF(ASSOCIATED(self%src_kr))DEALLOCATE(self%src_kr)
IF(ASSOCIATED(self%src_lc))DEALLOCATE(self%src_lc)
IF(ASSOCIATED(self%sol_loc))DEALLOCATE(self%sol_loc)
IF(ASSOCIATED(self%isol_loc))DEALLOCATE(self%isol_loc)
CALL lu_dist_sol_map_delete(self)
self%nsol=0
self%ng=0
self%nloc=0
self%offset=0
self%nnz=0
self%nproc_con=0
END SUBROUTINE lu_dist_delete
!------------------------------------------------------------------------------
!> Solve a linear system using a direct solve
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE ilusolver_apply(self,u,g)
CLASS(oft_ilusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: mode,nrhs,ldb,ierr,i,j,k,info
INTEGER(i4), POINTER :: csr_map(:),kr_tmp(:)
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: mat_vals
REAL(r8), POINTER, DIMENSION(:) :: csc_vals,vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vals,b
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: mytimer
TYPE(oft_graph) :: csr_graph
TYPE(oft_graph) :: csc_graph
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting ILU solver: ',self%package,self%refactor
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'ilusolver_apply',__FILE__)
IF(ASSOCIATED(A_native%Mfull))THEN
  mat_vals=>A_native%Mfull
ELSE
  mat_vals=>A_native%M
END IF
!---Initialize solvers
IF(.NOT.self%initialized)THEN
  !---
  SELECT CASE(TRIM(self%package))
    CASE DEFAULT
      CALL oft_abort('Unknown ILU package','ilusolver_apply', __FILE__)
    CASE("native")
      ALLOCATE(self%native_struct%ju(A_native%nr))
      ALLOCATE(self%native_struct%jlu(A_native%nnz+1))
      ALLOCATE(self%native_struct%alu(A_native%nnz+1))
!     CASE("super")
! #ifdef HAVE_SUPERLU
!       ALLOCATE(self%superlu_struct%kr(A_native%nr+1))
!       ALLOCATE(self%superlu_struct%lc(A_native%nnz))
!       self%superlu_struct%kr=A_native%kr-1
!       self%superlu_struct%lc=A_native%lc-1
! #else
!       CALL oft_abort('OFT not compiled with SUPERLU','ilusolver_apply',__FILE__)
! #endif
    CASE("mkl")
#if !defined(HAVE_MKL)
      CALL oft_abort('OFT not compiled with MKL-PARDISO','ilusolver_apply',__FILE__)
#endif
  END SELECT
  self%initialized=.TRUE.
END IF
!---Call factorization
ALLOCATE(vals(u%n,self%nrhs))
nrhs=1
vtmp=>vals(:,1)
CALL g%get_local(vtmp)
IF(self%nrhs>1)THEN
  nrhs=self%nrhs
  DO i=2,self%nrhs
    vals(:,i) = self%sec_rhs(:,i-1)
  END DO
END IF
ldb=A_native%nr
! NULLIFY(vals)
! CALL g%get_local(vals)
SELECT CASE(TRIM(self%package))
  CASE("mkl")
#ifdef HAVE_MKL
    IF(self%refactor)THEN
      !---Update matrix
      CALL oft_mkl_ilu(1,A_native%nr,A_native%nnz,mat_vals,A_native%lc, &
        A_native%kr,vals,self%pardiso_struct%f_factors,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed','ilusolver_apply',__FILE__)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    CALL oft_mkl_ilu(2,A_native%nr,A_native%nnz,mat_vals,A_native%lc, &
      A_native%kr,vals,self%pardiso_struct%f_factors,ierr)
#else
    CALL oft_abort('OFT not compiled with MKL-PARDISO','ilusolver_apply',__FILE__)
#endif
!   CASE("super")
! #ifdef HAVE_SUPERLU
!     IF(self%refactor)THEN
!       mode=3
!       IF(self%update_graph)mode=1
! #if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
!       !$omp critical (superlu_solve)
! #endif
!       CALL oft_superlu_dgsisx(mode,A_native%nr,A_native%nnz,nrhs, &
!         mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
!         self%superlu_struct%f_factors,self%superlu_struct%col_perm,1.d-1,ierr)
!       IF(ierr/=0)THEN
!         WRITE(*,*)ierr,A_native%nr
!         IF(.NOT.((ierr>=0).AND.(ierr<=A_native%nr)))THEN
!           CALL oft_abort('Factorization failed','ilusolver_apply',__FILE__)
!         END IF
!       END IF
! #if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
!       !$omp end critical (superlu_solve)
! #endif
!       self%refactor=.FALSE.
!       self%update_graph=.FALSE.
!     END IF
!     mode=2
! #if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
!     !$omp critical (superlu_solve)
! #endif
!     CALL oft_superlu_dgsisx(mode,A_native%nr,A_native%nnz,nrhs, &
!       mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
!       self%superlu_struct%f_factors,self%superlu_struct%col_perm,1.d3,ierr)
!       IF(ierr/=0)THEN
!         WRITE(*,*)ierr
!         CALL oft_abort('Solve failed','ilusolver_apply',__FILE__)
!       END IF
! #if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
!     !$omp end critical (superlu_solve)
! #endif
! #else
!   CALL oft_abort('OFT not compiled with SuperLU','ilusolver_apply',__FILE__)
! #endif
  CASE("native")
    IF(self%refactor)THEN
      CALL ilu0(A_native%nr,mat_vals,A_native%lc,A_native%kr,self%native_struct%alu, &
        self%native_struct%jlu,self%native_struct%ju, ierr)
      IF(ierr/=0)THEN
        WRITE(*,*)ierr
        CALL oft_abort('Factorization failed','ilusolver_apply',__FILE__)
      END IF
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    CALL lusol(A_native%nr,nrhs,vals,self%native_struct%alu,self%native_struct%jlu, &
      self%native_struct%ju)
  CASE DEFAULT
    CALL oft_abort('Unknown factorization package','ilusolver_apply',__FILE__)
END SELECT
CALL u%restore_local(vals(:,1))
IF(self%nrhs>1)THEN
  DO i=2,self%nrhs
    self%sec_rhs(:,i-1) = vals(:,i)
  END DO
END IF
DEALLOCATE(vals)
DEBUG_STACK_POP
END SUBROUTINE ilusolver_apply
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine ilusolver_update(self,new_pattern)
class(oft_ilusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))self%update_graph=new_pattern
end subroutine ilusolver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition
!------------------------------------------------------------------------------
subroutine ilusolver_setup_xml(self,solver_node,level)
CLASS(oft_ilusolver), INTENT(inout) :: self
TYPE(xml_node), INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
!---
CHARACTER(LEN=:), ALLOCATABLE :: factor_package
CHARACTER(LEN=3) :: fac_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---
CALL xml_get_element(solver_node,"package",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,factor_package,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `package` node","ilusolver_setup_xml",__FILE__)
  IF(ALLOCATED(factor_package))THEN
    IF(LEN(factor_package)>7)CALL oft_abort('Factorization package name too long','ilusolver_setup_xml',__FILE__)
    self%package=factor_package
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'LU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:  ',self%package
END IF
DEBUG_STACK_POP
end subroutine ilusolver_setup_xml
!------------------------------------------------------------------------------
!> Check for thread safety
!------------------------------------------------------------------------------
recursive function ilusolver_check_thread(self) result(thread_safe)
class(oft_ilusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
! IF(self%package(1:5)=='mumps')thread_safe=.FALSE.
! IF(self%package(1:6)=='superd')thread_safe=.FALSE.
end function ilusolver_check_thread
!------------------------------------------------------------------------------
!> Destroy direct solver and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine ilusolver_delete(self,propogate)
class(oft_ilusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: mode,nrhs,ierr,ldb,nr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ivals
REAL(r8), ALLOCATABLE, DIMENSION(:) :: rvals
self%refactor=.TRUE.
self%update_graph=.TRUE.
IF(.NOT.self%initialized)THEN
  NULLIFY(self%A)
  RETURN
END IF
ALLOCATE(ivals(1),rvals(1))
SELECT CASE(TRIM(self%package))
! #ifdef HAVE_SUPERLU
!   CASE("super")
!     mode=4
!     CALL oft_superlu_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
!       rvals,ldb,self%superlu_struct%f_factors,nrhs,ierr)
! #endif
#ifdef HAVE_MKL
  CASE("mkl")
    CALL oft_mkl_ilu(3,1,1,rvals,ivals,ivals,rvals,self%pardiso_struct%f_factors,ierr)
#endif
  CASE("native")
    IF(ASSOCIATED(self%native_struct%alu))THEN
      DEALLOCATE(self%native_struct%alu,self%native_struct%jlu,self%native_struct%ju)
    END IF
END SELECT
DEALLOCATE(ivals,rvals)
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine ilusolver_delete
!---------------------------------------------------------------------------------
!> Invert a general real matrix in place using LAPACK (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_real(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(8), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,LWORK,info
REAL(8) :: elapsed_time
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ipiv
REAL(8), ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL dgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'DGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL dgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'DGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(rwork(1),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL dgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'DGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_matinv_real
!---------------------------------------------------------------------------------
!> Invert a general complex matrix in place using LAPACK (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_complex(nrows,Amat,error)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows/columns
DOUBLE COMPLEX, INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(i4) :: N,LWORK,info
REAL(r8) :: elapsed_timer
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipiv
DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting complex matrix'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL zgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'ZGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'ZGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(REAL(rwork(1),r8),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'ZGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_timer=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_timer
END IF
END SUBROUTINE lapack_matinv_complex
!---------------------------------------------------------------------------------
!> Invert a general complex matrix in place using LAPACK (single precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_matinv_complex4(nrows,Amat,error)
INTEGER(i4), INTENT(in) :: nrows !< Number of rows/columns
COMPLEX(c4), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(i4) :: N,LWORK,info
REAL(r8) :: elapsed_timer
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipiv
COMPLEX(c4), ALLOCATABLE, DIMENSION(:) :: rwork
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting complex matrix (single precision)'
  CALL stimer%tick
END IF
N = nrows
ALLOCATE(ipiv(N),rwork(1))
CALL cgetrf(N,N,Amat,N,ipiv,info)
IF(info/=0)THEN
  WRITE(*,*)'CGETRF',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL cgetri(N,Amat,N,ipiv,rwork,lwork,info)
IF(info/=0)THEN
  WRITE(*,*)'CGETRI workspace query',info
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=MAX(INT(REAL(rwork(1),r4),4),N)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL cgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
IF(info/=0)WRITE(*,*)'CGETRI',info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_timer=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_timer
END IF
END SUBROUTINE lapack_matinv_complex4
!---------------------------------------------------------------------------------
!> Invert a symmetric positive-definite real matrix in place using a Cholesky
!! factorization (double precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_cholesky_real(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(8), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,info,i,j
REAL(8) :: elapsed_time
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix'
  CALL stimer%tick
END IF
N = nrows
CALL dpotrf('U',N,Amat,N,info)
IF(info/=0)THEN
  error=info
  RETURN
END IF
CALL dpotri('U',N,Amat,N,info)
DO i=1,N
  DO j=1,i-1
    Amat(i,j)=Amat(j,i)
  END DO
END DO
error=info
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_cholesky_real
!---------------------------------------------------------------------------------
!> Invert a symmetric positive-definite real matrix in place using a Cholesky
!! factorization (single precision)
!---------------------------------------------------------------------------------
SUBROUTINE lapack_cholesky_real4(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(r4), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,info,i,j
REAL(8) :: elapsed_time
TYPE(oft_timer) :: stimer
!---Invert Mmat
IF(oft_env%pm)THEN
  WRITE(*,*)'Inverting real matrix (single precision)'
  CALL stimer%tick
END IF
N = nrows
CALL spotrf('U',N,Amat,N,info)
IF(info/=0)THEN
  error=info
  RETURN
END IF
CALL spotri('U',N,Amat,N,info)
DO i=1,N
  DO j=1,i-1
    Amat(i,j)=Amat(j,i)
  END DO
END DO
error=info
IF(oft_env%pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
END SUBROUTINE lapack_cholesky_real4
!---------------------------------------------------------------------------------
!> Native Incomplete LU factorization with no fill
!!
!! Based on ILU++ by Mayer (https://dx.doi.org/10.1002/pamm.200700911)
!!
!! @note Requires CSR graph with sorted column list and diagonal entries
!---------------------------------------------------------------------------------
subroutine ilu0(n,a,lc,kr,alu,jlu,ju,ierr)
integer(4), intent(in) :: n !< Size of matrix
real(8), contiguous, dimension(:), intent(in) :: a !< Matrix entries
real(8), contiguous, dimension(:), intent(inout) :: alu !< Entries for L/U factors
integer(4), contiguous, dimension(:), intent(inout) :: lc !< Column indices
integer(4), contiguous, dimension(:), intent(inout) :: kr !< Row pointer
integer(4), contiguous, dimension(:), intent(inout) :: jlu !< Column indices for alu
integer(4), contiguous, dimension(:), intent(inout) :: ju !< Diagonal indices for alu
integer(4), intent(out) :: ierr !< Error flag
integer(4), allocatable, dimension(:) :: iw
integer(4) :: i,ii,j,jj,jcol,jf,jm,jrow,js,ju0,jw
real(8) :: tl
!
ju0 = n+2
jlu(1) = ju0
! initialize work vector to zero's
ALLOCATE(iw(n))
do i=1, n
    iw(i) = 0
end do
!---
ierr = 0
do ii=1,n
  js=ju0
  !---Generating row ii of L and U
  do j=kr(ii),kr(ii+1)-1
    ! copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
    jcol=lc(j)
    if(jcol==ii)then
      alu(ii)=a(j)
      iw(jcol)=ii
      ju(ii)=ju0
    else
      alu(ju0)=a(j)
      jlu(ju0)=lc(j)
      iw(jcol)=ju0
      ju0=ju0+1
    endif
  end do
  jlu(ii+1)=ju0
  jf=ju0-1
  jm=ju(ii)-1
  !---Move up to diagonal
  do j=js,jm
    jrow=jlu(j)
    tl=alu(j)*alu(jrow)
    alu(j)=tl
    do jj=ju(jrow),jlu(jrow+1)-1
      jw=iw(jlu(jj))
      if(jw/=0)alu(jw)=alu(jw)-tl*alu(jj)
    end do
  end do
  if(alu(ii)==0.0d0)then
    ierr = ii
    exit
  endif
  alu(ii)=1.0d0/alu(ii)
  !---Reset iw pointer
  iw(ii)=0
  do i=js,jf
    iw(jlu(i))=0
  end do
end do
DEALLOCATE(iw)
END SUBROUTINE
!---------------------------------------------------------------------------------
!> Approximate solve using ILU(0) factorization produced by ilu0
!---------------------------------------------------------------------------------
subroutine lusol(n, nrhs, x, alu, jlu, ju)
integer(4), intent(in) :: n !< Size of matrix
integer(4), intent(in) :: nrhs !< Number of RHS to solve
real(8), intent(inout) :: x(n,nrhs) !< Input: RHS, Output: Solution
real(8), contiguous, dimension(:), intent(in) :: alu !< Entries for L/U factors
integer(4), contiguous, dimension(:), intent(in) :: jlu !< Column indices for alu
integer(4), contiguous, dimension(:), intent(in) :: ju !< Diagonal indices for alu
integer(4) :: i,k
! forward solve
do i = 1, n
  do k=jlu(i),ju(i)-1
    x(i,:)=x(i,:)-alu(k)*x(jlu(k),:)
  end do
end do
! backward solve
do i=n,1,-1
  do k=ju(i),jlu(i+1)-1
      x(i,:)=x(i,:)-alu(k)*x(jlu(k),:)
  end do
  x(i,:) = alu(i)*x(i,:)
end do
end subroutine
END MODULE oft_lu
