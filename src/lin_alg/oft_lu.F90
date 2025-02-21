!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_lu.F90
!
!> Native module for direct solves.
!! Available interfaces:
!! - LAPACK (dgetrf, dgetrs)
!! - SuperLU
!! - MUMPS
!! - PARDISO (MKL interface)
!!
!! @note Currently supports local solves only.
!!
!! @authors Chris Hansen
!! @date May 2014
!! @ingroup doxy_oft_lin_alg
!---------------------------------------------------------------------------
MODULE oft_lu
USE, INTRINSIC :: iso_c_binding, ONLY: c_bool, c_int, c_double, c_ptr, c_null_ptr
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
#if !defined(DEF_LU_PACK) && defined(HAVE_MUMPS)
#define DEF_LU_PACK "mumps"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU)
#define DEF_LU_PACK "super"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_SUPERLU_DIST)
#define DEF_LU_PACK "superd"
#endif
#if !defined(DEF_LU_PACK) && defined(HAVE_UMFPACK)
#define DEF_LU_PACK "umfpack"
#endif
#if !defined(DEF_LU_PACK)
#define DEF_LU_PACK "none"
#endif
PRIVATE
PUBLIC lapack_matinv, lapack_cholesky
!---------------------------------------------------------------------------
!> Native ILU(0) preconditioner information object
!---------------------------------------------------------------------------
TYPE :: native_ilu_struc
  INTEGER(i4), POINTER, DIMENSION(:) :: ju => NULL() !< Diagonal entry indices
  INTEGER(i4), POINTER, DIMENSION(:) :: jlu => NULL() !< Column values in iLU factors
  REAL(r8), POINTER, DIMENSION(:) :: alu => NULL() !< Values of iLU factors
END TYPE native_ilu_struc
!---------------------------------------------------------------------------
!> SuperLU solver information
!---------------------------------------------------------------------------
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
END TYPE superlu_struc
!---------------------------------------------------------------------------
!> Pardiso solver information
!---------------------------------------------------------------------------
TYPE :: pardiso_struc
  INTEGER(i4) :: mtype = 11 !< Matrix type (1 -> real/symmetric, 11 -> real/nonsymmetric)
  INTEGER(i4) :: msglvl = 0 !< Controls message output (0 -> no output, 1 -> print some stats)
  INTEGER(i8) :: pt(64) = 0 !< Handle to internal data structure
  INTEGER(i4) :: iparm(64) = 0 !< This array is used to pass various parameters to Intel MKL PARDISO
  INTEGER(i4), POINTER, DIMENSION(:) :: perm => NULL() !< Permutation vector
  TYPE(C_PTR) :: f_factors = C_NULL_PTR !< ILU0 factors pointer
END TYPE pardiso_struc
!---------------------------------------------------------------------------
!> LU solver class
!---------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_lusolver
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  LOGICAL :: update_graph = .TRUE. !< Perform full factorization including symbolic steps
  LOGICAL(kind=c_bool) :: iter_refine = .FALSE. !< Perform iterative refinement
  INTEGER(i4) :: nrhs = 1 !< Number of right hand sides
  INTEGER(i4), POINTER, DIMENSION(:) :: ipiv => NULL() !< Row pivot array ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: atmp => NULL() !< Local dense matrix ("lapack" only)
  REAL(r8), POINTER, DIMENSION(:,:) :: sec_rhs => NULL() !< Storage for additional RHS data
  CHARACTER(LEN=7) :: package = DEF_LU_PACK !< Factorization package
#if defined( HAVE_SUPERLU ) || defined( HAVE_SUPERLU_DIST ) || defined( HAVE_UMFPACK )
  TYPE(superlu_struc) :: superlu_struct !< SuperLU solver information
#endif
#ifdef HAVE_MKL
  TYPE(pardiso_struc) :: pardiso_struct !< MKL/PARDISO solver information
#endif
#ifdef HAVE_MUMPS
  TYPE(DMUMPS_STRUC) :: mumps_struct !< MUMPS solver information
#endif
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
!---------------------------------------------------------------------------
!> ILU solver class
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Interface to dgsisx (iLUT) from SuperLU
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Interface to dgssv from SuperLU-DIST (local only)
!---------------------------------------------------------------------------
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
  TYPE(c_ptr), INTENT(inout) :: grid_handle !< Pointer to SuperLU-DIST communication object
  TYPE(c_ptr), INTENT(inout) :: f_factors !< Pointer to SuperLU-DIST internal data storage
  LOGICAL(c_bool), VALUE, INTENT(in) :: iter_refine !< Perform iterative refinement?
  INTEGER(c_int), VALUE, INTENT(in) :: col_perm !< Column permutation method (0 -> natural ordering,
  !! 1 -> min degree on A'*A, 2 -> min degree on A'+A, 3 -> approx min degree
  !! for unsymmetric matrices)
  INTEGER(c_int), INTENT(inout) :: info !< Solver return status (not currently used)
  END SUBROUTINE oft_superlu_dist_dgssv
!---------------------------------------------------------------------------
!> Interface to dgssv from SuperLU
!---------------------------------------------------------------------------
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
#endif
!---------------------------------------------------------------------------
!> Interface to dgssv from UMFPACK
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Interface to MKL iLU0 solver
!---------------------------------------------------------------------------
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
END INTERFACE lapack_matinv
INTERFACE lapack_cholesky
  MODULE PROCEDURE lapack_cholesky_real
  ! MODULE PROCEDURE lapack_cholesky_complex
END INTERFACE lapack_cholesky
CONTAINS
!---------------------------------------------------------------------------
!> Convert csr graph to csc storage
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Solve a linear system using a direct solve
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE lusolver_apply(self,u,g)
CLASS(oft_lusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: mode,nrhs,ldb,ierr,i,j,k,info
INTEGER(i4), POINTER :: csr_map(:),kr_tmp(:)
REAL(r8), POINTER, DIMENSION(:) :: mat_vals,csc_vals,vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vals,b
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: mytimer
TYPE(oft_graph) :: csr_graph
TYPE(oft_graph) :: csc_graph
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting LU solver: ',self%package,self%refactor
IF(native_matrix_cast(A_native,self%A)<0)CALL oft_abort('Native matrix required', &
  'lusolver_apply',__FILE__)
!---Check call threading for MUMPS
IF(self%package(1:5)=='mumps')THEN
  IF(omp_get_num_threads()>1)CALL oft_abort('MUMPS used by multiple threads', &
    'lusolver_apply',__FILE__)
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
      ALLOCATE(self%superlu_struct%csc_map(A_native%nnz),self%superlu_struct%csc_vals(A_native%nnz))
      csr_graph%nr=A_native%nr; csr_graph%nc=A_native%nc; csr_graph%nnz=A_native%nnz
      csr_graph%kr=>A_native%kr; csr_graph%lc=>A_native%lc
      CALL convert_csr_to_csc(csr_graph,csc_graph,self%superlu_struct%csc_map)
      csc_graph%kr=csc_graph%kr-1; csc_graph%lc=csc_graph%lc-1
      self%superlu_struct%kr=>csc_graph%kr; self%superlu_struct%lc=>csc_graph%lc
      self%superlu_struct%csc_vals=mat_vals(self%superlu_struct%csc_map)
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
      IF(ierr/=0)CALL oft_abort('Factorization failed','lusolver_apply',__FILE__)
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
    IF(ierr/=0)CALL oft_abort('Solve failed','lusolver_apply',__FILE__)
#if !defined( SUPERLU_VER_MAJOR ) || SUPERLU_VER_MAJOR < 5
    !$omp end critical (superlu_solve)
#endif
#else
    CALL oft_abort('OFT not compiled with SuperLU','lusolver_apply',__FILE__)
#endif
CASE("superd")
#ifdef HAVE_SUPERLU_DIST
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)THEN
        DEALLOCATE(self%superlu_struct%csc_map,self%superlu_struct%csc_vals)
        DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
        csr_graph%nr=A_native%nr; csr_graph%nc=A_native%nc; csr_graph%nnz=A_native%nnz
        csr_graph%kr=>A_native%kr; csr_graph%lc=>A_native%lc
        ALLOCATE(self%superlu_struct%csc_map(A_native%nnz))
        ALLOCATE(self%superlu_struct%csc_vals(A_native%nnz))
        CALL convert_csr_to_csc(csr_graph,csc_graph,self%superlu_struct%csc_map)
        csc_graph%kr=csc_graph%kr-1; csc_graph%lc=csc_graph%lc-1
        self%superlu_struct%kr=>csc_graph%kr; self%superlu_struct%lc=>csc_graph%lc
        mode=1
      END IF
      self%superlu_struct%csc_vals=mat_vals(self%superlu_struct%csc_map)
      !!$omp critical (superlu_solve)
      CALL oft_superlu_dist_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        self%superlu_struct%csc_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%grid_handle,self%superlu_struct%f_factors, &
        self%superlu_struct%col_perm,self%iter_refine,ierr)
      !!$omp end critical (superlu_solve)
      IF(ierr/=0)CALL oft_abort('Factorization failed','lusolver_apply',__FILE__)
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
    IF(ierr/=0)CALL oft_abort('Solve failed','lusolver_apply',__FILE__)
#else
    CALL oft_abort('OFT not compiled with SuperLU-DIST','lusolver_apply',__FILE__)
#endif
  CASE("umfpack")
#ifdef HAVE_UMFPACK
    IF(self%refactor)THEN
      mode=3
      IF(self%update_graph)mode=1
      !$omp critical (umfpack_solve)
      CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
        mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
        self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
      IF(ierr/=0)CALL oft_abort('Factorization failed','lusolver_apply',__FILE__)
      !$omp end critical (umfpack_solve)
      self%refactor=.FALSE.
      self%update_graph=.FALSE.
    END IF
    mode=2
    !$omp critical (umfpack_solve)
    CALL oft_umfpack_dgssv(mode,A_native%nr,A_native%nnz,nrhs, &
      mat_vals,self%superlu_struct%lc,self%superlu_struct%kr,vals,ldb, &
      self%superlu_struct%f_factors,self%superlu_struct%col_perm,self%iter_refine,ierr)
    IF(ierr/=0)CALL oft_abort('Solve failed','lusolver_apply',__FILE__)
    !$omp end critical (umfpack_solve)
#else
    CALL oft_abort('OFT not compiled with UMFPACK','lusolver_apply',__FILE__)
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
      IF(ierr/=0)CALL oft_abort('Factorization failed','lusolver_apply',__FILE__)
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
    IF(ierr/=0)CALL oft_abort('Solve failed','lusolver_apply',__FILE__)
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
    self%mumps_struct%rhs => vals
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
      self%refactor=.FALSE.
    END IF
    CALL dgetrs('N',A_native%nr,nrhs,self%atmp,A_native%nr,self%ipiv,vals,ldb,info)
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
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
recursive subroutine lusolver_update(self,new_pattern)
class(oft_lusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))self%update_graph=new_pattern
end subroutine lusolver_update
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!---------------------------------------------------------------------------
subroutine lusolver_setup_xml(self,solver_node,level)
CLASS(oft_lusolver), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
!---
CHARACTER(LEN=7) :: factor_package
CHARACTER(LEN=3) :: fac_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---
current_nodes=>fox_getElementsByTagName(solver_node,"package")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,factor_package,num=nread,iostat=ierr)
  IF(nread==1)THEN
    self%package=factor_package
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'LU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:  ',self%package
END IF
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','lusolver_setup_xml',__FILE__)
#endif
end subroutine lusolver_setup_xml
!---------------------------------------------------------------------------
!> Check for thread safety
!---------------------------------------------------------------------------
recursive function lusolver_check_thread(self) result(thread_safe)
class(oft_lusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
IF(self%package(1:5)=='mumps')thread_safe=.FALSE.
IF(self%package(1:6)=='superd')thread_safe=.FALSE.
end function lusolver_check_thread
!---------------------------------------------------------------------------
!> Destroy direct solver and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine lusolver_delete(self)
class(oft_lusolver), intent(inout) :: self
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
#ifdef HAVE_SUPERLU
  CASE("super")
    mode=4
    CALL oft_superlu_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
#endif
#ifdef HAVE_SUPERLU_DIST
  CASE("superd")
    mode=4
    CALL oft_superlu_dist_dgssv(mode,self%A%nr,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%grid_handle,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    mode=2
    nrhs=1
    CALL oft_superlu_dist_slugrid(mode,self%superlu_struct%comm,nrhs,nrhs,self%superlu_struct%grid_handle,ierr)
    CALL MPI_COMM_FREE(self%superlu_struct%comm, ierr)
    DEALLOCATE(self%superlu_struct%csc_vals,self%superlu_struct%csc_map)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
    self%superlu_struct%comm=MPI_COMM_NULL
#endif
#ifdef HAVE_UMFPACK
  CASE("umfpack")
    mode=4
    CALL oft_umfpack_dgssv(mode,nrhs,nrhs,nrhs,rvals,ivals,ivals, &
      rvals,ldb,self%superlu_struct%f_factors,nrhs,self%iter_refine,ierr)
    DEALLOCATE(self%superlu_struct%kr,self%superlu_struct%lc)
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
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine lusolver_delete
!---------------------------------------------------------------------------
!> Solve a linear system using a direct solve
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE ilusolver_apply(self,u,g)
CLASS(oft_ilusolver), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(i4) :: mode,nrhs,ldb,ierr,i,j,k,info
INTEGER(i4), POINTER :: csr_map(:),kr_tmp(:)
REAL(r8), POINTER, DIMENSION(:) :: mat_vals,csc_vals,vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vals,b
CLASS(oft_native_matrix), POINTER :: A_native
TYPE(oft_timer) :: mytimer
TYPE(oft_graph) :: csr_graph
TYPE(oft_graph) :: csc_graph
DEBUG_STACK_PUSH
IF(TRIM(self%package)=='pardiso')self%package='mkl'
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting ILU solver: ',self%package,self%refactor
IF(native_matrix_cast(A_native,self%A)<0)CALL oft_abort('Native matrix required', &
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
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
recursive subroutine ilusolver_update(self,new_pattern)
class(oft_ilusolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
self%refactor=.TRUE.
IF(PRESENT(new_pattern))self%update_graph=new_pattern
end subroutine ilusolver_update
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!---------------------------------------------------------------------------
subroutine ilusolver_setup_xml(self,solver_node,level)
CLASS(oft_ilusolver), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node !< XML node containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
!---
CHARACTER(LEN=7) :: factor_package
CHARACTER(LEN=3) :: fac_type
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---
current_nodes=>fox_getElementsByTagName(solver_node,"package")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,factor_package,num=nread,iostat=ierr)
  IF(nread==1)THEN
    self%package=factor_package
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'LU solver setup:'
  WRITE(*,'(2X,2A)')'- Package:  ',self%package
END IF
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','lusolver_setup_xml',__FILE__)
#endif
end subroutine ilusolver_setup_xml
!---------------------------------------------------------------------------
!> Check for thread safety
!---------------------------------------------------------------------------
recursive function ilusolver_check_thread(self) result(thread_safe)
class(oft_ilusolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
! IF(self%package(1:5)=='mumps')thread_safe=.FALSE.
! IF(self%package(1:6)=='superd')thread_safe=.FALSE.
end function ilusolver_check_thread
!---------------------------------------------------------------------------
!> Destroy direct solver and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine ilusolver_delete(self)
class(oft_ilusolver), intent(inout) :: self
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
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
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
lwork=INT(rwork(1),4)
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
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
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
  error=info
  DEALLOCATE(ipiv,rwork)
  RETURN
END IF
lwork=-1
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
lwork=INT(rwork(1),4)
IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(rwork)
ALLOCATE(rwork(lwork))
CALL zgetri(N,Amat,N,ipiv,rwork,lwork,info)
error=info
DEALLOCATE(ipiv,rwork)
IF(oft_env%pm)THEN
  elapsed_timer=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_timer
END IF
END SUBROUTINE lapack_matinv_complex
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE lapack_cholesky_real(nrows,Amat,error)
INTEGER(4), INTENT(in) :: nrows !< Number of rows/columns
REAL(8), INTENT(inout) :: Amat(nrows,nrows) !< Matrix to invert
INTEGER(4), INTENT(out) :: error !< Error flag
!---
INTEGER(4) :: N,LWORK,info,i,j
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
!------------------------------------------------------------------------------
!> Native Incomplete LU factorization with no fill
!!
!! Based on ILU++ by Mayer (https://dx.doi.org/10.1002/pamm.200700911)
!!
!! @note Requires CSR graph with sorted column list and diagonal entries
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Approximate solve using ILU(0) factorization produced by ilu0
!------------------------------------------------------------------------------
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
