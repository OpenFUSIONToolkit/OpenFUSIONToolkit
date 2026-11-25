!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_native_solvers.F90
!
!> Abstract solver interfaces and select native implementations
!!
!! Abstract interface definitions
!! - Abstract solver class
!! - Abstract eigensolver class
!! - Abstract orthogonalization class
!! - Field creation prototype
!! - Boundary condition prototype
!!
!! Native solver implementations
!! - Newton iteration
!!
!! Preconditioner implementations
!! - Diagonal scaling
!! - Symmetric Point-Jacobi (Richardson iteration with diagonal scaling)
!! - Block-Jacobi
!! - Multi-Grid
!!
!! @sa oft_cg, oft_gmres, oft_petsc_solvers
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_native_la
!------------------------------------------------------------------------------
MODULE oft_native_solvers
USE oft_local
USE oft_base
USE oft_stitching, ONLY: oft_seam
USE oft_scatter, ONLY: oft_scatter_plan, scatter_plan_setup_local, scatter_plan_delete, scatter_local, &
  OFT_SCATTER_INSERT, OFT_SCATTER_ADD, oft_index_set, index_set_create, index_set_grow, index_set_delete, &
  oft_exchange_report, oft_exchange_nsent
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_cvector_ptr, &
  oft_matrix, oft_cmatrix, oft_graph, oft_ml_vecspace, OFT_ACCESS_READ, OFT_ACCESS_WRITE, &
  OFT_ACCESS_READWRITE, oft_overlap_dofs, overlap_dofs_delete, OFT_CELL_ADJ_POINT, OFT_CELL_ADJ_EDGE, &
  OFT_CELL_ADJ_FACE
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, &
  oft_native_matrix, native_matrix_cast, oft_native_submatrix, &
  partition_graph, native_matrix_setup_full, oft_native_dense_cmatrix
USE oft_native_overlap_mod, ONLY: oft_native_overlap, native_overlap_setup, native_overlap_update, &
  native_overlap_gather, native_overlap_scatter, native_overlap_delete, native_overlap_setup_dofs, &
  native_overlap_index, native_overlap_pu_setup
USE oft_native_kernels, ONLY: nk_jacobi_start, nk_jacobi_sweep, nk_jacobi_restart, nk_copy, nk_axpby
USE oft_lu, ONLY: lapack_matinv
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc, oft_solver_ptr, &
  oft_csolver, oft_csolver_bc, oft_eigsolver, solver_setup, csolver_setup, &
  eigsolver_setup, solver_delete, csolver_delete
IMPLICIT NONE
#include "local.h"
private
!> Re-orthogonalize GMRES Arnoldi vectors when orthogonalization reduces the norm
!! by more than this factor ("twice is enough" criterion)
REAL(r8), PARAMETER :: gmres_reorth_tol = 1.d0/SQRT(2.d0)
!> Solve the upper-triangular GMRES least-squares system
INTERFACE gmres_backsolve
  MODULE PROCEDURE gmres_backsolve_r8
  MODULE PROCEDURE gmres_backsolve_c8
END INTERFACE gmres_backsolve
PUBLIC gmres_backsolve
!------------------------------------------------------------------------------
!> CG solver class
!!
!! @note Matrix must be SPD, otherwise solver will fail.
!------------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_native_cg_solver
  logical :: precond=.FALSE. !< Solver as preconditioner
  class(oft_solver_bc), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_solver_bc), pointer :: cleaner => NULL() !< Null-Space cleaner
contains
  !> Solve system
  procedure :: apply => cg_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => cg_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => cg_delete
end type oft_native_cg_solver
!------------------------------------------------------------------------------
!> GMRES solver class
!------------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_native_gmres_solver
  logical :: precond = .FALSE. !< Flag to indicate use as a preconditioner (not used)
  integer(i4) :: nrits=10 !< Convergence specification
  integer(i4) :: min_nits=0 !< Minimum inner iterations per cycle before the \f$ |r|/|u| \f$ check (`its=-1,-2` only)
  class(oft_solver_bc), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_solver_bc), pointer :: cleaner => NULL() !< Null-Space cleaner
  class(oft_vector), private, pointer :: r => NULL() !< Temporary storage vector
  class(oft_vector), private, pointer :: w => NULL() !< Temporary storage vector
  type(oft_vector_ptr), private, pointer, dimension(:) :: v => NULL() !< Search directions
  type(oft_vector_ptr), private, pointer, dimension(:) :: z => NULL() !< Preconditioned directions
contains
  !> Solve system
  procedure :: apply => gmres_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => gmres_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => gmres_delete
end type oft_native_gmres_solver
!------------------------------------------------------------------------------
!> GMRES solver class (complex)
!------------------------------------------------------------------------------
type, public, extends(oft_csolver) :: oft_native_gmres_csolver
  logical :: precond = .FALSE. !< Flag to indicate use as a preconditioner (not used)
  integer(i4) :: nrits=10 !< Convergence specification
  integer(i4) :: min_nits=0 !< Minimum inner iterations per cycle before the \f$ |r|/|u| \f$ check (`its=-1,-2` only)
  class(oft_csolver_bc), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_csolver_bc), pointer :: cleaner => NULL() !< Null-Space cleaner
  class(oft_cvector), private, pointer :: r => NULL() !< Temporary storage vector
  class(oft_cvector), private, pointer :: w => NULL() !< Temporary storage vector
  type(oft_cvector_ptr), private, pointer, dimension(:) :: v => NULL() !< Search directions
  type(oft_cvector_ptr), private, pointer, dimension(:) :: z => NULL() !< Preconditioned directions
contains
  !> Solve system
  procedure :: apply => cgmres_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => cgmres_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => cgmres_delete
end type oft_native_gmres_csolver
!------------------------------------------------------------------------------
!> CG eigensolver class
!!
!! @note Matrix system must be SPD, otherwise solver will fail.
!------------------------------------------------------------------------------
type, public, extends(oft_eigsolver) :: oft_native_cg_eigsolver
  integer(i4) :: nrestarts = 2
  integer(i4) :: ninner = -1
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  class(oft_solver_bc), pointer :: orthog => NULL() !< Orthogonalization
contains
  !> Solve system
  procedure :: apply => cg_eigsolver_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => cg_eig_delete
end type oft_native_cg_eigsolver
!------------------------------------------------------------------------------
!> Native Newton solver
!!
!! Non-linear system solver based on Newton iteration for the system,
!! \f$ M_j(x_i) = y_j \f$, where the Jacobian is defined as
!! \f$ J_{ij} = \frac{\partial M_j}{\partial x_i} \f$. The non-linear update then
!! takes the form \f$ x^{n+1} = x^n - J^{-1}(x^n) \cdot M(x^n) \f$.
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_nksolver
  LOGICAL :: backtrack = .TRUE. !< Perform backtracking?
  INTEGER(i4) :: its = -1 !< Maximum iteration count
  INTEGER(i4) :: cits = 0 !< Number of iteractions to convergence
  INTEGER(i4) :: lits = 0 !< Number of linear solver iteractions
  INTEGER(i4) :: nlits = 0 !< Number of non-linear steps to convergence
  INTEGER(i4) :: up_freq = 1 !< Frequency of jacobian updates
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_vector), POINTER :: du => NULL() !< Storage for non-linear corrections
  CLASS(oft_matrix), POINTER :: A => NULL() !< Metric matrix
  CLASS(oft_solver), POINTER :: J_inv => NULL() !< Jacobian inverse solver
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  PROCEDURE(oft_update_jacobian), POINTER, NOPASS :: J_update => NULL() !< Jacobian update subroutine
CONTAINS
  !> Solve non-linear system
  PROCEDURE :: apply => nksolver_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => nksolver_delete
END TYPE oft_nksolver
!------------------------------------------------------------------------------
!> Identity matrix inversion
!------------------------------------------------------------------------------
TYPE, PUBLIC, extends(oft_solver) :: oft_identity_inv
CONTAINS
  !> Apply preconditioner
  PROCEDURE :: apply => identity_inv_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => identity_inv_delete
END TYPE oft_identity_inv
!------------------------------------------------------------------------------
!> Diagonal preconditioner
!------------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_diag_scale
contains
  !> Apply preconditioner
  procedure :: apply => diag_scale_apply
  !> Clean-up internal smoother storage
  procedure :: delete => diag_scale_delete
end type oft_diag_scale
!------------------------------------------------------------------------------
!> Diagonal preconditioner
!------------------------------------------------------------------------------
type, public, extends(oft_csolver) :: oft_diag_cscale
contains
  !> Apply preconditioner
  procedure :: apply => cdiag_scale_apply
  !> Clean-up internal smoother storage
  procedure :: delete => cdiag_scale_delete
end type oft_diag_cscale
!------------------------------------------------------------------------------
!> Symmetric Jacobi smoother
!------------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_jblock_precond
  logical, private :: down=.TRUE. !< Internal flag for smoother step
  logical, private :: warn_once = .TRUE. !< Internal flag for iteration check
  real(r8) :: df !< Damping factor
  real(r8), private, pointer, dimension(:) :: u_save => NULL() !< Internal storage
  real(r8), private, pointer, dimension(:) :: g_save => NULL() !< Internal storage
  class(oft_native_vector), private, pointer :: p => NULL() !< Internal storage
  class(oft_native_vector), private, pointer :: q => NULL() !< Internal storage
contains
  !> Apply 1-Step of the symmetric jacobi smoother
  procedure :: apply => jblock_precond_apply
  !> Setup solver from XML node
  procedure :: setup_from_xml => jblock_setup_xml
  !> Clean-up internal smoother storage
  procedure :: delete => jblock_precond_delete
end type oft_jblock_precond
!------------------------------------------------------------------------------
!> Multi-level preconditioner level
!!
!! @note This class acts as a driver/container for separate solver objects
!! as a result any linear algebra backend may be used with this type.
!------------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_ml_precond
  logical :: symmetric = .FALSE. !< Symmetric flag
  integer(i4) :: report_err = 0 !< Report errors after each cycle (0: no, 1: total, 2: block-wise)
  integer(i4) :: level !< Current level in ML context
  integer(i4) :: minlevel !< Lowest level
  real(r8) :: timings(4) = 0.d0 !< Timing for each stage
  class(oft_solver), pointer :: smooth_up => NULL() !< Smoother for up-cycle
  class(oft_solver), pointer :: smooth_down => NULL() !< Smoother for down-cycle
  class(oft_solver), pointer :: base_solve => NULL() !< Pointer to lower level solve
  class(oft_vector), pointer :: p => NULL() !< Temporary vector
  class(oft_vector), pointer :: r => NULL() !< Temporary vector
  class(oft_vector), pointer :: ucors => NULL() !< Temporary coarse vector
  class(oft_vector), pointer :: gcors => NULL() !< Temporary coarse vector
  class(oft_ml_vecspace), pointer :: ml_vecspace => NULL() !< Multi-level vector space
contains
  !> Apply smoother
  procedure :: apply => ml_precond_apply
  !> Update children with new settings/operators
  procedure :: update => ml_precond_update
  !> Print solver information
  PROCEDURE :: view => ml_precond_view
  !> Clean-up internal smoother storage
  procedure :: delete => ml_precond_delete
end type oft_ml_precond
!------------------------------------------------------------------------------
!> Multi-level transfer level
!------------------------------------------------------------------------------
type, public, extends(oft_ml_precond) :: oft_ml_trans
contains
  !> Transfer field between Base and MPI levels
  procedure :: apply => ml_trans_apply
  !> Clean-up internal smoother storage
  procedure :: delete => ml_trans_delete
end type oft_ml_trans
!------------------------------------------------------------------------------
!> Block-Jacobi preconditioner
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_bjprecond
  LOGICAL :: update_slice = .FALSE. !< Update local matrices on next application
  LOGICAL :: boundary_overlap = .FALSE. !< Use Additive-Schwarz with boundary elements
  LOGICAL :: loverlap = .FALSE. !< Use Additive-Schwarz with interior boundaries
  INTEGER(i4) :: nlocal = 1 !< Number of subdomains on each processor
  INTEGER(i4) :: slice_group(10) = -1 !< Slice groupings
  INTEGER(i4), POINTER, DIMENSION(:) :: part => NULL()
  TYPE(oft_native_submatrix), pointer, dimension(:) :: alocals => NULL() !< Local matrix block
  TYPE(oft_seam), pointer, dimension(:) :: local_stitch => NULL() !< Local vectors
  TYPE(oft_native_vector), POINTER, DIMENSION(:) :: ulocs => NULL() !< Local solution vector for each block
  TYPE(oft_native_vector), POINTER, DIMENSION(:) :: glocs => NULL() !< Local RHS vector for each block
  TYPE(oft_solver_ptr), pointer, dimension(:) :: solvers => NULL() !< Local block solvers
  TYPE(oft_1d_int), pointer, dimension(:) :: parts => NULL() !< Local partition indices
  !> Gather/scatter plan of each block (rows `parts(i)%v`; with overlap, weights 1 on
  !! rows assigned to the block and 0 on overlap rows, which are read but not written back,
  !! unless `additive`; partition of unity weights if `additive` and `weighted`)
  TYPE(oft_scatter_plan), POINTER, DIMENSION(:) :: bplans => NULL()
  !> Number of overlap layers added to each block within a process (-1: one layer if
  !! `loverlap`, none otherwise). Set by @ref oft_schwarz_precond.
  INTEGER(i4) :: overlap_levels = -1
  !> Combine overlapping blocks additively (sum of block solutions) instead of
  !! restricting the write-back to the rows assigned to each block. Set by
  !! @ref oft_schwarz_precond.
  LOGICAL :: additive = .FALSE.
  !> Weight the additive combination by a partition of unity: each block solution is
  !! scaled on each row by one over the number of blocks containing the row. Only used
  !! with `additive`. Set by @ref oft_schwarz_precond.
  LOGICAL :: weighted = .FALSE.
  INTEGER(i4), POINTER, DIMENSION(:) :: lunion => NULL() !< Rows covered by any block (additive combination)
CONTAINS
  !> Solve system
  PROCEDURE :: apply => bjprecond_apply
  !> Update solver with new settings/operators
  PROCEDURE :: update => bjprecond_update
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => bjprecond_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => bjprecond_delete
END TYPE oft_bjprecond
!------------------------------------------------------------------------------
!> Overlapping Schwarz preconditioner (restricted, additive or weighted additive)
!!
!! Subdomain problems are solved with the sub-solver `pre` (as for
!! @ref oft_bjprecond) on subdomains grown from a non-overlapping partition by
!! `overlap` layers along the matrix graph, and combined either restricted (each row
!! is taken from the subdomain it was assigned to; default), additively (sum of
!! the subdomain solutions, symmetric for symmetric sub-solves) or weighted additively
!! (`additive` and `weighted`: sum of the subdomain solutions scaled on each row by
!! one over the number of subdomains containing the row, a partition of unity).
!! - `nlocal=1` (default): one subdomain per process made of its owned rows, grown
!!   across processes (rows held only by other processes are transferred, see
!!   @ref oft_native_overlap_mod::native_overlap_setup "native_overlap_setup").
!!   `overlap=0` is block-Jacobi on the owned rows.
!! - `nlocal>1` (or negative for slice groups, as for @ref oft_bjprecond): blocks
!!   within each process, grown within the owned rows of the process (`overlap=1`
!!   restricted is the `loverlap` option of @ref oft_bjprecond).
!!
!! With `cell_mode` the blocks are instead grown by `overlap` layers of mesh cells
!! (adjacent by `adjacency`) by the cell overlap builder attached to the matrix
!! (`oft_matrix::cell_ovl`), which also sets the number of blocks per process. Blocks
!! may cross process boundaries for any number of blocks; their union on each process
!! is one cross-process subdomain (@ref oft_native_overlap_mod::native_overlap_setup_dofs
!! "native_overlap_setup_dofs") and each block is a submatrix of it (see
!! @ref schwarz_cell_apply). In this mode the block-Jacobi storage is reused per block:
!! `parts(i)%v` holds the (sorted) subdomain rows of block `i`, `alocals(i)` its
!! submatrix of `ovl%amat`, `bplans(i)` its plan between the subdomain and block
!! vectors (`ulocs(i)`, `glocs(i)`), `solvers(i)` its sub-solver and `nlocal` the
!! number of blocks.
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_bjprecond) :: oft_schwarz_precond
  INTEGER(i4) :: overlap = 1 !< Number of overlap layers
  LOGICAL :: cell_mode = .FALSE. !< Grow blocks by layers of cells (requires `A%cell_ovl`) instead of the matrix graph
  INTEGER(i4) :: adjacency = OFT_CELL_ADJ_POINT !< Cell adjacency used to grow blocks in cell mode (`OFT_CELL_ADJ_*`)
  INTEGER(i4) :: base_layers = -1 !< Base mesh layers bounding the halo in cell mode (<0: builder default)
  TYPE(oft_native_overlap), POINTER :: ovl => NULL() !< Subdomain across processes (`nlocal=1`, `overlap>0`; union of blocks in cell mode)
  TYPE(oft_seam), POINTER :: sub_seam => NULL() !< Seam of subdomain vectors (no stitching)
  TYPE(oft_native_vector), POINTER :: usub => NULL() !< Subdomain solution vector (initial guess in cell mode)
  TYPE(oft_native_vector), POINTER :: gsub => NULL() !< Subdomain right-hand side vector
  TYPE(oft_native_vector), POINTER :: usol => NULL() !< Combined block solutions on the subdomain (cell mode)
  TYPE(oft_overlap_dofs), POINTER :: cdofs => NULL() !< Block DoF sets from the cell overlap builder (cell mode)
CONTAINS
  !> Solve system
  PROCEDURE :: apply => schwarz_apply
  !> Update solver with new settings/operators
  PROCEDURE :: update => schwarz_update
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => schwarz_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => schwarz_delete
END TYPE oft_schwarz_precond
!------------------------------------------------------------------------------
!> Inverted diagonal block for @ref oft_native_dense_cbjprecond
!------------------------------------------------------------------------------
TYPE :: dense_cbjprecond_block
  INTEGER(i4) :: n = 0 !< Number of rows in block
  INTEGER(i4), POINTER, DIMENSION(:) :: ind => NULL() !< Row indices (negative for overlap rows)
  COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: inv => NULL() !< Inverse of diagonal block
END TYPE dense_cbjprecond_block
!------------------------------------------------------------------------------
!> Block-Jacobi preconditioner for native dense complex matrices
!!
!! Blocks are defined by partitioning a supplied graph (eg. the non-zero pattern
!! of a related sparse matrix) using METIS. Each diagonal block is extracted from
!! the dense matrix and inverted explicitly during @ref dense_cbjprecond_setup.
!! If `loverlap` is set each block is extended by one level of graph neighbors and
!! only rows owned by the block are updated (restricted Additive-Schwarz).
!!
!! @note Only supports local (non-distributed) matrices. For diagonal scaling use
!! @ref oft_diag_cscale instead.
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_csolver) :: oft_native_dense_cbjprecond
  LOGICAL :: loverlap = .FALSE. !< Extend blocks by one level of graph neighbors?
  INTEGER(i4) :: nlocal = 1 !< Number of blocks (must be >=1)
  INTEGER(i4), PRIVATE :: max_block_size = 0 !< Size of largest block
  TYPE(dense_cbjprecond_block), PRIVATE, POINTER, DIMENSION(:) :: blocks => NULL() !< Inverted diagonal blocks
CONTAINS
  !> Partition graph and factor diagonal blocks
  PROCEDURE :: setup => dense_cbjprecond_setup
  !> Apply preconditioner
  PROCEDURE :: apply => dense_cbjprecond_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => dense_cbjprecond_delete
END TYPE oft_native_dense_cbjprecond
!---Interfaces
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Needs docs
  !------------------------------------------------------------------------------
  SUBROUTINE oft_update_jacobian(a)
    IMPORT oft_vector
    CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Needs docs
  END SUBROUTINE oft_update_jacobian
END INTERFACE
!---Make classes and prototypes public
PUBLIC native_cg_solver_cast, native_gmres_solver_cast
PUBLIC jblock_precond_cast, ml_precond_cast
PUBLIC ml_trans_cast, diag_scale_cast!, oft_veccreate_proto, oft_interp_proto
CONTAINS
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_native_cg_solver
!!
!! The source matrix must be @ref oft_native_cg_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION native_cg_solver_cast(self,source) result(success)
type(oft_native_cg_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cg_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
END FUNCTION native_cg_solver_cast
!------------------------------------------------------------------------------
!> Solve a linear system using the Conjugate-Gradient method
!------------------------------------------------------------------------------
recursive subroutine cg_solver_apply(self,u,g)
class(oft_native_cg_solver), intent(inout) :: self !< Solver object
CLASS(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
logical :: pm_save,breakdown
integer(i4) :: k,its,max_its
real(r8) :: alfa,beta,uu,gg,ggold,f,fold,pdq,ggin,ggnew
CLASS(oft_matrix), pointer :: A
class(oft_vector), pointer :: p,q,s,r
DEBUG_STACK_PUSH
A=>self%A
its=self%its
self%cits=0
IF(.NOT.self%initialized)CALL solver_setup(self)
if((oft_env%pm.AND.oft_env%head_proc))write(*,'(A)')'Starting CG solver'
if(u%n/=g%n)call oft_abort('Vectors are not the same length','cgsolver_apply',__FILE__)
call u%new(p)
call u%new(q)
if(self%precond)then
  call u%new(r)
  CALL r%copy(g)
end if
if(associated(self%pre))call u%new(s)
if(associated(self%bc))call self%bc%apply(g)
if(associated(self%orthog))call self%orthog%apply(u)
uu = clamp_nonneg(u%dot(u))
if(uu==0.d0)then
  alfa=1.d0; f=0.d0
  call q%set(0.d0)
else
  call A%apply(u,q)
  if(associated(self%bc))call self%bc%apply(q)
  pdq = u%dot(q)
  gg = g%dot(u)
  alfa=1.d0; f=pdq/2.d0-gg
endif
call g%add(1.d0,-1.d0,q)
CALL p%copy(g)
if(associated(self%pre))then
  pm_save=oft_env%pm; oft_env%pm=self%pre%pm
  call self%pre%apply(s,p)
  oft_env%pm=pm_save
  if(associated(self%bc))call self%bc%apply(s)
  CALL p%copy(s)
end if
gg = g%dot(p)
ggin = clamp_nonneg(g%dot(g)) ! Clamp round-off (see GMRES Arnoldi step)
100 FORMAT (I6,3ES14.6)
110 FORMAT (I6,4ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,f,SQRT(uu),SQRT(ggin)
if(its==0.or.gg==0.d0)then
  call p%delete
  call q%delete
  deallocate(p,q)
  if(self%precond)then
    CALL g%copy(r)
    CALL r%delete
    DEALLOCATE(r)
  end if
  if(associated(self%pre))then
    call s%delete
    deallocate(s)
  end if
  DEBUG_STACK_POP
  return
endif
!---Begin CG iteration
max_its=INT(MIN(1000000_i8,u%ng+abs(its)),4)
breakdown=.FALSE.
do k=1,max_its
  if(associated(self%orthog))call self%orthog%apply(p)
  call A%apply(p,q)
  pdq = p%dot(q)
  if(pdq==0.d0)then ! Breakdown (zero search direction)
    breakdown=.TRUE.
    exit
  end if
  if(associated(self%bc))call self%bc%apply(q)
  alfa=gg/pdq
  fold=f; f=f-alfa*gg/2
  call u%add(1.d0,alfa,p)
  call g%add(1.d0,-alfa,q)
  if(associated(self%cleaner))call self%cleaner%apply(g)
  CALL q%copy(g)
  if(associated(self%pre))then
    pm_save=oft_env%pm; oft_env%pm=self%pre%pm
    call self%pre%apply(s,q)
    oft_env%pm=pm_save
    if(associated(self%bc))call self%bc%apply(s)
    CALL q%copy(s)
  end if
  if(associated(self%orthog))call self%orthog%apply(q)
  ggold=gg
  gg=g%dot(q)
  ! Without preconditioning or orthogonalization q==g, so skip the redundant reduction
  if(associated(self%pre).OR.associated(self%orthog))then
    ggnew = clamp_nonneg(g%dot(g))
  else
    ggnew = clamp_nonneg(gg)
  end if
  beta=gg/ggold
  if((k<=self%itplot.OR.MOD(k,self%itplot)==0).AND.oft_env%pm)THEN
    uu=clamp_nonneg(u%dot(u)) ! Only needed for output
    if(oft_env%head_proc)write(*,110)k,f,SQRT(uu),SQRT(ggnew),SQRT(ggnew/uu)
  END IF
  if(k==its)exit
  if(its==-1.and.REAL(f,4)>=REAL(fold,4))exit
  if(its==-2.and.f>=fold)exit
  if(SQRT(ggnew)<self%atol)exit
  if(SQRT(ggnew/ggin)<self%rtol)exit
  call p%add(beta,1.d0,q)
enddo
call p%delete
call q%delete
DEALLOCATE(p,q)
if(self%precond)then
  CALL g%copy(r)
  CALL r%delete
  DEALLOCATE(r)
end if
if(associated(self%pre))then
  call s%delete
  DEALLOCATE(s)
end if
if(k>=max_its.AND.oft_env%head_proc)write(*,*)'Full Its: cgsolver_apply'
self%cits=k
if(breakdown)then
  if(oft_env%head_proc)call oft_warn('CG breakdown (p.Ap == 0) in cg_solver_apply')
  self%cits=-k
end if
DEBUG_STACK_POP
end subroutine cg_solver_apply
!------------------------------------------------------------------------------
!> Setup CG solver from XML definition
!------------------------------------------------------------------------------
subroutine cg_setup_xml(self,solver_node,level)
CLASS(oft_native_cg_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_int
REAL(r8), POINTER, DIMENSION(:) :: read_real
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `its` node","cg_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `its` values specified","cg_setup_xml",__FILE__)
      self%its=read_int(val_level)
    ELSE
      self%its=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"atol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `atol` node","cg_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `atol` values specified","cg_setup_xml",__FILE__)
      self%atol=read_real(val_level)
    ELSE
      self%atol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"rtol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `rtol` node","cg_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `rtol` values specified","cg_setup_xml",__FILE__)
      self%rtol=read_real(val_level)
    ELSE
      self%rtol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,*)'CG solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Tolerance:   ',REAL(self%atol,4),REAL(self%rtol,4)
END IF
DEBUG_STACK_POP
end subroutine cg_setup_xml
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine cg_delete(self,propogate)
class(oft_native_cg_solver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
NULLIFY(self%A)
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
end subroutine cg_delete
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_native_gmres_solver
!!
!! The source matrix must be @ref oft_native_gmres_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION native_gmres_solver_cast(self,source) result(success)
type(oft_native_gmres_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_gmres_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION native_gmres_solver_cast
!------------------------------------------------------------------------------
!> Solve a linear system using the flexible GMRES method
!------------------------------------------------------------------------------
recursive subroutine gmres_solver_apply(self,u,g)
class(oft_native_gmres_solver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_matrix), pointer :: A
logical :: pm_save,pre_full_res,inner_conv
integer(i4) :: i,j,k,nits,its,nrits
real(r8) :: uu,uuold,gg,ggold,delta,hkmi,ggin,wnorm0,res_u_tol
real(r8), allocatable :: h(:,:),c(:),s(:),res(:),hdots(:),znorm(:)
class(oft_vector), pointer :: r,w
type(oft_vector_ptr), pointer, dimension(:) :: v,z
type(oft_vector_ptr), allocatable, dimension(:) :: vw
DEBUG_STACK_PUSH
!---
A=>self%A
its=self%its
nrits=self%nrits
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%v(nrits+1),self%z(nrits+1))
  DO i=1,nrits+1
    CALL u%new(self%v(i)%f)
    CALL u%new(self%z(i)%f)
  END DO
  CALL u%new(self%r)
  CALL u%new(self%w)
  CALL solver_setup(self)
END IF
if((oft_env%pm.AND.oft_env%head_proc))write(*,'(A)')'Starting GMRES solver'
if(u%n/=g%n)call oft_abort('Vectors are not the same length','gmres_solver_apply',__FILE__)
!
if(associated(self%bc))call self%bc%apply(g)
if(associated(self%cleaner))call self%cleaner%apply(g)
if(associated(self%orthog))call self%orthog%apply(u)
!
v=>self%v
z=>self%z
r=>self%r
w=>self%w
CALL w%set(0.d0)
allocate(c(nrits+1),s(nrits+1),res(nrits+1),h(nrits+1,nrits))
allocate(hdots(nrits+1),vw(nrits+1),znorm(nrits))
c=0.d0; s=0.d0; res=0.d0; h=0.d0
!
uu = clamp_nonneg(u%dot(u))
IF(uu>0.d0)call A%apply(u,w)
if(associated(self%bc))call self%bc%apply(w)
!if(associated(self%cleaner))call self%cleaner%apply(g)
call r%add(0.d0,1.d0,g,-1.d0,w)
!if(associated(self%cleaner))call self%cleaner%apply(r)
gg = clamp_nonneg(r%dot(r)) ! Clamp round-off (see Arnoldi step)
ggin=gg
100 FORMAT (I6,2ES14.6)
110 FORMAT (I6,3ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,SQRT(uu),SQRT(gg)
! Tolerance on |r|/|u| for its=-1/-2 (distinct from self%rtol, which is relative to |r0|)
res_u_tol=-1.d0
IF(its==-1)res_u_tol=1.d-7
IF(its==-2)res_u_tol=1.d-14
nits=0
if(self%its/=0.AND.self%nrits/=0)THEN
do j=1,INT(MIN(u%ng,INT(HUGE(j),i8)),4)
  IF(gg==0.d0)EXIT
  res=0.d0
  res(1)=sqrt(gg)
  call v(1)%f%add(0.d0,1.d0/res(1),r)
  inner_conv=.FALSE.
  do i=1,nrits
    ! Precondition search direction
    if(associated(self%pre))then
      CALL w%copy(v(i)%f)
      CALL z(i)%f%set(0.d0)
      pm_save=oft_env%pm; oft_env%pm=self%pre%pm
      pre_full_res=self%pre%full_residual
      self%pre%full_residual=.FALSE.
      call self%pre%apply(z(i)%f,w)
      self%pre%full_residual=pre_full_res
      oft_env%pm=pm_save
    else
      CALL z(i)%f%copy(v(i)%f)
    end if
    if(associated(self%bc))call self%bc%apply(z(i)%f)
    if(associated(self%orthog))call self%orthog%apply(z(i)%f)
    IF(res_u_tol>0.d0)znorm(i)=SQRT(clamp_nonneg(z(i)%f%dot(z(i)%f)))
    call A%apply(z(i)%f,w)
    if(associated(self%bc))call self%bc%apply(w)
    !---Arnoldi iteration (CGS with selective re-orthogonalization)
    do k=1,i
      vw(k)%f=>v(k)%f
    end do
    vw(i+1)%f=>w
    hdots(1:i+1)=w%mdot(vw(1:i+1),i+1) ! Fuse projections and initial norm into one reduction
    h(1:i,i)=hdots(1:i)
    wnorm0=SQRT(clamp_nonneg(hdots(i+1)))
    CALL w%maxpy(-h(1:i,i),v(1:i),i)
    !if(associated(self%cleaner))call self%cleaner%apply(w)
    ! Safety net: after near-complete cancellation (eg. an exact preconditioner) clamp any
    ! negative round-off in the self dot-product to zero, so this is treated as a (lucky)
    ! breakdown instead of producing NaN. NaN itself is preserved so failures still propagate.
    h(i+1,i)=SQRT(clamp_nonneg(w%dot(w)))
    IF(h(i+1,i)<gmres_reorth_tol*wnorm0)THEN ! Significant cancellation, re-orthogonalize
      hdots(1:i)=w%mdot(v(1:i),i)
      CALL w%maxpy(-hdots(1:i),v(1:i),i)
      h(1:i,i)=h(1:i,i)+hdots(1:i)
      h(i+1,i)=SQRT(clamp_nonneg(w%dot(w)))
    END IF
    IF(h(i+1,i)>0.d0)call v(i+1)%f%add(0.d0,1.d0/h(i+1,i),w)
    !if(associated(self%cleaner))call self%cleaner%apply(v(i+1)%f)
    !---Apply Givens rotation
    do k=2,i
      hkmi = h(k-1,i)
      h(k-1,i) = c(k-1)*hkmi + s(k-1)*h(k,i)
      h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
    enddo
    delta = SQRT( h(i,i)**2 + h(i+1,i)**2 )
    c(i) = h(i,i) / delta
    s(i) = h(i+1,i) / delta
    h(i,i) = c(i)*h(i,i) + s(i)*h(i+1,i)
    !---Update residual
    res(i+1) = -s(i)*res(i)
    res(i) = c(i)*res(i)
    nits=nits+1
    IF(i==nrits)EXIT
    IF(ABS(res(i+1))==0.d0)EXIT ! Lucky breakdown
    IF(nits==its)EXIT
    IF(ABS(res(i+1))<self%atol)EXIT
    IF(ABS(res(i+1))/SQRT(ggin)<self%rtol)EXIT
    !---Check |r|/|u| criterion for current iterate (after `min_nits` inner iterations
    !   to avoid stalling through many short restart cycles)
    IF(res_u_tol>0.d0.AND.i>=self%min_nits)THEN
      CALL gmres_backsolve(i,h,res,hdots)
      ! Necessary condition using |u0 + Z*y| <= |u0| + sum_k |y_k|*|z_k|
      IF(ABS(res(i+1))<res_u_tol*(SQRT(uu)+SUM(ABS(hdots(1:i))*znorm(1:i))))THEN
        CALL w%copy(u) ! Form iterate in w (used as scratch)
        CALL w%maxpy(hdots(1:i),z(1:i),i)
        inner_conv=(ABS(res(i+1))<res_u_tol*SQRT(clamp_nonneg(w%dot(w))))
        IF(inner_conv)EXIT
      END IF
    END IF
    IF(oft_env%head_proc.AND.oft_env%pm.AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
      WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
    END IF
  end do
  k=i
  ! Update iterate
  IF(inner_conv)THEN
    CALL u%copy(w) ! Already formed by inner convergence check
  ELSE
    CALL gmres_backsolve(k,h,res,hdots)
    CALL u%maxpy(hdots(1:k),z(1:k),k)
  END IF
  ! Early inner exits based on |r|/|u| must be confirmed with the true residual
  IF(((nits==its).OR.((k<nrits).AND.(.NOT.inner_conv))).AND.(.NOT.self%full_residual))EXIT
  call A%apply(u,w)
  IF(associated(self%bc))call self%bc%apply(w)
  call r%add(0.d0,1.d0,g,-1.d0,w)
  !if(associated(self%cleaner))call self%cleaner%apply(r)
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = clamp_nonneg(r%dot(r)) ! Clamp round-off (see Arnoldi step)
  uu = clamp_nonneg(u%dot(u))
  IF((oft_env%pm.AND.oft_env%head_proc).AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
    WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  END IF
  IF(res_u_tol>0.d0.AND.SQRT(gg/uu)<res_u_tol)EXIT
  IF(SQRT(gg)<self%atol)EXIT
  IF(SQRT(gg/ggin)<self%rtol)EXIT
  IF(nits==its)EXIT
  IF((k<nrits).AND.(.NOT.inner_conv))EXIT
  IF(uu==uuold)EXIT
end do
END IF
CALL g%copy(r)
deallocate(c,s,h,res,hdots,vw,znorm)
self%cits=nits
IF(nits==its)self%cits=-1
DEBUG_STACK_POP
end subroutine gmres_solver_apply
!------------------------------------------------------------------------------
!> Setup GMRES solver from XML definition
!------------------------------------------------------------------------------
subroutine gmres_setup_xml(self,solver_node,level)
CLASS(oft_native_gmres_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_int
REAL(r8), POINTER, DIMENSION(:) :: read_real
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `its` node","gmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `its` values specified","gmres_setup_xml",__FILE__)
      self%its=read_int(val_level)
    ELSE
      self%its=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"nrits",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `nrits` node","gmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `nrits` values specified","gmres_setup_xml",__FILE__)
      self%nrits=read_int(val_level)
    ELSE
      self%nrits=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"atol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `atol` node","gmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `atol` values specified","gmres_setup_xml",__FILE__)
      self%atol=read_real(val_level)
    ELSE
      self%atol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"rtol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `rtol` node","gmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `rtol` values specified","gmres_setup_xml",__FILE__)
      self%rtol=read_real(val_level)
    ELSE
      self%rtol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'GMRES solver setup:'
  WRITE(*,'(2X,A,I4)')    '- Iterations:  ',self%its
  WRITE(*,'(2X,A,I4)')    '- Restart:     ',self%nrits
  WRITE(*,'(2X,A,ES10.2)')'- Abs-Tol:     ',self%atol
  WRITE(*,'(2X,A,ES10.2)')'- Rel-Tol:     ',self%rtol
END IF
DEBUG_STACK_POP
end subroutine gmres_setup_xml
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine gmres_delete(self,propogate)
class(oft_native_gmres_solver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
integer(i4) :: i
NULLIFY(self%A)
IF(self%initialized)THEN
  DO i=1,self%nrits+1
    CALL self%v(i)%f%delete
    CALL self%z(i)%f%delete
    DEALLOCATE(self%v(i)%f,self%z(i)%f)
  END DO
  DEALLOCATE(self%v,self%z)
  CALL self%r%delete
  CALL self%w%delete
  DEALLOCATE(self%r,self%w)
END IF
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
end subroutine gmres_delete
!------------------------------------------------------------------------------
!> Solve a linear system using the flexible GMRES method (complex)
!------------------------------------------------------------------------------
recursive subroutine cgmres_solver_apply(self,u,g)
class(oft_native_gmres_csolver), intent(inout) :: self !< Solver object
class(oft_cvector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_cvector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_cmatrix), pointer :: A
logical :: pm_save,pre_full_res,inner_conv
integer(i4) :: i,j,k,nits,its,nrits
real(r8) :: uu,uuold,gg,ggold,ggin
real(r8) :: wnorm0,res_u_tol
real(r8), allocatable :: znorm(:)
complex(c8) :: delta,hkmi
complex(c8), allocatable :: h(:,:),c(:),s(:),res(:),hdots(:)
class(oft_cvector), pointer :: r,w
type(oft_cvector_ptr), pointer, dimension(:) :: v,z
type(oft_cvector_ptr), allocatable, dimension(:) :: vw
DEBUG_STACK_PUSH
!---
A=>self%A
its=self%its
nrits=self%nrits
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%v(nrits+1),self%z(nrits+1))
  DO i=1,nrits+1
    CALL u%new(self%v(i)%f)
    CALL u%new(self%z(i)%f)
  END DO
  CALL u%new(self%r)
  CALL u%new(self%w)
  CALL csolver_setup(self)
END IF
if((oft_env%pm.AND.oft_env%head_proc))write(*,'(A)')'Starting GMRES solver'
if(u%n/=g%n)call oft_abort('Vectors are not the same length','gmres_solver_apply',__FILE__)
!
if(associated(self%bc))call self%bc%apply(g)
if(associated(self%cleaner))call self%cleaner%apply(g)
if(associated(self%orthog))call self%orthog%apply(u)
!
v=>self%v
z=>self%z
r=>self%r
w=>self%w
CALL w%set((0.d0,0.d0))
allocate(c(nrits+1),s(nrits+1),res(nrits+1),h(nrits+1,nrits))
allocate(hdots(nrits+1),vw(nrits+1),znorm(nrits))
c=0.d0; s=0.d0; res=0.d0; h=0.d0
!
uu = clamp_nonneg(REAL(u%dot(u),r8))
IF(uu>0.d0)call A%apply(u,w)
if(associated(self%bc))call self%bc%apply(w)
!if(associated(self%cleaner))call self%cleaner%apply(g)
call r%add((0.d0,0.d0),(1.d0,0.d0),g,(-1.d0,0.d0),w)
!if(associated(self%cleaner))call self%cleaner%apply(r)
gg = clamp_nonneg(REAL(r%dot(r),r8)) ! Clamp round-off (see Arnoldi step)
ggin=gg
100 FORMAT (I6,2ES14.6)
110 FORMAT (I6,3ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,SQRT(uu),SQRT(gg)
! Tolerance on |r|/|u| for its=-1/-2 (distinct from self%rtol, which is relative to |r0|)
res_u_tol=-1.d0
IF(its==-1)res_u_tol=1.d-7
IF(its==-2)res_u_tol=1.d-14
nits=0
if(self%its/=0.AND.self%nrits/=0)THEN
do j=1,INT(MIN(u%ng,INT(HUGE(j),i8)),4)
  IF(gg==0.d0)EXIT
  res=0.d0
  res(1)=sqrt(gg)
  call v(1)%f%add((0.d0,0.d0),(1.d0,0.d0)/res(1),r)
  inner_conv=.FALSE.
  do i=1,nrits
    ! Precondition search direction
    if(associated(self%pre))then
      CALL w%copy(v(i)%f)
      CALL z(i)%f%set((0.d0,0.d0))
      pm_save=oft_env%pm; oft_env%pm=self%pre%pm
      pre_full_res=self%pre%full_residual
      self%pre%full_residual=.FALSE.
      call self%pre%apply(z(i)%f,w)
      self%pre%full_residual=pre_full_res
      oft_env%pm=pm_save
    else
      CALL z(i)%f%copy(v(i)%f)
    end if
    if(associated(self%bc))call self%bc%apply(z(i)%f)
    if(associated(self%orthog))call self%orthog%apply(z(i)%f)
    IF(res_u_tol>0.d0)znorm(i)=SQRT(clamp_nonneg(REAL(z(i)%f%dot(z(i)%f),r8)))
    call A%apply(z(i)%f,w)
    if(associated(self%bc))call self%bc%apply(w)
    !---Arnoldi iteration (CGS with selective re-orthogonalization)
    do k=1,i
      vw(k)%f=>v(k)%f
    end do
    vw(i+1)%f=>w
    hdots(1:i+1)=w%mdot(vw(1:i+1),i+1) ! Fuse projections and initial norm into one reduction
    h(1:i,i)=CONJG(hdots(1:i)) ! mdot gives w^H v_k
    wnorm0=SQRT(clamp_nonneg(REAL(hdots(i+1),r8)))
    CALL w%maxpy(-h(1:i,i),v(1:i),i)
    !if(associated(self%cleaner))call self%cleaner%apply(w)
    ! Safety net: after near-complete cancellation (eg. an exact preconditioner) clamp any
    ! negative round-off in the self dot-product to zero, so this is treated as a (lucky)
    ! breakdown instead of producing NaN. NaN itself is preserved so failures still propagate.
    h(i+1,i)=SQRT(clamp_nonneg(REAL(w%dot(w),r8)))
    IF(ABS(h(i+1,i))<gmres_reorth_tol*wnorm0)THEN ! Significant cancellation, re-orthogonalize
      hdots(1:i)=CONJG(w%mdot(v(1:i),i))
      CALL w%maxpy(-hdots(1:i),v(1:i),i)
      h(1:i,i)=h(1:i,i)+hdots(1:i)
      h(i+1,i)=SQRT(clamp_nonneg(REAL(w%dot(w),r8)))
    END IF
    IF(ABS(h(i+1,i))>0.d0)call v(i+1)%f%add((0.d0,0.d0),(1.d0,0.d0)/h(i+1,i),w)
    !if(associated(self%cleaner))call self%cleaner%apply(v(i+1)%f)
    !---Apply Givens rotation
    do k=2,i
      hkmi = h(k-1,i)
      h(k-1,i) = CONJG(c(k-1))*hkmi + CONJG(s(k-1))*h(k,i)
      h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
    enddo
    delta = SQRT( ABS(h(i,i))**2 + ABS(h(i+1,i))**2 )
    c(i) = h(i,i) / delta
    s(i) = h(i+1,i) / delta
    h(i,i) = CONJG(c(i))*h(i,i) + CONJG(s(i))*h(i+1,i)
    !---Update residual
    res(i+1) = -s(i)*res(i)
    res(i) = CONJG(c(i))*res(i)
    nits=nits+1
    IF(i==nrits)EXIT
    IF(ABS(res(i+1))==0.d0)EXIT ! Lucky breakdown
    IF(nits==its)EXIT
    IF(ABS(res(i+1))<self%atol)EXIT
    IF(ABS(res(i+1))/SQRT(ggin)<self%rtol)EXIT
    !---Check |r|/|u| criterion for current iterate (after `min_nits` inner iterations
    !   to avoid stalling through many short restart cycles)
    IF(res_u_tol>0.d0.AND.i>=self%min_nits)THEN
      CALL gmres_backsolve(i,h,res,hdots)
      ! Necessary condition using |u0 + Z*y| <= |u0| + sum_k |y_k|*|z_k|
      IF(ABS(res(i+1))<res_u_tol*(SQRT(uu)+SUM(ABS(hdots(1:i))*znorm(1:i))))THEN
        CALL w%copy(u) ! Form iterate in w (used as scratch)
        CALL w%maxpy(hdots(1:i),z(1:i),i)
        inner_conv=(ABS(res(i+1))<res_u_tol*SQRT(clamp_nonneg(REAL(w%dot(w),r8))))
        IF(inner_conv)EXIT
      END IF
    END IF
    IF(oft_env%head_proc.AND.oft_env%pm.AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
      WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
    END IF
  end do
  k=i
  ! Update iterate
  IF(inner_conv)THEN
    CALL u%copy(w) ! Already formed by inner convergence check
  ELSE
    CALL gmres_backsolve(k,h,res,hdots)
    CALL u%maxpy(hdots(1:k),z(1:k),k)
  END IF
  IF((nits==its).AND.(.NOT.self%full_residual))EXIT
  call A%apply(u,w)
  IF(associated(self%bc))call self%bc%apply(w)
  call r%add((0.d0,0.d0),(1.d0,0.d0),g,(-1.d0,0.d0),w)
  !if(associated(self%cleaner))call self%cleaner%apply(r)
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = clamp_nonneg(REAL(r%dot(r),r8)) ! Clamp round-off (see Arnoldi step)
  uu = clamp_nonneg(REAL(u%dot(u),r8))
  IF((oft_env%pm.AND.oft_env%head_proc).AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
    WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  END IF
  IF(res_u_tol>0.d0.AND.SQRT(gg/uu)<res_u_tol)EXIT
  IF(SQRT(gg)<self%atol)EXIT
  IF(SQRT(gg/ggin)<self%rtol)EXIT
  IF(nits==its)EXIT
  IF(uu==uuold)EXIT
end do
END IF
CALL g%copy(r)
deallocate(c,s,h,res,hdots,vw,znorm)
self%cits=nits
IF(nits==its)self%cits=-1
DEBUG_STACK_POP
end subroutine cgmres_solver_apply
!------------------------------------------------------------------------------
!> Setup GMRES solver from XML definition
!------------------------------------------------------------------------------
subroutine cgmres_setup_xml(self,solver_node,level)
CLASS(oft_native_gmres_csolver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_int
REAL(r8), POINTER, DIMENSION(:) :: read_real
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `its` node","cgmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `its` values specified","cgmres_setup_xml",__FILE__)
      self%its=read_int(val_level)
    ELSE
      self%its=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"nrits",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `nrits` node","cgmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `nrits` values specified","cgmres_setup_xml",__FILE__)
      self%nrits=read_int(val_level)
    ELSE
      self%nrits=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"atol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `atol` node","cgmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `atol` values specified","cgmres_setup_xml",__FILE__)
      self%atol=read_real(val_level)
    ELSE
      self%atol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"rtol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `rtol` node","cgmres_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `rtol` values specified","cgmres_setup_xml",__FILE__)
      self%rtol=read_real(val_level)
    ELSE
      self%rtol=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'GMRES solver setup:'
  WRITE(*,'(2X,A,I4)')    '- Iterations:  ',self%its
  WRITE(*,'(2X,A,I4)')    '- Restart:     ',self%nrits
  WRITE(*,'(2X,A,ES10.2)')'- Abs-Tol:     ',self%atol
  WRITE(*,'(2X,A,ES10.2)')'- Rel-Tol:     ',self%rtol
END IF
DEBUG_STACK_POP
end subroutine cgmres_setup_xml
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine cgmres_delete(self,propogate)
class(oft_native_gmres_csolver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
integer(i4) :: i
NULLIFY(self%A)
IF(self%initialized)THEN
  DO i=1,self%nrits+1
    CALL self%v(i)%f%delete
    CALL self%z(i)%f%delete
    DEALLOCATE(self%v(i)%f,self%z(i)%f)
  END DO
  DEALLOCATE(self%v,self%z)
  CALL self%r%delete
  CALL self%w%delete
  DEALLOCATE(self%r,self%w)
END IF
self%initialized=.FALSE.
CALL csolver_delete(self,propogate)
end subroutine cgmres_delete
!------------------------------------------------------------------------------
!> Solve the upper-triangular GMRES least-squares system \f$ R y = b \f$ (real)
!!
!! \f$ R \f$ is the leading `k x k` block of the Hessenberg matrix after Givens rotations
!------------------------------------------------------------------------------
subroutine gmres_backsolve_r8(k,h,b,y)
integer(i4), intent(in) :: k !< Size of system
real(r8), intent(in) :: h(:,:) !< Rotated Hessenberg matrix
real(r8), intent(in) :: b(:) !< Rotated residual vector
real(r8), intent(inout) :: y(:) !< Solution [k]
integer(i4) :: i
DEBUG_STACK_PUSH
y(1:k)=b(1:k)
DO i=k,1,-1
  y(i)=y(i)/h(i,i)
  y(1:i-1)=y(1:i-1)-y(i)*h(1:i-1,i)
END DO
DEBUG_STACK_POP
end subroutine gmres_backsolve_r8
!------------------------------------------------------------------------------
!> Solve the upper-triangular GMRES least-squares system \f$ R y = b \f$ (complex)
!!
!! \f$ R \f$ is the leading `k x k` block of the Hessenberg matrix after Givens rotations
!------------------------------------------------------------------------------
subroutine gmres_backsolve_c8(k,h,b,y)
integer(i4), intent(in) :: k !< Size of system
complex(c8), intent(in) :: h(:,:) !< Rotated Hessenberg matrix
complex(c8), intent(in) :: b(:) !< Rotated residual vector
complex(c8), intent(inout) :: y(:) !< Solution [k]
integer(i4) :: i
DEBUG_STACK_PUSH
y(1:k)=b(1:k)
DO i=k,1,-1
  y(i)=y(i)/h(i,i)
  y(1:i-1)=y(1:i-1)-y(i)*h(1:i-1,i)
END DO
DEBUG_STACK_POP
end subroutine gmres_backsolve_c8
!------------------------------------------------------------------------------
!> Cast an eigensolver object to a oft_native_cg_eigsolver
!!
!! The source matrix must be @ref oft_native_cg_eigsolver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION native_cg_eigsolver_cast(self,source) result(success)
type(oft_native_cg_eigsolver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_eigsolver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cg_eigsolver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION native_cg_eigsolver_cast
!------------------------------------------------------------------------------
!> Solve a general eigenvalue system using the Conjugate-Gradient method.
!!
!! This solver uses a Non-Linear Conjugate-Gradient method to minimize the
!! Rayleigh Quotient (\f$ R = \frac{x*A*x}{x*M*x} \f$).
!------------------------------------------------------------------------------
subroutine cg_eigsolver_apply(self,u,alam)
class(oft_native_cg_eigsolver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
real(r8), intent(inout) :: alam !< Eigenvalue
logical :: pm_save
integer(i4) :: n,ir,ninner,its
real(r8) :: e0,h0,ev,sf,gg,e1,h1,e2,h2,qa,qb,qc,qd,qq,alfa1,alfa2
real(r8) :: h01,e01,h02,e02,ev1,ev2,evlast,gglast,alfa,beta,alast,ggnew
class(oft_vector), pointer :: b,c,e,f,g,s
class(oft_matrix), pointer :: A,M
DEBUG_STACK_PUSH
A=>self%A
M=>self%M
its=self%its
self%cits=0
IF(.NOT.self%initialized)CALL eigsolver_setup(self)
if((oft_env%pm.AND.oft_env%head_proc))write(*,'(A)')'Starting CG eigensolver'
call u%new(b)
call u%new(c)
call u%new(e)
call u%new(f)
call u%new(g)
if(associated(self%pre))call u%new(s)
if(associated(self%bc))call self%bc%apply(u)
if(associated(self%orthog))call self%orthog%apply(u)
call A%apply(u,c)
call M%apply(u,b)
if(associated(self%bc))call self%bc%apply(b)
if(associated(self%bc))call self%bc%apply(c)
h0 = u%dot(b)
if(h0==0.d0)call oft_abort('Metric is Zero','cgeigsolver_apply',__FILE__)
e0 = u%dot(c)
if(e0<=0.d0)call oft_abort('Energy is Zero','cgeigsolver_apply',__FILE__)
ev=h0/e0
sf=1.d0/sqrt(e0)
e0=1.d0
h0=ev
gg=0.d0
call u%scale(sf)
call b%scale(sf)
call c%scale(sf)
call g%add(0.d0,1.d0,c,-1.d0/ev,b)
CALL f%copy(g)
if(associated(self%pre))then
  pm_save=oft_env%pm; oft_env%pm=self%pre%pm
  call self%pre%apply(s,f)
  oft_env%pm=pm_save
  if(associated(self%bc))call self%bc%apply(s)
  CALL f%copy(s)
end if
gg = g%dot(f)
ggnew = g%dot(g)
alam=1.d0/ev
100 FORMAT (I6,3ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,alam,ggnew
if(gg==0.d0.OR.self%its==0)then
  call b%delete
  call c%delete
  call e%delete
  call f%delete
  call g%delete
  deallocate(b,c,e,f,g)
  IF(associated(self%pre))THEN
    call s%delete
    deallocate(s)
  END IF
  DEBUG_STACK_POP
  return
end if
ninner = INT(u%ng,4)
if((self%ninner>0).AND.(self%ninner<u%ng))ninner=self%ninner
! iterate to find the solution
do ir=1,self%nrestarts ! restart loop
CALL e%copy(f)
do n=1,ninner
  if(associated(self%orthog))call self%orthog%apply(e)
  call A%apply(e,g)
  call M%apply(e,f)
  if(associated(self%bc))then
    call self%bc%apply(f)
    call self%bc%apply(g)
  end if
  h1 = u%dot(f)+e%dot(b)
  e1 = u%dot(g)+e%dot(c)
  h2 = e%dot(f)
  e2 = e%dot(g)
  qa=h2*e1-e2*h1; qb=h2*e0-e2*h0
  qc=h1*e0-e1*h0; qd=qb**2-qa*qc
  if(qd<0.d0)call oft_abort('Negative Discriminant','cgeigsolver_apply',__FILE__)
  qq=-qb-sqrt(qd)*sign(1.d0,qb)
  alfa1=qq/qa; alfa2=qc/qq
  h01=h0+alfa1*(h1+alfa1*h2)
  e01=e0+alfa1*(e1+alfa1*e2)
  h02=h0+alfa2*(h1+alfa2*h2)
  e02=e0+alfa2*(e1+alfa2*e2)
  ev1=h01/e01; ev2=h02/e02; evlast=ev
  if(ev1>ev2)then
    alfa=alfa1; ev=ev1; h0=h01; e0=e01
  else
    alfa=alfa2; ev=ev2; h0=h02; e0=e02
  endif
  sf=1.d0/sqrt(e0); e0=1.d0; h0=ev
  gglast=gg; gg=0.d0
  call u%add(sf,alfa*sf,e)
  call b%add(sf,alfa*sf,f)
  call c%add(sf,alfa*sf,g)
  call g%add(0.d0,1.d0,c,-1.d0/ev,b)
  CALL f%copy(g)
  if ( associated(self%pre) )then
    pm_save=oft_env%pm; oft_env%pm=self%pre%pm
    call self%pre%apply(s,f)
    oft_env%pm=pm_save
    if(associated(self%bc))call self%bc%apply(s)
    CALL f%copy(s)
  end if
  gg = g%dot(f)
  ggnew = g%dot(g)
  alam=1.d0/ev
  beta=gg/gglast
  if((n<=10.OR.MOD(n,self%itplot)==0).AND.(oft_env%pm.AND.oft_env%head_proc))write(*,100)n,alam,ggnew
  if(its==-1.AND.REAL(ev,4)==REAL(evlast,4))exit
  if(its==-2.AND.ev==evlast)exit
  if(gg==0.d0)exit
  call e%add(beta,1.d0,f)
end do ! end main loop
if(n<ninner)exit
end do ! end restart loop
call b%delete
call c%delete
call e%delete
call f%delete
call g%delete
DEALLOCATE(b,c,e,f,g)
IF(associated(self%pre))THEN
    call s%delete
    deallocate(s)
  END IF
IF(ir>self%nrestarts.AND.oft_env%head_proc)WRITE(*,*)'Full Its'
self%cits=n
DEBUG_STACK_POP
end subroutine cg_eigsolver_apply
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine cg_eig_delete(self)
class(oft_native_cg_eigsolver), intent(inout) :: self !< Solver object
NULLIFY(self%A,self%M)
self%initialized=.FALSE.
end subroutine cg_eig_delete
!------------------------------------------------------------------------------
!> Solve a nonlinear system Newton's method
!------------------------------------------------------------------------------
subroutine nksolver_apply(self,u,g)
class(oft_nksolver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_vector), pointer :: v
logical :: pm_save
integer(i4) :: i,its
real(r8) :: res,resp,uu,sback
real(r8) :: met_time,up_time,inv_time
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!---
met_time=0.d0
up_time=0.d0
inv_time=0.d0
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,'(A)')'Starting Newton solver'
IF(.NOT.ASSOCIATED(self%du))THEN
  CALL u%new(self%du)
END IF
CALL u%new(v)
self%lits=0
self%nlits=-1
self%cits=1
!IF(ASSOCIATED(self%bc))call self%bc%apply(g)
!---Initialize Newton loop
res=1.d99 ! Set residual to a large number
its=100
IF(self%its>0)its=self%its
outer: do i=0,its
  resp=res ! Retain old residual
  !---Initialize backtracking loop
  sback=.5d0 ! Set backtracking step size
  do while(.TRUE.)
    if(oft_env%head_proc)CALL mytimer%tick
    !---Compute current residual
    call self%A%apply(u,v)
    call v%add(1.d0,-1.d0,g)
    IF(ASSOCIATED(self%bc))call self%bc%apply(v)
    res=clamp_nonneg(v%dot(v))
    uu=clamp_nonneg(u%dot(u))
    if(oft_env%head_proc)met_time=met_time+mytimer%tock()
    if((res<resp).OR.(.NOT.self%backtrack))exit ! If residual has decreased exit
    !---If residual has reached an acceptable value exit
    IF(SQRT(res)<self%atol)EXIT
    IF(SQRT(res/uu)<self%rtol)EXIT
    !---Backtrack along Newton path to find minimum
    CALL u%add(1.d0,sback,self%du)
    !---Update backtrack step size and exit if step is too small
    sback=.5d0*sback
    if(sback<1.d-8)THEN
      self%cits=-1
      exit outer
    end if
  end do ! End of backtracking loop
  !---Print Newton residual
  if((oft_env%pm.AND.oft_env%head_proc))WRITE(*,'(I6,3ES14.6)')i,SQRT(uu),SQRT(res),SQRT(res/uu)
  IF(i==its)self%cits=-1
  !---If residual has reached an acceptable value exit
  IF(SQRT(res)<self%atol.AND.i>0)EXIT
  IF(SQRT(res/uu)<self%rtol.AND.i>0)EXIT
  !---Update jacobian
  IF(MOD(i,self%up_freq)==0.AND.ASSOCIATED(self%J_update))THEN
    if(oft_env%head_proc)CALL mytimer%tick
    CALL self%J_update(u)
    if(oft_env%head_proc)up_time=up_time+mytimer%tock()
  END IF
  !---Solve jacobian for new search direction
  if(oft_env%head_proc)CALL mytimer%tick
  call self%du%set(0.d0)
  pm_save=oft_env%pm; oft_env%pm=self%J_inv%pm;
  call self%J_inv%apply(self%du,v)
  oft_env%pm=pm_save
  IF(self%J_inv%cits<0)THEN
    self%cits=-1
    EXIT
  END IF
  self%lits=self%lits+self%J_inv%cits
  !---Update solution
  call u%add(1.d0,-1.d0,self%du)
  if(oft_env%head_proc)inv_time=inv_time+mytimer%tock()
end do outer ! End of Newton loop
IF((oft_env%pm.AND.oft_env%head_proc))THEN
  WRITE(*,'(2X,A)')'Timing:'
  WRITE(*,'(4X,A,ES11.3)')'Metric = ',met_time
  WRITE(*,'(4X,A,ES11.3)')'Update = ',up_time
  WRITE(*,'(4X,A,ES11.3)')'Invert = ',inv_time
END IF
self%nlits=i
CALL v%delete
DEALLOCATE(v)
DEBUG_STACK_POP
end subroutine nksolver_apply
!------------------------------------------------------------------------------
!> Destroy Newton solver and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine nksolver_delete(self)
class(oft_nksolver), intent(inout) :: self !< Solver object
IF(ASSOCIATED(self%du))THEN
  CALL self%du%delete
  DEALLOCATE(self%du)
END IF
IF(ASSOCIATED(self%A))NULLIFY(self%A)
IF(ASSOCIATED(self%J_inv))NULLIFY(self%J_inv)
IF(ASSOCIATED(self%J_update))NULLIFY(self%J_update)
IF(ASSOCIATED(self%bc))NULLIFY(self%bc)
end subroutine nksolver_delete
!------------------------------------------------------------------------------
!> Solve the trivial linear system for the identity matrix
!!
!! @note Typically used for simple cases when a solver wrapper is required but the trivial
!! case of \f$ u = g \f$ is the desired result
!------------------------------------------------------------------------------
subroutine identity_inv_apply(self,u,g)
class(oft_identity_inv), intent(inout) :: self !< Solver object
CLASS(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
DEBUG_STACK_PUSH
if(g%n/=u%n)call oft_abort('Size mismatch','identity_inv_apply',__FILE__)
CALL u%copy(g)
DEBUG_STACK_POP
end subroutine identity_inv_apply
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine identity_inv_delete(self,propogate)
class(oft_identity_inv), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
CALL solver_delete(self,propogate)
end subroutine identity_inv_delete
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_diag_scale
!!
!! The source matrix must be @ref oft_diag_scale or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION diag_scale_cast(self,source) result(success)
type(oft_diag_scale), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_diag_scale)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION diag_scale_cast
!------------------------------------------------------------------------------
!> Solve the linear system for a diagonal matrix
!------------------------------------------------------------------------------
subroutine diag_scale_apply(self,u,g)
class(oft_diag_scale), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
DEBUG_STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch.','diag_scale_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch.','diag_scale_apply',__FILE__)
if(.NOT.associated(self%A%D))call oft_abort('No diagonal entries.','diag_scale_apply',__FILE__)
!---
CALL u%copy(g)
call u%mult(self%A%D,div_flag=.TRUE.)
DEBUG_STACK_POP
end subroutine diag_scale_apply
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine diag_scale_delete(self,propogate)
class(oft_diag_scale), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
NULLIFY(self%A)
CALL solver_delete(self,propogate)
end subroutine diag_scale_delete
!------------------------------------------------------------------------------
!> Solve the linear system for a diagonal matrix (complex)
!------------------------------------------------------------------------------
subroutine cdiag_scale_apply(self,u,g)
class(oft_diag_cscale), intent(inout) :: self !< Solver object
class(oft_cvector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_cvector), intent(inout) :: g !< RHS (input), Residual (output)
DEBUG_STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch.','cdiag_scale_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch.','cdiag_scale_apply',__FILE__)
if(.NOT.associated(self%A%D))call oft_abort('No diagonal entries.','cdiag_scale_apply',__FILE__)
!---
call u%add((0.d0,0.d0),(1.d0,0.0),g)
call u%mult(self%A%D,div_flag=.TRUE.)
DEBUG_STACK_POP
end subroutine cdiag_scale_apply
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine cdiag_scale_delete(self,propogate)
class(oft_diag_cscale), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
NULLIFY(self%A)
CALL csolver_delete(self,propogate)
end subroutine cdiag_scale_delete
!------------------------------------------------------------------------------
!> Cast a solver object to a jblock_precond_cast
!!
!! The source matrix must be @ref jblock_precond_cast or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION jblock_precond_cast(self,source) result(success)
type(oft_jblock_precond), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_jblock_precond)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION jblock_precond_cast
!------------------------------------------------------------------------------
!> Apply 1-step of a symmetric Jacobi smoother with native CRS matrices
!------------------------------------------------------------------------------
subroutine jblock_precond_apply(self,u,g)
class(oft_jblock_precond), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_vector), pointer :: p,q
class(oft_native_vector), pointer :: uv,gv,md
INTEGER(i4) :: k
DEBUG_STACK_PUSH
IF(g%n/=self%A%nr)CALL oft_abort('Row mismatch','jblock_precond_apply',__FILE__)
IF(u%n/=self%A%nc)CALL oft_abort('Col mismatch','jblock_precond_apply',__FILE__)
IF(.NOT.native_vector_cast(uv,u))CALL oft_abort('"u" is not a vector object.','jblock_precond_apply',__FILE__)
IF(.NOT.native_vector_cast(gv,g))CALL oft_abort('"g" is not a vector object.','jblock_precond_apply',__FILE__)
IF(.NOT.native_vector_cast(md,self%A%D))CALL oft_abort('"A%D" is not a vector object.','jblock_precond_apply',__FILE__)
!---Initialize solver
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%u_save(u%n),self%g_save(u%n))
  CALL u%new(p)
  IF(.NOT.native_vector_cast(self%p,p))CALL oft_abort('"p" is not a vector object.','jblock_precond_apply',__FILE__)
  CALL u%new(q)
  IF(.NOT.native_vector_cast(self%q,q))CALL oft_abort('"q" is not a vector object.','jblock_precond_apply',__FILE__)
  NULLIFY(p,q)
  CALL solver_setup(self)
END IF
!---
IF(self%warn_once.AND.self%its<=0)THEN
  CALL oft_warn("SymJacobi smoother called with zero iterations [suppressed].")
  self%warn_once=.FALSE.
END IF
!---
if(self%down) then ! Down smoother
  CALL nk_jacobi_start(u%n,self%df,md%v,gv%v,uv%v,self%p%v)
  IF(ASSOCIATED(self%bc))call self%bc%apply(self%p)
  do k=1,self%its
    call self%A%apply(self%p,self%q)
    IF(ASSOCIATED(self%bc))call self%bc%apply(self%q)
    IF(k<self%its)THEN
      CALL nk_jacobi_sweep(u%n,self%df,md%v,self%q%v,gv%v,self%p%v,u=uv%v)
    ELSE ! Last iteration, also save state for up smoother (u and g are final here)
      CALL nk_jacobi_sweep(u%n,self%df,md%v,self%q%v,gv%v,self%p%v,u=uv%v, &
        usave=self%u_save,gsave=self%g_save)
    END IF
    IF(ASSOCIATED(self%bc))call self%bc%apply(self%p)
  end do
  IF(self%its<1)THEN
    CALL nk_copy(u%n,uv%v,self%u_save)
    CALL nk_copy(u%n,gv%v,self%g_save)
  END IF
  self%down=.FALSE.
else ! Up smoother
  CALL nk_jacobi_restart(u%n,uv%v,self%u_save,self%g_save,gv%v,self%p%v)
  IF(ASSOCIATED(self%bc))call self%bc%apply(self%p)
  !---The update of u with the (BC-applied) correction from step k is deferred into the
  !---sweep of step k+1, which is valid as neither A%apply nor bc%apply(q) touch u
  do k=1,self%its
    call self%A%apply(self%p,self%q)
    IF(ASSOCIATED(self%bc))call self%bc%apply(self%q)
    IF(k==1)THEN
      CALL nk_jacobi_sweep(u%n,self%df,md%v,self%q%v,gv%v,self%p%v)
    ELSE
      CALL nk_jacobi_sweep(u%n,self%df,md%v,self%q%v,gv%v,self%p%v,u=uv%v)
    END IF
    IF(ASSOCIATED(self%bc))call self%bc%apply(self%p)
  end do
  IF(self%its>0)CALL nk_axpby(u%n,1.d0,1.d0,self%p%v,uv%v)
  self%down=.TRUE.
end if
DEBUG_STACK_POP
end subroutine jblock_precond_apply
!------------------------------------------------------------------------------
!> Setup symmetric Jacobi smoother from XML definition
!------------------------------------------------------------------------------
subroutine jblock_setup_xml(self,solver_node,level)
CLASS(oft_jblock_precond), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_int
REAL(r8), POINTER, DIMENSION(:) :: read_real
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_int,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `its` node","jblock_setup_xml",__FILE__)
  IF(ASSOCIATED(read_int))THEN
    IF(SIZE(read_int)>1)THEN
      IF(SIZE(read_int)<val_level)CALL oft_abort("Not enough `its` values specified","jblock_setup_xml",__FILE__)
      self%its=read_int(val_level)
    ELSE
      self%its=read_int(1)
    END IF
    DEALLOCATE(read_int)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"df",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_real,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `df` node","jblock_setup_xml",__FILE__)
  IF(ASSOCIATED(read_real))THEN
    IF(SIZE(read_real)>1)THEN
      IF(SIZE(read_real)<val_level)CALL oft_abort("Not enough `df` values specified","jblock_setup_xml",__FILE__)
      self%df=read_real(val_level)
    ELSE
      self%df=read_real(1)
    END IF
    DEALLOCATE(read_real)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,*)'S-Jacobi solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Factor:      ',self%df
END IF
DEBUG_STACK_POP
end subroutine jblock_setup_xml
!------------------------------------------------------------------------------
!> Destroy symmetric Jacobi preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine jblock_precond_delete(self,propogate)
class(oft_jblock_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
DEBUG_STACK_PUSH
IF(self%initialized)THEN
  DEALLOCATE(self%u_save,self%g_save)
  CALL self%p%delete
  CALL self%q%delete
  DEALLOCATE(self%p,self%q)
END IF
NULLIFY(self%A)
NULLIFY(self%bc)
self%warn_once=.TRUE.
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
DEBUG_STACK_POP
end subroutine jblock_precond_delete
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_ml_precond
!!
!! The source matrix must be @ref oft_ml_precond or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION ml_precond_cast(self,source) result(success)
class(oft_ml_precond), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_ml_precond)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION ml_precond_cast
!------------------------------------------------------------------------------
!> Print solver configuration
!------------------------------------------------------------------------------
recursive subroutine ml_precond_view(self)
class(oft_ml_precond), intent(inout) :: self !< Solver object
!---
IF(ASSOCIATED(self%base_solve))THEN
  CALL self%base_solve%view
END IF
WRITE(*,*)'Multi-Grid level = ',self%level
WRITE(*,*)'  Timings: ',REAL(self%timings,4)
self%timings=0.d0
IF(ASSOCIATED(self%smooth_down))THEN
  WRITE(*,*)'Down Smoother:'
  CALL self%smooth_down%view
END IF
IF(ASSOCIATED(self%smooth_up))THEN
  WRITE(*,*)'Up Smoother:'
  CALL self%smooth_up%view
END IF
end subroutine ml_precond_view
!------------------------------------------------------------------------------
!> Apply 1-step of a multi-level preconditioner
!------------------------------------------------------------------------------
recursive subroutine ml_precond_apply(self,u,g)
class(oft_ml_precond), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!---Get RHS
IF(.NOT.self%initialized)THEN
  CALL u%new(self%r)
  IF(ASSOCIATED(self%base_solve))THEN
    CALL u%new(self%p)
    CALL self%ml_vecspace%vec_create(self%ucors,self%level-1)
    CALL self%ml_vecspace%vec_create(self%gcors,self%level-1)
  END IF
  CALL solver_setup(self)
END IF
CALL self%r%copy(g)
!---Report input residual
IF(self%report_err>0)CALL print_residual('In:')
!---Apply down smoother
if(associated(self%smooth_down))then
  CALL mytimer%tick
  call self%smooth_down%apply(u,g)
  self%timings(1)=self%timings(1)+mytimer%tock()
  !---Report down residual
  IF(self%report_err>0)CALL print_residual('Down:')
end if
if(associated(self%base_solve))then
  CALL mytimer%tick
  !---Inject RHS
  call self%ucors%set(0.d0)
  CALL self%ml_vecspace%inject(g,self%gcors)
  self%timings(2)=self%timings(2)+mytimer%tock()
  CALL mytimer%tick
  !---Apply lower level smoother
  self%base_solve%full_residual=.FALSE.
  call self%base_solve%apply(self%ucors,self%gcors)
  self%base_solve%full_residual=.TRUE.
  self%timings(3)=self%timings(3)+mytimer%tock()
  CALL mytimer%tick
  !---Interpolate Solution
  CALL self%ml_vecspace%interp(self%ucors,self%p)
  if(associated(self%bc))call self%bc%apply(self%p)
  !---Add coarse correction
  call u%add(1.d0,1.d0,self%p)
  self%timings(4)=self%timings(4)+mytimer%tock()
end if
!---Reset RHS
CALL g%copy(self%r)
!---Apply up smoother
if(associated(self%smooth_up))then
  !---Report up residual
  IF((self%report_err>0).AND.associated(self%base_solve))THEN
    CALL self%smooth_up%A%apply(u,g)
    CALL g%add(1.d0,-1.d0,self%r)
    CALL print_residual('Up:')
    CALL g%copy(self%r)
  END IF
  !
  CALL mytimer%tick
  IF(self%report_err==0)self%smooth_up%full_residual=.FALSE.
  call self%smooth_up%apply(u,g)
  self%smooth_up%full_residual=.TRUE.
  self%timings(1)=self%timings(1)+mytimer%tock()
end if
!---Report down residual
IF(self%report_err>0)CALL print_residual('Out:')
!---Complete
DEBUG_STACK_POP
contains
subroutine print_residual(stage)
character(LEN=*), intent(in) :: stage
class(oft_vector), pointer :: gtmp
INTEGER(i4) :: i
REAL(r8) :: uu,res
real(r8), pointer :: vals(:)
IF(oft_env%head_proc)WRITE(*,'(2X,A,A,I4)')'Residual ',stage,self%level
IF(self%report_err==1)THEN
  res = SQRT(clamp_nonneg(g%dot(g)))
  uu = SQRT(clamp_nonneg(u%dot(u)))
  IF(oft_env%head_proc)WRITE(*,'(4X,2ES14.6)')uu,res
ELSE IF(self%report_err==2)THEN
  CALL u%new(gtmp)
  ALLOCATE(vals(g%n))
  DO i=1,g%nblocks
    CALL gtmp%set(0.d0)
    vals=0.d0
    CALL u%get_local(vals,i)
    CALL gtmp%restore_local(vals,i)
    uu = SQRT(clamp_nonneg(gtmp%dot(gtmp)))
    CALL gtmp%set(0.d0)
    vals=0.d0
    CALL g%get_local(vals,i)
    CALL gtmp%restore_local(vals,i)
    res = SQRT(clamp_nonneg(gtmp%dot(gtmp)))
    IF(oft_env%head_proc)WRITE(*,'(4X,2ES14.6)')uu,res
  END DO
  CALL gtmp%delete()
  DEALLOCATE(gtmp,vals)
END IF
IF(oft_env%head_proc)WRITE(*,*)
end subroutine print_residual
end subroutine ml_precond_apply
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine ml_precond_update(self,new_pattern)
class(oft_ml_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern
!---
IF(ASSOCIATED(self%smooth_up))CALL self%smooth_up%update(new_pattern)
IF(ASSOCIATED(self%smooth_down))CALL self%smooth_down%update(new_pattern)
IF(ASSOCIATED(self%base_solve))CALL self%base_solve%update(new_pattern)
end subroutine ml_precond_update
!------------------------------------------------------------------------------
!> Destroy Multi-Level preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
recursive subroutine ml_precond_delete(self,propogate)
class(oft_ml_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
DEBUG_STACK_PUSH
if(ASSOCIATED(self%smooth_up))then
  call self%smooth_up%delete(propogate)
  DEALLOCATE(self%smooth_up)
end if
if(ASSOCIATED(self%smooth_down))then
  IF(self%symmetric)THEN
    NULLIFY(self%smooth_down)
  ELSE
    call self%smooth_down%delete(propogate)
    DEALLOCATE(self%smooth_down)
  END IF
end if
if(associated(self%base_solve))then
  call self%base_solve%delete(propogate)
  DEALLOCATE(self%base_solve)
end if
! IF(ASSOCIATED(self%interp))NULLIFY(self%interp)
! IF(ASSOCIATED(self%inject))NULLIFY(self%inject)
! IF(ASSOCIATED(self%vec_create))NULLIFY(self%vec_create)
IF(ASSOCIATED(self%ml_vecspace))NULLIFY(self%ml_vecspace)
IF(ASSOCIATED(self%bc))NULLIFY(self%bc)
IF(self%initialized)THEN
  CALL self%r%delete
  DEALLOCATE(self%r)
  IF(ASSOCIATED(self%p))THEN
    CALL self%p%delete
    CALL self%ucors%delete
    CALL self%gcors%delete
    DEALLOCATE(self%p,self%ucors,self%gcors)
  END IF
  self%initialized=.FALSE.
END IF
DEBUG_STACK_POP
end subroutine ml_precond_delete
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_ml_trans
!!
!! The source matrix must be @ref oft_ml_trans or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION ml_trans_cast(self,source) result(success)
class(oft_ml_trans), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_ml_trans)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION ml_trans_cast
!------------------------------------------------------------------------------
!> Transfer solution between distributed and shared levels as part of a ML
!! preconditioner
!------------------------------------------------------------------------------
recursive subroutine ml_trans_apply(self,u,g)
class(oft_ml_trans), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_vector), pointer :: ucors,gcors
DEBUG_STACK_PUSH
!---Create base level fields
CALL self%ml_vecspace%vec_create(ucors,self%level-1)
CALL self%ml_vecspace%vec_create(gcors,self%level-1)
!---Push RHS to base level
call self%ml_vecspace%inject(g,gcors)
!---Solve on base levels
call self%base_solve%apply(ucors,gcors)
!---Pop solution from base level
call self%ml_vecspace%interp(ucors,u)
!---Delete base level fields
call ucors%delete()
call gcors%delete()
!---Complete
DEBUG_STACK_POP
end subroutine ml_trans_apply
!------------------------------------------------------------------------------
!> Destroy Multi-Level preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine ml_trans_delete(self,propogate)
class(oft_ml_trans), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
DEBUG_STACK_PUSH
if(associated(self%smooth_up))then
  call self%smooth_up%delete(propogate)
  nullify(self%smooth_up)
end if
if(associated(self%smooth_down))then
  call self%smooth_down%delete(propogate)
  nullify(self%smooth_down)
end if
if(associated(self%base_solve))then
  call self%base_solve%delete(propogate)
  nullify(self%base_solve)
end if
! IF(ASSOCIATED(self%interp))NULLIFY(self%interp)
! IF(ASSOCIATED(self%inject))NULLIFY(self%inject)
! IF(ASSOCIATED(self%vec_create))NULLIFY(self%vec_create)
IF(ASSOCIATED(self%ml_vecspace))NULLIFY(self%ml_vecspace)
IF(ASSOCIATED(self%bc))NULLIFY(self%bc)
DEBUG_STACK_POP
end subroutine ml_trans_delete
!------------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE bjprecond_apply(self,u,g)
CLASS(oft_bjprecond), INTENT(inout) :: self !< Solver object
CLASS(oft_vector), INTENT(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), INTENT(inout) :: g !< RHS (input), Residual (output)
!---
LOGICAL :: thread_safe,color_avail,lcopy
LOGICAL, POINTER, DIMENSION(:) :: eflag
INTEGER(i4) :: i,j,jr,k,ngroups,ncolors
INTEGER(i4), POINTER, DIMENSION(:) :: part,part_tmp,sflag
TYPE(oft_native_vector) :: uloc
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: uvals,gvals,uinit
REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: nblk
LOGICAL, ALLOCATABLE, DIMENSION(:) :: allowed
INTEGER(i4) :: nlev
TYPE(oft_index_set) :: bset
TYPE(oft_graph) :: graph
CLASS(oft_native_matrix), POINTER :: A_native
CLASS(oft_solver), POINTER :: pretmp
TYPE(oft_1d_int), pointer, dimension(:) :: parts_tmp => NULL()
DEBUG_STACK_PUSH
!---
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'bjprecond_apply',__FILE__)
!---Initialize local matrix
IF(.NOT.self%initialized)THEN
  IF(.NOT.ASSOCIATED(self%pre))CALL oft_abort('No sub-solver (preconditioner defined)', &
    'bjprecond_apply',__FILE__)
  color_avail=.FALSE.
  CALL solver_setup(self)
  IF(self%nlocal==-1)THEN
    self%nlocal=-self%A%ni
    self%slice_group(1:self%A%ni)=(/(i,i=1,self%A%ni)/)
  END IF
  IF(ASSOCIATED(A_native%color))THEN
    ncolors=MAXVAL(A_native%color)
    IF((self%nlocal<1).AND.(ncolors>1))THEN
      self%nlocal=self%nlocal*ncolors
      color_avail=.TRUE.
    ELSE IF((self%nlocal>=1).AND.(ncolors>1))THEN
      IF(oft_env%head_proc.AND.self%nlocal/=ncolors)WRITE(*,'(A,X,I4)')' Updating number of local parts to match coloring',ncolors
      self%nlocal=ncolors
      color_avail=.TRUE.
    END IF
  END IF
  ALLOCATE(self%alocals(ABS(self%nlocal)),self%solvers(ABS(self%nlocal)))
  pretmp=>self%pre
  IF(self%nlocal==1)THEN
    IF(self%boundary_overlap)THEN
      CALL native_matrix_setup_full(A_native,u)
      CALL A_native%update_slice
    ELSE
      CALL self%alocals(1)%setup(self%A,u)
    END IF
    ALLOCATE(self%solvers(1)%s, source=pretmp)
    IF(ASSOCIATED(pretmp%pre))ALLOCATE(self%solvers(1)%s%pre, source=pretmp%pre)
  ELSE
    DEALLOCATE(self%alocals)
    ALLOCATE(self%alocals(0:ABS(self%nlocal)))
    ALLOCATE(part(A_native%nr))
    part=-1
    IF(self%boundary_overlap)THEN
      CALL native_matrix_setup_full(A_native,u)
      CALL A_native%update_slice
      IF(self%nlocal>0)THEN
        IF(color_avail)THEN
          part=A_native%color
        ELSE
          !---Partition local matrix
          graph%nr=A_native%nr
          graph%nnz=A_native%nnz
          graph%kr=>A_native%kr
          graph%lc=>A_native%lc
          CALL partition_graph(graph,self%nlocal,part)
        END IF
      ELSE
        !---Use slice partitioning
        self%nlocal=ABS(self%nlocal)
        IF(color_avail)THEN
          ngroups=MAXVAL(self%slice_group)
          j=0
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%n
              part(j+k)=self%slice_group(i)+(A_native%color(j+k)-1)*ngroups
            END DO
            j=j+self%A%i_map(i)%n
          END DO
        ELSE
          j=0
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%n
              part(j+k)=self%slice_group(i)
            END DO
            j=j+self%A%i_map(i)%n
          END DO
        END IF
      END IF
    ELSE
      IF(self%nlocal>0)THEN
        IF(color_avail)THEN
          CALL native_matrix_setup_full(A_native,u)
          CALL A_native%update_slice
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%nslice
              j=self%A%i_map(i)%offset+self%A%i_map(i)%slice(k)
              part(j)=A_native%color(j)
            END DO
          END DO
        ELSE
          !---Partition local matrix
          CALL self%alocals(0)%setup(self%A,u)
          graph%nr=self%alocals(0)%nr
          graph%nnz=self%alocals(0)%nnz
          graph%kr=>self%alocals(0)%kr
          graph%lc=>self%alocals(0)%lc
          ALLOCATE(part_tmp(graph%nr))
          CALL partition_graph(graph,self%nlocal,part_tmp)
          j=0
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%nslice
              part(self%A%i_map(i)%offset+self%A%i_map(i)%slice(k))=part_tmp(j+k)
            END DO
            j=j+self%A%i_map(i)%nslice
          END DO
          DEALLOCATE(part_tmp)
          CALL self%alocals(0)%delete
        END IF
      ELSE
        CALL native_matrix_setup_full(A_native,u)
        CALL A_native%update_slice
        !---Use slice partitioning
        self%nlocal=ABS(self%nlocal)
        IF(color_avail)THEN
          ngroups=MAXVAL(self%slice_group)
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%nslice
              j=self%A%i_map(i)%offset+self%A%i_map(i)%slice(k)
              part(j)=self%slice_group(i)+(A_native%color(j)-1)*ngroups
            END DO
          END DO
        ELSE
          DO i=1,self%A%ni
            DO k=1,self%A%i_map(i)%nslice
              j=self%A%i_map(i)%offset+self%A%i_map(i)%slice(k)
              part(j)=self%slice_group(i)
            END DO
          END DO
        END IF
      END IF
    END IF
    !---Do not include redundant rows
    IF(ASSOCIATED(A_native%redind))THEN
      DO i=1,A_native%nred
        part(A_native%redind(i))=-1
      END DO
    END IF
    !---Create partition mappings
    ALLOCATE(self%parts(self%nlocal))
    IF(ANY(part==0))CALL oft_abort("Unclaimed row after partitioning","bjprecond_apply",__FILE__)
    DO i=1,A_native%nr
      IF(part(i)<0)CYCLE
      self%parts(part(i))%n=self%parts(part(i))%n+1
    END DO
    DO i=1,self%nlocal
      ALLOCATE(self%parts(i)%v(self%parts(i)%n))
      self%parts(i)%n=0
    END DO
    DO i=1,A_native%nr
      IF(part(i)<0)CYCLE
      self%parts(part(i))%n=self%parts(part(i))%n+1
      self%parts(part(i))%v(self%parts(part(i))%n) = i
    END DO
    nlev=0
    IF(self%loverlap)nlev=1
    IF(self%overlap_levels>=0)nlev=self%overlap_levels
    IF(nlev>0.AND.self%nlocal>1)THEN
      !---Grow each block along the matrix graph, within rows of the same slice group
      ALLOCATE(sflag(A_native%nr),allowed(A_native%nr))
      sflag=-1
      DO i=1,self%A%ni
        IF(self%boundary_overlap)THEN
          DO k=1,self%A%i_map(i)%n
            sflag(self%A%i_map(i)%offset+k)=self%slice_group(i)
          END DO
        ELSE
          DO k=1,self%A%i_map(i)%nslice
            j=self%A%i_map(i)%offset+self%A%i_map(i)%slice(k)
            sflag(j)=self%slice_group(i)
          END DO
        END IF
      END DO
      allowed=(part>=0)
      !$omp parallel do private(bset)
      DO i=1,self%nlocal
        CALL index_set_create(bset,self%parts(i)%n,self%parts(i)%v)
        CALL index_set_grow(bset,A_native%nr,A_native%kr,A_native%lc,nlev,allowed=allowed,same=sflag)
        DEALLOCATE(self%parts(i)%v)
        self%parts(i)%n=bset%n
        ALLOCATE(self%parts(i)%v(bset%n))
        IF(bset%n>0)self%parts(i)%v=bset%ind
        CALL index_set_delete(bset)
      END DO
      DEALLOCATE(sflag,allowed)
    END IF
    !---Create dummy local vector
    uloc%n=A_native%nr
    ALLOCATE(uloc%stitch_info)
    uloc%stitch_info%skip=.TRUE.
    !---Create local sub-matrices
    !$omp parallel do
    DO i=1,self%nlocal
      CALL self%alocals(i)%setup(A_native,uloc,part=self%parts(i)%v)
      ALLOCATE(self%solvers(i)%s, source=pretmp)
      IF(ASSOCIATED(pretmp%pre))ALLOCATE(self%solvers(i)%s%pre, source=pretmp%pre)
    END DO
    !---Block gather/scatter plans (restricted write-back on overlap rows)
    ALLOCATE(self%bplans(self%nlocal))
    IF(nlev>0.AND.self%nlocal>1.AND.(.NOT.self%additive))THEN
      !$omp parallel do private(j,wt)
      DO i=1,self%nlocal
        ALLOCATE(wt(self%parts(i)%n))
        DO j=1,self%parts(i)%n
          wt(j)=MERGE(1.d0,0.d0,part(self%parts(i)%v(j))==i)
        END DO
        CALL scatter_plan_setup_local(self%bplans(i),self%parts(i)%n,self%parts(i)%v,wt=wt)
        DEALLOCATE(wt)
      END DO
    ELSE IF(nlev>0.AND.self%nlocal>1.AND.self%weighted)THEN
      !---Partition of unity: one over the number of blocks containing each row
      ALLOCATE(nblk(A_native%nr))
      nblk=0
      DO i=1,self%nlocal
        DO j=1,self%parts(i)%n
          nblk(self%parts(i)%v(j))=nblk(self%parts(i)%v(j))+1
        END DO
      END DO
      !$omp parallel do private(j,wt)
      DO i=1,self%nlocal
        ALLOCATE(wt(self%parts(i)%n))
        DO j=1,self%parts(i)%n
          wt(j)=1.d0/REAL(nblk(self%parts(i)%v(j)),r8)
        END DO
        CALL scatter_plan_setup_local(self%bplans(i),self%parts(i)%n,self%parts(i)%v,wt=wt)
        DEALLOCATE(wt)
      END DO
      DEALLOCATE(nblk)
    ELSE
      DO i=1,self%nlocal
        CALL scatter_plan_setup_local(self%bplans(i),self%parts(i)%n,self%parts(i)%v)
      END DO
    END IF
    !---Rows covered by the blocks (zeroed before an additive combination)
    IF(self%additive)THEN
      ALLOCATE(self%lunion(COUNT(part>=0)))
      j=0
      DO i=1,A_native%nr
        IF(part(i)<0)CYCLE
        j=j+1
        self%lunion(j)=i
      END DO
    END IF
    DEALLOCATE(part)
    DEALLOCATE(uloc%stitch_info)
  END IF
END IF
!---
thread_safe=self%pre%check_thread()
!---
IF(self%nlocal>1)THEN
  IF(.NOT.ASSOCIATED(self%local_stitch))THEN
    ALLOCATE(self%local_stitch(self%nlocal))
    DO i=1,self%nlocal
      self%local_stitch(i)%full=.TRUE.
      self%local_stitch(i)%nbe=0
      ALLOCATE(self%local_stitch(i)%be(self%parts(i)%n))
      self%local_stitch(i)%be=.FALSE.
    END DO
  END IF
  !---Persistent local vectors for each block
  IF(.NOT.ASSOCIATED(self%ulocs))ALLOCATE(self%ulocs(self%nlocal),self%glocs(self%nlocal))
  DO i=1,self%nlocal
    CALL bjprecond_local_vec(self%ulocs(i),self%parts(i)%n,self%local_stitch(i))
    CALL bjprecond_local_vec(self%glocs(i),self%parts(i)%n,self%local_stitch(i))
  END DO
  IF(self%update_slice)THEN
    CALL A_native%update_slice
    !$omp parallel do
    DO i=1,self%nlocal
      CALL self%alocals(i)%update_slice
    END DO
  END IF
  !---Borrow local values (no copy for native vectors)
  NULLIFY(uvals,gvals)
  CALL u%borrow_local(uvals,OFT_ACCESS_READWRITE)
  CALL g%borrow_local(gvals,OFT_ACCESS_READ)
  !---With overlap (restricted), blocks read rows written back by other blocks. Take
  !   initial guesses from a copy so they do not depend on the order blocks are solved in.
  lcopy=(.NOT.self%additive).AND.(self%overlap_levels>0.OR.(self%overlap_levels<0.AND.self%loverlap))
  IF(lcopy)THEN
    ALLOCATE(uinit(SIZE(uvals)))
    uinit=uvals
  ELSE
    uinit=>uvals
  END IF
  !---Only open a parallel region if the sub-solver is thread safe and we are
  !   not already in an active parallel region. Sub-solvers that are not thread
  !   safe (eg. SuperLU-DIST, MUMPS) may use OpenMP internally, and even an
  !   inactive outer region makes their teams nested. On targets where libgomp
  !   uses its POSIX barrier (eg. macOS) nested teams that create tasks can
  !   deadlock or release threads early. Each block (and its local vectors) is
  !   only accessed by one thread.
  IF(thread_safe.AND.(.NOT.omp_in_parallel()))THEN
    !$omp parallel do schedule(dynamic,1)
    DO i=1,self%nlocal
      CALL apply_block(i)
    END DO
  ELSE
    DO i=1,self%nlocal
      CALL apply_block(i)
    END DO
  END IF
  IF(lcopy)DEALLOCATE(uinit)
  !---Additive combination: sum of the (weighted) block solutions (blocks in order)
  IF(self%additive)THEN
    DO j=1,SIZE(self%lunion)
      uvals(self%lunion(j))=0.d0
    END DO
    DO i=1,self%nlocal
      CALL scatter_local(self%bplans(i),self%ulocs(i)%v,uvals,OFT_SCATTER_ADD,reverse=.TRUE.)
    END DO
  END IF
  !---Local solves on shared rows differ between processes, so synchronize from owner
  CALL u%release_local(uvals,OFT_ACCESS_READWRITE,sync=.TRUE.)
  CALL g%release_local(gvals,OFT_ACCESS_READ,sync=.FALSE.)
ELSE
  !---Persistent local vectors, with a local seam that skips all stitching/reductions
  IF(.NOT.ASSOCIATED(self%local_stitch))THEN
    ALLOCATE(self%local_stitch(1))
    self%local_stitch(1)%skip=.TRUE.
  END IF
  IF(.NOT.ASSOCIATED(self%ulocs))ALLOCATE(self%ulocs(1),self%glocs(1))
  IF(self%boundary_overlap)THEN
    IF(self%update_slice)CALL A_native%update_slice
    !---Setup local vectors
    CALL bjprecond_local_vec(self%glocs(1),g%n,self%local_stitch(1))
    CALL bjprecond_local_vec(self%ulocs(1),u%n,self%local_stitch(1))
    !---Get local slice for solve
    CALL g%get_local(self%glocs(1)%v)
    CALL u%get_local(self%ulocs(1)%v)
    !---Zero source on redundant rows
    IF(ASSOCIATED(A_native%redind))THEN
      DO i=1,A_native%nred
        self%glocs(1)%v(A_native%redind(i))=0.d0
      END DO
    END IF
    !---Solve sub-system
    self%solvers(1)%s%A=>self%A
    CALL self%solvers(1)%s%apply(self%ulocs(1),self%glocs(1))
    !---Replace local slice into solution
    CALL u%restore_local(self%ulocs(1)%v)
  ELSE
    IF(self%update_slice)CALL self%alocals(1)%update_slice
    !---Setup local vectors
    CALL bjprecond_local_vec(self%glocs(1),g%nslice,self%local_stitch(1))
    CALL bjprecond_local_vec(self%ulocs(1),u%nslice,self%local_stitch(1))
    !---Get local slice for solve
    CALL g%get_slice(self%glocs(1)%v)
    CALL u%get_slice(self%ulocs(1)%v)
    !---Solve sub-system
    self%solvers(1)%s%A=>self%alocals(1)
    CALL self%solvers(1)%s%apply(self%ulocs(1),self%glocs(1))
    !---Replace local slice into solution
    CALL u%restore_slice(self%ulocs(1)%v)
  END IF
END IF
self%update_slice=.FALSE.
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Solve the local sub-system for a single block
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE apply_block(ii)
INTEGER(i4), INTENT(in) :: ii !< Index of block to solve
!---Get local slice for solve
CALL scatter_local(self%bplans(ii),gvals,self%glocs(ii)%v,OFT_SCATTER_INSERT)
CALL scatter_local(self%bplans(ii),uinit,self%ulocs(ii)%v,OFT_SCATTER_INSERT)
!---Solve sub-system
self%solvers(ii)%s%A=>self%alocals(ii)
CALL self%solvers(ii)%s%apply(self%ulocs(ii),self%glocs(ii))
!---Replace local slice into solution (rows with non-zero weight only); additive
!   combination is done after all blocks are solved
IF(.NOT.self%additive)CALL scatter_local(self%bplans(ii),self%ulocs(ii)%v,uvals,OFT_SCATTER_INSERT,reverse=.TRUE.)
END SUBROUTINE apply_block
END SUBROUTINE bjprecond_apply
!------------------------------------------------------------------------------
!> Setup a persistent local vector for @ref oft_bjprecond, (re)allocating if needed
!------------------------------------------------------------------------------
SUBROUTINE bjprecond_local_vec(vec,n,seam)
TYPE(oft_native_vector), INTENT(inout) :: vec !< Local vector
INTEGER(i4), INTENT(in) :: n !< Number of local values
TYPE(oft_seam), TARGET, INTENT(inout) :: seam !< Local seam for vector
IF(ASSOCIATED(vec%v))THEN
  IF(SIZE(vec%v)/=n)DEALLOCATE(vec%v)
END IF
IF(.NOT.ASSOCIATED(vec%v))ALLOCATE(vec%v(n))
vec%n=n
vec%ng=INT(n,i8)
vec%nslice=n
vec%stitch_info=>seam
END SUBROUTINE bjprecond_local_vec
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine bjprecond_update(self,new_pattern)
class(oft_bjprecond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix non-zero pattern (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
self%update_slice=.TRUE.
IF(ASSOCIATED(self%solvers))THEN
  DO i=1,self%nlocal
    CALL self%solvers(i)%s%update(new_pattern)
  END DO
END IF
DEBUG_STACK_POP
end subroutine bjprecond_update
!------------------------------------------------------------------------------
!> Setup block-Jacobi smoother from XML definition
!------------------------------------------------------------------------------
subroutine bjprecond_setup_xml(self,solver_node,level)
CLASS(oft_bjprecond), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
INTEGER(i4) :: nnodes,nread
TYPE(xml_node) :: current_node
INTEGER(i4) :: val_level,ierr
INTEGER(i4), POINTER, DIMENSION(:) :: read_tmp
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---Read-in desired number of subdomains
CALL xml_get_element(solver_node,"nlocal",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,read_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `nlocal` node","bjprecond_setup_xml",__FILE__)
  IF(ASSOCIATED(read_tmp))THEN
    IF(SIZE(read_tmp)>1)THEN
      IF(SIZE(read_tmp)<val_level)CALL oft_abort("Not enough local sizes specified", &
        "bjprecond_setup_xml",__FILE__)
      self%nlocal=read_tmp(val_level)
    ELSE
      self%nlocal=read_tmp(1)
    END IF
    DEALLOCATE(read_tmp)
  END IF
END IF
IF(self%nlocal<-1)THEN
  !---Read-in desired number of subdomains
  CALL xml_get_element(solver_node,"groups",current_node,ierr)
  IF(ierr==0)THEN
    CALL xml_read_content(current_node,read_tmp,iostat=ierr)
    IF(ierr/=0)CALL oft_xml_abort("Error reading `groups` node","bjprecond_setup_xml",__FILE__)
    IF(ASSOCIATED(read_tmp))THEN
      IF(SIZE(read_tmp)>10)THEN
        CALL oft_abort("Too many slice groups specified","bjprecond_setup_xml",__FILE__)
      ELSE
        self%slice_group=read_tmp(1)
      END IF
      DEALLOCATE(read_tmp)
    END IF
  END IF
END IF
!---Read-in desired boundary overlap specification
CALL xml_get_element(solver_node,"boverlap",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%boundary_overlap,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `boverlap` node","bjprecond_setup_xml",__FILE__)
END IF
!---Read-in desired internal overlap specification
CALL xml_get_element(solver_node,"loverlap",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,self%loverlap,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `loverlap` node","bjprecond_setup_xml",__FILE__)
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'Block Jacobi solver setup:'
  WRITE(*,'(2X,A,I4)')  '- NLocal:    ',self%nlocal
  WRITE(*,'(2X,A,3X,L1)')'- Boverlap:  ',self%boundary_overlap
  WRITE(*,'(2X,A,3X,L1)')'- Loverlap:  ',self%loverlap
END IF
DEBUG_STACK_POP
end subroutine bjprecond_setup_xml
!------------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine bjprecond_delete(self,propogate)
class(oft_bjprecond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%solvers))THEN
  DO i=1,self%nlocal
    CALL self%solvers(i)%s%delete(propogate)
    DEALLOCATE(self%solvers(i)%s)
  END DO
  DEALLOCATE(self%solvers)
END IF
!---Destroy local matrices
IF(ASSOCIATED(self%alocals))THEN
  DO i=1,self%nlocal
    CALL self%alocals(i)%delete
  END DO
  DEALLOCATE(self%alocals)
END IF
IF(ASSOCIATED(self%lunion))DEALLOCATE(self%lunion)
!---Destroy block plans
IF(ASSOCIATED(self%bplans))THEN
  DO i=1,SIZE(self%bplans)
    CALL scatter_plan_delete(self%bplans(i))
  END DO
  DEALLOCATE(self%bplans)
END IF
!---Destory local partitions
IF(ASSOCIATED(self%parts))THEN
  DO i=1,self%nlocal
    IF(ASSOCIATED(self%parts(i)%v))DEALLOCATE(self%parts(i)%v)
  END DO
  DEALLOCATE(self%parts)
END IF
!---Destroy local vectors
IF(ASSOCIATED(self%ulocs))THEN
  DO i=1,SIZE(self%ulocs)
    IF(ASSOCIATED(self%ulocs(i)%v))DEALLOCATE(self%ulocs(i)%v)
    IF(ASSOCIATED(self%glocs(i)%v))DEALLOCATE(self%glocs(i)%v)
  END DO
  DEALLOCATE(self%ulocs,self%glocs)
END IF
!---Destroy local stitch information
IF(ASSOCIATED(self%local_stitch))THEN
  DO i=1,SIZE(self%local_stitch)
    IF(ASSOCIATED(self%local_stitch(i)%be))DEALLOCATE(self%local_stitch(i)%be)
  END DO
  DEALLOCATE(self%local_stitch)
END IF
!---Reset
self%nlocal=1
self%initialized=.FALSE.
DEBUG_STACK_POP
end subroutine bjprecond_delete
!------------------------------------------------------------------------------
!> Precondition a linear system using an overlapping Schwarz method
!!
!! With several blocks per process (`nlocal/=1`) or no overlap this is the
!! block-Jacobi method of @ref oft_bjprecond with blocks grown by `overlap` layers
!! within each process. Otherwise the subdomain of each process (owned rows plus
!! `overlap` layers, across processes) is solved with the sub-solver and combined
!! restricted (owned rows taken from the local subdomain, then synchronized) or
!! additively (sum over subdomains).
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE schwarz_apply(self,u,g)
CLASS(oft_schwarz_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_vector), INTENT(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), INTENT(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_native_matrix), POINTER :: A_native
CLASS(oft_solver), POINTER :: pretmp
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: uvals,gvals
TYPE(oft_timer) :: stimer
INTEGER(i8) :: nsent0
DEBUG_STACK_PUSH
IF(self%overlap<0)CALL oft_abort('Negative overlap','schwarz_apply',__FILE__)
!---Blocks grown by layers of cells (may cross processes)
IF(self%cell_mode)THEN
  IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
    'schwarz_apply',__FILE__)
  CALL schwarz_cell_apply(self,A_native,u,g)
  DEBUG_STACK_POP
  RETURN
END IF
!---Blocks within each process (or no overlap): block-Jacobi with grown blocks
IF((self%nlocal/=1).OR.(self%overlap==0).OR.self%boundary_overlap)THEN
  self%overlap_levels=self%overlap
  CALL bjprecond_apply(self,u,g)
  DEBUG_STACK_POP
  RETURN
END IF
!---One subdomain per process, overlapping across processes
IF(.NOT.native_matrix_cast(A_native,self%A))CALL oft_abort('Native matrix required', &
  'schwarz_apply',__FILE__)
IF(.NOT.self%initialized)THEN
  IF(.NOT.ASSOCIATED(self%pre))CALL oft_abort('No sub-solver (preconditioner defined)', &
    'schwarz_apply',__FILE__)
  CALL stimer%tick
  nsent0=oft_exchange_nsent
  CALL solver_setup(self)
  ALLOCATE(self%ovl)
  CALL native_overlap_setup(self%ovl,A_native,u,self%overlap)
  pretmp=>self%pre
  ALLOCATE(self%solvers(1))
  ALLOCATE(self%solvers(1)%s, source=pretmp)
  IF(ASSOCIATED(pretmp%pre))ALLOCATE(self%solvers(1)%s%pre, source=pretmp%pre)
  ALLOCATE(self%sub_seam,self%usub,self%gsub)
  self%sub_seam%skip=.TRUE.
  CALL bjprecond_local_vec(self%usub,self%ovl%n,self%sub_seam)
  CALL bjprecond_local_vec(self%gsub,self%ovl%n,self%sub_seam)
  self%update_slice=.FALSE.
  CALL oft_exchange_report('Schwarz graph mode setup (total)',stimer,nsent0)
ELSE IF(self%update_slice)THEN
  CALL native_overlap_update(self%ovl)
  self%update_slice=.FALSE.
END IF
!---Gather subdomain values, solve, combine
NULLIFY(uvals,gvals)
CALL u%borrow_local(uvals,OFT_ACCESS_READWRITE)
CALL g%borrow_local(gvals,OFT_ACCESS_READ)
CALL native_overlap_gather(self%ovl,gvals,self%gsub%v)
CALL native_overlap_gather(self%ovl,uvals,self%usub%v)
self%solvers(1)%s%A=>self%ovl%amat
CALL self%solvers(1)%s%apply(self%usub,self%gsub)
CALL native_overlap_scatter(self%ovl,self%usub%v,uvals,self%additive,weighted=self%weighted)
!---Restricted: copies of owned rows are synchronized; additive: already consistent
CALL u%release_local(uvals,OFT_ACCESS_READWRITE,sync=(.NOT.self%additive))
CALL g%release_local(gvals,OFT_ACCESS_READ,sync=.FALSE.)
DEBUG_STACK_POP
END SUBROUTINE schwarz_apply
!------------------------------------------------------------------------------
!> Set up the blocks of the Schwarz preconditioner in cell mode (collective)
!!
!! The cell overlap builder attached to the matrix (`A%cell_ovl`) gives the DoFs of
!! each block (core cells grown by `overlap` layers of cells adjacent by `adjacency`)
!! and the DoFs assigned to each block for the restricted write-back. Their union on
!! this process is set up as one subdomain across processes
!! (@ref oft_native_overlap_mod::native_overlap_setup_dofs "native_overlap_setup_dofs",
!! one row request round) and each block is the submatrix of the subdomain matrix
!! on its rows, with a copy of the sub-solver `pre` and a local plan between the
!! subdomain and block vectors. Plan weights are 1 on the rows assigned to the block
!! and 0 elsewhere (restricted) or 1 everywhere (additive). For the weighted additive
!! combination the partition of unity weights count the blocks containing each DoF
!! over all processes (@ref oft_native_overlap_mod::native_overlap_pu_setup
!! "native_overlap_pu_setup").
!------------------------------------------------------------------------------
SUBROUTINE schwarz_cell_setup(self,A_native,u)
CLASS(oft_schwarz_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), POINTER, INTENT(in) :: A_native !< Matrix (`self%A` as native matrix)
CLASS(oft_vector), INTENT(inout) :: u !< Vector of the matrix row space
INTEGER(i4) :: i,j,nb,n
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,basg
REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt,mult
TYPE(oft_native_vector) :: uskip
CLASS(oft_solver), POINTER :: pretmp
TYPE(oft_timer) :: stimer
INTEGER(i8) :: nsent0
DEBUG_STACK_PUSH
CALL stimer%tick
nsent0=oft_exchange_nsent
IF(.NOT.ASSOCIATED(self%pre))CALL oft_abort('No sub-solver (preconditioner defined)', &
  'schwarz_cell_setup',__FILE__)
IF(.NOT.ASSOCIATED(self%A%cell_ovl))CALL oft_abort('Cell mode requires a cell overlap builder '// &
  'attached to the matrix ("cell_ovl")','schwarz_cell_setup',__FILE__)
IF(self%boundary_overlap)CALL oft_abort('"boverlap" is not supported in cell mode', &
  'schwarz_cell_setup',__FILE__)
IF(self%nlocal<-1)CALL oft_abort('Slice groups are not supported in cell mode', &
  'schwarz_cell_setup',__FILE__)
IF((self%adjacency<OFT_CELL_ADJ_POINT).OR.(self%adjacency>OFT_CELL_ADJ_FACE))CALL oft_abort( &
  'Invalid cell adjacency','schwarz_cell_setup',__FILE__)
CALL solver_setup(self)
!---Number of blocks is set by the builder
nb=self%A%cell_ovl%nblocks()
IF(nb<1)CALL oft_abort('Cell overlap builder has no blocks','schwarz_cell_setup',__FILE__)
IF(oft_env%head_proc.AND.(self%nlocal>1).AND.(self%nlocal/=nb))WRITE(*,'(A,X,I4)') &
  ' Updating number of local parts to match cell overlap blocks',nb
self%nlocal=nb
!---Block DoF sets and their union as one subdomain across processes
ALLOCATE(self%cdofs)
CALL self%A%cell_ovl%build(self%A,u,self%overlap,self%adjacency,self%base_layers,self%cdofs)
IF(self%cdofs%nblocks/=nb)CALL oft_abort('Block count mismatch with cell overlap builder', &
  'schwarz_cell_setup',__FILE__)
ALLOCATE(self%ovl)
CALL native_overlap_setup_dofs(self%ovl,A_native,u,self%cdofs%n,self%cdofs%gdof,self%cdofs%odof)
IF(self%ovl%n/=self%cdofs%n)CALL oft_abort('Subdomain size mismatch','schwarz_cell_setup',__FILE__)
!---Subdomain index of each union DoF and block assigned to each subdomain DoF (0: none)
ALLOCATE(pmap(self%cdofs%n),basg(self%ovl%n))
CALL native_overlap_index(self%ovl,self%cdofs%n,self%cdofs%gdof,pmap)
IF(ANY(pmap<1))CALL oft_abort('Union DoF missing from subdomain','schwarz_cell_setup',__FILE__)
basg=0
DO i=1,nb
  DO j=1,self%cdofs%asg(i)%n
    IF(basg(pmap(self%cdofs%asg(i)%v(j)))/=0)CALL oft_abort('DoF assigned to several blocks', &
      'schwarz_cell_setup',__FILE__)
    basg(pmap(self%cdofs%asg(i)%v(j)))=i
  END DO
END DO
!---Subdomain vectors (no stitching)
ALLOCATE(self%sub_seam,self%usub,self%gsub,self%usol)
self%sub_seam%skip=.TRUE.
CALL bjprecond_local_vec(self%usub,self%ovl%n,self%sub_seam)
CALL bjprecond_local_vec(self%gsub,self%ovl%n,self%sub_seam)
CALL bjprecond_local_vec(self%usol,self%ovl%n,self%sub_seam)
uskip%n=self%ovl%n
uskip%stitch_info=>self%sub_seam
!---Blocks: subdomain rows, submatrix, sub-solver, plan and local vectors
ALLOCATE(self%local_stitch(1))
self%local_stitch(1)%skip=.TRUE.
ALLOCATE(self%parts(nb),self%alocals(nb),self%solvers(nb),self%bplans(nb))
ALLOCATE(self%ulocs(nb),self%glocs(nb))
pretmp=>self%pre
DO i=1,nb
  n=self%cdofs%blk(i)%n
  self%parts(i)%n=n
  ALLOCATE(self%parts(i)%v(n))
  DO j=1,n
    self%parts(i)%v(j)=pmap(self%cdofs%blk(i)%v(j))
  END DO
  CALL sort_array(self%parts(i)%v,n)
  ALLOCATE(self%solvers(i)%s, source=pretmp)
  IF(ASSOCIATED(pretmp%pre))ALLOCATE(self%solvers(i)%s%pre, source=pretmp%pre)
  CALL bjprecond_local_vec(self%ulocs(i),n,self%local_stitch(1))
  CALL bjprecond_local_vec(self%glocs(i),n,self%local_stitch(1))
  IF(n==0)CYCLE ! Empty block (no cells): skipped in apply
  CALL self%alocals(i)%setup(self%ovl%amat,uskip,rows=self%parts(i)%v)
  ALLOCATE(wt(n))
  IF(self%additive)THEN
    wt=1.d0
  ELSE
    DO j=1,n
      wt(j)=MERGE(1.d0,0.d0,basg(self%parts(i)%v(j))==i)
    END DO
  END IF
  CALL scatter_plan_setup_local(self%bplans(i),n,self%parts(i)%v,wt=wt)
  DEALLOCATE(wt)
END DO
DEALLOCATE(pmap,basg)
!---Partition of unity weights: blocks containing each DoF, over all processes
IF(self%additive.AND.self%weighted)THEN
  ALLOCATE(mult(MAX(1,self%ovl%n)))
  mult=0.d0
  DO i=1,nb
    DO j=1,self%parts(i)%n
      mult(self%parts(i)%v(j))=mult(self%parts(i)%v(j))+1.d0
    END DO
  END DO
  CALL native_overlap_pu_setup(self%ovl,mult)
  DEALLOCATE(mult)
END IF
self%update_slice=.FALSE.
CALL oft_exchange_report('Schwarz cell mode setup (total)',stimer,nsent0)
DEBUG_STACK_POP
END SUBROUTINE schwarz_cell_setup
!------------------------------------------------------------------------------
!> Precondition a linear system using the Schwarz method with cell-grown blocks
!!
!! Blocks are set up on first use (@ref schwarz_cell_setup) and their matrices
!! refreshed after an update (subdomain matrix, then block submatrices). Each
!! application gathers the right-hand side and guess on the subdomain (union of the
!! blocks of this process), solves every block from the subdomain values and
!! combines the block solutions on the subdomain:
!! - restricted: each block writes only the rows assigned to it (0/1 plan weights).
!!   These are disjoint between blocks, so blocks may be solved concurrently; other
!!   rows keep the gathered guess. Owned rows are then written back and copies
!!   synchronized from their owners.
!! - additive: block solutions are summed on the subdomain (serially, in block order,
!!   as blocks overlap), then summed over subdomains of all processes.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE schwarz_cell_apply(self,A_native,u,g)
CLASS(oft_schwarz_precond), INTENT(inout) :: self !< Solver object
CLASS(oft_native_matrix), POINTER, INTENT(in) :: A_native !< Matrix (`self%A` as native matrix)
CLASS(oft_vector), INTENT(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), INTENT(inout) :: g !< RHS (input), Residual (output)
LOGICAL :: thread_safe
INTEGER(i4) :: i
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: uvals,gvals
DEBUG_STACK_PUSH
IF(.NOT.self%initialized)THEN
  CALL schwarz_cell_setup(self,A_native,u)
ELSE IF(self%update_slice)THEN
  !---Subdomain matrix (and its full local matrix) first, so block updates only copy values
  CALL native_overlap_update(self%ovl)
  !$omp parallel do schedule(dynamic,1)
  DO i=1,self%nlocal
    IF(self%parts(i)%n>0)CALL self%alocals(i)%update_slice
  END DO
  self%update_slice=.FALSE.
END IF
!---Block solves run concurrently only for thread safe sub-solvers outside parallel regions
thread_safe=self%pre%check_thread()
!---Gather subdomain values (borrowed local values, no copy for native vectors)
NULLIFY(uvals,gvals)
CALL u%borrow_local(uvals,OFT_ACCESS_READWRITE)
CALL g%borrow_local(gvals,OFT_ACCESS_READ)
CALL native_overlap_gather(self%ovl,gvals,self%gsub%v)
CALL native_overlap_gather(self%ovl,uvals,self%usub%v)
IF(self%additive)THEN
  IF(self%ovl%n>0)self%usol%v=0.d0
ELSE
  IF(self%ovl%n>0)self%usol%v=self%usub%v
END IF
!---Blocks only read the gathered subdomain values and write disjoint rows of `usol`
!   (restricted), so the result does not depend on the order blocks are solved in.
!   Only open a parallel region if the sub-solver is thread safe and we are not already
!   in an active parallel region (see @ref bjprecond_apply). Each block (and its local
!   vectors) is only accessed by one thread.
IF(thread_safe.AND.(.NOT.omp_in_parallel()))THEN
  !$omp parallel do schedule(dynamic,1)
  DO i=1,self%nlocal
    CALL apply_cell_block(i)
  END DO
ELSE
  DO i=1,self%nlocal
    CALL apply_cell_block(i)
  END DO
END IF
!---Additive combination: sum of the block solutions (blocks in order)
IF(self%additive)THEN
  DO i=1,self%nlocal
    IF(self%parts(i)%n==0)CYCLE
    CALL scatter_local(self%bplans(i),self%ulocs(i)%v,self%usol%v,OFT_SCATTER_ADD,reverse=.TRUE.)
  END DO
END IF
CALL native_overlap_scatter(self%ovl,self%usol%v,uvals,self%additive,weighted=self%weighted)
!---Restricted: copies of owned rows are synchronized; additive: already consistent
CALL u%release_local(uvals,OFT_ACCESS_READWRITE,sync=(.NOT.self%additive))
CALL g%release_local(gvals,OFT_ACCESS_READ,sync=.FALSE.)
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Solve the local sub-system of a single cell block
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE apply_cell_block(ii)
INTEGER(i4), INTENT(in) :: ii !< Index of block to solve
IF(self%parts(ii)%n==0)RETURN
CALL scatter_local(self%bplans(ii),self%gsub%v,self%glocs(ii)%v,OFT_SCATTER_INSERT)
CALL scatter_local(self%bplans(ii),self%usub%v,self%ulocs(ii)%v,OFT_SCATTER_INSERT)
self%solvers(ii)%s%A=>self%alocals(ii)
CALL self%solvers(ii)%s%apply(self%ulocs(ii),self%glocs(ii))
IF(.NOT.self%additive)CALL scatter_local(self%bplans(ii),self%ulocs(ii)%v,self%usol%v, &
  OFT_SCATTER_INSERT,reverse=.TRUE.)
END SUBROUTINE apply_cell_block
END SUBROUTINE schwarz_cell_apply
!------------------------------------------------------------------------------
!> Update Schwarz preconditioner after changing settings/operators
!!
!! As @ref bjprecond_update (numeric update on next application), except that in
!! cell mode a new matrix pattern (`new_pattern=.TRUE.`) frees the blocks and
!! subdomain, which are rebuilt on the next application.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE schwarz_update(self,new_pattern)
CLASS(oft_schwarz_precond), INTENT(inout) :: self !< Solver object
LOGICAL, OPTIONAL, INTENT(in) :: new_pattern !< Update matrix non-zero pattern (optional)
LOGICAL :: rebuild
DEBUG_STACK_PUSH
rebuild=.FALSE.
IF(PRESENT(new_pattern))rebuild=new_pattern
IF(self%cell_mode.AND.rebuild.AND.ASSOCIATED(self%cdofs))THEN
  CALL schwarz_delete(self)
ELSE
  CALL bjprecond_update(self,new_pattern)
END IF
DEBUG_STACK_POP
END SUBROUTINE schwarz_update
!------------------------------------------------------------------------------
!> Setup Schwarz preconditioner from XML definition
!!
!! Reads the block-Jacobi options (`nlocal`, `groups`; see @ref bjprecond_setup_xml),
!! `overlap` (number of overlap layers, default 1), `type` (`restricted`, default,
!! `additive`, or `weighted` for additive with partition of unity weights),
!! `overlap_mode` (`graph`, default, or `cell`), and for cell mode `adjacency`
!! (`point`, default, `edge` or `face`) and `base_layers` (default -1, builder
!! default). `overlap` and `base_layers` take one value or one per level.
!------------------------------------------------------------------------------
subroutine schwarz_setup_xml(self,solver_node,level)
CLASS(oft_schwarz_precond), INTENT(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
TYPE(xml_node) :: current_node
INTEGER(i4) :: ierr,val_level
INTEGER(i4), POINTER, DIMENSION(:) :: read_tmp
CHARACTER(LEN=:), ALLOCATABLE :: stype
DEBUG_STACK_PUSH
val_level=1
IF(PRESENT(level))val_level=level
CALL bjprecond_setup_xml(self,solver_node,level)
IF(self%loverlap)CALL oft_abort('Use "overlap" instead of "loverlap" for the Schwarz preconditioner', &
  'schwarz_setup_xml',__FILE__)
!---Number of overlap layers (one value, or one per level)
CALL xml_get_element(solver_node,"overlap",current_node,ierr)
IF(ierr==0)THEN
  NULLIFY(read_tmp)
  CALL xml_read_content(current_node,read_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `overlap` node","schwarz_setup_xml",__FILE__)
  IF(ASSOCIATED(read_tmp))THEN
    IF(SIZE(read_tmp)>1)THEN
      IF(SIZE(read_tmp)<val_level)CALL oft_abort("Not enough overlap values specified", &
        "schwarz_setup_xml",__FILE__)
      self%overlap=read_tmp(val_level)
    ELSE
      self%overlap=read_tmp(1)
    END IF
    DEALLOCATE(read_tmp)
  END IF
END IF
!---Combination type
CALL xml_get_element(solver_node,"type",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,stype,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `type` node","schwarz_setup_xml",__FILE__)
  SELECT CASE(TRIM(stype))
  CASE("restricted")
    self%additive=.FALSE.
    self%weighted=.FALSE.
  CASE("additive")
    self%additive=.TRUE.
    self%weighted=.FALSE.
  CASE("weighted")
    self%additive=.TRUE.
    self%weighted=.TRUE.
  CASE DEFAULT
    CALL oft_abort('Unknown Schwarz type "'//TRIM(stype)//'"','schwarz_setup_xml',__FILE__)
  END SELECT
END IF
!---Overlap mode (graph or cell layers)
CALL xml_get_element(solver_node,"overlap_mode",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,stype,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `overlap_mode` node","schwarz_setup_xml",__FILE__)
  SELECT CASE(TRIM(stype))
  CASE("graph")
    self%cell_mode=.FALSE.
  CASE("cell")
    self%cell_mode=.TRUE.
  CASE DEFAULT
    CALL oft_abort('Unknown overlap mode "'//TRIM(stype)//'"','schwarz_setup_xml',__FILE__)
  END SELECT
END IF
!---Cell adjacency (cell mode)
CALL xml_get_element(solver_node,"adjacency",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_read_content(current_node,stype,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `adjacency` node","schwarz_setup_xml",__FILE__)
  SELECT CASE(TRIM(stype))
  CASE("point")
    self%adjacency=OFT_CELL_ADJ_POINT
  CASE("edge")
    self%adjacency=OFT_CELL_ADJ_EDGE
  CASE("face")
    self%adjacency=OFT_CELL_ADJ_FACE
  CASE DEFAULT
    CALL oft_abort('Unknown cell adjacency "'//TRIM(stype)//'"','schwarz_setup_xml',__FILE__)
  END SELECT
END IF
!---Base mesh layers bounding the halo (cell mode; one value, or one per level)
CALL xml_get_element(solver_node,"base_layers",current_node,ierr)
IF(ierr==0)THEN
  NULLIFY(read_tmp)
  CALL xml_read_content(current_node,read_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort("Error reading `base_layers` node","schwarz_setup_xml",__FILE__)
  IF(ASSOCIATED(read_tmp))THEN
    IF(SIZE(read_tmp)>1)THEN
      IF(SIZE(read_tmp)<val_level)CALL oft_abort("Not enough base_layers values specified", &
        "schwarz_setup_xml",__FILE__)
      self%base_layers=read_tmp(val_level)
    ELSE
      self%base_layers=read_tmp(1)
    END IF
    DEALLOCATE(read_tmp)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'Schwarz solver setup:'
  WRITE(*,'(2X,A,I4)')  '- Overlap:   ',self%overlap
  WRITE(*,'(2X,A,3X,L1)')'- Additive:  ',self%additive
  WRITE(*,'(2X,A,3X,L1)')'- Weighted:  ',self%weighted
  WRITE(*,'(2X,A,3X,L1)')'- Cell mode: ',self%cell_mode
  IF(self%cell_mode)THEN
    WRITE(*,'(2X,A,I4)')  '- Adjacency: ',self%adjacency
    WRITE(*,'(2X,A,I4)')  '- Base lays: ',self%base_layers
  END IF
END IF
DEBUG_STACK_POP
end subroutine schwarz_setup_xml
!------------------------------------------------------------------------------
!> Destroy Schwarz preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine schwarz_delete(self,propogate)
class(oft_schwarz_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
DEBUG_STACK_PUSH
!---Cell mode: blocks (submatrices of the subdomain matrix) before the subdomain
IF(ASSOCIATED(self%cdofs))THEN
  CALL bjprecond_delete(self,propogate)
  CALL overlap_dofs_delete(self%cdofs)
  DEALLOCATE(self%cdofs)
END IF
IF(ASSOCIATED(self%usol))THEN
  IF(ASSOCIATED(self%usol%v))DEALLOCATE(self%usol%v)
  DEALLOCATE(self%usol)
END IF
IF(ASSOCIATED(self%ovl))THEN
  CALL native_overlap_delete(self%ovl)
  DEALLOCATE(self%ovl)
  !---Single sub-solver of the cross-process subdomain
  IF(ASSOCIATED(self%solvers))THEN
    CALL self%solvers(1)%s%delete(propogate)
    DEALLOCATE(self%solvers(1)%s)
    DEALLOCATE(self%solvers)
  END IF
END IF
IF(ASSOCIATED(self%usub))THEN
  IF(ASSOCIATED(self%usub%v))DEALLOCATE(self%usub%v)
  IF(ASSOCIATED(self%gsub%v))DEALLOCATE(self%gsub%v)
  DEALLOCATE(self%usub,self%gsub,self%sub_seam)
END IF
CALL bjprecond_delete(self,propogate)
self%overlap_levels=-1
DEBUG_STACK_POP
end subroutine schwarz_delete
!------------------------------------------------------------------------------
!> Setup Block-Jacobi preconditioner for a native dense complex matrix
!!
!! Partitions `graph` into `self%nlocal` blocks using METIS, optionally extends
!! each block by one level of graph neighbors (`self%loverlap`), and explicitly
!! inverts the corresponding diagonal blocks of `A`. If `self%nlocal==1` a single
!! block containing the full matrix is used. Any existing factorization is replaced.
!------------------------------------------------------------------------------
SUBROUTINE dense_cbjprecond_setup(self,A,graph)
CLASS(oft_native_dense_cbjprecond), INTENT(inout) :: self !< Solver object
CLASS(oft_cmatrix), TARGET, INTENT(inout) :: A !< Matrix to precondition (must be @ref oft_native_dense_cmatrix)
TYPE(oft_graph), INTENT(inout) :: graph !< Graph used to partition the matrix into blocks [A%nr]
!---
INTEGER(i4) :: i,j,k,jr,n,info,nfail,max_size
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: part
LOGICAL, ALLOCATABLE, DIMENSION(:) :: eflag
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Amat
DEBUG_STACK_PUSH
SELECT TYPE(A)
CLASS IS(oft_native_dense_cmatrix)
  Amat=>A%M
CLASS DEFAULT
  CALL oft_abort('"A" is not a native dense complex matrix','dense_cbjprecond_setup',__FILE__)
END SELECT
IF(A%nr/=A%nc)CALL oft_abort('Matrix is not square','dense_cbjprecond_setup',__FILE__)
IF(graph%nr/=A%nr)CALL oft_abort('Graph and matrix sizes do not match','dense_cbjprecond_setup',__FILE__)
IF(self%nlocal<1)CALL oft_abort('"nlocal" must be >=1','dense_cbjprecond_setup',__FILE__)
CALL dense_cbjprecond_delete(self)
self%A=>A
IF(oft_env%pm.AND.oft_env%head_proc)WRITE(*,'(2A,I6)')oft_indent,'nLocal = ',self%nlocal
!---Partition graph
ALLOCATE(part(A%nr))
IF(self%nlocal>1)THEN
  CALL partition_graph(graph,self%nlocal,part)
  IF(ANY(part<1.OR.part>self%nlocal))CALL oft_abort('Unowned row after partitioning', &
    'dense_cbjprecond_setup',__FILE__)
ELSE
  part=1
END IF
ALLOCATE(self%blocks(self%nlocal))
IF(self%loverlap)THEN
  !$omp parallel private(j,jr,k,eflag)
  ALLOCATE(eflag(A%nr))
  !$omp do
  DO i=1,self%nlocal
    !---Flag owned rows and their graph neighbors
    eflag=(part==i)
    DO jr=1,A%nr
      IF(part(jr)/=i)CYCLE
      DO k=graph%kr(jr),graph%kr(jr+1)-1
        eflag(graph%lc(k))=.TRUE.
      END DO
    END DO
    !---Overlap rows are marked with negative indices (input only)
    self%blocks(i)%n=COUNT(eflag)
    ALLOCATE(self%blocks(i)%ind(self%blocks(i)%n))
    k=0
    DO j=1,A%nr
      IF(.NOT.eflag(j))CYCLE
      k=k+1
      IF(part(j)==i)THEN
        self%blocks(i)%ind(k)=j
      ELSE
        self%blocks(i)%ind(k)=-j
      END IF
    END DO
  END DO
  DEALLOCATE(eflag)
  !$omp end parallel
ELSE
  DO i=1,A%nr
    self%blocks(part(i))%n=self%blocks(part(i))%n+1
  END DO
  DO i=1,self%nlocal
    ALLOCATE(self%blocks(i)%ind(self%blocks(i)%n))
    self%blocks(i)%n=0
  END DO
  DO i=1,A%nr
    self%blocks(part(i))%n=self%blocks(part(i))%n+1
    self%blocks(part(i))%ind(self%blocks(part(i))%n)=i
  END DO
END IF
DEALLOCATE(part)
!---Extract and invert diagonal blocks
max_size=0
nfail=0
!$omp parallel do private(j,k,n,info) schedule(static) reduction(+:nfail) reduction(max:max_size)
DO i=1,self%nlocal
  n=self%blocks(i)%n
  max_size=MAX(max_size,n)
  IF(n==0)CYCLE
  ALLOCATE(self%blocks(i)%inv(n,n))
  DO k=1,n
    DO j=1,n
      self%blocks(i)%inv(j,k)=Amat(ABS(self%blocks(i)%ind(j)),ABS(self%blocks(i)%ind(k)))
    END DO
  END DO
  CALL lapack_matinv(n,self%blocks(i)%inv,info)
  IF(info/=0)nfail=nfail+1
END DO
IF(nfail>0)CALL oft_abort('Block inversion failed','dense_cbjprecond_setup',__FILE__)
self%max_block_size=max_size
self%initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE dense_cbjprecond_setup
!------------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method with native dense
!! complex matrices
!!
!! Computes \f$ u = M^{-1} g \f$, where \f$ M \f$ is the block-diagonal
!! approximation of the matrix. The input value of `u` is ignored and `g` is
!! unchanged.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE dense_cbjprecond_apply(self,u,g)
CLASS(oft_native_dense_cbjprecond), INTENT(inout) :: self !< Solver object
CLASS(oft_cvector), INTENT(inout) :: u !< Result \f$ M^{-1} g \f$ (output)
CLASS(oft_cvector), INTENT(inout) :: g !< Input vector (unchanged)
!---
INTEGER(i4) :: i,j,n
COMPLEX(c8), POINTER, CONTIGUOUS, DIMENSION(:) :: uvals,gvals
COMPLEX(c8), ALLOCATABLE, DIMENSION(:,:) :: vloc
DEBUG_STACK_PUSH
IF(.NOT.self%initialized)CALL oft_abort('Preconditioner has not been setup','dense_cbjprecond_apply',__FILE__)
IF(g%n/=self%A%nr)CALL oft_abort('Row mismatch','dense_cbjprecond_apply',__FILE__)
IF(u%n/=self%A%nc)CALL oft_abort('Col mismatch','dense_cbjprecond_apply',__FILE__)
IF(self%nlocal/=SIZE(self%blocks))CALL oft_abort('"nlocal" changed since setup','dense_cbjprecond_apply',__FILE__)
NULLIFY(uvals,gvals)
CALL u%borrow_local(uvals,OFT_ACCESS_WRITE)
CALL g%borrow_local(gvals,OFT_ACCESS_READ)
uvals=(0.d0,0.d0)
!$omp parallel private(vloc,i,j,n)
ALLOCATE(vloc(self%max_block_size,2))
!$omp do schedule(static)
DO i=1,self%nlocal
  n=self%blocks(i)%n
  IF(n==0)CYCLE
  DO j=1,n
    vloc(j,1)=gvals(ABS(self%blocks(i)%ind(j)))
  END DO
  CALL zgemv('N',n,n,(1.d0,0.d0),self%blocks(i)%inv,n,vloc(:,1),1,(0.d0,0.d0),vloc(:,2),1)
  DO j=1,n
    IF(self%blocks(i)%ind(j)>0)uvals(self%blocks(i)%ind(j))=vloc(j,2)
  END DO
END DO
!$omp end do
DEALLOCATE(vloc)
!$omp end parallel
!---Block solves on shared rows differ between processes, so synchronize from owner
CALL u%release_local(uvals,OFT_ACCESS_WRITE,sync=.TRUE.)
CALL g%release_local(gvals,OFT_ACCESS_READ,sync=.FALSE.)
DEBUG_STACK_POP
END SUBROUTINE dense_cbjprecond_apply
!------------------------------------------------------------------------------
!> Destroy dense Block-Jacobi preconditioner and deallocate all internal storage
!!
!! @note Settings (`nlocal`, `loverlap`) are retained
!------------------------------------------------------------------------------
SUBROUTINE dense_cbjprecond_delete(self,propogate)
CLASS(oft_native_dense_cbjprecond), INTENT(inout) :: self !< Solver object
LOGICAL, OPTIONAL, INTENT(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%blocks))THEN
  DO i=1,SIZE(self%blocks)
    IF(ASSOCIATED(self%blocks(i)%ind))DEALLOCATE(self%blocks(i)%ind)
    IF(ASSOCIATED(self%blocks(i)%inv))DEALLOCATE(self%blocks(i)%inv)
  END DO
  DEALLOCATE(self%blocks)
END IF
self%max_block_size=0
NULLIFY(self%A)
self%initialized=.FALSE.
CALL csolver_delete(self,propogate)
DEBUG_STACK_POP
END SUBROUTINE dense_cbjprecond_delete
end module oft_native_solvers
