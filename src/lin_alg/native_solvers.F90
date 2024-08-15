!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
MODULE oft_native_solvers
USE oft_local
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_cvector_ptr, &
  oft_matrix, oft_cmatrix, oft_graph
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, &
  oft_native_matrix, native_matrix_cast, oft_native_submatrix, &
  partition_graph, native_matrix_setup_full, oft_native_submatrix
USE oft_solver_base, ONLY: oft_solver, oft_bc_proto, oft_orthog, oft_solver_ptr, &
  oft_csolver, oft_cbc_proto, oft_corthog, oft_eigsolver, solver_setup, csolver_setup, &
  eigsolver_setup
IMPLICIT NONE
#include "local.h"
private
!---------------------------------------------------------------------------
! TYPE oft_native_cg_solver
!---------------------------------------------------------------------------
!> CG solver class
!!
!! @note Matrix must be SPD, otherwise solver will fail.
!---------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_native_cg_solver
  logical :: precond=.FALSE. !< Solver as preconditioner
  class(oft_orthog), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_orthog), pointer :: cleaner => NULL() !< Null-Space cleaner
contains
  !> Solve system
  procedure :: apply => cg_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => cg_setup_xml
  !> Clean-up internal storage
  PROCEDURE :: delete => cg_delete
end type oft_native_cg_solver
!---------------------------------------------------------------------------
! TYPE oft_native_gmres_solver
!---------------------------------------------------------------------------
!> GMRES solver class
!---------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_native_gmres_solver
  logical :: precond = .FALSE. !< Flag to indicate use as a preconditioner (not used)
  integer(i4) :: nrits=10 !< Convergence specification
  class(oft_orthog), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_orthog), pointer :: cleaner => NULL() !< Null-Space cleaner
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
!---------------------------------------------------------------------------
! TYPE oft_native_gmres_csolver
!---------------------------------------------------------------------------
!> GMRES solver class
!---------------------------------------------------------------------------
type, public, extends(oft_csolver) :: oft_native_gmres_csolver
  logical :: precond = .FALSE. !< Flag to indicate use as a preconditioner (not used)
  integer(i4) :: nrits=10 !< Convergence specification
  class(oft_corthog), pointer :: orthog => NULL() !< Orthogonalization
  class(oft_corthog), pointer :: cleaner => NULL() !< Null-Space cleaner
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
!---------------------------------------------------------------------------
! TYPE oft_native_cg_eigsolver
!---------------------------------------------------------------------------
!> CG eigensolver class
!!
!! @note Matrix system must be SPD, otherwise solver will fail.
!---------------------------------------------------------------------------
type, public, extends(oft_eigsolver) :: oft_native_cg_eigsolver
  integer(i4) :: nrestarts = 2
  integer(i4) :: ninner = -1
  !> Boundary condition
  procedure(oft_bc_proto), pointer, nopass :: bc => NULL()
  class(oft_orthog), pointer :: orthog => NULL() !< Orthogonalization
contains
  !> Solve system
  procedure :: apply => cg_eigsolver_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => cg_eig_delete
end type oft_native_cg_eigsolver
!---------------------------------------------------------------------------
! TYPE oft_nksolver
!---------------------------------------------------------------------------
!> Native Newton solver
!!
!! Non-linear system solver based on Newton iteration for the system,
!! \f$ M_j(x_i) = y_j \f$, where the Jacobian is defined as
!! \f$ J_{ij} = \frac{\partial M_j}{\partial x_i} \f$. The non-linear update then
!! takes the form \f$ x^{n+1} = x^n - J^{-1}(x^n) \cdot M(x^n) \f$.
!---------------------------------------------------------------------------
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
  PROCEDURE(oft_update_jacobian), POINTER, NOPASS :: J_update => NULL() !< Jacobian update subroutine
  PROCEDURE(oft_bc_proto), POINTER, NOPASS :: bc => NULL() !< Boundary condition
CONTAINS
  !> Solve non-linear system
  PROCEDURE :: apply => nksolver_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => nksolver_delete
END TYPE oft_nksolver
!---------------------------------------------------------------------------
! TYPE oft_identity_inv
!---------------------------------------------------------------------------
!> Identity matrix inversion
!---------------------------------------------------------------------------
TYPE, PUBLIC, extends(oft_solver) :: oft_identity_inv
CONTAINS
  !> Apply preconditioner
  PROCEDURE :: apply => identity_inv_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => identity_inv_delete
END TYPE oft_identity_inv
!---------------------------------------------------------------------------
! TYPE oft_diag_scale
!---------------------------------------------------------------------------
!> Diagonal preconditioner
!---------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_diag_scale
contains
  !> Apply preconditioner
  procedure :: apply => diag_scale_apply
  !> Clean-up internal smoother storage
  procedure :: delete => diag_scale_delete
end type oft_diag_scale
!---------------------------------------------------------------------------
! TYPE oft_diag_cscale
!---------------------------------------------------------------------------
!> Diagonal preconditioner
!---------------------------------------------------------------------------
type, public, extends(oft_csolver) :: oft_diag_cscale
contains
  !> Apply preconditioner
  procedure :: apply => cdiag_scale_apply
  !> Clean-up internal smoother storage
  procedure :: delete => cdiag_scale_delete
end type oft_diag_cscale
!---------------------------------------------------------------------------
! TYPE oft_jblock_precond
!---------------------------------------------------------------------------
!> Symmetric Jacobi smoother
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
! TYPE oft_ml_precond
!---------------------------------------------------------------------------
!> Multi-level preconditioner level
!!
!! @note This class acts as a driver/container for separate solver objects
!! as a result any linear algebra backend may be used with this type.
!---------------------------------------------------------------------------
type, public, extends(oft_solver) :: oft_ml_precond
  logical :: symmetric = .FALSE. !< Symmetric flag
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
  !> Interpolate field to high level
  procedure(oft_interp_proto), pointer, nopass :: interp => NULL()
  !> Inject field to lower level
  procedure(oft_interp_proto), pointer, nopass :: inject => NULL()
  !> Create a new field on a specified level
  procedure(oft_veccreate_proto), pointer, nopass :: vec_create => NULL()
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
!---------------------------------------------------------------------------
! TYPE oft_ml_trans
!---------------------------------------------------------------------------
!> Multi-level transfer level
!---------------------------------------------------------------------------
type, public, extends(oft_ml_precond) :: oft_ml_trans
contains
  !> Transfer field between Base and MPI levels
  procedure :: apply => ml_trans_apply
  !> Clean-up internal smoother storage
  procedure :: delete => ml_trans_delete
end type oft_ml_trans
!---------------------------------------------------------------------------
! TYPE oft_bjprecond
!---------------------------------------------------------------------------
!> Block-Jacobi preconditioner
!---------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_bjprecond
  LOGICAL :: update_slice = .FALSE. !< Update local matrices on next application
  LOGICAL :: boundary_overlap = .FALSE. !< Use Additive-Schwarz with boundary elements
  LOGICAL :: loverlap = .FALSE. !< Use Additive-Schwarz with interior boundaries
  INTEGER(i4) :: nlocal = 1 !< Number of subdomains on each processor
  INTEGER(i4) :: slice_group(10) = -1 !< Slice groupings
  INTEGER(i4), POINTER, DIMENSION(:) :: part => NULL()
  TYPE(oft_native_submatrix), pointer, dimension(:) :: alocals => NULL() !< Local matrix block
  TYPE(oft_solver_ptr), pointer, dimension(:) :: solvers => NULL() !< Local block solvers
  TYPE(oft_1d_int), pointer, dimension(:) :: parts => NULL() !< Local partition indices
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
!---------------------------------------------------------------------------
! TYPE oft_bc_ptr
!---------------------------------------------------------------------------
!> Boundary condition container
!---------------------------------------------------------------------------
type, public :: oft_bc_ptr
  !> Apply boundary condition to vector
  PROCEDURE(oft_bc_proto), POINTER, NOPASS :: bc => NULL()
end type oft_bc_ptr
!---Interfaces
ABSTRACT INTERFACE
!---------------------------------------------------------------------------
! INTERFACE oft_interp_proto
!---------------------------------------------------------------------------
!> Abstract interpolation prototype
!!
!! @param[in] a Temp doc
!! @param[in,out] b Temp doc
!---------------------------------------------------------------------------
  SUBROUTINE oft_interp_proto(a,b)
    IMPORT oft_vector
    CLASS(oft_vector), INTENT(inout) :: a
    CLASS(oft_vector), INTENT(inout) :: b
  END SUBROUTINE oft_interp_proto
!---------------------------------------------------------------------------
! INTERFACE oft_getop_proto
!---------------------------------------------------------------------------
!> Abstract operator retrieval prototype
!!
!! @param[in] a Temp doc
!---------------------------------------------------------------------------
  SUBROUTINE oft_getop_proto(a)
    IMPORT oft_matrix
    CLASS(oft_matrix), INTENT(in) :: a
  END SUBROUTINE oft_getop_proto
!---------------------------------------------------------------------------
! INTERFACE oft_veccreate_proto
!---------------------------------------------------------------------------
!> Abstract field creation prototype
!!
!! @param[out] new Temp doc
!! @param[in] level Temp doc
!---------------------------------------------------------------------------
  SUBROUTINE oft_veccreate_proto(new,level,cache,native)
    IMPORT oft_vector, i4
    CLASS(oft_vector), POINTER, INTENT(out) :: new
    INTEGER(i4), OPTIONAL, INTENT(in) :: level
    LOGICAL, OPTIONAL, INTENT(in) :: cache
    LOGICAL, OPTIONAL, INTENT(in) :: native
  END SUBROUTINE oft_veccreate_proto
!---------------------------------------------------------------------------
! INTERFACE oft_update_jacobian
!---------------------------------------------------------------------------
!> Abstract operator retrieval prototype
!!
!! @param[in] a Temp doc
!---------------------------------------------------------------------------
  SUBROUTINE oft_update_jacobian(a)
    IMPORT oft_vector
    CLASS(oft_vector), TARGET, INTENT(inout) :: a
  END SUBROUTINE oft_update_jacobian
END INTERFACE
!---Make classes and prototypes public
PUBLIC native_cg_solver_cast, native_gmres_solver_cast
PUBLIC jblock_precond_cast, ml_precond_cast
PUBLIC ml_trans_cast, diag_scale_cast, oft_veccreate_proto, oft_interp_proto
CONTAINS
!---------------------------------------------------------------------------
! FUNCTION: native_cg_solver_cast
!---------------------------------------------------------------------------
!> Cast oft_solver to oft_native_cg_solver
!!
!! The source solver must be cg_solver or a child class, otherwise an error will be thrown.
!!
!! @param[out] self Pointer to cast cg_solver
!! @param[in] source Source solver to cast
!---------------------------------------------------------------------------
FUNCTION native_cg_solver_cast(self,source) result(ierr)
type(oft_native_cg_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cg_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
END FUNCTION native_cg_solver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: cg_solver_apply
!---------------------------------------------------------------------------
!> Solve a linear system using the Conjugate-Gradient method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
recursive subroutine cg_solver_apply(self,u,g)
class(oft_native_cg_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
logical :: pm_save
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
  call r%add(0.d0,1.d0,g)
end if
if(associated(self%pre))call u%new(s)
if(associated(self%bc))call self%bc(g)
if(associated(self%orthog))call self%orthog%apply(u)
gg = g%dot(u)
if(gg==0.d0)then
  alfa=1.d0; f=0.d0
  call q%set(0.d0)
else
  call A%apply(u,q)
  if(associated(self%bc))call self%bc(q)
  pdq = u%dot(q)
  alfa=1.d0; f=pdq/2.d0-gg
endif
call g%add(1.d0,-1.d0,q)
call p%add(0.d0,1.d0,g)
if(associated(self%pre))then
  pm_save=oft_env%pm; oft_env%pm=self%pre%pm
  call self%pre%apply(s,p)
  oft_env%pm=pm_save
  if(associated(self%bc))call self%bc(s)
  call p%add(0.d0,1.d0,s)
end if
gg = g%dot(p)
ggin = g%dot(g)
uu=u%dot(u)
100 FORMAT (I6,3ES14.6)
110 FORMAT (I6,4ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,f,SQRT(uu),SQRT(ggin)
if(its==0.or.gg==0.d0)then
  call p%delete
  call q%delete
  deallocate(p,q)
  if(associated(self%pre))then
    call s%delete
    deallocate(s)
  end if
  DEBUG_STACK_POP
  return
endif
!---Begin CG iteration
max_its=INT(MIN(1000000_i8,u%ng+abs(its)),4)
do k=1,max_its
  if(associated(self%orthog))call self%orthog%apply(p)
  call A%apply(p,q)
  pdq = p%dot(q)
  if(associated(self%bc))call self%bc(q)
  alfa=gg/pdq
  fold=f; f=f-alfa*gg/2
  call u%add(1.d0,alfa,p)
  call g%add(1.d0,-alfa,q)
  if(associated(self%cleaner))call self%cleaner%apply(g)
  call q%add(0.d0,1.d0,g)
  if(associated(self%pre))then
    pm_save=oft_env%pm; oft_env%pm=self%pre%pm
    call self%pre%apply(s,q)
    oft_env%pm=pm_save
    if(associated(self%bc))call self%bc(s)
    call q%add(0.d0,1.d0,s)
  end if
  if(associated(self%orthog))call self%orthog%apply(q)
  ggold=gg
  gg=g%dot(q)
  ggnew = g%dot(g)
  uu=u%dot(u)
  beta=gg/ggold
  if((k<=self%itplot.OR.MOD(k,self%itplot)==0).AND.(oft_env%pm.AND.oft_env%head_proc))THEN
    write(*,110)k,f,SQRT(uu),SQRT(ggnew),SQRT(ggnew/uu)
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
  call g%add(0.d0,1.d0,r)
  DEALLOCATE(r)
end if
if(associated(self%pre))then
  call s%delete
  DEALLOCATE(s)
end if
if(k>=max_its)write(*,*)'Full Its: cgsolver_apply'
self%cits=k
DEBUG_STACK_POP
end subroutine cg_solver_apply
!---------------------------------------------------------------------------
! SUBROUTINE: cg_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine cg_setup_xml(self,solver_node,level)
CLASS(oft_native_cg_solver), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),atol(val_level),rtol(val_level))
!---
current_nodes=>fox_getElementsByTagName(solver_node,"its")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","cg_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"atol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,atol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough atol values specified","cg_setup_xml",__FILE__)
    self%atol=atol(val_level)
  ELSE
    self%atol=atol(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"rtol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,rtol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough rtol values specified","cg_setup_xml",__FILE__)
    self%rtol=rtol(val_level)
  ELSE
    self%rtol=rtol(1)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,*)'CG solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Tolerance:   ',REAL(self%atol,4),REAL(self%rtol,4)
END IF
!---
DEALLOCATE(its,atol,rtol)
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','cg_setup_xml',__FILE__)
#endif
end subroutine cg_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: cg_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine cg_delete(self)
class(oft_native_cg_solver), intent(inout) :: self
NULLIFY(self%A)
self%initialized=.FALSE.
end subroutine cg_delete
!---------------------------------------------------------------------------
! FUNCTION: native_gmres_solver_cast
!---------------------------------------------------------------------------
!> Cast oft_solver to oft_native_gmres_solver
!!
!! The source solver must be gmres_solver or a child class, otherwise an error will be thrown.
!!
!! @param[out] self Pointer to cast gmres_solver
!! @param[in] source Source solver to cast
!---------------------------------------------------------------------------
FUNCTION native_gmres_solver_cast(self,source) result(ierr)
type(oft_native_gmres_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_gmres_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION native_gmres_solver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: gmres_solver_apply
!---------------------------------------------------------------------------
!> Temp doc
!!
!! @param[in,out] self Temp doc
!! @param[in,out] u Temp doc
!! @param[in,out] g Temp doc
!---------------------------------------------------------------------------
recursive subroutine gmres_solver_apply(self,u,g)
class(oft_native_gmres_solver), intent(inout) :: self
class(oft_vector), intent(inout) :: u,g
class(oft_matrix), pointer :: A
logical :: pm_save
integer(i4) :: i,j,k,nits,its,nrits
real(r8) :: uu,uuold,gg,ggold,delta,hkmi,ggin
real(r8), allocatable :: h(:,:),c(:),s(:),res(:)
class(oft_vector), pointer :: r,w
type(oft_vector_ptr), pointer, dimension(:) :: v,z
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
if(associated(self%bc))call self%bc(g)
if(associated(self%cleaner))call self%cleaner%apply(g)
if(associated(self%orthog))call self%orthog%apply(u)
!
v=>self%v
z=>self%z
r=>self%r
w=>self%w
CALL w%set(0.d0)
allocate(c(nrits+1),s(nrits+1),res(nrits+1),h(nrits+1,nrits))
c=0.d0; s=0.d0; res=0.d0; h=0.d0
!
uu = u%dot(u)
IF(uu>0.d0)call A%apply(u,w)
if(associated(self%bc))call self%bc(w)
!if(associated(self%cleaner))call self%cleaner%apply(g)
call r%add(0.d0,1.d0,g,-1.d0,w)
!if(associated(self%cleaner))call self%cleaner%apply(r)
gg = r%dot(r)
ggin=gg
100 FORMAT (I6,2ES14.6)
110 FORMAT (I6,3ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,SQRT(uu),SQRT(gg)
if(self%its/=0.AND.self%nrits/=0)THEN
nits=0
do j=1,u%ng
  IF(gg==0.d0)EXIT
  res=0.d0
  res(1)=sqrt(gg)
  call v(1)%f%add(0.d0,1.d0/res(1),r)
  do i=1,nrits
    ! Precondition search direction
    if(associated(self%pre))then
      call w%add(0.d0,1.d0,v(i)%f)
      CALL z(i)%f%set(0.d0)
      pm_save=oft_env%pm; oft_env%pm=self%pre%pm
      self%pre%full_residual=.FALSE.
      call self%pre%apply(z(i)%f,w)
      self%pre%full_residual=.TRUE.
      oft_env%pm=pm_save
    else
      call z(i)%f%add(0.d0,1.d0,v(i)%f)
    end if
    if(associated(self%bc))call self%bc(z(i)%f)
    if(associated(self%orthog))call self%orthog%apply(z(i)%f)
    call A%apply(z(i)%f,w)
    if(associated(self%bc))call self%bc(w)
    !---Arnoldi iteration
    h(1:i,i)=w%mdot(v(1:i),i)
    do k=1,i
      call w%add(1.d0,-h(k,i),v(k)%f)
    end do
    !if(associated(self%cleaner))call self%cleaner%apply(w)
    h(i+1,i)=SQRT(w%dot(w))
    call v(i+1)%f%add(0.d0,1.d0/h(i+1,i),w)
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
    IF(nits==its)EXIT
    IF(ABS(res(i+1))<self%atol)EXIT
    IF(ABS(res(i+1))/SQRT(ggin)<self%rtol)EXIT
    IF(oft_env%head_proc.AND.oft_env%pm.AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
      WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
    END IF
  end do
  k=i
  do i=k,2,-1
    res(i) = res(i) / h(i,i)
    res(1:i-1) = res(1:i-1) - res(i)*h(1:i-1,i)
  end do
  res(1) = res(1) / h(1,1)
  ! Update iterate
  do i=1,k
    call u%add(1.d0,res(i),z(i)%f)
  enddo
  IF(((nits==its).OR.(k<nrits)).AND.(.NOT.self%full_residual))EXIT
  call A%apply(u,w)
  IF(associated(self%bc))call self%bc(w)
  call r%add(0.d0,1.d0,g,-1.d0,w)
  !if(associated(self%cleaner))call self%cleaner%apply(r)
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = r%dot(r)
  uu = u%dot(u)
  IF((oft_env%pm.AND.oft_env%head_proc).AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
    WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  END IF
  IF(its==-1.AND.sqrt(gg/uu)<1.d-7)EXIT
  IF(its==-2.AND.sqrt(gg/uu)<1.d-14)EXIT
  IF(SQRT(gg)<self%atol)EXIT
  IF(SQRT(gg/ggin)<self%rtol)EXIT
  IF(nits==its)EXIT
  IF(k<nrits)EXIT
  IF(uu==uuold)EXIT
end do
END IF
call g%add(0.d0,1.d0,r)
deallocate(c,s,h,res)
self%cits=nits
IF(nits==its)self%cits=-1
DEBUG_STACK_POP
end subroutine gmres_solver_apply
!---------------------------------------------------------------------------
! SUBROUTINE: gmres_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine gmres_setup_xml(self,solver_node,level)
CLASS(oft_native_gmres_solver), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its,nrits
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),nrits(val_level),atol(val_level),rtol(val_level))
!---
current_nodes=>fox_getElementsByTagName(solver_node,"its")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","gmres_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"nrits")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,nrits,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough nrits values specified","gmres_setup_xml",__FILE__)
    self%nrits=nrits(val_level)
  ELSE
    self%nrits=nrits(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"atol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,atol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough atol values specified","gmres_setup_xml",__FILE__)
    self%atol=atol(val_level)
  ELSE
    self%atol=atol(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"rtol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,rtol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough rtol values specified","gmres_setup_xml",__FILE__)
    self%rtol=rtol(val_level)
  ELSE
    self%rtol=rtol(1)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'GMRES solver setup:'
  WRITE(*,'(2X,A,I4)')    '- Iterations:  ',self%its
  WRITE(*,'(2X,A,I4)')    '- Restart:     ',self%nrits
  WRITE(*,'(2X,A,ES10.2)')'- Abs-Tol:     ',self%atol
  WRITE(*,'(2X,A,ES10.2)')'- Rel-Tol:     ',self%rtol
END IF
!---
DEALLOCATE(its,nrits,atol,rtol)
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','gmres_setup_xml',__FILE__)
#endif
end subroutine gmres_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: gmres_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine gmres_delete(self)
class(oft_native_gmres_solver), intent(inout) :: self
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
end subroutine gmres_delete
!---------------------------------------------------------------------------
! SUBROUTINE: cgmres_solver_apply
!---------------------------------------------------------------------------
!> Temp doc
!!
!! @param[in,out] self Temp doc
!! @param[in,out] u Temp doc
!! @param[in,out] g Temp doc
!---------------------------------------------------------------------------
recursive subroutine cgmres_solver_apply(self,u,g)
class(oft_native_gmres_csolver), intent(inout) :: self
class(oft_cvector), intent(inout) :: u,g
class(oft_cmatrix), pointer :: A
logical :: pm_save
integer(i4) :: i,j,k,nits,its,nrits
real(r8) :: uu,uuold,gg,ggold,ggin
complex(c8) :: delta,hkmi
complex(c8), allocatable :: h(:,:),c(:),s(:),res(:)
class(oft_cvector), pointer :: r,w
type(oft_cvector_ptr), pointer, dimension(:) :: v,z
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
if(associated(self%bc))call self%bc(g)
if(associated(self%cleaner))call self%cleaner%apply(g)
if(associated(self%orthog))call self%orthog%apply(u)
!
v=>self%v
z=>self%z
r=>self%r
w=>self%w
CALL w%set((0.d0,0.d0))
allocate(c(nrits+1),s(nrits+1),res(nrits+1),h(nrits+1,nrits))
c=0.d0; s=0.d0; res=0.d0; h=0.d0
!
uu = REAL(u%dot(u))
IF(uu>0.d0)call A%apply(u,w)
if(associated(self%bc))call self%bc(w)
!if(associated(self%cleaner))call self%cleaner%apply(g)
call r%add((0.d0,0.d0),(1.d0,0.d0),g,(-1.d0,0.d0),w)
!if(associated(self%cleaner))call self%cleaner%apply(r)
gg = REAL(r%dot(r))
ggin=gg
100 FORMAT (I6,2ES14.6)
110 FORMAT (I6,3ES14.6)
if((oft_env%pm.AND.oft_env%head_proc))write(*,100)0,SQRT(uu),SQRT(gg)
if(self%its/=0.AND.self%nrits/=0)THEN
nits=0
do j=1,u%ng
  IF(gg==0.d0)EXIT
  res=0.d0
  res(1)=sqrt(gg)
  call v(1)%f%add((0.d0,0.d0),(1.d0,0.d0)/res(1),r)
  do i=1,nrits
    ! Precondition search direction
    if(associated(self%pre))then
      call w%add((0.d0,0.d0),(1.d0,0.d0),v(i)%f)
      CALL z(i)%f%set((0.d0,0.d0))
      pm_save=oft_env%pm; oft_env%pm=self%pre%pm
      self%pre%full_residual=.FALSE.
      call self%pre%apply(z(i)%f,w)
      self%pre%full_residual=.TRUE.
      oft_env%pm=pm_save
    else
      call z(i)%f%add((0.d0,0.d0),(1.d0,0.d0),v(i)%f)
    end if
    if(associated(self%bc))call self%bc(z(i)%f)
    if(associated(self%orthog))call self%orthog%apply(z(i)%f)
    call A%apply(z(i)%f,w)
    if(associated(self%bc))call self%bc(w)
    !---Arnoldi iteration
    ! h(1:i,i)=w%mdot(v,i)
    do k=1,i
      h(k,i)=w%dot(v(k)%f)
      call w%add((1.d0,0.d0),-h(k,i),v(k)%f)
    end do
    !if(associated(self%cleaner))call self%cleaner%apply(w)
    h(i+1,i)=SQRT(w%dot(w))
    call v(i+1)%f%add((0.d0,0.d0),(1.d0,0.d0)/h(i+1,i),w)
    !if(associated(self%cleaner))call self%cleaner%apply(v(i+1)%f)
    !---Apply Givens rotation
    do k=2,i
      hkmi = h(k-1,i)
      h(k-1,i) = c(k-1)*hkmi + s(k-1)*h(k,i)
      h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
    enddo
    delta = SQRT( ABS(h(i,i))**2 + ABS(h(i+1,i))**2 )
    c(i) = h(i,i) / delta
    s(i) = h(i+1,i) / delta
    h(i,i) = c(i)*h(i,i) + s(i)*h(i+1,i)
    !---Update residual
    res(i+1) = -s(i)*res(i)
    res(i) = c(i)*res(i)
    nits=nits+1
    IF(i==nrits)EXIT
    IF(nits==its)EXIT
    IF(ABS(res(i+1))<self%atol)EXIT
    IF(ABS(res(i+1))/SQRT(ggin)<self%rtol)EXIT
    IF(oft_env%head_proc.AND.oft_env%pm.AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
      WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
    END IF
  end do
  k=i
  do i=k,2,-1
    res(i) = res(i) / h(i,i)
    res(1:i-1) = res(1:i-1) - res(i)*h(1:i-1,i)
  end do
  res(1) = res(1) / h(1,1)
  ! Update iterate
  do i=1,k
    call u%add((1.d0,0.d0),res(i),z(i)%f)
  enddo
  IF((nits==its).AND.(.NOT.self%full_residual))EXIT
  call A%apply(u,w)
  IF(associated(self%bc))call self%bc(w)
  call r%add((0.d0,0.d0),(1.d0,0.d0),g,(-1.d0,0.d0),w)
  !if(associated(self%cleaner))call self%cleaner%apply(r)
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = REAL(r%dot(r))
  uu = REAL(u%dot(u))
  IF((oft_env%pm.AND.oft_env%head_proc).AND.(nits<=self%itplot.OR.MOD(nits,self%itplot)==0))THEN
    WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  END IF
  IF(its==-1.AND.sqrt(gg/uu)<1.d-7)EXIT
  IF(its==-2.AND.sqrt(gg/uu)<1.d-14)EXIT
  IF(SQRT(gg)<self%atol)EXIT
  IF(SQRT(gg/ggin)<self%rtol)EXIT
  IF(nits==its)EXIT
  IF(uu==uuold)EXIT
end do
END IF
call g%add((0.d0,0.d0),(1.d0,0.d0),r)
deallocate(c,s,h,res)
self%cits=nits
IF(nits==its)self%cits=-1
DEBUG_STACK_POP
end subroutine cgmres_solver_apply
!---------------------------------------------------------------------------
! SUBROUTINE: cgmres_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine cgmres_setup_xml(self,solver_node,level)
CLASS(oft_native_gmres_csolver), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its,nrits
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),nrits(val_level),atol(val_level),rtol(val_level))
!---
current_nodes=>fox_getElementsByTagName(solver_node,"its")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","gmres_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"nrits")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,nrits,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough nrits values specified","gmres_setup_xml",__FILE__)
    self%nrits=nrits(val_level)
  ELSE
    self%nrits=nrits(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"atol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,atol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough atol values specified","gmres_setup_xml",__FILE__)
    self%atol=atol(val_level)
  ELSE
    self%atol=atol(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"rtol")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,rtol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough rtol values specified","gmres_setup_xml",__FILE__)
    self%rtol=rtol(val_level)
  ELSE
    self%rtol=rtol(1)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'GMRES solver setup:'
  WRITE(*,'(2X,A,I4)')    '- Iterations:  ',self%its
  WRITE(*,'(2X,A,I4)')    '- Restart:     ',self%nrits
  WRITE(*,'(2X,A,ES10.2)')'- Abs-Tol:     ',self%atol
  WRITE(*,'(2X,A,ES10.2)')'- Rel-Tol:     ',self%rtol
END IF
!---
DEALLOCATE(its,nrits,atol,rtol)
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','cgmres_setup_xml',__FILE__)
#endif
end subroutine cgmres_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: cgmres_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine cgmres_delete(self)
class(oft_native_gmres_csolver), intent(inout) :: self
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
end subroutine cgmres_delete
!---------------------------------------------------------------------------
! FUNCTION: native_cg_eigsolver_cast
!---------------------------------------------------------------------------
!> Cast oft_solver to oft_native_cg_solver
!!
!! The source eigsolver must be cg_eigsolver or a child class, otherwise an error will be thrown.
!!
!! @param[out] self Pointer to cast cg_eigsolver
!! @param[in] source Source eigsolver to cast
!---------------------------------------------------------------------------
FUNCTION native_cg_eigsolver_cast(self,source) result(ierr)
type(oft_native_cg_eigsolver), pointer, intent(out) :: self
class(oft_eigsolver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_native_cg_eigsolver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION native_cg_eigsolver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: cg_eigsolver_apply
!---------------------------------------------------------------------------
!> Solve a general eigenvalue system using the Conjugate-Gradient method.
!!
!! This solver uses a Non-Linear Conjugate-Gradient method to minimize the
!! Rayleigh Quotient (\f$ R = \frac{x*A*x}{x*M*x} \f$).
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] alam Eigenvalue
!---------------------------------------------------------------------------
subroutine cg_eigsolver_apply(self,u,alam)
class(oft_native_cg_eigsolver), intent(inout) :: self
class(oft_vector), intent(inout) :: u
real(r8), intent(inout) :: alam
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
if(associated(self%bc))call self%bc(u)
if(associated(self%orthog))call self%orthog%apply(u)
call A%apply(u,c)
call M%apply(u,b)
if(associated(self%bc))call self%bc(b)
if(associated(self%bc))call self%bc(c)
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
call f%add(0.d0,1.d0,g)
if(associated(self%pre))then
  pm_save=oft_env%pm; oft_env%pm=self%pre%pm
  call self%pre%apply(s,f)
  oft_env%pm=pm_save
  if(associated(self%bc))call self%bc(s)
  call f%add(0.d0,1.d0,s)
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
ninner = u%ng
if((self%ninner>0).AND.(self%ninner<u%ng))ninner=self%ninner
! iterate to find the solution
do ir=1,self%nrestarts ! restart loop
call e%add(0.d0,1.d0,f)
do n=1,ninner
  if(associated(self%orthog))call self%orthog%apply(e)
  call A%apply(e,g)
  call M%apply(e,f)
  if(associated(self%bc))then
    call self%bc(f)
    call self%bc(g)
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
  call f%add(0.d0,1.d0,g)
  if ( associated(self%pre) )then
    pm_save=oft_env%pm; oft_env%pm=self%pre%pm
    call self%pre%apply(s,f)
    oft_env%pm=pm_save
    if(associated(self%bc))call self%bc(s)
    call f%add(0.d0,1.d0,s)
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
!---------------------------------------------------------------------------
! SUBROUTINE: cg_eig_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine cg_eig_delete(self)
class(oft_native_cg_eigsolver), intent(inout) :: self
NULLIFY(self%A,self%M)
self%initialized=.FALSE.
end subroutine cg_eig_delete
!---------------------------------------------------------------------------
! SUBROUTINE: nksolver_apply
!---------------------------------------------------------------------------
!> Apply Newton's method to compute \f$ F(u) = g\f$
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine nksolver_apply(self,u,g)
class(oft_nksolver), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
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
!IF(ASSOCIATED(self%bc))call self%bc(g)
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
    IF(ASSOCIATED(self%bc))call self%bc(v)
    res=v%dot(v)
    uu=u%dot(u)
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
!---------------------------------------------------------------------------
! SUBROUTINE: nksolver_delete
!---------------------------------------------------------------------------
!> Destroy Newton solver and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine nksolver_delete(self)
class(oft_nksolver), intent(inout) :: self
IF(ASSOCIATED(self%du))THEN
  CALL self%du%delete
  DEALLOCATE(self%du)
END IF
IF(ASSOCIATED(self%A))NULLIFY(self%A)
IF(ASSOCIATED(self%J_inv))NULLIFY(self%J_inv)
IF(ASSOCIATED(self%J_update))NULLIFY(self%J_update)
IF(ASSOCIATED(self%bc))NULLIFY(self%bc)
end subroutine nksolver_delete
!---------------------------------------------------------------------------
! SUBROUTINE: identity_inv_apply
!---------------------------------------------------------------------------
!> Solver container for trivial inverse of the identity matrix
!!
!! @note Used for simple cases when a solver wrapper is required but the trivial
!! case of \f$ u = g \f$ is the desired result
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine identity_inv_apply(self,u,g)
class(oft_identity_inv), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
DEBUG_STACK_PUSH
if(g%n/=u%n)call oft_abort('Size mismatch','identity_inv_apply',__FILE__)
call u%add(0.d0,1.d0,g)
DEBUG_STACK_POP
end subroutine identity_inv_apply
!---------------------------------------------------------------------------
! SUBROUTINE: identity_inv_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine identity_inv_delete(self)
class(oft_identity_inv), intent(inout) :: self
end subroutine identity_inv_delete
!---------------------------------------------------------------------------
! FUNCTION: diag_scale_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_native_solvers::oft_solver "oft_solver" to @ref oft_native_solvers::oft_diag_scale
!! "oft_diag_scale"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION diag_scale_cast(self,source) result(ierr)
type(oft_diag_scale), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  type is(oft_diag_scale)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION diag_scale_cast
!---------------------------------------------------------------------------
! SUBROUTINE: diag_scale_apply
!---------------------------------------------------------------------------
!> Apply diagonal scaling
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine diag_scale_apply(self,u,g)
class(oft_diag_scale), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
DEBUG_STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch.','diag_scale_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch.','diag_scale_apply',__FILE__)
if(.NOT.associated(self%A%D))call oft_abort('No diagonal entries.','diag_scale_apply',__FILE__)
!---
call u%add(0.d0,1.d0,g)
call u%mult(self%A%D,div_flag=.TRUE.)
DEBUG_STACK_POP
end subroutine diag_scale_apply
!---------------------------------------------------------------------------
! SUBROUTINE: diag_scale_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine diag_scale_delete(self)
class(oft_diag_scale), intent(inout) :: self
NULLIFY(self%A)
end subroutine diag_scale_delete
!---------------------------------------------------------------------------
! SUBROUTINE: cdiag_scale_apply
!---------------------------------------------------------------------------
!> Apply diagonal scaling
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine cdiag_scale_apply(self,u,g)
class(oft_diag_cscale), intent(inout) :: self
class(oft_cvector), intent(inout) :: u
class(oft_cvector), intent(inout) :: g
DEBUG_STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch.','cdiag_scale_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch.','cdiag_scale_apply',__FILE__)
if(.NOT.associated(self%A%D))call oft_abort('No diagonal entries.','cdiag_scale_apply',__FILE__)
!---
call u%add((0.d0,0.d0),(1.d0,0.0),g)
call u%mult(self%A%D,div_flag=.TRUE.)
DEBUG_STACK_POP
end subroutine cdiag_scale_apply
!---------------------------------------------------------------------------
! SUBROUTINE: cdiag_scale_delete
!---------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine cdiag_scale_delete(self)
class(oft_diag_cscale), intent(inout) :: self
NULLIFY(self%A)
end subroutine cdiag_scale_delete
!---------------------------------------------------------------------------
! FUNCTION: jblock_precond_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_native_solvers::oft_solver "oft_solver" to @ref oft_native_solvers::oft_jblock_precond
!! "oft_jblock_precond"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION jblock_precond_cast(self,source) result(ierr)
type(oft_jblock_precond), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  type is(oft_jblock_precond)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION jblock_precond_cast
!---------------------------------------------------------------------------
! SUBROUTINE: jblock_precond_apply
!---------------------------------------------------------------------------
!> Apply 1-step of a symmetric Jacobi smoother with native CRS matrices
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine jblock_precond_apply(self,u,g)
class(oft_jblock_precond), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
class(oft_vector), pointer :: p,q
class(oft_native_vector), pointer :: uv,gv,md
integer(i4) :: i,k,ierr
DEBUG_STACK_PUSH
IF(g%n/=self%A%nr)CALL oft_abort('Row mismatch','jblock_precond_apply',__FILE__)
IF(u%n/=self%A%nc)CALL oft_abort('Col mismatch','jblock_precond_apply',__FILE__)
IF(native_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a vector object.','jblock_precond_apply',__FILE__)
IF(native_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a vector object.','jblock_precond_apply',__FILE__)
IF(native_vector_cast(md,self%A%D)<0)CALL oft_abort('"A%D" is not a vector object.','jblock_precond_apply',__FILE__)
!---Initialize solver
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%u_save(u%n),self%g_save(u%n))
  CALL u%new(p)
  ierr=native_vector_cast(self%p,p)
  CALL u%new(q)
  ierr=native_vector_cast(self%q,q)
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
  !$omp parallel do if(u%n>OFT_OMP_VTHRESH/2)
  do i=1,u%n
    uv%v(i)=0.d0
    self%p%v(i)=self%df*gv%v(i)/md%v(i)
  end do
  IF(ASSOCIATED(self%bc))call self%bc(self%p)
  do k=1,self%its
    call self%A%apply(self%p,self%q)
    IF(ASSOCIATED(self%bc))call self%bc(self%q)
    !$omp parallel do if(u%n>OFT_OMP_VTHRESH/2)
    do i=1,u%n
      uv%v(i)=uv%v(i)+self%p%v(i)
      gv%v(i)=gv%v(i)-self%q%v(i)
      self%p%v(i)=self%df*gv%v(i)/md%v(i)
    end do
    IF(ASSOCIATED(self%bc))call self%bc(self%p)
  end do
  !$omp parallel do if(u%n>OFT_OMP_VTHRESH)
  do i=1,u%n
    self%u_save(i)=uv%v(i)
    self%g_save(i)=gv%v(i)
  end do
  self%down=.FALSE.
else ! Up smoother
  !$omp parallel do if(u%n>OFT_OMP_VTHRESH)
  do i=1,u%n
    self%p%v(i)=uv%v(i)-self%u_save(i)
    gv%v(i)=self%g_save(i)
  end do
  IF(ASSOCIATED(self%bc))call self%bc(self%p)
  do k=1,self%its
    call self%A%apply(self%p,self%q)
    IF(ASSOCIATED(self%bc))call self%bc(self%q)
    !$omp parallel do if(u%n>OFT_OMP_VTHRESH/2)
    do i=1,u%n
      gv%v(i)=gv%v(i)-self%q%v(i)
      self%p%v(i)=self%df*gv%v(i)/md%v(i)
    end do
    IF(ASSOCIATED(self%bc))call self%bc(self%p)
    !$omp parallel do if(u%n>OFT_OMP_VTHRESH)
    do i=1,u%n
      uv%v(i)=uv%v(i)+self%p%v(i)
    end do
  end do
  self%down=.TRUE.
end if
DEBUG_STACK_POP
end subroutine jblock_precond_apply
!---------------------------------------------------------------------------
! SUBROUTINE: jblock_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine jblock_setup_xml(self,solver_node,level)
CLASS(oft_jblock_precond), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node
TYPE(fox_nodelist), POINTER :: current_nodes
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),df(val_level))
!---
current_nodes=>fox_getElementsByTagName(solver_node,"its")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","jblock_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
current_nodes=>fox_getElementsByTagName(solver_node,"df")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,df,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough df values specified","jblock_setup_xml",__FILE__)
    self%df=df(val_level)
  ELSE
    self%df=df(1)
  END IF
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,*)'S-Jacobi solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Factor:      ',self%df
END IF
!---
DEALLOCATE(its,df)
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','jblock_setup_xml',__FILE__)
#endif
end subroutine jblock_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: jblock_precond_delete
!---------------------------------------------------------------------------
!> Destroy symmetric Jacobi preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine jblock_precond_delete(self)
class(oft_jblock_precond), intent(inout) :: self
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
DEBUG_STACK_POP
end subroutine jblock_precond_delete
!---------------------------------------------------------------------------
! FUNCTION: ml_precond_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_native_solvers::oft_solver "oft_solver" to @ref oft_native_solvers::oft_ml_precond
!! "oft_ml_precond"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION ml_precond_cast(self,source) result(ierr)
class(oft_ml_precond), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_ml_precond)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION ml_precond_cast
!---------------------------------------------------------------------------
! SUBROUTINE: ml_precond_view
!---------------------------------------------------------------------------
!> Print solver configuration
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!---------------------------------------------------------------------------
recursive subroutine ml_precond_view(self)
class(oft_ml_precond), intent(inout) :: self
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
!---------------------------------------------------------------------------
! SUBROUTINE: ml_precond_apply
!---------------------------------------------------------------------------
!> Apply 1-step of a Multi-Level preconditioner
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
recursive subroutine ml_precond_apply(self,u,g)
class(oft_ml_precond), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!---Get RHS
IF(.NOT.self%initialized)THEN
  CALL u%new(self%r)
  IF(ASSOCIATED(self%base_solve))THEN
    CALL u%new(self%p)
    CALL self%vec_create(self%ucors,self%level-1)
    CALL self%vec_create(self%gcors,self%level-1)
  END IF
  CALL solver_setup(self)
END IF
call self%r%add(0.d0,1.d0,g)
!---Apply down smoother
if(associated(self%smooth_down))then
  CALL mytimer%tick
  call self%smooth_down%apply(u,g)
  self%timings(1)=self%timings(1)+mytimer%tock()
end if
if(associated(self%base_solve))then
  CALL mytimer%tick
  !---Inject RHS
  call self%ucors%set(0.d0)
  call self%inject(g,self%gcors)
  self%timings(2)=self%timings(2)+mytimer%tock()
  CALL mytimer%tick
  !---Apply lower level smoother
  self%base_solve%full_residual=.FALSE.
  call self%base_solve%apply(self%ucors,self%gcors)
  self%base_solve%full_residual=.TRUE.
  self%timings(3)=self%timings(3)+mytimer%tock()
  CALL mytimer%tick
  !---Interpolate Solution
  call self%interp(self%ucors,self%p)
  if(associated(self%bc))call self%bc(self%p)
  !---Add coarse correction
  call u%add(1.d0,1.d0,self%p)
  self%timings(4)=self%timings(4)+mytimer%tock()
end if
!---Reset RHS
call g%add(0.d0,1.d0,self%r)
!---Apply up smoother
if(associated(self%smooth_up))then
  CALL mytimer%tick
  self%smooth_up%full_residual=.FALSE.
  call self%smooth_up%apply(u,g)
  self%smooth_up%full_residual=.TRUE.
  self%timings(1)=self%timings(1)+mytimer%tock()
end if
!---Complete
DEBUG_STACK_POP
end subroutine ml_precond_apply
!---------------------------------------------------------------------------
! SUBROUTINE: ml_precond_update
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
recursive subroutine ml_precond_update(self,new_pattern)
class(oft_ml_precond), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
!---
IF(ASSOCIATED(self%smooth_up))CALL self%smooth_up%update(new_pattern)
IF(ASSOCIATED(self%smooth_down))CALL self%smooth_down%update(new_pattern)
IF(ASSOCIATED(self%base_solve))CALL self%base_solve%update(new_pattern)
end subroutine ml_precond_update
!---------------------------------------------------------------------------
! SUBROUTINE: ml_precond_delete
!---------------------------------------------------------------------------
!> Destroy Multi-Level preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
recursive subroutine ml_precond_delete(self)
class(oft_ml_precond), intent(inout) :: self
DEBUG_STACK_PUSH
if(ASSOCIATED(self%smooth_up))then
  call self%smooth_up%delete
  DEALLOCATE(self%smooth_up)
end if
if(ASSOCIATED(self%smooth_down))then
  IF(self%symmetric)THEN
    NULLIFY(self%smooth_down)
  ELSE
    call self%smooth_down%delete
    DEALLOCATE(self%smooth_down)
  END IF
end if
if(associated(self%base_solve))then
  call self%base_solve%delete
  DEALLOCATE(self%base_solve)
end if
IF(ASSOCIATED(self%interp))NULLIFY(self%interp)
IF(ASSOCIATED(self%inject))NULLIFY(self%inject)
IF(ASSOCIATED(self%vec_create))NULLIFY(self%vec_create)
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
!---------------------------------------------------------------------------
! FUNCTION: ml_trans_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_native_solvers::oft_solver "oft_solver" to @ref oft_native_solvers::oft_ml_trans
!! "oft_ml_trans"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION ml_trans_cast(self,source) result(ierr)
class(oft_ml_trans), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_ml_trans)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION ml_trans_cast
!---------------------------------------------------------------------------
! SUBROUTINE: ml_trans_apply
!---------------------------------------------------------------------------
!> Transfer solution between distributed and shared levels as part of a ML
!! preconditioner
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
recursive subroutine ml_trans_apply(self,u,g)
class(oft_ml_trans), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
class(oft_vector), pointer :: ucors,gcors
DEBUG_STACK_PUSH
!---Create base level fields
call self%vec_create(ucors,self%level-1)
call self%vec_create(gcors,self%level-1)
!---Push RHS to base level
call self%inject(g,gcors)
!---Solve on base levels
call self%base_solve%apply(ucors,gcors)
!---Pop solution from base level
call self%interp(ucors,u)
!---Delete base level fields
call ucors%delete()
call gcors%delete()
!---Complete
DEBUG_STACK_POP
end subroutine ml_trans_apply
!---------------------------------------------------------------------------
! SUBROUTINE: ml_trans_delete
!---------------------------------------------------------------------------
!> Destroy Multi-Level preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine ml_trans_delete(self)
class(oft_ml_trans), intent(inout) :: self
DEBUG_STACK_PUSH
if(associated(self%smooth_up))then
  call self%smooth_up%delete
  nullify(self%smooth_up)
end if
if(associated(self%smooth_down))then
  call self%smooth_down%delete
  nullify(self%smooth_down)
end if
if(associated(self%base_solve))then
  call self%base_solve%delete
  nullify(self%base_solve)
end if
IF(ASSOCIATED(self%interp))NULLIFY(self%interp)
IF(ASSOCIATED(self%inject))NULLIFY(self%inject)
IF(ASSOCIATED(self%vec_create))NULLIFY(self%vec_create)
IF(ASSOCIATED(self%bc))NULLIFY(self%bc)
DEBUG_STACK_POP
end subroutine ml_trans_delete
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_apply
!---------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE bjprecond_apply(self,u,g)
CLASS(oft_bjprecond), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u,g
!---
LOGICAL :: thread_safe,color_avail
LOGICAL, POINTER, DIMENSION(:) :: eflag
INTEGER(i4) :: i,j,jr,k,ngroups,ncolors
INTEGER(i4), POINTER, DIMENSION(:) :: part,part_tmp,sflag
TYPE(oft_native_vector) :: uloc,gloc
REAL(r8), POINTER, DIMENSION(:) :: vals,uvals,gvals
TYPE(oft_graph) :: graph
CLASS(oft_native_matrix), POINTER :: A_native
CLASS(oft_solver), POINTER :: pretmp
TYPE(oft_1d_int), pointer, dimension(:) :: parts_tmp => NULL()
DEBUG_STACK_PUSH
!---
IF(native_matrix_cast(A_native,self%A)<0)CALL oft_abort('Native matrix required', &
  'bjprecond_apply',__FILE__)
!---Initialize local matrix
IF(.NOT.self%initialized)THEN
  color_avail=.FALSE.
  CALL solver_setup(self)
  IF(self%nlocal==-1)THEN
    self%nlocal=-self%A%ni
    self%slice_group(1:self%A%ni)=(/(i,i=1,self%A%ni)/)
  END IF
  IF(ASSOCIATED(A_native%color))THEN
    ncolors=MAXVAL(A_native%color)
    IF(self%nlocal<=1.AND.ncolors>1)THEN
      self%nlocal=self%nlocal*ncolors
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
      !---Do not include redundant rows
      IF(ASSOCIATED(A_native%redind))THEN
        DO i=1,A_native%nred
          part(A_native%redind(i))=-1
        END DO
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
    !---Create partition mappings
    ALLOCATE(self%parts(self%nlocal))
    IF(ANY(part==0))CALL oft_abort("Unclaimed partition","bjprecond_apply",__FILE__)
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
    IF(self%loverlap.AND.self%nlocal>1)THEN
      ALLOCATE(sflag(A_native%nr))
      DO i=1,self%A%ni
        DO k=1,self%A%i_map(i)%nslice
          j=self%A%i_map(i)%offset+self%A%i_map(i)%slice(k)
          sflag(j)=self%slice_group(i)
        END DO
      END DO
      parts_tmp=>self%parts
      ALLOCATE(self%parts(self%nlocal))
      !$omp parallel private(j,jr,k,eflag)
      ALLOCATE(eflag(A_native%nr))
      !$omp do
      DO i=1,self%nlocal
        eflag=(part==i)
        DO j=1,parts_tmp(i)%n
          jr=parts_tmp(i)%v(j)
          DO k=A_native%kr(jr),A_native%kr(jr+1)-1
            IF(sflag(A_native%lc(k))==sflag(jr))eflag(A_native%lc(k))=.TRUE.
          END DO
        END DO
        DO j=1,A_native%nr
          IF(part(j)<0)eflag(j)=.FALSE.
        END DO
        self%parts(i)%n=COUNT(eflag)
        ALLOCATE(self%parts(i)%v(self%parts(i)%n))
        k=0
        DO j=1,A_native%nr
          IF(.NOT.eflag(j))CYCLE
          k=k+1
          self%parts(i)%v(k)=j
        END DO
        DEALLOCATE(parts_tmp(i)%v)
      END DO
      DEALLOCATE(eflag)
      !$omp end parallel
      DEALLOCATE(parts_tmp,sflag)
    END IF
    !---Create dummy local vector
    uloc%n=A_native%nr
    ALLOCATE(uloc%stitch_info)
    uloc%stitch_info%skip=.TRUE.
    !---Create local sub-matrices
    DO i=1,self%nlocal
      CALL self%alocals(i)%setup(A_native,uloc,part=self%parts(i)%v)
      ALLOCATE(self%solvers(i)%s, source=pretmp)
    END DO
    IF(self%loverlap.AND.self%nlocal>1)THEN
      !$omp parallel do private(j)
      DO i=1,self%nlocal
        DO j=1,self%parts(i)%n
          IF(part(self%parts(i)%v(j))/=i)self%parts(i)%v(j)=-self%parts(i)%v(j)
        END DO
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
  IF(self%update_slice)THEN
    CALL A_native%update_slice
    DO i=1,self%nlocal
      CALL self%alocals(i)%update_slice
    END DO
  END IF
  NULLIFY(uvals,gvals)
  CALL u%get_local(uvals)
  CALL g%get_local(gvals)
  !$omp parallel private(j,uloc,gloc) if(thread_safe)
  NULLIFY(gloc%local_tmp,uloc%local_tmp)
  ALLOCATE(gloc%stitch_info,uloc%stitch_info)
  gloc%stitch_info%skip=.TRUE.
  uloc%stitch_info%skip=.TRUE.
  !$omp do schedule(dynamic,1)
  DO i=1,self%nlocal
    !---Setup local vectors
    gloc%n=self%parts(i)%n; uloc%n=self%parts(i)%n
    ALLOCATE(gloc%v(gloc%n),uloc%v(uloc%n))
    !---Get local slice for solve
    !$omp simd
    DO j=1,self%parts(i)%n
      gloc%v(j)=gvals(ABS(self%parts(i)%v(j)))
      uloc%v(j)=uvals(ABS(self%parts(i)%v(j)))
    END DO
    !---Solve sub-system
    self%solvers(i)%s%A=>self%alocals(i)
    CALL self%solvers(i)%s%apply(uloc,gloc)
    !---Replace local slice into solution
    !$omp simd
    DO j=1,self%parts(i)%n
      IF(self%parts(i)%v(j)>0)uvals(self%parts(i)%v(j))=uloc%v(j)
    END DO
    DEALLOCATE(gloc%v,uloc%v)
  END DO
  DEALLOCATE(gloc%stitch_info,uloc%stitch_info)
  !$omp end parallel
  CALL u%restore_local(uvals)
  DEALLOCATE(uvals,gvals)
ELSE
  IF(self%boundary_overlap)THEN
    IF(self%update_slice)CALL A_native%update_slice
    !---Setup local vectors
    gloc%n=g%n; uloc%n=u%n
    NULLIFY(gloc%local_tmp,uloc%local_tmp)
    ALLOCATE(gloc%v(gloc%n),uloc%v(uloc%n))
    ALLOCATE(gloc%stitch_info,uloc%stitch_info)
    gloc%stitch_info%skip=.TRUE.
    uloc%stitch_info%skip=.TRUE.
    !---Get local slice for solve
    CALL g%get_local(gloc%v)
    CALL u%get_local(uloc%v)
    !---Zero source on redundant rows
    IF(ASSOCIATED(A_native%redind))THEN
      DO i=1,A_native%nred
        gloc%v(A_native%redind(i))=0.d0
      END DO
    END IF
    !---Solve sub-system
    self%solvers(1)%s%A=>self%A
    CALL self%solvers(1)%s%apply(uloc,gloc)
    !---Replace local slice into solution
    CALL u%restore_local(uloc%v)
  ELSE
    IF(self%update_slice)CALL self%alocals(1)%update_slice
    !---Setup local vectors
    gloc%n=g%nslice; uloc%n=u%nslice
    NULLIFY(gloc%local_tmp,uloc%local_tmp)
    ALLOCATE(gloc%v(gloc%n),uloc%v(uloc%n))
    ALLOCATE(gloc%stitch_info,uloc%stitch_info)
    gloc%stitch_info%skip=.TRUE.
    uloc%stitch_info%skip=.TRUE.
    !---Get local slice for solve
    CALL g%get_slice(gloc%v)
    CALL u%get_slice(uloc%v)
    !---Solve sub-system
    self%solvers(1)%s%A=>self%alocals(1)
    CALL self%solvers(1)%s%apply(uloc,gloc)
    !---Replace local slice into solution
    CALL u%restore_slice(uloc%v)
  END IF
  DEALLOCATE(gloc%v,uloc%v,gloc%stitch_info,uloc%stitch_info)
END IF
self%update_slice=.FALSE.
DEBUG_STACK_POP
END SUBROUTINE bjprecond_apply
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_update
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! @param[in] new_pattern Update matrix pattern (optional)
!---------------------------------------------------------------------------
recursive subroutine bjprecond_update(self,new_pattern)
class(oft_bjprecond), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
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
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine bjprecond_setup_xml(self,solver_node,level)
CLASS(oft_bjprecond), INTENT(inout) :: self
TYPE(fox_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
INTEGER(i4) :: nnodes,nread
TYPE(fox_node), POINTER :: current_node,sub_node
TYPE(fox_nodelist), POINTER :: current_nodes
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: nlocals
DEBUG_STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(nlocals(val_level))
!---Read-in desired number of subdomains
current_nodes=>fox_getElementsByTagName(solver_node,"nlocal")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,nlocals,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough local sizes specified", &
    "bjprecond_setup_xml",__FILE__)
    self%nlocal=nlocals(val_level)
  ELSE
    self%nlocal=nlocals(1)
  END IF
END IF
IF(self%nlocal<-1)THEN
  !---Read-in desired number of subdomains
  current_nodes=>fox_getElementsByTagName(solver_node,"groups")
  nnodes=fox_getLength(current_nodes)
  IF(nnodes==1)THEN
    current_node=>fox_item(current_nodes,0)
    CALL fox_extractDataContent(current_node,self%slice_group,num=nread,iostat=ierr)
  END IF
END IF
!---Read-in desired boundary overlap specification
current_nodes=>fox_getElementsByTagName(solver_node,"boverlap")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,self%boundary_overlap,num=nread,iostat=ierr)
  IF(nread>1)CALL oft_abort("boverlap must be single value","bjprecond_setup_xml",__FILE__)
END IF
!---Read-in desired internal overlap specification
current_nodes=>fox_getElementsByTagName(solver_node,"loverlap")
nnodes=fox_getLength(current_nodes)
IF(nnodes==1)THEN
  current_node=>fox_item(current_nodes,0)
  CALL fox_extractDataContent(current_node,self%loverlap,num=nread,iostat=ierr)
  IF(nread>1)CALL oft_abort("loverlap must be single value","bjprecond_setup_xml",__FILE__)
END IF
IF(oft_debug_print(1))THEN
  WRITE(*,'(A)')'Block Jacobi solver setup:'
  WRITE(*,'(2X,A,I4)')  '- NLocal:    ',self%nlocal
  WRITE(*,'(2X,A,3X,L)')'- Boverlap:  ',self%boundary_overlap
  WRITE(*,'(2X,A,3X,L)')'- Loverlap:  ',self%loverlap
END IF
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','bjprecond_setup_xml',__FILE__)
#endif
end subroutine bjprecond_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_delete
!---------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine bjprecond_delete(self)
class(oft_bjprecond), intent(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%solvers))THEN
  DO i=1,self%nlocal
    CALL self%solvers(i)%s%delete
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
!---Destory local partitions
IF(ASSOCIATED(self%parts))THEN
  DO i=1,self%nlocal
    IF(ASSOCIATED(self%parts(i)%v))DEALLOCATE(self%parts(i)%v)
  END DO
  DEALLOCATE(self%parts)
END IF
!---Reset
self%nlocal=1
self%initialized=.FALSE.
DEBUG_STACK_POP
end subroutine bjprecond_delete
end module oft_native_solvers
