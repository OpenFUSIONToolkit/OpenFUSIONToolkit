!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_petsc_solvers.F90
!
!> PETSc solver implementations
!!
!! Wrappers for PETSc solvers implementing the common @ref oft_solvers::oft_solver
!! "oft_solver" interface.
!! - Conjugate-Gradient
!! - Generalized Minimum Residual
!! - LU Factorization (SUPERLU, MUMPS, etc.)
!!
!! Wrappers for PETSc preconditioners implementing the common @ref oft_solvers::oft_solver
!! "oft_solver" interface.
!! - Diagonal scaling
!! - Additive Schwarz
!! - Block Jacobi
!! - Multi-Grid
!! - LU Factorization (SUPERLU, MUMPS, etc.)
!!
!! @authors Chris Hansen
!! @date January 2013
!! @ingroup doxy_oft_petsc_la
!------------------------------------------------------------------------------
MODULE oft_petsc_solvers
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver, oft_bc_proto
#ifdef HAVE_PETSC
USE oft_petsc_la, ONLY: oft_petsc_vector, oft_petsc_vector_cast, &
  oft_petsc_matrix, oft_petsc_matrix_cast
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>5)
#if PETSC_VERSION_MINOR<8
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petsckspdef.h"
#include "petsc/finclude/petscpcdef.h"
#define PETSC_NULL_KSP PETSC_NULL_OBJECT
#define PETSC_NULL_PC PETSC_NULL_OBJECT
#define PETSC_NULL_FUNCTION PETSC_NULL_OBJECT
#else
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#endif
#else
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#endif
USE petscmat
USE petscksp
#if !(PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
#define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION
#endif
USE petscpc
#undef Vec
#undef Mat
#undef IS
#if !defined(OFT_LU)
#define OFT_LU 2
#endif
#if OFT_LU == 1
#define OFT_PETSC_LU MATSOLVERSUPERLU
#elif OFT_LU == 2
#define OFT_PETSC_LU MATSOLVERSUPERLU_DIST
#elif OFT_LU == 3
#define OFT_PETSC_LU MATSOLVERMUMPS
#endif
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! TYPE oft_petsc_factordef
!---------------------------------------------------------------------------
!> PETSc abstract preconditioner class
!---------------------------------------------------------------------------
type :: oft_petsc_factordef
  CHARACTER(LEN=6) :: type = 'lu' !< Factor type
  CHARACTER(LEN=6) :: package = 'superd' !< Factor package
end type oft_petsc_factordef
!---------------------------------------------------------------------------
! TYPE oft_petsc_precond
!---------------------------------------------------------------------------
!> PETSc abstract preconditioner class
!!
!! @note If called from a native solver a solver wrapper will be created
!---------------------------------------------------------------------------
type, abstract, public, extends(oft_solver) :: oft_petsc_precond
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
  INTEGER(petsc_addr) :: pre_obj !< PETSc preconditioner object
  INTEGER(petsc_addr) :: solver !< PETSc solver object
#else
  TYPE(tpc) :: pre_obj !< PETSc preconditioner object
  TYPE(tksp) :: solver !< PETSc solver object
#endif
contains
  !> Update preconditioner with new settings/operators
  PROCEDURE :: update => precond_update
  !> Solve system
  PROCEDURE :: apply => precond_apply
  !> Setup preconditioner
  PROCEDURE :: setup => precond_dummy_setup
  !> Destory wrapper solver
  PROCEDURE :: delete => precond_delete
end type oft_petsc_precond
!---------------------------------------------------------------------------
! TYPE oft_petsc_mgprecond
!---------------------------------------------------------------------------
!> PETSc MG preconditioner
!!
!! Use PETSc MG framework for preconditioning. \b V and \b W cycle types are
!! supported. \b Galerkin construction of sub-level matrices is also available.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_mgprecond
  TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: I => NULL() !< Interpolation matrices
  TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: Mats => NULL() !< Operator matrices
  INTEGER(i4) :: nlevels = 0 !< Number of levels
  INTEGER(i4), POINTER, DIMENSION(:) :: nu => NULL() !< Number of smoother applications
  REAL(r8), POINTER, DIMENSION(:) :: df => NULL() !< Damping factors
  CHARACTER(LEN=1) :: cycle_type = "V" !< MG cycle type
  CHARACTER(LEN=8) :: smooth_type = "jacobi" !< Smoother type
  LOGICAL :: direct_base = .TRUE. !< Use direct solver for base level
  LOGICAL :: galerkin = .FALSE. !< Use Galerkin method to construct operator matrices
contains
   !> Setup preconditioner
  procedure :: setup => mgprecond_setup
end type oft_petsc_mgprecond
!---------------------------------------------------------------------------
! TYPE oft_petsc_diagprecond
!---------------------------------------------------------------------------
!> PETSc diagonal scaling
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_diagprecond
contains
  !> Setup preconditioner
  procedure :: setup => diagprecond_setup
  !> Print solver information
  PROCEDURE :: view => diagprecond_view
end type oft_petsc_diagprecond
!---------------------------------------------------------------------------
! TYPE oft_petsc_luprecond
!---------------------------------------------------------------------------
!> PETSc LU preconditioner
!!
!! Use PETSc direct solver, LU decomposition, for preconditioning.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_luprecond
contains
  !> Setup preconditioner
  procedure :: setup => luprecond_setup
end type oft_petsc_luprecond
!---------------------------------------------------------------------------
! TYPE oft_petsc_asprecond
!---------------------------------------------------------------------------
!> PETSc Additive-Schwarz preconditioner
!!
!! Use PETSc AS framework for preconditioning.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_asprecond
  INTEGER(i4) :: overlap = 1 !< Size of domain overlap
  INTEGER(i4) :: n_local = 0 !< Use local field partitioning?
  TYPE(oft_petsc_factordef) :: sub_factor !< PETSc factor definition
contains
  !> Setup preconditioner
  procedure :: setup => asprecond_setup
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => asprecond_setup_xml
end type oft_petsc_asprecond
!---------------------------------------------------------------------------
! TYPE oft_petsc_solver
!---------------------------------------------------------------------------
!> PETSc abstract solver class
!---------------------------------------------------------------------------
type, abstract, public, extends(oft_solver) :: oft_petsc_solver
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
  INTEGER(petsc_addr) :: solver !< PETSc solver object
#else
  TYPE(tksp) :: solver !< PETSc solver object
#endif
  LOGICAL :: dist = .FALSE. !< Matrix/Vector are distributed
CONTAINS
  !> Update solver with new settings/operators
  PROCEDURE :: update => petsc_solver_update
  !> Setup PETSc solver object
  PROCEDURE :: setup_ksp => petsc_solver_setup_ksp
  !> Clean-up internal storage
  PROCEDURE :: delete => petsc_solver_delete
end type oft_petsc_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_pre_solver
!---------------------------------------------------------------------------
!> PETSc preconditioner as solver
!!
!! Apply PETSc preconditioner as a solver. For use with native solvers where
!! the preconditioner must be wrapped as a stand alone solver.
!!
!! @deprecated This implementation is obsolete, preconditioner objects should
!! now be refenced by native solvers directly.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_pre_solver
contains
  !> Solve system
  procedure :: apply => pre_solver_apply
end type oft_petsc_pre_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_cg_solver
!---------------------------------------------------------------------------
!> PETSc CG solver
!!
!! Apply PETSc Conjugate-Gradient algorithm.
!!
!! @note Matrix must be SPD, otherwise solver will fail.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_cg_solver
contains
  !> Solve system
  procedure :: apply => cg_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => cg_setup_xml
  !> Print solver information
  PROCEDURE :: view => cg_view
  !> Setup PETSc solver object
  PROCEDURE :: setup_ksp => cg_setup_ksp
end type oft_petsc_cg_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_gmres_solver
!---------------------------------------------------------------------------
!> PETSc GMRES solver
!!
!! Apply PETSc GMRES algorithm.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_gmres_solver
  integer(i4) :: nrits=10 !< Number of iterations before restart
  integer(i4) :: algorithm=1 !< Algorithm type
contains
  !> Solve system
  procedure :: apply => gmres_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => gmres_setup_xml
  !> Print solver information
  PROCEDURE :: view => gmres_view
  !> Setup PETSc solver object
  PROCEDURE :: setup_ksp => gmres_setup_ksp
end type oft_petsc_gmres_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_jacobi_solver
!---------------------------------------------------------------------------
!> PETSc Point-Jacobi solver
!!
!! Apply PETSc Point-Jacobi algorithm.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_jacobi_solver
  REAL(r8) :: df = .9d0 !< Damping factor
contains
  !> Solve system
  procedure :: apply => jacobi_solver_apply
end type oft_petsc_jacobi_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_sjacobi_solver
!---------------------------------------------------------------------------
!> PETSc symmetric Point-Jacobi solver
!!
!! Apply symmetric version of Point-Jacobi algorithm for PETSc matrices
!!
!! @sa oft_solvers::oft_jblock_precond
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_sjacobi_solver
  REAL(r8) :: df = .9d0 !< Damping factor
  logical, private :: down=.TRUE. !< Internal flag for smoother step
  logical, private :: warn_once = .TRUE. !< Internal flag for iteration check
  class(oft_vector), private, pointer :: D => NULL() !< Internal storage
  class(oft_vector), private, pointer :: u_save => NULL() !< Internal storage
  class(oft_vector), private, pointer :: g_save => NULL() !< Internal storage
contains
  !> Solve system
  procedure :: apply => sjacobi_solver_apply
  !> Setup solver from XML node
  procedure :: setup_from_xml => sjacobi_setup_xml
end type oft_petsc_sjacobi_solver
!---------------------------------------------------------------------------
! TYPE oft_petsc_direct_solver
!---------------------------------------------------------------------------
!> PETSc direct solver
!!
!! Use PETSc direct solver, LU decomposition, to invert system. Actual solve
!! employs a GMRES method preconditioned by an LU factorization.
!---------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_direct_solver
  TYPE(oft_petsc_factordef) :: pre_factor !< PETSc factor definition
contains
  !> Solve system
  procedure :: apply => direct_solver_apply
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => direct_setup_xml
  !> Print solver information
  PROCEDURE :: view => direct_view
  !> Setup PETSc solver object
  PROCEDURE :: setup_ksp => direct_setup_ksp
end type oft_petsc_direct_solver
!---Start subroutines
contains
!---------------------------------------------------------------------------
! FUNCTION: petsc_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_solver "oft_petsc_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_solver_cast(self,source) result(ierr)
CLASS(oft_petsc_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  class is(oft_petsc_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_solver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: petsc_solver_setup
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
subroutine petsc_solver_setup(self)
class(oft_petsc_solver), intent(inout) :: self
!---
CLASS(oft_petsc_matrix), pointer :: mat,pmat
CLASS(oft_petsc_precond), POINTER :: pre
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr) :: pc
#else
TYPE(tpc) :: pc
#endif
INTEGER(i4) :: ierr
STACK_PUSH
self%solver=PETSC_NULL_KSP
IF(oft_petsc_matrix_cast(mat,self%A)<0) &
  CALL oft_abort('"A" is not a PETSc matrix object.','petsc_solver_setup',__FILE__)
!---
pmat=>mat
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(petsc_precond_cast(pre,self%pre)<0) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','petsc_solver_setup',__FILE__)
  IF(ASSOCIATED(pre%A))THEN
    IF(oft_petsc_matrix_cast(pmat,pre%A)<0) &
      CALL oft_abort('"pre%A" is not a PETSc matrix object.','petsc_solver_setup',__FILE__)
  ELSE
    pre%A=>self%A
  END IF
END IF
IF(mat%nr==mat%nrg.AND.mat%nc==mat%ncg)THEN
  CALL KSPCreate(PETSC_COMM_SELF,self%solver,ierr)
  self%dist=.FALSE.
ELSE
  CALL KSPCreate(oft_env%COMM,self%solver,ierr)
  self%dist=.TRUE.
END IF
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
#else
CALL KSPSetOperators(self%solver,mat%M,pmat%M,SAME_NONZERO_PATTERN,ierr)
#endif
CALL KSPSetInitialGuessNonzero(self%solver,.TRUE.,ierr) ! Default to non-zero guess
!---
CALL KSPGetPC(self%solver,pc,ierr)
CALL PCSetType(pc,PCNONE,ierr)
STACK_POP
end subroutine petsc_solver_setup
!---------------------------------------------------------------------------
! SUBROUTINE: petsc_solver_update
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! @param[in] new_pattern Update matrix pattern (optional)
!---------------------------------------------------------------------------
recursive subroutine petsc_solver_update(self,new_pattern)
class(oft_petsc_solver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
!---
CLASS(oft_petsc_matrix), pointer :: mat,pmat
INTEGER(i4) :: ierr
LOGICAL :: update_pattern
STACK_PUSH
IF(self%initialized)THEN
  !---Update matrix references
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','petsc_solver_update',__FILE__)
  !---Preconditioner
  pmat=>mat
  IF(ASSOCIATED(self%pre))THEN
    IF(ASSOCIATED(self%pre%A))THEN
      IF(oft_petsc_matrix_cast(pmat,self%pre%A)<0) &
        CALL oft_abort('"pre%A" is not a PETSc matrix object.','petsc_solver_update',__FILE__)
    END IF
  END IF
  update_pattern=.FALSE.
  IF(PRESENT(new_pattern))update_pattern=new_pattern
  IF(update_pattern)THEN
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
#else
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,DIFFERENT_NONZERO_PATTERN,ierr)
#endif
  ELSE
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
#else
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,SAME_NONZERO_PATTERN,ierr)
#endif
  END IF
END IF
STACK_POP
end subroutine petsc_solver_update
!---------------------------------------------------------------------------
! SUBROUTINE: petsc_solver_setup_ksp
!---------------------------------------------------------------------------
!> Setup PETSc ksp object
!---------------------------------------------------------------------------
subroutine petsc_solver_setup_ksp(self,ksp)
class(oft_petsc_solver), intent(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), intent(inout) :: ksp
#else
TYPE(tksp), intent(inout) :: ksp
#endif
end subroutine petsc_solver_setup_ksp
!------------------------------------------------------------------------------
! SUBROUTINE: petsc_solver_delete
!------------------------------------------------------------------------------
!> Delete PETSc solver
!------------------------------------------------------------------------------
subroutine petsc_solver_delete(self)
class(oft_petsc_solver), intent(inout) :: self
integer(i4) :: ierr
STACK_PUSH
IF(self%initialized)CALL KSPDestroy(self%solver,ierr)
NULLIFY(self%A)
self%initialized=.FALSE.
STACK_POP
end subroutine petsc_solver_delete
!---------------------------------------------------------------------------
! FUNCTION: petsc_pre_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_pre_solver "oft_petsc_pre_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_pre_solver_cast(self,source) result(ierr)
type(oft_petsc_pre_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_pre_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_pre_solver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: petsc_solver_setup_pm
!---------------------------------------------------------------------------
!> Setup PETSc ksp object
!---------------------------------------------------------------------------
subroutine petsc_solver_setup_pm(ksp)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: ksp
#else
TYPE(tksp), INTENT(inout) :: ksp
#endif
INTEGER(i4) :: ierr
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<7)
  CALL KSPMonitorSet(ksp,KSPMonitorDefault,PETSC_NULL_FUNCTION,PETSC_NULL_FUNCTION,ierr)
#else
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr) :: vf
#else
TYPE(tPetscViewer) :: vf
#endif
  CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr)
  CALL KSPMonitorSet(ksp,KSPMonitorDefault,vf,PetscViewerAndFormatDestroy,ierr)
#endif
end subroutine petsc_solver_setup_pm
!------------------------------------------------------------------------------
! SUBROUTINE: pre_solver_apply
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine pre_solver_apply(self,u,g)
class(oft_petsc_pre_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_petsc_vector), pointer :: uv,gv
CLASS(oft_petsc_matrix), pointer :: mat
CLASS(oft_petsc_precond), POINTER :: pre
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr) :: pc
#else
TYPE(tpc) :: pc
#endif
LOGICAL :: dist
STACK_PUSH
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','pre_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','pre_solver_apply',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  self%solver=PETSC_NULL_KSP
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','pre_solver_apply',__FILE__)
  IF(mat%nr==mat%nrg.AND.mat%nc==mat%ncg)THEN
    CALL KSPCreate(PETSC_COMM_SELF,self%solver,ierr)
    dist=.FALSE.
  ELSE
    CALL KSPCreate(oft_env%COMM,self%solver,ierr)
    dist=.TRUE.
  END IF
  CALL KSPSetType(self%solver,KSPPREONLY,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
  CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
#else
  CALL KSPSetOperators(self%solver,mat%M,mat%M,SAME_PRECONDITIONER,ierr)
#endif
  !---
  IF(ASSOCIATED(self%pre))THEN
    IF(petsc_precond_cast(pre,self%pre)<0) &
      CALL oft_abort('"pre" is not a PETSc preconditioner object.','pre_solver_apply',__FILE__)
    CALL KSPSetUp(self%solver,ierr)
    CALL KSPGetPC(self%solver,pc,ierr)
    CALL pre%setup(pc,dist)
  END IF
  self%initialized=.TRUE.
END IF
!---
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
STACK_POP
end subroutine pre_solver_apply
!---------------------------------------------------------------------------
! FUNCTION: petsc_cg_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_cg_solver "oft_petsc_cg_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_cg_solver_cast(self,source) result(ierr)
type(oft_petsc_cg_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_cg_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_cg_solver_cast
!------------------------------------------------------------------------------
! SUBROUTINE: cg_solver_apply
!------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the Conjugate-Gradient
!! method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!------------------------------------------------------------------------------
subroutine cg_solver_apply(self,u,g)
class(oft_petsc_cg_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), POINTER :: uv,gv
STACK_PUSH
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','cg_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','cg_solver_apply',__FILE__)
self%cits=0
!---
IF(.NOT.self%initialized)THEN
  CALL petsc_solver_setup(self)
  CALL self%setup_ksp(self%solver)
  self%initialized=.TRUE.
END IF
!---Apply solver
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting PETSc CG solver'
IF(oft_env%pm)CALL petsc_solver_setup_pm(self%solver)
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
CALL KSPGetIterationNumber(self%solver,self%cits,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
#endif
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
!---Get residual
IF(self%full_residual)THEN
  CALL g%new(tmp)
  CALL self%A%apply(u,tmp)
  CALL g%add(1.d0,-1.d0,tmp)
  CALL tmp%delete
  DEALLOCATE(tmp)
END IF
STACK_POP
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
CLASS(oft_petsc_cg_solver), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),atol(val_level),rtol(val_level))
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","cg_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"atol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,atol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough atol values specified","cg_setup_xml",__FILE__)
    self%atol=atol(val_level)
  ELSE
    self%atol=atol(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"rtol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,rtol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough rtol values specified","cg_setup_xml",__FILE__)
    self%rtol=rtol(val_level)
  ELSE
    self%rtol=rtol(1)
  END IF
END IF
!---
IF(oft_debug_print(1))THEN
  WRITE(*,*)'CG solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Tolerance:   ',REAL(self%atol,4),REAL(self%rtol,4)
END IF
DEALLOCATE(its,atol,rtol)
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','cg_setup_xml',__FILE__)
#endif
end subroutine cg_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: cg_view
!---------------------------------------------------------------------------
!> Print solver configuration
!---------------------------------------------------------------------------
subroutine cg_view(self)
class(oft_petsc_cg_solver), intent(inout) :: self
!---
WRITE(*,*)'PETSc CG:'
WRITE(*,*)'  Max its = ',self%its
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine cg_view
!---------------------------------------------------------------------------
! SUBROUTINE: cg_setup_ksp
!---------------------------------------------------------------------------
!> Setup PETSc ksp object
!---------------------------------------------------------------------------
subroutine cg_setup_ksp(self,ksp)
CLASS(oft_petsc_cg_solver), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: ksp
INTEGER(petsc_addr) :: pc
#else
TYPE(tksp), INTENT(inout) :: ksp
TYPE(tpc) :: pc
#endif
INTEGER(i4) :: its,ierr
CLASS(oft_petsc_precond), POINTER :: pre
STACK_PUSH
!---
CALL KSPSetType(ksp,KSPCG,ierr)
its=PETSC_DEFAULT_INTEGER
IF(self%its>0)its=self%its
CALL KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED,ierr)
CALL KSPSetTolerances(ksp,self%rtol,self%atol, &
PETSC_DEFAULT_REAL,its,ierr)
!---
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(petsc_precond_cast(pre,self%pre)<0) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','cg_solver_apply',__FILE__)
END IF
!---
IF(ASSOCIATED(pre))THEN
  CALL KSPSetUp(ksp,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL pre%setup(pc,self%dist)
END IF
STACK_POP
end subroutine cg_setup_ksp
!---------------------------------------------------------------------------
! FUNCTION: petsc_gmres_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_gmres_solver "oft_petsc_gmres_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_gmres_solver_cast(self,source) result(ierr)
type(oft_petsc_gmres_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_gmres_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_gmres_solver_cast
!------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the FGMRES method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!------------------------------------------------------------------------------
subroutine gmres_solver_apply(self,u,g)
class(oft_petsc_gmres_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), pointer :: uv,gv
STACK_PUSH
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','gmres_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','gmres_solver_apply',__FILE__)
self%cits=0
IF(self%its==0.OR.self%nrits==0)THEN
  STACK_POP
  RETURN
END IF
!---
IF(.NOT.self%initialized)THEN
  CALL petsc_solver_setup(self)
  CALL self%setup_ksp(self%solver)
  self%initialized=.TRUE.
END IF
!---
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting PETSc GMRES solver'
IF(oft_env%pm)CALL petsc_solver_setup_pm(self%solver)
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
CALL KSPGetIterationNumber(self%solver,self%cits,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
#endif
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
!---Get residual
IF(self%full_residual)THEN
  CALL g%new(tmp)
  CALL self%A%apply(u,tmp)
  CALL g%add(1.d0,-1.d0,tmp)
  CALL tmp%delete
  DEALLOCATE(tmp)
END IF
STACK_POP
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
CLASS(oft_petsc_gmres_solver), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its,nrits
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),nrits(val_level),atol(val_level),rtol(val_level))
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","gmres_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"nrits",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,nrits,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough nrits values specified","gmres_setup_xml",__FILE__)
    self%nrits=nrits(val_level)
  ELSE
    self%nrits=nrits(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"atol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,atol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough atol values specified","gmres_setup_xml",__FILE__)
    self%atol=atol(val_level)
  ELSE
    self%atol=atol(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"rtol",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,rtol,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough rtol values specified","gmres_setup_xml",__FILE__)
    self%rtol=rtol(val_level)
  ELSE
    self%rtol=rtol(1)
  END IF
END IF
!---
IF(oft_debug_print(1))THEN
  WRITE(*,*)'GMRES solver setup:'
  WRITE(*,*)' - Iterations:  ',self%its
  WRITE(*,*)' - Restart:     ',self%nrits
  WRITE(*,*)' - Tolerance:   ',REAL(self%atol,4),REAL(self%rtol,4)
END IF
DEALLOCATE(its,nrits,atol,rtol)
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','gmres_setup_xml',__FILE__)
#endif
end subroutine gmres_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: gmres_view
!---------------------------------------------------------------------------
!> Print solver configuration
!---------------------------------------------------------------------------
subroutine gmres_view(self)
class(oft_petsc_gmres_solver), intent(inout) :: self
!---
WRITE(*,*)'PETSc GMRES:'
WRITE(*,*)'  Max its = ',self%its
WRITE(*,*)'  Restart = ',self%nrits
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine gmres_view
!---------------------------------------------------------------------------
! SUBROUTINE: gmres_setup_ksp
!---------------------------------------------------------------------------
!> Setup PETSc ksp object
!---------------------------------------------------------------------------
subroutine gmres_setup_ksp(self,ksp)
class(oft_petsc_gmres_solver), intent(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: ksp
INTEGER(petsc_addr) :: pc
#else
TYPE(tksp), INTENT(inout) :: ksp
TYPE(tpc) :: pc
#endif
!---
INTEGER(i4) :: its,ierr
CLASS(oft_petsc_precond), POINTER :: pre
STACK_PUSH
!---
SELECT CASE(self%algorithm)
  CASE(1)
    CALL KSPSetType(ksp,KSPFGMRES,ierr)
  CASE(2)
    CALL KSPSetType(ksp,KSPGMRES,ierr)
  CASE(3)
    CALL KSPSetType(ksp,KSPPGMRES,ierr)
  CASE DEFAULT
    CALL oft_abort('Invalid algorithm type','gmres_setup_ksp',__FILE__)
END SELECT
CALL KSPGMRESSetRestart(ksp,self%nrits,ierr)
its=PETSC_DEFAULT_INTEGER
IF(self%its>0)its=self%its
CALL KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED,ierr)
CALL KSPSetTolerances(ksp,self%rtol,self%atol, &
PETSC_DEFAULT_REAL,its,ierr)
!---
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(petsc_precond_cast(pre,self%pre)<0) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','gmres_solver_apply',__FILE__)
END IF
!---
IF(ASSOCIATED(pre))THEN
  CALL KSPSetUp(ksp,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL pre%setup(pc,self%dist)
END IF
STACK_POP
end subroutine gmres_setup_ksp
!---------------------------------------------------------------------------
! FUNCTION: petsc_jacobi_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_jacobi_solver "oft_petsc_jacobi_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_jacobi_solver_cast(self,source) result(ierr)
type(oft_petsc_jacobi_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_jacobi_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_jacobi_solver_cast
!------------------------------------------------------------------------------
! SUBROUTINE: jacobi_solver_apply
!------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the Point-Jacobi method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!------------------------------------------------------------------------------
subroutine jacobi_solver_apply(self,u,g)
class(oft_petsc_jacobi_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_petsc_vector), pointer :: uv,gv
CLASS(oft_petsc_matrix), pointer :: mat,pmat
CLASS(oft_petsc_precond), POINTER :: pre
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr) :: pc
#else
TYPE(tpc) :: pc
#endif
LOGICAL :: dist
STACK_PUSH
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','jacobi_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','jacobi_solver_apply',__FILE__)
self%cits=0
!---
IF(.NOT.self%initialized)THEN
  self%solver=PETSC_NULL_KSP
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','jacobi_solver_apply',__FILE__)
  !---
  pmat=>mat
  NULLIFY(pre)
  IF(ASSOCIATED(self%pre))THEN
    IF(petsc_precond_cast(pre,self%pre)<0) &
      CALL oft_abort('"pre" is not a PETSc preconditioner object.','jacobi_solver_apply',__FILE__)
    IF(ASSOCIATED(pre%A))THEN
      IF(oft_petsc_matrix_cast(pmat,pre%A)<0) &
        CALL oft_abort('"pre%A" is not a PETSc matrix object.','jacobi_solver_apply',__FILE__)
    END IF
  END IF
  !---
  IF(mat%nr==mat%nrg.AND.mat%nc==mat%ncg)THEN
    CALL KSPCreate(PETSC_COMM_SELF,self%solver,ierr)
    dist=.FALSE.
  ELSE
    CALL KSPCreate(oft_env%COMM,self%solver,ierr)
    dist=.TRUE.
  END IF
  CALL KSPSetType(self%solver,KSPRICHARDSON,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
  CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
#else
  CALL KSPSetOperators(self%solver,mat%M,pmat%M,SAME_NONZERO_PATTERN,ierr)
#endif
  !---
  IF(ASSOCIATED(pre))THEN
    CALL KSPSetUp(self%solver,ierr)
    CALL KSPGetPC(self%solver,pc,ierr)
    CALL pre%setup(pc,dist)
  ELSE
    CALL KSPGetPC(self%solver,pc,ierr)
    CALL PCSetType(pc,PCJACOBI,ierr)
  END IF
  !---
  its=PETSC_DEFAULT_INTEGER
  IF(self%its>0)its=self%its
  CALL KSPSetTolerances(self%solver,self%rtol,self%atol, &
  PETSC_DEFAULT_REAL,its,ierr)
  self%initialized=.TRUE.
END IF
!---
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting PETSc Jacobi solver'
IF(oft_env%pm)CALL petsc_solver_setup_pm(self%solver)
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
CALL KSPGetIterationNumber(self%solver,self%cits,ierr)
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
STACK_POP
end subroutine jacobi_solver_apply
!---------------------------------------------------------------------------
! FUNCTION: petsc_sjacobi_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_sjacobi_solver "oft_petsc_sjacobi_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_sjacobi_solver_cast(self,source) result(ierr)
type(oft_petsc_sjacobi_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_sjacobi_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_sjacobi_solver_cast
!---------------------------------------------------------------------------
! SUBROUTINE: sjacobi_solver_apply
!---------------------------------------------------------------------------
!> Apply 1-step of a symmetric Jacobi smoother with PETSc matrices
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
subroutine sjacobi_solver_apply(self,u,g)
class(oft_petsc_sjacobi_solver), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), intent(inout) :: g
class(oft_vector), pointer :: p,q
CLASS(oft_petsc_vector), pointer :: Vec
CLASS(oft_petsc_matrix), pointer :: Mat
integer(i4) :: i,k,ierr
STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch','sjacobi_solver_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch','sjacobi_solver_apply',__FILE__)
if(.NOT.associated(self%u_save))call u%new(self%u_save)
if(.NOT.associated(self%g_save))call u%new(self%g_save)
call u%new(p)
call u%new(q)
!---
IF(oft_petsc_matrix_cast(Mat,self%A)<0)CALL oft_abort('"Mat%D" is not a vector object.', &
'jblock_precond_apply',__FILE__)
IF(.NOT.ASSOCIATED(self%D))THEN
  CALL u%new(self%D)
  IF(oft_petsc_vector_cast(Vec,self%D)<0)CALL oft_abort('"Mat%D" is not a vector object.', &
  'jblock_precond_apply',__FILE__)
  CALL MatGetDiagonal(Mat%M,Vec%v,ierr)
  CALL VecReciprocal(Vec%v,ierr)
END IF
!---
IF(self%warn_once.AND.self%its<=0)THEN
  CALL oft_warn("SymJacobi smoother called with zero iterations [suppressed].")
  self%warn_once=.FALSE.
END IF
if(self%down) then ! Down smoother
  CALL u%set(0.d0)
  CALL p%add(0.d0,self%df,g)
  CALL p%mult(self%D)
  IF(ASSOCIATED(self%bc))call self%bc(p)
  do k=1,self%its
    call self%A%apply(p,q)
    IF(ASSOCIATED(self%bc))call self%bc(q)
    CALL u%add(1.d0,1.d0,p)
    CALL g%add(1.d0,-1.d0,q)
    CALL p%add(0.d0,self%df,g)
    CALL p%mult(self%D)
    IF(ASSOCIATED(self%bc))call self%bc(p)
  end do
  CALL self%u_save%add(0.d0,1.d0,u)
  CALL self%g_save%add(0.d0,1.d0,g)
  self%down=.FALSE.
else ! Up smoother
  CALL p%add(0.d0,1.d0,u,-1.d0,self%u_save)
  CALL g%add(0.d0,1.d0,self%g_save)
  IF(ASSOCIATED(self%bc))call self%bc(p)
  do k=1,self%its
    call self%A%apply(p,q)
    IF(ASSOCIATED(self%bc))call self%bc(q)
    CALL g%add(1.d0,-1.d0,q)
    CALL p%add(0.d0,self%df,g)
    CALL p%mult(self%D)
    IF(ASSOCIATED(self%bc))call self%bc(p)
    CALL u%add(1.d0,1.d0,p)
  end do
  self%down=.TRUE.
end if
call p%delete
call q%delete
nullify(p,q)
STACK_POP
end subroutine sjacobi_solver_apply
!---------------------------------------------------------------------------
! SUBROUTINE: sjacobi_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine sjacobi_setup_xml(self,solver_node,level)
CLASS(oft_petsc_sjacobi_solver), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(its(val_level),df(val_level))
!---
CALL xml_get_element(solver_node,"its",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,its,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough its values specified","jblock_setup_xml",__FILE__)
    self%its=its(val_level)
  ELSE
    self%its=its(1)
  END IF
END IF
!---
CALL xml_get_element(solver_node,"df",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,df,num=nread,iostat=ierr)
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
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','sjacobi_setup_xml',__FILE__)
#endif
end subroutine sjacobi_setup_xml
!---------------------------------------------------------------------------
! FUNCTION: petsc_direct_solver_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_direct_solver "oft_petsc_direct_solver"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_direct_solver_cast(self,source) result(ierr)
type(oft_petsc_direct_solver), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_direct_solver)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_direct_solver_cast
!------------------------------------------------------------------------------
! SUBROUTINE: direct_solver_apply
!------------------------------------------------------------------------------
!> Solve a linear system using a direct solver through PETSc
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!------------------------------------------------------------------------------
subroutine direct_solver_apply(self,u,g)
class(oft_petsc_direct_solver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---Local variables
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), pointer :: uv,gv
STACK_PUSH
!---Get solver fields
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','direct_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','direct_solver_apply',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  CALL petsc_solver_setup(self)
  CALL self%setup_ksp(self%solver)
  self%initialized=.TRUE.
END IF
!---Apply solver
IF((oft_env%pm.AND.oft_env%head_proc))WRITE(*,*)'Starting PETSc LU solver'
IF(oft_env%pm)CALL petsc_solver_setup_pm(self%solver)
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
#endif
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
!---Get residual
IF(self%full_residual)THEN
  CALL g%new(tmp)
  CALL self%A%apply(u,tmp)
  CALL g%add(1.d0,-1.d0,tmp)
  CALL tmp%delete
  DEALLOCATE(tmp)
END IF
STACK_POP
end subroutine direct_solver_apply
!---------------------------------------------------------------------------
! SUBROUTINE: direct_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine direct_setup_xml(self,solver_node,level)
CLASS(oft_petsc_direct_solver), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node,sub_node
CHARACTER(LEN=20) :: sub_type
CHARACTER(LEN=6) :: factor_type,factor_package
!---
INTEGER(i4) :: val_level,ierr
STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL lu_pc_load_xml(self%pre_factor,solver_node,val_level)
IF(oft_debug_print(1))THEN
  WRITE(*,*)'LU solver setup:'
  WRITE(*,*)' - Factor type:    ',self%pre_factor%type
  WRITE(*,*)' - Factor package: ',self%pre_factor%package
END IF
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','direct_setup_xml',__FILE__)
#endif
end subroutine direct_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: direct_view
!---------------------------------------------------------------------------
!> Print solver configuration
!---------------------------------------------------------------------------
subroutine direct_view(self)
class(oft_petsc_direct_solver), intent(inout) :: self
!---
WRITE(*,*)'PETSc LU:'
WRITE(*,*)'  Factor Type    = ',self%pre_factor%type
WRITE(*,*)'  Factor Package = ',self%pre_factor%package
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine direct_view
!---------------------------------------------------------------------------
! SUBROUTINE: direct_setup_ksp
!---------------------------------------------------------------------------
!> Setup PETSc ksp object
!---------------------------------------------------------------------------
subroutine direct_setup_ksp(self,ksp)
class(oft_petsc_direct_solver), intent(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: ksp
INTEGER(petsc_addr) :: pc
#else
TYPE(tksp), INTENT(inout) :: ksp
TYPE(tpc) :: pc
#endif
INTEGER(i4) :: its,ierr
STACK_PUSH
!---
CALL KSPSetType(ksp,KSPPREONLY,ierr)
CALL KSPGetPC(ksp,pc,ierr)
CALL create_lu_pc(self%pre_factor,pc)
its=PETSC_DEFAULT_INTEGER
IF(self%its>0)its=self%its
CALL KSPSetTolerances(ksp,self%rtol,self%atol, &
PETSC_DEFAULT_REAL,its,ierr)
STACK_POP
end subroutine direct_setup_ksp
!---------------------------------------------------------------------------
! FUNCTION: petsc_precond_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_precond "oft_petsc_precond"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_precond_cast(self,source) result(ierr)
class(oft_petsc_precond), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  class is(oft_petsc_precond)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_precond_cast
!---------------------------------------------------------------------------
! SUBROUTINE: precond_update
!---------------------------------------------------------------------------
!> Update preconditioner after changing settings/operators
!---------------------------------------------------------------------------
recursive subroutine precond_update(self,new_pattern)
class(oft_petsc_precond), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
!---
CLASS(oft_petsc_matrix), pointer :: mat
INTEGER(i4) :: ierr
LOGICAL :: update_pattern
STACK_PUSH
IF(self%initialized)THEN
  !---Update matrix references
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','precond_update',__FILE__)
    update_pattern=.FALSE.
  IF(PRESENT(new_pattern))update_pattern=new_pattern
  IF(update_pattern)THEN
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
#else
    CALL KSPSetOperators(self%solver,mat%M,mat%M,DIFFERENT_NONZERO_PATTERN,ierr)
#endif
  ELSE
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
#else
    CALL KSPSetOperators(self%solver,mat%M,mat%M,SAME_NONZERO_PATTERN,ierr)
#endif
  END IF
END IF
STACK_POP
end subroutine precond_update
!------------------------------------------------------------------------------
! SUBROUTINE: precond_apply
!------------------------------------------------------------------------------
!> Precondition a linear system using one iteration of a PETSc preconditioner
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!------------------------------------------------------------------------------
subroutine precond_apply(self,u,g)
class(oft_petsc_precond), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u,g
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_petsc_vector), pointer :: uv,gv
CLASS(oft_petsc_matrix), pointer :: mat
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr) :: pc
#else
TYPE(tpc) :: pc
#endif
LOGICAL :: dist
STACK_PUSH
IF(oft_petsc_vector_cast(uv,u)<0)CALL oft_abort('"u" is not a PETSc vector object.','pre_solver_apply',__FILE__)
IF(oft_petsc_vector_cast(gv,g)<0)CALL oft_abort('"g" is not a PETSc vector object.','pre_solver_apply',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  self%solver=PETSC_NULL_KSP
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','pre_solver_apply',__FILE__)
  IF(mat%nr==mat%nrg.AND.mat%nc==mat%ncg)THEN
    CALL KSPCreate(PETSC_COMM_SELF,self%solver,ierr)
    dist=.FALSE.
  ELSE
    CALL KSPCreate(oft_env%COMM,self%solver,ierr)
    dist=.TRUE.
  END IF
  CALL KSPSetType(self%solver,KSPPREONLY,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
  CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
#else
  CALL KSPSetOperators(self%solver,mat%M,mat%M,SAME_NONZERO_PATTERN,ierr)
#endif
  CALL KSPSetUp(self%solver,ierr)
  CALL KSPGetPC(self%solver,pc,ierr)
  CALL self%setup(pc,dist)
  self%initialized=.TRUE.
END IF
!---
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
#endif
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
STACK_POP
end subroutine precond_apply
!------------------------------------------------------------------------------
! SUBROUTINE: precond_delete
!------------------------------------------------------------------------------
!> Delete PETSc solver
!------------------------------------------------------------------------------
subroutine precond_delete(self)
class(oft_petsc_precond), intent(inout) :: self
integer(i4) :: ierr
STACK_PUSH
IF(self%initialized)CALL KSPDestroy(self%solver,ierr)
NULLIFY(self%A)
self%initialized=.FALSE.
STACK_POP
end subroutine precond_delete
!------------------------------------------------------------------------------
! SUBROUTINE: precond_dummy_setup
!------------------------------------------------------------------------------
!> Setup preconditioner
!!
!! @param[in,out] pc PETSc preconditioner object from parent KSP
!! @param[in] dist Flag for local vs distributed solve
!------------------------------------------------------------------------------
subroutine precond_dummy_setup(self,pc,dist)
CLASS(oft_petsc_precond), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
#else
TYPE(tpc), INTENT(inout) :: pc
#endif
LOGICAL, INTENT(in) :: dist
CALL oft_abort('No preconditioner type set.','precond_dummy_setup',__FILE__)
end subroutine precond_dummy_setup
!---------------------------------------------------------------------------
! FUNCTION: petsc_diagprecond_cast
!---------------------------------------------------------------------------
!> Cast @ref oft_solvers::oft_solver "oft_solver" to
!! @ref oft_petsc_solvers::oft_petsc_diagprecond "oft_petsc_diagprecond"
!!
!! @param[out] self Object of desired type, unassociated if cast fails
!! @param[in] source Source object to cast
!! @result Error flag
!---------------------------------------------------------------------------
FUNCTION petsc_diagprecond_cast(self,source) result(ierr)
type(oft_petsc_diagprecond), pointer, intent(out) :: self
class(oft_solver), target, intent(in) :: source
integer(i4) :: ierr
STACK_PUSH
select type(source)
  type is(oft_petsc_diagprecond)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
STACK_POP
end FUNCTION petsc_diagprecond_cast
!------------------------------------------------------------------------------
! SUBROUTINE: diagprecond_setup
!------------------------------------------------------------------------------
!> @copydoc oft_petsc_solvers::precond_dummy_setup
!------------------------------------------------------------------------------
subroutine diagprecond_setup(self,pc,dist)
CLASS(oft_petsc_diagprecond), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
#else
TYPE(tpc), INTENT(inout) :: pc
#endif
LOGICAL, INTENT(in) :: dist
!---Local variables
INTEGER(i4) :: ierr
STACK_PUSH
!---Basic preconditioner setup
CALL PCSetType(pc,PCJACOBI,ierr)
STACK_POP
end subroutine diagprecond_setup
!---------------------------------------------------------------------------
! SUBROUTINE: diagprecond_view
!---------------------------------------------------------------------------
!> Print solver configuration
!---------------------------------------------------------------------------
subroutine diagprecond_view(self)
class(oft_petsc_diagprecond), intent(inout) :: self
!---
WRITE(*,*)'PETSc Point-Jacobi'
end subroutine diagprecond_view
!------------------------------------------------------------------------------
! SUBROUTINE: luprecond_setup
!------------------------------------------------------------------------------
!> @copydoc oft_petsc_solvers::precond_dummy_setup
!------------------------------------------------------------------------------
subroutine luprecond_setup(self,pc,dist)
CLASS(oft_petsc_luprecond), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
#else
TYPE(tpc), INTENT(inout) :: pc
#endif
LOGICAL, INTENT(in) :: dist
!---Local variables
INTEGER(i4) :: ierr
STACK_PUSH
!---Basic preconditioner setup
CALL PCSetType(pc,PCLU,ierr)
!---Set solver package
CALL set_solver_package(pc,OFT_PETSC_LU,ierr)
STACK_POP
end subroutine luprecond_setup
!------------------------------------------------------------------------------
! SUBROUTINE: mgprecond_setup
!------------------------------------------------------------------------------
!> @copydoc oft_petsc_solvers::precond_dummy_setup
!------------------------------------------------------------------------------
subroutine mgprecond_setup(self,pc,dist)
CLASS(oft_petsc_mgprecond), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
INTEGER(petsc_addr) :: ksp,ksps(4),pc_local,sub_pc
#else
TYPE(tpc), INTENT(inout) :: pc
TYPE(tksp) :: ksp,ksps(4)
TYPE(tpc) :: pc_local,sub_pc
#endif
LOGICAL, INTENT(in) :: dist
!---Local variables
INTEGER(i4) :: i,j,n_local,ierr
CHARACTER(LEN=20) :: prefix
STACK_PUSH
!------------------------------------------------------------------------------
! Basic preconditioner setup
!------------------------------------------------------------------------------
CALL PCSetType(pc,PCMG,ierr)
CALL PCMGSetLevels(pc,self%nlevels,PETSC_NULL_INTEGER,ierr)
CALL PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
!---Set galerkin
IF(self%galerkin)THEN
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<10)
  CALL PCMGSetGalerkin(pc,PETSC_TRUE,ierr)
#else
  CALL PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH,ierr)
#endif
END IF
!---Set cycle type
IF(self%cycle_type=="v".OR.self%cycle_type=="V")THEN
  CALL PCMGSetCycleType(pc,PC_MG_CYCLE_V,ierr)
ELSE IF(self%cycle_type=="w".OR.self%cycle_type=="W")THEN
  CALL PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
ELSE
  CALL oft_abort('Invalid cycle type.','mgprecond_setup',__FILE__)
END IF
!------------------------------------------------------------------------------
! Loop over levels
!------------------------------------------------------------------------------
DO i=1,self%nlevels
  IF(i>1)THEN
    !---Set interpolation matrices
    SELECT TYPE(this=>self%I(i-1)%M)
      CLASS is(oft_petsc_matrix)
        CALL PCMGSetInterpolation(pc,i-1,this%M,ierr)
    END SELECT
  END IF
  !---Setup level smoother
  CALL PCMGGetSmoother(pc,i-1,ksp,ierr)
  !---Set level matrices
  IF(.NOT.self%galerkin)THEN
    SELECT TYPE(this=>self%Mats(i)%M)
      CLASS IS(oft_petsc_matrix)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>4)
        CALL KSPSetOperators(ksp,this%M,this%M,ierr)
#else
        CALL KSPSetOperators(ksp,this%M,this%M,SAME_NONZERO_PATTERN,ierr)
#endif
    END SELECT
  END IF
  !---
  IF(i==1.AND.self%direct_base)THEN
    CALL KSPSetType(ksp,KSPGMRES,ierr)
    CALL KSPGetPC(ksp,pc_local,ierr)
    CALL PCSetType(pc_local,PCLU,ierr)
    CALL set_solver_package(pc_local,OFT_PETSC_LU,ierr)
  ELSE
    SELECT CASE(TRIM(self%smooth_type))
      CASE("jacobi")
        CALL KSPSetType(ksp,KSPRICHARDSON,ierr)
        CALL KSPRichardsonSetScale(ksp,self%df(i),ierr)
        CALL KSPSetTolerances(ksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
        PETSC_DEFAULT_REAL,self%nu(i),ierr)
        CALL KSPGetPC(ksp,pc_local,ierr)
        CALL PCSetType(pc_local,PCJACOBI,ierr)
        ! CALL PCMGSetCyclesOnLevel(pc,i-1,1,ierr)
      CASE("gmres")
        CALL KSPSetType(ksp,KSPFGMRES,ierr)
        CALL KSPGMRESSetRestart(ksp,self%nu(i),ierr)
        CALL KSPSetTolerances(ksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
        PETSC_DEFAULT_REAL,self%nu(i),ierr)
        CALL KSPGetPC(ksp,pc_local,ierr)
        CALL PCSetType(pc_local,PCJACOBI,ierr)
        ! CALL PCMGSetCyclesOnLevel(pc,i-1,1,ierr)
      CASE("gmres-as")
        CALL KSPSetType(ksp,KSPGMRES,ierr)
        CALL KSPGMRESSetRestart(ksp,self%nu(i),ierr)
        CALL KSPSetTolerances(ksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
        PETSC_DEFAULT_REAL,self%nu(i),ierr)
        CALL KSPGetPC(ksp,pc_local,ierr)
        !---
        CALL PCSetType(pc_local,PCASM,ierr)
        CALL PCSetUp(pc_local,ierr)
        CALL PCASMGetSubKSP(pc_local,n_local,PETSC_NULL_INTEGER,ksps,ierr)
        DO j=1,n_local
          CALL KSPSetType(ksps(j),KSPPREONLY,ierr)
          CALL KSPGetPC(ksps(j),sub_pc,ierr)
          CALL PCSetType(sub_pc,PCLU,ierr)
          CALL set_solver_package(sub_pc,OFT_PETSC_LU,ierr)
        END DO
        ! CALL PCMGSetCyclesOnLevel(pc,i-1,1,ierr)
      CASE DEFAULT
        CALL oft_abort('Invalid smoother type.','mgprecond_setup',__FILE__)
    END SELECT
  END IF
END DO
CALL PCMGGetCoarseSolve(pc,ksp,ierr)
!CALL PCView(pc,PETSC_VIEWER_STDOUT_SELF,ierr)
STACK_POP
end subroutine mgprecond_setup
!------------------------------------------------------------------------------
! SUBROUTINE: asprecond_setup
!------------------------------------------------------------------------------
!> @copydoc oft_petsc_solvers::precond_dummy_setup
!------------------------------------------------------------------------------
subroutine asprecond_setup(self,pc,dist)
CLASS(oft_petsc_asprecond), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
INTEGER(petsc_addr) :: pc_local
INTEGER(petsc_addr), ALLOCATABLE, DIMENSION(:) :: ksp,is_parts
#else
TYPE(tpc), INTENT(inout) :: pc
TYPE(tpc) :: pc_local
TYPE(tksp), ALLOCATABLE, DIMENSION(:) :: ksp
TYPE(tis), ALLOCATABLE, DIMENSION(:) :: is_parts
#endif
LOGICAL, INTENT(in) :: dist
!---
INTEGER(i4) :: i,ierr,n_local,loc_start,loc_end
CLASS(oft_petsc_matrix), pointer :: mat
CLASS(oft_petsc_solver), POINTER :: pre
STACK_PUSH
!---
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(petsc_solver_cast(pre,self%pre)<0) &
    CALL oft_abort('"pre" is not a PETSc solver object.','asprecond_setup',__FILE__)
END IF
!---
CALL PCSetType(pc,PCASM,ierr)
CALL PCASMSetOverlap(pc,self%overlap,ierr)
!---
IF(self%n_local==-1)THEN
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','asprecond_setup',__FILE__)
  n_local=mat%ni
  ALLOCATE(is_parts(n_local))
  CALL MatGetOwnershipRange(mat%M,loc_start,loc_end,ierr)
  DO i=1,mat%ni
    loc_end = mat%i_map(i)%nslice
    CALL ISCreateStride(PETSC_COMM_SELF,loc_end,loc_start,1,is_parts(i),ierr)
    loc_start=loc_start+loc_end
  END DO
  CALL PCASMSetLocalSubdomains(pc,n_local,is_parts,is_parts,ierr)
  DEALLOCATE(is_parts)
ELSE IF(self%n_local>1)THEN
  IF(oft_petsc_matrix_cast(mat,self%A)<0) &
    CALL oft_abort('"A" is not a PETSc matrix object.','asprecond_setup',__FILE__)
  n_local=self%n_local
  ALLOCATE(is_parts(n_local))
  CALL MatGetOwnershipRange(mat%M,loc_start,loc_end,ierr)
  DO i=1,self%n_local
    loc_end = FLOOR(mat%nrslice/REAL(self%n_local,8))
    IF(i==self%n_local)loc_end = loc_end + MOD(mat%nrslice,self%n_local)
    CALL ISCreateStride(PETSC_COMM_SELF,loc_end,loc_start,1,is_parts(i),ierr)
    loc_start=loc_start+loc_end
  END DO
  CALL PCASMSetLocalSubdomains(pc,n_local,is_parts,is_parts,ierr)
  DEALLOCATE(is_parts)
ELSE
  n_local=1
END IF
!---
CALL PCSetUp(pc,ierr)
ALLOCATE(ksp(n_local))
CALL PCASMGetSubKSP(pc,n_local,PETSC_NULL_INTEGER,ksp,ierr)
DO i=1,n_local
  IF(ASSOCIATED(pre))THEN
    CALL pre%setup_ksp(ksp(i))
  ELSE
    CALL KSPSetType(ksp(i),KSPPREONLY,ierr)
    CALL KSPGetPC(ksp(i),pc_local,ierr)
    CALL create_lu_pc(self%sub_factor,pc_local)
  END IF
END DO
DEALLOCATE(ksp)
STACK_POP
end subroutine asprecond_setup
!---------------------------------------------------------------------------
! SUBROUTINE: asprecond_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
subroutine asprecond_setup_xml(self,solver_node,level)
CLASS(oft_petsc_asprecond), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node,sub_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: overlaps,nlocals
STACK_PUSH
!---
val_level=1
IF(PRESENT(level))val_level=level
ALLOCATE(overlaps(val_level),nlocals(val_level))
!---Read in overlap size
CALL xml_get_element(solver_node,"overlap",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,overlaps,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough overlap values specified","asprecond_setup_xml",__FILE__)
    self%overlap=overlaps(val_level)
  ELSE
    self%overlap=overlaps(1)
  END IF
END IF
!---Read in local field splitting flag
CALL xml_get_element(solver_node,"nlocal",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,nlocals,num=nread,iostat=ierr)
  IF(nread>1)THEN
    IF(ierr<0)CALL oft_abort("Not enough local flags specified","asprecond_setup_xml",__FILE__)
    self%n_local=nlocals(val_level)
  ELSE
    self%n_local=nlocals(1)
  END IF
END IF
IF(self%n_local<-1)CALL oft_abort("Slice grouping not supported with PETSc","asprecond_setup_xml",__FILE__)
IF(self%n_local==-1)self%overlap=0
!---Read-in desired overlap specification
CALL xml_get_element(solver_node,"boverlap",current_node,ierr)
IF(ierr==0)CALL oft_abort("Boundary overlap not supported with PETSc","asprecond_setup_xml",__FILE__)
IF(oft_debug_print(1))THEN
  WRITE(*,*)'Additive Schwartz solver setup:'
  WRITE(*,*)' - Overlap:    ',self%overlap
  WRITE(*,*)' - NLocal:     ',self%n_local
END IF
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','asprecond_setup_xml',__FILE__)
#endif
end subroutine asprecond_setup_xml
!------------------------------------------------------------------------------
! SUBROUTINE: set_solver_package
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine set_solver_package(pc,sol_type,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
#else
TYPE(tpc), INTENT(inout) :: pc
#endif
CHARACTER(LEN=*), INTENT(IN) :: sol_type
INTEGER(i4), INTENT(out) :: ierr
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<9)
CALL PCFactorSetMatSolverPackage(pc,sol_type,ierr)
#else
CALL PCFactorSetMatSolverType(pc,sol_type,ierr)
#endif
end subroutine set_solver_package
!------------------------------------------------------------------------------
! SUBROUTINE: create_lu_pc
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine create_lu_pc(self,pc)
CLASS(oft_petsc_factordef), INTENT(inout) :: self
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<8)
INTEGER(petsc_addr), INTENT(inout) :: pc
#else
TYPE(tpc), INTENT(inout) :: pc
#endif
INTEGER(i4) :: i,ierr
LOGICAL :: factor_req
STACK_PUSH
!---
SELECT CASE(TRIM(self%type))
  CASE('lu')
    CALL PCSetType(pc,PCLU,ierr)
    factor_req=.TRUE.
  CASE('ilu')
    CALL PCSetType(pc,PCILU,ierr)
    factor_req=.TRUE.
  CASE DEFAULT
    CALL oft_abort("Unknown sub solver.","create_lu_pc",__FILE__)
END SELECT
!---
IF(factor_req)THEN
  SELECT CASE(TRIM(self%package))
    CASE('petsc')
      CALL set_solver_package(pc,MATSOLVERPETSC,ierr)
    CASE('super')
      CALL set_solver_package(pc,MATSOLVERSUPERLU,ierr)
    CASE('superd')
      CALL set_solver_package(pc,MATSOLVERSUPERLU_DIST,ierr)
    CASE('mumps')
      CALL set_solver_package(pc,MATSOLVERMUMPS,ierr)
    CASE DEFAULT
      CALL oft_abort("Unknown solver package.","create_lu_pc",__FILE__)
  END SELECT
END IF
STACK_POP
end subroutine create_lu_pc
!------------------------------------------------------------------------------
! SUBROUTINE: lu_pc_load_xml
!------------------------------------------------------------------------------
!> Setup LU preconditioner object from XML definition
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!------------------------------------------------------------------------------
subroutine lu_pc_load_xml(self,solver_node,level)
CLASS(oft_petsc_factordef), INTENT(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: i,ierr,nnodes,nread
TYPE(xml_node), POINTER :: current_node
CHARACTER(LEN=6) :: factor_type,factor_package
STACK_PUSH
!---
CALL xml_get_element(solver_node,"type",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,factor_type,num=nread,iostat=ierr)
  IF(nread==1)THEN
    self%type=factor_type
  END IF
END IF
!---
CALL xml_get_element(solver_node,"package",current_node,ierr)
IF(ierr==0)THEN
  CALL xml_extractDataContent(current_node,factor_package,num=nread,iostat=ierr)
  IF(nread==1)THEN
    self%package=factor_package
  END IF
END IF
STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','lu_pc_load_xml',__FILE__)
#endif
end subroutine lu_pc_load_xml
#endif
end module oft_petsc_solvers
