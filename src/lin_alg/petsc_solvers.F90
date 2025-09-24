!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
MODULE oft_petsc_solvers
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc, solver_delete
#ifdef HAVE_PETSC
USE oft_petsc_la, ONLY: oft_petsc_vector, oft_petsc_vector_cast, &
  oft_petsc_matrix, oft_petsc_matrix_cast
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
USE petscmat
USE petscksp
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
!------------------------------------------------------------------------------
!> PETSc factorization package information
!------------------------------------------------------------------------------
type :: oft_petsc_factordef
  CHARACTER(LEN=6) :: type = 'lu' !< Factor type
  CHARACTER(LEN=6) :: package = 'superd' !< Factor package
end type oft_petsc_factordef
!------------------------------------------------------------------------------
!> PETSc abstract preconditioner class
!!
!! @note If called from a native solver a solver wrapper will be created
!------------------------------------------------------------------------------
type, abstract, public, extends(oft_solver) :: oft_petsc_precond
  TYPE(tpc) :: pre_obj !< PETSc preconditioner object
  TYPE(tksp) :: solver !< PETSc solver object
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
!------------------------------------------------------------------------------
!> PETSc multi-grid preconditioner
!!
!! Use PETSc multi-grid framework for preconditioning. V and W cycle types are
!! supported. Galerkin construction of sub-level matrices is also available.
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> PETSc diagonal preconditioning
!------------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_diagprecond
contains
  !> Setup preconditioner
  procedure :: setup => diagprecond_setup
  !> Print solver information
  PROCEDURE :: view => diagprecond_view
end type oft_petsc_diagprecond
!------------------------------------------------------------------------------
!> PETSc LU factorization/solve
!------------------------------------------------------------------------------
type, public, extends(oft_petsc_precond) :: oft_petsc_luprecond
contains
  !> Setup preconditioner
  procedure :: setup => luprecond_setup
end type oft_petsc_luprecond
!------------------------------------------------------------------------------
!> PETSc Additive-Schwarz preconditioner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> PETSc abstract solver class
!------------------------------------------------------------------------------
type, abstract, public, extends(oft_solver) :: oft_petsc_solver
  TYPE(tksp) :: solver !< PETSc solver object
  LOGICAL :: dist = .FALSE. !< Matrix/Vector are distributed
CONTAINS
  !> Update solver with new settings/operators
  PROCEDURE :: update => petsc_solver_update
  !> Setup PETSc solver object
  PROCEDURE :: setup_ksp => petsc_solver_setup_ksp
  !> Clean-up internal storage
  PROCEDURE :: delete => petsc_solver_delete
end type oft_petsc_solver
!------------------------------------------------------------------------------
!> PETSc Conjugate-Gradient solver
!!
!! @note Matrix must be SPD, otherwise solver will fail.
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> PETSc GMRES solver
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> PETSc Point-Jacobi solver
!------------------------------------------------------------------------------
type, public, extends(oft_petsc_solver) :: oft_petsc_jacobi_solver
  REAL(r8) :: df = .9d0 !< Damping factor
contains
  !> Solve system
  procedure :: apply => jacobi_solver_apply
end type oft_petsc_jacobi_solver
!------------------------------------------------------------------------------
!> PETSc symmetric Point-Jacobi solver
!!
!! @sa oft_native_solvers::oft_jblock_precond
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> PETSc direct solver (LU factorization)
!!
!! Use PETSc direct solver, LU decomposition, to invert system. Actual solve
!! employs a GMRES method preconditioned by an LU factorization.
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Cast a vector object to a oft_petsc_vector
!!
!! The source vector must be @ref oft_petsc_vector or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_solver_cast(self,source) result(success)
CLASS(oft_petsc_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_petsc_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_solver_cast
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!------------------------------------------------------------------------------
subroutine petsc_solver_setup(self)
class(oft_petsc_solver), intent(inout) :: self !< Solver object
CLASS(oft_petsc_matrix), pointer :: mat,pmat
CLASS(oft_petsc_precond), POINTER :: pre
TYPE(tpc) :: pc
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
self%solver=PETSC_NULL_KSP
IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
  CALL oft_abort('"A" is not a PETSc matrix object.','petsc_solver_setup',__FILE__)
!---
pmat=>mat
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.petsc_precond_cast(pre,self%pre)) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','petsc_solver_setup',__FILE__)
  IF(ASSOCIATED(pre%A))THEN
    IF(.NOT.oft_petsc_matrix_cast(pmat,pre%A)) &
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
CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
CALL KSPSetInitialGuessNonzero(self%solver,.TRUE.,ierr) ! Default to non-zero guess
!---
CALL KSPGetPC(self%solver,pc,ierr)
CALL PCSetType(pc,PCNONE,ierr)
DEBUG_STACK_POP
end subroutine petsc_solver_setup
!------------------------------------------------------------------------------
!> Update solver after updating settings/matrix
!------------------------------------------------------------------------------
recursive subroutine petsc_solver_update(self,new_pattern)
class(oft_petsc_solver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
CLASS(oft_petsc_matrix), pointer :: mat,pmat
INTEGER(i4) :: ierr
LOGICAL :: update_pattern
DEBUG_STACK_PUSH
IF(self%initialized)THEN
  !---Update matrix references
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
    CALL oft_abort('"A" is not a PETSc matrix object.','petsc_solver_update',__FILE__)
  !---Preconditioner
  pmat=>mat
  IF(ASSOCIATED(self%pre))THEN
    IF(ASSOCIATED(self%pre%A))THEN
      IF(.NOT.oft_petsc_matrix_cast(pmat,self%pre%A)) &
        CALL oft_abort('"pre%A" is not a PETSc matrix object.','petsc_solver_update',__FILE__)
    END IF
  END IF
  update_pattern=.FALSE.
  IF(PRESENT(new_pattern))update_pattern=new_pattern
  IF(update_pattern)THEN
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
  ELSE
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
  END IF
END IF
DEBUG_STACK_POP
end subroutine petsc_solver_update
!------------------------------------------------------------------------------
!> Setup PETSc ksp object
!------------------------------------------------------------------------------
subroutine petsc_solver_setup_ksp(self,ksp)
class(oft_petsc_solver), intent(inout) :: self !< Solver object
TYPE(tksp), intent(inout) :: ksp
end subroutine petsc_solver_setup_ksp
!---------------------------------------------------------------------------------
!> Delete PETSc solver
!---------------------------------------------------------------------------------
subroutine petsc_solver_delete(self,propogate)
class(oft_petsc_solver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
integer(i4) :: ierr
DEBUG_STACK_PUSH
IF(self%initialized)CALL KSPDestroy(self%solver,ierr)
NULLIFY(self%A)
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
DEBUG_STACK_POP
end subroutine petsc_solver_delete
!------------------------------------------------------------------------------
!> Setup performance monitor for PETSc solver
!------------------------------------------------------------------------------
subroutine petsc_solver_setup_pm(ksp)
TYPE(tksp), INTENT(inout) :: ksp
INTEGER(i4) :: ierr
TYPE(tPetscViewer) :: vf
CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr)
CALL KSPMonitorSet(ksp,KSPMonitorResidual,vf,PetscViewerAndFormatDestroy,ierr)
end subroutine petsc_solver_setup_pm
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_cg_solver
!!
!! The source solver must be @ref oft_petsc_cg_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_cg_solver_cast(self,source) result(success)
type(oft_petsc_cg_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_cg_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_cg_solver_cast
!---------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the Conjugate-Gradient
!! method
!---------------------------------------------------------------------------------
subroutine cg_solver_apply(self,u,g)
class(oft_petsc_cg_solver), intent(inout) :: self !< Solver object
CLASS(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), POINTER :: uv,gv
DEBUG_STACK_PUSH
IF(.NOT.oft_petsc_vector_cast(uv,u))CALL oft_abort('"u" is not a PETSc vector object.','cg_solver_apply',__FILE__)
IF(.NOT.oft_petsc_vector_cast(gv,g))CALL oft_abort('"g" is not a PETSc vector object.','cg_solver_apply',__FILE__)
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
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
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
DEBUG_STACK_POP
end subroutine cg_solver_apply
!------------------------------------------------------------------------------
!> Setup CG solver from XML definition
!------------------------------------------------------------------------------
subroutine cg_setup_xml(self,solver_node,level)
CLASS(oft_petsc_cg_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: atol,rtol
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','cg_setup_xml',__FILE__)
#endif
end subroutine cg_setup_xml
!------------------------------------------------------------------------------
!> Print solver configuration
!------------------------------------------------------------------------------
subroutine cg_view(self)
class(oft_petsc_cg_solver), intent(inout) :: self !< Solver object
!---
WRITE(*,*)'PETSc CG:'
WRITE(*,*)'  Max its = ',self%its
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine cg_view
!------------------------------------------------------------------------------
!> Setup PETSc ksp object for CG
!------------------------------------------------------------------------------
subroutine cg_setup_ksp(self,ksp)
CLASS(oft_petsc_cg_solver), INTENT(inout) :: self !< Solver object
TYPE(tksp), INTENT(inout) :: ksp !< PETSc KSP object
TYPE(tpc) :: pc
INTEGER(i4) :: its,ierr
CLASS(oft_petsc_precond), POINTER :: pre
DEBUG_STACK_PUSH
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
  IF(.NOT.petsc_precond_cast(pre,self%pre)) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','cg_solver_apply',__FILE__)
END IF
!---
IF(ASSOCIATED(pre))THEN
  CALL KSPSetUp(ksp,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL pre%setup(pc,self%dist)
END IF
DEBUG_STACK_POP
end subroutine cg_setup_ksp
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_gmres_solver
!!
!! The source solver must be @ref oft_petsc_gmres_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_gmres_solver_cast(self,source) result(success)
type(oft_petsc_gmres_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_gmres_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_gmres_solver_cast
!---------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the FGMRES method
!---------------------------------------------------------------------------------
subroutine gmres_solver_apply(self,u,g)
class(oft_petsc_gmres_solver), intent(inout) :: self !< Solver object
CLASS(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
CLASS(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), pointer :: uv,gv
DEBUG_STACK_PUSH
IF(.NOT.oft_petsc_vector_cast(uv,u))CALL oft_abort('"u" is not a PETSc vector object.','gmres_solver_apply',__FILE__)
IF(.NOT.oft_petsc_vector_cast(gv,g))CALL oft_abort('"g" is not a PETSc vector object.','gmres_solver_apply',__FILE__)
self%cits=0
IF(self%its==0.OR.self%nrits==0)THEN
  DEBUG_STACK_POP
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
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
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
DEBUG_STACK_POP
end subroutine gmres_solver_apply
!------------------------------------------------------------------------------
!> Setup FGMRES solver from XML definition
!------------------------------------------------------------------------------
subroutine gmres_setup_xml(self,solver_node,level)
CLASS(oft_petsc_gmres_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','gmres_setup_xml',__FILE__)
#endif
end subroutine gmres_setup_xml
!------------------------------------------------------------------------------
!> Print solver configuration
!------------------------------------------------------------------------------
subroutine gmres_view(self)
class(oft_petsc_gmres_solver), intent(inout) :: self !< Solver object
!---
WRITE(*,*)'PETSc GMRES:'
WRITE(*,*)'  Max its = ',self%its
WRITE(*,*)'  Restart = ',self%nrits
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine gmres_view
!------------------------------------------------------------------------------
!> Setup PETSc ksp object for FGMRES
!------------------------------------------------------------------------------
subroutine gmres_setup_ksp(self,ksp)
class(oft_petsc_gmres_solver), intent(inout) :: self !< Solver object
TYPE(tksp), INTENT(inout) :: ksp !< PETSc KSP object
TYPE(tpc) :: pc
INTEGER(i4) :: its,ierr
CLASS(oft_petsc_precond), POINTER :: pre
DEBUG_STACK_PUSH
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
  IF(.NOT.petsc_precond_cast(pre,self%pre)) &
    CALL oft_abort('"pre" is not a PETSc preconditioner object.','gmres_solver_apply',__FILE__)
END IF
!---
IF(ASSOCIATED(pre))THEN
  CALL KSPSetUp(ksp,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL pre%setup(pc,self%dist)
END IF
DEBUG_STACK_POP
end subroutine gmres_setup_ksp
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_jacobi_solver
!!
!! The source solver must be @ref oft_petsc_jacobi_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_jacobi_solver_cast(self,source) result(success)
type(oft_petsc_jacobi_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_jacobi_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_jacobi_solver_cast
!---------------------------------------------------------------------------------
!> Solve a linear system using PETSc's implementation of the Point-Jacobi method
!---------------------------------------------------------------------------------
subroutine jacobi_solver_apply(self,u,g)
class(oft_petsc_jacobi_solver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_petsc_vector), pointer :: uv,gv
CLASS(oft_petsc_matrix), pointer :: mat,pmat
CLASS(oft_petsc_precond), POINTER :: pre
TYPE(tpc) :: pc
LOGICAL :: dist
DEBUG_STACK_PUSH
IF(.NOT.oft_petsc_vector_cast(uv,u))CALL oft_abort('"u" is not a PETSc vector object.','jacobi_solver_apply',__FILE__)
IF(.NOT.oft_petsc_vector_cast(gv,g))CALL oft_abort('"g" is not a PETSc vector object.','jacobi_solver_apply',__FILE__)
self%cits=0
!---
IF(.NOT.self%initialized)THEN
  self%solver=PETSC_NULL_KSP
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
    CALL oft_abort('"A" is not a PETSc matrix object.','jacobi_solver_apply',__FILE__)
  !---
  pmat=>mat
  NULLIFY(pre)
  IF(ASSOCIATED(self%pre))THEN
    IF(.NOT.petsc_precond_cast(pre,self%pre)) &
      CALL oft_abort('"pre" is not a PETSc preconditioner object.','jacobi_solver_apply',__FILE__)
    IF(ASSOCIATED(pre%A))THEN
      IF(.NOT.oft_petsc_matrix_cast(pmat,pre%A)) &
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
  CALL KSPSetOperators(self%solver,mat%M,pmat%M,ierr)
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
DEBUG_STACK_POP
end subroutine jacobi_solver_apply
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_sjacobi_solver
!!
!! The source solver must be @ref oft_petsc_sjacobi_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_sjacobi_solver_cast(self,source) result(success)
type(oft_petsc_sjacobi_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_sjacobi_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_sjacobi_solver_cast
!------------------------------------------------------------------------------
!> Apply 1-step of a symmetric Jacobi smoother with PETSc matrices
!------------------------------------------------------------------------------
subroutine sjacobi_solver_apply(self,u,g)
class(oft_petsc_sjacobi_solver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
class(oft_vector), pointer :: p,q
CLASS(oft_petsc_vector), pointer :: Vec
CLASS(oft_petsc_matrix), pointer :: Mat
integer(i4) :: i,k,ierr
DEBUG_STACK_PUSH
if(g%n/=self%A%nr)call oft_abort('Row mismatch','sjacobi_solver_apply',__FILE__)
if(u%n/=self%A%nc)call oft_abort('Col mismatch','sjacobi_solver_apply',__FILE__)
if(.NOT.associated(self%u_save))call u%new(self%u_save)
if(.NOT.associated(self%g_save))call u%new(self%g_save)
call u%new(p)
call u%new(q)
!---
IF(.NOT.oft_petsc_matrix_cast(Mat,self%A))CALL oft_abort('"Mat%D" is not a vector object.', &
'jblock_precond_apply',__FILE__)
IF(.NOT.ASSOCIATED(self%D))THEN
  CALL u%new(self%D)
  IF(.NOT.oft_petsc_vector_cast(Vec,self%D))CALL oft_abort('"Mat%D" is not a vector object.', &
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
  IF(ASSOCIATED(self%bc))call self%bc%apply(p)
  do k=1,self%its
    call self%A%apply(p,q)
    IF(ASSOCIATED(self%bc))call self%bc%apply(q)
    CALL u%add(1.d0,1.d0,p)
    CALL g%add(1.d0,-1.d0,q)
    CALL p%add(0.d0,self%df,g)
    CALL p%mult(self%D)
    IF(ASSOCIATED(self%bc))call self%bc%apply(p)
  end do
  CALL self%u_save%add(0.d0,1.d0,u)
  CALL self%g_save%add(0.d0,1.d0,g)
  self%down=.FALSE.
else ! Up smoother
  CALL p%add(0.d0,1.d0,u,-1.d0,self%u_save)
  CALL g%add(0.d0,1.d0,self%g_save)
  IF(ASSOCIATED(self%bc))call self%bc%apply(p)
  do k=1,self%its
    call self%A%apply(p,q)
    IF(ASSOCIATED(self%bc))call self%bc%apply(q)
    CALL g%add(1.d0,-1.d0,q)
    CALL p%add(0.d0,self%df,g)
    CALL p%mult(self%D)
    IF(ASSOCIATED(self%bc))call self%bc%apply(p)
    CALL u%add(1.d0,1.d0,p)
  end do
  self%down=.TRUE.
end if
call p%delete
call q%delete
nullify(p,q)
DEBUG_STACK_POP
end subroutine sjacobi_solver_apply
!------------------------------------------------------------------------------
!> Setup symmetric Jacobi solver from XML definition
!------------------------------------------------------------------------------
subroutine sjacobi_setup_xml(self,solver_node,level)
CLASS(oft_petsc_sjacobi_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: its
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','sjacobi_setup_xml',__FILE__)
#endif
end subroutine sjacobi_setup_xml
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_direct_solver
!!
!! The source solver must be @ref oft_petsc_direct_solver or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_direct_solver_cast(self,source) result(success)
type(oft_petsc_direct_solver), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_direct_solver)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_direct_solver_cast
!---------------------------------------------------------------------------------
!> Solve a linear system using a direct solver through PETSc (LU factorization)
!---------------------------------------------------------------------------------
subroutine direct_solver_apply(self,u,g)
class(oft_petsc_direct_solver), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_matrix), pointer :: A
!---Local variables
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_vector), pointer :: tmp
CLASS(oft_petsc_vector), pointer :: uv,gv
DEBUG_STACK_PUSH
!---Get solver fields
IF(.NOT.oft_petsc_vector_cast(uv,u))CALL oft_abort('"u" is not a PETSc vector object.','direct_solver_apply',__FILE__)
IF(.NOT.oft_petsc_vector_cast(gv,g))CALL oft_abort('"g" is not a PETSc vector object.','direct_solver_apply',__FILE__)
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
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
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
DEBUG_STACK_POP
end subroutine direct_solver_apply
!------------------------------------------------------------------------------
!> Setup direct solver from XML definition
!------------------------------------------------------------------------------
subroutine direct_setup_xml(self,solver_node,level)
CLASS(oft_petsc_direct_solver), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node,sub_node
CHARACTER(LEN=20) :: sub_type
CHARACTER(LEN=6) :: factor_type,factor_package
!---
INTEGER(i4) :: val_level,ierr
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','direct_setup_xml',__FILE__)
#endif
end subroutine direct_setup_xml
!------------------------------------------------------------------------------
!> Print solver configuration
!------------------------------------------------------------------------------
subroutine direct_view(self)
class(oft_petsc_direct_solver), intent(inout) :: self !< Solver object
!---
WRITE(*,*)'PETSc LU:'
WRITE(*,*)'  Factor Type    = ',self%pre_factor%type
WRITE(*,*)'  Factor Package = ',self%pre_factor%package
IF(ASSOCIATED(self%pre))THEN
  WRITE(*,*)'Preconditioner:'
  CALL self%pre%view
END IF
end subroutine direct_view
!------------------------------------------------------------------------------
!> Setup PETSc ksp object for direct solver
!------------------------------------------------------------------------------
subroutine direct_setup_ksp(self,ksp)
class(oft_petsc_direct_solver), intent(inout) :: self !< Solver object
TYPE(tksp), INTENT(inout) :: ksp !< PETSc KSP solver object
TYPE(tpc) :: pc
INTEGER(i4) :: its,ierr
DEBUG_STACK_PUSH
!---
CALL KSPSetType(ksp,KSPPREONLY,ierr)
CALL KSPGetPC(ksp,pc,ierr)
CALL create_lu_pc(self%pre_factor,pc)
its=PETSC_DEFAULT_INTEGER
IF(self%its>0)its=self%its
CALL KSPSetTolerances(ksp,self%rtol,self%atol, &
PETSC_DEFAULT_REAL,its,ierr)
DEBUG_STACK_POP
end subroutine direct_setup_ksp
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_precond
!!
!! The source solver must be @ref oft_petsc_precond or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_precond_cast(self,source) result(success)
class(oft_petsc_precond), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_petsc_precond)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_precond_cast
!------------------------------------------------------------------------------
!> Update preconditioner after changing settings/operators
!------------------------------------------------------------------------------
recursive subroutine precond_update(self,new_pattern)
class(oft_petsc_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix pattern (optional)
!---
CLASS(oft_petsc_matrix), pointer :: mat
INTEGER(i4) :: ierr
LOGICAL :: update_pattern
DEBUG_STACK_PUSH
IF(self%initialized)THEN
  !---Update matrix references
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
    CALL oft_abort('"A" is not a PETSc matrix object.','precond_update',__FILE__)
    update_pattern=.FALSE.
  IF(PRESENT(new_pattern))update_pattern=new_pattern
  IF(update_pattern)THEN
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
  ELSE
    CALL KSPSetReusePreconditioner(self%solver,PETSC_FALSE,ierr)
    CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
  END IF
END IF
DEBUG_STACK_POP
end subroutine precond_update
!---------------------------------------------------------------------------------
!> Precondition a linear system using one iteration of a PETSc preconditioner
!---------------------------------------------------------------------------------
subroutine precond_apply(self,u,g)
class(oft_petsc_precond), intent(inout) :: self !< Solver object
class(oft_vector), intent(inout) :: u !< Guess (input), Solution (output)
class(oft_vector), intent(inout) :: g !< RHS (input), Residual (output)
CLASS(oft_matrix), pointer :: A
!---
INTEGER(i4) :: ierr,its,nlocal,first
REAL(r8) :: norm
CLASS(oft_petsc_vector), pointer :: uv,gv
CLASS(oft_petsc_matrix), pointer :: mat
TYPE(tpc) :: pc
LOGICAL :: dist
DEBUG_STACK_PUSH
IF(.NOT.oft_petsc_vector_cast(uv,u))CALL oft_abort('"u" is not a PETSc vector object.','pre_solver_apply',__FILE__)
IF(.NOT.oft_petsc_vector_cast(gv,g))CALL oft_abort('"g" is not a PETSc vector object.','pre_solver_apply',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  self%solver=PETSC_NULL_KSP
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
    CALL oft_abort('"A" is not a PETSc matrix object.','pre_solver_apply',__FILE__)
  IF(mat%nr==mat%nrg.AND.mat%nc==mat%ncg)THEN
    CALL KSPCreate(PETSC_COMM_SELF,self%solver,ierr)
    dist=.FALSE.
  ELSE
    CALL KSPCreate(oft_env%COMM,self%solver,ierr)
    dist=.TRUE.
  END IF
  CALL KSPSetType(self%solver,KSPPREONLY,ierr)
  CALL KSPSetOperators(self%solver,mat%M,mat%M,ierr)
  CALL KSPSetUp(self%solver,ierr)
  CALL KSPGetPC(self%solver,pc,ierr)
  CALL self%setup(pc,dist)
  self%initialized=.TRUE.
END IF
!---
CALL KSPSolve(self%solver,gv%v,uv%v,ierr)
CALL KSPSetReusePreconditioner(self%solver,PETSC_TRUE,ierr)
uv%loc_current=.FALSE.
gv%loc_current=.FALSE.
DEBUG_STACK_POP
end subroutine precond_apply
!---------------------------------------------------------------------------------
!> Delete PETSc solver
!---------------------------------------------------------------------------------
subroutine precond_delete(self,propogate)
class(oft_petsc_precond), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
integer(i4) :: ierr
DEBUG_STACK_PUSH
IF(self%initialized)CALL KSPDestroy(self%solver,ierr)
NULLIFY(self%A)
self%initialized=.FALSE.
CALL solver_delete(self,propogate)
DEBUG_STACK_POP
end subroutine precond_delete
!---------------------------------------------------------------------------------
!> Setup preconditioner
!---------------------------------------------------------------------------------
subroutine precond_dummy_setup(self,pc,dist)
CLASS(oft_petsc_precond), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
LOGICAL, INTENT(in) :: dist !< Flag for local vs distributed solve?
CALL oft_abort('No preconditioner type set.','precond_dummy_setup',__FILE__)
end subroutine precond_dummy_setup
!------------------------------------------------------------------------------
!> Cast a solver object to a oft_petsc_diagprecond
!!
!! The source solver must be @ref oft_petsc_diagprecond or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION petsc_diagprecond_cast(self,source) result(success)
type(oft_petsc_diagprecond), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_solver), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  type is(oft_petsc_diagprecond)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION petsc_diagprecond_cast
!---------------------------------------------------------------------------------
!> Setup diagonal preconditioner
!---------------------------------------------------------------------------------
subroutine diagprecond_setup(self,pc,dist)
CLASS(oft_petsc_diagprecond), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
LOGICAL, INTENT(in) :: dist !< Flag for local vs distributed solve?
!---Local variables
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---Basic preconditioner setup
CALL PCSetType(pc,PCJACOBI,ierr)
DEBUG_STACK_POP
end subroutine diagprecond_setup
!------------------------------------------------------------------------------
!> Print solver configuration
!------------------------------------------------------------------------------
subroutine diagprecond_view(self)
class(oft_petsc_diagprecond), intent(inout) :: self !< Solver object
WRITE(*,*)'PETSc Point-Jacobi'
end subroutine diagprecond_view
!---------------------------------------------------------------------------------
!> Setup LU factorization preconditioner
!---------------------------------------------------------------------------------
subroutine luprecond_setup(self,pc,dist)
CLASS(oft_petsc_luprecond), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
LOGICAL, INTENT(in) :: dist !< Flag for local vs distributed solve?
!---Local variables
INTEGER(i4) :: ierr
DEBUG_STACK_PUSH
!---Basic preconditioner setup
CALL PCSetType(pc,PCLU,ierr)
!---Set solver package
CALL set_solver_package(pc,OFT_PETSC_LU,ierr)
DEBUG_STACK_POP
end subroutine luprecond_setup
!---------------------------------------------------------------------------------
!> Setup multi-grid preconditioner
!---------------------------------------------------------------------------------
subroutine mgprecond_setup(self,pc,dist)
CLASS(oft_petsc_mgprecond), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
LOGICAL, INTENT(in) :: dist !< Flag for local vs distributed solve?
TYPE(tksp) :: ksp,ksps(4)
TYPE(tpc) :: pc_local,sub_pc
!---Local variables
INTEGER(i4) :: i,j,n_local,ierr
CHARACTER(LEN=20) :: prefix
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------------
! Basic preconditioner setup
!---------------------------------------------------------------------------------
CALL PCSetType(pc,PCMG,ierr)
CALL PCMGSetLevels(pc,self%nlevels,PETSC_NULL_INTEGER,ierr)
CALL PCMGSetType(pc,PC_MG_MULTIPLICATIVE,ierr)
!---Set galerkin
IF(self%galerkin)THEN
  CALL PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH,ierr)
END IF
!---Set cycle type
IF(self%cycle_type=="v".OR.self%cycle_type=="V")THEN
  CALL PCMGSetCycleType(pc,PC_MG_CYCLE_V,ierr)
ELSE IF(self%cycle_type=="w".OR.self%cycle_type=="W")THEN
  CALL PCMGSetCycleType(pc,PC_MG_CYCLE_W,ierr)
ELSE
  CALL oft_abort('Invalid cycle type.','mgprecond_setup',__FILE__)
END IF
!---------------------------------------------------------------------------------
! Loop over levels
!---------------------------------------------------------------------------------
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
        CALL KSPSetOperators(ksp,this%M,this%M,ierr)
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
DEBUG_STACK_POP
end subroutine mgprecond_setup
!---------------------------------------------------------------------------------
!> Setup additiv Schwarz preconditioner
!---------------------------------------------------------------------------------
subroutine asprecond_setup(self,pc,dist)
CLASS(oft_petsc_asprecond), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
LOGICAL, INTENT(in) :: dist !< Flag for local vs distributed solve?
TYPE(tpc) :: pc_local
TYPE(tksp), ALLOCATABLE, DIMENSION(:) :: ksp
TYPE(tis), ALLOCATABLE, DIMENSION(:) :: is_parts
!---
INTEGER(i4) :: i,ierr,n_local,loc_start,loc_end
CLASS(oft_petsc_matrix), pointer :: mat
CLASS(oft_petsc_solver), POINTER :: pre
DEBUG_STACK_PUSH
!---
NULLIFY(pre)
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.petsc_solver_cast(pre,self%pre)) &
    CALL oft_abort('"pre" is not a PETSc solver object.','asprecond_setup',__FILE__)
END IF
!---
CALL PCSetType(pc,PCASM,ierr)
CALL PCASMSetOverlap(pc,self%overlap,ierr)
!---
IF(self%n_local==-1)THEN
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
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
  IF(.NOT.oft_petsc_matrix_cast(mat,self%A)) &
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
DEBUG_STACK_POP
end subroutine asprecond_setup
!------------------------------------------------------------------------------
!> Setup additive Schwarz preconditioner from XML definition
!------------------------------------------------------------------------------
subroutine asprecond_setup_xml(self,solver_node,level)
CLASS(oft_petsc_asprecond), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: nnodes,nread
TYPE(xml_node), POINTER :: current_node,sub_node
!---
INTEGER(i4) :: val_level,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: overlaps,nlocals
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','asprecond_setup_xml',__FILE__)
#endif
end subroutine asprecond_setup_xml
!---------------------------------------------------------------------------------
!> Set factorization package for PETSc PC object
!---------------------------------------------------------------------------------
subroutine set_solver_package(pc,sol_type,ierr)
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
CHARACTER(LEN=*), INTENT(IN) :: sol_type !< Solver type
INTEGER(i4), INTENT(out) :: ierr !< Error flag
CALL PCFactorSetMatSolverType(pc,sol_type,ierr)
end subroutine set_solver_package
!---------------------------------------------------------------------------------
!> Create PETSc LU PC object
!---------------------------------------------------------------------------------
subroutine create_lu_pc(self,pc)
CLASS(oft_petsc_factordef), INTENT(inout) :: self !< Solver object
TYPE(tpc), INTENT(inout) :: pc !< PETSc PC object
INTEGER(i4) :: i,ierr
LOGICAL :: factor_req
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
end subroutine create_lu_pc
!---------------------------------------------------------------------------------
!> Setup LU preconditioner object from XML definition
!---------------------------------------------------------------------------------
subroutine lu_pc_load_xml(self,solver_node,level)
CLASS(oft_petsc_factordef), INTENT(inout) :: self !< Solver object
TYPE(xml_node), POINTER, INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
#ifdef HAVE_XML
!---
INTEGER(i4) :: i,ierr,nnodes,nread
TYPE(xml_node), POINTER :: current_node
CHARACTER(LEN=6) :: factor_type,factor_package
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','lu_pc_load_xml',__FILE__)
#endif
end subroutine lu_pc_load_xml
#endif
end module oft_petsc_solvers
