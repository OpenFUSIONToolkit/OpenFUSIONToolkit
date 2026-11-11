!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file solver_base.F90
!
!> Abstract solver interfaces
!!
!! - Linear solver classes (@ref oft_solver, @ref oft_csolver)
!! - Eigenvalue solver class (@ref oft_eigsolver)
!! - Boundary condition classes (@ref oft_solver_bc, @ref oft_csolver_bc)
!!
!! See @ref doc_solvers for available implementations and their options.
!!
!! @sa oft_native_solvers, oft_petsc_solvers, oft_lu, oft_solver_utils
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_native_la
!------------------------------------------------------------------------------
MODULE oft_solver_base
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix, oft_cmatrix
IMPLICIT NONE
#include "local.h"
private
!------------------------------------------------------------------------------
!> Base class for solver boundary conditions
!!
!! Applied by solvers to the RHS, residual and search directions to enforce
!! boundary conditions (eg. zero Dirichlet) or remove null-space components.
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_solver_bc
CONTAINS
  !> Apply boundary condition to field
  PROCEDURE(oft_bc_proto), DEFERRED :: apply
  !> Clean-up internal storage
  procedure(oft_bc_delete), deferred :: delete
END TYPE oft_solver_bc
!------------------------------------------------------------------------------
!> Base class for solver boundary conditions (complex)
!!
!! @sa oft_solver_bc
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_csolver_bc
CONTAINS
  !> Apply boundary condition to field
  PROCEDURE(oft_cbc_proto), DEFERRED :: apply
  !> Clean-up internal storage
  procedure(oft_cbc_delete), deferred :: delete
END TYPE oft_csolver_bc
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Abstract boundary condition prototype
  !------------------------------------------------------------------------------
  SUBROUTINE oft_bc_proto(self,a)
    IMPORT oft_solver_bc, oft_vector
    CLASS(oft_solver_bc), INTENT(inout) :: self !< Boundary condition object
    CLASS(oft_vector), INTENT(inout) :: a !< Field to apply BC to
  END SUBROUTINE oft_bc_proto
  !------------------------------------------------------------------------------
  !> Destroy boundary condition structure and deallocate all internal storage
  !------------------------------------------------------------------------------
  SUBROUTINE oft_bc_delete(self)
    IMPORT oft_solver_bc
    CLASS(oft_solver_bc), INTENT(inout) :: self !< Boundary condition object
  END SUBROUTINE oft_bc_delete
  !------------------------------------------------------------------------------
  !> Abstract boundary condition prototype (complex)
  !------------------------------------------------------------------------------
  SUBROUTINE oft_cbc_proto(self,a)
    IMPORT oft_csolver_bc, oft_cvector
    CLASS(oft_csolver_bc), INTENT(inout) :: self !< Boundary condition object
    CLASS(oft_cvector), INTENT(inout) :: a !< Field to apply BC to
  END SUBROUTINE oft_cbc_proto
  !------------------------------------------------------------------------------
  !> Destroy boundary condition structure and deallocate all internal storage (complex)
  !------------------------------------------------------------------------------
  SUBROUTINE oft_cbc_delete(self)
    IMPORT oft_csolver_bc
    CLASS(oft_csolver_bc), INTENT(inout) :: self !< Boundary condition object
  END SUBROUTINE oft_cbc_delete
END INTERFACE
!------------------------------------------------------------------------------
!> Base class for linear solvers and preconditioners
!!
!! Defines the common interface and options shared by all solvers. See
!! @ref doc_solvers_usage for usage and @ref doc_solvers_common for the meaning
!! of the common options.
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_solver
  LOGICAL :: pm = .FALSE. !< Print convergence information
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  LOGICAL :: full_residual = .TRUE. !< Compute true residual on exit
  INTEGER(i4) :: its = -1 !< Iteration limit (see @ref doc_solvers_common)
  INTEGER(i4) :: cits = 0 !< Number of iterations taken in last solve (<0 indicates failure)
  INTEGER(i4) :: itplot=10 !< Output frequency for iterative solvers when pm=.TRUE.
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_matrix), POINTER :: A => NULL() !< Matrix to be inverted
  CLASS(oft_solver), POINTER :: pre => NULL() !< Preconditioner (or sub-solver for block preconditioners)
  CLASS(oft_solver_bc), POINTER :: bc => NULL() !< Boundary condition
CONTAINS
  !> Solve system
  PROCEDURE(solver_apply), DEFERRED :: apply
  !> Update solver with new settings/operators
  PROCEDURE :: update => solver_update
  !> Setup solver from XML node (see @ref doc_solvers_xml)
  PROCEDURE :: setup_from_xml => solver_setup_xml
  !> Print solver information
  PROCEDURE :: view => solver_view
  !> Check thread safety
  PROCEDURE :: check_thread => solver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => solver_delete
END TYPE oft_solver
!------------------------------------------------------------------------------
!> Solver pointer container for building arrays of solvers
!------------------------------------------------------------------------------
type, public :: oft_solver_ptr
  CLASS(oft_solver), POINTER :: s => NULL() !< Solver object
end type oft_solver_ptr
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Apply linear solver to compute \f$ A u = g\f$
  !------------------------------------------------------------------------------
  SUBROUTINE solver_apply(self,u,g)
  IMPORT oft_solver, oft_vector
  CLASS(oft_solver), INTENT(inout) :: self !< Solver object
  CLASS(oft_vector), INTENT(inout) :: u !< Guess (input), Solution (output)
  CLASS(oft_vector), INTENT(inout) :: g !< RHS (input), overwritten with residual for iterative solvers (output)
  END SUBROUTINE solver_apply
END INTERFACE
!------------------------------------------------------------------------------
!> Base class for linear solvers and preconditioners (complex)
!!
!! @sa oft_solver
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_csolver
  LOGICAL :: pm = .FALSE. !< Print convergence information
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  LOGICAL :: full_residual = .TRUE. !< Compute true residual on exit
  INTEGER(i4) :: its = -1 !< Iteration limit (see @ref doc_solvers_common)
  INTEGER(i4) :: cits = 0 !< Number of iterations taken in last solve (<0 indicates failure)
  INTEGER(i4) :: itplot=10 !< Output frequency for iterative solvers when pm=.TRUE.
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_cmatrix), POINTER :: A => NULL() !< Matrix to be inverted
  CLASS(oft_csolver), POINTER :: pre => NULL() !< Preconditioner (or sub-solver for block preconditioners)
  CLASS(oft_csolver_bc), POINTER :: bc => NULL() !< Boundary condition
CONTAINS
  !> Solve system
  PROCEDURE(csolver_apply), DEFERRED :: apply
  !> Update solver with new settings/operators
  PROCEDURE :: update => csolver_update
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => csolver_setup_xml
  !> Print solver information
  PROCEDURE :: view => csolver_view
  !> Check thread safety
  PROCEDURE :: check_thread => csolver_check_thread
  !> Clean-up internal storage
  PROCEDURE :: delete => csolver_delete
END TYPE oft_csolver
!------------------------------------------------------------------------------
!> Solver pointer container for building arrays of solvers (complex)
!------------------------------------------------------------------------------
type, public :: oft_csolver_ptr
  CLASS(oft_csolver), POINTER :: s => NULL() !< Solver object
end type oft_csolver_ptr
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Apply linear solver to compute \f$ A u = g\f$ (complex)
  !------------------------------------------------------------------------------
  SUBROUTINE csolver_apply(self,u,g)
  IMPORT oft_csolver, oft_cvector
  CLASS(oft_csolver), INTENT(inout) :: self !< Solver object
  CLASS(oft_cvector), INTENT(inout) :: u !< Guess (input), Solution (output)
  CLASS(oft_cvector), INTENT(inout) :: g !< RHS (input), overwritten with residual for iterative solvers (output)
  END SUBROUTINE csolver_apply
END INTERFACE
!------------------------------------------------------------------------------
!> Base class for eigenvalue solvers of \f$ A u = \lambda M u \f$
!!
!! See @ref doc_solvers_eig for available implementations.
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_eigsolver
  LOGICAL :: pm = .FALSE. !< Print convergence information
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  INTEGER(i4) :: its = -1 !< Iteration limit (meaning depends on implementation)
  INTEGER(i4) :: cits = 0 !< Number of iterations taken in last solve
  INTEGER(i4) :: itplot=10 !< Output frequency for iterative solvers when pm=.TRUE.
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_matrix), POINTER :: A => NULL() !< LHS matrix
  CLASS(oft_matrix), POINTER :: M => NULL() !< RHS matrix
  CLASS(oft_solver), POINTER :: pre => NULL() !< Preconditioner
CONTAINS
  !> Solve Eigen-system
  PROCEDURE(eigsolver_apply), DEFERRED :: apply
  !> Clean-up internal storage
  PROCEDURE(eigsolver_delete), DEFERRED :: delete
END TYPE oft_eigsolver
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Apply eigensystem solver to compute eigenvalues/eigenvectors of the system
  !! \f$ A u_i = \lambda_i M u_i \f$
  !------------------------------------------------------------------------------
  SUBROUTINE eigsolver_apply(self,u,alam)
  IMPORT oft_eigsolver, oft_vector, r8
  CLASS(oft_eigsolver), INTENT(inout) :: self !< Solver object
  CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
  REAL(r8), INTENT(inout) :: alam !< Eigenvalue
  END SUBROUTINE eigsolver_apply
  !------------------------------------------------------------------------------
  !> Destroy eigensystem solver and deallocate all internal storage
  !------------------------------------------------------------------------------
  SUBROUTINE eigsolver_delete(self)
  IMPORT oft_eigsolver
  CLASS(oft_eigsolver), INTENT(inout) :: self !< Solver object
  END SUBROUTINE eigsolver_delete
END INTERFACE
!---Make classes and prototypes public
PUBLIC solver_setup, solver_delete, csolver_setup, csolver_delete, eigsolver_setup
CONTAINS
!------------------------------------------------------------------------------
!> Common solver initialization, called by implementations on first application
!!
!! Sets the preconditioner matrix to the solver matrix if not already set
!------------------------------------------------------------------------------
SUBROUTINE solver_setup(self)
CLASS(oft_solver), INTENT(inout) :: self !< Solver object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE solver_setup
!------------------------------------------------------------------------------
!> Destroy linear solver and deallocate all internal storage
!!
!! Base implementation only handles the preconditioner, called by
!! implementations after freeing their own storage
!------------------------------------------------------------------------------
SUBROUTINE solver_delete(self,propogate)
CLASS(oft_solver), INTENT(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Also delete preconditioner and sub-solvers? (optional)
IF(PRESENT(propogate))THEN
  IF(propogate.AND.ASSOCIATED(self%pre))THEN
    CALL self%pre%delete(propogate)
    DEALLOCATE(self%pre)
  END IF
END IF
END SUBROUTINE solver_delete
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! Base implementation only propagates the update to the preconditioner
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE solver_update(self,new_pattern)
CLASS(oft_solver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix non-zero pattern? (optional)
IF(ASSOCIATED(self%pre))CALL self%pre%update(new_pattern)
END SUBROUTINE solver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! Base implementation reads no options. See @ref doc_solvers_xml for the
!! options supported by each implementation.
!------------------------------------------------------------------------------
SUBROUTINE solver_setup_xml(self,solver_node,level)
class(oft_solver), intent(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
end SUBROUTINE solver_setup_xml
!------------------------------------------------------------------------------
!> Print solver configuration (base implementation prints nothing)
!------------------------------------------------------------------------------
SUBROUTINE solver_view(self)
class(oft_solver), intent(inout) :: self !< Solver object
end SUBROUTINE solver_view
!------------------------------------------------------------------------------
!> Check if solver (and preconditioner) can be applied concurrently from multiple threads
!------------------------------------------------------------------------------
recursive function solver_check_thread(self) result(thread_safe)
class(oft_solver), intent(inout) :: self !< Solver object
logical :: thread_safe !< Solver is thread safe?
thread_safe=.TRUE.
IF(ASSOCIATED(self%pre))thread_safe=(thread_safe.AND.self%pre%check_thread())
end function solver_check_thread
!------------------------------------------------------------------------------
!> Common solver initialization, called by implementations on first application (complex)
!------------------------------------------------------------------------------
SUBROUTINE csolver_setup(self)
CLASS(oft_csolver), INTENT(inout) :: self !< Solver object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE csolver_setup
!------------------------------------------------------------------------------
!> Destroy linear solver and deallocate all internal storage (complex)
!------------------------------------------------------------------------------
SUBROUTINE csolver_delete(self,propogate)
CLASS(oft_csolver), INTENT(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: propogate !< Also delete preconditioner and sub-solvers? (optional)
IF(PRESENT(propogate))THEN
  IF(propogate.AND.ASSOCIATED(self%pre))THEN
    CALL self%pre%delete(propogate)
    DEALLOCATE(self%pre)
  END IF
END IF
END SUBROUTINE csolver_delete
!------------------------------------------------------------------------------
!> Update solver after changing settings/operators (complex)
!------------------------------------------------------------------------------
recursive SUBROUTINE csolver_update(self,new_pattern)
CLASS(oft_csolver), intent(inout) :: self !< Solver object
LOGICAL, optional, intent(in) :: new_pattern !< Update matrix non-zero pattern? (optional)
IF(ASSOCIATED(self%pre))CALL self%pre%update(new_pattern)
END SUBROUTINE csolver_update
!------------------------------------------------------------------------------
!> Setup solver from XML definition (complex, base implementation reads no options)
!------------------------------------------------------------------------------
SUBROUTINE csolver_setup_xml(self,solver_node,level)
class(oft_csolver), intent(inout) :: self !< Solver object
TYPE(xml_node), INTENT(in) :: solver_node !< XML element containing solver definition
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Level in MG hierarchy (optional)
end SUBROUTINE csolver_setup_xml
!------------------------------------------------------------------------------
!> Print solver configuration (complex, base implementation prints nothing)
!------------------------------------------------------------------------------
SUBROUTINE csolver_view(self)
class(oft_csolver), intent(inout) :: self !< Solver object
end SUBROUTINE csolver_view
!------------------------------------------------------------------------------
!> Check if solver (and preconditioner) can be applied concurrently from multiple threads (complex)
!------------------------------------------------------------------------------
recursive function csolver_check_thread(self) result(thread_safe)
class(oft_csolver), intent(inout) :: self !< Solver object
logical :: thread_safe !< Solver is thread safe?
thread_safe=.TRUE.
IF(ASSOCIATED(self%pre))thread_safe=(thread_safe.AND.self%pre%check_thread())
end function csolver_check_thread
!------------------------------------------------------------------------------
!> Common eigensolver initialization, called by implementations on first application
!------------------------------------------------------------------------------
subroutine eigsolver_setup(self)
class(oft_eigsolver), intent(inout) :: self !< Solver object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
end subroutine eigsolver_setup
end module oft_solver_base
