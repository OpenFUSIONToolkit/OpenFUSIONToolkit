!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_solver_base.F90
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
MODULE oft_solver_base
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix, oft_cmatrix
IMPLICIT NONE
#include "local.h"
private
!---------------------------------------------------------------------------
! TYPE oft_solver
!---------------------------------------------------------------------------
!> Base class for OFT solvers
!---------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_solver
  LOGICAL :: pm = .FALSE. !< Performance monitor override
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  LOGICAL :: full_residual = .TRUE. !< Output true residual on exit
  INTEGER(i4) :: its = -1 !< Maximum iteration count
  INTEGER(i4) :: cits = 0 !< Number of iteractions to convergence
  INTEGER(i4) :: itplot=10 !< Output frequency for iterative solvers when pm=.TRUE.
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_matrix), POINTER :: A => NULL() !< Matrix to be inverted
  CLASS(oft_solver), POINTER :: pre => NULL() !< Preconditioner
  !> Boundary condition
  PROCEDURE(oft_bc_proto), POINTER, NOPASS :: bc => NULL()
CONTAINS
  !> Solve system
  PROCEDURE(solver_apply), DEFERRED :: apply
  !> Update solver with new settings/operators
  PROCEDURE :: update => solver_update
  !> Setup solver from XML node
  PROCEDURE :: setup_from_xml => solver_setup_xml
  !> Print solver information
  PROCEDURE :: view => solver_view
  !> Check thread safety
  PROCEDURE :: check_thread => solver_check_thread
  !> Clean-up internal storage
  PROCEDURE(solver_delete), DEFERRED :: delete
END TYPE oft_solver
!---------------------------------------------------------------------------
! TYPE oft_solver_ptr
!---------------------------------------------------------------------------
!> Solver container
!---------------------------------------------------------------------------
type, public :: oft_solver_ptr
  CLASS(oft_solver), POINTER :: s => NULL()
end type oft_solver_ptr
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------
  ! SUBROUTINE: solver_apply
  !---------------------------------------------------------------------------
  !> Apply linear solver to compute \f$ A u = g\f$
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !!
  !! @param[in,out] u Guess/Solution field
  !! @param[in,out] g RHS/Residual field
  !---------------------------------------------------------------------------
  SUBROUTINE solver_apply(self,u,g)
  IMPORT oft_solver, oft_vector
  CLASS(oft_solver), INTENT(inout) :: self
  CLASS(oft_vector), INTENT(inout) :: u
  CLASS(oft_vector), INTENT(inout) :: g
  END SUBROUTINE solver_apply
  !---------------------------------------------------------------------------
  ! SUBROUTINE: solver_delete
  !---------------------------------------------------------------------------
  !> Destroy linear solver and deallocate all internal storage
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !---------------------------------------------------------------------------
  SUBROUTINE solver_delete(self)
  IMPORT oft_solver
  CLASS(oft_solver), INTENT(inout) :: self
  END SUBROUTINE solver_delete
  !---------------------------------------------------------------------------
  ! INTERFACE oft_bc_proto
  !---------------------------------------------------------------------------
  !> Abstract boundary condition prototype
  !!
  !! @param[in,out] a Temp doc
  !---------------------------------------------------------------------------
  SUBROUTINE oft_bc_proto(a)
    IMPORT oft_vector
    CLASS(oft_vector), INTENT(inout) :: a
  END SUBROUTINE oft_bc_proto
END INTERFACE
!---------------------------------------------------------------------------
! TYPE oft_csolver
!---------------------------------------------------------------------------
!> Base class for OFT solvers
!---------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_csolver
  LOGICAL :: pm = .FALSE. !< Performance monitor override
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  LOGICAL :: full_residual = .TRUE. !< Output true residual on exit
  INTEGER(i4) :: its = -1 !< Maximum iteration count
  INTEGER(i4) :: cits = 0 !< Number of iteractions to convergence
  INTEGER(i4) :: itplot=10 !< Output frequency for iterative solvers when pm=.TRUE.
  REAL(r8) :: atol = 1.d-14 !< Absolute convergence tolerance \f$ |res| < atol \f$
  REAL(r8) :: rtol = 1.d-14 !< Relative convergence tolerance \f$ |res|/|res_0| < rtol \f$
  CLASS(oft_cmatrix), POINTER :: A => NULL() !< Matrix to be inverted
  CLASS(oft_csolver), POINTER :: pre => NULL() !< Preconditioner
  !> Boundary condition
  PROCEDURE(oft_cbc_proto), POINTER, NOPASS :: bc => NULL()
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
  PROCEDURE(csolver_delete), DEFERRED :: delete
END TYPE oft_csolver
!---------------------------------------------------------------------------
! TYPE oft_solver_ptr
!---------------------------------------------------------------------------
!> Solver container
!---------------------------------------------------------------------------
type, public :: oft_csolver_ptr
  CLASS(oft_csolver), POINTER :: s => NULL()
end type oft_csolver_ptr
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------
  ! SUBROUTINE: csolver_apply
  !---------------------------------------------------------------------------
  !> Apply linear solver to compute \f$ A u = g\f$
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !!
  !! @param[in,out] u Guess/Solution field
  !! @param[in,out] g RHS/Residual field
  !---------------------------------------------------------------------------
  SUBROUTINE csolver_apply(self,u,g)
  IMPORT oft_csolver, oft_cvector
  CLASS(oft_csolver), INTENT(inout) :: self
  CLASS(oft_cvector), INTENT(inout) :: u
  CLASS(oft_cvector), INTENT(inout) :: g
  END SUBROUTINE csolver_apply
  !---------------------------------------------------------------------------
  ! SUBROUTINE: csolver_delete
  !---------------------------------------------------------------------------
  !> Destroy linear solver and deallocate all internal storage
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !---------------------------------------------------------------------------
  SUBROUTINE csolver_delete(self)
  IMPORT oft_csolver
  CLASS(oft_csolver), INTENT(inout) :: self
  END SUBROUTINE csolver_delete
  !---------------------------------------------------------------------------
  ! INTERFACE oft_cbc_proto
  !---------------------------------------------------------------------------
  !> Abstract boundary condition prototype
  !!
  !! @param[in,out] a Temp doc
  !---------------------------------------------------------------------------
  SUBROUTINE oft_cbc_proto(a)
    IMPORT oft_cvector
    CLASS(oft_cvector), INTENT(inout) :: a
  END SUBROUTINE oft_cbc_proto
END INTERFACE
!---------------------------------------------------------------------------
! TYPE oft_eigsolver
!---------------------------------------------------------------------------
!> Base class for OFT eigenvalue solvers
!---------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_eigsolver
  LOGICAL :: pm = .FALSE. !< Performance monitor override
  LOGICAL :: initialized = .FALSE. !< Solver has been constructed
  INTEGER(i4) :: its = -1 !< Maximum iteration count
  INTEGER(i4) :: cits = 0 !< Number of iteractions to convergence
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
  !---------------------------------------------------------------------------
  ! SUBROUTINE: eigsolver_apply
  !---------------------------------------------------------------------------
  !> Apply eigensystem solver to compute eigenvalues/eigenvectors of the system
  !! \f$ A u_i = \lambda_i u_i \f$
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !!
  !! @param[in,out] u Guess/Solution field
  !! @param[in,out] alam Eigenvalue
  !---------------------------------------------------------------------------
  SUBROUTINE eigsolver_apply(self,u,alam)
  IMPORT oft_eigsolver, oft_vector, r8
  CLASS(oft_eigsolver), INTENT(inout) :: self
  CLASS(oft_vector), INTENT(inout) :: u
  REAL(r8), INTENT(inout) :: alam
  END SUBROUTINE eigsolver_apply
  !---------------------------------------------------------------------------
  ! SUBROUTINE: eigsolver_delete
  !---------------------------------------------------------------------------
  !> Destroy eigensystem solver and deallocate all internal storage
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized solvers
  !---------------------------------------------------------------------------
  SUBROUTINE eigsolver_delete(self)
  IMPORT oft_eigsolver
  CLASS(oft_eigsolver), INTENT(inout) :: self
  END SUBROUTINE eigsolver_delete
END INTERFACE
!---------------------------------------------------------------------------
! TYPE oft_orthog
!---------------------------------------------------------------------------
!> Base class for field orthogonalization
!---------------------------------------------------------------------------
type, abstract, public :: oft_orthog
contains
  !> Orthogonalize field
  procedure(orthog_apply), deferred :: apply
  !> Clean-up internal storage
  procedure(orthog_delete), deferred :: delete
end type oft_orthog
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------
  ! SUBROUTINE: orthog_apply
  !---------------------------------------------------------------------------
  !> Orthogonalize field to a specified subspace
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized orthog types
  !!
  !! @param[in,out] u Guess/Solution field
  !---------------------------------------------------------------------------
  SUBROUTINE orthog_apply(self,u)
  IMPORT oft_orthog, oft_vector
  CLASS(oft_orthog), INTENT(inout) :: self
  CLASS(oft_vector), INTENT(inout) :: u
  END SUBROUTINE orthog_apply
  !---------------------------------------------------------------------------
  ! SUBROUTINE: orthog_delete
  !---------------------------------------------------------------------------
  !> Destroy orthogonalization structure and deallocate all internal storage
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized orthog types
  !---------------------------------------------------------------------------
  SUBROUTINE orthog_delete(self)
  IMPORT oft_orthog
  CLASS(oft_orthog), INTENT(inout) :: self
  END SUBROUTINE orthog_delete
END INTERFACE
!---------------------------------------------------------------------------
! TYPE oft_corthog
!---------------------------------------------------------------------------
!> Base class for field orthogonalization
!---------------------------------------------------------------------------
type, abstract, public :: oft_corthog
contains
  !> Orthogonalize field
  procedure(corthog_apply), deferred :: apply
  !> Clean-up internal storage
  procedure(corthog_delete), deferred :: delete
end type oft_corthog
ABSTRACT INTERFACE
  !---------------------------------------------------------------------------
  ! SUBROUTINE: corthog_apply
  !---------------------------------------------------------------------------
  !> Orthogonalize field to a specified subspace
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized orthog types
  !!
  !! @param[in,out] u Guess/Solution field
  !---------------------------------------------------------------------------
  SUBROUTINE corthog_apply(self,u)
  IMPORT oft_corthog, oft_cvector
  CLASS(oft_corthog), INTENT(inout) :: self
  CLASS(oft_cvector), INTENT(inout) :: u
  END SUBROUTINE corthog_apply
  !---------------------------------------------------------------------------
  ! SUBROUTINE: corthog_delete
  !---------------------------------------------------------------------------
  !> Destroy orthogonalization structure and deallocate all internal storage
  !!
  !! @note This subroutine is a dummy routine used to specify the interface
  !! of the member function and catch errors in uninitialized orthog types
  !---------------------------------------------------------------------------
  SUBROUTINE corthog_delete(self)
  IMPORT oft_corthog
  CLASS(oft_corthog), INTENT(inout) :: self
  END SUBROUTINE corthog_delete
END INTERFACE
!---Make classes and prototypes public
PUBLIC solver_setup, csolver_setup, eigsolver_setup, oft_bc_proto, oft_cbc_proto
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: solver_setup
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
SUBROUTINE solver_setup(self)
CLASS(oft_solver), INTENT(inout) :: self
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE solver_setup
!---------------------------------------------------------------------------
! SUBROUTINE: solver_update
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!!
!! @param[in] new_pattern Update matrix pattern (optional)
!---------------------------------------------------------------------------
recursive SUBROUTINE solver_update(self,new_pattern)
CLASS(oft_solver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
IF(ASSOCIATED(self%pre))CALL self%pre%update(new_pattern)
END SUBROUTINE solver_update
!---------------------------------------------------------------------------
! SUBROUTINE: solver_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
SUBROUTINE solver_setup_xml(self,solver_node,level)
class(oft_solver), intent(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
end SUBROUTINE solver_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: solver_view
!---------------------------------------------------------------------------
!> Print solver configuration
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!---------------------------------------------------------------------------
SUBROUTINE solver_view(self)
class(oft_solver), intent(inout) :: self
end SUBROUTINE solver_view
!---------------------------------------------------------------------------
! FUNCTION: solver_check_thread
!---------------------------------------------------------------------------
!> Check thread safety
!---------------------------------------------------------------------------
recursive function solver_check_thread(self) result(thread_safe)
class(oft_solver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
IF(ASSOCIATED(self%pre))thread_safe=(thread_safe.AND.self%pre%check_thread())
end function solver_check_thread
!---------------------------------------------------------------------------
! SUBROUTINE: csolver_setup
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
SUBROUTINE csolver_setup(self)
CLASS(oft_csolver), INTENT(inout) :: self
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE csolver_setup
!---------------------------------------------------------------------------
! SUBROUTINE: csolver_update
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!!
!! @param[in] new_pattern Update matrix pattern (optional)
!---------------------------------------------------------------------------
recursive SUBROUTINE csolver_update(self,new_pattern)
CLASS(oft_csolver), intent(inout) :: self
LOGICAL, optional, intent(in) :: new_pattern
IF(ASSOCIATED(self%pre))CALL self%pre%update(new_pattern)
END SUBROUTINE csolver_update
!---------------------------------------------------------------------------
! SUBROUTINE: csolver_setup_xml
!---------------------------------------------------------------------------
!> Setup solver from XML definition
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!!
!! @param[in] solver_node XML node containing solver definition
!! @param[in] level Level in MG hierarchy (optional)
!---------------------------------------------------------------------------
SUBROUTINE csolver_setup_xml(self,solver_node,level)
class(oft_csolver), intent(inout) :: self
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
end SUBROUTINE csolver_setup_xml
!---------------------------------------------------------------------------
! SUBROUTINE: csolver_view
!---------------------------------------------------------------------------
!> Print solver configuration
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized solvers
!---------------------------------------------------------------------------
SUBROUTINE csolver_view(self)
class(oft_csolver), intent(inout) :: self
end SUBROUTINE csolver_view
!---------------------------------------------------------------------------
! FUNCTION: csolver_check_thread
!---------------------------------------------------------------------------
!> Check thread safety
!---------------------------------------------------------------------------
recursive function csolver_check_thread(self) result(thread_safe)
class(oft_csolver), intent(inout) :: self
logical :: thread_safe
thread_safe=.TRUE.
IF(ASSOCIATED(self%pre))thread_safe=(thread_safe.AND.self%pre%check_thread())
end function csolver_check_thread
!---------------------------------------------------------------------------
! SUBROUTINE: eigsolver_setup
!---------------------------------------------------------------------------
!> Update solver after changing settings/operators
!---------------------------------------------------------------------------
subroutine eigsolver_setup(self)
class(oft_eigsolver), intent(inout) :: self
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%pre))THEN
  IF(.NOT.ASSOCIATED(self%pre%A))self%pre%A=>self%A
END IF
self%initialized=.TRUE.
DEBUG_STACK_POP
end subroutine eigsolver_setup
end module oft_solver_base
