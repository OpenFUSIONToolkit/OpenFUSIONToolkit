!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_cvode.F90
!
!> Minimal interface to CVODE (SUNDIALS) for non-stiff ODE integration
!!
!! Provides LSODE-style usage (right-hand-side callback `f(neq,t,y,ydot)`, scalar
!! tolerances, one internal step per call) for the Adams-Moulton method with
!! functional iteration. Each \ref oft_cvode::oft_cvode_ode "oft_cvode_ode" object owns an
!! independent integrator, so no state is shared between objects or threads.
!!
!! This interface requires OFT to be built with CVODE (`OFT_CVODE_ROOT`). Without
!! it @ref oft_cvode::oft_cvode_enabled is `.FALSE.` and all methods return an error.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_base
!--------------------------------------------------------------------------------
MODULE oft_cvode
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_double, c_ptr, c_null_ptr, c_funptr, &
  c_funloc, c_loc, c_f_pointer, c_associated
USE oft_local
IMPLICIT NONE
#define MOD_NAME oft_cvode
#include "local.h"
PRIVATE
!> Error code returned when OFT was built without CVODE
INTEGER(i4), PUBLIC, PARAMETER :: OFT_CVODE_UNAVAILABLE = -100
!> Flag indicating CVODE is available in this build of OFT
#ifdef HAVE_CVODE
LOGICAL, PUBLIC, PARAMETER :: oft_cvode_enabled = .TRUE.
#else
LOGICAL, PUBLIC, PARAMETER :: oft_cvode_enabled = .FALSE.
#endif
!------------------------------------------------------------------------------
!> Abstract interface for the ODE right-hand-side function
!------------------------------------------------------------------------------
ABSTRACT INTERFACE
  !> Evaluate `ydot = f(t,y)`
  SUBROUTINE cvode_rhs(neq,t,y,ydot)
    IMPORT i4, r8
    INTEGER(i4), INTENT(in) :: neq !< Number of equations
    REAL(r8), INTENT(in) :: t !< Current time
    REAL(r8), INTENT(in) :: y(neq) !< Current solution
    REAL(r8), INTENT(out) :: ydot(neq) !< Time derivative of solution
  END SUBROUTINE cvode_rhs
END INTERFACE
PUBLIC :: cvode_rhs
!> Callback storage passed through CVODE (non-polymorphic so it can be referenced from C)
TYPE :: cvode_rhs_holder
  PROCEDURE(cvode_rhs), POINTER, NOPASS :: f => NULL()
END TYPE cvode_rhs_holder
!------------------------------------------------------------------------------
!> Non-stiff ODE integrator using CVODE (Adams method, functional iteration)
!!
!! Usage:
!! ```
!! CALL ode%setup(neq,rhs,rtol,atol,t0,y0,ierr)
!! DO
!!   CALL ode%step(y,t,tout,ierr)  ! One internal step, direction/scale set by `tout`
!!   ...
!! END DO
!! CALL ode%delete
!! ```
!! Step return codes are `0` for success, `1` if `tstop` was reached, and `<0` for failure.
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_cvode_ode
  INTEGER(i4) :: neq = 0 !< Number of equations
  REAL(r8) :: h = 0.d0 !< Size of last successful step
  TYPE(c_ptr), PRIVATE :: handle = c_null_ptr !< Handle to C integrator object
  TYPE(cvode_rhs_holder), POINTER, PRIVATE :: holder => NULL() !< Right-hand-side callback
CONTAINS
  !> Create integrator and set initial condition
  PROCEDURE :: setup => cvode_setup
  !> Restart integration from a new initial condition
  PROCEDURE :: reinit => cvode_reinit
  !> Advance by one internal step
  PROCEDURE :: step => cvode_step
  !> Check if integrator is active
  PROCEDURE :: active => cvode_active
  !> Free integrator
  PROCEDURE :: delete => cvode_delete
END TYPE oft_cvode_ode
!> @cond
#ifdef HAVE_CVODE
INTERFACE
  FUNCTION oft_cvode_create_c(neq,t0,y0,rtol,atol,f,udata,handle) RESULT(ierr) BIND(C)
    IMPORT c_int, c_double, c_funptr, c_ptr
    INTEGER(c_int), VALUE, INTENT(in) :: neq
    REAL(c_double), VALUE, INTENT(in) :: t0
    REAL(c_double), INTENT(in) :: y0(*)
    REAL(c_double), VALUE, INTENT(in) :: rtol
    REAL(c_double), VALUE, INTENT(in) :: atol
    TYPE(c_funptr), VALUE, INTENT(in) :: f
    TYPE(c_ptr), VALUE, INTENT(in) :: udata
    TYPE(c_ptr), INTENT(out) :: handle
    INTEGER(c_int) :: ierr
  END FUNCTION oft_cvode_create_c
  FUNCTION oft_cvode_reinit_c(handle,t0,y0,h0) RESULT(ierr) BIND(C)
    IMPORT c_int, c_double, c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: handle
    REAL(c_double), VALUE, INTENT(in) :: t0
    REAL(c_double), INTENT(in) :: y0(*)
    REAL(c_double), VALUE, INTENT(in) :: h0
    INTEGER(c_int) :: ierr
  END FUNCTION oft_cvode_reinit_c
  FUNCTION oft_cvode_step_c(handle,tout,use_tstop,tstop,y,t,hlast) RESULT(ierr) BIND(C)
    IMPORT c_int, c_double, c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: handle
    REAL(c_double), VALUE, INTENT(in) :: tout
    INTEGER(c_int), VALUE, INTENT(in) :: use_tstop
    REAL(c_double), VALUE, INTENT(in) :: tstop
    REAL(c_double), INTENT(inout) :: y(*)
    REAL(c_double), INTENT(inout) :: t
    REAL(c_double), INTENT(inout) :: hlast
    INTEGER(c_int) :: ierr
  END FUNCTION oft_cvode_step_c
  SUBROUTINE oft_cvode_free_c(handle) BIND(C)
    IMPORT c_ptr
    TYPE(c_ptr), INTENT(inout) :: handle
  END SUBROUTINE oft_cvode_free_c
END INTERFACE
#endif
!> @endcond
CONTAINS
#ifdef HAVE_CVODE
!------------------------------------------------------------------------------
!> Right-hand-side adapter called from C, forwards to the user function
!------------------------------------------------------------------------------
FUNCTION cvode_rhs_adapter(neq,t,y,ydot,udata) RESULT(ierr) BIND(C,NAME='oft_cvode_rhs_adapter_f')
INTEGER(c_int), VALUE, INTENT(in) :: neq !< Number of equations
REAL(c_double), VALUE, INTENT(in) :: t !< Current time
REAL(c_double), INTENT(in) :: y(neq) !< Current solution
REAL(c_double), INTENT(out) :: ydot(neq) !< Time derivative of solution
TYPE(c_ptr), VALUE, INTENT(in) :: udata !< Pointer to callback storage
INTEGER(c_int) :: ierr !< Status flag (always 0)
TYPE(cvode_rhs_holder), POINTER :: holder
CALL c_f_pointer(udata,holder)
CALL holder%f(neq,t,y,ydot)
ierr=0
END FUNCTION cvode_rhs_adapter
#endif
!------------------------------------------------------------------------------
!> Create integrator and set initial condition
!------------------------------------------------------------------------------
SUBROUTINE cvode_setup(self,neq,rhs,rtol,atol,t0,y0,ierr)
CLASS(oft_cvode_ode), INTENT(inout) :: self !< Integrator object
INTEGER(i4), INTENT(in) :: neq !< Number of equations
PROCEDURE(cvode_rhs) :: rhs !< Right-hand-side function
REAL(r8), INTENT(in) :: rtol !< Relative tolerance
REAL(r8), INTENT(in) :: atol !< Absolute tolerance
REAL(r8), INTENT(in) :: t0 !< Initial time
REAL(r8), INTENT(in) :: y0(neq) !< Initial solution
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
#ifdef HAVE_CVODE
IF(self%active())CALL self%delete()
ALLOCATE(self%holder)
self%holder%f=>rhs
self%neq=neq
self%h=0.d0
ierr=oft_cvode_create_c(neq,t0,y0,rtol,atol,c_funloc(cvode_rhs_adapter),c_loc(self%holder),self%handle)
IF(ierr/=0)CALL self%delete()
#else
ierr=OFT_CVODE_UNAVAILABLE
#endif
END SUBROUTINE cvode_setup
!------------------------------------------------------------------------------
!> Restart integration from a new initial condition
!!
!! Resets the method order and step size history. An initial step size can
!! optionally be supplied (eg. the last step before a restart).
!------------------------------------------------------------------------------
SUBROUTINE cvode_reinit(self,t0,y0,ierr,h0)
CLASS(oft_cvode_ode), INTENT(inout) :: self !< Integrator object
REAL(r8), INTENT(in) :: t0 !< Initial time
REAL(r8), INTENT(in) :: y0(:) !< Initial solution
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
REAL(r8), OPTIONAL, INTENT(in) :: h0 !< Initial step size (optional, default is estimated)
#ifdef HAVE_CVODE
REAL(r8) :: h0_loc
h0_loc=0.d0
IF(PRESENT(h0))h0_loc=h0
ierr=oft_cvode_reinit_c(self%handle,t0,y0,h0_loc)
#else
ierr=OFT_CVODE_UNAVAILABLE
#endif
END SUBROUTINE cvode_reinit
!------------------------------------------------------------------------------
!> Advance by one internal step
!!
!! As with LSODE `ITASK=2`, only the sign and scale of `tout` matter. If `tstop`
!! is given the integrator will not step past it (LSODE `ITASK=5`) and returns
!! `ierr=1` when it is reached.
!------------------------------------------------------------------------------
SUBROUTINE cvode_step(self,y,t,tout,ierr,tstop)
CLASS(oft_cvode_ode), INTENT(inout) :: self !< Integrator object
REAL(r8), INTENT(inout) :: y(:) !< Solution at time `t`
REAL(r8), INTENT(inout) :: t !< Current time (updated to the time reached)
REAL(r8), INTENT(in) :: tout !< Target time (sets direction and time scale)
INTEGER(i4), INTENT(out) :: ierr !< Return code (0: ok, 1: `tstop` reached, <0: failure)
REAL(r8), OPTIONAL, INTENT(in) :: tstop !< Time that must not be passed (optional)
#ifdef HAVE_CVODE
INTEGER(c_int) :: use_tstop
REAL(r8) :: tstop_loc
use_tstop=0
tstop_loc=0.d0
IF(PRESENT(tstop))THEN
  use_tstop=1
  tstop_loc=tstop
END IF
ierr=oft_cvode_step_c(self%handle,tout,use_tstop,tstop_loc,y,t,self%h)
#else
ierr=OFT_CVODE_UNAVAILABLE
#endif
END SUBROUTINE cvode_step
!------------------------------------------------------------------------------
!> Check if integrator has been created and not deleted
!------------------------------------------------------------------------------
FUNCTION cvode_active(self) RESULT(is_active)
CLASS(oft_cvode_ode), INTENT(in) :: self !< Integrator object
LOGICAL :: is_active !< Integrator is active?
is_active=c_associated(self%handle)
END FUNCTION cvode_active
!------------------------------------------------------------------------------
!> Free integrator
!------------------------------------------------------------------------------
SUBROUTINE cvode_delete(self)
CLASS(oft_cvode_ode), INTENT(inout) :: self !< Integrator object
#ifdef HAVE_CVODE
IF(c_associated(self%handle))CALL oft_cvode_free_c(self%handle)
#endif
self%handle=c_null_ptr
IF(ASSOCIATED(self%holder))DEALLOCATE(self%holder)
self%neq=0
self%h=0.d0
END SUBROUTINE cvode_delete
END MODULE oft_cvode
