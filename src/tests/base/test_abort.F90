!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_abort.F90
!
!> Regression tests for @ref oft_base::oft_abort "oft_abort"
!!
!! Each test ends the run through oft_abort, so test_abort.py checks the exit
!! code and the number of error blocks printed instead of a result.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
MODULE test_abort_helpers
USE oft_base
IMPLICIT NONE
CONTAINS
! Abort callback that re-enters oft_abort
SUBROUTINE reentrant_abort_cb() BIND(C)
CALL oft_abort('Abort called from abort callback','reentrant_abort_cb',__FILE__)
END SUBROUTINE reentrant_abort_cb
END MODULE test_abort_helpers

PROGRAM test_abort
USE oft_base
USE test_abort_helpers
IMPLICIT NONE
INTEGER(i4) :: i,iounit,ierr
REAL(r8) :: x
INTEGER(i4) :: test_id=1
NAMELIST/test_abort_options/test_id
!------------------------------------------------------------------------------
! Initialize environment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_abort_options,IOSTAT=ierr)
CLOSE(iounit)
!---Run test
SELECT CASE(test_id)
  CASE(1) ! Abort from every thread of a parallel loop
    x=0.d0
    !$omp parallel do reduction(+:x)
    DO i=1,64
      x=x+REAL(i,8)
      CALL oft_abort('Parallel abort test','test_abort',__FILE__)
    END DO
  CASE(2) ! Abort outside a parallel region
    CALL oft_abort('Serial abort test','test_abort',__FILE__)
  CASE(3) ! Abort callback that calls oft_abort again
    oft_abort_cb=>reentrant_abort_cb
    CALL oft_abort('Re-entrant abort test','test_abort',__FILE__)
  CASE DEFAULT
    CALL oft_abort('Unknown test','test_abort',__FILE__)
END SELECT
!---Should never be reached
WRITE(*,*)'Run continued after oft_abort'
CALL oft_finalize
END PROGRAM test_abort
