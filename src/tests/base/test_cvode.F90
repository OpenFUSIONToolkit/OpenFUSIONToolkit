!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_cvode.F90
!
!> Regression tests for the CVODE ODE integrator interface
!!
!! - Test 1: Exponential decay stepped to a stop time
!! - Test 2: Harmonic oscillator over one period with restart
!! - Test 3: Free and re-create integrator, with step size carried over a restart
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
PROGRAM test_cvode
USE oft_base
USE oft_cvode, ONLY: oft_cvode_ode, oft_cvode_enabled
IMPLICIT NONE
INTEGER(i4) :: iounit,ierr,nsteps,ierr_step
REAL(r8) :: t,y1(1),y2(2),h_last
TYPE(oft_cvode_ode) :: ode
INTEGER(i4) :: test_id=1
NAMELIST/test_cvode_options/test_id
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_cvode_options,IOSTAT=ierr)
CLOSE(iounit)
IF(.NOT.oft_cvode_enabled)CALL oft_abort('OFT built without CVODE','test_cvode',__FILE__)

SELECT CASE(test_id)
CASE(1)
  !---y' = -y, y(0) = 1: step until t = 1 (must stop exactly at tstop)
  WRITE(*,*)'Test 1'
  t=0.d0
  y1=1.d0
  CALL ode%setup(1_i4,decay_rhs,1.d-10,1.d-10,t,y1,ierr)
  IF(ierr/=0)CALL oft_abort('Setup failed','test_cvode',__FILE__)
  nsteps=0
  DO
    CALL ode%step(y1,t,1.d0,ierr_step,tstop=1.d0)
    IF(ierr_step<0)CALL oft_abort('Step failed','test_cvode',__FILE__)
    nsteps=nsteps+1
    IF(ierr_step==1)EXIT
    IF(nsteps>1000)CALL oft_abort('Too many steps','test_cvode',__FILE__)
  END DO
  IF(ABS(t-1.d0)>1.d-14)CALL oft_abort('Did not stop at tstop','test_cvode',__FILE__)
  IF(ABS(y1(1)-EXP(-1.d0))>1.d-8)CALL oft_abort('Incorrect solution','test_cvode',__FILE__)
  IF(nsteps<2)CALL oft_abort('Expected multiple internal steps','test_cvode',__FILE__)
  IF(ode%h<=0.d0)CALL oft_abort('Last step size not reported','test_cvode',__FILE__)
  CALL ode%delete()
  IF(ode%active())CALL oft_abort('Integrator still active after delete','test_cvode',__FILE__)
CASE(2)
  !---y1' = y2, y2' = -y1 with (1,0): one full period, then a restart
  WRITE(*,*)'Test 2'
  t=0.d0
  y2=[1.d0,0.d0]
  CALL ode%setup(2_i4,oscillator_rhs,1.d-10,1.d-10,t,y2,ierr)
  IF(ierr/=0)CALL oft_abort('Setup failed','test_cvode',__FILE__)
  DO
    CALL ode%step(y2,t,2.d0*pi,ierr_step,tstop=2.d0*pi)
    IF(ierr_step<0)CALL oft_abort('Step failed','test_cvode',__FILE__)
    IF(ierr_step==1)EXIT
  END DO
  IF(ABS(y2(1)-1.d0)>1.d-7.OR.ABS(y2(2))>1.d-7)CALL oft_abort('Incorrect solution after one period','test_cvode',__FILE__)
  !---Restart from a new state and integrate a quarter period
  t=0.d0
  y2=[0.d0,1.d0]
  CALL ode%reinit(t,y2,ierr)
  IF(ierr/=0)CALL oft_abort('Reinit failed','test_cvode',__FILE__)
  DO
    CALL ode%step(y2,t,pi/2.d0,ierr_step,tstop=pi/2.d0)
    IF(ierr_step<0)CALL oft_abort('Step failed after reinit','test_cvode',__FILE__)
    IF(ierr_step==1)EXIT
  END DO
  IF(ABS(y2(1)-1.d0)>1.d-7.OR.ABS(y2(2))>1.d-7)CALL oft_abort('Incorrect solution after reinit','test_cvode',__FILE__)
  CALL ode%delete()
CASE(3)
  !---Repeat Test 1 in "one step" mode with a restart that reuses the last step size
  WRITE(*,*)'Test 3'
  t=0.d0
  y1=1.d0
  CALL ode%setup(1_i4,decay_rhs,1.d-10,1.d-10,t,y1,ierr)
  IF(ierr/=0)CALL oft_abort('Setup failed','test_cvode',__FILE__)
  DO
    CALL ode%step(y1,t,t+1.d0,ierr_step)
    IF(ierr_step<0)CALL oft_abort('Step failed','test_cvode',__FILE__)
    IF(t>=0.5d0)EXIT
  END DO
  !---Restart at the current state (as done after a domain transfer)
  h_last=ode%h
  IF(h_last<=0.d0)CALL oft_abort('Last step size not reported','test_cvode',__FILE__)
  CALL ode%reinit(t,y1,ierr,h0=h_last)
  IF(ierr/=0)CALL oft_abort('Reinit failed','test_cvode',__FILE__)
  DO
    CALL ode%step(y1,t,t+1.d0,ierr_step)
    IF(ierr_step<0)CALL oft_abort('Step failed after reinit','test_cvode',__FILE__)
    IF(t>=1.d0)EXIT
  END DO
  IF(ABS(y1(1)-EXP(-t))>1.d-8)CALL oft_abort('Incorrect solution after restart','test_cvode',__FILE__)
  CALL ode%delete()
END SELECT

!---Finalize enviroment
CALL oft_finalize
CONTAINS
!> Exponential decay
SUBROUTINE decay_rhs(neq,t,y,ydot)
INTEGER(i4), INTENT(in) :: neq
REAL(r8), INTENT(in) :: t,y(neq)
REAL(r8), INTENT(out) :: ydot(neq)
ydot(1)=-y(1)
END SUBROUTINE decay_rhs
!> Harmonic oscillator
SUBROUTINE oscillator_rhs(neq,t,y,ydot)
INTEGER(i4), INTENT(in) :: neq
REAL(r8), INTENT(in) :: t,y(neq)
REAL(r8), INTENT(out) :: ydot(neq)
ydot(1)=y(2)
ydot(2)=-y(1)
END SUBROUTINE oscillator_rhs
END PROGRAM test_cvode
