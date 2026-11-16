!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_tracing.F90
!
!> Regression tests for 3D field line tracing
!!
!! Traces field lines of an analytic helical field, B = (-y, x, Bz), through a cube
!! mesh (distributed over processors for parallel runs) and compares the results to
!! the exact solution. Covers tracer types 1 (LSODE) and 4 (Runge-Kutta with steps limited
!! at cell boundaries), including transfers of the integrator state between domains.
!!
!! Three modes are available:
!! - Single field line (`tracing_line`): the end point is compared to the exact solution
!! - Poincare section (`tracing_poincare`): several field lines are traced in parallel
!!   and the crossings of the plane x=0 are compared to the exact helices
!! - Field with jumps between cells (`field_jump>0`): the helical field is scaled by a factor
!!   that varies from cell to cell, so the field jumps at every cell face as for finite element
!!   fields while its normal component keeps its sign. The field line still follows the exact
!!   helix, which is checked here, while the time along it (traced as a fourth equation) is
!!   compared by the driver to a reference trace with a tight tolerance
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
MODULE test_tracing_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
REAL(r8) :: swirl_bz = 0.3d0 !< Axial field (pitch) of the helical field
REAL(r8) :: swirl_jump = 0.d0 !< Amplitude of the variation of the field strength between cells
!> Analytic helical field
TYPE, EXTENDS(fem_interp) :: swirl_interp
CONTAINS
  PROCEDURE :: interp => swirl_interp_apply
END TYPE swirl_interp
CONTAINS
!> Evaluate field at a logical position in a cell
SUBROUTINE swirl_interp_apply(self,cell,f,gop,val)
CLASS(swirl_interp), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3),pc(3)
pt=self%mesh%log2phys(cell,f)
val(1)=-pt(2)
val(2)=pt(1)
val(3)=swirl_bz
IF(swirl_jump>0.d0)THEN
  !---Term set by the cell center, so it is independent of the cell numbering and decomposition
  IF(self%mesh%type==3)THEN
    pc=self%mesh%log2phys(cell,[0.5d0,0.5d0,0.5d0,0.d0])
  ELSE
    pc=self%mesh%log2phys(cell,[0.25d0,0.25d0,0.25d0,0.25d0])
  END IF
  val(1:3)=val(1:3)*(1.d0+swirl_jump*SIN(37.d0*pc(1)+53.d0*pc(2)+71.d0*pc(3)))
END IF
END SUBROUTINE swirl_interp_apply
!> Field line ODE with time as a fourth equation
SUBROUTINE time_ydot(t,y,B,n,ydot)
REAL(r8), INTENT(in) :: t
REAL(r8), INTENT(in) :: y(n)
REAL(r8), INTENT(in) :: B(3)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(out) :: ydot(n)
ydot(1:3)=B
ydot(4)=1.d0
END SUBROUTINE time_ydot
END MODULE test_tracing_helpers
!
PROGRAM test_tracing
USE oft_base
USE oft_io, ONLY: hdf5_create_file
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE tracing, ONLY: oft_tracer, create_tracer, tracing_line, tracing_poincare
USE test_tracing_helpers, ONLY: swirl_interp, swirl_bz, swirl_jump, time_ydot
IMPLICIT NONE
INTEGER(i4) :: iounit,ierr
REAL(r8) :: pt(3),phase,err
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(swirl_interp), TARGET :: Bfield
CLASS(oft_tracer), POINTER :: tracer
INTEGER(i4) :: tracer_type=1
LOGICAL :: test_poincare=.FALSE.
REAL(r8) :: field_jump=0.d0
REAL(r8) :: tol=1.d-10
NAMELIST/test_tracing_options/tracer_type,test_poincare,field_jump,tol
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=iounit,FILE=oft_env%ifile)
READ(iounit,test_tracing_options,IOSTAT=ierr)
CLOSE(iounit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','test_tracing',__FILE__)
!---Setup tracer
Bfield%mesh=>mg_mesh%mesh
Bfield%dim=3
CALL create_tracer(tracer,tracer_type)
tracer%tol=tol
tracer%maxsteps=100000
tracer%maxtrans=1000
tracer%B=>Bfield
tracer%ml_mesh=>mg_mesh
IF(test_poincare)THEN
  CALL run_poincare_test
  CALL tracer%delete()
  CALL oft_finalize
END IF
IF(field_jump>0.d0)THEN
  CALL run_jump_test
  CALL tracer%delete()
  CALL oft_finalize
END IF
!---Trace from (0.5,0,-0.9): circle of radius 0.5 advancing in z until exiting at z=1
pt=[0.5d0,0.d0,-0.9d0]
IF(oft_env%head_proc)CALL hdf5_create_file('trace_test')
CALL oft_mpi_barrier(ierr)
CALL tracing_line(tracer,pt,'trace_test')
!---Check result against exact solution x=0.5*cos(t), y=0.5*sin(t), z=z0+0.3*t
IF(tracer%status/=-1)CALL oft_abort('Tracer did not exit mesh as expected','test_tracing',__FILE__)
IF(tracer%nsteps<10)CALL oft_abort('Too few tracer steps','test_tracing',__FILE__)
IF(oft_env%nprocs>1.AND.tracer%ntrans<1)CALL oft_abort('No domain transfers occurred','test_tracing',__FILE__)
IF(tracer%y(3)<0.9d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
phase=(tracer%y(3)-pt(3))/0.3d0
err=SQRT((tracer%y(1)-0.5d0*COS(phase))**2+(tracer%y(2)-0.5d0*SIN(phase))**2)
IF(oft_env%head_proc)WRITE(*,'(A,ES11.3,2I8)')'  Final error, steps, transfers = ',err,tracer%nsteps,tracer%ntrans
IF(err>1.d-6)CALL oft_abort('Tracer end point differs from exact solution','test_tracing',__FILE__)
CALL tracer%delete()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---Poincare section of several field lines, compared to the exact helices
SUBROUTINE run_poincare_test
INTEGER(i4), PARAMETER :: nlines=3
INTEGER(i4), PARAMETER :: ncross=6
REAL(r8), PARAMETER :: z0=-0.9d0
INTEGER(i4) :: i,k,iline,ncount(nlines),nfound(ncross,nlines),io_stat
REAL(r8) :: pts(3,nlines),radius(nlines),pt_cross(3),z_exact,err,err_max
LOGICAL :: have_omp
!---Poincare tracing requires OpenMP
have_omp=.FALSE.
!$ have_omp=.TRUE.
IF(.NOT.have_omp)THEN
  IF(oft_env%head_proc)WRITE(*,*)'SKIP TEST'
  CALL oft_finalize
END IF
!---Lower pitch so each line makes several turns before exiting at z=1
swirl_bz=0.05d0
radius=[0.25d0,0.5d0,0.75d0]
DO i=1,nlines
  pts(:,i)=[radius(i),0.d0,z0]
END DO
CALL tracing_poincare(tracer,pts,nlines,'poincare_test.dat')
CALL oft_mpi_barrier(ierr)
!---Exact crossings of x=0 (x increasing): theta=3*pi/2+2*pi*k -> (0,-r,z0+Bz*theta)
!---Lines exit at z=1 (theta=38), so k=0,...,5 are crossed and the next at z=1.22 is not
IF(oft_env%head_proc)THEN
  ncount=0
  nfound=0
  err_max=0.d0
  OPEN(NEWUNIT=iounit,FILE='poincare_test.dat',STATUS='OLD')
  DO
    READ(iounit,*,IOSTAT=io_stat)iline,pt_cross
    IF(io_stat/=0)EXIT
    IF(iline<1.OR.iline>nlines)CALL oft_abort('Invalid field line index','test_tracing',__FILE__)
    ncount(iline)=ncount(iline)+1
    k=NINT(((pt_cross(3)-z0)/swirl_bz-1.5d0*pi)/(2.d0*pi))
    IF(k<0.OR.k>=ncross)CALL oft_abort('Unexpected crossing','test_tracing',__FILE__)
    nfound(k+1,iline)=nfound(k+1,iline)+1
    z_exact=z0+swirl_bz*(1.5d0*pi+2.d0*pi*k)
    err=SQRT(pt_cross(1)**2+(pt_cross(2)+radius(iline))**2+(pt_cross(3)-z_exact)**2)
    err_max=MAX(err_max,err)
  END DO
  CLOSE(iounit)
  WRITE(*,'(A,ES11.3,3I4)')'  Max crossing error, crossings per line = ',err_max,ncount
  IF(ANY(nfound/=1))CALL oft_abort('Wrong number of crossings','test_tracing',__FILE__)
  IF(err_max>1.d-3)CALL oft_abort('Crossings differ from exact solution','test_tracing',__FILE__)
END IF
END SUBROUTINE run_poincare_test
!---Single field line with a field that jumps at every cell face
SUBROUTINE run_jump_test
REAL(r8) :: pt4(4)
swirl_jump=field_jump
tracer%neq=4
tracer%ydot=>time_ydot
pt4=[0.5d0,0.d0,-0.9d0,0.d0]
IF(oft_env%head_proc)CALL hdf5_create_file('trace_test')
CALL oft_mpi_barrier(ierr)
CALL tracing_line(tracer,pt4,'trace_test')
IF(tracer%status/=-1)CALL oft_abort('Tracer did not exit mesh as expected','test_tracing',__FILE__)
IF(oft_env%nprocs>1.AND.tracer%ntrans<1)CALL oft_abort('No domain transfers occurred','test_tracing',__FILE__)
IF(tracer%y(3)<0.9d0)CALL oft_abort('Tracer exited early','test_tracing',__FILE__)
!---The jumps only change the field strength, so the field line is still the exact helix
phase=(tracer%y(3)-pt4(3))/swirl_bz
err=SQRT((tracer%y(1)-0.5d0*COS(phase))**2+(tracer%y(2)-0.5d0*SIN(phase))**2)
IF(oft_env%head_proc)WRITE(*,'(A,ES11.3,2I8)')'  Final error, steps, transfers = ',err,tracer%nsteps,tracer%ntrans
IF(err>MAX(1.d2*tol,1.d-8))CALL oft_abort('Tracer end point is not on the exact helix','test_tracing',__FILE__)
END SUBROUTINE run_jump_test
END PROGRAM test_tracing
