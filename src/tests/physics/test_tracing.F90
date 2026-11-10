!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_tracing.F90
!
!> Regression tests for 3D field line tracing
!!
!! Traces a field line of an analytic helical field, B = (-y, x, 0.3), through a cube
!! mesh (distributed over processors for parallel runs) and compares the end point to
!! the exact solution. Covers tracer types 1 (CVODE if available, otherwise LSODE) and
!! 4 (LSODE), including restart of the integrator after transfers between domains.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
MODULE test_tracing_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
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
REAL(r8) :: pt(3)
pt=self%mesh%log2phys(cell,f)
val(1)=-pt(2)
val(2)=pt(1)
val(3)=0.3d0
END SUBROUTINE swirl_interp_apply
END MODULE test_tracing_helpers
!
PROGRAM test_tracing
USE oft_base
USE oft_io, ONLY: hdf5_create_file
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE tracing, ONLY: oft_tracer, create_tracer, tracing_line
USE test_tracing_helpers, ONLY: swirl_interp
IMPLICIT NONE
INTEGER(i4) :: iounit,ierr
REAL(r8) :: pt(3),phase,err
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(swirl_interp), TARGET :: Bfield
CLASS(oft_tracer), POINTER :: tracer
INTEGER(i4) :: tracer_type=1
NAMELIST/test_tracing_options/tracer_type
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
tracer%tol=1.d-10
tracer%maxsteps=100000
tracer%maxtrans=1000
tracer%B=>Bfield
tracer%ml_mesh=>mg_mesh
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
END PROGRAM test_tracing
