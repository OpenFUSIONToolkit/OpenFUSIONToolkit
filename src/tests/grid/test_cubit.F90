!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_cubit.F90
!
!> Regression tests for the CUBIT mesh interface. Test perform refinements
!! and compare the resulting mesh volume and surface area with results from
!! refinement within CUBIT.
!!
!! The current test cases are:
!! - Unit sphere
!! - Unit sphere with cuts
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_cubit
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh
USE oft_mesh_native, ONLY: mesh_native_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct, multigrid_construct_surf
IMPLICIT NONE
#include "local.h"
INTEGER(i4) :: io_unit
INTEGER(i4) :: ierr
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
LOGICAL :: test_surf = .FALSE.
NAMELIST/cubit_test_options/test_surf
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,cubit_test_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check inputs
IF(ierr<0)CALL oft_abort('No "cubit_test_options" found in input file.', &
  'test_cubit',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "cubit_test_options" in input file.', &
  'test_cubit',__FILE__)
IF(test_surf)THEN
  CALL multigrid_construct_surf(mg_mesh)
  IF(mg_mesh%smesh%cad_type/=mesh_native_id)CALL oft_abort('Wrong mesh type.','test_cubit',__FILE__)
  CALL plot_file%setup("Test")
  CALL mg_mesh%smesh%setup_io(plot_file,1)
  IF(oft_env%head_proc)THEN
    OPEN(NEWUNIT=io_unit,FILE='cubit.results')
    WRITE(io_unit,*)0.0_r8
  END IF
  CALL compute_area
  IF(oft_env%head_proc)THEN
    WRITE(io_unit,*)mg_mesh%smesh%global%nbp
    WRITE(io_unit,*)mg_mesh%smesh%global%nbe
    CLOSE(io_unit)
  END IF
ELSE
  !---Setup grid
  CALL multigrid_construct(mg_mesh)
  IF(mg_mesh%mesh%cad_type/=mesh_native_id)CALL oft_abort('Wrong mesh type.','test_cubit',__FILE__)
  CALL check_bmesh_ho
  CALL plot_file%setup("Test")
  CALL mg_mesh%mesh%setup_io(plot_file,1)
  IF(oft_env%head_proc)OPEN(NEWUNIT=io_unit,FILE='cubit.results')
  CALL compute_volume
  CALL compute_area
  IF(oft_env%head_proc)THEN
    WRITE(io_unit,*)mg_mesh%mesh%global%nbp
    WRITE(io_unit,*)mg_mesh%mesh%global%nbe
    CLOSE(io_unit)
  END IF
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
! Check that the boundary mesh high-order nodes match the volume mesh on every level
SUBROUTINE check_bmesh_ho
INTEGER(i4) :: i,j,level
REAL(r8) :: max_err
CLASS(oft_mesh), POINTER :: mesh
max_err=0._r8
DO level=1,mg_mesh%mgmax
  mesh=>mg_mesh%meshes(level)
  IF(.NOT.ASSOCIATED(mesh%ho_info%r))CYCLE
  DO i=1,mesh%bmesh%ne
    j=ABS(mesh%bmesh%parent%le(i))
    max_err=MAX(max_err,MAXVAL(ABS(mesh%bmesh%ho_info%r(:,mesh%bmesh%ho_info%lep(1,i)) &
      - mesh%ho_info%r(:,mesh%ho_info%lep(1,j)))))
  END DO
  IF(mesh%ho_info%nfp==1)THEN
    DO i=1,mesh%bmesh%nc
      j=ABS(mesh%bmesh%parent%lf(i))
      max_err=MAX(max_err,MAXVAL(ABS(mesh%bmesh%ho_info%r(:,mesh%bmesh%ho_info%lcp(1,i)) &
        - mesh%ho_info%r(:,mesh%ho_info%lfp(1,j)))))
    END DO
  END IF
END DO
max_err=oft_mpi_max(max_err)
IF(oft_env%head_proc)WRITE(*,*)'Max boundary high-order node error =',max_err
IF(max_err>1.E-12_r8)CALL oft_abort('Boundary high-order nodes do not match volume mesh', &
  'check_bmesh_ho',__FILE__)
END SUBROUTINE check_bmesh_ho
!------------------------------------------------------------------------------
!> Compute volume of the current mesh and output
!------------------------------------------------------------------------------
SUBROUTINE compute_volume
INTEGER(i4) :: i,m
REAL(r8) :: v,det,goptmp(3,4),volume
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
volume=0._r8
CALL mg_mesh%mesh%quad_rule(8,quad)
!---
DO i=1,mg_mesh%mesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL mg_mesh%mesh%jacobian(i,quad%pts(:,m),goptmp,v)
    det=v*quad%wts(m)
    volume=volume+det
  END DO
END DO
CALL quad%delete
volume=oft_mpi_sum(volume)
IF(oft_env%head_proc)THEN
  WRITE(*,*)'Mesh Volume =',volume
  WRITE(io_unit,*)volume
END IF
END SUBROUTINE compute_volume
!------------------------------------------------------------------------------
!> Compute surface area of the current mesh and output
!------------------------------------------------------------------------------
SUBROUTINE compute_area
INTEGER(i4) :: i,j,m
REAL(r8) :: a,det,goptmp(3,4),area
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
area=0._r8
CALL mg_mesh%smesh%quad_rule(8,quad)
!---
DO i=1,mg_mesh%smesh%nc
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    CALL mg_mesh%smesh%jacobian(i,quad%pts(:,m),goptmp,a)
    det=a*quad%wts(m)
    area=area+det
  END DO
END DO
CALL quad%delete
area=oft_mpi_sum(area)
IF(oft_env%head_proc)THEN
  WRITE(*,*)'Mesh Area   =',area
  WRITE(io_unit,*)area
END IF
END SUBROUTINE compute_area
END PROGRAM test_cubit
