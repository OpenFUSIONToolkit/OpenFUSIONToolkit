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
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
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
  CALL check_smesh_ho
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
! Check that each level's boundary mesh matches the volume mesh: the boundary quadratic mapping
! at each boundary cell's nodes must equal the volume quadratic mapping at the same points of the
! volume cell behind that face. Corners are matched through vertex indices only, so this does not
! reuse the edge/face lookup that fills the boundary nodes.
SUBROUTINE check_bmesh_ho
INTEGER(i4) :: i,j,k,c,v,level,morder,border
REAL(r8) :: max_err,fb(4,4),fv(4,4),gb(4),gv(4)
LOGICAL, ALLOCATABLE :: mcurved(:),bcurved(:)
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_bmesh), POINTER :: bmesh
max_err=0._r8
DO level=1,mg_mesh%mgmax
  mesh=>mg_mesh%meshes(level)
  bmesh=>mesh%bmesh
  IF(.NOT.ASSOCIATED(mesh%ho_info%r))CYCLE
  IF(.NOT.ASSOCIATED(bmesh%ho_info%r))CALL oft_abort('Boundary high-order nodes missing', &
    'check_bmesh_ho',__FILE__)
  morder=mesh%order
  border=bmesh%order
  CALL mesh%set_order(2)
  CALL bmesh%set_order(2)
  ! Use the quadratic map on every cell of both meshes
  mcurved=mesh%ho_info%is_curved
  bcurved=bmesh%ho_info%is_curved
  mesh%ho_info%is_curved=.TRUE.
  bmesh%ho_info%is_curved=.TRUE.
  DO i=1,bmesh%nc
    c=mesh%lfc(1,ABS(bmesh%parent%lf(i)))
    ! Logical position of each boundary cell corner in the boundary (fb) and volume (fv) cells
    fb=0._r8
    fv=-1._r8
    DO k=1,bmesh%cell_np
      CALL bmesh%vlog(k,fb(:,k))
      v=bmesh%parent%lp(bmesh%lc(k,i))
      DO j=1,mesh%cell_np
        IF(mesh%lc(j,c)==v)THEN
          fv(:,k)=0._r8
          CALL mesh%vlog(j,fv(:,k))
        END IF
      END DO
      IF(ALL(fv(:,k)<0._r8))CALL oft_abort('Boundary vertex not found in volume cell', &
        'check_bmesh_ho',__FILE__)
    END DO
    ! Corners
    DO k=1,bmesh%cell_np
      max_err=MAX(max_err,MAXVAL(ABS(bmesh%log2phys(i,fb(:,k)) - mesh%log2phys(c,fv(:,k)))))
    END DO
    ! Edge midpoints
    DO k=1,bmesh%cell_ne
      gb=(fb(:,bmesh%cell_ed(1,k))+fb(:,bmesh%cell_ed(2,k)))/2._r8
      gv=(fv(:,bmesh%cell_ed(1,k))+fv(:,bmesh%cell_ed(2,k)))/2._r8
      max_err=MAX(max_err,MAXVAL(ABS(bmesh%log2phys(i,gb) - mesh%log2phys(c,gv))))
    END DO
    ! Cell (face) centre
    IF(bmesh%ho_info%ncp==1)THEN
      gb=SUM(fb(:,1:bmesh%cell_np),DIM=2)/REAL(bmesh%cell_np,8)
      gv=SUM(fv(:,1:bmesh%cell_np),DIM=2)/REAL(bmesh%cell_np,8)
      max_err=MAX(max_err,MAXVAL(ABS(bmesh%log2phys(i,gb) - mesh%log2phys(c,gv))))
    END IF
  END DO
  mesh%ho_info%is_curved=mcurved
  bmesh%ho_info%is_curved=bcurved
  CALL mesh%set_order(morder)
  CALL bmesh%set_order(border)
END DO
max_err=oft_mpi_max(max_err)
IF(oft_env%head_proc)WRITE(*,*)'Max boundary high-order node error =',max_err
IF(max_err>1.E-12_r8)CALL oft_abort('Boundary high-order nodes do not match volume mesh', &
  'check_bmesh_ho',__FILE__)
END SUBROUTINE check_bmesh_ho
! Check that each surface level's high-order nodes match its parent level: copied nodes on the
! transfer level, and the parent's quadratic mapping at each child cell's nodes on refined levels
SUBROUTINE check_smesh_ho
INTEGER(i4) :: i,j,k,l,v,cc,level,porder,sorder
REAL(r8) :: max_err,f1(4),fc(4,4),fp(4,4),gc(4),gp(4)
LOGICAL, ALLOCATABLE :: pcurved(:),scurved(:)
CLASS(oft_bmesh), POINTER :: smesh,pmesh
max_err=0._r8
DO level=2,mg_mesh%mgmax
  smesh=>mg_mesh%smeshes(level)
  pmesh=>mg_mesh%smeshes(level-1)
  ! Skip unless both levels have high-order nodes (with a linear input mesh only the top level does)
  IF(.NOT.ASSOCIATED(smesh%ho_info%r))CYCLE
  IF(.NOT.ASSOCIATED(pmesh%ho_info%r))CYCLE
  IF(pmesh%fullmesh.AND.(.NOT.smesh%fullmesh))THEN
    !---Transfer level: nodes are copied from the parent level
    DO i=1,smesh%ne
      j=ABS(smesh%base%le(i))
      max_err=MAX(max_err,MAXVAL(ABS(smesh%ho_info%r(:,smesh%ho_info%lep(1,i)) &
        - pmesh%ho_info%r(:,pmesh%ho_info%lep(1,j)))))
    END DO
    IF(smesh%ho_info%ncp==1)THEN
      DO i=1,smesh%nc
        j=ABS(smesh%base%lc(i))
        max_err=MAX(max_err,MAXVAL(ABS(smesh%ho_info%r(:,smesh%ho_info%lcp(1,i)) &
          - pmesh%ho_info%r(:,pmesh%ho_info%lcp(1,j)))))
      END DO
    END IF
    CYCLE
  END IF
  !---Refined level: compare the child and parent quadratic mappings at the child's nodes
  porder=pmesh%order
  sorder=smesh%order
  CALL pmesh%set_order(2)
  CALL smesh%set_order(2)
  ! Use the quadratic map on every cell of both levels
  pcurved=pmesh%ho_info%is_curved
  scurved=smesh%ho_info%is_curved
  pmesh%ho_info%is_curved=.TRUE.
  smesh%ho_info%is_curved=.TRUE.
  DO i=1,pmesh%nc
    DO l=1,4
      cc=(i-1)*4+l
      !---Logical position of each child corner in the child (fc) and parent (fp) cells
      fc=0._r8
      fp=0._r8
      DO k=1,smesh%cell_np
        CALL smesh%vlog(k,fc(:,k))
        v=smesh%lc(k,cc)
        DO j=1,pmesh%cell_np
          f1=0._r8
          CALL pmesh%vlog(j,f1)
          IF(v<=pmesh%np)THEN
            ! Parent vertex
            IF(pmesh%lc(j,i)==v)fp(:,k)=f1
          ELSE IF(v<=pmesh%np+pmesh%ne)THEN
            ! Parent edge midpoint
            IF(ANY(pmesh%lc(j,i)==pmesh%le(:,v-pmesh%np)))fp(:,k)=fp(:,k)+f1/2._r8
          ELSE
            ! Parent cell centre
            fp(:,k)=fp(:,k)+f1/REAL(pmesh%cell_np,8)
          END IF
        END DO
      END DO
      !---Corners
      DO k=1,smesh%cell_np
        max_err=MAX(max_err,MAXVAL(ABS(smesh%log2phys(cc,fc(:,k)) - pmesh%log2phys(i,fp(:,k)))))
      END DO
      !---Edge midpoints
      DO k=1,smesh%cell_ne
        gc=(fc(:,smesh%cell_ed(1,k))+fc(:,smesh%cell_ed(2,k)))/2._r8
        gp=(fp(:,smesh%cell_ed(1,k))+fp(:,smesh%cell_ed(2,k)))/2._r8
        max_err=MAX(max_err,MAXVAL(ABS(smesh%log2phys(cc,gc) - pmesh%log2phys(i,gp))))
      END DO
      !---Cell centre
      IF(smesh%ho_info%ncp==1)THEN
        gc=SUM(fc,DIM=2)/REAL(smesh%cell_np,8)
        gp=SUM(fp,DIM=2)/REAL(smesh%cell_np,8)
        max_err=MAX(max_err,MAXVAL(ABS(smesh%log2phys(cc,gc) - pmesh%log2phys(i,gp))))
      END IF
    END DO
  END DO
  pmesh%ho_info%is_curved=pcurved
  smesh%ho_info%is_curved=scurved
  CALL pmesh%set_order(porder)
  CALL smesh%set_order(sorder)
END DO
max_err=oft_mpi_max(max_err)
IF(oft_env%head_proc)WRITE(*,*)'Max surface high-order node error =',max_err
IF(max_err>1.E-12_r8)CALL oft_abort('Surface high-order nodes do not match parent level', &
  'check_smesh_ho',__FILE__)
END SUBROUTINE check_smesh_ho
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
