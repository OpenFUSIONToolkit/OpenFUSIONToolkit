!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_mesh_sphere.F90
!
!> Mesh handling for a spherical test mesh.
!!
!! Functions to create and refine a spherical mesh
!! - Mesh setup
!! - Boundary point location and refinement
!!
!! @authors George Marklin and Chris Hansen
!! @date April 2008 - Present
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_mesh_sphere
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
USE oft_mesh_local_util, ONLY: mesh_local_findface
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_hexmesh_type, ONLY: oft_hexmesh
USE oft_quadmesh_type, ONLY: oft_quadmesh
USE multigrid, ONLY: multigrid_mesh, multigrid_level
IMPLICIT NONE
#include "local.h"
private
INTEGER(i4), PARAMETER, PUBLIC :: mesh_sphere_id = 91
public mesh_sphere_load, mesh_sphere_cadlink, mesh_sphere_reffix
public mesh_sphere_add_quad, smesh_circle_load, smesh_circle_cadlink
public smesh_circle_reffix, smesh_circle_add_quad
contains
!---------------------------------------------------------------------------------
!> Setup a spherical test mesh
!! The mesh is initialized with a minimal "wheel" of cells
!! - 7 Points
!! - 8 Cells
!---------------------------------------------------------------------------------
subroutine mesh_sphere_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,ierr,io_unit,mesh_type
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
namelist/sphere_options/mesh_type
DEBUG_STACK_PUSH
mesh_type = 1
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,sphere_options,IOSTAT=ierr)
  CLOSE(io_unit)
  !IF(ierr<0)CALL oft_abort('No "sphere_options" found in input file.','mesh_sphere_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "sphere_options" in input file.','mesh_sphere_load',__FILE__)
  WRITE(*,'(2A)')oft_indent,'Sphere volume mesh:'
  CALL oft_increase_indent
  WRITE(*,'(2A,I4)')oft_indent,'Mesh Type = ',mesh_type
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(mesh_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_sphere_load',__FILE__)
#endif
!
IF(mesh_type==1)THEN
  allocate(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
  allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%meshes(i)%setup(mesh_sphere_id)
    CALL mg_mesh%smeshes(i)%setup(mesh_sphere_id,.TRUE.)
    mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
  END DO
  CALL multigrid_level(mg_mesh,1)
  mesh=>mg_mesh%meshes(1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !
    !---Setup points
    !
    mesh%np=7
    allocate(mesh%r(3,mesh%np))
    !---Center point
    mesh%r(:,1)=(/0.d0,0.d0,0.d0/)
    !---Midplane points
    mesh%r(:,2)=(/1.d0,0.d0,0.d0/)
    mesh%r(:,3)=(/0.d0,1.d0,0.d0/)
    mesh%r(:,4)=(/-1.d0,0.d0,0.d0/)
    mesh%r(:,5)=(/0.d0,-1.d0,0.d0/)
    !---Pole points
    mesh%r(:,6)=(/0.d0,0.d0,1.d0/)
    mesh%r(:,7)=(/0.d0,0.d0,-1.d0/)
    !
    !---Setup cells
    !
    mesh%nc=8
    allocate(mesh%lc(4,mesh%nc),mesh%reg(mesh%nc))
    mesh%reg=1
    !---Top wheel
    mesh%lc(:,1)=(/1,2,3,6/)
    mesh%lc(:,2)=(/1,3,4,6/)
    mesh%lc(:,3)=(/1,4,5,6/)
    mesh%lc(:,4)=(/1,5,2,6/)
    !---Bottom wheel
    mesh%lc(:,5)=(/1,2,3,7/)
    mesh%lc(:,6)=(/1,3,4,7/)
    mesh%lc(:,7)=(/1,4,5,7/)
    mesh%lc(:,8)=(/1,5,2,7/)
  END IF
ELSE
  allocate(oft_hexmesh::mg_mesh%meshes(mg_mesh%mgdim))
  allocate(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%meshes(i)%setup(mesh_sphere_id)
    CALL mg_mesh%smeshes(i)%setup(mesh_sphere_id,.TRUE.)
    mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
  END DO
  CALL multigrid_level(mg_mesh,1)
  mesh=>mg_mesh%meshes(1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !---Setup points
    mesh%np=16
    allocate(mesh%r(3,mesh%np))
    !
    mesh%r(:,1)=(/-0.5d0,0.5d0,0.5d0/)*0.5d0
    mesh%r(:,2)=(/ 0.5d0,0.5d0,0.5d0/)*0.5d0
    mesh%r(:,3)=(/0.5d0,-0.5d0,0.5d0/)*0.5d0
    mesh%r(:,4)=(/-0.5d0,-0.5d0,0.5d0/)*0.5d0
    !
    mesh%r(:,5)=(/-0.5d0,0.5d0,-0.5d0/)*0.5d0
    mesh%r(:,6)=(/ 0.5d0,0.5d0,-0.5d0/)*0.5d0
    mesh%r(:,7)=(/0.5d0,-0.5d0,-0.5d0/)*0.5d0
    mesh%r(:,8)=(/-0.5d0,-0.5d0,-0.5d0/)*0.5d0
    !
    mesh%r(:,9)= (/-0.5d0,0.5d0,0.5d0/)/SQRT(0.75d0)
    mesh%r(:,10)=(/ 0.5d0,0.5d0,0.5d0/)/SQRT(0.75d0)
    mesh%r(:,11)=(/0.5d0,-0.5d0,0.5d0/)/SQRT(0.75d0)
    mesh%r(:,12)=(/-0.5d0,-0.5d0,0.5d0/)/SQRT(0.75d0)
    !
    mesh%r(:,13)=(/-0.5d0,0.5d0,-0.5d0/)/SQRT(0.75d0)
    mesh%r(:,14)=(/ 0.5d0,0.5d0,-0.5d0/)/SQRT(0.75d0)
    mesh%r(:,15)=(/0.5d0,-0.5d0,-0.5d0/)/SQRT(0.75d0)
    mesh%r(:,16)=(/-0.5d0,-0.5d0,-0.5d0/)/SQRT(0.75d0)
    !---Setup cells
    mesh%nc=7
    allocate(mesh%lc(8,mesh%nc),mesh%reg(mesh%nc))
    mesh%reg=1
    !
    mesh%lc(:,1)=(/1,2,3,4, 5,6,7,8/)
    mesh%lc(:,2)=(/1,2,3,4, 9,10,11,12/)
    mesh%lc(:,3)=(/1,5,6,2, 9,13,14,10/)
    mesh%lc(:,4)=(/2,6,7,3, 10,14,15,11/)
    mesh%lc(:,5)=(/3,7,8,4, 11,15,16,12/)
    mesh%lc(:,6)=(/1,4,8,5, 9,12,16,13/)
    mesh%lc(:,7)=(/5,8,7,6, 13,16,15,14/)
    !   1,2,3,4, 1,5,6,2, 2,6,7,3, 3,7,8,4, 1,4,8,5, 5,8,7,6/),(/4,6/))
  END IF
END IF
call mesh_global_resolution(mesh)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_sphere_load
!---------------------------------------------------------------------------------
!> Setup surface IDs
!---------------------------------------------------------------------------------
subroutine mesh_sphere_cadlink(mesh)
class(oft_mesh), intent(inout) :: mesh
mesh%bfs=1
end subroutine mesh_sphere_cadlink
!---------------------------------------------------------------------------------
!> Refine boundary points onto the sphere
!---------------------------------------------------------------------------------
subroutine mesh_sphere_reffix(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j
real(r8) :: u1,v1,u2,v2,u,v,pt(3),r
class(oft_mesh), pointer :: pmesh,mesh
DEBUG_STACK_PUSH
!---Get parent mesh
mesh=>mg_mesh%mesh
pmesh=>mg_mesh%meshes(mg_mesh%level-1)
IF(pmesh%fullmesh.AND.(.NOT.mesh%fullmesh))THEN
  ! Do nothing
ELSE
  if(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Adjusting points to sphere boundary'
  !---Locate edge end points and place daughter point
  !$omp parallel do private(pt)
  do i=1,mesh%nbp
    if(mesh%global%gbp(mesh%lbp(i)))then
      pt=mesh%r(:,mesh%lbp(i))
      mesh%r(:,mesh%lbp(i))=pt/sqrt(sum(pt**2))
    endif
  enddo
  ! !---Locate edge end points and place daughter node
  ! !$omp parallel do private(j)
  ! do i=1,mesh%bmesh%np
  !   IF(i<=pmesh%bmesh%np)CYCLE
  !   mesh%bmesh%r(:,i)=mesh%r(:,mesh%bmesh%parent%lp(i))
  ! enddo
END IF
DEBUG_STACK_POP
end subroutine mesh_sphere_reffix
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points
!---------------------------------------------------------------------------------
subroutine mesh_sphere_add_quad(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,k
real(r8) :: pt(3)
DEBUG_STACK_PUSH
if(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting Sphere Quadratic Nodes'
!---Setup quadratic mesh
! CALL mesh%set_order(2)
! CALL mesh_global_set_curved(mesh,2)
!---Locate edge end points and place daughter point
!!$omp parallel do private(j,u1,v1,u2,v2,pt)
do i=1,mesh%ne
  mesh%ho_info%lep(1,i)=i
  mesh%ho_info%r(:,i)=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
  if(mesh%global%gbe(i))then
    pt=mesh%ho_info%r(:,i)
    mesh%ho_info%r(:,i)=pt/sqrt(sum(pt**2))
    ! DO k=mesh%kec(i),mesh%kec(i+1)-1
    !   mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
    ! END DO
  endif
enddo
IF(mesh%ho_info%nfp==1)THEN
  do i=1,mesh%nf
    mesh%ho_info%lfp(1,i)=i+mesh%ne
    pt=0.d0
    DO j=1,mesh%face_np
      pt=pt+mesh%r(:,mesh%lf(j,i))
    END DO
    mesh%ho_info%r(:,i+mesh%ne)=pt/REAL(mesh%face_np,8)
    if(mesh%global%gbf(i))then
      pt=mesh%ho_info%r(:,i+mesh%ne)
      mesh%ho_info%r(:,i+mesh%ne)=pt/sqrt(sum(pt**2))
      ! mesh%ho_info%is_curved(mesh%lfc(1,i))=.TRUE.
    endif
  end do
END IF
IF(mesh%ho_info%ncp==1)THEN
  do i=1,mesh%nc
    mesh%ho_info%lcp(1,i)=i+mesh%ne+mesh%nf
    pt=0.d0
    DO j=1,mesh%cell_np
      pt=pt+mesh%r(:,mesh%lc(j,i))
    END DO
    mesh%ho_info%r(:,i+mesh%ne+mesh%nf)=pt/REAL(mesh%cell_np,8)
  end do
END IF
DEBUG_STACK_POP
end subroutine mesh_sphere_add_quad
!---------------------------------------------------------------------------------
!> Setup a spherical test mesh
!! The mesh is initialized with a minimal "wheel" of cells
!! - 7 Points
!! - 8 Cells
!---------------------------------------------------------------------------------
subroutine smesh_circle_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,ierr,io_unit,mesh_type
class(oft_bmesh), pointer :: smesh
namelist/sphere_options/mesh_type
DEBUG_STACK_PUSH
mesh_type = 1
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,sphere_options,IOSTAT=ierr)
  CLOSE(io_unit)
  !IF(ierr<0)CALL oft_abort('No "sphere_options" found in input file.','smesh_circle_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "sphere_options" in input file.','smesh_circle_load',__FILE__)
  WRITE(*,'(2A)')oft_indent,'Circle surface mesh:'
  CALL oft_increase_indent
  WRITE(*,'(2A,I4)')oft_indent,'Mesh Type = ',mesh_type
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(mesh_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_circle_load',__FILE__)
#endif
!
IF(mesh_type==1)THEN
  allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%smeshes(i)%setup(mesh_sphere_id,.FALSE.)
  END DO
  CALL multigrid_level(mg_mesh,1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !---Setup points
    smesh%np=5
    allocate(smesh%r(3,smesh%np))
    smesh%r(:,1)=[0.d0,0.d0,0.d0]
    smesh%r(:,2)=[1.d0,0.d0,0.d0]
    smesh%r(:,3)=[0.d0,1.d0,0.d0]
    smesh%r(:,4)=[-1.d0,0.d0,0.d0]
    smesh%r(:,5)=[0.d0,-1.d0,0.d0]
    !---Setup cells
    smesh%nc=4
    allocate(smesh%lc(3,smesh%nc),smesh%reg(smesh%nc))
    smesh%reg=1
    smesh%lc(:,1)=[2,3,1]
    smesh%lc(:,2)=[3,4,1]
    smesh%lc(:,3)=[4,5,1]
    smesh%lc(:,4)=[5,2,1]
  END IF
ELSE
  allocate(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%smeshes(i)%setup(mesh_sphere_id,.FALSE.)
  END DO
  CALL multigrid_level(mg_mesh,1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !---Setup points
    smesh%np=8
    allocate(smesh%r(3,smesh%np))
    smesh%r(:,1)=[1.d0,1.d0,0.d0]*0.5d0/SQRT(2.d0)
    smesh%r(:,2)=[-1.d0,1.d0,0.d0]*0.5d0/SQRT(2.d0)
    smesh%r(:,3)=[-1.d0,-1.d0,0.d0]*0.5d0/SQRT(2.d0)
    smesh%r(:,4)=[1.d0,-1.d0,0.d0]*0.5d0/SQRT(2.d0)
    smesh%r(:,5)=[1.d0,1.d0,0.d0]/SQRT(2.d0)
    smesh%r(:,6)=[-1.d0,1.d0,0.d0]/SQRT(2.d0)
    smesh%r(:,7)=[-1.d0,-1.d0,0.d0]/SQRT(2.d0)
    smesh%r(:,8)=[1.d0,-1.d0,0.d0]/SQRT(2.d0)
    !---Setup cells
    smesh%nc=5
    allocate(smesh%lc(4,smesh%nc),smesh%reg(smesh%nc))
    smesh%reg=1
    !---
    smesh%lc(:,1)=[1,2,3,4]
    smesh%lc(:,2)=[1,2,6,5]
    smesh%lc(:,3)=[2,3,7,6]
    smesh%lc(:,4)=[3,4,8,7]
    smesh%lc(:,5)=[4,1,5,8]
  END IF
END IF
call mesh_global_resolution(smesh)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine smesh_circle_load
!---------------------------------------------------------------------------------
!> Setup surface IDs
!---------------------------------------------------------------------------------
subroutine smesh_circle_cadlink(smesh)
class(oft_bmesh), intent(inout) :: smesh
smesh%bes=1
end subroutine smesh_circle_cadlink
!---------------------------------------------------------------------------------
!> Refine boundary points onto the sphere
!---------------------------------------------------------------------------------
subroutine smesh_circle_reffix(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j
real(r8) :: u1,v1,u2,v2,u,v,pt(3),r
class(oft_bmesh), pointer :: pmesh,smesh
DEBUG_STACK_PUSH
!---Get parent mesh
smesh=>mg_mesh%smesh
pmesh=>mg_mesh%smeshes(mg_mesh%level-1)
IF(pmesh%fullmesh.AND.(.NOT.smesh%fullmesh))THEN
  ! Do nothing
ELSE
  if(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Adjusting points to circle boundary'
  !---Locate edge end points and place daughter point
  !$omp parallel do private(pt)
  do i=1,smesh%nbp
    if(smesh%global%gbp(smesh%lbp(i)))then
      pt=smesh%r(:,smesh%lbp(i))
      smesh%r(:,smesh%lbp(i))=pt/sqrt(sum(pt**2))
    endif
  enddo
END IF
DEBUG_STACK_POP
end subroutine smesh_circle_reffix
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points
!---------------------------------------------------------------------------------
subroutine smesh_circle_add_quad(smesh)
class(oft_bmesh), intent(inout) :: smesh
integer(i4) :: i,j,k
real(r8) :: pt(3)
DEBUG_STACK_PUSH
if(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting circle quadratic nodes'
!---Setup quadratic mesh
CALL smesh%set_order(2)
!---Locate edge end points and place daughter point
!!$omp parallel do private(j,u1,v1,u2,v2,pt)
do i=1,smesh%ne
  smesh%ho_info%lep(1,i)=i
  smesh%ho_info%r(:,i)=(smesh%r(:,smesh%le(1,i))+smesh%r(:,smesh%le(2,i)))/2.d0
  if(smesh%global%gbe(i))then
    pt=smesh%ho_info%r(:,i)
    smesh%ho_info%r(:,i)=pt/sqrt(sum(pt**2))
    DO k=smesh%kec(i),smesh%kec(i+1)-1
      smesh%ho_info%is_curved(smesh%lec(k))=.TRUE.
    END DO
  endif
enddo
IF(smesh%ho_info%ncp>0)THEN
  do i=1,smesh%nc
    smesh%ho_info%lcp(1,i)=i+smesh%ne
    pt=0.d0
    DO j=1,smesh%cell_np
      pt=pt+smesh%r(:,smesh%lc(j,i))
    END DO
    smesh%ho_info%r(:,i+smesh%ne)=pt/REAL(smesh%cell_np,8)
    ! if(smesh%global%gbc(i))then
    !   pt=smesh%ho_info%r(:,i+smesh%ne)
    !   smesh%ho_info%r(:,i+smesh%ne)=pt/sqrt(sum(pt**2))
    !   smesh%ho_info%is_curved(i)=.TRUE.
    ! endif
  end do
END IF
DEBUG_STACK_POP
end subroutine smesh_circle_add_quad
end module oft_mesh_sphere
