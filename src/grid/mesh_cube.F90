!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_mesh_cube.F90
!
!> Mesh handling for a 1x1x1 cube test mesh.
!!
!! Functions to create and refine a 1x1x1 cube mesh
!! - Mesh setup
!! - Boundary point location and refinement
!!
!! @authors Chris Hansen
!! @date August 2012
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_mesh_cube
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
PRIVATE
INTEGER(i4), PARAMETER, PUBLIC :: mesh_cube_id = 92
LOGICAL :: ref_per(3) = .FALSE. !< Character flag for periodic reflections
public mesh_cube_load, mesh_cube_cadlink, mesh_cube_set_periodic
public smesh_square_load, smesh_square_cadlink
public smesh_square_set_periodic
contains
!---------------------------------------------------------------------------------
!> Setup a 1x1x1 cube test mesh
!! The mesh is initialized with a basic set of cells
!! - 9 Points
!! - 12 Cells
!---------------------------------------------------------------------------------
subroutine mesh_cube_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,j,k,ierr,io_unit,mesh_type,ni(3)
INTEGER(i4), ALLOCATABLE :: pmap(:,:,:)
REAL(r8) :: rscale(3),shift(3),packing(3),alpha,beta,xtmp,ytmp,ztmp
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
namelist/cube_options/mesh_type,ni,rscale,shift,ref_per,packing
DEBUG_STACK_PUSH
mesh_type=1
ni=1
rscale=1.d0
shift=0.d0
packing=1.d0
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,cube_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "cube_options" found in input file.','mesh_cube_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "cube_options" in input file.','mesh_cube_load',__FILE__)
  WRITE(*,'(2A)')oft_indent,'Cube volume mesh:'
  CALL oft_increase_indent
  WRITE(*,'(2A,I4)')oft_indent,     'Mesh Type   = ',mesh_type
  WRITE(*,'(2A,3I4)')oft_indent,    'nx, ny, nz  = ',ni
  WRITE(*,'(2A,3ES11.3)')oft_indent,'Scale facs  = ',rscale
  WRITE(*,'(2A,3L2)')oft_indent,    'Periodicity = ',ref_per
  WRITE(*,'(2A,3ES11.3)')oft_indent,'Packing     = ',packing
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(mesh_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cube_load',__FILE__)
CALL MPI_Bcast(ni,3,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cube_load',__FILE__)
CALL MPI_Bcast(rscale,3,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cube_load',__FILE__)
CALL MPI_Bcast(ref_per,  3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cube_load',__FILE__)
#endif
!--Build mesh
IF(mesh_type==1)THEN
  allocate(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
  allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%meshes(i)%setup(mesh_cube_id)
    CALL mg_mesh%smeshes(i)%setup(mesh_cube_id,.TRUE.)
    mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
  END DO
  CALL multigrid_level(mg_mesh,1)
  mesh=>mg_mesh%meshes(1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !---Setup points
    mesh%np=9
    allocate(mesh%r(3,mesh%np))
    !---
    mesh%r(:,1)=(/1.d0,0.d0,1.d0/)
    mesh%r(:,2)=(/1.d0,1.d0,1.d0/)
    mesh%r(:,3)=(/0.d0,1.d0,1.d0/)
    mesh%r(:,4)=(/0.d0,0.d0,1.d0/)
    !---
    mesh%r(:,5)=(/1.d0,1.d0,0.d0/)
    mesh%r(:,6)=(/1.d0,0.d0,0.d0/)
    mesh%r(:,7)=(/0.d0,0.d0,0.d0/)
    mesh%r(:,8)=(/0.d0,1.d0,0.d0/)
    mesh%r(:,9)=(/.5d0,.5d0,.5d0/)
    !---Setup cells
    mesh%nc=12
    allocate(mesh%lc(4,mesh%nc))
    !---
    mesh%lc(:,1)=(/7,3,9,8/)
    mesh%lc(:,2)=(/4,3,9,7/)
    mesh%lc(:,3)=(/6,4,9,7/)
    mesh%lc(:,4)=(/2,3,9,1/)
    mesh%lc(:,5)=(/1,4,9,6/)
    mesh%lc(:,6)=(/3,8,2,9/)
    mesh%lc(:,7)=(/1,3,9,4/)
    mesh%lc(:,8)=(/7,9,5,8/)
    mesh%lc(:,9)=(/5,2,8,9/)
    mesh%lc(:,10)=(/2,9,5,1/)
    mesh%lc(:,11)=(/9,6,5,1/)
    mesh%lc(:,12)=(/7,6,5,9/)
    !---
    DO i=1,3
      IF(ref_per(i))THEN
        CALL reflect_tet(i)
        mesh%r(i,:)=(mesh%r(i,:)+1.d0)/2.d0
      END IF
    END DO
    !---
    mesh%r(1,:)=mesh%r(1,:)*rscale(1)+shift(1)
    mesh%r(2,:)=mesh%r(2,:)*rscale(2)+shift(2)
    mesh%r(3,:)=mesh%r(3,:)*rscale(3)+shift(3)
    allocate(mesh%reg(mesh%nc))
    mesh%reg=1
  END IF
ELSE
  allocate(oft_hexmesh::mg_mesh%meshes(mg_mesh%mgdim))
  allocate(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
  DO i=1,mg_mesh%mgdim
    CALL mg_mesh%meshes(i)%setup(mesh_cube_id)
    CALL mg_mesh%smeshes(i)%setup(mesh_cube_id,.TRUE.)
    mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
  END DO
  CALL multigrid_level(mg_mesh,1)
  mesh=>mg_mesh%meshes(1)
  smesh=>mg_mesh%smeshes(1)
  IF(oft_env%rank==0)THEN
    !---Setup points
    mesh%np=PRODUCT(ni+1)
    allocate(mesh%r(3,mesh%np))
    ALLOCATE(pmap(ni(1)+1,ni(2)+1,ni(3)+1))
    pmap=-1
    mesh%np=0
    ! Adjust scales for packing
    beta = 2.0*(packing(1)-1.0); alpha = -2.0*beta/3.0
    rscale(1) = rscale(1)/(alpha + beta + 1.d0)
    beta = 2.0*(packing(2)-1.0); alpha = -2.0*beta/3.0
    rscale(2) = rscale(2)/(alpha + beta + 1.d0)
    beta = 2.0*(packing(3)-1.0); alpha = -2.0*beta/3.0
    rscale(3) = rscale(3)/(alpha + beta + 1.d0)
    DO i=1,ni(1)+1
      beta = 2.0*(packing(1)-1.0); alpha = -2.0*beta/3.0
      xtmp = (i-1)/REAL(ni(1),8)
      xtmp = rscale(1)*(alpha*(xtmp**3) + beta*(xtmp**2) + xtmp)
      DO j=1,ni(2)+1
        beta = 2.0*(packing(2)-1.0); alpha = -2.0*beta/3.0
        ytmp = (j-1)/REAL(ni(2),8)
        ytmp = rscale(2)*(alpha*(ytmp**3) + beta*(ytmp**2) + ytmp)
        DO k=1,ni(3)+1
          beta = 2.0*(packing(3)-1.0); alpha = -2.0*beta/3.0
          ztmp = (k-1)/REAL(ni(3),8)
          ztmp = rscale(3)*(alpha*(ztmp**3) + beta*(ztmp**2) + ztmp)
          !
          mesh%np=mesh%np+1
          pmap(i,j,k)=mesh%np
          mesh%r(:,mesh%np)=(/xtmp,ytmp,ztmp/)+shift
        END DO
      END DO
    END DO
    !---Setup cells
    mesh%nc=PRODUCT(ni)
    allocate(mesh%lc(8,mesh%nc),mesh%reg(mesh%nc))
    mesh%reg=1
    mesh%nc=0
    DO i=1,ni(1)
      DO j=1,ni(2)
        DO k=1,ni(3)
          mesh%nc=mesh%nc+1
          mesh%lc(:,mesh%nc)=(/pmap(i,j,k),pmap(i+1,j,k),pmap(i+1,j+1,k),pmap(i,j+1,k), &
            pmap(i,j,k+1),pmap(i+1,j,k+1),pmap(i+1,j+1,k+1),pmap(i,j+1,k+1)/)
        END DO
      END DO
    END DO
    DEALLOCATE(pmap)
  END IF
END IF
call mesh_global_resolution(mesh)
CALL oft_decrease_indent
DEBUG_STACK_POP
contains
subroutine reflect_tet(ref_index)
integer(i4), intent(in) :: ref_index
integer(i4) :: npold,ncold,i,j
integer(i4), allocatable :: newindex(:),ltemp(:,:)
real(r8), allocatable :: rtemp(:,:)
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Reflecting to support periodicity'
npold=mesh%np
allocate(newindex(2*mesh%np),rtemp(3,2*mesh%np))
rtemp(:,1:mesh%np)=mesh%r
deallocate(mesh%r)
do i=1,npold
    IF(ABS(rtemp(ref_index,i))<=1.d-1)THEN
        rtemp(ref_index,i)=0.d0
        newindex(i)=i
    ELSE
        mesh%np=mesh%np+1
        rtemp(:,mesh%np) = rtemp(:,i)
        rtemp(ref_index,mesh%np) =-rtemp(ref_index,i)
        newindex(i)=mesh%np
    ENDIF
enddo
allocate(mesh%r(3,mesh%np))
mesh%r=rtemp(:,1:mesh%np)
deallocate(rtemp)
!---Reflect cells
ncold=mesh%nc
allocate(ltemp(mesh%cell_np,2*mesh%nc))
ltemp(:,1:mesh%nc)=mesh%lc
deallocate(mesh%lc)
do i=1,ncold
    mesh%nc=mesh%nc+1
    DO j=1,mesh%cell_np
        ltemp(j,mesh%nc)=newindex(ltemp(j,i))
    END DO
enddo
allocate(mesh%lc(mesh%cell_np,mesh%nc))
mesh%lc=ltemp(:,1:mesh%nc)
deallocate(ltemp)
end subroutine reflect_tet
end subroutine mesh_cube_load
!---------------------------------------------------------------------------------
!> Setup surface IDs
!---------------------------------------------------------------------------------
subroutine mesh_cube_cadlink(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j
DO i=1,mesh%nbf
  j=mesh%lbf(i)
  IF(ALL(mesh%r(1,mesh%lf(:,j))==1.d0))THEN
    mesh%bfs(i)=1
  ELSE IF(ALL(mesh%r(1,mesh%lf(:,j))==0.d0))THEN
    mesh%bfs(i)=2
  ELSE IF(ALL(mesh%r(2,mesh%lf(:,j))==1.d0))THEN
    mesh%bfs(i)=3
  ELSE IF(ALL(mesh%r(2,mesh%lf(:,j))==0.d0))THEN
    mesh%bfs(i)=4
  ELSE IF(ALL(mesh%r(3,mesh%lf(:,j))==1.d0))THEN
    mesh%bfs(i)=5
  ELSE IF(ALL(mesh%r(3,mesh%lf(:,j))==0.d0))THEN
    mesh%bfs(i)=6
  END IF
END DO
end subroutine mesh_cube_cadlink
!---------------------------------------------------------------------------------
!>
!---------------------------------------------------------------------------------
subroutine mesh_cube_set_periodic(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,jj,k,kk,l,m,n,iper
real(r8) :: pt_i(3),pt_j(3),d_plane,per_dir(3),i_cc(3),j_cc(3)
real(r8), parameter :: tol=1.d-6
logical :: flag(4)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Setting Cube periodicity'
!---Find periodic faces
mesh%periodic%nper=COUNT(ref_per)
ALLOCATE(mesh%periodic%lf(mesh%nf))
mesh%periodic%lf=-1
DO iper=1,3
  IF(.NOT.ref_per(iper))CYCLE
  DO jj=1,mesh%nbf
    j=mesh%lbf(jj)
    IF(mesh%periodic%lf(j)>0)CYCLE
    pt_i=cross_product(mesh%r(:,mesh%lf(2,j))-mesh%r(:,mesh%lf(1,j)), &
                       mesh%r(:,mesh%lf(3,j))-mesh%r(:,mesh%lf(1,j)))
    IF(MAXLOC(ABS(pt_i),DIM=1)/=iper)CYCLE
    i_cc=0.d0
    DO l=1,mesh%face_np
      i_cc=i_cc+mesh%r(:,mesh%lf(l,j))
    END DO
    i_cc=i_cc/REAL(mesh%face_np,8)
    i_cc(iper)=0.d0
    DO kk=1,mesh%nbf
      k=mesh%lbf(kk)
      IF(k<=j)CYCLE
      IF(mesh%periodic%lf(k)>0)CYCLE
      pt_j=cross_product(mesh%r(:,mesh%lf(2,k))-mesh%r(:,mesh%lf(1,k)), &
                         mesh%r(:,mesh%lf(3,k))-mesh%r(:,mesh%lf(1,k)))
      IF(MAXLOC(ABS(pt_j),DIM=1)==iper)THEN
        j_cc=0.d0
        DO l=1,mesh%face_np
          j_cc=j_cc+mesh%r(:,mesh%lf(l,k))
        END DO
        j_cc=j_cc/REAL(mesh%face_np,8)
        j_cc(iper)=0.d0
        IF(magnitude(i_cc-j_cc)<1.d-3)THEN
          mesh%periodic%lf(k)=j
          EXIT
        END IF
      END IF
    END DO
  END DO
END DO
!---Set periodic faces
ALLOCATE(mesh%periodic%lp(mesh%np))
ALLOCATE(mesh%periodic%le(mesh%ne))
mesh%periodic%lp=-1
mesh%periodic%le=-1
DO iper=1,3
  IF(.NOT.ref_per(iper))CYCLE
  per_dir=0.d0
  per_dir(iper)=1.d0
  !$omp parallel do private(j,i,l,k,n,m,pt_i,pt_j,d_plane,flag)
  DO jj=1,mesh%nbf
    j=mesh%lbf(jj)
    i=mesh%periodic%lf(j)
    IF(i<=0)CYCLE
    !---Check edges
    flag=.FALSE.
    DO n=1,mesh%face_np
      m=ABS(mesh%lfe(n,j))
      IF(mesh%periodic%le(m)>0)CYCLE
      pt_j=(mesh%r(:,mesh%le(1,m))+mesh%r(:,mesh%le(2,m)))/2.d0
      pt_j=pt_j-DOT_PRODUCT(pt_j,per_dir)*per_dir
      DO l=1,mesh%face_np
        IF(flag(l))CYCLE
        k=ABS(mesh%lfe(l,i))
        pt_i=(mesh%r(:,mesh%le(1,k))+mesh%r(:,mesh%le(2,k)))/2.d0
        pt_i=pt_i-DOT_PRODUCT(pt_i,per_dir)*per_dir
        !---
        pt_i=(pt_i-pt_j)
        d_plane=DOT_PRODUCT(pt_i,pt_i)
        IF(d_plane<tol)THEN
          flag(l)=.TRUE.
          mesh%periodic%le(m)=k
        END IF
      END DO
    END DO
    !---Check points
    flag=.FALSE.
    DO n=1,mesh%face_np
      m=ABS(mesh%lf(n,j))
      IF(mesh%periodic%lp(m)>0)CYCLE
      pt_j=mesh%r(:,m)
      pt_j=pt_j-DOT_PRODUCT(pt_j,per_dir)*per_dir
      DO l=1,mesh%face_np
        IF(flag(l))CYCLE
        k=ABS(mesh%lf(l,i))
        pt_i=mesh%r(:,k)
        pt_i=pt_i-DOT_PRODUCT(pt_i,per_dir)*per_dir
        !---
        pt_i=(pt_i-pt_j)
        d_plane=DOT_PRODUCT(pt_i,pt_i)
        IF(d_plane<tol)THEN
          flag(l)=.TRUE.
          mesh%periodic%lp(m)=k
        END IF
      END DO
    END DO
  END DO
END DO
DEBUG_STACK_POP
end subroutine mesh_cube_set_periodic
!---------------------------------------------------------------------------------
!> Setup a 1x1x1 cube test mesh
!! The mesh is initialized with a basic set of cells
!! - 9 Points
!! - 12 Cells
!---------------------------------------------------------------------------------
subroutine smesh_square_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,j,k,ierr,io_unit,nptmp,nctmp,mesh_type,ni(3)
INTEGER(i4), ALLOCATABLE :: pmap(:,:),lctmp(:,:)
REAL(r8) :: rscale(3),shift(3)
REAL(r8), ALLOCATABLE :: rtmp(:,:)
class(oft_bmesh), pointer :: smesh
namelist/cube_options/mesh_type,ni,rscale,shift,ref_per
DEBUG_STACK_PUSH
mesh_type=1
ni=1
rscale=1.d0
shift=0.d0
ref_per=.FALSE.
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,cube_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "cube_options" found in input file.','smesh_square_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "cube_options" in input file.','smesh_square_load',__FILE__)
  ni(3)=1
  WRITE(*,'(2A)')oft_indent,'Square surface mesh:'
  CALL oft_increase_indent
  WRITE(*,'(2A,I4)')oft_indent,     'Mesh Type   = ',mesh_type
  WRITE(*,'(2A,2I4)')oft_indent,    'nx, ny      = ',ni(1:2)
  WRITE(*,'(2A,2ES11.3)')oft_indent,'Scale facs  = ',rscale(1:2)
  WRITE(*,'(2A,2L2)')oft_indent,    'Periodicity = ',ref_per(1:2)
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(mesh_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_square_load',__FILE__)
CALL MPI_Bcast(ni,3,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_square_load',__FILE__)
CALL MPI_Bcast(rscale,3,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_square_load',__FILE__)
CALL MPI_Bcast(ref_per,3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_square_load',__FILE__)
#endif
!--Allocate mesh
IF(mesh_type==1)THEN
  allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
ELSE
  allocate(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
END IF
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%smeshes(i)%setup(mesh_cube_id,.FALSE.)
END DO
CALL multigrid_level(mg_mesh,1)
smesh=>mg_mesh%smeshes(1)
smesh%dim=2
!---Create mesh as quads
IF(oft_env%rank==0)THEN
  !---Setup points
  nptmp=PRODUCT(ni+1)
  allocate(rtmp(2,nptmp))
  ALLOCATE(pmap(ni(1)+1,ni(2)+1))
  pmap=-1
  nptmp=0
  DO i=1,ni(1)+1
    DO j=1,ni(2)+1
      nptmp=nptmp+1
      pmap(i,j)=nptmp
      rtmp(:,nptmp)=[(i-1)*rscale(1)/REAL(ni(1),8), &
        (j-1)*rscale(2)/REAL(ni(2),8)]+shift(1:2)
    END DO
  END DO
  !---Setup cells
  nctmp=PRODUCT(ni(1:2))
  allocate(lctmp(4,nctmp))
  nctmp=0
  DO i=1,ni(1)
    DO j=1,ni(2)
      nctmp=nctmp+1
      lctmp(:,nctmp)=[pmap(i,j),pmap(i+1,j),pmap(i+1,j+1),pmap(i,j+1)]
    END DO
  END DO
END IF
!--Build mesh
IF(mesh_type==1)THEN
  IF(oft_env%rank==0)THEN
    !---Setup points
    smesh%np=nptmp+nctmp
    allocate(smesh%r(3,smesh%np))
    smesh%r=0.d0
    DO i=1,nptmp
      smesh%r(1:2,i)=rtmp(:,i)
    END DO
    DO i=1,nctmp
      DO j=1,4
        smesh%r(1:2,nptmp+i)=smesh%r(1:2,nptmp+i)+rtmp(:,lctmp(j,i))/4.d0
      END DO
    END DO
    !---Setup cells
    smesh%nc=4*nctmp
    allocate(smesh%lc(3,smesh%nc),smesh%reg(smesh%nc))
    smesh%reg=1
    DO i=1,nctmp
      smesh%lc(:,(i-1)*4+1)=[lctmp(1,i),lctmp(2,i),nptmp+i]
      smesh%lc(:,(i-1)*4+2)=[lctmp(2,i),lctmp(3,i),nptmp+i]
      smesh%lc(:,(i-1)*4+3)=[lctmp(3,i),lctmp(4,i),nptmp+i]
      smesh%lc(:,(i-1)*4+4)=[lctmp(4,i),lctmp(1,i),nptmp+i]
    END DO
  END IF
ELSE
  IF(oft_env%rank==0)THEN
    !---Copy points
    smesh%np=nptmp
    allocate(smesh%r(3,smesh%np))
    smesh%r=0.d0
    DO i=1,nptmp
      smesh%r(1:2,i)=rtmp(:,i)
    END DO
    !---Setup cells
    smesh%nc=nctmp
    allocate(smesh%lc(4,smesh%nc),smesh%reg(smesh%nc))
    smesh%reg=1
    DO i=1,nctmp
      smesh%lc(:,i)=lctmp(:,i)
    END DO
  END IF
END IF
IF(oft_env%rank==0)DEALLOCATE(rtmp,lctmp)
call mesh_global_resolution(smesh)
DEBUG_STACK_POP
end subroutine smesh_square_load
!---------------------------------------------------------------------------------
!> Setup surface IDs
!---------------------------------------------------------------------------------
subroutine smesh_square_cadlink(smesh)
class(oft_bmesh), intent(inout) :: smesh
integer(i4) :: i,j
DO i=1,smesh%nbe
  j=smesh%lbe(i)
  IF(ALL(smesh%r(1,smesh%le(:,j))==1.d0))THEN
    smesh%bes(i)=1
  ELSE IF(ALL(smesh%r(1,smesh%le(:,j))==0.d0))THEN
    smesh%bes(i)=2
  ELSE IF(ALL(smesh%r(2,smesh%le(:,j))==1.d0))THEN
    smesh%bes(i)=3
  ELSE IF(ALL(smesh%r(2,smesh%le(:,j))==0.d0))THEN
    smesh%bes(i)=4
  END IF
END DO
end subroutine smesh_square_cadlink
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points
!! @note All edges are straight so construction is trivial
!---------------------------------------------------------------------------------
subroutine smesh_square_add_quad(smesh)
class(oft_bmesh), intent(inout) :: smesh
integer(i4) :: i,j
real(r8) :: pt(3)
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,*)'Setting square quadratic nodes'
!---Setup quadratic mesh
CALL smesh%set_order(2)
!---Locate edge end points and place daughter node
!$omp parallel do
do i=1,smesh%ne
  smesh%ho_info%lep(1,i)=i
  smesh%ho_info%r(:,i)=(smesh%r(:,smesh%le(1,i))+smesh%r(:,smesh%le(2,i)))/2.d0
enddo
IF(smesh%ho_info%ncp>0)THEN
  !---Locate cell vertices and place daughter node
  !$omp parallel do
  do i=1,smesh%nc
    smesh%ho_info%lcp(1,i)=i+smesh%ne
    pt=0.d0
    DO j=1,smesh%cell_np
      pt=pt+smesh%r(:,smesh%lc(j,i))
    END DO
    smesh%ho_info%r(:,i+smesh%ne)=pt/REAL(smesh%cell_np,8)
  enddo
END IF
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine smesh_square_add_quad
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine smesh_square_set_periodic(smesh)
class(oft_bmesh), intent(inout) :: smesh
integer(i4) :: i,j,jj,k,kk,l,m,n,iper
real(r8) :: pt_i(3),pt_j(3),d_plane,per_dir(3),i_cc(3),j_cc(3)
real(r8), parameter :: tol=1.d-6
logical :: flag(4)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Setting square periodicity'
!---Find periodic faces
smesh%periodic%nper=COUNT(ref_per(1:2))
ALLOCATE(smesh%periodic%le(smesh%ne))
smesh%periodic%le=-1
DO iper=1,2
  IF(.NOT.ref_per(iper))CYCLE
  DO jj=1,smesh%nbe
    j=smesh%lbe(jj)
    IF(smesh%periodic%le(j)>0)CYCLE
    pt_i=[-smesh%r(2,smesh%le(2,j))+smesh%r(2,smesh%le(1,j)), &
        smesh%r(1,smesh%le(2,j))-smesh%r(1,smesh%le(1,j)),0.d0]
    IF(MAXLOC(ABS(pt_i),DIM=1)/=iper)CYCLE
    i_cc=0.d0
    DO l=1,2
      i_cc=i_cc+smesh%r(:,smesh%le(l,j))
    END DO
    i_cc=i_cc/2.d0
    i_cc(iper)=0.d0
    DO kk=1,smesh%nbe
      k=smesh%lbe(kk)
      IF(k<=j)CYCLE
      IF(smesh%periodic%le(k)>0)CYCLE
      pt_j=[-smesh%r(2,smesh%le(2,k))+smesh%r(2,smesh%le(1,k)), &
        smesh%r(1,smesh%le(2,k))-smesh%r(1,smesh%le(1,k)),0.d0]
      IF(MAXLOC(ABS(pt_j),DIM=1)==iper)THEN
        j_cc=0.d0
        DO l=1,2
          j_cc=j_cc+smesh%r(:,smesh%le(l,k))
        END DO
        j_cc=j_cc/2.d0
        j_cc(iper)=0.d0
        IF(magnitude(i_cc-j_cc)<1.d-3)THEN
          smesh%periodic%le(k)=j
          EXIT
        END IF
      END IF
    END DO
  END DO
END DO
!---Set periodic points
ALLOCATE(smesh%periodic%lp(smesh%np))
smesh%periodic%lp=-1
DO iper=1,2
  IF(.NOT.ref_per(iper))CYCLE
  per_dir=0.d0
  per_dir(iper)=1.d0
  !$omp parallel do private(j,i,l,k,n,m,pt_i,pt_j,d_plane,flag)
  DO jj=1,smesh%nbe
    j=smesh%lbe(jj)
    i=smesh%periodic%le(j)
    IF(i<=0)CYCLE
    !---Check points
    flag=.FALSE.
    DO n=1,2
      m=ABS(smesh%le(n,j))
      IF(smesh%periodic%lp(m)>0)CYCLE
      pt_j=smesh%r(:,m)
      pt_j=pt_j-DOT_PRODUCT(pt_j,per_dir)*per_dir
      DO l=1,2
        IF(flag(l))CYCLE
        k=ABS(smesh%le(l,i))
        pt_i=smesh%r(:,k)
        pt_i=pt_i-DOT_PRODUCT(pt_i,per_dir)*per_dir
        !---
        pt_i=(pt_i-pt_j)
        d_plane=DOT_PRODUCT(pt_i,pt_i)
        IF(d_plane<tol)THEN
          flag(l)=.TRUE.
          smesh%periodic%lp(m)=k
        END IF
      END DO
    END DO
  END DO
END DO
DEBUG_STACK_POP
end subroutine smesh_square_set_periodic
end module oft_mesh_cube
