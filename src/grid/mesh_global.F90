!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_global.F90
!
!> MPI constructs and subroutines for global operations
!!
!! MPI types for global mesh information
!! - Global Mesh indices
!! - Proc to proc linkage information
!!
!! Functions to intialize the MPI structure
!! - Preform mesh decomposition and scatter
!! - Construct MPI linkage information
!!
!! @author Chris Hansen
!! @date Spring 2010
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_global
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, oft_amesh, oft_bmesh
USE oft_mesh_local, ONLY: mesh_local_init, mesh_local_init, oft_metis_partmesh, &
  bmesh_local_init
USE oft_stitching, ONLY: oft_seam, oft_stitch_check
IMPLICIT NONE
#include "local.h"
contains
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_init
!------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!------------------------------------------------------------------------------
subroutine mesh_global_init(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i
DEBUG_STACK_PUSH
call mesh_global_sync(mesh) ! Sync mesh information
call mesh_local_init(mesh)
!---
allocate(mesh%global%lp(mesh%np))
!$omp parallel do
do i=1,mesh%np
  mesh%global%lp(i)=i
end do
!---
allocate(mesh%global%le(mesh%ne))
!$omp parallel do
do i=1,mesh%ne
  mesh%global%le(i)=i
end do
!---
allocate(mesh%global%lf(mesh%nf))
!$omp parallel do
do i=1,mesh%nf
  mesh%global%lf(i)=i
end do
!---
ALLOCATE(mesh%global%lc(mesh%nc))
!$omp parallel do
do i=1,mesh%nc
  mesh%global%lc(i)=i
end do
!---
mesh%global%np=mesh%np
mesh%global%ne=mesh%ne
mesh%global%nf=mesh%nf
mesh%global%nc=mesh%nc
ALLOCATE(mesh%global%gbp(mesh%np))
ALLOCATE(mesh%global%gbe(mesh%ne))
ALLOCATE(mesh%global%gbf(mesh%nf))
ALLOCATE(mesh%global%gbc(mesh%nc))
mesh%global%gbp=mesh%bp
mesh%global%gbe=mesh%be
mesh%global%gbf=mesh%bf
mesh%global%gbc=mesh%bc
mesh%save%nbf=mesh%nbf
DEBUG_STACK_POP
end subroutine mesh_global_init
!------------------------------------------------------------------------------
! SUBROUTINE: tetmesh_global_init
!------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!------------------------------------------------------------------------------
subroutine bmesh_global_init(mesh)
class(oft_mesh), intent(inout) :: mesh
INTEGER(i4) :: i,j,np,nc
INTEGER(i4), ALLOCATABLE :: lptmp(:)
DEBUG_STACK_PUSH
!---Create surface mesh from boundary
nc=0
DO i=1,mesh%nbf
  j=mesh%lbf(i)
  IF(mesh%global%gbf(j))nc=nc+1
END DO
mesh%bmesh%nc=nc
!---
np=0
DO i=1,mesh%nbp
  j=mesh%lbp(i)
  IF(mesh%global%gbp(j))np=np+1
END DO
IF(np==0)mesh%bmesh%skip=.TRUE.
mesh%bmesh%np=np
ALLOCATE(mesh%bmesh%r(3,np))
ALLOCATE(mesh%bmesh%parent%lp(np))
ALLOCATE(lptmp(mesh%np))
lptmp=0
np=0
DO i=1,mesh%nbp
  j=mesh%lbp(i)
  IF(mesh%global%gbp(j))THEN
    np=np+1
    mesh%bmesh%r(:,np)=mesh%r(:,j)
    mesh%bmesh%parent%lp(np)=j
    lptmp(j)=np
  END IF
END DO
!---
ALLOCATE(mesh%bmesh%lc(mesh%bmesh%cell_np,nc))
ALLOCATE(mesh%bmesh%reg(nc))
ALLOCATE(mesh%bmesh%parent%lf(nc))
nc=0
DO i=1,mesh%nbf
  j=mesh%lbf(i)
  IF(mesh%global%gbf(j))THEN
    nc=nc+1
    mesh%bmesh%lc(:,nc)=lptmp(mesh%lf(:,j))
    mesh%bmesh%parent%lf(nc)=j
    mesh%bmesh%reg(nc)=mesh%reg(mesh%lfc(1,j))
  END IF
END DO
DEALLOCATE(lptmp)
!---
CALL bmesh_local_init(mesh%bmesh,mesh)
!---Copy high-order mapping if available
IF(ASSOCIATED(mesh%ho_info%r))THEN
  mesh%bmesh%order=1
  mesh%bmesh%ho_info%nep=mesh%ho_info%nep
  mesh%bmesh%ho_info%ncp=mesh%ho_info%nfp
  ALLOCATE(mesh%bmesh%ho_info%r(3,mesh%bmesh%ne+mesh%bmesh%nc))
  ALLOCATE(mesh%bmesh%ho_info%lep(mesh%bmesh%ho_info%nep,mesh%bmesh%ne))
  ALLOCATE(mesh%bmesh%ho_info%lcp(mesh%bmesh%ho_info%ncp,mesh%bmesh%nc))
  !!$omp parallel do private(j)
  do i=1,mesh%bmesh%ne
    mesh%bmesh%ho_info%lep(1,i)=i
    j=INT(ABS(mesh%bmesh%parent%le(i)),4)
    mesh%bmesh%ho_info%r(:,i)=mesh%ho_info%r(:,mesh%ho_info%lep(1,j))
  enddo
  IF(mesh%bmesh%ho_info%ncp==1)THEN
    !!$omp parallel do private(j)
    do i=1,mesh%bmesh%nc
      j=INT(ABS(mesh%bmesh%parent%lf(i)),4)
      mesh%bmesh%ho_info%lcp(1,i)=i+mesh%bmesh%ne
      mesh%bmesh%ho_info%r(:,i+mesh%bmesh%ne)=mesh%ho_info%r(:,mesh%ho_info%lfp(1,j))
    enddo
  END IF
END IF
DEBUG_STACK_POP
end subroutine bmesh_global_init
!------------------------------------------------------------------------------
! SUBROUTINE: smesh_global_init
!------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!------------------------------------------------------------------------------
subroutine smesh_global_init(mesh)
class(oft_bmesh), intent(inout) :: mesh
integer(i4) :: i
DEBUG_STACK_PUSH
call mesh_global_sync(mesh) ! Sync mesh information
call bmesh_local_init(mesh)
!---
mesh%global%np=mesh%np
ALLOCATE(mesh%global%gbp(mesh%np))
mesh%global%gbp=mesh%bp
allocate(mesh%global%lp(mesh%np))
!$omp parallel do
do i=1,mesh%np
  mesh%global%lp(i)=i
end do
!---
mesh%global%ne=mesh%ne
ALLOCATE(mesh%global%gbe(mesh%ne))
mesh%global%gbe=mesh%be
allocate(mesh%global%le(mesh%ne))
!$omp parallel do
do i=1,mesh%ne
  mesh%global%le(i)=i
end do
!---
mesh%global%nc=mesh%nc
ALLOCATE(mesh%global%gbc(mesh%nc))
mesh%global%gbc=mesh%bc
ALLOCATE(mesh%global%lc(mesh%nc))
!$omp parallel do
do i=1,mesh%nc
  mesh%global%lc(i)=i
end do
!---
ALLOCATE(mesh%lco(mesh%nc))
mesh%lco=1
ALLOCATE(mesh%pstitch%leo(mesh%nbp))
ALLOCATE(mesh%estitch%leo(mesh%nbe))
mesh%pstitch%leo=.TRUE.
mesh%estitch%leo=.TRUE.
DEBUG_STACK_POP
end subroutine smesh_global_init
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_link
!------------------------------------------------------------------------------
!> Driver for global seam linkage construction
!!
!! Construct inter-processor linkage information
!! - Link seam points, edges and faces
!------------------------------------------------------------------------------
subroutine mesh_global_link(mesh)
class(oft_mesh), intent(inout) :: mesh
DEBUG_STACK_PUSH
IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'Generating domain linkage'
CALL oft_increase_indent
call mesh_global_plinkage(mesh) ! Link seam points
call mesh_global_elinkage(mesh) ! Link seam edges
call mesh_global_flinkage(mesh) ! Link seam faces
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_global_link
!------------------------------------------------------------------------------
! SUBROUTINE: bmesh_global_link
!------------------------------------------------------------------------------
!> Driver for global seam linkage construction
!!
!! Construct inter-processor linkage information
!! - Link seam points, edges and faces
!------------------------------------------------------------------------------
subroutine bmesh_global_link(mesh,parent)
class(oft_bmesh), intent(inout) :: mesh
class(oft_mesh), optional, intent(inout) :: parent
DEBUG_STACK_PUSH
IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'Generating boundary domain linkage'
CALL oft_increase_indent
IF(PRESENT(parent))THEN
  call plinkage_from_parent(mesh,parent) ! Link seam points
  call elinkage_from_parent(mesh,parent) ! Link seam edges
ELSE
  call mesh_global_plinkage(mesh)
  call mesh_global_elinkage(mesh)
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine bmesh_global_link
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_sync
!------------------------------------------------------------------------------
!> Scatters base mesh information to all processors.
!!
!! Scatters mesh information from head task to all other tasks using MPI_BCAST calls.
!! - Communicates base mesh information (np,nc,lc,r)
!------------------------------------------------------------------------------
subroutine mesh_global_sync(mesh)
class(oft_amesh), intent(inout) :: mesh
#ifdef HAVE_MPI
integer(i4) :: ierr
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Syncing mesh'
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(mesh%cad_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast mesh type
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
!---
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(mesh%nc,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell count
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(mesh%np,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast point count
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(mesh%meshname,20,OFT_MPI_CHAR,0,oft_env%COMM,ierr) ! Broadcast mesh name
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
!---
if(oft_env%rank>0)allocate(mesh%r(3,mesh%np),mesh%lc(mesh%cell_np,mesh%nc),mesh%reg(mesh%nc)) ! Allocate point and cell arrays
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(mesh%lc,mesh%cell_np*mesh%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(mesh%r,3*mesh%np,OFT_MPI_R8,0,oft_env%COMM,ierr) ! Broadcast point list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
call MPI_Bcast(mesh%reg,mesh%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast region list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
DEBUG_STACK_POP
#endif
end subroutine mesh_global_sync
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_partition
!------------------------------------------------------------------------------
!> Perform mesh decomposition (METIS) and local construction.
!! - Decompose domain of head task using METIS library.
!! - Scatter decomposition to all tasks.
!------------------------------------------------------------------------------
subroutine mesh_global_partition(mesh,meshpart,part_meth)
class(oft_amesh), intent(inout) :: mesh
integer(i4), intent(inout) :: meshpart(:)
integer(i4), intent(in) :: part_meth
integer(i4) :: i,k,m,nctmp,nptmp,ierr
integer(i4), allocatable, dimension(:) :: lptmp,lctmp,pflag,ncells,cpart,isort
integer(i4), allocatable, dimension(:,:) :: lcctmp
real(r8) :: pt(3),zmin,zmax,zrange
real(r8), allocatable, dimension(:) :: part_sort
DEBUG_STACK_PUSH
!---Synchronize before deomposition
CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Decompose domain using METIS
if(oft_env%rank==0)then
  !---
  allocate(cpart(mesh%nc))
  cpart=0
  if(part_meth==1)THEN ! Metis partitioning
    CALL oft_metis_partMesh(mesh%nc,mesh%np,mesh%cell_np,mesh%lc,oft_env%nnodes,cpart,ierr)
    IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_global_partition',__FILE__)
  elseif((part_meth>1).AND.(part_meth<5))THEN ! Axial spatial partitioning
    zmin = MINVAL(mesh%r(part_meth-1,:))
    zmax = MAXVAL(mesh%r(part_meth-1,:))
    zrange = (zmax-zmin)/REAL(oft_env%nnodes,8)
    ! IF(MOD(mesh%nc,oft_env%nnodes)/=0)CALL oft_abort( &
    !   "Partition count invalid for spatial partitioning",'mesh_global_partition',__FILE__)
    ALLOCATE(pflag(mesh%np),part_sort(mesh%np),isort(mesh%np))
    pflag=oft_env%nnodes+1
    isort=[(i,i=1,mesh%np)]
    part_sort=mesh%r(part_meth-1,:)
    CALL sort_array(part_sort,isort,mesh%np)
    k=1
    DO i=1,mesh%np
      IF(i>REAL(k*mesh%np,8)/REAL(oft_env%nnodes,8))THEN
        k=k+1
      END IF
      pflag(isort(i))=k
    END DO
    WRITE(*,*)'CHK',COUNT(pflag<=0),COUNT(pflag>oft_env%nnodes)
    DO i=1,mesh%nc
      cpart(i)=oft_env%nnodes+1
      DO k=1,mesh%cell_np
        cpart(i)=MIN(cpart(i),pflag(mesh%lc(k,i)))
      END DO
      IF(cpart(i)<1.OR.cpart(i)>oft_env%nnodes)WRITE(*,*)'cBAD',i,cpart(i)
    END DO
  elseif(part_meth==5)THEN ! Toroidal spatial partitioning
    zmin = -pi
    zmax = pi
    zrange = (zmax-zmin)/REAL(oft_env%nnodes,8)
    ! IF(MOD(mesh%nc,oft_env%nnodes)/=0)CALL oft_abort( &
    ! 'Partition count invalid for spatial partitioning','mesh_global_partition',__FILE__)
    ALLOCATE(pflag(mesh%np),part_sort(mesh%np),isort(mesh%np))
    pflag=oft_env%nnodes+1
    isort=[(i,i=1,mesh%np)]
    part_sort=ATAN2(mesh%r(2,:),mesh%r(1,:))
    CALL sort_array(part_sort,isort,mesh%np)
    k=1
    DO i=1,mesh%np
      IF(i>REAL(k*mesh%np,8)/REAL(oft_env%nnodes,8))THEN
        k=k+1
      END IF
      pflag(isort(i))=k
    END DO
    WRITE(*,*)'CHK',COUNT(pflag<=0),COUNT(pflag>oft_env%nnodes)
    DO i=1,mesh%nc
      cpart(i)=oft_env%nnodes+1
      DO k=1,mesh%cell_np
        cpart(i)=MIN(cpart(i),pflag(mesh%lc(k,i)))
      END DO
      IF(cpart(i)<1.OR.cpart(i)>oft_env%nnodes)WRITE(*,*)'cBAD',i,cpart(i)
    END DO
  else
    CALL oft_abort('Invalid partitioning method', 'mesh_global_partition', __FILE__)
  endif
  ! CALL oft_metis_partMesh(mesh%nc,mesh%np,mesh%cell_np,mesh%lc,oft_env%nnodes,cpart,ierr)
  ! IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_global_partition',__FILE__)
  meshpart=cpart
  deallocate(cpart)
  if(oft_env%ppn>1)then
    !---
    allocate(lptmp(mesh%np),lctmp(mesh%nc),lcctmp(mesh%cell_np,mesh%nc))
    allocate(pflag(mesh%np))
    do m=1,oft_env%nnodes
      !---Initialize counters and flags
      pflag=0
      lctmp=0
      lptmp=0
      lcctmp=0
      nctmp=0
      nptmp=0
      do i=1,mesh%nc ! Loop over global cells
        if(meshpart(i)==m)then ! Processor owns this cell
          nctmp = nctmp + 1 ! Increment cell counter
          lctmp(nctmp) = i ! Link cell to global index
          do k=1,mesh%cell_np ! Loop over points
            if(pflag(mesh%lc(k,i))==0)then ! Point is not claimed yet
              nptmp=nptmp+1 ! Increment point counter
              pflag(mesh%lc(k,i))=nptmp ! Link global point to new index
              lptmp(nptmp)=mesh%lc(k,i) ! Link point to global index
              lcctmp(k,nctmp)=nptmp ! Insert local point index into local cell list
            else ! Already have this point
              lcctmp(k,nctmp)=pflag(mesh%lc(k,i)) ! Insert local point index into local cell list
            endif
          enddo
        endif
      enddo
      !---
      allocate(cpart(nctmp))
      CALL oft_metis_partMesh(nctmp,nptmp,mesh%cell_np,lcctmp,oft_env%ppn,cpart,ierr)
      IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_global_partition',__FILE__)
      !$omp parallel do
      do i=1,nctmp
        meshpart(lctmp(i))=-((m-1)*oft_env%ppn+cpart(i))
      end do
      deallocate(cpart)
    end do
    !$omp parallel do
    do i=1,mesh%nc
      meshpart(i)=-meshpart(i)
    end do
  end if
  allocate(ncells(oft_env%nprocs))
  ncells=0
  do i=1,mesh%nc
    ncells(meshpart(i))=ncells(meshpart(i))+1
  end do
  write(*,'(2A)')oft_indent,'Mesh Partitioned:'
  CALL oft_increase_indent
  write(*,'(2A,I8)')oft_indent,'Min # of cells/domain = ',MINVAL(ncells)
  write(*,'(2A,I8)')oft_indent,'Max # of cells/domain = ',MAXVAL(ncells)
  CALL oft_decrease_indent
endif
!---Scatter decomposition to all processors
CALL oft_mpi_barrier(ierr) ! Wait for all processes
#ifdef HAVE_MPI
call MPI_Bcast(meshpart,mesh%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell partition information
#endif
DEBUG_STACK_POP
end subroutine mesh_global_partition
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_decomp
!------------------------------------------------------------------------------
!> Perform mesh decomposition and local construction.
!! - Decompose global mesh
!! - Perform mesh construction of local domain
!------------------------------------------------------------------------------
subroutine mesh_global_decomp(mesh,part_meth)
class(oft_amesh), intent(inout) :: mesh
integer(i4), intent(in) :: part_meth
integer(i4), allocatable :: lptmp(:),lctmp(:),lcctmp(:,:),isort(:),isorti(:)
integer(i4), allocatable :: pflag(:),conflag(:),proc_contmp(:)
integer(i4) :: i,k,npcors
real(r8), allocatable :: rtmp(:,:)
DEBUG_STACK_PUSH
!---Synchronize before deomposition
allocate(mesh%base%lcpart(mesh%nc))
call mesh_global_partition(mesh,mesh%base%lcpart,part_meth)
!---
npcors=mesh%np ! Get coarse point count
!---Create local mesh indexing for points and cells
allocate(lptmp(mesh%np),lctmp(mesh%nc),lcctmp(mesh%cell_np,mesh%nc),rtmp(3,mesh%np))
allocate(pflag(mesh%np),conflag(npcors),proc_contmp(oft_env%nprocs))
!---Initialize counters and flags
pflag=0
conflag=0
proc_contmp=0
lctmp=0
lptmp=0
lcctmp=0
rtmp=0.d0
mesh%nc=0
mesh%np=0
do i=1,INT(mesh%global%nc,4) ! Loop over global cells
  if(mesh%base%lcpart(i)==(oft_env%rank+1))then ! Processor owns this cell
    mesh%nc = mesh%nc + 1 ! Increment cell counter
    lctmp(mesh%nc) = i ! Link cell to global index
    do k=1,mesh%cell_np ! Loop over points
      if(pflag(mesh%lc(k,i))==0)then ! Point is not claimed yet
        mesh%np=mesh%np+1 ! Increment point counter
        pflag(mesh%lc(k,i))=mesh%np ! Link global point to new index
        conflag(mesh%global%lp(mesh%lc(k,i)))=mesh%np ! Mark global point list
        lptmp(mesh%np)=mesh%lc(k,i) ! Link point to global index
        rtmp(:,mesh%np)=mesh%r(:,mesh%lc(k,i)) ! Insert point into local list
        lcctmp(k,mesh%nc)=mesh%np ! Insert local point index into local cell list
      else ! Already have this point
        lcctmp(k,mesh%nc)=pflag(mesh%lc(k,i)) ! Insert local point index into local cell list
      endif
    enddo
  endif
enddo
call mesh_global_proccon(conflag,npcors) ! Determine neighbor processors
allocate(oft_env%send(oft_env%nproc_con),oft_env%recv(oft_env%nproc_con))
!---Sort point list
allocate(isort(mesh%np))
isort=(/(i,i=1,mesh%np)/)
CALl sort_array(lptmp, isort, mesh%np)
allocate(isorti(mesh%np))
!---Allocate local mesh arrays
allocate(mesh%r(3,mesh%np),mesh%lc(mesh%cell_np,mesh%nc))
DO i=1,mesh%np
  isorti(isort(i))=i
  mesh%r(:,i)=rtmp(:,isort(i))
END DO
DO i=1,mesh%nc
  mesh%lc(:,i)=isorti(lcctmp(:,i))
END DO
! mesh%r=rtmp(:,1:mesh%np) ! Set local point list
! mesh%lc=lcctmp(:,1:mesh%nc) ! Set local cell list
DEALLOCATE(isort,isorti)
!---Get indexing on base mesh
allocate(mesh%base%lp(mesh%np),mesh%base%lc(mesh%nc))
!$omp parallel do
do i=1,mesh%np
  mesh%base%lp(i) = lptmp(i) ! Set base point index
  lptmp(i) = INT(mesh%global%lp(lptmp(i)),4) ! Set global point index
end do
!$omp parallel do
do i=1,mesh%nc
  mesh%base%lc(i) = lctmp(i) ! Set base cell index
  lctmp(i) = INT(mesh%global%lc(lctmp(i)),4) ! Set global cell index
end do
!---Assign global indices
allocate(mesh%global%lp(mesh%np),mesh%global%lc(mesh%nc))
mesh%global%lp = INT(lptmp(1:mesh%np),8) ! Set global point index
mesh%global%lc = INT(lctmp(1:mesh%nc),8) ! Set global cell index
!---Kill temporary work arrays
deallocate(lptmp,lctmp,lcctmp,rtmp)
deallocate(pflag,conflag)
DEBUG_STACK_POP
end subroutine mesh_global_decomp
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_proccon
!------------------------------------------------------------------------------
!> Determine neighboring processors for MPI linkage.
!!
!! @param[in] lptmp List of global point indices from initialization
!------------------------------------------------------------------------------
subroutine mesh_global_proccon(lptmp,np)
integer(i4), intent(in) :: lptmp(:)
integer(i4), intent(in) :: np
integer(i4), allocatable :: c(:),b(:),a(:)
integer(i4) :: i,j,k,l,ierr
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
type(mpi_request) :: req
type(mpi_status) :: stat
#else
integer(i4) :: req,stat(MPI_STATUS_SIZE)
#endif
#endif
DEBUG_STACK_PUSH
!---Initialize counters for ring update
j=oft_env%rank+1
if(j>oft_env%nprocs-1)j = 0
k=oft_env%rank - 1
if(k<0)k=oft_env%nprocs-1
l=k
!
allocate(c(np),a(np),b(oft_env%nprocs)) ! Allocate temporary work arrays
b=0
a=lptmp
CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Loop over processors to find shared points
#ifdef HAVE_MPI
do while(l.NE.oft_env%rank)
  CALL MPI_ISEND(a,np,OFT_MPI_I4,j,1,oft_env%COMM,req,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_proccon',__FILE__)
  CALL MPI_RECV(c,np,OFT_MPI_I4,k,1,oft_env%COMM,stat,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_RECV','mesh_global_proccon',__FILE__)
  CALL MPI_WAIT(req,stat,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAIT','mesh_global_proccon',__FILE__)
  a=c
  !---Test if point is present on processor
  do i=1,np
    if(c(i)>0.AND.lptmp(i)>0)then
      b(l+1)=1 ! Processor shares at least 1 point
      exit
    endif
  enddo
  !---Update processor index
  l=l-1
  if(l<0)l=oft_env%nprocs-1
enddo
#endif
DEALLOCATE(a,c)
oft_env%nproc_con=sum(b) ! Number of processor connections
oft_env%proc_split=oft_env%nproc_con
allocate(oft_env%proc_con(oft_env%nproc_con)) ! Allocate processor linkage list
k=1
do i=1,oft_env%nprocs
  if(b(i)==1)then ! Processor links to this domain
    oft_env%proc_con(k)=i-1 ! Populate processor linkage list
    IF((i-1>oft_env%rank).AND.(k-1<oft_env%proc_split))oft_env%proc_split=k-1
    k=k+1
  endif
enddo
DEALLOCATE(b)
DEBUG_STACK_POP
end subroutine mesh_global_proccon
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_plinkage
!------------------------------------------------------------------------------
!> Construct processor to processor point linkage for stitching operations.
!! - Create linkage of boundary points to other processors
!! - Determine ownership for shared points
!------------------------------------------------------------------------------
SUBROUTINE mesh_global_plinkage(mesh)
class(oft_amesh), intent(inout) :: mesh
TYPE :: pout
  INTEGER(8), POINTER, DIMENSION(:) :: lp => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: rp => NULL()
END TYPE pout
TYPE(pout), ALLOCATABLE, DIMENSION(:) :: lprecv,lpsend
INTEGER(8), POINTER, DIMENSION(:) :: lptmp,lpout
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: lsort
INTEGER(i4), ALLOCATABLE :: linktmp(:,:,:),isort(:),ncon(:),bpi(:),child_list(:,:)
INTEGER(i4) :: nppl
INTEGER(i4) :: m,mm,nptmp,nbpmax
INTEGER(i4) :: ierr,i,j
! TYPE(oft_seam) :: pstitch
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Point linkage'
CALL oft_increase_indent
IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---------------------------------------------------------------------------
! Determine maximum boundary point count
!---------------------------------------------------------------------------
nptmp=mesh%nbp
nbpmax=mesh%nbp
IF(.NOT.mesh%fullmesh)nbpmax=oft_mpi_max(nptmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam points =',nbpmax
mesh%pstitch%nbemax=nbpmax
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,(mesh%periodic%nper+2)*mesh%nbp,0:oft_env%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:oft_env%nproc_con))
ncon=0
!---Pointer into linkage array
ALLOCATE(mesh%pstitch%kle(0:oft_env%nproc_con+1))
mesh%pstitch%kle = 0
!---Global boundary point tag
IF(ASSOCIATED(mesh%global%gbp))DEALLOCATE(mesh%global%gbp)
ALLOCATE(mesh%global%gbp(mesh%np))
mesh%global%gbp = .FALSE.
!---Allocate temporary Send/Recv arrays
ALLOCATE(lpsend(oft_env%nproc_con+1),lprecv(oft_env%nproc_con+1))
ALLOCATE(lpsend(1)%lp(nbpmax))
ALLOCATE(lpsend(1)%rp(nbpmax))
IF(.NOT.mesh%fullmesh)THEN
  DO i=1,oft_env%nproc_con
    ALLOCATE(lprecv(i)%lp(nbpmax))
    ALLOCATE(lprecv(i)%rp(nbpmax))
  END DO
END IF
!---Count point to point interactions for load balancing
CALL oft_random_number(lpsend(1)%rp,nbpmax)
DO i=1,mesh%nbp
  lpsend(1)%rp(i)=lpsend(1)%rp(i)+REAL(mesh%kpp(mesh%lbp(i)+1)-mesh%kpp(mesh%lbp(i)),8)
END DO
!---Point dummy output array to Send array
lpout=>lpsend(1)%lp
lpout=0
!---Populate with global point index
!$omp parallel do
DO i=1,mesh%nbp
  lpout(i) = mesh%global%lp(mesh%lbp(i))
END DO
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  !---Point dummy Send arrays to main Send array
  DO i=2,oft_env%nproc_con
    lpsend(i)%lp=>lpsend(1)%lp
    lpsend(i)%rp=>lpsend(1)%rp
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!---------------------------------------------------------------------------
! Create Send and Recv calls
!---------------------------------------------------------------------------
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lpsend(j)%lp,nbpmax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%send(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_plinkage',__FILE__)
    CALL MPI_IRECV(lprecv(j)%lp,nbpmax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%recv(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_plinkage',__FILE__)
  END DO
!---------------------------------------------------------------------------
! Loop over each connected processor
!---------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    !---All recieves have been processed
    IF(oft_mpi_check_reqs(oft_env%nproc_con,oft_env%recv))EXIT
    !---Wait for completed recieve
    CALL oft_mpi_waitany(oft_env%nproc_con,oft_env%recv,j,ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_plinkage',__FILE__)
    !---Point dummy input array to current Recv array
    lptmp=>lprecv(j)%lp
    !---Determine location of boundary points on other processors
    nppl=0
    !$omp parallel do private(mm) reduction(+:nppl)
    DO m=1,mesh%nbp ! Loop over boundary points
      DO mm=1,nbpmax ! Loop over input points
        IF(lpout(m)==lptmp(mm))THEN ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          nppl=nppl+1
        END IF
        IF(lptmp(mm)==0)EXIT ! End of input points
      END DO
    END DO
    mesh%pstitch%kle(j)=nppl
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_plinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Check internal connections
!---------------------------------------------------------------------------
ALLOCATE(mesh%pstitch%leo(mesh%nbp))
mesh%pstitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic points on current processor
nppl=0
IF(mesh%periodic%nper>0)THEN
  IF(ASSOCIATED(mesh%periodic%lp))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bpi(mesh%np))
    CALL get_inverse_map(mesh%lbp,mesh%nbp,bpi,mesh%np)
    !---Construct child point list
    ALLOCATE(child_list(8,mesh%nbp))
    child_list=0
    DO m=1,mesh%nbp
      mm=mesh%lbp(m)
      DO
        IF(mesh%periodic%lp(mm)<=0)EXIT
        mm=mesh%periodic%lp(mm)
      END DO
      IF(mm/=mesh%lbp(m))THEN
        mm=bpi(mm)
        DO i=1,8
          IF(child_list(i,mm)==0)THEN
            child_list(i,mm)=m
            EXIT
          END IF
        END DO
      END IF
    END DO
    !---Find point pairs by linking child points
    DO m=1,mesh%nbp
      DO i=1,8
        !---Link parent point to child
        IF(child_list(i,m)>0)THEN
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/child_list(i,m),m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,child_list(i,m)/)
          mesh%pstitch%leo(child_list(i,m))=.FALSE.
          lpsend(1)%rp(child_list(i,m))=-1.d0
          nppl=nppl+2
          !---Link other child points to eachother
          DO j=i+1,8
            IF(child_list(j,m)>0)THEN
              ncon(0)=ncon(0)+1
              linktmp(:,ncon(0),0)=(/child_list(j,m),child_list(i,m)/)
              ncon(0)=ncon(0)+1
              linktmp(:,ncon(0),0)=(/child_list(i,m),child_list(j,m)/)
              nppl=nppl+2
            END IF
          END DO
        END IF
      END DO
    END DO
    DEALLOCATE(bpi,child_list)
  ELSE
    !$omp parallel do private(mm) reduction(+:nppl)
    DO m=1,mesh%nbp ! Loop over boundary points
      DO mm=m+1,mesh%nbp ! Loop over boundary points
        IF(lpout(m)==lpout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          nppl=nppl+2
          mesh%pstitch%leo(mm) = .FALSE.
          lpsend(1)%rp(mm)=-1.d0
        END IF
        IF(lpout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
mesh%pstitch%kle(0)=nppl
!---Transfer load balancing arrays
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lpsend(j)%rp,nbpmax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%send(j),ierr)
    CALL MPI_IRECV(lprecv(j)%rp,nbpmax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%recv(j),ierr)
  END DO
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%recv,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_plinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Condense linkage to sparse rep
!---------------------------------------------------------------------------
mesh%pstitch%nle=SUM(mesh%pstitch%kle)
mesh%pstitch%kle(oft_env%nproc_con+1)=mesh%pstitch%nle+1
!---Cumulative unique point linkage count
DO i=oft_env%nproc_con,0,-1
  mesh%pstitch%kle(i)=mesh%pstitch%kle(i+1)-mesh%pstitch%kle(i)
END DO
IF(mesh%pstitch%kle(0)/=1)CALL oft_abort('Bad point linkage count','mesh_global_plinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(mesh%pstitch%lle(2,mesh%pstitch%nle))
!---Create seam object
mesh%pstitch%full=mesh%fullmesh
mesh%pstitch%nbe=mesh%nbp
mesh%pstitch%lbe=>mesh%lbp
ALLOCATE(mesh%pstitch%send(0:oft_env%nproc_con),mesh%pstitch%recv(0:oft_env%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,oft_env%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(oft_env%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lprecv(i)%rp(linktmp(1,j,i))>lpsend(1)%rp(m))THEN
        mesh%pstitch%leo(m) = .FALSE.
      ELSE IF(lprecv(i)%rp(linktmp(1,j,i))==lpsend(1)%rp(m).AND.oft_env%proc_con(i)<oft_env%rank)THEN
        mesh%pstitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    mesh%pstitch%lle(:,m+mesh%pstitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Create preallocated send/recv arrays
  mesh%pstitch%send(i)%n=ncon(i)
  ALLOCATE(mesh%pstitch%send(i)%v(mesh%pstitch%send(i)%n))
  mesh%pstitch%recv(i)%n=ncon(i)
  ALLOCATE(mesh%pstitch%recv(i)%v(mesh%pstitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!---------------------------------------------------------------------------
! Clean-up transfers
!---------------------------------------------------------------------------
DEALLOCATE(lpsend(1)%lp)
DEALLOCATE(lpsend(1)%rp)
IF(.NOT.mesh%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,oft_env%nproc_con
    DEALLOCATE(lprecv(i)%lp)
    DEALLOCATE(lprecv(i)%rp)
  END DO
END IF
!---Check stitching information
CALL oft_stitch_check(mesh%pstitch)
DEALLOCATE(linktmp,ncon)
DEALLOCATE(lpsend,lprecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_plinkage
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_elinkage
!------------------------------------------------------------------------------
!> Construct processor to processor edge linkage for stitching operations.
!! - Create linkage of boundary edges to other processors with orientation.
!! - Determine ownership for shared edges
!------------------------------------------------------------------------------
SUBROUTINE mesh_global_elinkage(mesh)
class(oft_amesh), intent(inout) :: mesh
TYPE :: eout
  INTEGER(8), POINTER, DIMENSION(:) :: le => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: re => NULL()
END TYPE eout
TYPE(eout), ALLOCATABLE, DIMENSION(:) :: lerecv,lesend
INTEGER(8), POINTER, DIMENSION(:) :: letmp,leout
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: lsort
INTEGER(i4), ALLOCATABLE :: bpi(:),bei(:),linktmp(:,:,:),isort(:),ncon(:),child_list(:,:)
INTEGER(i4) :: m,mm,netmp,nbemax,ll(2),neel
INTEGER(i4) :: ierr,i,j,js,jn
LOGICAL :: etest,set_gbe
LOGICAL, ALLOCATABLE :: echeck(:,:)
! TYPE(oft_seam) :: estitch
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Edge linkage'
CALL oft_increase_indent
IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---------------------------------------------------------------------------
! Determine maximum boundary edge count
!---------------------------------------------------------------------------
netmp=mesh%nbe
nbemax=mesh%nbe
IF(.NOT.mesh%fullmesh)nbemax=oft_mpi_max(netmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam edges =',nbemax
mesh%estitch%nbemax=nbemax
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,(mesh%periodic%nper+2)*mesh%nbe,0:oft_env%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:oft_env%nproc_con),echeck(0:oft_env%nproc_con,mesh%nbe))
ncon=0
echeck=.FALSE.
!---Pointer into linkage array
ALLOCATE(mesh%estitch%kle(0:oft_env%nproc_con+1))
mesh%estitch%kle = 0
!---Global boundary edge tag (setup for surface meshes only)
IF(ASSOCIATED(mesh%global%gbe))DEALLOCATE(mesh%global%gbe)
ALLOCATE(mesh%global%gbe(mesh%ne))
mesh%global%gbe = .FALSE.
set_gbe=.FALSE.
SELECT TYPE(mesh)
CLASS IS(oft_bmesh)
  IF(.NOT.ASSOCIATED(mesh%parent))set_gbe=.TRUE.
END SELECT
IF(set_gbe)THEN
  DO i=1,mesh%ne
    IF(mesh%be(i))mesh%global%gbe(i)=.TRUE.
  END DO
END IF
!---Allocate temporary Send/Recv arrays
ALLOCATE(lesend(oft_env%nproc_con+1),lerecv(oft_env%nproc_con+1))
ALLOCATE(lesend(1)%le(nbemax))
ALLOCATE(lesend(1)%re(nbemax))
IF(.NOT.mesh%fullmesh)THEN
  DO i=1,oft_env%nproc_con
    ALLOCATE(lerecv(i)%le(nbemax))
    ALLOCATE(lerecv(i)%re(nbemax))
  END DO
ENDIF
!---Count edge to edge interactions for load balancing
CALL oft_random_number(lesend(1)%re,nbemax)
DO i=1,mesh%nbe
  lesend(1)%re(i)=lesend(1)%re(i)+REAL(mesh%kee(mesh%lbe(i)+1)-mesh%kee(mesh%lbe(i)),8)
END DO
!---Point dummy output array to Send array
leout=>lesend(1)%le
leout=0
!---Create link from local index to boundary index
ALLOCATE(bpi(mesh%np))
CALL get_inverse_map(mesh%lbp,mesh%nbp,bpi,mesh%np)
!---
!$omp parallel do private(m,js,jn,ll,j,etest)
DO i=1,mesh%nbe
  DO m=1,2
    ll(m)=bpi(mesh%le(m,mesh%lbe(i))) ! Get endpoints in boundary index
  END DO
  leout(i)=ABS(mesh%global%le(mesh%lbe(i))) ! Populate linkage array
  DO m=0,oft_env%nproc_con ! Flag processors to be checked
    js=mesh%pstitch%kle(m)
    jn=mesh%pstitch%kle(m+1)-1
    etest=.TRUE.
    DO j=1,2
      etest=etest.AND.ANY(mesh%pstitch%lle(1,js:jn)==ll(j))
    END DO
    IF(etest)echeck(m,i)=.TRUE.
  END DO
END DO
DEALLOCATE(bpi)
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  !---Point dummy Send arrays to main Send array
  DO i=2,oft_env%nproc_con
    lesend(i)%le=>lesend(1)%le
    lesend(i)%re=>lesend(1)%re
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!---------------------------------------------------------------------------
! Create Send and Recv calls
!---------------------------------------------------------------------------
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lesend(j)%le,nbemax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%send(j), ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_elinkage',__FILE__)
    CALL MPI_IRECV(lerecv(j)%le,nbemax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%recv(j), ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_elinkage',__FILE__)
  END DO
!---------------------------------------------------------------------------
! Loop over each connected processor
!---------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    !---All recieves have been processed
    IF(oft_mpi_check_reqs(oft_env%nproc_con,oft_env%recv))EXIT
    !---Wait for completed recieve
    CALL oft_mpi_waitany(oft_env%nproc_con,oft_env%recv,j,ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_elinkage',__FILE__)
    !---Point dummy input array to current Recv array
    letmp=>lerecv(j)%le
    !---
    neel=0
    !$omp parallel do private(mm) reduction(+:neel)
    DO m=1,mesh%nbe ! Loop over boundary edges
      IF(.NOT.echeck(j,m))CYCLE
      DO mm=1,nbemax ! Loop over input edges
        IF(leout(m)==letmp(mm))THEN ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          IF(set_gbe)mesh%global%gbe(mesh%lbe(m))=.FALSE. ! Edge is not on global boundary
          neel=neel+1
        END IF
      END DO
    END DO
    mesh%estitch%kle(j)=neel
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_elinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Check internal connections
!---------------------------------------------------------------------------
ALLOCATE(mesh%estitch%leo(mesh%nbe))
mesh%estitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic edges on current processor
neel=0
IF(mesh%periodic%nper>0)THEN
  IF(ASSOCIATED(mesh%periodic%le))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bei(mesh%ne))
    CALL get_inverse_map(mesh%lbe,mesh%nbe,bei,mesh%ne)
    !---Construct child edge list
    ALLOCATE(child_list(4,mesh%nbe))
    child_list=0
    DO m=1,mesh%nbe
      mm=mesh%lbe(m)
      DO
        IF(mesh%periodic%le(mm)<=0)EXIT
        mm=mesh%periodic%le(mm)
      END DO
      IF(mm/=mesh%lbe(m))THEN
        mm=bei(mm)
        DO i=1,4
          IF(child_list(i,mm)==0)THEN
            child_list(i,mm)=m
            EXIT
          END IF
        END DO
      END IF
    END DO
    !---Find edge pairs by linking child edges
    DO m=1,mesh%nbe
      DO i=1,4
        !---Link parent edge to child
        IF(child_list(i,m)>0)THEN
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/child_list(i,m),m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,child_list(i,m)/)
          mesh%estitch%leo(child_list(i,m))=.FALSE.
          lesend(1)%re(child_list(i,m))=-1.d0
          neel=neel+2
          IF(set_gbe)THEN ! Edge is not on global boundary
            mesh%global%gbe(mm)=.FALSE.
            mesh%global%gbe(mesh%periodic%le(mm))=.FALSE.
          END IF
          !---Link other child edges to eachother
          DO j=i+1,4
            IF(child_list(j,m)>0)THEN
              ncon(0)=ncon(0)+1
              linktmp(:,ncon(0),0)=(/child_list(j,m),child_list(i,m)/)
              ncon(0)=ncon(0)+1
              linktmp(:,ncon(0),0)=(/child_list(i,m),child_list(j,m)/)
              neel=neel+2
            END IF
          END DO
        END IF
      END DO
    END DO
    DEALLOCATE(bei,child_list)
  ELSE
    !$omp parallel do private(mm) reduction(+:neel)
    DO m=1,mesh%nbe ! Loop over boundary points
      IF(.NOT.echeck(0,m))CYCLE
      DO mm=m+1,mesh%nbe ! Loop over boundary points
        IF(leout(m)==leout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          neel=neel+2
          IF(set_gbe)THEN ! Edge is not on global boundary
            mesh%global%gbe(mesh%lbe(m))=.FALSE. 
            mesh%global%gbe(mesh%lbe(mm))=.FALSE.
          END IF
          mesh%estitch%leo(mm) = .FALSE.
          lesend(1)%re(mm)=-1.d0
        END IF
        IF(leout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
mesh%estitch%kle(0)=neel
!---Transfer load balancing arrays
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lesend(j)%re,nbemax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%send(j),ierr)
    CALL MPI_IRECV(lerecv(j)%re,nbemax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%recv(j),ierr)
  END DO
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%recv,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_elinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Condense linkage to sparse rep
!---------------------------------------------------------------------------
mesh%estitch%nle=SUM(mesh%estitch%kle)
mesh%estitch%kle(oft_env%nproc_con+1)=mesh%estitch%nle+1
DO i=oft_env%nproc_con,0,-1 ! cumulative unique edge linkage count
  mesh%estitch%kle(i)=mesh%estitch%kle(i+1)-mesh%estitch%kle(i)
END DO
IF(mesh%estitch%kle(0)/=1)CALL oft_abort('Bad edge linkage count','mesh_global_elinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(mesh%estitch%lle(2,mesh%estitch%nle))
!---Create seam object
mesh%estitch%full=mesh%fullmesh
mesh%estitch%nbe=mesh%nbe
mesh%estitch%lbe=>mesh%lbe
ALLOCATE(mesh%estitch%send(0:oft_env%nproc_con),mesh%estitch%recv(0:oft_env%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,oft_env%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(oft_env%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lerecv(i)%re(linktmp(1,j,i))>lesend(1)%re(m))THEN
        mesh%estitch%leo(m) = .FALSE.
      ELSE IF(lerecv(i)%re(linktmp(1,j,i))==lesend(1)%re(m).AND.oft_env%proc_con(i)<oft_env%rank)THEN
        mesh%estitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    mesh%estitch%lle(:,m+mesh%estitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  mesh%estitch%send(i)%n=ncon(i)
  ALLOCATE(mesh%estitch%send(i)%v(mesh%estitch%send(i)%n))
  mesh%estitch%recv(i)%n=ncon(i)
  ALLOCATE(mesh%estitch%recv(i)%v(mesh%estitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!---------------------------------------------------------------------------
! Clean-up transfers
!---------------------------------------------------------------------------
DEALLOCATE(lesend(1)%le)
DEALLOCATE(lesend(1)%re)
IF(.NOT.mesh%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,oft_env%nproc_con
    DEALLOCATE(lerecv(i)%le)
    DEALLOCATE(lerecv(i)%re)
  END DO
END IF
!---Check stitching information
CALL oft_stitch_check(mesh%estitch)
DEALLOCATE(linktmp,ncon,echeck)
DEALLOCATE(lesend,lerecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_elinkage
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_flinkage
!------------------------------------------------------------------------------
!> Construct processor to processor face linkage for stitching operations.
!! - Create linkage of boundary faces to other processors.
!! - Determine ownership for shared faces
!! - Set global boundary face flag
!------------------------------------------------------------------------------
SUBROUTINE mesh_global_flinkage(mesh)
class(oft_mesh), intent(inout) :: mesh
TYPE :: fout
  INTEGER(8), POINTER, DIMENSION(:) :: lf => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: rf => NULL()
END TYPE fout
TYPE(fout), ALLOCATABLE, DIMENSION(:) :: lfrecv,lfsend
INTEGER(8), POINTER, DIMENSION(:) :: lftmp,lfout
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: lsort
INTEGER(i4), ALLOCATABLE :: bpi(:),bfi(:),linktmp(:,:,:),isort(:),ncon(:)
INTEGER(i4) :: m,mm,nftmp,nbfmax,ll(4)
INTEGER(i4) :: ierr,i,j,js,jn,nffl
LOGICAL :: ftest
LOGICAL, ALLOCATABLE :: fcheck(:,:)
! TYPE(oft_seam) :: fstitch
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Face linkage'
CALL oft_increase_indent
IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---------------------------------------------------------------------------
! Determine maximum boundary face count
!---------------------------------------------------------------------------
nftmp=mesh%nbf
nbfmax=mesh%nbf
IF(.NOT.mesh%fullmesh)nbfmax=oft_mpi_max(nftmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam faces =',nbfmax
mesh%fstitch%nbemax=nbfmax
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,mesh%nbf,0:oft_env%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:oft_env%nproc_con),fcheck(0:oft_env%nproc_con,mesh%nbf))
ncon=0
fcheck=.FALSE.
!---Pointer into linkage array
ALLOCATE(mesh%fstitch%kle(0:oft_env%nproc_con+1))
mesh%fstitch%kle = 0
!---Global boundary face tag
IF(ASSOCIATED(mesh%global%gbf))DEALLOCATE(mesh%global%gbf)
ALLOCATE(mesh%global%gbf(mesh%nf))
mesh%global%gbf = .FALSE.
DO i=1,mesh%nf
  IF(mesh%bf(i))mesh%global%gbf(i)=.TRUE.
END DO
!---Allocate temporary Send/Recv arrays
ALLOCATE(lfsend(oft_env%nproc_con+1),lfrecv(oft_env%nproc_con+1))
ALLOCATE(lfsend(1)%lf(nbfmax))
ALLOCATE(lfsend(1)%rf(nbfmax))
IF(.NOT.mesh%fullmesh)THEN
  DO i=1,oft_env%nproc_con
    ALLOCATE(lfrecv(i)%lf(nbfmax))
    ALLOCATE(lfrecv(i)%rf(nbfmax))
  END DO
ENDIF
!---Create link from local index to boundary index
ALLOCATE(bpi(mesh%np))
CALL get_inverse_map(mesh%lbp,mesh%nbp,bpi,mesh%np)
!---Count corner ownership for load balancing
CALL oft_random_number(lfsend(1)%rf,nbfmax)
DO i=1,mesh%nbf
  DO m=1,mesh%face_np
    IF(mesh%pstitch%leo(bpi(mesh%lf(m,mesh%lbf(i)))))lfsend(1)%rf(i) = lfsend(1)%rf(i) + 1.d0
  END DO
END DO
!---Point dummy output array to Send array
lfout=>lfsend(1)%lf
lfout=0
!$omp parallel do private(m,ll,js,jn,j,ftest)
DO i=1,mesh%nbf
  DO m=1,mesh%face_np
    ll(m)=bpi(mesh%lf(m,mesh%lbf(i))) ! Get endpoints in boundary index
  END DO
  lfout(i)=ABS(mesh%global%lf(mesh%lbf(i))) ! Populate output array with global face indices
  DO m=0,oft_env%nproc_con  ! Flag processors to be checked
    js=mesh%pstitch%kle(m)
    jn=mesh%pstitch%kle(m+1)-1
    ftest=.TRUE.
    DO j=1,mesh%face_np
      ftest=ftest.AND.ANY(mesh%pstitch%lle(1,js:jn)==ll(j))
    END DO
    IF(ftest)fcheck(m,i)=.TRUE.
  END DO
END DO
DEALLOCATE(bpi)
!---Point dummy Send arrays to main Send array
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  DO i=2,oft_env%nproc_con
    lfsend(i)%lf=>lfsend(1)%lf
    lfsend(i)%rf=>lfsend(1)%rf
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!---------------------------------------------------------------------------
! Create Send and Recv calls
!---------------------------------------------------------------------------
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lfsend(j)%lf,nbfmax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%send(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_flinkage',__FILE__)
    CALL MPI_IRECV(lfrecv(j)%lf,nbfmax,OFT_MPI_I8,oft_env%proc_con(j),1,oft_env%COMM,oft_env%recv(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_flinkage',__FILE__)
  END DO
!---------------------------------------------------------------------------
! Loop over each connected processor
!---------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    IF(oft_mpi_check_reqs(oft_env%nproc_con,oft_env%recv))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(oft_env%nproc_con,oft_env%recv,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_flinkage',__FILE__)
    lftmp=>lfrecv(j)%lf ! Point dummy input array to current Recv array
    !---
    nffl=0
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,mesh%nbf ! Loop over boundary faces
      IF(.NOT.fcheck(j,m))CYCLE
      DO mm=1,nbfmax ! Loop over input faces
        IF(lfout(m)==lftmp(mm))THEN ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          mesh%global%gbf(mesh%lbf(m))=.FALSE. ! Face is not on global boundary
          nffl=nffl+1
          EXIT
        END IF
      END DO
    END DO
    mesh%fstitch%kle(j)=nffl
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_flinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Check internal connections
!---------------------------------------------------------------------------
ALLOCATE(mesh%fstitch%leo(mesh%nbf))
mesh%fstitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic faces on current processor
nffl=0
IF(mesh%periodic%nper>0)THEN
  IF(ASSOCIATED(mesh%periodic%lf))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bfi(mesh%nf))
    CALL get_inverse_map(mesh%lbf,mesh%nbf,bfi,mesh%nf)
    !---Find face pairs from periodic parents
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,mesh%nbf
      mm=mesh%lbf(m)
      IF(mesh%periodic%lf(mm)>0)THEN
        !$omp critical
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/bfi(mesh%periodic%lf(mm)),m/)
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/m,bfi(mesh%periodic%lf(mm))/)
        !$omp end critical
        nffl=nffl+2
        mesh%global%gbf(mm)=.FALSE.
        mesh%global%gbf(mesh%periodic%lf(mm))=.FALSE.
        IF(bfi(mesh%periodic%lf(mm))<m)THEN
          mesh%fstitch%leo(m) = .FALSE.
          lfsend(1)%rf(m)=-1.d0
        END IF
      END IF
    END DO
    DEALLOCATE(bfi)
  ELSE
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,mesh%nbf ! Loop over boundary faces
      IF(.NOT.fcheck(0,m))CYCLE
      DO mm=m+1,mesh%nbf ! Loop over boundary faces
        IF(lfout(m)==lfout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          nffl=nffl+2
          mesh%global%gbf(mesh%lbf(m))=.FALSE. ! Face is not on global boundary
          mesh%global%gbf(mesh%lbf(mm))=.FALSE. ! Face is not on global boundary
          mesh%fstitch%leo(mm) = .FALSE.
          lfsend(1)%rf(mm)=-1.d0
          EXIT
        END IF
        IF(lfout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
mesh%fstitch%kle(0)=nffl
!---Transfer load balancing arrays
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
  DO j=1,oft_env%nproc_con
    CALL MPI_ISEND(lfsend(j)%rf,nbfmax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%send(j),ierr)
    CALL MPI_IRECV(lfrecv(j)%rf,nbfmax,OFT_MPI_R8,oft_env%proc_con(j), &
                   1,oft_env%COMM,oft_env%recv(j),ierr)
  END DO
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%recv,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
  CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_flinkage",__FILE__)
#endif
END IF
!---------------------------------------------------------------------------
! Condense linkage to sparse rep
!---------------------------------------------------------------------------
mesh%fstitch%nle=sum(mesh%fstitch%kle)
mesh%fstitch%kle(oft_env%nproc_con+1)=mesh%fstitch%nle+1
DO i=oft_env%nproc_con,0,-1 ! cumulative unique edge linkage count
  mesh%fstitch%kle(i)=mesh%fstitch%kle(i+1)-mesh%fstitch%kle(i)
END DO
IF(mesh%fstitch%kle(0)/=1)CALL oft_abort('Bad face linkage count','mesh_global_flinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(mesh%fstitch%lle(2,mesh%fstitch%nle))
!---
mesh%fstitch%full=mesh%fullmesh
! mesh%fstitch%nbemax=mesh%linkage%nbfmax
mesh%fstitch%nbe=mesh%nbf
mesh%fstitch%lbe=>mesh%lbf
! mesh%fstitch%nle=mesh%linkage%nlf
! mesh%fstitch%kle=>mesh%linkage%klf
! mesh%fstitch%lle=>mesh%linkage%llf
! mesh%fstitch%leo=>mesh%linkage%lfo
ALLOCATE(mesh%fstitch%send(0:oft_env%nproc_con),mesh%fstitch%recv(0:oft_env%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,oft_env%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(oft_env%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lfrecv(i)%rf(linktmp(1,j,i))>lfsend(1)%rf(m))THEN
        mesh%fstitch%leo(m) = .FALSE.
      ELSE IF(lfrecv(i)%rf(linktmp(1,j,i))==lfsend(1)%rf(m).AND.oft_env%proc_con(i)<oft_env%rank)THEN
        mesh%fstitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    mesh%fstitch%lle(:,m+mesh%fstitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  mesh%fstitch%send(i)%n=ncon(i)
  ALLOCATE(mesh%fstitch%send(i)%v(mesh%fstitch%send(i)%n))
  mesh%fstitch%recv(i)%n=ncon(i)
  ALLOCATE(mesh%fstitch%recv(i)%v(mesh%fstitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!---------------------------------------------------------------------------
! Clean-up transfers
!---------------------------------------------------------------------------
DEALLOCATE(lfsend(1)%lf)
DEALLOCATE(lfsend(1)%rf)
IF(.NOT.mesh%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,oft_env%nproc_con
    DEALLOCATE(lfrecv(i)%lf)
    DEALLOCATE(lfrecv(i)%rf)
  END DO
END IF
!---Check stitching information
! fstitch%full=mesh%fullmesh
! fstitch%nbemax=mesh%linkage%nbfmax
! fstitch%nbe=mesh%nbf
! fstitch%lbe=>mesh%lbf
! fstitch%nle=mesh%linkage%nlf
! fstitch%kle=>mesh%linkage%klf
! fstitch%lle=>mesh%linkage%llf
! fstitch%leo=>mesh%linkage%lfo
! fstitch%send=>mesh%fsend
! fstitch%recv=>mesh%frecv
CALL oft_stitch_check(mesh%fstitch)
DEALLOCATE(linktmp,ncon,fcheck)
DEALLOCATE(lfsend,lfrecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_flinkage
!------------------------------------------------------------------------------
! SUBROUTINE: plinkage_from_parent
!------------------------------------------------------------------------------
!> Construct processor to processor point linkage for stitching operations.
!! - Create linkage of boundary points to other processors
!! - Determine ownership for shared points
!------------------------------------------------------------------------------
SUBROUTINE plinkage_from_parent(smesh,parent)
class(oft_bmesh), intent(inout) :: smesh
class(oft_mesh), intent(inout) :: parent
INTEGER(i4) :: i,j,k,m
INTEGER(i4), POINTER :: lpbound(:),lploc(:)
IF(smesh%np==0)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Point linkage'
ALLOCATE(lpbound(smesh%np))
CALL get_inverse_map(smesh%lbp,smesh%nbp,lpbound,smesh%np)
ALLOCATE(lploc(parent%np))
CALL get_inverse_map(smesh%parent%lp,smesh%np,lploc,smesh%parent%np)
!---
ALLOCATE(smesh%pstitch%kle(oft_env%nproc_con+1))
smesh%pstitch%kle=0
DO i=1,oft_env%nproc_con
  m=0
  DO j=parent%pstitch%kle(i),parent%pstitch%kle(i+1)-1
    k=parent%lbp(parent%pstitch%lle(1,j))
    IF(lploc(k)/=0)THEN
      m=m+1
    END IF
  END DO
  smesh%pstitch%kle(i)=m
END DO
!---------------------------------------------------------------------------
! Condense linkage to sparse rep
!---------------------------------------------------------------------------
smesh%pstitch%nle=SUM(smesh%pstitch%kle)
smesh%pstitch%kle(oft_env%nproc_con+1)=smesh%pstitch%nle+1
!---Cumulative unique point linkage count
DO i=oft_env%nproc_con,1,-1
  smesh%pstitch%kle(i)=smesh%pstitch%kle(i+1)-smesh%pstitch%kle(i)
END DO
!---
IF(smesh%pstitch%kle(1)/=1)CALL oft_abort('Bad point linkage count','tetmesh_global_plinkage',__FILE__)
!---
ALLOCATE(smesh%pstitch%lle(2,smesh%pstitch%nle),smesh%pstitch%leo(parent%nbp))
ALLOCATE(smesh%pstitch%send(oft_env%nproc_con),smesh%pstitch%recv(oft_env%nproc_con))
!---
!!$omp do
DO i=1,oft_env%nproc_con
  m=0
  DO j=parent%pstitch%kle(i),parent%pstitch%kle(i+1)-1
    k=parent%lbp(parent%pstitch%lle(1,j))
    IF(lploc(k)/=0)THEN
      smesh%pstitch%lle(1,m+smesh%pstitch%kle(i))=lpbound(lploc(k))
      smesh%pstitch%lle(2,m+smesh%pstitch%kle(i))=m+1
      m=m+1
    END IF
  END DO
  !---Allocate permanent stitching arrays
  smesh%pstitch%send(i)%n=m; smesh%pstitch%recv(i)%n=m
  ALLOCATE(smesh%pstitch%send(i)%v(smesh%pstitch%send(i)%n))
  ALLOCATE(smesh%pstitch%recv(i)%v(smesh%pstitch%recv(i)%n))
END DO
!!$omp end parallel
!---------------------------------------------------------------------------
! Assemble ownership tag of local boundary points
!---------------------------------------------------------------------------
DEALLOCATE(lpbound)
ALLOCATE(lpbound(parent%np))
CALL get_inverse_map(parent%lbp,parent%nbp,lpbound,parent%np)
!$omp parallel do private(k)
DO i=1,smesh%nbp
  k=INT(smesh%parent%lp(smesh%lbp(i)),4)
  IF(lploc(k)/=0)THEN
    smesh%pstitch%leo(i) = parent%pstitch%leo(lpbound(k))
  END IF
END DO
!---------------------------------------------------------------------------
! Clean-up
!---------------------------------------------------------------------------
!---Deallocate temporary work arrays
DEALLOCATE(lpbound,lploc)
DEBUG_STACK_POP
END SUBROUTINE plinkage_from_parent
!------------------------------------------------------------------------------
! SUBROUTINE: elinkage_from_parent
!------------------------------------------------------------------------------
!> Construct processor to processor edge linkage for stitching operations.
!! - Create linkage of boundary edges to other processors with orientation.
!! - Determine ownership for shared edges
!------------------------------------------------------------------------------
SUBROUTINE elinkage_from_parent(smesh,parent)
class(oft_bmesh), intent(inout) :: smesh
class(oft_mesh), intent(inout) :: parent
INTEGER(i4) :: i,j,k,m
INTEGER(i4), POINTER :: leloc(:),lebound(:)
IF(smesh%ne==0)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Edge linkage'
ALLOCATE(lebound(smesh%ne))
CALL get_inverse_map(smesh%lbe,smesh%nbe,lebound,smesh%ne)
ALLOCATE(leloc(parent%ne))
CALL get_inverse_map(smesh%parent%le,smesh%ne,leloc,smesh%parent%ne)
!---
ALLOCATE(smesh%estitch%kle(oft_env%nproc_con+1))
smesh%estitch%kle=0
DO i=1,oft_env%nproc_con
  m=0
  DO j=parent%estitch%kle(i),parent%estitch%kle(i+1)-1
    k=parent%lbe(parent%estitch%lle(1,j))
    IF(leloc(k)/=0)THEN
      m=m+1
    END IF
  END DO
  smesh%estitch%kle(i)=m
END DO
!---------------------------------------------------------------------------
! Condense linkage to sparse rep
!---------------------------------------------------------------------------
smesh%estitch%nle=SUM(smesh%estitch%kle)
smesh%estitch%kle(oft_env%nproc_con+1)=smesh%estitch%nle+1
!---Cumulative unique point linkage count
DO i=oft_env%nproc_con,1,-1
  smesh%estitch%kle(i)=smesh%estitch%kle(i+1)-smesh%estitch%kle(i)
END DO
!---
IF(smesh%estitch%kle(1)/=1)CALL oft_abort('Bad point linkage count','tetmesh_global_plinkage',__FILE__)
!---
ALLOCATE(smesh%estitch%lle(2,smesh%estitch%nle),smesh%estitch%leo(parent%nbe))
ALLOCATE(smesh%estitch%send(oft_env%nproc_con),smesh%estitch%recv(oft_env%nproc_con))
!---
!!$omp do
DO i=1,oft_env%nproc_con
  m=0
  DO j=parent%estitch%kle(i),parent%estitch%kle(i+1)-1
    k=parent%lbe(parent%estitch%lle(1,j))
    IF(leloc(k)/=0)THEN
      smesh%estitch%lle(1,m+smesh%estitch%kle(i))=lebound(leloc(k))
      smesh%estitch%lle(2,m+smesh%estitch%kle(i))=m+1
      m=m+1
    END IF
  END DO
  !---Allocate permanent stitching arrays
  smesh%estitch%send(i)%n=m; smesh%estitch%recv(i)%n=m
  ALLOCATE(smesh%estitch%send(i)%v(smesh%estitch%send(i)%n))
  ALLOCATE(smesh%estitch%recv(i)%v(smesh%estitch%recv(i)%n))
END DO
!!$omp end parallel
!---------------------------------------------------------------------------
! Assemble ownership tag of local boundary points
!---------------------------------------------------------------------------
DEALLOCATE(lebound)
ALLOCATE(lebound(parent%ne))
CALL get_inverse_map(parent%lbe,parent%nbe,lebound,parent%ne)
!$omp parallel do private(k)
DO i=1,smesh%nbe
  k=INT(ABS(smesh%parent%le(smesh%lbe(i))),4)
  IF(leloc(k)/=0)THEN
    smesh%estitch%leo(i) = parent%estitch%leo(lebound(k))
  END IF
END DO
!---------------------------------------------------------------------------
! Clean-up
!---------------------------------------------------------------------------
!---Deallocate temporary work arrays
DEALLOCATE(lebound,leloc)
DEBUG_STACK_POP
END SUBROUTINE elinkage_from_parent
end module oft_mesh_global
