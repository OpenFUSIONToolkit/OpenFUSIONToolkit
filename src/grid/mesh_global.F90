!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
MODULE oft_mesh_global
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, oft_amesh, oft_bmesh, oft_init_seam
USE oft_mesh_local, ONLY: mesh_local_init, mesh_local_init, oft_metis_partmesh, &
  bmesh_local_init
USE oft_stitching, ONLY: oft_seam, oft_stitch_check
IMPLICIT NONE
#include "local.h"
contains
!---------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!---------------------------------------------------------------------------------
subroutine mesh_global_init(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i4) :: i
DEBUG_STACK_PUSH
call mesh_global_sync(self) ! Sync mesh information
call mesh_local_init(self)
!---
allocate(self%global%lp(self%np))
!$omp parallel do
do i=1,self%np
  self%global%lp(i)=i
end do
!---
allocate(self%global%le(self%ne))
!$omp parallel do
do i=1,self%ne
  self%global%le(i)=i
end do
!---
allocate(self%global%lf(self%nf))
!$omp parallel do
do i=1,self%nf
  self%global%lf(i)=i
end do
!---
ALLOCATE(self%global%lc(self%nc))
!$omp parallel do
do i=1,self%nc
  self%global%lc(i)=i
end do
!---
self%global%np=self%np
self%global%ne=self%ne
self%global%nf=self%nf
self%global%nc=self%nc
ALLOCATE(self%global%gbp(self%np))
ALLOCATE(self%global%gbe(self%ne))
ALLOCATE(self%global%gbf(self%nf))
ALLOCATE(self%global%gbc(self%nc))
self%global%gbp=self%bp
self%global%gbe=self%be
self%global%gbf=self%bf
self%global%gbc=self%bc
self%save%nbf=self%nbf
DEBUG_STACK_POP
end subroutine mesh_global_init
!---------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!---------------------------------------------------------------------------------
subroutine bmesh_global_init(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
INTEGER(i4) :: i,j,np,nc
INTEGER(i4), ALLOCATABLE :: lptmp(:)
DEBUG_STACK_PUSH
!---Create surface mesh from boundary
nc=0
DO i=1,self%nbf
  j=self%lbf(i)
  IF(self%global%gbf(j))nc=nc+1
END DO
self%bmesh%nc=nc
!---
np=0
DO i=1,self%nbp
  j=self%lbp(i)
  IF(self%global%gbp(j))np=np+1
END DO
IF(np==0)self%bmesh%skip=.TRUE.
self%bmesh%np=np
ALLOCATE(self%bmesh%r(3,np))
ALLOCATE(self%bmesh%parent%lp(np))
ALLOCATE(lptmp(self%np))
lptmp=0
np=0
DO i=1,self%nbp
  j=self%lbp(i)
  IF(self%global%gbp(j))THEN
    np=np+1
    self%bmesh%r(:,np)=self%r(:,j)
    self%bmesh%parent%lp(np)=j
    lptmp(j)=np
  END IF
END DO
!---
ALLOCATE(self%bmesh%lc(self%bmesh%cell_np,nc))
ALLOCATE(self%bmesh%reg(nc))
ALLOCATE(self%bmesh%parent%lf(nc))
nc=0
DO i=1,self%nbf
  j=self%lbf(i)
  IF(self%global%gbf(j))THEN
    nc=nc+1
    self%bmesh%lc(:,nc)=lptmp(self%lf(:,j))
    self%bmesh%parent%lf(nc)=j
    self%bmesh%reg(nc)=self%reg(self%lfc(1,j))
  END IF
END DO
DEALLOCATE(lptmp)
!---
CALL bmesh_local_init(self%bmesh,self)
!---Copy high-order mapping if available
IF(ASSOCIATED(self%ho_info%r))THEN
  self%bmesh%order=1
  self%bmesh%ho_info%nep=self%ho_info%nep
  self%bmesh%ho_info%ncp=self%ho_info%nfp
  ALLOCATE(self%bmesh%ho_info%r(3,self%bmesh%ne+self%bmesh%nc))
  ALLOCATE(self%bmesh%ho_info%lep(self%bmesh%ho_info%nep,self%bmesh%ne))
  ALLOCATE(self%bmesh%ho_info%lcp(self%bmesh%ho_info%ncp,self%bmesh%nc))
  !!$omp parallel do private(j)
  do i=1,self%bmesh%ne
    self%bmesh%ho_info%lep(1,i)=i
    j=INT(ABS(self%bmesh%parent%le(i)),4)
    self%bmesh%ho_info%r(:,i)=self%ho_info%r(:,self%ho_info%lep(1,j))
  enddo
  IF(self%bmesh%ho_info%ncp==1)THEN
    !!$omp parallel do private(j)
    do i=1,self%bmesh%nc
      j=INT(ABS(self%bmesh%parent%lf(i)),4)
      self%bmesh%ho_info%lcp(1,i)=i+self%bmesh%ne
      self%bmesh%ho_info%r(:,i+self%bmesh%ne)=self%ho_info%r(:,self%ho_info%lfp(1,j))
    enddo
  END IF
END IF
DEBUG_STACK_POP
end subroutine bmesh_global_init
!---------------------------------------------------------------------------------
!> Driver for global mesh initialization
!!
!! Initialize global mesh enviroment
!! - Sync mesh information from proc 0 to all
!! - Set base global geometry counts, indices and boundary info
!! - Construct lowest level mesh
!---------------------------------------------------------------------------------
subroutine smesh_global_init(self)
class(oft_bmesh), intent(inout) :: self !< Mesh object
integer(i4) :: i
DEBUG_STACK_PUSH
call mesh_global_sync(self) ! Sync mesh information
call bmesh_local_init(self)
!---
self%global%np=self%np
ALLOCATE(self%global%gbp(self%np))
self%global%gbp=self%bp
allocate(self%global%lp(self%np))
!$omp parallel do
do i=1,self%np
  self%global%lp(i)=i
end do
!---
self%global%ne=self%ne
ALLOCATE(self%global%gbe(self%ne))
self%global%gbe=self%be
allocate(self%global%le(self%ne))
!$omp parallel do
do i=1,self%ne
  self%global%le(i)=i
end do
!---
self%global%nc=self%nc
ALLOCATE(self%global%gbc(self%nc))
self%global%gbc=self%bc
ALLOCATE(self%global%lc(self%nc))
!$omp parallel do
do i=1,self%nc
  self%global%lc(i)=i
end do
!---
ALLOCATE(self%lco(self%nc))
self%lco=1
ALLOCATE(self%pstitch%leo(self%nbp))
ALLOCATE(self%estitch%leo(self%nbe))
self%pstitch%leo=.TRUE.
self%estitch%leo=.TRUE.
DEBUG_STACK_POP
end subroutine smesh_global_init
!---------------------------------------------------------------------------------
!> Driver for global seam linkage construction
!!
!! Construct inter-processor linkage information
!! - Link seam points, edges and faces
!---------------------------------------------------------------------------------
subroutine mesh_global_link(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
DEBUG_STACK_PUSH
IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'Generating domain linkage'
CALL oft_increase_indent
call mesh_global_plinkage(self) ! Link seam points
call mesh_global_elinkage(self) ! Link seam edges
call mesh_global_flinkage(self) ! Link seam faces
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_global_link
!---------------------------------------------------------------------------------
!> Driver for global seam linkage construction
!!
!! Construct inter-processor linkage information
!! - Link seam points, edges and faces
!---------------------------------------------------------------------------------
subroutine bmesh_global_link(self,parent)
class(oft_bmesh), intent(inout) :: self !< Mesh object
class(oft_mesh), optional, intent(inout) :: parent !< Parent volume mesh (if present)
DEBUG_STACK_PUSH
IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'Generating boundary domain linkage'
CALL oft_increase_indent
IF(PRESENT(parent))THEN
  call plinkage_from_parent(self,parent) ! Link seam points
  call elinkage_from_parent(self,parent) ! Link seam edges
ELSE
  call mesh_global_plinkage(self)
  call mesh_global_elinkage(self)
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine bmesh_global_link
!---------------------------------------------------------------------------------
!> Scatters base mesh information to all processors.
!!
!! Scatters mesh information from head task to all other tasks using MPI_BCAST calls.
!! - Communicates base mesh information (np,nc,lc,r)
!---------------------------------------------------------------------------------
subroutine mesh_global_sync(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
#ifdef HAVE_MPI
integer(i4) :: ierr
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Syncing mesh'
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(self%cad_type,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast mesh type
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
!---
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(self%nc,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell count
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(self%np,1,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast point count
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(self%meshname,20,OFT_MPI_CHAR,0,oft_env%COMM,ierr) ! Broadcast mesh name
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
!---
if(oft_env%rank>0)allocate(self%r(3,self%np),self%lc(self%cell_np,self%nc),self%reg(self%nc)) ! Allocate point and cell arrays
CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL MPI_Bcast(self%lc,self%cell_np*self%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
CALL MPI_Bcast(self%r,3*self%np,OFT_MPI_R8,0,oft_env%COMM,ierr) ! Broadcast point list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
call MPI_Bcast(self%reg,self%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast region list
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_global_sync',__FILE__)
DEBUG_STACK_POP
#endif
end subroutine mesh_global_sync
!---------------------------------------------------------------------------------
!> Perform mesh decomposition (METIS) and local construction.
!! - Decompose domain of head task using METIS library.
!! - Scatter decomposition to all tasks.
!!
!! Supported partition methods are:
!! - `part_meth==1`: Partition using METIS library
!! - `2 <= part_meth <= 4`: Axial spatial partitioning along coordinate `part_meth-1`
!! - `part_meth==5`: Cylindrical spatial partitioning in azimuthal direction
!---------------------------------------------------------------------------------
subroutine mesh_global_partition(self,meshpart,part_meth)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), intent(inout) :: meshpart(:) !< Partition flag [self%nc]
integer(i4), intent(in) :: part_meth !< Method to use for partitioning
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
  allocate(cpart(self%nc))
  cpart=0
  if(part_meth==1)THEN ! Metis partitioning
    CALL oft_metis_partMesh(self%nc,self%np,self%cell_np,self%lc,oft_env%nnodes,cpart,ierr)
    IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_global_partition',__FILE__)
  elseif((part_meth>1).AND.(part_meth<5))THEN ! Axial spatial partitioning
    zmin = MINVAL(self%r(part_meth-1,:))
    zmax = MAXVAL(self%r(part_meth-1,:))
    zrange = (zmax-zmin)/REAL(oft_env%nnodes,8)
    ! IF(MOD(mesh%nc,oft_env%nnodes)/=0)CALL oft_abort( &
    !   "Partition count invalid for spatial partitioning",'mesh_global_partition',__FILE__)
    ALLOCATE(pflag(self%np),part_sort(self%np),isort(self%np))
    pflag=oft_env%nnodes+1
    isort=[(i,i=1,self%np)]
    part_sort=self%r(part_meth-1,:)
    CALL sort_array(part_sort,isort,self%np)
    k=1
    DO i=1,self%np
      IF(i>REAL(k*self%np,8)/REAL(oft_env%nnodes,8))THEN
        k=k+1
      END IF
      pflag(isort(i))=k
    END DO
    DO i=1,self%nc
      cpart(i)=oft_env%nnodes+1
      DO k=1,self%cell_np
        cpart(i)=MIN(cpart(i),pflag(self%lc(k,i)))
      END DO
      IF(cpart(i)<1.OR.cpart(i)>oft_env%nnodes)WRITE(*,*)'cBAD',i,cpart(i)
    END DO
  elseif(part_meth==5)THEN ! Toroidal spatial partitioning
    zmin = -pi
    zmax = pi
    zrange = (zmax-zmin)/REAL(oft_env%nnodes,8)
    ! IF(MOD(mesh%nc,oft_env%nnodes)/=0)CALL oft_abort( &
    ! 'Partition count invalid for spatial partitioning','mesh_global_partition',__FILE__)
    ALLOCATE(pflag(self%np),part_sort(self%np),isort(self%np))
    pflag=oft_env%nnodes+1
    isort=[(i,i=1,self%np)]
    part_sort=ATAN2(self%r(2,:),self%r(1,:))
    CALL sort_array(part_sort,isort,self%np)
    k=1
    DO i=1,self%np
      IF(i>REAL(k*self%np,8)/REAL(oft_env%nnodes,8))THEN
        k=k+1
      END IF
      pflag(isort(i))=k
    END DO
    DO i=1,self%nc
      cpart(i)=oft_env%nnodes+1
      DO k=1,self%cell_np
        cpart(i)=MIN(cpart(i),pflag(self%lc(k,i)))
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
    allocate(lptmp(self%np),lctmp(self%nc),lcctmp(self%cell_np,self%nc))
    allocate(pflag(self%np))
    do m=1,oft_env%nnodes
      !---Initialize counters and flags
      pflag=0
      lctmp=0
      lptmp=0
      lcctmp=0
      nctmp=0
      nptmp=0
      do i=1,self%nc ! Loop over global cells
        if(meshpart(i)==m)then ! Processor owns this cell
          nctmp = nctmp + 1 ! Increment cell counter
          lctmp(nctmp) = i ! Link cell to global index
          do k=1,self%cell_np ! Loop over points
            if(pflag(self%lc(k,i))==0)then ! Point is not claimed yet
              nptmp=nptmp+1 ! Increment point counter
              pflag(self%lc(k,i))=nptmp ! Link global point to new index
              lptmp(nptmp)=self%lc(k,i) ! Link point to global index
              lcctmp(k,nctmp)=nptmp ! Insert local point index into local cell list
            else ! Already have this point
              lcctmp(k,nctmp)=pflag(self%lc(k,i)) ! Insert local point index into local cell list
            endif
          enddo
        endif
      enddo
      !---
      allocate(cpart(nctmp))
      CALL oft_metis_partMesh(nctmp,nptmp,self%cell_np,lcctmp,oft_env%ppn,cpart,ierr)
      IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_global_partition',__FILE__)
      !$omp parallel do
      do i=1,nctmp
        meshpart(lctmp(i))=-((m-1)*oft_env%ppn+cpart(i))
      end do
      deallocate(cpart)
    end do
    !$omp parallel do
    do i=1,self%nc
      meshpart(i)=-meshpart(i)
    end do
  end if
  allocate(ncells(oft_env%nprocs))
  ncells=0
  do i=1,self%nc
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
call MPI_Bcast(meshpart,self%nc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast cell partition information
#endif
DEBUG_STACK_POP
end subroutine mesh_global_partition
!---------------------------------------------------------------------------------
!> Perform mesh decomposition and local construction.
!! - Decompose global mesh
!! - Perform mesh construction of local domain
!---------------------------------------------------------------------------------
subroutine mesh_global_decomp(self,part_meth)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), intent(in) :: part_meth !< Method to use for partitioning (see @ref mesh_global_partition)
integer(i4), allocatable :: lptmp(:),lctmp(:),lcctmp(:,:),isort(:),isorti(:)
integer(i4), allocatable :: pflag(:),conflag(:),proc_contmp(:)
integer(i4) :: i,k,npcors
real(r8), allocatable :: rtmp(:,:)
DEBUG_STACK_PUSH
!---Synchronize before deomposition
allocate(self%base%lcpart(self%nc))
call mesh_global_partition(self,self%base%lcpart,part_meth)
!---
npcors=self%np ! Get coarse point count
!---Create local mesh indexing for points and cells
allocate(lptmp(self%np),lctmp(self%nc),lcctmp(self%cell_np,self%nc),rtmp(3,self%np))
allocate(pflag(self%np),conflag(npcors),proc_contmp(oft_env%nprocs))
!---Initialize counters and flags
pflag=0
conflag=0
proc_contmp=0
lctmp=0
lptmp=0
lcctmp=0
rtmp=0.d0
self%nc=0
self%np=0
do i=1,INT(self%global%nc,4) ! Loop over global cells
  if(self%base%lcpart(i)==(oft_env%rank+1))then ! Processor owns this cell
    self%nc = self%nc + 1 ! Increment cell counter
    lctmp(self%nc) = i ! Link cell to global index
    do k=1,self%cell_np ! Loop over points
      if(pflag(self%lc(k,i))==0)then ! Point is not claimed yet
        self%np=self%np+1 ! Increment point counter
        pflag(self%lc(k,i))=self%np ! Link global point to new index
        conflag(self%global%lp(self%lc(k,i)))=self%np ! Mark global point list
        lptmp(self%np)=self%lc(k,i) ! Link point to global index
        rtmp(:,self%np)=self%r(:,self%lc(k,i)) ! Insert point into local list
        lcctmp(k,self%nc)=self%np ! Insert local point index into local cell list
      else ! Already have this point
        lcctmp(k,self%nc)=pflag(self%lc(k,i)) ! Insert local point index into local cell list
      endif
    enddo
  endif
enddo
ALLOCATE(self%global%seam) ! Allocate grid seam structure
call mesh_global_proccon(self,conflag,npcors) ! Determine neighbor processors
allocate(self%global%seam%send_reqs(self%global%seam%nproc_con),self%global%seam%recv_reqs(self%global%seam%nproc_con))
!---Sort point list
allocate(isort(self%np))
isort=(/(i,i=1,self%np)/)
CALl sort_array(lptmp, isort, self%np)
allocate(isorti(self%np))
!---Allocate local mesh arrays
allocate(self%r(3,self%np),self%lc(self%cell_np,self%nc))
DO i=1,self%np
  isorti(isort(i))=i
  self%r(:,i)=rtmp(:,isort(i))
END DO
DO i=1,self%nc
  self%lc(:,i)=isorti(lcctmp(:,i))
END DO
! mesh%r=rtmp(:,1:mesh%np) ! Set local point list
! mesh%lc=lcctmp(:,1:mesh%nc) ! Set local cell list
DEALLOCATE(isort,isorti)
!---Get indexing on base mesh
allocate(self%base%lp(self%np),self%base%lc(self%nc))
!$omp parallel do
do i=1,self%np
  self%base%lp(i) = lptmp(i) ! Set base point index
  lptmp(i) = INT(self%global%lp(lptmp(i)),4) ! Set global point index
end do
!$omp parallel do
do i=1,self%nc
  self%base%lc(i) = lctmp(i) ! Set base cell index
  lctmp(i) = INT(self%global%lc(lctmp(i)),4) ! Set global cell index
end do
!---Assign global indices
allocate(self%global%lp(self%np),self%global%lc(self%nc))
self%global%lp = INT(lptmp(1:self%np),8) ! Set global point index
self%global%lc = INT(lctmp(1:self%nc),8) ! Set global cell index
!---Kill temporary work arrays
deallocate(lptmp,lctmp,lcctmp,rtmp)
deallocate(pflag,conflag)
DEBUG_STACK_POP
end subroutine mesh_global_decomp
!---------------------------------------------------------------------------------
!> Determine neighboring processors for MPI linkage
!---------------------------------------------------------------------------------
subroutine mesh_global_proccon(self,lptmp,np)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), intent(in) :: np !< number of points
integer(i4), intent(in) :: lptmp(np) !< List of global point indices from initialization
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
self%global%seam%nproc_con=sum(b) ! Number of processor connections
self%global%seam%proc_split=self%global%seam%nproc_con
allocate(self%global%seam%proc_con(self%global%seam%nproc_con)) ! Allocate processor linkage list
k=1
do i=1,oft_env%nprocs
  if(b(i)==1)then ! Processor links to this domain
    self%global%seam%proc_con(k)=i-1 ! Populate processor linkage list
    IF((i-1>oft_env%rank).AND.(k-1<self%global%seam%proc_split))self%global%seam%proc_split=k-1
    k=k+1
  endif
enddo
DEALLOCATE(b)
DEBUG_STACK_POP
end subroutine mesh_global_proccon
!---------------------------------------------------------------------------------
!> Construct processor to processor point linkage for stitching operations.
!! - Create linkage of boundary points to other processors
!! - Determine ownership for shared points
!---------------------------------------------------------------------------------
SUBROUTINE mesh_global_plinkage(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
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
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Point linkage'
CALL oft_increase_indent
IF(.NOT.self%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!------------------------------------------------------------------------------
! Determine maximum boundary point count
!------------------------------------------------------------------------------
nptmp=self%nbp
nbpmax=self%nbp
IF(.NOT.self%fullmesh)nbpmax=oft_mpi_max(nptmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam points =',nbpmax
self%pstitch%nbemax=nbpmax
!---Initialize seam structure from mesh
CALL oft_init_seam(self,self%pstitch)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,(self%periodic%nper+2)*self%nbp,0:self%pstitch%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:self%pstitch%nproc_con))
ncon=0
!---Pointer into linkage array
ALLOCATE(self%pstitch%kle(0:self%pstitch%nproc_con+1))
self%pstitch%kle = 0
!---Global boundary point tag
IF(ASSOCIATED(self%global%gbp))DEALLOCATE(self%global%gbp)
ALLOCATE(self%global%gbp(self%np))
self%global%gbp = .FALSE.
!---Allocate temporary Send/Recv arrays
ALLOCATE(lpsend(self%pstitch%nproc_con+1),lprecv(self%pstitch%nproc_con+1))
ALLOCATE(lpsend(1)%lp(nbpmax))
ALLOCATE(lpsend(1)%rp(nbpmax))
IF(.NOT.self%fullmesh)THEN
  DO i=1,self%pstitch%nproc_con
    ALLOCATE(lprecv(i)%lp(nbpmax))
    ALLOCATE(lprecv(i)%rp(nbpmax))
  END DO
END IF
!---Count point to point interactions for load balancing
CALL oft_random_number(lpsend(1)%rp,nbpmax)
DO i=1,self%nbp
  lpsend(1)%rp(i)=lpsend(1)%rp(i)+REAL(self%kpp(self%lbp(i)+1)-self%kpp(self%lbp(i)),8)
END DO
!---Point dummy output array to Send array
lpout=>lpsend(1)%lp
lpout=0
!---Populate with global point index
!$omp parallel do
DO i=1,self%nbp
  lpout(i) = self%global%lp(self%lbp(i))
END DO
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  !---Point dummy Send arrays to main Send array
  DO i=2,self%pstitch%nproc_con
    lpsend(i)%lp=>lpsend(1)%lp
    lpsend(i)%rp=>lpsend(1)%rp
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!------------------------------------------------------------------------------
! Create Send and Recv calls
!------------------------------------------------------------------------------
  DO j=1,self%pstitch%nproc_con
    CALL MPI_ISEND(lpsend(j)%lp,nbpmax,OFT_MPI_I8,self%pstitch%proc_con(j),1,oft_env%COMM,self%pstitch%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_plinkage',__FILE__)
    CALL MPI_IRECV(lprecv(j)%lp,nbpmax,OFT_MPI_I8,self%pstitch%proc_con(j),1,oft_env%COMM,self%pstitch%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_plinkage',__FILE__)
  END DO
!------------------------------------------------------------------------------
! Loop over each connected processor
!------------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    !---All recieves have been processed
    IF(oft_mpi_check_reqs(self%pstitch%nproc_con,self%pstitch%recv_reqs))EXIT
    !---Wait for completed recieve
    CALL oft_mpi_waitany(self%pstitch%nproc_con,self%pstitch%recv_reqs,j,ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_plinkage',__FILE__)
    !---Point dummy input array to current Recv array
    lptmp=>lprecv(j)%lp
    !---Determine location of boundary points on other processors
    nppl=0
    !$omp parallel do private(mm) reduction(+:nppl)
    DO m=1,self%nbp ! Loop over boundary points
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
    self%pstitch%kle(j)=nppl
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_plinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Check internal connections
!------------------------------------------------------------------------------
ALLOCATE(self%pstitch%leo(self%nbp))
self%pstitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic points on current processor
nppl=0
IF(self%periodic%nper>0)THEN
  IF(ASSOCIATED(self%periodic%lp))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bpi(self%np))
    CALL get_inverse_map(self%lbp,self%nbp,bpi,self%np)
    !---Construct child point list
    ALLOCATE(child_list(8,self%nbp))
    child_list=0
    DO m=1,self%nbp
      mm=self%lbp(m)
      IF(self%periodic%lp(mm)>0)THEN
        mm=bpi(self%periodic%lp(mm))
        DO i=1,8
          IF(child_list(i,mm)==0)THEN
            child_list(i,mm)=m
            EXIT
          END IF
        END DO
      END IF
    END DO
    !---Find point pairs by linking child points
    DO m=1,self%nbp
      DO i=1,8
        !---Link parent point to child
        IF(child_list(i,m)>0)THEN
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/child_list(i,m),m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,child_list(i,m)/)
          self%pstitch%leo(child_list(i,m))=.FALSE.
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
    DO m=1,self%nbp ! Loop over boundary points
      DO mm=m+1,self%nbp ! Loop over boundary points
        IF(lpout(m)==lpout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          nppl=nppl+2
          self%pstitch%leo(mm) = .FALSE.
          lpsend(1)%rp(mm)=-1.d0
        END IF
        IF(lpout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
self%pstitch%kle(0)=nppl
!---Transfer load balancing arrays
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(self%pstitch%nproc_con,self%pstitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
  DO j=1,self%pstitch%nproc_con
    CALL MPI_ISEND(lpsend(j)%rp,nbpmax,OFT_MPI_R8,self%pstitch%proc_con(j), &
                   1,oft_env%COMM,self%pstitch%send_reqs(j),ierr)
    CALL MPI_IRECV(lprecv(j)%rp,nbpmax,OFT_MPI_R8,self%pstitch%proc_con(j), &
                   1,oft_env%COMM,self%pstitch%recv_reqs(j),ierr)
  END DO
  CALL oft_mpi_waitall(self%pstitch%nproc_con,self%pstitch%recv_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
  CALL oft_mpi_waitall(self%pstitch%nproc_con,self%pstitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_plinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_plinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Condense linkage to sparse rep
!------------------------------------------------------------------------------
self%pstitch%nle=SUM(self%pstitch%kle)
self%pstitch%kle(self%pstitch%nproc_con+1)=self%pstitch%nle+1
!---Cumulative unique point linkage count
DO i=self%pstitch%nproc_con,0,-1
  self%pstitch%kle(i)=self%pstitch%kle(i+1)-self%pstitch%kle(i)
END DO
IF(self%pstitch%kle(0)/=1)CALL oft_abort('Bad point linkage count','mesh_global_plinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(self%pstitch%lle(2,self%pstitch%nle))
!---Create seam object
self%pstitch%full=self%fullmesh
self%pstitch%nbe=self%nbp
self%pstitch%lbe=>self%lbp
ALLOCATE(self%pstitch%send(0:self%pstitch%nproc_con),self%pstitch%recv(0:self%pstitch%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,self%pstitch%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%pstitch%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lprecv(i)%rp(linktmp(1,j,i))>lpsend(1)%rp(m))THEN
        self%pstitch%leo(m) = .FALSE.
      ELSE IF(lprecv(i)%rp(linktmp(1,j,i))==lpsend(1)%rp(m).AND.self%pstitch%proc_con(i)<oft_env%rank)THEN
        self%pstitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    self%pstitch%lle(:,m+self%pstitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Create preallocated send/recv arrays
  self%pstitch%send(i)%n=ncon(i)
  ALLOCATE(self%pstitch%send(i)%v(self%pstitch%send(i)%n))
  self%pstitch%recv(i)%n=ncon(i)
  ALLOCATE(self%pstitch%recv(i)%v(self%pstitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!------------------------------------------------------------------------------
! Clean-up transfers
!------------------------------------------------------------------------------
DEALLOCATE(lpsend(1)%lp)
DEALLOCATE(lpsend(1)%rp)
IF(.NOT.self%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,self%pstitch%nproc_con
    DEALLOCATE(lprecv(i)%lp)
    DEALLOCATE(lprecv(i)%rp)
  END DO
END IF
!---Check stitching information
CALL oft_stitch_check(self%pstitch)
DEALLOCATE(linktmp,ncon)
DEALLOCATE(lpsend,lprecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_plinkage
!---------------------------------------------------------------------------------
!> Construct processor to processor edge linkage for stitching operations.
!! - Create linkage of boundary edges to other processors with orientation.
!! - Determine ownership for shared edges
!---------------------------------------------------------------------------------
SUBROUTINE mesh_global_elinkage(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
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
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Edge linkage'
CALL oft_increase_indent
IF(.NOT.self%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!------------------------------------------------------------------------------
! Determine maximum boundary edge count
!------------------------------------------------------------------------------
netmp=self%nbe
nbemax=self%nbe
IF(.NOT.self%fullmesh)nbemax=oft_mpi_max(netmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam edges =',nbemax
self%estitch%nbemax=nbemax
!---Initialize seam structure from mesh
CALL oft_init_seam(self,self%estitch)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,(self%periodic%nper+2)*self%nbe,0:self%estitch%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:self%estitch%nproc_con),echeck(0:self%estitch%nproc_con,self%nbe))
ncon=0
echeck=.FALSE.
!---Pointer into linkage array
ALLOCATE(self%estitch%kle(0:self%estitch%nproc_con+1))
self%estitch%kle = 0
!---Global boundary edge tag (setup for surface meshes only)
IF(ASSOCIATED(self%global%gbe))DEALLOCATE(self%global%gbe)
ALLOCATE(self%global%gbe(self%ne))
self%global%gbe = .FALSE.
set_gbe=.FALSE.
SELECT TYPE(self)
CLASS IS(oft_bmesh)
  IF(.NOT.ASSOCIATED(self%parent))set_gbe=.TRUE.
END SELECT
IF(set_gbe)THEN
  DO i=1,self%ne
    IF(self%be(i))self%global%gbe(i)=.TRUE.
  END DO
END IF
!---Allocate temporary Send/Recv arrays
ALLOCATE(lesend(self%estitch%nproc_con+1),lerecv(self%estitch%nproc_con+1))
ALLOCATE(lesend(1)%le(nbemax))
ALLOCATE(lesend(1)%re(nbemax))
IF(.NOT.self%fullmesh)THEN
  DO i=1,self%estitch%nproc_con
    ALLOCATE(lerecv(i)%le(nbemax))
    ALLOCATE(lerecv(i)%re(nbemax))
  END DO
ENDIF
!---Count edge to edge interactions for load balancing
CALL oft_random_number(lesend(1)%re,nbemax)
DO i=1,self%nbe
  lesend(1)%re(i)=lesend(1)%re(i)+REAL(self%kee(self%lbe(i)+1)-self%kee(self%lbe(i)),8)
END DO
!---Point dummy output array to Send array
leout=>lesend(1)%le
leout=0
!---Create link from local index to boundary index
ALLOCATE(bpi(self%np))
CALL get_inverse_map(self%lbp,self%nbp,bpi,self%np)
!---
!$omp parallel do private(m,js,jn,ll,j,etest)
DO i=1,self%nbe
  DO m=1,2
    ll(m)=bpi(self%le(m,self%lbe(i))) ! Get endpoints in boundary index
  END DO
  leout(i)=ABS(self%global%le(self%lbe(i))) ! Populate linkage array
  DO m=0,self%estitch%nproc_con ! Flag processors to be checked
    js=self%pstitch%kle(m)
    jn=self%pstitch%kle(m+1)-1
    etest=.TRUE.
    DO j=1,2
      etest=etest.AND.ANY(self%pstitch%lle(1,js:jn)==ll(j))
    END DO
    IF(etest)echeck(m,i)=.TRUE.
  END DO
END DO
DEALLOCATE(bpi)
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  !---Point dummy Send arrays to main Send array
  DO i=2,self%estitch%nproc_con
    lesend(i)%le=>lesend(1)%le
    lesend(i)%re=>lesend(1)%re
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!------------------------------------------------------------------------------
! Create Send and Recv calls
!------------------------------------------------------------------------------
  DO j=1,self%estitch%nproc_con
    CALL MPI_ISEND(lesend(j)%le,nbemax,OFT_MPI_I8,self%estitch%proc_con(j),1,oft_env%COMM,self%estitch%send_reqs(j), ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_elinkage',__FILE__)
    CALL MPI_IRECV(lerecv(j)%le,nbemax,OFT_MPI_I8,self%estitch%proc_con(j),1,oft_env%COMM,self%estitch%recv_reqs(j), ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_elinkage',__FILE__)
  END DO
!------------------------------------------------------------------------------
! Loop over each connected processor
!------------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    !---All recieves have been processed
    IF(oft_mpi_check_reqs(self%estitch%nproc_con,self%estitch%recv_reqs))EXIT
    !---Wait for completed recieve
    CALL oft_mpi_waitany(self%estitch%nproc_con,self%estitch%recv_reqs,j,ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_elinkage',__FILE__)
    !---Point dummy input array to current Recv array
    letmp=>lerecv(j)%le
    !---
    neel=0
    !$omp parallel do private(mm) reduction(+:neel)
    DO m=1,self%nbe ! Loop over boundary edges
      IF(.NOT.echeck(j,m))CYCLE
      DO mm=1,nbemax ! Loop over input edges
        IF(leout(m)==letmp(mm))THEN ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          IF(set_gbe)self%global%gbe(self%lbe(m))=.FALSE. ! Edge is not on global boundary
          neel=neel+1
        END IF
      END DO
    END DO
    self%estitch%kle(j)=neel
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_elinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Check internal connections
!------------------------------------------------------------------------------
ALLOCATE(self%estitch%leo(self%nbe))
self%estitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic edges on current processor
neel=0
IF(self%periodic%nper>0)THEN
  IF(ASSOCIATED(self%periodic%le))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bei(self%ne))
    CALL get_inverse_map(self%lbe,self%nbe,bei,self%ne)
    !---Construct child edge list
    ALLOCATE(child_list(4,self%nbe))
    child_list=0
    DO m=1,self%nbe
      mm=self%lbe(m)
      IF(self%periodic%le(mm)>0)THEN
        mm=bei(self%periodic%le(mm))
        DO i=1,4
          IF(child_list(i,mm)==0)THEN
            child_list(i,mm)=m
            EXIT
          END IF
        END DO
      END IF
    END DO
    !---Find edge pairs by linking child edges
    DO m=1,self%nbe
      DO i=1,4
        !---Link parent edge to child
        IF(child_list(i,m)>0)THEN
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/child_list(i,m),m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,child_list(i,m)/)
          self%estitch%leo(child_list(i,m))=.FALSE.
          lesend(1)%re(child_list(i,m))=-1.d0
          neel=neel+2
          IF(set_gbe)THEN ! Edge is not on global boundary
            self%global%gbe(self%lbe(m))=.FALSE.
            self%global%gbe(self%lbe(child_list(i,m)))=.FALSE.
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
    DO m=1,self%nbe ! Loop over boundary points
      IF(.NOT.echeck(0,m))CYCLE
      DO mm=m+1,self%nbe ! Loop over boundary points
        IF(leout(m)==leout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          neel=neel+2
          IF(set_gbe)THEN ! Edge is not on global boundary
            self%global%gbe(self%lbe(m))=.FALSE. 
            self%global%gbe(self%lbe(mm))=.FALSE.
          END IF
          self%estitch%leo(mm) = .FALSE.
          lesend(1)%re(mm)=-1.d0
        END IF
        IF(leout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
self%estitch%kle(0)=neel
!---Transfer load balancing arrays
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(self%estitch%nproc_con,self%estitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
  DO j=1,self%estitch%nproc_con
    CALL MPI_ISEND(lesend(j)%re,nbemax,OFT_MPI_R8,self%estitch%proc_con(j), &
                   1,oft_env%COMM,self%estitch%send_reqs(j),ierr)
    CALL MPI_IRECV(lerecv(j)%re,nbemax,OFT_MPI_R8,self%estitch%proc_con(j), &
                   1,oft_env%COMM,self%estitch%recv_reqs(j),ierr)
  END DO
  CALL oft_mpi_waitall(self%estitch%nproc_con,self%estitch%recv_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
  CALL oft_mpi_waitall(self%estitch%nproc_con,self%estitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_elinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_elinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Condense linkage to sparse rep
!------------------------------------------------------------------------------
self%estitch%nle=SUM(self%estitch%kle)
self%estitch%kle(self%estitch%nproc_con+1)=self%estitch%nle+1
DO i=self%estitch%nproc_con,0,-1 ! cumulative unique edge linkage count
  self%estitch%kle(i)=self%estitch%kle(i+1)-self%estitch%kle(i)
END DO
IF(self%estitch%kle(0)/=1)CALL oft_abort('Bad edge linkage count','mesh_global_elinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(self%estitch%lle(2,self%estitch%nle))
!---Create seam object
self%estitch%full=self%fullmesh
self%estitch%nbe=self%nbe
self%estitch%lbe=>self%lbe
ALLOCATE(self%estitch%send(0:self%estitch%nproc_con),self%estitch%recv(0:self%estitch%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,self%estitch%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%estitch%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lerecv(i)%re(linktmp(1,j,i))>lesend(1)%re(m))THEN
        self%estitch%leo(m) = .FALSE.
      ELSE IF(lerecv(i)%re(linktmp(1,j,i))==lesend(1)%re(m).AND.self%estitch%proc_con(i)<oft_env%rank)THEN
        self%estitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    self%estitch%lle(:,m+self%estitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  self%estitch%send(i)%n=ncon(i)
  ALLOCATE(self%estitch%send(i)%v(self%estitch%send(i)%n))
  self%estitch%recv(i)%n=ncon(i)
  ALLOCATE(self%estitch%recv(i)%v(self%estitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!------------------------------------------------------------------------------
! Clean-up transfers
!------------------------------------------------------------------------------
DEALLOCATE(lesend(1)%le)
DEALLOCATE(lesend(1)%re)
IF(.NOT.self%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,self%estitch%nproc_con
    DEALLOCATE(lerecv(i)%le)
    DEALLOCATE(lerecv(i)%re)
  END DO
END IF
!---Check stitching information
CALL oft_stitch_check(self%estitch)
DEALLOCATE(linktmp,ncon,echeck)
DEALLOCATE(lesend,lerecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_elinkage
!---------------------------------------------------------------------------------
!> Construct processor to processor face linkage for stitching operations.
!! - Create linkage of boundary faces to other processors.
!! - Determine ownership for shared faces
!! - Set global boundary face flag
!---------------------------------------------------------------------------------
SUBROUTINE mesh_global_flinkage(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
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
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Face linkage'
CALL oft_increase_indent
IF(.NOT.self%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!------------------------------------------------------------------------------
! Determine maximum boundary face count
!------------------------------------------------------------------------------
nftmp=self%nbf
nbfmax=self%nbf
IF(.NOT.self%fullmesh)nbfmax=oft_mpi_max(nftmp)
IF(oft_debug_print(2))WRITE(*,'(2A,I8)')oft_indent,'Max # of seam faces =',nbfmax
self%fstitch%nbemax=nbfmax
!---Initialize seam structure from mesh
CALL oft_init_seam(self,self%fstitch)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Temporary linkage array
ALLOCATE(linktmp(2,self%nbf,0:self%fstitch%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:self%fstitch%nproc_con),fcheck(0:self%fstitch%nproc_con,self%nbf))
ncon=0
fcheck=.FALSE.
!---Pointer into linkage array
ALLOCATE(self%fstitch%kle(0:self%fstitch%nproc_con+1))
self%fstitch%kle = 0
!---Global boundary face tag
IF(ASSOCIATED(self%global%gbf))DEALLOCATE(self%global%gbf)
ALLOCATE(self%global%gbf(self%nf))
self%global%gbf = .FALSE.
DO i=1,self%nf
  IF(self%bf(i))self%global%gbf(i)=.TRUE.
END DO
!---Allocate temporary Send/Recv arrays
ALLOCATE(lfsend(self%fstitch%nproc_con+1),lfrecv(self%fstitch%nproc_con+1))
ALLOCATE(lfsend(1)%lf(nbfmax))
ALLOCATE(lfsend(1)%rf(nbfmax))
IF(.NOT.self%fullmesh)THEN
  DO i=1,self%fstitch%nproc_con
    ALLOCATE(lfrecv(i)%lf(nbfmax))
    ALLOCATE(lfrecv(i)%rf(nbfmax))
  END DO
ENDIF
!---Create link from local index to boundary index
ALLOCATE(bpi(self%np))
CALL get_inverse_map(self%lbp,self%nbp,bpi,self%np)
!---Count corner ownership for load balancing
CALL oft_random_number(lfsend(1)%rf,nbfmax)
DO i=1,self%nbf
  DO m=1,self%face_np
    IF(self%pstitch%leo(bpi(self%lf(m,self%lbf(i)))))lfsend(1)%rf(i) = lfsend(1)%rf(i) + 1.d0
  END DO
END DO
!---Point dummy output array to Send array
lfout=>lfsend(1)%lf
lfout=0
!$omp parallel do private(m,ll,js,jn,j,ftest)
DO i=1,self%nbf
  DO m=1,self%face_np
    ll(m)=bpi(self%lf(m,self%lbf(i))) ! Get endpoints in boundary index
  END DO
  lfout(i)=ABS(self%global%lf(self%lbf(i))) ! Populate output array with global face indices
  DO m=0,self%fstitch%nproc_con  ! Flag processors to be checked
    js=self%pstitch%kle(m)
    jn=self%pstitch%kle(m+1)-1
    ftest=.TRUE.
    DO j=1,self%face_np
      ftest=ftest.AND.ANY(self%pstitch%lle(1,js:jn)==ll(j))
    END DO
    IF(ftest)fcheck(m,i)=.TRUE.
  END DO
END DO
DEALLOCATE(bpi)
!---Point dummy Send arrays to main Send array
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  DO i=2,self%fstitch%nproc_con
    lfsend(i)%lf=>lfsend(1)%lf
    lfsend(i)%rf=>lfsend(1)%rf
  END DO
  !---Wait for all processes
  CALL oft_mpi_barrier(ierr)
!------------------------------------------------------------------------------
! Create Send and Recv calls
!------------------------------------------------------------------------------
  DO j=1,self%fstitch%nproc_con
    CALL MPI_ISEND(lfsend(j)%lf,nbfmax,OFT_MPI_I8,self%fstitch%proc_con(j),1,oft_env%COMM,self%fstitch%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','mesh_global_flinkage',__FILE__)
    CALL MPI_IRECV(lfrecv(j)%lf,nbfmax,OFT_MPI_I8,self%fstitch%proc_con(j),1,oft_env%COMM,self%fstitch%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','mesh_global_flinkage',__FILE__)
  END DO
!------------------------------------------------------------------------------
! Loop over each connected processor
!------------------------------------------------------------------------------
  DO WHILE(.TRUE.)
    IF(oft_mpi_check_reqs(self%fstitch%nproc_con,self%fstitch%recv_reqs))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(self%fstitch%nproc_con,self%fstitch%recv_reqs,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','mesh_global_flinkage',__FILE__)
    lftmp=>lfrecv(j)%lf ! Point dummy input array to current Recv array
    !---
    nffl=0
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,self%nbf ! Loop over boundary faces
      IF(.NOT.fcheck(j,m))CYCLE
      DO mm=1,nbfmax ! Loop over input faces
        IF(lfout(m)==lftmp(mm))THEN ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          self%global%gbf(self%lbf(m))=.FALSE. ! Face is not on global boundary
          nffl=nffl+1
          EXIT
        END IF
      END DO
    END DO
    self%fstitch%kle(j)=nffl
  END DO
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_flinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Check internal connections
!------------------------------------------------------------------------------
ALLOCATE(self%fstitch%leo(self%nbf))
self%fstitch%leo = .TRUE. ! Default to ownership
!---Determine location of periodic faces on current processor
nffl=0
IF(self%periodic%nper>0)THEN
  IF(ASSOCIATED(self%periodic%lf))THEN
    !---Create link from local index to boundary index
    ALLOCATE(bfi(self%nf))
    CALL get_inverse_map(self%lbf,self%nbf,bfi,self%nf)
    !---Find face pairs from periodic parents
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,self%nbf
      mm=self%lbf(m)
      IF(self%periodic%lf(mm)>0)THEN
        !$omp critical
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/bfi(self%periodic%lf(mm)),m/)
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/m,bfi(self%periodic%lf(mm))/)
        !$omp end critical
        nffl=nffl+2
        self%global%gbf(mm)=.FALSE.
        self%global%gbf(self%periodic%lf(mm))=.FALSE.
        IF(bfi(self%periodic%lf(mm))<m)THEN
          self%fstitch%leo(m) = .FALSE.
          lfsend(1)%rf(m)=-1.d0
        END IF
      END IF
    END DO
    DEALLOCATE(bfi)
  ELSE
    !$omp parallel do private(mm) reduction(+:nffl)
    DO m=1,self%nbf ! Loop over boundary faces
      IF(.NOT.fcheck(0,m))CYCLE
      DO mm=m+1,self%nbf ! Loop over boundary faces
        IF(lfout(m)==lfout(mm))THEN ! Found match
          !$omp critical
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/mm,m/)
          ncon(0)=ncon(0)+1
          linktmp(:,ncon(0),0)=(/m,mm/)
          !$omp end critical
          nffl=nffl+2
          self%global%gbf(self%lbf(m))=.FALSE. ! Face is not on global boundary
          self%global%gbf(self%lbf(mm))=.FALSE. ! Face is not on global boundary
          self%fstitch%leo(mm) = .FALSE.
          lfsend(1)%rf(mm)=-1.d0
          EXIT
        END IF
        IF(lfout(mm)==0)EXIT ! End of input points
      END DO
    END DO
  END IF
END IF
self%fstitch%kle(0)=nffl
!---Transfer load balancing arrays
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(self%fstitch%nproc_con,self%fstitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
  DO j=1,self%fstitch%nproc_con
    CALL MPI_ISEND(lfsend(j)%rf,nbfmax,OFT_MPI_R8,self%fstitch%proc_con(j), &
                   1,oft_env%COMM,self%fstitch%send_reqs(j),ierr)
    CALL MPI_IRECV(lfrecv(j)%rf,nbfmax,OFT_MPI_R8,self%fstitch%proc_con(j), &
                   1,oft_env%COMM,self%fstitch%recv_reqs(j),ierr)
  END DO
  CALL oft_mpi_waitall(self%fstitch%nproc_con,self%fstitch%recv_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
  CALL oft_mpi_waitall(self%fstitch%nproc_con,self%fstitch%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','mesh_global_flinkage',__FILE__)
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_flinkage",__FILE__)
#endif
END IF
!------------------------------------------------------------------------------
! Condense linkage to sparse rep
!------------------------------------------------------------------------------
self%fstitch%nle=sum(self%fstitch%kle)
self%fstitch%kle(self%fstitch%nproc_con+1)=self%fstitch%nle+1
DO i=self%fstitch%nproc_con,0,-1 ! cumulative unique edge linkage count
  self%fstitch%kle(i)=self%fstitch%kle(i+1)-self%fstitch%kle(i)
END DO
IF(self%fstitch%kle(0)/=1)CALL oft_abort('Bad face linkage count','mesh_global_flinkage',__FILE__)
!---Generate linkage lists
ALLOCATE(self%fstitch%lle(2,self%fstitch%nle))
!---
self%fstitch%full=self%fullmesh
! mesh%fstitch%nbemax=mesh%linkage%nbfmax
self%fstitch%nbe=self%nbf
self%fstitch%lbe=>self%lbf
! mesh%fstitch%nle=mesh%linkage%nlf
! mesh%fstitch%kle=>mesh%linkage%klf
! mesh%fstitch%lle=>mesh%linkage%llf
! mesh%fstitch%leo=>mesh%linkage%lfo
ALLOCATE(self%fstitch%send(0:self%fstitch%nproc_con),self%fstitch%recv(0:self%fstitch%nproc_con))
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
DO i=0,self%fstitch%nproc_con
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%fstitch%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
      IF(lfrecv(i)%rf(linktmp(1,j,i))>lfsend(1)%rf(m))THEN
        self%fstitch%leo(m) = .FALSE.
      ELSE IF(lfrecv(i)%rf(linktmp(1,j,i))==lfsend(1)%rf(m).AND.self%fstitch%proc_con(i)<oft_env%rank)THEN
        self%fstitch%leo(m) = .FALSE.
      END IF
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    self%fstitch%lle(:,m+self%fstitch%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  self%fstitch%send(i)%n=ncon(i)
  ALLOCATE(self%fstitch%send(i)%v(self%fstitch%send(i)%n))
  self%fstitch%recv(i)%n=ncon(i)
  ALLOCATE(self%fstitch%recv(i)%v(self%fstitch%recv(i)%n))
END DO
DEALLOCATE(lsort,isort)
!$omp end parallel
!------------------------------------------------------------------------------
! Clean-up transfers
!------------------------------------------------------------------------------
DEALLOCATE(lfsend(1)%lf)
DEALLOCATE(lfsend(1)%rf)
IF(.NOT.self%fullmesh)THEN
  !---Deallocate temporary work arrays
  DO i=1,self%fstitch%nproc_con
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
CALL oft_stitch_check(self%fstitch)
DEALLOCATE(linktmp,ncon,fcheck)
DEALLOCATE(lfsend,lfrecv)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_global_flinkage
!---------------------------------------------------------------------------------
!> Construct processor to processor point linkage for stitching operations.
!! - Create linkage of boundary points to other processors
!! - Determine ownership for shared points
!---------------------------------------------------------------------------------
SUBROUTINE plinkage_from_parent(self,parent)
class(oft_bmesh), intent(inout) :: self !< Mesh object
class(oft_mesh), intent(inout) :: parent !< Parent volume mesh
INTEGER(i4) :: i,j,k,m
INTEGER(i4), POINTER :: lpbound(:),lploc(:)
IF(self%np==0)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Point linkage'
ALLOCATE(lpbound(self%np))
CALL get_inverse_map(self%lbp,self%nbp,lpbound,self%np)
ALLOCATE(lploc(parent%np))
CALL get_inverse_map(self%parent%lp,self%np,lploc,self%parent%np)
!---
CALL oft_init_seam(self,self%pstitch)
ALLOCATE(self%pstitch%kle(self%pstitch%nproc_con+1))
self%pstitch%kle=0
DO i=1,self%pstitch%nproc_con
  m=0
  DO j=parent%pstitch%kle(i),parent%pstitch%kle(i+1)-1
    k=parent%lbp(parent%pstitch%lle(1,j))
    IF(lploc(k)/=0)THEN
      m=m+1
    END IF
  END DO
  self%pstitch%kle(i)=m
END DO
!------------------------------------------------------------------------------
! Condense linkage to sparse rep
!------------------------------------------------------------------------------
self%pstitch%nle=SUM(self%pstitch%kle)
self%pstitch%kle(self%pstitch%nproc_con+1)=self%pstitch%nle+1
!---Cumulative unique point linkage count
DO i=self%pstitch%nproc_con,1,-1
  self%pstitch%kle(i)=self%pstitch%kle(i+1)-self%pstitch%kle(i)
END DO
!---
IF(self%pstitch%kle(1)/=1)CALL oft_abort('Bad point linkage count','tetmesh_global_plinkage',__FILE__)
!---
ALLOCATE(self%pstitch%lle(2,self%pstitch%nle),self%pstitch%leo(parent%nbp))
ALLOCATE(self%pstitch%send(0:self%pstitch%nproc_con),self%pstitch%recv(0:self%pstitch%nproc_con))
!---
!!$omp do
DO i=1,self%pstitch%nproc_con
  m=0
  DO j=parent%pstitch%kle(i),parent%pstitch%kle(i+1)-1
    k=parent%lbp(parent%pstitch%lle(1,j))
    IF(lploc(k)/=0)THEN
      self%pstitch%lle(1,m+self%pstitch%kle(i))=lpbound(lploc(k))
      self%pstitch%lle(2,m+self%pstitch%kle(i))=m+1
      m=m+1
    END IF
  END DO
  !---Allocate permanent stitching arrays
  self%pstitch%send(i)%n=m; self%pstitch%recv(i)%n=m
  ALLOCATE(self%pstitch%send(i)%v(self%pstitch%send(i)%n))
  ALLOCATE(self%pstitch%recv(i)%v(self%pstitch%recv(i)%n))
END DO
!!$omp end parallel
!------------------------------------------------------------------------------
! Assemble ownership tag of local boundary points
!------------------------------------------------------------------------------
DEALLOCATE(lpbound)
ALLOCATE(lpbound(parent%np))
CALL get_inverse_map(parent%lbp,parent%nbp,lpbound,parent%np)
!$omp parallel do private(k)
DO i=1,self%nbp
  k=INT(self%parent%lp(self%lbp(i)),4)
  IF(lploc(k)/=0)THEN
    self%pstitch%leo(i) = parent%pstitch%leo(lpbound(k))
  END IF
END DO
!------------------------------------------------------------------------------
! Clean-up
!------------------------------------------------------------------------------
!---Deallocate temporary work arrays
DEALLOCATE(lpbound,lploc)
DEBUG_STACK_POP
END SUBROUTINE plinkage_from_parent
!---------------------------------------------------------------------------------
!> Construct processor to processor edge linkage for stitching operations.
!! - Create linkage of boundary edges to other processors with orientation.
!! - Determine ownership for shared edges
!---------------------------------------------------------------------------------
SUBROUTINE elinkage_from_parent(self,parent)
class(oft_bmesh), intent(inout) :: self !< Mesh object
class(oft_mesh), intent(inout) :: parent !< Parent volume mesh
INTEGER(i4) :: i,j,k,m
INTEGER(i4), POINTER :: leloc(:),lebound(:)
IF(self%ne==0)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Edge linkage'
ALLOCATE(lebound(self%ne))
CALL get_inverse_map(self%lbe,self%nbe,lebound,self%ne)
ALLOCATE(leloc(parent%ne))
CALL get_inverse_map(self%parent%le,self%ne,leloc,self%parent%ne)
!---
CALL oft_init_seam(self,self%estitch)
ALLOCATE(self%estitch%kle(self%estitch%nproc_con+1))
self%estitch%kle=0
DO i=1,self%estitch%nproc_con
  m=0
  DO j=parent%estitch%kle(i),parent%estitch%kle(i+1)-1
    k=parent%lbe(parent%estitch%lle(1,j))
    IF(leloc(k)/=0)THEN
      m=m+1
    END IF
  END DO
  self%estitch%kle(i)=m
END DO
!------------------------------------------------------------------------------
! Condense linkage to sparse rep
!------------------------------------------------------------------------------
self%estitch%nle=SUM(self%estitch%kle)
self%estitch%kle(self%estitch%nproc_con+1)=self%estitch%nle+1
!---Cumulative unique point linkage count
DO i=self%estitch%nproc_con,1,-1
  self%estitch%kle(i)=self%estitch%kle(i+1)-self%estitch%kle(i)
END DO
!---
IF(self%estitch%kle(1)/=1)CALL oft_abort('Bad point linkage count','tetmesh_global_plinkage',__FILE__)
!---
ALLOCATE(self%estitch%lle(2,self%estitch%nle),self%estitch%leo(parent%nbe))
ALLOCATE(self%estitch%send(0:self%estitch%nproc_con),self%estitch%recv(0:self%estitch%nproc_con))
!---
!!$omp do
DO i=1,self%estitch%nproc_con
  m=0
  DO j=parent%estitch%kle(i),parent%estitch%kle(i+1)-1
    k=parent%lbe(parent%estitch%lle(1,j))
    IF(leloc(k)/=0)THEN
      self%estitch%lle(1,m+self%estitch%kle(i))=lebound(leloc(k))
      self%estitch%lle(2,m+self%estitch%kle(i))=m+1
      m=m+1
    END IF
  END DO
  !---Allocate permanent stitching arrays
  self%estitch%send(i)%n=m; self%estitch%recv(i)%n=m
  ALLOCATE(self%estitch%send(i)%v(self%estitch%send(i)%n))
  ALLOCATE(self%estitch%recv(i)%v(self%estitch%recv(i)%n))
END DO
!!$omp end parallel
!------------------------------------------------------------------------------
! Assemble ownership tag of local boundary points
!------------------------------------------------------------------------------
DEALLOCATE(lebound)
ALLOCATE(lebound(parent%ne))
CALL get_inverse_map(parent%lbe,parent%nbe,lebound,parent%ne)
!$omp parallel do private(k)
DO i=1,self%nbe
  k=INT(ABS(self%parent%le(self%lbe(i))),4)
  IF(leloc(k)/=0)THEN
    self%estitch%leo(i) = parent%estitch%leo(lebound(k))
  END IF
END DO
!------------------------------------------------------------------------------
! Clean-up
!------------------------------------------------------------------------------
!---Deallocate temporary work arrays
DEALLOCATE(lebound,leloc)
DEBUG_STACK_POP
END SUBROUTINE elinkage_from_parent
end module oft_mesh_global
