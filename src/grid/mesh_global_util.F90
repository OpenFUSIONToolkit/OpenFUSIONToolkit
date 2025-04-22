!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_mesh_global_util.F90
!
!> Global context creation and utility functions for distributed meshes.
!! - Global boundary tagging
!! - Global element orientation
!!   - Compute global orientations
!!   - Orient element
!! - Global I/O indexing
!! - Ground point location
!! - Global mesh resolution and statistics
!!
!! @author Chris Hansen
!! @date Spring 2010
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_mesh_global_util
USE oft_base
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh, oft_bmesh
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_stitching, ONLY: oft_global_stitch
IMPLICIT NONE
#include "local.h"
contains
!---------------------------------------------------------------------------------
!> Set global boundary flags for distributed meshes.
!---------------------------------------------------------------------------------
subroutine mesh_global_boundary(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,k
real(r8), allocatable, dimension(:) :: btrans
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Constructing global boundary'
!------------------------------------------------------------------------------
! Set global boundary points and cells
!------------------------------------------------------------------------------
IF(ASSOCIATED(mesh%global%gbc))DEALLOCATE(mesh%global%gbc)
allocate(mesh%global%gbc(mesh%nc)) ! Allocate boundary cell flag
mesh%global%gbc = .FALSE. ! Initiliaze boundary flag
allocate(btrans(mesh%np)) ! Allocate Send/Recv refernce arrays
btrans=0.d0
!---Mark global boundary points and cells on local processor
do i=1,mesh%nbf
  if(mesh%global%gbf(mesh%lbf(i)))then ! Face is on global boundary (Root Element)
    do j=1,mesh%face_np ! Loop over corner points
      k=mesh%lf(j,mesh%lbf(i))
      mesh%global%gbp(k)=.TRUE. ! Point is a boundary point
      btrans(k)=1.d0
    enddo
    mesh%global%gbc(mesh%lfc(1,mesh%lbf(i)))=.TRUE. ! Shared cell is a boundary cell
  endif
enddo
!---Synchronize across seams
CALL oft_global_stitch(mesh%pstitch,btrans,1) ! Sync boundary points across processors
do i=1,mesh%nbp
  if(btrans(mesh%lbp(i))>0.d0)mesh%global%gbp(mesh%lbp(i))=.TRUE. ! Point is a global boundary point
enddo
deallocate(btrans)
!------------------------------------------------------------------------------
! Set global boundary edges
!------------------------------------------------------------------------------
allocate(btrans(mesh%ne))
btrans=0.d0
!---Mark global boundary edges on local processor
do i=1,mesh%nbf
  if(mesh%global%gbf(mesh%lbf(i)))then ! Face is on global boundary (Root Element)
    do j=1,mesh%face_np
      k=abs(mesh%lfe(j,mesh%lbf(i)))
      mesh%global%gbe(k)=.TRUE. ! Edge is a boundary edge
      btrans(k)=1.d0
    enddo
  endif
enddo
!---Synchronize across seams
CALL oft_global_stitch(mesh%estitch,btrans,1) ! Sync boundary edges across processors
do i=1,mesh%nbe
  if(btrans(mesh%lbe(i))>0.d0)mesh%global%gbe(mesh%lbe(i))=.TRUE. ! Edge is a global boundary edge
enddo
deallocate(btrans)
DEBUG_STACK_POP
end subroutine mesh_global_boundary
!---------------------------------------------------------------------------------
!> Set global indexing from mesh periodicity
!---------------------------------------------------------------------------------
subroutine mesh_global_periodic(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,iper
INTEGER(8) :: npp,nep,nfp
INTEGER(8), ALLOCATABLE :: gtmp(:)
IF(mesh%periodic%nper==0)RETURN
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Setting Up Periodic Mesh'
!---Faces can only have one periodic parent
DO i=1,mesh%nbf
  j=mesh%lbf(i)
  IF(mesh%periodic%lf(j)>0)THEN
    mesh%global%lf(j)=mesh%global%lf(mesh%periodic%lf(j))
  END IF
END DO
!---Propogate edge and point parents through periodicity
DO iper=1,mesh%periodic%nper
  DO i=1,mesh%nbp
    j=mesh%lbp(i)
    IF(mesh%periodic%lp(j)>0)THEN
      mesh%global%lp(j)=mesh%global%lp(mesh%periodic%lp(j))
    END IF
  END DO
  !---
  DO i=1,mesh%nbe
    j=mesh%lbe(i)
    IF(mesh%periodic%le(j)>0)THEN
      mesh%global%le(j)=mesh%global%le(mesh%periodic%le(j))
    END IF
  END DO
END DO
!---Reindex globals
ALLOCATE(gtmp(mesh%np))
gtmp=-1
j=0
DO i=1,mesh%np
  IF(gtmp(mesh%global%lp(i))==-1)THEN
    j=j+1
    gtmp(mesh%global%lp(i))=j
  END IF
END DO
npp=mesh%global%np-j
mesh%global%np=j
DO i=1,mesh%np
  mesh%global%lp(i)=gtmp(mesh%global%lp(i))
END DO
DEALLOCATE(gtmp)
!---
ALLOCATE(gtmp(mesh%ne))
gtmp=-1
j=0
DO i=1,mesh%ne
  IF(gtmp(mesh%global%le(i))==-1)THEN
    j=j+1
    gtmp(mesh%global%le(i))=j
  END IF
END DO
nep=mesh%global%ne-j
mesh%global%ne=j
DO i=1,mesh%ne
  mesh%global%le(i)=gtmp(mesh%global%le(i))
END DO
DEALLOCATE(gtmp)
!---
ALLOCATE(gtmp(mesh%nf))
gtmp=-1
j=0
DO i=1,mesh%nf
  IF(gtmp(mesh%global%lf(i))==-1)THEN
    j=j+1
    gtmp(mesh%global%lf(i))=j
  END IF
END DO
nfp=mesh%global%nf-j
mesh%global%nf=j
DO i=1,mesh%nf
  mesh%global%lf(i)=gtmp(mesh%global%lf(i))
END DO
DEALLOCATE(gtmp)
!---
if(oft_debug_print(1))then
  CALL oft_increase_indent
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic directions =',mesh%periodic%nper
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Points     =',npp
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Edges      =',nep
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Faces      =',nfp
  CALL oft_decrease_indent
end if
DEBUG_STACK_POP
end subroutine mesh_global_periodic
!---------------------------------------------------------------------------------
!> Set global orientations for edges and faces
!---------------------------------------------------------------------------------
subroutine mesh_global_orient(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,k,ind,ed(2)
integer(i4), ALLOCATABLE, DIMENSION(:) :: ltmp
integer(8) :: etmp(2)
integer(8), ALLOCATABLE, DIMENSION(:) :: kc,kf
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Orienting geometry'
!$omp parallel do private(etmp,ed)
do i=1,mesh%ne
  ed=mesh%le(:,i)
  etmp=mesh%global%lp(ed)
  if(mesh%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
    IF(ed(1)>ed(2))THEN
      etmp(1)=etmp(1)+mesh%global%np
    ELSE
      etmp(2)=etmp(2)+mesh%global%np
    END IF
  end if
  if(etmp(1)>etmp(2))then
    mesh%global%le(i)=-ABS(mesh%global%le(i))
  end if
end do
allocate(mesh%lfo(mesh%nf))
! !$omp parallel do private(kf,j,etmp,ftmp)
allocate(kf(mesh%face_np))
allocate(ltmp(mesh%face_np))
do i=1,mesh%nf
  kf=mesh%global%lp(mesh%lf(:,i))
  do j=1,mesh%face_np
    ed=mesh%face_ed(:,j)
    etmp=kf(ed) !(/kf(mesh%bmesh%face_ed(1,j)),kf(mesh%bmesh%face_ed(2,j))/)
    if(mesh%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
      IF(mesh%lf(ed(1),i)>mesh%lf(ed(2),i))THEN
        etmp(1)=etmp(1)+mesh%global%np
      ELSE
        etmp(2)=etmp(2)+mesh%global%np
      END IF
    end if
    if(etmp(1)>etmp(2))then
      mesh%lfe(j,i)=-ABS(mesh%lfe(j,i))
    else
      mesh%lfe(j,i)=ABS(mesh%lfe(j,i))
    end if
  end do
  ltmp=INT(kf,4)
  CALL find_orient_listn(mesh%lfo(i),ltmp,mesh%face_np)
end do
allocate(kc(mesh%cell_np))
! !$omp parallel do private(kc,j,etmp)
do i=1,mesh%nc
  kc=mesh%global%lp(mesh%lc(:,i))
  do j=1,mesh%cell_ne
    ed=mesh%cell_ed(:,j)
    etmp=kc(ed) !(/kc(mesh%cell_ed(1,j)),kc(mesh%cell_ed(2,j))/)
    if(mesh%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
      IF(mesh%lc(ed(1),i)>mesh%lc(ed(2),i))THEN
        etmp(1)=etmp(1)+mesh%global%np
      ELSE
        etmp(2)=etmp(2)+mesh%global%np
      END IF
    end if
    if(etmp(1)>etmp(2))then
      mesh%lce(j,i)=-ABS(mesh%lce(j,i))
    else
      mesh%lce(j,i)=ABS(mesh%lce(j,i))
    end if
  end do
end do
allocate(mesh%lcfo(mesh%cell_nf,mesh%nc))
! !$omp parallel do private(kc,j,ftmp)
do i=1,mesh%nc
  kc=mesh%global%lp(mesh%lc(:,i))
  ! do j=1,mesh%cell_ne
  !   ed=mesh%cell_ed(:,j)
  !   if((mesh%periodic%nper>0).AND.(kc(ed(1))==kc(ed(2))))then
  !     IF(mesh%lc(ed(1),i)>mesh%lc(ed(2),i))THEN
  !       kc(ed(1))=kc(ed(1))+mesh%global%np
  !     ELSE
  !       kc(ed(2))=kc(ed(2))+mesh%global%np
  !     END IF
  !   end if
  ! end do
  do j=1,mesh%cell_nf
    kf=kc(mesh%cell_fc(:,j))
    do k=1,mesh%face_np
      ed=(/k,k+1/)
      IF(k==mesh%face_np)ed(2)=1
      if((mesh%periodic%nper>0).AND.(kf(ed(1))==kf(ed(2))))then
        IF(mesh%lc(mesh%cell_fc(ed(1),j),i)>mesh%lc(mesh%cell_fc(ed(2),j),i))THEN
          kf(ed(1))=kf(ed(1))+mesh%global%np
        ELSE
          kf(ed(2))=kf(ed(2))+mesh%global%np
        END IF
      end if
    end do
    ltmp=INT(kf,4)
    CALL find_orient_listn(mesh%lcfo(j,i),ltmp,mesh%face_np)
  end do
end do
deallocate(kc,kf,ltmp)
DEBUG_STACK_POP
end subroutine mesh_global_orient
!---------------------------------------------------------------------------------
!> Set cell curvature flag for common cases
!---------------------------------------------------------------------------------
subroutine mesh_global_set_curved(self,flag)
CLASS(oft_amesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: flag !< How to set curvature (1 -> all cells, 2 -> contains boundary edge)
INTEGER(i4) :: i,k
SELECT CASE(flag)
CASE(1)
  self%ho_info%is_curved=.TRUE.
CASE(2)
  self%ho_info%is_curved=.FALSE.
  DO i=1,self%nbe
    DO k=self%kec(self%lbe(i)),self%kec(self%lbe(i)+1)-1
      self%ho_info%is_curved(self%lec(k))=.TRUE.
    END DO
  END DO
CASE DEFAULT
  CALL oft_abort("Unknown setup flag","mesh_global_set_curved",__FILE__)
END SELECT
end subroutine mesh_global_set_curved
!---------------------------------------------------------------------------------
!> Get I/O scope for current mesh
!! - Element offsets for each geometric type
!! - Element span for each geometric type
!---------------------------------------------------------------------------------
subroutine mesh_global_save(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(8) :: nloc(6),ntmp(4),nfb,nfl,nfg,nfmax
integer(i4) :: i,j,k,l
#ifdef HAVE_MPI
integer(i4) :: ierr
#ifdef OFT_MPI_F08
type(mpi_status) :: stat
#else
integer(i4) :: stat(MPI_STATUS_SIZE)
#endif
#endif
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting I/O information'
!---Initialize counters for ring update
j=oft_env%rank+1
if(j>oft_env%nprocs-1)j=0
k=oft_env%rank-1
if(k<0)k=oft_env%nprocs-1
l=k
!---Initialize counts
nloc=0
!---Local count of owned boundary points
!$omp parallel do reduction(+:nloc)
do i=1,self%nbp ! Loop over boundary points
  if(self%pstitch%leo(i))nloc(4) = nloc(4) + 1 ! Processor owns this point
enddo
!$omp parallel do reduction(+:nloc)
do i=1,self%np
  if(self%bp(i))cycle
  nloc(4) = nloc(4) + 1 ! Interior points (all owned)
enddo
!---Local count of owned boundary edges
!$omp parallel do reduction(+:nloc)
do i=1,self%nbe ! Loop over boundary edges
  if(self%estitch%leo(i))THEN
    nloc(1) = nloc(1) + 1 ! Processor owns this edge
    IF(self%global%gbe(self%lbe(i)))nloc(6)=nloc(6)+1
  END IF
enddo
!$omp parallel do reduction(+:nloc)
do i=1,self%ne
  if(self%be(i))cycle
  nloc(1) = nloc(1) + 1 ! Interior edges (all owned)
enddo
nloc(3)=INT(self%nc,8) ! All cells are owned
!---Handle faces for volume meshes
SELECT TYPE(self)
CLASS IS(oft_mesh)
  !---Local count of owned boundary faces
  !$omp parallel do reduction(+:nloc)
  do i=1,self%nbf ! Loop over boundary faces
    if(self%fstitch%leo(i))nloc(2)=nloc(2)+1 ! Processor owns this face
  enddo
  !$omp parallel do reduction(+:nloc)
  do i=1,self%nf
    if(self%bf(i))cycle
    nloc(2)=nloc(2)+1 ! Interior face (all owned)
  enddo
  !---Local count of global boundary faces
  !$omp parallel do reduction(+:nloc)
  do i=1,self%nbf ! Loop over boundary faces
    if(self%global%gbf(self%lbf(i)))nloc(5)=nloc(5)+1 ! Face is global boundary
  enddo
  self%save%nbf=nloc(5)
END SELECT
!---Get maximum element counts over every processor
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(nloc(1:4),ntmp,4,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
#else
ntmp=nloc(1:4)
#endif
self%save%npmax=ntmp(4)
self%save%nemax=ntmp(1)
! mesh%save%nfmax=ntmp(2)
nfmax=ntmp(2)
self%save%ncmax=ntmp(3)
!---Initialize ouput index counters
self%save%npb=0
self%save%neb=0
! mesh%save%nfb=0
nfb=0
self%save%ncb=0
self%save%npl=nloc(4)
self%save%nel=nloc(1)
! mesh%save%nfl=nloc(2)
nfl=nloc(2)
self%global%ne=nloc(1)
! mesh%global%nf=nloc(2)
nfg=nloc(2)
! mesh%bmesh%global%ne=nloc(6)
!---Loop over processors to determine output index counts
IF(.NOT.self%fullmesh)THEN
#ifdef HAVE_MPI
  do while(l.NE.oft_env%rank)
    call MPI_SEND(nloc,6,OFT_MPI_I8,j,1,oft_env%COMM,ierr)
    call MPI_RECV(nloc,6,OFT_MPI_I8,k,1,oft_env%COMM,stat,ierr)

    if(l<oft_env%rank)then
      self%save%neb = self%save%neb + nloc(1) ! Update first edge index
      nfb = nfb + nloc(2) ! Update first face index
      self%save%ncb = self%save%ncb + nloc(3) ! Update first cell index
      self%save%npb = self%save%npb + nloc(4) ! Update first point index
    endif
    self%global%ne = self%global%ne + nloc(1) ! Update global edge count
    nfg = nfg + nloc(2) ! Update global face count
    ! mesh%bmesh%global%ne = mesh%bmesh%global%ne + nloc(6)
    !---Increment ring comm
    l=l-1
    if(l<0)l=oft_env%nprocs-1
  enddo
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_save",__FILE__)
#endif
END IF
self%save%npl = self%save%npl + self%save%npb ! Set last point index
self%save%nel = self%save%nel + self%save%neb ! Set last edge index
!---Handle faces for volume meshes
SELECT TYPE(self)
CLASS IS(oft_mesh)
  self%save%nfmax = nfmax
  self%save%nfb = nfb       
  self%save%nfl = nfl + nfb ! Set last face index
  self%global%nf = nfg  
END SELECT
DEBUG_STACK_POP
end subroutine mesh_global_save
!---------------------------------------------------------------------------------
!> Determine location of grounding point on mesh.
!---------------------------------------------------------------------------------
SUBROUTINE mesh_global_igrnd(self,grnd_pt)
class(oft_amesh), intent(inout) :: self !< Mesh object
real(r8), optional, intent(in) :: grnd_pt(3)
REAL(r8) :: c,r2min,r2,rgrnd(3)
INTEGER(i4) :: i,j,k,l,ierr
INTEGER(i8) :: grnd_global
LOGICAL :: owned
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
type(mpi_status) :: stat
#else
integer(i4) :: stat(MPI_STATUS_SIZE)
#endif
#endif
DEBUG_STACK_PUSH
rgrnd=[1.d0,0.d0,0.d0]
IF(PRESENT(grnd_pt))rgrnd=grnd_pt
!---Locate local surface ground point
self%igrnd=-1
r2min=1.d99
DO i=1,self%nbp
  IF(.NOT.self%global%gbp(self%lbp(i)))CYCLE
  r2=sum((self%r(:,self%lbp(i))-rgrnd)**2)
  IF(r2<r2min)THEN
    r2min=r2
    self%igrnd(1)=self%lbp(i)
  END IF
END DO
IF(self%fullmesh)THEN
  !---Check periodicity
  IF(self%igrnd(1)==-1)THEN
    IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'No grounding point found: self is fully periodic!'
  ELSE
    grnd_global=self%global%lp(self%igrnd(1))
    self%igrnd=-1
    DO i=1,self%nbp
      IF(self%global%lp(self%lbp(i))==grnd_global)THEN
        IF(ALL(self%igrnd/=-1))CALL oft_abort('Ground point repeated on same domain', &
          'mesh_global_igrnd',__FILE__)
        IF(self%igrnd(1)==-1)THEN
          self%igrnd(1)=self%lbp(i)
        ELSE
          self%igrnd(2)=self%lbp(i)
        END IF
      END IF
    END DO
    IF(oft_env%head_proc)WRITE(*,'(2A,I8)')oft_indent,'Surface grounded at vertex',self%igrnd(1)
  END IF
  DEBUG_STACK_POP
  RETURN
END IF
!---Initialize counters for ring update
j=oft_env%rank+1
IF(j>oft_env%nprocs-1)j=0
k=oft_env%rank-1
IF(k<0)k=oft_env%nprocs-1
l=k
CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Loop over processors to find global igrnd
c=r2min
#ifdef HAVE_MPI
DO WHILE(l.NE.oft_env%rank)
  CALL MPI_SEND(c,1,OFT_MPI_R8,j,1,oft_env%COMM,ierr)
  CALL MPI_RECV(c,1,OFT_MPI_R8,k,1,oft_env%COMM,stat,ierr)
  !---Test if igrnd is on local mesh and update
  IF(c<r2min)self%igrnd=-1
  !---Increment ring comm
  l=l-1
  IF(l<0)l=oft_env%nprocs-1
END DO
#endif
!---Handle periodicity
grnd_global=-1
IF(self%igrnd(1)/=-1)grnd_global=self%global%lp(self%igrnd(1))
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,grnd_global,1,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
#endif
self%igrnd=-1
IF(grnd_global==-1)THEN
  IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'No grounding point found: self is fully periodic!'
ELSE
  owned=.FALSE.
  DO i=1,self%nbp
    IF(self%global%lp(self%lbp(i))==grnd_global)THEN
      IF(ALL(self%igrnd/=-1))CALL oft_abort('Ground point repeated on same domain', &
        'mesh_global_igrnd',__FILE__)
      IF(self%igrnd(1)==-1)THEN
        self%igrnd(1)=self%lbp(i)
        owned=self%pstitch%leo(i)
      ELSE
        self%igrnd(2)=self%lbp(i)
        owned=owned.OR.self%pstitch%leo(i)
      END IF
    END IF
  END DO
  !---Write out grounding information
  IF(oft_env%head_proc)WRITE(*,'(2A,I8)')oft_indent,'Surface grounded at vertex ',grnd_global
END IF
DEBUG_STACK_POP
END SUBROUTINE mesh_global_igrnd
!---------------------------------------------------------------------------------
!> Compute mesh resolution statistics.
!---------------------------------------------------------------------------------
subroutine mesh_global_resolution(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
real(r8) :: hmax,hrms,a,b,c,hmin
real(r8), ALLOCATABLE, DIMENSION(:) :: ed_lens
integer(i4) :: i,j,ierr
DEBUG_STACK_PUSH
hmin=huge(hmin)
hmax=0.d0
hrms=0.d0
!$omp parallel private(j,ed_lens) &
!$omp reduction(min:hmin) reduction(max:hmax) reduction(+:hrms)
ALLOCATE(ed_lens(self%cell_ne))
!$omp do
do i=1,self%nc ! loop over all cells
  do j=1,self%cell_ne
    ed_lens(j)=sum((self%r(:,self%lc(self%cell_ed(2,j),i)) &
      - self%r(:,self%lc(self%cell_ed(1,j),i)))**2)
  end do
  hmin = min(hmin,minval(ed_lens))
  hmax = max(hmax,maxval(ed_lens))
  hrms = hrms+sum(ed_lens)
enddo
deallocate(ed_lens)
!$omp end parallel

if(self%fullmesh)then
  self%hmin=sqrt(hmin)
  self%hmax=sqrt(hmax)
  self%hrms=sqrt(hrms/REAL(self%cell_ne*self%nc,8))
else
#ifdef HAVE_MPI
  call MPI_ALLREDUCE(hmin,a,1,OFT_MPI_R8,MPI_MIN,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(hmax,b,1,OFT_MPI_R8,MPI_MAX,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(hrms,c,1,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  self%hmin=sqrt(a)
  self%hmax=sqrt(b)
SELECT TYPE(self)
CLASS IS(oft_mesh)
  self%hrms=sqrt(c/REAL(self%cell_ne*self%global%nc,8))
CLASS IS(oft_bmesh)
  self%hrms=sqrt(c/REAL(self%cell_ne*self%global%nc,8))
END SELECT
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_resolution",__FILE__)
#endif
endif
DEBUG_STACK_POP
end subroutine mesh_global_resolution
!---------------------------------------------------------------------------------
!> Display mesh statistics and counts.
!---------------------------------------------------------------------------------
subroutine mesh_global_stats(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i8) :: a,i,tmp(4)
integer(i4) :: nbtmp,ierr
real(r8) :: tmpr(2),vol,area
DEBUG_STACK_PUSH
call mesh_global_resolution(self)
vol = self%volume()
area = self%bmesh%area()
if(self%fullmesh)then
  IF(self%periodic%nper>0)THEN
    self%global%nbp=COUNT(self%global%gbp(self%lbp).AND.self%pstitch%leo)
    self%global%nbe=COUNT(self%global%gbe(self%lbe).AND.self%estitch%leo)
    self%global%nbf=COUNT(self%global%gbf(self%lbf).AND.self%fstitch%leo)
    self%global%nbc=COUNT(self%global%gbc)
  ELSE
    self%global%nbp=COUNT(self%global%gbp)
    self%global%nbe=COUNT(self%global%gbe)
    self%global%nbf=COUNT(self%global%gbf)
    self%global%nbc=COUNT(self%global%gbc)
  END IF
  if(oft_env%head_proc)then
    write(*,'(2A)')oft_indent,'Mesh statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'Volume          =',vol
    write(*,'(2A,ES11.3)')oft_indent,'Surface area    =',area
    write(*,'(2A,I8)')    oft_indent,'# of points     =',self%global%np
    write(*,'(2A,I8)')    oft_indent,'# of edges      =',self%global%ne
    write(*,'(2A,I8)')    oft_indent,'# of faces      =',self%global%nf
    write(*,'(2A,I8)')    oft_indent,'# of cells      =',self%global%nc
    write(*,'(2A,I8)')    oft_indent,'# of boundary points =',self%global%nbp
    write(*,'(2A,I8)')    oft_indent,'# of boundary edges  =',self%global%nbe
    write(*,'(2A,I8)')    oft_indent,'# of boundary faces  =',self%global%nbf
    write(*,'(2A,I8)')    oft_indent,'# of boundary cells  =',self%global%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',self%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',self%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',self%hmax
    CALL oft_decrease_indent
  endif
else
#ifdef HAVE_MPI
  if(oft_env%head_proc)then
    write(*,'(2A)')oft_indent,'Mesh statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'Volume          =',vol
    write(*,'(2A,ES11.3)')oft_indent,'Surface area    =',area
    write(*,'(2A,I8)')    oft_indent,'# of points     =',self%global%np
    write(*,'(2A,I8)')    oft_indent,'# of edges      =',self%global%ne
    write(*,'(2A,I8)')    oft_indent,'# of faces      =',self%global%nf
    write(*,'(2A,I8)')    oft_indent,'# of cells      =',self%global%nc
  endif
  !---Count boundary points
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbp
    if(self%global%gbp(self%lbp(i)).AND.self%pstitch%leo(i))a=a+1
  enddo
  tmp(1)=a
  !---Count boundary edges
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbe
    if(self%global%gbe(self%lbe(i)).AND.self%estitch%leo(i))a=a+1
  enddo
  tmp(2)=a
  !---Count boundary faces
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbf
    if(self%global%gbf(self%lbf(i)).AND.self%fstitch%leo(i))a=a+1
  enddo
  tmp(3)=a
  !---Count boundary cells
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbc
    if(self%global%gbc(self%lbc(i)))a=a+1
  enddo
  tmp(4)=a
  !---Check volumes
  tmpr(1)=SUM(self%cv)
  tmpr(2)=SUM(self%vv)
  !---Get global sums
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,4,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmpr,2,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  self%global%nbp=tmp(1)
  self%global%nbe=tmp(2)
  self%global%nbf=tmp(3)
  self%global%nbc=tmp(4)
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',self%global%nbp
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',self%global%nbe
    write(*,'(2A,I8)')oft_indent,'# of boundary faces  =',self%global%nbf
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',self%global%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',self%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',self%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',self%hmax
    CALL oft_decrease_indent
  endif
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_stats",__FILE__)
#endif
endif
DEBUG_STACK_POP
end subroutine mesh_global_stats
!---------------------------------------------------------------------------------
!> Set global boundary flags for distributed meshes.
!---------------------------------------------------------------------------------
subroutine bmesh_global_boundary(self)
class(oft_bmesh), intent(inout) :: self !< Mesh object
integer(i4) :: i,j,k
real(r8), allocatable, dimension(:) :: btrans
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Constructing global boundary'
!------------------------------------------------------------------------------
! Set global boundary points and cells
!------------------------------------------------------------------------------
IF(ASSOCIATED(self%global%gbc))DEALLOCATE(self%global%gbc)
allocate(self%global%gbc(self%nc)) ! Allocate boundary cell flag
self%global%gbc = .FALSE. ! Initiliaze boundary flag
allocate(btrans(self%np)) ! Allocate Send/Recv refernce arrays
btrans=0.d0
!---Mark global boundary points and cells on local processor
do i=1,self%nbe
  if(self%global%gbe(self%lbe(i)))then ! Face is on global boundary (Root Element)
    do j=1,2 ! Loop over edge end points
      k=self%le(j,self%lbe(i))
      self%global%gbp(k)=.TRUE. ! Point is a boundary point
      btrans(k)=1.d0
    enddo
    self%global%gbc(self%lec(self%kec(self%lbe(i))))=.TRUE. ! Shared cell is a boundary cell
  endif
enddo
!---Synchronize across seams
CALL oft_global_stitch(self%pstitch,btrans,1) ! Sync boundary points across processors
do i=1,self%nbp
  if(btrans(self%lbp(i))>0.d0)self%global%gbp(self%lbp(i))=.TRUE. ! Point is a global boundary point
enddo
deallocate(btrans)
DEBUG_STACK_POP
end subroutine bmesh_global_boundary
!---------------------------------------------------------------------------------
!> Set global indexing from mesh periodicity
!---------------------------------------------------------------------------------
subroutine bmesh_global_periodic(self)
class(oft_bmesh), intent(inout) :: self !< Mesh object
integer(i4) :: i,j,iper
INTEGER(8) :: npp,nep,nfp
INTEGER(8), ALLOCATABLE :: gtmp(:)
IF(self%periodic%nper==0)RETURN
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Setting Up Periodic Mesh'
!---Edges can only have one periodic parent
DO i=1,self%nbe
  j=self%lbe(i)
  IF(self%periodic%le(j)>0)THEN
    self%global%le(j)=self%global%le(self%periodic%le(j))
  END IF
END DO
!---Propogate point parents through periodicity
DO iper=1,self%periodic%nper
  DO i=1,self%nbp
    j=self%lbp(i)
    IF(self%periodic%lp(j)>0)THEN
      self%global%lp(j)=self%global%lp(self%periodic%lp(j))
    END IF
  END DO
END DO
!---Reindex globals
ALLOCATE(gtmp(self%np))
gtmp=-1
j=0
DO i=1,self%np
  IF(gtmp(self%global%lp(i))==-1)THEN
    j=j+1
    gtmp(self%global%lp(i))=j
  END IF
END DO
npp=self%global%np-j
self%global%np=j
DO i=1,self%np
  self%global%lp(i)=gtmp(self%global%lp(i))
END DO
DEALLOCATE(gtmp)
!---
ALLOCATE(gtmp(self%ne))
gtmp=-1
j=0
DO i=1,self%ne
  IF(gtmp(self%global%le(i))==-1)THEN
    j=j+1
    gtmp(self%global%le(i))=j
  END IF
END DO
nep=self%global%ne-j
self%global%ne=j
DO i=1,self%ne
  self%global%le(i)=gtmp(self%global%le(i))
END DO
DEALLOCATE(gtmp)
!---
if(oft_debug_print(1))then
  CALL oft_increase_indent
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic directions =',self%periodic%nper
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Points     =',npp
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Edges      =',nep
  CALL oft_decrease_indent
end if
DEBUG_STACK_POP
end subroutine bmesh_global_periodic
!---------------------------------------------------------------------------------
!> Set global orientations for edges and faces
!---------------------------------------------------------------------------------
subroutine bmesh_global_orient(self,parent)
CLASS(oft_bmesh), intent(inout) :: self !< Mesh object
CLASS(oft_mesh), optional, intent(in) :: parent !< Parent volume mesh (if present)
integer(i4) :: i,j,ind,ed(2)
integer(i4), ALLOCATABLE, DIMENSION(:) :: kf,forder,ltmp
integer(8) :: etmp(2)
integer(8), ALLOCATABLE, DIMENSION(:) :: kfg
DEBUG_STACK_PUSH
IF(PRESENT(parent))THEN
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Copying orientation to surface'
  !$omp parallel do private(etmp)
  do i=1,self%ne
    etmp=self%le(:,i)
    etmp=self%parent%lp(etmp)
    etmp=parent%global%lp(etmp)
    if(etmp(1)>etmp(2))then
      self%global%le(i)=-ABS(self%global%le(i))
    end if
  end do
  allocate(self%lco(self%nc))
  allocate(kf(self%cell_np),kfg(self%cell_np),forder(self%cell_np))
  ! !$omp parallel do private(kf,j,etmp,forder)
  do i=1,self%nc
    kfg=self%lc(:,i)
    kfg=self%parent%lp(kfg)
    kfg=parent%global%lp(kfg)
    do j=1,self%cell_np
      etmp=(/kfg(self%cell_ed(1,j)),kfg(self%cell_ed(2,j))/)
      if(etmp(1)>etmp(2))then
        self%lce(j,i)=-ABS(self%lce(j,i))
      else
        self%lce(j,i)=ABS(self%lce(j,i))
      end if
    end do
    !---
    kf=INT(kfg,4)
    CALL find_orient_listn(self%lco(i),kf,self%cell_np)
  end do
  deallocate(kf,kfg,forder)
ELSE
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Orienting surface geometry'
  !$omp parallel do private(etmp,ed)
  do i=1,self%ne
    ed=self%le(:,i)
    etmp=self%global%lp(ed)
    if(self%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
      IF(ed(1)>ed(2))THEN
        etmp(1)=etmp(1)+self%global%np
      ELSE
        etmp(2)=etmp(2)+self%global%np
      END IF
    end if
    if(etmp(1)>etmp(2))then
      self%global%le(i)=-ABS(self%global%le(i))
    end if
  end do
  allocate(self%lco(self%nc))
  ! !$omp parallel do private(kf,j,etmp,ftmp)
  allocate(kf(self%cell_np),kfg(self%cell_np))
  do i=1,self%nc
    kfg=self%global%lp(self%lc(:,i))
    do j=1,self%cell_ne
      ed=self%cell_ed(:,j)
      etmp=kfg(ed)
      if(self%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
        IF(self%lc(ed(1),i)>self%lc(ed(2),i))THEN
          etmp(1)=etmp(1)+self%global%np
        ELSE
          etmp(2)=etmp(2)+self%global%np
        END IF
      end if
      if(etmp(1)>etmp(2))then
        self%lce(j,i)=-ABS(self%lce(j,i))
      else
        self%lce(j,i)=ABS(self%lce(j,i))
      end if
    end do
    kf=INT(kfg,4)
    CALL find_orient_listn(self%lco(i),kf,self%cell_np)
  end do
  deallocate(kf,kfg)
END IF
DEBUG_STACK_POP
end subroutine bmesh_global_orient
!---------------------------------------------------------------------------------
!> Display mesh statistics and counts.
!---------------------------------------------------------------------------------
subroutine bmesh_global_stats(self)
class(oft_bmesh), intent(inout) :: self !< Mesh object
integer(i8) :: a,i,tmp(3)
integer(i4) :: nbtmp,ierr
real(r8) :: tmpr(2),area
DEBUG_STACK_PUSH
call mesh_global_resolution(self)
area = self%area()
if(oft_env%head_proc)then
  write(*,'(2A)')oft_indent,'Mesh statistics:'
  CALL oft_increase_indent
  write(*,'(2A,ES11.3)')oft_indent,'Area         =',area
  write(*,'(2A,I8)')oft_indent,'# of points  =',self%global%np
  write(*,'(2A,I8)')oft_indent,'# of edges   =',self%global%ne
  write(*,'(2A,I8)')oft_indent,'# of cells   =',self%global%nc
endif
if(self%fullmesh)then
  IF(self%periodic%nper>0)THEN
    self%global%nbp=COUNT(self%global%gbp(self%lbp).AND.self%pstitch%leo)
    self%global%nbe=COUNT(self%global%gbe(self%lbe).AND.self%estitch%leo)
    self%global%nbc=COUNT(self%global%gbc)
  ELSE
    self%global%nbp=COUNT(self%global%gbp)
    self%global%nbe=COUNT(self%global%gbe)
    self%global%nbc=COUNT(self%global%gbc)
  END IF
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',self%global%nbp
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',self%global%nbe
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',self%global%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',self%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',self%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',self%hmax
    CALL oft_decrease_indent
  endif
else
#ifdef HAVE_MPI
  !---Count boundary points
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbp
    if(self%global%gbp(self%lbp(i)).AND.self%pstitch%leo(i))a=a+1
  enddo
  tmp(1)=a
  !---Count boundary edges
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbe
    if(self%global%gbe(self%lbe(i)).AND.self%estitch%leo(i))a=a+1
  enddo
  tmp(2)=a
  !---Count boundary cells
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,self%nbc
    if(self%global%gbc(self%lbc(i)))a=a+1
  enddo
  tmp(3)=a
  !---Check volumes
  tmpr(1)=SUM(self%ca)
  tmpr(2)=SUM(self%va)
  !---Get global sums
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,3,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmpr,2,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  self%global%nbp=tmp(1)
  self%global%nbe=tmp(2)
  self%global%nbc=tmp(3)
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',self%global%nbp
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',self%global%nbe
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',self%global%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',self%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',self%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',self%hmax
    CALL oft_decrease_indent
  endif
#else
  CALL oft_abort("Distributed mesh requires MPI","bmesh_global_stats",__FILE__)
#endif
endif
DEBUG_STACK_POP
end subroutine bmesh_global_stats
end module oft_mesh_global_util
