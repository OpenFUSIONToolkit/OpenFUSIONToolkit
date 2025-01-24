!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
MODULE oft_mesh_global_util
USE oft_base
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh, oft_bmesh, rgrnd
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_stitching, ONLY: oft_global_stitch !mesh_global_pstitch, mesh_global_estitch
IMPLICIT NONE
#include "local.h"
contains
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_boundary
!------------------------------------------------------------------------------
!> Set global boundary flags for distributed meshes.
!------------------------------------------------------------------------------
subroutine mesh_global_boundary(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,k
real(r8), allocatable, dimension(:) :: btrans
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Constructing global boundary'
!---------------------------------------------------------------------------
! Set global boundary points and cells
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
! Set global boundary edges
!---------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_periodic
!------------------------------------------------------------------------------
!> Set global indexing from mesh periodicity
!------------------------------------------------------------------------------
subroutine mesh_global_periodic(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,iper
INTEGER(8) :: npp,nep,nfp
INTEGER(8), ALLOCATABLE :: gtmp(:)
IF(mesh%periodic%nper==0)RETURN
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2X,A)')'Setting Up Periodic Mesh'
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
  WRITE(*,'(4X,A,I8)')'  # of Periodic directions =',mesh%periodic%nper
  WRITE(*,'(4X,A,I8)')'  # of Periodic Points     =',npp
  WRITE(*,'(4X,A,I8)')'  # of Periodic Edges      =',nep
  WRITE(*,'(4X,A,I8)')'  # of Periodic Faces      =',nfp
end if
DEBUG_STACK_POP
end subroutine mesh_global_periodic
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_orient
!------------------------------------------------------------------------------
!> Set global orientations for edges and faces
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Set cell curvature flag for common cases
!!
!! @param[in,out] self Mesh object to set curvature on
!! @param[in,out] flag How to set curvature (1 -> all cells, 2 -> contains boundary edge)
!------------------------------------------------------------------------------
subroutine mesh_global_set_curved(self,flag)
CLASS(oft_amesh), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: flag
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
! !------------------------------------------------------------------------------
! ! SUBROUTINE: mesh_global_orient_face_edge
! !------------------------------------------------------------------------------
! !> Apply global orientation to face edges
! !!
! !! @param[in,out] fc Face orientation index
! !! @param[in,out] ftmp Face edges in local rep, transposed and oriented on output to global orientation
! !------------------------------------------------------------------------------
! subroutine mesh_global_orient_face_edge(fc,ftmp)
! integer(i4), intent(in) :: fc
! integer(i4), intent(inout) :: ftmp(3)
! DEBUG_STACK_PUSH
! CALL orient_listn(fc,ftmp,3_i4)
! IF(fc<0)ftmp=-ftmp
! DEBUG_STACK_POP
! end subroutine mesh_global_orient_face_edge
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_save
!------------------------------------------------------------------------------
!> Get I/O scope for current mesh
!! - Element offsets for each geometric type
!! - Element span for each geometric type
!------------------------------------------------------------------------------
subroutine mesh_global_save(mesh)
class(oft_amesh), intent(inout) :: mesh
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
do i=1,mesh%nbp ! Loop over boundary points
  if(mesh%pstitch%leo(i))nloc(4) = nloc(4) + 1 ! Processor owns this point
enddo
!$omp parallel do reduction(+:nloc)
do i=1,mesh%np
  if(mesh%bp(i))cycle
  nloc(4) = nloc(4) + 1 ! Interior points (all owned)
enddo
!---Local count of owned boundary edges
!$omp parallel do reduction(+:nloc)
do i=1,mesh%nbe ! Loop over boundary edges
  if(mesh%estitch%leo(i))THEN
    nloc(1) = nloc(1) + 1 ! Processor owns this edge
    IF(mesh%global%gbe(mesh%lbe(i)))nloc(6)=nloc(6)+1
  END IF
enddo
!$omp parallel do reduction(+:nloc)
do i=1,mesh%ne
  if(mesh%be(i))cycle
  nloc(1) = nloc(1) + 1 ! Interior edges (all owned)
enddo
nloc(3)=INT(mesh%nc,8) ! All cells are owned
!---Handle faces for volume meshes
SELECT TYPE(mesh)
CLASS IS(oft_mesh)
  !---Local count of owned boundary faces
  !$omp parallel do reduction(+:nloc)
  do i=1,mesh%nbf ! Loop over boundary faces
    if(mesh%fstitch%leo(i))nloc(2)=nloc(2)+1 ! Processor owns this face
  enddo
  !$omp parallel do reduction(+:nloc)
  do i=1,mesh%nf
    if(mesh%bf(i))cycle
    nloc(2)=nloc(2)+1 ! Interior face (all owned)
  enddo
  !---Local count of global boundary faces
  !$omp parallel do reduction(+:nloc)
  do i=1,mesh%nbf ! Loop over boundary faces
    if(mesh%global%gbf(mesh%lbf(i)))nloc(5)=nloc(5)+1 ! Face is global boundary
  enddo
  mesh%save%nbf=nloc(5)
END SELECT
!---Get maximum element counts over every processor
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(nloc(1:4),ntmp,4,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
#else
ntmp=nloc(1:4)
#endif
mesh%save%npmax=ntmp(4)
mesh%save%nemax=ntmp(1)
! mesh%save%nfmax=ntmp(2)
nfmax=ntmp(2)
mesh%save%ncmax=ntmp(3)
!---Initialize ouput index counters
mesh%save%npb=0
mesh%save%neb=0
! mesh%save%nfb=0
nfb=0
mesh%save%ncb=0
mesh%save%npl=nloc(4)
mesh%save%nel=nloc(1)
! mesh%save%nfl=nloc(2)
nfl=nloc(2)
mesh%global%ne=nloc(1)
! mesh%global%nf=nloc(2)
nfg=nloc(2)
! mesh%bmesh%global%ne=nloc(6)
!---Loop over processors to determine output index counts
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  do while(l.NE.oft_env%rank)
    call MPI_SEND(nloc,6,OFT_MPI_I8,j,1,oft_env%COMM,ierr)
    call MPI_RECV(nloc,6,OFT_MPI_I8,k,1,oft_env%COMM,stat,ierr)

    if(l<oft_env%rank)then
      mesh%save%neb = mesh%save%neb + nloc(1) ! Update first edge index
      nfb = nfb + nloc(2) ! Update first face index
      mesh%save%ncb = mesh%save%ncb + nloc(3) ! Update first cell index
      mesh%save%npb = mesh%save%npb + nloc(4) ! Update first point index
    endif
    mesh%global%ne = mesh%global%ne + nloc(1) ! Update global edge count
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
mesh%save%npl = mesh%save%npl + mesh%save%npb ! Set last point index
mesh%save%nel = mesh%save%nel + mesh%save%neb ! Set last edge index
!---Handle faces for volume meshes
SELECT TYPE(mesh)
CLASS IS(oft_mesh)
  mesh%save%nfmax = nfmax
  mesh%save%nfb = nfb       
  mesh%save%nfl = nfl + nfb ! Set last face index
  mesh%global%nf = nfg  
END SELECT
DEBUG_STACK_POP
end subroutine mesh_global_save
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_igrnd
!------------------------------------------------------------------------------
!> Determine location of grounding point on mesh.
!------------------------------------------------------------------------------
SUBROUTINE mesh_global_igrnd(mesh)
class(oft_amesh), intent(inout) :: mesh
REAL(r8) :: c,r2min,r2
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
!---Locate local surface ground point
mesh%igrnd=-1
r2min=1.d99
DO i=1,mesh%nbp
  IF(.NOT.mesh%global%gbp(mesh%lbp(i)))CYCLE
  r2=sum((mesh%r(:,mesh%lbp(i))-rgrnd)**2)
  IF(r2<r2min)THEN
    r2min=r2
    mesh%igrnd(1)=mesh%lbp(i)
  END IF
END DO
IF(mesh%fullmesh)THEN
  !---Check periodicity
  IF(mesh%igrnd(1)==-1)THEN
    IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'No grounding point found: Mesh is fully periodic!'
  ELSE
    grnd_global=mesh%global%lp(mesh%igrnd(1))
    mesh%igrnd=-1
    DO i=1,mesh%nbp
      IF(mesh%global%lp(mesh%lbp(i))==grnd_global)THEN
        IF(ALL(mesh%igrnd/=-1))CALL oft_abort('Ground point repeated on same domain', &
          'mesh_global_igrnd',__FILE__)
        IF(mesh%igrnd(1)==-1)THEN
          mesh%igrnd(1)=mesh%lbp(i)
        ELSE
          mesh%igrnd(2)=mesh%lbp(i)
        END IF
      END IF
    END DO
    IF(oft_env%head_proc)WRITE(*,'(2A,I8)')oft_indent,'Surface grounded at vertex',mesh%igrnd(1)
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
  IF(c<r2min)mesh%igrnd=-1
  !---Increment ring comm
  l=l-1
  IF(l<0)l=oft_env%nprocs-1
END DO
#endif
!---Handle periodicity
grnd_global=-1
IF(mesh%igrnd(1)/=-1)grnd_global=mesh%global%lp(mesh%igrnd(1))
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,grnd_global,1,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
#endif
mesh%igrnd=-1
IF(grnd_global==-1)THEN
  IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'No grounding point found: Mesh is fully periodic!'
ELSE
  owned=.FALSE.
  DO i=1,mesh%nbp
    IF(mesh%global%lp(mesh%lbp(i))==grnd_global)THEN
      IF(ALL(mesh%igrnd/=-1))CALL oft_abort('Ground point repeated on same domain', &
        'mesh_global_igrnd',__FILE__)
      IF(mesh%igrnd(1)==-1)THEN
        mesh%igrnd(1)=mesh%lbp(i)
        owned=mesh%pstitch%leo(i)
      ELSE
        mesh%igrnd(2)=mesh%lbp(i)
        owned=owned.OR.mesh%pstitch%leo(i)
      END IF
    END IF
  END DO
  !---Write out grounding information
  IF(oft_env%head_proc)WRITE(*,'(2A,I8)')oft_indent,'Surface grounded at vertex ',grnd_global
END IF
DEBUG_STACK_POP
END SUBROUTINE mesh_global_igrnd
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_resolution
!------------------------------------------------------------------------------
!> Compute mesh resolution statistics.
!------------------------------------------------------------------------------
subroutine mesh_global_resolution(self)
class(oft_amesh), intent(inout) :: self
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
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_global_stats
!------------------------------------------------------------------------------
!> Display mesh statistics and counts.
!------------------------------------------------------------------------------
subroutine mesh_global_stats(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i8) :: a,i,tmp(4)
integer(i4) :: ierr
real(r8) :: tmpr(2),vol,area
DEBUG_STACK_PUSH
call mesh_global_resolution(mesh)
vol = mesh%volume()
area = mesh%bmesh%area()
if(mesh%fullmesh)then
  if(oft_env%head_proc)then
    write(*,'(2A)')oft_indent,'Mesh statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'Volume          =',vol
    write(*,'(2A,ES11.3)')oft_indent,'Surface area    =',area
    write(*,'(2A,I8)')    oft_indent,'# of points     =',mesh%global%np
    write(*,'(2A,I8)')    oft_indent,'# of edges      =',mesh%global%ne
    write(*,'(2A,I8)')    oft_indent,'# of faces      =',mesh%global%nf
    write(*,'(2A,I8)')    oft_indent,'# of cells      =',mesh%global%nc
    write(*,'(2A,I8)')    oft_indent,'# of boundary points =',mesh%nbp
    write(*,'(2A,I8)')    oft_indent,'# of boundary edges  =',mesh%nbe
    write(*,'(2A,I8)')    oft_indent,'# of boundary faces  =',mesh%nbf
    write(*,'(2A,I8)')    oft_indent,'# of boundary cells  =',mesh%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',mesh%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',mesh%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',mesh%hmax
    CALL oft_decrease_indent
  endif
else
#ifdef HAVE_MPI
  if(oft_env%head_proc)then
    write(*,'(2A)')oft_indent,'Mesh statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'Volume          =',vol
    write(*,'(2A,ES11.3)')oft_indent,'Surface area    =',area
    write(*,'(2A,I8)')    oft_indent,'# of points     =',mesh%global%np
    write(*,'(2A,I8)')    oft_indent,'# of edges      =',mesh%global%ne
    write(*,'(2A,I8)')    oft_indent,'# of faces      =',mesh%global%nf
    write(*,'(2A,I8)')    oft_indent,'# of cells      =',mesh%global%nc
  endif
  !---Count boundary points
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbp
    if(mesh%global%gbp(mesh%lbp(i)).AND.mesh%pstitch%leo(i))a=a+1
  enddo
  tmp(1)=a
  !---Count boundary edges
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbe
    if(mesh%global%gbe(mesh%lbe(i)).AND.mesh%estitch%leo(i))a=a+1
  enddo
  tmp(2)=a
  !---Count boundary faces
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbf
    if(mesh%global%gbf(mesh%lbf(i)).AND.mesh%fstitch%leo(i))a=a+1
  enddo
  tmp(3)=a
  !---Count boundary cells
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbc
    if(mesh%global%gbc(mesh%lbc(i)))a=a+1
  enddo
  tmp(4)=a
  !---Check volumes
  tmpr(1)=SUM(mesh%cv)
  tmpr(2)=SUM(mesh%vv)
  !---Get global sums
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,4,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmpr,2,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',tmp(1)
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',tmp(2)
    write(*,'(2A,I8)')oft_indent,'# of boundary faces  =',tmp(3)
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',tmp(4)
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',mesh%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',mesh%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',mesh%hmax
    CALL oft_decrease_indent
  endif
#else
CALL oft_abort("Distributed mesh requires MPI","mesh_global_stats",__FILE__)
#endif
endif
DEBUG_STACK_POP
end subroutine mesh_global_stats
!------------------------------------------------------------------------------
! SUBROUTINE: bmesh_global_boundary
!------------------------------------------------------------------------------
!> Set global boundary flags for distributed meshes.
!------------------------------------------------------------------------------
subroutine bmesh_global_boundary(mesh)
class(oft_bmesh), intent(inout) :: mesh
integer(i4) :: i,j,k
real(r8), allocatable, dimension(:) :: btrans
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Constructing global boundary'
!---------------------------------------------------------------------------
! Set global boundary points and cells
!---------------------------------------------------------------------------
IF(ASSOCIATED(mesh%global%gbc))DEALLOCATE(mesh%global%gbc)
allocate(mesh%global%gbc(mesh%nc)) ! Allocate boundary cell flag
mesh%global%gbc = .FALSE. ! Initiliaze boundary flag
allocate(btrans(mesh%np)) ! Allocate Send/Recv refernce arrays
btrans=0.d0
!---Mark global boundary points and cells on local processor
do i=1,mesh%nbe
  if(mesh%global%gbe(mesh%lbe(i)))then ! Face is on global boundary (Root Element)
    do j=1,2 ! Loop over edge end points
      k=mesh%le(j,mesh%lbe(i))
      mesh%global%gbp(k)=.TRUE. ! Point is a boundary point
      btrans(k)=1.d0
    enddo
    mesh%global%gbc(mesh%lec(mesh%kec(mesh%lbe(i))))=.TRUE. ! Shared cell is a boundary cell
  endif
enddo
!---Synchronize across seams
CALL oft_global_stitch(mesh%pstitch,btrans,1) ! Sync boundary points across processors
do i=1,mesh%nbp
  if(btrans(mesh%lbp(i))>0.d0)mesh%global%gbp(mesh%lbp(i))=.TRUE. ! Point is a global boundary point
enddo
deallocate(btrans)
DEBUG_STACK_POP
end subroutine bmesh_global_boundary
!------------------------------------------------------------------------------
! SUBROUTINE: bmesh_global_periodic
!------------------------------------------------------------------------------
!> Set global indexing from mesh periodicity
!------------------------------------------------------------------------------
subroutine bmesh_global_periodic(mesh)
class(oft_bmesh), intent(inout) :: mesh
integer(i4) :: i,j,iper
INTEGER(8) :: npp,nep,nfp
INTEGER(8), ALLOCATABLE :: gtmp(:)
IF(mesh%periodic%nper==0)RETURN
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Setting Up Periodic Mesh'
!---Edges can only have one periodic parent
DO i=1,mesh%nbe
  j=mesh%lbe(i)
  IF(mesh%periodic%le(j)>0)THEN
    mesh%global%le(j)=mesh%global%le(mesh%periodic%le(j))
  END IF
END DO
!---Propogate point parents through periodicity
DO iper=1,mesh%periodic%nper
  DO i=1,mesh%nbp
    j=mesh%lbp(i)
    IF(mesh%periodic%lp(j)>0)THEN
      mesh%global%lp(j)=mesh%global%lp(mesh%periodic%lp(j))
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
if(oft_debug_print(1))then
  CALL oft_increase_indent
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic directions =',mesh%periodic%nper
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Points     =',npp
  WRITE(*,'(2A,I8)')oft_indent,'# of Periodic Edges      =',nep
  CALL oft_decrease_indent
end if
DEBUG_STACK_POP
end subroutine bmesh_global_periodic
!------------------------------------------------------------------------------
! SUBROUTINE: bmesh_global_orient
!------------------------------------------------------------------------------
!> Set global orientations for edges and faces
!------------------------------------------------------------------------------
subroutine bmesh_global_orient(smesh,parent)
CLASS(oft_bmesh), intent(inout) :: smesh
CLASS(oft_mesh), optional, intent(in) :: parent
integer(i4) :: i,j,ind,ed(2)
integer(i4), ALLOCATABLE, DIMENSION(:) :: kf,forder,ltmp
integer(8) :: etmp(2)
integer(8), ALLOCATABLE, DIMENSION(:) :: kfg
DEBUG_STACK_PUSH
IF(PRESENT(parent))THEN
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Copying orientation to surface'
  !$omp parallel do private(etmp)
  do i=1,smesh%ne
    etmp=smesh%le(:,i)
    etmp=smesh%parent%lp(etmp)
    etmp=parent%global%lp(etmp)
    if(etmp(1)>etmp(2))then
      smesh%global%le(i)=-ABS(smesh%global%le(i))
    end if
  end do
  allocate(smesh%lco(smesh%nc))
  allocate(kf(smesh%cell_np),kfg(smesh%cell_np),forder(smesh%cell_np))
  ! !$omp parallel do private(kf,j,etmp,forder)
  do i=1,smesh%nc
    kfg=smesh%lc(:,i)
    kfg=smesh%parent%lp(kfg)
    kfg=parent%global%lp(kfg)
    do j=1,smesh%cell_np
      etmp=(/kfg(smesh%cell_ed(1,j)),kfg(smesh%cell_ed(2,j))/)
      if(etmp(1)>etmp(2))then
        smesh%lce(j,i)=-ABS(smesh%lce(j,i))
      else
        smesh%lce(j,i)=ABS(smesh%lce(j,i))
      end if
    end do
    !---
    kf=INT(kfg,4)
    CALL find_orient_listn(smesh%lco(i),kf,smesh%cell_np)
  end do
  deallocate(kf,kfg,forder)
ELSE
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Orienting surface geometry'
  !$omp parallel do private(etmp,ed)
  do i=1,smesh%ne
    ed=smesh%le(:,i)
    etmp=smesh%global%lp(ed)
    if(smesh%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
      IF(ed(1)>ed(2))THEN
        etmp(1)=etmp(1)+smesh%global%np
      ELSE
        etmp(2)=etmp(2)+smesh%global%np
      END IF
    end if
    if(etmp(1)>etmp(2))then
      smesh%global%le(i)=-ABS(smesh%global%le(i))
    end if
  end do
  allocate(smesh%lco(smesh%nc))
  ! !$omp parallel do private(kf,j,etmp,ftmp)
  allocate(kf(smesh%cell_np),kfg(smesh%cell_np))
  do i=1,smesh%nc
    kfg=smesh%global%lp(smesh%lc(:,i))
    do j=1,smesh%cell_ne
      ed=smesh%cell_ed(:,j)
      etmp=kfg(ed)
      if(smesh%periodic%nper>0.AND.(etmp(1)==etmp(2)))then
        IF(smesh%lc(ed(1),i)>smesh%lc(ed(2),i))THEN
          etmp(1)=etmp(1)+smesh%global%np
        ELSE
          etmp(2)=etmp(2)+smesh%global%np
        END IF
      end if
      if(etmp(1)>etmp(2))then
        smesh%lce(j,i)=-ABS(smesh%lce(j,i))
      else
        smesh%lce(j,i)=ABS(smesh%lce(j,i))
      end if
    end do
    kf=INT(kfg,4)
    CALL find_orient_listn(smesh%lco(i),kf,smesh%cell_np)
  end do
  deallocate(kf,kfg)
END IF
DEBUG_STACK_POP
end subroutine bmesh_global_orient
!------------------------------------------------------------------------------
! SUBROUTINE: bmesh_global_stats
!------------------------------------------------------------------------------
!> Display mesh statistics and counts.
!------------------------------------------------------------------------------
subroutine bmesh_global_stats(mesh)
class(oft_bmesh), intent(inout) :: mesh
integer(i8) :: a,i,tmp(3)
integer(i4) :: ierr
real(r8) :: tmpr(2),area
DEBUG_STACK_PUSH
call mesh_global_resolution(mesh)
area = mesh%area()
if(oft_env%head_proc)then
  write(*,'(2A)')oft_indent,'Mesh statistics:'
  CALL oft_increase_indent
  write(*,'(2A,ES11.3)')oft_indent,'Area         =',area
  write(*,'(2A,I8)')oft_indent,'# of points  =',mesh%global%np
  write(*,'(2A,I8)')oft_indent,'# of edges   =',mesh%global%ne
  write(*,'(2A,I8)')oft_indent,'# of cells   =',mesh%global%nc
endif
if(mesh%fullmesh)then
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',mesh%nbp
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',mesh%nbe
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',mesh%nbc
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',mesh%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',mesh%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',mesh%hmax
    CALL oft_decrease_indent
  endif
else
#ifdef HAVE_MPI
  !---Count boundary points
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbp
    if(mesh%global%gbp(mesh%lbp(i)).AND.mesh%pstitch%leo(i))a=a+1
  enddo
  tmp(1)=a
  !---Count boundary edges
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbe
    if(mesh%global%gbe(mesh%lbe(i)).AND.mesh%estitch%leo(i))a=a+1
  enddo
  tmp(2)=a
  !---Count boundary cells
  a=0
  !$omp parallel do reduction(+:a)
  do i=1,mesh%nbc
    if(mesh%global%gbc(mesh%lbc(i)))a=a+1
  enddo
  tmp(3)=a
  !---Check volumes
  tmpr(1)=SUM(mesh%ca)
  tmpr(2)=SUM(mesh%va)
  !---Get global sums
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,3,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,tmpr,2,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  if(oft_env%head_proc)then
    write(*,'(2A,I8)')oft_indent,'# of boundary points =',tmp(1)
    write(*,'(2A,I8)')oft_indent,'# of boundary edges  =',tmp(2)
    write(*,'(2A,I8)')oft_indent,'# of boundary cells  =',tmp(3)
    CALL oft_decrease_indent
    write(*,'(2A)')oft_indent,'Resolution statistics:'
    CALL oft_increase_indent
    write(*,'(2A,ES11.3)')oft_indent,'hmin =',mesh%hmin
    write(*,'(2A,ES11.3)')oft_indent,'hrms =',mesh%hrms
    write(*,'(2A,ES11.3)')oft_indent,'hmax =',mesh%hmax
    CALL oft_decrease_indent
  endif
#else
  CALL oft_abort("Distributed mesh requires MPI","bmesh_global_stats",__FILE__)
#endif
endif
DEBUG_STACK_POP
end subroutine bmesh_global_stats
end module oft_mesh_global_util
