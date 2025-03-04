!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_local.F90
!
!> Local grid construction.
!! - Geometric link lists
!! - Local boundary information
!!
!! @authors George Marklin and Chris Hansen
!! @date April 2008 - Present
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_local
USE oft_base
USE oft_sort, ONLY: search_array, sort_array
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh, oft_bmesh
USE oft_mesh_local_util, ONLY: mesh_local_findedge, mesh_local_findface
IMPLICIT NONE
#include "local.h"
INTERFACE
!---------------------------------------------------------------------------
!> Parition a tetrahedral mesh using METIS
!---------------------------------------------------------------------------
  SUBROUTINE oft_metis_partmesh(nc,np,ncp,lc,npart,cpart,info) BIND(C,NAME="oft_metis_partMesh")
  IMPORT c_int
  INTEGER(c_int), INTENT(in) :: nc !< Number of cells
  INTEGER(c_int), INTENT(in) :: np !< Number of vertices
  INTEGER(c_int), INTENT(in) :: ncp !< Number of points per cell
  INTEGER(c_int), DIMENSION(ncp,nc), INTENT(inout) :: lc !< Cell list [ncp,nc]
  INTEGER(c_int), INTENT(in) :: npart !< Number of partitions
  INTEGER(c_int), DIMENSION(npart), INTENT(inout) :: cpart !< Partition for each cell [nc]
  INTEGER(c_int), INTENT(inout) :: info !< Parition return status
  END SUBROUTINE oft_metis_partmesh
END INTERFACE
contains
!------------------------------------------------------------------------------
!> Driver for grid construction
!------------------------------------------------------------------------------
subroutine mesh_local_init(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
DEBUG_STACK_PUSH
ALLOCATE(self%ho_info%is_curved(self%nc))
self%ho_info%is_curved=.FALSE.
call mesh_volumes(self)
call amesh_edges(self)
call amesh_to_cell(self)
call mesh_faces(self)
call mesh_neighbors(self)
call mesh_boundary(self)
call amesh_interactions(self)
IF(.NOT.ASSOCIATED(self%tloc_c))THEN
  self%nparts=oft_env%nthreads*oft_env%nparts
  ALLOCATE(self%tloc_p(self%nparts))
  ALLOCATE(self%tloc_e(self%nparts))
  ALLOCATE(self%tloc_f(self%nparts))
  ALLOCATE(self%tloc_c(self%nparts))
  call mesh_local_partition(self,self%tloc_p,self%tloc_e,self%tloc_f,self%tloc_c,self%nparts)
END IF
DEBUG_STACK_POP
end subroutine mesh_local_init
!------------------------------------------------------------------------------
!> Compute cell and vertex volumes, rectifying negative volume cells.
!!
!! Zero volume cells or vertices are also caught for mesh validation.
!------------------------------------------------------------------------------
subroutine mesh_volumes(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i4) :: i,j
real(r8) :: vol,f(4),gop(3,4)
DEBUG_STACK_PUSH
if(associated(self%cv))then
  if(size(self%cv).NE.self%nc)allocate(self%cv(self%nc),self%vv(self%np))
else
  allocate(self%cv(self%nc),self%vv(self%np))
endif
self%vv=0.d0 ! initialize vertex volumes
f=1.d0/4.d0 ! NEED TO UPDATE
do i=1,self%nc ! loop over cells
  CALL self%jacobian(i,f,gop,self%cv(i))
  IF(self%cv(i)<0.d0)CALL self%invert_cell(i)
  self%cv(i)=ABS(self%cv(i))
enddo
do i=1,self%nc ! accumulate vertex volumes
  vol=self%cv(i)/REAL(self%cell_np,8) ! corner volume
  do j=1,self%cell_np
    self%vv(self%lc(j,i))=self%vv(self%lc(j,i))+vol
  enddo
enddo
!---Check for zero volume cells
if(any(self%cv==0))call oft_abort('Zero volume cell detected','mesh_volumes',__FILE__)
!---Check for zero volume vertices
if(any(self%vv==0))call oft_abort('Floating vertex detected','mesh_volumes',__FILE__)
CALL bmesh_areas(self%bmesh)
DEBUG_STACK_POP
end subroutine mesh_volumes
!------------------------------------------------------------------------------
!> Identify global edges and link to points and cells.
!! - le List construction
!! - lce List construction
!! - klpe, llpe Linkage construction
!------------------------------------------------------------------------------
subroutine amesh_edges(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), allocatable :: nr(:),ir(:),jr(:)
integer(i4) :: i,j,js,je,jp,jn,k(8)
DEBUG_STACK_PUSH
!---Allocate local counters
allocate(nr(self%np+1),ir(self%np+1))
allocate(jr(self%cell_ne*self%nc))
!---Loop over cells & count raw edges
nr=0 ! initialize raw edge counter
do i=1,self%nc
  k(1:self%cell_np)=self%lc(:,i) ! cell corners
  do j=1,self%cell_ne ! loop over edges
    js=min(k(self%cell_ed(1,j)),k(self%cell_ed(2,j))) ! low points
    nr(js)=nr(js)+1 ! edge counter
  enddo
enddo
!---Build low point raw pointer list
ir(self%np+1)=self%cell_ne*self%nc+1
do i=self%np,1,-1 ! cumulative raw edge count
  ir(i)=ir(i+1)-nr(i)
enddo
if(ir(1)/=1)call oft_abort('Bad raw edge count','mesh_edges',__FILE__)
nr=0 ! reset raw edge counter
!---Loop over cells & index raw edges
do i=1,self%nc
  k(1:self%cell_np)=self%lc(:,i) ! cell corners
  do j=1,self%cell_ne  ! loop over edges
    js=min(k(self%cell_ed(1,j)),k(self%cell_ed(2,j))) ! low points
    je=max(k(self%cell_ed(1,j)),k(self%cell_ed(2,j))) ! high points
    jr(ir(js)+nr(js))=je      ! high point list
    nr(js)=nr(js)+1           ! edge counter
  enddo
enddo
!---Sort raw high point list
!$omp parallel do
do i=1,self%np
  if(nr(i)>1)call sort_array(jr(ir(i):ir(i)+nr(i)-1),nr(i))
enddo
deallocate(nr)
allocate(self%klpe(self%np+1))
self%klpe=0
!---Loop over raw edges & count unique edges
!$omp parallel do private(j,js,je)
do i=1,self%np
  js=ir(i)
  je=ir(i+1)-1
  IF(je<js)CYCLE
  self%klpe(i)=1
  do j=js+1,je
    if(jr(j)>jr(j-1))self%klpe(i)=self%klpe(i)+1
  enddo
enddo
self%ne=sum(self%klpe) ! total number of unique edges
self%klpe(self%np+1)=self%ne+1
do i=self%np,1,-1 ! cumulative unique edge count
  self%klpe(i)=self%klpe(i+1)-self%klpe(i)
enddo
if(self%klpe(1)/=1)call oft_abort('Bad unique edge count','mesh_edges',__FILE__)
allocate(self%llpe(self%ne))
!---Loop over raw edges & index unique edges
!$omp parallel do private(j,js,je,jn)
do i=1,self%np
  js=ir(i)
  je=ir(i+1)-1
  IF(je<js)CYCLE
  jn=self%klpe(i)
  self%llpe(jn)=jr(js)
  do j=js+1,je
    if(jr(j)>jr(j-1))then
      jn=jn+1
      self%llpe(jn)=jr(j)
    endif
  enddo
enddo
deallocate(ir,jr)
allocate(self%le(2,self%ne))
!$omp parallel do
do i=1,self%np ! construct global edge list
  self%le(1,self%klpe(i):self%klpe(i+1)-1)=i
enddo
!$omp parallel do
do i=1,self%ne
  self%le(2,i)=self%llpe(i)
enddo
allocate(self%lce(self%cell_ne,self%nc))
!---Loop over cells & index cells to edges
!$omp parallel do private(k,j,js,je,jp,jn)
do i=1,self%nc
  k(1:self%cell_np)=self%lc(:,i) ! cell corners
  do j=1,self%cell_ne  ! loop over edges
    js=min(k(self%cell_ed(1,j)),k(self%cell_ed(2,j))) ! low points
    je=max(k(self%cell_ed(1,j)),k(self%cell_ed(2,j))) ! high points
    jp=self%klpe(js)      ! pointer into low point list
    jn=self%klpe(js+1)-jp ! number of shared edges
    self%lce(j,i)=search_array(je,self%llpe(jp:jp+jn-1),jn)+jp-1
    self%lce(j,i)=sign(self%lce(j,i),k(self%cell_ed(2,j))-k(self%cell_ed(1,j))) ! apply orientation
  enddo
enddo
DEBUG_STACK_POP
end subroutine amesh_edges
!------------------------------------------------------------------------------
!> Locate point, edge, face and cell neighbor cells.
!! - lfc List construction
!! - lcc List construction
!! - kpc, lpc Linkage construction
!! - kec, lec Linkage construction
!------------------------------------------------------------------------------
subroutine amesh_to_cell(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), allocatable :: nr(:)
integer(i4) :: i,j,k
DEBUG_STACK_PUSH
allocate(nr(self%np+1),self%kpc(self%np+1))
nr=0 ! initialize cell counter
do i=1,self%nc ! loop over cells
  do j=1,self%cell_np  ! loop over corners
    k=self%lc(j,i) ! corner number
    nr(k)=nr(k)+1 ! count cell to corner
  enddo
enddo
self%npc=sum(nr)
self%kpc(self%np+1)=self%npc+1
do i=self%np,1,-1 ! cumulative point to cell count
  self%kpc(i)=self%kpc(i+1)-nr(i)
enddo
if(self%kpc(1)/=1)call oft_abort('Bad point to cell count','mesh_to_cell',__FILE__)
allocate(self%lpc(self%npc))
nr=0 ! reset cell counter
do i=1,self%nc ! loop over cells
  do j=1,self%cell_np  ! loop over corners
    k=self%lc(j,i) ! corner number
    self%lpc(self%kpc(k)+nr(k))=i ! index cell to corner
    nr(k)=nr(k)+1       ! count cell to corner
  enddo
enddo
deallocate(nr)
allocate(nr(self%ne+1),self%kec(self%ne+1))
nr=0 ! initialize cell counter
do i=1,self%nc ! loop over cells
  do j=1,self%cell_ne  ! loop over edges
    k=abs(self%lce(j,i)) ! edge number
    nr(k)=nr(k)+1 ! count cell to edge
  enddo
enddo
self%nec=sum(nr)
self%kec(self%ne+1)=self%nec+1
do i=self%ne,1,-1 ! cumulative edge to cell count
  self%kec(i)=self%kec(i+1)-nr(i)
enddo
if(self%kec(1)/=1)call oft_abort('Bad edge to cell count','mesh_to_cell',__FILE__)
allocate(self%lec(self%nec))
nr=0 ! reset cell counter
do i=1,self%nc ! loop over cells
  do j=1,self%cell_ne ! loop over edges
    k=abs(self%lce(j,i)) ! edge number
    self%lec(self%kec(k)+nr(k))=i ! index cell to edge
    nr(k)=nr(k)+1       ! count cell to edge
  enddo
enddo
deallocate(nr)
DEBUG_STACK_POP
end subroutine amesh_to_cell
!------------------------------------------------------------------------------
!> Index point to point, point to edge and edge to edge interactions.
!! - kpp, lpp Linkage construction
!! - kpe, lpe Linkage construction
!! - kee, lee Linkage construction
!------------------------------------------------------------------------------
subroutine amesh_interactions(self)
class(oft_amesh), intent(inout) :: self !< Mesh object
integer(i4), allocatable :: lcx(:,:),nr(:)
integer(i4) :: i,j,js,je,jn,jp,k
DEBUG_STACK_PUSH
! first do point to point interactions
allocate(self%kpp(self%np+1),nr(self%np+1))
allocate(lcx(self%cell_np,self%npc))
nr=0 ! initialize interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%np ! loop over points
  js=self%kpc(i)
  je=self%kpc(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    lcx(:,j)=self%lc(:,self%lpc(j)) ! get corners
  enddo
  jn=(je-js+1)
  call sort_array(lcx(:,js:je),self%cell_np,jn)
  !
  nr(i)=1
  jn=lcx(1,js)
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_np ! loop over corners
      if(lcx(k,j)/=jn)then
        nr(i)=nr(i)+1
        jn=lcx(k,j)
      endif
    enddo
  enddo
enddo
self%npp=sum(nr)
self%kpp(self%np+1)=self%npp+1
do i=self%np,1,-1 ! cumulative point to point count
  self%kpp(i)=self%kpp(i+1)-nr(i)
enddo
if(self%kpp(1)/=1)call oft_abort('Bad point to point count','mesh_interactions',__FILE__)
allocate(self%lpp(self%npp))
nr=0 ! reset interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%np ! loop over points
  js=self%kpc(i)
  je=self%kpc(i+1)-1
  IF(je<js)CYCLE
  jn=lcx(1,js)
  self%lpp(self%kpp(i))=jn
  nr(i)=1
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_np ! loop over corners
      if(lcx(k,j)/=jn)then
        jn=lcx(k,j)
        self%lpp(self%kpp(i)+nr(i))=jn
        nr(i)=nr(i)+1
      endif
    enddo
  enddo
enddo
deallocate(nr,lcx)
! next do point to edge interactions
allocate(self%kpe(self%np+1),nr(self%np+1),lcx(self%cell_ne,self%npc))
nr=0 ! initialize interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%np ! loop over points
  js=self%kpc(i)
  je=self%kpc(i+1)-1
  IF(je<js)CYCLE
  do j=js,je ! loop over neighbor cells
    lcx(:,j)=ABS(self%lce(:,self%lpc(j))) ! get corners
  enddo
  jn=(je-js+1)
  call sort_array(lcx(:,js:je),self%cell_ne,jn)
  !
  nr(i)=0
  jn=-1
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_ne ! loop over edges
      if(ALL(self%le(:,lcx(k,j))/=i))cycle
      if(lcx(k,j)/=jn)then
        nr(i)=nr(i)+1
        jn=lcx(k,j)
      endif
    enddo
  enddo
enddo
self%npe=sum(nr)
self%kpe(self%np+1)=self%npe+1
do i=self%np,1,-1 ! cumulative edge to edge count
  self%kpe(i)=self%kpe(i+1)-nr(i)
enddo
if(self%kpe(1)/=1)call oft_abort('Bad point to edge count','mesh_interactions',__FILE__)
allocate(self%lpe(self%npe))
nr=0 ! reset interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%np ! loop over points
  js=self%kpc(i)
  je=self%kpc(i+1)-1
  IF(je<js)CYCLE
  nr(i)=0
  jn=-1 !lcx(1,js)
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_ne ! loop over edges
      if(ALL(self%le(:,lcx(k,j))/=i))cycle
      if(lcx(k,j)/=jn)then
        self%lpe(self%kpe(i)+nr(i))=lcx(k,j)
        nr(i)=nr(i)+1
        jn=lcx(k,j)
      endif
    enddo
  enddo
enddo
deallocate(nr,lcx)
! next do edge to edge interactions
allocate(self%kee(self%ne+1),nr(self%ne+1),lcx(self%cell_ne,self%nec))
nr=0 ! initialize interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%ne ! loop over edges
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    lcx(:,j)=abs(self%lce(:,self%lec(j))) ! get edges
  enddo
  jn=(je-js+1)
  call sort_array(lcx(:,js:je),self%cell_ne,jn)
  !
  nr(i)=1
  jn=lcx(1,js)
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_ne ! loop over edges
      if(lcx(k,j)/=jn)then
        nr(i)=nr(i)+1
        jn=lcx(k,j)
      endif
    enddo
  enddo
enddo
self%nee=sum(nr)
self%kee(self%ne+1)=self%nee+1
do i=self%ne,1,-1 ! cumulative edge to edge count
  self%kee(i)=self%kee(i+1)-nr(i)
enddo
if(self%kee(1)/=1)call oft_abort('Bad edge to edge count','mesh_interactions',__FILE__)
allocate(self%lee(self%nee))
nr=0 ! reset interaction counter
!$omp parallel do private(j,js,je,jn,k)
do i=1,self%ne ! loop over edges
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  self%lee(self%kee(i))=lcx(1,js)
  nr(i)=1
  jn=lcx(1,js)
  do j=js,je ! loop over neighbor cells
    do k=1,self%cell_ne ! loop over edges
      if(lcx(k,j)/=jn)then
        self%lee(self%kee(i)+nr(i))=lcx(k,j)
        nr(i)=nr(i)+1
        jn=lcx(k,j)
      endif
    enddo
  enddo
enddo
deallocate(nr,lcx)
DEBUG_STACK_POP
end subroutine amesh_interactions
!------------------------------------------------------------------------------
!> Identify global faces and link to points, edges and cells.
!! - lf List construction
!! - lcf List construction
!! - lfe List construction
!! - klef, llef Linkage construction
!------------------------------------------------------------------------------
subroutine mesh_faces(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i4), allocatable :: nr(:),ir(:),jr(:)
integer(i4) :: i,j,js,je,jp,jn,oflag,k(12)
DEBUG_STACK_PUSH
!---Allocate local counters
allocate(nr(self%ne+1),ir(self%ne+1))
allocate(jr(self%cell_nf*self%nc))
!---Loop over cells & count raw faces
nr=0 ! initialize raw face counter
do i=1,self%nc
  k(1:self%cell_ne)=abs(self%lce(:,i)) ! cell edges
  do j=1,self%cell_nf ! loop over faces
    js=minval(k(self%cell_fe(:,j))) ! low edge
    nr(js)=nr(js)+1 ! face counter
  enddo
enddo
ir(self%ne+1)=self%cell_nf*self%nc+1
do i=self%ne,1,-1 ! cumulative raw face count
  ir(i)=ir(i+1)-nr(i)
enddo
if(ir(1)/=1)call oft_abort('Bad raw face count','mesh_faces',__FILE__)
!---Loop over cells & index raw faces
nr=0 ! reset raw face counter
do i=1,self%nc
  k(1:self%cell_ne)=abs(self%lce(:,i)) ! cell edges
  do j=1,self%cell_nf ! loop over faces
    js=minval(k(self%cell_fe(:,j))) ! low edge
    je=maxval(k(self%cell_fe(:,j))) ! high edge
    jr(ir(js)+nr(js))=je ! high edge list
    nr(js)=nr(js)+1      ! face counter
  enddo
enddo
!$omp parallel do
do i=1,self%ne ! sort raw face list
  if(nr(i)>1)call sort_array(jr(ir(i):ir(i)+nr(i)-1),nr(i))
enddo
deallocate(nr)
allocate(self%klef(self%ne+1))
self%klef=0 ! initialize unique face counter
!---Loop over raw faces & count unique faces
!$omp parallel do private(j,js,je)
do i=1,self%ne
  js=ir(i)
  je=ir(i+1)-1
  IF(je<js)CYCLE
  self%klef(i)=1
  do j=js+1,je
    if(jr(j)>jr(j-1))self%klef(i)=self%klef(i)+1
  enddo
enddo
self%nf=sum(self%klef) ! total number of unique faces
self%klef(self%ne+1)=self%nf+1
do i=self%ne,1,-1 ! cumulative unique face count
  self%klef(i)=self%klef(i+1)-self%klef(i)
enddo
if(self%klef(1)/=1)call oft_abort('Bad unique edge count','mesh_faces',__FILE__)
allocate(self%llef(self%nf))
self%llef=0
!---Loop over raw faces & index unique faces
!$omp parallel do private(j,js,je,jn)
do i=1,self%ne
  js=ir(i)
  je=ir(i+1)-1
  IF(je<js)CYCLE
  jn=self%klef(i)
  self%llef(jn)=jr(js)
  do j=js+1,je
    if(jr(j)>jr(j-1))then
      jn=jn+1
      self%llef(jn)=jr(j)
    endif
  enddo
enddo
deallocate(ir,jr)
allocate(self%lf(self%face_np,self%nf),self%lfe(self%face_np,self%nf))
allocate(self%lcf(self%cell_nf,self%nc))
self%lf=0
!---Loop over cells & index cells to faces
!$omp parallel do private(k,j,js,je,jp,jn)
do i=1,self%nc
  k(1:self%cell_ne)=abs(self%lce(:,i)) ! cell edges
  do j=1,self%cell_nf ! loop over faces
    js=minval(k(self%cell_fe(:,j))) ! low edge
    je=maxval(k(self%cell_fe(:,j))) ! high edge
    jp=self%klef(js)      ! pointer into high edge list
    jn=self%klef(js+1)-jp ! number of high edges
    self%lcf(j,i)=search_array(je,self%llef(jp:jp+jn-1),jn)+jp-1
    ! Construct face
    !$omp critical
    IF(self%lf(1,self%lcf(j,i))==0)THEN
      DO jp=1,self%face_np
        self%lf(jp,self%lcf(j,i)) = self%lc(self%cell_fc(jp,j),i)
      END DO
      CALL find_orient_listn(oflag,self%lf(:,self%lcf(j,i)),self%face_np)
      CALL orient_listn(oflag,self%lf(:,self%lcf(j,i)),self%face_np)
    END IF
    !$omp end critical
  enddo
enddo
!---Loop over cells & index cells to edges
!$omp parallel do private(k,j,js,je,jp,jn)
do i=1,self%nf
  k(1:self%face_np)=self%lf(:,i) ! cell corners
  do j=1,self%face_np ! loop over edges
    js=MIN(k(self%face_ed(1,j)),k(self%face_ed(2,j))) ! low points
    je=MAX(k(self%face_ed(1,j)),k(self%face_ed(2,j))) ! high points
    jp=self%klpe(js)      ! pointer into low point list
    jn=self%klpe(js+1)-jp ! number of shared edges
    self%lfe(j,i)=search_array(je,self%llpe(jp:jp+jn-1),jn)+jp-1
    self%lfe(j,i)=SIGN(self%lfe(j,i),k(self%face_ed(2,j))-k(self%face_ed(1,j))) ! apply orientation
  enddo
enddo
DEBUG_STACK_POP
end subroutine mesh_faces
!------------------------------------------------------------------------------
!> Locate point, edge, face and cell neighbor cells.
!! - lfc List construction
!! - lcc List construction
!! - kpc, lpc Linkage construction
!! - kec, lec Linkage construction
!------------------------------------------------------------------------------
subroutine mesh_neighbors(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i4), allocatable :: nr(:)
integer(i4) :: i,j,k
DEBUG_STACK_PUSH
allocate(self%lfc(2,self%nf),self%lcc(self%cell_nf,self%nc))
self%lcc=0 ! initialize cell list
self%lfc=0
do i=1,self%nc ! loop over cells & index to faces
  do j=1,self%cell_nf  ! loop over faces
    k=abs(self%lcf(j,i)) ! face numbers
    if(self%lfc(1,k)==0)then
      self%lfc(1,k)=i ! record first cell
    else
      self%lfc(2,k)=i ! record second cell
    endif
  enddo ! **** may need to apply oriented ****
enddo
!$omp parallel do private(j,k)
do i=1,self%nc ! loop over cells & locate neighbors
  do j=1,self%cell_nf  ! loop over faces
    k=abs(self%lcf(j,i)) ! face numbers
    self%lcc(j,i)=sum(self%lfc(:,k))-i
  enddo
enddo
DEBUG_STACK_POP
end subroutine mesh_neighbors
!------------------------------------------------------------------------------
!> Locate and index boundary points, edges, faces and cells.
!! - nbp, nbe, nbf, nbc Counts
!! - bp, be, bf, bc Flag construction
!! - lbp, lbe, lbf, lbc List construction
!------------------------------------------------------------------------------
subroutine mesh_boundary(self)
class(oft_mesh), intent(inout) :: self !< Mesh object
integer(i4) :: i,j
DEBUG_STACK_PUSH
allocate(self%bp(self%np),self%be(self%ne),self%bf(self%nf),self%bc(self%nc))
self%bp=.false.
self%be=.false.
! boundary cells have 1 or more faces on boundary
DO i=1,self%nf
  self%bf(i)=ANY(self%lfc(:,i)==0)
END DO
DO i=1,self%nc
  self%bc(i)=ANY(self%lcc(:,i)==0)
END DO
self%nbf=count(self%bf)
self%nbc=count(self%bc)
allocate(self%bfs(self%nbf))
self%bfs=-1
allocate(self%lbf(self%nbf),self%lbc(self%nbc))
j=0
do i=1,self%nf
  if (self%bf(i)) then
    j=j+1
    self%lbf(j)=i
  endif
enddo
j=0
do i=1,self%nc
  if (self%bc(i)) then
    j=j+1
    self%lbc(j)=i
  endif
enddo
!$omp parallel do private(j)
do i=1,self%nbf
  do j=1,self%face_np
    self%bp(self%lf(j,self%lbf(i)))=.true.
  end do
  do j=1,self%face_np
    self%be(ABS(self%lfe(j,self%lbf(i))))=.true.
  end do
enddo
self%nbp=count(self%bp)
self%nbe=count(self%be)
allocate(self%lbp(self%nbp),self%lbe(self%nbe))
j=0
do i=1,self%np
  if (self%bp(i)) then
    j=j+1
    self%lbp(j)=i
  endif
enddo
j=0
do i=1,self%ne
  if (self%be(i)) then
    j=j+1
    self%lbe(j)=i
  endif
enddo
DEBUG_STACK_POP
end subroutine mesh_boundary
!------------------------------------------------------------------------------
!> Perform local mesh decomposition (METIS)
!------------------------------------------------------------------------------
SUBROUTINE mesh_local_partition(self,tloc_p,tloc_e,tloc_f,tloc_c,nparts)
CLASS(oft_mesh), INTENT(inout) :: self !< Mesh to partition
TYPE(oft_1d_int), DIMENSION(:), INTENT(inout) :: tloc_p !< Point partitioning [self%np]
TYPE(oft_1d_int), DIMENSION(:), INTENT(inout) :: tloc_e !< Edge partitioning [self%ne]
TYPE(oft_1d_int), DIMENSION(:), INTENT(inout) :: tloc_f !< Face partitioning [self%nf]
TYPE(oft_1d_int), DIMENSION(:), INTENT(inout) :: tloc_c !< Cell partitioning [self%nc]
INTEGER(i4), INTENT(in) :: nparts !< Number of partitions
INTEGER(i4) :: i,j,k,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cpart,ppart,epart,fpart
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: town
REAL(r4), ALLOCATABLE, DIMENSION(:) :: trand,own_tmp
DEBUG_STACK_PUSH
!---Decompose domain using METIS
ALLOCATE(cpart(self%nc))
IF(nparts>1)THEN
  CALL oft_metis_partMesh(self%nc,self%np,self%cell_np,self%lc,nparts,cpart,ierr)
  IF(ierr<0)CALL oft_abort('Mesh partitioning failed','mesh_local_partition',__FILE__)
ELSE
  cpart=1
END IF
ALLOCATE(trand(nparts),own_tmp(nparts))
CALL RANDOM_NUMBER(trand)
!---Get cells
DO i=1,nparts
  tloc_c(i)%n=0
END DO
DO i=1,self%nc
  tloc_c(cpart(i))%n=tloc_c(cpart(i))%n+1
END DO
DO i=1,nparts
  ALLOCATE(tloc_c(i)%v(tloc_c(i)%n))
  tloc_c(i)%n=0
  DO j=1,self%nc
    IF(cpart(j)==i)THEN
      tloc_c(i)%n=tloc_c(i)%n+1
      tloc_c(i)%v(tloc_c(i)%n)=j
    END IF
  END DO
END DO
!---Get points
ALLOCATE(ppart(self%np),town(nparts,self%np))
ppart=0
town=.FALSE.
DO i=1,self%np
  DO j=self%kpc(i),self%kpc(i+1)-1
    town(cpart(self%lpc(j)),i)=.TRUE.
  END DO
END DO
DO i=1,self%np
  own_tmp=trand
  DO j=self%kpp(i),self%kpp(i+1)-1
    DO k=1,nparts
      IF(town(k,self%lpp(j)))own_tmp(k)=own_tmp(k)+1.
    END DO
  END DO
  ppart(i)=MAXLOC(own_tmp,DIM=1)
END DO
DEALLOCATE(town)
DO i=1,nparts
  tloc_p(i)%n=0
END DO
DO i=1,self%np
  tloc_p(ppart(i))%n=tloc_p(ppart(i))%n+1
END DO
DO i=1,nparts
  ALLOCATE(tloc_p(i)%v(tloc_p(i)%n))
  tloc_p(i)%n=0
  DO j=1,self%np
    IF(ppart(j)==i)THEN
      tloc_p(i)%n=tloc_p(i)%n+1
      tloc_p(i)%v(tloc_p(i)%n)=j
    END IF
  END DO
END DO
!---Get edges
ALLOCATE(epart(self%ne),town(nparts,self%ne))
epart=0
town=.FALSE.
DO i=1,self%ne
  DO j=self%kec(i),self%kec(i+1)-1
    town(cpart(self%lec(j)),i)=.TRUE.
  END DO
END DO
DO i=1,self%ne
  own_tmp=trand
  DO j=self%kee(i),self%kee(i+1)-1
    DO k=1,nparts
      IF(town(k,self%lee(j)))own_tmp(k)=own_tmp(k)+1.
    END DO
  END DO
  epart(i)=MAXLOC(own_tmp,DIM=1)
END DO
DEALLOCATE(town)
DO i=1,nparts
  tloc_e(i)%n=0
END DO
DO i=1,self%ne
  tloc_e(epart(i))%n=tloc_e(epart(i))%n+1
END DO
DO i=1,nparts
  ALLOCATE(tloc_e(i)%v(tloc_e(i)%n))
  tloc_e(i)%n=0
  DO j=1,self%ne
    IF(epart(j)==i)THEN
      tloc_e(i)%n=tloc_e(i)%n+1
      tloc_e(i)%v(tloc_e(i)%n)=j
    END IF
  END DO
END DO
DEALLOCATE(epart)
!---Get faces
ALLOCATE(fpart(self%nf))
fpart=0
DO i=1,self%nf
  own_tmp=trand
  own_tmp(cpart(self%lfc(1,i)))=own_tmp(cpart(self%lfc(1,i)))+1.
  IF(self%lfc(2,i)>0)own_tmp(cpart(self%lfc(2,i)))=own_tmp(cpart(self%lfc(2,i)))+1.
  DO j=1,self%face_np
    own_tmp(ppart(self%lf(j,i)))=own_tmp(ppart(self%lf(j,i)))+1.
  END DO
  fpart(i)=MAXLOC(own_tmp,DIM=1)
END DO
DO i=1,nparts
  tloc_f(i)%n=0
END DO
DO i=1,self%nf
  tloc_f(fpart(i))%n=tloc_f(fpart(i))%n+1
END DO
DO i=1,nparts
  ALLOCATE(tloc_f(i)%v(tloc_f(i)%n))
  tloc_f(i)%n=0
  DO j=1,self%nf
    IF(fpart(j)==i)THEN
      tloc_f(i)%n=tloc_f(i)%n+1
      tloc_f(i)%v(tloc_f(i)%n)=j
    END IF
  END DO
END DO
DEALLOCATE(ppart,fpart,cpart,trand,own_tmp)
DEBUG_STACK_POP
END SUBROUTINE mesh_local_partition
!------------------------------------------------------------------------------
!> Driver for grid construction
!------------------------------------------------------------------------------
subroutine bmesh_local_init(self,parent,sync_normals)
class(oft_bmesh), intent(inout) :: self !< Mesh object
class(oft_mesh), optional, intent(in) :: parent !< Parent volume mesh (if present)
LOGICAL, OPTIONAL, INTENT(IN) :: sync_normals !< Sync unit normal directions between faces?
REAL(r8) :: f(4),ftmp(4),gop(3,4),vol
REAL(r8), ALLOCATABLE :: fn(:,:)
integer(i4) :: i,j,k,l,ind,ed(2)
integer(i4), ALLOCATABLE, DIMENSION(:) :: fc
! IF(PRESENT(parent))THEN
!   IF(.NOT.ASSOCIATED(self%parent))ALLOCATE(self%parent)
! END IF
IF(self%skip)RETURN
DEBUG_STACK_PUSH
ALLOCATE(self%ho_info%is_curved(self%nc))
self%ho_info%is_curved=.FALSE.
IF(PRESENT(parent))THEN
  self%parent%np=parent%np
  self%parent%ne=parent%ne
  self%parent%nf=parent%nf
  IF(.NOT.ASSOCIATED(self%parent%lf))THEN
    ALLOCATE(fc(self%cell_np))
    ALLOCATE(self%parent%lf(self%nc))
    DO i=1,self%nc
      fc=INT(self%parent%lp(self%lc(:,i)),4)
      self%parent%lf(i)=ABS(mesh_local_findface(parent,fc))
      IF(self%parent%lf(i)==0)THEN
        WRITE(*,*)i,fc
        CALL oft_abort("Bad face link", "", __FILE__)
      END IF
    END DO
    DEALLOCATE(fc)
  END IF
  !---
  ALLOCATE(fn(3,self%nc))
  DO i=1,self%nc
    j=INT(self%parent%lf(i),4)
    k=parent%lfc(1,j)
    DO ind=1,parent%cell_nf
      IF(parent%lcf(ind,k)==j)EXIT
    END DO
    f=0.d0
    DO l=1,self%cell_np
      CALL parent%vlog(parent%cell_fc(l,ind), ftmp)
      f=f+ftmp
    END DO
    f=f/REAL(self%cell_np,8)
    CALL parent%snormal(k, ind, f, fn(:,i))
    fn(:,i)=-fn(:,i)
  END DO
  !
  call bmesh_areas(self,fn)
  DEALLOCATE(fn)
ELSE
  call bmesh_areas(self)
END IF
call amesh_edges(self)
call amesh_to_cell(self)
call bmesh_neighbors(self)
IF(PRESENT(sync_normals))THEN
  IF(sync_normals)CALL sync_face_normals(self)
END IF
call bmesh_boundary(self)
call amesh_interactions(self)
!---
IF(PRESENT(parent))THEN
  !---
  ALLOCATE(self%parent%le(self%ne))
  DO i=1,self%ne
    ed=INT(self%parent%lp(self%le(:,i)),4)
    self%parent%le(i)=ABS(mesh_local_findedge(parent,ed))
  END DO
! ELSE
!   self%global%np=self%np
!   ALLOCATE(self%global%lp(self%np))
!   DO i=1,self%np
!     self%global%lp(i)=INT(i,8)
!   END DO
!   self%global%ne=self%ne
!   ALLOCATE(self%global%le(self%ne))
!   DO i=1,self%ne
!     self%global%le(i)=INT(i,8)
!   END DO
!   self%global%nc=self%nc
!   ALLOCATE(self%global%lc(self%nc))
!   DO i=1,self%nc
!     self%global%lc(i)=INT(i,8)
!   END DO
!   ALLOCATE(self%lco(self%nc))
!   self%lco=1
!   ALLOCATE(self%pstitch%leo(self%nbp))
!   ALLOCATE(self%estitch%leo(self%nbe))
!   self%pstitch%leo=.TRUE.
!   self%estitch%leo=.TRUE.
END IF
! IF(.NOT.ASSOCIATED(self%tloc_c))THEN
!   self%nparts=oft_env%nthreads*oft_env%nparts
!   ALLOCATE(self%tloc_p(self%nparts))
!   ALLOCATE(self%tloc_e(self%nparts))
!   ALLOCATE(self%tloc_f(self%nparts))
!   ALLOCATE(self%tloc_c(self%nparts))
!   call mesh_local_partition(self,self%tloc_p,self%tloc_e,self%tloc_f,self%tloc_c,self%nparts)
! END IF
DEBUG_STACK_POP
end subroutine bmesh_local_init
!------------------------------------------------------------------------------
!> Compute cell and vertex volumes, rectifying negative volume cells.
!!
!! Zero volume cells or vertices are also caught for mesh validation.
!------------------------------------------------------------------------------
subroutine bmesh_areas(self,fn)
class(oft_bmesh), INTENT(inout) :: self !< Mesh object
REAL(r8), OPTIONAL, INTENT(in) :: fn(:,:) !< Unit normal directions for orientation matching [3,nc]
integer(i4) :: i,j
real(r8) :: area,f(3),gop(3,4),norm(3)
IF(self%nc==0)RETURN
DEBUG_STACK_PUSH
if(ASSOCIATED(self%ca))then
  if(SIZE(self%ca).NE.self%nc)ALLOCATE(self%ca(self%nc),self%va(self%np))
else
  allocate(self%ca(self%nc),self%va(self%np))
endif
self%va=0.d0 ! initialize vertex volumes
! !$omp parallel do private(k,e1,e2,n,area)
f = 1.d0/3.d0 ! NEED TO UPDATE
do i=1,self%nc ! loop over cells
  CALL self%jacobian(i,f,gop,self%ca(i))
  CALL self%norm(i,f,norm)
  !---Orient cell normals to outward direction
  IF(PRESENT(fn))THEN
    IF(DOT_PRODUCT(fn(:,i),norm)<0.d0)CALL self%invert_face(i)
  ELSE IF(self%dim==2)THEN
    IF(norm(3)<0.d0)CALL self%invert_face(i)
  END IF
  self%ca(i)=ABS(self%ca(i))
enddo
do i=1,self%nc ! accumulate vertex volumes
  area=self%ca(i)/REAL(self%cell_np,8) ! corner volume
  do j=1,self%cell_np
    self%va(self%lc(j,i))=self%va(self%lc(j,i))+area
  enddo
enddo
!---Check for zero volume cells
if(ANY(self%ca==0))call oft_abort('Zero volume cell detected','bmesh_areas',__FILE__)
!---Check for zero volume vertices
if(.NOT.ASSOCIATED(self%parent).AND.ANY(self%va==0))THEN
  call oft_abort('Floating vertex detected','bmesh_areas',__FILE__)
endif
DEBUG_STACK_POP
end subroutine bmesh_areas
!------------------------------------------------------------------------------
!> Compute cell and vertex volumes, rectifying negative volume cells.
!!
!! Zero volume cells or vertices are also caught for mesh validation.
!------------------------------------------------------------------------------
subroutine sync_face_normals(self)
class(oft_bmesh), INTENT(inout) :: self !< Mesh object
integer(i4) :: i,max_depth=1
logical, ALLOCATABLE, DIMENSION(:) :: oriented
IF(self%nc==0)RETURN
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Ensuring surface normal orientations'
CALL oft_increase_indent
ALLOCATE(oriented(self%nc))
oriented=.FALSE.
do i=1,self%nc ! loop over cells
  IF(.NOT.oriented(i))THEN
    oriented(i)=.TRUE.
    CALL orient_neighbors(i,1)
    IF(oft_debug_print(2))WRITE(*,'(2A,I8,I8)')oft_indent,'Chunk oriented',i,COUNT(oriented)
  END IF
enddo
WRITE(*,*)'Orientation depth =',max_depth
CALL oft_decrease_indent
IF(COUNT(oriented)/=self%nc)CALL oft_abort("Orientation failed", &
  "sync_face_normals",__FILE__)
DEALLOCATE(oriented)
CONTAINS
recursive subroutine orient_neighbors(face1,depth)
integer(i4), intent(in) :: face1,depth
integer(i4) :: j,face2,k,ed1(2),ed2(2)
! logical, allocatable :: mark(:)
max_depth=MAX(depth,max_depth)
! allocate(mark(self%cell_ne))
! mark=.FALSE.
DO j=1,self%cell_ne
  face2=self%lcc(j,face1)
  IF(face2==0)CYCLE
  IF(oriented(face2))CYCLE
  ed1=self%lc(self%cell_ed(:,j),face1)
  !---Ensure same sense as neighbor (opposite direction of shared edge)
  DO k=1,self%cell_ne
    IF(self%lcc(k,face2)==face1)EXIT
  END DO
  ed2=self%lc(self%cell_ed(:,k),face2)
  IF(ALL(ed1==ed2))THEN
    CALL self%invert_face(face2)
  END IF
  oriented(face2)=.TRUE.
  CALL orient_neighbors(face2,depth+1)
  ! mark(j)=.TRUE.
END DO
! DO j=1,3
!   IF(mark(j))CALL orient_neighbors(self%lcc(j,face1),depth+1)
! END DO
! DEALLOCATE(mark)
end subroutine orient_neighbors
end subroutine sync_face_normals
!------------------------------------------------------------------------------
!> Locate point, edge, face and cell neighbor cells.
!! - lfc List construction
!! - lcc List construction
!! - kpc, lpc Linkage construction
!! - kec, lec Linkage construction
!------------------------------------------------------------------------------
subroutine bmesh_neighbors(self)
class(oft_bmesh), INTENT(inout) :: self !< Mesh object
integer(i4), allocatable :: nr(:)
integer(i4) :: i,j,k
DEBUG_STACK_PUSH
allocate(self%lcc(self%cell_np,self%nc))
self%lcc=0 ! initialize cell list
!$omp parallel do private(j,k)
do i=1,self%nc ! loop over faces & locate neighbors
  do j=1,self%cell_np  ! loop over edges
    k=ABS(self%lce(j,i)) ! edge numbers
    self%lcc(j,i)=SUM(self%lec(self%kec(k):self%kec(k+1)-1))-i
  enddo
enddo
DEBUG_STACK_POP
end subroutine bmesh_neighbors
!------------------------------------------------------------------------------
!> Locate and index boundary points, edges, faces and cells.
!! - nbp, nbe, nbf, nbc Counts
!! - bp, be, bf, bc Flag construction
!! - lbp, lbe, lbf, lbc List construction
!------------------------------------------------------------------------------
subroutine bmesh_boundary(self)
class(oft_bmesh), INTENT(inout) :: self !< Mesh object
integer(i4) :: i,j
DEBUG_STACK_PUSH
allocate(self%bp(self%np),self%be(self%ne),self%bc(self%nc))
self%bp=.false.
! boundary faces have 1 or more edge on boundary
DO i=1,self%ne
  self%be(i)=(self%kec(i+1)-self%kec(i)==1)
END DO
DO i=1,self%nc
  self%bc(i)=ANY(self%lcc(:,i)==0)
END DO
self%nbe=COUNT(self%be)
self%nbc=COUNT(self%bc)
allocate(self%bes(self%nbe))
self%bes=-1
allocate(self%lbe(self%nbe),self%lbc(self%nbc))
j=0
do i=1,self%ne
  if (self%be(i)) then
    j=j+1
    self%lbe(j)=i
  endif
enddo
j=0
do i=1,self%nc
  if (self%bc(i)) then
    j=j+1
    self%lbc(j)=i
  endif
enddo
!$omp parallel do private(j)
do i=1,self%nbe
  do j=1,2
    self%bp(self%le(j,self%lbe(i)))=.TRUE.
  enddo
enddo
self%nbp=COUNT(self%bp)
allocate(self%lbp(self%nbp))
j=0
do i=1,self%np
  if (self%bp(i)) then
    j=j+1
    self%lbp(j)=i
  endif
enddo
DEBUG_STACK_POP
end subroutine bmesh_boundary
end module oft_mesh_local
