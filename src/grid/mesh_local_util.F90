!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_local_util.F90
!
!> Mesh utility functions, entity location and mappings.
!!
!! Utility functions.
!! - Edge, Face location
!! - Local to boundary index mapping
!! - Mesh copying
!!
!! @authors George Marklin and Chris Hansen
!! @date April 2008 - Present
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_local_util
USE oft_base
USE oft_sort, ONLY: search_array
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh
IMPLICIT NONE
#include "local.h"
INTEGER(i4) :: oriented_cell=-1
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: oriented_edges
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: oriented_faces
!$omp threadprivate(oriented_cell,oriented_edges,oriented_faces)
contains
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine mesh_local_orient(self,cell)
class(oft_mesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell
integer(i4) :: i
DEBUG_STACK_PUSH
IF(.NOT.ALLOCATED(oriented_edges))THEN
  ALLOCATE(oriented_edges(2,self%cell_ne))
  ALLOCATE(oriented_faces(self%face_np,self%cell_nf))
END IF
DO i=1,self%cell_ne
  oriented_edges(:,i)=self%cell_ed(:,i)
  call orient_list2(self%lce(i,cell),oriented_edges(:,i))
END DO
DO i=1,self%cell_nf
  oriented_faces(:,i)=self%cell_fc(:,i)
  call orient_listn(self%lcfo(i,cell), oriented_faces(:,i), self%face_np)
END DO
oriented_cell=cell
DEBUG_STACK_POP
end subroutine mesh_local_orient
!------------------------------------------------------------------------------
!> Find oriented edge connecting points i1 and i2.
!!
!! Orientation is (+,-) if [imin,imax] is [i1,i2] or [i2,i1]
!------------------------------------------------------------------------------
function mesh_local_findedge(self,inds)
class(oft_amesh), intent(in) :: self !< Mesh to search
integer(i4), intent(in) :: inds(2) !< Edge points (ordered)
integer(i4) :: mesh_local_findedge !< Oriented edge linking inds(1) -> inds(2), zero if no edge exists
integer(i4) :: js,je,jp,jn
DEBUG_STACK_PUSH
js=minval(inds) ! low point
je=maxval(inds) ! high point
jp=self%klpe(js)   ! pointer into low point list
jn=self%klpe(js+1)-jp ! number of shared edges
js=search_array(je,self%llpe(jp:jp+jn-1),jn)
if(js/=0)js=sign(js+jp-1,inds(2)-inds(1))
mesh_local_findedge=js
DEBUG_STACK_POP
end function mesh_local_findedge
!------------------------------------------------------------------------------
!> Find face with corner points `inds` (no orientation is applied)
!------------------------------------------------------------------------------
function mesh_local_findface(self,inds)
class(oft_mesh), intent(in) :: self !< Mesh to search
integer(i4), intent(in) :: inds(:) !< Face points (ordered)
integer(i4) :: mesh_local_findface !< Face composed of `inds`, zero if no face exists
integer(i4) :: ilo,imi,ihi,js,je,jp,jn,etmp(2),jsh(4)
DEBUG_STACK_PUSH
mesh_local_findface=0
IF(self%face_np==3)THEN
  ilo=minval(inds) ! low point
  ihi=maxval(inds) ! high point
  imi=SUM(inds)-ilo-ihi ! middle point
  etmp=(/ilo,imi/)
  js=ABS(mesh_local_findedge(self,etmp)) ! low edge
  if(js==0)THEN
    DEBUG_STACK_POP
    RETURN
  END IF
  etmp=(/imi,ihi/)
  je=ABS(mesh_local_findedge(self,etmp)) ! high edge
  if(je==0)THEN
    DEBUG_STACK_POP
    RETURN
  END IF
ELSE
  jsh(1)=ABS(mesh_local_findedge(self,(/inds(1), inds(2)/)))
  jsh(2)=ABS(mesh_local_findedge(self,(/inds(2), inds(3)/)))
  jsh(3)=ABS(mesh_local_findedge(self,(/inds(3), inds(4)/)))
  jsh(4)=ABS(mesh_local_findedge(self,(/inds(4), inds(1)/)))
  IF(ANY(jsh==0))THEN
    DEBUG_STACK_POP
    RETURN
  END IF
  js=MINVAL(jsh); je=MAXVAL(jsh)
END IF
jp=self%klef(js)      ! pointer into low edge list
jn=self%klef(js+1)-jp ! number of shared faces
js=search_array(je,self%llef(jp:jp+jn-1),jn)
if(js/=0)js=js+jp-1
mesh_local_findface=js
DEBUG_STACK_POP
end function mesh_local_findface
end module oft_mesh_local_util
