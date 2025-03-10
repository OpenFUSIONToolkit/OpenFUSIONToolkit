!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file trimesh_tessellation.F90
!
!> Subroutines for nested tesselation of a triangular mesh.
!!
!! @author Chris Hansen
!! @date Summer 2012
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------
MODULE trimesh_tessellation
USE oft_base
USE oft_mesh_type, ONLY: oft_bmesh
IMPLICIT NONE
#include "local.h"
PRIVATE
!---Quadratic tessellation
INTEGER(i4), PARAMETER :: tess2(3,4)=RESHAPE((/5,3,1,4,3,2,4,5,0,4,5,3/),(/3,4/))
!---Cubic tessellation
INTEGER(i4), PARAMETER :: tess3(3,9)=RESHAPE((/6,4,2,6,4,9,5,8,9,7,4,9,7,5,0,7,5,9,3,6, &
9,3,8,1,3,8,9/),(/3,9/))
!---Quartic tessellation
INTEGER(i4), PARAMETER :: tess4(3,16) = RESHAPE((/4,7,12,9,6,12,9,4,12,9,4,2,14,5,8,14,7, &
12,13,6,12,13,14,12,13,11,8,13,14,8,10,14,7,10,5,0,10,14,5,3,11,1,3,13,11,3,13,6/),(/3,16/))
!---
PUBLIC tessellate1, tessellate2, tessellate3, tessellate4
CONTAINS
!------------------------------------------------------------------------------
!> Construct point and face lists for a single tessellation
!!
!! This corresponds to the input triangulation with no refinement
!! - np_tess = np
!! - nc_tess = nc
!------------------------------------------------------------------------------
SUBROUTINE tessellate1(self,rtmp,lctmp)
CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
INTEGER(i4) :: i,np,nc
DEBUG_STACK_PUSH
!---
np = self%np
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!---
nc = self%nc
ALLOCATE(lctmp(3,nc))
!$omp parallel do
DO i=1,self%nc
  lctmp(:,i)=self%lc(:,i)-1
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate1
!------------------------------------------------------------------------------
!> Construct point and face lists for a 2 level tessellation
!!
!! This corresponds to the input triangulation with no refinement
!! - np_tess = np + ne
!! - nc_tess = 4*nc
!------------------------------------------------------------------------------
SUBROUTINE tessellate2(self,rtmp,lctmp)
CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,1)=RESHAPE((/.5d0,.5d0/),(/2,1/))
INTEGER(i4), PARAMETER :: tess(3,4)=RESHAPE((/6,4,2,5,4,3,5,6,1,5,6,4/),(/3,4/))
REAL(r8) :: f(3),pt(3)
INTEGER(i4) :: i,j,ind,np,nc,ed(2),dofs(6)
DEBUG_STACK_PUSH
!---
np = self%np+self%ne
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!$omp parallel do private(j,ind,ed,f,pt)
DO i=1,self%ne
  j=self%lec(self%kec(i))
  DO ind=1,3
    IF(ABS(self%lce(ind,j))==i)EXIT
  END DO
  ed=self%cell_ed(:,ind)
  f=0.d0
  f(ed(1))=ed_dofs(1,1)
  f(ed(2))=ed_dofs(2,1)
  pt=self%log2phys(j,f)
  rtmp(:,i+self%np)=pt
END DO
!---
nc = self%nc*4
ALLOCATE(lctmp(3,nc))
!$omp parallel do private(dofs,j)
DO i=1,self%nc
  dofs(1:3)=self%lc(:,i)
  DO j=1,3
    dofs(j+3)=ABS(self%lce(j,i))+self%np
  END DO
  DO j=1,4
    lctmp(:,(i-1)*4+j)=dofs(tess2(:,j)+1)-1
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate2
!------------------------------------------------------------------------------
!> Construct point and face lists for a 3 level tessellation
!!
!! This corresponds to the input triangulation with no refinement
!! - np_tess = np + 2*ne + nc
!! - nc_tess = 9*nc
!------------------------------------------------------------------------------
SUBROUTINE tessellate3(self,rtmp,lctmp)
CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,2)=RESHAPE((/1.d0/3.d0,2.d0/3.d0,2.d0/3.d0,1.d0/3.d0/),(/2,2/))
REAL(r8), PARAMETER :: fc_dofs(3,1)=RESHAPE((/1.d0/3.d0,1.d0/3.d0,1.d0/3.d0/),(/3,1/))
REAL(r8) :: f(3),pt(3)
INTEGER(i4) :: i,j,ind,np,nc,ed(2),fc(3),dofs(10)
DEBUG_STACK_PUSH
!---
np = self%np+2*self%ne+self%nc
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!$omp parallel do private(j,ind,ed,f,pt)
DO i=1,self%ne
  j=self%lec(self%kec(i))
  DO ind=1,3
    IF(ABS(self%lce(ind,j))==i)EXIT
  END DO
  ed=self%cell_ed(:,ind)
  CALL orient_list2(self%lce(ind,j),ed)
  f=0.d0
  f(ed(1))=ed_dofs(1,1)
  f(ed(2))=ed_dofs(2,1)
  pt=self%log2phys(j,f)
  rtmp(:,1+(i-1)*2+self%np)=pt
  f(ed(1))=ed_dofs(1,2)
  f(ed(2))=ed_dofs(2,2)
  pt=self%log2phys(j,f)
  rtmp(:,i*2+self%np)=pt
END DO
!$omp parallel do private(fc,f,pt)
DO i=1,self%nc
  fc=(/1,2,3/)
  CALL orient_listn(self%lco(i),fc,3_i4)
  f=0.d0
  f(fc(1))=fc_dofs(1,1)
  f(fc(2))=fc_dofs(2,1)
  f(fc(3))=fc_dofs(3,1)
  pt=self%log2phys(i,f)
  rtmp(:,i+self%np+2*self%ne)=pt
END DO
!---
nc = self%nc*9
ALLOCATE(lctmp(3,nc))
!$omp parallel do private(dofs,j)
DO i=1,self%nc
  dofs(1:3)=self%lc(:,i)
  DO j=1,3
    IF(self%lce(j,i)>0)THEN
      dofs(j+3)=(ABS(self%lce(j,i))-1)*2+self%np+1
      dofs(j+6)=ABS(self%lce(j,i))*2+self%np
    ELSE
      dofs(j+3)=ABS(self%lce(j,i))*2+self%np
      dofs(j+6)=(ABS(self%lce(j,i))-1)*2+self%np+1
    END IF
  END DO
  dofs(10)=i+self%np+2*self%ne
  DO j=1,9
    lctmp(:,(i-1)*9+j)=dofs(tess3(:,j)+1)-1
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate3
!------------------------------------------------------------------------------
!> Construct point and face lists for a 4 level tessellation
!!
!! This corresponds to the input triangulation with no refinement
!! - np_tess = np + 3*ne + 3*nc
!! - nc_tess = 16*nc
!------------------------------------------------------------------------------
SUBROUTINE tessellate4(self,rtmp,lctmp)
CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,3)=RESHAPE((/1.d0/4,3.d0/4,1.d0/2,1.d0/2,3.d0/4,1.d0/4/),(/2,3/))
REAL(r8), PARAMETER :: fc_dofs(3,3)=RESHAPE((/1.d0/4,1.d0/4,1.d0/2,1.d0/4,1.d0/2,1.d0/4, &
1.d0/2,1.d0/4,1.d0/4/),(/3,3/))
REAL(r8) :: f(3),pt(3)
INTEGER(i4) :: i,j,m,mm,ind,np,nc,ed(2),fc(3),dofs(15)
DEBUG_STACK_PUSH
!---
np = self%np+3*self%ne+3*self%nc
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!$omp parallel do private(j,ind,ed,f,pt,m)
DO i=1,self%ne
  j=self%lec(self%kec(i))
  DO ind=1,3
    IF(ABS(self%lce(ind,j))==i)EXIT
  END DO
  ed=self%cell_ed(:,ind)
  CALL orient_list2(self%lce(ind,j),ed)
  f=0.d0
  DO m=1,3
    f(ed(1))=ed_dofs(1,m)
    f(ed(2))=ed_dofs(2,m)
    pt=self%log2phys(j,f)
    rtmp(:,m+(i-1)*3+self%np)=pt
  END DO
END DO
!$omp parallel do private(fc,f,pt,m)
DO i=1,self%nc
  fc=(/1,2,3/)
  CALL orient_listn(self%lco(i),fc,3_i4)
  f=0.d0
  DO m=1,3
    f(fc(1))=fc_dofs(1,m)
    f(fc(2))=fc_dofs(2,m)
    f(fc(3))=fc_dofs(3,m)
    pt=self%log2phys(i,f)
    rtmp(:,m+(i-1)*3+3*self%ne+self%np)=pt
  END DO
END DO
!---
nc = self%nc*16
ALLOCATE(lctmp(3,nc))
!$omp parallel do private(dofs,j,m,mm,fc)
DO i=1,self%nc
  dofs(1:3)=self%lc(:,i)
  DO j=1,3
    DO m=1,3
      IF(self%lce(j,i)<0)THEN
        dofs(j+(3-m)*3+3)=m+(ABS(self%lce(j,i))-1)*3+self%np
      ELSE
        dofs(j+(m-1)*3+3)=m+(ABS(self%lce(j,i))-1)*3+self%np
      END IF
    END DO
  END DO
  fc=(/1,2,3/)
  CALL orient_listn(self%lco(i),fc,3_i4)
  DO m=1,3
    DO mm=1,3
      IF(ALL(fc_dofs(fc,m)==fc_dofs((/1,2,3/),mm)))EXIT
    END DO
    dofs(m+3*3+3)=mm+self%np+3*self%ne+(i-1)*3
  END DO
  DO j=1,16
    lctmp(:,(i-1)*16+j)=dofs(tess4(:,j)+1)-1
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate4
END MODULE trimesh_tessellation
