!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file tetmesh_tessellation.F90
!
!> Subroutines for nested tesselation of a tetrahedral mesh.
!!
!! @author Chris Hansen
!! @date June 2012
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE tetmesh_tessellation
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
IMPLICIT NONE
#include "local.h"
PRIVATE
!---Quadratic tessellation
INTEGER(i4), PARAMETER :: tess2(4,8) = RESHAPE((/6,5,4,3,7,6,8,2,9,8,4,0,9,7,5,1,7,6,8,4, &
7,6,5,4,9,7,5,8,9,5,8,4/),(/4,8/))
!---Cubic tessellation
INTEGER(i4), PARAMETER :: tess3(4,27) = RESHAPE((/16,17,19,18,9,15,19,18,10,4,17,18,7,12, &
14,2,11,5,16,18,8,10,15,0,8,14,17,19,6,12,16,17,6,5,4,3,13,11,9,1,13,7,16,19,5,4,17,18,5, &
16,17,18,7,14,17,19,7,16,17,19,7,12,14,17,7,12,16,17,8,10,15,18,8,15,19,18,8,17,19,18,8,10, &
17,18,6,16,4,17,6,5,16,4,13,11,9,18,13,9,19,18,13,16,19,18,13,11,16,18/),(/4,27/))
!---Quartic tessellation
INTEGER(i4), PARAMETER :: tess4(4,64) = RESHAPE((/7,26,13,25,22,23,24,34,22,26,30,34,22,5, &
11,24,31,26,25,34,31,23,27,34,32,33,27,34,8,33,14,27,10,4,23,24,10,32,16,27,6,4,5,3,6,22,12,23, &
20,31,14,25,18,31,12,26,18,20,7,2,17,19,9,1,21,32,33,15,21,8,16,0,29,33,25,34,29,19,13,30,28,29, &
30,34,28,32,24,34,28,17,11,30,28,29,9,15,26,30,25,34,26,13,30,25,33,27,25,34,33,14,27,25,22,4,5, &
24,22,4,23,24,22,30,24,34,22,11,30,24,31,22,23,34,31,22,26,34,31,22,12,23,31,22,12,26,31,27,25, &
34,31,14,27,25,8,32,16,27,8,32,33,27,10,23,24,34,10,23,27,34,10,32,24,34,10,32,27,34,6,22,5,23,6, &
4,5,23,20,7,26,25,20,31,26,25,18,20,31,7,18,31,7,26,21,8,33,16,21,32,33,16,29,13,30,25,29,30,25, &
34,28,32,33,34,28,29,33,34,28,17,19,30,28,29,19,30,28,29,33,15,28,32,33,15,28,29,19,9,28,17,19,9, &
28,30,24,34,28,11,30,24/),(/4,64/))
!---
PUBLIC tessellate1, tessellate2, tessellate3, tessellate4
CONTAINS
!---------------------------------------------------------------------------------
!> Construct point and cell lists for a single tessellation
!!
!! This corresponds to the input triangulation with no refinement
!! - np_tess = np
!! - nc_tess = nc
!---------------------------------------------------------------------------------
SUBROUTINE tessellate1(self,rtmp,lctmp)
CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
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
ALLOCATE(lctmp(4,nc))
!$omp parallel do
DO i=1,self%nc
  lctmp(:,i)=self%lc(:,i)-1
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate1
!---------------------------------------------------------------------------------
!> Construct point and cell lists for a 2 level tessellation
!!
!! This corresponds to the input triangulation with quadratic nodes added
!! - np_tess = np + ne
!! - nc_tess = nc*8
!---------------------------------------------------------------------------------
SUBROUTINE tessellate2(self,rtmp,lctmp)
CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,1)=RESHAPE((/.5d0,.5d0/),(/2,1/))
REAL(r8) :: f(4),pt(3),diag(3)
INTEGER(i4) :: i,j,k,ind,np,nc,ed(2),dofs(10),lccors(4),lcecors(6)
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
  DO ind=1,6
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
nc = self%nc*8
ALLOCATE(lctmp(4,nc))
! $omp parallel do private(lccors,lcecors,diag,k,nfskip,ncskip,j)
do i=1,self%nc
  lccors=self%lc(:,i) ! corner points
  lcecors=ABS(self%lce(:,i))+self%np ! edge points
  lctmp(:,(i-1)*8+1)=(/ lccors(1),lcecors(6),lcecors(5),lcecors(1)/)
  lctmp(:,(i-1)*8+2)=(/lcecors(6), lccors(2),lcecors(4),lcecors(2)/)
  lctmp(:,(i-1)*8+3)=(/lcecors(5),lcecors(4), lccors(3),lcecors(3)/)
  lctmp(:,(i-1)*8+4)=(/lcecors(1),lcecors(2),lcecors(3), lccors(4)/)
  diag(1)=magnitude(rtmp(:,lcecors(1))-rtmp(:,lcecors(4)))
  diag(2)=magnitude(rtmp(:,lcecors(2))-rtmp(:,lcecors(5)))
  diag(3)=magnitude(rtmp(:,lcecors(3))-rtmp(:,lcecors(6)))
  k=MINLOC(diag, DIM=1)
  select case(k)
    case(1) ! place diagonal from edge 1 --> 4
      lctmp(:,(i-1)*8+5)=(/lcecors(5),lcecors(6),lcecors(4),lcecors(1)/)
      lctmp(:,(i-1)*8+6)=(/lcecors(6),lcecors(2),lcecors(4),lcecors(1)/)
      lctmp(:,(i-1)*8+7)=(/lcecors(3),lcecors(5),lcecors(4),lcecors(1)/)
      lctmp(:,(i-1)*8+8)=(/lcecors(2),lcecors(3),lcecors(4),lcecors(1)/)
    case(2) ! place diagonal from edge 2 --> 5
      lctmp(:,(i-1)*8+5)=(/lcecors(1),lcecors(6),lcecors(5),lcecors(2)/)
      lctmp(:,(i-1)*8+6)=(/lcecors(6),lcecors(4),lcecors(5),lcecors(2)/)
      lctmp(:,(i-1)*8+7)=(/lcecors(4),lcecors(3),lcecors(5),lcecors(2)/)
      lctmp(:,(i-1)*8+8)=(/lcecors(3),lcecors(1),lcecors(5),lcecors(2)/)
    case(3) ! place diagonal from edge 3 --> 6
      lctmp(:,(i-1)*8+5)=(/lcecors(5),lcecors(1),lcecors(6),lcecors(3)/)
      lctmp(:,(i-1)*8+6)=(/lcecors(2),lcecors(4),lcecors(6),lcecors(3)/)
      lctmp(:,(i-1)*8+7)=(/lcecors(4),lcecors(5),lcecors(6),lcecors(3)/)
      lctmp(:,(i-1)*8+8)=(/lcecors(1),lcecors(2),lcecors(6),lcecors(3)/)
    case default
      call oft_abort('Invalid cell refinement','tessellate2',__FILE__)
  end select
enddo
lctmp=lctmp-1
! !$omp parallel do private(dofs,j)
! DO i=1,self%nc
!   dofs(1:4)=self%lc(:,i)
!   DO j=1,6
!     dofs(j+4)=ABS(self%lce(j,i))+self%np
!   END DO
!   DO j=1,8
!     lctmp(:,(i-1)*8+j)=dofs(tess2(:,j)+1)-1
!   END DO
! END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate2
!---------------------------------------------------------------------------------
!> Construct point and cell lists for a 3 level tessellation
!!
!! This corresponds to the input triangulation with cubic nodes added
!! - np_tess = np + 2*ne + nf
!! - nc_tess = nc*27
!---------------------------------------------------------------------------------
SUBROUTINE tessellate3(self,rtmp,lctmp)
CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,2)=RESHAPE((/1.d0/3.d0,2.d0/3.d0,2.d0/3.d0,1.d0/3.d0/),(/2,2/))
REAL(r8), PARAMETER :: fc_dofs(3,1)=RESHAPE((/1.d0/3.d0,1.d0/3.d0,1.d0/3.d0/),(/3,1/))
REAL(r8) :: f(4),pt(3)
INTEGER(i4) :: i,j,m,ind,np,nc,ed(2),fc(3),dofs(20)
DEBUG_STACK_PUSH
!---
np = self%np+2*self%ne+self%nf
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!$omp parallel do private(j,ind,ed,f,pt)
DO i=1,self%ne
  j=self%lec(self%kec(i))
  DO ind=1,6
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
!$omp parallel do private(j,ind,fc,f,pt)
DO i=1,self%nf
  j=self%lfc(1,i)
  DO ind=1,4
    IF(ABS(self%lcf(ind,j))==i)EXIT
  END DO
  fc=self%cell_fc(:,ind)
  CALL orient_listn(self%lcfo(ind,j),fc,3_i4)
  f=0.d0
  f(fc(1))=fc_dofs(1,1)
  f(fc(2))=fc_dofs(2,1)
  f(fc(3))=fc_dofs(3,1)
  pt=self%log2phys(j,f)
  rtmp(:,i+self%np+2*self%ne)=pt
END DO
!---
nc = self%nc*27
ALLOCATE(lctmp(4,nc))
!$omp parallel do private(dofs,j)
DO i=1,self%nc
  dofs(1:4)=self%lc(:,i)
  DO j=1,6
    DO m=1,2
      IF(self%lce(j,i)>0)THEN
        dofs(j+(m-1)*6+4)=(ABS(self%lce(j,i))-1)*2+self%np+m
      ELSE
        dofs(j+(2-m)*6+4)=(ABS(self%lce(j,i))-1)*2+self%np+m
      END IF
    END DO
  END DO
  DO j=1,4
    dofs(j+16)=ABS(self%lcf(j,i))+self%np+2*self%ne
  END DO
  DO j=1,27
    lctmp(:,(i-1)*27+j)=dofs(tess3(:,j)+1)-1
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate3
!---------------------------------------------------------------------------------
!> Construct point and cell lists for a 3 level tessellation
!!
!! This corresponds to the input triangulation with cubic nodes added
!! - np_tess = np + 3*ne + 3*nf + nc
!! - nc_tess = nc*64
!---------------------------------------------------------------------------------
SUBROUTINE tessellate4(self,rtmp,lctmp)
CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,np_tess]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,nc_tess]
REAL(r8), PARAMETER :: ed_dofs(2,3)=RESHAPE((/1.d0/4,3.d0/4,1.d0/2,1.d0/2,3.d0/4,1.d0/4/),(/2,3/))
REAL(r8), PARAMETER :: fc_dofs(3,3)=RESHAPE((/1.d0/4,1.d0/4,1.d0/2,1.d0/4,1.d0/2,1.d0/4, &
1.d0/2,1.d0/4,1.d0/4/),(/3,3/))
REAL(r8), PARAMETER :: c_dofs(4,1)=RESHAPE((/1.d0/4,1.d0/4,1.d0/4,1.d0/4/),(/4,1/))
REAL(r8) :: f(4),pt(3)
INTEGER(i4) :: i,j,m,mm,ind,np,nc,ed(2),fc(3),dofs(35)
DEBUG_STACK_PUSH
!---
np = self%np+3*self%ne+3*self%nf+self%nc
ALLOCATE(rtmp(3,np))
!$omp parallel do
DO i=1,self%np
  rtmp(:,i)=self%r(:,i)
END DO
!$omp parallel do private(j,ind,ed,f,pt,m)
DO i=1,self%ne
  j=self%lec(self%kec(i))
  DO ind=1,6
    IF(ABS(self%lce(ind,j))==i)EXIT
  END DO
  ed=self%cell_ed(:,ind)
  CALL orient_list2(self%lce(ind,j),ed)
  f=0.d0
  DO m=1,3
    f(ed(1))=ed_dofs(1,m)
    f(ed(2))=ed_dofs(2,m)
    pt=self%log2phys(j,f)
    rtmp(:,(i-1)*3+m+self%np)=pt
  END DO
END DO
!$omp parallel do private(j,ind,fc,f,pt,m)
DO i=1,self%nf
  j=self%lfc(1,i)
  DO ind=1,4
    IF(ABS(self%lcf(ind,j))==i)EXIT
  END DO
  fc=self%cell_fc(:,ind)
  CALL orient_listn(self%lcfo(ind,j),fc,3_i4)
  f=0.d0
  DO m=1,3
    f(fc(1))=fc_dofs(1,m)
    f(fc(2))=fc_dofs(2,m)
    f(fc(3))=fc_dofs(3,m)
    pt=self%log2phys(j,f)
    rtmp(:,(i-1)*3+m+3*self%ne+self%np)=pt
  END DO
END DO
!$omp parallel do private(j,ind,fc,f,pt)
DO i=1,self%nc
  f=c_dofs(:,1)
  pt=self%log2phys(i,f)
  rtmp(:,i+3*self%nf+3*self%ne+self%np)=pt
END DO
!---
nc = self%nc*64
ALLOCATE(lctmp(4,nc))
!$omp parallel do private(dofs,j,m,mm,fc)
DO i=1,self%nc
  dofs(1:4)=self%lc(:,i)
  DO j=1,6
    DO m=1,3
      IF(self%lce(j,i)<0)THEN
        dofs(j+(3-m)*6+4)=(ABS(self%lce(j,i))-1)*3+self%np+m
      ELSE
        dofs(j+(m-1)*6+4)=(ABS(self%lce(j,i))-1)*3+self%np+m
      END IF
    END DO
  END DO
  DO j=1,4
    fc=(/1,2,3/)
    CALL orient_listn(self%lcfo(j,i),fc,3_i4)
    DO m=1,3
      DO mm=1,3
        IF(ALL(fc_dofs(fc,m)==fc_dofs((/1,2,3/),mm)))EXIT
      END DO
      dofs(j+(m-1)*4+3*6+4)=(ABS(self%lcf(j,i))-1)*3+self%np+3*self%ne+mm
    END DO
  END DO
  dofs(1+3*4+3*6+4)=i+self%np+3*self%ne+3*self%nf
  DO j=1,64
    lctmp(:,(i-1)*64+j)=dofs(tess4(:,j)+1)-1
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tessellate4
END MODULE tetmesh_tessellation
