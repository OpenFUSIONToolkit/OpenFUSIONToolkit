!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_tetmesh_type.F90
!
!> Triangular boundary mesh definitions
!! - MPI seam
!! - MPI I/O index
!! - MPI global context
!! - Base mesh linkage
!! - MPI preallocated Send/Recv
!! - High order tet representation
!! - Mesh container
!!
!! Global Tet variables
!! - Tetrahedra edge and face lists
!! - Grounding point position
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_trimesh_type
USE oft_base
USE oft_lag_poly
USE oft_tet_quadrature, ONLY: oft_quad_type, set_quad_2d
USE oft_mesh_type, ONLY: oft_bmesh
USE trimesh_tessellation, ONLY: tessellate1, tessellate2, tessellate3, tessellate4
IMPLICIT NONE
#include "local.h"
INTEGER(i4), PARAMETER :: tri_ed(2,3)=RESHAPE((/3,2,1,3,2,1/),(/2,3/)) !< Triangle edge list
!------------------------------------------------------------------------------
! TYPE oft_tetmesh
!------------------------------------------------------------------------------
!> Triangular boundary mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_bmesh) :: oft_trimesh
  REAL(r8), POINTER, DIMENSION(:) :: xnodes => NULL()
CONTAINS
  PROCEDURE :: setup => trimesh_setup
  PROCEDURE :: load_from_file => trimesh_load
  PROCEDURE :: save_to_file => trimesh_save
  PROCEDURE :: set_order => trimesh_set_order
  PROCEDURE :: invert_face => trimesh_invert_face
  PROCEDURE :: log2phys => trimesh_log2phys
  PROCEDURE :: phys2log => trimesh_phys2log
  PROCEDURE :: jacobian => trimesh_jacobian
  PROCEDURE :: hessian => trimesh_hessian
  PROCEDURE :: norm => trimesh_norm
  PROCEDURE :: tang => trimesh_tang
  PROCEDURE :: vlog => trimesh_vlog
  PROCEDURE :: in_cell => trimesh_in_cell
  PROCEDURE :: quad_rule => trimesh_quad_rule
  PROCEDURE :: tessellate => trimesh_tessellate
  PROCEDURE :: get_io_sizes => trimesh_get_io_sizes
END TYPE oft_trimesh
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_load
!---------------------------------------------------------------------------
!> Load trimesh from transfer file
!---------------------------------------------------------------------------
SUBROUTINE trimesh_load(self,filename)
CLASS(oft_trimesh), INTENT(inout) :: self
CHARACTER(LEN=*), INTENT(in) :: filename
INTEGER(i4) :: i,io_unit
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
!---Read in points
READ(io_unit,'(I8)')self%np
ALLOCATE(self%r(3,self%np))
self%r=0.d0
DO i=1,self%np
  READ(io_unit,'(3E25.16)')self%r(:,i)
END DO
IF(ALL(ABS(self%r(3,:))<1.d-10))self%dim=2
!---Read in faces
READ(io_unit,'(I8)')self%nc
ALLOCATE(self%lc(3,self%nc))
DO i=1,self%nc
  READ(io_unit,*)self%lc(:,i)
END DO
CLOSE(io_unit)
! self%stand_alone=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE trimesh_load
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_save
!---------------------------------------------------------------------------
!> Load trimesh from transfer file
!---------------------------------------------------------------------------
SUBROUTINE trimesh_save(self,filename)
CLASS(oft_trimesh), INTENT(in) :: self
CHARACTER(LEN=*), INTENT(in) :: filename
INTEGER(i4) :: i,io_unit
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
!---Write out points
WRITE(io_unit,'(I8)')self%np
DO i=1,self%np
  WRITE(io_unit,'(3E25.16)')self%r(:,i)
END DO
!---Write out points
WRITE(io_unit,'(I8)')self%nc
DO i=1,self%nc
  WRITE(io_unit,'(3I8)')self%lc(:,i)
END DO
CLOSE(io_unit)
DEBUG_STACK_POP
END SUBROUTINE trimesh_save
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_setup
!---------------------------------------------------------------------------
!> Load trimesh from transfer file
!---------------------------------------------------------------------------
SUBROUTINE trimesh_setup(self,cad_type,has_parent)
CLASS(oft_trimesh), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cad_type
LOGICAL, INTENT(in) :: has_parent
self%cell_np=3
self%cell_ne=self%cell_np
self%type=1
self%cad_type=cad_type
ALLOCATE(self%cell_ed(2,3))
self%cell_ed=tri_ed
CALL self%set_order(1)
IF(has_parent)ALLOCATE(self%parent)
END SUBROUTINE trimesh_setup
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_set_order
!---------------------------------------------------------------------------
!> Load trimesh from transfer file
!---------------------------------------------------------------------------
SUBROUTINE trimesh_set_order(self,order)
CLASS(oft_trimesh), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: order
INTEGER(i4) :: i,counts(3)
IF(order<1.OR.order>3)CALL oft_abort("Invalid grid order", "trimesh_set_order", __FILE__)
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
self%order=order
select case(self%order)
case(1)
  counts=0
case(2)
  self%ho_info%nep=1
  counts=(/self%ne,0,0/)
end select
IF(.NOT.ASSOCIATED(self%ho_info%r))THEN
  IF(SUM(counts)>0)ALLOCATE(self%ho_info%r(3,SUM(counts)))
  IF(self%ho_info%nep>0)THEN
    ALLOCATE(self%ho_info%lep(self%ho_info%nep,self%ne))
    self%ho_info%lep=-1
  END IF
  IF(self%ho_info%ncp>0)THEN
    ALLOCATE(self%ho_info%lcp(self%ho_info%ncp,self%nc))
    self%ho_info%lcp=-1
  END IF
END IF
ALLOCATE(self%xnodes(order+1))
DO i=1,order+1
  self%xnodes(i)=REAL(i-1,8)/REAL(order,8)
END DO
END SUBROUTINE trimesh_set_order
!---------------------------------------------------------------------------
! SUBROUTINE tri_2d_grid
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE tri_2d_grid(order,xnodes,inodesf)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: inodesf
REAL(r8), POINTER, DIMENSION(:), INTENT(out) :: xnodes
INTEGER(i4) :: i,j,np
!---Points
ALLOCATE(xnodes(order+1))
DO i=1,order+1
  xnodes(i)=REAL(i-1,8)/REAL(order,8)
END DO
!---Faces
IF(order>2)THEN
  np=0
  DO i=1,order-1
    DO j=1,order-1-i
      np=np+1
    END DO
  END DO
  ALLOCATE(inodesf(2,np))
  np=0
  DO i=1,order-1
    DO j=1,order-1-i
      np=np+1
      inodesf(:,np)=(/i,j/)
    END DO
  END DO
END IF
END SUBROUTINE tri_2d_grid
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_invert_face
!---------------------------------------------------------------------------
!> Invert cell to produce positive volume
!---------------------------------------------------------------------------
SUBROUTINE trimesh_invert_face(self,i)
CLASS(oft_trimesh), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: i
self%lc(2:3,i)=self%lc(3:2:-1,i)
IF(ASSOCIATED(self%lce))self%lce(:,i)=-self%lce([1,3,2],i)
IF(ASSOCIATED(self%lcc))self%lcc(:,i)=self%lcc([1,3,2],i)
END SUBROUTINE trimesh_invert_face
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_quad_rule
!------------------------------------------------------------------------------
!> Create quadrature rule for tetrahedra
!!
!! @param[in] self Mesh containing cell
!! @param[in] i Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @result Physical position [3]
!------------------------------------------------------------------------------
SUBROUTINE trimesh_quad_rule(self,order,quad_rule)
CLASS(oft_trimesh), INTENT(in) :: self
INTEGER(i4), INTENT(in) :: order
TYPE(oft_quad_type), INTENT(out) :: quad_rule
CALL set_quad_2d(quad_rule, order)
END SUBROUTINE trimesh_quad_rule
!------------------------------------------------------------------------------
!> Logical locations of vertices
!!
!! @param[in] self Mesh containing cell
!! @param[in] i Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @result Physical position [3]
!------------------------------------------------------------------------------
SUBROUTINE trimesh_vlog(self,i,f)
CLASS(oft_trimesh), INTENT(in) :: self
INTEGER(i4), INTENT(in) :: i
REAL(r8), INTENT(out) :: f(:)
f=0.d0
f(i)=1.d0
END SUBROUTINE trimesh_vlog
!------------------------------------------------------------------------------
!> Logical locations of vertices
!!
!! @param[in] self Mesh containing cell
!! @param[in] i Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @result Physical position [3]
!------------------------------------------------------------------------------
FUNCTION trimesh_in_cell(self,f,tol) RESULT(eedge)
CLASS(oft_trimesh), INTENT(in) :: self
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: tol
INTEGER(i4) :: eedge
REAL(r8) :: fmin,fmax
fmin=MINVAL(f(1:3))
fmax=MAXVAL(f(1:3))
IF((fmax<=1.d0+tol).AND.(fmin>=-tol))THEN
  eedge=0
ELSE
  eedge=MINLOC(f(1:3), DIM=1)
END IF
END FUNCTION trimesh_in_cell
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_tessellate
!------------------------------------------------------------------------------
!> Driver for order specific tessellation subroutines
!!
!! @note Cell lists are returned with zero based indexing.
!!
!! @warning The maximum tessellation order currently supported is 4.
!!
!! @param[in] self Mesh to tessellate
!! @param[out] rtmp Tessellated mesh points [3,np_tess]
!! @param[out] lctmp Tessellated face list [3,nf_tess]
!! @param[in] order Desired tessellation order
!------------------------------------------------------------------------------
SUBROUTINE trimesh_tessellate(self,rtmp,lctmp,order)
CLASS(oft_trimesh), INTENT(in) :: self
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp
INTEGER(i4), INTENT(in) :: order
DEBUG_STACK_PUSH
IF(self%skip)THEN
  rtmp=>NULL()
  lctmp=>NULL()
  DEBUG_STACK_POP
  RETURN
END IF
!---
IF(order==1)THEN
  CALL tessellate1(self,rtmp,lctmp)
ELSE IF(order==2)THEN
  CALL tessellate2(self,rtmp,lctmp)
ELSE IF(order==3)THEN
  CALL tessellate3(self,rtmp,lctmp)
ELSE IF(order==4)THEN
  CALL tessellate4(self,rtmp,lctmp)
ELSE
  CALL oft_abort('Invalid tessellation order','trimesh_tessellate',__FILE__)
END IF
DEBUG_STACK_POP
END SUBROUTINE trimesh_tessellate
!---------------------------------------------------------------------------
! FUNCTION: trimesh_get_io_sizes
!---------------------------------------------------------------------------
!> Get variable sizes following tessellation
!---------------------------------------------------------------------------
FUNCTION trimesh_get_io_sizes(self) result(sizes)
CLASS(oft_trimesh), INTENT(in) :: self
INTEGER(i4) :: sizes(2)
SELECT CASE(self%tess_order)
CASE(1)
  sizes=[self%np, self%nc]
CASE(2)
  sizes=[self%np+self%ne, self%nc*4]
CASE(3)
  sizes=[self%np+2*self%ne+self%nc, self%nc*9]
CASE(4)
  sizes=[self%np+3*self%ne+3*self%nc, self%nc*16]
CASE DEFAULT
  CALL oft_abort("Unknown tessellation size","trimesh_get_io_sizes",__FILE__)
END SELECT
END FUNCTION trimesh_get_io_sizes
!---------------------------------------------------------------------------
! SUBROUTINE trimesh_phys2log
!---------------------------------------------------------------------------
!> Map from physical to logical coordinates for a linear element
!!
!! @param[in] self Mesh containing cell
!! @param[in] cell Index of cell for evaulation
!! @param[in] pt Physical position [3]
!! @result f Logical coordinates within the cell [4]
!---------------------------------------------------------------------------
SUBROUTINE trimesh_phys2log(self,cell,pt,f)
class(oft_trimesh), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: pt(3)
real(r8), intent(out) :: f(:)
real(r8) :: gop(3,3),v
integer(i4) :: k
call trimesh_jacl(self,cell,f,gop,v)
do k=1,3
  f(k)=1.d0+dot_product(pt-self%r(:,self%lc(k,cell)),gop(:,k))
end do
end SUBROUTINE trimesh_phys2log
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_log2phys
!------------------------------------------------------------------------------
!> Map from logical to physical coordinates
!!
!! Driver function calls mapping specific function depending on mesh order.
!!
!! @param[in] self Mesh containing cell
!! @param[in] cell Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @result pt Physical position [3]
!------------------------------------------------------------------------------
function trimesh_log2phys(self,cell,f) result(pt)
class(oft_trimesh), intent(in) :: self
integer, intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8) :: pt(3)
integer(i4) :: i,j,k,dof,ed
DEBUG_STACK_PUSH
if(self%order>3)call oft_abort('Invalid mesh order','tetmesh_log2phys',__FILE__)
pt=0.d0
IF(self%order>1)THEN
  do i=1,3
    pt = pt + lag_1d(self%order+1,f(i),self%xnodes,self%order+1) &
          *self%r(:,self%lc(i,cell))
  end do
  DO i=1,3
    ed=self%lce(i,cell)
    j=abs(ed)
    DO k=1,self%ho_info%nep
      dof=k
      IF(ed<0)dof=self%ho_info%nep+1-k
      pt = pt + lag_1d_bary(k,f(tri_ed(:,i)),self%xnodes,self%order+1) &
            *self%ho_info%r(:,self%ho_info%lep(dof,j))
    END DO
  END DO
  IF(self%ho_info%ncp>0)THEN
    pt = pt + lag_2d_bary([1,1],f,self%xnodes,self%order+1) &
          *self%ho_info%r(:,self%ho_info%lcp(1,cell))
  END IF
ELSE
  do i=1,3
    pt = pt + f(i)*self%r(:,self%lc(i,cell))
  end do
END IF
DEBUG_STACK_POP
end function trimesh_log2phys
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_jacobian
!------------------------------------------------------------------------------
!> Compute the jacobian matrix and its determinant for a grid cell
!!
!! Driver function calls mapping specific function depending on mesh order.
!!
!! @param[in] self Mesh containing cell
!! @param[in] cell Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @param[out] gop Jacobian matrix ( \f$ \frac{\partial x_i}{\partial \lambda_j} \f$ ) [3,4] (optional)
!! @param[out] j Jacobian of transformation from logical to physical coordinates (optional)
!------------------------------------------------------------------------------
subroutine trimesh_jacobian(self,cell,f,gop,j)
class(oft_trimesh), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: gop(:,:)
real(r8), intent(out) :: j
real(r8) :: jfull(2,3),pt(3),pt_proj(2),getmp(2),gftmp(3),A(2,2),C(2,2),t(3,2)
integer(i4) :: k,l,ed,etmp(2),dof
DEBUG_STACK_PUSH
if(self%order>3)call oft_abort('Invalid mesh order','trimesh_jacobian',__FILE__)
IF(self%dim==2)THEN
  jfull=0.d0
  IF(self%order>1)THEN
    do k=1,3 ! Get corner nodes
      pt=self%r(:,self%lc(k,cell))
      jfull(:,k)=jfull(:,k) &
        + dlag_1d(self%order+1,f(k),self%xnodes,self%order+1)*pt(1:2)
    end do
    do k=1,3 ! Get edge nodes
      ed=self%lce(k,cell)
      do dof=1,self%ho_info%nep
        IF(ed<0)THEN
          pt=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
        ELSE
          pt=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
        END IF
        getmp=dlag_1d_bary(dof,f(tri_ed(:,k)),self%xnodes,self%order+1)
        do l=1,2
          jfull(:,tri_ed(l,k))=jfull(:,tri_ed(l,k))+getmp(l)*pt(1:2)
        end do
      end do
    end do
    ! IF(self%ho_info%ncp==1)THEN
    !   ! Add loop and orient if add higher order
    !   pt=self%ho_info%r(:,self%ho_info%lcp(1,cell))
    !   gftmp=dlag_2d_bary((/1,1/),f,self%xnodes,self%order+1)
    !   do l=1,3
    !     jfull(:,l)=jfull(:,l)+gftmp(l)*pt(1:2)
    !   end do
    ! END IF
  ELSE
    ! Get node points
    do k=1,3
      jfull(:,k)=self%r(1:2,self%lc(k,cell))
    end do
  END IF
ELSE
  CALL trimesh_tang(self,cell,f,t)
  jfull=0.d0
  IF(self%order>1)THEN
    do k=1,3 ! Get corner nodes
      pt=self%r(:,self%lc(k,cell))
      pt_proj=MATMUL(pt,t)
      jfull(:,k)=jfull(:,k) + dlag_1d(self%order+1,f(k),self%xnodes,self%order+1) &
                            *pt_proj
    end do
    do k=1,3 ! Get edge nodes
      ed=self%lce(k,cell)
      do dof=1,self%ho_info%nep
        IF(ed<0)THEN
          pt=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
        ELSE
          pt=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
        END IF
        pt_proj=MATMUL(pt,t)
        getmp=dlag_1d_bary(dof,f(tri_ed(:,k)),self%xnodes,self%order+1)
        do l=1,2
          jfull(:,tri_ed(l,k))=jfull(:,tri_ed(l,k))+getmp(l)*pt_proj
        end do
      end do
    end do
    ! IF(self%ho_info%ncp==1)THEN
    !   ! Add loop and orient if add higher order
    !   pt=self%ho_info%r(:,self%ho_info%lcp(1,cell))
    !   pt_proj=MATMUL(pt,t)
    !   gftmp=dlag_2d_bary((/1,1/),f,self%xnodes,self%order+1)
    !   do l=1,3
    !     jfull(:,l)=jfull(:,l)+gftmp(l)*pt_proj
    !   end do
    ! END IF
  ELSE
    ! Get node points
    do k=1,3
      pt=self%r(:,self%lc(k,cell))
      pt_proj=MATMUL(pt,t)
      jfull(:,k)=pt_proj
    end do
  END IF
END IF
!---Condense Jacobian
do k=1,2
  A(k,:)=jfull(:,k+1)-jfull(:,1)
end do
!---Comute inverse
call trimesh_jacinv(A,C,j)
!---Expand gradient matrix
IF(self%dim==2)THEN
  gop(:,2)=[C(1,1),C(2,1),0.d0]
  gop(:,3)=[C(1,2),C(2,2),0.d0]
ELSE
  gop(:,2)=C(1,1)*t(:,1)+C(2,1)*t(:,2)
  gop(:,3)=C(1,2)*t(:,1)+C(2,2)*t(:,2)
END IF
gop(:,1)=-(gop(:,2)+gop(:,3))
j=j/2.d0
DEBUG_STACK_POP
end subroutine trimesh_jacobian
!------------------------------------------------------------------------------
! SUBROUTINE trimesh_hessian
!------------------------------------------------------------------------------
!> Compute the second order jacobians for a grid cell
!!
!! @param[in] self Mesh containing cell
!! @param[in] cell Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @param[out] g2op Second order Jacobian matrix
!! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
!! [6,6]
!! @param[out] K Gradient correction matrix
!! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [6,3]
!------------------------------------------------------------------------------
subroutine trimesh_hessian(self,cell,f,g2op,K)
class(oft_trimesh), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: g2op(:,:)
real(r8), intent(out) :: K(:,:)
real(r8) :: jfull(2,3),getmp(2),gftmp(3),d2etmp(3),d2ftmp(6),pt(3)
real(r8) :: ri(2,3),jac(2,2),A(3,3),C(3,3),v
integer(i4) :: j,m,l,etmp(2),dof,ed
integer(i4), parameter :: pmap(4)=(/1,5,8,10/)
integer(i4), parameter :: emap(6)=(/4,7,9,6,3,2/)
integer(i4), parameter :: fmap(4,4)=RESHAPE((/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/),(/4,4/))
if(self%order>3)call oft_abort('Invalid mesh order','trimesh_hessian',__FILE__)
IF(self%dim==3)call oft_abort('Only supported for dim=2','trimesh_hessian',__FILE__)
jfull=0.d0
K=0.d0
! IF(cell_is_curved(self, i))THEN
!   do m=1,4 ! Get corner nodes
!     jfull(:,m)=jfull(:,m) + dlag_1d(self%order+1,f(m),self%xnodes,self%order+1)*self%r(:,self%lc(m,i))
!     !
!     K(pmap(m),:) = K(pmap(m),:) + d2lag_1d(self%order+1,f(m),self%xnodes,self%order+1)*self%r(:,self%lc(m,i))
!   end do
!   do m=1,6 ! Get edge nodes
!     ed=self%lce(m,i)
!     do dof=1,self%ho_info%nep
!       IF(ed<0)THEN
!         pt=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
!       ELSE
!         pt=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
!       END IF
!       getmp=dlag_1d_bary(dof,f(tet_ed(:,m)),self%xnodes,self%order+1)
!       do l=1,2
!         jfull(:,tet_ed(l,m))=jfull(:,tet_ed(l,m))+getmp(l)*pt
!       end do
!       !
!       d2etmp=d2lag_1d_bary(dof,f(tet_ed(:,m)),self%xnodes,self%order+1)
!       K(pmap(tet_ed(1,m)),:)=K(pmap(tet_ed(1,m)),:)+d2etmp(1)*pt
!       K(emap(m),:)=K(emap(m),:)+d2etmp(2)*pt
!       K(pmap(tet_ed(2,m)),:)=K(pmap(tet_ed(2,m)),:)+d2etmp(3)*pt
!     end do
!   end do
!   IF(self%ho_info%ncp>0)THEN
!     do m=1,4 ! Get face nodes
!       pt=self%ho_info%r(:,self%ho_info%lcp(1,ABS(self%lcf(m,i))))
!       gftmp=dlag_2d_bary((/1,1/),f(tet_fc(:,m)),self%xnodes,self%order+1)
!       do l=1,3
!         jfull(:,tet_fc(l,m))=jfull(:,tet_fc(l,m))+gftmp(l)*pt
!       end do
!       !
!       d2ftmp=d2lag_2d_bary((/1,1/),f(tet_fc(:,m)),self%xnodes,self%order+1)
!       K(pmap(tet_fc(1,m)),:)=K(pmap(tet_fc(1,m)),:)+d2ftmp(1)*pt
!       K(fmap(tet_fc(1,m),tet_fc(2,m)),:)=K(fmap(tet_fc(1,m),tet_fc(2,m)),:)+d2ftmp(2)*pt
!       K(fmap(tet_fc(1,m),tet_fc(3,m)),:)=K(fmap(tet_fc(1,m),tet_fc(3,m)),:)+d2ftmp(3)*pt
!       K(pmap(tet_fc(2,m)),:)=K(pmap(tet_fc(2,m)),:)+d2ftmp(4)*pt
!       K(fmap(tet_fc(2,m),tet_fc(3,m)),:)=K(fmap(tet_fc(2,m),tet_fc(3,m)),:)+d2ftmp(5)*pt
!       K(pmap(tet_fc(3,m)),:)=K(pmap(tet_fc(3,m)),:)+d2ftmp(6)*pt
!     end do
!   END IF
! ELSE
  do j=1,3 ! Get node points
    ri(:,j)=self%r(1:2,self%lc(j,cell))
  end do
  !---Get Jacobian entries
  jfull=ri
  do j=1,2
    jac(j,:)=jfull(:,j+1)-jfull(:,1)
  end do
  !---Get 2nd order mapping
  !---Row 1
  A(1,1)=jac(1,1)**2
  A(1,2)=2*jac(1,1)*jac(1,2)
  A(1,3)=jac(1,2)**2
  !---Row 2
  A(2,1)=jac(1,1)*jac(2,1)
  A(2,2)=jac(1,1)*jac(2,2)+jac(2,1)*jac(1,2)
  A(2,3)=jac(1,2)*jac(2,2)
  !---Row 3
  A(3,1)=jac(2,1)**2
  A(3,2)=2*jac(2,1)*jac(2,2)
  A(3,3)=jac(2,2)**2
  !---Invert
  call trimesh_m3inv(A,C,v)
  !---Map back to g2op
  g2op=0.d0
  g2op([1,2,4],1)=C(:,1)
  g2op([1,2,4],2)=C(:,2)
  g2op([1,2,4],4)=C(:,3)
  !---Scatter to redundant dimension
  g2op([1,2,4],3)=-2*C(:,1)-C(:,2)
  g2op([1,2,4],5)=-C(:,2)-2*C(:,3)
  g2op([1,2,4],6)=C(:,1)+C(:,2)+C(:,3)
! END IF
end subroutine trimesh_hessian
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_jacl
!------------------------------------------------------------------------------
!> Compute the jacobian matrix and its determinant for a linear element
!!
!! @param[in] self Mesh containing cell
!! @param[in] i Index of cell for evaulation
!! @param[in] f Logical coordinate in cell [4]
!! @param[out] gop Jacobian matrix ( \f$ \frac{\partial x_i}{\partial \lambda_j} \f$ ) [3,4] (optional)
!! @param[out] j Jacobian of transformation from logical to physical coordinates (optional)
!------------------------------------------------------------------------------
subroutine trimesh_jacl(self,i,f,gop,j)
class(oft_trimesh), intent(in) :: self
integer(i4), intent(in) :: i
real(r8), intent(in) :: f(3)
real(r8), intent(out) :: gop(3,3)
real(r8), intent(out) :: j
real(r8) :: ri(3,3),jfull(2,3),A(2,2),C(2,2),t(3,2)
integer(i4) :: k
DEBUG_STACK_PUSH
IF(self%dim==2)THEN
  !---Get Jacobian entries
  DO k=1,3
    jfull(:,k)=self%r(1:2,self%lc(k,i))
  END DO
ELSE
  !---Get node points
  do k=1,3
    ri(:,k)=self%r(:,self%lc(k,i))
  end do
  CALL trimesh_tang(self,i,f,t)
  !---Get Jacobian entries
  DO k=1,3
    jfull(1,k)=DOT_PRODUCT(ri(:,k),t(:,1))
    jfull(2,k)=DOT_PRODUCT(ri(:,k),t(:,2))
  END DO
END IF
do k=1,2
  A(k,:)=jfull(:,k+1)-jfull(:,1)
end do
call trimesh_jacinv(A,C,j)
IF(self%dim==2)THEN
  gop(:,2)=[C(1,1),C(2,1),0.d0]
  gop(:,3)=[C(1,2),C(2,2),0.d0]
ELSE
  gop(:,2)=C(1,1)*t(:,1)+C(2,1)*t(:,2)
  gop(:,3)=C(1,2)*t(:,1)+C(2,2)*t(:,2)
END IF
gop(:,1)=-(gop(:,2)+gop(:,3))
j=j/2.d0
DEBUG_STACK_POP
end subroutine trimesh_jacl
!------------------------------------------------------------------------------
! FUNCTION: trimesh_norm
!------------------------------------------------------------------------------
!> Compute the unit normal vector to a face
!!
!! @param[in] self Mesh containing face
!! @param[in] i Index of face for evaulation
!! @param[in] f Logical coordinate on face [3]
!! @param[out] n Normal vector [3]
!------------------------------------------------------------------------------
subroutine trimesh_norm(self,i,f,n)
class(oft_trimesh), target, intent(in) :: self
integer(i4), intent(in) :: i
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: n(3)
real(r8) :: t(3,2)
DEBUG_STACK_PUSH
CALL trimesh_tang(self,i,f,t)
n = cross_product(t(:,1),t(:,2))
DEBUG_STACK_POP
end subroutine trimesh_norm
!------------------------------------------------------------------------------
! FUNCTION: trimesh_tang
!------------------------------------------------------------------------------
!> Compute a orthonormal set of axis tangential to a face
!!
!! @param[in] self Mesh containing face
!! @param[in] i Index of face for evaulation
!! @param[in] f Logical coordinate on face [3]
!! @param[out] t Orthonormal vectors [3,2]
!------------------------------------------------------------------------------
subroutine trimesh_tang(self,i,f,t)
class(oft_trimesh), target, intent(in) :: self
integer(i4), intent(in) :: i
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: t(3,2)
real(r8) :: pt(3)
DEBUG_STACK_PUSH
pt = trimesh_glogphys(self,i,1,f)
!---Get tangent direction 1
t(:,1) = trimesh_glogphys(self,i,2,f) - pt
t(:,1) = t(:,1)/magnitude(t(:,1))
!---Get tangent direction 2
t(:,2) = trimesh_glogphys(self,i,3,f) - pt
t(:,2) = t(:,2) - DOT_PRODUCT(t(:,2),t(:,1))*t(:,1)
t(:,2) = t(:,2)/magnitude(t(:,2))
DEBUG_STACK_POP
end subroutine trimesh_tang
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_glogphys
!------------------------------------------------------------------------------
!> Compute the partial derivative of the physical coordinates with a specific logical coordinate
!!
!! Driver function calls mapping specific function depending on mesh order.
!!
!! @param[in] self Mesh containing cell
!! @param[in] i Index of cell for evaulation
!! @param[in] k Logical coordinate for differentiation
!! @param[in] f Logical coordinate in cell [4]
!! @result pt \f$ \frac{\partial r}{\partial f_k} \f$ [3]
!------------------------------------------------------------------------------
function trimesh_glogphys(self,face,j,f) result(pt)
class(oft_trimesh), intent(in) :: self
integer, intent(in) :: face
integer, intent(in) :: j
real(r8), intent(in) :: f(3)
real(r8) :: pt(3)
real(r8) :: getmp(2),gftmp(3),pttmp(3)
integer(i4) :: k,l,ed,etmp(2),dof
DEBUG_STACK_PUSH
if(self%order>3)call oft_abort('Invalid mesh order','trimesh_jacobian',__FILE__)
IF(self%order>1)THEN
  pt=dlag_1d(self%order+1,f(j),self%xnodes,self%order+1)*self%r(:,self%lc(j,face))
  do k=1,3 ! Get edge nodes
    ed=self%lce(k,face)
    do dof=1,self%ho_info%nep
      IF(ed<0)THEN
        pttmp=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
      ELSE
        pttmp=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
      END IF
      getmp=dlag_1d_bary(dof,f(tri_ed(:,k)),self%xnodes,self%order+1)
      do l=1,2
        IF(tri_ed(l,k)==j)pt=pt+getmp(l)*pttmp
      end do
    end do
  end do
  IF(self%ho_info%ncp>0)THEN
    ! Add loop and orient if add higher order
    pttmp=self%ho_info%r(:,self%ho_info%lcp(1,face))
    gftmp=dlag_2d_bary((/1,1/),f,self%xnodes,self%order+1)
    pt=pt+gftmp(j)*pttmp
  END IF
ELSE
  pt=self%r(:,self%lc(j,face))
END IF
DEBUG_STACK_POP
end function trimesh_glogphys
!------------------------------------------------------------------------------
! SUBROUTINE: trimesh_jacinv
!------------------------------------------------------------------------------
!> Invert a 2x2 matrix.
!!
!! @param[in] A Matrix to invert
!! @param[out] C \f$ A^{-1} \f$
!! @param[out] j |A|
!------------------------------------------------------------------------------
subroutine trimesh_jacinv(A,C,j)
real(r8), intent(in) :: A(2,2)
real(r8), intent(out) :: C(2,2)
real(r8), intent(out) :: j
DEBUG_STACK_PUSH
! Compute Det(A)
j=A(1,1)*A(2,2)-A(1,2)*A(2,1)
! Compute cofactor matrix
! Diagonal terms
C(1,1)=A(2,2)
C(2,2)=A(1,1)
! Odd terms
C(1,2)=-A(1,2)
C(2,1)=-A(2,1)
C=C/j
DEBUG_STACK_POP
end subroutine trimesh_jacinv
!------------------------------------------------------------------------------
! SUBROUTINE trimesh_m3inv
!------------------------------------------------------------------------------
!> Invert a 3x3 matrix.
!!
!! @param[in] A Matrix to invert
!! @param[out] C \f$ A^{-1} \f$
!! @param[out] j |A|
!------------------------------------------------------------------------------
subroutine trimesh_m3inv(A,C,j)
real(8), intent(in) :: A(3,3)
real(8), intent(out) :: C(3,3)
real(8), intent(out) :: j
real(8) :: t1,t2,t3
!---Compute resusables
t1=A(2,2)*A(3,3)-A(3,2)*A(2,3)
t2=A(3,2)*A(1,3)-A(1,2)*A(3,3)
t3=A(1,2)*A(2,3)-A(2,2)*A(1,3)
!---Compute Det(A)
j=A(1,1)*t1 + A(2,1)*t2 + A(3,1)*t3
!---Compute cofactor matrix
!---Diagonal terms
C(1,1)=t1
C(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
C(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
!---Odd terms
C(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
C(1,2)=t2
C(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
C(2,3)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
!---Even terms
C(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
C(1,3)=t3
!---Scale
C=C/j
end subroutine trimesh_m3inv
END MODULE oft_trimesh_type
