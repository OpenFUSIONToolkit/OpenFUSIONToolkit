!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_quadmesh_type.F90
!
!> Quadralateral boundary mesh definitions
!!
!! @author Chris Hansen
!! @date October 2018
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_quadmesh_type
USE oft_base
USE oft_gauss_quadrature, ONLY: oft_quad_type, set_quad_2d
USE oft_lag_poly
USE oft_mesh_type, ONLY: oft_bmesh
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> Quadralateral surface mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!---------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_bmesh) :: oft_quadmesh
  INTEGER(i4) :: inodesp(2,4)
  INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodese => NULL()
  INTEGER(i4), POINTER, DIMENSION(:,:) :: inodesf => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: xnodes => NULL()
CONTAINS
  PROCEDURE :: setup => quadmesh_setup
  PROCEDURE :: load_from_file => quadmesh_load
  PROCEDURE :: save_to_file => quadmesh_save
  PROCEDURE :: set_order => quadmesh_set_order
  PROCEDURE :: invert_face => quadmesh_invert_face
  PROCEDURE :: log2phys => quadmesh_log2phys
  PROCEDURE :: phys2log => quadmesh_phys2log
  PROCEDURE :: jacobian => quadmesh_jacobian
  PROCEDURE :: hessian => quadmesh_hessian
  PROCEDURE :: norm => quadmesh_norm
  PROCEDURE :: tang => quadmesh_tang
  PROCEDURE :: vlog => quadmesh_vlog
  PROCEDURE :: in_cell => quadmesh_in_cell
  PROCEDURE :: quad_rule => quadmesh_quad_rule
  PROCEDURE :: tessellate => quadmesh_tessellate
  PROCEDURE :: tessellated_sizes => quadmesh_tessellated_sizes
END TYPE oft_quadmesh
INTEGER(i4), PARAMETER :: quad_ed(2,4)=RESHAPE((/1,2, 2,3, 3,4, 4,1/),(/2,4/)) !< Quad edge list
INTEGER(i4), PARAMETER :: quad_bary_map(4) = [-2,1,2,-1]
INTEGER(i4), PARAMETER :: inodesp_base(2,4) = RESHAPE((/ &
  0,0, 1,0, 1,1, 0,1/),(/2,4/))
INTEGER(i4), PARAMETER :: inodes1p(2,4) = RESHAPE((/ &
  1,1, 2,1, 2,2, 1,2/),(/2,4/))
INTEGER(i4), PARAMETER :: inodes2p(2,4) = RESHAPE((/ &
  1,1, 3,1, 3,3, 1,3/),(/2,4/))
INTEGER(i4), PARAMETER :: inodes2e(2,4) = RESHAPE((/ &
  2,1, 3,2, 2,3, 1,2/),(/2,4/))
INTEGER(i4), PARAMETER :: inodes2f(2) = (/2,2/)
!
PUBLIC quad_2d_grid, quad_grid_orient
CONTAINS
!------------------------------------------------------------------------------
!> Load trimesh from transfer file
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_load(self,filename)
CLASS(oft_quadmesh), INTENT(inout) :: self !< Mesh object
CHARACTER(LEN=*), INTENT(in) :: filename !< File to load mesh from
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
ALLOCATE(self%lc(4,self%nc))
DO i=1,self%nc
  READ(io_unit,'(4I8)')self%lc(:,i)
END DO
CLOSE(io_unit)
! self%stand_alone=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE quadmesh_load
!------------------------------------------------------------------------------
!> Save trimesh from transfer file
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_save(self,filename)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
CHARACTER(LEN=*), INTENT(in) :: filename !< File to save mesh to
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
  WRITE(io_unit,'(4I8)')self%lc(:,i)
END DO
CLOSE(io_unit)
DEBUG_STACK_POP
END SUBROUTINE quadmesh_save
!------------------------------------------------------------------------------
!> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_setup(self,cad_type,has_parent)
CLASS(oft_quadmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
LOGICAL, INTENT(in) :: has_parent !< Is this mesh the/a surface of a volume mesh?
self%cell_np=4
self%cell_ne=self%cell_np
self%type=3
self%cad_type=cad_type
ALLOCATE(self%cell_ed(2,4))
self%cell_ed=quad_ed
CALL self%set_order(1)
IF(has_parent)ALLOCATE(self%parent)
END SUBROUTINE quadmesh_setup
!------------------------------------------------------------------------------
!> Set maximum order of spatial mapping
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_set_order(self,order)
CLASS(oft_quadmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
INTEGER(i4) :: counts(2)
IF(order<1.OR.order>2)CALL oft_abort("Invalid grid order", "quadmesh_setup", __FILE__)
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
IF(ASSOCIATED(self%inodese))DEALLOCATE(self%inodese,self%inodesf)
self%order=order
!---Setup high order gemetry
select case(self%order)
case(1)
  counts=0
case(2)
  self%ho_info%nep=1
  self%ho_info%ncp=1
  counts=(/self%ne,self%nc/)
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
!---Points
CALL quad_2d_grid(order,self%xnodes,self%inodesp,self%inodese,self%inodesf)
END SUBROUTINE quadmesh_set_order
!------------------------------------------------------------------------------
!> Map from orthogonal logical coordinates to barycentric logical coordinates
!------------------------------------------------------------------------------
PURE FUNCTION quad_get_bary(flog) RESULT(fbary)
REAL(r8), INTENT(in) :: flog(:) !< Position in orthogonal logical coordinates [2]
REAL(r8) :: fbary(4) !< Position in barycentric logical coordinates [4]
INTEGER(i4) :: i
DO i=1,4
  IF(quad_bary_map(i)<0)THEN
    fbary(i)=1.d0-flog(ABS(quad_bary_map(i)))
  ELSE
    fbary(i)=flog(ABS(quad_bary_map(i)))
  END IF
END DO
END FUNCTION quad_get_bary
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE quad_2d_grid(order,xnodes,inodesp,inodese,inodesf)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4), INTENT(out) :: inodesp(2,4)
INTEGER(i4), POINTER, DIMENSION(:,:,:), INTENT(out) :: inodese
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: inodesf
REAL(r8), POINTER, DIMENSION(:), INTENT(out) :: xnodes
INTEGER(i4) :: i,j,np,inc1,ind1,inode(2),counts(3)
!---Points
inodesp = inodesp_base*order + 1
ALLOCATE(xnodes(order+1))
DO i=1,order+1
  xnodes(i)=REAL(i-1,8)/REAL(order,8)
END DO
!---Edges
np = order - 1
IF(np>0)THEN
  ALLOCATE(inodese(2,np,4))
  DO i=1,4
    ind1 = MAXLOC(ABS(inodesp(:,quad_ed(2,i))-inodesp(:,quad_ed(1,i))), DIM=1)
    inc1 = SIGN(1_i4, inodesp(ind1,quad_ed(2,i))-inodesp(ind1,quad_ed(1,i)))
    inode = inodesp(:,quad_ed(1,i))
    DO j=1,order-1
      inode(ind1) = inode(ind1) + inc1
      inodese(:,j,i) = inode
      ! WRITE(*,*)i,j,inode
    END DO
  END DO
  !---Cell
  ALLOCATE(inodesf(2,np*np))
  ! WRITE(*,*)
  DO i=1,np
    DO j=1,np
      inodesf(:,(i-1)*np+j) = (/1+i,1+j/)
    END DO
  END DO
END IF
END SUBROUTINE quad_2d_grid
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE quad_grid_orient(oflag,order,finds)
INTEGER(i4), INTENT(in) :: oflag,order
INTEGER(i4), DIMENSION(:), INTENT(inout) :: finds
INTEGER(i4), PARAMETER :: inodesp(2,4) = RESHAPE((/0,0, 1,0, 1,1, 0,1/), (/2,4/))
INTEGER(i4) :: i,j,inc1(2),inc2(2),np,ptmp(4),inode(2),p1(2),p2(2),p3(2)
REAL(r8) :: o8
!
ptmp=(/1,2,3,4/)
CALL orient_listn(oflag, ptmp, 4_i4)
p1=inodesp(:,ptmp(1))*order+1
p2=inodesp(:,ptmp(2))*order+1
p3=inodesp(:,ptmp(3))*order+1
np = order-1
! Directions
inc1 = (p2-p1)
inc2 = (p3-p2)
!
DO i=1,np
  DO j=1,np
    inode = p1 + inc1*i/(order) + inc2*j/(order) - 1
    finds((i-1)*np+j) = (inode(1)-1)*np + inode(2)
  END DO
END DO
END SUBROUTINE quad_grid_orient
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE quad_grid_orient_inv(oflag,order,finds)
INTEGER(i4), INTENT(in) :: oflag,order
INTEGER(i4), DIMENSION(:), INTENT(inout) :: finds
INTEGER(i4), PARAMETER :: inodesp(2,4) = RESHAPE((/0,0, 1,0, 1,1, 0,1/), (/2,4/))
INTEGER(i4) :: i,j,inc1(2),inc2(2),np,ptmp(4),inode(2),p1(2),p2(2),p3(2)
REAL(r8) :: o8
!
ptmp=(/1,2,3,4/)
CALL orient_listn_inv(oflag, ptmp, 4_i4)
p1=inodesp(:,ptmp(1))*order+1
p2=inodesp(:,ptmp(2))*order+1
p3=inodesp(:,ptmp(3))*order+1
np = order-1
! Directions
inc1 = (p2-p1)
inc2 = (p3-p2)
!
DO i=1,np
  DO j=1,np
    inode = p1 + inc1*i/(order) + inc2*j/(order) - 1
    finds((i-1)*np+j) = (inode(1)-1)*np + inode(2)
  END DO
END DO
END SUBROUTINE quad_grid_orient_inv
!------------------------------------------------------------------------------
!> Turn cell "inside out", used to ensure consistent orientations
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_invert_face(self,cell)
CLASS(oft_quadmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Maximum order of spatial mapping
self%lc((/2,4/),cell)=self%lc((/4,2/),cell)
IF(ASSOCIATED(self%lce))self%lce(:,cell)=-self%lce([4,3,2,1],cell)
IF(ASSOCIATED(self%lcc))self%lcc(:,cell)=self%lcc([4,3,2,1],cell)
END SUBROUTINE quadmesh_invert_face
!---------------------------------------------------------------------------------
!> Retrieve suitable quadrature rule for triangular mesh with given order
!---------------------------------------------------------------------------------
SUBROUTINE quadmesh_quad_rule(self,order,quad_rule)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
CALL set_quad_2d(quad_rule, order)
END SUBROUTINE quadmesh_quad_rule
!---------------------------------------------------------------------------------
!> Get position in logical space of cell vertex `i`
!---------------------------------------------------------------------------------
SUBROUTINE quadmesh_vlog(self,i,f)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: i !< Vertex to locate
REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
f(1:2)=self%xnodes(self%inodesp(:,i))
END SUBROUTINE quadmesh_vlog
!---------------------------------------------------------------------------------
!> Test if logical position lies within the base cell
!!
!! @returns Position `f` is inside the base cell?
!---------------------------------------------------------------------------------
FUNCTION quadmesh_in_cell(self,f,tol) RESULT(eedge)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
REAL(r8), INTENT(in) :: tol !< Tolerance for test
INTEGER(i4) :: eedge
REAL(r8) :: fmin,fmax,fbary(4)
fmin=MINVAL(f(1:2))
fmax=MAXVAL(f(1:2))
IF((fmax<=1.d0+tol).AND.(fmin>=-tol))THEN
  eedge=0
ELSE
  fbary=quad_get_bary(f(1:2))
  eedge=MINLOC(fbary, DIM=1)
END IF
END FUNCTION quadmesh_in_cell
!---------------------------------------------------------------------------------
!> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
!!
!! @note The maximum tessellation order currently supported is 4
!! (may be lower for certain mesh types).
!!
!! @warning Cell lists are returned with zero based indexing
!---------------------------------------------------------------------------------
SUBROUTINE quadmesh_tessellate(self,rtmp,lctmp,order)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,:]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [self%ncp,:]
INTEGER(i4), INTENT(in) :: order !< Tessellation order
INTEGER(i4) :: i,j,k,cell,np,inodesp(2,4),i1,i2,inc
INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodese
INTEGER(i4), POINTER, DIMENSION(:,:) :: inodesf,pmap
REAL(r8), POINTER, DIMENSION(:) :: xnodes
LOGICAL, ALLOCATABLE, DIMENSION(:) :: pflag
DEBUG_STACK_PUSH
!---Call desired driver
IF(order==1)THEN
  ALLOCATE(rtmp(3,self%np))
  rtmp=self%r
  ALLOCATE(lctmp(4,self%nc))
  lctmp=self%lc-1
ELSE
  CALL quad_2d_grid(order,xnodes,inodesp,inodese,inodesf)
  np=order-1
  ALLOCATE(rtmp(3,self%np + self%ne*np + self%nc*np**2))
  ALLOCATE(pflag( self%np + self%ne*np + self%nc*np**2))
  ALLOCATE(pmap(order+1,order+1))
  ALLOCATE(lctmp(4, self%nc*order**2))
  lctmp=-1
  pflag=.FALSE.
  rtmp(:,1:self%np)=self%r
  pflag(1:self%np)=.TRUE.
  DO cell=1,self%nc
    pmap=-1
    !---Point
    DO i=1,4
      j=self%lc(i,cell)
      pmap(inodesp(1,i),inodesp(2,i))=j
    END DO
    !---Edges
    DO i=1,4
      DO k=1,np
        i1=k
        IF(self%lce(i,cell)<0)i1=np+1-k
        j=(ABS(self%lce(i,cell))-1)*np + self%np + i1
        pmap(inodese(1,k,i),inodese(2,k,i))=j
        IF(.NOT.pflag(j))THEN
          rtmp(:,j)=self%log2phys(cell, xnodes(inodese(:,k,i)))
          pflag(j)=.TRUE.
        END IF
      END DO
    END DO
    !---Face interiors
    DO k=1,np**2
      j=(cell-1)*np**2 + self%ne*np + self%np + k
      pmap(inodesf(1,k),inodesf(2,k))=j
      IF(.NOT.pflag(j))THEN
        rtmp(:,j)=self%log2phys(cell, xnodes(inodesf(:,k)))
        pflag(j)=.TRUE.
      END IF
    END DO
    !---Create faces
    IF(ANY(pmap<0))THEN
      DO i=1,order+1
        DO j=1,order+1
          WRITE(*,*)i,j,pmap(i,j)
        END DO
      END DO
      CALL oft_abort("Bad pmap","",__FILE__)
    END IF
    DO i=1,order
      DO j=1,order
        lctmp(:,(cell-1)*order**2 + (i-1)*order + j) = &
        (/pmap(i,j),pmap(i+1,j),pmap(i+1,j+1),pmap(i,j+1)/)-1
      END DO
    END DO
  END DO
  DEALLOCATE(xnodes,inodese,inodesf)
  DEALLOCATE(pmap)
  !CALL oft_abort('Invalid tessellation order','quadmesh_tessellate',__FILE__)
END IF
DEBUG_STACK_POP
END SUBROUTINE quadmesh_tessellate
!------------------------------------------------------------------------------
!> Get sizes of arrays returned by @ref trimesh_tessellate
!------------------------------------------------------------------------------
FUNCTION quadmesh_tessellated_sizes(self) result(sizes)
CLASS(oft_quadmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
sizes(1)=self%np + self%ne*(self%tess_order-1) + self%nc*((self%tess_order-1)**2)
sizes(2)=self%nc*(self%tess_order**2)
END FUNCTION quadmesh_tessellated_sizes
!------------------------------------------------------------------------------
!> Map from physical to logical coordinates in a given cell
!------------------------------------------------------------------------------
SUBROUTINE quadmesh_phys2log(self,cell,pt,f)
class(oft_quadmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: pt(3) !< Physical position [3]
real(r8), intent(out) :: f(:) !< Logical coordinates within the cell [4]
real(r8) :: gop(3,3),v
integer(i4) :: k
CALL oft_abort("Operation not supported", "quadmesh_phys2log", __FILE__)
! call quadmesh_jacl(self,i,f,gop,v)
! do k=1,3
!   f(k)=1.d0+dot_product(pt-self%r(:,self%lc(k,i)),gop(:,k))
! end do
end SUBROUTINE quadmesh_phys2log
!---------------------------------------------------------------------------------
!> Map from logical to physical coordinates in a given cell
!---------------------------------------------------------------------------------
function quadmesh_log2phys(self,cell,f) result(pt)
class(oft_quadmesh), intent(in) :: self !< Mesh object
integer, intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8) :: pt(3) !< Physical position [3]
integer(i4) :: i,j,k
DEBUG_STACK_PUSH
pt=0.d0
! Get node points
do i=1,4
  pt=pt+lag_2d(self%inodesp(:,i),f,self%xnodes,self%order+1)*self%r(:,self%lc(i,cell))
end do
if(self%order>1)then
  ! Get edge points
  do i=1,4
    j=ABS(self%lce(i,cell))
    pt=pt+lag_2d(self%inodese(:,1,i),f,self%xnodes,self%order+1)*self%ho_info%r(:,self%ho_info%lep(1,j))
  end do
  ! Get face point
  DO k=1,self%ho_info%ncp
    pt=pt+lag_2d(self%inodesf(:,k),f,self%xnodes,self%order+1)*self%ho_info%r(:,self%ho_info%lcp(1,cell))
  END DO
end if
DEBUG_STACK_POP
end function quadmesh_log2phys
!---------------------------------------------------------------------------------
!> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
!---------------------------------------------------------------------------------
subroutine quadmesh_jacobian(self,cell,f,gop,j)
class(oft_quadmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [3]
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
integer(i4) :: i,k,l,m,inode(2)
real(r8) :: A(2,2),C(2,2),t(3,2),pt(3),ftmp
DEBUG_STACK_PUSH
if(self%order>2)call oft_abort('Invalid mesh order','quadmesh_jacobian',__FILE__)
A=0.d0
IF(self%dim==2)THEN
  !---Get node points
  do l=1,4
    pt=self%r(:,self%lc(l,cell))
    DO k=1,2
      A(k,:)=A(k,:)+dlag_2d(self%inodesp(:,l),k,f,self%xnodes,self%order+1)*pt(1:2)
    END DO
  end do
  if(self%order>1)then
    !---Get edge points
    do l=1,4
      m=ABS(self%lce(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lep(1,m))
      DO k=1,2
        A(k,:)=A(k,:)+dlag_2d(self%inodese(:,1,l),k,f,self%xnodes,self%order+1)*pt(1:2)
      END DO
    end do
    !---Get face point
    DO i=1,self%ho_info%ncp
      pt=self%ho_info%r(:,self%ho_info%lcp(i,cell))
      DO k=1,2
        A(k,:)=A(k,:)+dlag_2d(self%inodesf(:,i),k,f,self%xnodes,self%order+1)*pt(1:2)
      END DO
    END DO
  end if
ELSE
  CALL quadmesh_tang(self,cell,f,t)
  !---Get node points
  do l=1,4
    pt=self%r(:,self%lc(l,cell))
    DO k=1,2
      ftmp=dlag_2d(self%inodesp(:,l),k,f,self%xnodes,self%order+1)
      A(k,1)=A(k,1)+ftmp*DOT_PRODUCT(pt,t(:,1))
      A(k,2)=A(k,2)+ftmp*DOT_PRODUCT(pt,t(:,2))
    END DO
  end do
  if(self%order>1)then
    !---Get edge points
    do l=1,4
      m=ABS(self%lce(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lep(1,m))
      DO k=1,2
        ftmp=dlag_2d(self%inodese(:,1,l),k,f,self%xnodes,self%order+1)
        A(k,1)=A(k,1)+ftmp*DOT_PRODUCT(pt,t(:,1))
        A(k,2)=A(k,2)+ftmp*DOT_PRODUCT(pt,t(:,2))
      END DO
    end do
    !---Get face point
    DO i=1,self%ho_info%ncp
      pt=self%ho_info%r(:,self%ho_info%lcp(i,cell))
      DO k=1,2
        ftmp=dlag_2d(self%inodesf(:,i),k,f,self%xnodes,self%order+1)
        A(k,1)=A(k,1)+ftmp*DOT_PRODUCT(pt,t(:,1))
        A(k,2)=A(k,2)+ftmp*DOT_PRODUCT(pt,t(:,2))
      END DO
    END DO
  end if
END IF
call quadmesh_jacinv(A,C,j)
IF(self%dim==2)THEN
  gop(:,1)=[C(1,1),C(2,1),0.d0]
  gop(:,2)=[C(1,2),C(2,2),0.d0]
ELSE
  gop(:,1)=C(1,1)*t(:,1)+C(2,1)*t(:,2)
  gop(:,2)=C(1,2)*t(:,1)+C(2,2)*t(:,2)
END IF
DEBUG_STACK_POP
end subroutine quadmesh_jacobian
!---------------------------------------------------------------------------------
!> Compute the spatial hessian matrices for a given cell at a given logical position
!!
!! @warning Not presently supported for quadrilateral meshes
!---------------------------------------------------------------------------------
subroutine quadmesh_hessian(self,cell,f,g2op,K)
class(oft_quadmesh), intent(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
REAL(r8), INTENT(out) :: g2op(:,:) !< Second order Jacobian matrix
!! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
REAL(r8), INTENT(out) :: K(:,:) !< Gradient correction matrix
!! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [10,3]
call oft_abort('Hessian not yet supported for quadmesh','quadmesh_hessian',__FILE__)
end subroutine quadmesh_hessian
!---------------------------------------------------------------------------------
!> Get unit normal for surface at a given point in a given cell
!---------------------------------------------------------------------------------
subroutine quadmesh_norm(self,cell,f,n)
class(oft_quadmesh), target, intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Cell containing point
real(r8), intent(in) :: f(:) !< Logical coordinates in cell
real(r8), intent(out) :: n(3) !< Unit normal [3]
real(r8) :: t(3,2)
DEBUG_STACK_PUSH
CALL quadmesh_tang(self,cell,f,t)
n = cross_product(t(:,1),t(:,2))
DEBUG_STACK_POP
end subroutine quadmesh_norm
!---------------------------------------------------------------------------------
!> Get tangent basis set for surface at a given point in a given cell
!---------------------------------------------------------------------------------
subroutine quadmesh_tang(self,cell,f,t)
class(oft_quadmesh), target, intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Cell containing point
real(r8), intent(in) :: f(:) !< Logical coordinates in cell
real(r8), intent(out) :: t(3,2) !< Unit tangent basis set [3,2]
real(r8) :: pt(3)
DEBUG_STACK_PUSH
!---Get tangent direction 1
t(:,1) = quadmesh_glogphys(self,cell,1,f)
t(:,1) = t(:,1)/magnitude(t(:,1))
!---Get tangent direction 2
t(:,2) = quadmesh_glogphys(self,cell,2,f)
t(:,2) = t(:,2) - DOT_PRODUCT(t(:,2),t(:,1))*t(:,1)
t(:,2) = t(:,2)/magnitude(t(:,2))
DEBUG_STACK_POP
end subroutine quadmesh_tang
!---------------------------------------------------------------------------------
!> Compute the partial derivative of the physical coordinates with a specific logical coordinate
!!
!! Driver function calls mapping specific function depending on mesh order
!---------------------------------------------------------------------------------
function quadmesh_glogphys(self,cell,j,f) result(pt)
class(oft_quadmesh), intent(in) :: self !< Mesh object
integer, intent(in) :: cell !< Index of cell for evaulation
integer, intent(in) :: j !< Needs docs
real(r8), intent(in) :: f(3) !< Logical coordinate in cell [4]
real(r8) :: pt(3) !< \f$ \frac{\partial r}{\partial f_k} \f$ [3]
real(r8) :: gtmp(2),pttmp(3)
integer(i4) :: k,l,ed,etmp(2),dof
DEBUG_STACK_PUSH
if(self%order>2)call oft_abort('Invalid mesh order','quadmesh_glogphys',__FILE__)
pt=0.d0
do k=1,4 ! Get edge nodes
  pt=pt+dlag_2d(self%inodesp(:,k),j,f,self%xnodes,self%order+1)*self%r(:,self%lc(k,cell))
end do
!---Get edge points
IF(self%order>1)THEN
  do k=1,4 ! Get edge nodes
    ed=ABS(self%lce(k,cell))
    pttmp=self%ho_info%r(:,self%ho_info%lep(1,ed))
    pt=pt+dlag_2d(self%inodese(:,1,k),j,f,self%xnodes,self%order+1)*pttmp
  end do
  IF(self%ho_info%ncp>0)THEN
    ! Add loop and orient if add higher order
    pttmp=self%ho_info%r(:,self%ho_info%lcp(1,cell))
    pt=pt+dlag_2d(self%inodesf(:,1),j,f,self%xnodes,self%order+1)*pttmp
  END IF
END IF
DEBUG_STACK_POP
end function quadmesh_glogphys
!---------------------------------------------------------------------------------
!> Invert a 2x2 matrix
!---------------------------------------------------------------------------------
subroutine quadmesh_jacinv(A,C,j)
real(r8), intent(in) :: A(2,2) !< Matrix to invert
real(r8), intent(out) :: C(2,2) !< \f$ A^{-1} \f$
real(r8), intent(out) :: j !< |A|
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
end subroutine quadmesh_jacinv
END MODULE oft_quadmesh_type
