!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
MODULE oft_trimesh_type
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_field_get_sizes, hdf5_write, hdf5_create_file, hdf5_create_group
USE oft_lag_poly
USE oft_tet_quadrature, ONLY: oft_quad_type, set_quad_2d
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_destroy, cell_is_curved
USE trimesh_tessellation, ONLY: tessellate1, tessellate2, tessellate3, tessellate4
IMPLICIT NONE
#include "local.h"
INTEGER(i4), PARAMETER :: tri_ed(2,3)=RESHAPE((/3,2,1,3,2,1/),(/2,3/)) !< Triangle edge list
!---------------------------------------------------------------------------------
!> Triangular surface mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!---------------------------------------------------------------------------------
TYPE, EXTENDS(oft_bmesh) :: oft_trimesh
  REAL(r8), POINTER, DIMENSION(:) :: xnodes => NULL()
CONTAINS
  PROCEDURE :: setup => trimesh_setup
  PROCEDURE :: load_from_file => trimesh_load
  PROCEDURE :: save_to_file => trimesh_save
  PROCEDURE :: set_order => trimesh_set_order
  PROCEDURE :: invert_face => trimesh_invert_cell
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
  PROCEDURE :: tessellated_sizes => trimesh_tessellated_sizes
  !> Destroy mesh object
  PROCEDURE :: delete => trimesh_destroy
END TYPE oft_trimesh
CONTAINS
!------------------------------------------------------------------------------
!> Destroy mesh object
!------------------------------------------------------------------------------
SUBROUTINE trimesh_destroy(self)
CLASS(oft_trimesh), INTENT(inout) :: self
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
CALL bmesh_destroy(self)
DEBUG_STACK_POP
END SUBROUTINE trimesh_destroy
!------------------------------------------------------------------------------
!> Load trimesh from transfer file
!------------------------------------------------------------------------------
SUBROUTINE trimesh_load(self,filename,ascii)
CLASS(oft_trimesh), INTENT(inout) :: self !< Mesh object
CHARACTER(LEN=*), INTENT(in) :: filename !< File to load mesh from
LOGICAL, OPTIONAL, INTENT(in) :: ascii !< Mesh stored in ASCII format?
INTEGER(i4) :: i,io_unit,ndims
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: dim_sizes
LOGICAL :: read_ascii,success
DEBUG_STACK_PUSH
read_ascii=.FALSE.
IF(PRESENT(ascii))read_ascii=ascii
IF(read_ascii)THEN
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
ELSE
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/R",ndims,dim_sizes)
  IF(ndims==-1)CALL oft_abort('"mesh/R" field does not exist', 'trimesh_load', __FILE__)
  self%np=dim_sizes(2)
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/LC",ndims,dim_sizes)
  IF(ndims==-1)CALL oft_abort('"mesh/LC" field does not exist', 'trimesh_load', __FILE__)
  IF(dim_sizes(1)/=3)CALL oft_abort('"mesh/LC" has wrong first dimension (must be 3)', 'trimesh_load', __FILE__)
  self%nc=dim_sizes(2)
  ALLOCATE(self%r(3,self%np))
  CALL hdf5_read(self%r,TRIM(filename),"mesh/R",success)
  IF(.NOT.success)CALL oft_abort('Error reading points','trimesh_load',__FILE__)
  ALLOCATE(self%lc(3,self%nc))
  CALL hdf5_read(self%lc,TRIM(filename),"mesh/LC",success)
  IF(.NOT.success)CALL oft_abort('Error reading cells','trimesh_load',__FILE__)
END IF
ALLOCATE(self%reg(self%nc))
self%reg=1
! self%stand_alone=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE trimesh_load
!------------------------------------------------------------------------------
!> Save trimesh from transfer file
!------------------------------------------------------------------------------
SUBROUTINE trimesh_save(self,filename,ascii)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
CHARACTER(LEN=*), INTENT(in) :: filename !< File to save mesh to
LOGICAL, OPTIONAL, INTENT(in) :: ascii !< Save file in ASCII format?
INTEGER(i4) :: i,io_unit
LOGICAL :: write_ascii
DEBUG_STACK_PUSH
write_ascii=.FALSE.
IF(PRESENT(ascii))write_ascii=ascii
IF(write_ascii)THEN
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
ELSE
  CALL hdf5_create_file(TRIM(filename))
  CALL hdf5_create_group(TRIM(filename),'mesh')
  CALL hdf5_write(self%r,TRIM(filename),"mesh/R")
  CALL hdf5_write(self%lc,TRIM(filename),"mesh/LC")
END IF
DEBUG_STACK_POP
END SUBROUTINE trimesh_save
!------------------------------------------------------------------------------
!> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
!------------------------------------------------------------------------------
SUBROUTINE trimesh_setup(self,cad_type,has_parent)
CLASS(oft_trimesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
LOGICAL, INTENT(in) :: has_parent !< Is this mesh the/a surface of a volume mesh?
self%cell_np=3
self%cell_ne=self%cell_np
self%type=1
self%cad_type=cad_type
ALLOCATE(self%cell_ed(2,3))
self%cell_ed=tri_ed
CALL self%set_order(1)
IF(has_parent)ALLOCATE(self%parent)
END SUBROUTINE trimesh_setup
!------------------------------------------------------------------------------
!> Set maximum order of spatial mapping
!------------------------------------------------------------------------------
SUBROUTINE trimesh_set_order(self,order)
CLASS(oft_trimesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
INTEGER(i4) :: i,counts(3)
! Cubic triangles would also need nep=2, ncp=1 and the matching terms in trimesh_ho_map
IF(order<1.OR.order>2)CALL oft_abort("Invalid grid order", "trimesh_set_order", __FILE__)
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Turn cell "inside out", used to ensure consistent orientations
!------------------------------------------------------------------------------
SUBROUTINE trimesh_invert_cell(self,cell)
CLASS(oft_trimesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Maximum order of spatial mapping
self%lc(2:3,cell)=self%lc(3:2:-1,cell)
IF(ASSOCIATED(self%lce))self%lce(:,cell)=-self%lce([1,3,2],cell)
IF(ASSOCIATED(self%lcc))self%lcc(:,cell)=self%lcc([1,3,2],cell)
END SUBROUTINE trimesh_invert_cell
!---------------------------------------------------------------------------------
!> Retrieve suitable quadrature rule for triangular mesh with given order
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_quad_rule(self,order,quad_rule)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
CALL set_quad_2d(quad_rule, order)
END SUBROUTINE trimesh_quad_rule
!---------------------------------------------------------------------------------
!> Get position in logical space of vertex `i`
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_vlog(self,i,f)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: i !< Vertex to locate
REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
f=0.d0
f(i)=1.d0
END SUBROUTINE trimesh_vlog
!---------------------------------------------------------------------------------
!> Test if logical position lies within the base cell
!!
!! @returns Position `f` is inside the base cell?
!---------------------------------------------------------------------------------
FUNCTION trimesh_in_cell(self,f,tol) RESULT(eedge)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
REAL(r8), INTENT(in) :: tol !< Tolerance for test
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
!---------------------------------------------------------------------------------
!> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
!!
!! @note The maximum tessellation order currently supported is 4
!! (may be lower for certain mesh types).
!!
!! @warning Cell lists are returned with zero based indexing
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_tessellate(self,rtmp,lctmp,order)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,:]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [self%ncp,:]
INTEGER(i4), INTENT(in) :: order !< Tessellation order
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
!------------------------------------------------------------------------------
!> Get sizes of arrays returned by @ref trimesh_tessellate
!------------------------------------------------------------------------------
FUNCTION trimesh_tessellated_sizes(self) result(sizes)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
INTEGER(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
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
  CALL oft_abort("Unknown tessellation size","trimesh_tessellated_sizes",__FILE__)
END SELECT
END FUNCTION trimesh_tessellated_sizes
!---------------------------------------------------------------------------------
!> Map from physical to logical coordinates in a given cell
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_phys2log(self,cell,pt,f)
class(oft_trimesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: pt(3) !< Physical position [3]
real(r8), intent(out) :: f(:) !< Logical coordinates within the cell [4]
real(r8) :: gop(3,3),v
integer(i4) :: k
call trimesh_jacl(self,cell,f,gop,v)
do k=1,3
  f(k)=1.d0+dot_product(pt-self%r(:,self%lc(k,cell)),gop(:,k))
end do
end SUBROUTINE trimesh_phys2log
!---------------------------------------------------------------------------------
!> Map from logical to physical coordinates in a given cell
!---------------------------------------------------------------------------------
function trimesh_log2phys(self,cell,f) result(pt)
class(oft_trimesh), intent(in) :: self !< Mesh object
integer, intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [3]
real(r8) :: pt(3) !< Physical position [3]
INTEGER(i4) :: k
DEBUG_STACK_PUSH
IF(self%order>2)CALL oft_abort('Invalid mesh order','trimesh_log2phys',__FILE__)
IF(cell_is_curved(self,cell))THEN
  CALL trimesh_ho_map(self,cell,f,pt=pt)
ELSE
  pt=0.d0
  DO k=1,3
    pt=pt+f(k)*self%r(:,self%lc(k,cell))
  END DO
END IF
DEBUG_STACK_POP
end function trimesh_log2phys
!---------------------------------------------------------------------------------
!> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
!---------------------------------------------------------------------------------
subroutine trimesh_jacobian(self,cell,f,gop,j)
class(oft_trimesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [3]
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
REAL(r8) :: jfull(3,3),A(2,2),C(2,2),t(3,2)
INTEGER(i4) :: k
DEBUG_STACK_PUSH
IF(self%order>2)CALL oft_abort('Invalid mesh order','trimesh_jacobian',__FILE__)
CALL trimesh_dmap(self,cell,f,jfull)
!---Condense Jacobian (lambda_1 eliminated)
IF(self%dim==2)THEN
  DO k=1,2
    A(k,:)=jfull(1:2,k+1)-jfull(1:2,1)
  END DO
ELSE
  !---Project onto the local tangent plane
  CALL trimesh_tang_basis(jfull,t)
  DO k=1,2
    A(k,:)=MATMUL(jfull(:,k+1)-jfull(:,1),t)
  END DO
END IF
!---Comute inverse
CALL trimesh_jacinv(A,C,j)
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
!---------------------------------------------------------------------------------
!> Compute the spatial hessian matrices for a given cell at a given logical position
!!
!! Logical second derivatives are mapped with \f$ \lambda_3 \f$ eliminated. For curved
!! cells the Jacobian is evaluated at `f` and `K` holds the second derivatives of the
!! high-order mapping, so the physical Hessian of a field \f$ u \f$ is
!! `MATMUL(g2op, d2u - MATMUL(K, grad_u))`, where `d2u` are the logical second derivatives
!! and `grad_u` is the physical gradient (`K=0` on straight-sided cells).
!---------------------------------------------------------------------------------
subroutine trimesh_hessian(self,cell,f,g2op,K)
class(oft_trimesh), intent(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [3]
REAL(r8), INTENT(out) :: g2op(:,:) !< Second order Jacobian matrix
!! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$ [6,6]
REAL(r8), INTENT(out) :: K(:,:) !< Gradient correction matrix
!! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [6,3]
REAL(r8) :: jfull(3,3),Kfull(6,3),jac(2,2),A(3,3),C(3,3),v
INTEGER(i4) :: m
DEBUG_STACK_PUSH
IF(self%order>2)CALL oft_abort('Invalid mesh order','trimesh_hessian',__FILE__)
! Uses only the in-plane (x,y) components of the mapping; surfaces in 3D would also need
! the tangent-plane projection from trimesh_jacobian, and its variation for K
IF(self%dim==3)CALL oft_abort('Only supported for dim=2','trimesh_hessian',__FILE__)
CALL trimesh_dmap(self,cell,f,jfull,Kfull)
K=0.d0
K(1:6,1:2)=Kfull(:,1:2)
!---Get Jacobian entries (lambda_3 eliminated)
DO m=1,2
  jac(m,:)=jfull(1:2,m)-jfull(1:2,3)
END DO
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
CALL trimesh_m3inv(A,C,v)
!---Map back to g2op
g2op=0.d0
g2op([1,2,4],1)=C(:,1)
g2op([1,2,4],2)=C(:,2)
g2op([1,2,4],4)=C(:,3)
!---Scatter to redundant dimension
g2op([1,2,4],3)=-2*C(:,1)-C(:,2)
g2op([1,2,4],5)=-C(:,2)-2*C(:,3)
g2op([1,2,4],6)=C(:,1)+C(:,2)+C(:,3)
DEBUG_STACK_POP
end subroutine trimesh_hessian
!---------------------------------------------------------------------------------
!> Get the derivatives of the physical coordinates with respect to the logical
!! coordinates for a given cell at a given logical position
!!
!! All three logical coordinates are treated as independent, so callers must
!! eliminate the redundant one. Straight-sided cells use the vertex positions directly,
!! for which `jfull` is constant and `K` is zero.
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_dmap(self,cell,f,jfull,K)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [3]
REAL(r8), INTENT(out) :: jfull(3,3) !< \f$ \frac{\partial x_i}{\partial \lambda_j} \f$ [3,3]
REAL(r8), OPTIONAL, INTENT(out) :: K(6,3) !< \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$,
!! with rows ordered (11,12,13,22,23,33) [6,3]
INTEGER(i4) :: m
DEBUG_STACK_PUSH
IF(cell_is_curved(self,cell))THEN
  CALL trimesh_ho_map(self,cell,f,jfull=jfull,K=K)
ELSE
  DO m=1,3
    jfull(:,m)=self%r(:,self%lc(m,cell))
  END DO
  IF(PRESENT(K))K=0.d0
END IF
DEBUG_STACK_POP
END SUBROUTINE trimesh_dmap
!---------------------------------------------------------------------------------
!> Evaluate the quadratic mapping of a curved cell and its logical derivatives
!!
!! Uses the quadratic Lagrange shape functions \f$ \lambda_i (2 \lambda_i - 1) \f$ for
!! vertices and \f$ 4 \lambda_a \lambda_b \f$ for edge nodes. All three logical coordinates
!! are treated as independent, so callers must eliminate the redundant one.
!!
!! @note Only quadratic triangles are supported (see `trimesh_set_order`).
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_ho_map(self,cell,f,pt,jfull,K)
CLASS(oft_trimesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [3]
REAL(r8), OPTIONAL, INTENT(out) :: pt(3) !< Physical position [3]
REAL(r8), OPTIONAL, INTENT(out) :: jfull(3,3) !< \f$ \frac{\partial x_i}{\partial \lambda_j} \f$ [3,3]
REAL(r8), OPTIONAL, INTENT(out) :: K(6,3) !< \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$,
!! with rows ordered (11,12,13,22,23,33) [6,3]
REAL(r8) :: r(3)
INTEGER(i4) :: m,a,b
INTEGER(i4), PARAMETER :: pmap(3)=(/1,4,6/)
INTEGER(i4), PARAMETER :: fmap(3,3)=RESHAPE((/1,2,3,2,4,5,3,5,6/),(/3,3/))
DEBUG_STACK_PUSH
IF(PRESENT(pt))pt=0.d0
IF(PRESENT(jfull))jfull=0.d0
IF(PRESENT(K))K=0.d0
!---Vertex nodes
DO m=1,3
  r=self%r(:,self%lc(m,cell))
  IF(PRESENT(pt))pt=pt+f(m)*(2.d0*f(m)-1.d0)*r
  IF(PRESENT(jfull))jfull(:,m)=jfull(:,m)+(4.d0*f(m)-1.d0)*r
  IF(PRESENT(K))K(pmap(m),:)=K(pmap(m),:)+4.d0*r
END DO
!---Edge nodes (cubic triangles would need two oriented nodes per edge here)
DO m=1,3
  a=tri_ed(1,m)
  b=tri_ed(2,m)
  r=self%ho_info%r(:,self%ho_info%lep(1,ABS(self%lce(m,cell))))
  IF(PRESENT(pt))pt=pt+4.d0*f(a)*f(b)*r
  IF(PRESENT(jfull))THEN
    jfull(:,a)=jfull(:,a)+4.d0*f(b)*r
    jfull(:,b)=jfull(:,b)+4.d0*f(a)*r
  END IF
  IF(PRESENT(K))K(fmap(a,b),:)=K(fmap(a,b),:)+4.d0*r
END DO
!---Cubic triangles would also need the cell node (ho_info%lcp) here
DEBUG_STACK_POP
END SUBROUTINE trimesh_ho_map
!---------------------------------------------------------------------------------
!> Linear implementation of @trimesh_jacobian
!---------------------------------------------------------------------------------
subroutine trimesh_jacl(self,cell,f,gop,j)
class(oft_trimesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [3]
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
real(r8) :: ri(3,3),jfull(2,3),A(2,2),C(2,2),t(3,2)
integer(i4) :: k
DEBUG_STACK_PUSH
IF(self%dim==2)THEN
  !---Get Jacobian entries
  DO k=1,3
    jfull(:,k)=self%r(1:2,self%lc(k,cell))
  END DO
ELSE
  !---Get node points
  do k=1,3
    ri(:,k)=self%r(:,self%lc(k,cell))
  end do
  CALL trimesh_tang(self,cell,f,t)
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
!---------------------------------------------------------------------------------
!> Get unit normal for surface at a given point in a given cell
!---------------------------------------------------------------------------------
subroutine trimesh_norm(self,cell,f,n)
class(oft_trimesh), target, intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Cell containing point
real(r8), intent(in) :: f(:) !< Logical coordinates in cell
real(r8), intent(out) :: n(3) !< Unit normal [3]
real(r8) :: t(3,2)
DEBUG_STACK_PUSH
CALL trimesh_tang(self,cell,f,t)
n = cross_product(t(:,1),t(:,2))
DEBUG_STACK_POP
end subroutine trimesh_norm
!---------------------------------------------------------------------------------
!> Get tangent basis set for surface at a given point in a given cell
!---------------------------------------------------------------------------------
subroutine trimesh_tang(self,cell,f,t)
class(oft_trimesh), target, intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Cell containing point
real(r8), intent(in) :: f(:) !< Logical coordinates in cell
real(r8), intent(out) :: t(3,2) !< Unit tangent basis set [3,2]
REAL(r8) :: jfull(3,3)
DEBUG_STACK_PUSH
CALL trimesh_dmap(self,cell,f,jfull)
CALL trimesh_tang_basis(jfull,t)
DEBUG_STACK_POP
end subroutine trimesh_tang
!---------------------------------------------------------------------------------
!> Build an orthonormal tangent basis from the logical derivatives of the mapping
!!
!! For straight-sided cells `jfull` holds the vertex positions, whose differences
!! are the edge vectors.
!---------------------------------------------------------------------------------
SUBROUTINE trimesh_tang_basis(jfull,t)
REAL(r8), INTENT(in) :: jfull(3,3) !< \f$ \frac{\partial x_i}{\partial \lambda_j} \f$ (see `trimesh_dmap`) [3,3]
REAL(r8), INTENT(out) :: t(3,2) !< Unit tangent basis set [3,2]
DEBUG_STACK_PUSH
!---Get tangent direction 1
t(:,1) = jfull(:,2) - jfull(:,1)
t(:,1) = t(:,1)/magnitude(t(:,1))
!---Get tangent direction 2
t(:,2) = jfull(:,3) - jfull(:,1)
t(:,2) = t(:,2) - DOT_PRODUCT(t(:,2),t(:,1))*t(:,1)
t(:,2) = t(:,2)/magnitude(t(:,2))
DEBUG_STACK_POP
END SUBROUTINE trimesh_tang_basis
!---------------------------------------------------------------------------------
!> Invert a 2x2 matrix
!---------------------------------------------------------------------------------
subroutine trimesh_jacinv(A,C,j)
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
end subroutine trimesh_jacinv
!---------------------------------------------------------------------------------
!> Invert a 3x3 matrix
!---------------------------------------------------------------------------------
subroutine trimesh_m3inv(A,C,j)
real(8), intent(in) :: A(3,3) !< Matrix to invert
real(8), intent(out) :: C(3,3) !< \f$ A^{-1} \f$
real(8), intent(out) :: j !< |A|
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
