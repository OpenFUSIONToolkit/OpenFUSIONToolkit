!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_tetmesh_type.F90
!
!> Tetrahedral mesh structure definitions
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
MODULE oft_tetmesh_type
USE oft_base
USE oft_lag_poly
USE oft_tet_quadrature, ONLY: set_quad_3d, oft_quad_type
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
USE tetmesh_tessellation, ONLY: tessellate1, tessellate2, tessellate3, tessellate4
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Tetrahedral volume mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_mesh) :: oft_tetmesh
  REAL(r8), POINTER, DIMENSION(:) :: xnodes => NULL()
CONTAINS
  PROCEDURE :: setup => tetmesh_setup
  PROCEDURE :: set_order => tetmesh_set_order
  PROCEDURE :: invert_cell => tetmesh_invert_cell
  PROCEDURE :: log2phys => tetmesh_log2phys
  PROCEDURE :: phys2log => tetmesh_phys2log
  PROCEDURE :: jacobian => tetmesh_jacobian
  PROCEDURE :: hessian => tetmesh_hessian
  PROCEDURE :: snormal => tetmesh_snormal
  PROCEDURE :: ctang => tetmesh_ctang
  PROCEDURE :: get_surf_map => tetmesh_get_surf_map
  PROCEDURE :: surf_to_vol => tetmesh_surf_to_vol
  PROCEDURE :: vlog => tetmesh_vlog
  PROCEDURE :: in_cell => tetmesh_in_cell
  PROCEDURE :: quad_rule => tetmesh_quad_rule
  PROCEDURE :: tessellate => tetmesh_tessellate
  PROCEDURE :: tessellated_sizes => tetmesh_tessellated_sizes
END TYPE oft_tetmesh
!---
INTEGER(i4), PARAMETER :: tet_ed(2,6)=RESHAPE((/1,4, 2,4, 3,4, 2,3, 3,1, 1,2/),(/2,6/)) !< Tetrahedron edge list
INTEGER(i4), PARAMETER :: tet_fc(3,4)=RESHAPE((/2,3,4,3,1,4,1,2,4,1,2,3/),(/3,4/)) !< Tetrahedron face list
INTEGER(i4), PARAMETER :: tet_fe(3,4)=RESHAPE((/2,3,4, 1,3,5, 1,2,6, 4,5,6/),(/3,4/)) !< Tetrahedron face edge list
INTEGER(i4), PARAMETER :: tri_ed(2,3)=RESHAPE((/3,2,1,3,2,1/),(/2,3/)) !< Triangle edge list
!
INTEGER(i4), PRIVATE, PARAMETER :: ho_find_nsteps=100 !< Maximum number of steps during high order find_cell
REAL(r8), PRIVATE, PARAMETER :: ho_find_du=1.d-6 !< Step size used for jacobian eval during high order find_cell
REAL(r8), PRIVATE, PARAMETER :: ho_find_tol=1.d-6 !< Convergence tolerance for high order find_cell
CLASS(oft_tetmesh), PRIVATE, POINTER :: active_mesh => NULL() !< Active mesh for high order find_cell
INTEGER(i4), PRIVATE :: active_cell = 0 !< Active cell for high order find_cell
REAL(r8), PRIVATE :: active_pt(3) = 0.d0 !< Active point for high order find_cell
!$omp threadprivate(active_mesh,active_cell,active_pt)
CONTAINS
!---------------------------------------------------------------------------
!> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
!---------------------------------------------------------------------------
SUBROUTINE tetmesh_setup(self,cad_type)
CLASS(oft_tetmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
self%type=1
self%cad_type=cad_type
! self%cell_geom=(/4,6,4/)
self%cell_np=4
self%cell_ne=6
self%cell_nf=4
self%face_np=3
ALLOCATE(self%cell_ed(2,6),self%cell_fc(3,4),self%cell_fe(3,4))
self%cell_ed=tet_ed
self%cell_fc=tet_fc
self%cell_fe=tet_fe
ALLOCATE(self%face_ed(2,3))
self%face_ed=tri_ed
CALL self%set_order(1)
END SUBROUTINE tetmesh_setup
!---------------------------------------------------------------------------
!> Set maximum order of spatial mapping
!---------------------------------------------------------------------------
SUBROUTINE tetmesh_set_order(self,order)
CLASS(oft_tetmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
INTEGER(i4) :: i,counts(3)
IF(order<1.OR.order>3)CALL oft_abort("Invalid grid order", "tetmesh_set_order", __FILE__)
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
  IF(self%ho_info%nfp>0)THEN
    ALLOCATE(self%ho_info%lfp(self%ho_info%nfp,self%nf))
    self%ho_info%lfp=-1
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
END SUBROUTINE tetmesh_set_order
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
SUBROUTINE tet_3d_grid(order,xnodes,inodesf,inodesc)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4), POINTER, DIMENSION(:,:,:), INTENT(out) :: inodesf
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: inodesc
REAL(r8), POINTER, DIMENSION(:), INTENT(out) :: xnodes
INTEGER(i4) :: i,j,k,np,inc1,ind1,inc2,ind2,inode(3),counts(3)
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
  ALLOCATE(inodesf(2,np,1))
  np=0
  DO i=1,order-1
    DO j=1,order-1-i
      np=np+1
      inodesf(:,np,1)=(/i,j/)
    END DO
  END DO
END IF
!---Cell
IF(order>3)THEN
  np=0
  DO i=1,order-1
    DO j=1,order-1-i
      DO k=1,order-1-i-j
        np=np+1
      END DO
    END DO
  END DO
  ALLOCATE(inodesc(3,np))
  np=0
  DO i=1,order-1
    DO j=1,order-1-i
      DO k=1,order-1-i-j
        np=np+1
        inodesc(:,np)=(/i,j,k/)
      END DO
    END DO
  END DO
END IF
END SUBROUTINE tet_3d_grid
!---------------------------------------------------------------------------
!> Turn cell "inside out", used to ensure consistent orientations
!---------------------------------------------------------------------------
SUBROUTINE tetmesh_invert_cell(self,cell)
CLASS(oft_tetmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell to invert
self%lc(3:4,cell)=self%lc(4:3:-1,cell)
END SUBROUTINE tetmesh_invert_cell
!------------------------------------------------------------------------------
!> Retrieve suitable quadrature rule for mesh with given order
!------------------------------------------------------------------------------
SUBROUTINE tetmesh_quad_rule(self,order,quad_rule)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
CALL set_quad_3d(quad_rule, order)
END SUBROUTINE tetmesh_quad_rule
!------------------------------------------------------------------------------
!> Get position in logical space of vertex `i`
!------------------------------------------------------------------------------
SUBROUTINE tetmesh_vlog(self,i,f)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: i !< Vertex to locate
REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
f=0.d0
f(i)=1.d0
END SUBROUTINE tetmesh_vlog
!------------------------------------------------------------------------------
!> Test if logical position lies within the base cell
!!
!! @returns Position `f` is inside the base cell?
!------------------------------------------------------------------------------
FUNCTION tetmesh_in_cell(self,f,tol) RESULT(eface)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
REAL(r8), INTENT(in) :: tol !< Tolerance for test
INTEGER(i4) :: eface
REAL(r8) :: fmin,fmax
fmin=MINVAL(f(1:4))
fmax=MAXVAL(f(1:4))
IF((fmax<=1.d0+tol).AND.(fmin>=-tol))THEN
  eface=0
ELSE
  eface=MINLOC(f(1:4), DIM=1)
END IF
END FUNCTION tetmesh_in_cell
!------------------------------------------------------------------------------
!> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
!!
!! @note The maximum tessellation order currently supported is 4
!! (may be lower for certain mesh types).
!!
!! @warning Cell lists are returned with zero based indexing
!------------------------------------------------------------------------------
SUBROUTINE tetmesh_tessellate(self,rtmp,lctmp,order)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,:]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [4,:]
INTEGER(i4), INTENT(in) :: order !< Tessellation order
DEBUG_STACK_PUSH
!---Call desired driver
IF(order==1)THEN
  CALL tessellate1(self,rtmp,lctmp)
ELSE IF(order==2)THEN
  CALL tessellate2(self,rtmp,lctmp)
ELSE IF(order==3)THEN
  CALL tessellate3(self,rtmp,lctmp)
ELSE IF(order==4)THEN
  CALL tessellate4(self,rtmp,lctmp)
ELSE
  CALL oft_abort('Invalid tessellation order','tetmesh_tesselate',__FILE__)
END IF
DEBUG_STACK_POP
END SUBROUTINE tetmesh_tessellate
!---------------------------------------------------------------------------
!> Get sizes of arrays returned by @ref tetmesh_tessellate
!---------------------------------------------------------------------------
FUNCTION tetmesh_tessellated_sizes(self) result(sizes)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
SELECT CASE(self%tess_order)
CASE(1)
  sizes=[self%np, self%nc]
CASE(2)
  sizes=[self%np+self%ne, self%nc*8]
CASE(3)
  sizes=[self%np+2*self%ne+self%nf, self%nc*27]
CASE(4)
  sizes=[self%np+3*self%ne+3*self%nf+self%nc, self%nc*64]
CASE DEFAULT
  CALL oft_abort("Unkown tessellation size","tetmesh_tessellated_sizes",__FILE__)
END SELECT
END FUNCTION tetmesh_tessellated_sizes
!------------------------------------------------------------------------------
!> Map from logical to physical coordinates in a given cell
!------------------------------------------------------------------------------
function tetmesh_log2phys(self,cell,f) result(pt)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer, intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8) :: pt(3) !< Physical position [3]
integer(i4) :: i
DEBUG_STACK_PUSH
if(self%order>3)call oft_abort('Invalid mesh order','tetmesh_log2phys',__FILE__)
pt=0.d0
IF(cell_is_curved(self, cell))THEN
  SELECT CASE(self%order)
    CASE(2)
      CALL log2phys_quad()
    CASE DEFAULT
      CALL log2phys_gen()
  END SELECT
ELSE
  do i=1,4
    pt=pt + f(i)*self%r(:,self%lc(i,cell))
  end do
END IF
DEBUG_STACK_POP
contains
!
subroutine log2phys_quad()
integer(i4) :: k,ed
! Vertex nodes
do k=1,4
  pt = pt + f(k)*(f(k) - 0.5d0)*2.d0*self%r(:,self%lc(k,cell))
end do
! Edge nodes
do k=1,6
  ed=ABS(self%lce(k,cell))
  pt = pt + 4.d0*f(tet_ed(1,k))*f(tet_ed(2,k))*self%ho_info%r(:,self%ho_info%lep(1,ed))
end do
end subroutine log2phys_quad
!
subroutine log2phys_gen()
integer(i4) :: i,j,k,l,ed,dof,etmp(2)
! Vertex nodes
do i=1,4
  pt = pt + lag_1d(self%order+1,f(i),self%xnodes,self%order+1) &
    *self%r(:,self%lc(i,cell))
end do
! Edge nodes
DO i=1,6
  ed=self%lce(i,cell)
  j=abs(ed)
  DO k=1,self%ho_info%nep
    dof=k
    IF(ed<0)dof=self%ho_info%nep+1-k
    pt = pt + lag_1d_bary(k,f(tet_ed(:,i)),self%xnodes,self%order+1) &
      *self%ho_info%r(:,self%ho_info%lep(dof,j))
  END DO
END DO
! Face nodes
IF(self%ho_info%nfp>0)THEN
  do i=1,4
    j=ABS(self%lcf(i,cell))
    pt = pt + lag_2d_bary((/1,1/),f(tet_fc(:,i)),self%xnodes,self%order+1) &
      *self%ho_info%r(:,self%ho_info%lfp(1,j))
  end do
END IF
end subroutine log2phys_gen
end function tetmesh_log2phys
!------------------------------------------------------------------------------
!> Map from physical to logical coordinates in a given cell
!------------------------------------------------------------------------------
subroutine tetmesh_phys2log(self,cell,pt,f)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: pt(3) !< Physical position [3]
real(r8), intent(out) :: f(:) !< Logical coordinates within the cell [4]
DEBUG_STACK_PUSH
if(cell_is_curved(self,cell))then
  f=tetmesh_phys2logho(self,cell,pt)
else
  f=tetmesh_phys2logl(self,cell,pt)
end if
DEBUG_STACK_POP
contains
!------------------------------------------------------------------------------
! Linear element implementation of @ref tetmesh_phys2log
!------------------------------------------------------------------------------
function tetmesh_phys2logl(self,i,pt) result(f)
class(oft_tetmesh), intent(in) :: self
integer(i4), intent(in) :: i
real(r8), intent(in) :: pt(3)
real(r8) :: f(4),gop(3,4),v
integer(i4) :: k
DEBUG_STACK_PUSH
call tetmesh_jacl(self,i,gop,v)
do k=1,4
  f(k)=1.d0+dot_product(pt-self%r(:,self%lc(k,i)),gop(:,k))
end do
DEBUG_STACK_POP
end function tetmesh_phys2logl
!------------------------------------------------------------------------------
! General high-order implementation of @ref tetmesh_phys2log
!
! The MINPACK package is used with step size given by
! @ref tetmesh_mapping::ho_find_du "ho_find_du". The convergence tolerance is
! set by the variable @ref tetmesh_mapping::ho_find_tol "ho_find_tol".
!
! @note The final location may be outside the cell being searched. This is correct
! if the point is outside the cell, however it may also indicate a problem in the
! mapping, most likely due to a badly shaped cell
!------------------------------------------------------------------------------
function tetmesh_phys2logho(self,i,pt) result(f)
class(oft_tetmesh), target, intent(in) :: self
integer(i4), intent(in) :: i
real(r8), intent(in) :: pt(3)
real(r8) :: f(4)
!---
integer(i4), parameter :: nerr=3
integer(i4), parameter :: neq=3
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor,dmin
integer :: info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---
mode = 1
factor = 1.d0
maxfev = ho_find_nsteps
ftol = ho_find_tol
xtol = ho_find_tol
gtol = ho_find_tol
epsfcn = ho_find_du
nprint = 0
ldfjac = nerr
!---
active_mesh=>self
active_pt=pt
active_cell=i
uv=1.d0/4.d0
!---
call lmdif(tm_findcell_error,nerr,neq,uv,error, &
            ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
            nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!IF(info>4)WRITE(*,*)'High-order find failed',i,info,nfev
f(1:3)=uv; f(4)=1.d0-SUM(uv)
DEBUG_STACK_POP
end function tetmesh_phys2logho
!---------------------------------------------------------------------------
!> Evalute the error between a logical point and the current active point
!!
!! @note Designed to be used as the error function for minimization in
!! @ref tetmesh_mapping::tetmesh_phys2logho "tetmesh_phys2logho"
!---------------------------------------------------------------------------
subroutine tm_findcell_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions (3)
integer(i4), intent(in) :: n !< Number of parametric dimensions (3)
real(r8), intent(in) :: uv(n) !< Parametric possition [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3),f(4)
DEBUG_STACK_PUSH
f(1:3)=uv; f(4)=1.d0-SUM(uv)
pt=tetmesh_log2phys(active_mesh,active_cell,f)
err=active_pt-pt
DEBUG_STACK_POP
end subroutine tm_findcell_error
end subroutine tetmesh_phys2log
!------------------------------------------------------------------------------
!> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
!------------------------------------------------------------------------------
subroutine tetmesh_jacobian(self,cell,f,gop,j)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
real(r8) :: jfull(3,4)
integer(i4) :: k
DEBUG_STACK_PUSH
if(self%order>3)call oft_abort('Invalid mesh order','tetmesh_jacobian',__FILE__)
jfull=0.d0
IF(cell_is_curved(self, cell))THEN
  SELECT CASE(self%order)
    CASE(2)
      CALL jacobian_quad()
    CASE DEFAULT
      CALL jacobian_gen()
  END SELECT
ELSE
  ! Get node points
  do k=1,4
    jfull(:,k)=self%r(:,self%lc(k,cell))
  end do
END IF
call tetmesh_jacinv(jfull,gop,j)
DEBUG_STACK_POP
contains
!------------------------------------------------------------------------------
! Quadratic element implementation
!------------------------------------------------------------------------------
subroutine jacobian_quad()
real(r8) :: pt(3)
integer(i4) :: k,l,ed
! Vertex nodes
do k=1,4
  jfull(:,k) = jfull(:,k) + (4.d0*f(k) - 1.d0)*self%r(:,self%lc(k,cell))
end do
! Edge nodes
do k=1,6
  ed=ABS(self%lce(k,cell))
  pt=self%ho_info%r(:,self%ho_info%lep(1,ed))
  jfull(:,tet_ed(1,k)) = jfull(:,tet_ed(1,k)) + 4.d0*f(tet_ed(2,k))*pt
  jfull(:,tet_ed(2,k)) = jfull(:,tet_ed(2,k)) + 4.d0*f(tet_ed(1,k))*pt
end do
end subroutine jacobian_quad
!------------------------------------------------------------------------------
! General order element implementation
!------------------------------------------------------------------------------
subroutine jacobian_gen()
real(r8) :: pt(3),getmp(2),gftmp(3)
integer(i4) :: k,l,ed,etmp(2),dof
! Vertex nodes
do k=1,4
  jfull(:,k)=jfull(:,k) + dlag_1d(self%order+1,f(k),self%xnodes,self%order+1) &
    *self%r(:,self%lc(k,cell))
end do
! Edge nodes
do k=1,6
  ed=self%lce(k,cell)
  do dof=1,self%ho_info%nep
    IF(ed<0)THEN
      pt=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
    ELSE
      pt=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
    END IF
    getmp=dlag_1d_bary(dof,f(tet_ed(:,k)),self%xnodes,self%order+1)
    do l=1,2
      jfull(:,tet_ed(l,k)) = jfull(:,tet_ed(l,k)) + getmp(l)*pt
    end do
  end do
end do
! Face nodes
IF(self%ho_info%nfp>0)THEN
  do k=1,4
    ! Add loop and orient if add higher order
    pt=self%ho_info%r(:,self%ho_info%lfp(1,ABS(self%lcf(k,cell))))
    gftmp=dlag_2d_bary((/1,1/),f(tet_fc(:,k)),self%xnodes,self%order+1)
    do l=1,3
      jfull(:,tet_fc(l,k)) = jfull(:,tet_fc(l,k)) + gftmp(l)*pt
    end do
  end do
END IF
end subroutine jacobian_gen
end subroutine tetmesh_jacobian
!------------------------------------------------------------------------------
!> Linear implementation of @tetmesh_jacobian
!------------------------------------------------------------------------------
subroutine tetmesh_jacl(self,cell,gop,j)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
real(r8) :: jfull(3,4),A(3,3),C(3,3)
integer(i4) :: k
DEBUG_STACK_PUSH
! Get node points
do k=1,4
  jfull(:,k)=self%r(:,self%lc(k,cell))
end do
call tetmesh_jacinv(jfull,gop,j)
DEBUG_STACK_POP
end subroutine tetmesh_jacl
!------------------------------------------------------------------------------
!> Compute the spatial hessian matrices for a given cell at a given logical position
!------------------------------------------------------------------------------
subroutine tetmesh_hessian(self,cell,f,g2op,K)
class(oft_tetmesh), intent(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
REAL(r8), INTENT(out) :: g2op(:,:) !< Second order Jacobian matrix
!! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
REAL(r8), INTENT(out) :: K(:,:) !< Gradient correction matrix
!! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [10,3]
real(r8) :: jfull(3,4),getmp(2),gftmp(3),d2etmp(3),d2ftmp(6),pt(3)
integer(i4) :: j,m,l,etmp(2),dof,ed
integer(i4), parameter :: pmap(4)=(/1,5,8,10/)
integer(i4), parameter :: emap(6)=(/4,7,9,6,3,2/)
integer(i4), parameter :: fmap(4,4)=RESHAPE((/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/),(/4,4/))
if(self%order>3)call oft_abort('Invalid mesh order','tetmesh_hessian',__FILE__)
jfull=0.d0
K=0.d0
IF(cell_is_curved(self, cell))THEN
  do m=1,4 ! Get corner nodes
    jfull(:,m)=jfull(:,m) + dlag_1d(self%order+1,f(m),self%xnodes,self%order+1)*self%r(:,self%lc(m,cell))
    !
    K(pmap(m),:) = K(pmap(m),:) + d2lag_1d(self%order+1,f(m),self%xnodes,self%order+1)*self%r(:,self%lc(m,cell))
  end do
  do m=1,6 ! Get edge nodes
    ed=self%lce(m,cell)
    do dof=1,self%ho_info%nep
      IF(ed<0)THEN
        pt=self%ho_info%r(:,self%ho_info%lep(self%ho_info%nep+1-dof,ABS(ed)))
      ELSE
        pt=self%ho_info%r(:,self%ho_info%lep(dof,ABS(ed)))
      END IF
      getmp=dlag_1d_bary(dof,f(tet_ed(:,m)),self%xnodes,self%order+1)
      do l=1,2
        jfull(:,tet_ed(l,m))=jfull(:,tet_ed(l,m))+getmp(l)*pt
      end do
      !
      d2etmp=d2lag_1d_bary(dof,f(tet_ed(:,m)),self%xnodes,self%order+1)
      K(pmap(tet_ed(1,m)),:)=K(pmap(tet_ed(1,m)),:)+d2etmp(1)*pt
      K(emap(m),:)=K(emap(m),:)+d2etmp(2)*pt
      K(pmap(tet_ed(2,m)),:)=K(pmap(tet_ed(2,m)),:)+d2etmp(3)*pt
    end do
  end do
  IF(self%ho_info%nfp>0)THEN
    do m=1,4 ! Get face nodes
      pt=self%ho_info%r(:,self%ho_info%lfp(1,ABS(self%lcf(m,cell))))
      gftmp=dlag_2d_bary((/1,1/),f(tet_fc(:,m)),self%xnodes,self%order+1)
      do l=1,3
        jfull(:,tet_fc(l,m))=jfull(:,tet_fc(l,m))+gftmp(l)*pt
      end do
      !
      d2ftmp=d2lag_2d_bary((/1,1/),f(tet_fc(:,m)),self%xnodes,self%order+1)
      K(pmap(tet_fc(1,m)),:)=K(pmap(tet_fc(1,m)),:)+d2ftmp(1)*pt
      K(fmap(tet_fc(1,m),tet_fc(2,m)),:)=K(fmap(tet_fc(1,m),tet_fc(2,m)),:)+d2ftmp(2)*pt
      K(fmap(tet_fc(1,m),tet_fc(3,m)),:)=K(fmap(tet_fc(1,m),tet_fc(3,m)),:)+d2ftmp(3)*pt
      K(pmap(tet_fc(2,m)),:)=K(pmap(tet_fc(2,m)),:)+d2ftmp(4)*pt
      K(fmap(tet_fc(2,m),tet_fc(3,m)),:)=K(fmap(tet_fc(2,m),tet_fc(3,m)),:)+d2ftmp(5)*pt
      K(pmap(tet_fc(3,m)),:)=K(pmap(tet_fc(3,m)),:)+d2ftmp(6)*pt
    end do
  END IF
ELSE
  do m=1,4 ! Get corner nodes
    jfull(:,m)=self%r(:,self%lc(m,cell))
  end do
END IF
CALL tetmesh_g2inv(jfull,g2op)
end subroutine tetmesh_hessian
!------------------------------------------------------------------------------
!> Compute the surface normal vector for a given face on a cell
!!
!! If face is not a global boundary face the function returns with `norm = 0`
!!
!! @note The logical position in the cell must be on the chosen face for this
!! subroutine, else an error will be thrown
!------------------------------------------------------------------------------
subroutine tetmesh_snormal(self,cell,ind,f,norm)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell
integer(i4), intent(in) :: ind !< Index of edge within cell
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: norm(3) !< Unit vector normal to the face [3]
REAL(r8) :: goptmp(3,4),v
DEBUG_STACK_PUSH
!---
if(cell<0.OR.cell>self%nc)call oft_abort('Invalid cell index.','tetmesh_snormal',__FILE__)
if(ind<0.OR.ind>4)call oft_abort('Invalid face index.','tetmesh_snormal',__FILE__)
IF(f(ind)/=0._r8)call oft_abort('Invalid cell position.','tetmesh_snormal',__FILE__)
!---
norm=0._r8
! IF(.NOT.self%global%gbc(i))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
! IF(.NOT.self%global%gbf(self%lcf(ind,i)))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
!---
CALL tetmesh_jacobian(self,cell,f,goptmp,v)
norm=-goptmp(:,ind)/sqrt(sum(goptmp(:,ind)**2))
DEBUG_STACK_POP
end subroutine tetmesh_snormal
!------------------------------------------------------------------------------
!> Compute the curve tangent vector for a given edge on a cell
!!
!! If edge is not a global boundary edge the function returns with `tang = 0`
!!
!! @note The logical position in the cell must be on the chosen edge for this
!! subroutine to return a meaningful result
!------------------------------------------------------------------------------
subroutine tetmesh_ctang(self,cell,ind,f,tang)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell
integer(i4), intent(in) :: ind !< Index of edge within cell
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: tang(3) !< Unit vector tangent to the edge [3]
INTEGER(i4) :: j,k
REAL(r8) :: goptmp(3,4),v,norm(3,2),e(3)
DEBUG_STACK_PUSH
!---
IF(cell<0.OR.cell>self%nc)CALL oft_abort('Invalid cell index.','tetmesh_ctang',__FILE__)
IF(ind<0.OR.ind>6)CALL oft_abort('Invalid edge index.','tetmesh_ctang',__FILE__)
! IF(ANY(f(tet_ed(:,ind))/=0._r8))CALL oft_abort('Invalid cell position.','tetmesh_ctang',__FILE__)
! !---
! norm=0._r8
! if(.NOT.self%global%gbc(i))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
! if(.NOT.self%global%gbe(ABS(self%lce(ind,i))))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
!---
k=1
CALL tetmesh_jacobian(self,cell,f,goptmp,v)
DO j=1,4
  IF(ALL(tet_ed(:,ind)/=j))THEN
    norm(:,k)=goptmp(:,j)
    k=k+1
  END IF
END DO
e=self%r(:,self%lc(tet_ed(2,ind),cell))-self%r(:,self%lc(tet_ed(1,ind),cell))
tang=cross_product(norm(:,1),norm(:,2))
tang=tang/SQRT(SUM(tang**2))
tang=tang*SIGN(1_i4,self%lce(ind,cell))*SIGN(1._r8,DOT_PRODUCT(tang,e))
DEBUG_STACK_POP
end subroutine tetmesh_ctang
!------------------------------------------------------------------------------
!> Get mapping between boundary and volume logical coordinates
!------------------------------------------------------------------------------
subroutine tetmesh_get_surf_map(self,face,cell,lmap)
class(oft_tetmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: face !< Index of face on boundary mesh
integer(i4), intent(out) :: cell !< Cell containing face
integer(i4), intent(out) :: lmap(3) !< Coordinate mapping
integer(i4) :: j,ftmp,pmap(3)
!---Find parent cell and face index in cell
ftmp=ABS(self%bmesh%parent%lf(face))
cell=self%lfc(1,ftmp)
DO j=1,4
  IF(ABS(self%lcf(j,cell))==ftmp)EXIT
END DO
!---Map surface logical coordinates to cell logical
pmap=(/1,2,3/)
CALL orient_listn(self%bmesh%lco(face),pmap,3_i4)
CALL orient_listn_inv(self%lcfo(j,cell),pmap,3_i4)
lmap(pmap)=self%cell_fc(:,j)
end subroutine tetmesh_get_surf_map
!------------------------------------------------------------------------------
!> Map between surface and volume logical coordinates
!------------------------------------------------------------------------------
subroutine tetmesh_surf_to_vol(self,fsurf,lmap,fvol)
CLASS(oft_tetmesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: fsurf(:) !< Surface coordinates [3]
INTEGER(i4), INTENT(in) :: lmap(3) !< Coordinate mapping
REAL(r8), INTENT(out) :: fvol(:) !< Volume coordinates [4]
fvol=0.d0
fvol(lmap)=fsurf(1:3)
end subroutine tetmesh_surf_to_vol
!------------------------------------------------------------------------------
!> Invert a 3x3 matrix
!------------------------------------------------------------------------------
subroutine tetmesh_jacinv(jfull,gop,jac)
real(r8), intent(in) :: jfull(3,4) !< Matrix to invert
real(r8), intent(out) :: gop(3,4) !< \f$ A^{-1} \f$
real(r8), intent(out) :: jac !< |A|
real(r8) :: t1,t2,t3,A(3,3),C(3,3)
DEBUG_STACK_PUSH
!---
A(1,:)=jfull(:,2)-jfull(:,1)
A(2,:)=jfull(:,3)-jfull(:,1)
A(3,:)=jfull(:,4)-jfull(:,1)
!---Compute resusables
t1=A(2,2)*A(3,3)-A(3,2)*A(2,3)
t2=A(3,2)*A(1,3)-A(1,2)*A(3,3)
t3=A(1,2)*A(2,3)-A(2,2)*A(1,3)
!---Compute Det(A)
jac=A(1,1)*t1 + A(2,1)*t2 + A(3,1)*t3
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
C=C/jac
!---
gop(:,2:4)=C
gop(:,1)=-(gop(:,2)+gop(:,3)+gop(:,4))
jac=jac/6.d0
DEBUG_STACK_POP
end subroutine tetmesh_jacinv
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine tetmesh_g2inv(jfull,g2op)
real(r8), intent(in) :: jfull(3,4)
real(r8), intent(out) :: g2op(6,10)
REAL(r8) :: jac(3,3),A(6,6),work(12)
INTEGER(i4) :: info,ipiv(6)
DEBUG_STACK_PUSH
!---
jac(1,:)=jfull(:,2)-jfull(:,1)
jac(2,:)=jfull(:,3)-jfull(:,1)
jac(3,:)=jfull(:,4)-jfull(:,1)
!---Get 2nd order mapping
!---Row 1 (x_1/l_j)*(x_1/l_k)
A(1,1)=jac(1,1)**2
A(1,2)=2*jac(1,1)*jac(1,2)
A(1,3)=2*jac(1,1)*jac(1,3)
A(1,4)=jac(1,2)**2
A(1,5)=2*jac(1,2)*jac(1,3)
A(1,6)=jac(1,3)**2
!---Row 2 (x_1/l_j)*(x_2/l_k)
A(2,1)=jac(1,1)*jac(2,1)
A(2,2)=jac(1,1)*jac(2,2)+jac(2,1)*jac(1,2)
A(2,3)=jac(1,1)*jac(2,3)+jac(2,1)*jac(1,3)
A(2,4)=jac(1,2)*jac(2,2)
A(2,5)=jac(1,2)*jac(2,3)+jac(2,2)*jac(1,3)
A(2,6)=jac(1,3)*jac(2,3)
!---Row 3 (x_1/l_j)*(x_3/l_k)
A(3,1)=jac(1,1)*jac(3,1)
A(3,2)=jac(1,1)*jac(3,2)+jac(3,1)*jac(1,2)
A(3,3)=jac(1,1)*jac(3,3)+jac(3,1)*jac(1,3)
A(3,4)=jac(1,2)*jac(3,2)
A(3,5)=jac(1,2)*jac(3,3)+jac(3,2)*jac(1,3)
A(3,6)=jac(1,3)*jac(3,3)
!---Row 4 (x_2/l_j)*(x_2/l_k)
A(4,1)=jac(2,1)**2
A(4,2)=2*jac(2,1)*jac(2,2)
A(4,3)=2*jac(2,1)*jac(2,3)
A(4,4)=jac(2,2)**2
A(4,5)=2*jac(2,2)*jac(2,3)
A(4,6)=jac(2,3)**2
!---Row 5 (x_2/l_j)*(x_3/l_k)
A(5,1)=jac(2,1)*jac(3,1)
A(5,2)=jac(2,1)*jac(3,2)+jac(3,1)*jac(2,2)
A(5,3)=jac(2,1)*jac(3,3)+jac(3,1)*jac(2,3)
A(5,4)=jac(2,2)*jac(3,2)
A(5,5)=jac(2,2)*jac(3,3)+jac(3,2)*jac(2,3)
A(5,6)=jac(2,3)*jac(3,3)
!---Row 6 (x_3/l_j)*(x_3/l_k)
A(6,1)=jac(3,1)**2
A(6,2)=2*jac(3,1)*jac(3,2)
A(6,3)=2*jac(3,1)*jac(3,3)
A(6,4)=jac(3,2)**2
A(6,5)=2*jac(3,2)*jac(3,3)
A(6,6)=jac(3,3)**2
!---Invert
CALL dgetrf(6,6,A,6,ipiv,info)
IF(info/=0)WRITE(*,*)info
CALL dgetri(6,A,6,ipiv,work,12,info)
IF(info/=0)WRITE(*,*)info
!---Map back to g2op
g2op(:,5)=A(:,1)
g2op(:,6)=A(:,2)
g2op(:,7)=A(:,3)
g2op(:,8)=A(:,4)
g2op(:,9)=A(:,5)
g2op(:,10)=A(:,6)
!---Scatter to redundant dimension
g2op(:,2)=-2*A(:,1)-A(:,2)-A(:,3)
g2op(:,3)=-A(:,2)-2*A(:,4)-A(:,5)
g2op(:,4)=-A(:,3)-A(:,5)-2*A(:,6)
g2op(:,1)=SUM(A,DIM=2)
DEBUG_STACK_POP
end subroutine tetmesh_g2inv
END MODULE oft_tetmesh_type
