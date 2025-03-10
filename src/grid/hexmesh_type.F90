!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_hexmesh_type.F90
!
!> Hexhedral mesh definitions
!!
!! @author Chris Hansen
!! @date October 2018
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_hexmesh_type
USE oft_base
USE oft_gauss_quadrature, ONLY: set_quad_3d, oft_quad_type
USE oft_lag_poly
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> Hexahedral volume mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_mesh) :: oft_hexmesh
  INTEGER(i4) :: inodesp(3,8) !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:) :: inodesc => NULL() !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodese => NULL() !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodesf => NULL() !< Needs docs
  REAL(r8), POINTER, DIMENSION(:) :: xnodes => NULL() !< Needs docs
CONTAINS
  PROCEDURE :: setup => hexmesh_setup
  PROCEDURE :: set_order => hexmesh_set_order
  PROCEDURE :: invert_cell => hexmesh_invert_cell
  PROCEDURE :: log2phys => hexmesh_log2phys
  PROCEDURE :: phys2log => hexmesh_phys2log
  PROCEDURE :: jacobian => hexmesh_jacobian
  PROCEDURE :: hessian => hexmesh_hessian
  PROCEDURE :: snormal => hexmesh_snormal
  PROCEDURE :: ctang => hexmesh_ctang
  PROCEDURE :: get_surf_map => hexmesh_get_surf_map
  PROCEDURE :: surf_to_vol => hexmesh_surf_to_vol
  PROCEDURE :: vlog => hexmesh_vlog
  PROCEDURE :: in_cell => hexmesh_in_cell
  PROCEDURE :: quad_rule => hexmesh_quad_rule
  PROCEDURE :: tessellate => hexmesh_tessellate
  PROCEDURE :: tessellated_sizes => hexmesh_tessellated_sizes
END TYPE oft_hexmesh
!---
INTEGER(i4), PARAMETER :: hex_ed(2,12)=RESHAPE((/ & !< Hexahedron edge list
  1,2, 2,3, 3,4, 4,1, &
  1,5, 2,6, 3,7, 4,8, &
  5,6, 6,7, 7,8, 8,5/),(/2,12/))
INTEGER(i4), PARAMETER :: hex_fc(4,6)=RESHAPE((/ & !< Hexahedron face list
  1,2,3,4, 1,5,6,2, 2,6,7,3, &
  3,7,8,4, 1,4,8,5, 5,8,7,6/),(/4,6/))
INTEGER(i4), PARAMETER :: hex_fe(4,6)=RESHAPE((/ & !< Hexahedron face edge list
  1,2,3,4, 5,9,-6,-1, 6,10,-7,-2, &
  7,11,-8,-3, -4,8,12,-5, -12,-11,-10,-9/),(/4,6/))
INTEGER(i4), PARAMETER :: quad_ed(2,4)=RESHAPE((/1,2, 2,3, 3,4, 4,1/),(/2,4/)) !< Quad edge list
!
INTEGER(i4), PARAMETER :: hex_bary_map(6) = (/-3,-2,1, 2,-1,3/)
INTEGER(i4), PARAMETER :: hex_bary_pfcoords(3,8) = RESHAPE((/ & ! Do not vary
  1,2,5, 1,2,3, 1,3,4, 1,4,5, 2,5,6, 2,3,6, 3,4,6, 4,5,6/), (/3,8/))
INTEGER(i4), PARAMETER :: hex_bary_efcoords(2,12) = RESHAPE((/ & ! Do not vary
  1,2, 1,3, 1,4, 1,5, 2,5, 2,3, 3,4, 4,5, 2,6, 3,6, 4,6, 5,6/), (/2,12/))
INTEGER(i4), PARAMETER :: hex_bary_ecoords(2,12) = RESHAPE((/ & ! Do vary
  5,3, 2,4, 3,5, 4,2, &
  1,6, 1,6, 1,6, 1,6, &
  5,3, 2,4, 3,5, 4,2/), (/2,12/))
INTEGER(i4), PARAMETER :: hex_bary_fcoords(4,6) = RESHAPE((/ & ! Do vary
  5,2,3,4, 1,5,6,3, 1,2,6,4, &
  1,3,6,5, 2,1,4,6, 2,5,4,3/), (/4,6/))
!
INTEGER(i4), PARAMETER :: inodesp_base(3,8) = RESHAPE((/ &
  0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1/),(/3,8/))
INTEGER(i4), PARAMETER :: inodes1p(3,8) = RESHAPE((/ &
  1,1,1, 2,1,1, 2,2,1, 1,2,1, 1,1,2, 2,1,2, 2,2,2, 1,2,2/),(/3,8/))
INTEGER(i4), PARAMETER :: inodes2p(3,8) = RESHAPE((/ &
  1,1,1, 3,1,1, 3,3,1, 1,3,1, 1,1,3, 3,1,3, 3,3,3, 1,3,3/),(/3,8/))
INTEGER(i4), PARAMETER :: inodes2e(3,12) = RESHAPE((/ &
  2,1,1, 3,2,1, 2,3,1, 1,2,1, 1,1,2, 3,1,2, 3,3,2, 1,3,2, &
  2,1,3, 3,2,3, 2,3,3, 1,2,3/),(/3,12/))
INTEGER(i4), PARAMETER :: inodes2f(3,6) = RESHAPE((/ &
  2,2,1, 2,1,2, 3,2,2, 2,3,2, 1,2,2, 2,2,3/),(/3,6/))
!
INTEGER(i4), PRIVATE, PARAMETER :: ho_find_nsteps=100 !< Maximum number of steps during high order find_cell
REAL(r8), PRIVATE, PARAMETER :: ho_find_du=1.d-6 !< Step size used for jacobian eval during high order find_cell
REAL(r8), PRIVATE, PARAMETER :: ho_find_tol=1.d-6 !< Convergence tolerance for high order find_cell
CLASS(oft_hexmesh), PRIVATE, POINTER :: active_mesh => NULL() !< Active mesh for high order find_cell
INTEGER(i4), PRIVATE :: active_cell = 0 !< Active cell for high order find_cell
REAL(r8), PRIVATE :: active_pt(3) = 0.d0 !< Active point for high order find_cell
!$omp threadprivate(active_mesh,active_cell,active_pt)
!
PUBLIC hex_3d_grid, hex_grid_forient, hex_get_bary, hex_get_bary_gop, hex_get_bary_cgop
PUBLIC hex_bary_pfcoords, hex_bary_efcoords, hex_bary_ecoords, hex_bary_fcoords
CONTAINS
!---------------------------------------------------------------------------
!> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
!---------------------------------------------------------------------------
SUBROUTINE hexmesh_setup(self,cad_type)
CLASS(oft_hexmesh), INTENT(inout) :: self!< Mesh object
INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
self%type=3
self%cad_type=cad_type
! self%cell_geom=(/8,12,6/)
self%face_np=4
self%cell_np=8
self%cell_ne=12
self%cell_nf=6
ALLOCATE(self%cell_ed(2,12),self%cell_fc(4,6),self%cell_fe(4,6))
self%cell_ed=hex_ed
self%cell_fc=hex_fc
self%cell_fe=ABS(hex_fe)
ALLOCATE(self%face_ed(2,4))
self%face_ed=quad_ed  
CALL self%set_order(1)
END SUBROUTINE hexmesh_setup
!---------------------------------------------------------------------------
!> Set maximum order of spatial mapping
!---------------------------------------------------------------------------
SUBROUTINE hexmesh_set_order(self,order)
CLASS(oft_hexmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
INTEGER(i4) :: i,j,k,np,inc1,ind1,inc2,ind2,inode(3),counts(3)
IF(order<1.OR.order>2)CALL oft_abort("Invalid grid order", "hexmesh_set_order", __FILE__)
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
IF(ASSOCIATED(self%inodese))DEALLOCATE(self%inodese,self%inodesf,self%inodesc)
self%order=order
!---Setup high order gemetry
select case(self%order)
case(1)
  counts=0
case(2)
  self%ho_info%nep=1
  self%ho_info%nfp=1
  self%ho_info%ncp=1
  counts=(/self%ne,self%nf,self%nc/)
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
!---Points
CALL hex_3d_grid(order,self%xnodes,self%inodesp,self%inodese,self%inodesf,self%inodesc)
END SUBROUTINE hexmesh_set_order
!---------------------------------------------------------------------------
!> Map from orthogonal logical coordinates to barycentric logical coordinates
!---------------------------------------------------------------------------
PURE FUNCTION hex_get_bary(flog) RESULT(fbary)
REAL(r8), INTENT(in) :: flog(:) !< Position in orthogonal logical coordinates [3]
REAL(r8) :: fbary(6) !< Position in barycentric logical coordinates [6]
INTEGER(i4) :: i
DO i=1,6
  IF(hex_bary_map(i)<0)THEN
    fbary(i)=1.d0-flog(ABS(hex_bary_map(i)))
  ELSE
    fbary(i)=flog(ABS(hex_bary_map(i)))
  END IF
END DO
END FUNCTION  hex_get_bary
!---------------------------------------------------------------------------
!> Map gradients from orthogonal logical coordinates to barycentric logical coordinates
!---------------------------------------------------------------------------
PURE FUNCTION hex_get_bary_gop(glog) RESULT(gbary)
REAL(r8), INTENT(in) :: glog(:,:) !< Gradients in orthogonal logical coordinates [3,3]
REAL(r8) :: gbary(3,6) !< Gradients in barycentric logical coordinates [3,6]
INTEGER(i4) :: i
DO i=1,6
  IF(hex_bary_map(i)<0)THEN
    gbary(:,i)=-glog(:,ABS(hex_bary_map(i)))
  ELSE
    gbary(:,i)=glog(:,ABS(hex_bary_map(i)))
  END IF
END DO
END FUNCTION  hex_get_bary_gop
!---------------------------------------------------------------------------
!> Map gradient cross-products from orthogonal logical coordinates to barycentric logical coordinates
!---------------------------------------------------------------------------
PURE FUNCTION hex_get_bary_cgop(cglog,i,j) RESULT(cgbary)
REAL(r8), INTENT(in) :: cglog(3,3) !< Crossed-gradients in orthogonal logical coordinates [3,3]
INTEGER(i4), INTENT(in) :: i !< First barycentric coordinate
INTEGER(i4), INTENT(in) :: j !< Second barycentric coordinate
REAL(r8) :: cgbary(3) !< Crossed-gradient between barycentric coordinates `i` and `j` [3]
INTEGER(i4) :: ii,jj,k
INTEGER(i4), PARAMETER :: bary_cmap(3,3) = &
  RESHAPE((/0,-1,-2, 1,0,-3, 2,3,0/), (/3,3/))
ii=hex_bary_map(i)
jj=hex_bary_map(j)
k=bary_cmap(ABS(ii), ABS(jj))
IF(k==0)THEN
  cgbary=0.d0
ELSE
  cgbary=cglog(:, ABS(k))*SIGN(1, ii*jj*k)
END IF
END FUNCTION  hex_get_bary_cgop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
SUBROUTINE hex_3d_grid(order,xnodes,inodesp,inodese,inodesf,inodesc)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4), INTENT(out) :: inodesp(3,8)
INTEGER(i4), POINTER, DIMENSION(:,:,:), INTENT(out) :: inodese
INTEGER(i4), POINTER, DIMENSION(:,:,:), INTENT(out) :: inodesf
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: inodesc
REAL(r8), POINTER, DIMENSION(:), INTENT(out) :: xnodes
INTEGER(i4) :: i,j,k,np,inc1,ind1,inc2,ind2,inode(3),counts(3)
!---Points
inodesp = inodesp_base*order + 1
ALLOCATE(xnodes(order+1))
DO i=1,order+1
  xnodes(i)=REAL(i-1,8)/REAL(order,8)
END DO
!---Edges
np = order - 1
IF(np>0)THEN
  ALLOCATE(inodese(3,np,12))
  DO i=1,12
    ind1 = MAXLOC(ABS(inodesp(:,hex_ed(2,i))-inodesp(:,hex_ed(1,i))), DIM=1)
    inc1 = SIGN(1_i4, inodesp(ind1,hex_ed(2,i))-inodesp(ind1,hex_ed(1,i)))
    inode = inodesp(:,hex_ed(1,i))
    DO j=1,order-1
      inode(ind1) = inode(ind1) + inc1
      inodese(:,j,i) = inode
      ! WRITE(*,*)i,j,inode
    END DO
  END DO
  !---Faces
  ALLOCATE(inodesf(3,np**2,6))
  ! WRITE(*,*)
  DO i=1,6
    ! Direction 1
    ind1 = MAXLOC(ABS(inodesp(:,hex_fc(2,i))-inodesp(:,hex_fc(1,i))), DIM=1)
    inc1 = SIGN(1_i4, inodesp(ind1,hex_fc(2,i))-inodesp(ind1,hex_fc(1,i)))
    ! Direction 2
    ind2 = MAXLOC(ABS(inodesp(:,hex_fc(3,i))-inodesp(:,hex_fc(2,i))), DIM=1)
    inc2 = SIGN(1_i4, inodesp(ind2,hex_fc(3,i))-inodesp(ind2,hex_fc(2,i)))
    inode = inodesp(:,hex_fc(1,i))
    DO j=1,np
      inode(ind1) = inode(ind1) + inc1
      DO k=1,np
        inode(ind2) = inode(ind2) + inc2
        inodesf(:,(j-1)*np+k,i) = inode
        ! WRITE(*,*)i,j,k,inode
      END DO
      inode(ind2) = inodesp(ind2,hex_fc(1,i))
    END DO
  END DO
  !---Cell
  ALLOCATE(inodesc(3,np**3))
  ! WRITE(*,*)
  DO i=1,np
    DO j=1,np
      DO k=1,np
        inodesc(:,(i-1)*np*np+(j-1)*np+k) = (/1+i,1+j,1+k/)
        ! WRITE(*,*)i,j,k,(/1+i,1+j,1+k/)
      END DO
    END DO
  END DO
END IF
END SUBROUTINE hex_3d_grid
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
SUBROUTINE hex_grid_forient(oflag,order,finds)
INTEGER(i4), INTENT(in) :: oflag,order
INTEGER(i4), DIMENSION(:), INTENT(inout) :: finds
INTEGER(i4), PARAMETER :: inodesp(2,4) = RESHAPE((/0,0, 1,0, 1,1, 0,1/), (/2,4/))
INTEGER(i4) :: i,j,inc1(2),inc2(2),np,ptmp(4),inode(2),p1(2),p2(2),p3(2)
REAL(r8) :: o8
!
ptmp=(/1,2,3,4/)
CALL orient_listn(oflag, ptmp, 4_i4)
! CALL orient_listn_inv(oflag, ptmp, 4_i4)
p1=inodesp(:,ptmp(1))*order+1
p2=inodesp(:,ptmp(2))*order+1
p3=inodesp(:,ptmp(3))*order+1
np = order-1
! Direction 1
inc1 = (p2-p1)
! Direction 2
inc2 = (p3-p2)
!
! WRITE(*,*)p1
! WRITE(*,*)p2
! WRITE(*,*)p3
! WRITE(*,*)inc1
! WRITE(*,*)inc2
DO i=1,np
  DO j=1,np
    inode = p1 + inc1*i/(order) + inc2*j/(order) - 1
    ! WRITE(*,*)i,j,inode
    finds((i-1)*np+j) = (inode(1)-1)*np + inode(2)
  END DO
END DO
END SUBROUTINE hex_grid_forient
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
SUBROUTINE hex_grid_forient_inv(oflag,order,finds)
INTEGER(i4), INTENT(in) :: oflag,order
INTEGER(i4), DIMENSION(:), INTENT(inout) :: finds
INTEGER(i4), PARAMETER :: inodesp(2,4) = RESHAPE((/0,0, 1,0, 1,1, 0,1/), (/2,4/))
INTEGER(i4) :: i,j,inc1(2),inc2(2),np,ptmp(4),inode(2),p1(2),p2(2),p3(2)
REAL(r8) :: o8
!
ptmp=(/1,2,3,4/)
! CALL orient_listn(oflag, ptmp, 4_i4)
CALL orient_listn_inv(oflag, ptmp, 4_i4)
p1=inodesp(:,ptmp(1))*order+1
p2=inodesp(:,ptmp(2))*order+1
p3=inodesp(:,ptmp(3))*order+1
np = order-1
! Direction 1
inc1 = (p2-p1)
! Direction 2
inc2 = (p3-p2)
!
! WRITE(*,*)p1
! WRITE(*,*)p2
! WRITE(*,*)p3
! WRITE(*,*)inc1
! WRITE(*,*)inc2
DO i=1,np
  DO j=1,np
    inode = p1 + inc1*i/(order) + inc2*j/(order) - 1
    ! WRITE(*,*)i,j,inode
    finds((i-1)*np+j) = (inode(1)-1)*np + inode(2)
  END DO
END DO
END SUBROUTINE hex_grid_forient_inv
!---------------------------------------------------------------------------
!> Turn cell "inside out", used to ensure consistent orientations
!---------------------------------------------------------------------------
SUBROUTINE hexmesh_invert_cell(self,cell)
CLASS(oft_hexmesh), INTENT(inout) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: cell !< Index of cell to invert
INTEGER(i4) :: lctmp(8)
lctmp=self%lc(:,cell)
self%lc(1:4,cell)=lctmp(5:8)
self%lc(5:8,cell)=lctmp(1:4)
END SUBROUTINE hexmesh_invert_cell
!------------------------------------------------------------------------------
!> Retrieve suitable quadrature rule for mesh with given order
!------------------------------------------------------------------------------
SUBROUTINE hexmesh_quad_rule(self,order,quad_rule)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
CALL set_quad_3d(quad_rule, order)
END SUBROUTINE hexmesh_quad_rule
!------------------------------------------------------------------------------
!> Get position in logical space of vertex `i`
!------------------------------------------------------------------------------
SUBROUTINE hexmesh_vlog(self,i,f)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4), INTENT(in) :: i !< Vertex to locate
REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
f(1:3)=self%xnodes(self%inodesp(:,i))
END SUBROUTINE hexmesh_vlog
!------------------------------------------------------------------------------
!> Test if logical position lies within the base cell
!!
!! @returns Position `f` is inside the base cell?
!------------------------------------------------------------------------------
FUNCTION hexmesh_in_cell(self,f,tol) RESULT(eface)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
REAL(r8), INTENT(in) :: tol !< Tolerance for test
INTEGER(i4) :: eface
REAL(r8) :: fmin,fmax,fbary(6)
fmin=MINVAL(f(1:3))
fmax=MAXVAL(f(1:3))
IF((fmax<=1.d0+tol).AND.(fmin>=-tol))THEN
  eface=0
ELSE
  fbary=hex_get_bary(f(1:3))
  eface=MINLOC(fbary, DIM=1)
END IF
END FUNCTION hexmesh_in_cell
!------------------------------------------------------------------------------
!> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
!!
!! @note The maximum tessellation order currently supported is 4
!! (may be lower for certain mesh types).
!!
!! @warning Cell lists are returned with zero based indexing
!------------------------------------------------------------------------------
SUBROUTINE hexmesh_tessellate(self,rtmp,lctmp,order)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh to tessellate
REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,:]
INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [8,:]
INTEGER(i4), INTENT(in) :: order !< Tessellation order
INTEGER(i4) :: i,j,k,cell,np,inodesp(3,8),inodespf(3,4),i1,i2,inc
INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodese,inodesf,pmap
INTEGER(i4), POINTER, DIMENSION(:,:) :: inodesc
INTEGER(i4), POINTER, DIMENSION(:) :: finds
REAL(r8), POINTER, DIMENSION(:) :: xnodes
LOGICAL, ALLOCATABLE, DIMENSION(:) :: pflag
DEBUG_STACK_PUSH
!---Call desired driver
IF(order==1)THEN
  ALLOCATE(rtmp(3,self%np))
  rtmp=self%r
  ALLOCATE(lctmp(8,self%nc))
  lctmp=self%lc-1
ELSE
  CALL hex_3d_grid(order,xnodes,inodesp,inodese,inodesf,inodesc)
  np=order-1
  ALLOCATE(rtmp(3,self%np + self%ne*np + self%nf*np**2 + self%nc*np**3))
  ALLOCATE(pflag( self%np + self%ne*np + self%nf*np**2 + self%nc*np**3))
  ALLOCATE(pmap(order+1,order+1,order+1),finds(np**2))
  ALLOCATE(lctmp(8, self%nc*order**3))
  pflag=.FALSE.
  rtmp(:,1:self%np)=self%r
  pflag(1:self%np)=.TRUE.
  DO cell=1,self%nc
    !---Point
    DO i=1,8
      j=self%lc(i,cell)
      pmap(inodesp(1,i), inodesp(2,i), inodesp(3,i))=j
    END DO
    !---Edges
    DO i=1,12
      DO k=1,np
        i1=k
        IF(self%lce(i,cell)<0)i1=np+1-k
        j=(ABS(self%lce(i,cell))-1)*np + self%np + i1
        pmap(inodese(1,k,i), inodese(2,k,i), inodese(3,k,i))=j
        IF(.NOT.pflag(j))THEN
          rtmp(:,j)=self%log2phys(cell, xnodes(inodese(:,k,i)))
          pflag(j)=.TRUE.
        END IF
      END DO
    END DO
    !---Faces
    DO i=1,6
      CALL hex_grid_forient_inv(self%lcfo(i,cell),order,finds)
      DO k=1,np**2
        j=(ABS(self%lcf(i,cell))-1)*np**2 + self%ne*np + self%np + finds(k)
        pmap(inodesf(1,k,i), inodesf(2,k,i), inodesf(3,k,i))=j
        IF(.NOT.pflag(j))THEN
          rtmp(:,j)=self%log2phys(cell, xnodes(inodesf(:,k,i)))
          pflag(j)=.TRUE.
        END IF
      END DO
    END DO
    !---Cell interiors
    DO k=1,np**3
      j=(cell-1)*np**3 + self%nf*np**2 + self%ne*np + self%np + k
      pmap(inodesc(1,k), inodesc(2,k), inodesc(3,k))=j
      rtmp(:,j)=self%log2phys(cell, xnodes(inodesc(:,k)))
    END DO
    !---Create cells
    DO i=1,order
      DO j=1,order
        DO k=1,order
          lctmp(:,(cell-1)*order**3 + (i-1)*order**2 + (j-1)*order + k) = &
          (/pmap(i,j,k),pmap(i+1,j,k),pmap(i+1,j+1,k),pmap(i,j+1,k), &
            pmap(i,j,k+1),pmap(i+1,j,k+1),pmap(i+1,j+1,k+1),pmap(i,j+1,k+1)/)-1
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(xnodes,inodese,inodesf,inodesc)
  DEALLOCATE(pmap,finds)
  !CALL oft_abort('Invalid tessellation order','hexmesh_tessellate',__FILE__)
END IF
DEBUG_STACK_POP
END SUBROUTINE hexmesh_tessellate
!---------------------------------------------------------------------------
!> Get sizes of arrays returned by @ref tetmesh_tessellate
!---------------------------------------------------------------------------
FUNCTION hexmesh_tessellated_sizes(self) result(sizes)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh object
INTEGER(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
sizes(1)=self%np + self%ne*(self%tess_order-1) + self%nf*((self%tess_order-1)**2) &
  + self%nc*((self%tess_order-1)**3)
sizes(2)=self%nc*(self%tess_order**3)
END FUNCTION hexmesh_tessellated_sizes
!------------------------------------------------------------------------------
!> Map from logical to physical coordinates in a given cell
!------------------------------------------------------------------------------
function hexmesh_log2phys(self,cell,f) result(pt)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer, intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8) :: pt(3) !< Physical position [3]
integer(i4) :: k,i,j,inode(3)
DEBUG_STACK_PUSH
if(self%order>2)call oft_abort('Invalid mesh order','hexmesh_log2phys',__FILE__)
pt=0.d0
do i=1,8
  pt=pt+lag_3d(self%inodesp(:,i),f,self%xnodes,self%order+1)*self%r(:,self%lc(i,cell))
end do
IF(self%order>1)THEN
  DO k=1,self%ho_info%nep
    do i=1,12
      j=ABS(self%lce(i,cell))
      pt=pt+lag_3d(self%inodese(:,k,i),f,self%xnodes,self%order+1)*self%ho_info%r(:,self%ho_info%lep(k,j))
    end do
  END DO
  DO k=1,self%ho_info%nfp
    do i=1,6
      j=ABS(self%lcf(i,cell))
      pt=pt+lag_3d(self%inodesf(:,k,i),f,self%xnodes,self%order+1)*self%ho_info%r(:,self%ho_info%lfp(k,j))
    end do
  END DO
  DO k=1,self%ho_info%ncp
    pt=pt+lag_3d(self%inodesc(:,k),f,self%xnodes,self%order+1)*self%ho_info%r(:,self%ho_info%lcp(k,cell))
  END DO
END IF
DEBUG_STACK_POP
end function hexmesh_log2phys
!------------------------------------------------------------------------------
!> Map from physical to logical coordinates in a given cell
!------------------------------------------------------------------------------
subroutine hexmesh_phys2log(self,cell,pt,f)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: pt(3) !< Physical position [3]
real(r8), intent(out) :: f(:) !< Logical coordinates within the cell [4]
f=hexmesh_phys2logho(self,cell,pt)
end subroutine hexmesh_phys2log
!------------------------------------------------------------------------------
!> Implementation of @ref hexmesh_phys2log
!!
!! The MINPACK package is used with step size given by
!! @ref hexmesh_mapping::ho_find_du "ho_find_du". The convergence tolerance is
!! set by the variable @ref hexmesh_mapping::ho_find_tol "ho_find_tol".
!!
!! @note The final location may be outside the cell being searched. This is correct
!! if the point is outside the cell, however it may also indicate a problem in the
!! mapping, most likely due to a badly shaped cell
!------------------------------------------------------------------------------
function hexmesh_phys2logho(self,cell,pt) result(f)
class(oft_hexmesh), target, intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: pt(3) !< Physical position [3]
real(r8) :: f(4) !< Logical coordinates within the cell [4]
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
active_cell=cell
uv=1.d0/2.d0
!---
call lmdif(tm_findcell_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!IF(info>4)WRITE(*,*)'High-order find failed',i,info,nfev
f(1:3)=uv
DEBUG_STACK_POP
end function hexmesh_phys2logho
!---------------------------------------------------------------------------
!> Evalute the error between a logical point and the current active point
!!
!! @note Designed to be used as the error function for minimization in
!! @ref hexmesh_mapping::hexmesh_phys2logho "hexmesh_phys2logho"
!---------------------------------------------------------------------------
subroutine tm_findcell_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions (3)
integer(i4), intent(in) :: n !< Number of parametric dimensions (3)
real(r8), intent(in) :: uv(n) !< Parametric position [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
pt=hexmesh_log2phys(active_mesh,active_cell,uv)
err=active_pt-pt
DEBUG_STACK_POP
end subroutine tm_findcell_error
!------------------------------------------------------------------------------
!> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
!------------------------------------------------------------------------------
subroutine hexmesh_jacobian(self,cell,f,gop,j)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
real(r8), intent(out) :: j !< Jacobian of transformation from logical to physical coordinates
integer(i4) :: i,k,l,m
real(r8) :: A(3,3),pt(3),vtmp(3)
DEBUG_STACK_PUSH
if(self%order>2)call oft_abort('Invalid mesh order','hexmesh_jacobian',__FILE__)
A=0.d0
!---Get node points
do l=1,8
  pt = self%r(:,self%lc(l,cell))
  vtmp=dlag_3d(self%inodesp(:,l),f,self%xnodes,self%order+1)
  DO k=1,3
    A(k,:)=A(k,:)+vtmp(k)*pt
  END DO
end do
IF(self%order>1)THEN
  !---Get edge points
  DO i=1,self%ho_info%nep
    do l=1,12
      m=ABS(self%lce(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lep(1,m))
      vtmp=dlag_3d(self%inodese(:,i,l),f,self%xnodes,self%order+1)
      DO k=1,3
        A(k,:)=A(k,:)+vtmp(k)*pt
      END DO
    end do
  END DO
  !---Get face points
  DO i=1,self%ho_info%nfp
    do l=1,6
      m=ABS(self%lcf(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lfp(i,m))
      vtmp=dlag_3d(self%inodesf(:,i,l),f,self%xnodes,self%order+1)
      DO k=1,3
        A(k,:)=A(k,:)+vtmp(k)*pt
      END DO
    end do
  END DO
  !---Get cell points
  DO i=1,self%ho_info%ncp
    pt=self%ho_info%r(:,self%ho_info%lcp(i,cell))
    vtmp=dlag_3d(self%inodesc(:,i),f,self%xnodes,self%order+1)
    DO k=1,3
      A(k,:)=A(k,:)+vtmp(k)*pt
    END DO
  END DO
END IF
call hexmesh_jacinv(A,gop,j)
DEBUG_STACK_POP
end subroutine hexmesh_jacobian
!---------------------------------------------------------------------------
!> Compute the spatial hessian matrices for a given cell at a given logical position
!---------------------------------------------------------------------------
subroutine hexmesh_hessian(self,cell,f,g2op,K)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell for evaulation
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: g2op(:,:) !< Second order Jacobian matrix
!! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
real(r8), intent(out) :: K(:,:) !< Gradient correction matrix
!! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [10,3]
real(r8) :: ri(3,4),jac(3,3),A(6,6)
integer(i4) :: ii,l,m,j
REAL(r8) :: work(12),pt(3),vtmp(3)
INTEGER(i4) :: info,ipiv(6)
INTEGER(i4), PARAMETER :: dmap(2,6) = RESHAPE((/1,1, 1,2, 1,3, 2,2, 2,3, 3,3/), (/2,6/))
if(self%order>3)call oft_abort('Invalid mesh order','tetmesh_hessian',__FILE__)
jac=0.d0
K=0.d0
!---Get node points
do l=1,8
  pt = self%r(:,self%lc(l,cell))
  vtmp=dlag_3d(self%inodesp(:,l),f,self%xnodes,self%order+1)
  DO j=1,3
    jac(j,:)=jac(j,:)+vtmp(j)*pt
  END DO
  DO j=1,6
    K(j,:) = K(j,:) + d2lag_3d(self%inodesp(:,l),dmap(1,j), &
      dmap(2,j),f,self%xnodes,self%order+1)*pt
  END DO
end do
IF(self%order>1)THEN
  !---Get edge points
  DO ii=1,self%ho_info%nep
    do l=1,12
      m=ABS(self%lce(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lep(ii,m))
      vtmp=dlag_3d(self%inodese(:,ii,l),f,self%xnodes,self%order+1)
      DO j=1,3
        jac(j,:)=jac(j,:)+vtmp(j)*pt
      END DO
      DO j=1,6
        K(j,:) = K(j,:) + d2lag_3d(self%inodese(:,ii,l),dmap(1,j), &
          dmap(2,j),f,self%xnodes,self%order+1)*pt
      END DO
    end do
  END DO
  !---Get face points
  DO ii=1,self%ho_info%nfp
    do l=1,6
      m=ABS(self%lcf(l,cell))
      pt=self%ho_info%r(:,self%ho_info%lfp(ii,m))
      vtmp=dlag_3d(self%inodesf(:,ii,l),f,self%xnodes,self%order+1)
      DO m=1,3
        jac(m,:)=jac(m,:)+vtmp(m)*pt
      END DO
      DO j=1,6
        K(j,:) = K(j,:) + d2lag_3d(self%inodesf(:,ii,l),dmap(1,j), &
          dmap(2,j),f,self%xnodes,self%order+1)*pt
      END DO
    end do
  END DO
  !---Get cell points
  DO ii=1,self%ho_info%ncp
    pt=self%ho_info%r(:,self%ho_info%lcp(ii,cell))
    vtmp=dlag_3d(self%inodesc(:,ii),f,self%xnodes,self%order+1)
    DO j=1,3
      jac(j,:)=jac(j,:)+vtmp(j)*pt
    END DO
    DO j=1,6
      K(j,:) = K(j,:) + d2lag_3d(self%inodesc(:,ii),dmap(1,j), &
        dmap(2,j),f,self%xnodes,self%order+1)*pt
    END DO
  END DO
END IF
CALL hexmesh_g2inv(jac,g2op)
end subroutine hexmesh_hessian
!------------------------------------------------------------------------------
!> Compute the surface normal vector for a given face on a cell
!!
!! If face is not a global boundary face the function returns with `norm = 0`
!!
!! @note The logical position in the cell must be on the chosen face for this
!! subroutine, else an error will be thrown
!------------------------------------------------------------------------------
subroutine hexmesh_snormal(self,cell,ind,f,norm)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell
integer(i4), intent(in) :: ind !< Index of face within cell
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: norm(3) !< Unit vector normal to the face [3]
REAL(r8) :: goptmp(3,3),v
integer(i4), parameter :: gmap(6) = (/-3,-2,1,2,-1,3/)
DEBUG_STACK_PUSH
!---
if(cell<0.OR.cell>self%nc)call oft_abort('Invalid cell index.','hexmesh_snormal',__FILE__)
if(ind<0.OR.ind>6)call oft_abort('Invalid face index.','hexmesh_snormal',__FILE__)
! IF(f(ind)/=0._r8)call oft_abort('Invalid cell position.','hexmesh_snormal',__FILE__)
!---
! norm=0._r8
! IF(.NOT.self%global%gbc(i))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
! IF(.NOT.self%global%gbf(self%lcf(ind,i)))THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
!---
CALL hexmesh_jacobian(self,cell,f,goptmp,v)
norm = SIGN(1_i4,gmap(ind))*goptmp(:,ABS(gmap(ind)))
norm = norm/sqrt(sum(norm**2))
DEBUG_STACK_POP
end subroutine hexmesh_snormal
!------------------------------------------------------------------------------
!> Compute the curve tangent vector for a given edge on a cell
!!
!! If edge is not a global boundary edge the function returns with `tang = 0`
!!
!! @note The logical position in the cell must be on the chosen edge for this
!! subroutine to return a meaningful result
!------------------------------------------------------------------------------
subroutine hexmesh_ctang(self,cell,ind,f,tang)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: cell !< Index of cell
integer(i4), intent(in) :: ind !< Index of edge within cell
real(r8), intent(in) :: f(:) !< Logical coordinate in cell [4]
real(r8), intent(out) :: tang(3) !< Unit vector tangent to the edge [3]
INTEGER(i4) :: j,k
REAL(r8) :: goptmp(3,3),v,norm(3,2),e(3)
DEBUG_STACK_PUSH
!---
IF(cell<0.OR.cell>self%nc)CALL oft_abort('Invalid cell index.','hexmesh_ctang',__FILE__)
IF(ind<0.OR.ind>12)CALL oft_abort('Invalid edge index.','hexmesh_ctang',__FILE__)
!IF(f(ed1(ind))/=0._r8.AND.f(ed2(ind))/=0._r8)CALL oft_abort('Invalid cell position.','hexmesh_ctang',__FILE__)
!---
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
DO j=1,6
  IF(ANY(ABS(hex_fe(:,j))==ind))THEN
    CALL hexmesh_snormal(self, cell, j, f, norm(:,k))
    k=k+1
  END IF
END DO
e=self%r(:,self%lc(hex_ed(2,ind),cell))-self%r(:,self%lc(hex_ed(1,ind),cell))
tang=cross_product(norm(:,1),norm(:,2))
tang=tang/SQRT(SUM(tang**2))
tang=tang*SIGN(1_i4,self%lce(ind,cell))*SIGN(1._r8,DOT_PRODUCT(tang,e))
DEBUG_STACK_POP
end subroutine hexmesh_ctang
!------------------------------------------------------------------------------
!> Get mapping between boundary and volume logical coordinates
!------------------------------------------------------------------------------
subroutine hexmesh_get_surf_map(self,face,cell,lmap)
class(oft_hexmesh), intent(in) :: self !< Mesh object
integer(i4), intent(in) :: face !< Index of face on boundary mesh
integer(i4), intent(out) :: cell !< Cell containing face
integer(i4), intent(out) :: lmap(3) !< Coordinate mapping
integer(i4) :: i,j,k,iskip,find
real(r8) :: sgop(3,2),vgop(3,3),v,test(3),ftmp(3)
!---Find parent cell and face index in cell
find=ABS(self%bmesh%parent%lf(face))
cell=self%lfc(1,find)
DO j=1,6
  IF(ABS(self%lcf(j,cell))==find)EXIT
END DO
iskip = MAXLOC(ABS(inodes2f(:,j)-2), DIM=1)
!---Map surface logical coordinates to cell logical
call self%bmesh%jacobian(face,(/1.d0,1.d0/)/2.d0,sgop,v)
ftmp=1.d0/2.d0
IF(inodes2f(iskip,j)>2)THEN
  ftmp(iskip)=1.d0
  lmap(3)=-iskip
ELSE
  ftmp(iskip)=0.d0
  lmap(3)=iskip
END IF
call self%jacobian(cell,ftmp,vgop,v)
!---Find matching coordinates by gradient direction
DO k=1,2
  sgop(:,k)=sgop(:,k)/magnitude(sgop(:,k))
  test=0.d0
  DO i=1,3
    IF(k==1)vgop(:,i)=vgop(:,i)/magnitude(vgop(:,i))
    IF(i/=iskip)test(i)=DOT_PRODUCT(sgop(:,k),vgop(:,i))
  END DO
  i = MAXLOC(ABS(test), DIM=1)
  lmap(k)=INT(SIGN(REAL(i,8), test(i)),4)
END DO
end subroutine hexmesh_get_surf_map
!------------------------------------------------------------------------------
!> Map between surface and volume logical coordinates
!------------------------------------------------------------------------------
subroutine hexmesh_surf_to_vol(self,fsurf,lmap,fvol)
CLASS(oft_hexmesh), INTENT(in) :: self !< Mesh object
REAL(r8), INTENT(in) :: fsurf(:) !< Surface coordinates [2]
INTEGER(i4), INTENT(in) :: lmap(3) !< Coordinate mapping
REAL(r8), INTENT(out) :: fvol(:) !< Volume coordinates [3]
INTEGER(i4) :: k
DO k=1,2
  IF(lmap(k)<0)THEN
    fvol(ABS(lmap(k)))=1.d0-fsurf(k)
  ELSE IF(lmap(k)>0)THEN
    fvol(lmap(k))=fsurf(k)
  END IF
END DO
!---Static coordinate on face
IF(lmap(3)<0)THEN
  fvol(ABS(lmap(3)))=1.d0
ELSE IF(lmap(k)>0)THEN
  fvol(lmap(3))=0.d0
END IF
end subroutine hexmesh_surf_to_vol
!------------------------------------------------------------------------------
!> Invert a 3x3 matrix
!------------------------------------------------------------------------------
subroutine hexmesh_jacinv(A,C,j)
real(r8), intent(in) :: A(3,3) !< Matrix to invert
real(r8), intent(out) :: C(3,3) !< \f$ A^{-1} \f$
real(r8), intent(out) :: j !< |A|
real(r8) :: t1,t2,t3
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
end subroutine hexmesh_jacinv
!------------------------------------------------------------------------------
!> Expand and invert the matrix for the grid Hessian
!------------------------------------------------------------------------------
subroutine hexmesh_g2inv(jac,A)
real(r8), intent(in) :: jac(3,3)
real(r8), intent(out) :: A(6,6)
REAL(r8) :: work(12)
INTEGER(i4) :: info,ipiv(6)
DEBUG_STACK_PUSH
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
DEBUG_STACK_POP
end subroutine hexmesh_g2inv
END MODULE oft_hexmesh_type
