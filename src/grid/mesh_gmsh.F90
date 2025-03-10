!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_gmsh.F90
!
!> Mesh handling for GMSH meshes
!!
!! Functions to read in, intialize and refine GMSH meshes
!! - Mesh file read-in
!! - Boundary point location and refinement
!!
!! The CAD interface to [GMSH](https://geuz.org/gmsh/) uses a high-order, currently
!! quadratic, triangular tessellation of the boundary. This representation is used
!! to provide refinement and high order node placement. Note that this is not
!! equivalent to the CAD interface for @ref oft_mesh_cubit "CUBIT" or @ref oft_mesh_t3d
!! "T3D" as the form cannot analytically represent most model surfaces.
!!
!! @authors Chris Hansen
!! @date October 2013
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_gmsh
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_mesh_local, ONLY: bmesh_local_init
USE oft_mesh_local_util, ONLY: mesh_local_findedge, mesh_local_findface
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE multigrid, ONLY: multigrid_mesh, multigrid_level
!---End include modules
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> GMSH CAD linkage structure
!! - Linkage of mesh entities to CAD model
!------------------------------------------------------------------------------
type :: GMSH_cadlink
  integer(i4), pointer, dimension(:) :: lbfg => NULL() !< Linkage of mesh faces to CAD entities
end type GMSH_cadlink
!---Global variables
character(LEN=OFT_PATH_SLEN) :: filename = 'none' !< Name of GMSH mesh file
integer(i4) :: order = 1 !< Order of base mesh
type(oft_trimesh) :: cad_mesh !< Surface mesh representation of CAD geometry
type(GMSH_cadlink), pointer :: cad_link => NULL() !< Linkage of mesh to CAD geometry
type(GMSH_cadlink), pointer, dimension(:) :: ML_cad_link => NULL() !< ML CAD linkage
integer(i4) :: active_face = 0 !< Active face for MINPACK fitting
real(r8) :: active_endpts(3,3) = 0. !< Active constraint points for MINPACK fitting
real(r8) :: active_wts(3) = 0. !< Active constraint weights for MINPACK fitting
!$omp threadprivate(active_endpts,active_wts,active_face)
INTEGER(i4), PARAMETER, PUBLIC :: mesh_gmsh_id = 3
public mesh_gmsh_load, mesh_gmsh_cadlink, mesh_gmsh_reffix
public mesh_gmsh_add_quad, gmsh_finalize_setup
contains
!------------------------------------------------------------------------------
!> Finalize setup/load-in of GMSH mesh and destroy temporaries created
!! for grid construction (eg. high-order input nodes, in-memory data)
!------------------------------------------------------------------------------
subroutine gmsh_finalize_setup
integer(i4) :: i,n
CALL cad_mesh%delete()
!
IF(ASSOCIATED(ML_cad_link))THEN
  n=SIZE(ML_cad_link)
  DO i=1,n
    IF(ASSOCIATED(ML_cad_link(i)%lbfg))DEALLOCATE(ML_cad_link(i)%lbfg)
  END DO
  DEALLOCATE(ML_cad_link)
END IF
NULLIFY(cad_link)
end subroutine gmsh_finalize_setup
!------------------------------------------------------------------------------
!> Read in GMSH mesh file from file "filename"
!! - Read in GMSH options from input file
!! - Read in mesh points and cells
!------------------------------------------------------------------------------
subroutine mesh_gmsh_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j,ed,nptmp,id,ierr,io_unit
integer(i4) :: np_cad,ne_cad,nf_cad,nep,nfp,ncp
integer(i4), allocatable, dimension(:) :: cad_ptmp,ptmp,ctmp
integer(i4), allocatable, dimension(:,:) :: lftmp
real(r8), allocatable, dimension(:,:) :: rtmp
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
!---Read in mesh options
namelist/gmsh_options/filename,order
DEBUG_STACK_PUSH
filename = 'none'
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,gmsh_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "gmsh_options" found in input file.','mesh_gmsh_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "gmsh_options" in input file.','mesh_gmsh_load',__FILE__)
  IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','mesh_gmsh_load',__FILE__)
  WRITE(*,*)
  WRITE(*,'(A)')'**** Loading GMSH mesh'
  WRITE(*,'(2X,2A)')  'Mesh File = ',TRIM(filename)
  WRITE(*,'(2X,A,I4)')'Order     = ',order
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(filename,20,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_gmsh_load',__FILE__)
CALL MPI_Bcast(order,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_gmsh_load',__FILE__)
#endif
!---Allocate mesh
allocate(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%meshes(i)%setup(mesh_gmsh_id)
  CALL mg_mesh%smeshes(i)%setup(mesh_gmsh_id,.TRUE.)
  mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
END DO
CALL multigrid_level(mg_mesh,1)
mesh=>mg_mesh%meshes(1)
smesh=>mg_mesh%smeshes(1)
!---Setup geometry information
allocate(ML_cad_link(mg_mesh%mgdim))
cad_link=>ML_cad_link(mg_mesh%level)
nep=0
nfp=0
ncp=0
IF(order==1)THEN
  nep=2
  nfp=3
  ncp=4
ELSE IF(order==2)THEN
  nep=3
  nfp=6
  ncp=10
ELSE
  CALL oft_abort('Invalid GMSH order','mesh_gmsh_load',__FILE__)
END IF
!---Read in T3D mesh
open(NEWUNIT=io_unit,FILE=trim(filename))
DO i=1,4 ! Skip header lines
  read(io_unit,*)
END DO
read(io_unit,*)nptmp
!---Read in points
if(nptmp==0)call oft_abort('No points in GMSH file','mesh_gmsh_load',__FILE__)
allocate(rtmp(3,nptmp))
do i=1,nptmp
  read(io_unit,*)rtmp(:,i),id
end do
!---Skip boundary edges
read(io_unit,*)
read(io_unit,*)ne_cad
DO i=1,ne_cad
  read(io_unit,*)
END DO
!---Read-in boundary faces
read(io_unit,*)
read(io_unit,*)nf_cad
ALLOCATE(lftmp(nfp,nf_cad))
DO i=1,nf_cad
  read(io_unit,*)lftmp(:,i),id
END DO
!---Read in cells
read(io_unit,*)
read(io_unit,*)mesh%nc
!---Read in cells
if(mesh%nc==0)call oft_abort('No cells in GMSH file','mesh_gmsh_load',__FILE__)
allocate(mesh%lc(4,mesh%nc),ctmp(ncp),mesh%reg(mesh%nc))
mesh%reg=1
ALLOCATE(ptmp(nptmp))
ptmp=0; mesh%np=0
do i=1,mesh%nc
  read(io_unit,*)ctmp,id
  DO j=1,4
    IF(ptmp(ctmp(j))==0)THEN
      mesh%np=mesh%np+1
      ptmp(ctmp(j))=mesh%np
    END IF
    mesh%lc(j,i)=ptmp(ctmp(j))
  END DO
end do
deallocate(ctmp)
ALLOCATE(mesh%r(3,mesh%np))
DO i=1,nptmp
  IF(ptmp(i)/=0)mesh%r(:,ptmp(i))=rtmp(:,i)
END DO
close(io_unit)
!---Setup surface mesh
IF(oft_debug_print(2))WRITE(*,*)'  Creating CAD surface mesh'
ALLOCATE(cad_ptmp(nptmp))
cad_ptmp=0
np_cad=0
DO i=1,nf_cad
  DO j=1,3
    IF(cad_ptmp(lftmp(j,i))==0)THEN
      np_cad=np_cad+1
      cad_ptmp(lftmp(j,i))=np_cad
    END IF
  END DO
END DO
!---
CALL cad_mesh%setup(-1,.TRUE.)
cad_mesh%np=np_cad
ALLOCATE(cad_mesh%r(3,cad_mesh%np),cad_mesh%parent%lp(cad_mesh%np))
DO i=1,nptmp
  IF(cad_ptmp(i)/=0)THEN
    cad_mesh%r(:,cad_ptmp(i))=rtmp(:,i)
    cad_mesh%parent%lp(cad_ptmp(i))=ptmp(i)
  END IF
END DO
cad_mesh%nc=nf_cad
ALLOCATE(cad_mesh%lc(3,cad_mesh%nc))
DO i=1,cad_mesh%nc
  cad_mesh%lc(:,i)=cad_ptmp(lftmp(1:3,i))
END DO
CALL bmesh_local_init(cad_mesh)
!---Set quadratic nodes
IF(order==2)THEN
  CALL cad_mesh%set_order(2)
  DO i=1,cad_mesh%nc
    DO j=1,3
      ed=ABS(cad_mesh%lce(j,i))
      IF(cad_mesh%ho_info%lep(1,ed)<=0)THEN
        cad_mesh%ho_info%lep(1,ed)=ed
        cad_mesh%ho_info%r(:,ed)=rtmp(:,lftmp(4+MOD(j,3),i))
      END IF
    END DO
  END DO
END IF
IF(oft_debug_print(2))WRITE(*,*)'  Complete'
!---
call mesh_global_resolution(mesh)
DEBUG_STACK_POP
end subroutine mesh_gmsh_load
!------------------------------------------------------------------------------
!> Link GMSH CAD objects to mesh entities for use in refinement.
!------------------------------------------------------------------------------
subroutine mesh_gmsh_cadlink(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,ind,fp(3)
integer(i4), allocatable :: fmap(:)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,*)'Linking GMSH geometry to mesh'
!---Get face boundary mapping
allocate(fmap(mesh%nf))
allocate(cad_link%lbfg(mesh%nbf)) ! Create face linkage
cad_link%lbfg=0
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf) ! Get face map
!---Search CAD model faces for local faces
do i=1,cad_mesh%nc
  do j=1,3 ! Get corner points for CAD face
    fp(j)=INT(cad_mesh%parent%lp(cad_mesh%lc(j,i)),4)
  enddo
  ind=abs(mesh_local_findface(mesh,fp)) ! Find face on mesh
  if(ind==0)cycle ! Not on my domain
  if(fmap(ind)==0)cycle
  cad_link%lbfg(fmap(ind))=i
end do
!---Destroy mapping arrays
deallocate(fmap)
CALL mesh_gmsh_hobase(mesh)
if(oft_debug_print(2))write(*,*)'  Complete'
DEBUG_STACK_POP
end subroutine mesh_gmsh_cadlink
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from high order import
!------------------------------------------------------------------------------
subroutine mesh_gmsh_hobase(mesh)
class(oft_mesh), intent(inout) :: mesh
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2),ind
DEBUG_STACK_PUSH
IF(order==1)RETURN
if(oft_debug_print(1))write(*,*)'Quadratic mesh nodes imported'
!---Setup quadratic mesh
mesh%order=1
mesh%ho_info%nep=1
ALLOCATE(mesh%ho_info%r(3,mesh%ne),mesh%ho_info%lep(mesh%ho_info%nep,mesh%ne))
!---Initialize high order points with straight edges
DO i=1,mesh%ne
  mesh%ho_info%lep(1,i)=i
  mesh%ho_info%r(:,i)=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
END DO
!---Search CAD model faces for local faces
do i=1,cad_mesh%ne
  etmp(1)=cad_mesh%parent%lp(cad_mesh%le(1,i))
  etmp(2)=cad_mesh%parent%lp(cad_mesh%le(2,i))
  !
  ind=abs(mesh_local_findedge(mesh,etmp)) ! Find face on mesh
  if(ind==0)CALL oft_abort('Unlinked mesh edge','mesh_gmsh_hobase',__FILE__)
  mesh%ho_info%r(:,ind) = cad_mesh%ho_info%r(:,i)
end do
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine mesh_gmsh_hobase
!------------------------------------------------------------------------------
!> Adjust boundary points to CAD boundary.
!------------------------------------------------------------------------------
subroutine mesh_gmsh_reffix(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
real(r8) :: pt(3)
integer(i4) :: i,ierr,j,k,ind,ed,fp(3),npcors,edge,face
integer(i4), pointer :: tmp(:)
integer(i4), allocatable :: fmap(:)
class(oft_mesh), pointer :: pmesh,mesh
type(GMSH_cadlink), pointer :: pmesh_cad_link
DEBUG_STACK_PUSH
!---Get parent mesh
mesh=>mg_mesh%mesh
pmesh=>mg_mesh%meshes(mg_mesh%level-1)
!---If only one level do nothing
if(mg_mesh%level==1)THEN
  DEBUG_STACK_POP
  return
END IF
if(pmesh%fullmesh.AND.(.NOT.mesh%fullmesh))then ! Current level is transfer level
  !---Synchronize T3D linkage to distributed mesh
  if(oft_debug_print(1))write(*,*)'Copying geometry linkage to distributed mesh'
  !---Get global point mapping
  allocate(tmp(mesh%global%np))
  tmp=0
  !$omp parallel do
  do i=1,mesh%np
    tmp(mesh%base%lp(i))=i
  end do
  !---Get CAD representation aliases
  pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
  cad_link=>ML_cad_link(mg_mesh%level)
  !---Get face boundary mapping
  allocate(fmap(mesh%nf))
  CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf)
  !---Copy parent face linkage to new mesh
  allocate(cad_link%lbfg(mesh%nbf))
  cad_link%lbfg=0
  do i=1,pmesh%nbf
    fp=tmp(pmesh%lf(:,pmesh%lbf(i)))
    if(any(fp==0).OR.any(fp>mesh%np))cycle
    ind=abs(mesh_local_findface(mesh,fp))
    if(ind==0)cycle
    if(fmap(ind)==0)call oft_abort('Bad face link','mesh_gmsh_reffix',__FILE__)
    if(pmesh_cad_link%lbfg(i)/=0)cad_link%lbfg(fmap(ind))=pmesh_cad_link%lbfg(i)
  enddo
  !---Destroy mapping arrays
  deallocate(tmp,fmap)
  if(oft_debug_print(1))write(*,*)'Complete'
  DEBUG_STACK_POP
  return
endif
!---Refine new boundary points using CAD model
if(oft_debug_print(1))write(*,*)'Adjusting points to GMSH boundary'
!---Get CAD representation aliases
pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
cad_link=>ML_cad_link(mg_mesh%level)
!---Locate edge end points and place daughter point
!!$omp parallel do private(j,pt,ierr,obj_curve,obj_surf)
do i=1,pmesh%nbf
  face=pmesh%lbf(i)
  IF(pmesh_cad_link%lbfg(i)==0)cycle
  do j=1,3
    edge=ABS(pmesh%lfe(j,face))
    call gmsh_surf_midpoint(pmesh_cad_link%lbfg(i),pt,pmesh%r(:,pmesh%le(1,edge)), &
      pmesh%r(:,pmesh%le(2,edge)),.5d0,.5d0,ierr)
    mesh%r(:,pmesh%np+edge)=pt
  end do
enddo
!---Get edge and face boundary mapping
allocate(fmap(mesh%nf))
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf) ! Get face map
allocate(cad_link%lbfg(mesh%nbf))
cad_link%lbfg=0
npcors=pmesh%np
!---Copy face linkage to daughter edges and faces
!$omp parallel do private(i,k,ed)
do j=1,pmesh%nbf
  i=pmesh%lbf(j)
  do k=1,4
    ed=fmap(mg_mesh%inter(mg_mesh%level-1)%lfdf(k,i))
    if(pmesh_cad_link%lbfg(j)/=0)cad_link%lbfg(ed)=pmesh_cad_link%lbfg(j)
  end do
enddo
!---Destroy mapping arrays
deallocate(fmap)
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine mesh_gmsh_reffix
!------------------------------------------------------------------------------
!> Add quadratic mesh node points using CAD model
!------------------------------------------------------------------------------
subroutine mesh_gmsh_add_quad(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
class(oft_mesh), pointer :: mesh
real(r8) :: pt(3)
integer(i4) :: i,ierr,j,k,edge,face
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,*)'Setting GMSH quadratic nodes'
!---Get CAD representation alias
mesh=>mg_mesh%mesh
cad_link=>ML_cad_link(mg_mesh%level)
!---Locate edge end points and place daughter node
do i=1,mesh%ne
  mesh%ho_info%lep(1,i)=i
  mesh%ho_info%r(:,i)=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
end do
do i=1,mesh%nbf
  face=mesh%lbf(i)
  IF(cad_link%lbfg(i)==0)cycle
  do j=1,3
    edge=ABS(mesh%lfe(j,face))
    call gmsh_surf_midpoint(cad_link%lbfg(i),pt,mesh%r(:,mesh%le(1,edge)), &
    mesh%r(:,mesh%le(2,edge)),.5d0,.5d0,ierr)
    mesh%ho_info%r(:,edge)=pt
    ! !
    ! DO k=mesh%kec(edge),mesh%kec(edge+1)-1
    !   mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
    ! END DO
  end do
enddo
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine mesh_gmsh_add_quad
!---------------------------------------------------------------------------
!> Compute the weighted midpoint of a surface edge
!!
!! Locates the point on a given face of the imported GMSH boundary mesh by
!! minimizing the weighted sum of distances to 2 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!---------------------------------------------------------------------------
subroutine gmsh_surf_midpoint(face,pt,pt1,pt2,wt1,wt2,ierr)
integer(i4), intent(in) :: face !< Face index
real(r8), intent(inout) :: pt(3) !< Solution point
real(r8), intent(in) :: pt1(3) !< Constraint point 1
real(r8), intent(in) :: pt2(3) !< Constraint point 2
real(r8), intent(in) :: wt1 !< Constraint weight 1
real(r8), intent(in) :: wt2 !< Constraint weight 2
integer(i4), intent(out) :: ierr !< Error flag
!---
integer(i4), parameter :: nerr=3
integer(i4), parameter :: neq=2
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor,f(3)
real(r8) :: rt(3),uv1(2),uv2(2),u,v
integer :: i,info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---Set initial guess
u=.3d0
v=.7d0
ierr=0
!---Setup MINPACK parameters
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-6
xtol = 1.d-6
gtol = 1.d-3
epsfcn = 1.d-6
nprint = 0
ldfjac = nerr
!---Setup active points
active_face=face
!---Find first point
uv=(/u,v/)
active_endpts(1,:)=pt1; active_wts(1)=wt1
call lmdif(gmsh_spt_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
uv1=uv
IF(info<0)WRITE(*,*)'Surface midpoint failed',face
!---Find second point
uv=(/u,v/)
active_endpts(1,:)=pt2; active_wts(1)=wt1
call lmdif(gmsh_spt_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
uv2=uv
IF(info<0)WRITE(*,*)'Surface midpoint failed',face
!---Evaluate midpoint
uv=(uv1+uv2)/2.d0
f(1)=uv(1)
f(2)=uv(2)
f(3)=1.d0-f(1)-f(2)
pt=cad_mesh%log2phys(active_face,f)
DEBUG_STACK_POP
end subroutine gmsh_surf_midpoint
!---------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active point
!! used in a 1 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref mesh_gmsh::gmsh_surf_midpoint "gmsh_surf_midpoint"
!---------------------------------------------------------------------------
subroutine gmsh_spt_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions [3]
integer(i4), intent(in) :: n !< Number of parametric dimensions [2]
real(r8), intent(in) :: uv(n) !< Parametric possition [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3),f(3)
DEBUG_STACK_PUSH
f(1)=uv(1)
f(2)=uv(2)
f(3)=1.d0-f(1)-f(2)
pt=cad_mesh%log2phys(active_face,f)
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
DEBUG_STACK_POP
end subroutine gmsh_spt_error
end module oft_mesh_gmsh
