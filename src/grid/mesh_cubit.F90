!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_cubit.F90
!
!> Mesh handling for CUBIT meshes
!!
!! @author Chris Hansen
!! @date July 2012
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
module oft_mesh_cubit
#ifdef HAVE_NCDF
USE, INTRINSIC :: iso_c_binding, ONLY: C_LOC, C_CHAR, C_NULL_CHAR
USE netcdf
USE oft_base
#ifdef HAVE_ONURBS
USE nurbs_cad, ONLY: nurbs_curve, nurbs_surf, nurbs_entity_ptr, nurbs_init, &
  nurbs_read_in, nurbs_get_count, nurbs_curve_domain, nurbs_curve_periodic, &
  nurbs_curve_name, nurbs_surf_domain, nurbs_surf_periodic, nurbs_surf_singular, &
  nurbs_surf_name, nurbs_surf_avg, nurbs_surf_midpoint, nurbs_curve_midpoint, &
  nurbs_surf_center, nurbs_curve_linear, nurbs_surf_planar
#endif
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh, mesh, oft_bmesh, smesh
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_hexmesh_type, ONLY: oft_hexmesh
USE oft_quadmesh_type, ONLY: oft_quadmesh
USE oft_mesh_local_util, ONLY: mesh_local_findedge, mesh_local_findface
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE multigrid, ONLY: mg_mesh
!---End include modules
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
! TYPE Exodus_curve
!------------------------------------------------------------------------------
!> T3D CAD boundary structure
!! - CAD entity counts
!! - CAD wireframe entities
!------------------------------------------------------------------------------
TYPE :: Exodus_curve
  INTEGER(i4) :: cid = -1 !< Curve ID in OFT ordering
  INTEGER(i4) :: top_cid = -1 !< Curve ID from geometry file
  INTEGER(i4) :: ne = 0 !< Number of geometry model edges
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lgme => NULL() !< List of geometry model edges
  CHARACTER(LEN=OFT_PATH_SLEN) :: name = ''
END TYPE Exodus_curve
!------------------------------------------------------------------------------
! TYPE Exodus_surf
!------------------------------------------------------------------------------
!> T3D CAD boundary structure
!! - CAD entity counts
!! - CAD wireframe entities
!------------------------------------------------------------------------------
TYPE :: Exodus_surf
  INTEGER(i4) :: sid = -1 !< Surface ID in OFT ordering
  INTEGER(i4) :: top_sid = -1 !< Surface ID in OFT ordering
  INTEGER(i4) :: nf = 0 !< Number of geometry model faces
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lgmf => NULL() !< List of geometry model faces
  CHARACTER(LEN=OFT_PATH_SLEN) :: name = ''
END TYPE Exodus_surf
!------------------------------------------------------------------------------
! TYPE Exodus_cadlink
!------------------------------------------------------------------------------
!> T3D CAD linkage structure
!! - Linkage of mesh entities to CAD model
!------------------------------------------------------------------------------
#ifdef HAVE_ONURBS
TYPE :: Exodus_cadlink
  INTEGER(i4) :: nce = 0 !< Number of CAD linked edges
  INTEGER(i4) :: ncf = 0 !< Number of CAD linked faces
  INTEGER(i4), POINTER, DIMENSION(:) :: lce => NULL() !< List of CAD linked edges
  INTEGER(i4), POINTER, DIMENSION(:) :: lcf => NULL() !< List of CAD linked faces
  TYPE(nurbs_entity_ptr), POINTER, DIMENSION(:) :: lbeg => NULL() !< Linkage of mesh edges to CAD entities
  TYPE(nurbs_entity_ptr), POINTER, DIMENSION(:) :: lbfg => NULL() !< Linkage of mesh faces to CAD entities
END TYPE Exodus_cadlink
#endif
PRIVATE
INTEGER(i4), PARAMETER, PUBLIC :: mesh_cubit_id = 2
PUBLIC mesh_cubit_load, mesh_cubit_reffix, mesh_cubit_add_quad
PUBLIC mesh_cubit_hobase, mesh_cubit_cadlink
PUBLIC mesh_cubit_test_edge, mesh_cubit_set_periodic
PUBLIC smesh_cubit_load, cubit_read_nodesets, cubit_read_sidesets
!---Cubit "exodus" indexing maps
integer(i4), parameter :: exodus_tri_emap(2,3) = RESHAPE([1,2, 2,3, 3,1], [2,3])
integer(i4), parameter :: exodus_tet_emap(2,6) = RESHAPE((/1,2, 2,3, 3,1, 1,4, 2,4, 3,4/), (/2,6/))
integer(i4), parameter :: exodus_quad_emap(2,4) = RESHAPE([1,2, 2,3, 3,4, 4,1], [2,4])
integer(i4), parameter :: exodus_hex_emap(2,12) = RESHAPE((/1,2, 2,3, 3,4, 1,4, 1,5, 2,6, &
  3,7, 4,8, 5,6, 6,7, 7,8, 5,8/),(/2,12/))
integer(i4), parameter :: exodus_hex_fmap(4,6) = RESHAPE((/1,2,3,4, 5,6,7,8, 1,4,8,5, &
  2,3,7,6, 1,2,6,5, 3,4,8,7/),(/4,6/))
!---Global variables
INTEGER(i4), PARAMETER :: cubit_soffset=1E4
CHARACTER(LEN=OFT_PATH_SLEN), PUBLIC :: inpname = 'none' !< Name of Cubit input file for geometry (Used to retrieve CAD objects)
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'none' !< Name of Cubit input file for mesh
LOGICAL :: lf_file = .TRUE. !< Large format file flag
LOGICAL :: tor_mesh = .FALSE. !< Curve grid to toroidal shaping
LOGICAL :: reflect = .FALSE. !< Logical flag for mesh reflection (z-direction)
INTEGER(i4) :: per_ns = -1 !< Integer index of periodic nodeset
REAL(r8) :: zstretch = 1.d0 !< Scale for z-coordinates (useful for cylindrical pinch studies)
REAL(r8) :: tor_rmin = 0.d0
!---
#ifdef HAVE_ONURBS
INTEGER(i4) :: ngmc = 0 !< Number of geometry model curves
INTEGER(i4) :: ngms = 0 !< Number of geometry model surfaces
INTEGER(i4) :: ngwc = 0 !< Number of geometry wireframe curves
INTEGER(i4) :: ngws = 0 !< Number of geometry wireframe surfaces
TYPE(Exodus_curve), POINTER, DIMENSION(:) :: model_curves => NULL() !< List of model curves
TYPE(Exodus_surf), POINTER, DIMENSION(:) :: model_surfaces => NULL() !< List of model surfaces
TYPE(nurbs_curve), POINTER, DIMENSION(:) :: wf_curves => NULL() !< List of CAD wireframe curves
TYPE(nurbs_surf), POINTER, DIMENSION(:) :: wf_surfs => NULL() !< List of CAD wireframe surfaces
!---
TYPE(Exodus_cadlink), POINTER :: cad_link => NULL() !< Linkage of mesh to CAD geometry
TYPE(Exodus_cadlink), POINTER, DIMENSION(:) :: ML_cad_link => NULL() !< ML CAD linkage
!---
REAL(r4) :: cubit_version = 0.
INTEGER(i4), PARAMETER :: ex_topc_len=5
INTEGER(i4), PARAMETER :: ex_tops_len=5
INTEGER(i4), PARAMETER :: exodus_string_len=33
#endif
!---
INTEGER(i4) :: ncid = 0
INTEGER(i4) :: nblks = 0
INTEGER(i4) :: nregions = 0
INTEGER(i4) :: np_per = 0
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: per_nodes
!---
LOGICAL :: have_ho = .FALSE.
INTEGER(i4) :: np_ho = 0
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: lc_ho
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: r_ho
!
! TYPE(oft_1d_int), allocatable, target :: cubit_nsets(:)
! TYPE(oft_1d_int), allocatable, target :: cubit_ssets(:)
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_load
!------------------------------------------------------------------------------
!> Read in Exodus mesh and geometry information
!! - Read in Cubit options from input file
!! - Read in mesh points and cells
!! - Read in surface mesh
!! - Load and initialize OpenNURBS CAD representation
!------------------------------------------------------------------------------
subroutine mesh_cubit_load
real(r8), allocatable :: rtmp(:,:)
integer(i4) :: blkID,nodesID,att_len
integer(i4) :: i,j,id,it,lenreflag,ierr,io_unit
integer(i4), allocatable :: lptmp(:)
logical :: ltrans(3)
character(40) :: eltype
!---Read in mesh options
namelist/cubit_options/filename,inpname,lf_file,tor_mesh, &
  reflect,per_ns,zstretch,tor_rmin
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,cubit_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "cubit_options" found in input file.','mesh_cubit_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "cubit_options" in input file.','mesh_cubit_load',__FILE__)
  IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','mesh_cubit_load',__FILE__)
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'**** Loading Cubit mesh'
  CALL oft_increase_indent
  WRITE(*,'(3A)')oft_indent,'Mesh File = ',TRIM(filename)
  WRITE(*,'(3A)')oft_indent,'Geom File = ',TRIM(inpname)
ELSE
  CALL oft_increase_indent
END IF
!---Broadcast input information
#ifdef HAVE_MPI
ltrans=(/lf_file,tor_mesh,reflect/)
CALL MPI_Bcast(filename,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cubit_load',__FILE__)
CALL MPI_Bcast(inpname,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cubit_load',__FILE__)
CALL MPI_Bcast(ltrans,3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cubit_load',__FILE__)
CALL MPI_Bcast(per_ns,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cubit_load',__FILE__)
CALL MPI_Bcast(zstretch,1,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_cubit_load',__FILE__)
lf_file=ltrans(1)
tor_mesh=ltrans(2)
reflect=ltrans(3)
#endif
!---
call mesh_cubit_error(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
!---Check mesh type
call mesh_cubit_error(NF90_inq_varid(ncid,"connect1",blkID))
call mesh_cubit_error(NF90_INQUIRE_ATTRIBUTE(ncid,blkID,"elem_type",len=att_len))
call mesh_cubit_error(NF90_GET_ATT(ncid,blkID,"elem_type",eltype))
!---Allocate mesh
IF(eltype(1:3)=="TET")THEN
  ALLOCATE(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
  ALLOCATE(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
ELSE IF(eltype(1:3)=="HEX")THEN
  ALLOCATE(oft_hexmesh::mg_mesh%meshes(mg_mesh%mgdim))
  ALLOCATE(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
END IF
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%meshes(i)%setup(mesh_cubit_id)
  CALL mg_mesh%smeshes(i)%setup(mesh_cubit_id,.TRUE.)
  mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
END DO
mesh=>mg_mesh%meshes(1)
smesh=>mg_mesh%smeshes(1)
!---
call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_blk",blkID))
call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_nodes",nodesID))
call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,blkID,len = nblks))
call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nodesID,len = mesh%np))
!---
call mesh_cubit_read_regions
#ifdef HAVE_ONURBS
!---Setup geometry information
allocate(ML_cad_link(mg_mesh%mgdim))
cad_link=>ML_cad_link(mg_mesh%level)
!---
IF(TRIM(inpname)/='none')THEN
  call mesh_cubit_read_surface
  call mesh_cubit_geom
ELSE
  IF(oft_debug_print(2))WRITE(*,'(A)')oft_indent,'No Cubit geometry file specified'
END IF
#endif
!---Load mesh vertices
IF(lf_file)THEN
  allocate(mesh%r(3,mesh%np),rtmp(mesh%np,1))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coordx",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp(:,1)))
  mesh%r(1,:)=rtmp(:,1)
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coordy",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp(:,1)))
  mesh%r(2,:)=rtmp(:,1)
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coordz",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp(:,1)))
  mesh%r(3,:)=rtmp(:,1)
  deallocate(rtmp)
ELSE
  allocate(mesh%r(3,mesh%np),rtmp(mesh%np,3))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coord",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp))
  mesh%r(1,:)=rtmp(:,1)
  mesh%r(2,:)=rtmp(:,2)
  mesh%r(3,:)=rtmp(:,3)
  deallocate(rtmp)
END IF
call mesh_cubit_error(NF90_CLOSE(ncid))
!---Remove high order points from vertex list
IF(have_ho)THEN
  ALLOCATE(r_ho(3,mesh%np))
  r_ho=mesh%r
  np_ho=mesh%np
  !---Index points used in cell vertices
  ALLOCATE(lptmp(mesh%np))
  lptmp=0
  mesh%np=0
  DO i=1,mesh%nc
    DO j=1,mesh%cell_np
      IF(lptmp(mesh%lc(j,i))==0)THEN
        mesh%np = mesh%np + 1
        lptmp(mesh%lc(j,i)) = mesh%np
      END IF
    END DO
  END DO
  !---Reset point list
  DEALLOCATE(mesh%r)
  ALLOCATE(mesh%r(3,mesh%np))
  DO i=1,np_ho
    IF(lptmp(i)/=0)mesh%r(:,lptmp(i)) = r_ho(:,i)
  END DO
  !---Reindex cells
  DO i=1,mesh%nc
    DO j=1,mesh%cell_np
      mesh%lc(j,i) = lptmp(mesh%lc(j,i))
      lc_ho(j,i) = lptmp(lc_ho(j,i))
    END DO
  END DO
#ifdef HAVE_ONURBS
  !---Reindex CAD edges
  DO i=1,ngmc
    DO j=1,model_curves(i)%ne
      model_curves(i)%lgme(:,j) = lptmp(model_curves(i)%lgme(:,j))
    END DO
  END DO
  !---Reindex CAD faces
  DO i=1,ngms
    DO j=1,model_surfaces(i)%nf
      model_surfaces(i)%lgmf(:,j) = lptmp(model_surfaces(i)%lgmf(:,j))
    END DO
  END DO
#endif
  DEALLOCATE(lptmp)
  r_ho(3,:)=r_ho(3,:)*zstretch
END IF
mesh%r(3,:)=mesh%r(3,:)*zstretch
IF(reflect)THEN
  call mesh_global_resolution(mesh)
  call mesh_cubit_reflect(.1d0*mesh%hmin)
END IF
IF(oft_env%rank/=0)DEALLOCATE(mesh%r,mesh%lc,mesh%reg)
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,I8)')oft_indent,'# of Regions      = ',nblks
#ifdef HAVE_ONURBS
  WRITE(*,'(2A,2I8)')oft_indent,'# of CAD curves   = ',ngwc,ngmc
  WRITE(*,'(2A,2I8)')oft_indent,'# of CAD surfaces = ',ngws,ngms
#else
  WRITE(*,'(2A)')oft_indent,'Skipping CAD information: Not compiled with OpenNURBS'
#endif
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
contains
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_read_regions
!------------------------------------------------------------------------------
!> Read in Exodus volume mesh
!! - Merge all blocks
!------------------------------------------------------------------------------
subroutine mesh_cubit_read_regions
integer(i4) :: blkID,elemID,lcID,att_len
integer(i4) :: i,j,nc
integer(i4), allocatable :: rmark(:),lctmp(:,:)
character(3) :: blknm
character(40) :: eltype
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Loading Cubit regions'
CALL oft_increase_indent
!---
allocate(rmark(nblks))
rmark=0
nregions=0
mesh%nc=0
!---
do i=1,nblks
  !---
  if(i<10)then
    write(blknm,'(I1)')i
  else if(i>9.AND.i<100)then
    write(blknm,'(I2)')i
  else if(i>99.AND.i<1000)then
    write(blknm,'(I3)')i
  else if(i>=1000)then
    CALL oft_abort('Invalid Exodus variable index (connectXXX)','mesh_cubit_read_regions',__FILE__)
  end if
  IF(oft_debug_print(2))WRITE(*,'(3A)')oft_indent,'Loading cell count for region ',TRIM(blknm)
  !---
  call mesh_cubit_error(NF90_inq_varid(ncid,"connect"//TRIM(blknm),blkID))
  call mesh_cubit_error(NF90_INQUIRE_ATTRIBUTE(ncid,blkID,"elem_type",len=att_len))
  !---
  call mesh_cubit_error(NF90_GET_ATT(ncid,blkID,"elem_type",eltype))
  !---
  if(eltype(1:6)=="TETRA4")then
    rmark(i)=1
  else if(eltype(1:7)=="TETRA10")then
    rmark(i)=2
    have_ho=.TRUE.
  else if(eltype(1:4)=="HEX8")then
    rmark(i)=3
  else if(eltype(1:5)=="HEX27")then
    rmark(i)=4
    have_ho=.TRUE.
  ELSE
    !---Check for legacy definition
    IF(eltype(1:5)=="TETRA")THEN
      rmark(i)=1
    ELSE
      CYCLE
    END IF
  end if
  nregions=nregions+1
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_in_blk"//TRIM(blknm),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=nc))
  mesh%nc=mesh%nc+nc
end do
IF(mesh%nc==0)CALL oft_abort('No volume regions found!','mesh_cubit_read_regions',__FILE__)
IF(ANY(rmark==1).AND.ANY(rmark==2))CALL oft_abort('Two conflicting element types found', &
  'mesh_cubit_read_regions',__FILE__)
IF(have_ho)THEN
  IF(MAXVAL(rmark)==2)THEN
    ALLOCATE(lc_ho(10,mesh%nc))
  ELSE IF(MAXVAL(rmark)==4)THEN
    ALLOCATE(lc_ho(27,mesh%nc))
  END IF
  lc_ho=0
END IF
!---
allocate(mesh%lc(mesh%cell_np,mesh%nc),mesh%reg(mesh%nc))
mesh%nc=0
!---
do i=1,nblks
  if(rmark(i)==0)cycle
  !---
  if(i<10)then
    write(blknm,'(I1)')i
  else if(i>9.AND.i<100)then
    write(blknm,'(I2)')i
  else if(i>99.AND.i<1000)then
    write(blknm,'(I3)')i
  end if
  IF(oft_debug_print(2))WRITE(*,'(3A)')oft_indent,'Loading cell lists for region ',TRIM(blknm)
  !---
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_in_blk"//TRIM(blknm),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=nc))
  !---
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"connect"//TRIM(blknm),lcID))
  IF(rmark(i)==1)THEN
    ALLOCATE(lctmp(4,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    mesh%lc(:,mesh%nc+1:mesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==2)THEN
    ALLOCATE(lctmp(10,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    mesh%lc(:,mesh%nc+1:mesh%nc+nc)=lctmp(1:4,:)
    lc_ho(:,mesh%nc+1:mesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==3)THEN
    ALLOCATE(lctmp(8,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    mesh%lc(:,mesh%nc+1:mesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==4)THEN
    ALLOCATE(lctmp(27,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    mesh%lc(:,mesh%nc+1:mesh%nc+nc)=lctmp(1:8,:)
    lc_ho(:,mesh%nc+1:mesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  END IF
  mesh%reg(mesh%nc+1:mesh%nc+nc)=i
  mesh%nc=mesh%nc+nc
end do
IF(per_ns>0)THEN
  !---
  if(per_ns<10)then
    write(blknm,'(I1)')per_ns
  else if(per_ns>9.AND.per_ns<100)then
    write(blknm,'(I2)')per_ns
  else if(per_ns>=100)then
    CALL oft_abort('Invalid Exodus variable index (num_nod_nsXX)','mesh_cubit_read_regions',__FILE__)
  end if
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_nod_ns"//TRIM(blknm),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=np_per))
  ALLOCATE(per_nodes(np_per))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"node_ns"//TRIM(blknm),lcID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,per_nodes))
END IF
CALL oft_decrease_indent
end subroutine mesh_cubit_read_regions
end subroutine mesh_cubit_load
!------------------------------------------------------------------------------
! SUBROUTINE smesh_cubit_load
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine smesh_cubit_load()
logical :: read_2d = .FALSE.
logical :: ltrans(3)
real(8), allocatable :: rtmp(:),r2dtmp(:,:)
integer(4) :: i,j,ncid,blkID,nodesID,nblks,elemID,lcID,ierr
integer(4) :: nf,np,io_unit,nsID,num_nsets,num_sidesets,att_len
integer(i4), allocatable :: lptmp(:)
character(LEN=3) :: blknm
character(LEN=10) :: eltype
!---Read in mesh options
namelist/cubit_options/filename,inpname,lf_file,tor_mesh, &
  reflect,per_ns,zstretch,tor_rmin
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,cubit_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "cubit_options" found in input file.','smesh_cubit_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "cubit_options" in input file.','smesh_cubit_load',__FILE__)
  IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','smesh_cubit_load',__FILE__)
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'**** Loading CUBIT surface mesh'
  CALL oft_increase_indent
  WRITE(*,'(3A)')oft_indent,'Mesh File = ',TRIM(filename)
ELSE
  CALL oft_increase_indent
END IF
!---Broadcast input information
#ifdef HAVE_MPI
ltrans=(/lf_file,tor_mesh,reflect/)
CALL MPI_Bcast(filename,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_cubit_load',__FILE__)
CALL MPI_Bcast(inpname,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_cubit_load',__FILE__)
CALL MPI_Bcast(ltrans,3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_cubit_load',__FILE__)
CALL MPI_Bcast(per_ns,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_cubit_load',__FILE__)
CALL MPI_Bcast(zstretch,1,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','smesh_cubit_load',__FILE__)
lf_file=ltrans(1)
tor_mesh=ltrans(2)
reflect=ltrans(3)
#endif
!---Open mesh file
call mesh_cubit_error(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
!---Check mesh type
call mesh_cubit_error(NF90_inq_varid(ncid,"connect1",blkID))
call mesh_cubit_error(NF90_INQUIRE_ATTRIBUTE(ncid,blkID,"elem_type",len=att_len))
call mesh_cubit_error(NF90_GET_ATT(ncid,blkID,"elem_type",eltype))
!---Allocate mesh
IF(eltype(1:3)=="TRI")THEN
  ALLOCATE(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
ELSE IF(eltype(1:4)=="QUAD")THEN
  ALLOCATE(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
END IF
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%smeshes(i)%setup(mesh_cubit_id,.FALSE.)
END DO
smesh=>mg_mesh%smeshes(1)
!---Read vertex count
call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_nodes",nodesID))
call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nodesID,len = smesh%np))
call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_blk",blkID))
call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,blkID,len = nblks))
!---Read vertices
allocate(smesh%r(3,smesh%np))
ierr=NF90_INQ_VARID(ncid,"coordx",nodesID)
lf_file=(ierr==nf90_noerr)
IF(lf_file)THEN
  allocate(rtmp(smesh%np))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coordx",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp))
  smesh%r(1,:)=rtmp
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coordy",nodesID))
  call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp))
  smesh%r(2,:)=rtmp
  ierr=NF90_INQ_VARID(ncid,"coordz",nodesID)
  IF(ierr==nf90_noerr)THEN
    call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,rtmp))
    smesh%r(3,:)=rtmp
  ELSE
    smesh%r(3,:)=0.d0
  END IF
  deallocate(rtmp)
ELSE
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_dim",blkID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,blkID,len = np))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"coord",nodesID))
  IF(np==2)THEN
    allocate(r2dtmp(smesh%np,2))
    call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,r2dtmp))
    DO i=1,smesh%np
      smesh%r(1:2,i)=r2dtmp(i,:)
      smesh%r(3,i)=0.d0
    END DO
  ELSE IF(np==3)THEN
    allocate(r2dtmp(smesh%np,3))
    call mesh_cubit_error(NF90_GET_VAR(ncid,nodesID,r2dtmp))
    DO i=1,smesh%np
      smesh%r(:,i)=r2dtmp(i,:)
    END DO
  ELSE
    CALL oft_abort("Could not determine vertex list size","",__FILE__)
  END IF
  deallocate(r2dtmp)
END IF
IF(ALL(ABS(smesh%r(3,:))<1.d-10))smesh%dim=2
!---Read cell lists
CALL read_regions
!---Close file
call mesh_cubit_error(NF90_CLOSE(ncid))
!---Remove high order points from vertex list
IF(have_ho)THEN
  ALLOCATE(r_ho(3,smesh%np))
  r_ho=smesh%r
  np_ho=smesh%np
  !---Index points used in cell vertices
  ALLOCATE(lptmp(smesh%np))
  lptmp=0
  smesh%np=0
  DO i=1,smesh%nc
    DO j=1,smesh%cell_np
      IF(lptmp(smesh%lc(j,i))==0)THEN
        smesh%np = smesh%np + 1
        lptmp(smesh%lc(j,i)) = smesh%np
      END IF
    END DO
  END DO
  !---Reset point list
  DEALLOCATE(smesh%r)
  ALLOCATE(smesh%r(3,smesh%np))
  DO i=1,np_ho
    IF(lptmp(i)/=0)smesh%r(:,lptmp(i)) = r_ho(:,i)
  END DO
  !---Reindex cells
  DO i=1,smesh%nc
    DO j=1,smesh%cell_np
      smesh%lc(j,i) = lptmp(smesh%lc(j,i))
      lc_ho(j,i) = lptmp(lc_ho(j,i))
    END DO
  END DO
! #ifdef HAVE_ONURBS
!   !---Reindex CAD edges
!   DO i=1,ngmc
!     DO j=1,model_curves(i)%ne
!       model_curves(i)%lgme(:,j) = lptmp(model_curves(i)%lgme(:,j))
!     END DO
!   END DO
! #endif
  DEALLOCATE(lptmp)
  r_ho(3,:)=r_ho(3,:)*zstretch
END IF
smesh%r(3,:)=smesh%r(3,:)*zstretch
! IF(reflect)THEN
!   call mesh_global_resolution(mesh)
!   call mesh_cubit_reflect(.1d0*mesh%hmin)
! END IF
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,I8)')oft_indent,'# of Regions      = ',nblks
#ifdef HAVE_ONURBS
  WRITE(*,'(2A,2I8)')oft_indent,'# of CAD curves   = ',ngwc,ngmc
#else
  WRITE(*,'(2A)')oft_indent,'Skipping CAD information: Not compiled with OpenNURBS'
#endif
END IF
CALL oft_decrease_indent
call mesh_global_resolution(smesh)
IF(oft_env%rank/=0)DEALLOCATE(smesh%r,smesh%lc,smesh%reg)
contains
!---Formatting sub-function
function format_blknum(ind) result(blknum)
integer(4), intent(in) :: ind
character(LEN=3) :: blknum
if(ind<10)then
  write(blknum,'(I1)')ind
else if((ind>9).AND.(ind<100))then
  write(blknum,'(I2)')ind
else if((ind>99).AND.(ind<1000))then
  write(blknum,'(I3)')ind
end if
end function format_blknum
!------------------------------------------------------------------------------
! SUBROUTINE: read_regions
!------------------------------------------------------------------------------
!> Read in Exodus surface mesh
!! - Merge all blocks
!------------------------------------------------------------------------------
subroutine read_regions
integer(i4) :: blkID,elemID,lcID,att_len
integer(i4) :: i,j,nc
integer(i4), allocatable :: rmark(:),lctmp(:,:)
character(3) :: blknm
character(40) :: eltype
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Loading Cubit regions'
CALL oft_increase_indent
!---
allocate(rmark(nblks))
rmark=0
nregions=0
smesh%nc=0
!---
do i=1,nblks
  IF(oft_debug_print(2))WRITE(*,'(3A)')oft_indent,'Loading cell count for region ',TRIM(format_blknum(i))
  !---
  call mesh_cubit_error(NF90_inq_varid(ncid,"connect"//TRIM(format_blknum(i)),blkID))
  call mesh_cubit_error(NF90_INQUIRE_ATTRIBUTE(ncid,blkID,"elem_type",len=att_len))
  !---
  call mesh_cubit_error(NF90_GET_ATT(ncid,blkID,"elem_type",eltype))
  !---
  if(eltype(1:4)=="TRI3")then
    rmark(i)=1
  else if(eltype(1:4)=="TRI6")then
    rmark(i)=2
    have_ho=.TRUE.
  else if(eltype(1:5)=="QUAD4")then
    rmark(i)=3
  else if(eltype(1:5)=="QUAD9")then
    rmark(i)=4
    have_ho=.TRUE.
  ELSE
    if(eltype(1:3)=="TRI")then
      rmark(i)=1
    else if(eltype(1:4)=="QUAD")then
      rmark(i)=3
    else
      CYCLE
    endif
  end if
  nregions=nregions+1
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_in_blk"//TRIM(format_blknum(i)),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=nc))
  smesh%nc=smesh%nc+nc
end do
IF(smesh%nc==0)CALL oft_abort('No volume regions found!','smesh::read_regions',__FILE__)
IF(ANY(rmark==1).AND.ANY(rmark==3))CALL oft_abort('Two conflicting element types found', &
  'smesh::read_regions',__FILE__)
!---
allocate(smesh%lc(smesh%cell_np,smesh%nc),smesh%reg(smesh%nc))
IF(have_ho)THEN
  IF(MAXVAL(rmark)==2)THEN
    ALLOCATE(lc_ho(6,smesh%nc))
  ELSE IF(MAXVAL(rmark)==4)THEN
    ALLOCATE(lc_ho(9,smesh%nc))
  END IF
  lc_ho=0
END IF
!---
smesh%nc=0
do i=1,nblks
  if(rmark(i)==0)cycle
  IF(oft_debug_print(2))WRITE(*,'(3A)')oft_indent,'Loading cell lists for region ',TRIM(format_blknum(i))
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_in_blk"//TRIM(format_blknum(i)),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=nc))
  call mesh_cubit_error(NF90_INQ_VARID(ncid,"connect"//TRIM(format_blknum(i)),lcID))
  IF(rmark(i)==1)THEN
    ALLOCATE(lctmp(3,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    smesh%lc(:,smesh%nc+1:smesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==2)THEN
    ALLOCATE(lctmp(6,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    smesh%lc(:,smesh%nc+1:smesh%nc+nc)=lctmp(1:3,:)
    lc_ho(:,smesh%nc+1:smesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==3)THEN
    ALLOCATE(lctmp(4,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    smesh%lc(:,smesh%nc+1:smesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  ELSE IF(rmark(i)==4)THEN
    ALLOCATE(lctmp(9,nc))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,lctmp))
    smesh%lc(:,smesh%nc+1:smesh%nc+nc)=lctmp(1:4,:)
    lc_ho(:,smesh%nc+1:smesh%nc+nc)=lctmp
    DEALLOCATE(lctmp)
  END IF
  smesh%reg(smesh%nc+1:smesh%nc+nc)=i
  smesh%nc=smesh%nc+nc
end do
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine read_regions
end subroutine smesh_cubit_load
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine cubit_read_nodesets(nsets,cubit_filename)
TYPE(oft_1d_int), pointer, intent(inout) :: nsets(:)
CHARACTER(LEN=OFT_PATH_SLEN), OPTIONAL, INTENT(in) :: cubit_filename
integer(4) :: i,j,ncid,nsID,num_nsets
!---Open mesh file
IF(PRESENT(cubit_filename))THEN
  call mesh_cubit_error(NF90_OPEN(TRIM(cubit_filename),NF90_NOWRITE,ncid))
ELSE
  call mesh_cubit_error(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
END IF
!---Read nodesets
i=NF90_INQ_DIMID(ncid,"num_node_sets",nsID)
IF(i==nf90_noerr)THEN
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nsID,len = num_nsets))
  ALLOCATE(nsets(num_nsets))
  DO j=1,num_nsets
    i=NF90_INQ_DIMID(ncid,"num_nod_ns"//TRIM(format_blknum(j)),nsID)
    IF(i==nf90_noerr)THEN
      call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nsID,len = nsets(j)%n))
      ALLOCATE(nsets(j)%v(nsets(j)%n))
      call mesh_cubit_error(NF90_INQ_VARID(ncid,"node_ns"//TRIM(format_blknum(j)),nsID))
      call mesh_cubit_error(NF90_GET_VAR(ncid,nsID,nsets(j)%v))
    END IF
  END DO
END IF
!---Close file
call mesh_cubit_error(NF90_CLOSE(ncid))
CALL oft_decrease_indent
contains
!---Formatting sub-function
function format_blknum(ind) result(blknum)
integer(4), intent(in) :: ind
character(LEN=3) :: blknum
if(ind<10)then
  write(blknum,'(I1)')ind
else if((ind>9).AND.(ind<100))then
  write(blknum,'(I2)')ind
else if((ind>99).AND.(ind<1000))then
  write(blknum,'(I3)')ind
end if
end function format_blknum
END SUBROUTINE cubit_read_nodesets
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine cubit_read_sidesets(ssets,cubit_filename)
TYPE(oft_1d_int), pointer, intent(inout) :: ssets(:)
CHARACTER(LEN=OFT_PATH_SLEN), OPTIONAL, INTENT(in) :: cubit_filename
integer(4) :: i,j,ncid,nsID,num_sidesets
!---Open mesh file
IF(PRESENT(cubit_filename))THEN
  call mesh_cubit_error(NF90_OPEN(TRIM(cubit_filename),NF90_NOWRITE,ncid))
ELSE
  call mesh_cubit_error(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
END IF
!---Read sidesets
i=NF90_INQ_DIMID(ncid,"num_side_sets",nsID)
IF(i==nf90_noerr)THEN
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nsID,len = num_sidesets))
  ALLOCATE(ssets(num_sidesets))
  DO j=1,num_sidesets
    i=NF90_INQ_DIMID(ncid,"num_side_ss"//TRIM(format_blknum(j)),nsID)
    IF(i==nf90_noerr)THEN
      call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,nsID,len = ssets(j)%n))
      ALLOCATE(ssets(j)%v(ssets(j)%n))
      call mesh_cubit_error(NF90_INQ_VARID(ncid,"elem_ss"//TRIM(format_blknum(j)),nsID))
      call mesh_cubit_error(NF90_GET_VAR(ncid,nsID,ssets(j)%v))
    END IF
  END DO
END IF
!---Close file
call mesh_cubit_error(NF90_CLOSE(ncid))
CALL oft_decrease_indent
contains
!---Formatting sub-function
function format_blknum(ind) result(blknum)
integer(4), intent(in) :: ind
character(LEN=3) :: blknum
if(ind<10)then
  write(blknum,'(I1)')ind
else if((ind>9).AND.(ind<100))then
  write(blknum,'(I2)')ind
else if((ind>99).AND.(ind<1000))then
  write(blknum,'(I3)')ind
end if
end function format_blknum
END SUBROUTINE cubit_read_sidesets
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_read_surface
!------------------------------------------------------------------------------
!> Read in Exodus surface mesh
!! - Surface mesh used for indexing CAD geometry
!------------------------------------------------------------------------------
#ifdef HAVE_ONURBS
subroutine mesh_cubit_read_surface
integer(i4) :: blkID,elemID,lcID,att_len,ebID,qaID
integer(i4) :: i,j,nc
integer(i4), allocatable :: rmark(:)
logical :: valid_link
character(3) :: blknm
character(40) :: eltype
character(LEN=5) :: id
CHARACTER(LEN=exodus_string_len) :: ntmp
CHARACTER(LEN=exodus_string_len), ALLOCATABLE :: block_names(:),qa_records(:,:)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Loading Cubit surfaces'
!---
ALLOCATE(block_names(nblks))
call mesh_cubit_error(NF90_inq_varid(ncid,"eb_names",ebID))
call mesh_cubit_error(NF90_GET_VAR(ncid,ebID,block_names))
!---Read CUBIT version
ALLOCATE(qa_records(4,1))
call mesh_cubit_error(NF90_inq_varid(ncid,"qa_records",qaID))
call mesh_cubit_error(NF90_GET_VAR(ncid,qaID,qa_records))
READ(qa_records(2,1),'(F4.1)')cubit_version
IF(oft_debug_print(1))WRITE(*,'(4X,A,F4.1)')'Generated by CUBIT v',cubit_version
!---
allocate(rmark(nblks))
rmark=0
ngmc=0
ngms=0
!---
do i=1,nblks
  !---
  if(i<10)then
    write(blknm,'(I1)')i
  else if(i>9.AND.i<100)then
    write(blknm,'(I2)')i
  else if(i>99.AND.i<1000)then
    write(blknm,'(I3)')i
  else if(i>=1000)then
    CALL oft_abort('Invalid Exodus variable index (connectXXX)','mesh_cubit_read_surface',__FILE__)
  end if
  !---
  call mesh_cubit_error(NF90_inq_varid(ncid,"connect"//TRIM(blknm),blkID))
  call mesh_cubit_error(NF90_INQUIRE_ATTRIBUTE(ncid,blkID,"elem_type",len=att_len))
  !---
  call mesh_cubit_error(NF90_GET_ATT(ncid,blkID,"elem_type",eltype))
  !---
  ntmp=block_names(i)
  if(eltype(1:4)=="BAR2")then
    rmark(i)=1
    ngmc=ngmc+1
    !---
    valid_link=(ntmp(1:4)=='TBC_')
    IF(.NOT.valid_link)THEN
      WRITE(*,*)'======================'
      WRITE(*,*)'ERROR: Non-curve geometry linked to "BAR" region'
      WRITE(*,*)'  Block: ',i
      WRITE(*,*)'  Geom:  ',ntmp
      CALL oft_abort('Invalid Cubit geometry','mesh_cubit_read_surface',__FILE__)
    END IF
  else if(eltype(1:4)=="TRI3")then
    rmark(i)=2
    ngms=ngms+1
    !---
    IF(cubit_version>14.05)THEN
      valid_link=(ntmp(1:5)=='TBST_')
    ELSE
      valid_link=(ntmp(1:4)=='TBS_')
    END IF
    IF(.NOT.valid_link)THEN
      WRITE(*,*)'======================'
      WRITE(*,*)'ERROR: Non-surface geometry linked to "TRI" region'
      WRITE(*,*)'  Block: ',i
      WRITE(*,*)'  Geom:  ',ntmp
      CALL oft_abort('Invalid Cubit geometry','mesh_cubit_read_surface',__FILE__)
    END IF
  else if(eltype(1:5)=="QUAD4")then
    rmark(i)=3
    ngms=ngms+1
    !---
    IF(cubit_version>14.05)THEN
      valid_link=(ntmp(1:5)=='TBST_')
    ELSE
      valid_link=(ntmp(1:4)=='TBS_')
    END IF
    IF(.NOT.valid_link)THEN
      WRITE(*,*)'======================'
      WRITE(*,*)'ERROR: Non-surface geometry linked to "QUAD" region'
      WRITE(*,*)'  Block: ',i
      WRITE(*,*)'  Geom:  ',ntmp
      CALL oft_abort('Invalid Cubit geometry','mesh_cubit_read_surface',__FILE__)
    END IF
  end if
end do
!---
IF(reflect)THEN
  allocate(model_curves(2*ngmc),model_surfaces(2*ngms))
ELSE
  allocate(model_curves(ngmc),model_surfaces(ngms))
END IF
ngmc=0
ngms=0
!---
do i=1,nblks
  if(rmark(i)==0)cycle
  !---
  if(i<10)then
    write(blknm,'(I1)')i
  else if(i>9.AND.i<100)then
    write(blknm,'(I2)')i
  else if(i>99.AND.i<1000)then
    write(blknm,'(I3)')i
  end if
  !---
  call mesh_cubit_error(NF90_INQ_DIMID(ncid,"num_el_in_blk"//TRIM(blknm),elemID))
  call mesh_cubit_error(NF90_INQUIRE_DIMENSION(ncid,elemID,len=nc))
  if(rmark(i)==1)then
    ngmc=ngmc+1
    model_curves(ngmc)%ne=nc
    allocate(model_curves(ngmc)%lgme(2,model_curves(ngmc)%ne))
    !---
    call mesh_cubit_error(NF90_INQ_VARID(ncid,"connect"//TRIM(blknm),lcID))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,model_curves(ngmc)%lgme))
    model_curves(ngmc)%name=block_names(i)
    id=model_curves(ngmc)%name(ex_topc_len:ex_topc_len+4)
    !---Replace NULL characters
    DO j=1,5
      IF(id(j:j)==c_null_char)id(j:j)=""
    END DO
    read(id,'(I5)')j
    model_curves(ngmc)%top_cid=j
  else if((rmark(i)==2).OR.(rmark(i)==3))then
    ngms=ngms+1
    model_surfaces(ngms)%nf=nc
    allocate(model_surfaces(ngms)%lgmf(rmark(i)+1,model_surfaces(ngms)%nf))
    !---
    call mesh_cubit_error(NF90_INQ_VARID(ncid,"connect"//TRIM(blknm),lcID))
    call mesh_cubit_error(NF90_GET_VAR(ncid,lcID,model_surfaces(ngms)%lgmf))
    model_surfaces(ngms)%name=block_names(i)
    IF(cubit_version>14.05)THEN
      id=model_surfaces(ngms)%name(ex_tops_len+1:ex_tops_len+5)
    ELSE
      id=model_surfaces(ngms)%name(ex_tops_len:ex_tops_len+4)
    END IF
    !---Replace NULL characters
    DO j=1,5
      IF(id(j:j)==c_null_char)id(j:j)=""
    END DO
    read(id,'(I5)')j
    model_surfaces(ngms)%top_sid=j
  end if
end do
DEBUG_STACK_POP
end subroutine mesh_cubit_read_surface
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_geom
!------------------------------------------------------------------------------
!> Load and initialize OpenNURBS geometry for current Exodus mesh
!! - Initialize OpenNURBS library
!! - Construct FORTRAN objects for each CAD object
!------------------------------------------------------------------------------
subroutine mesh_cubit_geom
integer(i4) :: i,j,ierr,p1,p2,s(4)
real(r8) :: pt1(3),pt2(3)
character(len=5) :: id
character(kind=c_char,len=21) :: tmp_string
character(kind=c_char,len=31) :: obj_name
DEBUG_STACK_PUSH
IF(cubit_version<13.95)CALL oft_abort('CUBIT version too old, v14.0+ required', &
'mesh_cubit_geom',__FILE__)
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Loading OpenNURBS geometry:'
CALL oft_increase_indent
!---
call nurbs_init(ierr)
!---
tmp_string=TRIM(inpname)//C_NULL_CHAR
call nurbs_read_in(tmp_string,ierr)
!---
call nurbs_get_count(ngwc,ngws,ierr)
!---
IF(reflect)THEN
  allocate(wf_curves(2*ngwc),wf_surfs(2*ngws))
ELSE
  allocate(wf_curves(ngwc),wf_surfs(ngws))
END IF
do i=1,ngwc
  wf_curves(i)%coord_scales(3)=zstretch
  wf_curves(i)%id=i-1
  call nurbs_curve_domain(wf_curves(i)%id,wf_curves(i)%domain,ierr)
  IF(ierr<0)CALL oft_abort('Error getting curve domain','mesh_cubit_geom',__FILE__)
  !---Test if curve is periodic
  call nurbs_curve_periodic(wf_curves(i)%id,p1,ierr)
  IF(ierr<0)CALL oft_abort('Error getting curve periodicity','mesh_cubit_geom',__FILE__)
  wf_curves(i)%periodic=(p1==1)
  !---Get physical end points
  CALL wf_curves(i)%eval(pt1,wf_curves(i)%domain(1),0.d0)
  CALL wf_curves(i)%eval(pt2,wf_curves(i)%domain(2),0.d0)
  !---Force periodic for complex curves with closed geometry
  IF(.NOT.wf_curves(i)%periodic)THEN
    IF(SQRT(SUM((pt1-pt2)**2))<1.d-8)wf_curves(i)%periodic=.TRUE.
  END IF
  !---Test if curve is straight
  call nurbs_curve_linear(wf_curves(i)%id,p1,ierr)
  IF(ierr<0)CALL oft_abort('Error getting curve linearity','mesh_cubit_geom',__FILE__)
  wf_curves(i)%linear=(p1==1)
  call wf_curves(i)%grid
  !---
  obj_name=' '
  obj_name(31:31)=C_NULL_CHAR
  call nurbs_curve_name(wf_curves(i)%id,obj_name,ierr)
  IF(ierr<0)CALL oft_abort('Error getting curve name','mesh_cubit_geom',__FILE__)
  !wf_curves(i)%name=obj_name(1:ex_topc_len+5)
  id=obj_name(ex_topc_len:ex_topc_len+5)
  read(id,'(I5)')j
  wf_curves(i)%cid=j
  IF(oft_debug_print(3))THEN
    WRITE(*,'(2A)')oft_indent,'============================'
    WRITE(*,'(3A)')oft_indent,'Found CAD Curve: id = ',TRIM(obj_name)
    WRITE(*,'(2A,L)')oft_indent,'Linear: ',wf_curves(i)%linear
    WRITE(*,'(2A,L)')oft_indent,'Periodic: ',wf_curves(i)%periodic
    WRITE(*,'(2A,2ES12.4)')oft_indent,'Domain: ',wf_curves(i)%domain
    WRITE(*,'(2A,3ES12.4)')oft_indent,'Start Point: ',pt1
    WRITE(*,'(2A,3ES12.4)')oft_indent,'End Point:   ',pt2
    WRITE(*,'(2A)')oft_indent,'============================'
  END IF
end do
do i=1,ngws
  wf_surfs(i)%coord_scales(3)=zstretch
  wf_surfs(i)%id=i-1
  call nurbs_surf_domain(wf_surfs(i)%id,wf_surfs(i)%domain(:,1),wf_surfs(i)%domain(:,2),ierr)
  call nurbs_surf_periodic(wf_surfs(i)%id,p1,p2,ierr)
  wf_surfs(i)%periodic(1)=(p1==1)
  wf_surfs(i)%periodic(2)=(p2==1)
  call nurbs_surf_singular(wf_surfs(i)%id,s,ierr)
  wf_surfs(i)%singular(1,1)=(s(1)==1)
  wf_surfs(i)%singular(2,1)=(s(3)==1)
  wf_surfs(i)%singular(1,2)=(s(2)==1)
  wf_surfs(i)%singular(2,2)=(s(4)==1)
  call nurbs_surf_planar(wf_surfs(i)%id,p1,ierr)
  wf_surfs(i)%planar=(p1==1)
  call wf_surfs(i)%grid
  !---
  obj_name=' '
  obj_name(31:31)=C_NULL_CHAR
  call nurbs_surf_name(wf_surfs(i)%id,obj_name,ierr)
  !wf_surfs(i)%name=obj_name
  IF(cubit_version>14.05)THEN
    id=obj_name(ex_tops_len+1:ex_tops_len+6)
  ELSE
    id=obj_name(ex_tops_len:ex_tops_len+5)
  END IF
  read(id,'(I5)')j
  wf_surfs(i)%sid=j
  IF(oft_debug_print(3))THEN
    WRITE(*,'(2A)')oft_indent,'============================'
    WRITE(*,'(3A)')oft_indent,'Found CAD surface: id = ',TRIM(obj_name)
    WRITE(*,'(2A,L)')oft_indent,'Planar: ',wf_surfs(i)%planar
    WRITE(*,'(2A,2L)')oft_indent,'Periodic: ',wf_surfs(i)%periodic
    WRITE(*,'(2A,4L)')oft_indent,'Singular: ',wf_surfs(i)%singular(:,1),wf_surfs(i)%singular(:,2)
    WRITE(*,'(2A,4ES12.4)')oft_indent,'Domain: ',wf_surfs(i)%domain(:,1),wf_surfs(i)%domain(:,2)
    WRITE(*,'(2A)')oft_indent,'============================'
  END IF
end do
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_cubit_geom
#endif
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_cadlink
!------------------------------------------------------------------------------
!> Link OpenNURBS CAD objects to Exodus mesh entities for use in refinement.
!! - Map Exodus surface mesh to CAD objects using name attributes
!------------------------------------------------------------------------------
subroutine mesh_cubit_cadlink
#ifdef HAVE_ONURBS
real(r8) :: pt(3),u,v
integer(i4) :: i,j,m,ind,k,ep(2),fp(3),etmp(2),ftmp(3),ierr(4),chk
integer(i4), allocatable, dimension(:) :: fmap,lce_tmp,lcf_tmp
IF(TRIM(inpname)=='none')RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2A)')oft_indent,'Linking OpenNURBS geometry to mesh'
CALL oft_increase_indent
!---Allocate preliminary linkages
ALLOCATE(lce_tmp(mesh%ne),lcf_tmp(mesh%nf))
lce_tmp=0
lcf_tmp=0
cad_link%nce=0
cad_link%ncf=0
!---Loop over bar sets and link to CAD curves
DO i=1,ngmc
  chk=0
  DO j=1,ngwc
    IF(wf_curves(j)%cid/=model_curves(i)%top_cid)CYCLE
    IF(oft_debug_print(3))WRITE(*,*)model_curves(i)%top_cid,wf_curves(j)%cid
    model_curves(i)%cid=j
    !---Attempt to find edge vertices on current CAD curve
    etmp=model_curves(i)%lgme(:,1)
    pt=mesh%r(:,etmp(1))
    CALL wf_curves(j)%locate(pt,u,v,ierr(1))
    pt=mesh%r(:,etmp(2))
    CALL wf_curves(j)%locate(pt,u,v,ierr(2))
    IF(oft_debug_print(3))WRITE(*,*)i,j,ierr(1:2)
    IF(ALL(ierr(1:2)==0))THEN
      IF(oft_debug_print(3))WRITE(*,*)'Bar set ',INT(i,2),' linked with curve ',INT(j,2)
      IF(.NOT.wf_curves(j)%linear)THEN
        DO m=1,model_curves(i)%ne
          etmp=model_curves(i)%lgme(:,m)
          ind=ABS(mesh_local_findedge(mesh,etmp)) ! Find edge on mesh'
          IF(ind==0)CYCLE ! Not on my domain
          lce_tmp(ind)=j
          cad_link%nce=cad_link%nce+1
        END DO
      END IF
      chk=1
    END IF
  END DO
  IF(chk==0.AND.oft_debug_print(3))WRITE(*,'(2A,I6,A)')oft_indent,'Bar set ',i,' unlinked'
END DO
!---Get face boundary mapping
ALLOCATE(fmap(mesh%nf))
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf)
!---Loop over triangle/quad sets and link to CAD surfaces
DO i=1,ngms
  chk=0
  !---Check CAD surfaces for match
  DO j=1,ngws
    IF(wf_surfs(j)%sid/=model_surfaces(i)%top_sid)CYCLE
    model_surfaces(i)%sid=j
    !---Attempt to find face vertices on current CAD surface
    DO m=1,mesh%face_np
      pt=mesh%r(:,model_surfaces(i)%lgmf(m,1))
      CALL wf_surfs(j)%locate(pt,u,v,ierr(m))
    END DO
    IF(oft_debug_print(3))WRITE(*,*)i,j,ierr(1:mesh%face_np)
    !---If matched link remaining faces
    IF(ALL(ierr(1:mesh%face_np)==0))THEN
      IF(oft_debug_print(3))WRITE(*,'(2A,I6,A,I6)')oft_indent,'Tri set ',i, &
        ' linked with surface ',wf_surfs(j)%id
      DO m=1,model_surfaces(i)%nf
        ind=abs(mesh_local_findface(mesh,model_surfaces(i)%lgmf(:,m))) ! Find face on mesh
        if(ind==0)cycle ! Not on my domain
        lcf_tmp(ind) = cubit_soffset + j
        cad_link%ncf = cad_link%ncf + 1
        if(fmap(ind)==0)cycle
        mesh%bfs(fmap(ind))=j
      END DO
      chk=1
      EXIT
    END IF
  END DO
  IF(chk==0.AND.oft_debug_print(3))WRITE(*,'(2A,I6,A)')oft_indent,'Tri set ',i,' unlinked'
END DO
DEALLOCATE(fmap) ! Destroy mapping array
!---Link faces to CAD objects and propogate to edges
ALLOCATE(cad_link%lcf(cad_link%ncf),cad_link%lbfg(cad_link%ncf))
cad_link%ncf=0
do i=1,mesh%nf
  if(lcf_tmp(i)==0)cycle ! Face is not linked
  cad_link%ncf = cad_link%ncf + 1
  cad_link%lcf(cad_link%ncf) = i
  cad_link%lbfg(cad_link%ncf)%wo => wf_surfs(lcf_tmp(i)-cubit_soffset)
  !---Propogate to edges
  do j=1,mesh%face_np
    k=ABS(mesh%lfe(j,i))
    IF(lce_tmp(k)==0)THEN
      lce_tmp(k) = lcf_tmp(i)
      cad_link%nce = cad_link%nce+1
    END IF
  enddo
enddo
DEALLOCATE(lcf_tmp)
!---Link edges to CAD objects
ALLOCATE(cad_link%lce(cad_link%nce),cad_link%lbeg(cad_link%nce))
cad_link%nce=0
DO i=1,mesh%ne
  if(lce_tmp(i)==0)cycle ! Edge is not linked
  cad_link%nce = cad_link%nce + 1
  cad_link%lce(cad_link%nce) = i
  IF(lce_tmp(i)<cubit_soffset)THEN
    cad_link%lbeg(cad_link%nce)%wo => wf_curves(lce_tmp(i))
  ELSE
    cad_link%lbeg(cad_link%nce)%wo => wf_surfs(lce_tmp(i)-cubit_soffset)
  END IF
END DO
DEALLOCATE(lce_tmp)
CALL oft_decrease_indent
DEBUG_STACK_POP
#endif
end subroutine mesh_cubit_cadlink
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_hobase
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from high order import
!------------------------------------------------------------------------------
subroutine mesh_cubit_hobase(self)
CLASS(oft_amesh), INTENT(inout) :: self
IF(.NOT.have_ho)RETURN
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Importing quadratic mesh nodes'
SELECT CASE(self%type)
  CASE(1) ! Tet/Tri
    SELECT TYPE(self)
    CLASS IS(oft_bmesh)
      CALL mesh_cubit_hobase_tri(self)
    CLASS IS(oft_mesh)
      CALL mesh_cubit_hobase_tet(self)
    END SELECT
  CASE(3) ! Hex/Quad
    SELECT TYPE(self)
    CLASS IS(oft_bmesh)
      CALL mesh_cubit_hobase_quad(self)
    CLASS IS(oft_mesh)
      CALL mesh_cubit_hobase_hex(self)
    END SELECT
  CASE DEFAULT
    CALL oft_abort("Unknown element type", "mesh_cubit_hobase", __FILE__)
END SELECT
end subroutine mesh_cubit_hobase
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_hobase_tet
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from TRI6 import
!------------------------------------------------------------------------------
subroutine mesh_cubit_hobase_tri(self)
CLASS(oft_bmesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2)
integer(i4), allocatable :: emark(:)
DEBUG_STACK_PUSH
!---Setup quadratic mesh
self%order=1
self%ho_info%nep=1
ALLOCATE(self%ho_info%r(3,self%ne),self%ho_info%lep(self%ho_info%nep,self%ne))
!---Initialize high order points with straight edges
DO i=1,self%ne
  self%ho_info%lep(1,i)=i
  self%ho_info%r(:,i)=(self%r(:,self%le(1,i))+self%r(:,self%le(2,i)))/2.d0
END DO
!---Set midpoints from imported list
ALLOCATE(emark(self%ne))
emark=0
DO i=1,self%nc
  DO j=1,self%cell_ne
    etmp=lc_ho(exodus_tri_emap(:,j),i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','mesh_cubit_hobase_tri',__FILE__)
    IF(emark(k)==0)THEN
      self%ho_info%r(:,k) = r_ho(:,lc_ho(3+j,i))
      emark(k)=1
    END IF
  END DO
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,lc_ho,emark)
DEBUG_STACK_POP
end subroutine mesh_cubit_hobase_tri
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_hobase_tet
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from TETRA10 import
!------------------------------------------------------------------------------
subroutine mesh_cubit_hobase_tet(self)
CLASS(oft_mesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2)
integer(i4), allocatable :: emark(:)
DEBUG_STACK_PUSH
!---Setup quadratic mesh
self%order=1
self%ho_info%nep=1
ALLOCATE(self%ho_info%r(3,self%ne),self%ho_info%lep(self%ho_info%nep,self%ne))
!---Initialize high order points with straight edges
DO i=1,self%ne
  self%ho_info%lep(1,i)=i
  self%ho_info%r(:,i)=(self%r(:,self%le(1,i))+self%r(:,self%le(2,i)))/2.d0
END DO
!---Set midpoints from imported list
ALLOCATE(emark(self%ne))
emark=0
DO i=1,self%nc
  DO j=1,self%cell_ne
    etmp=lc_ho(exodus_tet_emap(:,j),i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','mesh_cubit_hobase_tet',__FILE__)
    IF(emark(k)==0)THEN
      self%ho_info%r(:,k) = r_ho(:,lc_ho(4+j,i))
      emark(k)=1
    END IF
  END DO
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,lc_ho,emark)
DEBUG_STACK_POP
end subroutine mesh_cubit_hobase_tet
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_hobase_quad
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from QUAD9 import
!------------------------------------------------------------------------------
subroutine mesh_cubit_hobase_quad(self)
CLASS(oft_bmesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2)
integer(i4), allocatable :: emark(:)
DEBUG_STACK_PUSH
!---Setup quadratic mesh
self%order=1
self%ho_info%nep=1
self%ho_info%ncp=1
ALLOCATE(self%ho_info%r(3,self%ne+self%nc),self%ho_info%lep(self%ho_info%nep,self%ne))
ALLOCATE(self%ho_info%lcp(self%ho_info%ncp,self%nc))
!---Initialize high order points with straight edges
DO i=1,self%ne
  self%ho_info%lep(1,i)=i
  self%ho_info%r(:,i)=(self%r(:,self%le(1,i))+self%r(:,self%le(2,i)))/2.d0
END DO
!---Set midpoints from imported list
ALLOCATE(emark(self%ne))
emark=0
DO i=1,self%nc
  DO j=1,self%cell_ne
    etmp=lc_ho(exodus_quad_emap(:,j),i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','mesh_cubit_hobase_quad',__FILE__)
    IF(emark(k)==0)THEN
      self%ho_info%r(:,k) = r_ho(:,lc_ho(4+j,i))
      emark(k)=1
    END IF
  END DO
END DO
DO i=1,self%nc
  self%ho_info%lcp(1,i)=i+self%ne
  self%ho_info%r(:,i+self%ne)=(self%r(:,self%lc(1,i))+self%r(:,self%lc(2,i)) &
    + self%r(:,self%lc(3,i))+self%r(:,self%lc(4,i)))/4.d0
  self%ho_info%r(:,i+self%ne) = r_ho(:,lc_ho(9,i))
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,lc_ho,emark)
DEBUG_STACK_POP
end subroutine mesh_cubit_hobase_quad
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_hobase_hex
!------------------------------------------------------------------------------
!> Add quadratic mesh node points from HEX27 import
!------------------------------------------------------------------------------
subroutine mesh_cubit_hobase_hex(self)
CLASS(oft_mesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2),ftmp(4)
integer(i4), allocatable :: emark(:),fmark(:)
DEBUG_STACK_PUSH
!---Setup quadratic mesh
self%order=1
self%ho_info%nep=1
self%ho_info%nfp=1
self%ho_info%ncp=1
ALLOCATE(self%ho_info%r(3,self%ne+self%nf+self%nc))
ALLOCATE(self%ho_info%lep(self%ho_info%nep,self%ne))
ALLOCATE(self%ho_info%lfp(self%ho_info%nfp,self%nf))
ALLOCATE(self%ho_info%lcp(self%ho_info%ncp,self%nc))
!---Initialize high order points with straight edges
DO i=1,self%ne
  self%ho_info%lep(1,i)=i
  self%ho_info%r(:,i)=(self%r(:,self%le(1,i))+self%r(:,self%le(2,i)))/2.d0
END DO
!---Set edge midpoints from imported list
ALLOCATE(emark(self%ne))
emark=0
DO i=1,self%nc
  DO j=1,self%cell_ne
    etmp=lc_ho(exodus_hex_emap(:,j),i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','mesh_cubit_hobase_hex',__FILE__)
    IF(emark(k)==0)THEN
      self%ho_info%r(:,k) = r_ho(:,lc_ho(8+j,i))
      emark(k)=1
    END IF
  END DO
END DO
DEALLOCATE(emark)
!---Initialize high order points with flat faces
DO i=1,self%nf
  self%ho_info%lfp(1,i)=i+self%ne
  self%ho_info%r(:,i+self%ne)=(self%r(:,self%lf(1,i)) + self%r(:,self%lf(2,i)) &
    + self%r(:,self%lf(3,i)) + self%r(:,self%lf(4,i)))/4.d0
END DO
!---Set face centers from imported list
ALLOCATE(fmark(self%nf))
fmark=0
DO i=1,self%nc
  DO j=1,self%cell_nf
    ftmp=lc_ho(exodus_hex_fmap(:,j),i)
    k=ABS(mesh_local_findface(self,ftmp))
    if(k==0)CALL oft_abort('Unlinked mesh face','mesh_cubit_hobase_hex',__FILE__)
    IF(fmark(k)==0)THEN
      self%ho_info%r(:,k+self%ne) = r_ho(:,lc_ho(20+j,i))
      fmark(k)=1
    END IF
  END DO
END DO
DEALLOCATE(fmark)
!---Set cell points from imported list
DO i=1,self%nc
  self%ho_info%lcp(1,i)=i+self%ne+self%nf
  self%ho_info%r(:,i+self%ne+self%nf)= (self%r(:,self%lc(1,i)) + self%r(:,self%lc(2,i)) &
    + self%r(:,self%lc(3,i)) + self%r(:,self%lc(4,i)) + self%r(:,self%lc(5,i)) &
    + self%r(:,self%lc(6,i)) + self%r(:,self%lc(7,i)) + self%r(:,self%lc(8,i)))/8.d0
  self%ho_info%r(:,i+self%ne+self%nf) = r_ho(:,lc_ho(27,i))
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,lc_ho)
DEBUG_STACK_POP
end subroutine mesh_cubit_hobase_hex
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_reffix
!------------------------------------------------------------------------------
!> Adjust points to CAD boundary and propogate CAD linkage
!------------------------------------------------------------------------------
subroutine mesh_cubit_reffix
#ifdef HAVE_ONURBS
real(r8) :: u1,v1,u2,v2,rp(3),p(2),pt(3),rad,p1(2),p2(2),f(4)
integer(i4) :: i,j,k,l,ed,lecors(2),npcors,lfecors(3),ho_count
integer(i4) :: ierr,nerr,nbecors,ep(2),fp(3),ind,sid,cell
integer(i4), allocatable, dimension(:) :: lce_tmp,lcf_tmp
type(Exodus_cadlink), pointer :: pmesh_cad_link
CHARACTER(LEN=60) :: error_str
#else
real(r8) :: rp(3)
integer(i4) :: i,j,ho_count
#endif
class(oft_mesh), pointer :: pmesh
DEBUG_STACK_PUSH
!---Get parent mesh
pmesh=>mg_mesh%meshes(mg_mesh%level-1)
!---If only one level do nothing
IF(mg_mesh%level==1)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(pmesh%fullmesh.AND.(.NOT.mesh%fullmesh))THEN ! Current level is transfer level
#ifdef HAVE_ONURBS
  !---Synchronize CAD linkage to distributed mesh
  IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Copying geometry linkage to distributed mesh'
  CALL oft_increase_indent
  !---Get CAD representation aliases
  pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
  cad_link=>ML_cad_link(mg_mesh%level)
  !---Index to edges on new mesh
  IF(pmesh_cad_link%nce>0)THEN
    ALLOCATE(lce_tmp(pmesh%ne))
    lce_tmp=0
    DO i=1,mesh%ne
      lce_tmp(ABS(mesh%base%le(i)))=i
    END DO
    cad_link%nce=0
    DO i=1,pmesh_cad_link%nce
      ind=pmesh_cad_link%lce(i)
      IF(lce_tmp(ind)==0)CYCLE
      cad_link%nce = cad_link%nce + 1
    END DO
    ALLOCATE(cad_link%lce(cad_link%nce),cad_link%lbeg(cad_link%nce))
    cad_link%nce=0
    DO i=1,pmesh_cad_link%nce
      ind=pmesh_cad_link%lce(i)
      IF(lce_tmp(ind)==0)CYCLE
      cad_link%nce = cad_link%nce + 1
      cad_link%lce(cad_link%nce) = lce_tmp(ind)
      cad_link%lbeg(cad_link%nce)%wo => pmesh_cad_link%lbeg(i)%wo
    END DO
    DEALLOCATE(lce_tmp)
  END IF
  !---Index to faces on new mesh
  IF(pmesh_cad_link%ncf>0)THEN
    ALLOCATE(lce_tmp(pmesh%nf))
    lce_tmp=0
    DO i=1,mesh%nf
      lce_tmp(ABS(mesh%base%lf(i)))=i
    END DO
    cad_link%ncf=0
    DO i=1,pmesh_cad_link%ncf
      ind=pmesh_cad_link%lcf(i)
      if(lce_tmp(ind)==0)CYCLE
      cad_link%ncf = cad_link%ncf + 1
    END DO
    ALLOCATE(cad_link%lcf(cad_link%ncf),cad_link%lbfg(cad_link%ncf))
    cad_link%ncf=0
    DO i=1,pmesh_cad_link%ncf
      ind=pmesh_cad_link%lcf(i)
      if(lce_tmp(ind)==0)CYCLE
      cad_link%ncf = cad_link%ncf + 1
      cad_link%lcf(cad_link%ncf) = lce_tmp(ind)
      cad_link%lbfg(cad_link%ncf)%wo => pmesh_cad_link%lbfg(i)%wo
    END DO
    DEALLOCATE(lce_tmp)
  END IF
#endif
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Handle toroidal grid
IF(tor_mesh)THEN
  DO i=1,pmesh%ne
    rp(1)=SQRT(SUM(pmesh%r(1:2,pmesh%le(1,i))**2))
    rp(2)=SQRT(SUM(pmesh%r(1:2,pmesh%le(2,i))**2))
    rp(3)=SQRT(SUM(mesh%r(1:2,i+pmesh%np)**2))
    mesh%r(1:2,i+pmesh%np)=mesh%r(1:2,i+pmesh%np)*(rp(1)+rp(2))/(2.d0*rp(3))
  END DO
  DEBUG_STACK_POP
  RETURN
END IF
!---Refine new boundary points using CAD model
IF(TRIM(inpname)/='none')THEN
#ifdef HAVE_ONURBS
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Adjusting points to OpenNURBS CAD boundary'
CALL oft_increase_indent
!---Get CAD representation aliases
pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
cad_link=>ML_cad_link(mg_mesh%level)
!---Locate edge end points and place daughter point
nerr=0
!!$omp parallel do private(j,u1,v1,u2,v2,pt) reduction(+:nerr)
DO i=1,pmesh_cad_link%nce
  j=pmesh_cad_link%lce(i)+pmesh%np
  IF(.NOT.ASSOCIATED(pmesh_cad_link%lbeg(i)%wo))CYCLE
  SELECT TYPE(obj=>pmesh_cad_link%lbeg(i)%wo)
    TYPE IS(nurbs_curve)
      IF(obj%linear)CYCLE
      pt=mesh%r(:,j)
      CALL nurbs_curve_midpoint(obj,pt,pmesh%r(:,pmesh%le(1,pmesh_cad_link%lce(i))), &
        pmesh%r(:,pmesh%le(2,pmesh_cad_link%lce(i))),.5d0,.5d0,ierr)
      IF(ierr==0)THEN
        mesh%r(:,j)=pt
      ELSE
        nerr=nerr+1
      END IF
    TYPE IS(nurbs_surf)
      IF(obj%planar)CYCLE
      pt=mesh%r(:,j)
      CALL nurbs_surf_midpoint(obj,pt,pmesh%r(:,pmesh%le(1,pmesh_cad_link%lce(i))), &
        pmesh%r(:,pmesh%le(2,pmesh_cad_link%lce(i))),.5d0,.5d0,ierr)
      IF(ierr==0)THEN
        mesh%r(:,j)=pt
      ELSE
        nerr=nerr+1
      END IF
  END SELECT
END DO
!---Synchornize error count
nerr=oft_mpi_sum(nerr)
IF(oft_env%head_proc.AND.nerr>0)THEN
  WRITE(error_str,'(A,I4,A)')'Refinement node placement failed for ',nerr,' edges'
  CALL oft_warn(error_str)
END IF
!---
cad_link%nce = pmesh_cad_link%nce*2 + pmesh_cad_link%ncf*3
ALLOCATE(cad_link%lce(cad_link%nce),cad_link%lbeg(cad_link%nce))
npcors=pmesh%np
!---Copy edge linkage to daughter edges
!$omp parallel do private(i,k)
DO j=1,pmesh_cad_link%nce
  i=pmesh_cad_link%lce(j)
  IF(ASSOCIATED(pmesh_cad_link%lbeg(j)%wo))THEN
    DO k=1,2
      cad_link%lce((j-1)*2+k) = mg_mesh%inter(mg_mesh%level-1)%lede(k,i)
      cad_link%lbeg((j-1)*2+k)%wo => pmesh_cad_link%lbeg(j)%wo
    END DO
  END IF
END DO
cad_link%ncf = pmesh_cad_link%ncf*4
ALLOCATE(cad_link%lcf(cad_link%ncf),cad_link%lbfg(cad_link%ncf))
!---Copy face linkage to daughter edges and faces
!$omp parallel do private(i,k)
DO j=1,pmesh_cad_link%ncf
  i=pmesh_cad_link%lcf(j)
  IF(associated(pmesh_cad_link%lbfg(j)%wo))THEN
    DO k=1,3
      cad_link%lce(pmesh_cad_link%nce*2 + (j-1)*3 + k) = mg_mesh%inter(mg_mesh%level-1)%lfde(k,i)
      cad_link%lbeg(pmesh_cad_link%nce*2 + (j-1)*3 + k)%wo => pmesh_cad_link%lbfg(j)%wo
    END DO
    DO k=1,4
      cad_link%lcf((j-1)*4 + k) = mg_mesh%inter(mg_mesh%level-1)%lfdf(k,i)
      cad_link%lbfg((j-1)*4 + k)%wo => pmesh_cad_link%lbfg(j)%wo
    END DO
  END IF
END DO
#else
  IF(.NOT.oft_env%test_run)CALL oft_abort('OFT not compiled with OpenNURBS', &
    'mesh_cubit_reffix',__FILE__)
#endif
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_cubit_reffix
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_add_quad
!------------------------------------------------------------------------------
!> Add quadratic mesh node points using CAD model
!------------------------------------------------------------------------------
subroutine mesh_cubit_add_quad
real(r8) :: pt(3),r1,r2,r3,pts_tmp(3,10),wts_tmp(10)
integer(i4) :: i,j,k,ierr,nerr
CHARACTER(LEN=60) :: error_str
DEBUG_STACK_PUSH
IF(tor_mesh)THEN
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Setting toroidal quadratic nodes'
  wts_tmp(1:2)=1.d0/2.d0
  DO i=1,mesh%ne
    pts_tmp(:,1)=mesh%r(:,mesh%le(1,i))
    pts_tmp(:,2)=mesh%r(:,mesh%le(2,i))
    CALL cubit_ho_tor(pts_tmp,wts_tmp,2,tor_rmin,mesh%ho_info%r(:,i))
  END DO
  IF(mesh%ho_info%nfp>0)THEN
    wts_tmp(1:4)=1.d0/4.d0
    do i=1,mesh%nf
      DO j=1,mesh%face_np
        pts_tmp(:,j)=mesh%r(:,mesh%lf(j,i))
      END DO
      CALL cubit_ho_tor(pts_tmp,wts_tmp,4,tor_rmin,mesh%ho_info%r(:,i+mesh%ne))
    end do
  END IF
  IF(mesh%ho_info%ncp>0)THEN
    wts_tmp(1:8)=1.d0/8.d0
    do i=1,mesh%nc
      DO j=1,mesh%cell_np
        pts_tmp(:,j)=mesh%r(:,mesh%lc(j,i))
      END DO
      CALL cubit_ho_tor(pts_tmp,wts_tmp,8,tor_rmin,mesh%ho_info%r(:,i+mesh%ne+mesh%nf))
    end do
  END IF
  DEBUG_STACK_POP
  RETURN
END IF
!
IF(TRIM(inpname)/='none')THEN
#ifdef HAVE_ONURBS
  if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Setting Cubit quadratic nodes'
  !---Get CAD representation alias
  cad_link=>ML_cad_link(mg_mesh%level)
  !---Locate edge end points and place daughter node
  nerr=0
  DO i=1,cad_link%nce
    IF(.NOT.associated(cad_link%lbeg(i)%wo))CYCLE
    j=cad_link%lce(i)
    SELECT TYPE(obj=>cad_link%lbeg(i)%wo)
      TYPE IS(nurbs_curve)
        IF(obj%linear)CYCLE
        pt=mesh%ho_info%r(:,j)
        CALL nurbs_curve_midpoint(obj,pt,mesh%r(:,mesh%le(1,j)),mesh%r(:,mesh%le(2,j)), &
          0.5d0,0.5d0,ierr)
        IF(ierr==0)THEN
          mesh%ho_info%r(:,j)=pt
          DO k=mesh%kec(j),mesh%kec(j+1)-1
            mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
          END DO
        ELSE
          nerr=nerr+1
        END IF
      TYPE IS(nurbs_surf)
        IF(obj%planar)CYCLE
        pt=mesh%ho_info%r(:,j)
        CALL nurbs_surf_midpoint(obj,pt,mesh%r(:,mesh%le(1,j)),mesh%r(:,mesh%le(2,j)), &
          0.5d0,0.5d0,ierr)
        IF(ierr==0)THEN
          mesh%ho_info%r(:,j)=pt
          DO k=mesh%kec(j),mesh%kec(j+1)-1
            mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
          END DO
        ELSE
          nerr=nerr+1
        END IF
    END SELECT
  END DO
  !---Synchornize error count
  nerr=oft_mpi_sum(nerr)
  IF(oft_env%head_proc.AND.nerr>0)THEN
    WRITE(error_str,'(A,I4,A)')'Quadratic node placement failed for ',nerr,' edges'
    CALL oft_warn(error_str)
  END IF
#else
  IF(.NOT.oft_env%test_run)CALL oft_abort('OFT not compiled with OpenNURBS', &
    'mesh_cubit_add_quad',__FILE__)
#endif
END IF
DEBUG_STACK_POP
end subroutine mesh_cubit_add_quad
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_test_edge
!------------------------------------------------------------------------------
!> Add quadratic mesh node points using CAD model
!------------------------------------------------------------------------------
subroutine mesh_cubit_test_edge(i)
integer(i4), intent(in) :: i
#ifdef HAVE_ONURBS
real(r8) :: rp(3),pt(3)
integer(i4) :: j,ierr
integer(i4), allocatable :: emap(:)
DEBUG_STACK_PUSH
if(oft_env%head_proc)write(*,'(2A)')oft_indent,'Testing Cubit edge location'
call oft_increase_indent
!---Get CAD representation alias
cad_link=>ML_cad_link(mg_mesh%level)
!---Locate edge end points and place daughter node
rp=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
DO j=1,cad_link%nce
  IF(cad_link%lce(j)==i)EXIT
END DO
!---
IF(j<=cad_link%nce)THEN
  IF(.NOT.ASSOCIATED(cad_link%lbeg(j)%wo))THEN
    WRITE(*,'(2A)')oft_indent,'No CAD object linked'
    call oft_decrease_indent
    DEBUG_STACK_POP
    RETURN
  END IF
  !---
  SELECT TYPE(obj=>cad_link%lbeg(j)%wo)
    TYPE IS(nurbs_curve)
      WRITE(*,'(2A)')oft_indent,'Entity is NURBS curve'
      call oft_increase_indent
      WRITE(*,'(2A,L)')oft_indent,'Linear:   ',obj%linear
      WRITE(*,'(2A,L)')oft_indent,'Periodic: ',obj%periodic
      WRITE(*,'(2A,2ES14.6)')oft_indent,'Domain:   ',obj%domain
      WRITE(*,'(2A,6ES14.6)')oft_indent,'Extent:   ',obj%rgrid(:,1),obj%rgrid(:,20)
      call oft_decrease_indent
      !---
      pt=rp
      call nurbs_curve_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)),mesh%r(:,mesh%le(2,i)), &
      1.d0/3,2.d0/3,ierr)
      !---
      pt=rp
      call nurbs_curve_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)),mesh%r(:,mesh%le(2,i)), &
      2.d0/3,1.d0/3,ierr)
    TYPE IS(nurbs_surf)
      WRITE(*,'(2A)')oft_indent,'Entity is NURBS surface',obj%sid
      call oft_increase_indent
      WRITE(*,'(2A,L)')oft_indent,'Planar:   ',obj%planar
      WRITE(*,'(2A,L,1X,L)')oft_indent,'Periodic: ',obj%periodic
      WRITE(*,'(2A,L,1X,L,1X,L,1X,L)')oft_indent,'Singular: ',obj%singular(:,1),obj%singular(:,2)
      WRITE(*,'(2A,4ES14.6)')oft_indent,'Domain:   ',obj%domain(:,1),obj%domain(:,2)
      WRITE(*,'(2A,6ES14.6)')oft_indent,'Edge 1:   ',obj%rgrid(:,1,1),obj%rgrid(:,2,1)
      WRITE(*,'(2A,6ES14.6)')oft_indent,'Edge 2:   ',obj%rgrid(:,20,1),obj%rgrid(:,20,2)
      WRITE(*,'(2A,6ES14.6)')oft_indent,'Edge 3:   ',obj%rgrid(:,20,20),obj%rgrid(:,19,20)
      WRITE(*,'(2A,6ES14.6)')oft_indent,'Edge 4:   ',obj%rgrid(:,1,20),obj%rgrid(:,1,19)
      call oft_decrease_indent
      !---
      pt=rp
      call nurbs_surf_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)),mesh%r(:,mesh%le(2,i)), &
      1.d0/3,2.d0/3,ierr)
      !---
      pt=rp
      call nurbs_surf_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)),mesh%r(:,mesh%le(2,i)), &
      2.d0/3,1.d0/3,ierr)
    CLASS DEFAULT
      WRITE(*,'(2A)')oft_indent,'Could not determine NURBS entity'
  END SELECT
ELSE
  WRITE(*,'(2A)')oft_indent,'Edge is not on a CAD surface'
END IF
call oft_decrease_indent
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with OpenNURBS.','mesh_cubit_test_edge',__FILE__)
#endif
end subroutine mesh_cubit_test_edge
!------------------------------------------------------------------------------
! SUBROUTINE mesh_cubit_error
!------------------------------------------------------------------------------
!> Catch NETCDF errors
!------------------------------------------------------------------------------
subroutine mesh_cubit_error(status)
integer(i4), intent(in) :: status
if(status /= nf90_noerr) THEN
  WRITE(*,'(A,I8.8)')'NETCDF-ERROR_CODE: ',status
  call oft_abort(TRIM(nf90_strerror(status)),'mesh_cubit_error',__FILE__)
end if
end subroutine mesh_cubit_error
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_reflect
!------------------------------------------------------------------------------
!> Reflect an exodus mesh and CAD model across the xy-plane.
!!
!! @param[in] tol Tolerance for marking point as on the reflection plane
!------------------------------------------------------------------------------
subroutine mesh_cubit_reflect(tol)
real(r8), intent(in) :: tol
integer(i4) :: npold,ncold,i,j,ic,is,cid_max,sid_max,nreg,n_ho
integer(i4), allocatable :: newindex(:),hoindex(:),regtmp(:),ltemp(:,:)
real(r8), allocatable :: rtemp(:,:)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2A)')oft_indent,'Reflecting mesh -> z'
CALL oft_increase_indent
!---Reflect points that are not on reflection plane
npold=mesh%np
allocate(newindex(2*mesh%np),rtemp(3,2*mesh%np))
rtemp(:,1:mesh%np)=mesh%r
deallocate(mesh%r)
do i=1,npold
  IF(ABS(rtemp(3,i))<=tol)THEN
    rtemp(3,i)=0.d0
    newindex(i)=i
  ELSE
    mesh%np=mesh%np+1
    rtemp(1:2,mesh%np)= rtemp(1:2,i)
    rtemp(3,mesh%np)  =-rtemp(3,i)
    newindex(i)=mesh%np
  ENDIF
enddo
allocate(mesh%r(3,mesh%np))
mesh%r=rtemp(:,1:mesh%np)
deallocate(rtemp)
#ifdef HAVE_ONURBS
!---Reflect CAD model edges that are not on reflection plane
if(ngmc>0)then
  cid_max=0
  DO ic=1,ngmc
    cid_max=MAX(cid_max,model_curves(ic)%top_cid)
  END DO
  DO ic=1,ngmc
    allocate(model_curves(ic+ngmc)%lgme(2,model_curves(ic)%ne))
    model_curves(ic+ngmc)%ne=model_curves(ic)%ne
    model_curves(ic+ngmc)%cid=model_curves(ic)%cid
    model_curves(ic+ngmc)%top_cid=model_curves(ic)%top_cid + cid_max
    model_curves(ic+ngmc)%name=model_curves(ic)%name
    do i=1,model_curves(ic)%ne
      model_curves(ic+ngmc)%lgme(:,i)=newindex(model_curves(ic)%lgme(:,i))
    enddo
  END DO
  ngmc=2*ngmc
endif
!---Reflect CAD model faces that are not on reflection plane
if(ngms>0)then
  sid_max=0
  DO is=1,ngms
    sid_max=MAX(sid_max,model_surfaces(is)%top_sid)
  END DO
  DO is=1,ngms
    allocate(model_surfaces(is+ngms)%lgmf(mesh%face_np,model_surfaces(is)%nf))
    model_surfaces(is+ngms)%nf=model_surfaces(is)%nf
    model_surfaces(is+ngms)%sid=model_surfaces(is)%sid
    model_surfaces(is+ngms)%top_sid=model_surfaces(is)%top_sid + sid_max
    model_surfaces(is+ngms)%name=model_surfaces(is)%name
    do i=1,model_surfaces(is)%nf
      model_surfaces(is+ngms)%lgmf(:,i)=newindex(model_surfaces(is)%lgmf(:,i))
    enddo
  END DO
  ngms=2*ngms
endif
!---START CAD reflection
!---Reflect CAD curves
if(ngwc>0)then
  DO ic=1,ngwc
    wf_curves(ic+ngwc)%id=wf_curves(ic)%id
    wf_curves(ic+ngwc)%cid=wf_curves(ic)%cid + cid_max
    wf_curves(ic+ngwc)%domain=wf_curves(ic)%domain
    wf_curves(ic+ngwc)%periodic=wf_curves(ic)%periodic
    wf_curves(ic+ngwc)%linear=wf_curves(ic)%linear
    wf_curves(ic+ngwc)%reflect=3
    CALL wf_curves(ic+ngwc)%grid()
  END DO
  ngwc=2*ngwc
endif
!---Reflect CAD surfaces
if(ngws>0)then
  DO is=1,ngws
    wf_surfs(is+ngws)%id=wf_surfs(is)%id
    wf_surfs(is+ngws)%sid=wf_surfs(is)%sid + sid_max
    wf_surfs(is+ngws)%domain=wf_surfs(is)%domain
    wf_surfs(is+ngws)%periodic=wf_surfs(is)%periodic
    wf_surfs(is+ngws)%singular=wf_surfs(is)%singular
    wf_surfs(is+ngws)%planar=wf_surfs(is)%planar
    wf_surfs(is+ngws)%reflect=3
    CALL wf_surfs(is+ngws)%grid()
  END DO
  ngws=2*ngws
endif
!---END CAD Reflection
#endif
!---Reflect cells
ncold=mesh%nc
allocate(ltemp(mesh%cell_np,2*mesh%nc),regtmp(2*mesh%nc))
ltemp(:,1:mesh%nc)=mesh%lc
regtmp(1:mesh%nc)=mesh%reg
nreg = MAXVAL(mesh%reg)
deallocate(mesh%lc,mesh%reg)
do i=1,ncold
  mesh%nc=mesh%nc+1
  DO j=1,mesh%cell_np
    ltemp(j,mesh%nc)=newindex(ltemp(j,i))
  END DO
  regtmp(mesh%nc)=regtmp(i)+nreg
enddo
allocate(mesh%lc(mesh%cell_np,mesh%nc),mesh%reg(mesh%nc))
mesh%lc=ltemp(:,1:mesh%nc)
mesh%reg=regtmp
deallocate(ltemp,regtmp)
!---Reflect high-order nodes
IF(have_ho)THEN
  n_ho=6
  IF(mesh%type==3)n_ho=19
  npold=np_ho
  allocate(rtemp(3,2*np_ho),hoindex(2*np_ho))
  rtemp(:,1:np_ho)=r_ho
  deallocate(r_ho)
  DO i=1,npold
    IF(ABS(rtemp(3,i))<=tol)THEN
      rtemp(3,i)=0.d0
      hoindex(i)=i
    ELSE
      np_ho=np_ho+1
      rtemp(1:2,np_ho)= rtemp(1:2,i)
      rtemp(3,np_ho)  =-rtemp(3,i)
      hoindex(i)=np_ho
    ENDIF
  END DO
  allocate(r_ho(3,np_ho))
  r_ho=rtemp(:,1:np_ho)
  deallocate(rtemp)
  !---
  mesh%nc=ncold
  allocate(ltemp(mesh%cell_np+n_ho,2*mesh%nc))
  ltemp(:,1:mesh%nc)=lc_ho
  deallocate(lc_ho)
  DO i=1,ncold
    mesh%nc=mesh%nc+1
    ltemp(1:mesh%cell_np,mesh%nc)  = newindex(ltemp(1:mesh%cell_np,i))
    ltemp(mesh%cell_np+1:mesh%cell_np+n_ho,mesh%nc) = &
      hoindex(ltemp(mesh%cell_np+1:mesh%cell_np+n_ho,i))
  END DO
  allocate(lc_ho(mesh%cell_np+n_ho,mesh%nc))
  lc_ho=ltemp(:,1:mesh%nc)
  deallocate(ltemp,hoindex)
END IF
!---Flag periodic points
IF(np_per>0)THEN
  ALLOCATE(mesh%periodic%lp(mesh%np))
  mesh%periodic%lp=-1
  DO i=1,np_per
    j=per_nodes(i)
    mesh%periodic%lp(newindex(j))=j
  END DO
END IF
deallocate(newindex)
! IF(oft_debug_print(3))write(*,'(A)')oft_indent,'Done'
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_cubit_reflect
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_cubit_set_periodic
!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine mesh_cubit_set_periodic
integer(i4) :: i,j,pt_e(2),ind,k,kk
integer(i4), ALLOCATABLE :: pt_f(:)
IF(np_per==0)RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting Cubit periodicity'
CALL oft_increase_indent
IF(.NOT.ASSOCIATED(mesh%periodic%lp))THEN
  ALLOCATE(mesh%periodic%lp(mesh%np))
  mesh%periodic%lp=-1
  DO i=1,np_per
    j=per_nodes(i)
    DO k=1,mesh%nbp
      kk=mesh%lbp(k)
      IF(kk==j)CYCLE
      IF(ALL(ABS(mesh%r(1:2,j)-mesh%r(1:2,kk))<1.d-8).AND.ABS(mesh%r(3,kk))<1.d-8)THEN
        IF(kk>j)THEN
          mesh%periodic%lp(kk)=j
        ELSE
          mesh%periodic%lp(j)=kk
        END IF
        EXIT
      END IF
    END DO
  END DO
END IF
!---Set periodic faces
mesh%periodic%nper=1
ALLOCATE(mesh%periodic%le(mesh%ne))
ALLOCATE(mesh%periodic%lf(mesh%nf))
mesh%periodic%le=-1
mesh%periodic%lf=-1
!---Flag periodic edges
!$omp parallel private(j,pt_e,pt_f,ind)
allocate(pt_f(mesh%face_np))
!$omp do
DO i=1,mesh%nbe
  j = mesh%lbe(i)
  pt_e=mesh%periodic%lp(mesh%le(:,j))
  IF(ALL(pt_e>0))THEN
    ind=ABS(mesh_local_findedge(mesh,pt_e))
    IF(ind==0)WRITE(*,'(2A,2I8)')oft_indent,'Bad edge',i,ind
    mesh%periodic%le(j)=ind
  END IF
END DO
!---Flag periodic faces
!$omp do
DO i=1,mesh%nbf
  j = mesh%lbf(i)
  pt_f=mesh%periodic%lp(mesh%lf(:,j))
  IF(ALL(pt_f>0))THEN
    ind=ABS(mesh_local_findface(mesh,pt_f))
    IF(ind==0)WRITE(*,'(2A,2I8)')oft_indent,'Bad face',i,ind
    mesh%periodic%lf(j)=ind
  END IF
END DO
deallocate(pt_f)
!$omp end parallel
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine mesh_cubit_set_periodic
#else
USE oft_base
INTEGER(i4), PARAMETER, PUBLIC :: mesh_cubit_id = 2
CONTAINS
#endif
!
subroutine cubit_ho_tor(pts,wts,n,rmin,pt)
real(r8), intent(in) :: pts(3,n),wts(n),rmin
real(r8), intent(inout) :: pt(3)
integer(i4), intent(in) :: n
real(r8) :: rpos,rnew
integer(i4) :: i
pt=0.d0
rpos=0.d0
DO i=1,n
  rpos=rpos+SQRT(SUM(pts(1:2,i)**2))*wts(i)
  pt=pt+pts(:,i)*wts(i)
END DO
!
IF(rpos>rmin)THEN
  rnew=SQRT(SUM(pt(1:2)**2))
  pt(1:2)=pt(1:2)*rpos/rnew
END IF
end subroutine cubit_ho_tor
end module oft_mesh_cubit
