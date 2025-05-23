!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_mesh_native.F90
!
!> Mesh handling for OFT native format, generated by python pre-processing
!! scripts.
!!
!! @author Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
module oft_mesh_native
USE oft_base
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
USE oft_mesh_type, ONLY: oft_amesh, oft_mesh, oft_bmesh
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_hexmesh_type, ONLY: oft_hexmesh
USE oft_quadmesh_type, ONLY: oft_quadmesh
USE oft_mesh_local_util, ONLY: mesh_local_findedge, mesh_local_findface
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE multigrid, ONLY: multigrid_mesh, multigrid_level
IMPLICIT NONE
#include "local.h"
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'none' !< Name of input file for mesh
INTEGER(i4), PARAMETER, PUBLIC :: mesh_native_id = 0
REAL(r8), ALLOCATABLE, PUBLIC :: r_mem(:,:)
INTEGER(i4), ALLOCATABLE, PUBLIC :: lc_mem(:,:)
INTEGER(i4), ALLOCATABLE, PUBLIC :: reg_mem(:)
!---
INTEGER(i4) :: np_ho = 0
REAL(r8), ALLOCATABLE :: r_ho(:,:)
INTEGER(i4), ALLOCATABLE :: le_ho(:,:),lf_ho(:,:)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: per_nodes
CONTAINS
!---------------------------------------------------------------------------------
!> Finalize setup/load-in of native mesh and destroy temporaries created
!! for grid construction (eg. high-order input nodes, in-memory data)
!---------------------------------------------------------------------------------
subroutine native_finalize_setup
!
np_ho=0
IF(ALLOCATED(r_ho))DEALLOCATE(r_ho)
IF(ALLOCATED(le_ho))DEALLOCATE(le_ho)
IF(ALLOCATED(lf_ho))DEALLOCATE(lf_ho)
IF(ALLOCATED(per_nodes))DEALLOCATE(per_nodes)
!
IF(ALLOCATED(r_mem))DEALLOCATE(r_mem)
IF(ALLOCATED(lc_mem))DEALLOCATE(lc_mem)
IF(ALLOCATED(reg_mem))DEALLOCATE(reg_mem)
end subroutine native_finalize_setup
!---------------------------------------------------------------------------------
!> Read in t3d mesh file from file "filename"
!! - Read in T3D options from input file
!! - Read in mesh points and cells
!! - Read in surface IDs for CAD edges and faces
!---------------------------------------------------------------------------------
subroutine native_load_vmesh(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
logical :: success
integer(i4) :: i,id,ierr,io_unit,ndims,np_mem,mesh_order
integer(i4), allocatable, dimension(:) :: dim_sizes
LOGICAL :: reflect = .FALSE.
LOGICAL :: ref_periodic = .FALSE.
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
!---Read in mesh options
namelist/native_mesh_options/filename,reflect,ref_periodic!,zstretch
DEBUG_STACK_PUSH
np_mem=-1
filename='none'
IF(oft_env%head_proc)THEN
    IF(ALLOCATED(r_mem).AND.ALLOCATED(lc_mem).AND.ALLOCATED(reg_mem))THEN
        np_mem=SIZE(r_mem,DIM=2)
    ELSE
        OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
        READ(io_unit,native_mesh_options,IOSTAT=ierr)
        CLOSE(io_unit)
        IF(ierr<0)CALL oft_abort('No "native_mesh_options" found in input file.','native_load_vmesh',__FILE__)
        IF(ierr>0)CALL oft_abort('Error parsing "native_mesh_options" in input file.','native_load_vmesh',__FILE__)
        IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','native_load_vmesh',__FILE__)
        !
        WRITE(*,'(2A)')oft_indent,'Native volume mesh:'
        CALL oft_increase_indent
        WRITE(*,'(3A)')oft_indent,'Filename = ',TRIM(filename)
    END IF
ELSE
    CALL oft_increase_indent
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(np_mem,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
CALL MPI_Bcast(filename,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
CALL MPI_Bcast(reflect,1,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
CALL MPI_Bcast(ref_periodic,1,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
! CALL MPI_Bcast(zstretch, 1,OFT_MPI_R8,0,oft_env%COMM,ierr)
! IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
#endif
!---Read in mesh sizes
IF(np_mem>0)THEN
    ALLOCATE(dim_sizes(2))
    dim_sizes=SHAPE(lc_mem)
ELSE
    CALL hdf5_field_get_sizes(TRIM(filename),"mesh/LC",ndims,dim_sizes)
    IF(ndims==-1)CALL oft_abort('"mesh/LC" field does not exist in input file', 'native_load_vmesh', __FILE__)
END IF
IF(dim_sizes(2)==0)CALL oft_abort('Zero length cell list','native_load_vmesh',__FILE__)
!---Allocate mesh
SELECT CASE(dim_sizes(1))
    CASE(4)
        ALLOCATE(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
        ALLOCATE(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
    CASE(8)
        ALLOCATE(oft_hexmesh::mg_mesh%meshes(mg_mesh%mgdim))
        ALLOCATE(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
    CASE DEFAULT
        CALL oft_abort('Invalid cell type','native_load_vmesh',__FILE__)
END SELECT
DO i=1,mg_mesh%mgdim
    CALL mg_mesh%meshes(i)%setup(mesh_native_id)
    CALL mg_mesh%smeshes(i)%setup(mesh_native_id,.TRUE.)
    mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
END DO
CALL multigrid_level(mg_mesh,1)
mesh=>mg_mesh%meshes(1)
smesh=>mg_mesh%smeshes(1)
mesh%nc=dim_sizes(2)
!
IF(np_mem>0)THEN
    ALLOCATE(dim_sizes(2))
    dim_sizes=SHAPE(r_mem)
ELSE
    CALL hdf5_field_get_sizes(TRIM(filename),"mesh/R",ndims,dim_sizes)
    IF(ndims==-1)CALL oft_abort('"mesh/R" field does not exist in input file', 'native_load_vmesh', __FILE__)
END IF
mesh%np=dim_sizes(2)
IF(mesh%np==0)CALL oft_abort('Zero length point list','native_load_vmesh',__FILE__)
!---Read in points
ALLOCATE(mesh%r(3,mesh%np))
IF(np_mem>0)THEN
    mesh%r=r_mem
    DEALLOCATE(r_mem)
ELSE
    CALL hdf5_read(mesh%r,TRIM(filename),"mesh/R",success)
    IF(.NOT.success)CALL oft_abort('Error reading points','native_load_vmesh',__FILE__)
END IF
!---Read in cells
ALLOCATE(mesh%lc(mesh%cell_np,mesh%nc))
IF(np_mem>0)THEN
    mesh%lc=lc_mem
    DEALLOCATE(lc_mem)
ELSE
    CALL hdf5_read(mesh%lc,TRIM(filename),"mesh/LC",success)
    IF(.NOT.success)CALL oft_abort('Error reading cells','native_load_vmesh',__FILE__)
END IF
!---Read in region list and other information
mesh_order=1
ALLOCATE(mesh%reg(mesh%nc))
IF(np_mem>0)THEN
    mesh%reg=reg_mem
    DEALLOCATE(reg_mem)
ELSE
    CALL hdf5_read(mesh%reg,TRIM(filename),"mesh/REG",success)
    IF(.NOT.success)CALL oft_abort('Error reading region list','native_load_vmesh',__FILE__)
    !---Read high-order mesh information
    IF(hdf5_field_exist(TRIM(filename),"mesh/ho_info/LE"))THEN
        mesh_order=2
        CALL hdf5_field_get_sizes(TRIM(filename),"mesh/ho_info/R",ndims,dim_sizes)
        IF(ndims==-1)CALL oft_abort('"mesh/ho_info/R" field does not exist in input file', 'native_load_vmesh', __FILE__)
        np_ho=dim_sizes(2)
        ALLOCATE(r_ho(3,np_ho))
        CALL hdf5_read(r_ho,TRIM(filename),"mesh/ho_info/R",success)
        IF(.NOT.success)CALL oft_abort('Error reading quadratic points','native_load_vmesh',__FILE__)
        CALL hdf5_field_get_sizes(TRIM(filename),"mesh/ho_info/LE",ndims,dim_sizes)
        ALLOCATE(le_ho(2,dim_sizes(2)))
        CALL hdf5_read(le_ho,TRIM(filename),"mesh/ho_info/LE",success)
        IF(.NOT.success)CALL oft_abort('Error reading quadratic edges','native_load_vmesh',__FILE__)
        IF(hdf5_field_exist(TRIM(filename),"mesh/ho_info/LF"))THEN
            CALL hdf5_field_get_sizes(TRIM(filename),"mesh/ho_info/LF",ndims,dim_sizes)
            ALLOCATE(lf_ho(4,dim_sizes(2)))
            CALL hdf5_read(lf_ho,TRIM(filename),"mesh/ho_info/LF",success)
            IF(.NOT.success)CALL oft_abort('Error reading quadratic faces','native_load_vmesh',__FILE__)
        END IF
    END IF
    !---Read periodicity information
    IF(ref_periodic)THEN
        IF(.NOT.hdf5_field_exist(TRIM(filename),"mesh/periodicity/nodes"))CALL oft_abort( &
        "Periodic nodeset not found in file","native_load_vmesh",__FILE__)
        CALL hdf5_field_get_sizes(TRIM(filename),"mesh/periodicity/nodes",ndims,dim_sizes)
        ALLOCATE(per_nodes(dim_sizes(1)))
        CALL hdf5_read(per_nodes,TRIM(filename),"mesh/periodicity/nodes",success)
        WRITE(*,'(2A,I8)')oft_indent,'Found periodic points',dim_sizes(1)
    END IF
END IF
!---
IF(reflect)THEN
    call mesh_global_resolution(mesh)
    call native_reflect(mesh,.1d0*mesh%hmin)
END IF
IF(oft_env%rank/=0)DEALLOCATE(mesh%r,mesh%lc,mesh%reg)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine native_load_vmesh
!---------------------------------------------------------------------------------
!> Read in t3d mesh file from file "filename"
!! - Read in T3D options from input file
!! - Read in mesh points and cells
!! - Read in surface IDs for CAD edges and faces
!---------------------------------------------------------------------------------
subroutine native_load_smesh(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
logical :: is_2d,success
integer(i4) :: i,id,lenreflag,ierr,io_unit,ndims,np_mem,mesh_order
integer(i4), allocatable, dimension(:) :: dim_sizes
real(r8), allocatable, dimension(:,:) :: rtmp
LOGICAL :: reflect = .FALSE.
LOGICAL :: ref_periodic = .FALSE.
class(oft_bmesh), pointer :: smesh
!---Read in mesh options
namelist/native_mesh_options/filename,reflect,ref_periodic!,zstretch
DEBUG_STACK_PUSH
!---
np_mem=-1
filename='none'
IF(oft_env%head_proc)THEN
    IF(ALLOCATED(r_mem).AND.ALLOCATED(lc_mem).AND.ALLOCATED(reg_mem))THEN
        np_mem=SIZE(r_mem,DIM=2)
    ELSE
        OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
        READ(io_unit,native_mesh_options,IOSTAT=ierr)
        CLOSE(io_unit)
        IF(ierr<0)CALL oft_abort('No "native_mesh_options" found in input file.','native_load_smesh',__FILE__)
        IF(ierr>0)CALL oft_abort('Error parsing "native_mesh_options" in input file.','native_load_smesh',__FILE__)
        IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','native_load_smesh',__FILE__)
        !
        WRITE(*,'(2A)')oft_indent,'Native surface mesh:'
        CALL oft_increase_indent
        WRITE(*,'(3A)')oft_indent,'Filename = ',TRIM(filename)
    END IF
ELSE
    CALL oft_increase_indent
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(np_mem,1,OFT_MPI_I4,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(filename,80,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(reflect,1,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
CALL MPI_Bcast(ref_periodic,1,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
! CALL MPI_Bcast(zstretch, 1,OFT_MPI_R8,0,oft_env%COMM,ierr)
! IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','native_load_vmesh',__FILE__)
#endif
!---Read in mesh sizes
IF(np_mem>0)THEN
    ALLOCATE(dim_sizes(2))
    dim_sizes=SHAPE(lc_mem)
ELSE
    CALL hdf5_field_get_sizes(TRIM(filename),"mesh/LC",ndims,dim_sizes)
    IF(ndims==-1)CALL oft_abort('"mesh/LC" field does not exist in input file', 'native_load_smesh', __FILE__)
END IF
IF(dim_sizes(2)==0)CALL oft_abort('Zero length cell list','native_load_smesh',__FILE__)
!---Allocate mesh
SELECT CASE(dim_sizes(1))
    CASE(3)
        ALLOCATE(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
    CASE(4)
        ALLOCATE(oft_quadmesh::mg_mesh%smeshes(mg_mesh%mgdim))
    CASE DEFAULT
        CALL oft_abort('Invalid cell type','native_load_smesh',__FILE__)
END SELECT
DO i=1,mg_mesh%mgdim
    CALL mg_mesh%smeshes(i)%setup(mesh_native_id,.FALSE.)
END DO
CALL multigrid_level(mg_mesh,1)
smesh=>mg_mesh%smeshes(1)
smesh%nc=dim_sizes(2)
!
DEALLOCATE(dim_sizes)
IF(np_mem>0)THEN
    ALLOCATE(dim_sizes(2))
    dim_sizes=SHAPE(r_mem)
ELSE
    CALL hdf5_field_get_sizes(TRIM(filename),"mesh/R",ndims,dim_sizes)
    IF(ndims==-1)CALL oft_abort('"mesh/R" field does not exist in input file', 'native_load_smesh', __FILE__)
END IF
is_2d=(dim_sizes(1)==2)
smesh%np=dim_sizes(2)
IF(smesh%np==0)CALL oft_abort('Zero length point list','native_load_smesh',__FILE__)
DEALLOCATE(dim_sizes)
!---Read in points
ALLOCATE(smesh%r(3,smesh%np))
IF(is_2d)THEN
    smesh%dim=2
    smesh%r=0.d0
    IF(np_mem>0)THEN
        smesh%r(1:2,:)=r_mem
        DEALLOCATE(r_mem)
    ELSE
        ALLOCATE(rtmp(2,smesh%np))
        CALL hdf5_read(rtmp,TRIM(filename),"mesh/R",success)
        IF(.NOT.success)CALL oft_abort('Error reading points','native_load_smesh',__FILE__)
        smesh%r(1:2,:)=rtmp
        DEALLOCATE(rtmp)
    END IF
ELSE
    IF(np_mem>0)THEN
        smesh%r=r_mem
        DEALLOCATE(r_mem)
    ELSE
        CALL hdf5_read(smesh%r,TRIM(filename),"mesh/R",success)
        IF(.NOT.success)CALL oft_abort('Error reading points','native_load_smesh',__FILE__)
    END IF
END IF
!---Read in cells
ALLOCATE(smesh%lc(smesh%cell_np,smesh%nc))
IF(np_mem>0)THEN
    smesh%lc=lc_mem
    DEALLOCATE(lc_mem)
ELSE
    CALL hdf5_read(smesh%lc,TRIM(filename),"mesh/LC",success)
    IF(.NOT.success)CALL oft_abort('Error reading cells','native_load_smesh',__FILE__)
END IF
!---Read in region list and other information
ALLOCATE(smesh%reg(smesh%nc))
mesh_order=1
IF(np_mem>0)THEN
    smesh%reg=reg_mem
    DEALLOCATE(reg_mem)
ELSE
    CALL hdf5_read(smesh%reg,TRIM(filename),"mesh/REG",success)
    IF(.NOT.success)CALL oft_abort('Error reading region list','native_load_smesh',__FILE__)
    !---Read high-order mesh information
    IF(hdf5_field_exist(TRIM(filename),"mesh/ho_info/LE"))THEN
        mesh_order=2
        CALL hdf5_field_get_sizes(TRIM(filename),"mesh/ho_info/R",ndims,dim_sizes)
        IF(ndims==-1)CALL oft_abort('"mesh/ho_info/R" field does not exist in input file', 'native_load_smesh', __FILE__)
        np_ho=dim_sizes(2)
        ALLOCATE(r_ho(3,np_ho))
        IF(dim_sizes(1)==2)THEN
            ALLOCATE(rtmp(2,np_ho))
            CALL hdf5_read(rtmp,TRIM(filename),"mesh/ho_info/R",success)
            r_ho(1:2,:)=rtmp
            DEALLOCATE(rtmp)
        ELSE
            CALL hdf5_read(r_ho,TRIM(filename),"mesh/ho_info/R",success)
        END IF
        IF(.NOT.success)CALL oft_abort('Error reading quadratic points','native_load_smesh',__FILE__)
        ALLOCATE(le_ho(2,np_ho))
        CALL hdf5_read(le_ho,TRIM(filename),"mesh/ho_info/LE",success)
        IF(.NOT.success)CALL oft_abort('Error reading quadratic edges','native_load_smesh',__FILE__)
    END IF
    !---Read periodicity information
    IF(ref_periodic)THEN
        IF(.NOT.hdf5_field_exist(TRIM(filename),"mesh/periodicity/nodes"))CALL oft_abort( &
        "Periodic nodeset not found in file","native_load_smesh",__FILE__)
        CALL hdf5_field_get_sizes(TRIM(filename),"mesh/periodicity/nodes",ndims,dim_sizes)
        ALLOCATE(per_nodes(dim_sizes(1)))
        CALL hdf5_read(per_nodes,TRIM(filename),"mesh/periodicity/nodes",success)
        WRITE(*,'(2A,I8)')oft_indent,'Found periodic points',dim_sizes(1)
    END IF
END IF
!---
IF(reflect)THEN
    call mesh_global_resolution(smesh)
    call native_reflect(smesh,.1d0*smesh%hmin)
END IF
IF(oft_env%rank/=0)DEALLOCATE(smesh%r,smesh%lc,smesh%reg)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine native_load_smesh
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
subroutine native_read_nodesets(nsets,native_filename)
TYPE(oft_1d_int), pointer, intent(inout) :: nsets(:)
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: native_filename
INTEGER(4) :: j,num_nsets,ndims
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: dim_sizes
LOGICAL :: success
CHARACTER(LEN=4) :: blknum
CHARACTER(LEN=OFT_PATH_SLEN) :: local_filename
!---Open mesh file
local_filename=filename
IF(PRESENT(native_filename))local_filename=native_filename
! Look for sideset field
IF(.NOT.hdf5_field_exist(local_filename,"mesh/NUM_NODESETS"))RETURN
!---Read nodesets
CALL hdf5_read(num_nsets,local_filename,"mesh/NUM_NODESETS",success=success)
IF(.NOT.success)RETURN
ALLOCATE(nsets(num_nsets))
DO j=1,num_nsets
  WRITE(blknum,'(I4.4)')j
  CALL hdf5_field_get_sizes(local_filename,"mesh/NODESET"//blknum,ndims,dim_sizes)
  nsets(j)%n=dim_sizes(1)
  DEALLOCATE(dim_sizes)
  ALLOCATE(nsets(j)%v(nsets(j)%n))
  CALL hdf5_read(nsets(j)%v,local_filename,"mesh/NODESET"//blknum,success=success)
  IF(.NOT.success)nsets(j)%n=-1
END DO
CALL oft_decrease_indent
END SUBROUTINE native_read_nodesets
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
subroutine native_read_sidesets(ssets,native_filename)
TYPE(oft_1d_int), pointer, intent(inout) :: ssets(:)
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: native_filename
INTEGER(4) :: j,num_ssets,ndims
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: dim_sizes
LOGICAL :: success
CHARACTER(LEN=4) :: blknum
CHARACTER(LEN=OFT_PATH_SLEN) :: local_filename
!---Open mesh file
local_filename=filename
IF(PRESENT(native_filename))local_filename=native_filename
! Look for sideset field
IF(.NOT.hdf5_field_exist(local_filename,"mesh/NUM_SIDESETS"))RETURN
!---Read nodesets
CALL hdf5_read(num_ssets,local_filename,"mesh/NUM_SIDESETS",success=success)
IF(.NOT.success)RETURN
ALLOCATE(ssets(num_ssets))
DO j=1,num_ssets
    WRITE(blknum,'(I4.4)')j
    CALL hdf5_field_get_sizes(local_filename,"mesh/SIDESET"//blknum,ndims,dim_sizes)
    ssets(j)%n=dim_sizes(1)
    DEALLOCATE(dim_sizes)
    ALLOCATE(ssets(j)%v(ssets(j)%n))
    CALL hdf5_read(ssets(j)%v,local_filename,"mesh/SIDESET"//blknum,success=success)
    IF(.NOT.success)ssets(j)%n=-1
END DO
CALL oft_decrease_indent
END SUBROUTINE native_read_sidesets
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points from high order import
!---------------------------------------------------------------------------------
subroutine native_hobase(self)
CLASS(oft_amesh), INTENT(inout) :: self
IF(np_ho==0)RETURN
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Importing quadratic mesh nodes'
SELECT CASE(self%type)
CASE(1) ! Tet/Tri
    CALL native_hobase_simplex(self)
CASE(3) ! Hex/Quad
  SELECT TYPE(self)
    CLASS IS(oft_bmesh)
      CALL native_hobase_quad(self)
    CLASS IS(oft_mesh)
      CALL native_hobase_hex(self)
  END SELECT
CASE DEFAULT
  CALL oft_abort("Unknown element type", "native_hobase", __FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine native_hobase
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points to a simplex mesh (Tet/Tri) from high order import
!---------------------------------------------------------------------------------
subroutine native_hobase_simplex(self)
CLASS(oft_amesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2)
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
DO i=1,np_ho
    etmp=le_ho(:,i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','native_hobase_simplex',__FILE__)
    self%ho_info%r(:,k) = r_ho(:,i)
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,le_ho)
DEBUG_STACK_POP
end subroutine native_hobase_simplex
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points to a simplex mesh (Tet/Tri) from high order import
!---------------------------------------------------------------------------------
subroutine native_hobase_quad(self)
CLASS(oft_bmesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2)
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
!---Set edge midpoints from imported list
DO i=1,self%ne
    etmp=le_ho(:,i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','native_hobase_quad',__FILE__)
    self%ho_info%r(:,k) = r_ho(:,i)
END DO
!---Set cell centerpoints from imported list
DO i=1,self%nc
  self%ho_info%lcp(1,i)=i+self%ne
  self%ho_info%r(:,i+self%ne)=(self%r(:,self%lc(1,i))+self%r(:,self%lc(2,i)) &
    + self%r(:,self%lc(3,i))+self%r(:,self%lc(4,i)))/4.d0
  self%ho_info%r(:,i+self%ne) = r_ho(:,self%ne+i)
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,le_ho)
DEBUG_STACK_POP
end subroutine native_hobase_quad
!---------------------------------------------------------------------------------
!> Add quadratic mesh node points to a simplex mesh (Tet/Tri) from high order import
!---------------------------------------------------------------------------------
subroutine native_hobase_hex(self)
CLASS(oft_mesh), INTENT(inout) :: self
real(r8) :: pt(3)
integer(i4) :: i,j,k,etmp(2),ftmp(4)
DEBUG_STACK_PUSH
!---Setup quadratic mesh
self%order=1
self%ho_info%nep=1
self%ho_info%nfp=1
self%ho_info%ncp=1
ALLOCATE(self%ho_info%r(3,self%ne+self%nf+self%nc),self%ho_info%lep(self%ho_info%nep,self%ne))
ALLOCATE(self%ho_info%lfp(self%ho_info%nfp,self%nf),self%ho_info%lcp(self%ho_info%ncp,self%nc))
!---Initialize high order points with straight edges
DO i=1,self%ne
    self%ho_info%lep(1,i)=i
    self%ho_info%r(:,i)=(self%r(:,self%le(1,i))+self%r(:,self%le(2,i)))/2.d0
END DO
!---Set edge midpoints from imported list
DO i=1,self%ne
    etmp=le_ho(:,i)
    k=ABS(mesh_local_findedge(self,etmp))
    if(k==0)CALL oft_abort('Unlinked mesh edge','native_hobase_hex',__FILE__)
    self%ho_info%r(:,k) = r_ho(:,i)
END DO
!---Initialize high order points with flat faces
DO i=1,self%nf
    self%ho_info%lfp(1,i)=i+self%ne
    self%ho_info%r(:,i+self%ne)=(self%r(:,self%lf(1,i))+self%r(:,self%lf(2,i)) &
        + self%r(:,self%lf(3,i))+self%r(:,self%lf(4,i)))/4.d0
END DO
!---Set face centerpoints from imported list
DO i=1,self%nf
    ftmp=lf_ho(:,i)
    k=ABS(mesh_local_findface(self,ftmp))
    if(k==0)CALL oft_abort('Unlinked mesh face','native_hobase_hex',__FILE__)
    self%ho_info%r(:,self%ne+k) = r_ho(:,self%ne+i)
END DO
!---Set cell centerpoints from imported list
DO i=1,self%nc
    self%ho_info%lcp(1,i)=i+self%ne+self%nf
    self%ho_info%r(:,i+self%ne+self%nf)= (self%r(:,self%lc(1,i)) + self%r(:,self%lc(2,i)) &
        + self%r(:,self%lc(3,i)) + self%r(:,self%lc(4,i)) + self%r(:,self%lc(5,i)) &
        + self%r(:,self%lc(6,i)) + self%r(:,self%lc(7,i)) + self%r(:,self%lc(8,i)))/8.d0
    self%ho_info%r(:,i+self%ne+self%nf) = r_ho(:,i+self%ne+self%nf)
END DO
!---Destory temporary storage
DEALLOCATE(r_ho,le_ho,lf_ho)
DEBUG_STACK_POP
end subroutine native_hobase_hex
!---------------------------------------------------------------------------------
!> Reflect an native mesh across the xy-plane
!---------------------------------------------------------------------------------
subroutine native_reflect(self,tol)
CLASS(oft_amesh), intent(inout) :: self !< Mesh to reflect
real(r8), intent(in) :: tol !< tol Tolerance for marking point as on the reflection plane
integer(i4) :: npold,neold,nfold,ncold,i,j,ic,is,cid_max,sid_max,nreg,npold_ho,ne_ho,nf_ho,np_per,ref_index
integer(i4), allocatable :: newindex(:),hoindex(:),regtmp(:),ltemp(:,:)
real(r8), allocatable :: rtemp(:,:),rlftemp(:,:),rctemp(:,:)
DEBUG_STACK_PUSH
SELECT TYPE(self)
CLASS IS(oft_bmesh)
  IF(self%dim==2)THEN
    ref_index=1
    IF(oft_debug_print(1))write(*,'(2A)')oft_indent,'Reflecting 2D surface mesh -> x'
  ELSE
    ref_index=3
    IF(oft_debug_print(1))write(*,'(2A)')oft_indent,'Reflecting 3D surface mesh -> z'
  END IF
CLASS IS(oft_mesh)
  ref_index=3
  IF(oft_debug_print(1))write(*,'(2A)')oft_indent,'Reflecting 3D volume mesh -> z'
END SELECT
CALL oft_increase_indent
!---Reflect points that are not on reflection plane
npold=self%np
allocate(newindex(2*self%np),rtemp(3,2*self%np))
rtemp(:,1:self%np)=self%r
deallocate(self%r)
do i=1,npold
    IF(ABS(rtemp(ref_index,i))<=tol)THEN
        rtemp(ref_index,i)=0.d0
        newindex(i)=i
    ELSE
        self%np=self%np+1
        rtemp(:,self%np) = rtemp(:,i)
        rtemp(ref_index,self%np) =-rtemp(ref_index,i)
        newindex(i)=self%np
    ENDIF
enddo
allocate(self%r(3,self%np))
self%r=rtemp(:,1:self%np)
deallocate(rtemp)
!---Reflect cells
ncold=self%nc
allocate(ltemp(self%cell_np,2*self%nc),regtmp(2*self%nc))
ltemp(:,1:self%nc)=self%lc
regtmp(1:self%nc)=self%reg
nreg = MAXVAL(self%reg)
deallocate(self%lc,self%reg)
do i=1,ncold
    self%nc=self%nc+1
    DO j=1,self%cell_np
        ltemp(j,self%nc)=newindex(ltemp(j,i))
    END DO
    regtmp(self%nc)=regtmp(i)+nreg
enddo
allocate(self%lc(self%cell_np,self%nc),self%reg(self%nc))
self%lc=ltemp(:,1:self%nc)
self%reg=regtmp
deallocate(ltemp,regtmp)
!---Reflect high-order nodes
IF(np_ho>0)THEN
    npold=SIZE(r_ho,DIM=2)
    !---Reflect edges
    neold=SIZE(le_ho,DIM=2)
    allocate(ltemp(2,2*neold),rtemp(3,2*neold))
    ltemp(:,1:neold)=le_ho
    rtemp(:,1:neold)=r_ho(:,1:neold)
    deallocate(le_ho)
    ne_ho=neold
    DO i=1,neold
        IF(ABS(rtemp(ref_index,i))<=tol)THEN
            rtemp(ref_index,i)=0.d0
        ELSE
            ne_ho=ne_ho+1
            ltemp(:,ne_ho) = newindex(ltemp(:,i))
            rtemp(:,ne_ho) = rtemp(:,i)
            rtemp(ref_index,ne_ho) =-rtemp(ref_index,i)
        ENDIF
    END DO
    allocate(le_ho(2,ne_ho))
    le_ho=ltemp(:,1:ne_ho)
    deallocate(ltemp)
    np_ho = ne_ho
    SELECT TYPE(self)
      CLASS IS(oft_mesh)
        IF(ALLOCATED(lf_ho))THEN
            nfold=SIZE(lf_ho,DIM=2)
            allocate(ltemp(self%face_np,2*nfold),rlftemp(3,2*nfold))
            ltemp(:,1:nfold)=lf_ho
            rlftemp(:,1:nfold)=r_ho(:,neold+1:neold+nfold)
            deallocate(lf_ho)
            nf_ho=nfold
            DO i=1,nfold
                IF(ABS(rlftemp(3,i))<=tol)THEN
                    rlftemp(3,i)=0.d0
                ELSE
                    nf_ho=nf_ho+1
                    ltemp(:,nf_ho) = newindex(ltemp(:,i))
                    rlftemp(:,nf_ho) = rlftemp(:,i)
                    rlftemp(3,nf_ho) =-rlftemp(3,i)
                ENDIF
            END DO
            allocate(lf_ho(self%face_np,nf_ho))
            lf_ho=ltemp(:,1:nf_ho)
            deallocate(ltemp)
            np_ho = ne_ho+nf_ho+self%nc
            allocate(rctemp(3,np_ho))
            rctemp(:,1:ne_ho)=rtemp(:,1:ne_ho)
            rctemp(:,ne_ho+1:ne_ho+nf_ho) = rlftemp(:,1:nf_ho)
            deallocate(rtemp,rlftemp)
            rctemp(:,ne_ho+nf_ho+1:ne_ho+nf_ho+ncold) = r_ho(:,neold+nfold+1:neold+nfold+ncold)
            DO i=1,ncold
                rctemp(:,ne_ho+nf_ho+ncold+i) = rctemp(:,ne_ho+nf_ho+i)
                rctemp(3,ne_ho+nf_ho+ncold+i) =-rctemp(3,ne_ho+nf_ho+i)
            END DO
            allocate(rtemp(3,np_ho))
            rtemp=rctemp
            deallocate(rctemp)
        END IF
      CLASS IS(oft_bmesh)
        IF(npold>neold)THEN
            np_ho = ne_ho+self%nc
            allocate(rctemp(ref_index,np_ho))
            rctemp(:,1:ne_ho)=rtemp(:,1:ne_ho)
            deallocate(rtemp)
            rctemp(:,ne_ho+1:ne_ho+ncold) = r_ho(:,neold+1:neold+ncold)
            DO i=1,ncold
                rctemp(:,ne_ho+ncold+i) = rctemp(:,ne_ho+i)
                rctemp(ref_index,ne_ho+ncold+i) =-rctemp(ref_index,ne_ho+i)
            END DO
            allocate(rtemp(3,np_ho))
            rtemp=rctemp
            deallocate(rctemp)
        END IF
    END SELECT
    !---
    deallocate(r_ho)
    allocate(r_ho(3,np_ho))
    r_ho = rtemp(:,1:np_ho)
    deallocate(rtemp)
END IF
!---Flag periodic points
IF(ALLOCATED(per_nodes))THEN
    np_per=SIZE(per_nodes)
    ALLOCATE(self%periodic%lp(self%np))
    self%periodic%lp=-1
    DO i=1,np_per
        j=per_nodes(i)
        self%periodic%lp(newindex(j))=j
    END DO
END IF
deallocate(newindex)
! IF(oft_debug_print(3))write(*,'(A)')oft_indent,'Done'
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine native_reflect
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine native_set_periodic(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,pt_e(2),ind,k,kk,np_per
integer(i4), ALLOCATABLE :: pt_f(:)
IF(.NOT.ALLOCATED(per_nodes))RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting native volume mesh periodicity'
CALL oft_increase_indent
IF(.NOT.ASSOCIATED(mesh%periodic%lp))THEN
    ALLOCATE(mesh%periodic%lp(mesh%np))
    mesh%periodic%lp=-1
    np_per=SIZE(per_nodes)
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
end subroutine native_set_periodic
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine native_bset_periodic(mesh)
class(oft_bmesh), intent(inout) :: mesh
integer(i4) :: i,j,pt_e(2),ind,k,kk,np_per,ref_index
real(r8) :: pt1(3),pt2(3)
IF(.NOT.ALLOCATED(per_nodes))RETURN
DEBUG_STACK_PUSH
IF(mesh%dim==2)THEN
  ref_index=1
ELSE
  ref_index=3
END IF
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Setting native surface mesh periodicity'
CALL oft_increase_indent
IF(.NOT.ASSOCIATED(mesh%periodic%lp))THEN
  ALLOCATE(mesh%periodic%lp(mesh%np))
  mesh%periodic%lp=-1
  np_per=SIZE(per_nodes)
  DO i=1,np_per
    j=per_nodes(i)
    pt1=mesh%r(:,j)
    pt1(ref_index)=0.d0
    DO k=1,mesh%nbp
      kk=mesh%lbp(k)
      IF(kk==j)CYCLE
      pt2=mesh%r(:,kk)
      pt2(ref_index)=0.d0
      IF(ALL(ABS(pt1-pt2)<1.d-8).AND.ABS(mesh%r(ref_index,kk))<1.d-8)THEN
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
mesh%periodic%le=-1
!---Flag periodic edges
!$omp parallel private(j,pt_e,ind)
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
!$omp end parallel
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine native_bset_periodic
end module oft_mesh_native