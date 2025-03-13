!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_t3d.F90
!
!> Mesh handling for T3D meshes
!!
!! Functions to read in, intialize and refine T3D based geometry
!! - Mesh file read-in
!! - CAD file read-in
!! - Boundary point location and refinement
!!
!! @authors George Marklin and Chris Hansen
!! @date April 2008 - Present
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_t3d
USE oft_base
USE bezier_cad, ONLY: cad_vertex, cad_curve, cad_surf, cad_entity_ptr, &
  cad_curve_midpoint, cad_surf_midpoint, cad_surf_center, cad_curve_cast, cad_surf_cast
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_mesh_local_util, ONLY: mesh_local_findedge, mesh_local_findface
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE multigrid, ONLY: multigrid_mesh, multigrid_level
!---End include modules
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> T3D CAD boundary structure
!! - CAD entity counts
!! - CAD wireframe entities
!------------------------------------------------------------------------------
type :: T3D_cadgeom
  integer(i4) :: ngme = 0 !< Number of geometry model edges
  integer(i4) :: ngmf = 0 !< Number of geometry model faces
  integer(i4) :: ngwv = 0 !< number of CAD wireframe vertices
  integer(i4) :: ngwc = 0 !< number of CAD wireframe curves
  integer(i4) :: ngws = 0 !< number of CAD wireframe surfaces
  integer(i4), pointer, dimension(:,:) :: lgme => NULL() !< List of geometry model edges
  integer(i4), pointer, dimension(:,:) :: lgmf => NULL() !< List of geometry model faces
  integer(i4), pointer, dimension(:) :: lgmf_par => NULL() !< List of geometry model faces
  type(cad_vertex), pointer, dimension(:) :: vert => NULL() !< List of CAD wireframe vertices
  type(cad_curve), pointer, dimension(:) :: curve => NULL() !< List of CAD wireframe curves
  type(cad_surf), pointer, dimension(:) :: surf => NULL() !< List of CAD wireframe surfaces
end type T3D_cadgeom
!------------------------------------------------------------------------------
!> T3D CAD linkage structure
!! - Linkage of mesh entities to CAD model
!------------------------------------------------------------------------------
type :: T3D_cadlink
  type(cad_entity_ptr), pointer, dimension(:) :: lbeg => NULL() !< Linkage of mesh edges to CAD entities
  type(cad_entity_ptr), pointer, dimension(:) :: lbfg => NULL() !< Linkage of mesh faces to CAD entities
end type T3D_cadlink
!------------------------------------------------------------------------------
!> T3D raw CAD structure
!------------------------------------------------------------------------------
type :: T3D_cad
  integer(i4) :: ngwv = 0 !< number of CAD wireframe vertices
  integer(i4) :: ngww = 0 !< number of CAD wireframe weight points
  integer(i4) :: ngwc = 0 !< number of CAD wireframe curves
  integer(i4) :: ngws = 0 !< number of CAD wireframe surfaces
  integer(i4) :: ngwp = 0 !< number of CAD wireframe patches
  integer(i4) :: nv = 0 !< number of unique CAD vertices
  integer(i4) :: nc = 0 !< number of unique CAD curves
  integer(i4) :: ns = 0 !< number of unique CAD surfaces
  integer(i4) :: np = 0 !< number of unique CAD patches
  integer(i4), pointer, dimension(:,:) :: lgwc => NULL() !< List of CAD wireframe curves
  integer(i4), pointer, dimension(:,:) :: lgws => NULL() !< List of CAD wireframe surfaces
  integer(i4), pointer, dimension(:) :: vtmp => NULL() !< Vertex index mapping
  integer(i4), pointer, dimension(:) :: ctmp => NULL() !< Curve index mapping
  integer(i4), pointer, dimension(:) :: stmp => NULL() !< Surface index mapping
  integer(i4), pointer, dimension(:) :: ptmp => NULL() !< Patch index mapping
  real(r8), pointer, dimension(:,:) :: lgwv => NULL() !< List of CAD wireframe vertices
  real(r8), pointer, dimension(:,:) :: lgww => NULL() !< List of CAD wireframe weight points
end type T3D_cad
INTEGER(i4), PARAMETER, PUBLIC :: mesh_t3d_id = 1
PUBLIC mesh_t3d_load, mesh_t3d_reffix, mesh_t3d_cadsync, mesh_t3d_set_periodic
PUBLIC mesh_t3d_cadlink, mesh_t3d_add_quad, smesh_t3d_load
!---Global variables
character(LEN=OFT_PATH_SLEN) :: inpname = 'none' !< Name of T3D input file for geometry (Used to retrieve CAD objects)
character(LEN=OFT_PATH_SLEN) :: filename = 'none' !< Name of T3D input file for mesh
character(LEN=3) :: reflect = '' !< Character flag for mesh reflection (eg. 'xy')
LOGICAL :: ref_per(3) = .FALSE. !< Character flag for periodic reflections
real(r8) :: zstretch = 1.d0 !< Scale for z-coordinates (useful for cylindrical pinch studies)
type(T3D_cad) :: cad_tmp !< Raw CAD information from "inpname"
type(T3D_cadgeom), pointer :: cad_rep => NULL() !< CAD representation
type(T3D_cadlink), pointer :: cad_link => NULL() !< Linkage of mesh to CAD geometry
type(T3D_cadgeom), pointer, dimension(:) :: ML_cad_rep => NULL() !< ML CAD representation array
type(T3D_cadlink), pointer, dimension(:) :: ML_cad_link => NULL() !< ML CAD linkage
contains
!------------------------------------------------------------------------------
!> Read in t3d mesh file from file "filename"
!! - Read in T3D options from input file
!! - Read in mesh points and cells
!! - Read in surface IDs for CAD edges and faces
!------------------------------------------------------------------------------
subroutine mesh_t3d_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,id,lenreflag,ierr,io_unit
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
!---Read in mesh options
namelist/t3d_options/filename,inpname,reflect,ref_per,zstretch
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,t3d_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "t3d_options" found in input file.','mesh_t3d_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "t3d_options" in input file.','mesh_t3d_load',__FILE__)
  IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','mesh_t3d_load',__FILE__)
  IF(TRIM(inpname)=='none')CALL oft_abort('No T3D input file specified','mesh_t3d_load',__FILE__)
  lenreflag=lnblnk(reflect)
  !
  WRITE(*,*)
  WRITE(*,'(A)')'**** Loading T3D mesh'
  WRITE(*,'(2X,2A)')'Mesh File = ',TRIM(filename)
  WRITE(*,'(2X,2A)')'Geom File = ',TRIM(inpname)
  IF(lenreflag>0)THEN
    WRITE(*,'(2X,A)')'Reflection:'
    DO i=1,lenreflag
      SELECT CASE(reflect(i:i))
        CASE('x')
          WRITE(*,'(4X,A,L)')'YZ-plane, periodic = ',ref_per(1)
        CASE('y')
          WRITE(*,'(4X,A,L)')'XZ-plane, periodic = ',ref_per(2)
        CASE('z')
          WRITE(*,'(4X,A,L)')'XY-plane, periodic = ',ref_per(3)
      END SELECT
    END DO
  END IF
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(filename,40,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(inpname, 40,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(reflect,  3,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(ref_per,  3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(zstretch, 1,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
#endif
!---Allocate mesh
allocate(oft_tetmesh::mg_mesh%meshes(mg_mesh%mgdim))
allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%meshes(i)%setup(mesh_t3d_id)
  CALL mg_mesh%smeshes(i)%setup(mesh_t3d_id,.TRUE.)
  mg_mesh%meshes(i)%bmesh=>mg_mesh%smeshes(i)
END DO
CALL multigrid_level(mg_mesh,1)
mesh=>mg_mesh%meshes(1)
smesh=>mg_mesh%smeshes(1)
IF(.NOT.oft_env%head_proc)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---Setup geometry information
allocate(ML_cad_rep(mg_mesh%mgdim))
cad_rep=>ML_cad_rep(mg_mesh%level)
allocate(ML_cad_link(mg_mesh%mgdim))
cad_link=>ML_cad_link(mg_mesh%level)
!---Read in T3D mesh
open(NEWUNIT=io_unit,FILE=trim(filename))
read(io_unit,*)
read(io_unit,*)mesh%np,cad_rep%ngme,cad_rep%ngmf,mesh%nc
read(io_unit,*)
!---Read in points
if(mesh%np==0)call oft_abort('No points in T3D file','mesh_t3d_load',__FILE__)
allocate(mesh%r(3,mesh%np))
do i=1,mesh%np
  read(io_unit,*)id,mesh%r(1,i),mesh%r(2,i),mesh%r(3,i)
  mesh%r(3,i)=mesh%r(3,i)*zstretch
end do
read(io_unit,*)
!---Read in CAD model edges
if(cad_rep%ngme>0)then
  allocate(cad_rep%lgme(4,cad_rep%ngme))
  do i=1,cad_rep%ngme
    read(io_unit,*)id,cad_rep%lgme(:,i)
  end do
  read(io_unit,*)
end if
!---Read in CAD model faces
if(cad_rep%ngmf>0)then
  allocate(cad_rep%lgmf(5,cad_rep%ngmf),cad_rep%lgmf_par(cad_rep%ngmf))
  cad_rep%lgmf_par=0
  do i=1,cad_rep%ngmf
    read(io_unit,*)id,cad_rep%lgmf(:,i)
  end do
  read(io_unit,*)
end if
!---Read in cells
if(mesh%nc==0)call oft_abort('No cells in T3D file','mesh_t3d_load',__FILE__)
allocate(mesh%lc(4,mesh%nc),mesh%reg(mesh%nc))
mesh%reg=1
do i=1,mesh%nc
  read(io_unit,*)id,mesh%lc(1,i),mesh%lc(2,i),mesh%lc(3,i),mesh%lc(4,i)
end do
close(io_unit)
IF(oft_debug_print(2))WRITE(*,*)'  Complete'
!---
call mesh_t3d_geom
call mesh_global_resolution(mesh)
!---
do i=1,lenreflag
  if(reflect(i:i)=='x')call mesh_t3d_reflect(mesh,1,.1d0*mesh%hmin,ref_per(i))
  if(reflect(i:i)=='y')call mesh_t3d_reflect(mesh,2,.1d0*mesh%hmin,ref_per(i))
  if(reflect(i:i)=='z')call mesh_t3d_reflect(mesh,3,.1d0*mesh%hmin,ref_per(i))
end do
DEBUG_STACK_POP
end subroutine mesh_t3d_load
!------------------------------------------------------------------------------
!> Read in t3d mesh file from file "filename"
!! - Read in T3D options from input file
!! - Read in mesh points and cells
!! - Read in surface IDs for CAD edges and faces
!------------------------------------------------------------------------------
subroutine smesh_t3d_load(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,id,lenreflag,ierr,io_unit
class(oft_bmesh), pointer :: smesh
!---Read in mesh options
namelist/t3d_options/filename,inpname,reflect,ref_per,zstretch
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
  READ(io_unit,t3d_options,IOSTAT=ierr)
  CLOSE(io_unit)
  IF(ierr<0)CALL oft_abort('No "t3d_options" found in input file.','smesh_t3d_load',__FILE__)
  IF(ierr>0)CALL oft_abort('Error parsing "t3d_options" in input file.','smesh_t3d_load',__FILE__)
  IF(TRIM(filename)=='none')CALL oft_abort('No mesh file specified','smesh_t3d_load',__FILE__)
  IF(TRIM(inpname)=='none')CALL oft_abort('No T3D input file specified','smesh_t3d_load',__FILE__)
  lenreflag=lnblnk(reflect)
  !
  WRITE(*,*)
  WRITE(*,'(A)')'**** Loading T3D mesh'
  WRITE(*,'(2X,2A)')'Mesh File = ',TRIM(filename)
  WRITE(*,'(2X,2A)')'Geom File = ',TRIM(inpname)
  IF(lenreflag>0)THEN
    WRITE(*,'(2X,A)')'Reflection:'
    DO i=1,lenreflag
      SELECT CASE(reflect(i:i))
        CASE('x')
          WRITE(*,'(4X,A,L)')'YZ-plane, periodic = ',ref_per(1)
        CASE('y')
          WRITE(*,'(4X,A,L)')'XZ-plane, periodic = ',ref_per(2)
        CASE('z')
          WRITE(*,'(4X,A,L)')'XY-plane, periodic = ',ref_per(3)
      END SELECT
    END DO
  END IF
END IF
!---Broadcast input information
#ifdef HAVE_MPI
CALL MPI_Bcast(filename,40,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(inpname, 40,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(reflect,  3,OFT_MPI_CHAR,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(ref_per,  3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
CALL MPI_Bcast(zstretch, 1,OFT_MPI_R8,0,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_Bcast','mesh_t3d_load',__FILE__)
#endif
!---Allocate mesh
allocate(oft_trimesh::mg_mesh%smeshes(mg_mesh%mgdim))
DO i=1,mg_mesh%mgdim
  CALL mg_mesh%smeshes(i)%setup(mesh_t3d_id,.FALSE.)
END DO
CALL multigrid_level(mg_mesh,1)
smesh=>mg_mesh%smeshes(1)
IF(.NOT.oft_env%head_proc)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---Setup geometry information
allocate(ML_cad_rep(mg_mesh%mgdim))
cad_rep=>ML_cad_rep(mg_mesh%level)
allocate(ML_cad_link(mg_mesh%mgdim))
cad_link=>ML_cad_link(mg_mesh%level)
!---Read in T3D mesh
open(NEWUNIT=io_unit,FILE=trim(filename))
read(io_unit,*)
read(io_unit,*)smesh%np,cad_rep%ngme,smesh%nc
read(io_unit,*)
!---Read in points
if(smesh%np==0)call oft_abort('No points in T3D file','smesh_t3d_load',__FILE__)
allocate(smesh%r(3,smesh%np))
do i=1,smesh%np
  read(io_unit,*)id,smesh%r(1,i),smesh%r(2,i),smesh%r(3,i)
  smesh%r(3,i)=smesh%r(3,i)*zstretch
end do
read(io_unit,*)
!---Read in CAD model edges
if(cad_rep%ngme>0)then
  allocate(cad_rep%lgme(4,cad_rep%ngme))
  do i=1,cad_rep%ngme
    read(io_unit,*)id,cad_rep%lgme(:,i)
  end do
  read(io_unit,*)
end if
cad_rep%ngmf=0 ! No faces in surface model
! !---Read in CAD model faces
! if(cad_rep%ngmf>0)then
!   allocate(cad_rep%lgmf(5,cad_rep%ngmf),cad_rep%lgmf_par(cad_rep%ngmf))
!   cad_rep%lgmf_par=0
!   do i=1,cad_rep%ngmf
!     read(io_unit,*)id,cad_rep%lgmf(:,i)
!   end do
!   read(io_unit,*)
! end if
!---Read in cells
if(smesh%nc==0)call oft_abort('No cells in T3D file','smesh_t3d_load',__FILE__)
allocate(smesh%lc(3,smesh%nc),smesh%reg(smesh%nc))
smesh%reg=1
do i=1,smesh%nc
  read(io_unit,*)id,smesh%lc(1,i),smesh%lc(2,i),smesh%lc(3,i)
end do
close(io_unit)
IF(oft_debug_print(2))WRITE(*,*)'  Complete'
!---
! call mesh_t3d_geom
call mesh_global_resolution(smesh)
! !---
! do i=1,lenreflag
!   if(reflect(i:i)=='x')call mesh_t3d_reflect(1,.1d0*mesh%hmin,ref_per(i))
!   if(reflect(i:i)=='y')call mesh_t3d_reflect(2,.1d0*mesh%hmin,ref_per(i))
!   if(reflect(i:i)=='z')call mesh_t3d_reflect(3,.1d0*mesh%hmin,ref_per(i))
! end do
DEBUG_STACK_POP
end subroutine smesh_t3d_load
!------------------------------------------------------------------------------
!> Read in T3D geometry information from file 'inpname'
!! - CAD vertices
!! - Weight points
!! - Curves
!! - Surfaces
!!
!! @note Limited to quadratic curves and dual quadratic surfaces.
!! Cubic curves and surfaces are not implemented at this time.
!------------------------------------------------------------------------------
subroutine mesh_t3d_geom
integer(i4) :: i,k,id,eo,eid(2),sid(4),io_unit
real(r8) :: cords(3),weight
character(2) :: test1,test2
character(40) :: t1,t2,t3
DEBUG_STACK_PUSH
!---Only lead process reads file
if(oft_env%head_proc)then
  if(oft_debug_print(1))WRITE(*,'(2X,A)')'Loading T3D geometry:'
  !---Open File for parsing
  open(NEWUNIT=io_unit,FILE=TRIM(inpname),ACTION='read')
  !---Initialize CAD counters
  cad_tmp%ngwv=0
  cad_tmp%ngww=0
  cad_tmp%ngwc=0
  cad_tmp%ngws=0
  cad_tmp%ngwp=0
  !---Count CAD objects
  do i=1,1000
    read(io_unit,*)test1 ! Get next line in file
    if(TRIM(test1)=='ve')then ! Vertex found
      backspace(io_unit)
      read(io_unit,*)t1,id
      cad_tmp%ngwv=max(id,cad_tmp%ngwv) ! Update vertex counter
    elseif(TRIM(test1)=='cu')then ! Curve found
      backspace(io_unit)
      read(io_unit,*)t1,id,t2
      if(TRIM(t2)=='order')then ! Higher order curve
        read(io_unit,*)test2
        if(test2=='po')cad_tmp%ngww=cad_tmp%ngww+1 ! Update weight counter
      endif
      cad_tmp%ngwc=max(id,cad_tmp%ngwc) ! Update curve counter
    elseif(TRIM(test1)=='su')then ! Surface found
      backspace(io_unit)
      read(io_unit,*)t1,id
      read(io_unit,*)test2
      backspace(io_unit)
      if(test2=='po')then ! High order surface
        cad_tmp%ngww=cad_tmp%ngww+1 ! Update weight counter
      endif
      cad_tmp%ngws=max(id,cad_tmp%ngws) ! Update surface counter
    elseif(TRIM(test1)=='pa')then ! Patch found
      backspace(io_unit)
      read(io_unit,*)t1,id
      cad_tmp%ngwp=max(id,cad_tmp%ngwp)
    elseif(TRIM(test1)=='re')then ! Region found (EOF)
      exit
    endif
  enddo
  !---Return to beginning of file to read in
  rewind(io_unit)
  !---Create CAD object arrays
  allocate(cad_tmp%lgwv(3,cad_tmp%ngwv),cad_tmp%lgww(4,cad_tmp%ngww))
  allocate(cad_tmp%lgwc(3,cad_tmp%ngwc),cad_tmp%lgws(5,cad_tmp%ngws))
  allocate(cad_tmp%vtmp(cad_tmp%ngwv),cad_tmp%ctmp(cad_tmp%ngwc))
  allocate(cad_tmp%stmp(cad_tmp%ngws),cad_tmp%ptmp(cad_tmp%ngwp))
  cad_tmp%lgwc=0; cad_tmp%lgwv=0; cad_tmp%lgww=0; cad_tmp%lgws=0
  cad_tmp%vtmp=0; cad_tmp%ctmp=0; cad_tmp%stmp=0; cad_tmp%ptmp=0
  cad_tmp%nv=0; cad_tmp%nc=0; cad_tmp%ns=0; cad_tmp%np=0
  !---Read in CAD objects
  k=1 ! Initialize weight counter
  do i=1,1000 ! Parse file reading in objects
    read(io_unit,*)test1
    if(TRIM(test1)=='ve')then ! Vertex
      cad_tmp%nv=cad_tmp%nv+1
      backspace(io_unit)
      read(io_unit,*)t1,id,t2
      backspace(io_unit)
      if(trim(t2)=='fixed')then ! Fixed to other vertex
        read(io_unit,*)t1,id,t2,t3,eo
        cad_tmp%lgwv(:,id)=cad_tmp%lgwv(:,eo)
      else ! New vertex coordinates
        read(io_unit,*)t1,id,t2,cords
        cords(3)=cords(3)*zstretch
        cad_tmp%lgwv(:,id)=cords
      endif
      cad_tmp%vtmp(id)=cad_tmp%nv
    elseif(TRIM(test1)=='cu')then ! Curve
      cad_tmp%nc=cad_tmp%nc+1
      backspace(io_unit)
      read(io_unit,*)t1,id,t2
      cad_tmp%lgwc(:,id)=0
      cad_tmp%ctmp(id)=cad_tmp%nc
      backspace(io_unit)
      if(TRIM(t2)=='order')then ! High order
        read(io_unit,*)t1,id,t2,eo,t3,eid
        cad_tmp%lgwc(1:2,id)=eid ! Get vertex IDs
        read(io_unit,*)test2
        if(test2=='po')then ! Weight information
          backspace(io_unit)
          cad_tmp%lgwc(3,id)=k
          read(io_unit,*)t1,id,t2,cords,t3,weight
          cords(3)=cords(3)*zstretch
          cad_tmp%lgww(1:3,k)=cords ! Weight position
          cad_tmp%lgww(4,k)=weight ! Weight value
          k=k+1 ! Increment weight point counter
        endif
      else ! Straight line
        read(io_unit,*)t1,id,t2,eid
        cad_tmp%lgwc(1:2,id)=eid ! get vertex IDs
      endif
    elseif(TRIM(test1)=='su')then ! Surface
      cad_tmp%ns=cad_tmp%ns+1
      backspace(io_unit)
      read(io_unit,*)t1,id,t2,sid
      cad_tmp%lgws(1:4,id)=sid ! Get curve IDs
      cad_tmp%stmp(id)=cad_tmp%ns
      read(io_unit,*)test2
      backspace(io_unit)
      if(test2=='po')then ! Weight information
        cad_tmp%lgws(5,id)=k
        read(io_unit,*)t1,id,eo,t2,cords,t3,weight
        cords(3)=cords(3)*zstretch
        cad_tmp%lgww(1:3,k)=cords ! Weight position
        cad_tmp%lgww(4,k)=weight ! Weight value
        k=k+1 ! Increment weight point counter
      endif
    elseif(TRIM(test1)=='pa')then ! Patch
      cad_tmp%np=cad_tmp%np+1
      backspace(io_unit)
      read(io_unit,*)t1,id
      cad_tmp%ptmp(id)=cad_tmp%np
    elseif(TRIM(test1)=='re')then ! Region found (EOF)
      exit
    endif
  enddo
  !---Close file
  close(io_unit)
  !---Write out CAD stats
  IF(oft_debug_print(1))THEN
    WRITE(*,'(4X,A,I8)')'# of vertices          =',cad_tmp%ngwv
    WRITE(*,'(4X,A,I8)')'# of weight points     =',cad_tmp%ngww
    WRITE(*,'(4X,A,I8)')'# of curves            =',cad_tmp%ngwc
    WRITE(*,'(4X,A,I8)')'# of surfaces          =',cad_tmp%ngws
    WRITE(*,'(4X,A,I8)')'# of patches           =',cad_tmp%ngwp
  END IF
end if
DEBUG_STACK_POP
end subroutine mesh_t3d_geom
!------------------------------------------------------------------------------
!> Synchronize T3D geometry information.
!------------------------------------------------------------------------------
subroutine mesh_t3d_cadsync(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: tmp(11),ierr
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2X,A)')'Syncing T3D geometry'
CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Communicate T3D boundary information
if(oft_env%rank==0)then
  tmp(1)=cad_rep%ngme
  tmp(2)=cad_rep%ngmf
  tmp(3)=cad_tmp%ngwv
  tmp(4)=cad_tmp%ngww
  tmp(5)=cad_tmp%ngwc
  tmp(6)=cad_tmp%ngws
  tmp(7)=cad_tmp%ngwp
  tmp(8)=cad_tmp%nv
  tmp(9)=cad_tmp%nc
  tmp(10)=cad_tmp%ns
  tmp(11)=cad_tmp%np
endif
#ifdef HAVE_MPI
call MPI_Bcast(tmp,11,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD object counts
#endif
if(oft_env%rank>0)then
  allocate(ML_cad_rep(mg_mesh%mgdim))
  cad_rep=>ML_cad_rep(mg_mesh%level)
  allocate(ML_cad_link(mg_mesh%mgdim))
  cad_link=>ML_cad_link(mg_mesh%level)
  cad_rep%ngme=tmp(1)
  cad_rep%ngmf=tmp(2)
  cad_tmp%ngwv=tmp(3)
  cad_tmp%ngww=tmp(4)
  cad_tmp%ngwc=tmp(5)
  cad_tmp%ngws=tmp(6)
  cad_tmp%ngwp=tmp(7)
  cad_tmp%nv=tmp(8)
  cad_tmp%nc=tmp(9)
  cad_tmp%ns=tmp(10)
  cad_tmp%np=tmp(11)
endif
if(oft_env%rank>0)then ! Allocate CAD object arrays
  allocate(cad_rep%lgme(4,cad_rep%ngme),cad_rep%lgmf(5,cad_rep%ngmf))
  allocate(cad_rep%lgmf_par(cad_rep%ngmf))
  allocate(cad_tmp%lgwv(3,cad_tmp%ngwv),cad_tmp%lgww(4,cad_tmp%ngww))
  allocate(cad_tmp%lgwc(3,cad_tmp%ngwc),cad_tmp%lgws(5,cad_tmp%ngws))
  allocate(cad_tmp%vtmp(cad_tmp%ngwv),cad_tmp%ctmp(cad_tmp%ngwc))
  allocate(cad_tmp%stmp(cad_tmp%ngws),cad_tmp%ptmp(cad_tmp%ngwp))
endif
#ifdef HAVE_MPI
call MPI_Bcast(cad_rep%lgme,4*cad_rep%ngme,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD edge list
call MPI_Bcast(cad_rep%lgmf,5*cad_rep%ngmf,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD face list
call MPI_Bcast(cad_rep%lgmf_par,cad_rep%ngmf,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD face parent list
call MPI_Bcast(cad_tmp%lgwv,3*cad_tmp%ngwv,OFT_MPI_R8,0,oft_env%COMM,ierr) ! Broadcast CAD vertex list
call MPI_Bcast(cad_tmp%lgww,4*cad_tmp%ngww,OFT_MPI_R8,0,oft_env%COMM,ierr) ! Broadcast CAD weight point list
call MPI_Bcast(cad_tmp%lgwc,3*cad_tmp%ngwc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD curve list
call MPI_Bcast(cad_tmp%lgws,5*cad_tmp%ngws,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD surface list
call MPI_Bcast(cad_tmp%vtmp,cad_tmp%ngwv,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD vertex list
call MPI_Bcast(cad_tmp%ctmp,cad_tmp%ngwc,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD curve list
call MPI_Bcast(cad_tmp%stmp,cad_tmp%ngws,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD surface list
call MPI_Bcast(cad_tmp%ptmp,cad_tmp%ngwp,OFT_MPI_I4,0,oft_env%COMM,ierr) ! Broadcast CAD patch list
call MPI_Bcast(ref_per,3,OFT_MPI_LOGICAL,0,oft_env%COMM,ierr) ! Broadcast periodicity flags
#endif
call mesh_t3d_cadconv
DEBUG_STACK_POP
end subroutine mesh_t3d_cadsync
!------------------------------------------------------------------------------
!> Convert T3D CAD representation to code represenation.
!! - Convert CAD entities to bezier_cad objects.
!------------------------------------------------------------------------------
subroutine mesh_t3d_cadconv
integer(i4) :: i,k
DEBUG_STACK_PUSH
!---Copy CAD vertices to code objects
allocate(cad_rep%vert(cad_tmp%nv))
do i=1,cad_tmp%ngwv
  k=cad_tmp%vtmp(i) ! Get compact index
  if(k>0)then ! Vertex found
    allocate(cad_rep%vert(k)%pt(3))
    cad_rep%vert(k)%pt=cad_tmp%lgwv(:,i)
  end if
end do
!---Copy CAD curves to code objects
allocate(cad_rep%curve(cad_tmp%nc))
do i=1,cad_tmp%ngwc
  k=cad_tmp%ctmp(i) ! Get compact index
  if(k>0)then ! Curve found
    if(cad_tmp%lgwc(3,i)>0)then ! Quadratic curve
      cad_rep%curve(k)%id=i
      cad_rep%curve(k)%order=3
      allocate(cad_rep%curve(k)%pt(3,3),cad_rep%curve(k)%wt(3))
      cad_rep%curve(k)%pt(:,1)=cad_tmp%lgwv(:,cad_tmp%lgwc(1,i))
      cad_rep%curve(k)%pt(:,2)=cad_tmp%lgww(1:3,cad_tmp%lgwc(3,i))
      cad_rep%curve(k)%pt(:,3)=cad_tmp%lgwv(:,cad_tmp%lgwc(2,i))
      cad_rep%curve(k)%wt(1)=1.d0
      cad_rep%curve(k)%wt(2)=cad_tmp%lgww(4,cad_tmp%lgwc(3,i))
      cad_rep%curve(k)%wt(3)=1.d0
    else ! Linear curve
      cad_rep%curve(k)%id=i
      cad_rep%curve(k)%order=2
      allocate(cad_rep%curve(k)%pt(3,2),cad_rep%curve(k)%wt(2))
      cad_rep%curve(k)%pt(:,1)=cad_tmp%lgwv(:,cad_tmp%lgwc(1,i))
      cad_rep%curve(k)%pt(:,2)=cad_tmp%lgwv(:,cad_tmp%lgwc(2,i))
      cad_rep%curve(k)%wt(1)=1.d0
      cad_rep%curve(k)%wt(2)=1.d0
    end if
  end if
end do
!---Create CAD curve grids
do i=1,cad_tmp%ngwc
  k=cad_tmp%ctmp(i) ! Get compact index
  if(k>0)call cad_rep%curve(k)%grid
end do
!---Copy CAD surfaces to code objects
allocate(cad_rep%surf(cad_tmp%ns))
do i=1,cad_tmp%ngws
  k=cad_tmp%stmp(i) ! Get compact index
  if(k>0)then
    cad_rep%surf(k)%id=i
    call mesh_t3d_surfconst(i,cad_rep%surf(k))
    call cad_rep%surf(k)%grid
  end if
end do
!---Destroy temporary arrays
deallocate(cad_tmp%lgwc,cad_tmp%lgws,cad_tmp%lgww,cad_tmp%lgwv)
DEBUG_STACK_POP
end subroutine mesh_t3d_cadconv
!------------------------------------------------------------------------------
!> Link T3D CAD objects to mesh entities for use in refinement.
!------------------------------------------------------------------------------
subroutine mesh_t3d_cadlink(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,ind,k,ep(2),fp(3),ind_par
integer(i4), allocatable :: emap(:),fmap(:)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2X,A)')'Linking T3D geometry to mesh'
!---Get edge boundary mapping
allocate(emap(mesh%ne))
allocate(cad_link%lbeg(mesh%nbe)) ! Create edge linkage
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne) ! Get edge map
!---Search CAD model edges for local edges
outer1:do i=1,cad_rep%ngme
  do j=1,2 ! Get end points for CAD edge
    if(cad_rep%lgme(j,i)==0)then
      cycle outer1 ! Skip empty entries
    endif
    ep(j)=cad_rep%lgme(j,i)
  enddo
  ind=abs(mesh_local_findedge(mesh,ep)) ! Find edge on mesh
  if(ind==0)cycle ! Not on my domain
  !if(emap(ind)==0)call oft_abort('mesh_t3d_cadlink','Bad edge link')
  if(emap(ind)==0)cycle
  k=cad_tmp%ctmp(cad_rep%lgme(4,i)) ! Get CAD curve index
  cad_link%lbeg(emap(ind))%wo=>cad_rep%curve(k) ! Link edge to CAD curve
enddo outer1
!---Get face boundary mapping
allocate(fmap(mesh%nf))
allocate(cad_link%lbfg(mesh%nbf)) ! Create face linkage
ALLOCATE(mesh%periodic%lf(mesh%nf))
mesh%periodic%lf=-1
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf) ! Get face map
!---Search CAD model faces for local faces
outer2:do i=1,cad_rep%ngmf
  !if(cad_rep%lgmf(4,i)==3)then ! Cycle if on patch
    do j=1,3 ! Get corner points for CAD face
      if(cad_rep%lgmf(j,i)==0)then
        cycle outer2 ! Skip empty entries
      endif
      fp(j)=cad_rep%lgmf(j,i)
    enddo
    ind=abs(mesh_local_findface(mesh,(/fp(1),fp(2),fp(3)/))) ! Find face on mesh
    if(ind==0)cycle ! Not on my domain
    if(fmap(ind)==0)cycle!call oft_abort('Bad face link','mesh_t3d_cadlink',__FILE__)
  if(cad_rep%lgmf(4,i)==3)then
    k=cad_tmp%stmp(cad_rep%lgmf(5,i)) ! Get CAD surface index
    cad_link%lbfg(fmap(ind))%wo=>cad_rep%surf(k) ! Link face to CAD surface
  elseif(cad_rep%lgmf(4,i)==5)then
    k=cad_tmp%ns+cad_tmp%ptmp(cad_rep%lgmf(5,i)) ! Get CAD surface index
  else
    k=-1
  endif
  mesh%bfs(fmap(ind))=k
  !---Get parent
  IF(cad_rep%lgmf_par(i)/=0)THEN
    k=cad_rep%lgmf_par(i)
    do j=1,3 ! Get corner points for CAD face
      if(cad_rep%lgmf(j,k)==0)cycle outer2 ! Skip empty entries
      fp(j)=cad_rep%lgmf(j,k)
    enddo
    ind_par=abs(mesh_local_findface(mesh,(/fp(1),fp(2),fp(3)/))) ! Find face on mesh
    if(ind_par==0)cycle ! Not on my domain
    if(fmap(ind_par)==0)cycle
    mesh%periodic%lf(ind)=ind_par
  END IF
enddo outer2
!---Link edges on surfaces from faces
do i=1,mesh%nbf
  if(.NOT.associated(cad_link%lbfg(i)%wo))cycle ! Face is not on surface
  do j=1,3 ! Check each edge on face
    k=abs(mesh%lfe(j,mesh%lbf(i)))
    if(emap(k)==0)call oft_abort('Bad edge link','mesh_t3d_cadlink',__FILE__)
    if(associated(cad_link%lbeg(emap(k))%wo))cycle ! Edge is already linked
    cad_link%lbeg(emap(k))%wo=>cad_link%lbfg(i)%wo ! Link edge to CAD surface
  enddo
enddo
!---Destroy temporary arrays
deallocate(cad_tmp%vtmp,cad_tmp%ctmp,cad_tmp%stmp)
!---Destroy mapping arrays
deallocate(emap,fmap)
DEBUG_STACK_POP
end subroutine mesh_t3d_cadlink
!------------------------------------------------------------------------------
!> Construct CAD surface object
!------------------------------------------------------------------------------
subroutine mesh_t3d_surfconst(si,surf)
integer(i4), intent(in) :: si !< Index of CAD surface to use as source
type(cad_surf), intent(out) :: surf !< Surface object
real(r8) :: wcheck(2),tmp(3),p(3,3,3),w(3,3)
integer(i4) :: k,n,m,cind,wind,ind,ctemp(2),ccheck(2)
logical :: flip
DEBUG_STACK_PUSH
!---Check if "si" is a valid surface
if(cad_tmp%lgws(1,si)>0)then
  !---Determine order of boundary curves
  do k=1,2
    wcheck(k)=cad_tmp%lgwc(3,cad_tmp%lgws(k,si))
  enddo
  n=3; m=3
  if(wcheck(1)==0.d0)m=2
  if(wcheck(2)==0.d0)n=2
  surf%order=(/n,m/)
  allocate(surf%pt(3,n,m),surf%wt(n,m))
  !---Initialize point and weight matrices
  p=0.d0
  w=0.d0
  !---Select first curve
  cind=cad_tmp%lgws(1,si)
  ctemp=cad_tmp%lgwc(1:2,cind)
  !---Insert first boundary curve into matrices
  if(m==2)then ! Linear curve
    p(:,1,1)=cad_tmp%lgwv(:,ctemp(1))
    p(:,1,2)=cad_tmp%lgwv(:,ctemp(2))
    w(1,1:2)=1.d0
  else ! Quadratic curve
    wind=cad_tmp%lgwc(3,cind)
    p(:,1,1)=cad_tmp%lgwv(:,ctemp(1))
    p(:,1,2)=cad_tmp%lgww(1:3,wind)
    p(:,1,3)=cad_tmp%lgwv(:,ctemp(2))
    w(1,:)=(/1.d0,cad_tmp%lgww(4,wind), 1.d0/)
  endif
  !---Select second curve
  cind=cad_tmp%lgws(2,si)
  ctemp=cad_tmp%lgwc(1:2,cind)
  !---Check relative orientations of first and second curves
  ccheck=cad_tmp%lgwc(1:2,cad_tmp%lgws(1,si))
  ind=2
  !---If endpoint of first curve is not on second curve flip first curve
  flip=(cad_tmp%lgwv(1,ctemp(2)).NE.p(1,1,m)).OR.(cad_tmp%lgwv(2,ctemp(2)).NE.p(2,1,m)).OR.(cad_tmp%lgwv(3,ctemp(2)).NE.p(3,1,m))
  flip=flip.AND.((cad_tmp%lgwv(1,ctemp(1)).NE.p(1,1,m)) .OR. &
  (cad_tmp%lgwv(2,ctemp(1)).NE.p(2,1,m)).OR.(cad_tmp%lgwv(3,ctemp(1)).NE.p(3,1,m)))
  if(flip)then
    tmp=p(:,1,1)
    p(:,1,1)=p(:,1,m)
    p(:,1,m)=tmp
    ccheck(2)=ccheck(1)
    ccheck(1)=cad_tmp%lgwc(2,cad_tmp%lgws(1,si))
  endif
  !---Check curve 2 orientation and flip if necissary
  flip=(cad_tmp%lgwv(1,ctemp(2))==p(1,1,m)).AND.(cad_tmp%lgwv(2,ctemp(2))==p(2,1,m)).AND.(cad_tmp%lgwv(3,ctemp(2))==p(3,1,m))
  if(flip)ind=1
  !---Insert second boundary curve into matrices
  if(n==2)then ! Linear curve
    p(:,n,m)=cad_tmp%lgwv(:,ctemp(ind))
    w(n,m)=1.d0
  else ! Quadratic curve
    wind=cad_tmp%lgwc(3,cind)
    p(:,2,m)=cad_tmp%lgww(1:3,wind)
    p(:,3,m)=cad_tmp%lgwv(:,ctemp(ind))
    w(2:3,m)=(/cad_tmp%lgww(4,wind),1.d0/)
  endif
  !---Select third curve
  cind=cad_tmp%lgws(3,si)
  ctemp=cad_tmp%lgwc(1:2,cind)
  ind=2
  !---Check curve 3 orientation and flip if necissary
  flip=(cad_tmp%lgwv(1,ctemp(2))==p(1,n,m)).AND.(cad_tmp%lgwv(2,ctemp(2))==p(2,n,m)).AND.(cad_tmp%lgwv(3,ctemp(2))==p(3,n,m))
  if(flip)ind=1
  !---Insert third boundary curve into matrices
  if(m==2)then ! Linear curve
    p(:,n,1)=cad_tmp%lgwv(:,ctemp(ind))
    w(n,1)=1.d0
  else ! Quadratic curve
    wind=cad_tmp%lgwc(3,cind)
    p(:,n,2)=cad_tmp%lgww(1:3,wind)
    p(:,n,1)=cad_tmp%lgwv(:,ctemp(ind))
    w(n,1:2)=(/1.d0,cad_tmp%lgww(4,wind)/)
  endif
  !---Select fourth curve
  cind=cad_tmp%lgws(4,si)
  ctemp=cad_tmp%lgwc(1:2,cind)
  !---Insert fourth boundary curve into matrices
  if(n==3)then ! Only if quadratic
    wind=cad_tmp%lgwc(3,cind)
    p(:,2,1)=cad_tmp%lgww(1:3,wind)
    w(2,1)=cad_tmp%lgww(4,wind)
  endif
  !---If dual quadratic surface insert surface point and weight into matrices
  if(n==3.AND.m==3)then
    wind=cad_tmp%lgws(5,si)
    p(:,2,2)=cad_tmp%lgww(1:3,wind)
    w(2,2)=cad_tmp%lgww(4,wind)
  endif
  surf%pt=p(:,1:n,1:m)
  surf%wt=w(1:n,1:m)
endif
DEBUG_STACK_POP
end subroutine mesh_t3d_surfconst
!------------------------------------------------------------------------------
!> Adjust boundary points to CAD boundary.
!------------------------------------------------------------------------------
subroutine mesh_t3d_reffix(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
real(r8) :: pt(3)
integer(i4) :: i,ierr,j,k,ind,ed,ep(2),fp(3),npcors,nerr
integer(i4), pointer :: tmp(:)
integer(i4), allocatable :: emap(:),fmap(:)
class(oft_mesh), pointer :: pmesh,mesh
type(T3D_cadgeom), pointer :: pmesh_cad_rep
type(T3D_cadlink), pointer :: pmesh_cad_link
CHARACTER(LEN=60) :: error_str
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
  allocate(tmp(mesh%base%np))
  tmp=0
  !$omp parallel do
  do i=1,mesh%np
    !tmp(mesh%global%lp(i))=i
    tmp(mesh%base%lp(i))=i
  end do
  !---Get CAD representation aliases
  pmesh_cad_rep=>ML_cad_rep(mg_mesh%level-1)
  pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
  cad_link=>ML_cad_link(mg_mesh%level)
  cad_rep=>ML_cad_rep(mg_mesh%level)
  !---Get edge boundary mapping
  allocate(emap(mesh%ne))
  CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
  !---Copy parent edge linkage to new mesh
  allocate(cad_link%lbeg(mesh%nbe))
  do i=1,pmesh%nbe
    if(associated(pmesh_cad_link%lbeg(i)%wo))then
      ep=tmp(pmesh%le(:,pmesh%lbe(i)))
      if(any(ep==0).OR.any(ep>mesh%np))cycle
      ind=abs(mesh_local_findedge(mesh,(/ep(1),ep(2)/)))
      if(ind==0)cycle
      if(emap(ind)==0)call oft_abort('Bad edge link','mesh_t3d_reffix',__FILE__)
      cad_link%lbeg(emap(ind))%wo=>pmesh_cad_link%lbeg(i)%wo
    endif
  enddo
  !---Get face boundary mapping
  allocate(fmap(mesh%nf))
  CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf)
  !---Copy parent face linkage to new mesh
  allocate(cad_link%lbfg(mesh%nbf))
  do i=1,pmesh%nbf
    fp=tmp(pmesh%lf(:,pmesh%lbf(i)))
    if(any(fp==0).OR.any(fp>mesh%np))cycle
    ind=abs(mesh_local_findface(mesh,(/fp(1),fp(2),fp(3)/)))
    if(ind==0)cycle
    if(fmap(ind)==0)call oft_abort('Bad face link','mesh_t3d_reffix',__FILE__)
    if(associated(pmesh_cad_link%lbfg(i)%wo))cad_link%lbfg(fmap(ind))%wo=>pmesh_cad_link%lbfg(i)%wo
  enddo
  !---Destroy mapping arrays
  deallocate(tmp,emap,fmap)
  if(oft_debug_print(1))write(*,*)'Complete'
  DEBUG_STACK_POP
  return
endif
!---Refine new boundary points using CAD model
if(oft_debug_print(1))write(*,*)'Adjusting points to T3D boundary'
!---Get CAD representation aliases
pmesh_cad_rep=>ML_cad_rep(mg_mesh%level-1)
pmesh_cad_link=>ML_cad_link(mg_mesh%level-1)
cad_rep=>pmesh_cad_rep
cad_link=>ML_cad_link(mg_mesh%level)
!---Get parent edge boundary mapping
allocate(emap(pmesh%ne))
CALL get_inverse_map(pmesh%lbe,pmesh%nbe,emap,pmesh%ne)
!---Locate edge end points and place daughter point
nerr=0
!!$omp parallel do private(j,pt,ierr,obj_curve,obj_surf) reduction(+:nerr)
do i=1,mesh%nbp
  if(mesh%global%gbp(mesh%lbp(i)))then
    if(mesh%lbp(i)<=pmesh%np)cycle
    j=emap(mesh%lbp(i)-pmesh%np)
    if(.NOT.associated(pmesh_cad_link%lbeg(j)%wo))cycle
    pt=mesh%r(:,mesh%lbp(i))
    SELECT TYPE(obj=>pmesh_cad_link%lbeg(j)%wo)
      TYPE IS(cad_curve)
        call cad_curve_midpoint(obj,pt,pmesh%r(:,pmesh%le(1,pmesh%lbe(j))), &
                                pmesh%r(:,pmesh%le(2,pmesh%lbe(j))),.5d0,.5d0,ierr)
        IF(ierr==0)THEN
          mesh%r(:,mesh%lbp(i))=pt
        ELSE
          nerr=nerr+1
        END IF
      TYPE IS(cad_surf)
        call cad_surf_midpoint(obj,pt,pmesh%r(:,pmesh%le(1,pmesh%lbe(j))), &
                               pmesh%r(:,pmesh%le(2,pmesh%lbe(j))),.5d0,.5d0,ierr)
        IF(ierr==0)THEN
          mesh%r(:,mesh%lbp(i))=pt
        ELSE
          nerr=nerr+1
        END IF
    END SELECT
  endif
enddo
deallocate(emap)
!---Synchornize error count
nerr=oft_mpi_sum(nerr)
IF(oft_env%head_proc.AND.nerr>0)THEN
  WRITE(error_str,'(A,I4,A)')'Refinement node placement failed for ',INT(nerr,2),' edges'
  CALL oft_warn(error_str)
END IF
!---Get edge and face boundary mapping
allocate(emap(mesh%ne),fmap(mesh%nf))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne) ! Get edge map
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf) ! Get face map
allocate(cad_link%lbeg(mesh%nbe),cad_link%lbfg(mesh%nbf))
npcors=pmesh%np
!---Copy edge linkage to daughter edges
!$omp parallel do private(i,k,ed)
do j=1,pmesh%nbe
  i=pmesh%lbe(j)
  if(.NOT.associated(pmesh_cad_link%lbeg(j)%wo))cycle
  do k=1,2
    ed=emap(mg_mesh%inter(mg_mesh%level-1)%lede(k,i))
    cad_link%lbeg(ed)%wo=>pmesh_cad_link%lbeg(j)%wo
  end do
enddo
!---Copy face linkage to daughter edges and faces
!$omp parallel do private(i,k,ed)
do j=1,pmesh%nbf
  i=pmesh%lbf(j)
  if(associated(pmesh_cad_link%lbfg(j)%wo))then
    do k=1,3
      ed=emap(mg_mesh%inter(mg_mesh%level-1)%lfde(k,i))
      cad_link%lbeg(ed)%wo=>pmesh_cad_link%lbfg(j)%wo
    end do
  end if
  do k=1,4
    ed=fmap(mg_mesh%inter(mg_mesh%level-1)%lfdf(k,i))
    if(associated(pmesh_cad_link%lbfg(j)%wo))cad_link%lbfg(ed)%wo=>pmesh_cad_link%lbfg(j)%wo
  end do
enddo
!---Destroy mapping arrays
deallocate(emap,fmap)
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine mesh_t3d_reffix
!------------------------------------------------------------------------------
!> Add quadratic mesh node points using CAD model
!------------------------------------------------------------------------------
subroutine mesh_t3d_add_quad(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
real(r8) :: pt(3)
integer(i4) :: i,j,k,ierr,nerr
integer(i4), allocatable :: emap(:)
class(oft_mesh), pointer :: mesh
CHARACTER(LEN=60) :: error_str
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,*)'Setting T3D quadratic nodes'
!---Get CAD representation alias
mesh=>mg_mesh%mesh
cad_rep=>ML_cad_rep(mg_mesh%level)
cad_link=>ML_cad_link(mg_mesh%level)
!---Get edge boundary mapping
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
!---Locate edge end points and place daughter node
nerr=0
!!$omp parallel do private(j,pt,ierr) reduction(+:nerr)
do i=1,mesh%ne
  mesh%ho_info%lep(1,i)=i
  mesh%ho_info%r(:,i)=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
  if(mesh%global%gbe(i))then
    j=emap(i)
    if(.NOT.associated(cad_link%lbeg(j)%wo))cycle
    SELECT TYPE(obj=>cad_link%lbeg(j)%wo)
      TYPE IS(cad_curve)
        pt=mesh%ho_info%r(:,i)
        call cad_curve_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)), &
                                mesh%r(:,mesh%le(2,i)),.5d0,.5d0,ierr)
        IF(ierr==0)THEN
          mesh%ho_info%r(:,i)=pt
          DO k=mesh%kec(i),mesh%kec(i+1)-1
            mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
          END DO
        ELSE
          nerr=nerr+1
        END IF
      TYPE IS(cad_surf)
        pt=mesh%ho_info%r(:,i)
        call cad_surf_midpoint(obj,pt,mesh%r(:,mesh%le(1,i)), &
                               mesh%r(:,mesh%le(2,i)),.5d0,.5d0,ierr)
        IF(ierr==0)THEN
          mesh%ho_info%r(:,i)=pt
          DO k=mesh%kec(i),mesh%kec(i+1)-1
            mesh%ho_info%is_curved(mesh%lec(k))=.TRUE.
          END DO
        ELSE
          nerr=nerr+1
        END IF
    END SELECT
  endif
enddo
!---Destroy mapping array
deallocate(emap)
!---Synchornize error count
nerr=oft_mpi_sum(nerr)
IF(oft_env%head_proc.AND.nerr>0)THEN
  WRITE(error_str,'(A,I4,A)')'Quadratic node placement failed for ',INT(nerr,2),' edges'
  CALL oft_warn(error_str)
END IF
if(oft_debug_print(1))write(*,*)'Complete'
DEBUG_STACK_POP
end subroutine mesh_t3d_add_quad
!------------------------------------------------------------------------------
!> Reflect a T3D mesh and CAD model across a plane
!------------------------------------------------------------------------------
subroutine mesh_t3d_reflect(mesh,k,tol,per_flag)
class(oft_mesh), intent(inout) :: mesh
integer(i4), intent(in) :: k !< Index of plane normal coordinate ( eg. 1 -> y-z plane )
real(r8), intent(in) :: tol !< Tolerance for marking point as on the reflection plane
logical, intent(in) :: per_flag !< Flag for periodicity
integer(i4) :: k1,k2,npold,ncold,i
integer(i4) :: ngme,ngmf,ngwv,ngww,ngwc,ngws,ngwp
integer(i4) :: ngmeold,ngmfold,ngwvold,ngwwold,ngwcold,ngwsold,ngwpold
integer(i4), allocatable :: newindex(:),ltemp(:,:),lpartmp(:),newnmgf(:)
integer(i4), allocatable, dimension(:) :: vtmp,ctmp,stmp,ptmp
real(r8), allocatable :: rtemp(:,:)
real(r8) :: rnorm
character(LEN=1), PARAMETER :: coords(3)=(/'x','y','z'/)
DEBUG_STACK_PUSH
if(abs(k-2)>1)call oft_abort('Invalid coordinate index','mesh_t3d_reflect',__FILE__)
IF(oft_debug_print(1))write(*,'(2X,2A)')'Reflecting mesh -> ',coords(k)
k1=mod(k ,3)+1
k2=mod(k1,3)+1
!---Reflect points that are not on reflection plane
npold=mesh%np
allocate(newindex(2*mesh%np),rtemp(3,2*mesh%np))
rtemp(:,1:mesh%np)=mesh%r
deallocate(mesh%r)
do i=1,npold
  if (abs(rtemp(k,i))<=tol) then
    rtemp(k,i)=0.d0
    newindex(i)=i
  else
    mesh%np=mesh%np+1
    rtemp(k ,mesh%np)=-rtemp(k ,i)
    rtemp(k1,mesh%np)= rtemp(k1,i)
    rtemp(k2,mesh%np)= rtemp(k2,i)
    newindex(i)=mesh%np
  endif
enddo
allocate(mesh%r(3,mesh%np))
mesh%r=rtemp(:,1:mesh%np)
deallocate(rtemp)
!---Reflect CAD model edges that are not on reflection plane
ngme=cad_rep%ngme
ngwv=cad_tmp%ngwv
ngww=cad_tmp%ngww
ngwc=cad_tmp%ngwc
ngws=cad_tmp%ngws
ngwp=cad_tmp%ngwp
if(ngme>0)then
  ngmeold=ngme
  allocate(ltemp(4,2*ngme))
  ltemp(:,1:ngme)=cad_rep%lgme
  deallocate(cad_rep%lgme)
  do i=1,ngmeold
    if (all(mesh%r(k,ltemp(1:2,i))<=0.d0)) then
      ngme=ngme+1
      ltemp(1,ngme)=newindex(ltemp(2,i))
      ltemp(2,ngme)=newindex(ltemp(1,i))
      ltemp(3,ngme)=ltemp(3,i)
      if(ltemp(3,i)==2)then
        ltemp(4,ngme)=ltemp(4,i)+ngwc
      elseif(ltemp(3,i)==3)then
        ltemp(4,ngme)=ltemp(4,i)+ngws
      else
        ltemp(3:4,ngme)=0
      endif
    else
      ngme=ngme+1
      ltemp(1,ngme)=newindex(ltemp(2,i))
      ltemp(2,ngme)=newindex(ltemp(1,i))
      ltemp(3,ngme)=ltemp(3,i)
      if(ltemp(3,i)==2)then
        ltemp(4,ngme)=ltemp(4,i)+ngwc
      elseif(ltemp(3,i)==3)then
        ltemp(4,ngme)=ltemp(4,i)+ngws
      else
        ltemp(3:4,ngme)=0
      endif
    endif
  enddo
  allocate(cad_rep%lgme(4,ngme))
  cad_rep%lgme=ltemp
  cad_rep%ngme=ngme
  deallocate(ltemp)
endif
!---Reflect CAD model faces that are not on reflection plane
ngmf=cad_rep%ngmf
if(ngmf>0)then
  ngmfold=ngmf
  !allocate(ltemp(5,2*ngmf),lpartmp(2*ngmf),newnmgf(ngmf))
  allocate(ltemp(5,2*ngmf),lpartmp(2*ngmf),newnmgf(ngmf))
  lpartmp=0
  ltemp(:,1:ngmf)=cad_rep%lgmf
  lpartmp(1:ngmf)=cad_rep%lgmf_par
  !deallocate(lgmf,lgpf)
  deallocate(cad_rep%lgmf,cad_rep%lgmf_par)
  do i=1,ngmfold
    if (all(mesh%r(k,ltemp(1:3,i))<=tol)) then
      ltemp(4:5,i)=0
      newnmgf(i)=0
    else
      ngmf=ngmf+1
      ltemp(1,ngmf)=newindex(ltemp(1,i))
      ltemp(2,ngmf)=newindex(ltemp(3,i))
      ltemp(3,ngmf)=newindex(ltemp(2,i))
      newnmgf(i)=ngmf
      if(ltemp(4,i)==3)then
        ltemp(4,ngmf)=ltemp(4,i)
        ltemp(5,ngmf)=ltemp(5,i)+ngws
      elseif(ltemp(4,i)==5)then
        ltemp(4,ngmf)=ltemp(4,i)
        ltemp(5,ngmf)=ltemp(5,i)+ngwp
      else
        ltemp(4:5,ngmf)=0
      end if
      IF(lpartmp(i)>0)THEN
        IF(newnmgf(lpartmp(i))>0)lpartmp(ngmf)=newnmgf(lpartmp(i))
      END IF
      if(per_flag)then
        rnorm=abs(mesh%r(k,ltemp(2,ngmf))-mesh%r(k,ltemp(1,ngmf)))
        rnorm=rnorm+abs(mesh%r(k,ltemp(3,ngmf))-mesh%r(k,ltemp(1,ngmf)))
        if(rnorm<=tol)then
          lpartmp(ngmf)=i
        end if
      end if
    endif
  enddo
  !allocate(lgmf(5,ngmf),lgpf(ngmf))
  allocate(cad_rep%lgmf(5,ngmf),cad_rep%lgmf_par(ngmf))
  cad_rep%lgmf=ltemp(:,1:ngmf)
  cad_rep%ngmf=ngmf
  cad_rep%lgmf_par=lpartmp(1:ngmf)
  !lgpf=lpartmp(1:ngmf)
  !deallocate(ltemp,lpartmp,newnmgf)
  deallocate(ltemp,lpartmp,newnmgf)
endif
!---START CAD reflection
!---Reflect CAD vertices
if(ngwv>0)then
  ngwvold=ngwv
  allocate(rtemp(3,2*ngwv),vtmp(2*ngwv)) ! Create temporary array
  rtemp(:,1:ngwv)=cad_tmp%lgwv ! Popluate with current objects
  vtmp(1:ngwv)=cad_tmp%vtmp
  deallocate(cad_tmp%lgwv,cad_tmp%vtmp) ! Kill old arrays
  do i=1,ngwvold ! Reflect and copy to make new objects
    if(abs(rtemp(k,i))<=tol)then ! Object on Ref-Plane
      rtemp(:,i+ngwv)=rtemp(:,i)
      IF(vtmp(i)==0)THEN
        vtmp(i+ngwv)=0
      ELSE
        vtmp(i+ngwv)=vtmp(i)+cad_tmp%nv
      END IF
    else ! Object not on Ref-Plane
      rtemp(k ,i+ngwv)=-rtemp(k ,i)
      rtemp(k1,i+ngwv)= rtemp(k1,i)
      rtemp(k2,i+ngwv)= rtemp(k2,i)
      vtmp(i+ngwv)=vtmp(i)+cad_tmp%nv
    endif
  enddo
  allocate(cad_tmp%lgwv(3,2*ngwv),cad_tmp%vtmp(2*ngwv)) ! Create new object array
  cad_tmp%lgwv=rtemp ! Populate new object array
  cad_tmp%vtmp=vtmp
  deallocate(rtemp,vtmp) ! Kill temporary array
endif
!---Reflect CAD weighting points
if(ngww>0)then
  ngwwold=ngww
  allocate(rtemp(4,2*ngww)) ! Create temporary array
  rtemp(:,1:ngww)=cad_tmp%lgww ! Popluate with current objects
  deallocate(cad_tmp%lgww) ! Kill old array
  do i=1,ngwwold ! Reflect and copy to make new objects
    if(abs(rtemp(k,i))<=tol)then ! Object on Ref-Plane
      rtemp(:,i+ngww)=rtemp(:,i)
    else ! Object not on Ref-Plane
      rtemp(k ,i+ngww)=-rtemp(k ,i)
      rtemp(k1,i+ngww)= rtemp(k1,i)
      rtemp(k2,i+ngww)= rtemp(k2,i)
      rtemp(4,i+ngww)=rtemp(4,i)
    endif
  enddo
  allocate(cad_tmp%lgww(4,2*ngww)) ! Create new object array
  cad_tmp%lgww=rtemp ! Populate new object array
  deallocate(rtemp) ! Kill temporary array
endif
!---Reflect CAD curves
if(ngwc>0)then
  ngwcold=ngwc
  allocate(ltemp(3,2*ngwc),ctmp(2*ngwc)) ! Create temporary array
  ltemp(:,1:ngwc)=cad_tmp%lgwc ! Popluate with current objects
  ctmp(1:ngwc)=cad_tmp%ctmp
  deallocate(cad_tmp%lgwc,cad_tmp%ctmp) ! Kill old array
  do i=1,ngwcold ! Copy to make new objects
    if(ctmp(i)==0)then
      ltemp(1:3,i+ngwc)=0
      ctmp(i+ngwc)=0
    else
      ltemp(1:2,i+ngwc)=ltemp(1:2,i)+ngwv ! Reference copied vertices
      IF(ctmp(i)==0)THEN
        ctmp(i+ngwc)=0
      ELSE
        ctmp(i+ngwc)=ctmp(i)+cad_tmp%nc
      END IF
      if(ltemp(3,i)==0)then ! No weight point
        ltemp(3,i+ngwc)=0
      else ! Reference copied weight point
        ltemp(3,i+ngwc)=ltemp(3,i)+ngww
      end if
    end if
  enddo
  allocate(cad_tmp%lgwc(3,2*ngwc),cad_tmp%ctmp(2*ngwc)) ! Create new object array
  cad_tmp%lgwc=ltemp ! Populate new object array
  cad_tmp%ctmp=ctmp
  deallocate(ltemp,ctmp) ! Kill temporary array
endif
!---Reflect CAD surfaces
if(ngws>0)then
  ngwsold=ngws
  allocate(ltemp(5,2*ngws),stmp(2*ngws)) ! Create temporary array
  ltemp(:,1:ngws)=cad_tmp%lgws ! Popluate with current objects
  stmp(1:ngws)=cad_tmp%stmp
  deallocate(cad_tmp%lgws,cad_tmp%stmp) ! Kill old array
  do i=1,ngwsold ! Copy to make new objects
    if(stmp(i)==0)then
      ltemp(1:5,i+ngws)=0
      stmp(i+ngws)=0
    else
      ltemp(1:4,i+ngws)=ltemp(1:4,i)+ngwc ! Reference copied curves
      IF(stmp(i)==0)THEN
        stmp(i+ngws)=0
      ELSE
        stmp(i+ngws)=stmp(i)+cad_tmp%ns
      END IF
      if(ltemp(5,i)==0)then ! No weight point
        ltemp(5,i+ngws)=0
      else ! Reference copied weight point
        ltemp(5,i+ngws)=ltemp(5,i)+ngww
      end if
    end if
  enddo
  allocate(cad_tmp%lgws(5,2*ngws),cad_tmp%stmp(2*ngws)) ! Create new object array
  cad_tmp%lgws=ltemp ! Populate new object array
  cad_tmp%stmp=stmp
  deallocate(ltemp,stmp) ! Kill temporary array
endif
!---Reflect CAD patches
if(ngwp>0)then
  ngwpold=ngwp
  allocate(ptmp(2*ngwp)) ! Create temporary array
  ptmp(1:ngwp)=cad_tmp%ptmp ! Popluate with current objects
  deallocate(cad_tmp%ptmp) ! Kill old array
  do i=1,ngwpold ! Copy to make new objects
    if(ptmp(i)==0)then
      ptmp(i+ngwp)=0
    else
      ptmp(i+ngwp)=ptmp(i)+cad_tmp%np
    end if
  enddo
  allocate(cad_tmp%ptmp(2*ngwp)) ! Create new object array
  cad_tmp%ptmp=ptmp ! Populate new object array
  deallocate(ptmp) ! Kill temporary array
endif
!---Increment object counters
cad_tmp%ngwv=cad_tmp%ngwv*2
cad_tmp%ngww=cad_tmp%ngww*2
cad_tmp%ngwc=cad_tmp%ngwc*2
cad_tmp%ngws=cad_tmp%ngws*2
cad_tmp%ngwp=cad_tmp%ngwp*2
cad_tmp%nv=cad_tmp%nv*2
cad_tmp%nc=cad_tmp%nc*2
cad_tmp%ns=cad_tmp%ns*2
cad_tmp%np=cad_tmp%np*2
!---END CAD Reflection
!---Reflect cells
ncold=mesh%nc
allocate(ltemp(4,2*mesh%nc))
ltemp(:,1:mesh%nc)=mesh%lc
deallocate(mesh%lc,mesh%reg)
do i=1,ncold
  mesh%nc=mesh%nc+1
  ltemp(1,mesh%nc)=newindex(ltemp(1,i))
  ltemp(2,mesh%nc)=newindex(ltemp(2,i))
  ltemp(3,mesh%nc)=newindex(ltemp(4,i))
  ltemp(4,mesh%nc)=newindex(ltemp(3,i))
enddo
allocate(mesh%lc(4,mesh%nc),mesh%reg(mesh%nc))
mesh%lc=ltemp(:,1:mesh%nc)
mesh%reg=1
deallocate(ltemp,newindex)
DEBUG_STACK_POP
end subroutine mesh_t3d_reflect
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine mesh_t3d_set_periodic(mesh)
class(oft_mesh), intent(inout) :: mesh
integer(i4) :: i,j,jj,k,l,m,n,iper
real(r8) :: pt_i(3),pt_j(3),d_plane,per_dir(3)
real(r8), parameter :: tol=1.d-6
logical :: flag(3)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Setting T3D periodicity'
!---Set periodic faces
mesh%periodic%nper=COUNT(ref_per)
ALLOCATE(mesh%periodic%lp(mesh%np))
ALLOCATE(mesh%periodic%le(mesh%ne))
mesh%periodic%lp=-1
mesh%periodic%le=-1
DO iper=1,3
  IF(.NOT.ref_per(iper))CYCLE
  per_dir=0.d0
  per_dir(iper)=1.d0
  !$omp parallel do private(j,i,l,k,n,m,pt_i,pt_j,d_plane,flag)
  DO jj=1,mesh%nbf
    j=mesh%lbf(jj)
    i=mesh%periodic%lf(j)
    IF(i<=0)CYCLE
    !---Check edges
    flag=.FALSE.
    DO n=1,3
      m=ABS(mesh%lfe(n,j))
      IF(mesh%periodic%le(m)>0)CYCLE
      pt_j=(mesh%r(:,mesh%le(1,m))+mesh%r(:,mesh%le(2,m)))/2.d0
      pt_j=pt_j-DOT_PRODUCT(pt_j,per_dir)*per_dir
      DO l=1,3
        IF(flag(l))CYCLE
        k=ABS(mesh%lfe(l,i))
        pt_i=(mesh%r(:,mesh%le(1,k))+mesh%r(:,mesh%le(2,k)))/2.d0
        pt_i=pt_i-DOT_PRODUCT(pt_i,per_dir)*per_dir
        !---
        pt_i=(pt_i-pt_j)
        d_plane=DOT_PRODUCT(pt_i,pt_i)
        IF(d_plane<tol)THEN
          flag(l)=.TRUE.
          mesh%periodic%le(m)=k
        END IF
      END DO
    END DO
    !---Check points
    flag=.FALSE.
    DO n=1,3
      m=ABS(mesh%lf(n,j))
      IF(mesh%periodic%lp(m)>0)CYCLE
      pt_j=mesh%r(:,m)
      pt_j=pt_j-DOT_PRODUCT(pt_j,per_dir)*per_dir
      DO l=1,3
        IF(flag(l))CYCLE
        k=ABS(mesh%lf(l,i))
        pt_i=mesh%r(:,k)
        pt_i=pt_i-DOT_PRODUCT(pt_i,per_dir)*per_dir
        !---
        pt_i=(pt_i-pt_j)
        d_plane=DOT_PRODUCT(pt_i,pt_i)
        IF(d_plane<tol)THEN
          flag(l)=.TRUE.
          mesh%periodic%lp(m)=k
        END IF
      END DO
    END DO
  END DO
END DO
DEBUG_STACK_POP
end subroutine mesh_t3d_set_periodic
end module oft_mesh_t3d
