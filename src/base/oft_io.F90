!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_io.F90
!
!> HDF5 file manipulation for output and restart data.
!!
!! Functions to write and read HDF5 data files and XDMF files for visit
!! - Parallel restart file creation and read in
!! - Data output for plotting
!! - XDMF reference files for plot files
!!
!! @author Chris Hansen
!! @date Summer 2010
!! @ingroup doxy_oft_base
!------------------------------------------------------------------------------
MODULE oft_io
USE oft_base
USE hdf5
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Binary output file object
!------------------------------------------------------------------------------
TYPE :: oft_bin_file
  LOGICAL :: header_written = .FALSE. !< Was header written during this run
  INTEGER(i4) :: io_unit = 0 !< Output unit for file
  INTEGER(i4) :: ncomm = 0 !< Number of comment lines
  INTEGER(i4) :: nfields = 0 !< Number of fields
  INTEGER(i4) :: nbytes = 0 !< Number of bytes per line
  INTEGER(i4), POINTER, DIMENSION(:) :: field_size => NULL() !< Dimension of each field
  CHARACTER(LEN=OFT_PATH_SLEN) :: filename = '' !< Output filename
  CHARACTER(LEN=OFT_SLEN) :: filedesc = '' !< Description string
  CHARACTER(LEN=2), POINTER, DIMENSION(:) :: field_type => NULL() !< Field types
  CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: field_names => NULL() !< Field names
  CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: field_desc => NULL() !< Field descriptions
  CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: comm_lines => NULL() !< Comment lines
CONTAINS
  !> Create object and allocate storage
  PROCEDURE :: setup => bin_file_setup
  !> Open output file for writing
  PROCEDURE :: open => bin_file_open
  !> Close output file
  PROCEDURE :: close => bin_file_close
  !> Add field to output file
  PROCEDURE :: add_field => bin_file_add
  !> Add comment line to output file
  PROCEDURE :: add_comm => bin_file_add_comm
  !> Write header containing metadata to file
  PROCEDURE :: write_header => bin_file_header
  !>
  PROCEDURE :: write => bin_file_write
  !>
  PROCEDURE :: flush => bin_file_flush
END TYPE oft_bin_file
!------------------------------------------------------------------------------
!> HDF5 restart structure
!!
!! Contains definition and data for saving/loading to distributed restart
!! files.
!!
!! @note Local ownership in HDF5 files is different from linear algebra
!! ownership.
!------------------------------------------------------------------------------
type :: hdf5_rst
  logical :: full = .TRUE. !< Distributed data flag
  integer(i4) :: count = 0 !< Number of values owned by local process
  integer(i8) :: offset = 0 !< Offset in global array
  integer(i8) :: dim = 0 !< Length of global array
  real(r8), pointer, dimension(:) :: data => NULL() !< Array holding local data
end type hdf5_rst
!------------------------------------------------------------------------------
!> Information for XDMF plotting groups in HDF5 plot file
!------------------------------------------------------------------------------
type :: xdmf_plot_file
  integer(i4) :: n_ts = 0
  integer(i4) :: curr_ts = 0
  integer(i4) :: n_grids = 0
  character(LEN=OFT_PATH_SLEN) :: file_path = ''
  character(LEN=OFT_PATH_SLEN) :: group_name = ''
  character(LEN=OFT_PATH_SLEN) :: grid_names(10) = ''
CONTAINS
  PROCEDURE :: setup => xmdf_setup
  PROCEDURE :: add_mesh => xdmf_add_mesh
  PROCEDURE :: add_timestep => xdmf_add_timestep
  PROCEDURE :: clear_timesteps => xdmf_clear_timesteps
  PROCEDURE :: write_scalar => xdmf_write_scalar
  PROCEDURE :: write_vector => xdmf_write_vector
  GENERIC :: write => write_scalar, write_vector
end type xdmf_plot_file
!------------------------------------------------------------------------------
!> Write data to an HDF5 file
!------------------------------------------------------------------------------
INTERFACE hdf5_write
  MODULE PROCEDURE hdf5_write_scalar_r8
  MODULE PROCEDURE hdf5_write_scalar_i4
  MODULE PROCEDURE hdf5_write_1d_r8
  MODULE PROCEDURE hdf5_write_1d_i4
  MODULE PROCEDURE hdf5_write_2d_r8
  MODULE PROCEDURE hdf5_write_2d_i4
  MODULE PROCEDURE hdf5_write_rst
END INTERFACE hdf5_write
!------------------------------------------------------------------------------
!> Read data from an HDF5 file
!------------------------------------------------------------------------------
INTERFACE hdf5_read
  MODULE PROCEDURE hdf5_read_scalar_r8
  MODULE PROCEDURE hdf5_read_scalar_i4
  MODULE PROCEDURE hdf5_read_1d_r8
  MODULE PROCEDURE hdf5_read_1d_i4
  MODULE PROCEDURE hdf5_read_2d_r8
  MODULE PROCEDURE hdf5_read_2d_i4
  MODULE PROCEDURE hdf5_read_rst
END INTERFACE hdf5_read
!---Global variables
integer(i4) :: hdf5_fcount=0 !< Number of HDF5 fields
integer(i4) :: hdf5_ts=0 !< Number of HDF5 time steps
integer(i4), parameter :: hdf5_nl=800 !< Maximum number of HDF5 fields
integer(i4), parameter :: hdf5_flen=40 !< Maximum size of HDF5 filenames
contains
!---------------------------------------------------------------------------------
!> Setup Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_setup(self,filename,desc)
CLASS(oft_bin_file), INTENT(inout) :: self
CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: desc
IF(LEN(filename)>OFT_PATH_SLEN)CALL oft_abort("Filename too long","bin_file_setup",__FILE__)
self%filename = filename
self%nfields = 0
self%ncomm = 0
self%nbytes = 0
IF(PRESENT(desc))THEN
  IF(LEN(desc)>LEN(self%filedesc))CALL oft_abort("Description too long","bin_file_setup",__FILE__)
  self%filedesc=desc
END IF
END SUBROUTINE bin_file_setup
!---------------------------------------------------------------------------------
!> Open Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_open(self)
CLASS(oft_bin_file), INTENT(inout) :: self
IF(self%header_written)THEN
  OPEN(NEWUNIT=self%io_unit,FILE=TRIM(self%filename),POSITION="APPEND",STATUS="OLD")
  WRITE(self%io_unit,*)
  WRITE(self%io_unit,'(A)')"--- BEGIN DATA ---"
  CLOSE(self%io_unit)
  OPEN(NEWUNIT=self%io_unit,FILE=TRIM(self%filename),FORM='UNFORMATTED', &
    POSITION="APPEND",STATUS="OLD",ACCESS="STREAM")
ELSE
  OPEN(NEWUNIT=self%io_unit,FILE=TRIM(self%filename),FORM='UNFORMATTED', &
    POSITION="APPEND",STATUS="OLD",ACCESS="STREAM")
END IF
self%header_written=.FALSE.
END SUBROUTINE bin_file_open
!---------------------------------------------------------------------------------
!> Close Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_close(self)
CLASS(oft_bin_file), INTENT(inout) :: self
CLOSE(self%io_unit)
END SUBROUTINE bin_file_close
!---------------------------------------------------------------------------------
!> Add field to Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_add(self,fieldname,type_str,desc,fsize)
CLASS(oft_bin_file), INTENT(inout) :: self !< File object
CHARACTER(LEN=*), INTENT(in) :: fieldname !< Field name to add
CHARACTER(LEN=2), INTENT(in) :: type_str !< Type of field
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: desc !< Description of field
INTEGER(i4), OPTIONAL, INTENT(in) :: fsize !< Size of field
INTEGER(i4) :: fsize_tmp
INTEGER(i4), POINTER, DIMENSION(:) :: sizes_tmp
CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: names_tmp
CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: desc_tmp
CHARACTER(LEN=2), POINTER, DIMENSION(:) :: types_tmp
IF((type_str(1:1)/='i').AND.(type_str(1:1)/='r'))THEN
  CALL oft_abort("Invalid field type", "bin_file_add", __FILE__)
END IF
!
fsize_tmp = 1
IF(PRESENT(fsize))fsize_tmp = fsize
IF(type_str(2:2)=='4')THEN
  self%nbytes = self%nbytes + 4*fsize_tmp
ELSE IF(type_str(2:2)=='8')THEN
  self%nbytes = self%nbytes + 8*fsize_tmp
ELSE
  CALL oft_abort("Invalid field size", "bin_file_add", __FILE__)
END IF
IF(LEN(fieldname)>OFT_SLEN)CALL oft_abort('Name too long', &
  "bin_file_add",__FILE__)
IF(PRESENT(desc))THEN
  IF(LEN(desc)>OFT_SLEN)CALL oft_abort('Description too long', &
    "bin_file_add",__FILE__)
END IF
!
IF(self%nfields>0)THEN
  sizes_tmp=>self%field_size
  names_tmp=>self%field_names
  types_tmp=>self%field_type
  desc_tmp=>self%field_desc
  ALLOCATE(self%field_size(self%nfields+1))
  ALLOCATE(self%field_names(self%nfields+1))
  ALLOCATE(self%field_type(self%nfields+1))
  ALLOCATE(self%field_desc(self%nfields+1))
  self%field_size(1:self%nfields)=sizes_tmp
  self%field_names(1:self%nfields)=names_tmp
  self%field_type(1:self%nfields)=types_tmp
  self%field_desc(1:self%nfields)=desc_tmp
  DEALLOCATE(sizes_tmp,names_tmp,types_tmp,desc_tmp)
ELSE
  ALLOCATE(self%field_size(self%nfields+1))
  ALLOCATE(self%field_names(self%nfields+1))
  ALLOCATE(self%field_type(self%nfields+1))
  ALLOCATE(self%field_desc(self%nfields+1))
END IF
self%nfields=self%nfields+1
self%field_names(self%nfields)=fieldname
self%field_type(self%nfields)=type_str
self%field_size(self%nfields)=fsize_tmp
IF(PRESENT(desc))THEN
  self%field_desc(self%nfields)=desc
ELSE
  self%field_desc(self%nfields)="No description"
END IF
END SUBROUTINE bin_file_add
!---------------------------------------------------------------------------------
!> Add comment to Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_add_comm(self,comment)
CLASS(oft_bin_file), INTENT(inout) :: self !< File object
CHARACTER(LEN=*), INTENT(in) :: comment !< Comment to add
CHARACTER(LEN=OFT_SLEN), POINTER, DIMENSION(:) :: comm_tmp
IF(LEN(comment)>OFT_SLEN)CALL oft_abort('Comment too long', &
  "bin_file_add_comm",__FILE__)
IF(self%ncomm>0)THEN
  comm_tmp=>self%comm_lines
  ALLOCATE(self%comm_lines(self%ncomm+1))
  self%comm_lines(1:self%ncomm)=comm_tmp
  DEALLOCATE(comm_tmp)
ELSE
  ALLOCATE(self%comm_lines(1))
END IF
self%ncomm=self%ncomm+1
self%comm_lines(self%ncomm)=comment
END SUBROUTINE bin_file_add_comm
!---------------------------------------------------------------------------------
!> Write header for Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_header(self)
CLASS(oft_bin_file), INTENT(inout) :: self
INTEGER(i4) :: i
INTEGER(i8) :: value(8)
!---Create file with binary format description
CALL DATE_AND_TIME(VALUES=value)
OPEN(NEWUNIT=self%io_unit,FILE=TRIM(self%filename))
100 FORMAT("# Created: ",I2,':',I2.2,':'I2.2,' on ',I2,"-",I2.2,"-",I4)
WRITE(self%io_unit,'(A)')"# Open FUSION Toolkit binary output"
WRITE(self%io_unit,'(2A)')"# Description: ",TRIM(self%filedesc)
WRITE(self%io_unit,100)value(5:7),value(2:3),value(1)
IF(self%ncomm>0)THEN
  WRITE(self%io_unit,'(A)')"# "
  WRITE(self%io_unit,'(A)')"# Comments:"
  DO i=1,self%ncomm
    WRITE(self%io_unit,'(2A)')"#  ",self%comm_lines
  END DO
END IF
WRITE(self%io_unit,*)
WRITE(self%io_unit,'(A,I6)')"nfields: ",self%nfields
!---Write out field names
WRITE(self%io_unit,'(A)',ADVANCE="NO")"fields: "
DO i=1,self%nfields
  WRITE(self%io_unit,'(2A)',ADVANCE="NO")TRIM(self%field_names(i))," "
END DO
WRITE(self%io_unit,*)
!---Write out field types
WRITE(self%io_unit,'(A)',ADVANCE="NO")"field_types: "
DO i=1,self%nfields
  WRITE(self%io_unit,'(2A)',ADVANCE="NO")TRIM(self%field_type(i))," "
END DO
WRITE(self%io_unit,*)
!---Write out field sizes
WRITE(self%io_unit,'(A)',ADVANCE="NO")"field_sizes: "
DO i=1,self%nfields
  WRITE(self%io_unit,'(I6,A)',ADVANCE="NO")self%field_size(i)," "
END DO
WRITE(self%io_unit,*)
!---Write descriptions
WRITE(self%io_unit,'(A)')"descriptions:"
DO i=1,self%nfields
  WRITE(self%io_unit,'(A,A,A)',ADVANCE="NO")"  - ",TRIM(self%field_names(i)),":"
  WRITE(self%io_unit,'(A,A)',ADVANCE="NO")" ",TRIM(self%field_desc(i))
  WRITE(self%io_unit,*)
END DO
WRITE(self%io_unit,*)
CLOSE(self%io_unit)
self%header_written=.TRUE.
END SUBROUTINE bin_file_header
!---------------------------------------------------------------------------------
!> Write single set of data to Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_write(self,data_i4,data_i8,data_r4,data_r8)
CLASS(oft_bin_file), INTENT(inout) :: self !< File object
INTEGER(i4), OPTIONAL, INTENT(in) :: data_i4(:) !< integer(i4) data
INTEGER(i8), OPTIONAL, INTENT(in) :: data_i8(:) !< integer(i8) data
REAL(r4), OPTIONAL, INTENT(in) :: data_r4(:) !< real(r4) data
REAL(r8), OPTIONAL, INTENT(in) :: data_r8(:) !< real(r8) data
WRITE(self%io_unit)self%nbytes
IF(PRESENT(data_i4))WRITE(self%io_unit)data_i4
IF(PRESENT(data_i8))WRITE(self%io_unit)data_i8
IF(PRESENT(data_r4))WRITE(self%io_unit)data_r4
IF(PRESENT(data_r8))WRITE(self%io_unit)data_r8
WRITE(self%io_unit)self%nfields
END SUBROUTINE bin_file_write
!---------------------------------------------------------------------------------
!> Flush I/O buffer for Open FUSION Toolkit binary I/O file
!---------------------------------------------------------------------------------
SUBROUTINE bin_file_flush(self)
CLASS(oft_bin_file), INTENT(inout) :: self
FLUSH(self%io_unit)
END SUBROUTINE bin_file_flush
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine xmdf_setup(self,group_name,basepath,persistent_space_tracking)
class(xdmf_plot_file), intent(inout) :: self
CHARACTER(LEN=*), intent(in) :: group_name !< Path to mesh in HDF5 file
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
LOGICAL, OPTIONAL, INTENT(in) :: persistent_space_tracking
integer :: error
IF(PRESENT(basepath))THEN
  call execute_command_line('mkdir -p '//TRIM(basepath), exitstat=error)
  IF(error/=0)CALL oft_abort('Failed to create output directory: '//TRIM(basepath), &
    "xmdf_setup", __FILE__)
  self%file_path=TRIM(basepath)//"oft_xdmf."//hdf5_proc_str()//".h5"
ELSE
  self%file_path="oft_xdmf."//hdf5_proc_str()//".h5"
END IF
self%group_name=TRIM(group_name)
CALL string_to_lower(self%group_name)
CALL hdf5_create_file(TRIM(self%file_path),persistent_space_tracking)
CALL hdf5_create_group(TRIM(self%file_path),TRIM(self%group_name))
end subroutine xmdf_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine xdmf_add_mesh(self,mesh_type,pt_list,cell_list,grid_name)
class(xdmf_plot_file), intent(inout) :: self
integer(i4), intent(in) :: mesh_type !< Mesh type flag (Tet/Tri or Hex/Quad)
real(r8), intent(in) :: pt_list(:,:) !< Point list [3,np]
integer(i4), intent(in) :: cell_list(:,:) !< Cell list [:,nc]
CHARACTER(LEN=*), intent(in) :: grid_name !< Path to mesh in HDF5 file
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix,filename
CHARACTER(LEN=200) :: hdf5_path
integer(i4) :: i,ntrans(4),ierr,io_unit
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
type(mpi_status) :: mpi_stat
#else
integer(i4) :: mpi_stat(MPI_STATUS_SIZE)
#endif
#endif
DEBUG_STACK_PUSH
self%n_grids=self%n_grids+1
self%grid_names(self%n_grids)=TRIM(grid_name)
CALL string_to_lower(self%grid_names(self%n_grids))
IF(.NOT.oft_file_exist(TRIM(self%file_path)))CALL oft_abort("File does not exist", &
  "xdmf_add_mesh",__FILE__)
hdf5_path=TRIM(self%group_name)//"/"//TRIM(self%grid_names(self%n_grids))
CALL hdf5_create_group(TRIM(self%file_path),TRIM(hdf5_path))
CALL hdf5_write(mesh_type,TRIM(self%file_path),TRIM(hdf5_path)//"/TYPE")
CALL hdf5_write(pt_list,TRIM(self%file_path),TRIM(hdf5_path)//"/R",single_prec=.TRUE.)
CALL hdf5_write(cell_list,TRIM(self%file_path),TRIM(hdf5_path)//"/LC")
#ifdef HAVE_MPI
CALL hdf5_write(oft_env%nprocs,TRIM(self%file_path),TRIM(hdf5_path)//"/NBLOCKS")
#endif
!---Create static storage
CALL hdf5_create_group(TRIM(self%file_path),TRIM(hdf5_path)//"/"//hdf5_ts_str(0))
DEBUG_STACK_POP
end subroutine xdmf_add_mesh
!------------------------------------------------------------------------------
!> Adds a timestep to the dump metadata file.
!!
!! Subsequent output will be added to this timestep until another call
!! to this subroutine
!------------------------------------------------------------------------------
subroutine xdmf_add_timestep(self,t)
class(xdmf_plot_file), intent(inout) :: self
real(r8), intent(in) :: t !< Time value
integer(i4) :: i
CHARACTER(LEN=200) :: hdf5_path
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A,ES11.4)')oft_indent,'Creating plot time: ',t
IF(.NOT.oft_file_exist(TRIM(self%file_path)))CALL oft_abort("File does not exist", &
  "xdmf_add_timestep",__FILE__)
self%n_ts=self%n_ts+1
DO i=1,self%n_grids
  hdf5_path=TRIM(self%group_name)//"/"//TRIM(self%grid_names(i))//"/"//TRIM(hdf5_ts_str(self%n_ts))
  CALL hdf5_create_group(TRIM(self%file_path),TRIM(hdf5_path))
  CALL hdf5_write(t,TRIM(self%file_path),TRIM(hdf5_path)//"/TIME")
END DO
DEBUG_STACK_POP
end subroutine xdmf_add_timestep
!------------------------------------------------------------------------------
!> Clear existing timesteps and reset to static fields
!------------------------------------------------------------------------------
subroutine xdmf_clear_timesteps(self,clear_static)
class(xdmf_plot_file), intent(inout) :: self
logical, optional, intent(in) :: clear_static !< Clear static fields as well?
integer(i4) :: j,i,istart
CHARACTER(LEN=200) :: hdf5_path
integer(HSSIZE_T) :: file_sizes(2)
DEBUG_STACK_PUSH
istart=1
IF(PRESENT(clear_static))THEN
  IF(clear_static)istart=0
END IF
if(oft_debug_print(1))write(*,'(2A,I6,A)')oft_indent,'Clearing ',self%n_ts+(1-istart),' existing timesteps'
IF(.NOT.oft_file_exist(TRIM(self%file_path)))CALL oft_abort("File does not exist", &
  "xdmf_clear_timesteps",__FILE__)
DO i=istart,self%n_ts
  DO j=1,self%n_grids
    hdf5_path=TRIM(self%group_name)//"/"//TRIM(self%grid_names(j))//"/"//TRIM(hdf5_ts_str(i))
    CALL hdf5_delete_obj(TRIM(self%file_path),TRIM(hdf5_path))
  END DO
END DO
self%n_ts=0
IF(istart==0)THEN
  hdf5_path=TRIM(self%group_name)//"/"//TRIM(self%grid_names(j))
  CALL hdf5_create_group(TRIM(self%file_path),TRIM(hdf5_path)//"/"//hdf5_ts_str(0))
END IF
if(oft_debug_print(2))THEN
  file_sizes = hdf5_file_size(TRIM(self%file_path))
  write(*,'(2A,ES14.6,A,ES14.6,A)')oft_indent,'Freed ',REAL(file_sizes(2),8)/1.d6,' of ',REAL(file_sizes(1),8)/1.d6,' [MB]'
end if
DEBUG_STACK_POP
end subroutine xdmf_clear_timesteps
!------------------------------------------------------------------------------
!> Write scalar field to plot file
!------------------------------------------------------------------------------
subroutine xdmf_write_scalar(self,data,grid_name,path,centering,single_prec)
class(xdmf_plot_file), intent(in) :: self
real(r8), intent(in) :: data(:) !< Scalar data
CHARACTER(LEN=*), intent(in) :: grid_name !< Grid name
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4), intent(in) :: centering !< Centering of data (1-> vertex; 2-> cell)
logical, optional, intent(in) :: single_prec !< Save as single precision?
CHARACTER(LEN=OFT_PATH_SLEN) :: hdf5_path,grid_lower
CHARACTER(LEN=80) :: attr_data
IF(.NOT.oft_file_exist(TRIM(self%file_path)))CALL oft_abort("File does not exist", &
  "xdmf_write_scalar",__FILE__)
grid_lower = TRIM(grid_name)
CALL string_to_lower(grid_lower)
hdf5_path=TRIM(self%group_name)//"/"//TRIM(grid_lower)//"/"//TRIM(hdf5_ts_str(self%n_ts))
IF(.NOT.hdf5_field_exist(TRIM(self%file_path),TRIM(hdf5_path)))CALL oft_abort("Timestep does not exist", &
  "xdmf_write_scalar",__FILE__)
hdf5_path=TRIM(hdf5_path)//"/"//TRIM(path)
CALL hdf5_write(data,TRIM(self%file_path),TRIM(hdf5_path),single_prec)
attr_data='Scalar'
CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'TYPE',[attr_data])
SELECT CASE(centering)
  CASE(1)
    attr_data='Node'
    CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'CENTERING',[attr_data])
  CASE(2)
    attr_data='Cell'
    CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'CENTERING',[attr_data])
  CASE DEFAULT
    CALL oft_abort('Unknown field centering','xdmf_write_scalar',__FILE__)
END SELECT
end subroutine xdmf_write_scalar
!------------------------------------------------------------------------------
!> Write vector field to plot file
!------------------------------------------------------------------------------
subroutine xdmf_write_vector(self,data,grid_name,path,centering,single_prec)
class(xdmf_plot_file), intent(in) :: self
real(r8), intent(in) :: data(:,:) !< Vector data
CHARACTER(LEN=*), intent(in) :: grid_name !< Grid name
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4), intent(in) :: centering !< Centering of data (1-> vertex; 2-> cell)
logical, optional, intent(in) :: single_prec !< Save as single precision?
CHARACTER(LEN=OFT_PATH_SLEN) :: hdf5_path,grid_lower
CHARACTER(LEN=80) :: attr_data
IF(.NOT.oft_file_exist(TRIM(self%file_path)))CALL oft_abort("File does not exist", &
  "xdmf_write_vector",__FILE__)
grid_lower = TRIM(grid_name)
CALL string_to_lower(grid_lower)
hdf5_path=TRIM(self%group_name)//"/"//TRIM(grid_lower)//"/"//TRIM(hdf5_ts_str(self%n_ts))
IF(.NOT.hdf5_field_exist(TRIM(self%file_path),TRIM(hdf5_path)))CALL oft_abort("Timestep does not exist", &
  "xdmf_write_vector",__FILE__)
hdf5_path=TRIM(hdf5_path)//"/"//TRIM(path)
CALL hdf5_write(data,TRIM(self%file_path),TRIM(hdf5_path),single_prec)
attr_data='Vector'
CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'TYPE',[attr_data])
SELECT CASE(centering)
  CASE(1)
    attr_data='Node'
    CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'CENTERING',[attr_data])
  CASE(2)
    attr_data='Cell'
    CALL hdf5_add_string_attribute(TRIM(self%file_path),TRIM(hdf5_path),'CENTERING',[attr_data])
  CASE DEFAULT
    CALL oft_abort('Unknown field centering','xdmf_write_vector',__FILE__)
END SELECT
end subroutine xdmf_write_vector
!------------------------------------------------------------------------------
!> Get processor rank as string for HDF5 I/O
!------------------------------------------------------------------------------
function hdf5_proc_str(proc_ind) result(proc_str)
integer(i4), optional, intent(in) :: proc_ind
character(LEN=OFT_MPI_PLEN) :: proc_str
100 FORMAT (I OFT_MPI_PLEN.OFT_MPI_PLEN)
IF(PRESENT(proc_ind))THEN
  write(proc_str,100)proc_ind+1
ELSE
  write(proc_str,100)oft_env%rank+1
END IF
end function hdf5_proc_str
!------------------------------------------------------------------------------
!> Get timestep index as string for HDF5 I/O
!------------------------------------------------------------------------------
function hdf5_ts_str(ts_in) result(ts_str)
integer(i4), optional, intent(in) :: ts_in
character(LEN=OFT_MPI_PLEN) :: ts_str
100 FORMAT (I OFT_MPI_PLEN.OFT_MPI_PLEN)
IF(PRESENT(ts_in))THEN
  write(ts_str,100)ts_in
ELSE
  write(ts_str,100)hdf5_ts
END IF
end function hdf5_ts_str
!------------------------------------------------------------------------------
!> Test for exitence of a file
!!
!! @result Logical flag indicating existence
!------------------------------------------------------------------------------
function oft_file_exist(filepath) result(exists)
character(LEN=*), intent(in) :: filepath !< Path to file
logical :: exists
INQUIRE(FILE=TRIM(filepath),EXIST=exists)
end function oft_file_exist
!------------------------------------------------------------------------------
!> Test for exitence of a field in a HDF5 file
!!
!! @result Logical flag indicating existence of field and file
!------------------------------------------------------------------------------
function hdf5_file_size(filepath) result(sizes)
character(LEN=*), intent(in) :: filepath !< Path to file
integer :: error
logical :: exists
integer(HID_T) :: file_id
integer(HSSIZE_T) :: sizes(2)
exists=oft_file_exist(filepath)
sizes=-1
IF(.NOT.exists)RETURN
DEBUG_STACK_PUSH
!---Try to open file as HDF5 file
call h5open_f(error)
call h5fopen_f(TRIM(filepath), H5F_ACC_RDONLY_F, file_id, error)
IF(error/=0)THEN
  call h5close_f(error)
  DEBUG_STACK_POP
  RETURN
END IF
CALL H5Fget_filesize_f(file_id,sizes(1),error)
CALL H5Fget_freespace_f(file_id,sizes(2),error)
!---Close file and finalize HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end function hdf5_file_size
!------------------------------------------------------------------------------
!> Test for exitence of a field in a HDF5 file
!!
!! @result Logical flag indicating existence of field and file
!------------------------------------------------------------------------------
function hdf5_field_exist(filepath,path) result(exists)
character(LEN=*), intent(in) :: filepath !< Path to file
character(LEN=*), intent(in) :: path !< Path of field in file
integer :: access_flag,error,subpath
integer(HID_T) :: file_id,dset_id
logical :: exists
exists=oft_file_exist(filepath)
IF(.NOT.exists)RETURN
DEBUG_STACK_PUSH
!---Try to open file as HDF5 file
access_flag=H5F_ACC_RDONLY_F
call h5open_f(error)
call h5fopen_f(TRIM(filepath), access_flag, file_id, error)
IF(error/=0)exists=.FALSE.
!---Check for desired field
DO subpath=1,LEN_TRIM(path)
  IF(.NOT.exists)EXIT
  IF(subpath<LEN_TRIM(path).AND.path(subpath:subpath)/="/")CYCLE
  CALL h5lexists_f(file_id, "/"//path(1:subpath), exists, error)
END DO
!---Close file and finalize HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end function hdf5_field_exist
!------------------------------------------------------------------------------
!> Test for exitence of a field in a HDF5 file
!!
!! @result Logical flag indicating existence of field and file
!------------------------------------------------------------------------------
subroutine hdf5_field_get_sizes(filepath,path,ndims,dim_sizes)
character(LEN=*), intent(in) :: filepath !< Path to file
character(LEN=*), intent(in) :: path !< Path of field in file
integer(4), intent(out)  :: ndims !< Number of dimensions in field
integer(4), allocatable, dimension(:), intent(out) :: dim_sizes !< Size of each dimension
integer(4) :: access_flag,error
integer(hsize_t), allocatable, dimension(:) :: tmp_sizes,maxdims
integer(HID_T) :: file_id,dset_id,dspace_id
ndims=-1
IF(.NOT.hdf5_field_exist(filepath,path))RETURN
DEBUG_STACK_PUSH
!---Try to open file as HDF5 file
access_flag=H5F_ACC_RDONLY_F
call h5open_f(error)
call h5fopen_f(TRIM(filepath), access_flag, file_id, error)
IF(error==0)THEN
  !---Check for desired field
  CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
  IF(error==0)THEN
    CALL h5dget_space_f(dset_id, dspace_id, error)
    CALL h5sget_simple_extent_ndims_f(dspace_id, ndims, error)
    ALLOCATE(dim_sizes(ndims),tmp_sizes(ndims),maxdims(ndims))
    CALL h5sget_simple_extent_dims_f(dspace_id, tmp_sizes, maxdims, error)
    dim_sizes=INT(tmp_sizes,4)
    DEALLOCATE(tmp_sizes,maxdims)
    !---Close dataspace/set
    call h5sclose_f(dspace_id, error)
    call h5dclose_f(dset_id, error)
  END IF
END IF
!---Close file and finalize HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_field_get_sizes
!------------------------------------------------------------------------------
!> Create an empty HDF5 output file
!------------------------------------------------------------------------------
subroutine hdf5_create_file(filename,persistent_space_tracking)
character(LEN=*), intent(in) :: filename !< Name of file to be created
logical, optional, intent(in) :: persistent_space_tracking
integer :: error
integer(HID_T) :: file_id,plist_id
integer(HSIZE_T) :: zero = 0
logical :: track_free
DEBUG_STACK_PUSH
!---Get processor rank for file creation
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Creating HDF5 output file: ',TRIM(filename)
!---Initialize HDF5 system
call h5open_f(error)
!---Create HDF5 file
track_free=.FALSE.
IF(PRESENT(persistent_space_tracking))track_free=persistent_space_tracking
IF(track_free)THEN
  CALL H5Pcreate_f(H5P_FILE_CREATE_F,plist_id,error)
  CALL H5Pset_file_space_strategy_f(plist_id,H5F_FSPACE_STRATEGY_FSM_AGGR_F,.TRUE.,zero,error)
  CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error, creation_prp=plist_id)
  CALL h5pclose_f(plist_id, error)
ELSE
  call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error)
END IF
call h5fclose_f(file_id, error)
!---Finalize HDF5 system
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_create_file
!------------------------------------------------------------------------------
!> Create HDF5 group in existing file
!------------------------------------------------------------------------------
subroutine hdf5_create_group(filename,group_name)
character(LEN=*), intent(in) :: filename !< Name of HDF5 file
character(LEN=*), intent(in) :: group_name !< Group path
integer :: error
integer(HID_T) :: file_id,grp_id
DEBUG_STACK_PUSH
!---Initialize HDF5 and open vector dump file
call h5open_f(error)
CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_create_group',__FILE__)
CALL h5gcreate_f(file_id, "/"//TRIM(group_name), grp_id, error)
CALL h5gclose_f(grp_id, error)
!---Close vector dump file and finalize HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_create_group
!------------------------------------------------------------------------------
!> Delete HDF5 object in existing file
!------------------------------------------------------------------------------
subroutine hdf5_delete_obj(filename,obj_path)
character(LEN=*), intent(in) :: filename !< Name of HDF5 file
character(LEN=*), intent(in) :: obj_path !< Object path
integer :: error
integer(HID_T) :: file_id,grp_id
DEBUG_STACK_PUSH
!---Initialize HDF5 and open vector dump file
call h5open_f(error)
CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_delete_obj',__FILE__)
CALL h5ldelete_f(file_id, "/"//TRIM(obj_path), error)
IF(error/=0)CALL oft_abort('Error deleting object','hdf5_delete_obj',__FILE__)
!---Close vector dump file and finalize HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_delete_obj
!------------------------------------------------------------------------------
!> Add string attribute to existing object (group or dataset)
!------------------------------------------------------------------------------
subroutine hdf5_add_string_attribute(filename,objname,aname,attr_data)
character(LEN=*), intent(in) :: filename !< Name of HDF5 file
character(LEN=*), intent(in) :: objname !< Name of object (dataset or group)
character(LEN=*), intent(in) :: aname !< Attribute name
character(LEN=80), dimension(:), intent(in) :: attr_data !< Attribute data (80-character lines)
integer :: arank,error
integer(HID_T) :: file_id,obj_id,aspace_id,atype_id,attr_id
INTEGER(SIZE_T) :: attrlen
INTEGER(HSIZE_T), DIMENSION(1) :: adims
DEBUG_STACK_PUSH
!---Initialize HDF5 and open vector dump file
call h5open_f(error)
CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)THEN
  call h5close_f(error)
  CALL oft_abort('Error opening file','hdf5_add_string_attribute',__FILE__)
END IF
!---Open an existing object
CALL h5oopen_f(file_id, objname, obj_id, error)
IF(error/=0)THEN
  call h5fclose_f(file_id, error)
  call h5close_f(error)
  CALL oft_abort('Error opening object','hdf5_add_string_attribute',__FILE__)
END IF
!---Create scalar data space for the attribute
arank=1
adims = SIZE(attr_data)
CALL h5screate_simple_f(arank, adims, aspace_id, error)
!---Create datatype for the attribute
CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
attrlen = 80
CALL h5tset_size_f(atype_id, attrlen, error)
!---Create attribute and write data
CALL h5acreate_f(obj_id, aname, atype_id, aspace_id, attr_id, error)
CALL h5awrite_f(attr_id, atype_id, attr_data, adims, error)
!---Close HDF5 objects
CALL h5aclose_f(attr_id, error)
CALL h5tclose_f(atype_id, error)
CALL h5sclose_f(aspace_id, error)
CALL h5oclose_f(obj_id, error)
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_add_string_attribute
!------------------------------------------------------------------------------
!> real(r8) scalar implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_scalar_r8(val,filename,path,single_prec)
real(r8), intent(in) :: val !< Value to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(in) :: single_prec !< Save as single precision?
real(r8) :: tmpval(1)
tmpval(1)=val
CALL hdf5_write_1d_r8(tmpval,filename,path,single_prec)
end subroutine hdf5_write_scalar_r8
!------------------------------------------------------------------------------
!> integer(i4) scalar implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_scalar_i4(val,filename,path)
integer(i4), intent(in) :: val !< Value to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
integer(i4) :: tmpval(1)
DEBUG_STACK_PUSH
tmpval(1)=val
CALL hdf5_write_1d_i4(tmpval,filename,path)
DEBUG_STACK_POP
end subroutine hdf5_write_scalar_i4
!------------------------------------------------------------------------------
!> real(r8) 1D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_1d_r8(array,filename,path,single_prec)
real(r8), intent(in) :: array(:) !< Values to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(in) :: single_prec !< Save as single precision?
logical :: write_single
integer(i4) :: error
integer(i4), parameter :: one=1
integer(HID_T) :: file_id,dspace_id,dset_id
INTEGER(HSIZE_T), DIMENSION(1) :: dims
DEBUG_STACK_PUSH
!---Remove field if it already exists
IF(hdf5_field_exist(filename,path))THEN
  CALL oft_warn('Overwriting existing object "'//TRIM(path)//'" in file '//TRIM(filename))
  CALL hdf5_delete_obj(filename,path)
END IF
!---Write at lower precision?
write_single=.FALSE.
IF(PRESENT(single_prec))write_single=single_prec
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_write_1d_r8',__FILE__)
!---Create dataset and perform write
dims=SHAPE(array)
call h5screate_simple_f(one,dims,dspace_id,error)
IF(write_single)THEN
  call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_REAL, &
    dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, REAL(array,4), dims, error)
ELSE
  call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_DOUBLE, &
    dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, error)
END IF
!---Close dataset
call h5dclose_f(dset_id, error)
call h5sclose_f(dspace_id, error)
!---Close file and HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_write_1d_r8
!------------------------------------------------------------------------------
!> integer(i4) 1D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_1d_i4(array,filename,path)
integer(i4), intent(in) :: array(:) !< Values to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical :: write_single
integer(i4) :: error
integer(i4), parameter :: one=1
integer(HID_T) :: file_id,dspace_id,dset_id
INTEGER(HSIZE_T), DIMENSION(1) :: dims
DEBUG_STACK_PUSH
!---Remove field if it already exists
IF(hdf5_field_exist(filename,path))THEN
  CALL oft_warn('Overwriting existing object "'//TRIM(path)//'" in file '//TRIM(filename))
  CALL hdf5_delete_obj(filename,path)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_write_1d_i4',__FILE__)
!---Create dataset and perform write
dims=SHAPE(array)
call h5screate_simple_f(one,dims,dspace_id,error)
call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_INTEGER, &
  dspace_id, dset_id, error)
call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, error)
!---Close dataset
call h5dclose_f(dset_id, error)
call h5sclose_f(dspace_id, error)
!---Close file and HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_write_1d_i4
!------------------------------------------------------------------------------
!> real(r8) 2D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_2d_r8(array,filename,path,single_prec)
real(r8), intent(in) :: array(:,:) !< Values to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(in) :: single_prec !< Save as single precision?
logical :: write_single
integer(i4) :: error
integer(i4), parameter :: two=2
integer(HID_T) :: file_id,dspace_id,dset_id
INTEGER(HSIZE_T), DIMENSION(2) :: dims
DEBUG_STACK_PUSH
!---Remove field if it already exists
IF(hdf5_field_exist(filename,path))THEN
  CALL oft_warn('Overwriting existing object "'//TRIM(path)//'" in file '//TRIM(filename))
  CALL hdf5_delete_obj(filename,path)
END IF
!---Write at lower precision?
write_single=.FALSE.
IF(PRESENT(single_prec))write_single=single_prec
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_write_2d_r8',__FILE__)
!---Create dataset and perform write
dims=SHAPE(array)
call h5screate_simple_f(two,dims,dspace_id,error)
IF(write_single)THEN
  call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_REAL, &
    dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, REAL(array,4), dims, error)
ELSE
  call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_DOUBLE, &
    dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, error)
END IF
!---Close dataset
call h5dclose_f(dset_id, error)
call h5sclose_f(dspace_id, error)
!---Close file and HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_write_2d_r8
!------------------------------------------------------------------------------
!> integer(i4) 2D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_2d_i4(array,filename,path)
integer(i4), intent(in) :: array(:,:) !< Values to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical :: write_single
integer(i4) :: error
integer(i4), parameter :: two=2
integer(HID_T) :: file_id
integer(HID_T) :: dspace_id,dset_id
INTEGER(HSIZE_T), DIMENSION(2) :: dims
DEBUG_STACK_PUSH
!---Remove field if it already exists
IF(hdf5_field_exist(filename,path))THEN
  CALL oft_warn('Overwriting existing object "'//TRIM(path)//'" in file '//TRIM(filename))
  CALL hdf5_delete_obj(filename,path)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
IF(error/=0)CALL oft_abort('Error opening file','hdf5_write_2d_i4',__FILE__)
!---Create dataset and perform write
dims=SHAPE(array)
call h5screate_simple_f(two,dims,dspace_id,error)
call h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_INTEGER, &
  dspace_id, dset_id, error)
call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims, error)
!---Close dataset
call h5dclose_f(dset_id, error)
call h5sclose_f(dspace_id, error)
!---Close file and HDF5
call h5fclose_f(file_id, error)
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_write_2d_i4
!------------------------------------------------------------------------------
!> FE vector implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_write_rst(rst_info,filename,path)
type(hdf5_rst), intent(in) :: rst_info !< Restart data (structure containing data and mapping)
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
integer(i4) :: i,error,one=1
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: filespace,memspace
integer(HSIZE_T), DIMENSION(1) :: dims,count
integer(HSSIZE_T), DIMENSION(1) :: offset
real(r8), pointer :: data(:)
#ifdef HAVE_PHDF5
integer(i4) :: info
integer(HID_T) :: plist_id
#else
integer(i4) :: j
#endif
DEBUG_STACK_PUSH
if(oft_debug_print(2))write(*,'(3A)')oft_indent,'Writing FE vector to file: ',TRIM(filename)
!---Remove field if it already exists
IF(oft_env%head_proc)THEN
  IF(hdf5_field_exist(filename,path))THEN
    CALL oft_warn('Overwriting existing object "'//TRIM(path)//'" in file '//TRIM(filename))
    CALL hdf5_delete_obj(filename,path)
  END IF
END IF
CALL oft_mpi_barrier(error)
!---Setup local chunk
dims=rst_info%dim
count=rst_info%count
offset=rst_info%offset
data=>rst_info%data
!---Initialize HDF5 interface
call h5open_f(error)
!---Wait up to 30 seconds for file to be available
DO i=1,11
  IF(oft_file_exist(filename))EXIT
  error=oft_sleep(1)
END DO
IF(i>10)CALL oft_abort("Exceeded timeout waiting for output file", "hdf5_write_rst", __FILE__)
!---Save to file
if(rst_info%full)then
  IF(oft_env%rank==0)THEN
    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
    IF(error/=0)CALL oft_abort('Error opening file', 'hdf5_write_rst', __FILE__)
    !---Create dataset
    CALL h5screate_simple_f(one, dims, filespace, error)
    CALL h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_DOUBLE, filespace, dset_id, error)
    !---Write data to file
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)
    !---Close file and HDF5 internal repesentations
    CALL h5sclose_f(filespace, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
  END IF
else
#ifdef HAVE_MPI
#ifdef HAVE_PHDF5
  info = MPI_INFO_NULL
  !---Open file
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
  CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error, access_prp=plist_id)
  IF(error/=0)CALL oft_abort('Error opening file', 'hdf5_write_rst', __FILE__)
  CALL h5pclose_f(plist_id, error)
  !---Create dataset
  CALL h5screate_simple_f(one, dims, filespace, error)
  CALL h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_DOUBLE, filespace, dset_id, error)
  CALL h5sclose_f(filespace, error)
  !---Define local chunk
  CALL h5screate_simple_f(one, count, memspace, error)
  CALL h5dget_space_f(dset_id, filespace, error)
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
  !---Setup parallel I/O
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  !---Write data to file
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, &
                  file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
  !---Close file and HDF5 internal repesentations
  CALL h5sclose_f(filespace, error)
  CALL h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error)
  call h5pclose_f(plist_id, error)
  call h5fclose_f(file_id, error)
#else
  if(oft_env%rank==0)then
    !---Open file and create dataset
    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
    IF(error/=0)CALL oft_abort('Error opening file', 'hdf5_write_rst', __FILE__)
    CALL h5screate_simple_f(one, dims, filespace, error)
    CALL h5dcreate_f(file_id, "/"//TRIM(path), H5T_NATIVE_DOUBLE, filespace, dset_id, error)
    !---Close file and HDF5 internal repesentations
    CALL h5sclose_f(filespace, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
  end if
  !---Loop over processors outputing chunks
  do i=1,oft_env%nprocs
    CALL oft_mpi_barrier(error)
    if(oft_env%rank+1==i)then
      !---Open file and dataset
      CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, file_id, error)
      CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
      !---Define local chunk
      CALL h5screate_simple_f(one, count, memspace, error)
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
      !---Write data to file
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, &
        file_space_id=filespace, mem_space_id=memspace)
      !---Close file and HDF5 internal repesentations
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5fclose_f(file_id, error)
    end if
  end do
#endif
#else
  call oft_abort("Requested distributed write without MPI","hdf5_write_rst",__FILE__)
#endif
end if
!---Finalize HDF5 interface
call h5close_f(error)
!---Synchronize processes
CALL oft_mpi_barrier(error)
DEBUG_STACK_POP
end subroutine hdf5_write_rst
!------------------------------------------------------------------------------
!> real(r8) scalar implementation of \ref oft_io::hdf5_read
!------------------------------------------------------------------------------
subroutine hdf5_read_scalar_r8(val,filename,path,success)
real(r8), intent(out) :: val !< Value to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
real(r8) :: tmpval(1)
CALL hdf5_read_1d_r8(tmpval,filename,path,success)
val=tmpval(1)
end subroutine hdf5_read_scalar_r8
!------------------------------------------------------------------------------
!> integer(i4) scalar implementation of \ref oft_io::hdf5_read
!------------------------------------------------------------------------------
subroutine hdf5_read_scalar_i4(val,filename,path,success)
integer(i4), intent(out) :: val !< Value to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: tmpval(1)
CALL hdf5_read_1d_i4(tmpval,filename,path,success)
val=tmpval(1)
end subroutine hdf5_read_scalar_i4
!------------------------------------------------------------------------------
!> real(r8) 1D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_read_1d_r8(array,filename,path,success)
real(r8), intent(inout) :: array(:) !< Values to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: error
integer(i4), parameter :: one=1,zero=0
integer(HID_T) :: file_id,dset_id
INTEGER(HSIZE_T), DIMENSION(1) :: dims
! DEBUG_STACK_PUSH
IF(PRESENT(success))THEN
  success=.FALSE.
  CALL h5eset_auto_f(zero, error)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
IF(error/=0)GOTO 102
CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
IF(error/=0)GOTO 101
!---
dims=SHAPE(array)
call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, error)
IF(error/=0)GOTO 100
!---Close and finalize HDF5
call h5dclose_f(dset_id, error)
call h5fclose_f(file_id, error)
call h5close_f(error)
! DEBUG_STACK_POP
IF(PRESENT(success))THEN
  success=.TRUE.
  CALL h5eset_auto_f(one, error)
END IF
RETURN
100 CALL h5dclose_f(dset_id, error)
101 CALL h5fclose_f(file_id, error)
102 CALL h5close_f(error)
IF(PRESENT(success))CALL h5eset_auto_f(one, error)
end subroutine hdf5_read_1d_r8
!------------------------------------------------------------------------------
!> integer(i8) 1D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_read_1d_i4(array,filename,path,success)
integer(i4), intent(inout) :: array(:) !< Values to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: error
integer(i4), parameter :: one=1,zero=0
integer(HID_T) :: file_id,dset_id
INTEGER(HSIZE_T), DIMENSION(1) :: dims
! DEBUG_STACK_PUSH
IF(PRESENT(success))THEN
  success=.FALSE.
  CALL h5eset_auto_f(zero, error)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
IF(error/=0)GOTO 102
CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
IF(error/=0)GOTO 101
!---
dims=SHAPE(array)
call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, error)
IF(error/=0)GOTO 100
!---Close and finalize HDF5
call h5dclose_f(dset_id, error)
call h5fclose_f(file_id, error)
call h5close_f(error)
! DEBUG_STACK_POP
IF(PRESENT(success))THEN
  success=.TRUE.
  CALL h5eset_auto_f(one, error)
END IF
RETURN
100 CALL h5dclose_f(dset_id, error)
101 CALL h5fclose_f(file_id, error)
102 CALL h5close_f(error)
IF(PRESENT(success))CALL h5eset_auto_f(one, error)
end subroutine hdf5_read_1d_i4
!------------------------------------------------------------------------------
!> real(r8) 2D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_read_2d_r8(array,filename,path,success)
real(r8), intent(inout) :: array(:,:) !< Values to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: error
integer(i4), parameter :: two=2,one=1,zero=0
integer(HID_T) :: file_id,dset_id
INTEGER(HSIZE_T), DIMENSION(2) :: dims
! DEBUG_STACK_PUSH
IF(PRESENT(success))THEN
  success=.FALSE.
  CALL h5eset_auto_f(zero, error)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
IF(error/=0)GOTO 102
CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
IF(error/=0)GOTO 101
!---
dims=SHAPE(array)
call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, error)
IF(error/=0)GOTO 100
!---Close and finalize HDF5
call h5dclose_f(dset_id, error)
call h5fclose_f(file_id, error)
call h5close_f(error)
! DEBUG_STACK_POP
IF(PRESENT(success))THEN
  success=.TRUE.
  CALL h5eset_auto_f(one, error)
END IF
RETURN
100 CALL h5dclose_f(dset_id, error)
101 CALL h5fclose_f(file_id, error)
102 CALL h5close_f(error)
IF(PRESENT(success))CALL h5eset_auto_f(one, error)
end subroutine hdf5_read_2d_r8
!------------------------------------------------------------------------------
!> integer(i4) 2D array implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_read_2d_i4(array,filename,path,success)
integer(i4), intent(inout) :: array(:,:) !< Values to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: error
integer(i4), parameter :: two=2,one=1,zero=0
integer(HID_T) :: file_id,dset_id
INTEGER(HSIZE_T), DIMENSION(2) :: dims
! DEBUG_STACK_PUSH
IF(PRESENT(success))THEN
  success=.FALSE.
  CALL h5eset_auto_f(zero, error)
END IF
!---Initialize HDF5 and open file
call h5open_f(error)
call h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
IF(error/=0)GOTO 102
CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
IF(error/=0)GOTO 101
!---
dims=SHAPE(array)
call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, error)
IF(error/=0)GOTO 100
!---Close and finalize HDF5
call h5dclose_f(dset_id, error)
call h5fclose_f(file_id, error)
call h5close_f(error)
! DEBUG_STACK_POP
IF(PRESENT(success))THEN
  success=.TRUE.
  CALL h5eset_auto_f(one, error)
END IF
RETURN
100 CALL h5dclose_f(dset_id, error)
101 CALL h5fclose_f(file_id, error)
102 CALL h5close_f(error)
IF(PRESENT(success))CALL h5eset_auto_f(one, error)
end subroutine hdf5_read_2d_i4
!------------------------------------------------------------------------------
!> FE vector implementation of \ref oft_io::hdf5_write
!------------------------------------------------------------------------------
subroutine hdf5_read_rst(rst_info,filename,path,success)
type(hdf5_rst), intent(inout) :: rst_info !< Restart data (structure containing mapping and data holder)
character(*), intent(in) :: filename !< Path to file
character(*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: i,error
integer(i4), parameter :: two=2,one=1,zero=0
integer :: access_flag
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: filespace,memspace
integer(HSIZE_T) :: space_count
integer(HSIZE_T), DIMENSION(1) :: dims,count
integer(HSSIZE_T), DIMENSION(1) :: offset
real(r8), pointer :: data(:)
logical :: exists
#ifdef HAVE_PHDF5
integer(i4) :: info
integer(HID_T) :: plist_id
#endif
! DEBUG_STACK_PUSH
IF(PRESENT(success))THEN
  success=.FALSE.
  CALL h5eset_auto_f(zero, error)
END IF
!---Set filename and check for file
CALL oft_mpi_barrier(error)
!---Set access flag
access_flag=H5F_ACC_RDONLY_F
if(oft_debug_print(2))write(*,'(2X,2A)')'Reading field from restart file: ',TRIM(filename)
!---Wait up to 30 seconds for file to be available
DO i=1,31
  IF(oft_file_exist(filename))EXIT
  error=oft_sleep(1)
END DO
IF(i>30)CALL oft_abort("Exceeded timeout waiting for output file", "hdf5_read_rst", __FILE__)
!---Read data
if(rst_info%full)then
  !---Initialize HDF5 and open vector dump file
  call h5open_f(error)
  call h5fopen_f(TRIM(filename), access_flag, file_id, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 103
    CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
  END IF
  !---
  dims=rst_info%dim
  count=rst_info%count
  offset=0
  data=>rst_info%data
  !---
  call h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 102
    CALL oft_abort('Error opening dataset','hdf5_read_rst',__FILE__)
  END IF
  !---Check size
  CALL h5dget_space_f(dset_id, memspace, error)
  CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 101
    CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
  END IF
  IF(space_count/=rst_info%dim)THEN
    IF(PRESENT(success))GOTO 101
    CALL oft_abort('Dataset size does not match','hdf5_read_rst',__FILE__)
  END IF
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 100
    CALL oft_abort('Error reading dataset','hdf5_read_rst',__FILE__)
  END IF
  !---Close dataset
  CALL h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error)
  !---Close vector dump file and finalize HDF5
  call h5fclose_f(file_id, error)
  call h5close_f(error)
else
#ifdef HAVE_MPI
#ifdef HAVE_PHDF5
  info = MPI_INFO_NULL
  !---Initialize HDF5 and open restart file
  call h5open_f(error)
  !---Initliazie parallel I/O
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, int(oft_env%comm,4), info, error)
  !---
  CALL h5fopen_f(TRIM(filename), access_flag, file_id, error, access_prp=plist_id)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 103
    CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
  END IF
  CALL h5pclose_f(plist_id, error)
  !---
  dims=rst_info%dim
  count=rst_info%count
  offset=rst_info%offset
  data=>rst_info%data
  !---Write out cell data
  CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 102
    CALL oft_abort('Error opening dataset','hdf5_read_rst',__FILE__)
  END IF
  !---Check size
  CALL h5dget_space_f(dset_id, memspace, error)
  CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 101
    CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
  END IF
  IF(space_count/=rst_info%dim)THEN
    IF(PRESENT(success))GOTO 101
    CALL oft_abort('Dataset size does not match','hdf5_read_rst',__FILE__)
  END IF
  !---Write out cell data
  CALL h5screate_simple_f(one, count, memspace, error)
  !---Write out cell data
  CALL h5dget_space_f(dset_id, filespace, error)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
  !---Write out cell data
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, &
    file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
  IF(error/=0)THEN
    IF(PRESENT(success))GOTO 100
    CALL oft_abort('Error reading dataset','hdf5_read_rst',__FILE__)
  END IF
  !---Write out cell data
  CALL h5sclose_f(filespace, error)
  CALL h5sclose_f(memspace, error)
  !---Write out cell data
  call h5dclose_f(dset_id, error)
  call h5pclose_f(plist_id, error)
  !---Write out cell data
  call h5fclose_f(file_id, error)
  !---Synchronize processes
  CALL oft_mpi_barrier(error)
  call h5close_f(error)
#else
  !---Initialize HDF5 and open vector dump file
  call h5open_f(error)
  !---
  dims=rst_info%dim
  count=rst_info%count
  offset=rst_info%offset
  data=>rst_info%data
  do i=1,oft_env%nprocs
    if(oft_env%rank+1==i)then
      !---
      call h5fopen_f(TRIM(filename), access_flag, file_id, error)
      IF(error/=0)THEN
        IF(PRESENT(success))GOTO 103
        CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
      END IF
      !---
      CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
      IF(error/=0)THEN
        IF(PRESENT(success))GOTO 102
        CALL oft_abort('Error opening dataset','hdf5_read_rst',__FILE__)
      END IF
      !---Check size
      CALL h5dget_space_f(dset_id, memspace, error)
      CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
      IF(error/=0)THEN
        IF(PRESENT(success))GOTO 101
        CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
      END IF
      IF(space_count/=rst_info%dim)THEN
        IF(PRESENT(success))GOTO 101
        CALL oft_abort('Dataset size does not match','hdf5_read_rst',__FILE__)
      END IF
      !---
      CALL h5screate_simple_f(one, count, memspace, error)
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
      !---
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, &
        file_space_id=filespace, mem_space_id=memspace)
      IF(error/=0)THEN
        IF(PRESENT(success))GOTO 100
        CALL oft_abort('Error reading dataset','hdf5_read_rst',__FILE__)
      END IF
      !---Close dataset
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      !---Close dataset
      call h5dclose_f(dset_id, error)
      !---Close vector dump file and finalize HDF5
      call h5fclose_f(file_id, error)
    end if
    CALL oft_mpi_barrier(error)
  end do
  call h5close_f(error)
#endif
#else
  call oft_abort("Requested distributed read without MPI","hdf5_read_rst",__FILE__)
#endif
end if
IF(PRESENT(success))THEN
  success=.TRUE.
  CALL h5eset_auto_f(one, error)
END IF
! DEBUG_STACK_POP
RETURN
100 CALL h5sclose_f(memspace, error)
#ifdef HAVE_MPI
IF(.NOT.rst_info%full)THEN
#ifdef HAVE_PHDF5
  CALL h5pclose_f(plist_id, error)
#endif
  CALL h5sclose_f(filespace, error)
END IF
#endif
101 CALL h5dclose_f(dset_id, error)
102 CALL h5fclose_f(file_id, error)
103 CALL h5close_f(error)
IF(PRESENT(success))CALL h5eset_auto_f(one, error)
#ifdef HAVE_MPI
IF(.NOT.rst_info%full)THEN
#ifdef HAVE_PHDF5
  CALL oft_mpi_barrier(error)
#else
  do i=oft_env%rank+1,oft_env%nprocs
    CALL oft_mpi_barrier(error)
  end do
#endif
END IF
#endif
! DEBUG_STACK_POP
end subroutine hdf5_read_rst
!------------------------------------------------------------------------------
!> Deallocate internal storage fields created for HDF5 collective I/O
!------------------------------------------------------------------------------
subroutine hdf5_rst_destroy(self)
type(hdf5_rst), intent(inout) :: self
IF(ASSOCIATED(self%data))DEALLOCATE(self%data)
end subroutine hdf5_rst_destroy
end module oft_io
