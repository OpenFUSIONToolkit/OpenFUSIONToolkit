!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
MODULE oft_io
USE oft_base
USE hdf5
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Binary output file object
!---------------------------------------------------------------------------
TYPE :: oft_bin_file
  LOGICAL :: header_written = .FALSE. !< Was header written during this run
  INTEGER(i4) :: io_unit = 0 !< Output unit for file
  INTEGER(i4) :: ncomm = 0 !< Number of comment lines
  INTEGER(i4) :: nfields = 0 !< Number of fields
  INTEGER(i4) :: nbytes = 0 !< Number of bytes per line
  INTEGER(i4), POINTER, DIMENSION(:) :: field_size => NULL() !< Dimension of each field
  CHARACTER(LEN=OFT_PATH_SLEN) :: filename = '' !< Output filename
  CHARACTER(LEN=80) :: filedesc = '' !< Description string
  CHARACTER(LEN=2), POINTER, DIMENSION(:) :: field_type => NULL() !< Field types
  CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: field_names => NULL() !< Field names
  CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: field_desc => NULL() !< Field descriptions
  CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: comm_lines => NULL() !< Comment lines
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
!---------------------------------------------------------------------------
!> HDF5 restart structure
!!
!! Contains definition and data for saving/loading to distributed restart
!! files.
!!
!! @note Local ownership in HDF5 files is different from linear algebra
!! ownership.
!---------------------------------------------------------------------------
type :: hdf5_rst
  logical :: full = .TRUE. !< Distributed data flag
  integer(i4) :: count = 0 !< Number of values owned by local process
  integer(i8) :: offset = 0 !< Offset in global array
  integer(i8) :: dim = 0 !< Length of global array
  real(r8), pointer, dimension(:) :: data => NULL() !< Array holding local data
end type hdf5_rst
!---------------------------------------------------------------------------
!> Write data to an HDF5 file
!---------------------------------------------------------------------------
INTERFACE hdf5_write
  MODULE PROCEDURE hdf5_write_scalar_r8
  MODULE PROCEDURE hdf5_write_scalar_i4
  MODULE PROCEDURE hdf5_write_1d_r8
  MODULE PROCEDURE hdf5_write_1d_i4
  MODULE PROCEDURE hdf5_write_2d_r8
  MODULE PROCEDURE hdf5_write_2d_i4
  MODULE PROCEDURE hdf5_write_rst
END INTERFACE hdf5_write
!---------------------------------------------------------------------------
!> Read data from an HDF5 file
!---------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Setup Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Open Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Close Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
SUBROUTINE bin_file_close(self)
CLASS(oft_bin_file), INTENT(inout) :: self
CLOSE(self%io_unit)
END SUBROUTINE bin_file_close
!------------------------------------------------------------------------------
!> Add field to Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
SUBROUTINE bin_file_add(self,fieldname,type_str,desc,fsize)
CLASS(oft_bin_file), INTENT(inout) :: self !< File object
CHARACTER(LEN=*), INTENT(in) :: fieldname !< Field name to add
CHARACTER(LEN=2), INTENT(in) :: type_str !< Type of field
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: desc !< Description of field
INTEGER(i4), OPTIONAL, INTENT(in) :: fsize !< Size of field
INTEGER(i4) :: fsize_tmp
INTEGER(i4), POINTER, DIMENSION(:) :: sizes_tmp
CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: names_tmp
CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: desc_tmp
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
IF(LEN(fieldname)>OFT_HIST_SLEN)CALL oft_abort('Name too long', &
  "bin_file_add",__FILE__)
IF(PRESENT(desc))THEN
  IF(LEN(desc)>OFT_HIST_SLEN)CALL oft_abort('Description too long', &
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
!------------------------------------------------------------------------------
!> Add comment to Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
SUBROUTINE bin_file_add_comm(self,comment)
CLASS(oft_bin_file), INTENT(inout) :: self !< File object
CHARACTER(LEN=*), INTENT(in) :: comment !< Comment to add
CHARACTER(LEN=OFT_HIST_SLEN), POINTER, DIMENSION(:) :: comm_tmp
IF(LEN(comment)>OFT_HIST_SLEN)CALL oft_abort('Comment too long', &
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
!------------------------------------------------------------------------------
!> Write header for Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Write single set of data to Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Flush I/O buffer for Open FUSION Toolkit binary I/O file
!------------------------------------------------------------------------------
SUBROUTINE bin_file_flush(self)
CLASS(oft_bin_file), INTENT(inout) :: self
FLUSH(self%io_unit)
END SUBROUTINE bin_file_flush
!---------------------------------------------------------------------------
!> Creates HDF5 output files for the current mesh
!!
!! The following files are created
!! - mesh.[PROC_RANK].h5
!! - scalar_dump.[PROC_RANK].h5
!! - vector_dump.[PROC_RANK].h5
!!
!! @note One output file is created per MPI task
!---------------------------------------------------------------------------
subroutine hdf5_create_files(basepath)
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix
INTEGER(4) :: ierr
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Creating HDF5 plot files'
pathprefix=''
IF(PRESENT(basepath))THEN
  IF(LEN(basepath)>OFT_PATH_SLEN)CALL oft_abort("Basepath too long", &
    "hdf5_create_files", __FILE__)
  pathprefix=basepath
END IF
CALL oft_increase_indent
CALL hdf5_create_file(TRIM(pathprefix)//"scalar_dump."//hdf5_proc_str()//".h5")
CALL hdf5_create_file(TRIM(pathprefix)//"vector_dump."//hdf5_proc_str()//".h5")
CALL hdf5_create_file(TRIM(pathprefix)//"mesh."//hdf5_proc_str()//".h5")
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine hdf5_create_files
!---------------------------------------------------------------------------
!> Adds a timestep to the dump metadata file.
!!
!! Subsequent output will be added to this timestep until another call
!! to this subroutine
!---------------------------------------------------------------------------
subroutine hdf5_create_timestep(t,basepath)
real(r8), intent(in) :: t !< Time value
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix
integer(i4) :: io_unit
DEBUG_STACK_PUSH
if(oft_debug_print(1))write(*,'(2A,ES11.4)')oft_indent,'Creating plot time: ',t
hdf5_ts=hdf5_ts+1
if(oft_env%rank==0)then
  pathprefix=''
  IF(PRESENT(basepath))THEN
    IF(LEN(basepath)>OFT_PATH_SLEN)CALL oft_abort("Basepath too long", &
      "hdf5_create_timestep", __FILE__)
    pathprefix=basepath
  END IF
  OPEN(NEWUNIT=io_unit,FILE=TRIM(pathprefix)//'dump.dat',POSITION="APPEND",STATUS="OLD")
  WRITE(io_unit,*)
  WRITE(io_unit,*)'Time Step',REAL(t,4)
  WRITE(io_unit,*)'Field Data'
  CLOSE(io_unit)
end if
DEBUG_STACK_POP
end subroutine hdf5_create_timestep
!---------------------------------------------------------------------------
!> Get processor rank as string for HDF5 I/O
!---------------------------------------------------------------------------
function hdf5_proc_str(proc_ind) result(proc_str)
integer(i4), optional, intent(in) :: proc_ind
character(LEN=HDF5_TLEN) :: proc_str
100 FORMAT (I HDF5_TLEN.HDF5_TLEN)
IF(PRESENT(proc_ind))THEN
  write(proc_str,100)proc_ind+1
ELSE
  write(proc_str,100)oft_env%rank+1
END IF
end function hdf5_proc_str
!---------------------------------------------------------------------------
!> Get timestep index as string for HDF5 I/O
!---------------------------------------------------------------------------
function hdf5_ts_str() result(ts_str)
character(LEN=HDF5_TLEN) :: ts_str
100 FORMAT (I HDF5_TLEN.HDF5_TLEN)
write(ts_str,100)hdf5_ts
end function hdf5_ts_str
!---------------------------------------------------------------------------
!> Test for exitence of a file
!!
!! @result Logical flag indicating existence
!---------------------------------------------------------------------------
function oft_file_exist(filepath) result(exists)
character(LEN=*), intent(in) :: filepath !< Path to file
logical :: exists
INQUIRE(FILE=TRIM(filepath),EXIST=exists)
end function oft_file_exist
!---------------------------------------------------------------------------
!> Test for exitence of a field in a HDF5 file
!!
!! @result Logical flag indicating existence of field and file
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Test for exitence of a field in a HDF5 file
!!
!! @result Logical flag indicating existence of field and file
!---------------------------------------------------------------------------
subroutine hdf5_field_get_sizes(filepath,path,ndims,dim_sizes)
character(LEN=*), intent(in) :: filepath !< Path to file
character(LEN=*), intent(in) :: path !< Path of field in file
integer(4), intent(out)  :: ndims !< Number of dimensions in field
integer(4), allocatable, dimension(:), intent(out) :: dim_sizes !< Size of each dimension
integer(4) :: access_flag,error
integer(hsize_t), allocatable, dimension(:) :: tmp_sizes,maxdims
integer(HID_T) :: file_id,dset_id,dspace_id
DEBUG_STACK_PUSH
ndims=-1
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
!---------------------------------------------------------------------------
!> Create an empty HDF5 output file
!---------------------------------------------------------------------------
subroutine hdf5_create_file(filename)
character(LEN=*), intent(in) :: filename !< Name of file to be created
integer(4) :: error
integer(HID_T) :: file_id
DEBUG_STACK_PUSH
!---Get processor rank for file creation
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Creating HDF5 output file: ',TRIM(filename)
!---Initialize HDF5 system
call h5open_f(error)
!---Create HDF5 file
call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error)
call h5fclose_f(file_id, error)
!---Finalize HDF5 system
call h5close_f(error)
DEBUG_STACK_POP
end subroutine hdf5_create_file
!---------------------------------------------------------------------------
!> Create HDF5 group in existing file
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> Add string attribute to existing object (group or dataset)
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> real(r8) scalar implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
subroutine hdf5_write_scalar_r8(val,filename,path,single_prec)
real(r8), intent(in) :: val !< Value to write to file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(in) :: single_prec !< Save as single precision?
real(r8) :: tmpval(1)
tmpval(1)=val
CALL hdf5_write_1d_r8(tmpval,filename,path,single_prec)
end subroutine hdf5_write_scalar_r8
!---------------------------------------------------------------------------
!> integer(i4) scalar implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> real(r8) 1D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> integer(i4) 1D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> real(r8) 2D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> integer(i4) 2D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> FE vector implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
subroutine hdf5_write_rst(rst_info,filename,path,append)
type(hdf5_rst), intent(in) :: rst_info !< Restart data (structure containing data and mapping)
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(in) :: append !< Append to or create file? (optional)
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
CALL oft_mpi_barrier(error)
if(oft_debug_print(2))write(*,'(3A)')oft_indent,'Writing FE vector to file: ',TRIM(filename)
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
!---------------------------------------------------------------------------
!> real(r8) scalar implementation of \ref oft_io::hdf5_read
!---------------------------------------------------------------------------
subroutine hdf5_read_scalar_r8(val,filename,path,success)
real(r8), intent(out) :: val !< Value to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
real(r8) :: tmpval(1)
CALL hdf5_read_1d_r8(tmpval,filename,path,success)
val=tmpval(1)
end subroutine hdf5_read_scalar_r8
!---------------------------------------------------------------------------
!> integer(i4) scalar implementation of \ref oft_io::hdf5_read
!---------------------------------------------------------------------------
subroutine hdf5_read_scalar_i4(val,filename,path,success)
integer(i4), intent(out) :: val !< Value to read from file
character(LEN=*), intent(in) :: filename !< Path to file
character(LEN=*), intent(in) :: path !< Variable path in file
logical, optional, intent(out) :: success !< Successful read?
integer(i4) :: tmpval(1)
CALL hdf5_read_1d_i4(tmpval,filename,path,success)
val=tmpval(1)
end subroutine hdf5_read_scalar_i4
!---------------------------------------------------------------------------
!> real(r8) 1D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> integer(i8) 1D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> real(r8) 2D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> integer(i4) 2D array implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
!> FE vector implementation of \ref oft_io::hdf5_write
!---------------------------------------------------------------------------
subroutine hdf5_read_rst(rst_info,filename,path)
type(hdf5_rst), intent(inout) :: rst_info !< Restart data (structure containing mapping and data holder)
character(*), intent(in) :: filename !< Path to file
character(*), intent(in) :: path !< Variable path in file
integer(i4) :: i,error,one=1
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
DEBUG_STACK_PUSH
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
  IF(error/=0)CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
  !---
  dims=rst_info%dim
  count=rst_info%count
  offset=0
  data=>rst_info%data
  !---
  call h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
  IF(error/=0)CALL oft_abort('Dataset not found','hdf5_read_rst',__FILE__)
  !---Check size
  CALL h5dget_space_f(dset_id, memspace, error)
  CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
  IF(error/=0)CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
  IF(space_count/=rst_info%dim)CALL oft_abort('Dataset size does not match', &
    'hdf5_read_rst',__FILE__)
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error)
  IF(error/=0)CALL oft_abort('Error in HDF5 read','hdf5_read_rst',__FILE__)
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
  IF(error/=0)CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
  CALL h5pclose_f(plist_id, error)
  !---
  dims=rst_info%dim
  count=rst_info%count
  offset=rst_info%offset
  data=>rst_info%data
  !---Write out cell data
  CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
  IF(error/=0)CALL oft_abort('Dataset not found','hdf5_read_rst',__FILE__)
  !---Check size
  CALL h5dget_space_f(dset_id, memspace, error)
  CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
  IF(error/=0)CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
  IF(space_count/=rst_info%dim)CALL oft_abort('Dataset size does not match', &
    'hdf5_read_rst',__FILE__)
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
  IF(error/=0)CALL oft_abort('Error in HDF5 read','hdf5_read_rst',__FILE__)
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
      IF(error/=0)CALL oft_abort('Error opening file','hdf5_read_rst',__FILE__)
      !---
      CALL h5dopen_f(file_id, "/"//TRIM(path), dset_id, error)
      IF(error/=0)CALL oft_abort('Dataset not found','hdf5_read_rst',__FILE__)
      !---Check size
      CALL h5dget_space_f(dset_id, memspace, error)
      CALL h5sget_simple_extent_npoints_f(memspace, space_count, error)
      IF(error/=0)CALL oft_abort('Could not read dataset size','hdf5_read_rst',__FILE__)
      IF(space_count/=rst_info%dim)CALL oft_abort('Dataset size does not match', &
        'hdf5_read_rst',__FILE__)
      !---
      CALL h5screate_simple_f(one, count, memspace, error)
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
      !---
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, &
        file_space_id=filespace, mem_space_id=memspace)
      IF(error/=0)CALL oft_abort('Error in HDF5 read','hdf5_read_rst',__FILE__)
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
DEBUG_STACK_POP
end subroutine hdf5_read_rst
!---------------------------------------------------------------------------
!> Deallocate internal storage fields created for HDF5 collective I/O
!---------------------------------------------------------------------------
subroutine hdf5_rst_destroy(self)
type(hdf5_rst), intent(inout) :: self
IF(ASSOCIATED(self%data))DEALLOCATE(self%data)
end subroutine hdf5_rst_destroy
!---------------------------------------------------------------------------
!> Writes the metadata file referencing all of the mesh blocks
!!
!! - Writes out mesh sizes
!! - Initiliazes timestep 0
!---------------------------------------------------------------------------
subroutine oft_hdf5_write_dump(mesh_type,vol_sizes,surf_sizes,basepath)
integer(i4), intent(in) :: mesh_type !< Mesh type flag (Tet/Tri or Hex/Quad)
integer(i4), intent(in) :: vol_sizes(2) !< Volume mesh counts (np,nc)
integer(i4), intent(in) :: surf_sizes(2) !< Surface mesh counts (np,nc)
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix
integer(i4) :: i,ntrans(4),ierr,io_unit
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
type(mpi_status) :: mpi_stat
#else
integer(i4) :: mpi_stat(MPI_STATUS_SIZE)
#endif
#endif
DEBUG_STACK_PUSH
pathprefix=''
IF(PRESENT(basepath))THEN
  IF(LEN(basepath)>OFT_PATH_SLEN)CALL oft_abort("Basepath too long", &
    "oft_hdf5_write_dump", __FILE__)
  pathprefix=basepath
END IF
!---Get local mesh counts
ntrans(1:2)=vol_sizes
ntrans(3:4)=surf_sizes
!---Send to lead proc
#ifdef HAVE_MPI
if(oft_env%rank>0)call MPI_SEND(ntrans,4,OFT_MPI_I4,0,1,MPI_COMM_WORLD,ierr)
#endif
!---On lead proc setup metadata file
if(oft_env%rank==0)then
  !---Setup file
  open(NEWUNIT=io_unit,FILE=TRIM(pathprefix)//'dump.dat')
  write(io_unit,*)'Mesh Data'
  !---Local mesh size
  write(io_unit,*)mesh_type
  write(io_unit,*)hdf5_proc_str(),ntrans
  !---Mesh sizes of other tasks
#ifdef HAVE_MPI
  do i=1,oft_env%nprocs-1
    call MPI_RECV(ntrans,4,OFT_MPI_I4,i,1,MPI_COMM_WORLD,mpi_stat,ierr)
    write(io_unit,*)hdf5_proc_str(i),ntrans
  end do
#endif
  !---Initialize timestep 0
  write(io_unit,*)
  write(io_unit,*)'Time Step',real(0.d0)
  write(io_unit,*)'Field Data'
  close(io_unit)
end if
DEBUG_STACK_POP
end subroutine oft_hdf5_write_dump
!---------------------------------------------------------------------------
!> Adds an output field to the dump file for Xdmf construction
!---------------------------------------------------------------------------
subroutine oft_hdf5_add_dump(tag,type,basepath)
character(LEN=*), intent(in) :: tag !< Name of the field to add
integer(i4), intent(in) :: type !< Type of field being output
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix
integer(i4) :: io_unit
DEBUG_STACK_PUSH
if(oft_env%rank==0)then
  pathprefix=''
  IF(PRESENT(basepath))THEN
    IF(LEN(basepath)>OFT_PATH_SLEN)CALL oft_abort("Basepath too long", &
      "oft_hdf5_add_dump", __FILE__)
    pathprefix=basepath
  END IF
  open(NEWUNIT=io_unit,FILE=TRIM(pathprefix)//'dump.dat',POSITION="APPEND",STATUS="OLD")
  write(io_unit,*)tag,type
  close(io_unit)
end if
DEBUG_STACK_POP
end subroutine oft_hdf5_add_dump
end module oft_io
