!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_local.F90
!
!> @brief Machine and compiler specific settings.
!!
!! Machine and compiler specific settings and global constants.
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_base
!-----------------------------------------------------------------------------
MODULE oft_local
USE, INTRINSIC :: iso_c_binding, only: c_int, c_ptr, c_long
#ifdef __INTEL_COMPILER
USE ifport ! Intel fortran portability library
#endif
#ifdef HAVE_XML
USE fox_dom, ONLY: xml_node => node, xml_parsefile => parsefile, xml_hasAttribute => hasAttribute, &
  xml_extractDataAttribute => extractDataAttribute, xml_extractDataContent => extractDataContent
#endif
IMPLICIT NONE
!---Local types sizes
INTEGER, PARAMETER :: i4=SELECTED_INT_KIND(9)           !< 4-Byte integer spec
INTEGER, PARAMETER :: i8=SELECTED_INT_KIND(18)          !< 8-Byte integer spec
INTEGER, PARAMETER :: r4=SELECTED_REAL_KIND(6,37)       !< Single precision float spec (4 Bytes)
INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(13,307)     !< Double precision float spec (8 Bytes)
INTEGER, PARAMETER :: r10=SELECTED_REAL_KIND(18,4900)   !< Extended precision float spec (10 or 16 Bytes depending on platform)
INTEGER, PARAMETER :: c4=r4                             !< Single precision complex spec (4 Bytes)
INTEGER, PARAMETER :: c8=r8                             !< Double precision complex spec (8 Bytes)
REAL(r8), PARAMETER :: pi=3.141592653589793238462643_r8 !< \f$ \pi \f$
!------------------------------------------------------------------------------
! Define PETSc address type
! - This is the integer value of the C memory pointer to a given object
! and replaces Vec, Mat, IS, etc. I use this instead of the preprocessed
! definitions to keep documentation clean and make the real types obvious.
!------------------------------------------------------------------------------
#ifdef HAVE_PETSC
#include "petscconf.h"
#if (PETSC_SIZEOF_VOID_P == 8)
INTEGER, PARAMETER :: petsc_addr=i8 !< Size of address pointer (32 or 64) bits
#else
INTEGER, PARAMETER :: petsc_addr=i4
#endif
#else
INTEGER, PARAMETER :: petsc_addr=i4
#endif
INTERFACE
!---------------------------------------------------------------------------
!> Interface to C sleep function
!---------------------------------------------------------------------------
  FUNCTION oft_sleep(seconds)  BIND(C,name="sleep")
    IMPORT c_int
    INTEGER(c_int) :: oft_sleep !< Error code on return
    INTEGER(c_int), INTENT(in), VALUE :: seconds !< Length of time to pause in seconds
  END FUNCTION oft_sleep
!---------------------------------------------------------------------------
!> Simple in-memory hashing function for dataset checksumming
!---------------------------------------------------------------------------
  FUNCTION oft_simple_hash(key,length)  BIND(C)
    IMPORT c_int, c_long, c_ptr
    INTEGER(c_int) :: oft_simple_hash !< Hash of data
    TYPE(c_ptr), VALUE, INTENT(in) :: key !< Location of data
    INTEGER(c_long), VALUE, INTENT(in) :: length !< Length of data to hash in bytes
  END FUNCTION oft_simple_hash
END INTERFACE
!---------------------------------------------------------------------------
!> One dimensional integer set
!---------------------------------------------------------------------------
TYPE :: oft_1d_int
  INTEGER(i4) :: n = 0 !< Number of values in set
  INTEGER(i4), POINTER, DIMENSION(:) :: v => NULL() !< Values
END TYPE oft_1d_int
!---------------------------------------------------------------------------
!> One dimensional real set
!---------------------------------------------------------------------------
TYPE :: oft_1d_real
  INTEGER(i4) :: n = 0 !< Number of values in set
  REAL(r8), POINTER, DIMENSION(:) :: v => NULL() !< Values
END TYPE oft_1d_real
!---------------------------------------------------------------------------
!> One dimensional complex set
!---------------------------------------------------------------------------
TYPE :: oft_1d_comp
  INTEGER(i4) :: n = 0 !< Number of values in set
  COMPLEX(c8), POINTER, DIMENSION(:) :: v => NULL() !< Values
END TYPE oft_1d_comp
!---------------------------------------------------------------------------
!> Generate inverse of sparse indexing
!---------------------------------------------------------------------------
INTERFACE get_inverse_map
  MODULE PROCEDURE get_inverse_map_i4
  MODULE PROCEDURE get_inverse_map_i8
END INTERFACE get_inverse_map
!
ABSTRACT INTERFACE
!---------------------------------------------------------------------------
!> Generic interface for 1D function
!---------------------------------------------------------------------------
  FUNCTION oft_1d_func(x) result(f)
    IMPORT r8
    REAL(r8), INTENT(in) :: x !< Parameter 1
    REAL(r8) :: f !< Function value
  END FUNCTION oft_1d_func
!---------------------------------------------------------------------------
!> Generic interface for 2D function
!---------------------------------------------------------------------------
  FUNCTION oft_2d_func(x,y) result(f)
    IMPORT r8
    REAL(r8), INTENT(in) :: x !< Parameter 1
    REAL(r8), INTENT(in) :: y !< Parameter 2
    REAL(r8) :: f !< Function value
  END FUNCTION oft_2d_func
END INTERFACE
!------------------------------------------------------------------------------
!> Simple timer class
!------------------------------------------------------------------------------
TYPE :: oft_timer
  INTEGER(i8) :: count = 0 !< Integer value of system clock at last call
CONTAINS
  !> Start or reset timer
  procedure :: tick => oft_timer_start
  !> Set elapsed time since tick/tock
  procedure :: tock => oft_timer_elapsed
  !> Get elapsed time since tick/tock in integer counts
  procedure :: int_tock => oft_timer_intelapsed
  !> Check if time since tick/tock exceeds a limit
  procedure :: timeout => oft_timer_timeout
END TYPE oft_timer
PRIVATE oft_timer_start, oft_timer_elapsed, oft_timer_intelapsed, oft_timer_timeout
#ifdef HAVE_XML
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE :: xml_node_ptr
  TYPE(xml_node), POINTER :: this => NULL()
END TYPE xml_node_ptr
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE :: xml_nodelist
  INTEGER(i4) :: n = 0
  TYPE(xml_node_ptr), POINTER, DIMENSION(:) :: nodes => NULL()
END TYPE xml_nodelist
!---------------------------------------------------------------------------
!> Generate inverse of sparse indexing
!---------------------------------------------------------------------------
INTERFACE xml_get_element
  MODULE PROCEDURE xml_get_element_single
  MODULE PROCEDURE xml_get_element_list
END INTERFACE xml_get_element
#endif
CONTAINS
!------------------------------------------------------------------------------
!> Returns the corresponding lowercase letter, if `c` is an uppercase
!! ASCII character, otherwise `c` itself.
!!
!! Reproduced from the [Fortran stdlib](https://stdlib.fortran-lang.org/index.html)
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION char_to_lower(c) result(t)
CHARACTER(len=1), INTENT(in) :: c !< Input character
CHARACTER(len=1) :: t !< Lowercase version
INTEGER(i4), PARAMETER :: wp=iachar('a')-iachar('A')
INTEGER(i4), PARAMETER :: BA=iachar('A')
INTEGER(i4), PARAMETER :: BZ=iachar('Z')
INTEGER(i4) :: k
k = ichar(c) 
IF(k>=BA.and.k<=BZ)k = k + wp 
t = char(k)
END FUNCTION char_to_lower
!------------------------------------------------------------------------------
!> Converts a string to all lowercase characters
!------------------------------------------------------------------------------
SUBROUTINE string_to_lower(c)
CHARACTER(len=*), INTENT(inout) :: c !< Input/output string
INTEGER(i4) :: i
DO i=1,LEN(c)
  c(i:i)=char_to_lower(c(i:i))
END DO
END SUBROUTINE string_to_lower
!------------------------------------------------------------------------------
!> Returns the corresponding uppercase letter, if `c` is an lowercase
!! ASCII character, otherwise `c` itself.
!!
!! Reproduced from the [Fortran stdlib](https://stdlib.fortran-lang.org/index.html)
!------------------------------------------------------------------------------
ELEMENTAL FUNCTION char_to_upper(c) result(t)
CHARACTER(len=1), INTENT(in) :: c !< Input character
CHARACTER(len=1) :: t !< Uppercase version
INTEGER(i4), PARAMETER :: wp=iachar('A')-iachar('a')
INTEGER(i4), PARAMETER :: BA=iachar('a')
INTEGER(i4), PARAMETER :: BZ=iachar('z')
INTEGER(i4) :: k
k = ichar(c) 
IF(k>=BA.and.k<=BZ)k = k + wp 
t = char(k)
END FUNCTION char_to_upper
!------------------------------------------------------------------------------
!> Converts a string to all uppercase characters
!------------------------------------------------------------------------------
SUBROUTINE string_to_upper(c)
CHARACTER(len=*), INTENT(inout) :: c !< Input/output string
INTEGER(i4) :: i
DO i=1,LEN(c)
  c(i:i)=char_to_upper(c(i:i))
END DO
END SUBROUTINE string_to_upper
!------------------------------------------------------------------------------
!> Start or reset timer
!------------------------------------------------------------------------------
SUBROUTINE oft_timer_start(self)
CLASS(oft_timer), INTENT(inout) :: self !< Calling timer class
self%count=oft_time_i8()
END SUBROUTINE oft_timer_start
!------------------------------------------------------------------------------
!> Set elapsed time since last tick/tock
!------------------------------------------------------------------------------
FUNCTION oft_timer_elapsed(self) RESULT(time)
CLASS(oft_timer), INTENT(inout) :: self !< Calling timer class
REAL(r8) :: time !< Time since last tick/tock
INTEGER(i8) :: countnew,crate,cmax,dt
CALL system_clock(countnew,crate,cmax)
dt=countnew-self%count
IF(dt<0)dt=dt+cmax
time=dt/REAL(crate,8)
self%count=countnew
END FUNCTION oft_timer_elapsed
!------------------------------------------------------------------------------
!> Get elapsed time since last tick/tock in integer counts
!!
!! @param[in,out] self Calling timer class
!! @return Number of integer counts since last tick/tock
!------------------------------------------------------------------------------
FUNCTION oft_timer_intelapsed(self) result(dt)
CLASS(oft_timer), INTENT(inout) :: self !< Calling timer class
INTEGER(i8) :: dt !< Number of integer counts since last tick/tock
INTEGER(i8) :: countnew,crate,cmax
CALL system_clock(countnew,crate,cmax)
dt=countnew-self%count
IF(dt<0)dt=dt+cmax
self%count=countnew
END FUNCTION oft_timer_intelapsed
!------------------------------------------------------------------------------
!> Check if time since last tick/tock exceeds a limit
!------------------------------------------------------------------------------
FUNCTION oft_timer_timeout(self,timeout) result(test)
CLASS(oft_timer), INTENT(inout) :: self !< Calling timer class
REAL(r8), INTENT(in) :: timeout !< Length of timeout (seconds)
INTEGER(i8) :: countnew,crate,cmax,dt
REAL(r8) :: time
LOGICAL :: test
CALL system_clock(countnew,crate,cmax)
dt=countnew-self%count
IF(dt<0)dt=dt+cmax
time=dt/REAL(crate,8)
test=(time>timeout)
END FUNCTION oft_timer_timeout
!------------------------------------------------------------------------------
!> Get current system time in integer counts
!------------------------------------------------------------------------------
FUNCTION oft_time_i8() RESULT(time)
INTEGER(i8) :: time !< System time in integer counts
INTEGER(i8) :: crate,cmax
CALL system_clock(time,crate,cmax)
END FUNCTION oft_time_i8
!------------------------------------------------------------------------------
!> Get difference between timestamps, including wrapping
!------------------------------------------------------------------------------
FUNCTION oft_time_diff(timein) RESULT(dt)
INTEGER(i8), intent(in) :: timein !< Previous time in integer counts
INTEGER(i8) :: timenew,crate,cmax,dt
CALL system_clock(timenew,crate,cmax)
dt=timenew-timein
IF(dt<0)dt=dt+cmax
END FUNCTION oft_time_diff
!------------------------------------------------------------------------------
!> Skip comment lines in open file
!------------------------------------------------------------------------------
function skip_comment_lines(io_unit) result(status)
integer(i4), intent(in) :: io_unit !< I/O unit to advance
integer(i4) :: status !< IOSTAT from last read or -1 if io_unit is not open
integer(i4) :: i
logical :: io_open
CHARACTER(LEN=1) :: test_char
INQUIRE(unit=io_unit,opened=io_open)
IF(io_open)THEN
  status=0
  DO WHILE(status==0)
    READ(io_unit,"(A1)",IOSTAT=status,ERR=800)test_char
    BACKSPACE(io_unit,IOSTAT=status,ERR=800)
    IF(test_char=="#")THEN
      READ(io_unit,*,IOSTAT=status,ERR=800)
    ELSE
      RETURN
    END IF
  END DO
ELSE
  status=-1
END IF
800 RETURN
end function skip_comment_lines
!------------------------------------------------------------------------------
!> integer(i4) implementation of \ref oft_local::get_inverse_map
!------------------------------------------------------------------------------
subroutine get_inverse_map_i4(map,n1,imap,n2)
integer(i4), intent(inout) :: map(n1) !< Forward map [n1]
integer(i4), intent(inout) :: imap(n2) !< Inverse map [n2]
integer(i4), intent(in) :: n1 !< Length of forward map
integer(i4), intent(in) :: n2 !< Length of inverse map (n2<=MAX(map))
integer(i4) :: i
! if(size(imap)/=n2)call oft_abort('Invalid array size','get_inverse_map_i4',__FILE__)
imap=0
!$omp parallel do
do i=1,n1
  imap(map(i))=i
end do
end subroutine get_inverse_map_i4
!------------------------------------------------------------------------------
!> integer(i8) implementation of \ref oft_local::get_inverse_map
!------------------------------------------------------------------------------
subroutine get_inverse_map_i8(map,n1,imap,n2)
integer(i8), intent(inout) :: map(n1) !< Forward map [n1]
integer(i4), intent(inout) :: imap(n2) !< Inverse map [n2]
integer(i4), intent(in) :: n1 !< Length of forward map
integer(i8), intent(in) :: n2 !< Length of inverse map (n2<=MAX(map))
integer(i4) :: i
! if(size(imap)/=n2)call oft_abort('Invalid array size','get_inverse_map_i8',__FILE__)
imap=0
!$omp parallel do
do i=1,n1
  imap(map(i))=i
end do
end subroutine get_inverse_map_i8
!------------------------------------------------------------------------------
!> Compute the cross product of two 3 dimensional vectors
!------------------------------------------------------------------------------
PURE FUNCTION cross_product(a,b) RESULT(c)
REAL(r8), INTENT(in) :: a(3) !< Vector 1 [3]
REAL(r8), INTENT(in) :: b(3) !< Vector 2 [3]
REAL(r8) :: c(3) !< \f$ a \times b \f$ [3]
INTEGER(i4), PARAMETER :: i2(3)=[2,3,1],i3(3)=[3,1,2]
c=a(i2)*b(i3)-a(i3)*b(i2)
END FUNCTION cross_product
!------------------------------------------------------------------------------
!> Compute the 2-norm of an array
!------------------------------------------------------------------------------
PURE FUNCTION magnitude(a) RESULT(c)
REAL(r8), INTENT(in) :: a(:) !< Array
REAL(r8) :: c !< \f$ \sum_i a^2_i \f$
c=SQRT(SUM(a**2))
END FUNCTION magnitude
!------------------------------------------------------------------------------
!> Compute the 2-norm of an array
!------------------------------------------------------------------------------
PURE FUNCTION time_to_string(a) RESULT(c)
REAL(r8), INTENT(in) :: a !< Array
CHARACTER(LEN=13) :: c !< \f$ \sum_i a^2_i \f$
INTEGER(4) :: hours,minutes,seconds
hours = FLOOR(a/3600.d0)
minutes = FLOOR((a-hours*3600.d0)/60.d0)
seconds = FLOOR((a-hours*3600.d0-minutes*60.d0))
IF(hours>0)THEN
  WRITE(c,'(I4,A,I2,A,I2,A)')hours,'h ',minutes,'m ',seconds,'s'
ELSE IF(minutes>0)THEN
  WRITE(c,'(I2,A,I2,A,6X)')minutes,'m ',seconds,'s'
ELSE
  WRITE(c,'(I2,A,10X)')seconds,'s'
END IF
END FUNCTION time_to_string
#ifdef HAVE_XML
!------------------------------------------------------------------------------
!> Get child element with a specific name within a given XML node
!------------------------------------------------------------------------------
subroutine xml_get_element_single(parent,name,element,error_flag,index)
USE fox_dom, ONLY: nodelist, item, getLength, getChildNodes, getNodeName, DOMException, &
  getExceptionCode
TYPE(xml_node), POINTER, INTENT(in) :: parent !< Parent element
CHARACTER(LEN=*), INTENT(in) :: name !< Name of child element to find
TYPE(xml_node), POINTER, INTENT(inout) :: element !< Found element
INTEGER(i4), INTENT(out) :: error_flag !< Error flag (0 if successful)
INTEGER(i4), OPTIONAL, INTENT(in) :: index !< Optional index, defaults to first matching element
INTEGER(i4) :: i,req_index,nchildren,nelements
TYPE(xml_node), POINTER :: tmp_element
TYPE(nodelist), POINTER :: tmp_list
TYPE(DOMException) :: xml_ex
NULLIFY(element)
IF(.NOT.ASSOCIATED(parent))THEN
  error_flag=1
  RETURN
END IF
req_index=1
IF(PRESENT(index))req_index=index
IF(req_index<=0)THEN
  error_flag=3
  RETURN
END IF
tmp_list=>getChildNodes(parent,ex=xml_ex)
IF(getExceptionCode(xml_ex)/=0)THEN
  error_flag=2
  RETURN
END IF
nchildren=getLength(tmp_list)
nelements=0
DO i=1,nchildren
  tmp_element=>item(tmp_list,i-1)
  IF(getNodeName(tmp_element)==TRIM(name))THEN
    nelements=nelements+1
    IF(nelements==req_index)THEN
      element=>tmp_element
      EXIT
    END IF
  END IF
END DO
IF(nelements==0)THEN
  error_flag=4
  RETURN
END IF
IF(req_index>nelements)THEN
  error_flag=-nelements
  RETURN
END IF
error_flag=0
end subroutine xml_get_element_single
!------------------------------------------------------------------------------
!> Get all child elements with a specific name within a given XML node
!------------------------------------------------------------------------------
subroutine xml_get_element_list(parent,name,elements,error_flag)
USE fox_dom, ONLY: nodelist, item, getLength, getChildNodes, getNodeName, DOMException, &
  getExceptionCode
TYPE(xml_node), POINTER, INTENT(in) :: parent !< Parent element
CHARACTER(LEN=*), INTENT(in) :: name !< Name of child element to find
TYPE(xml_nodelist), INTENT(inout) :: elements !< Found elements
INTEGER(i4), INTENT(out) :: error_flag !< Error flag (0 if successful)
INTEGER(i4) :: i,nchildren
TYPE(xml_node), POINTER :: tmp_element
TYPE(nodelist), POINTER :: tmp_list
TYPE(DOMException) :: xml_ex
IF(ASSOCIATED(elements%nodes))DEALLOCATE(elements%nodes)
elements%n=0
IF(.NOT.ASSOCIATED(parent))THEN
  error_flag=1
  RETURN
END IF
tmp_list=>getChildNodes(parent,ex=xml_ex)
IF(getExceptionCode(xml_ex)/=0)THEN
  error_flag=2
  RETURN
END IF
nchildren=getLength(tmp_list)
DO i=1,nchildren
  tmp_element=>item(tmp_list,i-1)
  IF(getNodeName(tmp_element)==TRIM(name))THEN
    elements%n=elements%n+1
  END IF
END DO
IF(elements%n==0)THEN
  error_flag=4
  RETURN
END IF
ALLOCATE(elements%nodes(elements%n))
elements%n=0
DO i=1,nchildren
  tmp_element=>item(tmp_list,i-1)
  IF(getNodeName(tmp_element)==TRIM(name))THEN
    elements%n=elements%n+1
    elements%nodes(elements%n)%this=>tmp_element
  END IF
END DO
error_flag=0
end subroutine xml_get_element_list
#endif
END MODULE oft_local
