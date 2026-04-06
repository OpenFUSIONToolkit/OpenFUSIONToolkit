!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_xml.F90
!
!> @brief Libxml2 DOM interface for OFT
!!
!! Provides a limited interface between Libxml2 and Fortran using
!! iso_c_binding for DOM-style XML access, along with helper routines for
!! parsing string content into typed Fortran values. Data stored in XML
!! nodes/attributes can be either space-delimited or comma-delimited with
!! newlines to denote 2D arrays.
!!
!! @authors Chris Hansen
!! @date April 2026
!! @ingroup doxy_oft_base
!---------------------------------------------------------------------------------
MODULE oft_xml
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_ptr, c_char, c_null_char, &
  c_null_ptr, c_f_pointer, c_associated, c_double, c_loc
USE oft_local, ONLY: i4, r8, string_to_upper
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Libxml2 node wrapper
!------------------------------------------------------------------------------
TYPE :: xml_node
  TYPE(c_ptr) :: obj = c_null_ptr !< Opaque pointer to xmlNode C object
CONTAINS
  PROCEDURE :: associated => xml_node_associated
END TYPE xml_node
!------------------------------------------------------------------------------
!> List of XML nodes
!------------------------------------------------------------------------------
TYPE :: xml_nodelist
  INTEGER(i4) :: n = 0 !< Number of nodes in the list
  TYPE(xml_node), POINTER, DIMENSION(:) :: nodes => NULL() !< Pointer to nodes
END TYPE xml_nodelist
!------------------------------------------------------------------------------
!> XML document wrapper
!------------------------------------------------------------------------------
TYPE :: xml_doc
  TYPE(c_ptr) :: doc = c_null_ptr !< Opaque pointer to xmlDoc
  TYPE(xml_node), POINTER :: root => NULL() !< Pointer to root node
END TYPE xml_doc
INTERFACE
!------------------------------------------------------------------------------
!> Parse an XML file from a given file path, returning a pointer to the corresponding xmlDoc
!------------------------------------------------------------------------------
  FUNCTION oft_xml_parse_file_c(filepath,doc_ptr) BIND(C, NAME="oft_xml_parse_file") &
      RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: filepath(*) !< Path to XML file
    TYPE(c_ptr), INTENT(out) :: doc_ptr !< Pointer to parsed xmlDoc object
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_parse_file_c
!------------------------------------------------------------------------------
!> Retrieve a pointer to the root element of an xmlDoc object
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_root_c(doc_ptr,root_ptr) BIND(C, NAME="oft_xml_get_root") &
      RESULT(ierr)
    IMPORT c_int, c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: doc_ptr !< Pointer to xmlDoc
    TYPE(c_ptr), INTENT(out) :: root_ptr !< Pointer to root xmlNode
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_root_c
!------------------------------------------------------------------------------
!> Retrieve a pointer to the i-th xml node with a given name inside a parent node
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_element_c(parent_ptr,name,index,element_ptr) &
      BIND(C, NAME="oft_xml_get_element") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: parent_ptr !< Pointer to parent xmlNode
    CHARACTER(KIND=c_char), INTENT(in) :: name(*) !< Element name to search for
    INTEGER(c_int), VALUE, INTENT(in) :: index !< Desired index among matching children
    TYPE(c_ptr), INTENT(out) :: element_ptr !< Pointer to matching xmlNode
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_element_c
!------------------------------------------------------------------------------
!> Retrieve pointers to all xml nodes with a given name inside a parent node
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_elements_c(parent_ptr,name,n,elements_ptr) &
      BIND(C, NAME="oft_xml_get_elements") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: parent_ptr !< Pointer to parent xmlNode
    CHARACTER(KIND=c_char), INTENT(in) :: name(*) !< Element name to search for
    INTEGER(c_int), INTENT(out) :: n !< Number of matching children
    TYPE(c_ptr), INTENT(out) :: elements_ptr !< Pointer to array of matching nodes
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_elements_c
!------------------------------------------------------------------------------
!> Extract content from a given xml node as a string
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_content_c(node_ptr,content,content_len) &
      BIND(C, NAME="oft_xml_get_content") RESULT(ierr)
    IMPORT c_int, c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    TYPE(c_ptr), INTENT(out) :: content !< Output character buffer
    INTEGER(c_int), INTENT(out) :: content_len !< Length of content
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_content_c
!------------------------------------------------------------------------------
!> Extract content from an attribute on a given xml node as a string
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_attribute_c(node_ptr,attr_name,content,content_len) &
      BIND(C, NAME="oft_xml_get_attribute") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    TYPE(c_ptr), VALUE, INTENT(in) :: attr_name !< Attribute name
    TYPE(c_ptr), INTENT(out) :: content !< Output character buffer
    INTEGER(c_int), INTENT(out) :: content_len !< Length of content
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_attribute_c
!------------------------------------------------------------------------------
!> Test if a given xml node has a given attribute
!------------------------------------------------------------------------------
  FUNCTION oft_xml_has_attribute_c(node_ptr,attr_name) &
      BIND(C, NAME="oft_xml_has_attribute") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    TYPE(c_ptr), VALUE, INTENT(in) :: attr_name !< Attribute name
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_has_attribute_c
!------------------------------------------------------------------------------
!> Free c objects, for example arrays returned by \ref oft_xml_get_elements_c or
!! strings returned by \ref oft_xml_get_content_c
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_ptr_c(gen_ptr) BIND(C, NAME="oft_xml_free_ptr")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: gen_ptr !< Pointer to array to free
  END SUBROUTINE oft_xml_free_ptr_c
!------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load_file_c
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_doc_c(doc_ptr) BIND(C, NAME="oft_xml_free_doc")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: doc_ptr !< Pointer to xmlDoc to free
  END SUBROUTINE oft_xml_free_doc_c
END INTERFACE
!------------------------------------------------------------------------------
!> Get one or more elements by name from a parent node
!------------------------------------------------------------------------------
INTERFACE xml_get_element
  MODULE PROCEDURE xml_get_element_single
  MODULE PROCEDURE xml_get_element_list
END INTERFACE xml_get_element
!------------------------------------------------------------------------------
!> Extract content of a given type from an XML node
!------------------------------------------------------------------------------
INTERFACE xml_read_content
  MODULE PROCEDURE xml_extractDataContent_string
  MODULE PROCEDURE xml_extractDataContent_double
  MODULE PROCEDURE xml_extractDataContent_double1D
  MODULE PROCEDURE xml_extractDataContent_double2D
  MODULE PROCEDURE xml_extractDataContent_int
  MODULE PROCEDURE xml_extractDataContent_int1D
  MODULE PROCEDURE xml_extractDataContent_int2D
  MODULE PROCEDURE xml_extractDataContent_logical
  MODULE PROCEDURE xml_extractDataContent_logical1D
  MODULE PROCEDURE xml_extractDataContent_logical2D
END INTERFACE xml_read_content
!------------------------------------------------------------------------------
!> Extract content of a given type from an attribute on an XML node
!------------------------------------------------------------------------------
INTERFACE xml_read_attribute
  MODULE PROCEDURE xml_extractDataAttribute_string
  MODULE PROCEDURE xml_extractDataAttribute_double
  MODULE PROCEDURE xml_extractDataAttribute_double1D
  MODULE PROCEDURE xml_extractDataAttribute_double2D
  MODULE PROCEDURE xml_extractDataAttribute_int
  MODULE PROCEDURE xml_extractDataAttribute_int1D
  MODULE PROCEDURE xml_extractDataAttribute_int2D
  MODULE PROCEDURE xml_extractDataAttribute_logical
  MODULE PROCEDURE xml_extractDataAttribute_logical1D
  MODULE PROCEDURE xml_extractDataAttribute_logical2D
END INTERFACE xml_read_attribute
!---Exception handling
TYPE :: xml_exception
  CHARACTER(LEN=:), ALLOCATABLE :: msg !< Exception message
END TYPE xml_exception
INTEGER(i4) :: xml_exception_depth = 0 !< Exception stack depth (0 if no exceptions)
TYPE(xml_exception) :: xml_exceptions(10) !< Stack of excepton messages
!$omp threadprivate(xml_exception_depth,xml_exceptions)
CONTAINS
!---------------------------------------------------------------------------------
!> Check if an `xml_node` is associated with a C object
!---------------------------------------------------------------------------------
FUNCTION xml_node_associated(self) RESULT(is_assoc)
CLASS(xml_node), INTENT(in) :: self !< XML node
LOGICAL :: is_assoc
is_assoc = c_associated(self%obj)
END FUNCTION xml_node_associated
!---------------------------------------------------------------------------------
!> Add exception to the exception stack
!---------------------------------------------------------------------------------
SUBROUTINE xml_set_exception(message)
CHARACTER(LEN=*), INTENT(in) :: message !< Exception message
xml_exception_depth=xml_exception_depth+1
xml_exceptions(xml_exception_depth)%msg=TRIM(ADJUSTL(message))
END SUBROUTINE xml_set_exception
!---------------------------------------------------------------------------------
!> Clear messages from the exception stack
!---------------------------------------------------------------------------------
SUBROUTINE xml_clear_exceptions()
INTEGER(i4) :: i
DO i=1,xml_exception_depth
  xml_exceptions(i)%msg=''
END DO
xml_exception_depth=0
END SUBROUTINE xml_clear_exceptions
!---------------------------------------------------------------------------------
!> Parse an XML file
!---------------------------------------------------------------------------------
SUBROUTINE xml_parsefile(filepath,doc,ierr)
CHARACTER(LEN=*), INTENT(in) :: filepath !< Path to XML file
TYPE(xml_doc), INTENT(out) :: doc !< XML document object
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_filepath(LEN_TRIM(filepath)+1)
INTEGER(i4) :: i
doc%doc=c_null_ptr
DO i=1,LEN_TRIM(filepath)
  c_filepath(i) = filepath(i:i)
END DO
c_filepath(LEN_TRIM(filepath)+1) = c_null_char
ierr = INT(oft_xml_parse_file_c(c_filepath, doc%doc), i4)
IF(ierr/=0)THEN
    doc%doc=c_null_ptr
    doc%root=>NULL()
  RETURN
END IF
ALLOCATE(doc%root)
ierr = INT(oft_xml_get_root_c(doc%doc, doc%root%obj), i4)
IF(ierr/=0)THEN
    CALL oft_xml_free_doc_c(doc%doc)
    doc%doc=c_null_ptr
    doc%root=>NULL()
  RETURN
END IF
END SUBROUTINE xml_parsefile
!---------------------------------------------------------------------------------
!> Retrieve the i-th xml node with a given name contained within a parent node
!---------------------------------------------------------------------------------
SUBROUTINE xml_get_element_single(parent,name,element,error_flag,index)
TYPE(xml_node), INTENT(in) :: parent !< Parent element
CHARACTER(LEN=*), INTENT(in) :: name !< Element name to find
TYPE(xml_node), INTENT(out) :: element !< Found element
INTEGER(i4), INTENT(out) :: error_flag !< Error flag (0 if successful)
INTEGER(i4), OPTIONAL, INTENT(in) :: index !< Optional index, defaults to first matching element
CHARACTER(KIND=c_char) :: c_name(LEN_TRIM(name)+1)
INTEGER(c_int) :: req_index
INTEGER(i4) :: i
element%obj = c_null_ptr
DO i=1,LEN_TRIM(name)
  c_name(i)=name(i:i)
END DO
c_name(LEN_TRIM(name)+1)=c_null_char
req_index=1
IF(PRESENT(index))req_index=INT(index,c_int)
error_flag = INT(oft_xml_get_element_c(parent%obj,c_name,req_index,element%obj), i4)
END SUBROUTINE xml_get_element_single
!---------------------------------------------------------------------------------
!> Retrieve all xml nodes with a given name contained within a parent node
!---------------------------------------------------------------------------------
SUBROUTINE xml_get_element_list(parent,name,elements,error_flag)
TYPE(xml_node), INTENT(in) :: parent !< Parent element
CHARACTER(LEN=*), INTENT(in) :: name !< Element name to find
TYPE(xml_nodelist), INTENT(inout) :: elements !< Found elements
INTEGER(i4), INTENT(out) :: error_flag !< Error flag (0 if successful)
TYPE(c_ptr) :: elements_c
TYPE(c_ptr), POINTER :: elements_list(:)
CHARACTER(KIND=c_char) :: c_name(LEN_TRIM(name)+1)
INTEGER(c_int) :: n_c
INTEGER(i4) :: i
elements_c=c_null_ptr
elements%n=0
IF(ASSOCIATED(elements%nodes))DEALLOCATE(elements%nodes)
DO i=1,LEN_TRIM(name)
  c_name(i)=name(i:i)
END DO
c_name(LEN_TRIM(name)+1)=c_null_char
error_flag=INT(oft_xml_get_elements_c(parent%obj,c_name,n_c,elements_c),i4)
IF(error_flag/=0)RETURN
elements%n=INT(n_c,i4)
CALL c_f_pointer(elements_c,elements_list,[elements%n])
ALLOCATE(elements%nodes(elements%n))
DO i=1,elements%n
  elements%nodes(i)%obj=elements_list(i)
END DO
CALL oft_xml_free_ptr_c(elements_c)
END SUBROUTINE xml_get_element_list
!---------------------------------------------------------------------------------
!> Extract content from a given xml node into a Fortran string
!---------------------------------------------------------------------------------
SUBROUTINE xml_get_content(node,content,ierr)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: content !< Output string
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
TYPE(c_ptr) :: buffer_ptr
INTEGER(c_int) :: content_len
INTEGER(i4) :: i
CHARACTER(KIND=c_char), POINTER :: buffer(:)
CHARACTER(LEN=:), ALLOCATABLE :: content_tmp
ierr=INT(oft_xml_get_content_c(node%obj,buffer_ptr,content_len),i4)
IF(ierr/=0)RETURN
IF(content_len<=0)THEN
  content=''
  RETURN
END IF
CALL c_f_pointer(buffer_ptr, buffer, [content_len])
ALLOCATE(character(LEN=content_len-1)::content_tmp)
DO i=1,content_len-1
  IF(buffer(i)==c_null_char)EXIT
  content_tmp(i:i)=buffer(i)
END DO
CALL oft_xml_free_ptr_c(buffer_ptr)
content=TRIM(ADJUSTL(content_tmp))
DEALLOCATE(content_tmp)
END SUBROUTINE xml_get_content
!---------------------------------------------------------------------------------
!> Extract content from a specified attribute on a given xml node into a Fortran string
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_attr(node,attr_name,content,ierr)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: content !< Output string
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char), TARGET :: c_name(LEN_TRIM(attr_name)+1)
TYPE(c_ptr) :: buffer_ptr
INTEGER(c_int) :: content_len
INTEGER(i4) :: i
CHARACTER(KIND=c_char), POINTER :: buffer(:)
CHARACTER(LEN=:), ALLOCATABLE :: content_tmp
DO i=1,LEN_TRIM(attr_name)
  c_name(i)=attr_name(i:i)
END DO
c_name(LEN_TRIM(attr_name)+1)=c_null_char
ierr=INT(oft_xml_get_attribute_c(node%obj,C_LOC(c_name),buffer_ptr,content_len),i4)
IF(ierr/=0)RETURN
CALL c_f_pointer(buffer_ptr, buffer, [content_len])
ALLOCATE(character(LEN=content_len-1)::content_tmp)
DO i=1,content_len-1
  IF(buffer(i)==c_null_char)EXIT
  content_tmp(i:i)=buffer(i)
END DO
CALL oft_xml_free_ptr_c(buffer_ptr)
content=TRIM(ADJUSTL(content_tmp))
DEALLOCATE(content_tmp)
END SUBROUTINE oft_xml_get_attr
!---------------------------------------------------------------------------------
!> Test if a given xml node has a specified attribute
!---------------------------------------------------------------------------------
FUNCTION xml_hasAttribute(node,attr_name) RESULT(has_attr)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
CHARACTER(KIND=c_char), TARGET :: c_name(LEN_TRIM(attr_name)+1)
INTEGER(i4) :: i,ierr
LOGICAL :: has_attr
has_attr=.FALSE.
DO i=1,LEN_TRIM(attr_name)
  c_name(i)=attr_name(i:i)
END DO
c_name(LEN_TRIM(attr_name)+1)=c_null_char
ierr=INT(oft_xml_has_attribute_c(node%obj,C_LOC(c_name)),i4)
IF(ierr/=0)RETURN
has_attr=.TRUE.
END FUNCTION xml_hasAttribute
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_free(doc_ptr)
TYPE(c_ptr), INTENT(inout) :: doc_ptr !< C pointer to xmlDoc to free
IF (c_associated(doc_ptr)) THEN
  CALL oft_xml_free_doc_c(doc_ptr)
  doc_ptr=c_null_ptr
END IF
END SUBROUTINE oft_xml_free
!---------------------------------------------------------------------------------
!> Normalize a string by removing leading/trailing whitespace and blank/comment lines
!---------------------------------------------------------------------------------
SUBROUTINE normalize_string(content_in,content_out)
CHARACTER(LEN=*), INTENT(in) :: content_in !< Input string
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: content_out !< Output string
INTEGER(i4) :: i,j,k,ilast,nchars
CHARACTER(*), PARAMETER :: NL = new_line('A')
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
nchars=LEN(content_in)
content_out=content_in
!---Remove leading and trailing whitespace
DO i=1,nchars
  IF((content_in(i:i)/=' ').AND.(content_in(i:i)/=NL))EXIT
END DO
DO j=nchars,1,-1
  IF((content_in(j:j)/=' ').AND.(content_in(j:j)/=NL))EXIT
END DO
nchars=j
content_out=content_in(i:nchars)
!---Remove blank lines and comment lines, and trim whitespace from remaining lines
ilast=i
k=1
DO j=i,nchars
  IF(content_in(j:j)==NL)THEN
    line_strip = TRIM(ADJUSTL(content_in(ilast:j-1)))
    IF(LEN(line_strip)==0)THEN
      ilast=j+1
      CYCLE
    END IF
    IF(line_strip(1:1)=='#')THEN
      CYCLE
    ELSE IF(line_strip(1:1)==',')THEN
      line_strip = TRIM(ADJUSTL(line_strip(2:)))
    END IF
    IF(line_strip(LEN(line_strip):LEN(line_strip))==',')THEN
      line_strip = line_strip(:LEN(line_strip)-1)
    END IF
    content_out(k:k+LEN(line_strip))=line_strip//NL
    k=k+LEN(line_strip)+1
    ilast=j+1
  END IF
END DO
!---Add final line if needed
line_strip = TRIM(ADJUSTL(content_in(ilast:nchars)))
! WRITE(*,*)'Remaining: "',line_strip,'"'
IF(LEN(line_strip)>0)THEN
  IF(line_strip(1:1)/='#')THEN
    IF(line_strip(1:1)==',')THEN
      line_strip = TRIM(ADJUSTL(line_strip(2:)))
    END IF
    IF(line_strip(LEN(line_strip):LEN(line_strip))==',')THEN
      line_strip = line_strip(:LEN(line_strip)-1)
    END IF
    content_out(k:k+LEN(line_strip)-1)=line_strip
    k=k+LEN(line_strip)
  END IF
  content_out=content_out(:k-1)
ELSE
  content_out=content_out(:k-2)
END IF
END SUBROUTINE normalize_string
!---------------------------------------------------------------------------------
!> Find token boundaries in a string containing comma-delimited and/or
!! newline-delimited values to be parsed into typed Fortran values
!---------------------------------------------------------------------------------
SUBROUTINE tokenize_string(content,tokens,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Input string
INTEGER(i4), ALLOCATABLE, INTENT(out) :: tokens(:,:,:) !< String tokens as (start,end) indices for each column and row
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
LOGICAL :: is_delim
INTEGER(i4) :: i,j,nchars,nrows,ncols,istart,iend,nr_loc,iend2,offset
REAL(r8) :: val_tmp
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
character(len=256) :: msg
character(*), parameter :: NL = new_line('A')
ierr=0
!---Count number of lines/rows
nchars=LEN(content)
nrows=1
DO i=1,nchars
  IF(content(i:i)==NL)THEN
    nrows=nrows+1
  END IF
END DO
!---Convert from space-delimited to CSV if needed
IF(INDEX(content,',')==0)THEN
  is_delim=((content(1:1)==' ').OR.(content(1:1)==NL))
  DO i=2,nchars
    IF(content(i:i)==' ')THEN
      IF(.NOT.is_delim)THEN
        content(i:i)=','
        is_delim=.TRUE.
      END IF
    ELSE
      is_delim=.FALSE.
    END IF
  END DO
END IF
! WRITE(*,*)'Parsing: ncols=',ncols
!---Count number of columns
istart=1
iend=INDEX(content(istart:nchars),NL)
IF(iend==0)THEN
  iend=nchars+1
ELSE
  iend=istart+iend-1
END IF
ncols=1
DO j=istart,iend-1
  IF(content(j:j)==',')THEN
    ncols=ncols+1
  END IF
END DO
! WRITE(*,*)'Parsing: ncols=',ncols-offset,' nrows=',nrows
!---Generate actual tokens
ALLOCATE(tokens(2,ncols,nrows))
istart=1
outer_read: DO i=1,nrows
  iend=INDEX(content(istart:nchars),NL)
  IF(iend==0)THEN
    iend=nchars+1
  ELSE
    iend=istart+iend-1
  END IF
  DO j=1,ncols
    iend2=INDEX(content(istart:iend-1),',')-1
    IF(iend2<0)THEN
      IF(j<ncols)THEN
        ierr=-10000
        CALL xml_set_exception('Inconsistent number of columns when tokenizing comma/newline-delimited string')
        EXIT outer_read
      END IF
      iend2=MIN(iend-istart,nchars-istart+1)
    END IF
    ! WRITE(*,*)i,j,'"',content(istart:istart+iend2-1),'"',istart,istart+iend2-1,nchars
    tokens(:,j,i)=[istart,istart+iend2-1]
    istart=istart+iend2+1
  END DO
  istart=iend+1
END DO outer_read
! WRITE(*,*)'Parsing error code: ',ierr
IF(ierr/=0)THEN
  DEALLOCATE(tokens)
  RETURN
END IF
END SUBROUTINE tokenize_string
!---------------------------------------------------------------------------------
!> Parse string content into a 1D array of double precision values with shape information
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_doubles(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Input string
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,nrows,ncols
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
REAL(r8) :: val_tmp
CHARACTER(LEN=8) :: col_char,row_char
character(len=256) :: msg
CHARACTER(LEN=:), ALLOCATABLE :: content_tmp

CALL normalize_string(content,content_tmp)
CALL tokenize_string(content_tmp,tokens,ierr)
IF(ierr/=0)THEN
  ierr=-1
  DEALLOCATE(content_tmp)
  RETURN
END IF

nrows=SIZE(tokens,3)
ncols=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[ncols,nrows]
ALLOCATE(output(ncols*nrows))
outer_read: DO i=1,nrows
  DO j=1,ncols
    READ(content_tmp(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr,IOMSG=msg)val_tmp
    IF(ierr/=0)THEN
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing double value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      CALL xml_set_exception(msg)
      EXIT outer_read
    END IF
    output((i-1)*ncols+j)=val_tmp
  END DO
END DO outer_read
DEALLOCATE(tokens,content_tmp)
IF(ierr/=0)THEN
  ierr=-2
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_doubles
!---------------------------------------------------------------------------------
!> Parse string content into a 1D array of integer values with shape information
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_integers(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Input string
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,nrows,ncols
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
INTEGER(i4) :: val_tmp
CHARACTER(LEN=8) :: col_char,row_char
character(len=256) :: msg
CHARACTER(LEN=:), ALLOCATABLE :: content_tmp

CALL normalize_string(content,content_tmp)
CALL tokenize_string(content_tmp,tokens,ierr)
IF(ierr/=0)THEN
  ierr=-1
  DEALLOCATE(content_tmp)
  RETURN
END IF

nrows=SIZE(tokens,3)
ncols=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[ncols,nrows]
ALLOCATE(output(ncols*nrows))
outer_read: DO i=1,nrows
  DO j=1,ncols
    READ(content_tmp(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr,IOMSG=msg)val_tmp
    IF(ierr/=0)THEN
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing integer value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      CALL xml_set_exception(msg)
      EXIT outer_read
    END IF
    output((i-1)*ncols+j)=val_tmp
  END DO
END DO outer_read
DEALLOCATE(tokens,content_tmp)
IF(ierr/=0)THEN
  ierr=-2
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_integers
!---------------------------------------------------------------------------------
!> Parse string content into a 1D array of logical values with shape information
!!
!! Supported values are:
!! - True: "true", ".true.", "1", "t" (case-insensitive)
!! - False: "false", ".false.", "0", "f" (case-insensitive)
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_logicals(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Input string
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,nrows,ncols
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
CHARACTER(LEN=8) :: col_char,row_char
CHARACTER(LEN=:), ALLOCATABLE :: token_strip
CHARACTER(LEN=:), ALLOCATABLE :: content_tmp

CALL normalize_string(content,content_tmp)
CALL tokenize_string(content_tmp,tokens,ierr)
IF(ierr/=0)THEN
  ierr=-1
  DEALLOCATE(content_tmp)
  RETURN
END IF
ierr=0

nrows=SIZE(tokens,3)
ncols=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[ncols,nrows]
ALLOCATE(output(ncols*nrows))
outer_read: DO i=1,nrows
  DO j=1,ncols
    token_strip=TRIM(ADJUSTL(content_tmp(tokens(1,j,i):tokens(2,j,i))))
    CALL string_to_upper(token_strip)
    IF((token_strip=='TRUE').OR.(token_strip=='.TRUE.').OR.(token_strip=='1').OR.(token_strip=='T'))THEN
      output((i-1)*ncols+j)=.TRUE.
    ELSE IF((token_strip=='FALSE').OR.(token_strip=='.FALSE.').OR.(token_strip=='0').OR.(token_strip=='F'))THEN
      output((i-1)*ncols+j)=.FALSE.
    ELSE
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing logical value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      ierr=-2
      EXIT outer_read
    END IF
  END DO
END DO outer_read
DEALLOCATE(tokens,content_tmp)
IF(ierr/=0)THEN
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_logicals
!---------------------------------------------------------------------------------
!> Get content from an XML node as a string
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_string(node,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CALL xml_get_content(node,output,ierr)
IF(ierr/=0)CALL xml_set_exception('Error in child of xml_extractDataContent_string')
IF(PRESENT(iostat))iostat=ierr
END SUBROUTINE xml_extractDataContent_string
!---------------------------------------------------------------------------------
!> Parse XML content into a 1D array of double precision values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL xml_get_content(node,content,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_double1D')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
CALL parse_string_to_doubles(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_double1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double1D
!---------------------------------------------------------------------------------
!> Parse a single double precision value from an XML node's content
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
REAL(r8), INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataContent_double1D(node,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_double')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double
!---------------------------------------------------------------------------------
!> Parse XML content into a 2D array of double precision values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
REAL(r8), POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: shape_tmp(2),ierr
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataContent_double1D(node,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_double2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double2D
!---------------------------------------------------------------------------------
!> Parse XML content into a 1D array of integer values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL xml_get_content(node,content,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_int1D')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
CALL parse_string_to_integers(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_int1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_int1D
!---------------------------------------------------------------------------------
!> Parse XML content into a single integer value
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
INTEGER(i4), INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4), POINTER :: output_tmp(:)
CALL xml_extractDataContent_int1D(node,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_int')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_int
!---------------------------------------------------------------------------------
!> Parse XML content into a 2D array of integer values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
INTEGER(i4), POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: i,j,ierr,shape_tmp(2)
INTEGER(i4), POINTER :: output_tmp(:)
CALL xml_extractDataContent_int1D(node,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_int2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_int2D
!---------------------------------------------------------------------------------
!> Parse XML content into a 1D array of logical values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL xml_get_content(node,content,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_logical1D')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
CALL parse_string_to_logicals(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_logical1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical1D
!---------------------------------------------------------------------------------
!> Parse XML content into a single logical value
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
LOGICAL, INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
LOGICAL, POINTER :: output_tmp(:)
CALL xml_extractDataContent_logical1D(node,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_logical')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical
!---------------------------------------------------------------------------------
!> Parse XML content into a 2D array of logical values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
LOGICAL, POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
LOGICAL, POINTER :: output_tmp(:)
CALL xml_extractDataContent_logical1D(node,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataContent_logical2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical2D
!---------------------------------------------------------------------------------
!> Get content from an XML attribute as a string
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_string(node,attr_name,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
CHARACTER(LEN=:), ALLOCATABLE, INTENT(inout) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CALL oft_xml_get_attr(node,attr_name,output,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_string')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_string
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 1D array of double precision values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_double1D')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
CALL parse_string_to_doubles(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_double1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_double1D
!---------------------------------------------------------------------------------
!> Parse a single double precision value from an XML attribute's content
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_double1D(node,attr_name,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_double')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_double
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 2D array of double precision values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_double1D(node,attr_name,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_double2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_double2D
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 1D array of integer values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
 CALL xml_set_exception('Error in child of xml_extractDataAttribute_int1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
CALL parse_string_to_integers(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_int1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int1D
!---------------------------------------------------------------------------------
!> Parse a single integer value from an XML attribute's content
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_int1D(node,attr_name,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_int')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 2D array of integer values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
INTEGER(i4), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_int1D(node,attr_name,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_int2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int2D
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 1D array of logical values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_logical1D')
  IF(PRESENT(iostat))iostat=-1
  RETURN
END IF
CALL parse_string_to_logicals(content,output,output_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_logical1D')
  IF(PRESENT(iostat))iostat=-2
  RETURN
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical1D
!---------------------------------------------------------------------------------
!> Parse a single logical value from an XML attribute's content
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, INTENT(out) :: output !< Output value
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
LOGICAL, POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_logical1D(node,attr_name,output_tmp,output_shape,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_logical2')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
output=output_tmp(1)
DEALLOCATE(output_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical
!---------------------------------------------------------------------------------
!> Parse XML attribute content into a 2D array of logical values
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< XML node
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(out) :: output(:,:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (ncols, nrows)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
LOGICAL, POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_logical1D(node,attr_name,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  CALL xml_set_exception('Error in child of xml_extractDataAttribute_logical2D')
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(output_shape))output_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical2D
END MODULE oft_xml
