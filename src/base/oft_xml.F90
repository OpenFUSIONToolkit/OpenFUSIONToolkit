!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_xml.F90
!
!> @brief LIBXML2 DOM interface for Fortran.
!!
!! Provides a limited interface between LIBXML2 and Fortran using
!! iso_c_binding for DOM-style XML access, along with helper routines for
!! parsing string content into typed Fortran values.
!!
!! The XML document access functions are guarded by @c HAVE_LIBXML2 and are
!! only available when the library is found at build time.  The string parsing
!! routines (\ref oft_xml_parse_logical, \ref oft_xml_parse_int,
!! \ref oft_xml_parse_real and their @c _array variants) are always available.
!!
!! The @c _array variants accept comma-and-newline-delimited strings, reading
!! data into a flat array.  A 2-component @p shape (nrows, ncols) is returned,
!! where commas delimit columns within a row and newlines delimit rows.
!!
!! Typical usage (with LIBXML2):
!! @code{.f90}
!! TYPE(c_ptr) :: doc, root, node
!! TYPE(c_ptr) :: elements_c
!! TYPE(c_ptr), POINTER :: elements(:)
!! INTEGER(i4) :: arr_shape(2), ierr
!! REAL(r8) :: vals(64)
!! CHARACTER(LEN=256) :: content
!!
!! CALL oft_xml_load("input.xml", doc, ierr)
!! CALL oft_xml_get_root(doc, root, ierr)
!! CALL oft_xml_get_element(root, "section", node, ierr)
!! CALL oft_xml_get_content(node, content, ierr)
!! CALL oft_xml_parse_real_array(content, vals, arr_shape, ierr)
!!
!! CALL oft_xml_get_all_elements(root, "item", n, elements_c, ierr)
!! CALL c_f_pointer(elements_c, elements, [n])
!! ! ... use elements(i) ...
!! CALL oft_xml_free_elements(elements_c)
!! NULLIFY(elements)
!!
!! CALL oft_xml_free(doc)
!! @endcode
!!
!! @authors Chris Hansen
!! @date 2024
!! @ingroup doxy_oft_base
!---------------------------------------------------------------------------------
MODULE oft_xml
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_ptr, c_char, c_null_char, &
  c_null_ptr, c_f_pointer, c_associated, c_double, c_loc
USE oft_local, ONLY: i4, r8, string_to_upper
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE :: xml_node
  TYPE(c_ptr) :: obj = c_null_ptr !< Opaque pointer to xmlNode
CONTAINS
  PROCEDURE :: associated => xml_node_associated
END TYPE xml_node
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE :: xml_nodelist
  INTEGER(i4) :: n = 0 !< Number of nodes in the list
  TYPE(xml_node), POINTER, DIMENSION(:) :: nodes => NULL() !< Pointer to nodes
END TYPE xml_nodelist
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE :: xml_doc
  TYPE(c_ptr) :: doc = c_null_ptr !< Opaque pointer to xmlDoc
  TYPE(xml_node), POINTER :: root => NULL() !< Pointer to root node
END TYPE xml_doc
INTERFACE
!------------------------------------------------------------------------------
!> Parse an XML file from a given file path, returning a pointer to the xmlDoc.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_parse_file_c(filepath,doc_ptr) BIND(C, NAME="oft_xml_parse_file") &
      RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: filepath(*) !< Null-terminated file path
    TYPE(c_ptr), INTENT(out) :: doc_ptr !< Pointer to parsed xmlDoc on return
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_parse_file_c
!------------------------------------------------------------------------------
!> Retrieve a pointer to the root element of an XML document.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_root_c(doc_ptr,root_ptr) BIND(C, NAME="oft_xml_get_root") &
      RESULT(ierr)
    IMPORT c_int, c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: doc_ptr !< Pointer to xmlDoc
    TYPE(c_ptr), INTENT(out) :: root_ptr !< Pointer to root xmlNode on return
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_root_c
!------------------------------------------------------------------------------
!> Retrieve a pointer to the I-th xml node with a given name inside a parent
!! node (1-based index).
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_element_c(parent_ptr,name,index,element_ptr) &
      BIND(C, NAME="oft_xml_get_element") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: parent_ptr !< Pointer to parent xmlNode
    CHARACTER(KIND=c_char), INTENT(in) :: name(*) !< Null-terminated element name
    INTEGER(c_int), VALUE, INTENT(in) :: index !< 1-based index among matching children
    TYPE(c_ptr), INTENT(out) :: element_ptr !< Pointer to matching xmlNode on return
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_element_c
!------------------------------------------------------------------------------
!> Retrieve pointers to all xml nodes with a given name inside a parent node.
!!
!! On success @p elements_ptr points to a heap-allocated array of @c c_ptr
!! values (one per match); the caller must free it with
!! \ref oft_xml_free_elements.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_elements_c(parent_ptr,name,n,elements_ptr) &
      BIND(C, NAME="oft_xml_get_elements") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: parent_ptr !< Pointer to parent xmlNode
    CHARACTER(KIND=c_char), INTENT(in) :: name(*) !< Null-terminated element name
    INTEGER(c_int), INTENT(out) :: n !< Number of matching children on return
    TYPE(c_ptr), INTENT(out) :: elements_ptr !< Pointer to allocated array of c_ptr on return
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_elements_c
!------------------------------------------------------------------------------
!> Extract the string content from a given xml node.
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
!> Extract the string value of a given attribute on a given xml node.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_attribute_c(node_ptr,attr_name,content,content_len) &
      BIND(C, NAME="oft_xml_get_attribute") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    TYPE(c_ptr), VALUE, INTENT(in) :: attr_name !< Null-terminated attribute name
    TYPE(c_ptr), INTENT(out) :: content !< Output character buffer
    INTEGER(c_int), INTENT(out) :: content_len !< Length of content
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_attribute_c
!------------------------------------------------------------------------------
!> Extract the string value of a given attribute on a given xml node.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_has_attribute_c(node_ptr,attr_name) &
      BIND(C, NAME="oft_xml_has_attribute") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    TYPE(c_ptr), VALUE, INTENT(in) :: attr_name !< Null-terminated attribute name
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_has_attribute_c
!------------------------------------------------------------------------------
!> Free the array of node pointers allocated by \ref oft_xml_get_elements_c.
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_ptr_c(gen_ptr) BIND(C, NAME="oft_xml_free_ptr")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: gen_ptr !< Pointer to array to free
  END SUBROUTINE oft_xml_free_ptr_c
!------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load_file_c.
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_doc_c(doc_ptr) BIND(C, NAME="oft_xml_free_doc")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: doc_ptr !< Pointer to xmlDoc to free
  END SUBROUTINE oft_xml_free_doc_c
END INTERFACE
!------------------------------------------------------------------------------
!> Generate inverse of sparse indexing
!------------------------------------------------------------------------------
INTERFACE xml_get_element
  MODULE PROCEDURE xml_get_element_single
  MODULE PROCEDURE xml_get_element_list
END INTERFACE xml_get_element
!------------------------------------------------------------------------------
!> Generate inverse of sparse indexing
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
!> Generate inverse of sparse indexing
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
!---
TYPE :: xml_exception
  CHARACTER(LEN=:), ALLOCATABLE :: msg
END TYPE xml_exception
INTEGER(i4) :: xml_exception_depth = 0
TYPE(xml_exception) :: xml_exceptions(10)
!$omp threadprivate(xml_exception_depth,xml_exceptions)
CONTAINS
!---------------------------------------------------------------------------------
!> Check if an xml_node is associated with a valid xmlNode pointer.
!---------------------------------------------------------------------------------
FUNCTION xml_node_associated(self) RESULT(is_assoc)
CLASS(xml_node), INTENT(in) :: self
LOGICAL :: is_assoc
is_assoc = c_associated(self%obj)
END FUNCTION xml_node_associated
!---------------------------------------------------------------------------------
!> Check if an xml_node is associated with a valid xmlNode pointer.
!---------------------------------------------------------------------------------
SUBROUTINE xml_set_exception(message)
CHARACTER(LEN=*), INTENT(in) :: message
xml_exception_depth=xml_exception_depth+1
xml_exceptions(xml_exception_depth)%msg=TRIM(ADJUSTL(message))
END SUBROUTINE xml_set_exception
!---------------------------------------------------------------------------------
!> Check if an xml_node is associated with a valid xmlNode pointer.
!---------------------------------------------------------------------------------
SUBROUTINE xml_clear_exceptions()
INTEGER(i4) :: i
DO i=1,xml_exception_depth
  xml_exceptions(i)%msg=''
END DO
xml_exception_depth=0
END SUBROUTINE xml_clear_exceptions
!---------------------------------------------------------------------------------
!> Parse an XML file from a given file path.
!---------------------------------------------------------------------------------
SUBROUTINE xml_parsefile(filepath,doc,ierr)
CHARACTER(LEN=*), INTENT(in) :: filepath !< Path to XML file
TYPE(xml_doc), INTENT(out) :: doc !< C pointer to parsed xmlDoc
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
!> Retrieve a pointer to the I-th xml node with a given name contained within
!! a parent node (1-based index, defaults to 1).
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
!> Retrieve C pointers to all xml nodes with a given name contained within a
!! parent node.
!!
!! @p elements_c is set to a raw C pointer to a heap-allocated array of
!! @c c_ptr values (one per matching child).  Use @c c_f_pointer to obtain a
!! Fortran array over this memory, and call \ref oft_xml_free_elements with
!! @p elements_c when done.
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
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE xml_get_content(node,content,ierr)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
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
!> Extract the string value of a given attribute on a given xml node.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_attr(node,attr_name,content,ierr)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
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
!> Extract the string value of a given attribute on a given xml node.
!---------------------------------------------------------------------------------
FUNCTION xml_hasAttribute(node,attr_name) RESULT(has_attr)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_free(doc_ptr)
TYPE(c_ptr), INTENT(inout) :: doc_ptr !< C pointer to xmlDoc to free
IF (c_associated(doc_ptr)) THEN
  CALL oft_xml_free_doc_c(doc_ptr)
  doc_ptr=c_null_ptr
END IF
END SUBROUTINE oft_xml_free
!---------------------------------------------------------------------------------
!> Extract the string content from a given xml node into a Fortran string.
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
  IF((content_in(i:i)/=' ').OR.(content_in(i:i)/=NL))EXIT
END DO
DO j=nchars,1,-1
  IF((content_in(j:j)/=' ').OR.(content_in(j:j)/=NL))EXIT
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
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE tokenize_string(content,tokens,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Output string
INTEGER(i4), ALLOCATABLE, INTENT(out) :: tokens(:,:,:) !< Output array
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
LOGICAL :: is_delim
INTEGER(i4) :: i,j,nchars,nrows,ncols,istart,iend,nr_loc,iend2,offset
REAL(r8) :: val_tmp
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
character(len=256) :: msg
character(*), parameter :: NL = new_line('A')
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
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_doubles(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Output string
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
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

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    READ(content_tmp(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr,IOMSG=msg)val_tmp
    IF(ierr/=0)THEN
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing double value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      CALL xml_set_exception(msg)
      EXIT outer_read
    END IF
    output((i-1)*nrows+j)=val_tmp
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
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_integers(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Output string
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
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

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    READ(content_tmp(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr)val_tmp
    IF(ierr/=0)THEN
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing integer value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      CALL xml_set_exception(msg)
      EXIT outer_read
    END IF
    output((i-1)*nrows+j)=val_tmp
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
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_logicals(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(inout) :: content !< Output string
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
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

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    token_strip=TRIM(ADJUSTL(content_tmp(tokens(1,j,i):tokens(2,j,i))))
    CALL string_to_upper(token_strip)
    IF((token_strip=='TRUE').OR.(token_strip=='.TRUE.').OR.(token_strip=='1').OR.(token_strip=='T'))THEN
      output((i-1)*nrows+j)=.TRUE.
    ELSE IF((token_strip=='FALSE').OR.(token_strip=='.FALSE.').OR.(token_strip=='0').OR.(token_strip=='F'))THEN
      output((i-1)*nrows+j)=.FALSE.
    ELSE
      WRITE(col_char,'(I8)')i
      WRITE(row_char,'(I8)')j
      CALL xml_set_exception('Error parsing logical value at column '//TRIM(ADJUSTL(col_char))//' row '//TRIM(ADJUSTL(row_char)))
      EXIT outer_read
    END IF
  END DO
END DO outer_read
DEALLOCATE(tokens,content_tmp)
IF(ierr/=0)THEN
  ierr=-2
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_logicals
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_string(node,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CALL xml_get_content(node,output,ierr)
CALL xml_set_exception('Error in child of xml_extractDataContent_string')
IF(PRESENT(iostat))iostat=ierr
END SUBROUTINE xml_extractDataContent_string
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical1D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical2D(node,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_string(node,attr_name,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical1D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical2D(node,attr_name,output,output_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(out) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
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
