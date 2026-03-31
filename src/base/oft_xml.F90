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
!> Free the array of node pointers allocated by \ref oft_xml_get_elements_c.
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_elements_c(elements_ptr) BIND(C, NAME="oft_xml_free_elements")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: elements_ptr !< Pointer to array to free
  END SUBROUTINE oft_xml_free_elements_c
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
  MODULE PROCEDURE xml_extractDataContent_string2
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
  MODULE PROCEDURE xml_extractDataAttribute_string2
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
CONTAINS
!---------------------------------------------------------------------------------
!> Parse an XML file from a given file path.
!---------------------------------------------------------------------------------
SUBROUTINE xml_parsefile(filepath,doc,ierr)
CHARACTER(LEN=*), INTENT(in) :: filepath !< Path to XML file
TYPE(xml_doc), INTENT(out) :: doc !< C pointer to parsed xmlDoc
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_filepath(LEN_TRIM(filepath)+1)
INTEGER(i4) :: i
! WRITE(*,*)1,1
doc%doc=c_null_ptr
DO i=1,LEN_TRIM(filepath)
  c_filepath(i) = filepath(i:i)
END DO
c_filepath(LEN_TRIM(filepath)+1) = c_null_char
! WRITE(*,*)1,2
ierr = INT(oft_xml_parse_file_c(c_filepath, doc%doc), i4)
! WRITE(*,*)1,3
IF(ierr/=0)THEN
    doc%doc=c_null_ptr
    doc%root=>NULL()
  RETURN
END IF
! WRITE(*,*)1,4
ALLOCATE(doc%root)
ierr = INT(oft_xml_get_root_c(doc%doc, doc%root%obj), i4)
! WRITE(*,*)1,5
IF(ierr/=0)THEN
    CALL oft_xml_free_doc_c(doc%doc)
    doc%doc=c_null_ptr
    doc%root=>NULL()
  RETURN
END IF
! WRITE(*,*)1,6
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
CALL oft_xml_free_elements_c(elements_c)
END SUBROUTINE xml_get_element_list
!---------------------------------------------------------------------------------
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_content(node,content,ierr)
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
content=TRIM(ADJUSTL(content_tmp))
DEALLOCATE(content_tmp)
END SUBROUTINE oft_xml_get_content
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
! WRITE(*,*)5,1
has_attr=.FALSE.
DO i=1,LEN_TRIM(attr_name)
  c_name(i)=attr_name(i:i)
END DO
c_name(LEN_TRIM(attr_name)+1)=c_null_char
! WRITE(*,*)5,2,c_name
ierr=INT(oft_xml_has_attribute_c(node%obj,C_LOC(c_name)),i4)
IF(ierr/=0)RETURN
! WRITE(*,*)5,3
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
SUBROUTINE tokenize_string(content,tokens,ierr)
CHARACTER(LEN=*), INTENT(in) :: content !< Output string
INTEGER(i4), ALLOCATABLE, INTENT(out) :: tokens(:,:,:) !< Output array
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,nchars,ncols,nrows,istart,iend,nr_loc,iend2
REAL(r8) :: val_tmp
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
character(len=256) :: msg
character(*), parameter :: NL = new_line('A')
nchars=LEN(content)
ncols=0
DO i=1,nchars
  IF(content(i:i)==NL)THEN
    ncols=ncols+1
  END IF
END DO
ncols=ncols+1
WRITE(*,*)'Parsing: ncols=',ncols
istart=1
DO i=1,ncols
  iend = INDEX(content(istart:nchars),NL)
  IF(iend==0)THEN
    iend=nchars+1
  ELSE
    iend=istart+iend-1
  END IF
  line_strip = TRIM(ADJUSTL(content(istart:iend-1)))
  nr_loc=0
  DO j=1,LEN(line_strip)
    IF(line_strip(j:j)==',')THEN
      nr_loc=nr_loc+1
    END IF
  END DO
  WRITE(*,*)'"',line_strip,'" ',nr_loc
  IF(i==1)THEN
    nrows=nr_loc
  ELSE IF((i==2).AND.(line_strip(1:1)==','))THEN
    nrows=nrows+1
  ELSE IF((i==ncols).AND.(line_strip(1:1)/=','))THEN
    nr_loc=nr_loc+1
  END IF
  IF(nr_loc/=nrows)THEN
    ierr=-10000
    RETURN
  END IF
  istart=iend+1
END DO
IF(ncols==1)nrows=nrows+1
WRITE(*,*)'Parsing: ncols=',ncols,' nrows=',nrows
!
ALLOCATE(tokens(2,nrows,ncols))
istart=1
outer_read: DO i=1,ncols
  iend = INDEX(content(istart+1:nchars),NL)
  IF(iend==0)THEN
    iend=nchars+1
  ELSE
    iend=istart+iend
  END IF
  ! iend=istart+1+MIN(INDEX(content(istart+1:),NL),nchars-istart+1)
  ! ! WRITE(*,*)'"',content(istart:iend-1),'"'
  line_strip = TRIM(ADJUSTL(content(istart:iend-1)))
  ! IF(line_strip(1:1)==NL)istart=istart+INDEX(content(istart:iend-1),NL)
  IF(line_strip(1:1)==',')istart=istart+INDEX(content(istart:iend-1),',')
  DO j=1,nrows
    iend2=INDEX(content(istart:iend-1),',')-1
    IF(iend2<0)iend2=MIN(iend-istart,nchars-istart+1)
    WRITE(*,*)i,j,'"',content(istart:istart+iend2-1),'"',istart,istart+iend2-1,nchars
    tokens(:,j,i)=[istart,istart+iend2-1]
    ! READ(content(istart:istart+iend2-1),*,IOSTAT=ierr,iomsg=msg)val_tmp
    ! IF(ierr/=0)THEN
    !   ! WRITE(*,*)'Error parsing double value at column ',i,' row ',j
    !   ! WRITE(*,*)msg
    !   EXIT outer_read
    ! END IF
    ! output((i-1)*nrows+j)=val_tmp
    istart=istart+iend2+1
  END DO
  istart=istart+1
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
CHARACTER(LEN=*), INTENT(in) :: content !< Output string
REAL(r8), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
REAL(r8) :: val_tmp
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
character(len=256) :: msg

CALL tokenize_string(content,tokens,ierr)
IF(ierr/=0)THEN
  DEALLOCATE(tokens)
  RETURN
END IF

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    READ(content(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr)val_tmp
    IF(ierr/=0)THEN
      ! WRITE(*,*)'Error parsing double value at column ',i,' row ',j
      ! WRITE(*,*)msg
      EXIT outer_read
    END IF
    output((i-1)*nrows+j)=val_tmp
  END DO
END DO outer_read
DEALLOCATE(tokens)
IF(ierr/=0)THEN
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_doubles
!---------------------------------------------------------------------------------
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_integers(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(in) :: content !< Output string
INTEGER(i4), POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
INTEGER(i4) :: val_tmp
CHARACTER(LEN=:), ALLOCATABLE :: line_strip
character(len=256) :: msg

CALL tokenize_string(content,tokens,ierr)
IF(ierr/=0)THEN
  DEALLOCATE(tokens)
  RETURN
END IF

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    READ(content(tokens(1,j,i):tokens(2,j,i)),*,IOSTAT=ierr)val_tmp
    IF(ierr/=0)THEN
      ! WRITE(*,*)'Error parsing integer value at column ',i,' row ',j
      ! WRITE(*,*)msg
      EXIT outer_read
    END IF
    output((i-1)*nrows+j)=val_tmp
  END DO
END DO outer_read
DEALLOCATE(tokens)
IF(ierr/=0)THEN
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_integers
!---------------------------------------------------------------------------------
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE parse_string_to_logicals(content,output,output_shape,ierr)
CHARACTER(LEN=*), INTENT(in) :: content !< Output string
LOGICAL, POINTER, INTENT(out) :: output(:) !< Output array
INTEGER(i4), OPTIONAL, INTENT(out) :: output_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4) :: i,j,ncols,nrows
INTEGER(i4), ALLOCATABLE :: tokens(:,:,:)
CHARACTER(LEN=:), ALLOCATABLE :: token_strip

CALL tokenize_string(content,tokens,ierr)
IF(ierr/=0)THEN
  DEALLOCATE(tokens)
  RETURN
END IF

ncols=SIZE(tokens,3)
nrows=SIZE(tokens,2)
IF(PRESENT(output_shape))output_shape=[nrows,ncols]
ALLOCATE(output(nrows*ncols))
outer_read: DO i=1,ncols
  DO j=1,nrows
    token_strip=TRIM(ADJUSTL(content(tokens(1,j,i):tokens(2,j,i))))
    CALL string_to_upper(token_strip)
    IF((token_strip=='TRUE').OR.(token_strip=='.TRUE.').OR.(token_strip=='1').OR.(token_strip=='T'))THEN
      output((i-1)*nrows+j)=.TRUE.
    ELSE IF((token_strip=='FALSE').OR.(token_strip=='.FALSE.').OR.(token_strip=='0').OR.(token_strip=='F'))THEN
      output((i-1)*nrows+j)=.FALSE.
    ELSE
      ! WRITE(*,*)'Error parsing logical value at column ',i,' row ',j
      ! WRITE(*,*)'Unrecognized token: "',token_strip,'"'
      EXIT outer_read
    END IF
  END DO
END DO outer_read
IF(ierr/=0)THEN
  DEALLOCATE(output)
  RETURN
END IF
END SUBROUTINE parse_string_to_logicals
! !---------------------------------------------------------------------------------
! !> Free an XML document previously parsed with \ref oft_xml_load.
! !---------------------------------------------------------------------------------
! SUBROUTINE xml_extractDataContent_string1(node,data,iostat)
! TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
! CHARACTER(LEN=*), INTENT(out) :: data !< Output string
! INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
! INTEGER(i4) :: ierr
! CHARACTER(LEN=:), ALLOCATABLE :: data_tmp
! CALL oft_xml_get_content(node,data_tmp,ierr)
! IF(PRESENT(iostat))THEN
!   iostat=ierr
! ELSE
!   STOP 'Error extracting string content from XML node in xml_extractDataContent_string1'
!   ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataContent_string1', __FILE__)
! END IF
! IF(LEN(data_tmp)>LEN(data))THEN
!   STOP 'Extracted string from XML content is too long for the provided output variable in xml_extractDataContent_string1'
!   ! CALL oft_abort('Extracted string from XML content is too long for the provided output variable', 'xml_extractDataContent_string1', __FILE__)
! END IF
! data=data_tmp(1:LEN(data))
! DEALLOCATE(data_tmp)
! END SUBROUTINE xml_extractDataContent_string1
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_string2(node,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CALL oft_xml_get_content(node,output,ierr)
IF(PRESENT(iostat))THEN
  iostat=ierr
ELSE
  STOP 'Error extracting string content from XML node in xml_extractDataContent_string2'
  ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataContent_string2', __FILE__)
END IF
END SUBROUTINE xml_extractDataContent_string2
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double1D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), POINTER, INTENT(inout) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_content(node,content,ierr)
WRITE(*,*)'rCHK',ierr
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataContent_double1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataContent_double1D', __FILE__)
  END IF
END IF
CALL parse_string_to_doubles(content,output,data_shape,ierr)
WRITE(*,*)'pCHK',ierr
DEALLOCATE(content)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    RETURN
  ELSE
    STOP 'Error parsing double values from XML node in xml_extractDataContent_double1D'
    ! CALL oft_abort('Error parsing double values from XML node', 'xml_extractDataContent_double1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr
REAL(r8), POINTER :: data_tmp(:)
CALL xml_extractDataContent_double1D(node,data_tmp,data_shape,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
output=data_tmp(1)
DEALLOCATE(data_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_double2D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
REAL(r8), POINTER, INTENT(inout) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: shape_tmp(2),ierr
REAL(r8), POINTER :: data_tmp(:)
CALL xml_extractDataContent_double1D(node,data_tmp,shape_tmp,ierr)
WRITE(*,*)'CHK',ierr
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(data_tmp,shape_tmp)
DEALLOCATE(data_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_double2D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int1D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), POINTER, INTENT(inout) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_content(node,content,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    IF(ALLOCATED(content))DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataContent_int1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataContent_int1D', __FILE__)
  END IF
END IF
CALL parse_string_to_integers(content,output,data_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    IF(ASSOCIATED(output))DEALLOCATE(output)
    RETURN
  ELSE
    STOP 'Error parsing integer values from XML node in xml_extractDataContent_int1D'
    ! CALL oft_abort('Error parsing integer values from XML node', 'xml_extractDataContent_int1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_int1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4), POINTER :: data_tmp(:)
CALL xml_extractDataContent_int1D(node,data_tmp,data_shape,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
output=data_tmp(1)
DEALLOCATE(data_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_int
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_int2D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
INTEGER(i4), POINTER, INTENT(inout) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
INTEGER(i4), POINTER :: data_tmp(:)
CALL xml_extractDataContent_int1D(node,data_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
WRITE(*,*)'CHK',SIZE(data_tmp),shape_tmp,ierr
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
WRITE(*,*)'CHK2',SHAPE(output),shape_tmp
output=RESHAPE(data_tmp,shape_tmp)
IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
WRITE(*,*)'Done',ierr
END SUBROUTINE xml_extractDataContent_int2D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical1D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, POINTER, INTENT(inout) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_content(node,content,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataContent_logical1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataContent_logical1D', __FILE__)
  END IF
END IF
CALL parse_string_to_logicals(content,output,data_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    RETURN
  ELSE
    STOP 'Error parsing logical values from XML node in xml_extractDataContent_logical1D'
    ! CALL oft_abort('Error parsing logical values from XML node', 'xml_extractDataContent_logical1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
LOGICAL, POINTER :: data_tmp(:)
CALL xml_extractDataContent_logical1D(node,data_tmp,data_shape,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
output=data_tmp(1)
DEALLOCATE(data_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataContent_logical2D(node,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
LOGICAL, POINTER, INTENT(inout) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
LOGICAL, POINTER :: data_tmp(:)
CALL xml_extractDataContent_logical1D(node,data_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(data_tmp,shape_tmp)
DEALLOCATE(data_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataContent_logical2D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_string2(node,attr_name,output,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
CHARACTER(LEN=:), ALLOCATABLE, INTENT(inout) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CALL oft_xml_get_attr(node,attr_name,output,ierr)
IF(PRESENT(iostat))THEN
  iostat=ierr
ELSE
  STOP 'Error extracting string content from XML node in xml_extractDataAttribute_string2'
  ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataAttribute_string', __FILE__)
END IF
END SUBROUTINE xml_extractDataAttribute_string2
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double1D(node,attr_name,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(inout) :: output(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataAttribute_double1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataAttribute_double1D', __FILE__)
  END IF
END IF
CALL parse_string_to_doubles(content,output,data_shape,ierr)
DEALLOCATE(content)
WRITE(*,*)'CHK',ierr,ASSOCIATED(output)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    RETURN
  ELSE
    STOP 'Error parsing double values from XML node in xml_extractDataAttribute_double1D'
    ! CALL oft_abort('Error parsing double values from XML node', 'xml_extractDataAttribute_double1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_double1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_double(node,attr_name,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), INTENT(out) :: output !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_double1D(node,attr_name,output_tmp,data_shape,ierr)
WRITE(*,*)'CHK',ierr,ASSOCIATED(output_tmp)
IF(ierr/=0)THEN
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
SUBROUTINE xml_extractDataAttribute_double2D(node,attr_name,output,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
REAL(r8), POINTER, INTENT(inout) :: output(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
REAL(r8), POINTER :: output_tmp(:)
CALL xml_extractDataAttribute_double1D(node,attr_name,output_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(output_tmp))DEALLOCATE(output_tmp)
  RETURN
END IF
ALLOCATE(output(shape_tmp(1),shape_tmp(2)))
output=RESHAPE(output_tmp,shape_tmp)
DEALLOCATE(output_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_double2D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int1D(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(inout) :: data(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataAttribute_int1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataAttribute_int1D', __FILE__)
  END IF
END IF
CALL parse_string_to_integers(content,data,data_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    RETURN
  ELSE
    STOP 'Error parsing integer values from XML node in xml_extractDataAttribute_int1D'
    ! CALL oft_abort('Error parsing integer values from XML node', 'xml_extractDataAttribute_int1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), INTENT(out) :: data !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4), POINTER :: data_tmp(:)
CALL xml_extractDataAttribute_int1D(node,attr_name,data_tmp,data_shape,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
data=data_tmp(1)
DEALLOCATE(data_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_int2D(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
INTEGER(i4), POINTER, INTENT(inout) :: data(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
INTEGER(i4), POINTER :: data_tmp(:)
CALL xml_extractDataAttribute_int1D(node,attr_name,data_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
ALLOCATE(data(shape_tmp(1),shape_tmp(2)))
data=RESHAPE(data_tmp,shape_tmp)
DEALLOCATE(data_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_int2D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical1D(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(inout) :: data(:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
CALL oft_xml_get_attr(node,attr_name,content,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    DEALLOCATE(content)
    RETURN
  ELSE
    STOP 'Error extracting string content from XML node in xml_extractDataAttribute_logical1D'
    ! CALL oft_abort('Error extracting string content from XML node', 'xml_extractDataAttribute_logical1D', __FILE__)
  END IF
END IF
CALL parse_string_to_logicals(content,data,data_shape,ierr)
DEALLOCATE(content)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))THEN
    iostat=ierr
    RETURN
  ELSE
    STOP 'Error parsing logical values from XML node in xml_extractDataAttribute_logical1D'
    ! CALL oft_abort('Error parsing logical values from XML node', 'xml_extractDataAttribute_logical1D', __FILE__)
  END IF
END IF
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical1D
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, INTENT(out) :: data !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
INTEGER(i4) :: ierr
CHARACTER(LEN=:), ALLOCATABLE :: content
LOGICAL, POINTER :: data_tmp(:)
CALL xml_extractDataAttribute_logical1D(node,attr_name,data_tmp,data_shape,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
data=data_tmp(1)
DEALLOCATE(data_tmp)
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE xml_extractDataAttribute_logical2D(node,attr_name,data,data_shape,iostat)
TYPE(xml_node), INTENT(in) :: node !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
LOGICAL, POINTER, INTENT(inout) :: data(:,:) !< Output string
INTEGER(i4), OPTIONAL, INTENT(out) :: data_shape(2) !< Shape of the output array (nrows, ncols)
INTEGER(i4), OPTIONAL, INTENT(out) :: iostat !< I/O status flag (0 on success)
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: ierr,shape_tmp(2)
LOGICAL, POINTER :: data_tmp(:)
CALL xml_extractDataAttribute_logical1D(node,attr_name,data_tmp,shape_tmp,ierr)
IF(ierr/=0)THEN
  IF(PRESENT(iostat))iostat=ierr
  IF(ASSOCIATED(data_tmp))DEALLOCATE(data_tmp)
  RETURN
END IF
ALLOCATE(data(shape_tmp(1),shape_tmp(2)))
data=RESHAPE(data_tmp,shape_tmp)
DEALLOCATE(data_tmp)
IF(PRESENT(data_shape))data_shape=shape_tmp
IF(PRESENT(iostat))iostat=0
END SUBROUTINE xml_extractDataAttribute_logical2D
END MODULE oft_xml
