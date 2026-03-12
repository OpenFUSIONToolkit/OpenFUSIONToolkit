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
  c_null_ptr, c_f_pointer, c_associated, c_double
USE oft_local, ONLY: i4, r8, string_to_upper
IMPLICIT NONE
!> Maximum length (in characters) for XML content and attribute value buffers
INTEGER(i4), PARAMETER :: OFT_XML_SLEN = 2048
!------------------------------------------------------------------------------
! C-binding interfaces for string parsing functions in oft_xml_c.c
!------------------------------------------------------------------------------
INTERFACE
!------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of 32-bit
!! integers, returning a 2-component shape (nrows, ncols).
!------------------------------------------------------------------------------
  FUNCTION oft_xml_parse_int_array_c(str, arr, max_n, shape) &
      BIND(C, NAME="oft_xml_parse_int_array_c") RESULT(ierr)
    IMPORT c_int, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: str(*) !< Null-terminated input string
    INTEGER(c_int), INTENT(out) :: arr(*) !< Output flat integer array
    INTEGER(c_int), VALUE, INTENT(in) :: max_n !< Maximum number of elements
    INTEGER(c_int), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_parse_int_array_c
!------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of doubles,
!! returning a 2-component shape (nrows, ncols).
!------------------------------------------------------------------------------
  FUNCTION oft_xml_parse_real_array_c(str, arr, max_n, shape) &
      BIND(C, NAME="oft_xml_parse_real_array_c") RESULT(ierr)
    IMPORT c_int, c_char, c_double
    CHARACTER(KIND=c_char), INTENT(in) :: str(*) !< Null-terminated input string
    REAL(c_double), INTENT(out) :: arr(*) !< Output flat real array
    INTEGER(c_int), VALUE, INTENT(in) :: max_n !< Maximum number of elements
    INTEGER(c_int), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_parse_real_array_c
!------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of logical
!! values encoded as integers (1=true, 0=false), returning a 2-component shape.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_parse_logical_array_c(str, arr, max_n, shape) &
      BIND(C, NAME="oft_xml_parse_logical_array_c") RESULT(ierr)
    IMPORT c_int, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: str(*) !< Null-terminated input string
    INTEGER(c_int), INTENT(out) :: arr(*) !< Output flat logical-as-int array
    INTEGER(c_int), VALUE, INTENT(in) :: max_n !< Maximum number of elements
    INTEGER(c_int), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_parse_logical_array_c
END INTERFACE
#ifdef HAVE_LIBXML2
!------------------------------------------------------------------------------
! C-binding interfaces for oft_xml_c.c
!------------------------------------------------------------------------------
INTERFACE
!------------------------------------------------------------------------------
!> Parse an XML file from a given file path, returning a pointer to the xmlDoc.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_load_file_c(filepath, doc_ptr) BIND(C, NAME="oft_xml_load_file") &
      RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: filepath(*) !< Null-terminated file path
    TYPE(c_ptr), INTENT(out) :: doc_ptr !< Pointer to parsed xmlDoc on return
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_load_file_c
!------------------------------------------------------------------------------
!> Retrieve a pointer to the root element of an XML document.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_root_c(doc_ptr, root_ptr) BIND(C, NAME="oft_xml_get_root") &
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
  FUNCTION oft_xml_get_element_c(parent_ptr, name, index, element_ptr) &
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
  FUNCTION oft_xml_get_elements_c(parent_ptr, name, n, elements_ptr) &
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
  FUNCTION oft_xml_get_content_c(node_ptr, content, max_len) &
      BIND(C, NAME="oft_xml_get_content") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    CHARACTER(KIND=c_char), INTENT(out) :: content(*) !< Output character buffer
    INTEGER(c_int), VALUE, INTENT(in) :: max_len !< Buffer size including null terminator
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_content_c
!------------------------------------------------------------------------------
!> Extract the string value of a given attribute on a given xml node.
!------------------------------------------------------------------------------
  FUNCTION oft_xml_get_attribute_c(node_ptr, attr_name, value, max_len) &
      BIND(C, NAME="oft_xml_get_attribute") RESULT(ierr)
    IMPORT c_int, c_ptr, c_char
    TYPE(c_ptr), VALUE, INTENT(in) :: node_ptr !< Pointer to xmlNode
    CHARACTER(KIND=c_char), INTENT(in) :: attr_name(*) !< Null-terminated attribute name
    CHARACTER(KIND=c_char), INTENT(out) :: value(*) !< Output character buffer
    INTEGER(c_int), VALUE, INTENT(in) :: max_len !< Buffer size including null terminator
    INTEGER(c_int) :: ierr !< 0 on success, nonzero on error
  END FUNCTION oft_xml_get_attribute_c
!------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load_file_c.
!------------------------------------------------------------------------------
  SUBROUTINE oft_xml_free_doc_c(doc_ptr) BIND(C, NAME="oft_xml_free_doc")
    IMPORT c_ptr
    TYPE(c_ptr), VALUE, INTENT(in) :: doc_ptr !< Pointer to xmlDoc to free
  END SUBROUTINE oft_xml_free_doc_c
END INTERFACE
#endif
CONTAINS
#ifdef HAVE_LIBXML2
!---------------------------------------------------------------------------------
!> Parse an XML file from a given file path.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_load(filepath, doc_ptr, ierr)
CHARACTER(LEN=*), INTENT(in) :: filepath !< Path to XML file
TYPE(c_ptr), INTENT(out) :: doc_ptr !< C pointer to parsed xmlDoc
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_filepath(LEN_TRIM(filepath)+1)
INTEGER(i4) :: i
doc_ptr = c_null_ptr
DO i = 1, LEN_TRIM(filepath)
  c_filepath(i) = filepath(i:i)
END DO
c_filepath(LEN_TRIM(filepath)+1) = c_null_char
ierr = INT(oft_xml_load_file_c(c_filepath, doc_ptr), i4)
END SUBROUTINE oft_xml_load
!---------------------------------------------------------------------------------
!> Retrieve a pointer to the root element of an XML document.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_root(doc_ptr, root_ptr, ierr)
TYPE(c_ptr), INTENT(in) :: doc_ptr !< C pointer to xmlDoc
TYPE(c_ptr), INTENT(out) :: root_ptr !< C pointer to root xmlNode
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
root_ptr = c_null_ptr
ierr = INT(oft_xml_get_root_c(doc_ptr, root_ptr), i4)
END SUBROUTINE oft_xml_get_root
!---------------------------------------------------------------------------------
!> Retrieve a pointer to the I-th xml node with a given name contained within
!! a parent node (1-based index, defaults to 1).
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_element(parent_ptr, name, element_ptr, ierr, index)
TYPE(c_ptr), INTENT(in) :: parent_ptr !< C pointer to parent xmlNode
CHARACTER(LEN=*), INTENT(in) :: name !< Element name to find
TYPE(c_ptr), INTENT(out) :: element_ptr !< C pointer to matching xmlNode
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
INTEGER(i4), OPTIONAL, INTENT(in) :: index !< 1-based index (default 1)
CHARACTER(KIND=c_char) :: c_name(LEN_TRIM(name)+1)
INTEGER(c_int) :: req_index
INTEGER(i4) :: i
element_ptr = c_null_ptr
DO i = 1, LEN_TRIM(name)
  c_name(i) = name(i:i)
END DO
c_name(LEN_TRIM(name)+1) = c_null_char
req_index = 1_c_int
IF (PRESENT(index)) req_index = INT(index, c_int)
ierr = INT(oft_xml_get_element_c(parent_ptr, c_name, req_index, element_ptr), i4)
END SUBROUTINE oft_xml_get_element
!---------------------------------------------------------------------------------
!> Retrieve C pointers to all xml nodes with a given name contained within a
!! parent node.
!!
!! @p elements_c is set to a raw C pointer to a heap-allocated array of
!! @c c_ptr values (one per matching child).  Use @c c_f_pointer to obtain a
!! Fortran array over this memory, and call \ref oft_xml_free_elements with
!! @p elements_c when done.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_all_elements(parent_ptr, name, n, elements_c, ierr)
TYPE(c_ptr), INTENT(in) :: parent_ptr !< C pointer to parent xmlNode
CHARACTER(LEN=*), INTENT(in) :: name !< Element name to find
INTEGER(i4), INTENT(out) :: n !< Number of matching children
TYPE(c_ptr), INTENT(out) :: elements_c !< Raw C pointer to array of c_ptr
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_name(LEN_TRIM(name)+1)
INTEGER(c_int) :: n_c
INTEGER(i4) :: i
elements_c = c_null_ptr
n = 0
DO i = 1, LEN_TRIM(name)
  c_name(i) = name(i:i)
END DO
c_name(LEN_TRIM(name)+1) = c_null_char
ierr = INT(oft_xml_get_elements_c(parent_ptr, c_name, n_c, elements_c), i4)
IF (ierr == 0) n = INT(n_c, i4)
END SUBROUTINE oft_xml_get_all_elements
!---------------------------------------------------------------------------------
!> Free the array of node pointers allocated by \ref oft_xml_get_all_elements.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_free_elements(elements_c)
TYPE(c_ptr), INTENT(inout) :: elements_c !< Raw C pointer to free
IF (c_associated(elements_c)) THEN
  CALL oft_xml_free_elements_c(elements_c)
  elements_c = c_null_ptr
END IF
END SUBROUTINE oft_xml_free_elements
!---------------------------------------------------------------------------------
!> Extract the string content from a given xml node into a Fortran string.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_content(node_ptr, content, ierr)
TYPE(c_ptr), INTENT(in) :: node_ptr !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(out) :: content !< Output string
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: buf(OFT_XML_SLEN)
INTEGER(i4) :: i
content = ''
ierr = INT(oft_xml_get_content_c(node_ptr, buf, INT(OFT_XML_SLEN, c_int)), i4)
IF (ierr /= 0) RETURN
DO i = 1, MIN(LEN(content), OFT_XML_SLEN-1)
  IF (buf(i) == c_null_char) EXIT
  content(i:i) = buf(i)
END DO
END SUBROUTINE oft_xml_get_content
!---------------------------------------------------------------------------------
!> Extract the string value of a given attribute on a given xml node.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_get_attr(node_ptr, attr_name, value, ierr)
TYPE(c_ptr), INTENT(in) :: node_ptr !< C pointer to xmlNode
CHARACTER(LEN=*), INTENT(in) :: attr_name !< Attribute name
CHARACTER(LEN=*), INTENT(out) :: value !< Output string
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_name(LEN_TRIM(attr_name)+1)
CHARACTER(KIND=c_char) :: buf(OFT_XML_SLEN)
INTEGER(i4) :: i
value = ''
DO i = 1, LEN_TRIM(attr_name)
  c_name(i) = attr_name(i:i)
END DO
c_name(LEN_TRIM(attr_name)+1) = c_null_char
ierr = INT(oft_xml_get_attribute_c(node_ptr, c_name, buf, INT(OFT_XML_SLEN, c_int)), i4)
IF (ierr /= 0) RETURN
DO i = 1, MIN(LEN(value), OFT_XML_SLEN-1)
  IF (buf(i) == c_null_char) EXIT
  value(i:i) = buf(i)
END DO
END SUBROUTINE oft_xml_get_attr
!---------------------------------------------------------------------------------
!> Free an XML document previously parsed with \ref oft_xml_load.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_free(doc_ptr)
TYPE(c_ptr), INTENT(inout) :: doc_ptr !< C pointer to xmlDoc to free
IF (c_associated(doc_ptr)) THEN
  CALL oft_xml_free_doc_c(doc_ptr)
  doc_ptr = c_null_ptr
END IF
END SUBROUTINE oft_xml_free
#endif
!---------------------------------------------------------------------------------
!> Parse a string into a single logical value.
!!
!! Accepted representations (case-insensitive): @c T, @c F, @c TRUE,
!! @c FALSE, @c .TRUE., @c .FALSE., @c 1 (true), @c 0 (false).
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_logical(str, val, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input string
LOGICAL, INTENT(out) :: val !< Parsed value
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(LEN=LEN_TRIM(str)) :: tmp
INTEGER(i4) :: read_err
val = .FALSE.
ierr = 0
tmp = ADJUSTL(TRIM(str))
CALL string_to_upper(tmp)
SELECT CASE (TRIM(tmp))
CASE ('T', 'TRUE', '.TRUE.')
  val = .TRUE.
CASE ('F', 'FALSE', '.FALSE.')
  val = .FALSE.
CASE ('1')
  val = .TRUE.
CASE ('0')
  val = .FALSE.
CASE DEFAULT
  READ(tmp, *, IOSTAT=read_err) val
  IF (read_err /= 0) ierr = read_err
END SELECT
END SUBROUTINE oft_xml_parse_logical
!---------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of logical values.
!!
!! Commas separate columns; newlines separate rows.  The data is read into a
!! flat array and the 2-component @p shape (nrows, ncols) is returned.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_logical_array(str, vals, shape, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input comma-and-newline-delimited string
LOGICAL, INTENT(out) :: vals(:) !< Output flat logical array
INTEGER(i4), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_str(LEN_TRIM(str)+1)
INTEGER(c_int), ALLOCATABLE :: int_buf(:)
INTEGER(c_int) :: shape_c(2)
INTEGER(i4) :: i, ntotal
shape = 0
ierr = 0
DO i = 1, LEN_TRIM(str)
  c_str(i) = str(i:i)
END DO
c_str(LEN_TRIM(str)+1) = c_null_char
ALLOCATE(int_buf(SIZE(vals)))
ierr = INT(oft_xml_parse_logical_array_c(c_str, int_buf, INT(SIZE(vals), c_int), &
    shape_c), i4)
IF (ierr == 0) THEN
  shape(1) = INT(shape_c(1), i4)
  shape(2) = INT(shape_c(2), i4)
  ntotal = shape(1) * shape(2)
  DO i = 1, ntotal
    vals(i) = (int_buf(i) /= 0_c_int)
  END DO
END IF
DEALLOCATE(int_buf)
END SUBROUTINE oft_xml_parse_logical_array
!---------------------------------------------------------------------------------
!> Parse a string into a single integer value.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_int(str, val, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input string
INTEGER(i4), INTENT(out) :: val !< Parsed value
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(LEN=LEN_TRIM(str)) :: tmp
ierr = 0
tmp = ADJUSTL(TRIM(str))
READ(tmp, *, IOSTAT=ierr) val
END SUBROUTINE oft_xml_parse_int
!---------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of 32-bit integers.
!!
!! Commas separate columns; newlines separate rows.  The data is read into a
!! flat array and the 2-component @p shape (nrows, ncols) is returned.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_int_array(str, vals, shape, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input comma-and-newline-delimited string
INTEGER(i4), INTENT(out) :: vals(:) !< Output flat integer array
INTEGER(i4), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_str(LEN_TRIM(str)+1)
INTEGER(c_int) :: shape_c(2)
INTEGER(i4) :: i
shape = 0
ierr = 0
DO i = 1, LEN_TRIM(str)
  c_str(i) = str(i:i)
END DO
c_str(LEN_TRIM(str)+1) = c_null_char
ierr = INT(oft_xml_parse_int_array_c(c_str, vals, INT(SIZE(vals), c_int), shape_c), i4)
IF (ierr == 0) THEN
  shape(1) = INT(shape_c(1), i4)
  shape(2) = INT(shape_c(2), i4)
END IF
END SUBROUTINE oft_xml_parse_int_array
!---------------------------------------------------------------------------------
!> Parse a string into a single double-precision real value.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_real(str, val, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input string
REAL(r8), INTENT(out) :: val !< Parsed value
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(LEN=LEN_TRIM(str)) :: tmp
ierr = 0
tmp = ADJUSTL(TRIM(str))
READ(tmp, *, IOSTAT=ierr) val
END SUBROUTINE oft_xml_parse_real
!---------------------------------------------------------------------------------
!> Parse a comma-and-newline-delimited string into a flat array of
!! double-precision real values.
!!
!! Commas separate columns; newlines separate rows.  The data is read into a
!! flat array and the 2-component @p shape (nrows, ncols) is returned.
!---------------------------------------------------------------------------------
SUBROUTINE oft_xml_parse_real_array(str, vals, shape, ierr)
CHARACTER(LEN=*), INTENT(in) :: str !< Input comma-and-newline-delimited string
REAL(r8), INTENT(out) :: vals(:) !< Output flat real array
INTEGER(i4), INTENT(out) :: shape(2) !< shape(1)=nrows, shape(2)=ncols
INTEGER(i4), INTENT(out) :: ierr !< Error flag (0 on success)
CHARACTER(KIND=c_char) :: c_str(LEN_TRIM(str)+1)
INTEGER(c_int) :: shape_c(2)
INTEGER(i4) :: i
shape = 0
ierr = 0
DO i = 1, LEN_TRIM(str)
  c_str(i) = str(i:i)
END DO
c_str(LEN_TRIM(str)+1) = c_null_char
ierr = INT(oft_xml_parse_real_array_c(c_str, vals, INT(SIZE(vals), c_int), shape_c), i4)
IF (ierr == 0) THEN
  shape(1) = INT(shape_c(1), i4)
  shape(2) = INT(shape_c(2), i4)
END IF
END SUBROUTINE oft_xml_parse_real_array
END MODULE oft_xml
