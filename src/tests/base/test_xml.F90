!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_xml.F90
!
!> Regression tests for some XML I/O
!!
!! @authors Chris Hansen
!! @date March 2025
!! @ingroup testing
!--------------------------------------------------------------------------------
PROGRAM test_xml
USE oft_base
IMPLICIT NONE
LOGICAL :: success,bool_tmp
INTEGER(i4) :: iounit,ierr,data_shape(2),int_tmp
INTEGER(i4), POINTER :: data_int(:),data2D_int(:,:)
REAL(r8) :: real_tmp
REAL(r8), POINTER :: data(:), data2D(:,:)
TYPE(xml_doc) :: doc
TYPE(xml_node) :: root_node,current_node
TYPE(xml_nodelist) :: child_nodes
CHARACTER(LEN=:), ALLOCATABLE :: content
INTEGER(i4) :: test_id=1
NAMELIST/test_io_options/test_id
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init

!---Parse test file
CALL xml_parsefile('test_file.xml',doc,ierr)
IF(ierr/=0)CALL oft_abort('Error parsing XML input file','test_xml',__FILE__)

!---Look for element in root node
WRITE(*,*)'Test 1'
CALL xml_get_element(doc%root,'test1',current_node,ierr)
IF(ierr/=0)CALL oft_abort('Error finding XML element "test1"','test_xml',__FILE__)
!---Read content
CALL xml_read_content(current_node,content,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting string content from XML node','test_xml',__FILE__) 
!---Read attribute
success=xml_hasAttribute(current_node,'attr1')
IF(.NOT.success)CALL oft_abort('Error finding "attr1" in XML node "test1"','test_xml',__FILE__)
CALL xml_read_attribute(current_node,'attr1',content,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting string attribute from XML node','test_xml',__FILE__)
!---Confirm missing attribute is handled correctly
success=xml_hasAttribute(current_node,'missing_attr')
IF(success)CALL oft_abort('Unexpectedly found "missing_attr" in XML node "test1"','test_xml',__FILE__)

!---Test 1 for parsing of numeric content and attributes
WRITE(*,*)'Test 2'
CALL xml_get_element(doc%root,'test2',current_node,ierr)
IF(ierr/=0)CALL oft_abort('Error finding XML element','test_xml',__FILE__)
CALL xml_read_content(current_node,data,data_shape,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting double values from XML node','test_xml',__FILE__)
IF(data_shape(1) /= 4 .OR. data_shape(2) /= 1)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
CALL xml_read_content(current_node,data2D,data_shape,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting double values from XML node','test_xml',__FILE__)
IF(data_shape(1) /= 4 .OR. data_shape(2) /= 1)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
CALL xml_read_attribute(current_node,'attr1',real_tmp,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting double attribute from XML node','test_xml',__FILE__)
IF(ABS(real_tmp-2.d0)>1.d-8)CALL oft_abort('Unexpected value for attr1','test_xml',__FILE__)
CALL xml_read_attribute(current_node,'attr2',real_tmp,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting double attribute from XML node','test_xml',__FILE__)
IF(ABS(real_tmp-1.d-3)>1.d-8)CALL oft_abort('Unexpected value for attr2','test_xml',__FILE__)
CALL xml_read_attribute(current_node,'attr3',real_tmp,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting double attribute from XML node','test_xml',__FILE__)
IF(ABS(real_tmp+1.d-3)>1.d-8)CALL oft_abort('Unexpected value for attr3','test_xml',__FILE__)

!---Test 2 for parsing of numeric content and attributes
WRITE(*,*)'Test 3'
CALL xml_get_element(doc%root,'test3',current_node,ierr)
IF(ierr/=0)CALL oft_abort('Error finding XML element','test_xml',__FILE__)
CALL xml_read_content(current_node,data_int,data_shape,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting integer values from XML node','test_xml',__FILE__)
IF(data_shape(1) /= 3 .OR. data_shape(2) /= 3)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
CALL xml_read_content(current_node,data2D_int,data_shape,iostat=ierr)
IF(ierr/=0)CALL oft_abort('Error extracting integer values from XML node','test_xml',__FILE__)
IF(data_shape(1) /= 3 .OR. data_shape(2) /= 3)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr1',int_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting integer attribute from XML node','test_xml',__FILE__)
! IF(int_tmp /= 2)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr2',int_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting integer attribute from XML node','test_xml',__FILE__)
! IF(int_tmp /= -2)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)

! !---Test 3 for parsing of numeric content and attributes
! WRITE(*,*)'Test 4'
! CALL xml_get_element(doc%root,'test4',current_node,ierr)
! IF(ierr/=0)CALL oft_abort('Error finding XML element','test_xml',__FILE__)
! CALL xml_read_content(current_node,data,data_shape,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting double values from XML node','test_xml',__FILE__)
! IF(data_shape(1) /= 3 .OR. data_shape(2) /= 3)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
! CALL xml_read_content(current_node,data2D,data_shape,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting double values from XML node','test_xml',__FILE__)
! IF(data_shape(1) /= 3 .OR. data_shape(2) /= 3)CALL oft_abort('Unexpected data shape extracted from XML node','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr1',bool_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting boolean attribute from XML node','test_xml',__FILE__)
! IF(bool_tmp)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr2',bool_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting boolean attribute from XML node','test_xml',__FILE__)
! IF(.NOT.bool_tmp)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr3',bool_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting boolean attribute from XML node','test_xml',__FILE__)
! IF(.NOT.bool_tmp)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)
! CALL xml_read_attribute(current_node,'attr4',bool_tmp,iostat=ierr)
! IF(ierr/=0)CALL oft_abort('Error extracting boolean attribute from XML node','test_xml',__FILE__)
! IF(bool_tmp)CALL oft_abort('Unexpected value extracted from XML attribute','test_xml',__FILE__)

!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_xml