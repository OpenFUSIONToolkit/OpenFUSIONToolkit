!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (OpenFUSIONToolkit)
!------------------------------------------------------------------------------
!> @file oft_sort.F90
!
!> Portable array sorting and search functions
!!
!! @note On Intel compilers binary sort and search are used.
!! Otherwise C++ templated sort search functions defined in c_sort.cxx are used.
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_base
!-----------------------------------------------------------------------------
MODULE oft_sort
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_long, c_double
USE oft_local
IMPLICIT NONE
#define MOD_NAME oft_sort
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------
!> Sort an array in ascending order
!---------------------------------------------------------------------------
INTERFACE sort_array
  ! Sort 1D array
  MODULE PROCEDURE sort4
  MODULE PROCEDURE sort8
  ! Sort 1D array and companion array (done in C++)
  SUBROUTINE int_sort_array44(list,ind,n) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE int_sort_array44
  SUBROUTINE int_sort_array48(list,ind,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE int_sort_array48
  SUBROUTINE int_sort_array84(list,ind,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE int_sort_array84
  SUBROUTINE int_sort_array88(list,ind,n) BIND(C)
    IMPORT c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE int_sort_array88
  SUBROUTINE real_sort_array4(list,ind,n) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: n
    REAL(c_double), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE real_sort_array4
  SUBROUTINE real_sort_array8(list,ind,n) BIND(C)
    IMPORT c_long, c_double
    INTEGER(c_long), INTENT(in) :: n
    REAL(c_double), DIMENSION(n), INTENT(inout) :: list
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE real_sort_array8
  ! Sort 2D array as 1D
  MODULE PROCEDURE sorta4
  MODULE PROCEDURE sorta8
END INTERFACE sort_array
!---------------------------------------------------------------------------
!> Search a sorted integer array
!---------------------------------------------------------------------------
INTERFACE search_array
  MODULE PROCEDURE search4
  MODULE PROCEDURE search8
END INTERFACE search_array
!---------------------------------------------------------------------------
! External C++ subroutine definitions
!---------------------------------------------------------------------------
INTERFACE
#if !defined(__INTEL_COMPILER)
  SUBROUTINE int_sort44(list,n) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: list
  END SUBROUTINE
  SUBROUTINE int_sort48(list,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: list
  END SUBROUTINE
  SUBROUTINE int_sort84(list,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: list
  END SUBROUTINE
  SUBROUTINE int_sort88(list,n) BIND(C)
    IMPORT c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: list
  END SUBROUTINE
  PURE SUBROUTINE int_search44(list,n,item,result) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(in) :: list
    INTEGER(c_int), INTENT(in) :: item
    INTEGER(c_int), INTENT(in) :: result
  END SUBROUTINE
  PURE SUBROUTINE int_search48(list,n,item,result) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n), INTENT(in) :: list
    INTEGER(c_int), INTENT(in) :: item
    INTEGER(c_long), INTENT(in) :: result
  END SUBROUTINE
  PURE SUBROUTINE int_search84(list,n,item,result) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(in) :: list
    INTEGER(c_long), INTENT(in) :: item
    INTEGER(c_int), INTENT(in) :: result
  END SUBROUTINE
  PURE SUBROUTINE int_search88(list,n,item,result) BIND(C)
    IMPORT c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n), INTENT(in) :: list
    INTEGER(c_long), INTENT(in) :: item
    INTEGER(c_long), INTENT(in) :: result
  END SUBROUTINE
#endif
END INTERFACE
PUBLIC sort_array, search_array
!------------------------------------------------------------------------------
!> Sort rows of a matrix and associated index vector by the first column in ascending order
!!
!! @param[in,out] matrix Matrix to sort [2,n]
!! @param[in,out] ind Index vector [n]
!! @param[in] n Number of rows in matrix
!------------------------------------------------------------------------------
INTERFACE sort_matrix
  SUBROUTINE int_sort_matrix44(matrix,ind,n) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
  SUBROUTINE int_sort_matrix48(matrix,ind,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_int), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
  SUBROUTINE int_sort_matrix84(matrix,ind,n) BIND(C)
    IMPORT c_int, c_long
    INTEGER(c_int), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
  SUBROUTINE int_sort_matrix88(matrix,ind,n) BIND(C)
    IMPORT c_long
    INTEGER(c_long), INTENT(in) :: n
    INTEGER(c_long), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
  SUBROUTINE real_sort_matrix4(matrix,ind,n) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: n
    REAL(c_double), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_int), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
  SUBROUTINE real_sort_matrix8(matrix,ind,n) BIND(C)
    IMPORT c_long, c_double
    INTEGER(c_long), INTENT(in) :: n
    REAL(c_double), DIMENSION(n,2), INTENT(inout) :: matrix
    INTEGER(c_long), DIMENSION(n), INTENT(inout) :: ind
  END SUBROUTINE
END INTERFACE sort_matrix
PUBLIC sort_matrix
CONTAINS
!---------------------------------------------------------------------------
!> Sorts an int4 vector.
!!
!! @note This function should not be called directly,
!! use the oft_sort interface
!---------------------------------------------------------------------------
SUBROUTINE sort4(ia,n)
INTEGER(i4), INTENT(in) :: n
INTEGER(i4), INTENT(inout) :: ia(n)
IF(n==0)RETURN
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
CALL sortqq(loc(ia(1)),int8(n),srt$integer4)
#else
CALL int_sort44(ia,n)
#endif
! DEBUG_STACK_POP
END SUBROUTINE sort4
!---------------------------------------------------------------------------
!> Sorts an int8 vector.
!!
!! @note This function should not be called directly,
!! use the oft_sort interface
!---------------------------------------------------------------------------
SUBROUTINE sort8(ia,n)
INTEGER(i8), INTENT(in) :: n
INTEGER(i8), INTENT(inout) :: ia(n)
IF(n==0)RETURN
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
CALL sortqq(loc(ia(1)),int8(n),srt$integer8)
#else
CALL int_sort88(ia,n)
#endif
! DEBUG_STACK_POP
END SUBROUTINE sort8
!---------------------------------------------------------------------------
!> Sorts an int4 MxN array.
!!
!! @note This function should not be called directly,
!! use the oft_sort interface
!---------------------------------------------------------------------------
SUBROUTINE sorta4(ia,n,m)
INTEGER(i4), INTENT(in) :: n
INTEGER(i4), INTENT(in) :: m
INTEGER(i4), INTENT(inout) :: ia(n,m)
IF(n==0.OR.m==0)RETURN
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
CALL sortqq(loc(ia(1,1)),int8(n*m),srt$integer4)
#else
CALL int_sort44(ia,n*m)
#endif
! DEBUG_STACK_POP
END SUBROUTINE sorta4
!---------------------------------------------------------------------------
!> Sorts an int8 MxN array.
!!
!! @note This function should not be called directly,
!! use the oft_sort interface
!---------------------------------------------------------------------------
SUBROUTINE sorta8(ia,n,m)
INTEGER(i8), INTENT(in) :: n
INTEGER(i8), INTENT(in) :: m
INTEGER(i8), INTENT(inout) :: ia(n,m)
IF(n==0.OR.m==0)RETURN
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
CALL sortqq(loc(ia(1,1)),int8(n*m),srt$integer8)
#else
CALL int_sort88(ia,n*m)
#endif
! DEBUG_STACK_POP
END SUBROUTINE sorta8
!---------------------------------------------------------------------------
!> Search for a value in a sorted integer array
!!
!! @result Index of item in array
!---------------------------------------------------------------------------
FUNCTION search4(item,list,n)
INTEGER(i4), INTENT(in) :: item !< Value to search for
INTEGER(i4), INTENT(in) :: list(n) !< Array to search [n]
INTEGER(i4), INTENT(in) :: n !< Size of array
INTEGER(i4) :: search4
IF(n==0)THEN
  search4=0
  RETURN
END IF
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
search4=bsearchqq(loc(item),loc(list(1)),int8(n),srt$integer4)
#else
CALL int_search44(list,n,item,search4)
#endif
! DEBUG_STACK_POP
END FUNCTION search4
!---------------------------------------------------------------------------
!> Search for a value in a sorted integer array
!!
!! @result Index of item in array
!---------------------------------------------------------------------------
FUNCTION search8(item,list,n)
INTEGER(i8), intent(in) :: item !< Value to search for
INTEGER(i8), intent(in) :: list(n) !< Array to search [n]
INTEGER(i8), intent(in) :: n !< Size of array
INTEGER(i8) :: search8
IF(n==0)THEN
  search8=0
  RETURN
END IF
! DEBUG_STACK_PUSH
#ifdef __INTEL_COMPILER
search8=bsearchqq(loc(item),loc(list(1)),int8(n),srt$integer8)
#else
CALL int_search88(list,n,item,search8)
#endif
! DEBUG_STACK_POP
END FUNCTION search8
END MODULE oft_sort
