!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_deriv_matrices.F90
!
!> Abstract operator interfaces and native matrix implementations
!!
!! Abstract interface definitions
!! - Abstract matrix class
!!
!! Native vector implementations
!! - CSR matrix class
!! - Identity matrix class
!! - Diagonal matrix class
!! - JMLB matrix class
!!
!! @authors Chris Hansen
!! @date Feburary 2012
!! @ingroup la
!------------------------------------------------------------------------------
MODULE oft_deriv_matrices
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix
IMPLICIT NONE
#include "local.h"
private
!------------------------------------------------------------------------------
! CLASS oft_noop_matrix
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_noop_matrix
contains
  procedure :: apply_real => noop_apply_real
  procedure :: apply_complex => noop_apply_comp
  procedure :: applyt_real => noop_apply_real
  procedure :: applyt_complex => noop_apply_comp
  procedure :: set_values => noop_set_values
  procedure :: add_values => noop_add_values
  procedure :: atomic_add_values => noop_add_values
  procedure :: assemble => noop_assemble
  procedure :: zero => noop_zero
  procedure :: zero_rows => noop_zero_rows
  procedure :: delete => noop_delete
end type oft_noop_matrix
!------------------------------------------------------------------------------
! CLASS oft_diagmatrix
!------------------------------------------------------------------------------
!> Diagonal matrix class
!------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_diagmatrix
contains
  procedure :: apply_real => diagmat_apply_real
  procedure :: apply_complex => diagmat_apply_comp
  procedure :: applyt_real => diagmat_apply_real
  procedure :: applyt_complex => diagmat_apply_comp
  procedure :: assemble => diagmat_assemble
  procedure :: zero => diagmat_zero
  procedure :: delete => diagmat_delete
end type oft_diagmatrix
!------------------------------------------------------------------------------
! CLASS oft_sum_matrix
!------------------------------------------------------------------------------
!> Sum matrix class
!------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_sum_matrix
  real(r8) :: alam !< Lambda factor
  class(oft_matrix), pointer :: J => NULL() !< J matrix
  class(oft_matrix), pointer :: K => NULL() !< K matrix
contains
  procedure :: apply_real => sum_mat_apply_real
  procedure :: apply_complex => sum_mat_apply_comp
  procedure :: applyt_real => sum_mat_applyt_real
  procedure :: applyt_complex => sum_mat_applyt_comp
  !> Complete matrix assembly
  procedure :: assemble => sum_mat_assemble
end type oft_sum_matrix
!------------------------------------------------------------------------------
! CLASS oft_mf_matrix
!------------------------------------------------------------------------------
!> Wrapper matrix for matrix-free Jacobians
!!
!! Computes a Finite Difference approximation to the Jacobian of a function
!! \f$ F'(v) = [F(u0 + h * v) - F(u0)]/h \f$
!------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_mf_matrix
  real(r8) :: b0 = 1.d-6 !< Approximate relative error in function eval
  class(oft_matrix), pointer :: f => NULL() !< Non-linear function
  class(oft_vector), pointer :: u0 => NULL() !< Linearization point
  class(oft_vector), pointer :: f0 => NULL() !< Linearization value "F(u0)"
  class(oft_vector), pointer :: tmp => NULL() !< Internal work vector
  class(oft_vector), pointer :: utyp => NULL() !< Typical values vector
contains
  procedure :: apply_real => mf_mat_apply_real
  !> Setup MF jacobian operator
  procedure :: setup => mf_mat_setup
  !> Update linearization point
  procedure :: update => mf_mat_update
  !> Cleanup internal storage
  procedure :: delete => mf_mat_delete
end type oft_mf_matrix
public create_diagmatrix
contains
!------------------------------------------------------------------------------
! SUBROUTINE: noop_apply
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine noop_apply_real(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
call oft_abort('No matrix type specified','noop_apply_real',__FILE__)
end subroutine noop_apply_real
!------------------------------------------------------------------------------
! SUBROUTINE: noop_apply
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine noop_apply_comp(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
call oft_abort('No matrix type specified','noop_apply_comp',__FILE__)
end subroutine noop_apply_comp
!------------------------------------------------------------------------------
! SUBROUTINE: noop_applyt
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = selfT * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine noop_applyt(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
call oft_abort('No matrix type specified','noop_applyt',__FILE__)
end subroutine noop_applyt
!------------------------------------------------------------------------------
! SUBROUTINE: noop_set_values
!------------------------------------------------------------------------------
!> Set values of a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] i_inds Row indices of entries to set [n]
!! @param[in] j_inds Column indices of entries to set [m]
!! @param[in] b Values to set [n,m]
!! @param[in] n Number of rows in local matrix
!! @param[in] m Number of columns in local matrix
!! @param[in] iblock Row block (optional)
!! @param[in] jblock Column block (optional)
!------------------------------------------------------------------------------
subroutine noop_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n)
integer(i4), intent(in) :: j_inds(m)
real(r8), intent(in) :: b(n,m)
integer(i4), intent(in) :: n
integer(i4), intent(in) :: m
integer(i4), optional, intent(in) :: iblock
integer(i4), optional, intent(in) :: jblock
call oft_abort('Invalid operation for matrix type','noop_set_values',__FILE__)
end subroutine noop_set_values
!------------------------------------------------------------------------------
! SUBROUTINE: noop_add_values
!------------------------------------------------------------------------------
!> Add values to a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] i_inds Row indices of entries to set [n]
!! @param[in] j_inds Column indices of entries to set [m]
!! @param[in] b Values to add [n,m]
!! @param[in] n Number of rows in local matrix
!! @param[in] m Number of columns in local matrix
!! @param[in] iblock Row block (optional)
!! @param[in] jblock Column block (optional)
!------------------------------------------------------------------------------
subroutine noop_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n)
integer(i4), intent(in) :: j_inds(m)
real(r8), intent(in) :: b(n,m)
integer(i4), intent(in) :: n
integer(i4), intent(in) :: m
integer(i4), optional, intent(in) :: iblock
integer(i4), optional, intent(in) :: jblock
integer(i4), optional, intent(inout) :: loc_cache(n,m)
call oft_abort('Invalid operation for matrix type','noop_add_values',__FILE__)
end subroutine noop_add_values
!------------------------------------------------------------------------------
! SUBROUTINE: noop_assemble
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in,out] diag Diagonal entries of matrix [nr] (optional)
!------------------------------------------------------------------------------
subroutine noop_assemble(self,diag)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag
call oft_abort('No matrix type specified','noop_assemble',__FILE__)
end subroutine noop_assemble
!------------------------------------------------------------------------------
! SUBROUTINE: noop_zero
!------------------------------------------------------------------------------
!> Zero all entries in matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine noop_zero(self)
class(oft_noop_matrix), intent(inout) :: self
call oft_abort('No matrix type specified','noop_zero',__FILE__)
end subroutine noop_zero
!------------------------------------------------------------------------------
! SUBROUTINE: noop_zero_rows
!------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] nrows Number of rows to zero
!! @param[in] irows Indices of rows to zero [nrows]
!! @param[in] iblock Row block (optional)
!! @param[in] keep_diag Keep diagonal entries
!------------------------------------------------------------------------------
subroutine noop_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: nrows
integer(i4), intent(in) :: irows(nrows)
integer(i4), optional, intent(in) :: iblock
logical, optional, intent(in) :: keep_diag
call oft_abort('No matrix type specified','noop_zero_rows',__FILE__)
end subroutine noop_zero_rows
!------------------------------------------------------------------------------
! SUBROUTINE: noop_delete
!------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine noop_delete(self)
class(oft_noop_matrix), intent(inout) :: self
call oft_warn('Finalizing general matrix, this may indicate an error.')
end subroutine noop_delete
!------------------------------------------------------------------------------
! SUBROUTINE: create_diagmatrix
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine create_diagmatrix(self,a)
class(oft_matrix), pointer, intent(inout) :: self
class(oft_vector), intent(inout) :: a
DEBUG_STACK_PUSH
ALLOCATE(oft_diagmatrix::self)
CALL a%new(self%D)
CALL self%D%add(0.d0,1.d0,a)
self%nr=a%n; self%nrg=a%ng
self%nc=a%n; self%ncg=a%ng
DEBUG_STACK_POP
end subroutine create_diagmatrix
!------------------------------------------------------------------------------
! FUNCTION: diagmatrix_cast
!------------------------------------------------------------------------------
!> Cast a matrix object to a oft_diagmatrix.
!!
!! The source matrix must be oft_diagmatrix or a child class, otherwise an error will be thrown.
!!
!! @param[out] self Pointer to cast crsmatrix
!! @param[in] source Source matrix to cast
!------------------------------------------------------------------------------
FUNCTION diagmatrix_cast(self,source) result(ierr)
class(oft_diagmatrix), pointer, intent(out) :: self
class(oft_matrix), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_diagmatrix)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION diagmatrix_cast
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_apply_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine diagmat_apply_real(self,a,b)
class(oft_diagmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','diagmat_apply_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','diagmat_apply_real',__FILE__)
call b%add(0.d0,1.d0,a)
call b%mult(self%D)
DEBUG_STACK_POP
end subroutine diagmat_apply_real
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_apply_comp
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine diagmat_apply_comp(self,a,b)
class(oft_diagmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','diagmat_apply_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','diagmat_apply_comp',__FILE__)
call b%add((0.d0,0.d0),(1.d0,0.d0),a)
call b%mult(self%D)
DEBUG_STACK_POP
end subroutine diagmat_apply_comp
!------------------------------------------------------------------------------
! SUBROUTINE: noop_set_values
!------------------------------------------------------------------------------
!> Set values of a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] i_inds Row indices of entries to set [n]
!! @param[in] j_inds Column indices of entries to set [m]
!! @param[in] b Values to set [n,m]
!! @param[in] n Number of rows in local matrix
!! @param[in] m Number of columns in local matrix
!! @param[in] iblock Row block (optional)
!! @param[in] jblock Column block (optional)
!------------------------------------------------------------------------------
subroutine diagmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_diagmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n)
integer(i4), intent(in) :: j_inds(m)
real(r8), intent(in) :: b(n,m)
integer(i4), intent(in) :: n
integer(i4), intent(in) :: m
integer(i4), optional, intent(in) :: iblock
integer(i4), optional, intent(in) :: jblock
call oft_abort('Invalid operation for matrix type','diagmat_set_values',__FILE__)
end subroutine diagmat_set_values
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_add_values
!------------------------------------------------------------------------------
!> Add values to a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] i_inds Row indices of entries to set [n]
!! @param[in] j_inds Column indices of entries to set [m]
!! @param[in] b Values to add [n,m]
!! @param[in] n Number of rows in local matrix
!! @param[in] m Number of columns in local matrix
!! @param[in] iblock Row block (optional)
!! @param[in] jblock Column block (optional)
!------------------------------------------------------------------------------
subroutine diagmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_diagmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n)
integer(i4), intent(in) :: j_inds(m)
real(r8), intent(in) :: b(n,m)
integer(i4), intent(in) :: n
integer(i4), intent(in) :: m
integer(i4), optional, intent(in) :: iblock
integer(i4), optional, intent(in) :: jblock
integer(i4), optional, intent(inout) :: loc_cache(n,m)
call oft_abort('Invalid operation for matrix type','diagmat_add_values',__FILE__)
end subroutine diagmat_add_values
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_assemble
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!!
!! @param[in,out] diag Diagonal entries of matrix [nr] (optional)
!------------------------------------------------------------------------------
SUBROUTINE diagmat_assemble(self,diag)
CLASS(oft_diagmatrix), INTENT(inout) :: self
CLASS(oft_vector), OPTIONAL, TARGET, INTENT(inout) :: diag
DEBUG_STACK_PUSH
!---Setup diagonal scaling
IF(present(diag))THEN
  IF(associated(self%D))THEN
    CALL diag%add(0.d0,1.d0,self%D)
  ELSE
    CALL oft_abort('Matrix not initialized.','diagmat_assemble',__FILE__)
  END IF
END IF
DEBUG_STACK_POP
!---Common assembly tasks
END SUBROUTINE diagmat_assemble
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_zero
!------------------------------------------------------------------------------
!> Zero all entries in matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine diagmat_zero(self)
class(oft_diagmatrix), intent(inout) :: self
CALL self%D%set(0.d0)
end subroutine diagmat_zero
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_zero_rows
!------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] nrows Number of rows to zero
!! @param[in] irows Indices of rows to zero [nrows]
!! @param[in] iblock Row block (optional)
!! @param[in] keep_diag Keep diagonal entries
!------------------------------------------------------------------------------
subroutine diagmat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_diagmatrix), intent(inout) :: self
integer(i4), intent(in) :: nrows
integer(i4), intent(in) :: irows(nrows)
integer(i4), optional, intent(in) :: iblock
logical, optional, intent(in) :: keep_diag
call oft_abort('Invalid operation for matrix type','diagmat_zero_rows',__FILE__)
end subroutine diagmat_zero_rows
!------------------------------------------------------------------------------
! SUBROUTINE: diagmat_delete
!------------------------------------------------------------------------------
!> Delete matrix
!------------------------------------------------------------------------------
subroutine diagmat_delete(self)
class(oft_diagmatrix), intent(inout) :: self
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  NULLIFY(self%D)
END IF
self%nr=0; self%nrg=0
self%nc=0; self%ncg=0
end subroutine diagmat_delete
!------------------------------------------------------------------------------
! FUNCTION: sum_matrix_cast
!------------------------------------------------------------------------------
!> Cast a matrix object to a oft_sum_matrix
!!
!! The source matrix must be oft_sum_matrix or a child class, otherwise an error will be thrown.
!!
!! @param[out] self Pointer to cast sum_matrix
!! @param[in] source Source matrix to cast
!------------------------------------------------------------------------------
FUNCTION sum_matrix_cast(self,source) result(ierr)
class(oft_sum_matrix), pointer, intent(out) :: self
class(oft_matrix), target, intent(in) :: source
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(oft_sum_matrix)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION sum_matrix_cast
!------------------------------------------------------------------------------
! SUBROUTINE: sum_mat_apply_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine sum_mat_apply_real(self,a,b)
class(oft_sum_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
class(oft_vector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_real',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%apply(a,b)
CALL self%K%apply(a,vtmp)
CALL b%add(1.d0,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_apply_real
!------------------------------------------------------------------------------
! SUBROUTINE: sum_mat_apply_comp
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine sum_mat_apply_comp(self,a,b)
class(oft_sum_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%apply(a,b)
CALL self%K%apply(a,vtmp)
CALL b%add((1.d0,0.d0),CMPLX(self%alam,KIND=c8),vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_apply_comp
!------------------------------------------------------------------------------
! SUBROUTINE: sum_mat_applyt_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine sum_mat_applyt_real(self,a,b)
class(oft_sum_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
class(oft_vector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_real',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%applyt(a,b)
CALL self%K%applyt(a,vtmp)
CALL b%add(1.d0,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_applyt_real
!------------------------------------------------------------------------------
! SUBROUTINE: sum_mat_applyt_comp
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine sum_mat_applyt_comp(self,a,b)
class(oft_sum_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a
class(oft_cvector), intent(inout) :: b
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%applyt(a,b)
CALL self%K%applyt(a,vtmp)
CALL b%add((1.d0,0.d0),CMPLX(self%alam,KIND=c8),vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_applyt_comp
!------------------------------------------------------------------------------
! SUBROUTINE: sum_mat_assemble
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!!
!! @param[in,out] diag Diagonal entries of matrix [nr] (optional)
!------------------------------------------------------------------------------
subroutine sum_mat_assemble(self,diag)
class(oft_sum_matrix), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag
DEBUG_STACK_PUSH
!---
if(present(diag))then
  !if(associated(self%D))then
  !  call diag%add(0.d0,1.d0,self%D)
  !else
    if(associated(self%D))call self%D%delete
    call diag%new(self%D)
    !---Get diagonal matrix values
    CALL self%J%assemble(self%D)
    CALL self%K%assemble(self%D)
    !---
    CALL self%D%add(0.d0,1.d0,self%J%D,self%alam,self%K%D)
  !end if
else
  !---Assemble matrices
  CALL self%J%assemble
  CALL self%K%assemble
end if
self%nr=self%J%nr; self%nrg=self%J%nrg
self%nc=self%J%nc; self%ncg=self%J%ncg
DEBUG_STACK_POP
end subroutine sum_mat_assemble
!---------------------------------------------------------------------------
! SUBROUTINE: mf_mat_apply_real
!---------------------------------------------------------------------------
!> Compute Jacobian approximation
!!
!! @param[in] a Source field
!! @param[in,out] b F(a)
!---------------------------------------------------------------------------
subroutine mf_mat_apply_real(self,a,b)
class(oft_mf_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8) :: eps,bb,aa,a1
aa=a%dot(a)
IF(aa==0.d0)THEN
  CALL b%set(0.d0)
  RETURN
ELSE
  IF(ASSOCIATED(self%utyp))THEN
    bb=a%dot(self%u0)
    CALL self%tmp%add(0.d0,1.d0,a)
    CALL self%tmp%mult(self%utyp)
    a1 = self%tmp%norm(1)
    eps=self%b0*MAX(ABS(bb),a1)*SIGN(1.d0,bb)/aa
  ELSE
    bb=SQRT(self%u0%dot(self%u0))
    eps=self%b0*SQRT(1.d0+bb)/aa
  END IF
END IF
!---
CALL self%tmp%add(0.d0,1.d0,self%u0,eps,a)
CALL self%f%apply(self%tmp,b)
CALL b%add(1.d0/eps,-1.d0/eps,self%f0)
end subroutine mf_mat_apply_real
!---------------------------------------------------------------------------
! SUBROUTINE: mf_mat_setup
!---------------------------------------------------------------------------
!> Setup matrix-free Jacobian operator
!!
!! @param[in] a Field defining domain and range spaces
!! @param[in] f Non-linear function defining the Jacobian
!---------------------------------------------------------------------------
subroutine mf_mat_setup(self,a,f,utyp)
class(oft_mf_matrix), intent(inout) :: self
class(oft_vector), intent(inout) :: a
class(oft_matrix), target, intent(in) :: f
class(oft_vector), optional, intent(inout) :: utyp
self%f=>f
CALL a%new(self%u0)
CALL a%new(self%f0)
CALL a%new(self%tmp)
CALL a%new(self%utyp)
IF(PRESENT(utyp))THEN
  CALL self%utyp%add(0.d0,1.d0,utyp)
ELSE
  CALL self%utyp%set(1.d0)
END IF
end subroutine mf_mat_setup
!---------------------------------------------------------------------------
! SUBROUTINE mf_mat_update
!---------------------------------------------------------------------------
!> Update linearization point
!!
!! @param[in] a New linearization point
!---------------------------------------------------------------------------
subroutine mf_mat_update(self,a)
class(oft_mf_matrix), intent(inout) :: self
class(oft_vector), intent(inout) :: a
integer(i4) :: i,nslice,ierr
real(r8) :: umin
real(r8), pointer :: vals(:)
CALL self%u0%add(0.d0,1.d0,a)
CALL self%f%apply(self%u0,self%f0)
END SUBROUTINE mf_mat_update
!---------------------------------------------------------------------------
! SUBROUTINE: mf_mat_delete
!---------------------------------------------------------------------------
!> Cleanup internal storage and reset defaults
!---------------------------------------------------------------------------
subroutine mf_mat_delete(self)
class(oft_mf_matrix), intent(inout) :: self
IF(ASSOCIATED(self%u0))THEN
  CALL self%u0%delete
  CALL self%f0%delete
  CALL self%tmp%delete
  CALL self%utyp%delete
  DEALLOCATE(self%u0,self%f0,self%tmp,self%utyp)
END IF
NULLIFY(self%f)
end subroutine mf_mat_delete
end module oft_deriv_matrices
