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
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix, oft_cmatrix
IMPLICIT NONE
#include "local.h"
private
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_noop_matrix
contains
  procedure :: apply_real => noop_apply_real
  ! procedure :: apply_complex => noop_apply_comp
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
!> Needs docs
!------------------------------------------------------------------------------
type, public, extends(oft_cmatrix) :: oft_noop_cmatrix
contains
  procedure :: apply_real => cnoop_apply_real
  procedure :: apply_complex => cnoop_apply_comp
  procedure :: applyt_real => cnoop_apply_real
  procedure :: applyt_complex => cnoop_apply_comp
  procedure :: set_values => cnoop_set_values
  procedure :: add_values => cnoop_add_values
  procedure :: atomic_add_values => cnoop_add_values
  procedure :: assemble => cnoop_assemble
  procedure :: zero => cnoop_zero
  procedure :: zero_rows => cnoop_zero_rows
  procedure :: delete => cnoop_delete
end type oft_noop_cmatrix
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
!> Sum matrix class M = beta*J + alam*K
!------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_sum_matrix
  real(r8) :: alam = 1.d0 !< Lambda factor
  real(r8) :: beta = 1.d0 !< Beta factor
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
!> Sum matrix class M = beta*J + alam*K
!------------------------------------------------------------------------------
type, public, extends(oft_noop_cmatrix) :: oft_sum_cmatrix
  complex(c8) :: alam = (1.d0,0.d0) !< Lambda factor
  complex(c8) :: beta = (1.d0,0.d0) !< Beta factor
  class(oft_matrix), pointer :: rJ => NULL() !< J matrix
  class(oft_matrix), pointer :: rK => NULL() !< K matrix
  class(oft_cmatrix), pointer :: cJ => NULL() !< J matrix
  class(oft_cmatrix), pointer :: cK => NULL() !< K matrix
contains
  procedure :: apply_real => sum_cmat_apply_real
  procedure :: apply_complex => sum_cmat_apply_comp
  procedure :: applyt_real => sum_cmat_applyt_real
  procedure :: applyt_complex => sum_cmat_applyt_comp
  !> Complete matrix assembly
  procedure :: assemble => sum_cmat_assemble
end type oft_sum_cmatrix
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
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_apply_real(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','noop_apply_real',__FILE__)
end subroutine noop_apply_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_apply_comp(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','noop_apply_comp',__FILE__)
end subroutine noop_apply_comp
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self^T * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_applyt(self,a,b)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','noop_applyt',__FILE__)
end subroutine noop_applyt
!------------------------------------------------------------------------------
!> Set values of a matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
call oft_abort('Invalid operation for matrix type','noop_set_values',__FILE__)
end subroutine noop_set_values
!------------------------------------------------------------------------------
!> Add values to a matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
call oft_abort('Invalid operation for matrix type','noop_add_values',__FILE__)
end subroutine noop_add_values
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_assemble(self,diag)
class(oft_noop_matrix), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
call oft_abort('No matrix type specified','noop_assemble',__FILE__)
end subroutine noop_assemble
!------------------------------------------------------------------------------
!> Zero all entries in matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_zero(self)
class(oft_noop_matrix), intent(inout) :: self
call oft_abort('No matrix type specified','noop_zero',__FILE__)
end subroutine noop_zero
!------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine noop_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_noop_matrix), intent(inout) :: self
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
call oft_abort('No matrix type specified','noop_zero_rows',__FILE__)
end subroutine noop_zero_rows
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
!> Apply the matrix to a field
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_apply_real(self,a,b)
class(oft_noop_cmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','cnoop_apply_real',__FILE__)
end subroutine cnoop_apply_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_apply_comp(self,a,b)
class(oft_noop_cmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','cnoop_apply_comp',__FILE__)
end subroutine cnoop_apply_comp
!------------------------------------------------------------------------------
!> Apply the matrix to a field
!!
!! b = selfT * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_applyt(self,a,b)
class(oft_noop_cmatrix), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< Result of matrix product
call oft_abort('No matrix type specified','cnoop_applyt',__FILE__)
end subroutine cnoop_applyt
!------------------------------------------------------------------------------
!> Set values of a matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_noop_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
complex(c8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
call oft_abort('Invalid operation for matrix type','cnoop_set_values',__FILE__)
end subroutine cnoop_set_values
!------------------------------------------------------------------------------
!> Add values to a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_noop_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
complex(c8), intent(in) :: b(n,m) !< Values to add [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
call oft_abort('Invalid operation for matrix type','cnoop_add_values',__FILE__)
end subroutine cnoop_add_values
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_assemble(self,diag)
class(oft_noop_cmatrix), intent(inout) :: self
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
call oft_abort('No matrix type specified','cnoop_assemble',__FILE__)
end subroutine cnoop_assemble
!------------------------------------------------------------------------------
!> Zero all entries in matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine cnoop_zero(self)
class(oft_noop_cmatrix), intent(inout) :: self
call oft_abort('No matrix type specified','cnoop_zero',__FILE__)
end subroutine cnoop_zero
!------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine cnoop_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_noop_cmatrix), intent(inout) :: self
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
call oft_abort('No matrix type specified','cnoop_zero_rows',__FILE__)
end subroutine cnoop_zero_rows
!------------------------------------------------------------------------------
!> Delete matrix
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!------------------------------------------------------------------------------
subroutine cnoop_delete(self)
class(oft_noop_cmatrix), intent(inout) :: self
call oft_warn('Finalizing general matrix, this may indicate an error.')
end subroutine cnoop_delete
!------------------------------------------------------------------------------
!> Needs docs
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
!> Cast a matrix object to a oft_diagmatrix
!!
!! The source matrix must be @ref oft_diagmatrix or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION diagmatrix_cast(self,source) result(success)
class(oft_diagmatrix), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_matrix), target, intent(in) :: source !< Source matrix to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_diagmatrix)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION diagmatrix_cast
!------------------------------------------------------------------------------
!> Apply the matrix to a field
!!
!! b = self * a
!------------------------------------------------------------------------------
subroutine diagmat_apply_real(self,a,b)
class(oft_diagmatrix), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of matrix product
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','diagmat_apply_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','diagmat_apply_real',__FILE__)
call b%add(0.d0,1.d0,a)
call b%mult(self%D)
DEBUG_STACK_POP
end subroutine diagmat_apply_real
!------------------------------------------------------------------------------
!> Apply the matrix to a field
!!
!! b = self * a
!------------------------------------------------------------------------------
subroutine diagmat_apply_comp(self,a,b)
class(oft_diagmatrix), intent(inout) :: self
class(oft_cvector), target, intent(inout) :: a !< Source field
class(oft_cvector), intent(inout) :: b !< 
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','diagmat_apply_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','diagmat_apply_comp',__FILE__)
call b%add((0.d0,0.d0),(1.d0,0.d0),a)
call b%mult(self%D)
DEBUG_STACK_POP
end subroutine diagmat_apply_comp
!------------------------------------------------------------------------------
!> Set values of a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine diagmat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_diagmatrix), intent(inout) :: self
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
call oft_abort('Invalid operation for matrix type','diagmat_set_values',__FILE__)
end subroutine diagmat_set_values
!------------------------------------------------------------------------------
!> Add values to a matrix.
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices
!------------------------------------------------------------------------------
subroutine diagmat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_diagmatrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
call oft_abort('Invalid operation for matrix type','diagmat_add_values',__FILE__)
end subroutine diagmat_add_values
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!------------------------------------------------------------------------------
SUBROUTINE diagmat_assemble(self,diag)
CLASS(oft_diagmatrix), INTENT(inout) :: self !< Matrix object
CLASS(oft_vector), OPTIONAL, TARGET, INTENT(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
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
!> Zero all entries in matrix
!------------------------------------------------------------------------------
subroutine diagmat_zero(self)
class(oft_diagmatrix), intent(inout) :: self
CALL self%D%set(0.d0)
end subroutine diagmat_zero
!------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!------------------------------------------------------------------------------
subroutine diagmat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_diagmatrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
call oft_abort('Invalid operation for matrix type','diagmat_zero_rows',__FILE__)
end subroutine diagmat_zero_rows
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
!> Cast a matrix object to a sum_matrix_cast
!!
!! The source matrix must be @ref sum_matrix_cast or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION sum_matrix_cast(self,source) result(success)
class(oft_sum_matrix), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_matrix), target, intent(in) :: source !< Source matrix to cast
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_sum_matrix)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION sum_matrix_cast
!------------------------------------------------------------------------------
!> Compute matrix vector product
!!
!! b = (self%beta*self%J + self%alam*self%K) * a
!------------------------------------------------------------------------------
subroutine sum_mat_apply_real(self,a,b)
class(oft_sum_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_vector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_real',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%apply(a,b)
CALL self%K%apply(a,vtmp)
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_apply_real
!------------------------------------------------------------------------------
!> Compute matrix vector product (complex)
!!
!! b = (self%beta*self%J + self%alam*self%K) * a
!------------------------------------------------------------------------------
subroutine sum_mat_apply_comp(self,a,b)
class(oft_sum_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%apply(a,b)
CALL self%K%apply(a,vtmp)
CALL b%add(self%beta*(1.d0,0.d0),self%alam*(1.d0,0.d0),vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_apply_comp
!------------------------------------------------------------------------------
!> Compute matrix vector product for matrix transpose
!!
!! b = self%beta*self%J^T + self%alam*self%K^T * a
!------------------------------------------------------------------------------
subroutine sum_mat_applyt_real(self,a,b)
class(oft_sum_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_vector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_real',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%applyt(a,b)
CALL self%K%applyt(a,vtmp)
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_applyt_real
!------------------------------------------------------------------------------
!> Compute matrix vector product for matrix transpose (complex)
!!
!! b = self%beta*self%J^T + self%alam*self%K^T * a
!------------------------------------------------------------------------------
subroutine sum_mat_applyt_comp(self,a,b)
class(oft_sum_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
CALL self%J%applyt(a,b)
CALL self%K%applyt(a,vtmp)
CALL b%add(self%beta*(1.d0,0.d0),self%alam*(1.d0,0.d0),vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_mat_applyt_comp
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!------------------------------------------------------------------------------
subroutine sum_mat_assemble(self,diag)
class(oft_sum_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
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
    CALL self%D%add(0.d0,self%beta,self%J%D,self%alam,self%K%D)
    CALL diag%add(0.d0,1.d0,self%D)
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
!------------------------------------------------------------------------------
!> Compute matrix vector product (real)
!!
!! b = (self%beta*self%J + self%alam*self%K) * a
!------------------------------------------------------------------------------
subroutine sum_cmat_apply_real(self,a,b)
class(oft_sum_cmatrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: rtmp
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_real',__FILE__)
NULLIFY(rtmp,vtmp)
CALL a%new(vtmp)
IF(ASSOCIATED(self%rJ))THEN
  CALL a%new(rtmp)
  CALL self%rJ%apply(a,rtmp)
  CALL b%add((0.d0,0.d0),(1.d0,0.d0),rtmp)
  CALL rtmp%delete()
  DEALLOCATE(rtmp)
ELSE IF(ASSOCIATED(self%cJ))THEN
  CALL self%cJ%apply(a,b)
END IF
IF(ASSOCIATED(self%rK))THEN
  CALL a%new(rtmp)
  CALL self%rK%apply(a,rtmp)
  CALL vtmp%add((0.d0,0.d0),(1.d0,0.d0),rtmp)
  CALL rtmp%delete()
  DEALLOCATE(rtmp)
ELSE IF(ASSOCIATED(self%cK))THEN
  CALL self%cK%apply(a,vtmp)
END IF
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_cmat_apply_real
!------------------------------------------------------------------------------
!> Compute matrix vector product
!!
!! b = (self%beta*self%J + self%alam*self%K) * a
!------------------------------------------------------------------------------
subroutine sum_cmat_apply_comp(self,a,b)
class(oft_sum_cmatrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_apply_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_apply_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
IF(ASSOCIATED(self%rJ))THEN
  CALL self%rJ%apply(a,b)
ELSE IF(ASSOCIATED(self%cJ))THEN
  CALL self%cJ%apply(a,b)
END IF
IF(ASSOCIATED(self%rK))THEN
  CALL self%rK%apply(a,vtmp)
ELSE IF(ASSOCIATED(self%cK))THEN
  CALL self%cK%apply(a,vtmp)
END IF
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_cmat_apply_comp
!------------------------------------------------------------------------------
!> Compute matrix vector product for matrix transpose (real)
!!
!! b = self%beta*self%J^T + self%alam*self%K^T * a
!------------------------------------------------------------------------------
subroutine sum_cmat_applyt_real(self,a,b)
class(oft_sum_cmatrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_vector), pointer :: rtmp
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_real',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_real',__FILE__)
NULLIFY(vtmp,rtmp)
CALL a%new(vtmp)
IF(ASSOCIATED(self%rJ))THEN
  CALL a%new(rtmp)
  CALL self%rJ%apply(a,rtmp)
  CALL b%add((0.d0,0.d0),(1.d0,0.d0),rtmp)
  CALL rtmp%delete()
  DEALLOCATE(rtmp)
ELSE IF(ASSOCIATED(self%cJ))THEN
  CALL self%cJ%apply(a,b)
END IF
IF(ASSOCIATED(self%rK))THEN
  CALL a%new(rtmp)
  CALL self%rK%apply(a,rtmp)
  CALL vtmp%add((0.d0,0.d0),(1.d0,0.d0),rtmp)
  CALL rtmp%delete()
  DEALLOCATE(rtmp)
ELSE IF(ASSOCIATED(self%cK))THEN
  CALL self%cK%apply(a,vtmp)
END IF
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_cmat_applyt_real
!------------------------------------------------------------------------------
!> Compute matrix vector product for matrix transpose
!!
!! b = self%beta*self%J^T + self%alam*self%K^T * a
!------------------------------------------------------------------------------
subroutine sum_cmat_applyt_comp(self,a,b)
class(oft_sum_cmatrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Vector object
class(oft_cvector), intent(inout) :: b !< Result vector
class(oft_cvector), pointer :: vtmp
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','sum_mat_applyt_comp',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','sum_mat_applyt_comp',__FILE__)
NULLIFY(vtmp)
CALL a%new(vtmp)
IF(ASSOCIATED(self%rJ))THEN
  CALL self%rJ%apply(a,b)
ELSE IF(ASSOCIATED(self%cJ))THEN
  CALL self%cJ%apply(a,b)
END IF
IF(ASSOCIATED(self%rK))THEN
  CALL self%rK%apply(a,vtmp)
ELSE IF(ASSOCIATED(self%cK))THEN
  CALL self%cK%apply(a,vtmp)
END IF
CALL b%add(self%beta,self%alam,vtmp)
CALL vtmp%delete()
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine sum_cmat_applyt_comp
!------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!------------------------------------------------------------------------------
subroutine sum_cmat_assemble(self,diag)
class(oft_sum_cmatrix), intent(inout) :: self !< Matrix object
class(oft_cvector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
class(oft_vector), pointer :: rdiag
DEBUG_STACK_PUSH
!---
NULLIFY(rdiag)
if(present(diag))then
  if(associated(self%D))call self%D%delete
  call diag%new(self%D)
  !---Get diagonal matrix values
  IF(ASSOCIATED(self%cJ))THEN
    CALL self%cJ%assemble(diag)
    CALL self%D%add((0.d0,0.d0),self%beta,diag)
  ELSE IF(ASSOCIATED(self%rJ))THEN
    CALL diag%new(rdiag)
    CALL self%rJ%assemble(rdiag)
    CALL self%D%add((0.d0,0.d0),self%beta,rdiag)
    CALL rdiag%delete()
    DEALLOCATE(rdiag)
  END IF
  IF(ASSOCIATED(self%cK))THEN
    CALL self%cK%assemble(diag)
    CALL self%D%add((1.d0,0.d0),self%alam,diag)
  ELSE IF(ASSOCIATED(self%rK))THEN
    CALL diag%new(rdiag)
    CALL self%rK%assemble(rdiag)
    CALL self%D%add((1.d0,0.d0),self%alam,rdiag)
    CALL rdiag%delete()
    DEALLOCATE(rdiag)
  END IF
  !---
  CALL diag%add((0.d0,0.d0),(1.d0,0.d0),self%D)
else
  !---Assemble matrices
  IF(ASSOCIATED(self%rJ))THEN
    CALL self%rJ%assemble()
  ELSE IF(ASSOCIATED(self%cJ))THEN
    CALL self%cJ%assemble()
  END IF
  IF(ASSOCIATED(self%rK))THEN
    CALL self%rK%assemble()
  ELSE IF(ASSOCIATED(self%cK))THEN
    CALL self%cK%assemble()
  END IF
end if
IF(ASSOCIATED(self%rJ))THEN
  self%nr=self%rJ%nr; self%nrg=self%rJ%nrg
  self%nc=self%rJ%nc; self%ncg=self%rJ%ncg
ELSE IF(ASSOCIATED(self%cJ))THEN
  self%nr=self%cJ%nr; self%nrg=self%cJ%nrg
  self%nc=self%cJ%nc; self%ncg=self%cJ%ncg
END IF
DEBUG_STACK_POP
end subroutine sum_cmat_assemble
!---------------------------------------------------------------------------
!> Compute matrix vector product
!!
!! b = (self%f(a) - self%f0)/eps
!---------------------------------------------------------------------------
subroutine mf_mat_apply_real(self,a,b)
class(oft_mf_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Vector object
class(oft_vector), intent(inout) :: b !< Result vector
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
!> Setup matrix-free Jacobian operator
!---------------------------------------------------------------------------
subroutine mf_mat_setup(self,a,f,utyp)
class(oft_mf_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), intent(inout) :: a !< Vector defining domain and range spaces
class(oft_matrix), target, intent(in) :: f !< Non-linear function defining the Jacobian
class(oft_vector), optional, intent(inout) :: utyp !< Vector of "typical sizes" (optional)
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
!> Update linearization point for matrix-free Jacobian
!---------------------------------------------------------------------------
subroutine mf_mat_update(self,a)
class(oft_mf_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), intent(inout) :: a !< New linearization point
integer(i4) :: i,nslice,ierr
real(r8) :: umin
real(r8), pointer :: vals(:)
CALL self%u0%add(0.d0,1.d0,a)
CALL self%f%apply(self%u0,self%f0)
END SUBROUTINE mf_mat_update
!---------------------------------------------------------------------------
!> Cleanup internal storage and reset defaults
!---------------------------------------------------------------------------
subroutine mf_mat_delete(self)
class(oft_mf_matrix), intent(inout) :: self !< Matrix object
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
