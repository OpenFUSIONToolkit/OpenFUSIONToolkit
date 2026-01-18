!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_petsc_la.F90
!
!> @defgroup doxy_oft_petsc_la PETSc backend
!! PETSc vector and matrix backend
!! @ingroup doxy_oft_lin_alg
!
!> PETSc vector implementation
!!
!! @authors Chris Hansen
!! @date January 2013
!! @ingroup doxy_oft_petsc_la
!---------------------------------------------------------------------------------
MODULE oft_petsc_la
USE oft_base
USE oft_stitching, ONLY: oft_global_dp
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_cvector_ptr, &
  oft_map, oft_matrix, oft_cmatrix
#ifdef HAVE_PETSC
USE petscvec
USE petscmat
#undef Mat
#undef IS
#undef Vec
IMPLICIT NONE
#include "local.h"
private
!---------------------------------------------------------------------------------
!> PETSc vector implementation
!---------------------------------------------------------------------------------
type, public, extends(oft_vector) :: oft_petsc_vector
  TYPE(tvec) :: v !< PETSc vector object
  TYPE(tvec) :: vloc !< PETSc local vector object
  TYPE(tis) :: lis !< PETSc local to global mapping
  TYPE(tvecscatter) :: get !< PETSc local to global scatter context
  LOGICAL :: loc_store = .FALSE. !< Flag indicating pending local->global sync
  LOGICAL :: loc_current = .FALSE. !< Flag indicating local/global vectors are consistent
contains
  procedure :: new_real => vec_new_vec
  procedure :: new_complex => vec_new_cvec
  !> Set all elements to a scalar
  procedure :: set => vec_set
  !> Get local portion of field
  procedure :: get_local => vec_get_local
  !> Restore local portion of field
  procedure :: restore_local => vec_restore_local
  !> Get owned portion of field
  procedure :: get_slice => vec_get_slice
  !> Restore owned portion of field
  procedure :: restore_slice => vec_restore_slice
  !> Add vectors
  procedure :: add => vec_add
  !> Multiply fields element by element
  procedure :: mult => vec_mult
  !> Scale vector by a scalar
  procedure :: scale => vec_scale
  procedure :: dot_real => vec_dot_vec
  procedure :: dot_complex => vec_dot_cvec
  procedure :: mdot_real => vec_mdot_vec
  procedure :: mdot_complex => vec_mdot_cvec
  !> Sum reduction over vector
  procedure :: sum => vec_sum
  !> Norm of vector
  procedure :: norm => vec_norm
  !> Perform global stitching
  procedure :: stitch => vec_stitch
  !> Delete vector
  procedure :: delete => vec_delete
end type oft_petsc_vector
! !---------------------------------------------------------------------------------
! !> PETSc vector implementation (complex)
! !---------------------------------------------------------------------------------
! type, public, extends(oft_cvector) :: oft_petsc_cvector
!   TYPE(tvec) :: v !< PETSc vector object
!   TYPE(tvec) :: vloc !< PETSc local vector object
!   TYPE(tis) :: lis !< PETSc local to global mapping
!   TYPE(tvecscatter) :: get !< PETSc local to global scatter context
!   LOGICAL :: loc_store = .FALSE. !< Flag indicating pending local->global sync
!   LOGICAL :: loc_current = .FALSE. !< Flag indicating local/global vectors are consistent
! contains
!   procedure :: new_real => cvec_new_vec
!   procedure :: new_complex => cvec_new_cvec
!   !> Set all elements to a scalar
!   procedure :: set => cvec_set
!   !> Get local portion of field
!   procedure :: get_local => cvec_get_local
!   !> Restore local portion of field
!   procedure :: restore_local => cvec_restore_local
!   !> Get owned portion of field
!   procedure :: get_slice => cvec_get_slice
!   !> Restore owned portion of field
!   procedure :: restore_slice => cvec_restore_slice
!   procedure :: add_real => cvec_add_vec
!   procedure :: add_complex => cvec_add_cvec
!   procedure :: mult_real => cvec_mult_vec
!   procedure :: mult_complex => cvec_mult_cvec
!   procedure :: dot_real => cvec_dot_vec
!   procedure :: dot_complex => cvec_dot_cvec
!   procedure :: mdot_real => cvec_mdot_vec
!   procedure :: mdot_complex => cvec_mdot_cvec
!   !> Sum reduction over vector
!   procedure :: sum => cvec_sum
!   !> Norm of vector
!   procedure :: norm => cvec_norm
!   !> Perform global stitching
!   procedure :: stitch => cvec_stitch
!   !> Delete vector
!   procedure :: delete => cvec_delete
! end type oft_petsc_cvector
!---------------------------------------------------------------------------------
!> PETSc matrix implementation
!---------------------------------------------------------------------------------
type, public, extends(oft_matrix) :: oft_petsc_matrix
  TYPE(tmat) :: M !< PETSc matrix object
  TYPE(tvec) :: Md !< PETSc matrix diagonal object (unused)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<21)
  INTEGER(petsc_addr) :: r_lis !< PETSc Local to Global mapping for rows
  INTEGER(petsc_addr) :: c_lis !< PETSc Local to Global mapping for columns
#else
  TYPE(tislocaltoglobalmapping) :: r_lis !< PETSc Local to Global mapping for rows
  TYPE(tislocaltoglobalmapping) :: c_lis !< PETSc Local to Global mapping for columns
#endif
contains
  procedure :: apply_real => mat_apply_vec
  procedure :: apply_complex => mat_apply_cvec
  procedure :: applyt_real => mat_applyt_vec
  procedure :: applyt_complex => mat_applyt_cvec
  !> Set values of the matrix
  procedure :: set_values => mat_set_values
  !> Add values to the matrix
  procedure :: add_values => mat_add_values
  !> Add values atomically to the matrix
  procedure :: atomic_add_values => mat_add_values
  !> Complete matrix assembly
  procedure :: assemble => mat_assemble
  !> Zero all elements
  procedure :: zero => mat_zero
  !> Zero all elements in a given row
  procedure :: zero_rows => mat_zero_rows
  !> Delete matrix
  procedure :: delete => mat_delete
end type oft_petsc_matrix
! !---------------------------------------------------------------------------------
! !> PETSc matrix implementation (complex)
! !---------------------------------------------------------------------------------
! type, public, extends(oft_cmatrix) :: oft_petsc_cmatrix
!   TYPE(tmat) :: M !< PETSc matrix object
!   TYPE(tvec) :: Md !< PETSc matrix diagonal object (unused)
!   INTEGER(petsc_addr) :: r_lis !< PETSc Local to Global mapping for rows
!   INTEGER(petsc_addr) :: c_lis !< PETSc Local to Global mapping for columns
! contains
!   procedure :: apply_real => mat_apply_vec
!   procedure :: apply_complex => mat_apply_cvec
!   procedure :: applyt_real => mat_applyt_vec
!   procedure :: applyt_complex => mat_applyt_cvec
!   !> Set values of the matrix
!   procedure :: set_values => mat_set_values
!   !> Add values to the matrix
!   procedure :: add_values => mat_add_values
!   !> Complete matrix assembly
!   procedure :: assemble => mat_assemble
!   !> Zero all elements
!   procedure :: zero => mat_zero
!   !> Zero all elements in a given row
!   procedure :: zero_rows => mat_zero_rows
!   !> Delete matrix
!   procedure :: delete => mat_delete
! end type oft_petsc_cmatrix
!---Declare public entities
public oft_petsc_vector_cast, oft_petsc_matrix_cast
contains
!------------------------------------------------------------------------------
!> Cast a vector object to a oft_petsc_vector
!!
!! The source vector must be @ref oft_petsc_vector or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_petsc_vector_cast(self,source) result(success)
class(oft_petsc_vector), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_vector), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_petsc_vector)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION oft_petsc_vector_cast
!---------------------------------------------------------------------------------
!> Create a new vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine vec_new_vec(self,new)
class(oft_petsc_vector), intent(in) :: self !< Vector object
class(oft_vector), pointer, intent(out) :: new !< New vector
integer(i4) :: ierr
DEBUG_STACK_PUSH
ALLOCATE(oft_petsc_vector::new)
SELECT TYPE(this=>new)
TYPE IS(oft_petsc_vector)
  this%n=self%n
  this%ng=self%ng
  this%nblocks=self%nblocks
  this%map=>self%map
  this%stitch_info=>self%stitch_info
  this%lis=self%lis
  this%get=self%get
  !---
  CALL VecDuplicate(self%v,this%v,ierr)
  CALL VecSet(this%v,0.d0,ierr)
  CALL VecDuplicate(self%vloc,this%vloc,ierr)
  CALL VecSet(this%vloc,0.d0,ierr)
  this%loc_current=.TRUE.
  this%loc_store=.FALSE.
CLASS DEFAULT
  CALL oft_abort('Failure to allocate vector.','vec_new_vec',__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine vec_new_vec
!---------------------------------------------------------------------------------
!> Create a new complex vector as a bare copy of `self`
!---------------------------------------------------------------------------------
subroutine vec_new_cvec(self,new)
class(oft_petsc_vector), intent(in) :: self !< Vector object
class(oft_cvector), pointer, intent(out) :: new !< New vector
CALL oft_abort('Complex LA not yet supported with PETSc.','vec_new_cvec',__FILE__)
! DEBUG_STACK_PUSH
! ALLOCATE(oft_petsc_vector::new)
! SELECT TYPE(this=>new)
! TYPE IS(oft_petsc_vector)
!   this%n=self%n
!   this%ng=self%ng
!   this%nblocks=self%nblocks
!   this%map=>self%map
!   this%stitch_info=>self%stitch_info
!   this%lis=self%lis
!   this%get=self%get
!   !---
!   CALL VecDuplicate(self%v,this%v,ierr)
!   CALL VecSet(this%v,0.d0,ierr)
!   CALL VecDuplicate(self%vloc,this%vloc,ierr)
!   CALL VecSet(this%vloc,0.d0,ierr)
!   this%loc_current=.TRUE.
!   this%loc_store=.FALSE.
! CLASS DEFAULT
!   CALL oft_abort('Failure to allocate vector.','vec_new',__FILE__)
! END SELECT
! DEBUG_STACK_POP
end subroutine vec_new_cvec
!---------------------------------------------------------------------------------
!> Set all elements to a scalar
!---------------------------------------------------------------------------------
subroutine vec_set(self,alpha,iblock,random)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: alpha !< Updated vector value
integer(i4), optional, intent(in) :: iblock !< Vector sub-block to act on
logical, optional, intent(in) :: random !< Set to random number, if true alpha is ignored (optional)
logical :: random_flag
integer(i4) :: i,ierr
real(r8), pointer, dimension(:) :: array
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
random_flag=.FALSE.
IF(PRESENT(random))random_flag=random
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_set',__FILE__)
  !---
  IF(random_flag)THEN
    NULLIFY(array)
    CALL self%get_local(array,iblock)
    CALL oft_random_number(array,self%map(iblock)%n)
    CALL self%restore_local(array,iblock)
    DEALLOCATE(array)
  ELSE
    NULLIFY(array)
    CALL self%get_slice(array,iblock)
    array=alpha
    CALL self%restore_slice(array,iblock)
    DEALLOCATE(array)
  END IF
ELSE
  !---
  IF(random_flag)THEN
    NULLIFY(array)
    CALL self%get_local(array)
    CALL oft_random_number(array,self%n)
    CALL self%restore_local(array)
    DEALLOCATE(array)
  ELSE
    CALL VecSet(self%v,alpha,ierr)
  END IF
END IF
DEBUG_STACK_POP
end subroutine vec_set
!---------------------------------------------------------------------------------
!> Get values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine vec_get_slice(self,array,iblock)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), pointer, intent(inout) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,j,ib,n,nslice,ierr
real(r8), pointer :: x_array(:)
logical :: dist
DEBUG_STACK_PUSH
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requestd block exceeds nblocks.','vec_get_slice',__FILE__)
  !---
  nslice=self%map(iblock)%nslice
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(nslice))
  call VecGetArrayF90(self%v,x_array,ierr)
  DO i=1,nslice
    array(i)=x_array(self%map(iblock)%soffset+i)
  END DO
  call VecRestoreArrayF90(self%v,x_array,ierr)
  !CALL VecScatterBegin(self%get,self%v,self%vloc,INSERT_VALUES,SCATTER_FORWARD)
  !CALL VecGetValues(self%v,nslice,self%map(iblock)%slice,array,ierr)
  !call VecGetArray(self%vloc,x_array,i_x,ierr)
ELSE
  !---
  nslice=0
  DO i=1,self%nblocks
    nslice=nslice+self%map(i)%nslice
  END DO
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(nslice))
  call VecGetArrayF90(self%v,x_array,ierr)
  nslice=0
  do i=1,self%nblocks
    !CALL VecGetValues(self%v,self%map(i)%nslice, &
    !self%map(iblock)%slice,array(nslice+1:nslice+self%map(i)%nslice),ierr)
    DO j=1,self%map(i)%nslice
      array(nslice+j)=x_array(nslice+j)
    END DO
    nslice=nslice+self%map(i)%nslice
  end do
  call VecRestoreArrayF90(self%v,x_array,ierr)
END IF
DEBUG_STACK_POP
end subroutine vec_get_slice
!---------------------------------------------------------------------------------
!> Set/add values for locally-owned portion of vector (slice)
!---------------------------------------------------------------------------------
subroutine vec_restore_slice(self,array,iblock,wait)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: array(:) !< Slice values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: wait !< Wait to perform global sync?
integer(i4) :: i,j,n,nslice,ierr
real(r8), pointer :: x_array(:)
logical :: dist
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_slice',__FILE__)
  !---
  nslice=self%map(iblock)%nslice
  call VecGetArrayF90(self%v,x_array,ierr)
  DO i=1,nslice
    x_array(self%map(iblock)%soffset+i)=array(i)
  END DO
  call VecRestoreArrayF90(self%v,x_array,ierr)
ELSE
  !---
  nslice=0
  DO i=1,self%nblocks
    nslice=nslice+self%map(i)%nslice
  END DO
  call VecGetArrayF90(self%v,x_array,ierr)
  nslice=0
  do i=1,self%nblocks
    DO j=1,self%map(i)%nslice
      x_array(nslice+j)=array(nslice+j)
    END DO
    nslice=nslice+self%map(i)%nslice
  end do
  call VecRestoreArrayF90(self%v,x_array,ierr)
END IF
DEBUG_STACK_POP
end subroutine vec_restore_slice
!---------------------------------------------------------------------------------
!> Get local values from vector
!---------------------------------------------------------------------------------
subroutine vec_get_local(self,array,iblock)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), pointer, intent(inout) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to retrieve
integer(i4) :: i,nslice,ierr
real(r8), pointer :: x_array(:)
integer(i4), pointer :: inds(:)
logical :: dist
DEBUG_STACK_PUSH
!---
IF(.NOT.self%loc_current)THEN
  CALL VecScatterBegin(self%get,self%v,self%vloc,INSERT_VALUES,SCATTER_FORWARD,ierr)
  CALL VecScatterEnd(self%get,self%v,self%vloc,INSERT_VALUES,SCATTER_FORWARD,ierr)
  self%loc_current=.TRUE.
END IF
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_get_local',__FILE__)
  !---
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%map(iblock)%n))
  call VecGetArrayF90(self%vloc,x_array,ierr)
  DO i=1,self%map(iblock)%n
   array(i)=x_array(self%map(iblock)%offset+i)
  END DO
  call VecRestoreArrayF90(self%vloc,x_array,ierr)
ELSE
  !---
  IF(.NOT.ASSOCIATED(array))ALLOCATE(array(self%n))
  call VecGetArrayF90(self%vloc,x_array,ierr)
  DO i=1,self%n
    array(i)=x_array(i)
  END DO
  call VecRestoreArrayF90(self%vloc,x_array,ierr)
END IF
DEBUG_STACK_POP
end subroutine vec_get_local
!---------------------------------------------------------------------------------
!> Set/add local values to vector
!---------------------------------------------------------------------------------
subroutine vec_restore_local(self,array,iblock,add,wait)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: array(:) !< Local values
integer(i4), optional, intent(in) :: iblock !< Sub-block to restore
logical, optional, intent(in) :: add !< Add values instead of replace
logical, optional, intent(in) :: wait !< Wait to perform global sync
integer(i4) :: i,j,k,nslice,ierr
real(r8) :: val
real(r8), pointer :: x_array(:)
logical :: do_add,do_wait
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
do_add=.FALSE.
IF(PRESENT(add))do_add=add
do_wait=.FALSE.
IF(PRESENT(wait))do_wait=wait
!---
IF(PRESENT(iblock))THEN
  IF(iblock>self%nblocks)CALL oft_abort('Requested block exceeds nblocks.','vec_restore_local',__FILE__)
  !---
  !ALLOCATE(inds(self%map(iblock)%n))
  !inds=(/(self%map(iblock)%offset+i,i=1,self%map(iblock)%n)/)
  !CALL VecSetValues(self%vloc,self%map(iblock)%n,inds,array,ierr)
  IF(do_add)THEN
    call VecGetArrayF90(self%vloc,x_array,ierr)
    IF(.NOT.self%loc_store)x_array=0.d0
    DO i=1,self%map(iblock)%n
     x_array(self%map(iblock)%offset+i)=array(i)
    END DO
    call VecRestoreArrayF90(self%vloc,x_array,ierr)
    !DEALLOCATE(x_array)
    IF(.NOT.do_wait)THEN
      CALL VecScatterBegin(self%get,self%vloc,self%v,ADD_VALUES,SCATTER_REVERSE,ierr)
      CALL VecScatterEnd(self%get,self%vloc,self%v,ADD_VALUES,SCATTER_REVERSE,ierr)
      self%loc_store=.FALSE.
    ELSE
      self%loc_store=.TRUE.
    END IF
  ELSE
    call VecGetArrayF90(self%v,x_array,ierr)
    nslice=0
    DO k=1,iblock-1
      nslice=nslice+self%map(k)%nslice
    END DO
    DO i=1,self%map(iblock)%nslice
      j=self%map(iblock)%slice(i)
      x_array(nslice+i)=array(j)
    END DO
    call VecRestoreArrayF90(self%v,x_array,ierr)
  END IF
!  call VecRestoreArrayF90(self%vloc,x_array,ierr)
!  CALL VecScatterBegin(self%get,self%vloc,self%v,insert_flag,SCATTER_REVERSE,ierr)
!  CALL VecScatterEnd(self%get,self%vloc,self%v,insert_flag,SCATTER_REVERSE,ierr)
  !DEALLOCATE(inds)
ELSE
  !---
  !ALLOCATE(inds(self%n))
  !inds=(/(i,i=1,self%n)/)
  !CALL VecSetValues(self%vloc,self%n,inds,array,ierr)
  !call VecGetArrayF90(self%vloc,x_array,ierr)
  IF(do_add)THEN
    call VecGetArrayF90(self%vloc,x_array,ierr)
    x_array=0.d0
    DO i=1,self%n
      x_array(i)=array(i)
    END DO
    call VecRestoreArrayF90(self%vloc,x_array,ierr)
    CALL VecScatterBegin(self%get,self%vloc,self%v,ADD_VALUES,SCATTER_REVERSE,ierr)
    CALL VecScatterEnd(self%get,self%vloc,self%v,ADD_VALUES,SCATTER_REVERSE,ierr)
  ELSE
    call VecGetArrayF90(self%v,x_array,ierr)
    nslice=0
    DO k=1,self%nblocks
      DO i=1,self%map(k)%nslice
        j=self%map(k)%offset+self%map(k)%slice(i)
        x_array(nslice+i)=array(j)
      END DO
      nslice=nslice+self%map(k)%nslice
    END DO
    call VecRestoreArrayF90(self%v,x_array,ierr)
  END IF
!  call VecRestoreArrayF90(self%vloc,x_array,ierr)
!  CALL VecScatterBegin(self%get,self%vloc,self%v,insert_flag,SCATTER_REVERSE,ierr)
!  CALL VecScatterEnd(self%get,self%vloc,self%v,insert_flag,SCATTER_REVERSE,ierr)
  !DEALLOCATE(inds)
END IF
DEBUG_STACK_POP
end subroutine vec_restore_local
!---------------------------------------------------------------------------------
!> Add vectors
!!
!! self = \f$ \gamma \f$ self + \f$ \alpha \f$ a + \f$ \beta \f$ b
!---------------------------------------------------------------------------------
subroutine vec_add(self,gamma,alpha,a,beta,b)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: gamma !< Scale of source vector
real(r8), intent(in) :: alpha !< Scale of first vector
class(oft_vector), target, intent(inout) :: a !< First vector to add
real(r8), optional, intent(in) :: beta !< Scale of second vector (optional)
class(oft_vector), target, optional, intent(inout) :: b !< Second vector to add (optional)
class(oft_petsc_vector), pointer :: av,bv
integer(i4) :: i,ierr
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
if(present(beta).AND.present(b))then ! Add two vectors to the source
  IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a PETSc vector.','vec_add',__FILE__)
  IF(.NOT.oft_petsc_vector_cast(bv,b))CALL oft_abort('"b" is not a PETSc vector.','vec_add',__FILE__)
  CALL VecAXPBYPCZ(self%v,alpha,beta,gamma,av%v,bv%v,ierr)
else ! Add one vector to the source
  IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a PETSc vector.','vec_add',__FILE__)
  CALL VecAXPBY(self%v,alpha,gamma,av%v,ierr)
end if
DEBUG_STACK_POP
end subroutine vec_add
!---------------------------------------------------------------------------------
!> Elementwise multiplication with another vector
!!
!! \f$ self_i = self_i * a_i \f$
!---------------------------------------------------------------------------------
subroutine vec_mult(self,a,div_flag)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< vector for multiplication
logical, optional, intent(in) :: div_flag !< Divide instead of multiply?
class(oft_petsc_vector), pointer :: av
integer(i4) :: ierr
logical :: div
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
div=.FALSE.
if(present(div_flag))div=div_flag
!---
IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a PETSc vector.','vec_mult',__FILE__)
if(self%n/=av%n)call oft_abort('Vector sizes do not match.','vec_mult',__FILE__)
!---
if(div)then
  CALL VecPointwiseDivide(self%v,self%v,av%v,ierr)
else
  CALL VecPointwiseMult(self%v,self%v,av%v,ierr)
end if
DEBUG_STACK_POP
end subroutine vec_mult
!---------------------------------------------------------------------------------
!> Scale vector by a scalar
!---------------------------------------------------------------------------------
subroutine vec_scale(self,alpha)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8), intent(in) :: alpha !< Scale factor
integer(i4) :: i,ierr
DEBUG_STACK_PUSH
self%loc_current=.FALSE.
CALL VecScale(self%v,alpha,ierr)
DEBUG_STACK_POP
end subroutine vec_scale
!---------------------------------------------------------------------------------
!> Dot product with a vector
!---------------------------------------------------------------------------------
function vec_dot_vec(self,a) result(dot)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
class(oft_vector), target, intent(inout) :: a !< Second vector for dot product
real(r8) :: dot !< \f$ \sum_i self_i a_i \f$
class(oft_petsc_vector), pointer :: av
integer(i4) :: i,ierr
DEBUG_STACK_PUSH
IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a vector.','vec_dot_vec',__FILE__)
if(self%n/=av%n)call oft_abort('Vector lengths do not match.','vec_dot_vec',__FILE__)
CALL VecDot(self%v,av%v,dot,ierr)
DEBUG_STACK_POP
end function vec_dot_vec
!---------------------------------------------------------------------------------
!> Dot product with a complex vector
!---------------------------------------------------------------------------------
function vec_dot_cvec(self,a) result(dot)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
class(oft_cvector), target, intent(inout) :: a !< Second vector for dot product
complex(c8) :: dot !< \f$ \sum_i self_i a_i \f$
CALL oft_abort('Complex LA not yet supported with PETSc.','vec_dot_cvec',__FILE__)
end function vec_dot_cvec
!---------------------------------------------------------------------------------
!> Dot product with an array of vectors
!---------------------------------------------------------------------------------
function vec_mdot_vec(self,a,n) result(dots)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
type(oft_vector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
real(r8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
class(oft_petsc_vector), pointer :: av
integer(i4) :: i,ierr
TYPE(tvec), ALLOCATABLE, DIMENSION(:) :: vecs
DEBUG_STACK_PUSH
ALLOCATE(vecs(n))
DO i=1,n
  IF(.NOT.oft_petsc_vector_cast(av,a(i)%f))CALL oft_abort('"a" is not a vector.','vec_mdot_vec',__FILE__)
  IF(self%n/=av%n)call oft_abort('Vector lengths do not match.','vec_mdot_vec',__FILE__)
  vecs(i)=av%v
END DO
CALL VecMDot(self%v,n,vecs,dots,ierr)
DEALLOCATE(vecs)
DEBUG_STACK_POP
end function vec_mdot_vec
!---------------------------------------------------------------------------------
!> Dot product with an array of complex vectors
!---------------------------------------------------------------------------------
function vec_mdot_cvec(self,a,n) result(dots)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
type(oft_cvector_ptr), intent(inout) :: a(n) !< Array of vectors for dot product [n]
integer(i4), intent(in) :: n !< Length of vector array
complex(c8) :: dots(n) !< \f$ \sum_i self_i a(j)_i \f$
class(oft_petsc_vector), pointer :: av
CALL oft_abort('Complex LA not yet supported with PETSc.','vec_mdot_cvec',__FILE__)
end function vec_mdot_cvec
!---------------------------------------------------------------------------------
!> Sum reduction over vector
!---------------------------------------------------------------------------------
function vec_sum(self) result(sum)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
real(r8) :: sum !< Sum of vector elements
integer(i4) :: ierr
DEBUG_STACK_PUSH
CALL VecSum(self%v,sum,ierr)
DEBUG_STACK_POP
end function vec_sum
!---------------------------------------------------------------------------------
!> Compute norm of vector
!---------------------------------------------------------------------------------
function vec_norm(self,itype) result(norm)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: itype !< Type of norm (1-> 1-norm, 2-> 2-norm, 3-> Inf-norm)
real(r8) :: norm !< Specified norm of vector
integer(i4) :: ierr
DEBUG_STACK_PUSH
SELECT CASE(itype)
  CASE(1)
    CALL VecNorm(self%v,NORM_1,norm,ierr)
  CASE(2)
    CALL VecNorm(self%v,NORM_2,norm,ierr)
  CASE(3)
    CALL VecNorm(self%v,NORM_INFINITY,norm,ierr)
  CASE DEFAULT
    CALL oft_abort("Invalid norm type","vec_norm",__FILE__)
END SELECT
DEBUG_STACK_POP
end function vec_norm
!---------------------------------------------------------------------------------
!> Perform global stitching
!---------------------------------------------------------------------------------
subroutine vec_stitch(self,up_method)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
integer(i4), intent(in) :: up_method !< Type of stitching to perform
CALL oft_abort('Vector stitching is not supported for PETSc vectors.','vec_stitch',__FILE__)
end subroutine vec_stitch
!---------------------------------------------------------------------------------
!> Finalize vector
!---------------------------------------------------------------------------------
subroutine vec_delete(self)
class(oft_petsc_vector), intent(inout) :: self !< Vector object
integer(i4) :: ierr
DEBUG_STACK_PUSH
CALL VecDestroy(self%v,ierr)
CALL VecDestroy(self%vloc,ierr)
self%n=0
self%ng=0
self%loc_current=.FALSE.
self%loc_store=.FALSE.
DEBUG_STACK_POP
end subroutine vec_delete
!------------------------------------------------------------------------------
!> Cast a oft_matrix object to a oft_petsc_matrix
!!
!! The source matrix must be @ref oft_petsc_matrix or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_petsc_matrix_cast(self,source) result(success)
class(oft_petsc_matrix), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_matrix), target, intent(in) :: source !< Source solver to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_petsc_matrix)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
END FUNCTION oft_petsc_matrix_cast
!---------------------------------------------------------------------------------
!> Compute matrix-vector product
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine mat_apply_vec(self,a,b)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
class(oft_petsc_vector), pointer :: av,bv
integer(i4) :: ierr
DEBUG_STACK_PUSH
if(b%n/=self%nr)call oft_abort('Row mismatch','mat_apply_vec',__FILE__)
if(a%n/=self%nc)call oft_abort('Col mismatch','mat_apply_vec',__FILE__)
IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a PETSc vector object.','mat_apply_vec',__FILE__)
IF(.NOT.oft_petsc_vector_cast(bv,b))CALL oft_abort('"b" is not a PETSc vector object.','mat_apply_vec',__FILE__)
!---Apply operator
CALL MatMult(self%M,av%v,bv%v,ierr)
bv%loc_current=.FALSE.
DEBUG_STACK_POP
end subroutine mat_apply_vec
!---------------------------------------------------------------------------------
!> Compute matrix-vector product (complex)
!!
!! b = self * a
!---------------------------------------------------------------------------------
subroutine mat_apply_cvec(self,a,b)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source vector
class(oft_cvector), intent(inout) :: b !< Result of matrix product
CALL oft_abort('Complex LA not yet supported with PETSc.','mat_apply_cvec',__FILE__)
end subroutine mat_apply_cvec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine mat_applyt_vec(self,a,b)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), target, intent(inout) :: a !< Source vector
class(oft_vector), intent(inout) :: b !< Result of matrix product
class(oft_petsc_vector), pointer :: av,bv
integer(i4) :: ierr
DEBUG_STACK_PUSH
if(b%n/=self%nc)call oft_abort('Row mismatch','mat_applyt_vec',__FILE__)
if(a%n/=self%nr)call oft_abort('Col mismatch','mat_applyt_vec',__FILE__)
IF(.NOT.oft_petsc_vector_cast(av,a))CALL oft_abort('"a" is not a PETSc vector object.','mat_applyt_vec',__FILE__)
IF(.NOT.oft_petsc_vector_cast(bv,b))CALL oft_abort('"b" is not a PETSc vector object.','mat_applyt_vec',__FILE__)
!---Apply operator
CALL MatMultTranspose(self%M,av%v,bv%v,ierr)
bv%loc_current=.FALSE.
DEBUG_STACK_POP
end subroutine mat_applyt_vec
!---------------------------------------------------------------------------------
!> Apply matrix vector product for matrix transpose (complex vector)
!!
!! b = self^T * a
!---------------------------------------------------------------------------------
subroutine mat_applyt_cvec(self,a,b)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
class(oft_cvector), target, intent(inout) :: a !< Source vector
class(oft_cvector), intent(inout) :: b !< Result of matrix product
CALL oft_abort('Complex LA not yet supported with PETSc.','mat_applyt_cvec',__FILE__)
end subroutine mat_applyt_cvec
!---------------------------------------------------------------------------------
!> Set values of a matrix
!---------------------------------------------------------------------------------
subroutine mat_set_values(self,i_inds,j_inds,b,n,m,iblock,jblock)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to set [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to set [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4) :: j,ierr
DEBUG_STACK_PUSH
IF(XOR(PRESENT(iblock),PRESENT(jblock)))call oft_abort('Only one block index was supplied', &
  'mat_set_values',__FILE__)
!---
!$omp critical (petsc_mat)
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  DO j=1,m
    CALL MatSetValuesLocal(self%M,n,self%i_map(iblock)%offset+i_inds-1, &
      1,self%j_map(jblock)%offset+j_inds(j:j)-1,b(:,j),INSERT_VALUES,ierr)
  END DO
ELSE
  DO j=1,m
    CALL MatSetValuesLocal(self%M,n,i_inds-1,1,j_inds(j:j)-1,b(:,j),ADD_VALUES,ierr)
  END DO
END IF
!$omp end critical (petsc_mat)
DEBUG_STACK_POP
end subroutine mat_set_values
!---------------------------------------------------------------------------------
!> Add values to a matrix
!---------------------------------------------------------------------------------
subroutine mat_add_values(self,i_inds,j_inds,b,n,m,iblock,jblock,loc_cache)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: i_inds(n) !< Row indices of entries to add [n]
integer(i4), intent(in) :: j_inds(m) !< Column indices of entries to add [m]
real(r8), intent(in) :: b(n,m) !< Values to set [n,m]
integer(i4), intent(in) :: n !< Number of rows in local matrix
integer(i4), intent(in) :: m !< Number of columns in local matrix
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
integer(i4), optional, intent(in) :: jblock !< Column block (optional)
integer(i4), optional, intent(inout) :: loc_cache(n,m) !< Cache of entry locations
integer(i4) :: j,ierr
DEBUG_STACK_PUSH
IF(XOR(PRESENT(iblock),PRESENT(jblock)))call oft_abort('Only one block index was supplied', &
  'mat_add_values',__FILE__)
!---
!$omp critical (petsc_mat)
IF(PRESENT(iblock).AND.PRESENT(jblock))THEN
  DO j=1,m
    CALL MatSetValuesLocal(self%M,n,self%i_map(iblock)%offset+i_inds-1, &
      1,self%j_map(jblock)%offset+j_inds(j:j)-1,b(:,j),ADD_VALUES,ierr)
  END DO
ELSE
  DO j=1,m
    CALL MatSetValuesLocal(self%M,n,i_inds-1,1,j_inds(j:j)-1,b(:,j),ADD_VALUES,ierr)
  END DO
END IF
!$omp end critical (petsc_mat)
DEBUG_STACK_POP
end subroutine mat_add_values
!---------------------------------------------------------------------------------
!> Finish assembly of matrix and optionally extract diagonals
!---------------------------------------------------------------------------------
subroutine mat_assemble(self,diag)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
class(oft_vector), optional, target, intent(inout) :: diag !< Diagonal entries of matrix [nr] (optional)
integer(i4) :: i,ierr
INTEGER(petsc_addr) :: Md
real(r8) :: dmin,dmax
DEBUG_STACK_PUSH
!---
CALL MatAssemblyBegin(self%M,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(self%M,MAT_FINAL_ASSEMBLY,ierr)
!---Setup diagonal scaling
if(PRESENT(diag))then
  IF(associated(self%D))CALL self%D%delete
  CALL diag%new(self%D)
  SELECT TYPE(this=>self%D)
    CLASS IS(oft_petsc_vector)
      CALL MatGetDiagonal(self%M,this%v,ierr)
  END SELECT
  call diag%add(0.d0,1.d0,self%D)
end if
!---Common assembly tasks
DEBUG_STACK_POP
end subroutine mat_assemble
!---------------------------------------------------------------------------------
!> Zero all entries in matrix
!---------------------------------------------------------------------------------
subroutine mat_zero(self)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
integer(i4) :: ierr
DEBUG_STACK_PUSH
!---
CALL MatZeroEntries(self%m,ierr)
!---Common assembly tasks
DEBUG_STACK_POP
end subroutine mat_zero
!---------------------------------------------------------------------------------
!> Zero all entries in the specified rows
!---------------------------------------------------------------------------------
subroutine mat_zero_rows(self,nrows,irows,iblock,keep_diag)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
integer(i4), intent(in) :: nrows !< Number of rows to zero
integer(i4), intent(in) :: irows(nrows) !< Indices of rows to zero [nrows]
integer(i4), optional, intent(in) :: iblock !< Row block (optional)
logical, optional, intent(in) :: keep_diag !< Keep diagonal entries
integer(i4) :: ierr
real(r8), parameter :: one = 1.d0
real(r8), parameter :: zero = 0.d0
logical :: zero_diag
DEBUG_STACK_PUSH
zero_diag=.TRUE.
IF(PRESENT(keep_diag))zero_diag=keep_diag
IF(zero_diag)THEN
  CALL MatZeroRowsLocal(self%m,nrows,self%i_map(iblock)%offset+irows-1, &
    one,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
ELSE
  CALL MatZeroRowsLocal(self%m,nrows,self%i_map(iblock)%offset+irows-1, &
    zero,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
END IF
DEBUG_STACK_POP
end subroutine mat_zero_rows
!---------------------------------------------------------------------------------
!> Delete matrix
!---------------------------------------------------------------------------------
subroutine mat_delete(self)
class(oft_petsc_matrix), intent(inout) :: self !< Matrix object
integer(i4) :: ierr
DEBUG_STACK_PUSH
CALL MatDestroy(self%M,ierr)
! CALL VecDestroy(self%Md,ierr)
IF(ASSOCIATED(self%D))THEN
  CALL self%D%delete
  NULLIFY(self%D)
END IF
CALL ISLocalToGlobalMappingDestroy(self%r_lis,ierr)
CALL ISLocalToGlobalMappingDestroy(self%c_lis,ierr)
NULLIFY(self%i_map,self%j_map)
DEBUG_STACK_POP
end subroutine mat_delete
#endif
end module oft_petsc_la
