!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file fem_utils.F90
!
!> FEM utility classes and functions
!! - FEM interpolator classes
!! - Field averaging
!!
!! @authors Chris Hansen
!! @date March 2013
!! @ingroup doxy_oft_fem
!---------------------------------------------------------------------------
module fem_utils
USE oft_base
USE oft_quadrature
USE oft_mesh_type, ONLY: mesh, cell_is_curved
USE oft_mesh_local, ONLY: mesh_local_partition
USE oft_stitching, ONLY: oft_global_stitch
USE oft_la_base, ONLY: oft_matrix
USE fem_base, ONLY: oft_afem_type, oft_fem_type, oft_bfem_type
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! CLASS fem_interp
!---------------------------------------------------------------------------
!> Base class for interpolation of a FE field
!---------------------------------------------------------------------------
type, abstract :: fem_interp
  integer(i4) :: dim = 0 !< Dimension of field
  class(fem_interp), pointer :: parent => NULL() !< Parent interpolator
contains
  !> Reconstruct field
  procedure(oft_fem_interp), deferred :: interp
  !> Setup reconstruction object
  procedure :: setup => fem_interp_setup
  !> Delete reconstruction object
  procedure :: delete => fem_interp_delete
end type fem_interp
!---
abstract interface
!---------------------------------------------------------------------------
!> Protoype for FE interpolation method
!---------------------------------------------------------------------------
  subroutine oft_fem_interp(self,cell,f,gop,val)
    import fem_interp, i4, r8
    class(fem_interp), intent(inout) :: self
    integer(i4), intent(in) :: cell
    real(r8), intent(in) :: f(:)
    real(r8), intent(in) :: gop(3,4)
    real(r8), intent(out) :: val(:)
  end subroutine oft_fem_interp
end interface
!---------------------------------------------------------------------------
! CLASS bfem_interp
!---------------------------------------------------------------------------
!> Base class for interpolation of a FE field
!---------------------------------------------------------------------------
type, abstract :: bfem_interp
  integer(i4) :: dim = 0 !< Dimension of field
  class(bfem_interp), pointer :: parent => NULL() !< Parent interpolator
contains
  !> Reconstruct field
  procedure(oft_bfem_interp), deferred :: interp
  !> Setup reconstruction object
  procedure :: setup => bfem_interp_setup
  !> Delete reconstruction object
  procedure :: delete => bfem_interp_delete
end type bfem_interp
!---
abstract interface
!---------------------------------------------------------------------------
!> Protoype for FE interpolation method
!---------------------------------------------------------------------------
  subroutine oft_bfem_interp(self,cell,f,gop,val)
    import bfem_interp, i4, r8
    class(bfem_interp), intent(inout) :: self
    integer(i4), intent(in) :: cell
    real(r8), intent(in) :: f(:)
    real(r8), intent(in) :: gop(3,3)
    real(r8), intent(out) :: val(:)
  end subroutine oft_bfem_interp
end interface
!---------------------------------------------------------------------------
! CLASS cc_interp
!---------------------------------------------------------------------------
!> Interpolator for cell centered fields
!---------------------------------------------------------------------------
type, extends(fem_interp) :: cc_interp
  real(r8), contiguous, pointer :: bcc(:,:) => NULL() !< Field values in each cell
contains
  !> Reconstruct field
  procedure :: interp => cc_interp_apply
end type cc_interp
!---------------------------------------------------------------------------
! CLASS diff_interp
!---------------------------------------------------------------------------
!> Interpolator for difference between two fields
!---------------------------------------------------------------------------
type, extends(fem_interp) :: diff_interp
  class(fem_interp), pointer :: a => NULL() !< Field 1
  class(fem_interp), pointer :: b => NULL() !< Field 2
contains
  !> Reconstruct field
  procedure :: interp => diff_interp_apply
end type diff_interp
!---------------------------------------------------------------------------
! CLASS dot_interp
!---------------------------------------------------------------------------
!> Interpolator for dot-product of two fields
!---------------------------------------------------------------------------
type, extends(fem_interp) :: dot_interp
  class(fem_interp), pointer :: a => NULL() !< Field 1
  class(fem_interp), pointer :: b => NULL() !< Field 2
contains
  !> Reconstruct field
  procedure :: interp => dot_interp_apply
end type dot_interp
!---------------------------------------------------------------------------
! CLASS cross_interp
!---------------------------------------------------------------------------
!> Interpolator for cross-product of two fields
!---------------------------------------------------------------------------
type, extends(fem_interp) :: cross_interp
  class(fem_interp), pointer :: a => NULL() !< Field 1
  class(fem_interp), pointer :: b => NULL() !< Field 2
contains
  !> Reconstruct field
  procedure :: interp => cross_interp_apply
end type cross_interp
!---------------------------------------------------------------------------
! CLASS tensor_dot_interp
!---------------------------------------------------------------------------
!> Interpolator for the product of a vector and tensor field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: tensor_dot_interp
  INTEGER(i4) :: bshape = 1 !< Number of rows in tensor B
  REAL(r8), POINTER :: bvals(:) => NULL()  !< Internal storage for intepolated B
  class(fem_interp), pointer :: a => NULL() !< Vector field
  class(fem_interp), pointer :: b => NULL() !< Tensor field
contains
  !> Setup reconstruction object
  procedure :: setup => tensor_dot_interp_setup
  !> Reconstruct field
  procedure :: interp => tensor_dot_interp_apply
  !> Delete reconstruction object
  procedure :: delete => tensor_dot_interp_delete
end type tensor_dot_interp
contains
!---------------------------------------------------------------------------
! SUBROUTINE: cc_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct a cell centered field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine cc_interp_apply(self,cell,f,gop,val)
class(cc_interp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
DEBUG_STACK_PUSH
!---
val=self%bcc(:,cell)
DEBUG_STACK_POP
end subroutine cc_interp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: diff_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct the difference between two fields
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine diff_interp_apply(self,cell,f,gop,val)
class(diff_interp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8), allocatable, dimension(:) :: aval,bval
DEBUG_STACK_PUSH
IF(self%dim<=0)CALL oft_abort("Field dimension must be specified.", &
  "diff_interp_apply",__FILE__)
ALLOCATE(aval(self%dim),bval(self%dim))
CALL self%a%interp(cell,f,gop,aval)
CALL self%b%interp(cell,f,gop,bval)
val=aval-bval
DEALLOCATE(aval,bval)
DEBUG_STACK_POP
end subroutine diff_interp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: dot_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct the dot-product of two fields
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine dot_interp_apply(self,cell,f,gop,val)
class(dot_interp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: aval(3),bval(3)
DEBUG_STACK_PUSH
CALL self%a%interp(cell,f,gop,aval)
CALL self%b%interp(cell,f,gop,bval)
val(1)=DOT_PRODUCT(aval,bval)
DEBUG_STACK_POP
end subroutine dot_interp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: cross_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct the cross-product of two fields
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine cross_interp_apply(self,cell,f,gop,val)
class(cross_interp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: aval(3),bval(3)
DEBUG_STACK_PUSH
CALL self%a%interp(cell,f,gop,aval)
CALL self%b%interp(cell,f,gop,bval)
val=cross_product(aval,bval)
DEBUG_STACK_POP
end subroutine cross_interp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: tensor_dot_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct the product of a vector and tensor field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine tensor_dot_interp_apply(self,cell,f,gop,val)
class(tensor_dot_interp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: aval(3)
DEBUG_STACK_PUSH
CALL self%a%interp(cell,f,gop,aval)
CALL self%b%interp(cell,f,gop,self%bvals)
val=MATMUL(aval,RESHAPE(self%bvals,(/3,self%bshape/)))
DEBUG_STACK_POP
end subroutine tensor_dot_interp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: tensor_dot_interp_setup
!---------------------------------------------------------------------------
!> Setup composite interpolator for a matrix-vector product
!!
!! Allocates local interpolation objects. Setup of component fields must be
!! called separately before the interpolator may be used.
!!
!! @note Should only be used via class @ref tensor_dot_interp or children
!---------------------------------------------------------------------------
subroutine tensor_dot_interp_setup(self)
class(tensor_dot_interp), intent(inout) :: self
IF(.NOT.ASSOCIATED(self%bvals))ALLOCATE(self%bvals(3*self%bshape))
end subroutine tensor_dot_interp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: tensor_dot_interp_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class @ref tensor_dot_interp or children
!---------------------------------------------------------------------------
subroutine tensor_dot_interp_delete(self)
class(tensor_dot_interp), intent(inout) :: self
IF(ASSOCIATED(self%bvals))DEALLOCATE(self%bvals)
end subroutine tensor_dot_interp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: fem_interp_setup
!---------------------------------------------------------------------------
!> Dummy setup function
!---------------------------------------------------------------------------
subroutine fem_interp_setup(self)
class(fem_interp), intent(inout) :: self
call oft_warn('Setup called on general interpolator, this may indicate an error.')
end subroutine fem_interp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: fem_interp_delete
!---------------------------------------------------------------------------
!> Dummy destroy function
!---------------------------------------------------------------------------
subroutine fem_interp_delete(self)
class(fem_interp), intent(inout) :: self
call oft_warn('Finalize called on general interpolator, this may indicate an error.')
end subroutine fem_interp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: bfem_interp_setup
!---------------------------------------------------------------------------
!> Dummy setup function
!---------------------------------------------------------------------------
subroutine bfem_interp_setup(self)
class(bfem_interp), intent(inout) :: self
call oft_warn('Setup called on general interpolator, this may indicate an error.')
end subroutine bfem_interp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: bfem_interp_delete
!---------------------------------------------------------------------------
!> Dummy destroy function
!---------------------------------------------------------------------------
subroutine bfem_interp_delete(self)
class(bfem_interp), intent(inout) :: self
call oft_warn('Finalize called on general interpolator, this may indicate an error.')
end subroutine bfem_interp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: fem_avg_bcc
!---------------------------------------------------------------------------
!> Average a FE interpolator field to cell centers, by volume averaging.
!!
!! @param[in,out] field Source field intepolator
!! @param[in,out] bcc Averaged field over each cell
!! @param[in] order Desired integration order
!! @param[in] n Dimension of field (optional)
!---------------------------------------------------------------------------
subroutine fem_avg_bcc(field,bcc,order,n)
CLASS(fem_interp), INTENT(inout) :: field
REAL(r8), INTENT(inout) :: bcc(:,:)
INTEGER(i4), INTENT(in) :: order
INTEGER(i4), OPTIONAL, INTENT(in) :: n
INTEGER(i4) :: i,m,nf
REAL(r8) :: v,f(4),goptmp(3,4)
REAL(r8), ALLOCATABLE :: dets(:),bcctmp(:)
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
nf=1
IF(PRESENT(n))nf=n
CALL mesh%quad_rule(order,quad)
!---Construct operators
!$omp parallel private(dets,curved,goptmp,m,v,bcctmp)
ALLOCATE(dets(quad%np),bcctmp(nf))
!$omp do
DO i=1,mesh%nc
  bcc(:,i)=0.d0
  ! Get reconstructed operators
  curved=cell_is_curved(mesh,i)
  IF(.NOT.curved)CALL mesh%jacobian(i,quad%pts(:,1),goptmp,v)
  DO m=1,quad%np
    IF(curved)CALL mesh%jacobian(i,quad%pts(:,m),goptmp,v)
    dets(m)=v*quad%wts(m)
    CALL field%interp(i,quad%pts(:,m),goptmp,bcctmp)
    bcc(:,i)=bcc(:,i)+bcctmp*dets(m)
  END DO
  bcc(:,i)=bcc(:,i)/SUM(dets)
END DO
DEALLOCATE(dets,bcctmp)
!$omp end parallel
CALL quad%delete
DEBUG_STACK_POP
end subroutine fem_avg_bcc
!---------------------------------------------------------------------------
! SUBROUTINE: fem_partition
!---------------------------------------------------------------------------
!> Partition FE weights based on geometric connectivity
!!
!! @param[in,out] self Finite element structure
!! @param[in,out] part Weight partitioning [self%ne]
!! @param[in] nparts Number of partitions
!---------------------------------------------------------------------------
subroutine fem_partition(self,part,nparts)
class(oft_fem_type), intent(inout) :: self
INTEGER(i4), intent(inout) :: part(:)
INTEGER(i4), INTENT(in) :: nparts
INTEGER(i4) :: i,j,k,nloc,offset,mycounts(4)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: tloc
TYPE(oft_1d_int), ALLOCATABLE, DIMENSION(:) :: tloc_p,tloc_e,tloc_f,tloc_c
DEBUG_STACK_PUSH
IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Partitioning FE weights',nparts
!---Partition mesh
ALLOCATE(tloc_p(nparts),tloc_e(nparts))
ALLOCATE(tloc_f(nparts),tloc_c(nparts))
CALL mesh_local_partition(self%mesh,tloc_p,tloc_e,tloc_f,tloc_c,nparts)
!---Ownership
!$omp parallel do private(i,j,k,nloc,mycounts,offset,tloc)
DO k=1,nparts
  mycounts=(/tloc_p(k)%n,tloc_e(k)%n,tloc_f(k)%n,tloc_c(k)%n/)
  nloc=SUM(mycounts*self%gstruct)
  ALLOCATE(tloc(nloc))
  nloc=0
  offset=0
  do j=1,tloc_p(k)%n
    do i=1,self%gstruct(1)
      nloc=nloc+1
      tloc(nloc)=(tloc_p(k)%v(j)-1)*self%gstruct(1)+i
    end do
  end do
  offset=self%gstruct(1)*mesh%np
  do j=1,tloc_e(k)%n
    do i=1,self%gstruct(2)
      nloc=nloc+1
      tloc(nloc)=(tloc_e(k)%v(j)-1)*self%gstruct(2)+i+offset
    end do
  end do
  offset=offset+self%gstruct(2)*self%mesh%ne
  do j=1,tloc_f(k)%n
    do i=1,self%gstruct(3)
      nloc=nloc+1
      tloc(nloc)=(tloc_f(k)%v(j)-1)*self%gstruct(3)+i+offset
    end do
  end do
  offset=offset+self%gstruct(3)*self%mesh%nf
  do j=1,tloc_c(k)%n
    do i=1,self%gstruct(4)
      nloc=nloc+1
      tloc(nloc)=(tloc_c(k)%v(j)-1)*self%gstruct(4)+i+offset
    end do
  end do
  do i=1,nloc
    part(tloc(i))=k
  end do
  DEALLOCATE(tloc)
  DEALLOCATE(tloc_p(k)%v,tloc_e(k)%v,tloc_f(k)%v,tloc_c(k)%v)
END DO
DEALLOCATE(tloc_p,tloc_e,tloc_f,tloc_c)
DEBUG_STACK_POP
end subroutine fem_partition
!---------------------------------------------------------------------------
! SUBROUTINE: fem_map_flag
!---------------------------------------------------------------------------
!> Set diagonal elements to one on owned rows according to BC flag
!---------------------------------------------------------------------------
subroutine fem_map_flag(fem_obj,vert_flag,edge_flag,face_flag,fe_flag)
CLASS(oft_fem_type), INTENT(INOUT) :: fem_obj
LOGICAL, DIMENSION(:), INTENT(IN) :: vert_flag,edge_flag,face_flag
LOGICAL, DIMENSION(:), INTENT(INOUT) :: fe_flag
INTEGER(4) :: i,j,offset
REAL(8), ALLOCATABLE, DIMENSION(:) :: flag_tmp
ALLOCATE(flag_tmp(fem_obj%ne))
flag_tmp=0.d0
!---Point DOFs
do i=1,fem_obj%mesh%np
  IF(vert_flag(i))THEN
    offset=(i-1)*fem_obj%gstruct(1)
    do j=1,fem_obj%gstruct(1)
      flag_tmp(j+offset) = 1.d0
    end do
  END IF
end do
!---Edge DOFs
do i=1,fem_obj%mesh%ne
  IF(edge_flag(i))THEN
    offset=(i-1)*fem_obj%gstruct(2) + fem_obj%mesh%np*fem_obj%gstruct(1)
    do j=1,fem_obj%gstruct(2)
      flag_tmp(j+offset) = 1.d0
    end do
  END IF
end do
!---Face DOFs
do i=1,fem_obj%mesh%nf
  IF(face_flag(i))THEN
    offset=(i-1)*fem_obj%gstruct(3) + fem_obj%mesh%ne*fem_obj%gstruct(2) &
      + fem_obj%mesh%np*fem_obj%gstruct(1)
    do j=1,fem_obj%gstruct(3)
      flag_tmp(j+offset) = 1.d0
    end do
  END IF
end do
CALL oft_global_stitch(fem_obj%linkage,flag_tmp,1)
fe_flag=(flag_tmp>0.5d0)
DEALLOCATE(flag_tmp)
end subroutine fem_map_flag
!---------------------------------------------------------------------------
! SUBROUTINE: fem_map_flag
!---------------------------------------------------------------------------
!> Set diagonal elements to one on owned rows according to BC flag
!---------------------------------------------------------------------------
subroutine bfem_map_flag(fem_obj,vert_flag,edge_flag,fe_flag)
CLASS(oft_bfem_type), INTENT(INOUT) :: fem_obj
LOGICAL, DIMENSION(:), INTENT(IN) :: vert_flag,edge_flag
LOGICAL, DIMENSION(:), INTENT(INOUT) :: fe_flag
INTEGER(4) :: i,j,offset
REAL(8), ALLOCATABLE, DIMENSION(:) :: flag_tmp
ALLOCATE(flag_tmp(fem_obj%ne))
flag_tmp=0.d0
!---Point DOFs
do i=1,fem_obj%mesh%np
  IF(vert_flag(i))THEN
    offset=(i-1)*fem_obj%gstruct(1)
    do j=1,fem_obj%gstruct(1)
      flag_tmp(j+offset) = 1.d0
    end do
  END IF
end do
!---Edge DOFs
do i=1,fem_obj%mesh%ne
  IF(edge_flag(i))THEN
    offset=(i-1)*fem_obj%gstruct(2) + fem_obj%mesh%np*fem_obj%gstruct(1)
    do j=1,fem_obj%gstruct(2)
      flag_tmp(j+offset) = 1.d0
    end do
  END IF
end do
CALL oft_global_stitch(fem_obj%linkage,flag_tmp,1)
fe_flag=(flag_tmp>0.5d0)
DEALLOCATE(flag_tmp)
end subroutine bfem_map_flag
!---------------------------------------------------------------------------
! SUBROUTINE: fem_dirichlet_diag
!---------------------------------------------------------------------------
!> Set diagonal elements to one on owned rows according to BC flag
!---------------------------------------------------------------------------
subroutine fem_dirichlet_diag(fem_obj,mat,flag,iblock)
CLASS(oft_afem_type), INTENT(INOUT) :: fem_obj
CLASS(oft_matrix), INTENT(INOUT) :: mat
LOGICAL, DIMENSION(:), INTENT(IN) :: flag
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock
INTEGER(i4) :: i,j,jtmp(1)
REAL(r8) :: dtmp(1,1)
!---Set diagonal entries for dirichlet terms
dtmp=1.d0
DO i=1,fem_obj%ne
  IF(fem_obj%be(i).OR.(.NOT.flag(i)))CYCLE
  jtmp=i
  CALL mat%add_values(jtmp,jtmp,dtmp,1,1,iblock,iblock)
END DO
DO i=1,fem_obj%nbe
  IF(.NOT.fem_obj%linkage%leo(i))CYCLE
  j=fem_obj%lbe(i)
  IF(flag(j))THEN
    jtmp=j
    CALL mat%add_values(jtmp,jtmp,dtmp,1,1,iblock,iblock)
  END IF
END DO
end subroutine fem_dirichlet_diag
!---------------------------------------------------------------------------
! SUBROUTINE: fem_dirichlet_vec
!---------------------------------------------------------------------------
!> Replace values in local vector according to BC flag
!---------------------------------------------------------------------------
subroutine fem_dirichlet_vec(fem_obj,vecin,vecout,flag)
CLASS(oft_afem_type), INTENT(INOUT) :: fem_obj
REAL(r8), DIMENSION(:), INTENT(INOUT) :: vecin,vecout
LOGICAL, DIMENSION(:), INTENT(IN) :: flag
INTEGER(i4) :: i,j
!$omp parallel do if(fem_obj%ne>OFT_OMP_VTHRESH)
DO i=1,fem_obj%ne
  IF(flag(i))THEN
    vecout(i)=0.d0
    IF(.NOT.fem_obj%be(i))vecout(i)=vecin(i)
  END IF
END DO
DO i=1,fem_obj%nbe
  j=fem_obj%lbe(i)
  IF(.NOT.fem_obj%linkage%leo(i))CYCLE
  IF(flag(j))vecout(j)=vecin(j)
END DO
end subroutine fem_dirichlet_vec
end module fem_utils
