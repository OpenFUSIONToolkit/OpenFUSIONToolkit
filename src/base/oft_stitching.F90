!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_stitching.F90
!
!> Subroutines for general global stitching operations
!!
!! @author Chris Hansen
!! @date Spring 2010
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_stitching
USE, INTRINSIC :: iso_c_binding, ONLY: c_int8_t
USE oft_base
USE oft_scatter
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Global stitch structure
!! - Global linkage information
!! - Global ownership and boundary flags
!! - Preallocated seam arrays
!---------------------------------------------------------------------------------
TYPE :: oft_seam
  LOGICAL :: full = .FALSE. !< Flag indicating local mesh is complete
  LOGICAL :: skip = .FALSE. !< Flag for skipping stitching operations
  INTEGER(i4) :: nbemax = 0 !< Maximum number of boundary elements
  INTEGER(i4) :: nown = -1 !< Length of `own` (-1 if not yet built)
  INTEGER(i4) :: nbe = 0 !< Number of local boundary elements
  INTEGER(i4) :: nie = 0 !< Number of interior elements
  INTEGER(i4) :: nle = 0 !< Number of local boundary elements
  INTEGER(i4) :: nproc_con = 0 !< Number of connections to other processors
  INTEGER(i4) :: proc_split = 0 !< Location of self in processor list
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: leo => NULL() !< List of boundary element ownership
  INTEGER(c_int8_t), POINTER, CONTIGUOUS, DIMENSION(:) :: own => NULL() !< Ownership weight (1/0) for each local element (built on first reduction, which must occur outside of OpenMP parallel regions)
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: be => NULL() !< Boundary element flag
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lbe => NULL() !< List of boundary elements in full vector
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: lie => NULL() !< List of interior elements in full vector
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kle => NULL() !< Pointer to list of boundary element linkage
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: proc_con => NULL() !< Processor connectivity list
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:,:) :: lle => NULL() !< List of boundary element linkage
#ifdef OFT_MPI_F08
  TYPE(mpi_request), POINTER, DIMENSION(:) :: send_reqs => NULL() !< Asynchronous MPI Send requests
  TYPE(mpi_request), POINTER, DIMENSION(:) :: recv_reqs => NULL() !< Asynchronous MPI Recv requests
#else
  INTEGER(i4), POINTER, DIMENSION(:) :: send_reqs => NULL() !< Asynchronous MPI Send requests
  INTEGER(i4), POINTER, DIMENSION(:) :: recv_reqs => NULL() !< Asynchronous MPI Recv requests
#endif
  !> Scatter plan for stitching (built on first use by \ref oft_stitching::seam_setup_plan
  !! "seam_setup_plan"; owns the flat send/receive buffers and the MPI requests)
  !!
  !! A given seam must not be stitched concurrently from multiple threads. The plan
  !! and the arrays below are allocatable components, so a copied seam never shares
  !! them with the original.
  TYPE(oft_scatter_plan) :: plan
  INTEGER(i4) :: plan_nbe = -1 !< Value of `nbe` when the plan was built
  INTEGER(i4) :: plan_nle = -1 !< Value of `nle` when the plan was built
  LOGICAL :: plan_full = .FALSE. !< Value of `full` when the plan was built
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: smask !< Ownership of each packed plan entry
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: kcon !< Pointer to contributions of each boundary element [nbe+1]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: kloc !< First local-link contribution of each boundary element [nbe]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: icon !< Position of each contribution in the plan receive buffer
  INTEGER(i4) :: nsingle = 0 !< Number of boundary elements with exactly one contribution
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lsingle !< Location in full vector of elements with one contribution [nsingle]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: isingle !< Position of their contribution in the plan receive buffer [nsingle]
  !> Number of boundary elements with three or more values (own plus two or more
  !! contributions), summed in ascending order (-1 until the plan is built)
  INTEGER(i4) :: nmulti = -1
  INTEGER(i4) :: multi_max = 0 !< Maximum number of values of a boundary element
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lmulti !< Boundary elements with three or more values [nmulti]
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: rscratch !< Persistent scratch for real sorted sums [multi_max]
  COMPLEX(c8), ALLOCATABLE, DIMENSION(:) :: cscratch !< Persistent scratch for complex sorted sums [multi_max]
END TYPE oft_seam
!---------------------------------------------------------------------------------
!> Seam object pointer
!---------------------------------------------------------------------------------
TYPE :: seam_list
  TYPE(oft_seam), POINTER :: s => NULL() !< Seam information
END TYPE seam_list
!------------------------------------------------------------------------------
!> Perform a global dot product for vectors with a given seam structure
!------------------------------------------------------------------------------
INTERFACE oft_global_dp
  MODULE PROCEDURE global_dp_r8
  MODULE PROCEDURE global_dp_c8
  MODULE PROCEDURE global_dp_rc8
  MODULE PROCEDURE global_dp_cr8
END INTERFACE oft_global_dp
!---------------------------------------------------------------------------------
!> Perform a global 1-norm (sum of absolute values) for a vector with a given seam structure
!---------------------------------------------------------------------------------
INTERFACE oft_global_norm1
  MODULE PROCEDURE global_norm1_r8
  MODULE PROCEDURE global_norm1_c8
END INTERFACE oft_global_norm1
!---------------------------------------------------------------------------------
!> Perform a global sum reduction for a vector with a given seam structure
!---------------------------------------------------------------------------------
INTERFACE oft_global_reduction
  MODULE PROCEDURE global_reduction_r8
  MODULE PROCEDURE global_reduction_c8
END INTERFACE oft_global_reduction
!---------------------------------------------------------------------------------
!> Stitch values along domain boundaries
!!
!! General subroutine for domain stitching
!! Available stitching methods
!! - Synchronize global vector using element ownership [0] (initialization)
!! - Sum proc contributions [1] (stitching)
!!
!! Summation produces the same value, bit for bit, on every copy of a shared element.
!! With two contributions this holds because floating point addition is commutative.
!! Elements with three or more contributions (eg. periodic images on one process, or
!! elements shared by three or more processes) are summed in ascending order of value,
!! which every copy computes identically.
!---------------------------------------------------------------------------------
INTERFACE oft_global_stitch
  MODULE PROCEDURE global_stitch_r8
  MODULE PROCEDURE global_stitch_c8
END INTERFACE oft_global_stitch
CONTAINS
!---------------------------------------------------------------------------------
!> Build local ownership weights used by global reductions
!!
!! Reductions sum only owned entries directly, instead of summing all entries and
!! subtracting non-owned (shared/periodic) duplicates. The latter is subject to
!! catastrophic cancellation when duplicates dominate the sum, which can produce
!! negative values for self dot-products. A 0/1 integer weight is used instead of a
!! logical mask as the resulting multiply vectorizes with GCC, while a masked
!! selection (MERGE/IF) does not.
!!
!! @note Weights are built lazily on the first reduction with a given seam. This is
!! not thread safe, so it aborts if called from within an OpenMP parallel region.
!---------------------------------------------------------------------------------
SUBROUTINE seam_build_own(self,n)
TYPE(oft_seam), INTENT(inout) :: self !< Seam structure
INTEGER(i4), INTENT(in) :: n !< Length of local arrays
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Lazy initialization is not thread safe
IF(omp_in_parallel())CALL oft_abort('Ownership weights must be built outside of a parallel region', &
  'seam_build_own',__FILE__)
IF(ASSOCIATED(self%own))DEALLOCATE(self%own)
ALLOCATE(self%own(n))
self%own=1_c_int8_t
DO i=1,self%nbe
  IF(.NOT.self%leo(i))self%own(self%lbe(i))=0_c_int8_t
END DO
self%nown=n
DEBUG_STACK_POP
END SUBROUTINE seam_build_own
!---------------------------------------------------------------------------------
!> Local owned-only real(r8) dot product (explicit-shape arrays allow vectorization)
!---------------------------------------------------------------------------------
FUNCTION owned_dp_r8(a,b,w,n) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Length of local arrays
REAL(r8), INTENT(in) :: a(n) !< Local vector 1 for dot_product
REAL(r8), INTENT(in) :: b(n) !< Local vector 2 for dot_product
INTEGER(c_int8_t), INTENT(in) :: w(n) !< Ownership weights (1/0)
INTEGER(i4) :: i
REAL(r8) :: c
c=0.d0
!$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
DO i=1,n
  c=c+a(i)*b(i)*REAL(w(i),r8)
END DO
END FUNCTION owned_dp_r8
!---------------------------------------------------------------------------------
!> Local owned-only complex(c8) dot product (explicit-shape arrays allow vectorization)
!---------------------------------------------------------------------------------
FUNCTION owned_dp_c8(a,b,w,n) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Length of local arrays
COMPLEX(c8), INTENT(in) :: a(n) !< Local vector 1 for dot_product
COMPLEX(c8), INTENT(in) :: b(n) !< Local vector 2 for dot_product
INTEGER(c_int8_t), INTENT(in) :: w(n) !< Ownership weights (1/0)
INTEGER(i4) :: i
COMPLEX(c8) :: c
c=(0.d0,0.d0)
!$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
DO i=1,n
  c=c+CONJG(a(i))*b(i)*REAL(w(i),r8)
END DO
END FUNCTION owned_dp_c8
!---------------------------------------------------------------------------------
!> Local owned-only real(r8) sum (explicit-shape arrays allow vectorization)
!---------------------------------------------------------------------------------
FUNCTION owned_sum_r8(a,w,n) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Length of local array
REAL(r8), INTENT(in) :: a(n) !< Local data for reduction
INTEGER(c_int8_t), INTENT(in) :: w(n) !< Ownership weights (1/0)
INTEGER(i4) :: i
REAL(r8) :: c
c=0.d0
!$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
DO i=1,n
  c=c+a(i)*REAL(w(i),r8)
END DO
END FUNCTION owned_sum_r8
!---------------------------------------------------------------------------------
!> Local owned-only complex(c8) sum (explicit-shape arrays allow vectorization)
!---------------------------------------------------------------------------------
FUNCTION owned_sum_c8(a,w,n) RESULT(c)
INTEGER(i4), INTENT(in) :: n !< Length of local array
COMPLEX(c8), INTENT(in) :: a(n) !< Local data for reduction
INTEGER(c_int8_t), INTENT(in) :: w(n) !< Ownership weights (1/0)
INTEGER(i4) :: i
COMPLEX(c8) :: c
c=(0.d0,0.d0)
!$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
DO i=1,n
  c=c+a(i)*REAL(w(i),r8)
END DO
END FUNCTION owned_sum_c8
!---------------------------------------------------------------------------------
!> real(r8) dot_product implementation of \ref oft_stitching::oft_global_reduction
!---------------------------------------------------------------------------------
function global_dp_r8(self,a,b,n,no_reduce) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local arrays
real(r8), intent(in) :: a(n) !< Local vector 1 for dot_product
real(r8), intent(in) :: b(n) !< Local vector 2 for dot_product
logical, optional, intent(in) :: no_reduce !< Skip global reduction step
logical :: do_reduce
integer(i4) :: i
real(r8) :: c
DEBUG_STACK_PUSH
do_reduce=(.NOT.self%full)
IF(PRESENT(no_reduce))do_reduce=do_reduce.AND.(.NOT.no_reduce)
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  IF(self%nbe>0)THEN
    IF(self%nown/=n)CALL seam_build_own(self,n)
    c=owned_dp_r8(a,b,self%own,n)
  ELSE
    !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
    do i=1,n
      c=c+a(i)*b(i)
    end do
  END IF
END IF
IF(do_reduce)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_dp_r8
!---------------------------------------------------------------------------------
!> complex(c8) dot_product implementation of \ref oft_stitching::oft_global_reduction
!---------------------------------------------------------------------------------
function global_dp_c8(self,a,b,n,no_reduce) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local arrays
COMPLEX(c8), intent(in) :: a(n) !< Local vector 1 for dot_product
COMPLEX(c8), intent(in) :: b(n) !< Local vector 2 for dot_product
logical, optional, intent(in) :: no_reduce !< Skip global reduction step
logical :: do_reduce
integer(i4) :: i
COMPLEX(c8) :: c
DEBUG_STACK_PUSH
do_reduce=(.NOT.self%full)
IF(PRESENT(no_reduce))do_reduce=do_reduce.AND.(.NOT.no_reduce)
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  IF(self%nbe>0)THEN
    IF(self%nown/=n)CALL seam_build_own(self,n)
    c=owned_dp_c8(a,b,self%own,n)
  ELSE
    !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
    do i=1,n
      c=c+CONJG(a(i))*b(i)
    end do
  END IF
END IF
IF(do_reduce)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_dp_c8
!---------------------------------------------------------------------------------
!> Mixed real(r8)/complex(c8) dot_product implementation of \ref oft_stitching::oft_global_dp
!!
!! Computes \f$ \sum_i a_i b_i \f$ (no conjugation is needed for real `a`)
!---------------------------------------------------------------------------------
function global_dp_rc8(self,a,b,n,no_reduce) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local arrays
real(r8), intent(in) :: a(n) !< Local vector 1 for dot_product
COMPLEX(c8), intent(in) :: b(n) !< Local vector 2 for dot_product
logical, optional, intent(in) :: no_reduce !< Skip global reduction step
logical :: do_reduce
integer(i4) :: i
COMPLEX(c8) :: c
DEBUG_STACK_PUSH
do_reduce=(.NOT.self%full)
IF(PRESENT(no_reduce))do_reduce=do_reduce.AND.(.NOT.no_reduce)
!---Global reduction
c=(0.d0,0.d0)
IF(.NOT.self%skip)THEN
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+a(i)*b(i)
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-a(self%lbe(i))*b(self%lbe(i))
  end do
  !$omp end parallel
END IF
IF(do_reduce)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_dp_rc8
!---------------------------------------------------------------------------------
!> Mixed complex(c8)/real(r8) dot_product implementation of \ref oft_stitching::oft_global_dp
!!
!! Computes \f$ \sum_i \overline{a_i} b_i \f$
!---------------------------------------------------------------------------------
function global_dp_cr8(self,a,b,n,no_reduce) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local arrays
COMPLEX(c8), intent(in) :: a(n) !< Local vector 1 for dot_product
real(r8), intent(in) :: b(n) !< Local vector 2 for dot_product
logical, optional, intent(in) :: no_reduce !< Skip global reduction step
logical :: do_reduce
integer(i4) :: i
COMPLEX(c8) :: c
DEBUG_STACK_PUSH
do_reduce=(.NOT.self%full)
IF(PRESENT(no_reduce))do_reduce=do_reduce.AND.(.NOT.no_reduce)
!---Global reduction
c=(0.d0,0.d0)
IF(.NOT.self%skip)THEN
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+CONJG(a(i))*b(i)
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-CONJG(a(self%lbe(i)))*b(self%lbe(i))
  end do
  !$omp end parallel
END IF
IF(do_reduce)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_dp_cr8
!---------------------------------------------------------------------------------
!> real(r8) array implementation of \ref oft_stitching::oft_global_reduction
!---------------------------------------------------------------------------------
function global_reduction_r8(self,a,n) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local array for reduction
real(r8), intent(in) :: a(n) !< Local data for reduction
integer(i4) :: i
real(r8) :: c
DEBUG_STACK_PUSH
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  IF(self%nbe>0)THEN
    IF(self%nown/=n)CALL seam_build_own(self,n)
    c=owned_sum_r8(a,self%own,n)
  ELSE
    !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
    do i=1,n
      c=c+a(i)
    end do
  END IF
END IF
IF(.NOT.self%full)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_reduction_r8
!---------------------------------------------------------------------------------
!> complex(c8) array implementation of \ref oft_stitching::oft_global_reduction
!---------------------------------------------------------------------------------
function global_reduction_c8(self,a,n) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local array for reduction
COMPLEX(c8), intent(in) :: a(n) !< Local data for reduction
integer(i4) :: i
COMPLEX(c8) :: c
DEBUG_STACK_PUSH
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  IF(self%nbe>0)THEN
    IF(self%nown/=n)CALL seam_build_own(self,n)
    c=owned_sum_c8(a,self%own,n)
  ELSE
    !$omp parallel do simd reduction(+:c) if(n>OFT_OMP_VTHRESH)
    do i=1,n
      c=c+a(i)
    end do
  END IF
END IF
IF(.NOT.self%full)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_reduction_c8
!---------------------------------------------------------------------------------
!> real(r8) implementation of \ref oft_stitching::oft_global_norm1
!---------------------------------------------------------------------------------
function global_norm1_r8(self,a,n) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local array for reduction
real(r8), intent(in) :: a(n) !< Local data for reduction
integer(i4) :: i
real(r8) :: c !< \f$ \sum_i |a_i| \f$ over owned entries on all processes
DEBUG_STACK_PUSH
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+ABS(a(i))
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-ABS(a(self%lbe(i)))
  end do
  !$omp end parallel
END IF
IF(.NOT.self%full)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_norm1_r8
!---------------------------------------------------------------------------------
!> complex(c8) implementation of \ref oft_stitching::oft_global_norm1
!---------------------------------------------------------------------------------
function global_norm1_c8(self,a,n) result(c)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: n !< Length of local array for reduction
COMPLEX(c8), intent(in) :: a(n) !< Local data for reduction
integer(i4) :: i
real(r8) :: c !< \f$ \sum_i |a_i| \f$ over owned entries on all processes
DEBUG_STACK_PUSH
!---Global reduction
c=0.d0
IF(.NOT.self%skip)THEN
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+ABS(a(i))
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-ABS(a(self%lbe(i)))
  end do
  !$omp end parallel
END IF
IF(.NOT.self%full)c=oft_mpi_sum(c)
DEBUG_STACK_POP
end function global_norm1_c8
!---------------------------------------------------------------------------------
!> Build the scatter plan and combination lists used by global stitching
!!
!! The plan (tag @ref oft_scatter_base::oft_tag_stitch "OFT_TAG_STITCH") sends the
!! values of all linked boundary elements to each neighbour (`lsend`=`lrecv`, in
!! linkage order) and holds the local (periodic) links as its local part, so after
!! the exchange every contribution to a boundary element is in the plan's receive
!! buffer. The contributions to each boundary element are listed by element (`kcon`,
!! `icon`): remote contributions in neighbour order, then local links (from `kloc`).
!! Elements with exactly one contribution are also listed in `lsingle`/`isingle`
!! and elements with two or more (three or more values including their own) in
!! `lmulti`. Rebuilt automatically when `nbe`, `nle` or `full` change.
!---------------------------------------------------------------------------------
subroutine seam_setup_plan(self)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4) :: i,j,k,e,p,active_procs,nsend,nloc
integer(i4), allocatable :: ksend(:),lsend(:),lsrc(:),ldst(:),bpos(:),kfill(:),kfill_loc(:),nrem(:)
DEBUG_STACK_PUSH
CALL seam_free_plan(self)
active_procs=self%nproc_con
IF(self%full)active_procs=0
nsend=self%kle(active_procs+1)-self%kle(1)
nloc=self%kle(1)-self%kle(0)
!---Remote entries in linkage order, local links as local part
ALLOCATE(ksend(active_procs+1),lsend(MAX(1,nsend)),lsrc(MAX(1,nloc)),ldst(MAX(1,nloc)))
ALLOCATE(bpos(self%nle),self%smask(nsend+nloc))
DO j=1,active_procs+1
  ksend(j)=self%kle(j)-self%kle(1)+1
END DO
DO i=self%kle(1),self%kle(active_procs+1)-1
  p=i-self%kle(1)+1
  lsend(p)=self%lbe(self%lle(1,i))
  bpos(i)=p
END DO
DO i=self%kle(0),self%kle(1)-1
  k=i-self%kle(0)+1
  lsrc(k)=self%lbe(self%lle(1,i))
  ldst(k)=self%lbe(self%lle(2,i))
  bpos(i)=nsend+k
END DO
!---Ownership of packed entries (synchronization by ownership sends zero for unowned)
IF(ASSOCIATED(self%leo))THEN
  DO i=self%kle(0),self%kle(active_procs+1)-1
    self%smask(bpos(i))=self%leo(self%lle(1,i))
  END DO
ELSE
  self%smask=.TRUE.
END IF
CALL scatter_plan_setup(self%plan,active_procs,self%proc_con(1:active_procs),ksend,lsend,ksend,lsend, &
  nloc,lsrc,ldst,tag=OFT_TAG_STITCH)
!---Contributions to each boundary element (remote in neighbour order, then local links)
ALLOCATE(self%kcon(self%nbe+1),self%kloc(self%nbe),nrem(self%nbe))
self%kcon=0
nrem=0
DO i=self%kle(1),self%kle(active_procs+1)-1
  e=self%lle(1,i)
  nrem(e)=nrem(e)+1
  self%kcon(e+1)=self%kcon(e+1)+1
END DO
DO i=self%kle(0),self%kle(1)-1
  e=self%lle(2,i) ! Local linkage is applied in reverse
  self%kcon(e+1)=self%kcon(e+1)+1
END DO
self%kcon(1)=1
DO e=1,self%nbe
  self%kcon(e+1)=self%kcon(e+1)+self%kcon(e)
  self%kloc(e)=self%kcon(e)+nrem(e)
END DO
ALLOCATE(self%icon(MAX(1,self%kcon(self%nbe+1)-1)),kfill(self%nbe),kfill_loc(self%nbe))
kfill=self%kcon(1:self%nbe)
kfill_loc=self%kloc
DO i=self%kle(1),self%kle(active_procs+1)-1
  e=self%lle(1,i)
  self%icon(kfill(e))=bpos(i)
  kfill(e)=kfill(e)+1
END DO
DO i=self%kle(0),self%kle(1)-1
  e=self%lle(2,i)
  self%icon(kfill_loc(e))=bpos(i)
  kfill_loc(e)=kfill_loc(e)+1
END DO
!---Elements with one contribution and with two or more
self%nsingle=0
self%nmulti=0
self%multi_max=0
DO e=1,self%nbe
  k=self%kcon(e+1)-self%kcon(e)
  IF(k==1)self%nsingle=self%nsingle+1
  IF(k>=2)self%nmulti=self%nmulti+1
  self%multi_max=MAX(self%multi_max,k+1)
END DO
ALLOCATE(self%lsingle(MAX(1,self%nsingle)),self%isingle(MAX(1,self%nsingle)),self%lmulti(MAX(1,self%nmulti)))
self%nsingle=0
self%nmulti=0
DO e=1,self%nbe
  k=self%kcon(e+1)-self%kcon(e)
  IF(k==1)THEN
    self%nsingle=self%nsingle+1
    self%lsingle(self%nsingle)=self%lbe(e)
    self%isingle(self%nsingle)=self%icon(self%kcon(e))
  ELSE IF(k>=2)THEN
    self%nmulti=self%nmulti+1
    self%lmulti(self%nmulti)=e
  END IF
END DO
DEALLOCATE(ksend,lsend,lsrc,ldst,bpos,kfill,kfill_loc,nrem)
self%plan_nbe=self%nbe
self%plan_nle=self%nle
self%plan_full=self%full
DEBUG_STACK_POP
end subroutine seam_setup_plan
!---------------------------------------------------------------------------------
!> Release the scatter plan and combination lists of a seam (rebuilt on next stitch)
!---------------------------------------------------------------------------------
subroutine seam_free_plan(self)
type(oft_seam), intent(inout) :: self !< Seam structure
CALL scatter_plan_delete(self%plan)
IF(ALLOCATED(self%smask))DEALLOCATE(self%smask)
IF(ALLOCATED(self%kcon))DEALLOCATE(self%kcon)
IF(ALLOCATED(self%kloc))DEALLOCATE(self%kloc)
IF(ALLOCATED(self%icon))DEALLOCATE(self%icon)
IF(ALLOCATED(self%lsingle))DEALLOCATE(self%lsingle)
IF(ALLOCATED(self%isingle))DEALLOCATE(self%isingle)
IF(ALLOCATED(self%lmulti))DEALLOCATE(self%lmulti)
IF(ALLOCATED(self%rscratch))DEALLOCATE(self%rscratch)
IF(ALLOCATED(self%cscratch))DEALLOCATE(self%cscratch)
self%nsingle=0
self%nmulti=-1
self%multi_max=0
self%plan_nbe=-1
end subroutine seam_free_plan
!---------------------------------------------------------------------------------
!> Check that the stitching plan of a seam is current, building it if needed
!---------------------------------------------------------------------------------
subroutine seam_check_plan(self)
type(oft_seam), intent(inout) :: self !< Seam structure
IF((.NOT.self%plan%setup).OR.(self%plan_nbe/=self%nbe).OR.(self%plan_nle/=self%nle) &
  .OR.(self%plan_full.NEQV.self%full))CALL seam_setup_plan(self)
end subroutine seam_check_plan
!---------------------------------------------------------------------------------
!> real(r8) implementation of \ref oft_stitching::oft_global_stitch
!!
!! Values are exchanged with the seam's scatter plan (see @ref seam_setup_plan) and
!! combined in place by boundary element:
!! - method 0: unowned elements are set to the sum of the contributions, in which
!!   only the owner's value is non-zero;
!! - method 1: each element is incremented by its contribution if it has one, and
!!   elements with three or more values (own plus two or more contributions) are set
!!   to the sum of all values in ascending order, which every copy computes
!!   identically;
!! - method 2: each element is set to the sum of its remote contributions (local
!!   links and its own value are not included), in neighbour order.
!!
!! The plan, its buffers and the sorted-sum scratch belong to the seam, so a given
!! seam must not be stitched concurrently; MPI calls are made from the calling
!! thread only.
!---------------------------------------------------------------------------------
subroutine global_stitch_r8(self,a,up_method)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: up_method !< Stitching method (0, 1 or 2)
real(r8), intent(inout) :: a(:) !< Local data to be stitched
integer(i4) :: i,k,l,m,n,e
real(r8) :: tmp
IF(self%skip.OR.self%nbe==0)RETURN
DEBUG_STACK_PUSH
if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','global_stitch_r8',__FILE__)
CALL seam_check_plan(self)
IF(up_method==0)THEN
  IF(.NOT.ASSOCIATED(self%leo))CALL oft_abort('Synchronization requires ownership','global_stitch_r8',__FILE__)
  CALL scatter_begin(self%plan,a,mask=self%smask)
ELSE
  CALL scatter_begin(self%plan,a,nolocal=(up_method==2))
END IF
CALL scatter_wait_recv(self%plan,a)
ASSOCIATE(buf=>self%plan%rrbuf)
SELECT CASE(up_method)
CASE(0)
  !$omp parallel do private(k,tmp) if(self%nbe>OFT_OMP_VTHRESH)
  DO e=1,self%nbe
    IF(self%leo(e))CYCLE
    tmp=0.d0
    DO k=self%kcon(e),self%kcon(e+1)-1
      tmp=tmp+buf(self%icon(k))
    END DO
    a(self%lbe(e))=tmp
  END DO
CASE(1)
  !$omp parallel do simd if(self%nsingle>OFT_OMP_VTHRESH)
  DO i=1,self%nsingle
    a(self%lsingle(i))=a(self%lsingle(i))+buf(self%isingle(i))
  END DO
  !---Elements with three or more values: sum in ascending order
  IF(self%nmulti>0)THEN
    IF(.NOT.ALLOCATED(self%rscratch))ALLOCATE(self%rscratch(self%multi_max))
    DO m=1,self%nmulti
      e=self%lmulti(m)
      n=1
      self%rscratch(1)=a(self%lbe(e))
      DO k=self%kcon(e),self%kcon(e+1)-1
        n=n+1
        self%rscratch(n)=buf(self%icon(k))
      END DO
      !---Insertion sort (lists are short)
      DO k=2,n
        tmp=self%rscratch(k)
        DO l=k-1,1,-1
          IF(self%rscratch(l)<=tmp)EXIT
          self%rscratch(l+1)=self%rscratch(l)
        END DO
        self%rscratch(l+1)=tmp
      END DO
      tmp=0.d0
      DO k=1,n
        tmp=tmp+self%rscratch(k)
      END DO
      a(self%lbe(e))=tmp
    END DO
  END IF
CASE(2)
  !$omp parallel do private(k,tmp) if(self%nbe>OFT_OMP_VTHRESH)
  DO e=1,self%nbe
    tmp=0.d0
    DO k=self%kcon(e),self%kloc(e)-1
      tmp=tmp+buf(self%icon(k))
    END DO
    a(self%lbe(e))=tmp
  END DO
END SELECT
END ASSOCIATE
CALL scatter_wait_send(self%plan,a)
DEBUG_STACK_POP
end subroutine global_stitch_r8
!---------------------------------------------------------------------------------
!> complex(c8) implementation of \ref oft_stitching::oft_global_stitch
!!
!! See @ref global_stitch_r8. Values with three or more contributions are summed in
!! ascending order of real part, then imaginary part.
!---------------------------------------------------------------------------------
subroutine global_stitch_c8(self,a,up_method)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: up_method !< Stitching method (0, 1 or 2)
complex(c8), intent(inout) :: a(:) !< Local data to be stitched
integer(i4) :: i,k,l,m,n,e
complex(c8) :: tmp
IF(self%skip.OR.self%nbe==0)RETURN
DEBUG_STACK_PUSH
if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','global_stitch_c8',__FILE__)
CALL seam_check_plan(self)
IF(up_method==0)THEN
  IF(.NOT.ASSOCIATED(self%leo))CALL oft_abort('Synchronization requires ownership','global_stitch_c8',__FILE__)
  CALL scatter_begin(self%plan,a,mask=self%smask)
ELSE
  CALL scatter_begin(self%plan,a,nolocal=(up_method==2))
END IF
CALL scatter_wait_recv(self%plan,a)
ASSOCIATE(buf=>self%plan%crbuf)
SELECT CASE(up_method)
CASE(0)
  !$omp parallel do private(k,tmp) if(self%nbe>OFT_OMP_VTHRESH)
  DO e=1,self%nbe
    IF(self%leo(e))CYCLE
    tmp=(0.d0,0.d0)
    DO k=self%kcon(e),self%kcon(e+1)-1
      tmp=tmp+buf(self%icon(k))
    END DO
    a(self%lbe(e))=tmp
  END DO
CASE(1)
  !$omp parallel do simd if(self%nsingle>OFT_OMP_VTHRESH)
  DO i=1,self%nsingle
    a(self%lsingle(i))=a(self%lsingle(i))+buf(self%isingle(i))
  END DO
  !---Elements with three or more values: sum in ascending order
  IF(self%nmulti>0)THEN
    IF(.NOT.ALLOCATED(self%cscratch))ALLOCATE(self%cscratch(self%multi_max))
    DO m=1,self%nmulti
      e=self%lmulti(m)
      n=1
      self%cscratch(1)=a(self%lbe(e))
      DO k=self%kcon(e),self%kcon(e+1)-1
        n=n+1
        self%cscratch(n)=buf(self%icon(k))
      END DO
      !---Insertion sort (lists are short)
      DO k=2,n
        tmp=self%cscratch(k)
        DO l=k-1,1,-1
          IF(REAL(self%cscratch(l),r8)<REAL(tmp,r8))EXIT
          IF((REAL(self%cscratch(l),r8)==REAL(tmp,r8)).AND.(AIMAG(self%cscratch(l))<=AIMAG(tmp)))EXIT
          self%cscratch(l+1)=self%cscratch(l)
        END DO
        self%cscratch(l+1)=tmp
      END DO
      tmp=(0.d0,0.d0)
      DO k=1,n
        tmp=tmp+self%cscratch(k)
      END DO
      a(self%lbe(e))=tmp
    END DO
  END IF
CASE(2)
  !$omp parallel do private(k,tmp) if(self%nbe>OFT_OMP_VTHRESH)
  DO e=1,self%nbe
    tmp=(0.d0,0.d0)
    DO k=self%kcon(e),self%kloc(e)-1
      tmp=tmp+buf(self%icon(k))
    END DO
    a(self%lbe(e))=tmp
  END DO
END SELECT
END ASSOCIATE
CALL scatter_wait_send(self%plan,a)
DEBUG_STACK_POP
end subroutine global_stitch_c8
!---------------------------------------------------------------------------------
!> Validate stitching structure
!---------------------------------------------------------------------------------
subroutine oft_stitch_check(self)
type(oft_seam), intent(inout) :: self
integer(i4) :: j,ierr
integer(i4), allocatable :: le_size(:,:),le_out(:,:)
DEBUG_STACK_PUSH
CALL oft_mpi_barrier(ierr)
IF(self%full)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---Check linkage pointer
IF(ANY(self%kle>self%nle+1))THEN
  WRITE(*,*)'ERROR: ',oft_env%rank,MAXVAL(self%kle),MAXLOC(self%kle),self%nle
  CALL oft_abort('BAD LINKAGE PTR','oft_stitch_check',__FILE__)
END IF
!---Check matching transfer sizes
#ifdef HAVE_MPI
ALLOCATE(le_size(2,self%nproc_con),le_out(2,self%nproc_con))
DO j=1,self%nproc_con
  le_out(:,j)=self%kle(j+1)-self%kle(j)
  CALL MPI_ISEND(le_out(:,j),2,OFT_MPI_I4,self%proc_con(j),2,oft_env%comm,self%send_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','oft_stitch_check',__FILE__)
  CALL MPI_IRECV(le_size(:,j),2,OFT_MPI_I4,self%proc_con(j),2,oft_env%comm,self%recv_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','oft_stitch_check',__FILE__)
END DO
CALL oft_mpi_waitall(self%nproc_con,self%recv_reqs,ierr)
CALL oft_mpi_waitall(self%nproc_con,self%send_reqs,ierr)
DO j=1,self%nproc_con
  IF(ANY(le_size(:,j)/=le_out(:,j)))THEN
    WRITE(*,*)'STITCH ERROR: ',le_out(:,j),le_size(:,j)
    CALL oft_abort('BAD ARRAY SIZE MATCH','oft_stitch_check',__FILE__)
  END IF
END DO
DEALLOCATE(le_size,le_out)
#endif
DEBUG_STACK_POP
end subroutine oft_stitch_check
!---------------------------------------------------------------------------------
!> Release the index data and stitching plan of a seam
!---------------------------------------------------------------------------------
subroutine destory_seam(self)
type(oft_seam), intent(inout) :: self
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%leo))DEALLOCATE(self%leo)
IF(ASSOCIATED(self%own))DEALLOCATE(self%own)
self%nown=-1
IF(ASSOCIATED(self%be))DEALLOCATE(self%be)
IF(ASSOCIATED(self%lbe))DEALLOCATE(self%lbe)
IF(ASSOCIATED(self%lie))DEALLOCATE(self%lie)
IF(ASSOCIATED(self%kle))DEALLOCATE(self%kle)
IF(ASSOCIATED(self%lle))DEALLOCATE(self%lle)
CALL seam_free_plan(self)
DEBUG_STACK_POP
end subroutine destory_seam
end module oft_stitching
