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
USE, INTRINSIC :: iso_c_binding, ONLY: c_int8_t, c_loc, c_int, c_long, c_ptr
USE oft_base
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
  LOGICAL :: replicated = .FALSE. !< All entries are stored on every processor, each owned by exactly one processor
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
  TYPE(oft_1d_real), POINTER, DIMENSION(:) :: send => NULL() !< Preallocated MPI send arrays
  TYPE(oft_1d_real), POINTER, DIMENSION(:) :: recv => NULL() !< Preallocated MPI recv arrays
  TYPE(oft_1d_comp), POINTER, DIMENSION(:) :: csend => NULL() !< Preallocated MPI send arrays
  TYPE(oft_1d_comp), POINTER, DIMENSION(:) :: crecv => NULL() !< Preallocated MPI recv arrays
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: rep_counts => NULL() !< Number of entries owned by each processor (replicated only) [nprocs]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: rep_displs => NULL() !< Offset of each processor's entries in owner order (replicated only) [nprocs]
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: rep_perm => NULL() !< Local index of each entry in owner order (replicated only) [nbe]
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
END INTERFACE oft_global_dp
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
!! For replicated seams (see \ref oft_stitching::seam_setup_replicated) these operations
!! use collective communication over all processors instead of domain-to-domain linkage.
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
!> real(r8) implementation of \ref oft_stitching::oft_global_stitch
!---------------------------------------------------------------------------------
subroutine global_stitch_r8(self,a,up_method)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: up_method !< Stitching method (0 or 1)
real(r8), intent(inout) :: a(:) !< Local data to be stitched
real(r8), allocatable, dimension(:) :: atmp
integer(i4) :: i,j,js,jn,elocal,active_procs,ierr
IF(self%skip)RETURN
IF(self%replicated)THEN
  CALL replicated_stitch_r8(self,a,up_method)
  RETURN
END IF
IF(self%nbe==0)RETURN
DEBUG_STACK_PUSH
if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','global_stitch_r8',__FILE__)
active_procs=self%nproc_con
IF(self%full)active_procs=0
!---Create Recv calls
#ifdef HAVE_MPI
DO j=1,active_procs
  !---Skip processors with no elements
  IF(self%send(j)%n==0)THEN
    self%recv_reqs(j)=MPI_REQUEST_NULL
    CYCLE
  END IF
  CALL MPI_IRECV(self%recv(j)%v,self%recv(j)%n,OFT_MPI_R8,self%proc_con(j),1,oft_env%comm,self%recv_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','global_stitch_r8',__FILE__)
END DO
#endif
!---Extract boundary elements
ALLOCATE(atmp(self%nbe))
!$omp simd
DO i=1,self%nbe
  atmp(i)=a(self%lbe(i))
END DO
!---Zero unowned when synchronizing
IF(up_method==0)THEN
  !$omp simd
  DO i=1,self%nbe
    IF(.NOT.self%leo(i))atmp(i)=0.d0
  END DO
END IF
!---Create send arrays
!$omp parallel private(j,js,jn,i,elocal,ierr)
#ifdef HAVE_MPI
!$omp do
DO j=1,active_procs
  js=self%kle(j)
  jn=self%kle(j+1)-1
  !$omp simd private(elocal)
  do i=js,jn
    elocal=self%lle(1,i)
    self%send(j)%v(i-js+1)=atmp(elocal) ! Populate stitching array
  enddo
  !---Start MPI Send
  IF(self%send(j)%n==0)THEN ! Skip processors with no elements
    self%send_reqs(j)=MPI_REQUEST_NULL
  ELSE
    !$omp critical
    CALL MPI_ISEND(self%send(j)%v,self%send(j)%n,OFT_MPI_R8,self%proc_con(j),1,oft_env%comm,self%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','global_stitch_r8',__FILE__)
    !$omp end critical
  END IF
END DO
!$omp end do nowait
#endif
!$omp single
IF(up_method/=2)THEN
  js=self%kle(0)
  jn=self%kle(1)-1
  !$omp simd private(elocal)
  do i=js,jn
    elocal=self%lle(1,i)
    self%send(0)%v(i-js+1)=atmp(elocal) ! Populate stitching array
  enddo
END IF
!$omp end single
!$omp end parallel
atmp=0.d0
!---Add couplings (always same order)
#ifdef HAVE_MPI
IF((active_procs>0).AND.(self%proc_split>0))THEN
  !---Wait for all Recvs to complete
  CALL oft_mpi_waitall(self%proc_split,self%recv_reqs(1:self%proc_split),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_r8',__FILE__)
  DO j=1,self%proc_split
    js=self%kle(j)
    jn=self%kle(j+1)-1
    do i=js,jn
      elocal=self%lle(1,i)
      atmp(elocal)=atmp(elocal)+self%recv(j)%v(i-js+1) ! If linked sum values
    enddo
  END DO
END IF
#endif
IF(up_method==0)THEN
  !$omp simd
  DO i=1,self%nbe
    IF(self%leo(i))atmp(i)=atmp(i)+a(self%lbe(i))
  END DO
ELSE IF(up_method==1)THEN
  !$omp simd
  DO i=1,self%nbe
    atmp(i)=atmp(i)+a(self%lbe(i))
  END DO
END IF
!---Process local connectivity (Interaction is reverse)
IF(up_method/=2)THEN
  js=self%kle(0)
  jn=self%kle(1)-1
  do i=js,jn
    elocal=self%lle(2,i)
    atmp(elocal)=atmp(elocal)+self%send(0)%v(i-js+1) ! If linked sum values
  enddo
END IF
#ifdef HAVE_MPI
IF(active_procs>0)THEN
  CALL oft_mpi_waitall(self%nproc_con,self%recv_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_r8',__FILE__)
  DO j=self%proc_split+1,active_procs
    js=self%kle(j)
    jn=self%kle(j+1)-1
    do i=js,jn
      elocal=self%lle(1,i)
      atmp(elocal)=atmp(elocal)+self%recv(j)%v(i-js+1) ! If linked sum values
    enddo
  END DO
END IF
#endif
!---Move back to full array
!$omp simd
DO i=1,self%nbe
  a(self%lbe(i))=atmp(i)
END DO
DEALLOCATE(atmp)
!---Wait for all sends to complete
#ifdef HAVE_MPI
IF(active_procs>0)THEN
  CALL oft_mpi_waitall(self%nproc_con,self%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_r8',__FILE__)
END IF
#endif
DEBUG_STACK_POP
end subroutine global_stitch_r8
!---------------------------------------------------------------------------------
!> complex(c8) implementation of \ref oft_stitching::oft_global_stitch
!---------------------------------------------------------------------------------
subroutine global_stitch_c8(self,a,up_method)
type(oft_seam), intent(inout) :: self !< Seam structure
integer(i4), intent(in) :: up_method !< Stitching method (0 or 1)
COMPLEX(c8), intent(inout) :: a(:) !< Local data to be stitched
COMPLEX(c8), allocatable, dimension(:) :: atmp
integer(i4) :: i,j,js,jn,elocal,active_procs,ierr
IF(self%skip)RETURN
IF(self%replicated)THEN
  CALL replicated_stitch_c8(self,a,up_method)
  RETURN
END IF
IF(self%nbe==0)RETURN
DEBUG_STACK_PUSH
if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','global_stitch_c8',__FILE__)
active_procs=self%nproc_con
IF(self%full)active_procs=0
!---Create Recv calls
#ifdef HAVE_MPI
DO j=1,active_procs
  !---Skip processors with no elements
  IF(self%csend(j)%n==0)THEN
    self%recv_reqs(j)=MPI_REQUEST_NULL
    CYCLE
  END IF
  CALL MPI_IRECV(self%crecv(j)%v,self%crecv(j)%n,OFT_MPI_C8,self%proc_con(j),1,oft_env%comm,self%recv_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','global_stitch_c8',__FILE__)
END DO
#endif
!---Extract boundary elements
ALLOCATE(atmp(self%nbe))
!$omp simd
DO i=1,self%nbe
  atmp(i)=a(self%lbe(i))
END DO
!---Zero unowned when synchronizing
IF(up_method==0)THEN
  !$omp simd
  DO i=1,self%nbe
    IF(.NOT.self%leo(i))atmp(i)=0.d0
  END DO
END IF
!---Create send arrays
!$omp parallel private(j,js,jn,i,elocal,ierr)
#ifdef HAVE_MPI
!$omp do
DO j=1,active_procs
  js=self%kle(j)
  jn=self%kle(j+1)-1
  !$omp simd private(elocal)
  do i=js,jn
    elocal=self%lle(1,i)
    self%csend(j)%v(i-js+1)=atmp(elocal) ! Populate stitching array
  enddo
  !---Start MPI Send
  IF(self%csend(j)%n==0)THEN ! Skip processors with no elements
    self%send_reqs(j)=MPI_REQUEST_NULL
  ELSE
    !$omp critical
    CALL MPI_ISEND(self%csend(j)%v,self%csend(j)%n,OFT_MPI_C8,self%proc_con(j),1,oft_env%comm,self%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','global_stitch_c8',__FILE__)
    !$omp end critical
  END IF
END DO
!$omp end do nowait
#endif
!$omp single
IF(up_method/=2)THEN
  js=self%kle(0)
  jn=self%kle(1)-1
  !$omp simd private(elocal)
  do i=js,jn
    elocal=self%lle(1,i)
    self%csend(0)%v(i-js+1)=atmp(elocal) ! Populate stitching array
  enddo
END IF
!$omp end single
!$omp end parallel
atmp=0.d0
!---Add couplings (always same order)
#ifdef HAVE_MPI
IF((active_procs>0).AND.(self%proc_split>0))THEN
  !---Wait for all Recvs to complete
  CALL oft_mpi_waitall(self%proc_split,self%recv_reqs(1:self%proc_split),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_c8',__FILE__)
  DO j=1,self%proc_split
    js=self%kle(j)
    jn=self%kle(j+1)-1
    do i=js,jn
      elocal=self%lle(1,i)
      atmp(elocal)=atmp(elocal)+self%crecv(j)%v(i-js+1) ! If linked sum values
    enddo
  END DO
END IF
#endif
IF(up_method==0)THEN
  !$omp simd
  DO i=1,self%nbe
    IF(self%leo(i))atmp(i)=atmp(i)+a(self%lbe(i))
  END DO
ELSE IF(up_method==1)THEN
  !$omp simd
  DO i=1,self%nbe
    atmp(i)=atmp(i)+a(self%lbe(i))
  END DO
END IF
!---Process local connectivity (Interaction is reverse)
IF(up_method/=2)THEN
  js=self%kle(0)
  jn=self%kle(1)-1
  do i=js,jn
    elocal=self%lle(2,i)
    atmp(elocal)=atmp(elocal)+self%csend(0)%v(i-js+1) ! If linked sum values
  enddo
END IF
#ifdef HAVE_MPI
IF(active_procs>0)THEN
  CALL oft_mpi_waitall(self%nproc_con,self%recv_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_c8',__FILE__)
  DO j=self%proc_split+1,active_procs
    js=self%kle(j)
    jn=self%kle(j+1)-1
    do i=js,jn
      elocal=self%lle(1,i)
      atmp(elocal)=atmp(elocal)+self%crecv(j)%v(i-js+1) ! If linked sum values
    enddo
  END DO
END IF
#endif
!---Move back to full array
!$omp simd
DO i=1,self%nbe
  a(self%lbe(i))=atmp(i)
END DO
DEALLOCATE(atmp)
!---Wait for all sends to complete
#ifdef HAVE_MPI
IF(active_procs>0)THEN
  CALL oft_mpi_waitall(self%nproc_con,self%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','global_stitch_c8',__FILE__)
END IF
#endif
DEBUG_STACK_POP
end subroutine global_stitch_c8
! !---------------------------------------------------------------------------------
! !> Begin stitch operation and return while transfers complete
! !!
! !! @param[in,out] self Stitching information
! !! @param[in,out] a Local point data to be stitched
! !! @param[in] up_method Stitching method (0 or 1)
! !---------------------------------------------------------------------------------
! subroutine oft_global_stitch_begin(self,a,up_method)
! type(oft_seam), intent(inout) :: self
! integer(i4), intent(in) :: up_method
! real(r8), intent(inout) :: a(:)
! real(r8), pointer, dimension(:) :: leout
! integer(i4) :: i,j,js,jn,elocal,ierr
! DEBUG_STACK_PUSH
! IF(self%full.OR.self%skip)THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF

! if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','oft_global_stitch',__FILE__)
! !---
! if(up_method==1)then ! Sum values on each processor
!   DO j=1,oft_env%nproc_con
!     leout=>self%send(j)%v
!     js=self%kle(j)
!     jn=self%kle(j+1)-1
!     do i=js,jn
!       elocal=self%lle(1,i)
!       leout(i-js+1)=a(self%lbe(elocal)) ! Populate stitching array
!     enddo
!   END DO
! else if(up_method==0)then ! Sync values by ownership
!   DO j=1,oft_env%nproc_con
!     leout=>self%send(j)%v
!     js=self%kle(j)
!     jn=self%kle(j+1)-1
!     do i=js,jn
!       elocal=self%lle(1,i)
!       if(self%leo(elocal))then
!         leout(i-js+1)=a(self%lbe(elocal)) ! If point is owned populate stitching array
!       else
!         leout(i-js+1)=0.d0
!       end if
!     enddo
!   END DO
!   DO i=1,self%nbe
!     IF(self%leo(i))CYCLE
!     a(self%lbe(i))=0.d0
!   END DO
! else if(up_method==2)then ! Sync values by ownership
!   DO j=1,oft_env%nproc_con
!     leout=>self%send(j)%v
!     js=self%kle(j)
!     jn=self%kle(j+1)-1
!     do i=js,jn
!       elocal=self%lle(1,i)
!       leout(i-js+1)=a(self%lbe(elocal))
!       a(self%lbe(elocal))=0.d0
!     enddo
!   END DO
! endif
! !---Create Send and Recv calls
! #ifdef HAVE_MPI
! do j=1,oft_env%nproc_con
!   !---Skip processors with no elements
!   IF(self%send(j)%n==0)THEN
!     oft_env%send(j)=MPI_REQUEST_NULL
!     oft_env%recv(j)=MPI_REQUEST_NULL
!     CYCLE
!   END IF
!   call MPI_ISEND(self%send(j)%v,self%send(j)%n,OFT_MPI_R8,oft_env%proc_con(j),1,oft_env%comm,oft_env%send(j),ierr)
!   call MPI_IRECV(self%recv(j)%v,self%recv(j)%n,OFT_MPI_R8,oft_env%proc_con(j),1,oft_env%comm,oft_env%recv(j),ierr)
! enddo
! #endif
! DEBUG_STACK_POP
! end subroutine oft_global_stitch_begin
! !---------------------------------------------------------------------------------
! !> Complete stitch operation by processing transfers
! !!
! !! @param[in,out] self Stitching information
! !! @param[in,out] a Local point data to be stitched
! !! @param[in] up_method Stitching method (0 or 1)
! !---------------------------------------------------------------------------------
! subroutine oft_global_stitch_end(self,a,up_method)
! type(oft_seam), intent(inout) :: self
! integer(i4), intent(in) :: up_method
! real(r8), intent(inout) :: a(:)
! real(r8), pointer, dimension(:) :: letmp
! integer(i4) :: i,j,js,jn,elocal,eremote,ierr
! DEBUG_STACK_PUSH
! IF(self%full.OR.self%skip)THEN
!   DEBUG_STACK_POP
!   RETURN
! END IF
! if(up_method>2.OR.up_method<0)call oft_abort('Invalid stitch method','oft_global_stitch',__FILE__)
! !---Loop over each connected processor
! #ifdef HAVE_MPI
! do while(.TRUE.)
!   if(all(oft_env%recv==MPI_REQUEST_NULL))exit ! All recieves have been processed
!   CALL oft_mpi_waitany(oft_env%nproc_con,oft_env%recv,j,ierr) ! Wait for completed recieve
!   letmp=>self%recv(j)%v ! Point dummy input array to current Recv array
!   js=self%kle(j)
!   jn=self%kle(j+1)-1
!   if(up_method==1)then ! Sum values on each processor
!     do i=js,jn
!       elocal=self%lle(1,i)
!       eremote=i-js+1
!       a(self%lbe(elocal))=a(self%lbe(elocal))+letmp(eremote) ! If linked sum values
!     enddo
!   endif
!   if(up_method==0.OR.up_method==2)then ! Sync values by ownership
!     do i=js,jn
!       elocal=self%lle(1,i)
!       eremote=i-js+1
!       a(self%lbe(elocal))=a(self%lbe(elocal))+letmp(eremote) ! If contibution exists sync value
!     enddo
!   endif
! end do
! CALL oft_mpi_waitall(oft_env%nproc_con,oft_env%send,ierr) ! Wait for all sends to complete
! #endif
! DEBUG_STACK_POP
! end subroutine oft_global_stitch_end
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
IF(self%replicated)THEN
  CALL replicated_seam_check(self)
  DEBUG_STACK_POP
  RETURN
END IF
!---Check linkage pointer
IF(ANY(self%kle>self%nle+1))THEN
  WRITE(*,*)'ERROR: ',oft_env%rank,MAXVAL(self%kle),MAXLOC(self%kle),self%nle
  CALL oft_abort('BAD LINKAGE PTR','oft_stitch_check',__FILE__)
END IF
!---Check local sizes
DO j=1,self%nproc_con
  IF(self%send(j)%n==0.OR.self%recv(j)%n==0)CYCLE
  IF(SIZE(self%send(j)%v)/=self%send(j)%n)CALL oft_abort('BAD SEND ARRAY','oft_stitch_check',__FILE__)
  IF(SIZE(self%recv(j)%v)/=self%recv(j)%n)CALL oft_abort('BAD RECV ARRAY','oft_stitch_check',__FILE__)
END DO
!---Check matching transfer sizes
#ifdef HAVE_MPI
ALLOCATE(le_size(2,self%nproc_con),le_out(2,self%nproc_con))
DO j=1,self%nproc_con
  le_out(:,j)=[self%send(j)%n,self%recv(j)%n]
  CALL MPI_ISEND(le_out(:,j),2,OFT_MPI_I4,self%proc_con(j),2,oft_env%comm,self%send_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','oft_stitch_check',__FILE__)
  CALL MPI_IRECV(le_size(:,j),2,OFT_MPI_I4,self%proc_con(j),2,oft_env%comm,self%recv_reqs(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','oft_stitch_check',__FILE__)
END DO
CALL oft_mpi_waitall(self%nproc_con,self%recv_reqs,ierr)
CALL oft_mpi_waitall(self%nproc_con,self%send_reqs,ierr)
DO j=1,self%nproc_con
  IF((le_size(1,j)/=self%recv(j)%n).OR.(le_size(2,j)/=self%send(j)%n))THEN
    WRITE(*,*)'STITCH ERROR: ',[self%send(j)%n,self%recv(j)%n],le_size(:,j)
    CALL oft_abort('BAD ARRAY SIZE MATCH','oft_stitch_check',__FILE__)
  END IF
END DO
DEALLOCATE(le_size,le_out)
#endif
DEBUG_STACK_POP
end subroutine oft_stitch_check
!---------------------------------------------------------------------------------
!> Validate stitching structure
!---------------------------------------------------------------------------------
subroutine destory_seam(self)
type(oft_seam), intent(inout) :: self
integer(i4) :: i
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%leo))DEALLOCATE(self%leo)
IF(ASSOCIATED(self%own))DEALLOCATE(self%own)
self%nown=-1
IF(ASSOCIATED(self%be))DEALLOCATE(self%be)
IF(ASSOCIATED(self%lbe))DEALLOCATE(self%lbe)
IF(ASSOCIATED(self%lie))DEALLOCATE(self%lie)
IF(ASSOCIATED(self%kle))DEALLOCATE(self%kle)
IF(ASSOCIATED(self%lle))DEALLOCATE(self%lle)
IF(ASSOCIATED(self%rep_counts))DEALLOCATE(self%rep_counts)
IF(ASSOCIATED(self%rep_displs))DEALLOCATE(self%rep_displs)
IF(ASSOCIATED(self%rep_perm))DEALLOCATE(self%rep_perm)
self%replicated=.FALSE.
!---
IF(ASSOCIATED(self%send))THEN
  DO i=0,self%nproc_con
    IF(ASSOCIATED(self%send(i)%v))THEN
      DEALLOCATE(self%send(i)%v)
      DEALLOCATE(self%recv(i)%v)
    END IF
  END DO
  DEALLOCATE(self%send,self%recv)
END IF
IF(ASSOCIATED(self%csend))THEN
  DO i=0,self%nproc_con
    IF(ASSOCIATED(self%csend(i)%v))THEN
      DEALLOCATE(self%csend(i)%v)
      DEALLOCATE(self%crecv(i)%v)
    END IF
  END DO
  DEALLOCATE(self%csend,self%crecv)
END IF
DEBUG_STACK_POP
end subroutine destory_seam
!---------------------------------------------------------------------------------
!> Setup a replicated seam, where every processor stores all `n` entries of a vector and
!! each entry is owned by exactly one processor
!!
!! Replicated seams are used for problems where the full set of unknowns is stored on every
!! processor, but matrices and solvers are decomposed by row ownership (eg. ThinCurr). Instead
!! of domain-to-domain linkage, stitching uses collective communication over all processors:
!! - Method 0 (synchronize by ownership): owned entries are gathered to all processors
!!   (`MPI_Allgatherv`)
!! - Methods 1 and 2 (sum contributions): contributions from all processors are summed on the
!!   owner (`MPI_Reduce_scatter`) and then gathered to all processors (`MPI_Allgatherv`)
!!
!! In both cases every value is computed once, on its owner, so results are bitwise identical on
!! all processors. For compatibility with code that inspects ownership directly, every entry is
!! flagged as a boundary element (`nbe=n`, `be=.TRUE.`, `lbe(i)=i`) and `leo` marks owned entries,
!! so the usual owned-only reductions apply. No domain-to-domain linkage (`kle`, `lle`, `proc_con`)
!! is created.
!!
!! @note `owner` must be identical on all processors, which is checked here.
!---------------------------------------------------------------------------------
SUBROUTINE seam_setup_replicated(self,n,owner)
TYPE(oft_seam), INTENT(inout) :: self !< Seam structure (existing contents are destroyed)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), TARGET, CONTIGUOUS, INTENT(in) :: owner(:) !< Owning processor (0-based rank) of each entry [n]
INTEGER(i4) :: i,hash_loc
INTEGER(i4), ALLOCATABLE :: fill(:)
DEBUG_STACK_PUSH
IF(SIZE(owner)/=n)CALL oft_abort('Size of "owner" does not match "n"','seam_setup_replicated',__FILE__)
IF(n>0)THEN
  IF(MINVAL(owner)<0.OR.MAXVAL(owner)>=oft_env%nprocs)CALL oft_abort('Invalid owner rank', &
    'seam_setup_replicated',__FILE__)
END IF
!---Check that ownership is identical on all processors
hash_loc=INT(oft_simple_hash(C_LOC(owner),INT(4,c_long)*INT(n,c_long),starting_value=INT(n,c_int)),i4)
IF(.NOT.replicated_hash_match(hash_loc))CALL oft_abort('Ownership differs between processors', &
  'seam_setup_replicated',__FILE__)
!---Reset seam
CALL destory_seam(self)
self%full=.FALSE.
self%skip=.FALSE.
self%replicated=.TRUE.
self%nproc_con=0
self%proc_split=0
self%nie=0
self%nle=0
self%nbe=n
self%nbemax=n
ALLOCATE(self%be(n),self%lbe(n),self%leo(n),self%lie(0))
self%be=.TRUE.
DO i=1,n
  self%lbe(i)=i
  self%leo(i)=(owner(i)==oft_env%rank)
END DO
!---Owner ordering
ALLOCATE(self%rep_counts(oft_env%nprocs),self%rep_displs(oft_env%nprocs),self%rep_perm(n))
self%rep_counts=0
DO i=1,n
  self%rep_counts(owner(i)+1)=self%rep_counts(owner(i)+1)+1
END DO
self%rep_displs(1)=0
DO i=2,oft_env%nprocs
  self%rep_displs(i)=self%rep_displs(i-1)+self%rep_counts(i-1)
END DO
ALLOCATE(fill(oft_env%nprocs))
fill=0
DO i=1,n
  fill(owner(i)+1)=fill(owner(i)+1)+1
  self%rep_perm(self%rep_displs(owner(i)+1)+fill(owner(i)+1))=i
END DO
DEALLOCATE(fill)
DEBUG_STACK_POP
END SUBROUTINE seam_setup_replicated
!---------------------------------------------------------------------------------
!> Validate a replicated seam structure
!---------------------------------------------------------------------------------
SUBROUTINE replicated_seam_check(self)
TYPE(oft_seam), INTENT(inout) :: self !< Seam structure
INTEGER(i4) :: i,nowned
INTEGER(i4), ALLOCATABLE :: counts(:)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%rep_counts))CALL oft_abort('Replicated seam not initialized', &
  'replicated_seam_check',__FILE__)
IF(SUM(self%rep_counts)/=self%nbe)CALL oft_abort('Owner counts do not sum to vector size', &
  'replicated_seam_check',__FILE__)
nowned=COUNT(self%leo)
IF(nowned/=self%rep_counts(oft_env%rank+1))CALL oft_abort('Local ownership does not match owner counts', &
  'replicated_seam_check',__FILE__)
ALLOCATE(counts(oft_env%nprocs))
counts=0
counts(oft_env%rank+1)=nowned
counts=oft_mpi_sum(counts,oft_env%nprocs)
IF(ANY(counts/=self%rep_counts))CALL oft_abort('Owner counts differ between processors', &
  'replicated_seam_check',__FILE__)
DO i=1,self%rep_counts(oft_env%rank+1)
  IF(.NOT.self%leo(self%rep_perm(self%rep_displs(oft_env%rank+1)+i)))CALL oft_abort( &
    'Owner permutation does not match local ownership','replicated_seam_check',__FILE__)
END DO
DEALLOCATE(counts)
DEBUG_STACK_POP
END SUBROUTINE replicated_seam_check
!---------------------------------------------------------------------------------
!> real(r8) implementation of \ref oft_stitching::oft_global_stitch for replicated seams
!---------------------------------------------------------------------------------
SUBROUTINE replicated_stitch_r8(self,a,up_method)
TYPE(oft_seam), INTENT(inout) :: self !< Seam structure
REAL(r8), INTENT(inout) :: a(:) !< Local data to be stitched
INTEGER(i4), INTENT(in) :: up_method !< Stitching method (0, 1, or 2)
INTEGER(i4) :: i,n,nloc,ierr
INTEGER(i8) :: timein
REAL(r8), ALLOCATABLE, TARGET :: buf(:),loc(:)
IF(oft_env%nprocs==1)RETURN
DEBUG_STACK_PUSH
IF(up_method>2.OR.up_method<0)CALL oft_abort('Invalid stitch method','replicated_stitch_r8',__FILE__)
n=self%nbe
IF(SIZE(a)/=n)CALL oft_abort('Array size does not match seam','replicated_stitch_r8',__FILE__)
nloc=self%rep_counts(oft_env%rank+1)
#ifdef HAVE_MPI
timein=oft_time_i8()
ALLOCATE(buf(n),loc(MAX(1,nloc)))
IF(up_method==0)THEN
  !---Owned entries only
  !$omp parallel do simd if(nloc>OFT_OMP_VTHRESH)
  DO i=1,nloc
    loc(i)=a(self%rep_perm(self%rep_displs(oft_env%rank+1)+i))
  END DO
ELSE
  !---Sum contributions from all processors on the owner
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    buf(i)=a(self%rep_perm(i))
  END DO
  CALL MPI_REDUCE_SCATTER(buf,loc,self%rep_counts,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_REDUCE_SCATTER','replicated_stitch_r8',__FILE__)
END IF
!---Gather owned values to all processors
CALL MPI_ALLGATHERV(loc,nloc,OFT_MPI_R8,buf,self%rep_counts,self%rep_displs,OFT_MPI_R8,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLGATHERV','replicated_stitch_r8',__FILE__)
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  a(self%rep_perm(i))=buf(i)
END DO
!---Check consistency (debug only)
IF(oft_env%debug>0)CALL replicated_check_hash(C_LOC(buf),INT(8,c_long)*INT(n,c_long))
DEALLOCATE(buf,loc)
!$omp atomic
comm_times(1)=comm_times(1)+oft_time_diff(timein)
#endif
DEBUG_STACK_POP
END SUBROUTINE replicated_stitch_r8
!---------------------------------------------------------------------------------
!> complex(c8) implementation of \ref oft_stitching::oft_global_stitch for replicated seams
!---------------------------------------------------------------------------------
SUBROUTINE replicated_stitch_c8(self,a,up_method)
TYPE(oft_seam), INTENT(inout) :: self !< Seam structure
COMPLEX(c8), INTENT(inout) :: a(:) !< Local data to be stitched
INTEGER(i4), INTENT(in) :: up_method !< Stitching method (0, 1, or 2)
INTEGER(i4) :: i,n,nloc,ierr
INTEGER(i8) :: timein
COMPLEX(c8), ALLOCATABLE, TARGET :: buf(:),loc(:)
IF(oft_env%nprocs==1)RETURN
DEBUG_STACK_PUSH
IF(up_method>2.OR.up_method<0)CALL oft_abort('Invalid stitch method','replicated_stitch_c8',__FILE__)
n=self%nbe
IF(SIZE(a)/=n)CALL oft_abort('Array size does not match seam','replicated_stitch_c8',__FILE__)
nloc=self%rep_counts(oft_env%rank+1)
#ifdef HAVE_MPI
timein=oft_time_i8()
ALLOCATE(buf(n),loc(MAX(1,nloc)))
IF(up_method==0)THEN
  !---Owned entries only
  !$omp parallel do simd if(nloc>OFT_OMP_VTHRESH)
  DO i=1,nloc
    loc(i)=a(self%rep_perm(self%rep_displs(oft_env%rank+1)+i))
  END DO
ELSE
  !---Sum contributions from all processors on the owner
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO i=1,n
    buf(i)=a(self%rep_perm(i))
  END DO
  CALL MPI_REDUCE_SCATTER(buf,loc,self%rep_counts,OFT_MPI_C8,MPI_SUM,oft_env%COMM,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_REDUCE_SCATTER','replicated_stitch_c8',__FILE__)
END IF
!---Gather owned values to all processors
CALL MPI_ALLGATHERV(loc,nloc,OFT_MPI_C8,buf,self%rep_counts,self%rep_displs,OFT_MPI_C8,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLGATHERV','replicated_stitch_c8',__FILE__)
!$omp parallel do simd if(n>OFT_OMP_VTHRESH)
DO i=1,n
  a(self%rep_perm(i))=buf(i)
END DO
!---Check consistency (debug only)
IF(oft_env%debug>0)CALL replicated_check_hash(C_LOC(buf),INT(16,c_long)*INT(n,c_long))
DEALLOCATE(buf,loc)
!$omp atomic
comm_times(1)=comm_times(1)+oft_time_diff(timein)
#endif
DEBUG_STACK_POP
END SUBROUTINE replicated_stitch_c8
!---------------------------------------------------------------------------------
!> Abort if data differs between processors (used to check replicated stitching)
!---------------------------------------------------------------------------------
SUBROUTINE replicated_check_hash(data_ptr,nbytes)
TYPE(c_ptr), INTENT(in) :: data_ptr !< Location of data
INTEGER(c_long), INTENT(in) :: nbytes !< Size of data in bytes
IF(.NOT.replicated_hash_match(INT(oft_simple_hash(data_ptr,nbytes),i4)))CALL oft_abort( &
  'Replicated data differs between processors','replicated_check_hash',__FILE__)
END SUBROUTINE replicated_check_hash
!---------------------------------------------------------------------------------
!> Check if a hash value is the same on all processors
!!
!! @result Hash values match on all processors?
!---------------------------------------------------------------------------------
FUNCTION replicated_hash_match(hash_loc) RESULT(match)
INTEGER(i4), INTENT(in) :: hash_loc !< Local hash value
LOGICAL :: match
INTEGER(i8) :: hash_ext(2)
hash_ext=[INT(hash_loc,i8),-INT(hash_loc,i8)]
hash_ext=oft_mpi_max(hash_ext,2)
match=(hash_ext(1)==INT(hash_loc,i8)).AND.(-hash_ext(2)==INT(hash_loc,i8))
END FUNCTION replicated_hash_match
end module oft_stitching
