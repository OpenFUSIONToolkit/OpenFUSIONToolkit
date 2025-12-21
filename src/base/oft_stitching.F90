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
  INTEGER(i4) :: nbemax = 0 !< Maximum number of boundary elements
  INTEGER(i4) :: nbe = 0 !< Number of local boundary elements
  INTEGER(i4) :: nie = 0 !< Number of interior elements
  INTEGER(i4) :: nle = 0 !< Number of local boundary elements
  INTEGER(i4) :: nproc_con = 0 !< Number of connections to other processors
  INTEGER(i4) :: proc_split = 0 !< Location of self in processor list
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: leo => NULL() !< List of boundary element ownership
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
!---------------------------------------------------------------------------------
INTERFACE oft_global_stitch
  MODULE PROCEDURE global_stitch_r8
  MODULE PROCEDURE global_stitch_c8
END INTERFACE oft_global_stitch
CONTAINS
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
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+a(i)
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-a(self%lbe(i))
  end do
  !$omp end parallel
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
  !$omp parallel reduction(+:c) if(n>OFT_OMP_VTHRESH)
  !$omp do simd
  do i=1,n
    c=c+a(i)
  end do
  !$omp do simd
  do i=1,self%nbe
    if(.NOT.self%leo(i))c=c-a(self%lbe(i))
  end do
  !$omp end parallel
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
IF(self%skip.OR.self%nbe==0)RETURN
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
IF(self%skip.OR.self%nbe==0)RETURN
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
IF(ASSOCIATED(self%be))DEALLOCATE(self%be)
IF(ASSOCIATED(self%lbe))DEALLOCATE(self%lbe)
IF(ASSOCIATED(self%lie))DEALLOCATE(self%lie)
IF(ASSOCIATED(self%kle))DEALLOCATE(self%kle)
IF(ASSOCIATED(self%lle))DEALLOCATE(self%lle)
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
end module oft_stitching
