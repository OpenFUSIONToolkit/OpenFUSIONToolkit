!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Module body template for scatter plan operations, included by oft_scatter.F90
!
! Required macros (undefined by the including file after each instantiation):
!   ST     - Value type (eg. REAL(r8))
!   SZERO  - Zero of type ST
!   SKIND  - Value type code passed to scatter_plan_inflight (1: real, 2: complex)
!   SSBUF  - Send buffer component of oft_scatter_plan for type ST
!   SRBUF  - Receive buffer component of oft_scatter_plan for type ST
!   SMPI   - MPI datatype for type ST
!
! The including file provides the MODULE and END MODULE statements. Specific procedures
! are private; only the generic names are public, so generics of the same name from
! each instantiation merge in oft_scatter.
!---------------------------------------------------------------------------------
USE oft_local
USE oft_base
USE oft_scatter_base
IMPLICIT NONE
PRIVATE
PUBLIC :: scatter_begin, scatter_end, scatter_wait_recv, scatter_unpack, scatter_wait_send
PUBLIC :: scatter_local
INTERFACE scatter_begin
  MODULE PROCEDURE s_begin
END INTERFACE scatter_begin
INTERFACE scatter_end
  MODULE PROCEDURE s_end
END INTERFACE scatter_end
INTERFACE scatter_wait_recv
  MODULE PROCEDURE s_wait_recv
END INTERFACE scatter_wait_recv
INTERFACE scatter_unpack
  MODULE PROCEDURE s_unpack
END INTERFACE scatter_unpack
INTERFACE scatter_wait_send
  MODULE PROCEDURE s_wait_send
END INTERFACE scatter_wait_send
INTERFACE scatter_local
  MODULE PROCEDURE s_local
END INTERFACE scatter_local
CONTAINS
!---------------------------------------------------------------------------------
!> Start a scatter operation: post receives, pack values from `x` and post sends
!!
!! Forward (default): values `x(lsend)` are sent to the neighbours and `x(lsrc)` are
!! packed for the local part. Reverse: values `x(lrecv)` are sent back to the
!! neighbours they came from and `x(ldst)` are packed for the local part. Packed
!! entries with `mask` false are sent as zero. With `nolocal` the local part is not
!! packed (and must not be unpacked).
!!
!! Must be completed by @ref s_end "scatter_end" (or by @ref s_wait_recv
!! "scatter_wait_recv" and @ref s_wait_send "scatter_wait_send") before `x` is
!! released or another operation is started on the plan.
!---------------------------------------------------------------------------------
SUBROUTINE s_begin(plan,x,reverse,mask,nolocal)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
ST, INTENT(in) :: x(:) !< Source values
LOGICAL, OPTIONAL, INTENT(in) :: reverse !< Reverse direction? (default forward)
LOGICAL, OPTIONAL, INTENT(in) :: mask(:) !< Packed entries to send (zero otherwise) [nsend+nloc] or [nrecv+nloc] (optional)
LOGICAL, OPTIONAL, INTENT(in) :: nolocal !< Skip local part? (optional)
LOGICAL :: rev,doloc
INTEGER(i4) :: i,j,ierr,nout,nin
rev=.FALSE.
IF(PRESENT(reverse))rev=reverse
doloc=.TRUE.
IF(PRESENT(nolocal))doloc=(.NOT.nolocal)
IF(rev)THEN
  CALL scatter_plan_inflight(plan,-1,SKIND)
ELSE
  CALL scatter_plan_inflight(plan,1,SKIND)
END IF
IF(.NOT.ALLOCATED(plan%SSBUF))ALLOCATE(plan%SSBUF(MAX(1,plan%nsend+plan%nloc)))
IF(.NOT.ALLOCATED(plan%SRBUF))ALLOCATE(plan%SRBUF(MAX(1,plan%nrecv+plan%nloc)))
IF(rev)THEN
  nout=plan%nrecv
  nin=plan%nsend
ELSE
  nout=plan%nsend
  nin=plan%nrecv
END IF
!---Post receives
#ifdef HAVE_MPI
DO j=1,plan%nnbr
  IF(rev)THEN
    CALL post_recv(j,plan%ksend(j),plan%ksend(j+1)-plan%ksend(j),plan%SSBUF)
  ELSE
    CALL post_recv(j,plan%krecv(j),plan%krecv(j+1)-plan%krecv(j),plan%SRBUF)
  END IF
END DO
#endif
!---Pack values
IF(rev)THEN
  CALL pack_vals(nout,plan%lrecv,plan%SRBUF)
  IF(doloc)CALL pack_vals(plan%nloc,plan%ldst,plan%SSBUF(plan%nsend+1:),nout)
ELSE
  CALL pack_vals(nout,plan%lsend,plan%SSBUF)
  IF(doloc)CALL pack_vals(plan%nloc,plan%lsrc,plan%SRBUF(plan%nrecv+1:),nout)
END IF
!---Post sends from the calling thread (MPI is initialized with MPI_THREAD_FUNNELED)
#ifdef HAVE_MPI
DO j=1,plan%nnbr
  IF(rev)THEN
    CALL post_send(j,plan%krecv(j),plan%krecv(j+1)-plan%krecv(j),plan%SRBUF)
  ELSE
    CALL post_send(j,plan%ksend(j),plan%ksend(j+1)-plan%ksend(j),plan%SSBUF)
  END IF
END DO
#endif
CONTAINS
!> Pack values into a buffer, applying the mask (offset `moff` into `mask`)
SUBROUTINE pack_vals(n,idx,buf,moff)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: idx(:) !< Source index of each entry
ST, INTENT(inout) :: buf(:) !< Buffer
INTEGER(i4), OPTIONAL, INTENT(in) :: moff !< Offset of these entries in `mask`
INTEGER(i4) :: k,mo
mo=0
IF(PRESENT(moff))mo=moff
IF(PRESENT(mask))THEN
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO k=1,n
    IF(mask(mo+k))THEN
      buf(k)=x(idx(k))
    ELSE
      buf(k)=SZERO
    END IF
  END DO
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
  DO k=1,n
    buf(k)=x(idx(k))
  END DO
END IF
END SUBROUTINE pack_vals
#ifdef HAVE_MPI
!> Post a receive from neighbour `jn` into `buf(js:js+n-1)`
SUBROUTINE post_recv(jn,js,n,buf)
INTEGER(i4), INTENT(in) :: jn !< Neighbour
INTEGER(i4), INTENT(in) :: js !< First buffer entry
INTEGER(i4), INTENT(in) :: n !< Number of entries
ST, INTENT(inout) :: buf(:) !< Buffer
IF(n==0)THEN
  plan%rreq(jn)=MPI_REQUEST_NULL
  RETURN
END IF
CALL MPI_IRECV(buf(js:js+n-1),n,SMPI,plan%nbr(jn),plan%tag,oft_env%comm,plan%rreq(jn),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','scatter_begin',__FILE__)
END SUBROUTINE post_recv
!> Post a send of `buf(js:js+n-1)` to neighbour `jn`
SUBROUTINE post_send(jn,js,n,buf)
INTEGER(i4), INTENT(in) :: jn !< Neighbour
INTEGER(i4), INTENT(in) :: js !< First buffer entry
INTEGER(i4), INTENT(in) :: n !< Number of entries
ST, INTENT(inout) :: buf(:) !< Buffer
IF(n==0)THEN
  plan%sreq(jn)=MPI_REQUEST_NULL
  RETURN
END IF
CALL MPI_ISEND(buf(js:js+n-1),n,SMPI,plan%nbr(jn),plan%tag,oft_env%comm,plan%sreq(jn),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','scatter_begin',__FILE__)
END SUBROUTINE post_send
#endif
END SUBROUTINE s_begin
!---------------------------------------------------------------------------------
!> Wait for all receives of the operation in flight on a plan
!!
!! The received values are then available in the plan's buffers (see
!! @ref oft_scatter_base::oft_scatter_plan "oft_scatter_plan").
!---------------------------------------------------------------------------------
SUBROUTINE s_wait_recv(plan,x)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
ST, INTENT(in) :: x(:) !< Any array of the value type of the operation (selects the type, not accessed)
INTEGER(i4) :: ierr
IF(plan%inflight==0)CALL oft_abort('No operation in flight','scatter_wait_recv',__FILE__)
IF(plan%inflight_kind/=SKIND)CALL oft_abort('Value type differs from operation in flight', &
  'scatter_wait_recv',__FILE__)
#ifdef HAVE_MPI
IF(plan%nnbr>0)THEN
  CALL oft_mpi_waitall(plan%nnbr,plan%rreq,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','scatter_wait_recv',__FILE__)
END IF
#endif
END SUBROUTINE s_wait_recv
!---------------------------------------------------------------------------------
!> Unpack received (and local) values of the operation in flight into `y`
!!
!! Forward: `y(lrecv(i))` and `y(ldst(k))` are set or incremented. Reverse:
!! `y(lsend(i))` and `y(lsrc(k))` (times `wt(k)`, skipping local entries with zero
!! weight). Entries are processed in order (deterministic sums) unless
!! the destination indices are unique, in which case the loops are threaded. Must
!! follow @ref s_wait_recv "scatter_wait_recv".
!---------------------------------------------------------------------------------
SUBROUTINE s_unpack(plan,y,mode,nolocal)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
ST, INTENT(inout) :: y(:) !< Destination values
INTEGER(i4), INTENT(in) :: mode !< @ref oft_scatter_base::oft_scatter_insert "OFT_SCATTER_INSERT" or `OFT_SCATTER_ADD`
LOGICAL, OPTIONAL, INTENT(in) :: nolocal !< Skip local part? (must match `scatter_begin`)
LOGICAL :: doloc,par
doloc=.TRUE.
IF(PRESENT(nolocal))doloc=(.NOT.nolocal)
IF(plan%inflight_kind/=SKIND)CALL oft_abort('Value type differs from operation in flight', &
  'scatter_unpack',__FILE__)
IF(plan%inflight==1)THEN
  par=plan%recv_unique
  CALL unpack_vals(plan%nrecv,plan%lrecv,plan%SRBUF,.FALSE.)
  IF(doloc)CALL unpack_vals(plan%nloc,plan%ldst,plan%SRBUF(plan%nrecv+1:),.FALSE.)
ELSE IF(plan%inflight==-1)THEN
  par=plan%send_unique
  CALL unpack_vals(plan%nsend,plan%lsend,plan%SSBUF,.FALSE.)
  IF(doloc)CALL unpack_vals(plan%nloc,plan%lsrc,plan%SSBUF(plan%nsend+1:),.TRUE.)
ELSE
  CALL oft_abort('No operation in flight','scatter_unpack',__FILE__)
END IF
CONTAINS
!> Unpack buffer entries into `y`, with weights for local entries
SUBROUTINE unpack_vals(n,idx,buf,weighted)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: idx(:) !< Destination index of each entry
ST, INTENT(in) :: buf(:) !< Buffer
LOGICAL, INTENT(in) :: weighted !< Apply local weights (if allocated)?
INTEGER(i4) :: k
IF(weighted.AND.ALLOCATED(plan%wt))THEN
  CALL apply_weighted(n,idx,plan%wt,buf,y,mode,par)
ELSE IF(mode==OFT_SCATTER_ADD)THEN
  IF(par)THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
    DO k=1,n
      y(idx(k))=y(idx(k))+buf(k)
    END DO
  ELSE
    DO k=1,n
      y(idx(k))=y(idx(k))+buf(k)
    END DO
  END IF
ELSE
  IF(par)THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
    DO k=1,n
      y(idx(k))=buf(k)
    END DO
  ELSE
    DO k=1,n
      y(idx(k))=buf(k)
    END DO
  END IF
END IF
END SUBROUTINE unpack_vals
END SUBROUTINE s_unpack
!---------------------------------------------------------------------------------
!> Wait for all sends of the operation in flight on a plan and mark it complete
!---------------------------------------------------------------------------------
SUBROUTINE s_wait_send(plan,x)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
ST, INTENT(in) :: x(:) !< Any array of the value type of the operation (selects the type, not accessed)
INTEGER(i4) :: ierr
#ifdef HAVE_MPI
IF(plan%nnbr>0)THEN
  CALL oft_mpi_waitall(plan%nnbr,plan%sreq,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','scatter_wait_send',__FILE__)
END IF
#endif
CALL scatter_plan_inflight(plan,0,SKIND)
END SUBROUTINE s_wait_send
!---------------------------------------------------------------------------------
!> Complete a scatter operation: wait for receives, unpack into `y` and wait for sends
!!
!! See @ref s_unpack "scatter_unpack" for the destination entries.
!---------------------------------------------------------------------------------
SUBROUTINE s_end(plan,y,mode,nolocal)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
ST, INTENT(inout) :: y(:) !< Destination values
INTEGER(i4), INTENT(in) :: mode !< @ref oft_scatter_base::oft_scatter_insert "OFT_SCATTER_INSERT" or `OFT_SCATTER_ADD`
LOGICAL, OPTIONAL, INTENT(in) :: nolocal !< Skip local part? (must match `scatter_begin`)
CALL s_wait_recv(plan,y)
CALL s_unpack(plan,y,mode,nolocal)
CALL s_wait_send(plan,y)
END SUBROUTINE s_end
!---------------------------------------------------------------------------------
!> Apply only the local part of a plan, without communication or buffers
!!
!! Forward: \f$ y_{ldst(k)} \leftarrow x_{lsrc(k)} \f$. Reverse:
!! \f$ y_{lsrc(k)} \leftarrow w_k x_{ldst(k)} \f$, skipping entries with zero weight
!! (weights are optional). Insert or add. Threaded when the destination indices are
!! unique.
!---------------------------------------------------------------------------------
SUBROUTINE s_local(plan,x,y,mode,reverse)
TYPE(oft_scatter_plan), INTENT(in) :: plan !< Plan
ST, INTENT(in) :: x(:) !< Source values
ST, INTENT(inout) :: y(:) !< Destination values
INTEGER(i4), INTENT(in) :: mode !< @ref oft_scatter_base::oft_scatter_insert "OFT_SCATTER_INSERT" or `OFT_SCATTER_ADD`
LOGICAL, OPTIONAL, INTENT(in) :: reverse !< Reverse direction? (default forward)
LOGICAL :: rev
INTEGER(i4) :: k,n
rev=.FALSE.
IF(PRESENT(reverse))rev=reverse
n=plan%nloc
IF(rev.AND.ALLOCATED(plan%wt))THEN
  CALL apply_weighted_gather(n,plan%ldst,plan%lsrc,plan%wt,x,y,mode,plan%send_unique)
  RETURN
END IF
IF(rev)THEN
  IF(mode==OFT_SCATTER_ADD)THEN
    IF(plan%send_unique)THEN
      !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
      DO k=1,n
        y(plan%lsrc(k))=y(plan%lsrc(k))+x(plan%ldst(k))
      END DO
    ELSE
      DO k=1,n
        y(plan%lsrc(k))=y(plan%lsrc(k))+x(plan%ldst(k))
      END DO
    END IF
  ELSE
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH.AND.plan%send_unique)
    DO k=1,n
      y(plan%lsrc(k))=x(plan%ldst(k))
    END DO
  END IF
ELSE
  IF(mode==OFT_SCATTER_ADD)THEN
    IF(plan%recv_unique)THEN
      !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
      DO k=1,n
        y(plan%ldst(k))=y(plan%ldst(k))+x(plan%lsrc(k))
      END DO
    ELSE
      DO k=1,n
        y(plan%ldst(k))=y(plan%ldst(k))+x(plan%lsrc(k))
      END DO
    END IF
  ELSE
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH.AND.plan%recv_unique)
    DO k=1,n
      y(plan%ldst(k))=x(plan%lsrc(k))
    END DO
  END IF
END IF
END SUBROUTINE s_local
!---------------------------------------------------------------------------------
!> Weighted unpack of buffer entries, \f$ y_{idx(k)} \leftarrow w_k b_k \f$ (zero weights skipped)
!---------------------------------------------------------------------------------
SUBROUTINE apply_weighted(n,idx,wt,buf,y,mode,par)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: idx(:) !< Destination index of each entry [n]
REAL(r8), INTENT(in) :: wt(:) !< Weight of each entry [n]
ST, INTENT(in) :: buf(:) !< Values [n]
ST, INTENT(inout) :: y(:) !< Destination values
INTEGER(i4), INTENT(in) :: mode !< Insert or add
LOGICAL, INTENT(in) :: par !< Destination indices are unique (thread)?
INTEGER(i4) :: k
IF(mode==OFT_SCATTER_ADD)THEN
  IF(par)THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
    DO k=1,n
      IF(wt(k)/=0.d0)y(idx(k))=y(idx(k))+wt(k)*buf(k)
    END DO
  ELSE
    DO k=1,n
      IF(wt(k)/=0.d0)y(idx(k))=y(idx(k))+wt(k)*buf(k)
    END DO
  END IF
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH.AND.par)
  DO k=1,n
    IF(wt(k)/=0.d0)y(idx(k))=wt(k)*buf(k)
  END DO
END IF
END SUBROUTINE apply_weighted
!---------------------------------------------------------------------------------
!> Weighted local scatter, \f$ y_{idst(k)} \leftarrow w_k x_{isrc(k)} \f$ (zero weights skipped)
!---------------------------------------------------------------------------------
SUBROUTINE apply_weighted_gather(n,isrc,idst,wt,x,y,mode,par)
INTEGER(i4), INTENT(in) :: n !< Number of entries
INTEGER(i4), INTENT(in) :: isrc(:) !< Source index of each entry [n]
INTEGER(i4), INTENT(in) :: idst(:) !< Destination index of each entry [n]
REAL(r8), INTENT(in) :: wt(:) !< Weight of each entry [n]
ST, INTENT(in) :: x(:) !< Source values
ST, INTENT(inout) :: y(:) !< Destination values
INTEGER(i4), INTENT(in) :: mode !< Insert or add
LOGICAL, INTENT(in) :: par !< Destination indices are unique (thread)?
INTEGER(i4) :: k
IF(mode==OFT_SCATTER_ADD)THEN
  IF(par)THEN
    !$omp parallel do simd if(n>OFT_OMP_VTHRESH)
    DO k=1,n
      IF(wt(k)/=0.d0)y(idst(k))=y(idst(k))+wt(k)*x(isrc(k))
    END DO
  ELSE
    DO k=1,n
      IF(wt(k)/=0.d0)y(idst(k))=y(idst(k))+wt(k)*x(isrc(k))
    END DO
  END IF
ELSE
  !$omp parallel do simd if(n>OFT_OMP_VTHRESH.AND.par)
  DO k=1,n
    IF(wt(k)/=0.d0)y(idst(k))=wt(k)*x(isrc(k))
  END DO
END IF
END SUBROUTINE apply_weighted_gather
