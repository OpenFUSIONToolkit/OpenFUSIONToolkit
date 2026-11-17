!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_scatter.F90
!
!> Index sets and scatter plans: general gather/scatter of values between arrays on
!! one or several processes
!!
!! A scatter plan (@ref oft_scatter_base::oft_scatter_plan "oft_scatter_plan") moves
!! values from a source array `x` to a destination array `y`, which may be the same
!! array or arrays of different size (eg. a distributed vector and a local subdomain
!! vector). It is made of
!! - a remote part: for each neighbour process, a list of source indices whose values
!!   are sent (`lsend`, CSR by neighbour through `ksend`) and a list of destination
!!   indices that receive the values sent by that neighbour, in the same order as the
!!   neighbour's send list (`lrecv`, CSR through `krecv`);
!! - a local part: source indices `lsrc` and destination indices `ldst`, with
!!   optional weights `wt`, for values that stay on the process.
!!
!! Operations (generic over real and complex values, see @ref oft_scatter):
!! - forward: \f$ y_{lrecv} \leftarrow x_{lsend} \f$ (from neighbours) and
!!   \f$ y_{ldst} \leftarrow x_{lsrc} \f$, with insert or add (`OFT_SCATTER_INSERT`,
!!   `OFT_SCATTER_ADD`);
!! - reverse: \f$ x_{lsend} \leftarrow y_{lrecv} \f$ (to neighbours) and
!!   \f$ x_{lsrc} \leftarrow w y_{ldst} \f$, with insert or add (eg. accumulating
!!   copies on their owner).
!!
!! Weights apply to the local part in the reverse direction only (\f$ R^T W \f$ for a
!! restriction \f$ R \f$ given by the forward direction). Entries whose weight is zero
!! are skipped (a write mask), so 0/1 weights restrict which source entries are written
!! back (restricted Schwarz), unit weights everywhere give plain additive combinations
!! and fractional weights partition of unity combinations.
!!
!! Communication is split into `scatter_begin` (pack, post receives and sends) and
!! `scatter_end` (wait, unpack, wait for sends), so work can overlap transfers;
!! `scatter_wait_recv`, `scatter_unpack` and `scatter_wait_send` expose the steps of
!! `scatter_end` for callers that combine received values themselves (eg. summation
!! stitching, @ref oft_stitching::oft_global_stitch "oft_global_stitch").
!!
!! All index lists are flat arrays and all send/receive buffers are flat and
!! persistent (allocated on first use for each value type), so a plan does no
!! allocation per call once it has been used.
!!
!! MPI rules:
!! - MPI calls are only made from the calling thread (MPI is initialized with
!!   `MPI_THREAD_FUNNELED`); pack and unpack loops are threaded internally, so plan
!!   operations must be called from serial code.
!! - Each plan carries an MPI tag (`tag`, see `OFT_TAG_*`). At most one operation may
!!   be in flight on a plan at a time (checked). Operations on different plans that
!!   share a tag must be started in the same order on all processes, which holds for
!!   operations that complete before returning or that are called collectively in
!!   program order: MPI does not reorder messages with the same source, tag and
!!   communicator.
!!
!! An index set (@ref oft_scatter_base::oft_index_set "oft_index_set") is a sorted list
!! of unique local indices, with growth along the edges of a graph for overlapping
!! subdomains.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_base
!---------------------------------------------------------------------------------
#include "local.h"
!---------------------------------------------------------------------------------
!> Index set and scatter plan types, setup and type-independent operations
!---------------------------------------------------------------------------------
MODULE oft_scatter_base
USE oft_local
USE oft_base
USE oft_sort, ONLY: sort_array
IMPLICIT NONE
PRIVATE
!> Replace destination values
INTEGER(i4), PUBLIC, PARAMETER :: OFT_SCATTER_INSERT = 0
!> Add to destination values
INTEGER(i4), PUBLIC, PARAMETER :: OFT_SCATTER_ADD = 1
!> MPI tag for global stitching (@ref oft_stitching::oft_global_stitch "oft_global_stitch")
INTEGER(i4), PUBLIC, PARAMETER :: OFT_TAG_STITCH = 101
!> Default MPI tag for general scatter plans
INTEGER(i4), PUBLIC, PARAMETER :: OFT_TAG_SCATTER = 102
!> MPI tag for setup-time exchanges of index and matrix data between processes
INTEGER(i4), PUBLIC, PARAMETER :: OFT_TAG_SETUP = 103
!---------------------------------------------------------------------------------
!> Sorted set of unique local indices
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_index_set
  INTEGER(i4) :: n = 0 !< Number of indices in set
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ind !< Indices, in increasing order [n]
END TYPE oft_index_set
!---------------------------------------------------------------------------------
!> Plan for gathering/scattering values between arrays, locally and between processes
!!
!! See the description of oft_scatter.F90 for the operations. Set up with
!! @ref scatter_plan_setup, released with @ref scatter_plan_delete.
!!
!! Value buffers are laid out as `[remote entries (by neighbour) | local entries]`:
!! `rbuf(1:nrecv)` holds values received in a forward operation (entry `i` goes to
!! `lrecv(i)`) and `rbuf(nrecv+k)` the local value for entry `k`; `sbuf(1:nsend)`
!! holds values sent in a forward operation (entry `i` from `lsend(i)`) and receives
!! values in a reverse operation, with `sbuf(nsend+k)` the local value for reverse
!! entry `k`.
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_scatter_plan
  INTEGER(i4) :: tag = OFT_TAG_SCATTER !< MPI tag used by this plan
  INTEGER(i4) :: nnbr = 0 !< Number of neighbour processes
  INTEGER(i4) :: nsend = 0 !< Number of remote send entries (forward direction)
  INTEGER(i4) :: nrecv = 0 !< Number of remote receive entries (forward direction)
  INTEGER(i4) :: nloc = 0 !< Number of local entries
  INTEGER(i4) :: inflight = 0 !< Operation in flight (0: none, 1: forward, -1: reverse)
  INTEGER(i4) :: inflight_kind = 0 !< Value type of operation in flight (1: real, 2: complex)
  LOGICAL :: setup = .FALSE. !< Has the plan been set up?
  LOGICAL :: recv_unique = .TRUE. !< Destination indices (`lrecv` and `ldst`) are unique
  LOGICAL :: send_unique = .TRUE. !< Source indices (`lsend` and `lsrc`) are unique
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: nbr !< Neighbour ranks [nnbr]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ksend !< Pointer into `lsend` by neighbour [nnbr+1]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lsend !< Source index of each send entry [nsend]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: krecv !< Pointer into `lrecv` by neighbour [nnbr+1]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lrecv !< Destination index of each receive entry [nrecv]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lsrc !< Source index of each local entry [nloc]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ldst !< Destination index of each local entry [nloc]
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt !< Weight of each local entry in the reverse direction [nloc] (optional)
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: rsbuf !< Real send buffer [nsend+nloc] (allocated on first use)
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: rrbuf !< Real receive buffer [nrecv+nloc] (allocated on first use)
  COMPLEX(c8), ALLOCATABLE, DIMENSION(:) :: csbuf !< Complex send buffer [nsend+nloc] (allocated on first use)
  COMPLEX(c8), ALLOCATABLE, DIMENSION(:) :: crbuf !< Complex receive buffer [nrecv+nloc] (allocated on first use)
#ifdef OFT_MPI_F08
  TYPE(mpi_request), ALLOCATABLE, DIMENSION(:) :: sreq !< Send requests [nnbr]
  TYPE(mpi_request), ALLOCATABLE, DIMENSION(:) :: rreq !< Receive requests [nnbr]
#else
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: sreq !< Send requests [nnbr]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: rreq !< Receive requests [nnbr]
#endif
END TYPE oft_scatter_plan
PUBLIC scatter_plan_setup, scatter_plan_setup_local, scatter_plan_delete, scatter_plan_check, scatter_plan_inflight
PUBLIC index_set_create, index_set_delete, index_set_grow, index_set_map, index_set_contains
CONTAINS
!---------------------------------------------------------------------------------
!> Set up a scatter plan from its index lists
!!
!! Any previous content of `plan` is released. The remote part needs MPI and matching
!! lists on the neighbours: the number of entries `self` receives from neighbour `j`
!! must equal the number of entries `j` sends to `self` (checked by
!! @ref scatter_plan_check). Neighbours with no entries in either direction are
!! allowed and skipped by all operations.
!---------------------------------------------------------------------------------
SUBROUTINE scatter_plan_setup(plan,nnbr,nbr,ksend,lsend,krecv,lrecv,nloc,lsrc,ldst,wt,tag)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
INTEGER(i4), INTENT(in) :: nnbr !< Number of neighbour processes
INTEGER(i4), INTENT(in) :: nbr(nnbr) !< Neighbour ranks [nnbr]
INTEGER(i4), INTENT(in) :: ksend(nnbr+1) !< Pointer into `lsend` by neighbour (starting at 1) [nnbr+1]
INTEGER(i4), INTENT(in) :: lsend(:) !< Source index of each send entry [ksend(nnbr+1)-1]
INTEGER(i4), INTENT(in) :: krecv(nnbr+1) !< Pointer into `lrecv` by neighbour (starting at 1) [nnbr+1]
INTEGER(i4), INTENT(in) :: lrecv(:) !< Destination index of each receive entry [krecv(nnbr+1)-1]
INTEGER(i4), INTENT(in) :: nloc !< Number of local entries
INTEGER(i4), INTENT(in) :: lsrc(:) !< Source index of each local entry [nloc]
INTEGER(i4), OPTIONAL, INTENT(in) :: ldst(:) !< Destination index of each local entry [nloc] (default: 1..nloc)
REAL(r8), OPTIONAL, INTENT(in) :: wt(:) !< Weight of each local entry in the reverse direction [nloc] (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: tag !< MPI tag (default @ref OFT_TAG_SCATTER)
INTEGER(i4) :: i
CALL scatter_plan_delete(plan)
IF(PRESENT(tag))plan%tag=tag
plan%nnbr=nnbr
plan%nsend=ksend(nnbr+1)-1
plan%nrecv=krecv(nnbr+1)-1
plan%nloc=nloc
ALLOCATE(plan%nbr(nnbr),plan%ksend(nnbr+1),plan%krecv(nnbr+1))
ALLOCATE(plan%lsend(plan%nsend),plan%lrecv(plan%nrecv))
ALLOCATE(plan%lsrc(nloc),plan%ldst(nloc))
ALLOCATE(plan%sreq(MAX(1,nnbr)),plan%rreq(MAX(1,nnbr)))
plan%sreq=MPI_REQUEST_NULL
plan%rreq=MPI_REQUEST_NULL
IF(nnbr>0)THEN
  plan%nbr=nbr
  IF(ksend(1)/=1.OR.krecv(1)/=1)CALL oft_abort('Neighbour pointers must start at 1','scatter_plan_setup',__FILE__)
END IF
plan%ksend=ksend
plan%krecv=krecv
IF(plan%nsend>0)plan%lsend=lsend(1:plan%nsend)
IF(plan%nrecv>0)plan%lrecv=lrecv(1:plan%nrecv)
IF(nloc>0)plan%lsrc=lsrc(1:nloc)
IF(PRESENT(ldst))THEN
  IF(nloc>0)plan%ldst=ldst(1:nloc)
ELSE
  plan%ldst=[(i,i=1,nloc)]
END IF
IF(PRESENT(wt))THEN
  ALLOCATE(plan%wt(nloc))
  IF(nloc>0)plan%wt=wt(1:nloc)
END IF
plan%recv_unique=lists_unique(plan%lrecv,plan%ldst)
plan%send_unique=lists_unique(plan%lsend,plan%lsrc)
plan%setup=.TRUE.
CONTAINS
!> Check whether the entries of two index lists are unique (combined)
FUNCTION lists_unique(la,lb) RESULT(unique)
INTEGER(i4), INTENT(in) :: la(:) !< First list
INTEGER(i4), INTENT(in) :: lb(:) !< Second list
LOGICAL :: unique !< No index occurs twice?
INTEGER(i4) :: j,nmax
LOGICAL, ALLOCATABLE, DIMENSION(:) :: seen
unique=.TRUE.
nmax=0
IF(SIZE(la)>0)nmax=MAX(nmax,MAXVAL(la))
IF(SIZE(lb)>0)nmax=MAX(nmax,MAXVAL(lb))
IF(nmax==0)RETURN
ALLOCATE(seen(nmax))
seen=.FALSE.
DO j=1,SIZE(la)
  IF(seen(la(j)))unique=.FALSE.
  seen(la(j))=.TRUE.
END DO
DO j=1,SIZE(lb)
  IF(seen(lb(j)))unique=.FALSE.
  seen(lb(j))=.TRUE.
END DO
DEALLOCATE(seen)
END FUNCTION lists_unique
END SUBROUTINE scatter_plan_setup
!---------------------------------------------------------------------------------
!> Set up a purely local scatter plan (gather/scatter of a subset with optional weights)
!!
!! Equivalent to @ref scatter_plan_setup without neighbours. With `ldst` absent the
!! plan gathers `x(lsrc(k))` into `y(k)` (and scatters back in reverse).
!---------------------------------------------------------------------------------
SUBROUTINE scatter_plan_setup_local(plan,nloc,lsrc,ldst,wt)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
INTEGER(i4), INTENT(in) :: nloc !< Number of local entries
INTEGER(i4), INTENT(in) :: lsrc(:) !< Source index of each entry [nloc]
INTEGER(i4), OPTIONAL, INTENT(in) :: ldst(:) !< Destination index of each entry [nloc] (default: 1..nloc)
REAL(r8), OPTIONAL, INTENT(in) :: wt(:) !< Weight of each entry in the reverse direction [nloc] (optional)
INTEGER(i4) :: knull(1),lnull(1),nbr(1)
knull=1
lnull=0
nbr=0
CALL scatter_plan_setup(plan,0,nbr,knull,lnull,knull,lnull,nloc,lsrc,ldst,wt)
END SUBROUTINE scatter_plan_setup_local
!---------------------------------------------------------------------------------
!> Release all storage of a scatter plan
!---------------------------------------------------------------------------------
SUBROUTINE scatter_plan_delete(plan)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
IF(plan%inflight/=0)CALL oft_abort('Operation in flight','scatter_plan_delete',__FILE__)
IF(ALLOCATED(plan%nbr))DEALLOCATE(plan%nbr)
IF(ALLOCATED(plan%ksend))DEALLOCATE(plan%ksend)
IF(ALLOCATED(plan%lsend))DEALLOCATE(plan%lsend)
IF(ALLOCATED(plan%krecv))DEALLOCATE(plan%krecv)
IF(ALLOCATED(plan%lrecv))DEALLOCATE(plan%lrecv)
IF(ALLOCATED(plan%lsrc))DEALLOCATE(plan%lsrc)
IF(ALLOCATED(plan%ldst))DEALLOCATE(plan%ldst)
IF(ALLOCATED(plan%wt))DEALLOCATE(plan%wt)
IF(ALLOCATED(plan%rsbuf))DEALLOCATE(plan%rsbuf)
IF(ALLOCATED(plan%rrbuf))DEALLOCATE(plan%rrbuf)
IF(ALLOCATED(plan%csbuf))DEALLOCATE(plan%csbuf)
IF(ALLOCATED(plan%crbuf))DEALLOCATE(plan%crbuf)
IF(ALLOCATED(plan%sreq))DEALLOCATE(plan%sreq)
IF(ALLOCATED(plan%rreq))DEALLOCATE(plan%rreq)
plan%nnbr=0
plan%nsend=0
plan%nrecv=0
plan%nloc=0
plan%recv_unique=.TRUE.
plan%send_unique=.TRUE.
plan%setup=.FALSE.
END SUBROUTINE scatter_plan_delete
!---------------------------------------------------------------------------------
!> Check that the numbers of entries exchanged with each neighbour match on both sides
!!
!! Collective over the plan's neighbours; aborts on mismatch.
!---------------------------------------------------------------------------------
SUBROUTINE scatter_plan_check(plan)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
#ifdef HAVE_MPI
INTEGER(i4) :: j,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: nout,nin
IF(plan%nnbr==0)RETURN
ALLOCATE(nout(2,plan%nnbr),nin(2,plan%nnbr))
DO j=1,plan%nnbr
  nout(:,j)=[plan%ksend(j+1)-plan%ksend(j),plan%krecv(j+1)-plan%krecv(j)]
  CALL MPI_ISEND(nout(:,j),2,OFT_MPI_I4,plan%nbr(j),plan%tag,oft_env%comm,plan%sreq(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','scatter_plan_check',__FILE__)
  CALL MPI_IRECV(nin(:,j),2,OFT_MPI_I4,plan%nbr(j),plan%tag,oft_env%comm,plan%rreq(j),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','scatter_plan_check',__FILE__)
END DO
CALL oft_mpi_waitall(plan%nnbr,plan%rreq,ierr)
CALL oft_mpi_waitall(plan%nnbr,plan%sreq,ierr)
DO j=1,plan%nnbr
  IF((nin(1,j)/=nout(2,j)).OR.(nin(2,j)/=nout(1,j)))THEN
    WRITE(*,*)'SCATTER PLAN ERROR: ',oft_env%rank,plan%nbr(j),nout(:,j),nin(:,j)
    CALL oft_abort('Mismatched transfer sizes','scatter_plan_check',__FILE__)
  END IF
END DO
DEALLOCATE(nout,nin)
#endif
END SUBROUTINE scatter_plan_check
!---------------------------------------------------------------------------------
!> Mark an operation as started on a plan, checking that none is in flight
!---------------------------------------------------------------------------------
SUBROUTINE scatter_plan_inflight(plan,direction,kind)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
INTEGER(i4), INTENT(in) :: direction !< 1 for forward, -1 for reverse, 0 to mark completion
INTEGER(i4), INTENT(in) :: kind !< Value type (1: real, 2: complex)
IF(direction==0)THEN
  IF(plan%inflight==0)CALL oft_abort('No operation in flight','scatter_plan_inflight',__FILE__)
  IF(plan%inflight_kind/=kind)CALL oft_abort('Value type differs from operation in flight', &
    'scatter_plan_inflight',__FILE__)
  plan%inflight=0
  plan%inflight_kind=0
ELSE
  IF(.NOT.plan%setup)CALL oft_abort('Plan not set up','scatter_plan_inflight',__FILE__)
  IF(plan%inflight/=0)CALL oft_abort('Operation already in flight','scatter_plan_inflight',__FILE__)
  plan%inflight=direction
  plan%inflight_kind=kind
END IF
END SUBROUTINE scatter_plan_inflight
!---------------------------------------------------------------------------------
!> Create an index set from a list of indices (sorted, duplicates removed)
!---------------------------------------------------------------------------------
SUBROUTINE index_set_create(set,n,list)
TYPE(oft_index_set), INTENT(inout) :: set !< Index set
INTEGER(i4), INTENT(in) :: n !< Number of indices in list
INTEGER(i4), INTENT(in) :: list(n) !< Indices (any order, may repeat) [n]
INTEGER(i4) :: i,m
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: tmp
CALL index_set_delete(set)
ALLOCATE(tmp(MAX(1,n)))
IF(n>0)THEN
  tmp(1:n)=list
  CALL sort_array(tmp,n)
END IF
m=0
DO i=1,n
  IF(m>0)THEN
    IF(tmp(i)==tmp(m))CYCLE
  END IF
  m=m+1
  tmp(m)=tmp(i)
END DO
set%n=m
ALLOCATE(set%ind(m))
IF(m>0)set%ind=tmp(1:m)
DEALLOCATE(tmp)
END SUBROUTINE index_set_create
!---------------------------------------------------------------------------------
!> Release an index set
!---------------------------------------------------------------------------------
SUBROUTINE index_set_delete(set)
TYPE(oft_index_set), INTENT(inout) :: set !< Index set
IF(ALLOCATED(set%ind))DEALLOCATE(set%ind)
set%n=0
END SUBROUTINE index_set_delete
!---------------------------------------------------------------------------------
!> Grow an index set along the edges of a graph
!!
!! Each level adds every index `lc(k)`, `kr(i)<=k<kr(i+1)`, connected to an index `i`
!! in the set (one layer of overlap per level). If `allowed` is present only indices
!! with `allowed(j)` are added, and indices are only reached through rows `i` with
!! `same(i)==same(j)` when `same` is present (eg. to stay within a group of fields).
!---------------------------------------------------------------------------------
SUBROUTINE index_set_grow(set,nr,kr,lc,levels,allowed,same)
TYPE(oft_index_set), INTENT(inout) :: set !< Index set
INTEGER(i4), INTENT(in) :: nr !< Number of graph rows (indices range over 1..nr)
INTEGER(i4), INTENT(in) :: kr(nr+1) !< Graph row pointer [nr+1]
INTEGER(i4), INTENT(in) :: lc(:) !< Graph column list [kr(nr+1)-1]
INTEGER(i4), INTENT(in) :: levels !< Number of layers to add
LOGICAL, OPTIONAL, INTENT(in) :: allowed(nr) !< Indices that may be added [nr] (optional)
INTEGER(i4), OPTIONAL, INTENT(in) :: same(nr) !< Group of each index, edges must stay in a group [nr] (optional)
INTEGER(i4) :: i,j,k,l,nfront,nnew
LOGICAL, ALLOCATABLE, DIMENSION(:) :: inset
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: front,newf
IF(levels<=0)RETURN
ALLOCATE(inset(nr))
inset=.FALSE.
DO i=1,set%n
  inset(set%ind(i))=.TRUE.
END DO
ALLOCATE(front(MAX(1,set%n)))
nfront=set%n
IF(set%n>0)front(1:nfront)=set%ind
DO l=1,levels
  nnew=0
  ALLOCATE(newf(MAX(1,nr)))
  DO k=1,nfront
    i=front(k)
    DO j=kr(i),kr(i+1)-1
      IF(inset(lc(j)))CYCLE
      IF(PRESENT(allowed))THEN
        IF(.NOT.allowed(lc(j)))CYCLE
      END IF
      IF(PRESENT(same))THEN
        IF(same(lc(j))/=same(i))CYCLE
      END IF
      inset(lc(j))=.TRUE.
      nnew=nnew+1
      newf(nnew)=lc(j)
    END DO
  END DO
  DEALLOCATE(front)
  CALL MOVE_ALLOC(newf,front)
  nfront=nnew
  IF(nnew==0)EXIT
END DO
DEALLOCATE(front)
!---Rebuild sorted list
set%n=COUNT(inset)
IF(ALLOCATED(set%ind))DEALLOCATE(set%ind)
ALLOCATE(set%ind(set%n))
k=0
DO i=1,nr
  IF(.NOT.inset(i))CYCLE
  k=k+1
  set%ind(k)=i
END DO
DEALLOCATE(inset)
END SUBROUTINE index_set_grow
!---------------------------------------------------------------------------------
!> Position of each index in an index set (inverse map)
!!
!! `map(i)` is the position of `i` in `set%ind`, or 0 if `i` is not in the set.
!---------------------------------------------------------------------------------
SUBROUTINE index_set_map(set,n,map)
TYPE(oft_index_set), INTENT(in) :: set !< Index set
INTEGER(i4), INTENT(in) :: n !< Size of index range
INTEGER(i4), INTENT(out) :: map(n) !< Position in set (0 if absent) [n]
INTEGER(i4) :: i
map=0
DO i=1,set%n
  map(set%ind(i))=i
END DO
END SUBROUTINE index_set_map
!---------------------------------------------------------------------------------
!> Check whether an index is in an index set (binary search)
!---------------------------------------------------------------------------------
PURE FUNCTION index_set_contains(set,i) RESULT(pos)
TYPE(oft_index_set), INTENT(in) :: set !< Index set
INTEGER(i4), INTENT(in) :: i !< Index to look for
INTEGER(i4) :: pos !< Position in set (0 if absent)
INTEGER(i4) :: ilow,ihigh,imid
pos=0
ilow=1
ihigh=set%n
DO WHILE(ilow<=ihigh)
  imid=(ilow+ihigh)/2
  IF(set%ind(imid)==i)THEN
    pos=imid
    RETURN
  ELSE IF(set%ind(imid)<i)THEN
    ilow=imid+1
  ELSE
    ihigh=imid-1
  END IF
END DO
END FUNCTION index_set_contains
END MODULE oft_scatter_base
!---------------------------------------------------------------------------------
! Value type instantiations of the scatter operations
!---------------------------------------------------------------------------------
#define ST REAL(r8)
#define SZERO 0.d0
#define SKIND 1
#define SSBUF rsbuf
#define SRBUF rrbuf
#define SMPI OFT_MPI_R8
MODULE oft_scatter_r8
#include "oft_scatter_ops.inc"
END MODULE oft_scatter_r8
#undef ST
#undef SZERO
#undef SKIND
#undef SSBUF
#undef SRBUF
#undef SMPI
#define ST COMPLEX(c8)
#define SZERO (0.d0,0.d0)
#define SKIND 2
#define SSBUF csbuf
#define SRBUF crbuf
#define SMPI OFT_MPI_C8
MODULE oft_scatter_c8
#include "oft_scatter_ops.inc"
END MODULE oft_scatter_c8
#undef ST
#undef SZERO
#undef SKIND
#undef SSBUF
#undef SRBUF
#undef SMPI
!---------------------------------------------------------------------------------
!> Index sets and scatter plans for real and complex values
!!
!! Collects @ref oft_scatter_base and the generic operations of each value type
!! (`scatter_begin`, `scatter_end`, `scatter_wait_recv`, `scatter_unpack`,
!! `scatter_wait_send`, `scatter_local`). See the description of oft_scatter.F90.
!---------------------------------------------------------------------------------
MODULE oft_scatter
USE oft_scatter_base
USE oft_scatter_r8
USE oft_scatter_c8
IMPLICIT NONE
PUBLIC
END MODULE oft_scatter
