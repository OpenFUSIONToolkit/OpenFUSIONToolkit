!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file tracing.F90
!
!> Tracing implementation for the Open FUSION Toolkit
!!
!! The default tracing implementation uses LSODE as the ode solver for advancing
!! field lines.
!!
!! @authors Chris Hansen
!! @date Feburary 2012
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------------
module tracing
use oft_base
use oft_mesh_type, only: mesh_findcell
use multigrid, only: multigrid_mesh
USE oft_io, only: hdf5_write
!---
use fem_utils, only: fem_interp
implicit none
#include "local.h"
!---------------------------------------------------------------------------------
!> Abstract class for OFT tracers
!---------------------------------------------------------------------------------
type, abstract :: oft_tracer
  logical :: initialized = .FALSE. !< Flag to indicate tracer is ready
  integer(i4) :: neq = 3 !< Number of equations in ODE system (3)
  integer(i4) :: maxsteps = 1e6 !< Limit for the number of tracer steps
  integer(i4) :: nsteps = 0 !< Current number of steps
  integer(i4) :: maxtrans = 1e2 !< Limit for the number of domain transfers
  integer(i4) :: ntrans = 0 !< Current number of transfers
  integer(i4) :: cell = 0 !< Cell location of last/current trace point
  integer(i4) :: status = 0 !< Status of tracer
  integer(i4) :: estatus = 0 !< Error status of tracer
  integer(i4) :: proc = 0 !< Processor location of last/current trace point
  integer(i4) :: tids(2) = 0 !< Transfer tags for inter-processor communication
#ifdef OFT_MPI_F08
  type(mpi_request) :: rids(2) = MPI_REQUEST_NULL !< Active MPI recieve request IDs
#else
  integer(i4) :: rids(2) = MPI_REQUEST_NULL !< Active MPI recieve request IDs
#endif
  real(r8) :: tol = 1.d-4 !< Tolerance for ODE solver
  real(r8) :: dt = 0.d0 !< Time step for ODE solver
  real(r8), pointer, dimension(:) :: y !< Possition of current trace point
  real(r8), pointer, dimension(:) :: dy !< Change from previous trace point
  real(r8), pointer, dimension(:) :: dyp !< Previous derivative
  real(r8) :: f(4) = 0.d0 !< Logical cell position of last/current trace point
  type(multigrid_mesh), pointer :: ml_mesh => NULL() !< Muli-grid mesh
  class(fem_interp), pointer :: B => NULL() !< Interpolation operator for trace field
  procedure(tracer_ydot), pointer, nopass :: ydot => NULL() !< General ODE function
contains
  !> Initialize tracer at a given point
  procedure(tracer_setup), deferred :: setup
  !> Create a copy of this tracer
  procedure(tracer_copy), deferred :: copy
  !> Save tracer information for transfer
  procedure(tracer_step), deferred :: save
  !> Transfer tracer information to new domain
  procedure(tracer_send), deferred :: send
  !> Receive tracer information from other domain
  procedure(tracer_recv), deferred :: recv
  !> Advance tracer one step
  procedure(tracer_step), deferred :: step
  !> Clean-up internal storage
  procedure(tracer_step), deferred :: delete
end type oft_tracer
!---------------------------------------------------------------------------------
!> Tracer pointer
!---------------------------------------------------------------------------------
type, public :: oft_tracer_ptr
  class(oft_tracer), pointer :: t => NULL() !< Pointer to tracer object
end type oft_tracer_ptr
!---------------------------------------------------------------------------------
!> Tracer implementation using an Euler method
!---------------------------------------------------------------------------------
type, extends(oft_tracer) :: oft_tracer_euler
  integer(i4) :: isend(2) = 0 !< Integer MPI transfer array
  real(r8), pointer, dimension(:) :: rsend !< Real MPI transfer array
contains
  !> Initialize tracer at a given point
  procedure :: setup => tracer_euler_setup
  !> Create a copy of this tracer
  procedure :: copy => tracer_euler_copy
  !> Save tracer information for transfer
  procedure :: save => tracer_euler_save
  !> Transfer tracer information to new domain
  procedure :: send => tracer_euler_send
  !> Receive tracer information from other domain
  procedure :: recv => tracer_euler_recv
  !> Advance tracer one step
  procedure :: step => tracer_euler_step
  !> Clean-up internal storage
  procedure :: delete => tracer_euler_delete
end type oft_tracer_euler
!---------------------------------------------------------------------------------
!> Tracer implementation using LSODE as the ODE solver package
!---------------------------------------------------------------------------------
type, extends(oft_tracer) :: oft_tracer_lsode
  integer(i4) :: lrw = 0 !< Length of real work array for LSODE
  integer(i4) :: liw = 0 !< Length of integer work array for LSODE
  integer(i4) :: itask = 2 !< ODE advance method
  integer(i4) :: istate = 1 !< LSODE exit state
  integer(i4) :: mf = 10 !< Method flag (10 = Adams, 22 = BDF)
  integer(i4) :: nisend = 39 !< Number of values in isend
  integer(i4) :: nrsend = 222 !< Number of values in rsend
  integer(i4) :: isav(37) = 0 !< Integer common block array for LSODE
  integer(i4), pointer, dimension(:) :: iwork => NULL() !< Integer work array for LSODE
  integer(i4), pointer, dimension(:) :: isend => NULL() !< Integer MPI transfer array
  logical :: comm_load = .FALSE. !< Common block ready to load
  real(r8) :: t = 0.d0 !< Current time in ODE advance
  real(r8) :: tout = 1.d0 !< Desired completion time for ODE advance
  real(r8) :: rsav(218) = 0.d0 !< Real common block array for LSODE
  real(r8), pointer, dimension(:) :: rwork => NULL() !< Real work array for LSODE
  real(r8), pointer, dimension(:) :: rsend => NULL() !< Real MPI transfer array
contains
  !> Initialize tracer at a given point
  procedure :: setup => tracer_lsode_setup
  !> Create a copy of this tracer
  procedure :: copy => tracer_lsode_copy
  !> Save tracer information for transfer
  procedure :: save => tracer_lsode_save
  !> Transfer tracer information to new domain
  procedure :: send => tracer_lsode_send
  !> Receive tracer information from other domain
  procedure :: recv => tracer_lsode_recv
  !> Advance tracer one step
  procedure :: step => tracer_lsode_step
  !> Clean-up internal storage
  procedure :: delete => tracer_lsode_delete
end type oft_tracer_lsode
!---Abstract procedure prototypes
abstract interface
  !---------------------------------------------------------------------------------
  !> Abstract interface for tracer setup with a new start point
  !---------------------------------------------------------------------------------
  subroutine tracer_setup(self,y,cell,init)
  import oft_tracer, i4, r8
  class(oft_tracer), intent(inout) :: self !< Tracer object
  real(r8), intent(in) :: y(:) !< New start point [3]
  integer(i4), optional, intent(in) :: cell !< Guess cell for use in \ref oft_mesh_type::mesh_findcell
  !! "mesh_findcell" (optional)
  logical, optional, intent(in) :: init !< Flag indicating tracer is starting a new trace and all counts
  !! should be set to zero (optional)
  end subroutine tracer_setup
  !------------------------------------------------------------------------------
  !> Abstract interface for creating bare copies of a tracer object
  !------------------------------------------------------------------------------
  subroutine tracer_copy(self,new)
  import oft_tracer
  class(oft_tracer), intent(in) :: self !< Tracer object
  class(oft_tracer), pointer, intent(out) :: new !< Copy of tracer object with same type, tolerances and field
  !! interpolation object
  end subroutine tracer_copy
  !------------------------------------------------------------------------------
  !> Abstract interface for sending tracer information to a new domain
  !------------------------------------------------------------------------------
  function tracer_send(self,proc) result(sid)
  import oft_tracer, i4
#ifdef OFT_MPI_F08
  import mpi_request
#endif
  class(oft_tracer), intent(inout) :: self !< Tracer object
  integer(i4), intent(in) :: proc !< Processor ID to send tracer information to
#ifdef OFT_MPI_F08
  type(mpi_request) :: sid(2)
#else
  integer(i4) :: sid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
  end function tracer_send
  !------------------------------------------------------------------------------
  !> Abstract interface for recieving tracer information from a different domain
  !------------------------------------------------------------------------------
  function tracer_recv(self,proc) result(rid)
  import oft_tracer, i4
#ifdef OFT_MPI_F08
  import mpi_request
#endif
  class(oft_tracer), intent(inout) :: self !< Tracer object
  integer(i4), intent(in) :: proc !< Processor ID to recieve tracer information from
#ifdef OFT_MPI_F08
  type(mpi_request) :: rid(2)
#else
  integer(i4) :: rid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
  end function tracer_recv
  !---------------------------------------------------------------------------------
  !> Abstract interface for advancing a tracer by one step
  !---------------------------------------------------------------------------------
  subroutine tracer_step(self)
  import oft_tracer
  class(oft_tracer), intent(inout) :: self !< Tracer object
  end subroutine tracer_step
  !---------------------------------------------------------------------------------
  !> Abstract interface for general tracer ODE function
  !---------------------------------------------------------------------------------
  subroutine tracer_ydot(t,y,B,n,ydot)
  IMPORT i4, r8
  real(r8), intent(in) :: t
  real(r8), intent(in) :: y(n)
  real(r8), intent(in) :: B(3)
  integer(i4), intent(in) :: n
  real(r8), intent(out) :: ydot(n)
  end subroutine tracer_ydot
end interface
!---Global Variables
class(oft_tracer), pointer :: active_tracer => NULL() !< Active tracer for current thread
!$omp threadprivate(active_tracer)
INTEGER(i4), PRIVATE, PARAMETER :: nlocpts = 1E6 !< Number of points in processor buffers
REAL(r8), PRIVATE, PARAMETER :: offmesh_tol = 1.d-2 !< Tolerance for cell logical coordinate test
LOGICAL, PRIVATE :: test_timeout = .FALSE. !< Timeout has passed?
REAL(r8), PRIVATE :: timeout = 60.d0 !< Timeout for tracing calls (seconds)
TYPE(oft_timer), PRIVATE :: timeout_timer !< Timer for timeout tests
!---Tracer state definitions
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_INIT = 99 !< Tracer is awaiting setup
!---Error states
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_EXIT    = -1 !< Tracer exited mesh
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_TLIMIT  = -2 !< Tracer exceeded `maxsteps`
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_SLIMIT  = -3 !< Tracer exceeded `maxtrans`
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_BLIMIT  = -4 !< Tracer exceeded buffer size
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_TIME    = -5 !< Tracer exceeded timeout
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_ERROR_FAIL    = -6 !< Tracer failed
!---Local advance states
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_TRACE_READY   = 1 !< Tracer is ready to advance
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_TRACE_ACTIVE  = 2 !< Tracer is being advanced
!---Send states
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_SEND_READY   = 10 !< Tracer is ready to send
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_SEND_ACTIVE  = 11 !< Tracer is waiting for send to complete
!---Recieve state
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_RECV_READY   = 20 !< Tracer is ready to recv
INTEGER(i4), PRIVATE, PARAMETER :: TRACER_RECV_ACTIVE  = 21 !< Tracer is waiting for recv to complete
contains
!------------------------------------------------------------------------------
!> Set timeout for tracing calls
!!
!! @param[in] timeout New timeout (seconds)
!------------------------------------------------------------------------------
subroutine set_timeout(new_timeout)
real(r8), intent(in) :: new_timeout
timeout=new_timeout
end subroutine set_timeout
!------------------------------------------------------------------------------
!> Create tracer object
!!
!! The following tracer types are currently available:
!! - (1) \ref tracing::oft_tracer_lsode "LSODE tracer", using LSODE`s non-stiff
!! Adam's method
!! - (2) \ref tracing::oft_tracer_lsode "LSODE tracer", using LSODE`s stiff
!! BDF method
!! - (3) \ref tracing::oft_tracer_euler "EULER tracer", using an adaptive
!! step size
!------------------------------------------------------------------------------
subroutine create_tracer(tracer,type)
class(oft_tracer), pointer, intent(out) :: tracer !< Tracer object
integer(i4), intent(in) :: type !< Tracer type key
DEBUG_STACK_PUSH
SELECT CASE(type)
  CASE(1)
    ALLOCATE(oft_tracer_lsode::tracer)
  CASE(2)
    ALLOCATE(oft_tracer_lsode::tracer)
    SELECT TYPE(this=>tracer)
      CLASS IS(oft_tracer_lsode)
        this%mf=22
    END SELECT
  CASE(3)
    ALLOCATE(oft_tracer_euler::tracer)
  CASE DEFAULT
    CALL oft_abort('Invalid tracer type.','set_tracer',__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine create_tracer
!------------------------------------------------------------------------------
!> Get string describing tracer exit reason
!------------------------------------------------------------------------------
function termination_reason(status_flag) result(reason_str)
INTEGER(i4), INTENT(IN) :: status_flag !< Tracer exit status
CHARACTER(LEN=40) :: reason_str !< Tracer exit reason
SELECT CASE(status_flag)
  CASE(TRACER_ERROR_EXIT)
    reason_str = "Tracer exited mesh"
  CASE(TRACER_ERROR_TLIMIT)
    reason_str = "Tracer exceeded 'maxtrans'"
  CASE(TRACER_ERROR_SLIMIT)
    reason_str = "Tracer exceeded 'maxsteps'"
  CASE(TRACER_ERROR_BLIMIT)
    reason_str = "Tracer exceeded buffer size"
  CASE(TRACER_ERROR_TIME)
    reason_str = "Tracer exceeded timeout"
  CASE(TRACER_ERROR_FAIL)
    reason_str = "Tracer step failed"
  CASE DEFAULT
    reason_str = "Unknown tracer termination reason"
END SELECT
end function termination_reason
!------------------------------------------------------------------------------
!> Advance tracer from a single point and save field line to specified file
!------------------------------------------------------------------------------
subroutine tracing_line(tracer,pt,filename)
class(oft_tracer), TARGET, intent(in) :: tracer !< Tracer to use for field line trace
real(r8), intent(in) :: pt(:) !< Starting point for trace [3]
character(LEN=*), intent(in) :: filename !< Filename to save field line data to
integer(i4) :: i,j,ierr,cell,proc,nsteps,io_unit
#ifdef OFT_MPI_F08
type(mpi_request) :: tid(2)
type(mpi_request), allocatable, dimension(:) :: send_reqr,send_reqi
#else
integer(i4) :: tid(2)
integer(i4), allocatable, dimension(:) :: send_reqr,send_reqi
#endif
real(r8) :: f(4),fmin,fmax,ptmp(3)
real(r8), ALLOCATABLE, DIMENSION(:,:) :: pt_list,tmp_list
DEBUG_STACK_PUSH
ALLOCATE(pt_list(tracer%neq,tracer%maxsteps))
pt_list=0.d0
if(oft_env%head_proc)write(*,'(A)')'Starting fieldline trace'
CALL oft_mpi_barrier(ierr)
test_timeout=.FALSE.
CALL timeout_timer%tick
!---
active_tracer=>tracer
CALL active_tracer%setup(pt,init=.TRUE.)
!---Setup send/recv arrays
allocate(send_reqr(oft_env%nprocs),send_reqi(oft_env%nprocs))
send_reqr=MPI_REQUEST_NULL
send_reqi=MPI_REQUEST_NULL
!---
do while(active_tracer%ntrans<active_tracer%maxtrans)
    cell=0
    IF(active_tracer%ml_mesh%nbase<active_tracer%ml_mesh%mgdim)THEN
      call mesh_findcell(active_tracer%ml_mesh%meshes(active_tracer%ml_mesh%nbase),cell,active_tracer%y(1:3),f)
    ELSE
      call mesh_findcell(active_tracer%ml_mesh%mesh,cell,active_tracer%y(1:3),f)
    END IF
    fmin=MINVAL(f); fmax=MAXVAL(f)
    !---Disable point if off mesh
    if(( fmax>1.d0+offmesh_tol ).OR.( fmin<-offmesh_tol ))active_tracer%status=TRACER_ERROR_EXIT
    !---Disable point if over maxtrans
    if(active_tracer%ntrans>=active_tracer%maxtrans)active_tracer%status=TRACER_ERROR_TLIMIT
    !---Disable point if over maxsteps
    if(active_tracer%nsteps>=active_tracer%maxsteps)active_tracer%status=TRACER_ERROR_SLIMIT
    !---Check error status
    IF(active_tracer%estatus/=0)active_tracer%status=active_tracer%estatus
    !---Exit trace for fail status or timeout
    IF(active_tracer%status<0)EXIT
    IF(timeout_timer%timeout(timeout))EXIT
    !---Get new processor
    IF(active_tracer%ml_mesh%nbase<active_tracer%ml_mesh%mgdim)THEN
      active_tracer%proc=active_tracer%ml_mesh%meshes(active_tracer%ml_mesh%nbase+1)%base%lcpart(cell)-1
      active_tracer%cell=(cell-1)*8+1
    ELSE
      active_tracer%proc=oft_env%rank
      active_tracer%cell=cell
    END IF
    !---
    if(active_tracer%proc==oft_env%rank)then
      active_tracer%status=TRACER_TRACE_ACTIVE
      call active_tracer%setup(active_tracer%y)
      do while(active_tracer%nsteps<active_tracer%maxsteps)
        call active_tracer%step
        pt_list(:,active_tracer%nsteps)=active_tracer%y
        if((active_tracer%cell==0).OR.(active_tracer%estatus/=0))exit
      end do
      IF(active_tracer%nsteps>=active_tracer%maxsteps)active_tracer%estatus=TRACER_ERROR_SLIMIT
      CALL active_tracer%save
      !---
      IF(oft_env%nprocs>1)THEN
#ifdef HAVE_MPI
        !---Send point
        DO j=0,oft_env%nprocs-1
          IF(j==oft_env%rank)CYCLE
          tid=active_tracer%send(j)
          send_reqr(j+1)=tid(1)
          send_reqi(j+1)=tid(2)
        END DO
        !---Test for completion of send
        CALL MPI_WAITALL(oft_env%nprocs,send_reqr,MPI_STATUSES_IGNORE,ierr)
        CALL MPI_WAITALL(oft_env%nprocs,send_reqi,MPI_STATUSES_IGNORE,ierr)
        active_tracer%ntrans=active_tracer%ntrans+1
#else
        CALL oft_abort("Distributed mesh requires MPI","tracing_line",__FILE__)
#endif
      END IF
    else
      DO
        tid=active_tracer%recv(active_tracer%proc)
        IF(oft_mpi_check_reqs(2,tid))EXIT
      END DO
      active_tracer%ntrans=active_tracer%ntrans+1
    end if
    ! if(oft_env%head_proc.AND.MOD(active_tracer%ntrans,100)==0)write(*,*)active_tracer%ntrans,active_tracer%nsteps
end do
IF(oft_env%head_proc)THEN
  WRITE(*,'(2A)')  '  Tracing complete: ',termination_reason(active_tracer%status)
  WRITE(*,'(A,I8)')'    # of steps     = ',active_tracer%nsteps
  WRITE(*,'(A,I8)')'    # of transfers = ',active_tracer%ntrans
  WRITE(*,*)
  ALLOCATE(tmp_list(tracer%neq,tracer%nsteps))
#ifdef HAVE_MPI
  tmp_list=0.d0
  CALL MPI_Reduce(pt_list,tmp_list,tracer%neq*tracer%nsteps, &
    MPI_REAL8,MPI_SUM,0,oft_env%COMM,ierr)
#else
  tmp_list=pt_list(:,1:tracer%nsteps)
#endif
  IF(oft_env%rank==0)CALL hdf5_write(tmp_list,filename,'trace')
  DEALLOCATE(tmp_list)
ELSE
#ifdef HAVE_MPI
  CALL MPI_Reduce(pt_list,pt_list,tracer%neq*tracer%nsteps, &
    MPI_REAL8,MPI_SUM,0,oft_env%COMM,ierr)
#endif
END IF
deallocate(send_reqr,send_reqi,pt_list)
DEBUG_STACK_POP
end subroutine tracing_line
!------------------------------------------------------------------------------
!> Advance a group of tracers while accumulating crossing of the yz-plane into
!! a poincare section.
!!
!! Tracing parallelizes across domains and streamlines. Each domain advances
!! fieldlines in parallel using local worker threads. Communication between domains
!! is handled by a master thread.
!!
!! @note When the local buffer for point crossings is exceeded tracing will
!! exit. The size of the point buffer is set by the private module variable
!! @ref tracing::nlocpts "nlocpts"
!------------------------------------------------------------------------------
subroutine tracing_poincare(tracer,pts,n,filename,offset,pcoord,qfile)
class(oft_tracer), intent(in) :: tracer !< Base tracer for field line advance
integer(i4), intent(in) :: n !< Number of launch points
real(r8), intent(in) :: pts(3,n) !< Launch points for field lines [3,n]
character(LEN=*), intent(in) :: filename !< Filename for section data
real(r8), optional, intent(in) :: offset !< X-intercept of section plane (optional, default: 0.0)
integer(i4), optional, intent(in) :: pcoord !< Coordinate index orthogonal to section plane (optional, default: 1)
character(LEN=*), optional, intent(in) :: qfile !< Filename for approximate safety factor, must be z-axis oriented (optional)
integer(i4) :: i,j,ierr,nloc,nthread,io_unit,counts(6),icoord
integer(i4), allocatable, dimension(:) :: ibuff
real(r8) :: xcross,ntime,elapsed_time
real(r8), allocatable, dimension(:,:) :: outbuff,qbuff
type(oft_tracer_ptr), allocatable, dimension(:) :: tracers
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  CALL mytimer%tick
  WRITE(*,*)'Starting Poincare Section'
END IF
IF(tracer%neq/=3)CALL oft_abort("Poincare plot requires tracer with neq==3", &
  "tracing_poincare",__FILE__)
ALLOCATE(ibuff(nlocpts),outbuff(3,nlocpts))
ntime=10.d0
xcross=0.d0
icoord=1
IF(PRESENT(offset))xcross=offset
IF(PRESENT(pcoord))icoord=pcoord
IF(PRESENT(qfile))THEN
  ALLOCATE(qbuff(2,n))
  qbuff=0.d0
END IF
IF(oft_env%rank==0)THEN
  OPEN(NEWUNIT=io_unit,FILE=TRIM(filename),STATUS="REPLACE")
  CLOSE(io_unit)
END IF
ALLOCATE(tracers(n))
DO i=1,n
  CALL tracer%copy(tracers(i)%t)
  CALL tracers(i)%t%setup(pts(:,i),init=.TRUE.)
  tracers(i)%t%tids=(/2*(i-1)+1,2*i/)
END DO
CALL oft_mpi_barrier(ierr)
!$omp parallel
!$omp single
nthread=omp_get_num_threads()
!$omp end single
!$omp barrier
!$omp end parallel
IF(nthread==1)CALL omp_set_num_threads(2)
ibuff=-1
test_timeout=.FALSE.
CALL timeout_timer%tick
nloc=0
!$omp parallel
IF(.NOT.omp_in_parallel())CALL oft_abort('Poincare tracing requires >1 OpenMP threads.', &
  'tracing_poincare',__FILE__)
IF(omp_get_thread_num()==0)THEN
  !---Setup send/recv arrays
  CALL tracing_poincare_master(tracers,n)
ELSE
  !---Start loop
  IF(PRESENT(qfile))THEN
    CALL tracing_poincare_worker(tracers,n,ibuff,outbuff,nlocpts,nloc,xcross,icoord,qbuff)
  ELSE
    CALL tracing_poincare_worker(tracers,n,ibuff,outbuff,nlocpts,nloc,xcross,icoord)
  END IF
END IF
!$omp barrier
!$omp end parallel
IF(nloc>=nlocpts)WRITE(*,*)'Warning: point buffer limit reached on task ',oft_env%rank
IF(nthread==1)CALL omp_set_num_threads(1)
!---Save punctures
DO i=0,oft_env%nprocs-1
  CALL oft_mpi_barrier(ierr)
  IF(i/=oft_env%rank)CYCLE
  OPEN(NEWUNIT=io_unit,FILE=TRIM(filename),POSITION="APPEND",STATUS="OLD")
  DO j=1,nlocpts
    IF(ibuff(j)==-1)EXIT
    WRITE(io_unit,'(I6,3E14.6)')ibuff(j),outbuff(:,j)
  END DO
  CLOSE(io_unit)
END DO
DEALLOCATE(ibuff,outbuff)
!---Save safety factor
IF(PRESENT(qfile))THEN
  ALLOCATE(outbuff(2,n))
#ifdef HAVE_MPI
  outbuff=0.d0
  CALL MPI_Reduce(qbuff,outbuff,2*n,MPI_REAL8,MPI_SUM,0,oft_env%COMM,ierr)
#else
  outbuff=qbuff
#endif
  IF(oft_env%head_proc)THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(qfile))
    DO j=1,n
      WRITE(io_unit,'(5E14.6)')pts(:,j),outbuff(:,j)
    END DO
    CLOSE(io_unit)
  END IF
  DEALLOCATE(qbuff,outbuff)
END IF
!---Write out tracer stats
IF(oft_env%head_proc)THEN
  counts=0
  DO i=1,n
    SELECT CASE(tracers(i)%t%status)
      CASE(TRACER_ERROR_EXIT)
        counts(1)=counts(1)+1
      CASE(TRACER_ERROR_TLIMIT)
        counts(2)=counts(2)+1
      CASE(TRACER_ERROR_SLIMIT)
        counts(3)=counts(3)+1
      CASE(TRACER_ERROR_TIME)
        counts(4)=counts(4)+1
      CASE(TRACER_ERROR_BLIMIT)
        counts(5)=counts(5)+1
      CASE(TRACER_ERROR_FAIL)
        counts(6)=counts(6)+1
    END SELECT
  END DO
  elapsed_time=mytimer%tock()
  WRITE(*,'(1X,A,F6.1)')'Trace Complete: ',elapsed_time
  WRITE(*,*)'  # exited mesh              ',counts(1)
  WRITE(*,*)'  # exceeded transfer limit  ',counts(2)
  WRITE(*,*)'  # exceeded step limit      ',counts(3)
  WRITE(*,*)'  # exceeded timeout         ',counts(4)
  WRITE(*,*)'  # stopped by buffer size   ',counts(5)
  WRITE(*,*)'  # failure in tracer        ',counts(6)
END IF
!---Clean-up storage
DO i=1,n
  CALL tracers(i)%t%delete
  DEALLOCATE(tracers(i)%t)
END DO
DEALLOCATE(tracers)
DEBUG_STACK_POP
end subroutine tracing_poincare
!------------------------------------------------------------------------------
!> Manage tracer communication and status for parallel tracing
!------------------------------------------------------------------------------
subroutine tracing_poincare_master(tracers,n)
integer(i4), intent(in) :: n !< Number of active tracers
type(oft_tracer_ptr), target, intent(inout) :: tracers(n) !< Array of all active tracers [n]
integer(i4) :: i,j,ierr,cell
#ifdef OFT_MPI_F08
type(mpi_request) :: tid(2)
type(mpi_request), allocatable, dimension(:) :: recv_reqr,recv_reqi
type(mpi_request), allocatable, dimension(:,:) :: send_reqr,send_reqi
#else
integer(i4) :: tid(2)
integer(i4), allocatable, dimension(:) :: recv_reqr,recv_reqi
integer(i4), allocatable, dimension(:,:) :: send_reqr,send_reqi
#endif
logical :: testr,testi,testc
real(r8) :: f(4),fmin,fmax,tol_ratio
class(oft_tracer), pointer :: loc_tracer
type(multigrid_mesh), pointer :: mg_mesh
DEBUG_STACK_PUSH
mg_mesh=>tracers(1)%t%ml_mesh
!---Setup send/recv arrays
allocate(send_reqr(oft_env%nprocs,n),send_reqi(oft_env%nprocs,n))
allocate(recv_reqr(n),recv_reqi(n))
send_reqr=MPI_REQUEST_NULL
send_reqi=MPI_REQUEST_NULL
recv_reqr=MPI_REQUEST_NULL
recv_reqi=MPI_REQUEST_NULL
!---Start loop
DO
  testc=.TRUE.
  IF(.NOT.test_timeout)test_timeout=timeout_timer%timeout(timeout)
  do i=1,n
    loc_tracer=>tracers(i)%t
    testc=testc.AND.(loc_tracer%status<0)
    if(loc_tracer%status<10)cycle !---Point is tracing or disabled
    if(loc_tracer%status==TRACER_SEND_READY)then !---Point is ready for send
      IF(test_timeout)THEN
        loc_tracer%status=TRACER_ERROR_TIME
        CYCLE
      END IF
      !---Check if phony transfer
      IF(loc_tracer%estatus==0)THEN
        cell=0
        IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
          call mesh_findcell(mg_mesh%meshes(mg_mesh%nbase),cell,loc_tracer%y,f)
          tol_ratio=mg_mesh%mesh%hrms/mg_mesh%meshes(mg_mesh%nbase)%hrms*.1d0
        ELSE
          call mesh_findcell(mg_mesh%mesh,cell,loc_tracer%y,f)
          tol_ratio=.1d0
        END IF
        fmin=minval(f)
        fmax=maxval(f)
        !---Disable point if off mesh
        IF(( fmax<1.d0+offmesh_tol*tol_ratio ).AND.( fmin>-offmesh_tol*tol_ratio ))THEN
          !---Find new proc and setup recv
          IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
            IF(mg_mesh%meshes(mg_mesh%nbase+1)%base%lcpart(cell)-1==oft_env%rank)THEN
              loc_tracer%status=TRACER_TRACE_READY
              CYCLE
            END IF
          ELSE
            loc_tracer%status=TRACER_TRACE_READY
            CYCLE
          END IF
        END IF
      END IF
      !write(*,*)'Sending Point ',oft_env%rank,i,loc_tracer%y
      IF(oft_env%nprocs==1)THEN
        loc_tracer%status=TRACER_SEND_ACTIVE
      ELSE
        !---Send to all other procs
        do j=0,oft_env%nprocs-1
          if(j==oft_env%rank)cycle
            tid=loc_tracer%send(j)
            send_reqr(j+1,i)=tid(1)
            send_reqi(j+1,i)=tid(2)
        end do
      END IF
    else if(loc_tracer%status==TRACER_SEND_ACTIVE)then !---Point is waiting for send to complete
      IF(test_timeout)THEN
#ifdef HAVE_MPI
        DO j=1,oft_env%nprocs
          IF(send_reqr(j,i)/=MPI_REQUEST_NULL)CALL MPI_CANCEL(send_reqr(j,i),ierr)
          IF(send_reqi(j,i)/=MPI_REQUEST_NULL)CALL MPI_CANCEL(send_reqi(j,i),ierr)
        END DO
#endif
        loc_tracer%status=TRACER_ERROR_TIME
        CYCLE
      END IF
      !---Test for completion of send
#ifdef HAVE_MPI
      call MPI_TESTALL(oft_env%nprocs,send_reqr(:,i),testr,MPI_STATUSES_IGNORE,ierr)
      call MPI_TESTALL(oft_env%nprocs,send_reqi(:,i),testi,MPI_STATUSES_IGNORE,ierr)
#else
      testr=.TRUE.
      testi=.TRUE.
#endif
      !---Send is complete
      if(testr.AND.testi)then
        loc_tracer%ntrans=loc_tracer%ntrans+1
        !---Disable point if over maxtrans
        if(loc_tracer%ntrans>=loc_tracer%maxtrans)loc_tracer%status=TRACER_ERROR_TLIMIT
        !---Disable point if over maxsteps
        if(loc_tracer%nsteps>=loc_tracer%maxsteps)loc_tracer%status=TRACER_ERROR_SLIMIT
        !---Check error status
        IF(loc_tracer%estatus/=0)loc_tracer%status=loc_tracer%estatus
        if(loc_tracer%status<0)cycle
        cell=0
        IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
          call mesh_findcell(mg_mesh%meshes(mg_mesh%nbase),cell,loc_tracer%y,f)
          tol_ratio=mg_mesh%mesh%hrms/mg_mesh%meshes(mg_mesh%nbase)%hrms*.1d0
        ELSE
          call mesh_findcell(mg_mesh%mesh,cell,loc_tracer%y,f)
          tol_ratio=.1d0
        END IF
        fmin=minval(f)
        fmax=maxval(f)
        !---Disable point if off mesh
        if(( fmax>1.d0+offmesh_tol*tol_ratio ).OR.( fmin<-offmesh_tol*tol_ratio ))loc_tracer%status=TRACER_ERROR_EXIT
        if(loc_tracer%status<0)cycle
        !---Find new proc and setup recv
        IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
          loc_tracer%proc=mg_mesh%meshes(mg_mesh%nbase+1)%base%lcpart(cell)-1
          loc_tracer%cell=(cell-1)*8+1
        ELSE
          loc_tracer%proc=oft_env%rank
          loc_tracer%cell=cell
        END IF
        !---
        recv_reqr(i)=MPI_REQUEST_NULL
        recv_reqi(i)=MPI_REQUEST_NULL
        if(loc_tracer%proc==oft_env%rank)then
          loc_tracer%status=TRACER_TRACE_READY
        else ! Setup recv
          tid=loc_tracer%recv(loc_tracer%proc)
          recv_reqr(i)=tid(1)
          recv_reqi(i)=tid(2)
        end if
      end if
    else if(loc_tracer%status==TRACER_RECV_ACTIVE)then !---Point is waiting for recv to complete
      IF(test_timeout)THEN
#ifdef HAVE_MPI
        IF(recv_reqr(i)/=MPI_REQUEST_NULL)CALL MPI_CANCEL(recv_reqr(i),ierr)
        IF(recv_reqi(i)/=MPI_REQUEST_NULL)CALL MPI_CANCEL(recv_reqi(i),ierr)
#endif
        loc_tracer%status=TRACER_ERROR_TIME
        CYCLE
      END IF
      !---Test for completion of recv
      tid=loc_tracer%recv(loc_tracer%proc)
      IF(oft_mpi_check_reqs(2,tid))THEN
        !---Recv is complete
        loc_tracer%ntrans=loc_tracer%ntrans+1
        !---Disable point if over maxtrans
        if(loc_tracer%ntrans>=loc_tracer%maxtrans)loc_tracer%status=TRACER_ERROR_TLIMIT
        !---Disable point if over maxsteps
        if(loc_tracer%nsteps>=loc_tracer%maxsteps)loc_tracer%status=TRACER_ERROR_SLIMIT
        !---Check error status
        IF(loc_tracer%estatus/=0)loc_tracer%status=loc_tracer%estatus
        if(loc_tracer%status<0)CYCLE
        !---
        cell=0
        IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
          call mesh_findcell(mg_mesh%meshes(mg_mesh%nbase),cell,loc_tracer%y,f)
          tol_ratio=mg_mesh%mesh%hrms/mg_mesh%meshes(mg_mesh%nbase)%hrms*.1d0
        ELSE
          call mesh_findcell(mg_mesh%mesh,cell,loc_tracer%y,f)
          tol_ratio=.1d0
        END IF
        fmin=minval(f)
        fmax=maxval(f)
        !---Disable point if off mesh
        if(( fmax>1.d0+offmesh_tol*tol_ratio ).OR.( fmin<-offmesh_tol*tol_ratio ))loc_tracer%status=TRACER_ERROR_EXIT
        if(loc_tracer%status<0)CYCLE
        IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
          loc_tracer%proc=mg_mesh%meshes(mg_mesh%nbase+1)%base%lcpart(cell)-1
          loc_tracer%cell=(cell-1)*8+1
        ELSE
          loc_tracer%proc=oft_env%rank
          loc_tracer%cell=cell
        END IF
        recv_reqr(i)=MPI_REQUEST_NULL
        recv_reqi(i)=MPI_REQUEST_NULL
        if(loc_tracer%proc==oft_env%rank)then
          loc_tracer%status=TRACER_TRACE_READY
        else ! Setup recv
          tid=loc_tracer%recv(loc_tracer%proc)
          recv_reqr(i)=tid(1)
          recv_reqi(i)=tid(2)
        end if
      end if
    else if(loc_tracer%status==TRACER_INIT)then !---Point needs to be initialized
      IF(test_timeout)THEN
        loc_tracer%status=TRACER_ERROR_TIME
        CYCLE
      END IF
      !write(*,*)'Initializing Point ',oft_env%rank,i
      !---
      cell=0
      IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
        call mesh_findcell(mg_mesh%meshes(mg_mesh%nbase),cell,loc_tracer%y,f)
        tol_ratio=mg_mesh%mesh%hrms/mg_mesh%meshes(mg_mesh%nbase)%hrms*.1d0
      ELSE
        call mesh_findcell(mg_mesh%mesh,cell,loc_tracer%y,f)
        tol_ratio=.1d0
      END IF
      fmin=minval(f)
      fmax=maxval(f)
      !---Disable point if off mesh
      if(( fmax>=1.d0+offmesh_tol ).OR.( fmin<=-offmesh_tol ))loc_tracer%status=TRACER_ERROR_EXIT
      !---Disable point if over maxtrans
      if(loc_tracer%ntrans>=loc_tracer%maxtrans)loc_tracer%status=TRACER_ERROR_TLIMIT
      if(loc_tracer%status<0)exit
      IF(mg_mesh%nbase<mg_mesh%mgdim)THEN
        loc_tracer%proc=mg_mesh%meshes(mg_mesh%nbase+1)%base%lcpart(cell)-1
        loc_tracer%cell=(cell-1)*8+1
      ELSE
        loc_tracer%proc=oft_env%rank
        loc_tracer%cell=cell
      END IF
      recv_reqr(i)=MPI_REQUEST_NULL
      recv_reqi(i)=MPI_REQUEST_NULL
      if(loc_tracer%proc==oft_env%rank)then
        loc_tracer%status=TRACER_TRACE_READY
      else ! Setup recv
        tid=loc_tracer%recv(loc_tracer%proc)
        recv_reqr(i)=tid(1)
        recv_reqi(i)=tid(2)
      end if
    end if
  end do
  if(testc)exit !---Tracing complete
end do
!---Cleanup
deallocate(send_reqr,send_reqi)
deallocate(recv_reqr,recv_reqi)
DEBUG_STACK_POP
end subroutine tracing_poincare_master
!------------------------------------------------------------------------------
!> Advance a group of tracers on the local domain while accumulating crossings
!! of the yz-plane into a buffer array
!------------------------------------------------------------------------------
subroutine tracing_poincare_worker(tracers,n,ibuff,outbuff,nbuff,nloc,xcross,icoord,qbuff)
integer(i4), intent(in) :: n !< Number of active tracers
type(oft_tracer_ptr), target, intent(inout) :: tracers(n) !< Array of all active tracers [n]
integer(i4), intent(in) :: nbuff !< Size of local buffer
integer(i4), intent(inout) :: ibuff(nbuff) !< Local buffer of tracer indices [nbuff]
real(r8), intent(inout) :: outbuff(3,nbuff) !< Local buffer of plane crossings [3,nbuff]
integer(i4), intent(inout) :: nloc !< Number of crossings on local task
real(r8), optional, intent(in) :: xcross !< X-intercept of section plane (optional, default: 0.0)
integer(i4), optional, intent(in) :: icoord !< Coordinate index orthogonal to section plane (optional, default: 1)
real(r8), optional, intent(inout) :: qbuff(2,n) !< Buffer for q information
integer(i4) :: i,j,cell
logical :: testc
real(r8) :: yp(3),y_plane,yp_plane,dtheta,dphi
class(oft_tracer), pointer :: loc_tracer
DEBUG_STACK_PUSH
!---Start loop
do
  testc=.TRUE.
  do i=1,n
    j=0
    !$omp critical
    active_tracer=>tracers(i)%t
    loc_tracer=>tracers(i)%t
    testc=testc.AND.(loc_tracer%status<0)
    if(loc_tracer%status==TRACER_TRACE_READY)then
      IF(test_timeout)THEN
        loc_tracer%status=TRACER_ERROR_TIME
      ELSE
        loc_tracer%status=TRACER_TRACE_ACTIVE
        IF(nloc<nlocpts)THEN
          j=1
        ELSE
          loc_tracer%estatus=TRACER_ERROR_BLIMIT
          loc_tracer%status=TRACER_SEND_READY
        END IF
      END IF
    end if
    !$omp end critical
    if(j==0)cycle
    cell=loc_tracer%cell
    call loc_tracer%setup(loc_tracer%y)
    yp=loc_tracer%y
    yp_plane=loc_tracer%y(icoord)-xcross
    do while(loc_tracer%nsteps<loc_tracer%maxsteps.AND.loc_tracer%nsteps>=0)
      call loc_tracer%step
      y_plane=loc_tracer%y(icoord)-xcross
      if(yp_plane<=0.d0.AND.y_plane>0.d0)then
        !$omp critical
        IF(nloc<nlocpts)THEN
          nloc=nloc+1
          ibuff(nloc)=i
          outbuff(:,nloc)=loc_tracer%y &
            - (loc_tracer%y-yp)*(loc_tracer%y(icoord)-xcross)/(loc_tracer%y(icoord)-yp(icoord))
        END IF
        !$omp end critical
      end if
      if(PRESENT(qbuff))then
        IF(ALL(ABS(loc_tracer%dyp(1:3))>1.d-8))THEN
          CALL approx_angle_change(loc_tracer%dyp(1:3),loc_tracer%dy(1:3),yp,loc_tracer%y,dtheta,dphi)
          qbuff(:,i) = qbuff(:,i) + (/dtheta, dphi/)
        END IF
      end if
      if((loc_tracer%cell==0).OR.(loc_tracer%estatus/=0))exit
      if(test_timeout)exit
      yp=loc_tracer%y
      yp_plane=y_plane
    end do
    IF(test_timeout)THEN
      loc_tracer%status=TRACER_ERROR_TIME
    ELSE
      IF(loc_tracer%nsteps>=loc_tracer%maxsteps)loc_tracer%estatus=TRACER_ERROR_SLIMIT
      CALL loc_tracer%save
      loc_tracer%status=TRACER_SEND_READY
    END IF
  end do
  if(testc)exit !---Tracing complete
end do
DEBUG_STACK_POP
contains
!---Compute approzximate angle changes for step
subroutine approx_angle_change(b1,b2,r1,r2,dpol,dtor)
real(r8), intent(in) :: b1(3),b2(3),r1(3),r2(3)
real(r8), intent(out) :: dpol,dtor
real(r8) :: phi1,phi2,theta1,theta2
!
phi1 = ATAN2(r1(2), r1(1))
phi2 = ATAN2(r2(2), r2(1))
!
theta1 = ATAN2(b1(3), b1(1)*COS(phi1) + b1(2)*SIN(phi1))
theta2 = ATAN2(b2(3), b2(1)*COS(phi2) + b2(2)*SIN(phi2))
dpol = theta2 - theta1
IF(ABS(dpol)>pi)dpol = dpol - 2*pi*SIGN(1.d0,dpol)
!
theta1 = ATAN2(b1(2), b1(1))
theta2 = ATAN2(b2(2), b2(1))
dtor = theta2 - theta1
IF(ABS(dtor)>pi)dtor = dtor - 2*pi*SIGN(1.d0,dtor)
end subroutine approx_angle_change
end subroutine tracing_poincare_worker
!------------------------------------------------------------------------------
!> Setup a EULER tracer with a new starting point
!------------------------------------------------------------------------------
subroutine tracer_euler_setup(self,y,cell,init)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
real(r8), intent(in) :: y(:) !< New start point [3]
integer(i4), optional, intent(in) :: cell !< Guess cell for use in \ref oft_mesh_type::mesh_findcell
!! "mesh_findcell" (optional)
logical, optional, intent(in) :: init !< Flag indicating tracer is starting a new trace and all counts
!! should be set to zero (optional)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%y))THEN
  ALLOCATE(self%y(self%neq))
  ALLOCATE(self%dy(self%neq))
  ALLOCATE(self%dyp(self%neq))
  ALLOCATE(self%rsend(2*self%neq+1))
  self%dyp=0.d0
END IF
IF(PRESENT(init))THEN
  IF(init)THEN
    self%ntrans=0
    self%nsteps=0
    self%status=TRACER_INIT
    self%estatus=0
  END IF
END IF
!---Find starting point
self%cell=0
if(present(cell))self%cell=cell
self%y=y
call mesh_findcell(self%ml_mesh%mesh,self%cell,self%y(1:3),self%f)
self%initialized=.TRUE.
DEBUG_STACK_POP
end subroutine tracer_euler_setup
!------------------------------------------------------------------------------
!> Create a bare copy of the LSODE tracer object
!!
!! @param[out] new Copy of object with same tolerances and field interpolation
!! object
!------------------------------------------------------------------------------
subroutine tracer_euler_copy(self,new)
class(oft_tracer_euler), intent(in) :: self !< Tracer object
class(oft_tracer), pointer, intent(out) :: new !< Copy of object with same tolerances and referenced interpolator
DEBUG_STACK_PUSH
allocate(oft_tracer_euler::new)
new%B=>self%B
new%maxsteps=self%maxsteps
new%tol=self%tol
new%maxtrans=self%maxtrans
new%neq=self%neq
IF(ASSOCIATED(self%ydot))new%ydot=>self%ydot
DEBUG_STACK_POP
end subroutine tracer_euler_copy
!------------------------------------------------------------------------------
!> Advancing tracer by one step using LSODE
!------------------------------------------------------------------------------
subroutine tracer_euler_save(self)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
end subroutine tracer_euler_save
!------------------------------------------------------------------------------
!> Send tracer information to a new domain
!!
!! This subroutine transfers the number of steps and tracer position. Data
!! is packed into a single array for each real and integer data types and
!! non-blocking MPI sends are created to transfer each array
!------------------------------------------------------------------------------
function tracer_euler_send(self,proc) result(sid)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
integer(i4), intent(in) :: proc !< Processor ID to send tracer information to
#ifdef OFT_MPI_F08
type(mpi_request) :: sid(2)
#else
integer(i4) :: sid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
#ifdef HAVE_MPI
integer(i4) :: ierr
DEBUG_STACK_PUSH
!---Pack info
self%isend(1)=self%nsteps
self%isend(2)=self%estatus
self%rsend(1:self%neq)=self%y
self%rsend(self%neq+1:2*self%neq)=self%dy
self%rsend(2*self%neq+1)=self%dt
!---Set send flags
sid=MPI_REQUEST_NULL
!---Send to all other procs
CALL MPI_ISEND(self%rsend,self%neq+1,OFT_MPI_R8,proc,self%tids(1),oft_env%COMM,sid(1),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','tracer_euler_send',__FILE__)
CALL MPI_ISEND(self%isend,2,OFT_MPI_I4,proc,self%tids(2),oft_env%COMM,sid(2),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','tracer_euler_send',__FILE__)
self%status=TRACER_SEND_ACTIVE
DEBUG_STACK_POP
#else
sid=MPI_REQUEST_NULL
CALL oft_abort("Send requested without MPI","tracer_euler_send",__FILE__)
#endif
end function tracer_euler_send
!------------------------------------------------------------------------------
!> Recieve tracer information from a different domain
!!
!! This subroutine sets up non-blocking MPI recieve requests and unpacks the
!! resulting data sent by \ref tracing::tracer_euler_send "tracer_euler_send"
!------------------------------------------------------------------------------
function tracer_euler_recv(self,proc) result(rid)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
integer(i4), intent(in) :: proc !< Processor ID to receive tracer information from
#ifdef OFT_MPI_F08
type(mpi_request) :: rid(2)
#else
integer(i4) :: rid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
#ifdef HAVE_MPI
integer(i4) :: ierr
logical :: testr,testi
DEBUG_STACK_PUSH
IF(oft_mpi_check_reqs(2,self%rids))THEN
  !---Setup recv
  CALL MPI_IRECV(self%rsend,self%neq+1,OFT_MPI_R8,proc,self%tids(1),oft_env%COMM,self%rids(1),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','tracer_euler_recv',__FILE__)
  CALL MPI_IRECV(self%isend,2,OFT_MPI_I4,proc,self%tids(2),oft_env%COMM,self%rids(2),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','tracer_euler_recv',__FILE__)
  self%status=TRACER_RECV_ACTIVE
  rid=self%rids
ELSE
  !---Test for completion of recv
  call MPI_TEST(self%rids(1),testr,MPI_STATUS_IGNORE,ierr)
  call MPI_TEST(self%rids(2),testi,MPI_STATUS_IGNORE,ierr)
  rid=self%rids
  !---Set flag and restore if complete
  IF(testr.AND.testi)THEN
    self%nsteps=self%isend(1)
    self%estatus=self%isend(2)
    self%y=self%rsend(1:self%neq)
    self%dy=self%rsend(self%neq+1:2*self%neq)
    self%dt=self%rsend(2*self%neq+1)
    rid=MPI_REQUEST_NULL
    self%rids=MPI_REQUEST_NULL
  END IF
END IF
DEBUG_STACK_POP
#else
rid=MPI_REQUEST_NULL
CALL oft_abort("Recieve requested without MPI","tracer_euler_recv",__FILE__)
#endif
end function tracer_euler_recv
!------------------------------------------------------------------------------
!> Advancing tracer by one step using EULER
!------------------------------------------------------------------------------
subroutine tracer_euler_step(self)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
real(r8), ALLOCATABLE, DIMENSION(:) :: ydot
real(r8) :: B(3),goptmp(3,4),v,fmin,fmax
DEBUG_STACK_PUSH
ALLOCATE(ydot(self%neq))
!---Prediction
call mesh_findcell(self%ml_mesh%mesh,self%cell,self%y(1:3),self%f)
fmin=MINVAL(self%f); fmax=MAXVAL(self%f)
IF(( fmax>=1.d0+offmesh_tol ).OR.( fmin<=-offmesh_tol ))self%cell=0
IF(self%cell==0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
CALL self%ml_mesh%mesh%jacobian(self%cell,self%f,goptmp,v)
CALL self%B%interp(self%cell,self%f,goptmp,B)
self%dyp=self%dy
IF(ASSOCIATED(active_tracer%ydot))THEN
  CALL self%ydot(0.d0,self%y,B,self%neq,self%dy)
ELSE
  self%dy=B
END IF
IF(self%nsteps==0)self%dt=1.d-3/SQRT(SUM(self%dy**2))
!---Step
self%y=self%y+self%dt*self%dy
!---Compute second derivative for error approximation
call mesh_findcell(self%ml_mesh%mesh,self%cell,self%y(1:3),self%f)
fmin=MINVAL(self%f); fmax=MAXVAL(self%f)
IF(( fmax>=1.d0+offmesh_tol ).OR.( fmin<=-offmesh_tol ))self%cell=0
IF(self%cell==0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
CALL self%ml_mesh%mesh%jacobian(self%cell,self%f,goptmp,v)
CALL self%B%interp(self%cell,self%f,goptmp,B)
IF(ASSOCIATED(active_tracer%ydot))THEN
  CALL self%ydot(0.d0,self%y,B,self%neq,ydot)
ELSE
  ydot=B
END IF
!---Compute error and update time step
ydot=self%dt*(ydot - self%dy)/2.d0 ! LTE = h^2 * y'' / 2 = h^2 * (yd_1 - yd_0)/(2*h)
self%dt=MIN(self%tol/SQRT(SUM(ydot**2)),2.d0*self%dt)
!---
self%nsteps=self%nsteps+1
IF(ANY(self%y(1:3)>1.d20))THEN
  self%cell=0
  self%estatus=TRACER_ERROR_FAIL
END IF
DEALLOCATE(ydot)
DEBUG_STACK_POP
end subroutine tracer_euler_step
!------------------------------------------------------------------------------
!> Delete internal storage and reset counters
!------------------------------------------------------------------------------
subroutine tracer_euler_delete(self)
class(oft_tracer_euler), intent(inout) :: self !< Tracer object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%y))DEALLOCATE(self%y,self%dy,self%rsend)
self%nsteps=0
self%ntrans=0
self%dt=0.d0
self%initialized=.FALSE.
NULLIFY(self%ydot)
DEBUG_STACK_POP
end subroutine tracer_euler_delete
!---------------------------------------------------------------------------------
!> Cast a tracer object to a oft_tracer_lsode
!!
!! The source matrix must be @ref oft_tracer_lsode or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!---------------------------------------------------------------------------------
FUNCTION tracer_lsode_cast(self,source) result(success)
class(oft_tracer_lsode), pointer, intent(out) :: self !< Reference to source object with desired class
class(oft_tracer), target, intent(in) :: source !< Source tracer to cast
logical :: success !< Cast success flag
DEBUG_STACK_PUSH
select type(source)
  class is(oft_tracer_lsode)
    self=>source
    success=.TRUE.
  class default
    success=.FALSE.
end select
DEBUG_STACK_POP
end FUNCTION tracer_lsode_cast
!------------------------------------------------------------------------------
!> Setup a LSODE tracer with a new starting point
!------------------------------------------------------------------------------
subroutine tracer_lsode_setup(self,y,cell,init)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
real(r8), intent(in) :: y(:) !< New start point [3]
integer(i4), optional, intent(in) :: cell !< Guess cell for use in \ref oft_mesh_type::mesh_findcell
!! "mesh_findcell" (optional)
logical, optional, intent(in) :: init !< Flag indicating tracer is starting a new trace and all counts
!! should be set to zero (optional)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%rwork))THEN
  IF(self%mf==10)THEN
    self%lrw=20 + 16*self%neq !68
    self%liw=20
  ELSE IF(self%mf==22)THEN
    self%lrw=22 + 9*self%neq + self%neq**2
    self%liw=20 + self%neq
  END IF
  ALLOCATE(self%y(self%neq),self%dy(self%neq),self%dyp(self%neq))
  ALLOCATE(self%rwork(self%lrw),self%iwork(self%liw))
  self%dyp=0.d0
  self%rwork=0.d0
  self%iwork=0
  !---
  self%nisend=3+37+self%liw
  self%nrsend=2*self%neq+1+218+self%lrw
  ALLOCATE(self%isend(self%nisend),self%rsend(self%nrsend))
END IF
IF(PRESENT(init))THEN
  IF(init)THEN
    self%istate=1
    self%comm_load=.FALSE.
    !---
    self%ntrans=0
    self%nsteps=0
    self%status=TRACER_INIT
    self%estatus=0
  END IF
END IF
!---Find starting point
self%cell=0
if(present(cell))self%cell=cell
self%y=y
call mesh_findcell(self%ml_mesh%mesh,self%cell,self%y(1:3),self%f)
self%initialized=.TRUE.
IF(self%comm_load)THEN
  call dsrcom(self%rsav,self%isav,2)
  self%comm_load=.FALSE.
END IF
DEBUG_STACK_POP
end subroutine tracer_lsode_setup
!------------------------------------------------------------------------------
!> Create a bare copy of the LSODE tracer object
!------------------------------------------------------------------------------
subroutine tracer_lsode_copy(self,new)
class(oft_tracer_lsode), intent(in) :: self !< Tracer object
class(oft_tracer), pointer, intent(out) :: new !< Copy of object with same tolerances and referenced interpolator
class(oft_tracer_lsode), pointer :: newtmp
DEBUG_STACK_PUSH
allocate(oft_tracer_lsode::new)
IF(.NOT.tracer_lsode_cast(newtmp,new))CALL oft_abort('Failure to allocate LSODE.', &
  'tracer_lsode_copy',__FILE__)
new%B=>self%B
new%maxsteps=self%maxsteps
new%tol=self%tol
new%maxtrans=self%maxtrans
new%neq=self%neq
IF(ASSOCIATED(self%ydot))new%ydot=>self%ydot
newtmp%mf=self%mf
DEBUG_STACK_POP
end subroutine tracer_lsode_copy
!------------------------------------------------------------------------------
!> Save LSODE common block data to arrays
!------------------------------------------------------------------------------
subroutine tracer_lsode_save(self)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
DEBUG_STACK_PUSH
call dsrcom(self%rsav,self%isav,1)
self%comm_load=.TRUE.
DEBUG_STACK_POP
end subroutine tracer_lsode_save
!------------------------------------------------------------------------------
!> Send tracer information to a new domain
!!
!! This subroutine transfers the number of steps, tracer position and LSODE
!! working arrays. Data is packed into a single array for each real and integer
!! data types and non-blocking MPI sends are created to transfer each array
!------------------------------------------------------------------------------
function tracer_lsode_send(self,proc) result(sid)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
integer(i4), intent(in) :: proc !< Processor ID to send tracer information to
#ifdef OFT_MPI_F08
type(mpi_request) :: sid(2)
#else
integer(i4) :: sid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
#ifdef HAVE_MPI
integer(i4) :: ierr
DEBUG_STACK_PUSH
!---Pack info
self%isend(1)=self%nsteps
self%isend(2)=self%estatus
self%isend(3)=self%istate
self%isend(4:40)=self%isav
self%isend(41:self%nisend)=self%iwork
self%rsend(1:self%neq)=self%y
self%rsend(self%neq+1:2*self%neq)=self%dy
self%rsend(2*self%neq+1)=self%t
self%rsend(2*self%neq+2:2*self%neq+219)=self%rsav
self%rsend(2*self%neq+220:self%nrsend)=self%rwork
!---Set send flags
sid=MPI_REQUEST_NULL
!---Send to all other procs
CALL MPI_ISEND(self%rsend,self%nrsend,OFT_MPI_R8,proc,self%tids(1),oft_env%COMM,sid(1),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','tracer_lsode_send',__FILE__)
CALL MPI_ISEND(self%isend,self%nisend,OFT_MPI_I4,proc,self%tids(2),oft_env%COMM,sid(2),ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','tracer_lsode_send',__FILE__)
self%status=TRACER_SEND_ACTIVE
DEBUG_STACK_POP
#else
sid=MPI_REQUEST_NULL
CALL oft_abort("Send requested without MPI","tracer_lsode_send",__FILE__)
#endif
end function tracer_lsode_send
!------------------------------------------------------------------------------
!> Recieve tracer information from a different domain
!!
!! This subroutine sets up non-blocking MPI recieve requests and unpacks the
!! resulting data sent by \ref tracing::tracer_lsode_send "tracer_lsode_send"
!------------------------------------------------------------------------------
function tracer_lsode_recv(self,proc) result(rid)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
integer(i4), intent(in) :: proc !< Processor ID to receive tracer information from
#ifdef OFT_MPI_F08
type(mpi_request) :: rid(2)
#else
integer(i4) :: rid(2) !< MPI request IDs for transfer of real and integer data (2)
#endif
#ifdef HAVE_MPI
integer(i4) :: ierr
logical :: testr,testi
DEBUG_STACK_PUSH
IF(oft_mpi_check_reqs(2,self%rids))THEN
  !---Setup recv
  CALL MPI_IRECV(self%rsend,self%nrsend,OFT_MPI_R8,proc,self%tids(1),oft_env%COMM,self%rids(1),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','tracer_lsode_recv',__FILE__)
  CALL MPI_IRECV(self%isend,self%nisend,OFT_MPI_I4,proc,self%tids(2),oft_env%COMM,self%rids(2),ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','tracer_lsode_recv',__FILE__)
  self%status=TRACER_RECV_ACTIVE
  rid=self%rids
ELSE
  !---Test for completion of recv
  call MPI_TEST(self%rids(1),testr,MPI_STATUS_IGNORE,ierr)
  call MPI_TEST(self%rids(2),testi,MPI_STATUS_IGNORE,ierr)
  rid=self%rids
  !---Set flag and restore if complete
  IF(testr.AND.testi)THEN
    self%nsteps=self%isend(1)
    self%estatus=self%isend(2)
    self%istate=self%isend(3)
    self%isav=self%isend(4:40)
    self%iwork=self%isend(41:self%nisend)
    self%y=self%rsend(1:self%neq)
    self%dy=self%rsend(self%neq+1:2*self%neq)
    self%t=self%rsend(2*self%neq+1)
    self%rsav=self%rsend(2*self%neq+2:2*self%neq+219)
    self%rwork=self%rsend(2*self%neq+220:self%nrsend)
    self%comm_load=.TRUE.
    rid=MPI_REQUEST_NULL
    self%rids=MPI_REQUEST_NULL
  END IF
END IF
DEBUG_STACK_POP
#else
rid=MPI_REQUEST_NULL
CALL oft_abort("Receive requested without MPI","tracer_lsode_recv",__FILE__)
#endif
end function tracer_lsode_recv
!------------------------------------------------------------------------------
!> Advancing tracer by one step using LSODE
!------------------------------------------------------------------------------
subroutine tracer_lsode_step(self)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
real(r8) :: fmin,fmax
DEBUG_STACK_PUSH
self%dyp=self%dy
call dlsode(tracer_lsode_eval_B,self%neq,self%y,self%t,self%tout,1,self%tol,self%tol,self%itask, &
  self%istate,0,self%rwork,self%lrw,self%iwork,self%liw,tracer_lsode_eval_B,self%mf)
self%nsteps=self%nsteps+1
IF(self%istate<0.OR.ANY(self%y(1:3)>1.d20))THEN
  self%cell=0
  self%estatus=TRACER_ERROR_FAIL
END IF
IF(self%cell/=0)THEN ! Check for exit mesh
  call mesh_findcell(self%ml_mesh%mesh,self%cell,self%y(1:3),self%f)
  fmin=minval(self%f)
  fmax=maxval(self%f)
  IF(( fmax>=1.d0+offmesh_tol ).OR.( fmin<=-offmesh_tol ))self%cell=0
END IF
DEBUG_STACK_POP
end subroutine tracer_lsode_step
!------------------------------------------------------------------------------
!> Field evaluation callback used by `dlsode`
!!
!! See \ref tracing::tracer_lsode_step "tracer_lsode_step"
!------------------------------------------------------------------------------
subroutine tracer_lsode_eval_B(neq,t,y,ydot)
integer, intent(in) :: neq
real(r8), intent(in) :: t,y(neq)
real(r8), intent(out) :: ydot(neq)
real(r8) :: goptmp(3,4),v,B(3)
DEBUG_STACK_PUSH
ydot=active_tracer%dy
call mesh_findcell(active_tracer%ml_mesh%mesh,active_tracer%cell,y(1:3),active_tracer%f)
IF(active_tracer%cell==0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
CALL active_tracer%ml_mesh%mesh%jacobian(active_tracer%cell,active_tracer%f,goptmp,v)
CALL active_tracer%B%interp(active_tracer%cell,active_tracer%f,goptmp,B)
IF(ASSOCIATED(active_tracer%ydot))THEN
  CALL active_tracer%ydot(t,y,B,neq,ydot)
ELSE
  ydot=B
END IF
active_tracer%dy=ydot
DEBUG_STACK_POP
end subroutine tracer_lsode_eval_B
!------------------------------------------------------------------------------
!> Delete internal storage and reset counters
!------------------------------------------------------------------------------
subroutine tracer_lsode_delete(self)
class(oft_tracer_lsode), intent(inout) :: self !< Tracer object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%rwork))THEN
  DEALLOCATE(self%y,self%dy,self%dyp)
  DEALLOCATE(self%rwork,self%iwork,self%isend,self%rsend)
END IF
self%nsteps=0
self%ntrans=0
self%initialized=.FALSE.
self%comm_load=.FALSE.
NULLIFY(self%ydot)
DEBUG_STACK_POP
end subroutine tracer_lsode_delete
end module tracing
