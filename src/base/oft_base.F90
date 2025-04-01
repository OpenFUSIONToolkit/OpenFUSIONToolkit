!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_base.F90
!
!> @defgroup doxy_oft_base Open FUSION Toolkit Base
!! Enviroment functions and aliases for the Open FUSION Toolkit
!
!> Base environment functions and aliases for Open FUSION Toolkit.
!! - MPI aliases for Open FUSION Toolkit environment
!! - Runtime settings and identifiers
!! - Init, Abort, Finalize functions
!! oft_env_type class
!! - MPI context information
!! - MPI send/recv tags
!!
!! @authors Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_base
!--------------------------------------------------------------------------------
MODULE oft_base
USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit
USE omp_lib
USE oft_local
USE oft_sort, ONLY: sort_array
#ifdef HAVE_PETSC
USE petscsys
#else
#ifdef OFT_MPI_F08
USE mpi_f08
#endif
#endif
IMPLICIT NONE
#if defined(HAVE_MPI) && !defined(OFT_MPI_F08) && !defined(HAVE_PETSC)
#include "mpif.h"
#endif
#include "local.h"
#include "git_info.h"
!---MPI Type Aliases
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
TYPE(mpi_datatype) :: OFT_MPI_I8=MPI_INTEGER8 !< MPI_INT8 alias
TYPE(mpi_datatype) :: OFT_MPI_I4=MPI_INTEGER4 !< MPI_INT4 alias
TYPE(mpi_datatype) :: OFT_MPI_R4=MPI_REAL8 !< MPI_REAL4 alias
TYPE(mpi_datatype) :: OFT_MPI_R8=MPI_REAL8 !< MPI_REAL8 alias
TYPE(mpi_datatype) :: OFT_MPI_C4=MPI_COMPLEX8 !< MPI_COMPLEX8 alias
TYPE(mpi_datatype) :: OFT_MPI_C8=MPI_COMPLEX16 !< MPI_COMPLEX16 alias
TYPE(mpi_datatype) :: OFT_MPI_LOGICAL=MPI_LOGICAL !< MPI_LOGICAL alias
TYPE(mpi_datatype) :: OFT_MPI_CHAR=MPI_CHARACTER !< MPI_CHAR alias
#else
INTEGER(i4) :: OFT_MPI_I8=MPI_INTEGER8 !< MPI_INT8 alias
INTEGER(i4) :: OFT_MPI_I4=MPI_INTEGER4 !< MPI_INT4 alias
INTEGER(i4) :: OFT_MPI_R4=MPI_REAL8 !< MPI_REAL4 alias
INTEGER(i4) :: OFT_MPI_R8=MPI_REAL8 !< MPI_REAL8 alias
INTEGER(i4) :: OFT_MPI_C4=MPI_COMPLEX8 !< MPI_COMPLEX8 alias
INTEGER(i4) :: OFT_MPI_C8=MPI_COMPLEX16 !< MPI_COMPLEX16 alias
INTEGER(i4) :: OFT_MPI_LOGICAL=MPI_LOGICAL !< MPI_LOGICAL alias
INTEGER(i4) :: OFT_MPI_CHAR=MPI_CHARACTER !< MPI_CHAR alias
#endif
#else
INTEGER(i4), PARAMETER :: MPI_COMM_WORLD = -1 ! Dummy comm value for non-MPI runs
INTEGER(i4), PARAMETER :: MPI_COMM_NULL = -2 ! Dummy null comm value for non-MPI runs
INTEGER(i4), PARAMETER :: MPI_REQUEST_NULL = -3 ! Dummy null request value for non-MPI runs
INTEGER(i4), PARAMETER :: MPI_COMM_SELF = -4 ! Dummy self-comm value for non-MPI runs
#endif
!------------------------------------------------------------------------------
!> Perform a SUM/AND reduction over all processors
!------------------------------------------------------------------------------
INTERFACE oft_mpi_sum
  MODULE PROCEDURE oft_mpi_sumr
  MODULE PROCEDURE oft_mpi_sumra
  MODULE PROCEDURE oft_mpi_sumc
  MODULE PROCEDURE oft_mpi_sumca
  MODULE PROCEDURE oft_mpi_sumi4
  MODULE PROCEDURE oft_mpi_sumi4a
  MODULE PROCEDURE oft_mpi_sumi8
  MODULE PROCEDURE oft_mpi_sumi8a
END INTERFACE oft_mpi_sum
PRIVATE oft_mpi_sumr, oft_mpi_sumra, oft_mpi_sumc, oft_mpi_sumca, &
  oft_mpi_sumi4, oft_mpi_sumi4a, oft_mpi_sumi8, oft_mpi_sumi8a
!------------------------------------------------------------------------------
!> Perform a MAX reduction over all processors
!------------------------------------------------------------------------------
INTERFACE oft_mpi_max
  MODULE PROCEDURE oft_mpi_maxr
  MODULE PROCEDURE oft_mpi_maxra
  MODULE PROCEDURE oft_mpi_maxi
  MODULE PROCEDURE oft_mpi_maxia
  MODULE PROCEDURE oft_mpi_maxi8a
END INTERFACE oft_mpi_max
PRIVATE oft_mpi_maxr, oft_mpi_maxra, oft_mpi_maxi, oft_mpi_maxia
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
INTERFACE oft_random_number
  MODULE PROCEDURE oft_random_number_r8
  MODULE PROCEDURE oft_random_number_c8
END INTERFACE oft_random_number
PRIVATE oft_random_number_r8
!------------------------------------------------------------------------------
!> Dummy shadow type for Fox XML node
!------------------------------------------------------------------------------
#if !defined(HAVE_XML)
TYPE :: xml_node
  INTEGER(i4) :: dummy = 0
END TYPE
#endif
!------------------------------------------------------------------------------
!> Open FUSION Toolkit environment class
!!
!! Contains runtime enviroment information
!! - Global MPI COMM alias
!! - processor context
!------------------------------------------------------------------------------
TYPE :: oft_env_type
  INTEGER(i4) :: nbase = -1 !< Number of OpenMP base meshes
  INTEGER(i4) :: nparts = 1 !< Number of OpenMP paritions
#ifdef OFT_MPI_F08
  TYPE(mpi_comm) :: COMM = MPI_COMM_WORLD !< Open FUSION Toolkit MPI communicator
  TYPE(mpi_comm) :: NODE_COMM = MPI_COMM_SELF !< Open FUSION Toolkit MPI node-local communicator
#else
  INTEGER(i4) :: COMM = MPI_COMM_WORLD !< Open FUSION Toolkit MPI communicator
  INTEGER(i4) :: NODE_COMM = MPI_COMM_SELF !< Open FUSION Toolkit MPI node-local communicator
#endif
  INTEGER(i4) :: nnodes = -1 !< Number of MPI tasks
  INTEGER(i4) :: ppn = 1 !< Number of procs per NUMA node
  INTEGER(i4) :: nprocs = -1 !< Number of MPI tasks
  INTEGER(i4) :: nthreads = -1 !< Number of OpenMP threads
  INTEGER(i4) :: rank = -1 !< MPI rank
  INTEGER(i4) :: node_rank = -1 !< MPI node-local rank
  INTEGER(i4) :: debug = 0 !< Debug level (1-3)
  LOGICAL :: head_proc = .FALSE. !< Lead processor flag
  LOGICAL :: pm = .TRUE. !< Performance monitor (default T=on, F=off)
  LOGICAL :: test_run = .FALSE. !< Test run
  CHARACTER(LEN=4) :: crank = '' !< Processor rank in character form
  CHARACTER(LEN=OFT_PATH_SLEN) :: ifile = 'none' !< Name of input file
  CHARACTER(LEN=OFT_PATH_SLEN) :: xml_file = 'none' !< Name of XML input file
#ifdef HAVE_XML
  TYPE(xml_node), POINTER :: xml => NULL()
#endif
END TYPE oft_env_type
!---Global variables
TYPE(oft_env_type) :: oft_env !< Global container for environment information
character(len=:), allocatable :: oft_indent !< Indentation string for status messages
INTEGER(i4) :: oft_tid = 0 !< ID of current thread
REAL(r8), PRIVATE :: start_time !< Time that run began
INTEGER(i8) :: comm_times(4) = 0 !< Times for blocking communications
LOGICAL :: use_petsc = .FALSE. !< Use PETSc as linear algebra backend
INTEGER(i4), PARAMETER :: oft_test_seed(24) = [430470439, -303393496, -476850581, &
  -964913795, 995391627, 84909391, -586395292, -2070086573, -1010035798, 1012650827, &
  325297911, 701378007, 392909068, 379156631, 484729024, -292308758, -1669043581, &
  142231192, 708877466, -1255634259, 593274827, -561530186, -934579426, 900810854]
!$omp threadprivate(oft_tid)
!---Debugging stack information
LOGICAL :: oft_stack_disabled = .FALSE. !< Disable debug/profiling stack
LOGICAL, PRIVATE :: oft_prof_enabled = .FALSE. !< Enable profiling?
#ifdef OFT_STACK
INTEGER(i4), PRIVATE, PARAMETER :: stacksize = 40
INTEGER(i4), PRIVATE :: nstack = 0
INTEGER(i4), PRIVATE :: stack(2,stacksize) = 0
!$omp threadprivate(nstack,stack)
#include "stack_defs.h"
#ifdef OFT_PROFILE
TYPE(oft_timer), PRIVATE :: local_timer
INTEGER(i4), PRIVATE :: stack_nfun_c(stack_nfuns) = 0_i4
INTEGER(i8), PRIVATE :: stack_fun_time(0:stack_nfuns) = 0_i8
!$omp threadprivate(local_timer,stack_nfun_c,stack_fun_time)
#endif
#endif
INTERFACE
!------------------------------------------------------------------------------
!> Interface for setting UNIX signal handlers in "oft_local.c"
!------------------------------------------------------------------------------
  SUBROUTINE oft_set_signal_handlers()  BIND(C)
  END SUBROUTINE oft_set_signal_handlers
!------------------------------------------------------------------------------
!> Prototype for abort callback to override usual abort process
!------------------------------------------------------------------------------
  SUBROUTINE oft_abort_callback()  BIND(C)
  END SUBROUTINE oft_abort_callback
END INTERFACE
!> Abort callback for graceful abort in Python interface
PROCEDURE(oft_abort_callback), POINTER :: oft_abort_cb
CONTAINS
!------------------------------------------------------------------------------
!> Initializes Open FUSION Toolkit run environment
!!
!! Also calls MPI_INIT
!------------------------------------------------------------------------------
SUBROUTINE oft_init(nthreads)
INTEGER(i4), INTENT(in), OPTIONAL :: nthreads !< Number for threads to use (negative for default)
INTEGER(i4) :: ierr,thrdtype,nargs,io_unit
REAL(r8) :: elapsed_time
INTEGER(i4) :: ppn=1
INTEGER(i4) :: debug=0
INTEGER(i4) :: nparts=1
INTEGER(i4) :: omp_nthreads=-1
LOGICAL :: test_run=.FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: ifile
LOGICAL :: called_from_lib
#ifdef HAVE_XML
TYPE(xml_node), POINTER :: doc
#endif
LOGICAL :: rst
NAMELIST/runtime_options/ppn,omp_nthreads,debug,oft_stack_disabled,use_petsc,test_run,nparts
!---Initialize MPI
#ifdef HAVE_MPI
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,thrdtype,ierr)
#endif
#ifdef OFT_STACK
CALL oft_set_signal_handlers
#endif
!---Get MPI task information
#ifdef HAVE_MPI
CALL MPI_COMM_RANK(oft_env%comm,oft_env%rank,ierr)
CALL MPI_COMM_SIZE(oft_env%comm,oft_env%nprocs,ierr)
#else
oft_env%rank=0
oft_env%nprocs=1
#endif
IF(oft_env%rank==0)oft_env%head_proc=.TRUE.
WRITE(oft_env%crank,'(I4.4)')oft_env%rank
!---Read settings filename from command line
IF(PRESENT(nthreads))THEN
  omp_nthreads=nthreads
ELSE
  nargs=COMMAND_ARGUMENT_COUNT()
  oft_env%ifile=TRIM('oft.in') ! If none, use default
  oft_env%xml_file=TRIM('none') ! If none, specify
  IF(nargs>0)THEN
    CALL GET_COMMAND_ARGUMENT(1,ifile)
    IF(ifile(1:1)/="-")oft_env%ifile=TRIM(ifile)
    IF(nargs>1)THEN
      CALL GET_COMMAND_ARGUMENT(2,ifile)
      IF(ifile(1:1)/="-")oft_env%xml_file=TRIM(ifile)
    END IF
  END IF
END IF
!---Test for existence of settings file
INQUIRE(file=oft_env%ifile,exist=rst)
IF(.NOT.rst)CALL oft_abort('Input file does not exist.','oft_init',__FILE__)
!---Read in node options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,runtime_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No runtime options found in input file.','oft_init',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing runtime options.','oft_init',__FILE__)
!---Seed pRNG if test run for repeatability (does not work on < GCC 9)
#if !defined(__GNUC__) || (__GNUC__ > 9)
IF(test_run)THEN
  CALL random_seed(size=nargs)
  IF(nargs>SIZE(oft_test_seed))CALL oft_abort('pRNG seed size exceeds built in values', &
    'oft_init',__FILE__)
  CALL random_seed(put=oft_test_seed)
END IF
#endif
!---Initialize PETSc or exit if requested but not available
IF(use_petsc)THEN
#ifdef HAVE_PETSC
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#else
  IF(oft_env%head_proc)WRITE(*,*)'ERROR: Not compiled with PETSc'
  CALL oft_finalize()
#endif
END IF
!---Set run timer
start_time=omp_get_wtime()
#ifdef OFT_PROFILE
!$omp parallel
CALL local_timer%tick()
elapsed_time=local_timer%tock()
!$omp end parallel
#endif
!---Get number of OpenMP threads
IF(omp_nthreads>0)CALL omp_set_num_threads(omp_nthreads)
!$omp parallel
!$omp single
oft_env%nthreads=omp_get_num_threads()
!$omp end single
oft_tid=omp_get_thread_num()
!$omp end parallel
IF(PRESENT(nthreads))THEN
  IF((nthreads>0).AND.(oft_env%nthreads/=MAX(omp_nthreads,1)))CALL oft_abort( &
  'Error setting number of threads','oft_init',__FILE__)
END IF
!---Parse xml file if necessary
IF(oft_env%xml_file(1:4)/='none')THEN
#ifdef HAVE_XML
  !---Test for existence of XML file
  INQUIRE(FILE=TRIM(oft_env%xml_file),exist=rst)
  IF(.NOT.rst)CALL oft_abort('XML file specified but does not exist.','oft_init',__FILE__)
  doc=>xml_parseFile(TRIM(oft_env%xml_file),iostat=ierr)
  IF(ierr/=0)CALL oft_abort('Error parsing XML input file','oft_init',__FILE__)
  CALL xml_get_element(doc,"oft",oft_env%xml,ierr)
  IF(ierr/=0)CALL oft_abort('Error finding "oft" XML root element','oft_init',__FILE__)
#else
  CALL oft_warn("Open FUSION Toolkit not built wit xml support, ignoring xml input.")
#endif
END IF
!---
IF(MOD(oft_env%nprocs,ppn)/=0)CALL oft_abort('# of MPI tasks and Procs per node do not agree.','oft_init',__FILE__)
oft_env%nnodes=oft_env%nprocs/ppn
oft_env%ppn=ppn
oft_env%debug=debug
oft_env%test_run=test_run
oft_env%nparts=nparts
oft_indent=""
#ifdef HAVE_MPI
IF(oft_env%ppn>1)THEN
  CALL MPI_Comm_split(oft_env%comm,oft_env%rank/oft_env%ppn,oft_env%rank,oft_env%NODE_COMM,ierr)
  CALL MPI_Comm_rank(oft_env%NODE_COMM,oft_env%node_rank,ierr)
ELSE
  oft_env%node_rank=0
END IF
#endif
!---Print runtime information
IF(oft_env%rank==0)THEN
  WRITE(*,'(A)')    '#----------------------------------------------'
  WRITE(*,'(A)')    'Open FUSION Toolkit Initialized'
  WRITE(*,'(2A)')   'Development branch:   ',GITBRANCH
  WRITE(*,'(2A)')   'Revision id:          ',GITVERSION
  WRITE(*,'(2A)')   'Parallelization Info:'
#ifdef HAVE_MPI
  WRITE(*,'(A,I4)') '  # of MPI tasks      = ',oft_env%nprocs
  WRITE(*,'(A,I4)') '  # of NUMA nodes     = ',oft_env%nnodes
#else
  WRITE(*,'(A)')    '  Not compiled with MPI'
#endif
#if defined(_OPENMP)
  WRITE(*,'(A,I4)') '  # of OpenMP threads = ',oft_env%nthreads
#else
  WRITE(*,'(A)')    '  Not compiled with OpenMP'
#endif
  WRITE(*,'(2A)')   'Fortran input file    = ',TRIM(oft_env%ifile)
  WRITE(*,'(2A)')   'XML input file        = ',TRIM(oft_env%xml_file)
  WRITE(*,'(A,3I4)')'Integer Precisions    = ',i4,i8
  WRITE(*,'(A,3I4)')'Float Precisions      = ',r4,r8,r10
  WRITE(*,'(A,3I4)')'Complex Precisions    = ',c4,c8
IF(use_petsc)THEN
  WRITE(*,'(A)')    'LA backend            = PETSc'
ELSE
  WRITE(*,'(A)')    'LA backend            = native'
END IF
  WRITE(*,'(A)')    '#----------------------------------------------'
  WRITE(*,*)
END IF
END SUBROUTINE oft_init
!------------------------------------------------------------------------------
!> Finalize Open FUSION Toolkit environment
!!
!! Also calls PetscFinalize/MPI_FINALIZE
!------------------------------------------------------------------------------
SUBROUTINE oft_finalize() BIND(C)
INTEGER(i4) :: ierr
INTEGER(i8) :: comm_min(4),comm_max(4),comm_avg(4)
INTEGER(i8) :: countnew,crate,cmax
REAL(i8) :: tmp(3)
#ifdef HAVE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_ALLREDUCE(comm_times,comm_min,4,OFT_MPI_I8,MPI_MIN,oft_env%COMM,ierr)
CALL MPI_ALLREDUCE(comm_times,comm_max,4,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
CALL MPI_ALLREDUCE(comm_times,comm_avg,4,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
#endif
!---Print timing information
IF(oft_env%rank==0)THEN
  CALL system_clock(countnew,crate,cmax)
  WRITE(*,*)
  WRITE(*,'(A)')         '#----------------------------------------------'
  WRITE(*,'(A)')         'Run Complete'
  WRITE(*,'(A,ES11.3)')  'Run Duration = ',omp_get_wtime()-start_time
#ifdef HAVE_MPI
  WRITE(*,*)
  WRITE(*,'(A)')         '====== MPI Times (Avg, Min, Max) ======'
  tmp=(/comm_avg(1),comm_min(1),comm_max(1)/)/REAL(crate,8)
  WRITE(*,'(A,3ES11.3)') 'Barrier   = ',tmp(1)/oft_env%nprocs,tmp(2),tmp(3)
  tmp=(/comm_avg(2),comm_min(2),comm_max(2)/)/REAL(crate,8)
  WRITE(*,'(A,3ES11.3)') 'WaitAny   = ',tmp(1)/oft_env%nprocs,tmp(2),tmp(3)
  tmp=(/comm_avg(3),comm_min(3),comm_max(3)/)/REAL(crate,8)
  WRITE(*,'(A,3ES11.3)') 'WaitAll   = ',tmp(1)/oft_env%nprocs,tmp(2),tmp(3)
  tmp=(/comm_avg(4),comm_min(4),comm_max(4)/)/REAL(crate,8)
  WRITE(*,'(A,3ES11.3)') 'AllReduce = ',tmp(1)/oft_env%nprocs,tmp(2),tmp(3)
#endif
  WRITE(*,'(A)')         '#----------------------------------------------'
END IF
!---
CALL oft_prof_print
#ifdef HAVE_PETSC
IF(use_petsc)CALL PetscFinalize(ierr)
#endif
!---Finalize MPI environment
#ifdef HAVE_MPI
CALL MPI_FINALIZE(ierr)
#endif
STOP
END SUBROUTINE oft_finalize
!------------------------------------------------------------------------------
!> Graceful abort for Open FUSION Toolkit
!!
!! Also calls MPI_ABORT/STOP
!------------------------------------------------------------------------------
SUBROUTINE oft_abort(error_str,sname,fname)
CHARACTER(LEN=*), INTENT(in) :: error_str !< Error string
CHARACTER(LEN=*), INTENT(in) :: sname !< Source subroutine name
CHARACTER(LEN=*), INTENT(in) :: fname !< Source file name
INTEGER(i4) :: ierr,errcode,outunit
#ifdef OFT_ABORT_FILES
CHARACTER(LEN=5) :: proc
#endif
outunit=error_unit
#ifdef OFT_ABORT_FILES
WRITE(proc,'(I5.5)')oft_env%rank
OPEN(outunit,FILE='abort_'//proc//'.err')
#endif
!---Print error information
100 FORMAT (A,I5,2A)
101 FORMAT (2A)
WRITE(outunit,'(X)')
WRITE(outunit,'(A)')'#----------------------------------------------'
IF(oft_env%nprocs>1)THEN
  WRITE(outunit,100)  '[',oft_env%rank,'] ERROR: ',TRIM(error_str)
  WRITE(outunit,100)  '[',oft_env%rank,'] SUBROUTINE: ',TRIM(sname)
  WRITE(outunit,100)  '[',oft_env%rank,'] FILE: ',TRIM(fname)
ELSE
  WRITE(outunit,101)  'ERROR: ',TRIM(error_str)
  WRITE(outunit,101)  'SUBROUTINE: ',TRIM(sname)
  WRITE(outunit,101)  'FILE: ',TRIM(fname)
END IF
WRITE(outunit,'(A)')'#----------------------------------------------'
WRITE(outunit,'(X)')
#ifdef OFT_ABORT_FILES
CLOSE(outunit)
#endif
!---
CALL oft_stack_print
!---
IF(ASSOCIATED(oft_abort_cb))CALL oft_abort_cb
!---Abort run
errcode=99
#ifdef HAVE_MPI
CALL MPI_ABORT(oft_env%comm,errcode,ierr)
#else
ERROR STOP errcode
#endif
END SUBROUTINE oft_abort
!------------------------------------------------------------------------------
!> Graceful warning printing for Open FUSION Toolkit
!------------------------------------------------------------------------------
SUBROUTINE oft_warn(error_str)
CHARACTER(LEN=*) :: error_str
!---Print warning information
100 FORMAT (A,I5,2A)
101 FORMAT (2A)
IF(oft_env%nprocs>1)THEN
  WRITE(error_unit,100)'[',oft_env%rank,'] WARNING: ',TRIM(error_str)
ELSE
  WRITE(error_unit,101)'WARNING: ',TRIM(error_str)
END IF
END SUBROUTINE oft_warn
!------------------------------------------------------------------------------
!> Output control for performance messages
!!
!! @result oft_env\%pm.AND.oft_env\%head_proc
!------------------------------------------------------------------------------
PURE FUNCTION oft_pm_print() RESULT(pflag)
LOGICAL :: pflag
pflag=(oft_env%pm.AND.oft_env%head_proc)
END FUNCTION oft_pm_print
!------------------------------------------------------------------------------
!> Output control for debug messages
!!
!! @result (oft_env\%debug>=level).AND.oft_env\%head_proc
!------------------------------------------------------------------------------
PURE FUNCTION oft_debug_print(level) RESULT(pflag)
INTEGER(i4), INTENT(in) :: level !< Threshold debugging level
LOGICAL :: pflag
pflag=((oft_env%debug>=level).AND.oft_env%head_proc)
END FUNCTION oft_debug_print
!---------------------------------------------------------------------------------
!> Get thread ownership range for even spacing
!---------------------------------------------------------------------------------
PURE SUBROUTINE oft_thread_slice(tid,nthreads,length,i1,i2)
INTEGER(i4), INTENT(in) :: tid !< Thread index (0-indexed)
INTEGER(i4), INTENT(in) :: nthreads !< Total number of threads
INTEGER(i4), INTENT(in) :: length !< Length of range to partition
INTEGER(i4), INTENT(out) :: i1 !< Start of ownership range (1-indexed)
INTEGER(i4), INTENT(out) :: i2 !< End of ownership range (1-indexed)
INTEGER(i4) :: i,stride
stride = FLOOR(length/REAL(nthreads,8))
i1 = 1
i2 = stride
DO i=1,tid
  i1 = i1 + stride
  i2 = i2 + stride
END DO
IF(tid==(nthreads-1))i2=length
END SUBROUTINE oft_thread_slice
!------------------------------------------------------------------------------
!> Increase length of global indent string by 2
!------------------------------------------------------------------------------
SUBROUTINE oft_increase_indent
oft_indent=oft_indent//"  "
END SUBROUTINE oft_increase_indent
!------------------------------------------------------------------------------
!> Decrease length of global indent string by 2
!------------------------------------------------------------------------------
SUBROUTINE oft_decrease_indent
INTEGER(i4) :: nindent
nindent = LEN(oft_indent)
IF(nindent>2)THEN
  oft_indent=oft_indent(1:nindent-2)
ELSE
  oft_indent=""
END IF
END SUBROUTINE oft_decrease_indent
!------------------------------------------------------------------------------
!> Wrapper for MPI_BARRIER with MPI_COMM_WORLD to enable communication profiling
!------------------------------------------------------------------------------
SUBROUTINE oft_mpi_barrier(ierr)
INTEGER(i4), INTENT(out) :: ierr !< Error flag
INTEGER(i8) :: timein
#ifdef HAVE_MPI
timein=oft_time_i8()
CALL MPI_BARRIER(oft_env%comm,ierr)
!$omp atomic
comm_times(1)=comm_times(1)+oft_time_diff(timein)
#else
ierr=0
#endif
END SUBROUTINE oft_mpi_barrier
!------------------------------------------------------------------------------
!> Wrapper for MPI_WAITANY to enable communication profiling
!------------------------------------------------------------------------------
SUBROUTINE oft_mpi_waitany(n,req,j,ierr)
INTEGER(i4), INTENT(in) :: n !< Number of requests
#ifdef OFT_MPI_F08
TYPE(mpi_request), INTENT(inout) :: req(n) !< Array of requests
#else
INTEGER(i4), INTENT(inout) :: req(n) !< Array of requests
#endif
INTEGER(i4), INTENT(out) :: j !< Next completed request 
INTEGER(i4), INTENT(out) :: ierr !< Error flag
INTEGER(i8) :: timein
#ifdef HAVE_MPI
timein=oft_time_i8()
CALL MPI_WAITANY(n,req,j,MPI_STATUS_IGNORE,ierr)
!$omp atomic
comm_times(2)=comm_times(2)+oft_time_diff(timein)
#else
ierr=0
j=-1
#endif
END SUBROUTINE oft_mpi_waitany
!------------------------------------------------------------------------------
!> Wrapper for MPI_WAITALL to enable communication profiling
!------------------------------------------------------------------------------
SUBROUTINE oft_mpi_waitall(n,req,ierr)
INTEGER(i4), INTENT(in) :: n !< Number of requests
#ifdef OFT_MPI_F08
TYPE(mpi_request), INTENT(inout) :: req(n) !< Array of requests
#else
INTEGER(i4), INTENT(inout) :: req(n) !< Array of requests
#endif
INTEGER(i4), INTENT(out) :: ierr !< Error flag
INTEGER(i8) :: timein
#ifdef HAVE_MPI
timein=oft_time_i8()
CALL MPI_WAITALL(n,req,MPI_STATUSES_IGNORE,ierr)
!$omp atomic
comm_times(3)=comm_times(3)+oft_time_diff(timein)
#else
ierr=0
#endif
END SUBROUTINE oft_mpi_waitall
!------------------------------------------------------------------------------
!> Helper to check all requests from completion
!------------------------------------------------------------------------------
FUNCTION oft_mpi_check_reqs(n,req) result(all_null)
INTEGER(i4), INTENT(in) :: n !< Number of requests
#ifdef OFT_MPI_F08
TYPE(mpi_request), INTENT(inout) :: req(n) !< Array of requests
#else
INTEGER(i4), INTENT(inout) :: req(n) !< Array of requests
#endif
INTEGER(i4) :: i
LOGICAL :: all_null
all_null=.TRUE.
#ifdef HAVE_MPI
DO i=1,n
  all_null=all_null.AND.(req(i)==MPI_REQUEST_NULL)
END DO
#endif
END FUNCTION oft_mpi_check_reqs
!---------------------------------------------------------------------------------
!> real(r8) scalar implementation of global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumr(a) result(b)
REAL(r8), INTENT(in) :: a !< Local value for SUM
REAL(r8) :: b
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumr',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumr
!---------------------------------------------------------------------------------
!> real(r8) array implementation of element-wise global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumra(a,n) result(b)
REAL(r8), INTENT(in) :: a(n) !< Local values for SUM [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
REAL(r8) :: b(n)
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumra',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumra
!---------------------------------------------------------------------------------
!> complex(c8) scalar implementation of global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumc(a) result(b)
COMPLEX(r8), INTENT(in) :: a !< Local value for SUM
COMPLEX(r8) :: b
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_C8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumc',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumc
!---------------------------------------------------------------------------------
!> complex(c8) array implementation of element-wise global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumca(a,n) result(b)
COMPLEX(c8), INTENT(in) :: a(n) !< Local values for SUM [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
COMPLEX(c8) :: b(n)
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_C8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumca',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumca
!---------------------------------------------------------------------------------
!> integer(i4) scalar implementation of global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumi4(a) result(b)
INTEGER(i4), INTENT(in) :: a !< Local value for SUM
INTEGER(i4) :: b,ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_I4,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumi4',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumi4
!---------------------------------------------------------------------------------
!> integer(i4) array implementation of element-wise global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumi4a(a,n) result(b)
INTEGER(i4), INTENT(in) :: a(n) !< Local values for SUM [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
INTEGER(i4) :: b(n),ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_I4,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumi4a',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumi4a
!---------------------------------------------------------------------------------
!> integer(i8) scalar implementation of global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumi8(a) result(b)
INTEGER(i8), INTENT(in) :: a !< Local value for SUM
INTEGER(i8) :: b
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumi8',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumi8
!---------------------------------------------------------------------------------
!> integer(i8) array implementation of element-wise global SUM reduction
!!
!! @result \f$ sum(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_sumi8a(a,n) result(b)
INTEGER(i8), INTENT(in) :: a(n) !< Local values for SUM [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
INTEGER(i8) :: b(n)
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_sumi8a',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_sumi8a
!---------------------------------------------------------------------------------
!> logical scalar implementation of global AND (SUM) reduction
!!
!! @result \f$ ALL(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_and(a) result(b)
LOGICAL, INTENT(in) :: a !< Local value for AND
LOGICAL :: b
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_LOGICAL,MPI_LAND,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_and',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_and
!---------------------------------------------------------------------------------
!> real(r8) scalar implementation of global MAX reduction
!!
!! @result \f$ max(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_maxr(a) result(b)
REAL(i8), INTENT(in) :: a !< Local value for MAX
REAL(i8) :: b
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_R8,MPI_MAX,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_maxr',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_maxr
!---------------------------------------------------------------------------------
!> real(r8) array implementation of element-wise global MAX reduction
!!
!! @result \f$ max(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_maxra(a,n) result(b)
REAL(r8), INTENT(in) :: a(n) !< Local values for MAX [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
REAL(r8) :: b(n)
INTEGER(i4) :: ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_R8,MPI_MAX,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_maxra',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_maxra
!---------------------------------------------------------------------------------
!> integer(i4) scalar implementation of global MAX reduction
!!
!! @result \f$ max(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_maxi(a) result(b)
INTEGER(i4), INTENT(in) :: a !< Local value for MAX
INTEGER(i4) :: b,ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,1,OFT_MPI_I4,MPI_MAX,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_maxi',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_maxi
!---------------------------------------------------------------------------------
!> integer(i4) array implementation of element-wise global MAX reduction
!!
!! @result \f$ max(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_maxia(a,n) result(b)
INTEGER(i4), INTENT(in) :: a(n) !< Local values for MAX [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
INTEGER(i4) :: b(n),ierr
INTEGER(i8) :: timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_I4,MPI_MAX,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_maxia',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_maxia
!---------------------------------------------------------------------------------
!> integer(i8) array implementation of element-wise global MAX reduction
!!
!! @result \f$ max(a) \f$ over all processors
!---------------------------------------------------------------------------------
FUNCTION oft_mpi_maxi8a(a,n) result(b)
INTEGER(i8), INTENT(in) :: a(n) !< Local values for MAX [n]
INTEGER(i4), INTENT(in) :: n !< Length of array for reduction
INTEGER(i4) :: ierr
INTEGER(i8) :: b(n),timein
#ifdef HAVE_MPI
DEBUG_STACK_PUSH
timein=oft_time_i8()
CALL MPI_ALLREDUCE(a,b,n,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','oft_mpi_maxi8a',__FILE__)
!$omp atomic
comm_times(4)=comm_times(4)+oft_time_diff(timein)
DEBUG_STACK_POP
#else
b=a
#endif
END FUNCTION oft_mpi_maxi8a
!---------------------------------------------------------------------------------
!> Real implementation of toolkit random_number function
!---------------------------------------------------------------------------------
SUBROUTINE oft_random_number_r8(array,n)
REAL(r8), INTENT(out) :: array(n) !< Array to set with random values [n]
INTEGER(i4), INTENT(in) :: n !< Length of array
INTEGER(i4) :: i
IF(oft_env%test_run)THEN
  DO i=1,n
    array(i) = (1.d0+SIN(REAL(i+oft_env%rank,8)))/2.d0
  END DO
ELSE
  CALL RANDOM_NUMBER(array)
END IF
END SUBROUTINE oft_random_number_r8
!---------------------------------------------------------------------------------
!> Complex implementation of toolkit random_number function
!---------------------------------------------------------------------------------
SUBROUTINE oft_random_number_c8(array,n)
COMPLEX(c8), INTENT(out) :: array(n) !< Array to set with random values [n]
INTEGER(i4), INTENT(in) :: n !< Length of array
INTEGER(i4) :: i
REAL(r8), ALLOCATABLE, DIMENSION(:) :: harvest
IF(oft_env%test_run)THEN
  DO i=1,n
    array(i) = (1.d0,0.d0)*(1.d0+SIN(REAL(i+oft_env%rank,8)))/2.d0 &
      + (0.d0,1.d0)*(1.d0+COS(REAL(i+oft_env%rank,8)))/2.d0
  END DO
ELSE
  ALLOCATE(harvest(2*n))
  CALL RANDOM_NUMBER(harvest)
  DO i=1,n
    array(i)=(1.d0,0.d0)*harvest(i) + (0.d0,1.d0)*harvest(2*i)
  END DO
END IF
END SUBROUTINE oft_random_number_c8
!---------------------------------------------------------------------------------
!> Apply an orientation transform to a 2 value array
!!
!! list([1,2]) = list([1,2]) if oflag>0
!! list([1,2]) = list([2,1]) if oflag<0
!---------------------------------------------------------------------------------
PURE SUBROUTINE orient_list2(oflag,list)
INTEGER(i4), INTENT(in) :: oflag !< Orientation flag
INTEGER(i4), INTENT(inout) :: list(2) !< Array for orientation [2]
if(oflag<0)list([2,1])=list([1,2])
END SUBROUTINE orient_list2
!---------------------------------------------------------------------------------
!> Apply an orientation transform to a n-value array
!---------------------------------------------------------------------------------
PURE SUBROUTINE find_orient_listn(oflag,list,n)
INTEGER(i4), INTENT(out) :: oflag !< Orientation flag
INTEGER(i4), INTENT(in) :: list(n) !< Array for orientation [n]
INTEGER(i4), INTENT(in) :: n !< Length of array
INTEGER(i4) :: i,ind,comp(2)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: listtmp
! Nomial ordering
ind=MINLOC(list, DIM=1)
comp(1)=list(ind)
IF(ind==1)THEN
  comp(2)=MIN(list(n),list(2))
ELSEIF(ind==n)THEN
  comp(2)=MIN(list(n-1),list(1))
ELSE
  comp(2)=MIN(list(ind-1),list(ind+1))
END IF
!
oflag=0
ALLOCATE(listtmp(n))
listtmp=list
DO i=1,n
  IF(ALL(listtmp(1:2)==comp))THEN
    oflag=i
    RETURN
  END IF
  listtmp=CSHIFT(listtmp,1)
END DO
listtmp=list(n:1:-1) ! Reverse list
DO i=1,n
  IF(ALL(listtmp(1:2)==comp))THEN
    oflag=-i
    RETURN
  END IF
  listtmp=CSHIFT(listtmp,1)
END DO
END SUBROUTINE find_orient_listn
!---------------------------------------------------------------------------------
!> Apply an orientation transform to a n-value array
!---------------------------------------------------------------------------------
PURE SUBROUTINE orient_listn(oflag,list,n)
INTEGER(i4), INTENT(in) :: oflag !< Orientation flag
INTEGER(i4), INTENT(inout) :: list(n) !< Array for orientation [n]
INTEGER(i4), INTENT(in) :: n !< Length of array
INTEGER(i4) :: i
IF(oflag<0)list=list(n:1:-1) ! Reverse list
list=CSHIFT(list,ABS(oflag)-1) ! Shift list
END SUBROUTINE orient_listn
!---------------------------------------------------------------------------------
!> Apply an orientation transform to a n-value array
!---------------------------------------------------------------------------------
PURE SUBROUTINE orient_listn_inv(oflag,list,n)
INTEGER(i4), INTENT(in) :: oflag !< Orientation flag
INTEGER(i4), INTENT(inout) :: list(n) !< Array for orientation [n]
INTEGER(i4), INTENT(in) :: n !< Length of array
INTEGER(i4) :: i
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ltmp
ALLOCATE(ltmp(n))
ltmp=(/(i,i=1,n)/)
CALL orient_listn(oflag, ltmp, n)
list(ltmp)=list
DEALLOCATE(ltmp)
END SUBROUTINE orient_listn_inv
!---------------------------------------------------------------------------------
!> Compute coefficients for linear 1-D interpolation of function F(x)
!!
!! @warning This function requires `x` be sorted lowest to highest.
!! @note This function performs an interval search each time it is called.
!---------------------------------------------------------------------------------
SUBROUTINE linterp_facs(x,n,xx,inds,facs,extrap)
REAL(r8), INTENT(in) :: x(n) !< Paramaterizing array \f$ x_i \f$ [n]
REAL(r8), INTENT(in) :: xx !< Location to perform interpolation
INTEGER(i4), INTENT(in) :: n !< Length of function parameterization
INTEGER(i4), INTENT(out) :: inds(2) !< Indices of points bounding subinterval (-1 -> error)
REAL(r8), INTENT(out) :: facs(2) !< Interpolation factors for point in `inds(2)`
INTEGER(i4), OPTIONAL, INTENT(in) :: extrap !< Extrapolation mode (0: none, 1: constant, 2: linear)
INTEGER(i4) :: i
DO i=2,n
  IF(x(i-1)<=xx.AND.x(i)>=xx)EXIT
END DO
IF(i<=n)THEN
  inds=[i,i-1]
  facs=[1.d0,-1.d0]*(xx-x(i-1))/(x(i)-x(i-1)) + [0.d0,1.d0]
ELSE
  IF(PRESENT(extrap))THEN
    SELECT CASE(extrap)
    CASE(0)
      inds=[-1,-1]
      facs=[0.d0,0.d0]
    CASE(1)
      IF(xx<x(1))THEN
        inds=[1,1]
        facs=[1.d0,0.d0]
      ELSE IF(xx>x(n))THEN
        inds=[n,n]
        facs=[1.d0,0.d0]
      END IF
    CASE(2)
      IF(xx<x(1))THEN
        inds=[2,1]
        facs=[1.d0,-1.d0]*(xx-x(1))/(x(2)-x(1)) + [0.d0,1.d0]
      ELSE IF(xx>x(n))THEN
        inds=[n,n-1]
        facs=[1.d0,-1.d0]*(xx-x(n-1))/(x(n)-x(n-1)) + [0.d0,1.d0]
      END IF
    CASE DEFAULT
      inds=[-2,-2]
      facs=[0.d0,0.d0]
    END SELECT
  END IF
END IF
END SUBROUTINE linterp_facs
!---------------------------------------------------------------------------------
!> Perform linear 1-D interpolation of function F(x)
!!
!! @warning This function requires `x` be sorted lowest to highest.
!! @note This function performs an interval search each time it is called.
!!
!! @returns \f$ F(xx) \f$ (-1.E99 if outside domain and `extrap=0` or invalid value for `extrap`)
!---------------------------------------------------------------------------------
FUNCTION linterp(x,y,n,xx,extrap) result(yy)
REAL(r8), INTENT(in) :: x(n) !< Paramaterizing array \f$ x_i \f$ [n]
REAL(r8), INTENT(in) :: y(n) !< Function values \f$ F(x_i) \f$ [n]
REAL(r8), INTENT(in) :: xx !< Location to perform interpolation
INTEGER(i4), INTENT(in) :: n !< Length of function parameterization
INTEGER(i4), OPTIONAL, INTENT(in) :: extrap !< Extrapolation mode (0: none, 1: constant, 2: linear)
INTEGER(i4) :: inds(2)
REAL(r8) :: yy,facs(2)
CALL linterp_facs(x,n,xx,inds,facs,extrap)
IF(inds(1)>0)THEN
  yy=y(inds(1))*facs(1)+y(inds(2))*facs(2)
ELSE
  yy=-1.d99
END IF
END FUNCTION linterp
!------------------------------------------------------------------------------
!> Reset the stack
!!
!! @warning This should only be done for the main run program, otherwise the
!! stack may become corrupted.
!------------------------------------------------------------------------------
SUBROUTINE oft_stack_reset
#ifdef OFT_STACK
!$omp parallel
nstack = 0_i4
stack = 0_i4
!$omp end parallel
#endif
END SUBROUTINE oft_stack_reset
#ifdef OFT_STACK
!------------------------------------------------------------------------------
!> Add a subroutine call to the current stack
!------------------------------------------------------------------------------
SUBROUTINE oft_stack_push(mod_ind,sub_ind)
INTEGER(i4), INTENT(in) :: mod_ind !< Index of containing module
INTEGER(i4), INTENT(in) :: sub_ind !< Index of current subroutine
IF(oft_stack_disabled)RETURN
!---
#ifdef OFT_PROFILE
IF(oft_prof_enabled)THEN
  IF(nstack>0)THEN
    stack_fun_time(stack(2,nstack))=stack_fun_time(stack(2,nstack))+local_timer%int_tock()
  ELSE
    stack_fun_time(0)=stack_fun_time(0)+local_timer%int_tock()
  END IF
END IF
#endif
!---
nstack=nstack+1
IF(nstack>stacksize)nstack=stacksize
stack(:,nstack)=[mod_ind,sub_ind]
#ifdef OFT_PROFILE
IF(oft_prof_enabled)stack_nfun_c(sub_ind)=stack_nfun_c(sub_ind)+1
#endif
END SUBROUTINE oft_stack_push
!------------------------------------------------------------------------------
!> Pop a subroutine off of the current stack
!------------------------------------------------------------------------------
SUBROUTINE oft_stack_pop
IF(oft_stack_disabled)RETURN
!---
#ifdef OFT_PROFILE
IF(oft_prof_enabled)THEN
  IF(nstack>0)THEN
    stack_fun_time(stack(2,nstack))=stack_fun_time(stack(2,nstack))+local_timer%int_tock()
  ELSE
    stack_fun_time(0)=stack_fun_time(0)+local_timer%int_tock()
  END IF
END IF
#endif
!---
nstack=nstack-1
IF(nstack<0)nstack=0
END SUBROUTINE oft_stack_pop
#endif
!------------------------------------------------------------------------------
!> Print the current contents of the stack
!------------------------------------------------------------------------------
SUBROUTINE oft_stack_print() BIND(C)
#ifdef OFT_STACK
INTEGER(i4) :: i,outunit
#ifdef OFT_ABORT_FILES
CHARACTER(LEN=5) :: proc
#endif
outunit=error_unit
#ifdef OFT_ABORT_FILES
WRITE(proc,'(I5.5)')oft_env%rank
OPEN(outunit,FILE='stack_'//proc//'.err')
#endif
WRITE(outunit,*)'Dumping stack',nstack
IF(nstack==0)RETURN
!---
WRITE(outunit,'(A)')  ''
WRITE(outunit,'(A)')  'Stacktrace'
WRITE(outunit,'(A)')  '#----------------------------------------------'
100 FORMAT (A,I5,A,2X,I2,2X,A60)
101 FORMAT (2X,I2,2X,A60)
DO i=nstack,1,-1
  IF(oft_env%nprocs>0)THEN
    WRITE(outunit,100)  '[',oft_env%rank,']',i,ADJUSTR(TRIM(stack_mods(stack(1,i)))//"::"//TRIM(stack_funs(stack(2,i))))
  ELSE
    WRITE(outunit,101)  i,ADJUSTR(TRIM(stack_mods(stack(1,i)))//"::"//TRIM(stack_funs(stack(2,i))))
  END IF
END DO
WRITE(outunit,'(A)')  ''
#ifdef OFT_ABORT_FILES
CLOSE(outunit)
#endif
#endif
END SUBROUTINE oft_stack_print
!------------------------------------------------------------------------------
!> Activate profiling
!------------------------------------------------------------------------------
SUBROUTINE oft_profile_start()
#ifdef OFT_PROFILE
oft_prof_enabled=.TRUE.
#else
CALL oft_abort("OFT not built with profiling","oft_profile_start",__FILE__)
#endif
END SUBROUTINE oft_profile_start
!------------------------------------------------------------------------------
!> Reset all profiling counters
!------------------------------------------------------------------------------
SUBROUTINE oft_prof_reset
#ifdef OFT_PROFILE
!$omp parallel
stack_nfun_c=0_i4
stack_fun_time=0_i8
!$omp end parallel
#endif
END SUBROUTINE oft_prof_reset
!------------------------------------------------------------------------------
!> Print some basic profiling information
!------------------------------------------------------------------------------
SUBROUTINE oft_prof_print
#ifdef OFT_PROFILE
INTEGER(i4) :: i,j,k,ierr
INTEGER(i4), ALLOCATABLE :: ind(:)
INTEGER(i8), ALLOCATABLE :: lnfc(:),nfc(:)
INTEGER(i8), ALLOCATABLE :: lcft(:),cft(:)
REAL(r8) :: ntmp
INTEGER(i8) :: countnew,crate,cmax
CHARACTER(LEN=1) :: suffix
IF(.NOT.oft_prof_enabled)RETURN
call system_clock(countnew,crate,cmax)
!---
ALLOCATE(lnfc(stack_nfuns),lcft(0:stack_nfuns))
lnfc=0_i8
lcft=0_i8
!$omp parallel
!$omp atomic
lcft(0)=lcft(0)+stack_fun_time(0)
DO i=1,stack_nfuns
  !$omp atomic
  lnfc(i)=lnfc(i)+stack_nfun_c(i)
  !$omp atomic
  lcft(i)=lcft(i)+stack_fun_time(i)
END DO
!$omp end parallel
!---
IF(oft_env%rank==0)THEN
  ALLOCATE(nfc(stack_nfuns),cft(0:stack_nfuns))
ELSE
  ALLOCATE(nfc(1),cft(0:1))
END IF
nfc=oft_mpi_sum(lnfc, stack_nfuns)
cft=oft_mpi_sum(lcft, stack_nfuns+1)
! CALL MPI_Reduce(lnfc,nfc,stack_nfuns,OFT_MPI_I8,MPI_SUM,0,oft_env%comm,ierr)
! CALL MPI_Reduce(lcft,cft,stack_nfuns+1,OFT_MPI_I8,MPI_SUM,0,oft_env%comm,ierr)
IF(oft_env%rank==0)THEN
  !--- Sort by time used
  ALLOCATE(ind(0:stack_nfuns))
  ind=(/(i,i=0,stack_nfuns)/)
  ind(0)=1
  CALL sort_array(cft,ind,stack_nfuns+1)
  !---
  WRITE(*,'(A)') ""
  WRITE(*,'(A)') "Performance Summary:"
  WRITE(*,'(A)') " Rank  |        Module/Context        ::        Function/Region        |    NFC    |    CFT   "
  !                    1  tetmesh_mapping               ::             tetmesh_log_grad2    10239680         1.7
  WRITE(*,'(A)') "#----------------------------------------------------------------------------------------------"
  100 FORMAT (I6,2X,A40,A2,A40,2X,F9.2,A1,2X,ES12.2)
  DO i=0,20
    j=stack_nfuns-i
    k=ind(j)
    ntmp=REAL(nfc(k),8)
    suffix=" "
    IF(nfc(k)/1E3<=1._r8)THEN
      suffix=" "
    ELSE IF(nfc(k)/1E3>1._r8.AND.nfc(k)/1E6<=1._r8)THEN
      suffix="K"
      ntmp=REAL(nfc(k)/1E3,8)
    ELSE IF(nfc(k)/1E6>1._r8.AND.nfc(k)/1E9<=1._r8)THEN
      suffix="M"
      ntmp=REAL(nfc(k)/1E6,8)
    ELSE IF(nfc(k)/1E9>1._r8.AND.nfc(k)/1E12<=1._r8)THEN
      suffix="B"
      ntmp=REAL(nfc(k)/1E9,8)
    END IF
    WRITE(*,100)  i+1,ADJUSTR(stack_mods(stack_fun_mods(k))), &
      "::",ADJUSTL(stack_funs(k)),ntmp,suffix,REAL(cft(j),8)/crate
  END DO
  WRITE(*,'(A)')  ""
  DEALLOCATE(ind)
END IF
CALL oft_mpi_barrier(ierr)
DEALLOCATE(nfc,lnfc,cft,lcft)
#endif
END SUBROUTINE oft_prof_print
END MODULE oft_base
