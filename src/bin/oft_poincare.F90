!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file poincare_trace.F90
!
!> Driver to compute poincare section from fields stored in restart files
!!
!! **Option group:** `poincare_trace_options`
!! |  Option  |  Description  | Type [dim] |
!! |------------|------------|-----|
!! |  `order=1`              |  FE order  | int |
!! |  `type=1`               |  Field type (1 -> vLag, 2-> H1, 3->H1(Curl))  | int |
!! |  `fields=""`            |  Sub-field names in restart files  | str(10) [3] |
!! |  `pt_file="none"`       |  File containing launch point list  | str(40) |
!! |  `rst_list_file="none"` |  File containing restart file list  | str(40) |
!! |  `tracer_type=1`        |  Type of tracer to be used  | int |
!! |  `tracer_maxsteps=1E5`  |  Maximum number of tracer steps per field line  | int |
!! |  `tracer_maxtrans=1E4`  |  Maximum number of domain transfers per field line  | int |
!! |  `tracer_tol=1E-9`      |  Tracer tolerance  | real |
!! |  `tracer_timeout=60`    |  Timeout for each Poincare section  | real |
!! |  `compute_q=F`          |  Compute approximate safety factor  | bool |
!!
!! @authors Chris Hansen
!! @date December 2014
!! @ingroup doxy_oft_bin
!-----------------------------------------------------------------------------
PROGRAM poincare_trace
!---Base
USE ISO_FORTRAN_ENV, ONLY: IOSTAT_END
USE oft_base
!--Grid
USE oft_mesh_type, ONLY: mesh_findcell, mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_io, ONLY: oft_file_exist
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lagrange ,oft_lag_setup
USE oft_lag_fields, ONLY: oft_lag_create, oft_lag_vcreate!, oft_lag_load
USE oft_lag_operators, ONLY: oft_lag_vrinterp
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl, oft_hcurl_setup
USE oft_hcurl_fields, ONLY: oft_hcurl_create!, oft_hcurl_load
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
  hcurl_mloptions
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0, oft_h0_setup
USE oft_h0_fields, ONLY: oft_h0_create!, oft_h0_load
USE oft_h0_operators, ONLY: h0_mloptions, h0_setup_interp
!---H1 Full FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_rinterp
!---Tracing
USE tracing, ONLY: oft_tracer, create_tracer, tracing_poincare, set_timeout
IMPLICIT NONE
#include "local.h"
!---Local variables
INTEGER(i4) :: i,j,ind,npts,io_stat,ierr,rst_file_unit,nfinal,io_unit
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pt_flag
REAL(r8) :: f(4)
REAL(r8), POINTER, DIMENSION(:) :: valtmp => NULL()
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: pts,pts_final
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=OFT_PATH_SLEN) :: filename
CLASS(oft_vector), POINTER :: u => NULL()
CLASS(oft_vector), POINTER :: x1 => NULL()
CLASS(oft_vector), POINTER :: x2 => NULL()
TYPE(oft_lag_vrinterp), TARGET :: Bfield_lag
TYPE(oft_h1_rinterp), TARGET :: Bfield_H1
TYPE(oft_hcurl_cinterp), TARGET :: Bfield_HCurl
CLASS(oft_tracer), POINTER :: tracer
!---Input options
INTEGER(i4) :: order = 2
INTEGER(i4) :: type = 1
INTEGER(i4) :: tracer_type = 1
REAL(r8) :: tracer_maxsteps = 1.d5
REAL(r8) :: tracer_maxtrans = 1.d4
REAL(r8) :: tracer_tol = 1.d-9
REAL(r8) :: tracer_timeout = 60.d0
LOGICAL :: compute_q = .FALSE.
CHARACTER(LEN=10) :: fields(3) = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: pt_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_list_file = 'none'
NAMELIST/poincare_trace_options/order,type,fields,pt_file,rst_list_file, &
  tracer_type,tracer_maxsteps,tracer_maxtrans,tracer_tol,tracer_timeout, &
  compute_q
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,poincare_trace_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check inputs
IF(ierr<0)CALL oft_abort('No "poincare_trace_options" found in input file.', &
  'poincare_trace',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "poincare_trace_options" in input file.', &
  'poincare_trace',__FILE__)
IF(TRIM(pt_file)=='none')CALL oft_abort('No launch point file specified.',  &
  'poincare_trace',__FILE__)
!---Setup grid
CALL multigrid_construct
!---Setup tracer
CALL set_timeout(tracer_timeout)
CALL create_tracer(tracer,tracer_type)
tracer%tol=tracer_tol
tracer%maxsteps=INT(tracer_maxsteps)
tracer%maxtrans=INT(tracer_maxtrans)
!---Setup necessary FE space
SELECT CASE(type)
  CASE(1) ! Vector Lagrange field
    CALL oft_lag_setup(order, -1)
    !---Create field structure
    CALL oft_lag_create(x1)
    CALL oft_lag_vcreate(u)
    Bfield_lag%u=>u
    tracer%B=>Bfield_lag
  CASE(2) !  Nedelec H1 field
    CALL oft_hcurl_setup(order, -1)
    CALL oft_h0_setup(order+1, -1)
    CALL oft_h1_setup(order, -1)
    !---Create field structure
    CALL oft_hcurl_create(x1)
    CALL oft_h0_create(x2)
    CALL oft_h1_create(u)
    Bfield_H1%u=>u
    tracer%B=>Bfield_H1
  CASE(3) !  Nedelec HCurl field
    CALL oft_hcurl_setup(order, -1)
    !---Create field structure
    CALL oft_hcurl_create(u)
    Bfield_HCurl%u=>u
    tracer%B=>Bfield_HCurl
  CASE DEFAULT
    CALL oft_abort("Unknown field type", "poincare_trace", __FILE__)
END SELECT
!---Load launch point list
OPEN(NEWUNIT=io_unit,FILE=TRIM(pt_file))
READ(io_unit,*)npts
ALLOCATE(pts(3,npts))
DO i=1,npts
  READ(io_unit,*)pts(:,i)
END DO
CLOSE(io_unit)
!---Filter input points
ALLOCATE(pt_flag(npts))
nfinal=0
pt_flag=1
!$omp parallel do private(j,f)
DO i=1,npts
  j=0
  call mesh_findcell(mesh,j,pts(:,i),f)
  IF(( MAXVAL(f)>=1.d0+1.d-3 ).OR.( MINVAL(f)<=-1.d-3 ))pt_flag(i)=0
END DO
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,pt_flag,npts,OFT_MPI_I4,MPI_MAX,oft_env%COMM,ierr)
#endif
nfinal=SUM(pt_flag)
ALLOCATE(pts_final(3,nfinal))
nfinal=0
DO i=1,npts
  IF(pt_flag(i)==1)THEN
    nfinal=nfinal+1
    pts_final(:,nfinal)=pts(:,i)
  END IF
END DO
DEALLOCATE(pts,pt_flag)
IF(nfinal==0)CALL oft_abort('No launch points within mesh!','poincare_trace',__FILE__)
IF(oft_env%head_proc.AND.npts/=nfinal)THEN
  WRITE(*,*)'Removed ',INT(npts-nfinal,2),' offmesh points'
  WRITE(*,*)'Proceeding with remaining ',INT(nfinal,2),' points'
END IF
!---Loop over file list
ind=0
OPEN(newunit=rst_file_unit,FILE=TRIM(rst_list_file),IOSTAT=io_stat)
IF(io_stat/=0)CALL oft_abort("Error loading rst_file", "poincare_trace", __FILE__)
DO
  !---Read next filename
  READ(rst_file_unit,'(A)',IOSTAT=io_stat)filename
  IF(io_stat/=0)THEN
    IF(oft_env%head_proc.AND.io_stat/=IOSTAT_END)THEN
      WRITE(*,*)'Error reading rst_list_file at line ',INT(ind+1,2)
    END IF
    EXIT
  END IF
  IF(.NOT.oft_file_exist(filename))THEN
    IF(oft_env%head_proc)WRITE(*,*)'Restart file at line',INT(ind+1,2),' in rst_list_file does not exist'
    EXIT
  END IF
  !---Read field
  SELECT CASE(type)
    CASE(1) ! Vector Lagrange field
      !---Extract field components
      DO i=1,3
        CALL oft_lagrange%vec_load(x1,filename,fields(i))
        CALL x1%get_local(valtmp)
        CALL u%restore_local(valtmp,i)
      END DO
    CASE(2) !  Nedelec H1 field
      !---Extract H1(Curl) component
      CALL oft_hcurl%vec_load(x1,filename,fields(1))
      CALL x1%get_local(valtmp)
      CALL u%restore_local(valtmp,1)
      !---Extract H1(Grad) component
      CALL oft_h0%vec_load(x2,filename,fields(2))
      CALL x2%get_local(valtmp)
      CALL u%restore_local(valtmp,2)
    CASE(3) !  Nedelec HCurl field
      CALL oft_hcurl%vec_load(u,filename,fields(1))
  END SELECT
  CALL tracer%B%setup
  !---Perform tracing
  ind=ind+1
  WRITE(pltnum,'(I4.4)')ind
  IF(compute_q)THEN
    CALL tracing_poincare(tracer,pts_final,nfinal,'poin_'//pltnum//'.dat',qfile='q_'//pltnum//'.dat')
  ELSE
    CALL tracing_poincare(tracer,pts_final,nfinal,'poin_'//pltnum//'.dat')
  END IF
END DO
CLOSE(rst_file_unit)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM poincare_trace
