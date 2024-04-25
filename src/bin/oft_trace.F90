!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_trace.F90
!
!> Driver to trace field lines and particles from fields stored in a restart file
!!
!! **Option group:** `oft_trace_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `order=1`              |  FE order     | int |
!! |  `type=1`               |  Field type (1 -> vLag, 2-> H1, 3->H1(Curl))  | int |
!! |  `fields="","",""`      |  Sub-field names in restart files  | str(10) [3] |
!! |  `rst_file="none"`      |  Restart file containing fields  | str(40) |
!! |  `pt_file="none"`       |  File containing launch point list  | str(40) |
!! |  `adv_type=1`           |  Type of advance to be used (1=>b-field, 2=>lorentz)  | int |
!! |  `tracer_type=1`        |  Type of tracer to be used  | int |
!! |  `bscale=1`             |  Scale factor for magnetic field for `adv_type=2`  | real |
!! |  `mu_ion=1`             |  Ion mass in atomic units for `adv_type=2`  | real |
!! |  `tracer_maxsteps=1E5`  |  Maximum number of tracer steps per field line  | int |
!! |  `tracer_maxtrans=1E4`  |  Maximum number of domain transfers per field line  | int |
!! |  `tracer_tol=1E-9`      |  Tracer tolerance  | real |
!! |  `tracer_timeout=60`    |  Timeout for each Poincare section  | real |
!!
!! @authors Chris Hansen
!! @date June 2018
!! @ingroup doxy_oft_bin
!-----------------------------------------------------------------------------
PROGRAM oft_trace
!---Base
USE ISO_FORTRAN_ENV, ONLY: IOSTAT_END
USE oft_base
!--Grid
USE multigrid_build, ONLY: multigrid_construct
USE oft_io, ONLY: oft_file_exist
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lagrange, oft_lag_setup
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
USE mhd_utils, ONLY: elec_charge, proton_mass
USE tracing, ONLY: oft_tracer, create_tracer, tracing_line, set_timeout
IMPLICIT NONE
#include "local.h"
!---Local variables
INTEGER(i4) :: i,j,ind,npts,io_stat,ierr,pt_file_unit,nfinal,io_unit
REAL(r8) :: f(4)
REAL(r8), POINTER, DIMENSION(:) :: valtmp => NULL()
REAL(r8), ALLOCATABLE, DIMENSION(:) :: pt
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=OFT_PATH_SLEN) :: filename
CLASS(oft_vector), POINTER :: u => NULL()
CLASS(oft_vector), POINTER :: x1 => NULL()
CLASS(oft_vector), POINTER :: x2 => NULL()
TYPE(oft_lag_vrinterp), TARGET :: Bfield_lag
TYPE(oft_h1_rinterp), TARGET :: Bfield_H1
TYPE(oft_hcurl_cinterp), TARGET :: Bfield_HCurl
CLASS(oft_tracer), POINTER :: tracer
REAL(r8), PARAMETER :: vel_scale = 1.d3
!---Input options
INTEGER(i4) :: order = 2
INTEGER(i4) :: type = 1
INTEGER(i4) :: tracer_type = 1
INTEGER(i4) :: adv_type = 1
REAL(r8) :: tracer_maxsteps = 1.d5
REAL(r8) :: tracer_maxtrans = 1.d4
REAL(r8) :: tracer_tol = 1.d-9
REAL(r8) :: tracer_timeout = 60.d0
REAL(r8) :: bscale = 1.d0
REAL(r8) :: mu_ion = 1.d0
CHARACTER(LEN=10) :: fields(3) = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: pt_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_file = 'none'
NAMELIST/oft_trace_options/order,type,fields,bscale,mu_ion,pt_file,rst_file, &
  adv_type,tracer_type,tracer_maxsteps,tracer_maxtrans,tracer_tol,tracer_timeout
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,oft_trace_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check inputs
IF(ierr<0)CALL oft_abort('No "oft_trace_options" found in input file.', &
  'oft_trace',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "oft_trace_options" in input file.', &
  'oft_trace',__FILE__)
!---Setup grid
CALL multigrid_construct
!---Setup tracer
CALL set_timeout(tracer_timeout)
CALL create_tracer(tracer,tracer_type)
tracer%tol=tracer_tol
tracer%maxsteps=INT(tracer_maxsteps)
tracer%maxtrans=INT(tracer_maxtrans)
SELECT CASE(adv_type)
  CASE(1)
    ALLOCATE(pt(3))
  CASE(2)
    ALLOCATE(pt(6))
    tracer%neq=6
    tracer%ydot=>particle_lorentz
  CASE DEFAULT
    CALL oft_abort("Unknown advance type", "oft_trace", __FILE__)
END SELECT
!---Setup restart file
IF(.NOT.oft_file_exist(TRIM(rst_file)))THEN
  CALL oft_abort('Restart file does not exist','oft_trace',__FILE__)
END IF
!---Setup necessary FE space
SELECT CASE(type)
  CASE(1) ! Vector Lagrange field
    CALL oft_lag_setup(order, -1)
    !---Create field structure
    CALL oft_lag_create(x1)
    CALL oft_lag_vcreate(u)
    !---Load field components
    DO i=1,3
      CALL oft_lagrange%vec_load(x1,rst_file,fields(i))
      CALL x1%get_local(valtmp)
      CALL u%restore_local(valtmp,i)
    END DO
    !---Setup interplolation
    CALL u%scale(bscale)
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
    !---Load H1(Curl) component
    CALL oft_hcurl%vec_load(x1,rst_file,fields(1))
    CALL x1%get_local(valtmp)
    CALL u%restore_local(valtmp,1)
    !---Load H1(Grad) component
    CALL oft_h0%vec_load(x2,rst_file,fields(2))
    CALL x2%get_local(valtmp)
    CALL u%restore_local(valtmp,2)
    !---Setup interplolation
    CALL u%scale(bscale)
    Bfield_H1%u=>u
    tracer%B=>Bfield_H1
  CASE(3) !  Nedelec HCurl field
    CALL oft_hcurl_setup(order, -1)
    !---Create field structure
    CALL oft_hcurl_create(u)
    !---Load H1(Curl) field
    CALL oft_hcurl%vec_load(u,rst_file,fields(1))
    !---Setup interplolation
    CALL u%scale(bscale)
    Bfield_HCurl%u=>u
    tracer%B=>Bfield_HCurl
  CASE DEFAULT
    CALL oft_abort("Unknown field type", "oft_trace", __FILE__)
END SELECT
CALL tracer%B%setup()
!---Loop over launch points
ind=0
OPEN(newunit=pt_file_unit,FILE=TRIM(pt_file),IOSTAT=io_stat)
IF(io_stat/=0)CALL oft_abort("Error opening pt_file", "oft_trace", __FILE__)
READ(pt_file_unit,*,IOSTAT=io_stat)npts
DO i=1,npts
  !---Read next launch point
  READ(pt_file_unit,*,IOSTAT=io_stat)pt
  IF(io_stat/=0)THEN
    IF(oft_env%head_proc.AND.io_stat/=IOSTAT_END)THEN
      WRITE(*,*)'Error reading launch point at line ',INT(i+1,2)
    END IF
    EXIT
  END IF
  IF(adv_type==2)pt(4:6)=pt(4:6)/vel_scale
  !---Perform tracing
  WRITE(pltnum,'(I4.4)')i
  CALL tracing_line(tracer,pt,'trace_'//pltnum)
END DO
CLOSE(pt_file_unit)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: particle_lorentz
!---------------------------------------------------------------------------
!> ODE function for Lorentz force particle advance
!!
!! F = q * cross(V,B) / m_i
!---------------------------------------------------------------------------
SUBROUTINE particle_lorentz(t,y,B,n,ydot)
REAL(r8), INTENT(in) :: t
REAL(r8), INTENT(in) :: y(n)
REAL(r8), INTENT(in) :: B(3)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(out) :: ydot(n)
ydot(1:3)=y(4:6)*vel_scale
ydot(4:6)=elec_charge*cross_product(y(4:6),B)/(proton_mass*mu_ion)
END SUBROUTINE particle_lorentz
END PROGRAM oft_trace
