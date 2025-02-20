!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_fr.F90
!
!> Run thin wall frequency response simulations using ThinCurr
!!
!! **Option group:** `thincurr_fr_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `mesh_file="none"`     |  Surface mesh filename (Cubit) | str |
!! |  `plot_run=F`           |  Produce plot files from stored restart files | bool |
!! |  `direct=T`             |  Use direct solver | bool |
!! |  `force_f0=0.`          |  Toroidal field (R0*B0) for force calculation | float |
!! |  `freq=1.d3`            |  Frequency for FR run | float |
!! |  `mode_file="none"`     |  DCON mode surface file from "mode_to_tw" | str |
!! |  `fr_limit=0`           |  Frequency limit for FR run (1->Inf, 2->Zero) | str |
!!
!! @authors Chris Hansen
!! @date Feb 2022
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_fr
USE oft_base
USE oft_io, ONLY: oft_bin_file
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: mu0
USE thin_wall
USE thin_wall_hodlr
USE thin_wall_solvers, ONLY: frequency_response
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: nsensors = 0
TYPE(tw_type), TARGET :: tw_sim,mode_source
TYPE(tw_sensors) :: sensors
!
INTEGER(4) :: i,n,ierr,io_unit
REAL(8), POINTER, contiguous, DIMENSION(:,:) :: mode_driver,fr_driver,fr_sensor,driver_tmp
CHARACTER(LEN=2) :: eig_tag
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
TYPE(oft_tw_hodlr_op), TARGET :: tw_hodlr
!
INTEGER(4) :: fr_limit = 0
INTEGER(4) :: jumper_start = 0
LOGICAL :: direct = .TRUE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
REAL(8) :: freq = 1.d3
REAL(8) :: force_f0 = 0.d0
CHARACTER(LEN=OFT_PATH_SLEN) :: mode_file = 'none'
NAMELIST/thincurr_fr_options/direct,save_L,save_Mcoil,save_Msen,freq,mode_file, &
  fr_limit,jumper_start
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_fr_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)call oft_abort('No thin-wall options found in input file.', &
  'thincurr_fr',__FILE__)
if(ierr>0)call oft_abort('Error parsing thin-wall options in input file.', &
  'thincurr_fr',__FILE__)
!---Setup mesh
CALL multigrid_construct_surf
! ALLOCATE(mg_mesh)
! mg_mesh%mgmax=1
! mg_mesh%nbase=1
! oft_env%nbase=1
! mg_mesh%mgdim=mg_mesh%mgmax
! CALL smesh_cubit_load
SELECT CASE(smesh%cad_type)
CASE(0)
  CALL native_read_nodesets(mesh_nsets)
  CALL native_read_sidesets(mesh_ssets)
CASE(2)
#ifdef HAVE_NCDF
  CALL cubit_read_nodesets(mesh_nsets)
  CALL cubit_read_sidesets(mesh_ssets)
#endif
CASE DEFAULT
  CALL oft_abort("Unsupported mesh type","thincurr_fr",__FILE__)
END SELECT
IF(ASSOCIATED(mesh_ssets))THEN
  IF(mesh_ssets(1)%n>0)THEN
    tw_sim%nclosures=mesh_ssets(1)%n
    ALLOCATE(tw_sim%closures(tw_sim%nclosures))
    tw_sim%closures=mesh_ssets(1)%v
  END IF
END IF
tw_sim%mesh=>smesh
IF(jumper_start>0)THEN
  n=SIZE(mesh_nsets)
  hole_nsets=>mesh_nsets(1:jumper_start-1)
  ALLOCATE(tw_sim%jumper_nsets(n-jumper_start+1))
  DO i=jumper_start,n
    tw_sim%jumper_nsets(i-jumper_start+1)%n=mesh_nsets(i)%n
    ALLOCATE(tw_sim%jumper_nsets(i-jumper_start+1)%v(tw_sim%jumper_nsets(i-jumper_start+1)%n))
    tw_sim%jumper_nsets(i-jumper_start+1)%v=mesh_nsets(i)%v
  END DO
ELSE
  hole_nsets=>mesh_nsets
END IF
CALL tw_sim%setup(hole_nsets)
!---Setup I/0
CALL tw_sim%xdmf%setup("ThinCurr")
CALL smesh%setup_io(tw_sim%xdmf,1)
IF(oft_debug_print(1))CALL tw_sim%save_debug()
!---------------------------------------------------------------------------
! Frequency-response run
!---------------------------------------------------------------------------
!---Load sensors
CALL tw_load_sensors('floops.loc',tw_sim,sensors)
!---Setup voltage source
ALLOCATE(fr_driver(tw_sim%nelems,2),fr_sensor(sensors%nfloops,2))
fr_driver=0.d0
fr_sensor=0.d0
IF(TRIM(mode_file)/="none")THEN ! DCON mode source
  !---Load DCON mesh
  CALL tw_load_mode(mode_file,mode_source,mode_driver)
  !---Compute voltage source term from DCON mode
  CALL tw_compute_Lmat_MF(mode_source,tw_sim,2,mode_driver,fr_driver)
  !---Compute mode-sensor coupling
  CALL tw_compute_mutuals(mode_source,sensors%nfloops,sensors%floops)
  CALL dgemv('N',sensors%nfloops,mode_source%nelems,1.d0,mode_source%Ael2sen, &
    sensors%nfloops,mode_driver(:,1),1,0.d0,fr_sensor(:,1),1)
  CALL dgemv('N',sensors%nfloops,mode_source%nelems,1.d0,mode_source%Ael2sen, &
    sensors%nfloops,mode_driver(:,2),1,0.d0,fr_sensor(:,2),1)
  DEALLOCATE(mode_source%Adr2sen,mode_source%Ael2sen,mode_driver)
  !---Compute sensor mutuals for main model
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
ELSE ! Driver coil as source
  IF(tw_sim%n_icoils==0)CALL oft_abort('No coils specified for FR run', &
    'thincurr_fr', __FILE__)
  IF(tw_sim%n_icoils>2)CALL oft_abort('More than two coils specified for FR run', &
    'thincurr_fr', __FILE__)
  !---------------------------------------------------------------------------
  ! Load or build coil to element mutual matrix
  !---------------------------------------------------------------------------
  IF(save_Mcoil)THEN
    CALL tw_compute_Ael2dr(tw_sim,'Mcoil.save')
  ELSE
    CALL tw_compute_Ael2dr(tw_sim)
  END IF
  !---------------------------------------------------------------------------
  ! Load or build sensor mutual matrices
  !---------------------------------------------------------------------------
  IF(save_Msen)THEN
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops,'Msen.save')
  ELSE
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  END IF
  !---FR with coil 1
  fr_driver(:,1)=tw_sim%Ael2dr(:,1)/mu0
  fr_sensor(:,1)=tw_sim%Adr2sen(:,1)/mu0
  IF(tw_sim%n_icoils>1)THEN ! Second coil if present
    fr_driver(:,2)=tw_sim%Ael2dr(:,2)/mu0
    fr_sensor(:,2)=tw_sim%Adr2sen(:,2)/mu0
  END IF
END IF
!---------------------------------------------------------------------------
! Load or build element to element mutual matrix
!---------------------------------------------------------------------------
IF(fr_limit/=2)THEN
  tw_hodlr%tw_obj=>tw_sim
  CALL tw_hodlr%setup(.FALSE.)
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    IF(direct)CALL oft_abort('HODLR compression does not support "direct=T"','thincurr_fr',__FILE__)
    IF(save_L)THEN
      CALL tw_hodlr%compute_L(save_file='Lmat.save')
    ELSE
      CALL tw_hodlr%compute_L()
    END IF
  ELSE
    IF(save_L)THEN
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat,save_file='Lmat.save')
    ELSE
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat)
    END IF
  END IF
END IF
!---Setup resistivity matrix
CALL tw_compute_Rmat(tw_sim,.TRUE.)
!---Assume fixed current V = -i*omega*M*I
ALLOCATE(driver_tmp(tw_sim%nelems,1))
driver_tmp(:,1)=-fr_driver(:,1)
fr_driver(:,1)=fr_driver(:,2)
fr_driver(:,2)=driver_tmp(:,1)
IF(fr_limit==0)fr_driver=fr_driver*freq*2.d0*pi ! Scale driver by frequency
DEALLOCATE(driver_tmp)
!---Compute Frequency-response
IF(tw_hodlr%L_svd_tol>0.d0)THEN
  CALL frequency_response(tw_sim,direct,fr_limit,freq,fr_driver,hodlr_op=tw_hodlr)
ELSE
  CALL frequency_response(tw_sim,direct,fr_limit,freq,fr_driver)
END IF
CALL save_results
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE save_results()
INTEGER(i4) :: i,j,k,info
REAL(r8), ALLOCATABLE :: senin(:,:)
CLASS(oft_vector), POINTER :: uloc
TYPE(oft_bin_file) :: floop_hist
LOGICAL :: pm_save
!---Setup local vectors for I/O
CALL tw_sim%Uloc%new(uloc)
!---Save real part
CALL uloc%restore_local(fr_driver(:,1))
CALL tw_save_pfield(tw_sim,fr_driver(:,1),'JRe')
CALL tw_rst_save(tw_sim,uloc,'pThinCurr_frRe.rst','U')
!---Save imaginary part
CALL uloc%restore_local(fr_driver(:,2))
CALL tw_save_pfield(tw_sim,fr_driver(:,2),'JIm')
CALL tw_rst_save(tw_sim,uloc,'pThinCurr_frIm.rst','U')
!---Save probe signals
IF(sensors%nfloops>0)THEN
    ALLOCATE(senin(sensors%nfloops,2))
    senin=fr_sensor
    CALL dgemv('N',sensors%nfloops,tw_sim%nelems,1.d0,tw_sim%Ael2sen,sensors%nfloops,fr_driver(:,1),1,1.d0,fr_sensor(:,1),1)
    CALL dgemv('N',sensors%nfloops,tw_sim%nelems,1.d0,tw_sim%Ael2sen,sensors%nfloops,fr_driver(:,2),1,1.d0,fr_sensor(:,2),1)
    !---Setup history file
    IF(oft_env%head_proc)THEN
    floop_hist%filedesc = 'ThinCurr frequency-response flux loop signals (Re, Im, Re_vac, Im_vac)'
    CALL floop_hist%setup('thincurr_fr.hist')
    DO i=1,sensors%nfloops
        CALL floop_hist%add_field(sensors%floops(i)%name, 'r8')
    END DO
    CALL floop_hist%write_header
    CALL floop_hist%open
    CALL floop_hist%write(data_r8=fr_sensor(:,1))
    CALL floop_hist%write(data_r8=fr_sensor(:,2))
    CALL floop_hist%write(data_r8=senin(:,1))
    CALL floop_hist%write(data_r8=senin(:,2))
    CALL floop_hist%close
    END IF
    DEALLOCATE(senin)
END IF
!---Cleanup
CALL uloc%delete
DEALLOCATE(uloc)
END SUBROUTINE save_results
END PROGRAM thincurr_fr
