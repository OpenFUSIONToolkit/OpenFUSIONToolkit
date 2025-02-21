!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_td.F90
!
!> @defgroup doxy_thincurr ThinCurr
!! ThinCurr drivers for 3D thin-wall modeling
!! @ingroup doxy_oft_bin
!
!> Run time-dependent thin wall simulations using ThinCurr
!!
!! **Option group:** `thincurr_td_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `mesh_file="none"`     |  Surface mesh filename (Cubit) | str |
!! |  `dt=1.E-4`             |  Time step for time dependent run | float |
!! |  `nsteps=400`           |  Number of steps for time dependent run | int |
!! |  `timestep_cn=T`        |  Use Crank-Nicolson timestep | bool |
!! |  `cg_tol=1.E-6`         |  Convergence tolerance for `direct=F` | float |
!! |  `nplot=10`             |  Restart save frequency for time dependent run | int |
!! |  `plot_run=F`           |  Produce plot files from stored restart files | bool |
!! |  `direct=T`             |  Use direct solver | bool |
!! |  `save_mat=T`           |  Store inverted matrix for later use | bool |
!! |  `compute_b=F`          |  Compute magnetic fields on cell centers | bool |
!!
!! @authors Chris Hansen
!! @date Feb 2022
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_td
USE oft_base
USE oft_io, ONLY: oft_bin_file
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_lu, ONLY: lapack_cholesky
USE oft_native_la, ONLY: oft_native_dense_matrix
USE oft_deriv_matrices, ONLY: oft_sum_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE mhd_utils, ONLY: mu0
USE thin_wall
USE thin_wall_hodlr
USE thin_wall_solvers, ONLY: run_td_sim, plot_td_sim
IMPLICIT NONE
#include "local.h"
TYPE(tw_type), TARGET :: tw_sim
TYPE(tw_sensors) :: sensors
!
LOGICAL :: exists
INTEGER(4) :: i,j,k,n,ierr,io_unit,ncols,ntimes
REAL(8), ALLOCATABLE :: curr_ic(:)
REAL(8), POINTER, DIMENSION(:,:) :: curr_waveform,volt_waveform,sensor_waveform
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
TYPE(oft_tw_hodlr_op), TARGET :: tw_hodlr
!
REAL(8) :: dt = 1.d-4
REAL(8) :: cg_tol=1.d-6
INTEGER(4) :: nsteps = 400
INTEGER(4) :: nstatus = 10
INTEGER(4) :: nplot = 10
INTEGER(4) :: jumper_start = 0
LOGICAL :: timestep_cn=.TRUE.
LOGICAL :: direct = .TRUE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: compute_B = .FALSE.
LOGICAL :: plot_rebuild_sensors = .FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: curr_file="none"
CHARACTER(LEN=OFT_PATH_SLEN) :: volt_file="none"
NAMELIST/thincurr_td_options/curr_file,volt_file,dt,nsteps,nstatus,nplot,direct,save_L,save_Mcoil,save_Msen,  &
  plot_run,compute_B,plot_rebuild_sensors,timestep_cn,cg_tol,jumper_start
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_td_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)call oft_abort('No thin-wall options found in input file.','thincurr_td',__FILE__)
if(ierr>0)call oft_abort('Error parsing thin-wall options in input file.','thincurr_td',__FILE__)
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
  CALL oft_abort("Unsupported mesh type","thincurr_td",__FILE__)
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
IF((TRIM(curr_file)=="none").AND.(tw_sim%n_icoils>0))CALL oft_abort("No waveform filename specified", &
  "thincurr_td",__FILE__)
!---Setup I/0
CALL tw_sim%xdmf%setup("thincurr")
CALL smesh%setup_io(tw_sim%xdmf,1)
IF(oft_debug_print(1))CALL tw_sim%save_debug()
!---------------------------------------------------------------------------
! Time-dependent run
!---------------------------------------------------------------------------
!---Load drivers and sensors
CALL tw_load_sensors('floops.loc',tw_sim,sensors)
!---Compute inductances
WRITE(*,*)
IF(.NOT.plot_run)THEN
  IF(save_Mcoil)THEN
    CALL tw_compute_Ael2dr(tw_sim,'Mcoil.save')
  ELSE
    CALL tw_compute_Ael2dr(tw_sim)
  END IF
END IF
IF(save_Msen)THEN
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops,'Msen.save')
ELSE
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
END IF
!---------------------------------------------------------------------------
! Load or build element to element mutual matrix
!---------------------------------------------------------------------------
tw_hodlr%tw_obj=>tw_sim
CALL tw_hodlr%setup(.FALSE.)
IF(.NOT.plot_run)THEN
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    IF(direct)CALL oft_abort('HODLR compression does not support "direct=T"','thincurr_td',__FILE__)
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
!---------------------------------------------------------------------------
! Run main calculation or plots
!---------------------------------------------------------------------------
IF(plot_run)THEN
  NULLIFY(sensor_waveform)
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    CALL plot_td_sim(tw_sim,nsteps,nplot,sensors,compute_B,plot_rebuild_sensors,sensor_waveform,tw_hodlr)
  ELSE
    CALL plot_td_sim(tw_sim,nsteps,nplot,sensors,compute_B,plot_rebuild_sensors,sensor_waveform)
  END IF
ELSE
  !---Setup resistivity matrix
  CALL tw_compute_Rmat(tw_sim,.FALSE.)
  !---Load I-coil waveform
  IF(TRIM(curr_file)/="none")THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(curr_file))
    READ(io_unit,*)ncols,ntimes
    ALLOCATE(curr_waveform(ntimes,ncols))
    DO i=1,ntimes
      READ(io_unit,*)curr_waveform(i,:)
    END DO
    CLOSE(io_unit)
  ELSE
    NULLIFY(curr_waveform)
  END IF
  !---Load V-coil waveform
  IF(TRIM(volt_file)/="none")THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(volt_file))
    READ(io_unit,*)ncols,ntimes
    ALLOCATE(volt_waveform(ntimes,ncols))
    DO i=1,ntimes
      READ(io_unit,*)volt_waveform(i,:)
    END DO
    CLOSE(io_unit)
  ELSE
    NULLIFY(volt_waveform)
  END IF
  NULLIFY(sensor_waveform)
  !---Run time-dependent simulation
  ALLOCATE(curr_ic(tw_sim%nelems))
  curr_ic=0.d0
  oft_env%pm=.FALSE.
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    CALL run_td_sim(tw_sim,dt,nsteps,curr_ic,direct,cg_tol,timestep_cn,nstatus, &
      nplot,sensors,curr_waveform,volt_waveform,sensor_waveform,tw_hodlr)
  ELSE
    CALL run_td_sim(tw_sim,dt,nsteps,curr_ic,direct,cg_tol,timestep_cn,nstatus, &
      nplot,sensors,curr_waveform,volt_waveform,sensor_waveform)
  END IF
END IF
!---
CALL oft_finalize
END PROGRAM thincurr_td
