!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file tokamaker_fit.F90
!
!> Driver program for GS fitting
!!
!! @authors Chris Hansen
!! @date March 2014
!! @ingroup doxy_tokamaker
!---------------------------------------------------------------------------
PROGRAM tokamaker_fit
USE oft_base
USE oft_mesh_type, ONLY: smesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE fem_base, ONLY: oft_afem_type
USE oft_la_base, ONLY: oft_vector
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, oft_blagrange, &
  oft_lag_setup
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_eq, gs_save_fields, gs_save_fgrid, gs_setup_walls, gs_save_prof, &
  gs_fixed_vflux, gs_load_regions, oft_indent
USE oft_gs_util, ONLY: gs_save, gs_load, gs_analyze, gs_save_decon, gs_save_eqdsk, &
  gs_profile_load, gs_profile_save
USE oft_gs_fit, ONLY: fit_gs, fit_pm
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: i,ierr,io_unit,npts,iostat
LOGICAL :: file_exists
REAL(8), ALLOCATABLE :: pts(:,:)
TYPE(gs_eq) :: mygs
CLASS(oft_vector), POINTER :: xv
!---GS input options
INTEGER(4) :: order = 1
INTEGER(4) :: maxits = 30
INTEGER(4) :: ninner = 4
INTEGER(4) :: mode = 0
INTEGER(4) :: dcon_npsi = -1
INTEGER(4) :: dcon_ntheta = -1
INTEGER(4) :: eqdsk_nr = -1
INTEGER(4) :: eqdsk_nz = -1
LOGICAL :: pm = .FALSE.
REAL(8) :: urf = .3d0
REAL(8) :: nl_tol = 1.d-6
REAL(8) :: pnorm = 0.d0
REAL(8) :: alam = 0.d0
REAL(8) :: beta_mr = .3d0
REAL(8) :: f_offset = 0.d0
REAL(8) :: itor_target = -1.d0
REAL(8) :: estore_target = -1.d0
REAL(8) :: rmin = 0.d0
REAL(8) :: R0_target = -1.d0
REAL(8) :: V0_target = -1.d99
! REAL(8) :: rbounds(2) = (/-1.d99,1.d99/)
! REAL(8) :: zbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: lim_zmax = 1.d99
REAL(8) :: eqdsk_rbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: eqdsk_zbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: eqdsk_lcfs_pad = 0.01d0
REAL(8) :: init_r0(2)=[-1.d0,0.d0]
REAL(8) :: init_a=-1.d0
REAL(8) :: init_kappa=1.d0
REAL(8) :: init_delta=0.d0
LOGICAL :: free_boundary = .FALSE.
LOGICAL :: has_plasma = .TRUE.
LOGICAL :: save_mug = .FALSE.
LOGICAL :: fast_boundary = .TRUE.
LOGICAL :: limited_only = .FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: coil_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: limiter_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: eqdsk_filename = 'gTokaMaker'
CHARACTER(LEN=36) :: eqdsk_run_info = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: eqdsk_limiter_file = 'none'
!---Fit Input options
REAL(8) :: psinorm = 1.d0
LOGICAL :: adjust_pnorm = .FALSE.
LOGICAL :: adjust_alam = .FALSE.
LOGICAL :: adjust_R0 = .FALSE.
LOGICAL :: adjust_coils = .FALSE.
LOGICAL :: adjust_F0 = .FALSE.
LOGICAL :: adjust_V0 = .FALSE.
LOGICAL :: fixed_f = .FALSE.
LOGICAL :: fixed_p = .FALSE.
LOGICAL :: fixed_center = .FALSE.
NAMELIST/tokamaker_options/order,pm,mode,maxits,ninner,urf,nl_tol,itor_target,pnorm, &
alam,beta_mr,free_boundary,coil_file,limiter_file,f_offset,dcon_npsi,dcon_ntheta, &
has_plasma,rmin,R0_target,V0_target,save_mug,fast_boundary, &
limited_only,eqdsk_filename,eqdsk_nr,eqdsk_nz,eqdsk_rbounds,eqdsk_zbounds,eqdsk_run_info, &
eqdsk_limiter_file,eqdsk_lcfs_pad,init_r0,init_a,init_kappa,init_delta,lim_zmax,estore_target
NAMELIST/tokamaker_fit_options/psinorm,adjust_pnorm,adjust_alam,adjust_R0,adjust_V0, &
adjust_coils,adjust_F0,fixed_f,fixed_p,fixed_center
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---------------------------------------------------------------------------
! Load settings
!---------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,tokamaker_options,IOSTAT=ierr)
IF(ierr<0)CALL oft_abort('No "tokamaker_options" found in input file.', &
  'tokamaker_fit',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "tokamaker_options" in input file.', &
  'tokamaker_fit',__FILE__)
READ(io_unit,tokamaker_fit_options,IOSTAT=ierr)
IF(ierr<0)CALL oft_abort('No "tokamaker_fit_options" found in input file.', &
  'tokamaker_fit',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "tokamaker_fit_options" in input file.', &
  'tokamaker_fit',__FILE__)
CLOSE(io_unit)
oft_env%pm=pm
!---------------------------------------------------------------------------
! Check input files
!---------------------------------------------------------------------------
IF(TRIM(coil_file)/='none')THEN
  INQUIRE(EXIST=file_exists,FILE=TRIM(coil_file))
  IF(.NOT.file_exists)CALL oft_abort('Specified "coil_file" cannot be found', &
    'tokamaker_fit', __FILE__)
END IF
IF(TRIM(limiter_file)/='none')THEN
  INQUIRE(EXIST=file_exists,FILE=TRIM(limiter_file))
  IF(.NOT.file_exists)CALL oft_abort('Specified "limiter_file" cannot be found', &
    'tokamaker_fit', __FILE__)
END IF
IF(TRIM(eqdsk_limiter_file)/='none')THEN
  INQUIRE(EXIST=file_exists,FILE=TRIM(eqdsk_limiter_file))
  IF(.NOT.file_exists)CALL oft_abort('Specified "eqdsk_limiter_file" cannot be found', &
    'tokamaker_fit', __FILE__)
END IF
!---------------------------------------------------------------------------
! Setup Mesh
!---------------------------------------------------------------------------
CALL multigrid_construct_surf
CALL smesh%setup_io(order)
!---------------------------------------------------------------------------
! Setup Lagrange Elements
!---------------------------------------------------------------------------
CALL oft_lag_setup(order, -1)
!---------------------------------------------------------------------------
! Compute optimized smoother coefficients
!---------------------------------------------------------------------------
! CALL lag_mloptions
!---------------------------------------------------------------------------
! Setup experimental geometry
!---------------------------------------------------------------------------
! CALL exp_setup(mygs)
mygs%compute_chi=.FALSE.
mygs%free=free_boundary
IF(mygs%free)THEN
  mygs%itor_target=itor_target*mu0
  mygs%coil_file=coil_file
  CALL mygs%load_coils
ELSE
  IF(f_offset/=0.d0)mygs%itor_target=itor_target*mu0
END IF
CALL gs_setup_walls(mygs)
!---------------------------------------------------------------------------
! Setup profiles
!---------------------------------------------------------------------------
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'*** Loading flux and pressure profiles ***'
CALL gs_profile_load('f_prof.in',mygs%I)
mygs%I%f_offset=f_offset
IF(fixed_f)mygs%I%ncofs=0
CALL gs_profile_load('p_prof.in',mygs%P)
IF(fixed_p)mygs%P%ncofs=0
!---------------------------------------------------------------------------
! Initialize GS solution
!---------------------------------------------------------------------------
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'*** Initializing GS solution ***'
INQUIRE(EXIST=file_exists,FILE='tokamaker_psi_in.rst')
CALL mygs%init()!compute=(.NOT.file_exists),r0=init_r0,a=init_a,kappa=init_kappa,delta=init_delta)
IF(file_exists)THEN
  CALL oft_blagrange%vec_load(mygs%psi,'tokamaker_psi_in.rst','psi')
ELSE
  CALl mygs%init_psi(ierr,r0=init_r0,a=init_a,kappa=init_kappa,delta=init_delta)
  IF(ierr/=0)CALL oft_abort("Flux initialization failed","tokamaker_fit",__FILE__)
END IF
!---------------------------------------------------------------------------
! Compute GS fit
!---------------------------------------------------------------------------
mygs%mode=mode
mygs%urf=urf
mygs%maxits=maxits
mygs%nl_tol=nl_tol
mygs%pnorm=pnorm
mygs%psiscale=psinorm
mygs%limited_only=limited_only
mygs%R0_target=R0_target
mygs%V0_target=V0_target
mygs%limiter_file=limiter_file
CALL mygs%load_limiters
IF(mygs%free)THEN
  mygs%alam=alam
  mygs%boundary_limiter=.FALSE.
ELSE
  IF(f_offset/=0.d0)mygs%alam=alam
END IF
mygs%has_plasma=has_plasma
IF(.NOT.mygs%has_plasma)mygs%I%ncofs=0
! mygs%spatial_bounds(1,1)=MAX(MINVAL(smesh%r(1,:)),rbounds(1))
! mygs%spatial_bounds(2,1)=MIN(MAXVAL(smesh%r(1,:)),rbounds(2))
! mygs%spatial_bounds(1,2)=MAX(MINVAL(smesh%r(2,:)),zbounds(1))
! mygs%spatial_bounds(2,2)=MIN(MAXVAL(smesh%r(2,:)),zbounds(2))
mygs%lim_zmax=lim_zmax
mygs%rmin=rmin
mygs%estore_target=estore_target*mu0
oft_env%pm=.TRUE.
fit_pm=pm
mygs%plot_step=.FALSE.
mygs%plot_final=.TRUE.
!---Load input equilibrium
INQUIRE(EXIST=file_exists,FILE='tokamaker_fit_in.rst')
IF(file_exists)CALL gs_load(mygs,'tokamaker_fit_in.rst')
!---Solve
CALL fit_gs(mygs, fitPnorm=adjust_pnorm, fitAlam=adjust_alam, &
            fitR0=adjust_R0, fitV0=adjust_V0, fitCoils=adjust_coils, fitF0=adjust_F0, &
            fixedCentering=fixed_center)
!---------------------------------------------------------------------------
! Post-solution analysis
!---------------------------------------------------------------------------
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'*** Post-solution analysis ***'
!---Equilibrium information
CALL gs_analyze(mygs)
!---Save equilibrium and flux function
CALL gs_save(mygs,'tokamaker_fit.rst')
CALL oft_blagrange%vec_save(mygs%psi,'tokamaker_psi.rst','psi')
!---Save profile specifications
CALL gs_profile_save('f_prof.out',mygs%I)
CALL gs_profile_save('p_prof.out',mygs%P)
!---Save final flux profiles
CALL gs_save_prof(mygs,'fit.prof')
!---Save output grid
IF(save_mug)THEN
  CALL smesh%save_to_file('gs_trans_mesh.dat')
  CALL gs_save_fgrid(mygs,'gs_trans_fields.dat')
ELSE
  CALL gs_save_fgrid(mygs)
END IF
!---Sample fields at chosen locations
INQUIRE(EXIST=file_exists,FILE='tokamaker_fields.loc')
IF(file_exists)THEN
    OPEN(NEWUNIT=io_unit,FILE='tokamaker_fields.loc')
    READ(io_unit,*)npts
    ALLOCATE(pts(2,npts))
    DO i=1,npts
      READ(io_unit,*,IOSTAT=iostat)pts(:,i)
      IF(iostat<0)CALL oft_abort('EOF reached while reading "tokamaker_fields.loc"', &
                                 'tokamaker_fit',__FILE__)
    END DO
    CLOSE(io_unit)
    CALL gs_save_fields(mygs,pts,npts,'tokamaker_fields.dat')
    DEALLOCATE(pts)
ELSE
    WRITE(*,'(2A)')oft_indent,'No "tokamaker_fields.loc" file found, skipping field output'
END IF
!---Save DCON/EQDSK files
IF(has_plasma)THEN
  IF((dcon_npsi>0).AND.(dcon_ntheta>0))THEN
    CALL gs_save_decon(mygs,dcon_npsi,dcon_ntheta)
  END IF
  IF((eqdsk_nr>0).AND.(eqdsk_nz>0))THEN
    IF(ANY(eqdsk_rbounds<0.d0))CALL oft_abort('Invalid or unset EQDSK radial extents', &
                                              'tokamaker_fit',__FILE__)
    IF(ANY(ABS(eqdsk_zbounds)>1.d90))CALL oft_abort('Invalid or unset EQDSK vertical extents', &
                                              'tokamaker_fit',__FILE__)
    CALL gs_save_eqdsk(mygs,eqdsk_filename,eqdsk_nr,eqdsk_nz,eqdsk_rbounds,eqdsk_zbounds, &
      eqdsk_run_info,eqdsk_limiter_file,eqdsk_lcfs_pad)
  END IF
END IF
!---------------------------------------------------------------------------
! Terminate
!---------------------------------------------------------------------------
CALL oft_finalize
END PROGRAM tokamaker_fit
