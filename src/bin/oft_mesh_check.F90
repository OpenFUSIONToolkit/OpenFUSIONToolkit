!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_check.F90
!
!> @defgroup doxy_oft_bin Driver Programs
!! Utility driver programs
!
!> Driver to check grid generation and compute optimized smoother parameters
!!
!! **Option group:** `oft_mesh_check_options`
!! |  Option  |  Description  | Type [dim] |
!! |------------|------------|-----|
!! |  `order=1`   |  FE order  | int |
!! |  `minlev=1`  |  Lowest FE level to use  | int |
!!
!! @authors Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_bin
!-----------------------------------------------------------------------------
PROGRAM oft_mesh_check
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!--Grid
USE oft_mesh_type, ONLY: mesh
USE multigrid_build, ONLY: multigrid_construct
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_wop_eigs
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_lop_eigs
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_mop_eigs
IMPLICIT NONE
INTEGER(i4) :: ierr,io_unit
TYPE(xdmf_plot_file) :: plot_file
INTEGER(i4) :: order=1
INTEGER(i4) :: minlev=1
NAMELIST/oft_mesh_check_options/order,minlev
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,oft_mesh_check_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct
!---------------------------------------------------------------------------
! Output mesh
!---------------------------------------------------------------------------
CALL plot_file%setup("mesh_check")
CALL mesh%setup_io(plot_file,ABS(order))
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
IF(order>0)THEN
  !---Lagrange
  CALL oft_lag_setup(order,minlev)
  !---H1(Curl) subspace
  CALL oft_hcurl_setup(order,minlev)
  !---H1(Grad) subspace
  CALL oft_h0_setup(order+1,minlev)
  !---H1 space
  CALL oft_h1_setup(order,minlev)
!---------------------------------------------------------------------------
! Compute smoother coefficients
!---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'=============================================='
    WRITE(*,*)'Computing optimized smoother coefficients:'
    WRITE(*,*)'  NP = ',INT(order,2)
  END IF
  oft_env%pm=.FALSE.
  CALL lag_lop_eigs(minlev)
  CALL h0_lop_eigs(minlev)
  CALL hcurl_wop_eigs(minlev)
  CALL h1_mop_eigs(minlev)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM oft_mesh_check
