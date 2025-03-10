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
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs
!---H1 FE (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_lop_eigs
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_wop_eigs
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_mop_eigs
IMPLICIT NONE
INTEGER(i4) :: ierr,io_unit
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_ml_fem_type), TARGET :: ML_oft_h1,ML_h1grad
TYPE(oft_ml_fem_type), TARGET :: ML_oft_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad
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
CALL multigrid_construct(mg_mesh)
!---------------------------------------------------------------------------
! Output mesh
!---------------------------------------------------------------------------
CALL plot_file%setup("mesh_check")
CALL mg_mesh%mesh%setup_io(plot_file,ABS(order))
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
IF(order>0)THEN
  !---Lagrange
  CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=minlev)
  !---H^1 (Grad(H^1) subspace)
  CALL oft_h1_setup(mg_mesh,order+1,ML_oft_h1,minlev=minlev)
  !---H(Curl) subspace
  CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl,minlev=minlev)
  !---H(Curl) + Grad(H^1) space
  CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,minlev)
!---------------------------------------------------------------------------
! Compute smoother coefficients
!---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'=============================================='
    WRITE(*,*)'Computing optimized smoother coefficients:'
    WRITE(*,*)'  NP = ',INT(order,2)
  END IF
  oft_env%pm=.FALSE.
  CALL lag_lop_eigs(ML_oft_lagrange,minlev)
  CALL h1_lop_eigs(ML_oft_h1,minlev)
  CALL hcurl_wop_eigs(ML_oft_hcurl,minlev)
  CALL hcurl_grad_mop_eigs(ML_hcurl_grad,minlev)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM oft_mesh_check
