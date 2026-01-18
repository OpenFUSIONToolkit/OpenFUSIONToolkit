!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_taylor.F90
!
!> Regression tests for computing Taylor states in a spherical geometry.
!! States are computed at different polynomial orders and the resulting
!! eigenvalues and toroidal currents are compared to reference cases.
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_taylor
USE oft_base
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions
USE taylor, ONLY: taylor_hmodes, oft_taylor_hmodes
implicit none
INTEGER(i4) :: order=1,nm=1,ierr,io_unit
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_taylor_hmodes) :: taylor_states
LOGICAL :: mg_test=.FALSE.
NAMELIST/test_taylor_options/order,nm,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_taylor_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_sphere_id)CALL oft_abort('Wrong mesh type, test for SPHERE only.','main',__FILE__)
IF(mg_test)THEN
  taylor_states%minlev=2
  IF(oft_env%nprocs>1)taylor_states%minlev=mg_mesh%nbase+1
  IF(mg_mesh%mesh%type==3)taylor_states%minlev=mg_mesh%mgmax
ELSE
  taylor_states%minlev=mg_mesh%mgmax+order-1
END IF
ALLOCATE(taylor_states%ML_hcurl,taylor_states%ML_lagrange)
!---
CALL oft_hcurl_setup(mg_mesh,order,taylor_states%ML_hcurl,minlev=taylor_states%minlev)
IF(mg_test)THEN
  CALL hcurl_setup_interp(taylor_states%ML_hcurl)
  CALL hcurl_mloptions(taylor_states%ML_hcurl)
END IF
!---
CALL oft_lag_setup(mg_mesh,order,taylor_states%ML_lagrange,minlev=taylor_states%minlev)
IF(mg_test)THEN
  CALL lag_setup_interp(taylor_states%ML_lagrange)
  CALL lag_mloptions
END IF
!---Run tests
oft_env%pm=.FALSE.
CALL taylor_hmodes(taylor_states,nm,rst_filename='oft_taylor.rst')
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='taylor.results')
  WRITE(io_unit,*)taylor_states%hlam(:,taylor_states%ML_hcurl%nlevels)
  WRITE(io_unit,*)taylor_states%htor(:,taylor_states%ML_hcurl%nlevels)
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_taylor
