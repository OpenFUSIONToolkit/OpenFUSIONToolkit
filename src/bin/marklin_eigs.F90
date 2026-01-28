!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file marklin_eigs.F90
!
!> Driver to compute Taylor states (and higher force-free modes)
!!
!! **Option group:** `marklin_eigs_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `order=1`              |  FE order     | int |
!! |  `nmodes=1`             |  Number of modes to compute | int |
!! |  `minlev=-1`            |  Lowest level for MG (negative to disable) | int |
!!
!! @authors Chris Hansen
!! @date June 2018
!! @ingroup doxy_oft_bin
!--------------------------------------------------------------------------------
PROGRAM marklin_eigs
!---Base
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!--Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
  oft_lag_vgetmop, oft_lag_vproject
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
  hcurl_mloptions
!---Taylor state
USE taylor, ONLY: taylor_hmodes, oft_taylor_hmodes
IMPLICIT NONE
#include "local.h"
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
CHARACTER(LEN=3) :: pltnum
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_comp_type) :: ML_vlagrange
TYPE(oft_taylor_hmodes) :: hmodes
INTEGER(i4) :: order = 2
INTEGER(i4) :: nmodes = 1
INTEGER(i4) :: minlev = 1
NAMELIST/marklin_eigs_options/order,nmodes,minlev
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,marklin_eigs_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check inputs
IF(ierr<0)CALL oft_abort('No "marklin_eigs_options" found in input file.', &
  'marklin_eigs',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "marklin_eigs_options" in input file.', &
  'marklin_eigs',__FILE__)
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("marklin_eigs")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!
IF(minlev<0)THEN
  hmodes%minlev=mg_mesh%mgdim+(order-1)
ELSE
  hmodes%minlev=minlev
  IF(oft_env%nprocs>1)hmodes%minlev=MAX(oft_env%nbase+1,minlev)
END IF
ALLOCATE(hmodes%ML_hcurl,hmodes%ML_lagrange)
!---Lagrange
CALL oft_lag_setup(mg_mesh,order,hmodes%ML_lagrange,ML_vlag_obj=ML_vlagrange,minlev=hmodes%minlev)
CALL lag_setup_interp(hmodes%ML_lagrange)
CALL lag_mloptions
!---H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,hmodes%ML_hcurl,minlev=hmodes%minlev)
CALL hcurl_setup_interp(hmodes%ML_hcurl)
CALL hcurl_mloptions(hmodes%ML_hcurl)
oft_env%pm=.TRUE.
CALL taylor_hmodes(hmodes,nmodes)
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(ML_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL ML_vlagrange%vec_create(u)
CALL ML_vlagrange%vec_create(v)
ALLOCATE(bvout(3,u%n/3))
!---Save modes
DO i=1,nmodes
  WRITE(pltnum,'(I3.3)')i
  CALL hmodes%ML_hcurl%current_level%vec_save(hmodes%hffa(i,hmodes%ML_hcurl%level)%f, &
                          'marklin_eigs.rst','A_'//pltnum, append=(i/=1))
  !---Setup field interpolation
  Bfield%u=>hmodes%hffa(i,hmodes%ML_hcurl%level)%f
  CALL Bfield%setup(hmodes%ML_hcurl%current_level)
  !---Project field
  CALL oft_lag_vproject(hmodes%ML_lagrange%current_level,Bfield,v)
  CALL u%set(0.d0)
  CALL lminv%apply(u,v)
  !---Retrieve local values and save
  vals=>bvout(1,:)
  CALL u%get_local(vals,1)
  vals=>bvout(2,:)
  CALL u%get_local(vals,2)
  vals=>bvout(3,:)
  CALL u%get_local(vals,3)
  call mg_mesh%mesh%save_vertex_vector(bvout,plot_file,'B_'//pltnum)
END DO
!---Finalize enviroment
CALL oft_finalize
END PROGRAM marklin_eigs
