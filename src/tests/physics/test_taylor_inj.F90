!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_taylor_inj.F90
!
!> Regression tests for computing Taylor states in a spherical geometry.
!! States are computed at different polynomial orders and the resulting
!! eigenvalues and toroidal currents are compared to reference cases.
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_taylor_inj
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, oft_lag_vgetmop, &
  oft_lag_vproject
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_mloptions, h1_setup_interp
USE taylor, ONLY: taylor_vacuum, taylor_injectors, taylor_injector_single, &
  oft_taylor_rinterp, oft_taylor_hmodes, oft_taylor_ifield, taylor_tag_size
implicit none
INTEGER(i4) :: ierr,io_unit
REAL(r8) :: comps(3),diff_err
CLASS(oft_vector), POINTER :: gffa
INTEGER(i4), PARAMETER :: nh = 1
REAL(r8) :: fluxes(nh),hcpc(3,nh),hcpv(3,nh),energy(nh)
CHARACTER(LEN=taylor_tag_size) :: htags(nh)
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_comp_type) :: ML_vlagrange
TYPE(oft_taylor_hmodes) :: hmodes
TYPE(oft_taylor_ifield) :: ff_obj
INTEGER(i4) :: order=1
LOGICAL :: mg_test=.FALSE.
NAMELIST/test_taylor_options/order,mg_test
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_taylor_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("Test")
CALL mg_mesh%mesh%setup_io(plot_file,order)
IF(mg_test)THEN
  ff_obj%minlev=2
ELSE
  ff_obj%minlev=mg_mesh%mgmax+order-1
END IF
hmodes%minlev=ff_obj%minlev
!---
ALLOCATE(ff_obj%ML_hcurl,ff_obj%ML_h1,ff_obj%ML_hcurl_grad,ff_obj%ML_h1grad,ff_obj%ML_lagrange)
hmodes%ML_hcurl=>ff_obj%ML_hcurl
hmodes%ML_lagrange=>ff_obj%ML_lagrange
CALL oft_lag_setup(mg_mesh,order,ff_obj%ML_lagrange,ML_vlag_obj=ML_vlagrange,minlev=ff_obj%minlev)
CALL oft_hcurl_setup(mg_mesh,order,ff_obj%ML_hcurl,minlev=ff_obj%minlev)
CALL oft_h1_setup(mg_mesh,order+1,ff_obj%ML_h1,minlev=ff_obj%minlev)
CALL oft_hcurl_grad_setup(ff_obj%ML_hcurl,ff_obj%ML_h1,ff_obj%ML_hcurl_grad,ff_obj%ML_h1grad,ff_obj%minlev)
IF(mg_test)THEN
  CALL lag_setup_interp(ff_obj%ML_lagrange)
  CALL lag_mloptions
  CALL hcurl_setup_interp(ff_obj%ML_hcurl)
  CALL hcurl_mloptions(ff_obj%ML_hcurl)
  CALL h1_setup_interp(ff_obj%ML_h1)
  CALL h1_mloptions
END IF
!---Define jumps
hcpc(:,1)=(/1.d0, 0.d0, 0.d0/)
hcpv(:,1)=(/0.d0, 1.d0, 0.d0/)*1.5d0
htags(1)='Tinj'
CALL ff_obj%setup(nh,hcpc,hcpv,htags)
!---Run tests
ff_obj%jtol=1.d-4
oft_env%pm=.FALSE.
CALL taylor_vacuum(ff_obj,energy=energy,hmodes=hmodes,rst_filename='oft_taylor.rst')
CALL taylor_injectors(ff_obj,hmodes,5.d0,rst_filename='oft_taylor.rst')
comps(1) = ff_obj%hvac(1)%f%dot(ff_obj%hvac(1)%f)
comps(2) = ff_obj%hcur(1)%f%dot(ff_obj%hcur(1)%f)
comps(3) = ff_obj%gffa(1)%f%dot(ff_obj%gffa(1)%f)
CALL ff_obj%gffa(1)%f%new(gffa)
CALL taylor_injector_single(ff_obj,hmodes,5.d0,(/1.d0/),gffa)
CALL gffa%add(1.d0,-1.d0,ff_obj%gffa(1)%f)
diff_err = gffa%dot(gffa)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='taylor.results')
  WRITE(io_unit,*)energy
  WRITE(io_unit,*)comps
  WRITE(io_unit,*)diff_err
  CLOSE(io_unit)
END IF
CALL gffa%delete()
DEALLOCATE(gffa)
!---Test plotting Taylor fields
BLOCK
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), POINTER, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
TYPE(oft_taylor_rinterp), TARGET :: Bfield
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(ML_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL ML_vlagrange%current_level%vec_create(u)
CALL ML_vlagrange%current_level%vec_create(v)
!---Plot solution
Bfield%uvac=>ff_obj%hvac(1)%f
Bfield%ua=>ff_obj%gffa(1)%f
CALL Bfield%setup(ff_obj%ML_hcurl%current_level,ff_obj%ML_h1%current_level)
!---Project field
CALL oft_lag_vproject(ff_obj%ML_lagrange%current_level,Bfield,v)
CALL u%set(0.d0)
CALL lminv%apply(u,v)
!---Retrieve local values and save
ALLOCATE(bvout(3,u%n/3))
vals=>bvout(1,:)
CALL u%get_local(vals,1)
vals=>bvout(2,:)
CALL u%get_local(vals,2)
vals=>bvout(3,:)
CALL u%get_local(vals,3)
CALL mg_mesh%mesh%save_vertex_vector(bvout,plot_file,'B')
END BLOCK
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_taylor_inj
