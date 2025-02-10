!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_taylor_inj.F90
!
!> Regression tests for computing Taylor states in a spherical geometry.
!! States are computed at different polynomial orders and the resulting
!! eigenvalues and toroidal currents are compared to reference cases.
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_taylor_inj
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_mesh_type, ONLY: mesh
USE multigrid, ONLY: mg_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_lag_setup, oft_vlagrange
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, oft_lag_vgetmop, &
  oft_lag_vproject
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_nlevels
USE oft_hcurl_fields, ONLY: oft_hcurl_create
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_mloptions, h0_setup_interp
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_level
USE taylor, ONLY: taylor_vacuum, taylor_injectors, taylor_injector_single, &
  taylor_minlev, taylor_jtol, taylor_tag_size, taylor_hvac, taylor_hcur, &
  taylor_gffa, oft_taylor_rinterp
implicit none
INTEGER(i4) :: ierr,io_unit
REAL(r8) :: comps(3),diff_err
CLASS(oft_vector), POINTER :: gffa
INTEGER(i4), PARAMETER :: nh = 1
REAL(r8) :: fluxes(nh),hcpc(3,nh),hcpv(3,nh),energy(nh)
CHARACTER(LEN=taylor_tag_size) :: htags(nh)
TYPE(xdmf_plot_file) :: plot_file
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
CALL multigrid_construct
CALL plot_file%setup("Test")
CALL mesh%setup_io(plot_file,order)
IF(mg_test)THEN
  taylor_minlev=2
ELSE
  taylor_minlev=mg_mesh%mgmax+order-1
END IF
!---
CALL oft_lag_setup(order,taylor_minlev)
CALL oft_hcurl_setup(order,taylor_minlev)
CALL oft_h0_setup(order+1,taylor_minlev)
CALL oft_h1_setup(order,taylor_minlev)
IF(mg_test)THEN
  CALL lag_setup_interp
  CALL lag_mloptions
  CALL hcurl_setup_interp
  CALL hcurl_mloptions
  CALL h0_setup_interp
  CALL h0_mloptions
END IF
!---Define jumps
hcpc(:,1)=(/1.d0, 0.d0, 0.d0/)
hcpv(:,1)=(/0.d0, 1.d0, 0.d0/)*1.5d0
htags(1)='Tinj'
!---Run tests
taylor_jtol=1.d-4
oft_env%pm=.FALSE.
CALL taylor_vacuum(nh,hcpc,hcpv,htags,energy)
CALL taylor_injectors(5.d0)
comps(1) = taylor_hvac(1,oft_h1_level)%f%dot(taylor_hvac(1,oft_h1_level)%f)
comps(2) = taylor_hcur(1,oft_h1_level)%f%dot(taylor_hcur(1,oft_h1_level)%f)
comps(3) = taylor_gffa(1,oft_h1_level)%f%dot(taylor_gffa(1,oft_h1_level)%f)
CALL taylor_gffa(1,oft_h1_level)%f%new(gffa)
CALL taylor_injector_single(5.d0,(/1.d0/),gffa)
CALL gffa%add(1.d0,-1.d0,taylor_gffa(1,oft_h1_level)%f)
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
CALL oft_lag_vgetmop(lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL oft_vlagrange%vec_create(u)
CALL oft_vlagrange%vec_create(v)
!---Plot solution
Bfield%uvac=>taylor_hvac(1,oft_h1_level)%f
Bfield%ua=>taylor_gffa(1,oft_h1_level)%f
CALL Bfield%setup
!---Project field
CALL oft_lag_vproject(Bfield,v)
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
CALL mesh%save_vertex_vector(bvout,plot_file,'B')
END BLOCK
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_taylor_inj
