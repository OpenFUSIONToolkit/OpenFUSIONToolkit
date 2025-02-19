!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> Regression test for xMHD module. A traveling alfven wave is
!! initialized in a triply periodic box and advanced for one half period.
!! The resulting wave is then compared to the initial wave to confirm basic
!! operation of the xMHD module.
!!
!! @authors Chris Hansen
!! @date November 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_alfven
USE oft_base
!--Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE fem_utils, ONLY: diff_interp
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_create, oft_lag_vcreate
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, &
  oft_lag_vgetmop, oft_lag_vrinterp
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_level, oft_hcurl_nlevels
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions, hcurl_zerob
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1_setup_interp, h1_getmop, oft_h1_project, &
  h1grad_zerop, oft_h1_rinterp
!---Physics
USE oft_vector_inits, ONLY: uniform_field
USE diagnostic, ONLY: vec_energy
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, xmhd_lin_run, &
  xmhd_sub_fields
USE test_phys_helpers, ONLY: alfven_eig
IMPLICIT NONE
!---H1 metric solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Lagrange metric solver
CLASS(oft_solver), POINTER :: lag_minv => NULL()
CLASS(oft_matrix), POINTER :: lag_mop => NULL()
!---Local variables
CLASS(oft_vector), POINTER :: b,vel,db,dvel,u,v,lag_u,lag_v,bi,vi,be
CLASS(oft_vector), POINTER :: den,temp,dden,dtemp
TYPE(xmhd_sub_fields) :: equil_fields,pert_fields
TYPE(uniform_field) :: z_field
TYPE(alfven_eig), TARGET :: alf_field
TYPE(oft_h1_rinterp), TARGET :: bfield
TYPE(oft_lag_vrinterp), TARGET :: vfield
TYPE(diff_interp) :: err_field
INTEGER(i4) :: io_unit
REAL(r8) :: B0,verr,vierr,berr,bierr
REAL(r8), POINTER :: tmp(:),bvals(:),uvals(:)
INTEGER(i4) :: ierr
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: delta = 1.d-4
REAL(r8) :: v_alf = 1.d4
LOGICAL :: linear = .FALSE.
NAMELIST/test_alfven_options/order,minlev,delta,linear
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_alfven_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
rgrnd=(/2.d0,0.d0,0.d0/)
CALL multigrid_construct
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order,minlev)
CALL lag_setup_interp
!---H1(Curl) subspace
CALL oft_hcurl_setup(order,minlev)
CALL hcurl_setup_interp
!---H1(Grad) subspace
CALL oft_h0_setup(order+1,minlev)
CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order,minlev)
CALL h1_setup_interp
!---------------------------------------------------------------------------
! Create H1 metric solver
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL h1_getmop(mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL oft_h1_create(u)
CALL oft_h1_create(v)
CALL oft_h1_create(b)
CALL oft_h1_create(db)
CALL oft_h1_create(be)
CALL oft_h1_create(bi)
!---------------------------------------------------------------------------
! Set uniform B0 = zhat
!---------------------------------------------------------------------------
z_field%val=(/0.d0,0.d0,1.d0/)
CALL oft_h1_project(z_field,v)
CALL h1grad_zerop(v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL b%add(0.d0,1.d0,u)
CALL be%add(0.d0,1.d0,u)
!---------------------------------------------------------------------------
! Set dB from alfven wave init
!---------------------------------------------------------------------------
CALL oft_h1_project(alf_field,v)
CALL h1grad_zerop(v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL db%add(0.d0,delta,u)
B0=v_alf*SQRT(mu0*2.d19*proton_mass)
CALL bi%add(0.d0,1.d0,u)
!---------------------------------------------------------------------------
! Create Lagrange metric solver
!---------------------------------------------------------------------------
NULLIFY(lag_mop)
CALL oft_lag_vgetmop(lag_mop,"none")
CALL create_cg_solver(lag_minv)
lag_minv%A=>lag_mop
lag_minv%its=-3
lag_minv%atol=1.d-10
CALL create_diag_pre(lag_minv%pre)
!---
CALL oft_lag_vcreate(lag_u)
CALL oft_lag_vcreate(lag_v)
CALL oft_lag_vcreate(vel)
CALL oft_lag_vcreate(dvel)
CALL oft_lag_vcreate(vi)
!---------------------------------------------------------------------------
! Set dB from alfven wave init
!---------------------------------------------------------------------------
CALL oft_lag_vproject(alf_field,lag_v)
CALL lag_u%set(0.d0)
CALL lag_minv%apply(lag_u,lag_v)
CALL dvel%add(0.d0,delta,lag_u)
CALL vi%add(0.d0,1.d0,lag_u)
!---------------------------------------------------------------------------
! Run simulation and test result
!---------------------------------------------------------------------------
xmhd_minlev=minlev
xmhd_taxis=2
oft_env%pm=.FALSE.
IF(linear)THEN
  CALL b%scale(B0)
  CALL db%scale(B0)
  CALL dvel%scale(v_alf)
  !
  CALL oft_lag_create(den)
  CALL oft_lag_create(temp)
  CALL oft_lag_create(dden)
  CALL oft_lag_create(dtemp)
  CALL den%set(1.d19)
  CALL temp%set(1.d1)
  equil_fields%B=>b
  equil_fields%V=>vel
  equil_fields%Ne=>den
  equil_fields%Ti=>temp
  pert_fields%B=>db
  pert_fields%V=>dvel
  pert_fields%Ne=>dden
  pert_fields%Ti=>dtemp
  CALL xmhd_lin_run(equil_fields,pert_fields)
  CALL den%delete
  CALL dden%delete
  CALL temp%delete
  CALL dtemp%delete
  DEALLOCATE(den,dden,temp,dtemp)
  CALL b%add(1.d0,1.d0,db)
  CALL vel%add(1.d0,1.d0,dvel)
ELSE
  CALL b%add(1.d0,1.d0,db)
  CALL b%scale(B0)
  CALL vel%add(0.d0,1.d0,dvel)
  CALL vel%scale(v_alf)
  !
  CALL oft_lag_create(den)
  CALL oft_lag_create(temp)
  CALL den%set(1.d19)
  CALL temp%set(1.d1)
  equil_fields%B=>b
  equil_fields%V=>vel
  equil_fields%Ne=>den
  equil_fields%Ti=>temp
  CALL xmhd_run(equil_fields)
  CALL den%delete
  CALL temp%delete
  DEALLOCATE(den,temp)
END IF
!---Check magnetic field waveform
bierr=vec_energy(alf_field,order*2)
err_field%dim=3
err_field%a=>alf_field
err_field%b=>bfield
CALL b%scale(1.d0/B0)
CALL b%add(-1.d0,1.d0,be)
CALL b%scale(1.d0/delta)
bfield%u=>b
CALL bfield%setup
berr=vec_energy(err_field,order*2)
!---Check velocity field waveform
vierr=vec_energy(alf_field,order*2)
err_field%dim=3
err_field%a=>alf_field
err_field%b=>vfield
CALL vel%scale(-1.d0/(delta*v_alf))
vfield%u=>vel
CALL vfield%setup
verr=vec_energy(err_field,order*2)
!---Output wave comparisons
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='alfven.results')
  WRITE(io_unit,*)SQRT(berr/bierr)
  WRITE(io_unit,*)SQRT(verr/vierr)
  CLOSE(io_unit)
END IF
!---Test plotting routine
CALL xmhd_plot
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_alfven
