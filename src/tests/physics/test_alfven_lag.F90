!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> Regression test for xMHD_lag module. A traveling alfven wave is
!! initialized in a triply periodic box and advanced for one half period.
!! The resulting wave is then compared to the initial wave to confirm basic
!! operation of the xMHD_lag module.
!!
!! @authors Chris Hansen
!! @date November 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_alfven
USE oft_base
!--Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE fem_utils, ONLY: diff_interp
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, &
  oft_lag_vgetmop, oft_lag_vrinterp
!---Physics
USE oft_vector_inits, ONLY: uniform_field
USE diagnostic, ONLY: vec_energy
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_lag, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, xmhd_lin_run, &
  xmhd_sub_fields, xmhd_ML_lagrange, xmhd_ML_vlagrange
USE test_phys_helpers, ONLY: alfven_eig
IMPLICIT NONE
!---Lagrange mass matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Local variables
CLASS(oft_vector), POINTER :: b,vel,db,dvel,u,v,bi,vi,be
CLASS(oft_vector), POINTER :: den,temp,dden,dtemp
TYPE(xmhd_sub_fields) :: equil_fields,pert_fields
TYPE(uniform_field) :: z_field
TYPE(alfven_eig), TARGET :: alf_field
TYPE(oft_lag_vrinterp), TARGET :: Vfield,Bfield
TYPE(diff_interp) :: err_field
TYPE(multigrid_mesh) :: mg_mesh
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
CALL multigrid_construct(mg_mesh)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
!---Lagrange
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!---------------------------------------------------------------------------
! Create Lagrange metric solver
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(b)
CALL xmhd_ML_vlagrange%vec_create(db)
CALL xmhd_ML_vlagrange%vec_create(be)
CALL xmhd_ML_vlagrange%vec_create(bi)
CALL xmhd_ML_vlagrange%vec_create(vel)
CALL xmhd_ML_vlagrange%vec_create(dvel)
CALL xmhd_ML_vlagrange%vec_create(vi)
!---------------------------------------------------------------------------
! Set uniform B0 = zhat
!---------------------------------------------------------------------------
z_field%val=(/0.d0,0.d0,1.d0/)
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,z_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL b%add(0.d0,1.d0,u)
CALL be%add(0.d0,1.d0,u)
!---------------------------------------------------------------------------
! Set dB from alfven wave init
!---------------------------------------------------------------------------
alf_field%mesh=>mg_mesh%mesh
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,alf_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL db%add(0.d0,delta,u)
B0=v_alf*SQRT(mu0*2.d19*proton_mass)
CALL bi%add(0.d0,1.d0,u)
!---------------------------------------------------------------------------
! Set dV from alfven wave init
!---------------------------------------------------------------------------
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,alf_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL dvel%add(0.d0,delta,u)
CALL vi%add(0.d0,1.d0,u)
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
  CALL xmhd_ML_lagrange%vec_create(den)
  CALL xmhd_ML_lagrange%vec_create(temp)
  CALL xmhd_ML_lagrange%vec_create(dden)
  CALL xmhd_ML_lagrange%vec_create(dtemp)
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
  CALL xmhd_ML_lagrange%vec_create(den)
  CALL xmhd_ML_lagrange%vec_create(temp)
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
bierr=vec_energy(mg_mesh%mesh,alf_field,order*2)
err_field%dim=3
err_field%a=>alf_field
err_field%b=>bfield
CALL b%scale(1.d0/B0)
CALL b%add(-1.d0,1.d0,be)
CALL b%scale(1.d0/delta)
bfield%u=>b
CALL bfield%setup(xmhd_ML_lagrange%current_level)
berr=vec_energy(mg_mesh%mesh,err_field,order*2)
!---Check velocity field waveform
vierr=vec_energy(mg_mesh%mesh,alf_field,order*2)
err_field%dim=3
err_field%a=>alf_field
err_field%b=>vfield
CALL vel%scale(-1.d0/(delta*v_alf))
vfield%u=>vel
CALL vfield%setup(xmhd_ML_lagrange%current_level)
verr=vec_energy(mg_mesh%mesh,err_field,order*2)
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
