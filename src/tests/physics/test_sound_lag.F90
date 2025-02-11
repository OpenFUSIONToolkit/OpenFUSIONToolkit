!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> Regression test for xMHD_lag module. A traveling sound wave is initialized in
!! a triply periodic box and advanced for one period. The resulting wave is
!! then compared to the initial wave to confirm basic operation of the xMHD
!! module.
!!
!! @authors Chris Hansen
!! @date November 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_sound
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
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project, lag_zerob, oft_lag_rinterp, oft_lag_vrinterp
!---Physics
USE diagnostic, ONLY: scal_energy, vec_energy
USE mhd_utils, ONLY: elec_charge, proton_mass
USE xmhd_lag, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, temp_floor, &
  xmhd_lin_run, xmhd_adv_b, xmhd_sub_fields
USE test_phys_helpers, ONLY: sound_eig
IMPLICIT NONE
!---Lagrange Metric solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Local variables
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_vector), POINTER :: b,vel,n,temp,tempe
CLASS(oft_vector), POINTER :: db,dvel,dn,dtemp,dtempe
CLASS(oft_vector), POINTER :: ti,ni,vi
TYPE(xmhd_sub_fields) :: equil_fields,pert_fields
TYPE(oft_lag_rinterp), TARGET :: sfield
TYPE(oft_lag_vrinterp), TARGET :: vfield
TYPE(sound_eig), TARGET :: sound_field
TYPE(diff_interp) :: err_field
INTEGER(i4) :: io_unit
REAL(r8) :: T0,v_delta,nerr,nierr,terr,tierr,verr,vierr
REAL(r8), POINTER :: tmp(:),bvals(:),uvals(:),vals(:)
INTEGER(i4) :: ierr
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: v_sound = 2.d4
REAL(r8) :: N0 = 1.d19
REAL(r8) :: delta = 1.d-4
LOGICAL :: linear = .FALSE.
LOGICAL :: two_temp = .FALSE.
NAMELIST/test_sound_options/order,minlev,delta,linear,two_temp
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_sound_options,IOSTAT=ierr)
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
!---------------------------------------------------------------------------
! Create Lagrange metric solver
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_getmop(mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL oft_lag_create(u)
CALL oft_lag_create(v)
CALL oft_lag_create(n)
CALL oft_lag_create(dn)
CALL oft_lag_create(ni)
CALL oft_lag_create(temp)
CALL oft_lag_create(dtemp)
CALL oft_lag_create(ti)
!---------------------------------------------------------------------------
! Set dn from sound wave init
!---------------------------------------------------------------------------
sound_field%delta=delta
sound_field%field='n'
CALL oft_lag_project(sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL n%set(1.d0)
CALL dn%add(0.d0,1.d0,u,-1.d0,n)
CALL ni%add(0.d0,1.d0,dn)
!---------------------------------------------------------------------------
! Set dt from sound wave init
!---------------------------------------------------------------------------
sound_field%field='t'
CALL oft_lag_project(sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL temp%set(1.d0)
CALL dtemp%add(0.d0,1.d0,u,-1.d0,temp)
CALL ti%add(0.d0,1.d0,dtemp)
T0=(v_sound**2)*3.d0*proton_mass/(5.d0*elec_charge)
!---
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
CALL minv%delete
!---------------------------------------------------------------------------
! Create Lagrange vector metric solver
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_vgetmop(mop,"none")
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL oft_lag_vcreate(u)
CALL oft_lag_vcreate(v)
CALL oft_lag_vcreate(vel)
CALL oft_lag_vcreate(dvel)
CALL oft_lag_vcreate(vi)
!---------------------------------------------------------------------------
! Set dV from sound wave init
!---------------------------------------------------------------------------
sound_field%field='v'
CALL oft_lag_vproject(sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
v_delta=T0*elec_charge/(proton_mass*v_sound)
CALL dvel%add(0.d0,1.d0,u)
CALL vi%add(0.d0,1.d0,u)
!---
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
!---------------------------------------------------------------------------
! Run simulation and test result
!---------------------------------------------------------------------------
CALL oft_lag_vcreate(b)
CALL oft_lag_vcreate(db)
xmhd_minlev=minlev
xmhd_taxis=2
xmhd_adv_b=.FALSE.
oft_env%pm=.FALSE.
IF(linear)THEN
  CALL dvel%scale(v_delta)
  CALL n%scale(N0)
  CALL dn%scale(N0)
  CALL temp%scale(T0)
  CALL dtemp%scale(T0)
  IF(two_temp)THEN
    CALL oft_lag_create(tempe)
    CALL tempe%add(0.d0,1.d0,temp)
    CALL oft_lag_create(dtempe)
    CALL dtempe%add(0.d0,1.d0,dtemp)
    equil_fields%Te=>tempe
    pert_fields%Te=>dtempe
  END IF
  equil_fields%B=>b
  equil_fields%V=>vel
  equil_fields%Ne=>n
  equil_fields%Ti=>temp
  pert_fields%B=>db
  pert_fields%V=>dvel
  pert_fields%Ne=>dn
  pert_fields%Ti=>dtemp
  CALL xmhd_lin_run(equil_fields,pert_fields)
  CALL b%add(1.d0,1.d0,db)
  CALL vel%add(1.d0,1.d0,dvel)
  CALL n%add(1.d0,1.d0,dn)
  CALL temp%add(1.d0,1.d0,dtemp)
  IF(two_temp)CALL tempe%add(1.d0,1.d0,dtempe)
ELSE
  temp_floor=T0*1.e-2
  CALL vel%add(0.d0,v_delta,dvel)
  CALL n%add(1.d0,1.d0,dn)
  CALL n%scale(N0)
  CALL temp%add(1.d0,1.d0,dtemp)
  CALL temp%scale(T0)
  IF(two_temp)THEN
    CALL oft_lag_create(tempe)
    CALL tempe%add(0.d0,1.d0,temp)
    equil_fields%Te=>tempe
  END IF
  equil_fields%B=>b
  equil_fields%V=>vel
  equil_fields%Ne=>n
  equil_fields%Ti=>temp
  CALL xmhd_run(equil_fields)
END IF
CALL oft_lag_create(u)
CALL u%set(1.d0)
!---Compare density waveform
sound_field%field='n'
sound_field%diff=.TRUE.
nierr=scal_energy(sound_field,order*2)
err_field%dim=1
err_field%a=>sound_field
err_field%b=>sfield
CALL n%scale(1.d0/N0)
CALL n%add(1.d0,-1.d0,u)
sfield%u=>n
CALL sfield%setup
nerr=scal_energy(err_field,order*2)
!---Compare temperature waveform
sound_field%field='t'
sound_field%diff=.TRUE.
tierr=scal_energy(sound_field,order*2)
err_field%dim=1
err_field%a=>sound_field
err_field%b=>sfield
CALL temp%scale(1.d0/T0)
CALL temp%add(1.d0,-1.d0,u)
IF(two_temp)THEN
  CALL tempe%scale(1.d0/T0)
  CALL tempe%add(1.d0,-1.d0,u)
  CALL temp%add(5.d-1,5.d-1,tempe)
END IF
sfield%u=>temp
CALL sfield%setup
terr=scal_energy(err_field,order*2)
!---Compare velocity waveform
sound_field%field='v'
sound_field%diff=.FALSE.
vierr=vec_energy(sound_field,order*2)
err_field%dim=3
err_field%a=>sound_field
err_field%b=>vfield
CALL vel%scale(1.d0/v_delta)
vfield%u=>vel
CALL vfield%setup
verr=vec_energy(err_field,order*2)
!---Output wave comparisons
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='sound.results')
  WRITE(io_unit,*)SQRT(nerr/nierr)
  WRITE(io_unit,*)SQRT(terr/tierr)
  WRITE(io_unit,*)SQRT(verr/vierr)
  CLOSE(io_unit)
END IF
!---Test plotting routine
CALL xmhd_plot
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_sound
