!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> Regression test for xMHD module. A traveling sound wave is initialized in
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
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE fem_utils, ONLY: diff_interp
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project, oft_lag_rinterp, oft_lag_vrinterp
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, oft_hcurl_grad_project
!---Physics
USE diagnostic, ONLY: scal_energy, vec_energy
USE mhd_utils, ONLY: elec_charge, proton_mass
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, temp_floor, &
  xmhd_lin_run, xmhd_adv_b, xmhd_sub_fields, ML_oft_hcurl, ML_oft_h1, &
  ML_hcurl_grad, ML_h1grad, ML_oft_lagrange, ML_oft_vlagrange
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
TYPE(multigrid_mesh) :: mg_mesh
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
CALL multigrid_construct(mg_mesh)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!--- Lagrange
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange,minlev=minlev)
CALL lag_setup_interp(ML_oft_lagrange)
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_oft_h1,minlev=minlev)
CALL h1_setup_interp(ML_oft_h1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl,minlev=minlev)
CALL hcurl_setup_interp(ML_oft_hcurl)
!--- Full H(Curl) space
CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,minlev)
CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_oft_h1)
!---------------------------------------------------------------------------
! Create Lagrange metric solver
!---------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
CALL ML_oft_lagrange%vec_create(n)
CALL ML_oft_lagrange%vec_create(dn)
CALL ML_oft_lagrange%vec_create(ni)
CALL ML_oft_lagrange%vec_create(temp)
CALL ML_oft_lagrange%vec_create(dtemp)
CALL ML_oft_lagrange%vec_create(ti)
!---------------------------------------------------------------------------
! Set dn from sound wave init
!---------------------------------------------------------------------------
sound_field%mesh=>mg_mesh%mesh
sound_field%delta=delta
sound_field%field='n'
CALL oft_lag_project(ML_oft_lagrange%current_level,sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL n%set(1.d0)
CALL dn%add(0.d0,1.d0,u,-1.d0,n)
CALL ni%add(0.d0,1.d0,dn)
!---------------------------------------------------------------------------
! Set dt from sound wave init
!---------------------------------------------------------------------------
sound_field%field='t'
CALL oft_lag_project(ML_oft_lagrange%current_level,sound_field,v)
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
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,mop,"none")
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
CALL ML_oft_vlagrange%vec_create(vel)
CALL ML_oft_vlagrange%vec_create(dvel)
CALL ML_oft_vlagrange%vec_create(vi)
!---------------------------------------------------------------------------
! Set dV from sound wave init
!---------------------------------------------------------------------------
sound_field%field='v'
CALL oft_lag_vproject(ML_oft_lagrange%current_level,sound_field,v)
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
CALL ML_hcurl_grad%vec_create(b)
CALL ML_hcurl_grad%vec_create(db)
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
    CALL ML_oft_lagrange%vec_create(tempe)
    CALL tempe%add(0.d0,1.d0,temp)
    CALL ML_oft_lagrange%vec_create(dtempe)
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
    CALL ML_oft_lagrange%vec_create(tempe)
    CALL tempe%add(0.d0,1.d0,temp)
    equil_fields%Te=>tempe
  END IF
  equil_fields%B=>b
  equil_fields%V=>vel
  equil_fields%Ne=>n
  equil_fields%Ti=>temp
  CALL xmhd_run(equil_fields)
END IF
CALL ML_oft_lagrange%vec_create(u)
CALL u%set(1.d0)
!---Compare density waveform
sound_field%field='n'
sound_field%diff=.TRUE.
nierr=scal_energy(mg_mesh%mesh,sound_field,order*2)
err_field%dim=1
err_field%a=>sound_field
err_field%b=>sfield
CALL n%scale(1.d0/N0)
CALL n%add(1.d0,-1.d0,u)
sfield%u=>n
CALL sfield%setup(ML_oft_lagrange%current_level)
nerr=scal_energy(mg_mesh%mesh,err_field,order*2)
!---Compare temperature waveform
sound_field%field='t'
sound_field%diff=.TRUE.
tierr=scal_energy(mg_mesh%mesh,sound_field,order*2)
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
CALL sfield%setup(ML_oft_lagrange%current_level)
terr=scal_energy(mg_mesh%mesh,err_field,order*2)
!---Compare velocity waveform
sound_field%field='v'
sound_field%diff=.FALSE.
vierr=vec_energy(mg_mesh%mesh,sound_field,order*2)
err_field%dim=3
err_field%a=>sound_field
err_field%b=>vfield
CALL vel%scale(1.d0/v_delta)
vfield%u=>vel
CALL vfield%setup(ML_oft_lagrange%current_level)
verr=vec_energy(mg_mesh%mesh,err_field,order*2)
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
