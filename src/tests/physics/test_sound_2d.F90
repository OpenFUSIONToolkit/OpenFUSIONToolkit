PROGRAM test_sound_2d
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project, oft_lag_brinterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE diagnostic, ONLY: scal_energy_2d
USE fem_utils, ONLY: diff_interp_2d
USE xmhd_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_lag_brinterp), TARGET ::  ic_field, final_field
TYPE(diff_interp_2d) :: err_field
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,v_ic,n_ic,T_ic,tmp,uniform
!--Errors
REAL(r8) :: verr_i,nerr_i,Terr_i,verr,nerr,Terr,v_delta
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: dt = 5.d-7
REAL(r8) :: n0 = 1.d0
REAL(r8) :: velx0 = 1.d0
REAL(r8) :: vely0 = 1.d0
REAL(r8) :: velz0 = 1.d0
REAL(r8) :: t0 = 1.d0
REAL(r8) :: psi0 = 1.d0
REAL(r8) :: by0 = 1.d0
REAL(r8) :: bx0=0.d0
REAL(r8) :: bz0=0.d0
REAL(r8) :: chi=1.E-12 !< Needs docs
REAL(r8) :: eta=1.d0 !< Needs docs
REAL(r8) :: nu=1.E-12 !< Needs docs
REAL(r8) :: gamma=1.d0
REAL(r8) :: D_diff=1.E-12
REAL(r8) :: den_scale=1.d19
REAL(r8) :: k_dir(3) = (/1.d0,0.d0,0.d0/) !< Direction of wave propogation
REAL(r8) :: r0(3) = (/0.d0,0.d0,0.d0/)  !< Zero-phase position
REAL(r8) :: lam = 2.d0 !< Wavelength
REAL(r8) :: v_sound = 2.d4
REAL(r8) :: delta = 1.d-4 !< Relative size of perturbation (<<1)
REAL(r8) :: lin_tol = 1.d-8
REAL(r8) :: nl_tol = 1.d-5
LOGICAL :: pm=.FALSE.
LOGICAL :: linear=.FALSE.
LOGICAL :: cyl=.FALSE.
LOGICAL :: use_mfnk=.FALSE.

NAMELIST/xmhd_options/order,linear,cyl,chi,eta,nu,gamma,D_diff,dt,nsteps,rst_freq, &
use_mfnk,pm,n0,psi0,velx0,vely0,velz0,t0,by0,den_scale,bx0,bz0,lin_tol,nl_tol
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)

!---------------------------------------------------------------------------
! Set intial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop) ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Set constant values
t0=(v_sound**2)*3.d0*proton_mass/(5.d0*elec_charge*2.d0)
v_delta=2.d0*t0*elec_charge/(proton_mass*v_sound)
!v_delta = 2.d0/(v_sound*proton_mass*n0)


!---Project n initial condition onto scalar Lagrange basis
IF (linear)THEN
  CALL u%set(n0)
  CALL u%get_local(vec_vals)
  vec_vals = vec_vals / den_scale
  CALL mhd_sim%u0%restore_local(vec_vals,1)
END IF

field_init%func=>dens_sound
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
vec_vals = vec_vals / den_scale
mhd_sim%den_scale = den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)
CALL ML_oft_blagrange%vec_create(n_ic)
CALL n_ic%restore_local(vec_vals)


!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>velx_sound
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(v_delta)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,2)
CALL ML_oft_blagrange%vec_create(v_ic)
CALL v_ic%restore_local(vec_vals)

! !---Project v_z initial condition onto scalar Lagrange basis
! field_init%func=>velz_sound
! CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
! CALL u%set(0.d0)
! CALL minv%apply(u,v)
! CALL u%scale(v_delta)
! CALL u%get_local(vec_vals)
! CALL mhd_sim%u%restore_local(vec_vals,4)


!---Project T initial condition onto scalar Lagrange basis
! First save background constant temperature field
IF(linear)THEN
  CALL u%set(t0)
  CALL u%get_local(vec_vals)
  CALL mhd_sim%u0%restore_local(vec_vals,5)
END IF

field_init%func=> temp_sound
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(t0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,5)
CALL ML_oft_blagrange%vec_create(t_ic)
CALL t_ic%restore_local(vec_vals)

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%gamma=5.d0/3.d0 !gamma
mhd_sim%D_diff=D_diff
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
mhd_sim%linear = linear
mhd_sim%cyl_flag = cyl
mhd_sim%lin_tol = lin_tol
mhd_sim%nl_tol = nl_tol
mhd_sim%ittarget=300
mhd_sim%timestep_cn=.TRUE.
oft_env%pm=pm

IF (linear) THEN
  CALL mhd_sim%u%add(1.d0,-1.d0,mhd_sim%u0)
  CALL mhd_sim%run_lin_simulation()
  CALL mhd_sim%u%add(1.d0,1.d0,mhd_sim%u0)
ELSE 
  CALL mhd_sim%run_simulation()
END IF

CALL xmhd_2d_plot(mhd_sim)

CALL ML_oft_blagrange%vec_create(uniform)
CALL ML_oft_blagrange%vec_create(tmp)

!---Compare density waveform
CALL uniform%set(n0/den_scale)
CALL n_ic%add(1.d0,-1.d0,uniform)
ic_field%u=>n_ic
CALL ic_field%setup(ML_oft_blagrange%current_level)
nerr_i=scal_energy_2d(mg_mesh%smesh,ic_field,order*2)
err_field%dim=1
err_field%a=>ic_field
err_field%b=>final_field
CALL mhd_sim%u%get_local(vec_vals,1)
CALL tmp%restore_local(vec_vals)
CALL tmp%add(1.d0,-1.d0,uniform)
final_field%u=>tmp
CALL final_field%setup(ML_oft_blagrange%current_level)
nerr=scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!---Compare temperature waveform
CALL uniform%set(t0)
CALL t_ic%add(1.d0,-1.d0,uniform)
ic_field%u=>t_ic
CALL ic_field%setup(ML_oft_blagrange%current_level)
terr_i=scal_energy_2d(mg_mesh%smesh,ic_field,order*2)
err_field%dim=1
err_field%a=>ic_field
err_field%b=>final_field
CALL mhd_sim%u%get_local(vec_vals,5)
CALL tmp%restore_local(vec_vals)
CALL tmp%add(1.d0,-1.d0,uniform)
final_field%u=>tmp
CALL final_field%setup(ML_oft_blagrange%current_level)
terr=scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!---Compare velocity waveform
ic_field%u=>v_ic
CALL ic_field%setup(ML_oft_blagrange%current_level)
verr_i=scal_energy_2d(mg_mesh%smesh,ic_field,order*2)
err_field%dim=1
err_field%a=>ic_field
err_field%b=>final_field
CALL mhd_sim%u%get_local(vec_vals,2)
CALL tmp%restore_local(vec_vals)
final_field%u=>tmp
CALL final_field%setup(ML_oft_blagrange%current_level)
verr=scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!Write errors to file
OPEN(NEWUNIT=io_unit,FILE='sound_2d.results')
WRITE(io_unit,*)SQRT(verr/verr_i)
WRITE(io_unit,*)SQRT(nerr/nerr_i)
WRITE(io_unit,*)SQRT(Terr/Terr_i)
CLOSE(io_unit)

!---Finalize enviroment
CALL oft_finalize
CONTAINS
    
SUBROUTINE dens_sound(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val=(1.d0+delta*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam))**(3.d0/5.d0)
END SUBROUTINE dens_sound
    
SUBROUTINE velx_sound(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = delta*k_dir(1)*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)
END SUBROUTINE velx_sound
    
SUBROUTINE velz_sound(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = delta*k_dir(2)*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)
END SUBROUTINE velz_sound
    
SUBROUTINE temp_sound(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val=(1.d0+delta*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam))**(2.d0/5.d0)
END SUBROUTINE temp_sound

END PROGRAM test_sound_2d