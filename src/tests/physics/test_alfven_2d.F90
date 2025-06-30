PROGRAM test_alfven_2d
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
TYPE(oft_lag_brinterp), TARGET :: initial, final
TYPE(diff_interp_2d) :: err_field
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v, vx_ic, psi_ic, tmp
!--Errors
REAL(r8) :: verr_i, psierr_i, verr, psierr
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: dt = 5.d-7
REAL(r8) :: n0 = 1.d19
REAL(r8) :: velx0 = 0.d0
REAL(r8) :: vely0 = 0.d0
REAL(r8) :: velz0 = 0.d0
REAL(r8) :: t0 = 1.d3
REAL(r8) :: psi0 = 0.d0
REAL(r8) :: by0 = 0.d0
REAL(r8) :: bx0=0.d0
REAL(r8) :: bz0=0.d0
REAL(r8) :: chi=1.E-12 !< Needs docs
REAL(r8) :: eta=1.d12 !< Needs docs
REAL(r8) :: nu=1.E-12 !< Needs docs
REAL(r8) :: D_diff=1.E-12
REAL(r8) :: gamma=1.67d0
REAL(r8) :: den_scale=1.d19
REAL(r8) :: B_0(3)=0.d0
REAL(r8) :: k_dir(3) = (/0.d0,1.d0,0.d0/) !< Direction of wave propogation
Real(r8) :: v_dir(3) = (/1.d0,0.d0,0.d0/) !<Direction of velocity perturbation
REAL(r8) :: r0(3) = (/0.d0,0.d0,0.d0/)  !< Zero-phase position
REAL(r8) :: lam = 2.d0 !< Wavelength
REAL(r8) :: v_alf = 1.d4 !< Alfven speed
REAL(r8) :: v_delta = 1.d-4 !< Relative size of perturbation (<<1)
REAL(r8) :: B !<Background magnetic field magnitude
REAL(r8) :: B_delta !<Perturbed magnetic field magnitude
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/xmhd_options/order,nsteps, rst_freq, dt, n0, psi0, velx0,&
vely0,velz0, t0, by0, bx0, bz0,chi,eta,nu, D_diff, gamma, den_scale, use_mfnk,pm
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
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
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Set constant values
B=v_alf*SQRT(mu0*n0*proton_mass)
B_delta = v_delta*B/v_alf

!---Project n initial condition onto scalar Lagrange basis
field_init%func=>const_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'n0')
mhd_sim%den_scale = den_scale
vec_vals = vec_vals / den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)

!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>velx_alf
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(v_delta)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vx0')
CALL mhd_sim%u%restore_local(vec_vals,2)
CALL ML_oft_blagrange%vec_create(vx_ic)
CALL vx_ic%restore_local(vec_vals)

!---Project v_y initial condition onto scalar Lagrange basis
field_init%func=>vely_alf
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(v_delta)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vy0')
CALL mhd_sim%u%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
field_init%func=>velz_alf
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(v_delta)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vz0')
CALL mhd_sim%u%restore_local(vec_vals,4)

!---Project T initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(t0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'T0')
CALL mhd_sim%u%restore_local(vec_vals,5)

!---Project psi initial condition onto scalar Lagrange basis
field_init%func=>psi_alf
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(B_delta)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(vec_vals,6)
CALL ML_oft_blagrange%vec_create(psi_ic)
CALL psi_ic%restore_local(vec_vals)

!---Project by initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(by0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'by0')
CALL mhd_sim%u%restore_local(vec_vals,7)

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
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
B_0(3) = B
mhd_sim%B_0 = B_0
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
oft_env%pm=pm
CALL mhd_sim%run_simulation()

CALL ML_oft_blagrange%vec_create(tmp)

!---Compare vx waveform
initial%u => vx_ic
CALL initial%setup(ML_oft_blagrange%current_level)
verr_i=scal_energy_2d(mg_mesh%smesh,initial,order*2)
err_field%dim=1
err_field%a=>initial
err_field%b=>final
CALL mhd_sim%u%get_local(vec_vals,2)
CALL tmp%restore_local(vec_vals) !this line is the problem
final%u=>tmp
CALL final%setup(ML_oft_blagrange%current_level)
verr=scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!---Compare psi waveform
initial%u => psi_ic
CALL initial%setup(ML_oft_blagrange%current_level)
psierr_i=scal_energy_2d(mg_mesh%smesh,initial,order*2)
err_field%dim=1
err_field%a=>initial
err_field%b=>final
CALL mhd_sim%u%get_local(vec_vals,6)
CALL tmp%restore_local(vec_vals)
final%u=>tmp
CALL final%setup(ML_oft_blagrange%current_level)
psierr=scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!Write errors to file
OPEN(NEWUNIT=io_unit,FILE='alfven_2d.results')
WRITE(io_unit,*)SQRT(verr/verr_i)
WRITE(io_unit,*)SQRT(psierr/psierr_i)
CLOSE(io_unit)

!---Finalize enviroment
CALL oft_finalize
CONTAINS

SUBROUTINE const_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = 1.0
END SUBROUTINE const_init

SUBROUTINE psi_alf(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = COS(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)*lam/(2.d0*pi)
END SUBROUTINE psi_alf

SUBROUTINE velx_alf(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = v_dir(1)*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)
END SUBROUTINE velx_alf

SUBROUTINE vely_alf(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = v_dir(2)*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)
END SUBROUTINE vely_alf

SUBROUTINE velz_alf(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = v_dir(3)*SIN(DOT_PRODUCT(pt-r0,k_dir)*2.d0*pi/lam)
END SUBROUTINE velz_alf


END PROGRAM test_alfven_2d