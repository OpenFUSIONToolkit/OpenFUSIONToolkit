MODULE taylor_green_helpers
USE oft_local
IMPLICIT NONE
REAL(r8) :: rho_ana
CONTAINS


SUBROUTINE vx_init(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = SIN(2.d0*pi*pt(1))*COS(2.d0*pi*pt(2))
END SUBROUTINE vx_init

SUBROUTINE vz_init(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -COS(2.d0*pi*pt(1))*SIN(2.d0*pi*pt(2))
END SUBROUTINE vz_init

SUBROUTINE p_init(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = 0.25d0*rho_ana*(COS(4.d0*pi*pt(1)) + COS(4.d0*pi*pt(2)))
END SUBROUTINE p_init

END MODULE taylor_green_helpers

PROGRAM test_taylor_green_2d
!---------------------------------------------------------------------------
! Taylor-Green vortex test for the incompressible 2D xMHD solver.
!
! The doubly-periodic unit square is initialized with
!     v_x =  sin(2*pi*x) cos(2*pi*y)
!     v_z = -cos(2*pi*x) sin(2*pi*y)
! and no magnetic field. For an incompressible fluid, there is an exact
! solution of Navier-Stokes that decays in place:
!     v(t) = v(0) * exp(-8*pi^2*nu*t)
!     p(t) = (rho/4)*[cos(4*pi*x) + cos(4*pi*y)] * exp(-16*pi^2*nu*t)
!
! NOTE: the incompressible solver defaults to a direct LU preconditioner, which
! cannot factorize the periodic Jacobian . This test must therefore be run 
!  with an XML preconditioner definition wrapping the LU in block-Jacobi
!---------------------------------------------------------------------------
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
USE oft_blag_operators, ONLY: oft_blag_getmop, oft_blag_project, oft_lag_brinterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: proton_mass, mu0
USE diagnostic, ONLY: scal_energy_2d
USE fem_utils, ONLY: diff_interp_2d
USE xmhd_2d
USE taylor_green_helpers, ONLY: rho_ana, vx_init, vz_init, p_init
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_lag_brinterp), TARGET :: ana_field, final_field
TYPE(oft_lag_brinterp), TARGET :: pana_field, pfinal_field
TYPE(diff_interp_2d) :: err_field, perr_field
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL(), minv_p => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL(), mop_p => NULL()
CLASS(oft_vector), POINTER :: u,v,vx_ana,vz_ana,tmp
CLASS(oft_vector), POINTER :: up,vp,p_ana,p_num,ones
REAL(r8), POINTER :: pvec_vals(:)
!---Errors
REAL(r8) :: vxerr_i,vzerr_i,vxerr,vzerr,decay
REAL(r8) :: perr_i,perr,decay_p
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 20
INTEGER(i4) :: rst_freq = 1000
INTEGER(i4) :: ittarget = 1000
REAL(r8) :: n0 = 3.d24
REAL(r8) :: nu = 1.d-2
REAL(r8) :: eta = 1.d0
REAL(r8) :: chi = 1.d0
REAL(r8) :: gamma = 1.67d0
REAL(r8) :: D_diff = 1.d0
REAL(r8) :: m_i = 205.d0
REAL(r8) :: mu = mu0
REAL(r8) :: den_scale = 1.d24
REAL(r8) :: dt = 5.d-3
REAL(r8) :: lin_tol = 1.d-9
REAL(r8) :: nl_tol = 1.d-7
LOGICAL :: pm = .FALSE.
LOGICAL :: use_mfnk = .TRUE.
LOGICAL :: timestep_cn = .TRUE.

NAMELIST/xmhd_options/order,chi,eta,nu,gamma,D_diff,m_i,mu,dt,nsteps,rst_freq, &
use_mfnk,pm,n0,den_scale,lin_tol,nl_tol,timestep_cn,ittarget
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
mhd_sim%incomp = .TRUE. !Turn on incompressible flag
CALL mhd_sim%setup(mg_mesh, order)

!---------------------------------------------------------------------------
! Set initial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrices for both FE spaces (pressure space is lower order)
NULLIFY(u,v,mop,vec_vals)
NULLIFY(up,vp,mop_p,pvec_vals)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop)
CALL oft_blag_getmop(ML_oft_blagrange_p%current_level,mop_p)
!---Setup linear solvers
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-2
CALL create_diag_pre(minv%pre)
CALL create_cg_solver(minv_p)
minv_p%A=>mop_p
minv_p%its=-2
CALL create_diag_pre(minv_p%pre)
!---Create fields for solvers
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
CALL ML_oft_blagrange%vec_create(vx_ana)
CALL ML_oft_blagrange%vec_create(vz_ana)
CALL ML_oft_blagrange%vec_create(tmp)
CALL ML_oft_blagrange_p%vec_create(up)
CALL ML_oft_blagrange_p%vec_create(vp)
CALL ML_oft_blagrange_p%vec_create(p_ana)
CALL ML_oft_blagrange_p%vec_create(p_num)
CALL ML_oft_blagrange_p%vec_create(ones)

!---Density is uniform and held fixed
mhd_sim%den_scale = den_scale
CALL u%set(n0)
CALL u%get_local(vec_vals)
vec_vals = vec_vals / den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)

!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>vx_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,2)
CALL vx_ana%restore_local(vec_vals)

!---No out-of-plane flow
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
field_init%func=>vz_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,4)
CALL vz_ana%restore_local(vec_vals)

!---Project pressure initial condition onto scalar Lagrange pressure basis
rho_ana = m_i*proton_mass*n0
field_init%func=>p_init
CALL oft_blag_project(ML_oft_blagrange_p%current_level,field_init,vp)
CALL up%set(0.d0)
CALL minv_p%apply(up,vp)
CALL up%get_local(pvec_vals)
CALL mhd_sim%u%restore_local(pvec_vals,5)
CALL p_ana%restore_local(pvec_vals)

!---No magnetic field
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,6)
CALL mhd_sim%u%restore_local(vec_vals,7)

!---Cleanup objects used for projection
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
DEALLOCATE(minv%pre)
CALL minv%delete
DEALLOCATE(minv)
CALL up%delete
CALL vp%delete
CALL mop_p%delete
DEALLOCATE(up,vp,mop_p)
CALL minv_p%pre%delete
DEALLOCATE(minv_p%pre)
CALL minv_p%delete
DEALLOCATE(minv_p)

!---------------------------------------------------------------------------
! Boundary conditions
!---------------------------------------------------------------------------
!Turn off magnetic field evolution, defaults are fine for other fields
NULLIFY(mhd_sim%psi_bc,mhd_sim%by_bc)
ALLOCATE(mhd_sim%psi_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%by_bc(oft_blagrange%ne))
mhd_sim%psi_bc = .TRUE.
mhd_sim%by_bc  = .TRUE.

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%m_i=m_i*proton_mass
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
mhd_sim%lin_tol=lin_tol
mhd_sim%nl_tol=nl_tol
mhd_sim%timestep_cn=timestep_cn
mhd_sim%ittarget=ittarget
oft_env%pm=pm

CALL mhd_sim%run_simulation()

!---------------------------------------------------------------------------
! Compare to the analytic solution at the final time
!---------------------------------------------------------------------------
decay   = EXP( -8.d0*(pi**2)*nu*mhd_sim%t)
decay_p = EXP(-16.d0*(pi**2)*nu*mhd_sim%t)
CALL vx_ana%scale(decay)
CALL vz_ana%scale(decay)
CALL p_ana%scale(decay_p)
WRITE(*,'(A,ES14.6,A,ES14.6)')'Final time = ',mhd_sim%t,'   decay factor = ',decay

!---The velocity and pressure live on different spaces, so each needs its own
!   pair of interpolants and difference object
err_field%dim=1
err_field%a=>ana_field
err_field%b=>final_field
perr_field%dim=1
perr_field%a=>pana_field
perr_field%b=>pfinal_field

!---v_x error
CALL mhd_sim%u%get_local(vec_vals,2)
CALL tmp%restore_local(vec_vals)
ana_field%u=>vx_ana
final_field%u=>tmp
CALL ana_field%setup(ML_oft_blagrange%current_level)
CALL final_field%setup(ML_oft_blagrange%current_level)
vxerr_i=scal_energy_2d(mg_mesh%smesh,ana_field,order*2)
vxerr  =scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!---v_z error
CALL mhd_sim%u%get_local(vec_vals,4)
CALL tmp%restore_local(vec_vals)
ana_field%u=>vz_ana
final_field%u=>tmp
CALL ana_field%setup(ML_oft_blagrange%current_level)
CALL final_field%setup(ML_oft_blagrange%current_level)
vzerr_i=scal_energy_2d(mg_mesh%smesh,ana_field,order*2)
vzerr  =scal_energy_2d(mg_mesh%smesh,err_field,order*2)

!---Pressure error (both analytic and numerical pressures are shifted so their minimum is zero,
!--- because incompressible pressure only defined up to a constant)
CALL mhd_sim%u%get_local(pvec_vals,5)
CALL p_num%restore_local(pvec_vals)
CALL ones%set(1.d0)
CALL p_ana%get_local(pvec_vals)
CALL p_ana%add(1.d0,oft_mpi_max(-MINVAL(pvec_vals)),ones)
CALL p_num%get_local(pvec_vals)
CALL p_num%add(1.d0,oft_mpi_max(-MINVAL(pvec_vals)),ones)
pana_field%u=>p_ana
pfinal_field%u=>p_num
CALL pana_field%setup(ML_oft_blagrange_p%current_level)
CALL pfinal_field%setup(ML_oft_blagrange_p%current_level)
perr_i=scal_energy_2d(mg_mesh%smesh,pana_field,order*2)
perr  =scal_energy_2d(mg_mesh%smesh,perr_field,order*2)

!---Write errors to file
OPEN(NEWUNIT=io_unit,FILE='taylor_green_2d.results')
WRITE(io_unit,*)SQRT(vxerr/vxerr_i)
WRITE(io_unit,*)SQRT(vzerr/vzerr_i)
WRITE(io_unit,*)SQRT(perr/perr_i)
CLOSE(io_unit)

!---Finalize environment
CALL oft_finalize


END PROGRAM test_taylor_green_2d
