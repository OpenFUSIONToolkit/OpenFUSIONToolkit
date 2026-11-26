!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_rt_2d.F90
!
!> Regression test for gravity in the 2D xMHD module. An exponentially stratified plasma
!! (n0 = n_0*exp(z/L)) with a constant-Alfven-speed field along the periodic direction
!! is placed in a uniform gravitational field along -z between two conducting, free-slip
!! walls, and seeded with the incompressible m=1 Rayleigh-Taylor eigenmode in velocity.
!! The kinetic energy history (`xmhd_2d.hist`) is used by `test_rt_2d.py` to measure the
!! growth rate (Chandrasekhar, Hydrodynamic and Hydromagnetic Stability, 1961, sections
!! 92 and 97).
!!
!! @ingroup testing
!------------------------------------------------------------------------------
MODULE rt_2d_helpers
USE oft_local
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
IMPLICIT NONE
REAL(r8) :: n_0 = 1.d19 !< Density at z=0 [m^-3]
REAL(r8) :: T_0 = 1.d2 !< Temperature at z=0 [eV]
REAL(r8) :: L_n = 0.5d0 !< Density scale length [m] (> 0 is unstable)
REAL(r8) :: g = 1.d8 !< Gravitational acceleration along -z [m/s^2]
REAL(r8) :: B_0 = 0.d0 !< Field along x at z=0 [T]
REAL(r8) :: V1 = 1.d0 !< Seed amplitude [m/s]
REAL(r8) :: L_x = 1.d0 !< Periodic length in x [m]
REAL(r8) :: L_z = 1.d0 !< Distance between the walls [m]
CONTAINS
!> Equilibrium density n0(z), normalized by n_0
SUBROUTINE n_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = EXP(pt(2)/L_n)
END SUBROUTINE n_eq
!> Equilibrium pressure from total force balance d(p + Bx^2/(2 mu0))/dz = -m_i n0 g
FUNCTION p_eq(z) RESULT(p)
REAL(r8), INTENT(in) :: z
REAL(r8) :: p
p = 2.d0*elec_charge*n_0*T_0 - (proton_mass*n_0*g*L_n + B_0**2/(2.d0*mu0))*(EXP(z/L_n)-1.d0)
END FUNCTION p_eq
!> Equilibrium temperature [eV], with p = 2 n e T
SUBROUTINE T_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = p_eq(pt(2))/(2.d0*elec_charge*n_0*EXP(pt(2)/L_n))
END SUBROUTINE T_eq
!> Equilibrium flux, Bx = -dpsi/dz = B_0*exp(z/(2 L))
SUBROUTINE psi_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -2.d0*L_n*B_0*EXP(pt(2)/(2.d0*L_n))
END SUBROUTINE psi_eq
!> Seed vz: incompressible m=1 eigenmode
SUBROUTINE vz_seed(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
REAL(r8) :: k
k = 2.d0*pi/L_x
val = V1*EXP(-pt(2)/(2.d0*L_n))*COS(pi*pt(2)/L_z)*COS(k*pt(1))
END SUBROUTINE vz_seed
!> Seed vx, from div(v) = 0
SUBROUTINE vx_seed(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
REAL(r8) :: k
k = 2.d0*pi/L_x
val = (V1/k)*EXP(-pt(2)/(2.d0*L_n))*(COS(pi*pt(2)/L_z)/(2.d0*L_n) &
  + (pi/L_z)*SIN(pi*pt(2)/L_z))*SIN(k*pt(1))
END SUBROUTINE vx_seed
END MODULE rt_2d_helpers

PROGRAM test_rt_2d
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_blag_operators, ONLY: oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE xmhd_2d
USE rt_2d_helpers
IMPLICIT NONE
INTEGER(i4) :: io_unit, ierr, i
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u, v
REAL(r8) :: v_A_crit, B_crit, p_min
!---Runtime options
INTEGER(i4) :: order = 3
INTEGER(i4) :: nsteps = 100
LOGICAL :: linear = .FALSE.
LOGICAL :: use_mfnk = .FALSE.
REAL(r8) :: b_frac = 0.d0 !< B_0 as a fraction of the critical field for the m=1 mode
REAL(r8) :: dt = 2.d-6
REAL(r8) :: nu = 1.d0
REAL(r8) :: eta = 1.d0
REAL(r8) :: chi = 0.d0 ! Diffusion of n and T moves the equilibrium off hydrostatic balance
REAL(r8) :: D_diff = 0.d0
REAL(r8) :: lin_tol = 1.d-10
REAL(r8) :: nl_tol = 1.d-8
NAMELIST/rt_2d_options/order,nsteps,linear,use_mfnk,n_0,T_0,L_n,g,b_frac,V1,L_x,L_z, &
  dt,nu,eta,chi,D_diff,lin_tol,nl_tol
!------------------------------------------------------------------------------
! Initialize environment
!------------------------------------------------------------------------------
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,rt_2d_options,IOSTAT=ierr)
CLOSE(io_unit)
! Critical constant Alfven speed for the m=1 mode (Chandrasekhar 1961, section 97); ABS gives a
! reference field for stable (L_n < 0) stratifications
v_A_crit = SQRT(ABS(g/L_n)/((2.d0*pi/L_x)**2 + (pi/L_z)**2 + 1.d0/(4.d0*L_n**2)))
B_crit = v_A_crit*SQRT(mu0*proton_mass*n_0)
B_0 = b_frac*B_crit
p_min = MIN(p_eq(-L_z/2.d0), p_eq(L_z/2.d0))
IF(p_min<=0.d0)CALL oft_abort("Equilibrium pressure is not positive","test_rt_2d",__FILE__)
!---------------------------------------------------------------------------
! Setup grid and FE representation
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
! Free-slip walls: no Dirichlet condition on vx (the walls are the only boundary)
NULLIFY(mhd_sim%velx_bc)
ALLOCATE(mhd_sim%velx_bc(oft_blagrange%ne))
mhd_sim%velx_bc = .FALSE.
! Natural (zero-flux) conditions on n and T: with vz = 0 on a wall, dn/dt = -n dvz/dz is not zero
! there, and Dirichlet conditions on n and T in a stratified plasma drive a spurious instability
! growing at a rate ~c_s/L
NULLIFY(mhd_sim%n_bc,mhd_sim%T_bc)
ALLOCATE(mhd_sim%n_bc(oft_blagrange%ne))
mhd_sim%n_bc = .FALSE.
mhd_sim%T_bc => mhd_sim%n_bc
!---------------------------------------------------------------------------
! Project equilibrium and seed
!---------------------------------------------------------------------------
NULLIFY(mop,vec_vals)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop)
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-2
CALL create_diag_pre(minv%pre)
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
field_init%mesh=>mesh
CALL mhd_sim%u0%set(0.d0)
CALL mhd_sim%u%set(0.d0)
mhd_sim%den_scale = n_0
! Density (stored normalized by den_scale = n_0)
field_init%func=>n_eq
CALL project_field(1, mhd_sim%u0)
! Temperature
field_init%func=>T_eq
CALL project_field(5, mhd_sim%u0)
! Flux
field_init%func=>psi_eq
CALL project_field(6, mhd_sim%u0)
! Seed velocity, with vz exactly zero on the walls
field_init%func=>vx_seed
CALL project_field(2, mhd_sim%u)
field_init%func=>vz_seed
CALL project_field(4, mhd_sim%u, zero_bc=.TRUE.)
IF(.NOT.linear)CALL mhd_sim%u%add(1.d0,1.d0,mhd_sim%u0)
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
DEALLOCATE(minv%pre)
CALL minv%delete
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%g_accel = g
mhd_sim%chi = chi
mhd_sim%eta = eta
mhd_sim%nu = nu
mhd_sim%D_diff = D_diff
mhd_sim%gamma = 5.d0/3.d0
mhd_sim%dt = dt
mhd_sim%nsteps = nsteps
mhd_sim%rst_freq = nsteps
mhd_sim%mfnk = use_mfnk
mhd_sim%linear = linear
mhd_sim%lin_tol = lin_tol
mhd_sim%nl_tol = nl_tol
mhd_sim%ittarget = 1000 ! Keep dt fixed
mhd_sim%timestep_cn = .TRUE.
oft_env%pm = .FALSE.
IF(linear)THEN
  CALL mhd_sim%run_lin_simulation()
ELSE
  CALL mhd_sim%run_simulation()
END IF
CALL oft_finalize
CONTAINS
!> Project field_init onto the Lagrange basis and store it in field `ifield` of `uout`
SUBROUTINE project_field(ifield, uout, zero_bc)
INTEGER(i4), INTENT(in) :: ifield
CLASS(oft_vector), INTENT(inout) :: uout
LOGICAL, OPTIONAL, INTENT(in) :: zero_bc
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vec_vals)
IF(PRESENT(zero_bc))THEN
  IF(zero_bc)THEN
    DO i=1,oft_blagrange%ne
      IF(oft_blagrange%global%gbe(i))vec_vals(i)=0.d0
    END DO
  END IF
END IF
CALL uout%restore_local(vec_vals,ifield)
END SUBROUTINE project_field
END PROGRAM test_rt_2d
