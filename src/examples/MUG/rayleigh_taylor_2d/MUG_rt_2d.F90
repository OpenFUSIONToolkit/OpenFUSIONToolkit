!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: 2D Rayleigh-Taylor instability with and without a magnetic field    {#doc_mug2d_rt_ex}
!!============================
!!
!![TOC]
!!
!! This example uses the \ref xmhd_2d "2D extended MHD" module to simulate the Rayleigh-Taylor (RT)
!! instability of a plasma supported against gravity, with and without a magnetic field along the
!! interchange direction. It has two equilibria:
!!  - An exponential profile (`profile='exp'`), whose growth rate is known exactly, with and without a field
!!    ([S. Chandrasekhar, Hydrodynamic and Hydromagnetic Stability (1961)](https://store.doverpublications.com/products/9780486640716),
!!    sections 92 and 97). It validates the linearized model over a scan in field strength (\ref doc_mug2d_rt_ex_scan).
!!  - A transition layer (`profile='tanh'`) between two uniform plasmas, used for nonlinear runs. Its growth rate in
!!    the nonlinear model is checked against the linearized model on the same mesh, and the run is then followed
!!    into the bubble-and-spike phase (\ref doc_mug2d_rt_ex_nl).
!!
!! The exponential profile is not used for nonlinear runs (the driver stops): its largest gradients sit at the walls,
!! where they make the late nonlinear phase and any particle diffusion or thermal conduction difficult
!! (\ref doc_mug2d_rt_ex_exp_nl).
!!
!!\section doc_mug2d_rt_ex_setup Problem setup
!!
!! The 2D solver works in the plane of the mesh, with the first mesh coordinate as \f$ x \f$ and the second as
!! the vertical coordinate \f$ z \f$; \f$ y \f$ is the ignorable direction. Gravity acts along \f$ -\hat{z} \f$
!! (\ref xmhd_2d::oft_xmhd_2d_sim::g_accel "g_accel" \f$ = g > 0 \f$). The domain is
!! \f$ x \in [-L_x/2, L_x/2] \f$ (periodic) by \f$ z \in [-L_z/2, L_z/2] \f$, bounded by rigid, perfectly
!! conducting, free-slip walls. The pressure is \f$ p = 2 n e T \f$ (equal ion and electron temperatures, \f$ T \f$ in eV),
!! and in both cases it follows from total force balance, \f$ \frac{d}{dz} \left( p + B_x^2 / 2 \mu_0 \right) = -m_i n_0 g \f$,
!! with \f$ p_0 = 2 e n_0 T_0 \f$ at \f$ z = 0 \f$. The driver stops if \f$ p \le 0 \f$ anywhere in the domain.
!! A uniform out-of-plane field \f$ B_y = b_{y0} \f$ can be added to either case.
!!
!! The perturbation is applied to the velocity only, with a vertical structure \f$ f(z) \f$ that vanishes on the walls:
!! \f[ v_z = V_1 f(z) \cos(k x), \qquad v_x = -\frac{V_1}{k} f'(z) \sin(k x), \qquad k = 2 \pi / L_x \f]
!! It is divergence free and launches no sound waves. Since the initial displacement is zero, it is an equal mix of
!! growing and decaying modes.
!!
!!\subsection doc_mug2d_rt_ex_exp Exponential profile (linear benchmark)
!!
!! \f[ n_0(z) = n_0 e^{z/L}, \qquad B_x(z) = B_0 e^{z/2L}, \qquad \psi_0(z) = -2 L B_0 e^{z/2L} \f]
!! with the heavy plasma on top for \f$ L > 0 \f$, \f$ B_x = -\partial_z \psi \f$, and a constant Alfven speed
!! \f$ v_A = B_0 / \sqrt{\mu_0 m_i n_0} \f$. The pressure is
!! \f[ p(z) = p_0 - \left( m_i n_0 g L + \frac{B_0^2}{2 \mu_0} \right) \left( e^{z/L} - 1 \right) \f]
!! The seed is the \f$ m = 1 \f$ incompressible eigenmode, \f$ f(z) = e^{-z/2L} \cos(\pi z / L_z) \f$, so the velocity
!! grows as \f$ \cosh(\gamma t) \f$. For an incompressible, ideal plasma (Chandrasekhar 1961, sections 92 and 97) the
!! eigenfunctions are \f$ e^{-z/2L} \cos(m \pi z / L_z) \f$ (odd \f$ m \f$), with growth rate
!! \f[ \gamma^2 = k^2 \left[ \frac{g / L}{k^2 + m^2 \pi^2 / L_z^2 + 1 / 4L^2} - v_A^2 \right] \f]
!! The field along \f$ x \f$ stabilizes all modes once \f$ B_0 \f$ exceeds
!! \f[ B_c = \sqrt{\mu_0 m_i n_0} \, \sqrt{\frac{g / L}{k^2 + \pi^2 / L_z^2 + 1/4L^2}} \f]
!! while a field along \f$ y \f$ (\f$ \mathbf{k} \cdot \mathbf{B} = 0 \f$) does not change \f$ \gamma \f$.
!! With \f$ L_x = L_z = 1 \f$ m, \f$ L = 0.5 \f$ m, \f$ n_0 = 10^{19} \f$ m\f$^{-3}\f$, \f$ T_0 = 100 \f$ eV and \f$ g = 10^8 \f$ m/s\f$^2\f$, the
!! hydrodynamic growth rate is \f$ \gamma = 1.252 \times 10^4 \f$ s\f$^{-1}\f$, \f$ B_c = 2.89 \times 10^{-4} \f$ T, and
!! \f$ g L / c_s^2 = 1.6 \times 10^{-3} \f$ (\f$ c_s^2 = \Gamma p_0 / m_i n_0 \f$), so compressibility changes \f$ \gamma \f$ by much
!! less than 1%.
!!
!!\subsection doc_mug2d_rt_ex_tanh Transition layer (nonlinear runs)
!!
!! \f[ n_0(z) = n_0 \left[ 1 + A \tanh(z/\delta) \right], \qquad B_x = B_0, \qquad \psi_0(z) = -B_0 z \f]
!! with Atwood number \f$ A = (n_{top} - n_{bottom})/(n_{top} + n_{bottom}) \f$ and layer half-width \f$ \delta \f$. The uniform field
!! exerts no force, so
!! \f[ p(z) = p_0 - m_i n_0 g \left[ z + A \delta \ln \cosh(z/\delta) \right] \f]
!! The seed is localized on the layer, \f$ f(z) = \cos(\pi z / L_z) e^{-(z/\sigma)^2} \f$. Far from the layer both plasmas
!! are uniform, so the walls see (almost) no gradients. There is no closed form for \f$ \gamma \f$; the driver prints the
!! sharp-interface (\f$ \delta \to 0 \f$), incompressible estimate between walls,
!! \f[ \gamma^2 = A g k \tanh(k L_z / 2) - k^2 v_A^2 \f]
!! which is an upper bound: a layer of finite thickness grows more slowly. The defaults are \f$ L_x = 1 \f$ m, \f$ L_z = 2 \f$ m,
!! \f$ A = 0.5 \f$ (density ratio 3), \f$ \delta = 0.1 \f$ m, \f$ \sigma = 0.2 \f$ m, with \f$ n_0 \f$, \f$ T_0 \f$ and \f$ g \f$ as above
!! (\f$ T \f$ from 200 eV at the bottom to 66 eV at the top).
!!
!!\subsection doc_mug2d_rt_ex_bcs Boundary conditions
!!
!! The walls are free slip: \f$ v_z \f$ and \f$ v_y \f$ are zero (Dirichlet), and \f$ v_x \f$ has the natural
!! condition. The flux \f$ \psi \f$ and \f$ B_y \f$ keep their default (Dirichlet) wall conditions, so the walls
!! are perfect conductors. The density and temperature use natural (zero-flux) conditions. This matters:
!! with \f$ v_z = 0 \f$ on a wall, \f$ \partial_t n = -n \partial_z v_z \f$ is not zero there, and fixing
!! \f$ n \f$ and \f$ T \f$ on the walls of a stratified plasma drives a spurious mode growing at a rate
!! \f$ \sim c_s / L \f$, much faster than the RT instability (see \ref doc_mug_xmhd_2d_eqns_bc). Hydrostatic balance
!! of the discrete equilibrium also relies on the interpolated (group FE) pressure described in
!! \ref doc_mug_xmhd_2d_eqns_pres; the full 2D system, including gravity, is given in \ref doc_mug_xmhd_2d_eqns.
!!
!!\section doc_mug2d_rt_ex_code Code walk through
!!
!!\subsection doc_mug2d_rt_ex_code_helper Helper module
!!
!! The helper module holds the case parameters and the functions that evaluate the equilibrium profiles
!! and the seed. Each takes a position `pt`, with `pt(2)` the vertical coordinate \f$ z \f$, and is passed
!! to a \ref oft_scalar_inits::poss_scalar_bfield "poss_scalar_bfield" object for projection onto the
!! Lagrange basis.
! START SOURCE
MODULE rt_2d_helpers
USE oft_local
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
IMPLICIT NONE
CHARACTER(LEN=4) :: profile = 'exp' !< Equilibrium profile: 'exp' (linear benchmark) or 'tanh' (transition layer)
REAL(r8) :: n_0 = 1.d19 !< Density at z=0 [m^-3]
REAL(r8) :: T_0 = 1.d2 !< Temperature at z=0 [eV]
REAL(r8) :: L_n = 0.5d0 !< Density scale length [m] (> 0 is unstable, 'exp' only)
REAL(r8) :: atwood = 0.5d0 !< Atwood number (n_top-n_bot)/(n_top+n_bot) ('tanh' only, > 0 is unstable)
REAL(r8) :: delta = 0.04d0 !< Half-width of the transition layer [m] ('tanh' only)
REAL(r8) :: sigma = 0.15d0 !< Width of the seed about z=0 [m] ('tanh' only)
REAL(r8) :: g = 1.d8 !< Gravitational acceleration along -z [m/s^2]
REAL(r8) :: B_0 = 0.d0 !< In-plane field along x at z=0 [T]
REAL(r8) :: V1 = 1.d0 !< Seed amplitude [m/s]
REAL(r8) :: den_scale = 1.d19 !< Density normalization used by the solver [m^-3]
REAL(r8) :: L_x = 1.d0 !< Periodic length in x [m] (set from the mesh)
REAL(r8) :: L_z = 1.d0 !< Distance between the walls [m] (set from the mesh)
CONTAINS
!> Equilibrium density n0(z)/n_0
FUNCTION n_prof(z) RESULT(f)
REAL(r8), INTENT(in) :: z
REAL(r8) :: f
IF(profile=='tanh')THEN
  f = 1.d0 + atwood*TANH(z/delta)
ELSE
  f = EXP(z/L_n)
END IF
END FUNCTION n_prof
!> Equilibrium density n0(z), normalized by den_scale
SUBROUTINE n_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = (n_0/den_scale)*n_prof(pt(2))
END SUBROUTINE n_eq
!> Equilibrium pressure [Pa] from total force balance, d(p + Bx^2/(2 mu0))/dz = -m_i n0 g
FUNCTION p_eq(z) RESULT(p)
REAL(r8), INTENT(in) :: z
REAL(r8) :: p,x
IF(profile=='tanh')THEN
  ! int_0^z n0 dz' = n_0 (z + A delta ln cosh(z/delta)); the uniform field exerts no force
  x = ABS(z/delta)
  p = 2.d0*elec_charge*n_0*T_0 - proton_mass*n_0*g*(z + atwood*delta*(x + LOG(1.d0+EXP(-2.d0*x)) - LOG(2.d0)))
ELSE
  p = 2.d0*elec_charge*n_0*T_0 - (proton_mass*n_0*g*L_n + B_0**2/(2.d0*mu0))*(EXP(z/L_n)-1.d0)
END IF
END FUNCTION p_eq
!> Equilibrium temperature [eV], with p = 2 n e T
SUBROUTINE T_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = p_eq(pt(2))/(2.d0*elec_charge*n_0*n_prof(pt(2)))
END SUBROUTINE T_eq
!> Equilibrium flux, Bx = -dpsi/dz: B_0*exp(z/(2 L)) ('exp', constant Alfven speed) or B_0 ('tanh')
SUBROUTINE psi_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
IF(profile=='tanh')THEN
  val = -B_0*pt(2)
ELSE
  val = -2.d0*L_n*B_0*EXP(pt(2)/(2.d0*L_n))
END IF
END SUBROUTINE psi_eq
!> Vertical structure f(z) of the seed and its derivative: the m=1 eigenmode e^{-z/2L} cos(pi z/L_z) ('exp'),
!! or cos(pi z/L_z) e^{-(z/sigma)^2}, localized on the layer ('tanh')
SUBROUTINE seed_f(z, f, fp)
REAL(r8), INTENT(in) :: z
REAL(r8), INTENT(out) :: f, fp
REAL(r8) :: e
IF(profile=='tanh')THEN
  e = EXP(-(z/sigma)**2)
  f = COS(pi*z/L_z)*e
  fp = -(pi/L_z)*SIN(pi*z/L_z)*e - 2.d0*z/sigma**2*f
ELSE
  e = EXP(-z/(2.d0*L_n))
  f = e*COS(pi*z/L_z)
  fp = -e*(COS(pi*z/L_z)/(2.d0*L_n) + (pi/L_z)*SIN(pi*z/L_z))
END IF
END SUBROUTINE seed_f
!> Seed vz = V1 f(z) cos(kx)
SUBROUTINE vz_seed(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
REAL(r8) :: f,fp
CALL seed_f(pt(2),f,fp)
val = V1*f*COS(2.d0*pi*pt(1)/L_x)
END SUBROUTINE vz_seed
!> Seed vx = -(V1/k) f'(z) sin(kx), from div(v) = 0
SUBROUTINE vx_seed(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
REAL(r8) :: f,fp,k
k = 2.d0*pi/L_x
CALL seed_f(pt(2),f,fp)
val = -(V1/k)*fp*SIN(k*pt(1))
END SUBROUTINE vx_seed
END MODULE rt_2d_helpers
!!\subsection doc_mug2d_rt_ex_code_driver Driver program
!!
!! The driver reads the case parameters, builds the FE space, projects the equilibrium and the seed,
!! overrides the default boundary conditions, and runs a linear or nonlinear simulation (or, with
!! `plot_run=T`, converts the restart files to XDMF/HDF5 plot files).
PROGRAM MUG_rt_2d
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_blag_operators, ONLY: oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
!---Physics
USE xmhd_2d
USE rt_2d_helpers
IMPLICIT NONE
!!\subsubsection doc_mug2d_rt_ex_code_vars Local variables and input options
!!
!! The physical parameters (`profile`, `n_0`, `T_0`, `L_n`, `atwood`, `delta`, `sigma`, `g`, `B_0`, `V1`,
!! `den_scale`) live in the helper module and are read from the `rt_2d_options` group together with the run
!! options below. `L_n` is used only by the exponential profile, and `atwood`, `delta` and `sigma` only by the
!! transition layer. `B_0` is the in-plane field (at \f$ z = 0 \f$ for the exponential profile) [T], `by0` the
!! uniform out-of-plane field [T], and `V1` the seed amplitude [m/s]. The domain size is taken from the mesh.
INTEGER(i4) :: io_unit,ierr,i
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
REAL(r8) :: zmin,zmax,k,K2,v_A,v_A_crit,B_crit,gamma_rt,c_s2,p_min
INTEGER(i4) :: order = 2 !< FE order
INTEGER(i4) :: nsteps = 100 !< Number of time steps
INTEGER(i4) :: rst_freq = 10 !< Restart (and plot) file frequency
INTEGER(i4) :: ittarget = 1000 !< Target iteration count for time step adaptation (large keeps dt fixed)
INTEGER(i4) :: pre_freq = 1 !< Preconditioner update frequency in time steps (`use_mfnk=T` only)
LOGICAL :: linear = .FALSE. !< Run the linearized model?
LOGICAL :: plot_run = .FALSE. !< Generate plot files from existing restart files instead of running
LOGICAL :: pm = .FALSE. !< Print solver progress?
LOGICAL :: use_mfnk = .FALSE. !< Use the matrix-free Newton-Krylov method?
LOGICAL :: timestep_cn = .TRUE. !< Use second-order time stepping (CN nonlinear, BDF2 linear)
LOGICAL :: upwind = .FALSE. !< Use streamline upwinding for n, V and T (nonlinear runs only)
REAL(r8) :: by0 = 0.d0 !< Uniform out-of-plane field [T]
REAL(r8) :: nu = 1.d0 !< Kinematic viscosity [m^2/s]
REAL(r8) :: eta = 1.d0 !< Resistivity [m^2/s] (eta/mu0)
REAL(r8) :: chi = 0.d0 !< Thermal diffusivity [m^2/s]
REAL(r8) :: D_diff = 0.d0 !< Particle diffusivity [m^2/s]
REAL(r8) :: gamma = 5.d0/3.d0 !< Adiabatic index
REAL(r8) :: dt = 2.d-6 !< Time step [s]
REAL(r8) :: lin_tol = 1.d-10 !< Absolute linear solver tolerance
REAL(r8) :: nl_tol = 1.d-8 !< Absolute nonlinear solver tolerance
NAMELIST/rt_2d_options/order,linear,plot_run,pm,profile,n_0,T_0,L_n,atwood,delta,sigma,g,B_0,by0,V1,nu,eta,chi, &
  D_diff,gamma,dt,nsteps,rst_freq,use_mfnk,pre_freq,timestep_cn,upwind,ittarget,lin_tol,nl_tol,den_scale
!!\subsubsection doc_mug2d_rt_ex_code_setup Mesh and FE setup
!!
!! The mesh must be periodic in \f$ x \f$ and centered on \f$ z = 0 \f$, since the profiles are defined
!! relative to the midplane.
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,rt_2d_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading "rt_2d_options"','MUG_rt_2d',__FILE__)
CALL multigrid_construct_surf(mg_mesh)
L_x=oft_mpi_max(MAXVAL(mg_mesh%smesh%r(1,:)))+oft_mpi_max(-MINVAL(mg_mesh%smesh%r(1,:)))
zmin=-oft_mpi_max(-MINVAL(mg_mesh%smesh%r(2,:)))
zmax=oft_mpi_max(MAXVAL(mg_mesh%smesh%r(2,:)))
L_z=zmax-zmin
IF(ABS(zmax+zmin)>1.d-6*L_z)CALL oft_abort('Mesh must be centered in z','MUG_rt_2d',__FILE__)
CALL mhd_sim%setup(mg_mesh,order)
!!\subsubsection doc_mug2d_rt_ex_code_theory Theoretical growth rate
!!
!! For reference the driver prints the incompressible growth rate (exact for the exponential profile, the
!! sharp-interface estimate for the transition layer), the corresponding critical field, the compressibility
!! parameters, and the density and temperature on the walls. Nonlinear runs with the exponential profile are
!! stopped here (\ref doc_mug2d_rt_ex_exp_nl).
IF(profile/='exp'.AND.profile/='tanh')CALL oft_abort('Unknown profile','MUG_rt_2d',__FILE__)
IF(profile=='exp'.AND.(.NOT.linear).AND.(.NOT.plot_run))CALL oft_abort( &
  'Nonlinear runs require profile="tanh" (see the example docs)','MUG_rt_2d',__FILE__)
k=2.d0*pi/L_x
v_A=B_0/SQRT(mu0*proton_mass*n_0)
IF(profile=='tanh')THEN
  ! Sharp-interface estimate between walls: gamma^2 = A g k tanh(k L_z/2) - k^2 v_A^2
  gamma_rt=atwood*g*k*TANH(k*L_z/2.d0)
  v_A_crit=SQRT(ABS(gamma_rt))/k
  gamma_rt=gamma_rt-(k*v_A)**2
ELSE
  K2=k**2+(pi/L_z)**2+1.d0/(4.d0*L_n**2)
  v_A_crit=SQRT(ABS(g/L_n)/K2)
  gamma_rt=k**2*((g/L_n)/K2-v_A**2)
END IF
B_crit=v_A_crit*SQRT(mu0*proton_mass*n_0)
gamma_rt=SIGN(SQRT(ABS(gamma_rt)),gamma_rt)
c_s2=gamma*2.d0*elec_charge*T_0/proton_mass
p_min=MIN(p_eq(zmin),p_eq(zmax))
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')'Rayleigh-Taylor case, profile = ',TRIM(profile)
  WRITE(*,'(A,2ES12.4)')'  Lx, Lz [m]          = ',L_x,L_z
  WRITE(*,'(A,ES12.4)') '  B_0 [T]             = ',B_0
  WRITE(*,'(A,ES12.4)') '  B_c [T]             = ',B_crit
  WRITE(*,'(A,ES12.4)') '  B_0/B_c             = ',B_0/B_crit
  WRITE(*,'(A,ES12.4)') '  by0 [T]             = ',by0
  WRITE(*,'(A,ES12.4)') '  v_A [m/s]           = ',v_A
  WRITE(*,'(A,ES12.4)') '  v_A^2/c_s^2         = ',v_A**2/c_s2
  WRITE(*,'(A,2ES12.4)')'  n at walls [m^-3]   = ',n_0*n_prof(zmin),n_0*n_prof(zmax)
  WRITE(*,'(A,2ES12.4)')'  T at walls [eV]     = ',p_eq(zmin)/(2.d0*elec_charge*n_0*n_prof(zmin)), &
    p_eq(zmax)/(2.d0*elec_charge*n_0*n_prof(zmax))
  IF(profile=='tanh')THEN
    WRITE(*,'(A,ES12.4)') '  g delta/c_s^2       = ',g*delta/c_s2
    IF(gamma_rt>=0.d0)THEN
      WRITE(*,'(A,ES12.4)')'  gamma (sharp interface) [1/s] = ',gamma_rt
    ELSE
      WRITE(*,'(A,ES12.4)')'  Stable (sharp interface), omega [1/s] = ',-gamma_rt
    END IF
  ELSE
    WRITE(*,'(A,ES12.4)') '  g L/c_s^2           = ',g*L_n/c_s2
    IF(gamma_rt>=0.d0)THEN
      WRITE(*,'(A,ES12.4)')'  gamma (theory) [1/s] = ',gamma_rt
    ELSE
      WRITE(*,'(A,ES12.4)')'  Stable, omega (theory) [1/s] = ',-gamma_rt
    END IF
  END IF
  WRITE(*,*)
END IF
IF(p_min<=0.d0)CALL oft_abort('Equilibrium pressure is not positive (increase T_0)','MUG_rt_2d',__FILE__)
!!\subsubsection doc_mug2d_rt_ex_code_bcs Boundary conditions
!!
!! \ref xmhd_2d::oft_xmhd_2d_sim::setup "setup" makes all fields Dirichlet on the boundary by default.
!! Here \f$ v_x \f$, \f$ n \f$, and \f$ T \f$ are replaced with natural conditions (free-slip,
!! zero-flux walls); the walls are the only boundary, since \f$ x \f$ is periodic.
NULLIFY(mhd_sim%velx_bc,mhd_sim%n_bc,mhd_sim%T_bc)
ALLOCATE(mhd_sim%velx_bc(oft_blagrange%ne))
mhd_sim%velx_bc=.FALSE.
mhd_sim%n_bc=>mhd_sim%velx_bc
mhd_sim%T_bc=>mhd_sim%velx_bc
!!\subsubsection doc_mug2d_rt_ex_code_run_opts Simulation options
!!
!! In nonlinear runs the time step is reduced (never increased beyond `dt`) to keep the number of linear
!! iterations per step near `ittarget`; a large value keeps it fixed. Linear runs always use a fixed time step.
!! Time centering (`timestep_cn=T`, BDF2 in linear runs and Crank-Nicolson in nonlinear
!! runs) avoids the \f$ \gamma \Delta t / 2 \f$ bias of backward Euler in the measured growth rate.
!! `upwind=T` adds streamline upwinding (\ref doc_mug_xmhd_2d_eqns_upwind) to the advection of \f$ n \f$, \f$ V \f$
!! and \f$ T \f$ in nonlinear runs; it is ignored in linear runs.
mhd_sim%g_accel=g
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%D_diff=D_diff
mhd_sim%gamma=gamma
mhd_sim%den_scale=den_scale
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
mhd_sim%pre_freq=pre_freq
mhd_sim%linear=linear
mhd_sim%lin_tol=lin_tol
mhd_sim%nl_tol=nl_tol
mhd_sim%ittarget=ittarget
mhd_sim%timestep_cn=timestep_cn
mhd_sim%upwind=upwind
oft_env%pm=pm
!!\subsubsection doc_mug2d_rt_ex_code_plot Plotting
!!
!! With `plot_run=T` the restart files from a previous run (in the same directory and with the same
!! input file) are converted to XDMF/HDF5 plot files by \ref xmhd_2d::xmhd_2d_plot "xmhd_2d_plot".
!! In linear runs the equilibrium is added back, so the plotted fields are \f$ u_0 + \delta u \f$.
IF(plot_run)THEN
  CALL xmhd_2d_plot(mhd_sim)
  CALL oft_finalize
END IF
!!\subsubsection doc_mug2d_rt_ex_code_ic Initial conditions
!!
!! Each field is projected onto the Lagrange basis with a mass-matrix solve. The equilibrium (\f$ n \f$, \f$ T \f$,
!! \f$ \psi \f$, \f$ B_y \f$) goes into `u0`, and the seed (\f$ v_x \f$, \f$ v_z \f$) into `u`; \f$ v_z \f$ is set to
!! exactly zero on the walls. In nonlinear runs the equilibrium is then added to `u`, while linear runs
!! evolve the perturbation `u` about `u0`.
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
!---Equilibrium
field_init%func=>n_eq
CALL project_field(1,mhd_sim%u0)
field_init%func=>T_eq
CALL project_field(5,mhd_sim%u0)
field_init%func=>psi_eq
CALL project_field(6,mhd_sim%u0)
CALL u%set(by0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u0%restore_local(vec_vals,7)
!---Seed
field_init%func=>vx_seed
CALL project_field(2,mhd_sim%u)
field_init%func=>vz_seed
CALL project_field(4,mhd_sim%u,zero_bc=.TRUE.)
IF(.NOT.linear)CALL mhd_sim%u%add(1.d0,1.d0,mhd_sim%u0)
!---Cleanup
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
DEALLOCATE(minv%pre)
CALL minv%delete
DEALLOCATE(minv)
!!\subsubsection doc_mug2d_rt_ex_code_run Run simulation
!!
!! The kinetic energy (`ven`) in the history file `xmhd_2d.hist` gives the growth rate,
!! \f$ \gamma = \frac{1}{2} \partial_t \ln E_k \f$ (see `run_scan.py`). In linear runs it is computed with the
!! equilibrium density.
IF(linear)THEN
  CALL mhd_sim%run_lin_simulation()
ELSE
  CALL mhd_sim%run_simulation()
END IF
CALL oft_finalize
CONTAINS
!> Project field_init onto the Lagrange basis and store it in field `ifield` of `uout`
SUBROUTINE project_field(ifield,uout,zero_bc)
INTEGER(i4), INTENT(in) :: ifield !< Index of the field in `uout`
CLASS(oft_vector), INTENT(inout) :: uout !< Solution vector
LOGICAL, OPTIONAL, INTENT(in) :: zero_bc !< Zero the field on boundary nodes?
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
END PROGRAM MUG_rt_2d
!!\section doc_mug2d_rt_ex_input Input files
!!
!! The default input file runs the nonlinear transition-layer case without a field on a \f$ 48 \times 96 \f$ mesh at
!! second order (\ref doc_mug2d_rt_ex_nl), with an adaptive time step, streamline upwinding and \f$ D = \chi = 0.3 \f$ m\f$^2\f$/s
!! (\ref doc_mug2d_rt_ex_upwind), up to the time the bubble and spikes reach the walls (about an hour on two cores). `run_nonlinear.py` uses the same file for the linear/nonlinear comparison, and
!! `run_scan.py` sets up the exponential profile for the linear scan (\ref doc_mug2d_rt_ex_scan).
!!
! OFT_DOC_INCLUDE: oft.in
!!
!! The preconditioner is set in `oft_in.xml`: one GMRES iteration preconditioned by an LU factorization
!! of the local Jacobian (block Jacobi, one block per process). The default diagonal preconditioner
!! does not converge for this problem.
!!
! OFT_DOC_INCLUDE: oft_in.xml
!!
!! To plot a run, set `plot_run=T` in `oft.in` and run the executable again in the same directory. This writes
!! XDMF/HDF5 files (in `all/`) for VisIt or ParaView, which `plot_rt_2d.py` also reads to draw density snapshots
!! with matplotlib, eg.
!! \verbatim
!! python plot_rt_2d.py --runs . --times 6.E-4 7.E-4 7.4E-4 7.6E-4
!! \endverbatim
!!
!!\section doc_mug2d_rt_ex_results Results
!!
!!\subsection doc_mug2d_rt_ex_scan Linear growth rate scan (exponential profile)
!!
!! `run_scan.py` runs the linearized model for \f$ B_0 / B_c \in \{0, 0.25, 0.5, 0.75, 0.9, 1.1\} \f$, plus one case
!! with \f$ B_0 = 0 \f$ and a strong out-of-plane field \f$ b_{y0} = 10^{-2} \f$ T (\f$ 35 B_c \f$, \f$ \beta \approx 8 \f$).
!! It uses an \f$ 8 \times 8 \f$ mesh at third order, \f$ \Delta t = 2 \times 10^{-6} \f$ s (\f$ \gamma \Delta t \le 0.025 \f$),
!! \f$ \nu = 1 \f$ m\f$^2\f$/s, \f$ \eta = 0.1 \f$ m\f$^2\f$/s, \f$ D = \chi = 0 \f$, and \f$ V_1 = 1 \f$ m/s, and runs each case
!! to \f$ \gamma t = 6 \f$ (about 30 s each on one core). The growth rate is fit to the whole kinetic energy history
!! with \f$ E_k = A \cosh^2(\gamma t) \f$, and checked with the slope of \f$ \frac{1}{2} \ln E_k \f$ for \f$ \gamma t > 3 \f$:
!!
!! | \f$ B_0 / B_c \f$ | \f$ b_{y0} \f$ [T] | \f$ \gamma \f$ theory [s\f$^{-1}\f$] | \f$ \gamma \f$ MUG, \f$ \cosh^2 \f$ fit | Ratio | \f$ \gamma \f$ MUG, \f$ \gamma t > 3 \f$ | Ratio |
!! |------|------|------------|------------|--------|------------|--------|
!! | 0    | 0    | 1.2523e4   | 1.2508e4   | 0.9988 | 1.2502e4   | 0.9984 |
!! | 0.25 | 0    | 1.2125e4   | 1.2111e4   | 0.9988 | 1.2105e4   | 0.9984 |
!! | 0.5  | 0    | 1.0845e4   | 1.0832e4   | 0.9988 | 1.0827e4   | 0.9984 |
!! | 0.75 | 0    | 8.2831e3   | 8.2751e3   | 0.9990 | 8.2722e3   | 0.9987 |
!! | 0.9  | 0    | 5.4586e3   | 5.4649e3   | 1.0012 | 5.4645e3   | 1.0011 |
!! | 1.1  | 0    | stable     | stable     |        |            |        |
!! | 0    | 1e-2 | 1.2523e4   | 1.2507e4   | 0.9987 | 1.2501e4   | 0.9983 |
!!
!! \image html MUG_rt2d_ex-scan.png "Linear growth rate against in-plane field strength, compared with the incompressible, ideal result."
!!
!! The remaining ~0.1% comes mostly from viscosity and compressibility. The results do not change with a
!! \f$ 16 \times 16 \f$ mesh or a four times smaller time step. Above \f$ B_c \f$ the kinetic energy oscillates
!! (gravity-Alfven waves) without growing, and the out-of-plane field leaves \f$ \gamma \f$ unchanged, as expected
!! for \f$ \mathbf{k} \cdot \mathbf{B} = 0 \f$.
!!
!! \note Resistivity matters near marginal stability: the field slips through the plasma at a rate
!! \f$ \sim \eta k^2 \f$, which weakens field-line bending, and because \f$ \gamma^2 \f$ is a small difference near
!! \f$ B_c \f$ the relative effect on \f$ \gamma \f$ is amplified by \f$ 1 / (1 - B_0^2/B_c^2) \f$. With \f$ \eta = 1 \f$ m\f$^2\f$/s
!! the measured \f$ \gamma \f$ at \f$ 0.9 B_c \f$ is 1.6% above the ideal value (0.1% at \f$ \eta = 0.1 \f$), which is why the scan uses
!! the smaller resistivity.
!!
!!\subsection doc_mug2d_rt_ex_nl Transition layer: linear and nonlinear growth
!!
!! `run_nonlinear.py` runs the transition-layer case from the default input file, linearized and nonlinear, for
!! \f$ B_0 = 0 \f$ and \f$ B_0 = 2 \times 10^{-4} \f$ T (\f$ 0.49 B_c \f$ for a sharp interface), with \f$ V_1 = 1 \f$ m/s,
!! \f$ \nu = 10 \f$ m\f$^2\f$/s, \f$ \eta = 1 \f$ m\f$^2\f$/s, and \f$ D = \chi = 0 \f$ (upwinding, which is on in the input file, is
!! used in the nonlinear runs only; it changes the nonlinear growth rate by less than 0.03%). The viscosity damps the short-wavelength
!! modes of the layer, which grow faster than the seeded mode (up to \f$ \sim \sqrt{A g / \delta} \f$) and are poorly resolved.
!! Since the nonlinear runs start from the projected equilibrium, which is not in exact discrete balance, the layer
!! adjusts slightly at the start of the run. This x-independent motion adds a floor of about \f$ 4 \times 10^{-10} \f$ J/m to
!! the kinetic energy, not far below the energy of the seed (\f$ 1.7 \times 10^{-9} \f$ J/m), so the comparison uses the
!! amplitude of the \f$ \cos(kx) \f$ component of \f$ v_z \f$ instead, which the adjustment does not contribute to.
!! Both runs are fit over the same time window, in which the nonlinear amplitude is between 3 and 30 m/s:
!!
!! | \f$ B_0 \f$ [T] | Window [\f$ \mu \f$s] | \f$ \gamma \f$ linear [s\f$^{-1}\f$] | \f$ \gamma \f$ nonlinear [s\f$^{-1}\f$] | Ratio | \f$ \gamma \f$ sharp interface [s\f$^{-1}\f$] |
!! |--------|---------|----------|----------|--------|----------|
!! | 0      | 240-400 | 1.4010e4 | 1.4005e4 | 0.9996 | 1.7724e4 |
!! | 2e-4   | 310-520 | 1.0633e4 | 1.0629e4 | 0.9996 | 1.5461e4 |
!!
!! \image html MUG_rt2d_ex-modes.png "Amplitude of the cos(kx) component of v_z in the linearized (dashed) and nonlinear (solid) runs. The fit window is shaded."
!!
!! The two models agree to better than 0.05% until the amplitude reaches a few hundred m/s, where the nonlinear growth
!! starts to slow down. The linear growth rate changes by 0.09% between \f$ 32 \times 64 \f$ and \f$ 48 \times 96 \f$ meshes.
!! As expected, the layer grows more slowly than a sharp interface (\f$ k \delta = 0.63 \f$), and the field reduces the
!! growth rate more than for a sharp interface, since the finite layer weakens the buoyancy drive but not the
!! field-line bending.
!!
!!\subsection doc_mug2d_rt_ex_nl_evol Nonlinear evolution
!!
!! The nonlinear runs shown here use the default input file (\f$ \nu = 10 \f$, \f$ \eta = 1 \f$ m\f$^2\f$/s, \f$ D = \chi = 0.3 \f$
!! m\f$^2\f$/s, streamline upwinding, adaptive time step), run for longer. A single bubble of light plasma rises at \f$ x = 0 \f$ and the
!! heavy plasma falls in spikes at the periodic boundaries \f$ x = \pm L_x / 2 \f$. The shear along the bubble flanks rolls up into
!! Kelvin-Helmholtz vortices, giving the familiar mushroom shape, before the bubble and spikes reach the walls (\f$ t \approx 870 \, \mu \f$s
!! for \f$ B_0 = 0 \f$) and spread along them; the heavy plasma spreading along the lower wall from both sides then collides
!! below the bubble, at \f$ x = 0 \f$, and is pushed up into a thin jet. With \f$ B_0 = 2 \times 10^{-4} \f$ T the bubble grows more slowly
!! and the field, stretched along its flanks, suppresses the Kelvin-Helmholtz roll-up: the bubble keeps smooth flanks and a rounded head as
!! it approaches the upper wall (\f$ t \approx 1.3 \f$ ms), and its stem later pinches off.
!!
!! \image html MUG_rt2d_ex-density.png "Density for B_0 = 0 (top) and B_0 = 2e-4 T (bottom), with the default input file (upwinding, D = chi = 0.3 m^2/s, adaptive time step)."
!!
!! With these settings \f$ n \f$ stays within 5% of its initial range up to \f$ t = 930 \, \mu \f$s for \f$ B_0 = 0 \f$, and within
!! 6% below and 8% above it during the collision at the lower wall (to \f$ t = 1 \f$ ms), with \f$ T \f$ within 65-222 eV. With the field
!! \f$ n \f$ stays within 1% of its initial range to 1.3 ms. Without upwinding, particle diffusion and thermal conduction of about
!! 1 m\f$^2\f$/s are needed for the same result (\ref doc_mug2d_rt_ex_upwind).
!!
!! The time step must also adapt once the flow is nonlinear: without upwinding and with a fixed \f$ \Delta t = 2 \times 10^{-6} \f$ s,
!! the GMRES iteration count grows quickly after \f$ 650 \, \mu \f$s and the Newton solve fails at \f$ 710 \, \mu \f$s, whether the
!! preconditioner is rebuilt every 4 steps or every step. With `ittarget=60` the time step falls to \f$ 1-2 \times 10^{-6} \f$ s
!! around the wall contact with upwinding (\f$ 2-7 \times 10^{-7} \f$ s without it, from the roll-up on), and stays at
!! \f$ 2 \times 10^{-6} \f$ s in the field case.
!!
!!\subsection doc_mug2d_rt_ex_upwind Streamline upwinding
!!
!! The Galerkin discretization of advection has no damping at the grid scale, so without particle diffusion and thermal conduction
!! (\f$ D = \chi = 0 \f$) the density and temperature stay smooth and within their initial range only until the Kelvin-Helmholtz roll-ups
!! thin the layer to the grid scale (\f$ t \approx 730 \, \mu \f$s for \f$ B_0 = 0 \f$). Grid-scale oscillations then grow and fill the
!! domain: at \f$ 790 \, \mu \f$s \f$ n \f$ is already 9% below its initial minimum and 29% above its initial maximum, and by
!! \f$ 840 \, \mu \f$s the minimum is half its initial value. `upwind=T` adds streamline upwinding (an inconsistent SUPG method,
!! \ref doc_mug_xmhd_2d_eqns_upwind) to the advection of \f$ n \f$, \f$ V \f$ and \f$ T \f$. The table compares \f$ B_0 = 0 \f$ runs
!! of the default case (adaptive time step) up to \f$ t = 1 \f$ ms. The ranges are those of \f$ n \f$ and \f$ T \f$ over the whole run,
!! relative to the initial range for \f$ n \f$ (\f$ 5-15 \times 10^{18} \f$ m\f$^{-3}\f$; \f$ T \f$ is initially 66-201 eV). Part of the
!! departure is physical: the plasma is compressed where the bubble and spikes reach the walls, which raises \f$ n \f$ by a few percent.
!!
!! | \f$ D = \chi \f$ [m\f$^2\f$/s] | Upwinding | \f$ n \f$ leaves its initial range by 5% at | \f$ n \f$ range | \f$ T \f$ range [eV] | Run to |
!! |-----|-----|-----------------|-----------|--------|-----------------------|
!! | 0   | no  | 747 \f$ \mu \f$s (noise from 730 \f$ \mu \f$s) | 0.91-1.29 | 52-220 | 791 \f$ \mu \f$s (stopped, noise everywhere) |
!! | 0   | yes | 878 \f$ \mu \f$s (noise from 910 \f$ \mu \f$s) | 0.38-1.94 | 33-524 | 988 \f$ \mu \f$s (stopped) |
!! | 0.1 | yes | 912 \f$ \mu \f$s | 0.61-1.75 | 40-338 | 974 \f$ \mu \f$s (stopped) |
!! | 0.3 | yes | 931 \f$ \mu \f$s | 0.95-1.08 | 65-222 | 1 ms |
!! | 1   | yes | never           | 1.00-1.04 | 66-209 | 1 ms |
!! | 1   | no  | 892 \f$ \mu \f$s | 0.85-1.08 | 63-241 | 1 ms |
!!
!! \image html MUG_rt2d_ex-upwind.png "Density for B_0 = 0 with D = chi = 0 (first row), and with upwinding and D = chi = 0 (second row) or 0.3 m^2/s (third row), compared with D = chi = 1 m^2/s without upwinding (last row). The first run was stopped once grid-scale noise had filled the spikes."
!!
!! - Upwinding alone (\f$ D = \chi = 0 \f$) removes the grid-scale noise of the roll-up and of the spikes reaching the lower wall, about
!!   \f$ 190 \, \mu \f$s more than without it. It fails where the heavy plasma spreading along the lower wall collides with itself at
!!   \f$ x = 0 \f$, and in the lower corners: streamline diffusion acts only along the flow, so it does not damp the oscillations at
!!   the stagnation point, or the ripples across the layers that run along the walls. A small \f$ D \f$ and \f$ \chi \f$ are still
!!   needed there: \f$ 0.3 \f$ m\f$^2\f$/s with upwinding is cleaner than 1 m\f$^2\f$/s without it, and 1 m\f$^2\f$/s with upwinding
!!   keeps \f$ n \f$ within 4% of its initial range with no transients. The default input file uses 0.3 m\f$^2\f$/s, which changes the
!!   growth rate less (below).
!! - Upwinding does not change the growth of the instability: the growth rate of the seeded mode in the nonlinear run changes by
!!   +0.001% (\f$ B_0 = 0 \f$) and -0.03% (\f$ B_0 = 2 \times 10^{-4} \f$ T), against -0.20% for \f$ D = \chi = 0.3 \f$ m\f$^2\f$/s and -0.64% for
!!   \f$ D = \chi = 1 \f$ m\f$^2\f$/s. It is not used in linear runs, so `run_scan.py` and the linear rates are unaffected.
!! - With the field the flow stays smooth: with \f$ B_0 = 2 \times 10^{-4} \f$ T, upwinding and \f$ D = \chi = 0 \f$ the run is clean to
!!   2 ms (\f$ n \f$ within 2.6% of its initial range, no undershoot), through the bubble reaching the upper wall and the pinch-off of its stem,
!!   with the time step fixed at \f$ 2 \times 10^{-6} \f$ s throughout. It closely follows a run with \f$ D = \chi = 1 \f$ m\f$^2\f$/s
!!   without upwinding, apart from the sharper layer.
!!
!! \image html MUG_rt2d_ex-upwind-field.png "Density for B_0 = 2e-4 T with upwinding and D = chi = 0 (top), compared with D = chi = 1 m^2/s without upwinding (bottom)."
!!
!! Upwinding also helps the nonlinear solver: the extra streamline diffusion keeps the GMRES iteration count low for longer, so the
!! adaptive time step stays at \f$ 2 \times 10^{-6} \f$ s until \f$ t \approx 790 \, \mu \f$s, against \f$ 690-700 \, \mu \f$s without it.
!! It costs about 5% more per time step (\ref doc_mug2d_rt_ex_cost), but the runs reach a given time in fewer, cheaper steps: on one
!! core each, with four runs sharing the machine, \f$ t = 1 \f$ ms took 547 steps and 5800 s with upwinding (\f$ D = \chi = 0.3 \f$),
!! against 788 steps and 11200 s without it (\f$ D = \chi = 1 \f$).
!!
!!\subsection doc_mug2d_rt_ex_exp_nl Why the exponential profile is not used nonlinearly
!!
!! The exponential profile has its largest density and pressure gradients at the walls. In nonlinear runs:
!!  - The bubble and spikes reach the walls while still growing, at \f$ t \approx 440 \, \mu \f$s on a \f$ 32 \times 32 \f$ mesh
!!    (\f$ V_1 = 100 \f$ m/s). As the heavy plasma then piles up on the lower wall, grid-scale oscillations grow in
!!    \f$ n \f$ and \f$ T \f$ and the Newton solve fails, with fixed or adaptive time steps, and with \f$ \nu = \eta = 50 \f$ m\f$^2\f$/s.
!!  - Particle diffusion or thermal conduction, even \f$ D \f$ or \f$ \chi = 1 \f$ m\f$^2\f$/s, makes the first Newton solve fail at
!!    practical time steps (\f$ \Delta t \gtrsim 4 \times 10^{-7} \f$ s). The natural (zero-flux) walls are inconsistent with
!!    the equilibrium gradient \f$ dn/dz = n/L \f$, so the first step must form a thin boundary layer at each wall. The Jacobian
!!    itself is correct (it matches a finite-difference Jacobian), and the solve converges with a uniform density or with
!!    \f$ D \f$ acting on \f$ n - n_0 \f$ only.
!!  - Streamline upwinding (with \f$ D = \chi = 0 \f$, \f$ \nu = 10 \f$, \f$ \eta = 1 \f$ m\f$^2\f$/s, adaptive time step) removes the
!!    grid-scale noise that fills the domain without it by \f$ 500 \, \mu \f$s, but not the oscillations where the heavy plasma spreading
!!    along the lower wall collides at \f$ x = 0 \f$: by \f$ 530 \, \mu \f$s \f$ n \f$ ranges over 0.64-1.73 times its initial range,
!!    against 0.29-2.17 without upwinding at \f$ 516 \, \mu \f$s.
!!  - Fixing \f$ n \f$ and \f$ T \f$ on the walls instead drives a spurious, wall-localized mode (\f$ 3.4 \times 10^4 \f$ s\f$^{-1}\f$ on an
!!    \f$ 8 \times 8 \f$, third-order mesh, against \f$ 1.25 \times 10^4 \f$ s\f$^{-1}\f$ for RT; see \ref doc_mug2d_rt_ex_bcs).
!!
!! In the transition-layer case the plasma near the walls is uniform, so none of these problems arise: the
!! wall conditions are consistent with the equilibrium, and particle diffusion and thermal conduction can be used.
!!
!!\subsection doc_mug2d_rt_ex_cost Cost
!!
!! The nonlinear runs use the matrix-free Newton-Krylov method (`use_mfnk=T`), with the LU preconditioner (UMFPACK) rebuilt
!! every `pre_freq` steps. Rebuilding it every step roughly doubles the cost of the linear phase, while letting it age for more
!! than a few steps makes the GMRES iteration count grow quickly once the flow is nonlinear, so the default is `pre_freq=4`.
!! On the default \f$ 48 \times 96 \f$ mesh with two OpenMP threads this costs about 4-5 s per step while the time step is
!! fixed, and 15-20 s per step later; the default run (450 steps, to \f$ 880 \, \mu \f$s) takes about an hour. Upwinding adds about
!! 5% to the cost of a step (3.8 against 3.6 s on two threads in the linear phase, with about 15% more GMRES iterations), but
!! keeps the time step larger in the nonlinear phase (\ref doc_mug2d_rt_ex_upwind). The linear runs
!! cost about 0.5 s per step on this mesh, and the exponential-profile scan 0.01-0.04 s per step.
!!
!! \note `pre_freq` only takes effect with `use_mfnk=T`. With the assembled Jacobian as the Newton operator
!! (`use_mfnk=F`) the Jacobian and its factorization are rebuilt every step.
!!
!! \note UMFPACK fails to factorize the Jacobian on a \f$ 64 \times 128 \f$ second-order mesh in this configuration, so finer
!! meshes need another LU package or MPI (one block per process).
