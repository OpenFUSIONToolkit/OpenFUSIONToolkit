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
!! interchange direction. The linear growth rate is compared with the exact incompressible result of
!! [S. Chandrasekhar, Hydrodynamic and Hydromagnetic Stability (1961)](https://store.doverpublications.com/products/9780486640716),
!! sections 92 and 97, over a scan in field strength, and the nonlinear evolution is followed into
!! the bubble-and-spike phase.
!!
!!\section doc_mug2d_rt_ex_setup Problem setup
!!
!! The 2D solver works in the plane of the mesh, with the first mesh coordinate as \f$ x \f$ and the second as
!! the vertical coordinate \f$ z \f$; \f$ y \f$ is the ignorable direction. Gravity acts along \f$ -\hat{z} \f$
!! (\ref xmhd_2d::oft_xmhd_2d_sim::g_accel "g_accel" \f$ = g > 0 \f$). The domain is
!! \f$ x \in [-L_x/2, L_x/2] \f$ (periodic) by \f$ z \in [-L_z/2, L_z/2] \f$, bounded by rigid, perfectly
!! conducting, free-slip walls.
!!
!! The equilibrium is exponentially stratified, with the heavy plasma on top for \f$ L > 0 \f$:
!! \f[ n_0(z) = n_0 e^{z/L} \f]
!! The optional in-plane field lies along the periodic direction and has a constant Alfven speed,
!! \f$ v_A = B_0 / \sqrt{\mu_0 m_i n_0} \f$:
!! \f[ B_x(z) = B_0 e^{z/2L}, \qquad \psi_0(z) = -2 L B_0 e^{z/2L} \f]
!! where \f$ B_x = -\partial_z \psi \f$. A uniform out-of-plane field \f$ B_y = b_{y0} \f$ can also be added.
!! The pressure \f$ p = 2 n e T \f$ (equal ion and electron temperatures, \f$ T \f$ in eV) follows from total
!! force balance, \f$ \frac{d}{dz} \left( p + B_x^2 / 2 \mu_0 \right) = -m_i n_0 g \f$:
!! \f[ p(z) = p_0 - \left( m_i n_0 g L + \frac{B_0^2}{2 \mu_0} \right) \left( e^{z/L} - 1 \right),
!! \qquad p_0 = 2 e n_0 T_0 \f]
!! and the driver stops if \f$ p \le 0 \f$ anywhere in the domain.
!!
!! The perturbation is the \f$ m = 1 \f$ incompressible eigenmode, \f$ k = 2 \pi / L_x \f$, applied to the
!! velocity only:
!! \f[ v_z = V_1 e^{-z/2L} \cos(\pi z / L_z) \cos(k x) \f]
!! \f[ v_x = \frac{V_1}{k} e^{-z/2L} \left[ \frac{\cos(\pi z / L_z)}{2 L}
!!     + \frac{\pi}{L_z} \sin(\pi z / L_z) \right] \sin(k x) \f]
!! It is divergence free and launches no sound waves. Since the initial displacement is zero, it is an
!! equal mix of the growing and decaying modes, and the velocity grows as \f$ \cosh(\gamma t) \f$.
!!
!!\subsection doc_mug2d_rt_ex_theory Theory
!!
!! For an incompressible, ideal plasma (Chandrasekhar 1961, sections 92 and 97) the vertical
!! eigenfunctions are \f$ e^{-z/2L} \cos(m \pi z / L_z) \f$ (odd \f$ m \f$), with growth rate
!! \f[ \gamma^2 = k^2 \left[ \frac{g / L}{k^2 + m^2 \pi^2 / L_z^2 + 1 / 4L^2} - v_A^2 \right] \f]
!! The field along \f$ x \f$ stabilizes all modes once \f$ B_0 \f$ exceeds
!! \f[ B_c = \sqrt{\mu_0 m_i n_0} \, \sqrt{\frac{g / L}{k^2 + \pi^2 / L_z^2 + 1/4L^2}} \f]
!! while a field along \f$ y \f$ (\f$ \mathbf{k} \cdot \mathbf{B} = 0 \f$) does not change \f$ \gamma \f$.
!! The driver prints \f$ \gamma \f$, \f$ B_c \f$, and the compressibility parameters
!! \f$ g L / c_s^2 \f$ and \f$ v_A^2 / c_s^2 \f$ (\f$ c_s^2 = \Gamma p_0 / m_i n_0 \f$ at \f$ z = 0 \f$)
!! at start-up. With the defaults (\f$ L_x = L_z = 1 \f$ m, \f$ L = 0.5 \f$ m, \f$ n_0 = 10^{19} \f$ m\f$^{-3}\f$,
!! \f$ T_0 = 100 \f$ eV, \f$ g = 10^8 \f$ m/s\f$^2\f$) the hydrodynamic growth rate is
!! \f$ \gamma = 1.252 \times 10^4 \f$ s\f$^{-1}\f$, \f$ B_c = 2.89 \times 10^{-4} \f$ T, and
!! \f$ g L / c_s^2 = 1.6 \times 10^{-3} \f$, so compressibility changes \f$ \gamma \f$ by much less than 1%.
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
REAL(r8) :: n_0 = 1.d19 !< Density at z=0 [m^-3]
REAL(r8) :: T_0 = 1.d2 !< Temperature at z=0 [eV]
REAL(r8) :: L_n = 0.5d0 !< Density scale length [m] (> 0 is unstable)
REAL(r8) :: g = 1.d8 !< Gravitational acceleration along -z [m/s^2]
REAL(r8) :: B_0 = 0.d0 !< In-plane field along x at z=0 [T]
REAL(r8) :: V1 = 1.d0 !< Seed amplitude [m/s]
REAL(r8) :: den_scale = 1.d19 !< Density normalization used by the solver [m^-3]
REAL(r8) :: L_x = 1.d0 !< Periodic length in x [m] (set from the mesh)
REAL(r8) :: L_z = 1.d0 !< Distance between the walls [m] (set from the mesh)
CONTAINS
!> Equilibrium density n0(z), normalized by den_scale
SUBROUTINE n_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = (n_0/den_scale)*EXP(pt(2)/L_n)
END SUBROUTINE n_eq
!> Equilibrium pressure [Pa] from total force balance, d(p + Bx^2/(2 mu0))/dz = -m_i n0 g
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
!> Equilibrium flux, giving Bx = -dpsi/dz = B_0*exp(z/(2 L)) (constant Alfven speed)
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
!! The physical parameters (`n_0`, `T_0`, `L_n`, `g`, `B_0`, `V1`, `den_scale`) live in the helper module and
!! are read from the `rt_2d_options` group together with the run options below. `B_0` is the in-plane field
!! at \f$ z = 0 \f$ [T], `by0` the uniform out-of-plane field [T], and `V1` the seed amplitude [m/s]. The domain
!! size is taken from the mesh.
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
REAL(r8) :: by0 = 0.d0 !< Uniform out-of-plane field [T]
REAL(r8) :: nu = 1.d0 !< Kinematic viscosity [m^2/s]
REAL(r8) :: eta = 1.d0 !< Resistivity [m^2/s] (eta/mu0)
REAL(r8) :: chi = 0.d0 !< Thermal diffusivity [m^2/s]
REAL(r8) :: D_diff = 0.d0 !< Particle diffusivity [m^2/s]
REAL(r8) :: gamma = 5.d0/3.d0 !< Adiabatic index
REAL(r8) :: dt = 2.d-6 !< Time step [s]
REAL(r8) :: lin_tol = 1.d-10 !< Absolute linear solver tolerance
REAL(r8) :: nl_tol = 1.d-8 !< Absolute nonlinear solver tolerance
NAMELIST/rt_2d_options/order,linear,plot_run,pm,n_0,T_0,L_n,g,B_0,by0,V1,nu,eta,chi, &
  D_diff,gamma,dt,nsteps,rst_freq,use_mfnk,pre_freq,timestep_cn,ittarget,lin_tol,nl_tol,den_scale
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
!! The incompressible growth rate of the \f$ m = 1 \f$ mode, the critical field, and the compressibility
!! parameters are printed for reference.
k=2.d0*pi/L_x
K2=k**2+(pi/L_z)**2+1.d0/(4.d0*L_n**2)
v_A=B_0/SQRT(mu0*proton_mass*n_0)
v_A_crit=SQRT(ABS(g/L_n)/K2)
B_crit=v_A_crit*SQRT(mu0*proton_mass*n_0)
gamma_rt=k**2*((g/L_n)/K2-v_A**2)
gamma_rt=SIGN(SQRT(ABS(gamma_rt)),gamma_rt)
c_s2=gamma*2.d0*elec_charge*T_0/proton_mass
p_min=MIN(p_eq(zmin),p_eq(zmax))
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'Rayleigh-Taylor case'
  WRITE(*,'(A,2ES12.4)')'  Lx, Lz [m]          = ',L_x,L_z
  WRITE(*,'(A,ES12.4)') '  B_0 [T]             = ',B_0
  WRITE(*,'(A,ES12.4)') '  B_c [T]             = ',B_crit
  WRITE(*,'(A,ES12.4)') '  B_0/B_c             = ',B_0/B_crit
  WRITE(*,'(A,ES12.4)') '  by0 [T]             = ',by0
  WRITE(*,'(A,ES12.4)') '  v_A [m/s]           = ',v_A
  WRITE(*,'(A,ES12.4)') '  g L/c_s^2           = ',g*L_n/c_s2
  WRITE(*,'(A,ES12.4)') '  v_A^2/c_s^2         = ',v_A**2/c_s2
  WRITE(*,'(A,2ES12.4)')'  T at walls [eV]     = ',p_eq(zmin)/(2.d0*elec_charge*n_0*EXP(zmin/L_n)), &
    p_eq(zmax)/(2.d0*elec_charge*n_0*EXP(zmax/L_n))
  IF(gamma_rt>=0.d0)THEN
    WRITE(*,'(A,ES12.4)')'  gamma (theory) [1/s] = ',gamma_rt
  ELSE
    WRITE(*,'(A,ES12.4)')'  Stable, omega (theory) [1/s] = ',-gamma_rt
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
!! A large `ittarget` keeps the time step fixed in linear runs and lets it grow only slowly in
!! nonlinear runs. Time centering (`timestep_cn=T`, BDF2 in linear runs and Crank-Nicolson in nonlinear
!! runs) avoids the \f$ \gamma \Delta t / 2 \f$ bias of backward Euler in the measured growth rate.
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
!! The default input file runs the nonlinear case without a field on a \f$ 32 \times 32 \f$ mesh at
!! second order (see \ref doc_mug2d_rt_ex_nl). Linear runs and the field scan are set up by `run_scan.py`
!! (\ref doc_mug2d_rt_ex_scan) from this file.
!!
! OFT_DOC_INCLUDE: oft.in
!!
!! The preconditioner is set in `oft_in.xml`: one GMRES iteration preconditioned by an LU factorization
!! of the local Jacobian (block Jacobi, one block per process). The default diagonal preconditioner
!! does not converge for this problem.
!!
! OFT_DOC_INCLUDE: oft_in.xml
!!
!! To plot a run, set `plot_run=T` in `oft.in` and run the executable again in the same directory.
!!
!!\section doc_mug2d_rt_ex_results Results
!!
!! RESULTS_PLACEHOLDER
