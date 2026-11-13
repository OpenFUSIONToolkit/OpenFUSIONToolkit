!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Harris Sheet Reconnection    {#doc_mug2d_harris_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref xmhd_2d "2D MHD" module in OFT (see \ref doc_mug2d_main).
!! It simulates resistive reconnection of a Harris current sheet in a slab that is periodic in
!! \f$ x \f$, following the setup of the [GEM Challenge](https://doi.org/10.1029/1999JA900449).
!! A 3D version of a similar case is shown in \ref doc_mug_recon_ex.
!!
!!\section doc_mug2d_harris_ex_code Code Walk Through
!!
!! The code consists of the main run program and subroutines for field initialization.
!!
!!\subsection doc_mug2d_harris_ex_code_imports Module imports
!!
!! As usual, we start with a few module imports for required functionality.
! START SOURCE
PROGRAM harris_sheet
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
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE xmhd_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: free_bc => NULL() !< Neumann (no Dirichlet) BC flags for density and temperature
!!\subsection doc_mug2d_harris_ex_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve. Runtime options are read from the
!! `harris_sheet_options` group in the input file.
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
INTEGER(i4) :: pre_freq = 1 !< Preconditioner update frequency (steps, with use_mfnk=T)
INTEGER(i4) :: ittarget = 40 !< Target linear iterations per step for adaptive time step
REAL(r8) :: B_0(3)=0.d0  !< Background magnetic field
REAL(r8) :: n0 = 1.d0 !< Density scale
REAL(r8) :: velx0 = 1.d0 !< Velocity scale (x-direction)
REAL(r8) :: vely0 = 1.d0 !< Velocity scale (y-direction)
REAL(r8) :: velz0 = 1.d0 !< Velocity scale (z-direction)
REAL(r8) :: t0 = 1.d0 !< Temperature scale
REAL(r8) :: psi0 = 1.d0 !< Magnetic flux scale
REAL(r8) :: by0 = 1.d0 !< Magnetic field scale (y-direction)
REAL(r8) :: chi=1.d0 !< thermal diffusivity
REAL(r8) :: eta=1.d0 !< electrical resistivity
REAL(r8) :: nu=1.d0 !< viscosity
REAL(r8) :: gamma=2.d0/3.d0 !< adiabatic index
REAL(r8) :: D_diff=1.d0 !< Diffusion coefficient
REAL(r8) :: den_scale=1.d19 !< Density scale factor for normalization of equations (used for post-processing)
REAL(r8) :: den_inf = 2.d-1 !< Background density
REAL(r8) :: lam_n = 2.d0 !< eq density width
REAL(r8) :: lam_b = 2.d0 !< eq psi perturbation width
REAL(r8) :: L_x = 2.d0 !< Box length in x
REAL(r8) :: L_z = 2.d0 !< Box length in z
REAL(r8) :: delta = 1.d-1 !< Relative size of perturbation (<<1)
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
LOGICAL :: linear=.FALSE.
NAMELIST/harris_sheet_options/order, chi, eta, nu, gamma, D_diff, &
dt, nsteps, rst_freq, pre_freq, ittarget, use_mfnk, pm, linear, n0, psi0, velx0, vely0, velz0, &
t0, by0, den_scale, delta, lam_n, lam_b, L_x, L_z
!!\subsection doc_mug2d_harris_ex_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,harris_sheet_options,IOSTAT=ierr)
CLOSE(io_unit)
!!\subsection doc_mug2d_harris_ex_mesh Mesh and simulation setup
!! The mesh is a built-in rectangular grid of triangles (see `cube_options` in the
!! \ref doc_mug2d_harris_ex_input "input file"). Calling `setup` on the
!! \ref xmhd_2d::oft_xmhd_2d_sim "oft_xmhd_2d_sim" object builds the Lagrange finite element space,
!! the composite solution vector and the default boundary conditions.
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
blag_zerob%ML_lag_rep=>ML_oft_blagrange
!!\subsection doc_mug2d_harris_ex_bc Boundary conditions
!! Each field has a logical flag array (`n_bc`, `velx_bc`, `vely_bc`, `velz_bc`, `T_bc`, `psi_bc`, `by_bc`) with one entry per
!! degree of freedom, where `.TRUE.` holds that degree of freedom fixed at its initial value (Dirichlet). By default `setup`
!! points every flag array at the boundary nodes, so all fields are held fixed on the walls. Here the defaults for the density and
!! temperature are overridden by pointing `n_bc` and `T_bc` at an all-`.FALSE.` array after `setup`, which gives natural
!! (Neumann, zero-gradient) boundary conditions for those two fields. The velocity, \f$ \psi \f$ and \f$ B_y \f$ keep the
!! default Dirichlet conditions.
!---------------------------------------------------------------------------
! Set Neumann boundary conditions on density and temperature
!---------------------------------------------------------------------------
ALLOCATE(free_bc(SIZE(mhd_sim%n_bc)))
free_bc=.FALSE.
mhd_sim%n_bc=>free_bc
mhd_sim%T_bc=>free_bc
!!\subsection doc_mug2d_harris_ex_ic Initial conditions
!! Initial conditions are set field by field. Spatially varying fields (density and flux) are
!! \f$ L_2 \f$-projected onto the Lagrange basis using a mass matrix solve, while uniform fields
!! are set directly. Each field is then copied into the matching component of the solution vector
!! with `restore_local` (1 = density, 2-4 = velocity, 5 = temperature, 6 = \f$ \psi \f$, 7 = \f$ B_y \f$).
!! In linear mode the equilibrium is placed in `u0` and only the perturbation in `u`.
!! The density is divided by `den_scale`, because it is stored normalized.

!---------------------------------------------------------------------------
! Set initial conditions from analytic functions
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


!---Project n initial condition onto scalar Lagrange basis
field_init%func=>dens_harris
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
vec_vals = vec_vals / den_scale
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,1)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,1)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,1)
END IF

!---Set v_x initial condition
CALL u%set(velx0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,2)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,2)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,2)
END IF

!---Set v_y initial condition
CALL u%set(vely0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,3)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,3)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,3)
END IF

!---Set v_z initial condition
CALL u%set(velz0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,4)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,4)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,4)
END IF

!---Set T initial condition
CALL u%set(t0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,5)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,5)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,5)
END IF

!---Project psi initial condition onto scalar Lagrange basis
field_init%func=>psi_harris_eq
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(psi0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,6)
ELSE
    CALL mhd_sim%u%restore_local(vec_vals,6)
END IF
field_init%func=>psi_harris_pert
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(psi0)
IF(linear)THEN
    CALL u%get_local(vec_vals)
    CALL mhd_sim%u%restore_local(vec_vals,6)
ELSE
    CALL mhd_sim%u%get_local(vec_vals,6)
    CALL v%restore_local(vec_vals)
    CALL u%add(1.d0, 1.d0, v)
    CALL u%get_local(vec_vals)
    CALL mhd_sim%u%restore_local(vec_vals,6)
END IF

!---Set By initial condition
CALL u%set(by0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,7)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,7)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,7)
END IF

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)
!!\subsection doc_mug2d_harris_ex_run Run simulation
!! Finally, the physical coefficients and run parameters are copied into the simulation object.
!! Then either the nonlinear (`run_simulation`) or linear (`run_lin_simulation`) time advance is called.
!! All coefficients here are uniform, so assigning a scalar fills every mesh region.
!! After the run, \ref xmhd_2d::xmhd_2d_plot "xmhd_2d_plot" converts the restart files into XDMF/HDF5 plot files.
!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%B_0=B_0
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%dt=dt
mhd_sim%den_scale=den_scale
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%pre_freq=pre_freq
mhd_sim%ittarget=ittarget
mhd_sim%mfnk=use_mfnk
mhd_sim%linear=linear
oft_env%pm=pm

IF(linear)THEN
  CALL mhd_sim%run_lin_simulation()
ELSE
  CALL mhd_sim%run_simulation()
END IF

CALL xmhd_2d_plot(mhd_sim)
!---Finalize enviroment
CALL oft_finalize
CONTAINS

!!\subsection doc_mug2d_harris_ex_ic_funcs Initial condition functions
!! To set the initial conditions we define a set of functions, which are scaled by `n0` and `psi0` in the main program.
!! Here the mesh coordinates `pt(1)` and `pt(2)` correspond to \f$ x \f$ and \f$ z \f$.
!!
!! The equilibrium is a Harris sheet, with density
!! \f[ n = n_0 \left[ \mathrm{sech}^2(z / \lambda_n) + n_{\infty} \right] \f]
!! and flux \f$ \psi_{eq} = -\psi_0 \lambda_b \ln \cosh(z / \lambda_b) \f$, which gives
!! \f[ \mathbf{B}_{eq} = \psi_0 \tanh(z / \lambda_b) \hat{x} \f]
!! The temperature is uniform (`t0`). With \f$ \lambda_n = \lambda_b \f$, the equilibrium is in pressure balance when
!! \f$ 2 n_0 k T_0 = \psi_0^2 / 2 \mu_0 \f$, which the default input satisfies.
!!
!! Reconnection is seeded with the flux perturbation
!! \f$ \delta\psi = -\psi_0 \delta \cos(2 \pi x / L_x) \cos(\pi z / L_z) \f$, which gives
!! \f[ \delta \mathbf{B} = -\psi_0 \delta \frac{\pi}{L_z} \cos(2 \pi x / L_x) \sin(\pi z / L_z) \hat{x}
!! + \psi_0 \delta \frac{2 \pi}{L_x} \sin(2 \pi x / L_x) \cos(\pi z / L_z) \hat{z} \f]
!! This places the X-point at \f$ x = 0 \f$ and the O-points at \f$ x = \pm L_x/2 \f$.

SUBROUTINE dens_harris(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = den_inf + 1/((COSH(pt(2)/lam_n))**2)
END SUBROUTINE dens_harris

SUBROUTINE psi_harris_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -lam_b*LOG(COSH(pt(2)/lam_b))
END SUBROUTINE psi_harris_eq

SUBROUTINE psi_harris_pert(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = - delta*COS(2*pi*pt(1)/L_x)*COS(pi*pt(2)/L_z)
END SUBROUTINE psi_harris_pert

END PROGRAM harris_sheet
! STOP SOURCE
!!\section doc_mug2d_harris_ex_input Input file
!!
!! Below is an input file that can be used with this example. It uses a \f$ 64 \times 128 \f$ grid of quadratic
!! (`order=2`) elements on a \f$ 25.6 \times 12.8 \f$ m domain. The grid is periodic in \f$ x \f$ (`ref_per`), and
!! `packing` concentrates the cells in \f$ z \f$ around the current sheet (about 0.05 m at \f$ z = 0 \f$). This resolves the
!! thin reconnection layer, and keeps the density advection stable with a small particle diffusivity (`D_diff`). At the
!! walls (\f$ z = \pm 6.4 \f$ m) the density and temperature have Neumann (zero-gradient) boundary conditions, while the
!! velocity, \f$ \psi \f$ and \f$ B_y \f$ keep the default Dirichlet conditions. The run takes 4000 steps of \f$ 10^{-7} \f$ s
!! (0.4 ms), and writes a restart file every 50 steps. The `xmhd_plot_options` group sets the last restart file to plot
!! (`rst_end`, default 2000), so that the plot files cover the whole run.
!!
!!\verbatim
!!&runtime_options
!!  ppn=1
!!  debug=0
!!/
!!
!!&mesh_options
!!  meshname='test'
!!  cad_type=92
!!  nlevels=1
!!  nbase=1
!!  grid_order=2
!!  fix_boundary=T
!!/
!!
!!&cube_options
!! mesh_type=1
!! ref_per=T,F,F
!! ni=64,128,0
!! rscale=25.6,12.8,0.0
!! shift=-12.8,-6.4,0.0
!! packing=1.0,0.25,0.0
!!/
!!
!!&harris_sheet_options
!!  order=2
!!  n0=5.196374E16
!!  psi0=0.002045692328575E0,
!!  velx0=0,
!!  vely0=0,
!!  velz0=0,
!!  t0=1E2,
!!  by0=0,
!!  chi=2943.4,
!!  eta=742.59999855
!!  nu=7424.22118697
!!  delta = 1.E-1
!!  gamma=1.67
!!  D_diff=1.E1
!!  den_scale = 5.196374E16
!!  dt=1.E-7
!!  lam_b = 0.5E0
!!  lam_n = 0.5E0
!!  L_x = 25.6E0
!!  L_z = 12.8E0
!!  nsteps=4000
!!  rst_freq=50
!!  use_mfnk=T
!!  pm=F
!!  linear=F
!!/
!!
!!&xmhd_plot_options
!!  rst_end=4000
!!/
!!\endverbatim
!!
!!\subsection doc_mug2d_harris_ex_input_solver Solver specification
!!
!! The run uses matrix-free Newton-Krylov (`use_mfnk=T`). The approximate Jacobian is used only as the
!! preconditioner, and here it is applied as a few GMRES iterations with a block-Jacobi LU inner preconditioner. This is set in the
!! `xmhd2d` element of the XML input file:
!!
!!```xml
!!<oft>
!!  <xmhd2d>
!!    <pre type="gmres">
!!      <its>2</its>
!!      <nrits>2</nrits>
!!      <pre type="block_jacobi">
!!        <nlocal>-1</nlocal>
!!        <solver type="lu"></solver>
!!      </pre>
!!    </pre>
!!  </xmhd2d>
!!</oft>
!!```
!!
!! The XML file is passed as the second command-line argument:
!!
!!\verbatim
!!~$ ./harris_sheet oft.in oft_in.xml
!!\endverbatim
!!
!! If the XML file is omitted, the solver falls back to diagonal preconditioning. This roughly doubles the number of
!! linear iterations per step, and the adaptive time step then shrinks to keep them near `ittarget`.
!!
!!\section doc_mug2d_harris_ex_post Post-Processing
!!
!! During the run, restart files (`xmhd2d_NNNNN.rst`) are written every `rst_freq` steps, and global quantities
!! (energies, average density and temperature, iteration counts) are written to the history file `xmhd_2d.hist`.
!!
!! After the time advance, \ref xmhd_2d::xmhd_2d_plot "xmhd_2d_plot" reads the restart files and saves the fields for plotting
!! to `oft_xdmf.0001.h5` in the run directory. This file holds the density `n`, velocity `V`, flux `psi`, magnetic field `B` and,
!! for compressible runs, the temperature `T`. For incompressible runs, the pressure `p` is saved separately to
!! `pressure/oft_xdmf.0001.h5`, because it uses a lower-order mesh.
!!
!! To view the fields, run `bin/OFT_build_xdmf.py` in the run directory (or in `pressure/`). This writes `*.xmf` files next to the
!! `.h5` file, which can be loaded by [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!!
!!\verbatim
!!~$ python OFT_build_xdmf.py
!!\endverbatim
!!
!! The vector fields `B` and `V` are written in the coordinates of the mesh, whose first two coordinates are \f$ (x, z) \f$
!! and whose third is the out-of-plane direction. So `B` is ordered \f$ (B_x, B_z, B_y) \f$ and `V` is ordered
!! \f$ (u_x, u_z, u_y) \f$, and in a visualization program the vertical axis of the plot is \f$ z \f$.
!!
!! \note The plotted `B` and `psi` include only the evolved fields \f$ \tilde{\psi} \f$ and \f$ B_y \f$. The background field
!! `B_0` and the vacuum flux `psi_vac` are not added. Both are zero in this example, but for runs that use them they must be
!! added in post-processing to get the total field.
!!
!! The history file can be read in Python with `OpenFUSIONToolkit.io.histfile`.
!!
!! The figures below show contours of \f$ \psi \f$ (the in-plane field lines) at four times, and the in-plane flow and density at
!! the final time.
!!
!! \image html MUG2D_harris_ex-Psi.png "Flux contours at four times"
!! \image html MUG2D_harris_ex-Flow.png "In-plane flow and density at the final time"
