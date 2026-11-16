!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!ThickCurr Example: Eddy currents in a conducting slab    {#doc_thickcurr_slab_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref thick_wall "ThickCurr" module in OFT to model
!! eddy currents in volumetric ("thick") conducting structures. In this example a stainless steel
!! slab (\f$ \eta = 7 \times 10^{-7} \Omega \cdot m \f$, thickness 0.2 m) is placed at the
!! center of a 1 m cube with perfectly conducting walls, with the remainder of the cube treated
!! as vacuum. We will:
!!
!!  1. Compute the slowest decay times (L/R time constants) of the system with an eigenvalue solve
!!  2. Release an initial magnetic field and follow its decay in time
!!  3. Generate plot files of the eigenmodes and time-dependent fields for ParaView
!!
!! ThickCurr evolves the vector potential (\f$ B = \nabla \times A \f$, \f$ E = -\partial A / \partial t \f$)
!! following
!! \f[ \frac{1}{\eta} \frac{\partial A}{\partial t} + \nabla \times \left( \frac{1}{\mu_0}
!! \nabla \times A \right) = 0 \f]
!! where vacuum regions are treated as a pseudo-vacuum with a large resistivity.
!!
!!\section doc_thickcurr_slab_ex_analytic Analytic solution
!!
!! For this geometry the slowest mode varies as \f$ cos(\pi x) \f$ (or \f$ cos(\pi y) \f$), giving
!! a horizontal wavenumber \f$ \kappa = \pi \f$. Matching a potential field in each vacuum layer
!! (with \f$ B_z = 0 \f$ on the top and bottom walls) to the slab eigenmode
!! \f$ h(z) = cos(s(z-z_c)) \f$, where \f$ z_c = 0.5 \f$ is the center of the slab, gives
!! \f[ s \; tan(s d / 2) = \kappa \; coth(\kappa z_1) \f]
!! where \f$ d = 0.2 \f$ m is the slab thickness and \f$ z_1 = 0.4 \f$ m is the height of the
!! lower vacuum layer. The decay time is then \f$ \tau = \mu_0 / (\eta (s^2 + \kappa^2)) \f$,
!! which gives \f$ \tau \approx 42.05 \f$ ms for this case. As the mode can vary in either
!! the \f$ x \f$ or \f$ y \f$ direction, the slowest decay time is 2-fold degenerate.
!!
!!\section doc_thickcurr_slab_ex_code Code Walk Through
!!
!! The code consists of two components. First, a helper module is defined that includes
!! a class for generating the initial condition, followed by the main run program.
!!
!!\subsection doc_thickcurr_slab_ex_code_helper Helper module
!!
!! To set the initial condition we define an implementation of \ref fem_utils::fem_interp
!! "fem_interp" for the vector potential
!! \f[ A_0 = B_0 L / \pi \; sin(\pi x / L) sin(\pi z / L) \hat{y} \f]
!! whose tangential component vanishes on the boundary of the cube. This corresponds to a magnetic
!! field
!! \f[ B_0 = B_0 \left[ -sin(\pi x / L) cos(\pi z / L) \hat{x} + cos(\pi x / L) sin(\pi z / L) \hat{z} \right] \f]
!! that threads the slab vertically, with the same \f$ x \f$ variation as the slowest decaying mode.
! START SOURCE
MODULE slab_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
!---------------------------------------------------------------------------
! Field interpolation object for initial vector potential
!---------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: slab_init_interp
  REAL(8) :: B0 = 1.d-2 ! Peak magnetic field [T]
  REAL(8) :: L = 1.d0 ! Size of domain [m]
CONTAINS
  PROCEDURE :: interp => slab_init_apply ! Reconstruct field
END TYPE slab_init_interp
CONTAINS
SUBROUTINE slab_init_apply(self,cell,f,gop,val)
CLASS(slab_init_interp), INTENT(inout) :: self ! Interpolation object
INTEGER(4), INTENT(in) :: cell ! Cell for interpolation
REAL(8), INTENT(in) :: f(:) ! Position in cell in logical coordinates
REAL(8), INTENT(in) :: gop(3,4) ! Logical gradient vectors at f
REAL(8), INTENT(out) :: val(:) ! Vector potential at f [T-m]
REAL(8) :: pt(3)
! Map logical position to physical coordinates
pt = self%mesh%log2phys(cell,f)
! Evaluate vector potential
val(1) = 0.d0
val(2) = self%B0*self%L/pi*SIN(pi*pt(1)/self%L)*SIN(pi*pt(3)/self%L)
val(3) = 0.d0
END SUBROUTINE slab_init_apply
END MODULE slab_helpers
!!\subsection doc_thickcurr_slab_ex_driver Driver program
!!
!! As usual, we start with a few module imports for required functionality.
PROGRAM ThickCurr_slab
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---FE spaces
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_getmop, oft_hcurl_grad_project, &
  oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob
!---Physics
USE thick_wall, ONLY: thickcurr_ML_hcurl, thickcurr_ML_H1grad, thickcurr_ML_hcurl_grad, &
  thickcurr_ML_lagrange, thickcurr_ML_vlagrange, thickcurr_setup, thickcurr_eigs, &
  thickcurr_run, thickcurr_plot
!---Self
USE slab_helpers, ONLY: slab_init_interp
IMPLICIT NONE
!!\subsubsection doc_thickcurr_slab_ex_driver_vars Local Variables
!! Next we define the local variables needed to initialize our case, run the eigenvalue
!! and time-dependent solves, and perform post-processing.
!---FE spaces
TYPE(oft_ml_fem_type), TARGET :: ML_lagrange,ML_H1,ML_H1grad,ML_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_vlagrange,ML_hcurl_grad
!---Mass matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Local variables
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_hcurl_grad_gzerop) :: hcurl_grad_gzerop
TYPE(oft_hcurl_grad_zerob) :: hcurl_grad_zerob
TYPE(slab_init_interp) :: init_field
REAL(8), ALLOCATABLE :: tau(:)
INTEGER(4) :: i,ierr,io_unit,minlev
!---Input file options
INTEGER(4) :: order = 2
INTEGER(4) :: nev = 4
REAL(8) :: B0 = 1.d-2
LOGICAL :: run_eigs = .TRUE.
LOGICAL :: run_td = .TRUE.
LOGICAL :: plot_run = .FALSE.
NAMELIST/thickcurr_slab_options/order,nev,B0,run_eigs,run_td,plot_run
!!\subsubsection doc_thickcurr_slab_ex_driver_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thickcurr_slab_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
minlev=mg_mesh%mgmax+order-1 ! Only construct the highest FE level (no multigrid)
!!\subsubsection doc_thickcurr_slab_ex_driver_fe Build FE structures
!!
!! ThickCurr represents the vector potential using the full H(Curl) + Grad(H^1) space,
!! which is constructed in the same way as for \ref xmhd "MUG". The higher order Grad(H^1)
!! DOFs in vacuum regions are removed automatically by ThickCurr. Lagrange FE spaces are
!! also constructed, which are used to output fields at the mesh vertices for plotting (if
!! these are not set cell averaged values are output instead).
!--- Lagrange (for plotting)
CALL oft_lag_setup(mg_mesh,order,ML_lagrange,ML_vlag_obj=ML_vlagrange,minlev=minlev)
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_H1,minlev=minlev)
CALL h1_setup_interp(ML_H1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(ML_hcurl)
!--- Full H(Curl) + Grad(H^1) space
CALL oft_hcurl_grad_setup(ML_hcurl,ML_H1,ML_hcurl_grad,ML_H1grad,minlev=minlev)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
hcurl_grad_zerob%ML_hcurl_grad_rep=>ML_hcurl_grad
!--- Set FE spaces for ThickCurr
thickcurr_ML_hcurl=>ML_hcurl
thickcurr_ML_H1grad=>ML_H1grad
thickcurr_ML_hcurl_grad=>ML_hcurl_grad
thickcurr_ML_lagrange=>ML_lagrange
thickcurr_ML_vlagrange=>ML_vlagrange
!!\subsubsection doc_thickcurr_slab_ex_driver_tw Setup ThickCurr
!!
!! Next we setup ThickCurr, which reads the run settings from the `thickcurr_options` group
!! in the input file and the resistivity of each conducting region (and optionally a
!! preconditioner) from the `thickcurr` element in the XML input file, see
!! \ref doc_thickcurr_slab_ex_input. The operators and linear solver used for the eigenvalue
!! and time-dependent solves are then constructed. Resistivities may also be passed directly
!! using the optional `eta_reg` argument to \ref thick_wall::thickcurr_setup "thickcurr_setup".
CALL thickcurr_setup()
!!\subsubsection doc_thickcurr_slab_ex_driver_plot Perform post-processing
!!
!! To visualize the time-dependent solution once a simulation has completed the
!! \ref thick_wall::thickcurr_plot "thickcurr_plot" subroutine is used, which steps through
!! the restart files produced by \ref thick_wall::thickcurr_run "thickcurr_run" and saves the
!! vector potential (`A`), magnetic field (`B`), electric field (`E`), and current density
!! (`J`, which is zero in vacuum regions) to plot files. In
!! this example this is done automatically after the time-dependent run, but it can also be
!! run separately by setting `plot_run=T` in the `thickcurr_slab_options` input group.
IF(plot_run)THEN
  CALL thickcurr_plot()
  CALL oft_finalize
END IF
!!\subsubsection doc_thickcurr_slab_ex_driver_eigs Compute decay times
!!
!! The slowest decaying modes of the system are computed using
!! \ref thick_wall::thickcurr_eigs "thickcurr_eigs". The time step (`dt`) from the input file
!! is used as a shift for this solve and should be comparable to the time constants of interest.
!! The time constants are returned in order of decreasing decay time and, with `save_plot=.TRUE.`,
!! each mode is saved to the plot group `thickcurr_eigs` in the `eigs/` subdirectory, using the mode
!! index as the "time".
IF(run_eigs)THEN
  ALLOCATE(tau(nev))
  CALL thickcurr_eigs(nev,tau,save_plot=.TRUE.)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'Slowest decay times [ms]'
    DO i=1,nev
      WRITE(*,'(2X,I4,F12.4)')i,tau(i)*1.d3
    END DO
    WRITE(*,*)
  END IF
  DEALLOCATE(tau)
END IF
!!\subsubsection doc_thickcurr_slab_ex_driver_ic Initial conditions
!!
!! For the time-dependent run we use the usual projection times mass inverse method to map our
!! initial condition onto the H(Curl) + Grad(H^1) basis. The redundant vertex Grad(H^1) DOFs are
!! zeroed before solving. As the perfectly conducting walls hold the tangential components of
!! \f$ A \f$ on the boundary fixed in time, these components are also zeroed after the projection.
!! Otherwise small projection errors on the boundary would support a (small) static field.
IF(run_td)THEN
  !---Generate mass matrix
  NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
  CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,"none") ! Construct mass matrix with "none" BC
  !---Setup linear solver
  CALL create_cg_solver(minv)
  minv%A=>mop ! Set matrix to be solved
  minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
  CALL create_diag_pre(minv%pre) ! Setup Preconditioner
  !---Create fields for solver
  CALL ML_hcurl_grad%vec_create(u)
  CALL ML_hcurl_grad%vec_create(v)
  !---Project onto vector H(Curl) + Grad(H^1) basis
  init_field%mesh=>mg_mesh%mesh
  init_field%B0=B0
  CALL oft_hcurl_grad_project(ML_hcurl_grad%current_level,init_field,v)
  CALL hcurl_grad_gzerop%apply(v) ! Zero out redundant vertex degrees of freedom
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  CALL hcurl_grad_zerob%apply(u) ! Zero tangential components on boundary
  !---Cleanup objects used for projection
  CALL v%delete
  CALL mop%delete
  DEALLOCATE(v,mop)
  CALL minv%pre%delete
  DEALLOCATE(minv%pre)
  CALL minv%delete
  DEALLOCATE(minv)
!!\subsubsection doc_thickcurr_slab_ex_driver_run Run simulation
!!
!! Finally, the time-dependent simulation is run using
!! \ref thick_wall::thickcurr_run "thickcurr_run", which advances the vector potential with
!! the time step and number of steps given in the input file and produces restart files that
!! contain the solution at different times. On return `u` contains the solution at the final time.
!! The magnetic energy and Ohmic power dissipated in conductors are also recorded to the history
!! file `thickcurr.hist` during the simulation. Plot files are then generated from the restart
!! files as described above.
  CALL thickcurr_run(u)
  CALL u%delete
  DEALLOCATE(u)
  CALL thickcurr_plot()
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM ThickCurr_slab
! STOP SOURCE
!!
!!\section doc_thickcurr_slab_ex_input Input files
!!
!! Below is the input file used with this example. Some annotation of the options is provided
!! inline below, for more information on the available options in the `thickcurr_options` group
!! see \ref thick_wall::thickcurr_read_settings "thickcurr_read_settings".
!!
!! The mesh file (`thickcurr_slab.h5`) included with this example was generated using the script
!! `create_slab_mesh.py`, which can be modified to change the resolution or geometry.
!!
!!\verbatim
!!&runtime_options
!! ppn=1               ! Number of processors/tasks per node (used for heirarchical domain decomposition)
!! debug=0             ! Debug level (0-3)
!!/
!!
!!&mesh_options
!! meshname='slab'     ! Meshname (unused at present)
!! cad_type=0          ! Mesh format type (0 -> native)
!! nlevels=1           ! Number of total grid levels (see docs for more information)
!! nbase=1             ! Number grid levels before domain decomposition (see docs for more information)
!!/
!!
!!&native_mesh_options
!! filename='thickcurr_slab.h5' ! Mesh file
!!/
!!
!!&thickcurr_slab_options
!! order=2             ! FE order
!! nev=4               ! Number of eigenvalues to compute
!! B0=1.E-2            ! Peak initial magnetic field [T]
!! run_eigs=T          ! Compute decay times?
!! run_td=T            ! Run time-dependent simulation?
!! plot_run=F          ! Run plotting only (from existing restart files)
!!/
!!
!!&thickcurr_options
!! dt=2.E-3            ! Time step [s] (also used as the shift for eigenvalue solves)
!! theta=1.0           ! Time centering (1 -> backward Euler)
!! nsteps=150          ! Number of time steps to take
!! rst_freq=5          ! Restart file frequency
!! lin_tol=1.E-10      ! Linear solver tolerance
!!/
!!
!!&thickcurr_plot_options
!! rst_start=0         ! First restart file to plot
!! rst_end=150         ! Last restart file to plot
!!/
!!\endverbatim
!!
!!\subsection doc_thickcurr_slab_ex_input_xml XML input file
!!
!! The resistivity of each conducting region is specified using `region` elements, where `id`
!! gives the mesh region IDs and `eta` the resistivity in \f$ \Omega \cdot m \f$. All other
!! regions are treated as vacuum. A preconditioner for the linear solver is also specified, which
!! in this case uses a direct solve (LU factorization).
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./ThickCurr_slab oft.in oft_in.xml
!!\endverbatim
!!
! OFT_DOC_INCLUDE: oft_in.xml
!!
!! \note The pseudo-vacuum makes the linear system ill-conditioned for iterative methods, so a direct
!! solve is recommended. In parallel, the block-Jacobi preconditioner shown here only factors the local
!! portion of the matrix on each process and convergence degrades significantly with the number of
!! processes, so this example should be run in serial.
!!
!!\section doc_thickcurr_slab_ex_post Post-Processing
!!
!!\subsection doc_thickcurr_slab_ex_post_hist Energy history
!!
!! During the simulation the magnetic energy and Ohmic power in conductors are written to the history
!! file `thickcurr.hist`. The script `plot_thickcurr_hist.py` included with this example plots the
!! magnetic energy as a function of time, along with the expected decay of the slowest mode
!! (\f$ W \propto e^{-2 t / \tau} \f$). After a short initial transient, as the field in the vacuum
!! regions relaxes and faster modes decay, the energy follows the slowest decay time.
!!
!! \note With backward Euler time stepping (`theta=1`), a mode with decay time \f$ \tau \f$ is advanced
!! with an effective decay time of \f$ \Delta t / ln(1 + \Delta t / \tau) \f$, which is approximately
!! 43.0 ms (2% longer than the exact value) for the time step used here. Crank-Nicolson (`theta=0.5`)
!! is more accurate for the conducting regions, but only weakly damps the stiff pseudo-vacuum modes,
!! which can then oscillate from step to step.
!!
!!\verbatim
!!~$ python plot_thickcurr_hist.py
!!\endverbatim
!!
!!\subsection doc_thickcurr_slab_ex_post_plot Visualization in ParaView
!!
!! Once the run is complete, the `OFT_build_xdmf` script (see \ref build_xdmf) can be used to generate
!! `*.xmf` files that can be loaded by [ParaView](https://www.paraview.org/),
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), or other visualization programs.
!! Plot files are produced in two locations, which must be processed separately:
!!  - The run directory, containing the `thickcurr` group: time-dependent fields (`A`, `B`, `E`, and `J`)
!!    at each saved time
!!  - The `eigs` subdirectory, containing the `thickcurr_eigs` group: eigenmode fields, where the "time"
!!    index corresponds to the mode number (ordered from slowest to fastest decay)
!!
!!\verbatim
!!~$ OFT_build_xdmf
!!~$ cd eigs
!!~$ OFT_build_xdmf
!!\endverbatim
!!
!! The resulting `thickcurr_vmesh-*.xmf` (or `thickcurr_eigs_vmesh-*.xmf`) files can then be opened in
!! ParaView as a time series. The `*_smesh-*.xmf` files contain the boundary surface mesh.
!!
!! In ParaView the current density (`J`) is a good choice for visualizing the eddy currents in the slab,
!! for example using a "Slice" or "Glyph" filter. `J` is computed from the electric field in each
!! conducting region (\f$ J = E / \eta \f$) and is zero in vacuum regions. Surface currents on the
!! perfectly conducting walls are not included.
