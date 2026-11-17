!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!ThickCurr Example: Eddy currents in a conducting torus    {#doc_thickcurr_torus_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref thick_wall "ThickCurr" module in OFT to model
!! eddy currents in a solid conducting torus (major radius \f$ R_0 = 1 \f$ m, minor radius
!! \f$ a = 0.5 \f$ m, \f$ \eta = 7 \times 10^{-7} \Omega \cdot m \f$) at the center of a spherical
!! domain (radius 3 m) bounded by a perfectly conducting wall, with the remainder of the domain treated
!! as vacuum. In this example we will:
!!
!!  1. Compute the slowest decay times (L/R time constants) of the system with an eigenvalue solve
!!  2. Simulate the response of the torus to the sudden removal of a uniform vertical magnetic field,
!!     using a probe object to record the magnetic flux through the hole of the torus
!!  3. Generate plot files of the eigenmodes and time-dependent fields for ParaView
!!
!! This example follows a similar structure to \ref doc_thickcurr_slab_ex, see that example for
!! additional discussion of the setup and options.
!!
!!\section doc_thickcurr_torus_ex_estimate Approximate decay time
!!
!! The slowest decaying mode of this system corresponds to a toroidal current in the torus, which
!! links the hole of the torus. While there is no simple analytic solution for this case, the decay
!! time can be estimated from the inductance and resistance of a torus carrying a uniform toroidal
!! current
!! \f[ L \approx \mu_0 R_0 \left[ ln \left( \frac{8 R_0}{a} \right) - \frac{7}{4} \right],
!! \; R = \frac{2 \eta R_0}{a^2} \f]
!! which gives \f$ \tau = L/R \approx 230 \f$ ms for this case. As the eigenmode adjusts its current
!! distribution to maximize the decay time, the true value is expected to be somewhat larger than this
!! estimate, while the nearby conducting wall reduces the inductance (and decay time) slightly. Using the
!! included mesh ThickCurr gives \f$ \tau \approx 240.0 \f$ ms, which agrees with the result on a coarser
!! mesh (6337 cells) to within 0.1%. The next slowest modes, with \f$ \tau \approx 81 \f$ ms, are nearly
!! degenerate.
!!
!!\section doc_thickcurr_torus_ex_mesh Mesh generation
!!
!! The mesh for this example is generated using [GMSH](https://gmsh.info/) with the script
!! `create_torus_mesh.py`, which requires the `gmsh` python package. A second-order (curved) mesh
!! is generated to accurately represent the curved surfaces of the torus and outer boundary. The
!! script also embeds a disk of radius 0.45 m in the \f$ z = 0 \f$ plane within the hole of the torus,
!! which is used to measure the magnetic flux through the hole. The resulting mesh (`thickcurr_torus.h5`,
!! 12180 cells) is included with this example.
!!
!! \note With the p-multigrid preconditioner used in this example (see \ref doc_thickcurr_torus_ex_input_xml)
!! the example requires approximately 2.7 GB of memory and takes approximately 20 minutes to run in
!! serial with SuperLU, of which approximately 13 minutes are spent in the eigenvalue solve.
!!
!!\section doc_thickcurr_torus_ex_code Code Walk Through
!!
!!\subsection doc_thickcurr_torus_ex_code_helper Helper module
!!
!! As usual, we start with a few module imports for required functionality.
! START SOURCE
MODULE torus_helpers
USE oft_base
USE oft_io, ONLY: oft_bin_file
USE fem_utils, ONLY: fem_interp
USE oft_la_base, ONLY: oft_vector
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_cinterp
USE diagnostic, ONLY: flux_probe
USE thick_wall, ONLY: oft_thickcurr_probe, thickcurr_ML_hcurl_grad
IMPLICIT NONE
!!\subsubsection doc_thickcurr_torus_ex_code_helper_classes Class definitions
!!
!! We define a class for generating the initial condition through an implementation of
!! \ref fem_utils::fem_interp "fem_interp", and a class for measuring the magnetic flux through
!! the hole of the torus by extending the abstract \ref thick_wall::oft_thickcurr_probe
!! "oft_thickcurr_probe" class. Probe objects can be passed to \ref thick_wall::thickcurr_run
!! "thickcurr_run", where they are called after each time step, or to \ref thick_wall::thickcurr_plot
!! "thickcurr_plot", where they are called for each restart file.
!---------------------------------------------------------------------------
! Field interpolation object for initial vector potential
!---------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: torus_init_interp
  REAL(8) :: B0 = 1.d-1 ! Vertical magnetic field [T]
CONTAINS
  PROCEDURE :: interp => torus_init_apply ! Reconstruct field
END TYPE torus_init_interp
!---------------------------------------------------------------------------
! Probe object for flux through the hole of the torus
!---------------------------------------------------------------------------
TYPE, EXTENDS(oft_thickcurr_probe) :: torus_flux_probe
  LOGICAL :: initialized = .FALSE. ! Flag to indicate setup has been called
  REAL(8) :: r_probe = 0.45d0 ! Radius of probe disk [m]
  TYPE(oft_hcurl_grad_cinterp), POINTER :: Bfield => NULL() ! Magnetic field interpolation object
  TYPE(flux_probe) :: flux_probe ! Synthetic flux probe
  TYPE(oft_bin_file) :: flux_hist ! History file object
CONTAINS
  PROCEDURE :: apply => torus_probe_apply ! Sample probe signals
END TYPE torus_flux_probe
CONTAINS
!!\subsubsection doc_thickcurr_torus_ex_code_helper_interp Initial condition
!!
!! The initial condition is a uniform vertical magnetic field \f$ B_0 \hat{z} \f$, with vector potential
!! \f[ A_0 = \frac{B_0}{2} \left( -y \hat{x} + x \hat{y} \right) \f]
SUBROUTINE torus_init_apply(self,cell,f,gop,val)
CLASS(torus_init_interp), INTENT(inout) :: self ! Interpolation object
INTEGER(4), INTENT(in) :: cell ! Cell for interpolation
REAL(8), INTENT(in) :: f(:) ! Position in cell in logical coordinates
REAL(8), INTENT(in) :: gop(3,4) ! Logical gradient vectors at f
REAL(8), INTENT(out) :: val(:) ! Vector potential at f [T-m]
REAL(8) :: pt(3)
! Map logical position to physical coordinates
pt = self%mesh%log2phys(cell,f)
! Evaluate vector potential
val = self%B0*[-pt(2), pt(1), 0.d0]/2.d0
END SUBROUTINE torus_init_apply
!!\subsubsection doc_thickcurr_torus_ex_code_helper_probe Flux probe
!!
!! The flux through the hole of the torus is computed using a \ref diagnostic::flux_probe "flux_probe"
!! object, which integrates \f$ B \cdot \hat{n} \f$ over the mesh faces that lie within a disk (center
!! `hcpc`, with normal direction given by `hcpv` and radius \f$ 1/|hcpv| \f$). On the first call the
!! probe is initialized and on each call the magnetic field (\f$ B = \nabla \times A \f$) is sampled
!! and the result is saved to a history file.
SUBROUTINE torus_probe_apply(self,a,t)
CLASS(torus_flux_probe), INTENT(inout) :: self ! Probe object
CLASS(oft_vector), TARGET, INTENT(inout) :: a ! Vector potential at current time
REAL(8), INTENT(in) :: t ! Current time
REAL(8) :: tflux
!---------------------------------------------------------------------------
! Setup if necessary
!---------------------------------------------------------------------------
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%Bfield)
  !---Setup internal flux probe (disk in z=0 plane centered on the origin)
  self%flux_probe%hcpc=[0.d0, 0.d0, 0.d0]
  self%flux_probe%hcpv=[0.d0, 0.d0, 1.d0/self%r_probe]
  self%flux_probe%mesh=>thickcurr_ML_hcurl_grad%ml_fields(1)%ml%ml_mesh%mesh
  CALL self%flux_probe%setup
  self%flux_probe%B=>self%Bfield
  !---Setup history file I/O
  IF(oft_env%head_proc)THEN
    CALL self%flux_hist%setup('torus_flux.hist', desc='ThickCurr torus example flux through hole')
    CALL self%flux_hist%add_field('time', 'r8', desc="Simulation time [s]")
    CALL self%flux_hist%add_field('flux', 'r8', desc="Flux through hole of torus [Wb]")
    CALL self%flux_hist%write_header
  END IF
  self%initialized=.TRUE.
END IF
!---------------------------------------------------------------------------
! Sample signals and save to history file
!---------------------------------------------------------------------------
!---Setup interpolator for B = curl(A)
self%Bfield%u=>a
CALL self%Bfield%setup(thickcurr_ML_hcurl_grad%current_level)
!---Sample flux
CALL self%flux_probe%eval(tflux)
!---Save results
IF(oft_env%head_proc)THEN
  CALL self%flux_hist%open
  CALL self%flux_hist%write(data_r8=[t,tflux])
  CALL self%flux_hist%close
END IF
END SUBROUTINE torus_probe_apply
END MODULE torus_helpers
!!\subsection doc_thickcurr_torus_ex_driver Driver program
!!
!! With supporting infrastructure defined, the main run program follows the same structure as
!! \ref doc_thickcurr_slab_ex.
PROGRAM ThickCurr_torus
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
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, &
  oft_hcurl_grad_project, oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob
!---Physics
USE thick_wall, ONLY: thickcurr_ML_hcurl, thickcurr_ML_H1grad, thickcurr_ML_hcurl_grad, &
  thickcurr_ML_lagrange, thickcurr_ML_vlagrange, thickcurr_minlev, thickcurr_setup, &
  thickcurr_eigs, thickcurr_run, thickcurr_plot
!---Self
USE torus_helpers, ONLY: torus_init_interp, torus_flux_probe
IMPLICIT NONE
!!\subsubsection doc_thickcurr_torus_ex_driver_vars Local Variables
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
TYPE(torus_init_interp) :: init_field
TYPE(torus_flux_probe) :: flux_probe
REAL(8), ALLOCATABLE :: tau(:)
INTEGER(4) :: i,ierr,io_unit,minlev,toplev
!---Input file options
INTEGER(4) :: order = 2
INTEGER(4) :: nev = 4
REAL(8) :: B0 = 1.d-1
REAL(8) :: eig_shift = 1.d-1
LOGICAL :: run_eigs = .TRUE.
LOGICAL :: run_td = .TRUE.
LOGICAL :: plot_run = .FALSE.
NAMELIST/thickcurr_torus_options/order,nev,B0,eig_shift,run_eigs,run_td,plot_run
!!\subsubsection doc_thickcurr_torus_ex_driver_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thickcurr_torus_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
!!\subsubsection doc_thickcurr_torus_ex_driver_fe Build FE structures and setup ThickCurr
!!
!! See \ref doc_thickcurr_slab_ex_driver_fe and \ref doc_thickcurr_slab_ex_driver_tw for details.
!! Unlike \ref doc_thickcurr_slab_ex, this example uses a "p-multigrid" preconditioner, where the
!! linear system for the requested FE order (`order=2`) is preconditioned using a coarse level with
!! linear (`order=1`) elements on the same mesh. To support this, the FE spaces are constructed
!! starting from the linear level (`minlev`), the interpolation operators between levels are constructed,
!! and the lowest level for multigrid is passed to ThickCurr using \ref thick_wall::thickcurr_minlev
!! "thickcurr_minlev". The preconditioner itself is defined in the XML input file
!! (see \ref doc_thickcurr_torus_ex_input_xml).
toplev=mg_mesh%mgmax+order-1 ! Highest FE level (requested order on finest mesh)
minlev=mg_mesh%mgmax         ! Linear elements on finest mesh (coarse level for p-multigrid)
!--- Lagrange (for plotting, highest level only)
CALL oft_lag_setup(mg_mesh,order,ML_lagrange,ML_vlag_obj=ML_vlagrange,minlev=toplev)
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_H1,minlev=minlev)
CALL h1_setup_interp(ML_H1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(ML_hcurl)
!--- Full H(Curl) + Grad(H^1) space
CALL oft_hcurl_grad_setup(ML_hcurl,ML_H1,ML_hcurl_grad,ML_H1grad,minlev=minlev)
CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_H1,create_full=.TRUE.)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
hcurl_grad_zerob%ML_hcurl_grad_rep=>ML_hcurl_grad
!--- Set FE spaces for ThickCurr
thickcurr_ML_hcurl=>ML_hcurl
thickcurr_ML_H1grad=>ML_H1grad
thickcurr_ML_hcurl_grad=>ML_hcurl_grad
thickcurr_ML_lagrange=>ML_lagrange
thickcurr_ML_vlagrange=>ML_vlagrange
thickcurr_minlev=minlev ! Lowest level for multigrid preconditioner
!---Setup ThickCurr
CALL thickcurr_setup()
!!\subsubsection doc_thickcurr_torus_ex_driver_plot Perform post-processing
!!
!! If `plot_run=T` plot files are generated from existing restart files using
!! \ref thick_wall::thickcurr_plot "thickcurr_plot". As for the time-dependent run, the probe object is
!! passed to sample the flux through the hole of the torus at the time of each restart file.
IF(plot_run)THEN
  CALL thickcurr_plot(flux_probe)
  CALL oft_finalize
END IF
!!\subsubsection doc_thickcurr_torus_ex_driver_eigs Compute decay times
!!
!! The slowest decay times are computed using \ref thick_wall::thickcurr_eigs "thickcurr_eigs", with
!! the eigenmodes saved to the plot group `thickcurr_eigs` in the `eigs/` subdirectory. By default, the
!! eigenvalue solve uses the time step as a shift, which must be small enough to accurately resolve the
!! time-dependent evolution. As the convergence of the eigenvalue solve is much faster with a shift
!! comparable to the slowest decay time, a separate shift (`eig_shift`) is passed here. The LHS matrices
!! and preconditioner are rebuilt for the eigenvalue solve and restored afterward for the time-dependent
!! simulation. For this example, using `eig_shift=0.1` instead of the time step (`dt=1.E-2`) reduces the
!! cost of the eigenvalue solve by more than a factor of two.
IF(run_eigs)THEN
  ALLOCATE(tau(nev))
  CALL thickcurr_eigs(nev,tau,save_plot=.TRUE.,shift=eig_shift)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'Slowest decay times [ms]'
    DO i=1,nev
      WRITE(*,'(2X,I4,F12.4)')i,tau(i)*1.d3
    END DO
    WRITE(*,*)
  END IF
  DEALLOCATE(tau)
END IF
!!\subsubsection doc_thickcurr_torus_ex_driver_ic Initial conditions
!!
!! The initial uniform vertical field is projected onto the H(Curl) + Grad(H^1) basis. Note that
!! this field has a non-zero tangential vector potential (and normal magnetic field) on the outer
!! boundary, which would be held fixed by the perfectly conducting wall, representing an externally
!! applied field. Instead, these components are zeroed after the projection, which corresponds to the
!! external field being switched off at \f$ t = 0 \f$. The field in the vacuum region then relaxes
!! rapidly, while the flux linked by the torus is initially conserved by the eddy currents induced in
!! the torus, which then decay.
IF(run_td)THEN
  !---Generate mass matrix
  NULLIFY(mop)
  CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,"none")
  !---Setup linear solver
  CALL create_cg_solver(minv)
  minv%A=>mop
  minv%its=-2
  CALL create_diag_pre(minv%pre)
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
  CALL hcurl_grad_zerob%apply(u) ! Remove external field (zero tangential components on boundary)
  !---Cleanup objects used for projection
  CALL v%delete
  CALL mop%delete
  DEALLOCATE(v,mop)
  CALL minv%pre%delete
  DEALLOCATE(minv%pre)
  CALL minv%delete
  DEALLOCATE(minv)
!!\subsubsection doc_thickcurr_torus_ex_driver_run Run simulation
!!
!! The time-dependent simulation is run using \ref thick_wall::thickcurr_run "thickcurr_run", passing
!! the probe object so that the flux through the hole of the torus is recorded to the history file
!! `torus_flux.hist` after each time step. Plot files are then generated from the restart files.
  CALL thickcurr_run(u,probes=flux_probe)
  CALL u%delete
  DEALLOCATE(u)
  CALL thickcurr_plot()
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM ThickCurr_torus
! STOP SOURCE
!!
!!\section doc_thickcurr_torus_ex_input Input files
!!
!! Below is the input file used with this example. Note that `grid_order=2` is used to make use of the
!! curved (second-order) mesh. See \ref doc_thickcurr_slab_ex_input for additional discussion.
!!
!!\verbatim
!!&runtime_options
!! ppn=1               ! Number of processors/tasks per node (used for heirarchical domain decomposition)
!! debug=0             ! Debug level (0-3)
!!/
!!
!!&mesh_options
!! meshname='torus'    ! Meshname (unused at present)
!! cad_type=0          ! Mesh format type (0 -> native)
!! nlevels=1           ! Number of total grid levels (see docs for more information)
!! nbase=1             ! Number grid levels before domain decomposition (see docs for more information)
!! grid_order=2        ! Use curved (second-order) mesh
!!/
!!
!!&native_mesh_options
!! filename='thickcurr_torus.h5' ! Mesh file
!!/
!!
!!&thickcurr_torus_options
!! order=2             ! FE order
!! nev=4               ! Number of eigenvalues to compute
!! B0=0.1              ! Initial vertical magnetic field [T]
!! eig_shift=0.1       ! Shift for eigenvalue solve [s]
!! run_eigs=T          ! Compute decay times?
!! run_td=T            ! Run time-dependent simulation?
!! plot_run=F          ! Run plotting only (from existing restart files)
!!/
!!
!!&thickcurr_options
!! dt=1.E-2            ! Time step [s]
!! theta=1.0           ! Time centering (1 -> backward Euler)
!! nsteps=100          ! Number of time steps to take
!! rst_freq=5          ! Restart file frequency
!! lin_tol=1.E-10      ! Linear solver tolerance
!!/
!!
!!&thickcurr_plot_options
!! rst_start=0         ! First restart file to plot
!! rst_end=100         ! Last restart file to plot
!!/
!!\endverbatim
!!
!!\subsection doc_thickcurr_torus_ex_input_xml XML input file
!!
!! The torus is region 2 in the mesh (the vacuum region is region 1). The preconditioner is a p-multigrid
!! V-cycle (`type="mg"`) between the second-order (`order=2`) and linear FE levels on the same mesh. On the
!! second-order level the smoother is two iterations of GMRES, preconditioned with block-Jacobi using
!! 200 local blocks, each solved with a direct (LU) solver. On the linear (coarse) level a single block
!! (`nlocal=1`) is used, which corresponds to a direct solve of the full coarse system. Compared to a direct
!! solve of the full system, this approach greatly reduces the required memory and allows larger meshes
!! to be used (eg. a direct solve on a 6337 cell mesh required ~6 GB, while the included mesh with p-multigrid
!! requires ~2.7 GB). Typically, 7-15 GMRES iterations are required per time step.
!!
!!\verbatim
!!~$ ./ThickCurr_torus oft.in oft_in.xml
!!\endverbatim
!!
! OFT_DOC_INCLUDE: oft_in.xml
!!
!! \note As the p-multigrid preconditioner uses a direct solver on the coarse level, this example should be
!! run in serial.
!!
!!\section doc_thickcurr_torus_ex_post Post-Processing
!!
!!\subsection doc_thickcurr_torus_ex_post_hist Flux history
!!
!! The script `plot_torus_flux.py` included with this example plots the flux through the hole of the
!! torus (`torus_flux.hist`) and the magnetic energy (`thickcurr.hist`) as a function of time. Initially
!! the flux decays slowly, as faster decaying modes redistribute the induced currents, after which the flux
!! decays with the slowest decay time from the eigenvalue solve. With backward Euler time stepping the
!! observed decay time is \f$ \Delta t / ln(1 + \Delta t / \tau) \approx 244.9 \f$ ms for the time step used
!! here, see \ref doc_thickcurr_slab_ex_post_hist.
!!
!!\verbatim
!!~$ python plot_torus_flux.py
!!\endverbatim
!!
!!\subsection doc_thickcurr_torus_ex_post_plot Visualization in ParaView
!!
!! As for \ref doc_thickcurr_slab_ex_post_plot, plot files are produced for the time-dependent fields in
!! the run directory and for the eigenmodes in the `eigs` subdirectory, which can be converted for
!! visualization using `OFT_build_xdmf`.
!!
!!\verbatim
!!~$ OFT_build_xdmf
!!~$ cd eigs
!!~$ OFT_build_xdmf
!!\endverbatim
!!
!! The current density (`J`) shows the toroidal eddy currents induced in the torus by the removal of
!! the vertical field.
