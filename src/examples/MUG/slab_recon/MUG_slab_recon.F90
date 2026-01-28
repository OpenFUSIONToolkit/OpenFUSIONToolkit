!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Slab Reconnection    {#doc_mug_recon_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref xmhd "extended MHD" module in OFT. In
!! this example reconnection in a periodic-slab current sheet will be simulated following
!! the setup used to the [GEM Challange](https://doi.org/10.1029/1999JA900449).
!!
!!\section doc_mug_recon_ex_code_helper Code Walk Through
!!
!! The code consists of two components. First, a helper module is defined that includes
!! class implementations and other components for initialization and post-processing and 
!! the main run program.
!!
!!\section doc_mug_recon_ex_code_helper Helper module
!!
!! As usual, we start with a few module imports for required functionality.
! START SOURCE
MODULE GEM_helpers
USE oft_base
USE oft_io, ONLY: oft_bin_file
USE fem_utils, ONLY: fem_interp
USE mhd_utils, ONLY: mu0
USE oft_la_base, ONLY: oft_vector
USE oft_lag_operators, ONLY: oft_lag_rinterp
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp
USE diagnostic, ONLY: flux_probe
USE xmhd, ONLY: oft_xmhd_driver, oft_xmhd_probe, xmhd_sub_fields, xmhd_ML_hcurl_grad, xmhd_ML_hcurl
IMPLICIT NONE
!!\subsection doc_mug_recon_ex_code_helper_classes Class definitions
!!
!! Now we define a class for generating the initial condition through an implementation
!! of \ref fem_utils::fem_interp "fem_interp" for the current sheet and perturbation field.
!! Additionally, we define a class by extending the abstract \ref xmhd::oft_xmhd_probe "oft_xmhd_probe"
!! class, which will be used during post-processing (\ref xmhd::xmhd_plot "xmhd_plot") to evaluate
!! the reconnected flux as a function of time.
!---------------------------------------------------------------------------
! Field interpolation object for intial conditions
!---------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: GEM_interp
  REAL(8) :: dpsi = 1.d-1 ! Amplitude of flux perturbation
  REAL(8) :: den_inf = 2.d-1 ! Density at infinity
  REAL(8) :: lam = 0.5d0 ! Wavelength of perturbation
  REAL(8) :: Lx = 25.6 ! Length of domain in x-direction
  REAL(8) :: Lz = 12.8 ! Length of domain in z-direction
  CHARACTER(LEN=2) :: field = 'n' ! Field component to initialize
CONTAINS
  PROCEDURE :: interp => GEM_interp_apply ! Reconstruct field
END TYPE GEM_interp
!------------------------------------------------------------------------------
!> Probe object for evaluating reconnected flux
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_xmhd_probe) :: GEM_probe
  INTEGER(4) :: io_unit ! I/O unit for history file
  LOGICAL :: initialized = .FALSE. ! Flag to indicate setup has been called
  TYPE(oft_hcurl_grad_rinterp), POINTER :: Bfield => NULL() ! Magnetic field interpolation class
  TYPE(flux_probe) :: flux_probe ! Synthetic flux probe
  TYPE(oft_bin_file) :: flux_hist ! History file object
CONTAINS
  PROCEDURE :: apply => GEM_probe_apply ! Sample probe signals
END TYPE GEM_probe
CONTAINS
!!\subsection doc_mug_recon_ex_code_helper_interp Initial condition field interpolator
!!
!! To set the initial conditions we define a single interpolation object/function, with
!! a switch for each of the fields with non-uniform initial distributions. Uniform fields
!! are set in the main run program as they do not require FE projection.
!!
!! The non-uniform initial conditions for this case is given by
!! \f[ n_0 = (cosh(z / \lambda))^{-2} + n_{\infty} \f]
!! \f[ B_0 = tanh(z / \lambda) \hat{x} \f]
!!
!! The perturbation is then given by
!! \f[ \delta B = -\frac{\pi}{L_z} cos(2 \pi x / L_x) sin(\pi z / L_z) \hat{x} + \frac{2 \pi}{L_x} sin(2 \pi x / L_x) cos(\pi z / L_z) \hat{x} \f]
SUBROUTINE GEM_interp_apply(self,cell,f,gop,val)
CLASS(GEM_interp), INTENT(inout) :: self ! Needs docs
INTEGER(4), INTENT(in) :: cell ! Needs docs
REAL(8), INTENT(in) :: f(:) ! Needs docs
REAL(8), INTENT(in) :: gop(3,4) ! Needs docs
REAL(8), INTENT(out) :: val(:) ! Needs docs
REAL(8) :: pt(3),Beq(3),Bper(3)
! Map logical positionto physical coordinates
pt = self%mesh%log2phys(cell,f)
! Return requested field evaluated at "(cell,f) -> pt"
SELECT CASE(TRIM(self%field))
  CASE('n0') ! Density
    val = 1.d0/(COSH(pt(3)/self%lam))**2 + self%den_inf
  CASE('b0') ! Equilibrium magnetic field
    Beq = (/TANH(pt(3)/self%lam), 0.d0, 0.d0/)
    val = Beq 
  CASE('db') ! Perturbed magnetic field
    Bper(1) = -pi*COS(2.d0*pi*pt(1)/self%Lx)*SIN(pi*pt(3)/self%Lz)/self%Lz
    Bper(2) = 0.d0
    Bper(3) = 2.d0*pi*SIN(2.d0*pi*pt(1)/self%Lx)*COS(pi*pt(3)/self%Lz)/self%Lx
    val = Bper
  CASE DEFAULT
    CALL oft_abort('Unknown field component','GEM_interp_apply',__FILE__)
END SELECT
END SUBROUTINE GEM_interp_apply
!!\subsection doc_mug_recon_ex_code_helper_probe Reconnected flux probe
!!
!! Arbitrary signals can be computed for post-processing by passing an
!! implementation of \ref xmhd::oft_xmhd_probe "oft_xmhd_probe" to \ref xmhd::xmhd_plot "xmhd_plot".
!! In this case we compute the flux in the half plane (\f$ z>0; x=0 \f$), which can then be used
!! to compute the "reconnected flux" as \f$ \Psi(t_0) - Psi(t) \f$. On the first call a
!! \ref diagnostic::flux_probe "flux_probe" object is initialized and used on subsequent calls,
!! whose result it save to a history file.
SUBROUTINE GEM_probe_apply(self,sub_fields,t)
CLASS(GEM_probe), INTENT(inout) :: self ! Needs docs
type(xmhd_sub_fields), intent(inout) :: sub_fields ! Needs docs
REAL(8), INTENT(in) :: t
REAL(8) :: tflux
!---History file variables
REAL(4), DIMENSION(2) :: output
!---------------------------------------------------------------------------
! Setup if necessary
!---------------------------------------------------------------------------
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%Bfield)
  !---Setup internal flux probe
  self%flux_probe%hcpc=(/0.d0, 0.d0, 8.d0/)
  self%flux_probe%hcpv=(/0.125d0, 0.0d0, 0.0d0/)
  self%flux_probe%mesh=>xmhd_ML_hcurl%ml_mesh%mesh
  CALL self%flux_probe%setup
  self%flux_probe%B=>self%Bfield
  !---Setup history file I/O
  IF(oft_env%head_proc)THEN
    self%flux_hist%filedesc = 'GEM example reconnected flux'
    CALL self%flux_hist%setup('gem_flux.hist')
    CALL self%flux_hist%add_field('time', 'r4', desc="Simulation time [s]")
    CALL self%flux_hist%add_field('flux', 'r4', desc="Reconnected flux [Wb]")
    CALL self%flux_hist%write_header
  END IF
  self%initialized=.TRUE.
END IF
!---------------------------------------------------------------------------
! Sample signals and save to history file
!---------------------------------------------------------------------------
!---Setup interpolator
self%Bfield%u=>sub_fields%B
CALL self%Bfield%setup(xmhd_ML_hcurl_grad%current_level)
!---Sample flux
CALL self%flux_probe%eval(tflux)
!---Save results
output=REAL([t,tflux],4)
IF(oft_env%head_proc)THEN
  CALL self%flux_hist%open
  CALL self%flux_hist%write(data_r4=output)
  CALL self%flux_hist%close
END IF
END SUBROUTINE GEM_probe_apply
END MODULE GEM_helpers
!!\section doc_mug_recon_ex_driver Driver program
!!
!! With supporting infrasturcture defined, we now follow a similar structure to \ref doc_mug_sph_ex1 and \ref doc_mug_sph_ex2
!! for the main run program.
PROGRAM MUG_slab_recon
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, &
  create_ilu_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project
!---H1(Grad) FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, oft_hcurl_grad_project, &
  oft_hcurl_grad_gzerop, oft_hcurl_grad_rinterp, oft_hcurl_grad_cinterp
!---Physics
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, temp_floor, den_floor, den_scale, &
  xmhd_sub_fields, xmhd_lin_run, xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, xmhd_ML_lagrange, xmhd_ML_vlagrange
!---Self
USE GEM_helpers, ONLY: GEM_interp, GEM_probe
IMPLICIT NONE
!!\subsection doc_mug_sph_ex2_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve and post-processing. Compared to previous
!! examples, we now have specialized initial condition and post-processing
!! objects.
!---Fixed scaling parameters
REAL(8), PARAMETER :: N0 = 5.196374d16
REAL(8), PARAMETER :: T0 = 1.d2
REAL(8), PARAMETER :: B0 = 0.002045692328575d0
!---Mass matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Local variables
TYPE(multigrid_mesh) :: mg_mesh
TYPE(xdmf_plot_file) :: plot_file
TYPE(xmhd_sub_fields) :: ic_fields,pert_fields
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
TYPE(oft_hcurl_grad_rinterp), TARGET :: bfield
TYPE(oft_hcurl_grad_cinterp), TARGET :: jfield
TYPE(GEM_interp), TARGET :: GEM_field
TYPE(GEM_probe) :: GEM_probes
REAL(8), POINTER :: vals(:),vec_vals(:,:)
INTEGER(4) :: i,ierr,io_unit
!---Input file options
INTEGER(4) :: order = 2
INTEGER(4) :: minlev = 1
REAL(8) :: db = 1.d-1
LOGICAL :: pm = .FALSE.
LOGICAL :: linear = .FALSE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: view_ic = .FALSE.
NAMELIST/slab_recon_options/order,minlev,linear,db,plot_run,view_ic,pm
!!\subsection doc_mug_recon_ex_driver_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,slab_recon_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh,[0.d0,0.d0,1.d0])
!---Setup I/0
IF(view_ic)THEN
  CALL plot_file%setup("gem")
  CALL mg_mesh%mesh%setup_io(plot_file,order)
END IF
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!--- Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!--- Grad(H^1) subspace
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
!--- H(Curl) subspace
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
!--- Full H(Curl) space
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev=minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
!!\subsection doc_mug_recon_ex_driver_plot Perform post-processing
!!
!! To visualize the solution fields once a simulation has completed the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_run "xmhd_run" and generates plot files, and in our case
!! probe signals for the half-plane flux at evenly spaced points in time as specified in
!! the input file, see \ref xmhd::xmhd_plot "xmhd_plot" and \ref doc_mug_recon_ex_post.
IF(plot_run)THEN
  CALL xmhd_plot(GEM_probes)
  CALL oft_finalize
END IF
!!\subsection doc_mug_recon_ex_driver_ic Initial conditions
!!
!! For this simulation we set every field using either uniform fields or using our interpolation
!! object. As all non-zero uniform fields utilize a nodel basis, they can be set simply using `%set()`
!! on the relevant vector object. For non-uniform fields, we use the usual projection times mass inverse
!! method to map the field onto a given FE representation.
!---Set constant initial temperature
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
CALL ic_fields%Ti%set(T0)
!---Set zero initial velocity
CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
CALL ic_fields%V%set(0.d0)
!---------------------------------------------------------------------------
! Set intial density from analytic definition
!---------------------------------------------------------------------------
GEM_field%mesh=>mg_mesh%mesh
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL xmhd_ML_lagrange%vec_create(u)
CALL xmhd_ML_lagrange%vec_create(v)
!---Project onto scalar Lagrange basis
GEM_field%field='n0'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,GEM_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create density field and set values
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ne)
CALL ic_fields%Ne%add(0.d0,1.d0,u)
CALL ic_fields%Ne%scale(N0) ! Scale to desired value
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
! Set intial magnetic field from analytic definition
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
! CALL create_diag_pre(minv%pre) ! Setup Preconditioner
CALL create_bjacobi_pre(minv%pre,-1)
DEALLOCATE(minv%pre%pre)
CALL create_ilu_pre(minv%pre%pre)
!---Create fields for solver
CALL xmhd_ML_hcurl_grad%vec_create(u)
CALL xmhd_ML_hcurl_grad%vec_create(v)
!---Project onto vector H(Curl) basis
GEM_field%field='b0'
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,GEM_field,v)
CALL hcurl_grad_gzerop%apply(v) ! Zero out redundant vertex degrees of freedom
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create magnetic field and set values
CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
CALL ic_fields%B%add(0.d0,1.d0,u)
CALL ic_fields%B%scale(B0) ! Scale to desired value
!---Compute perturbed magnetic field
GEM_field%field='db'
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,GEM_field,v)
CALL hcurl_grad_gzerop%apply(v) ! Zero out redundant vertex degrees of freedom
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL xmhd_ML_hcurl_grad%vec_create(pert_fields%B)
CALL pert_fields%B%add(0.d0,1.d0,u)
CALL pert_fields%B%scale(B0*db) ! Scale to desired value
!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre%pre)
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)
!!\subsubsection doc_mug_recon_ex_driver_ic_plot Plot initial conditions
!!
!! We also add the ability to visualize the initial conditions before running the simulation
!! to verify expected distribution.
IF(view_ic)THEN
  !---Set up output arrays for plotting
  NULLIFY(vals)
  ALLOCATE(vec_vals(3,ic_fields%Ne%n))
  !---Plot density
  CALL ic_fields%Ne%get_local(vals) ! Fetch local values
  CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'N0') ! Add field to plotting file
  !---Plot temperature
  CALL ic_fields%Ti%get_local(vals) ! Fetch local values
  CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'T0') ! Add field to plotting file
  !---Plot velocity
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL ic_fields%V%get_local(vals,i)
  END DO
  CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'V0') ! Add field to plotting file
  !---------------------------------------------------------------------------
  ! Project B and J for plotting
  !---------------------------------------------------------------------------
  !---Generate mass matrix
  NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
  CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
  !---Setup linear solver
  CALL create_cg_solver(minv)
  minv%A=>mop ! Set matrix to be solved
  minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
  CALL create_diag_pre(minv%pre) ! Setup Preconditioner
  !---Create fields for solver
  CALL xmhd_ML_vlagrange%vec_create(u)
  CALL xmhd_ML_vlagrange%vec_create(v)
  !---Project B onto vector Lagrange basis
  bfield%u=>ic_fields%B
  CALL bfield%setup(xmhd_ML_hcurl_grad%current_level)
  CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected magnetic field
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'B0') ! Add field to plotting file
  !---Project B onto vector Lagrange basis
  bfield%u=>pert_fields%B
  CALL bfield%setup(xmhd_ML_hcurl_grad%current_level)
  CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected magnetic field
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'dB') ! Add field to plotting file
  !---Project J onto vector Lagrange basis
  jfield%u=>ic_fields%B
  CALL jfield%setup(xmhd_ML_hcurl_grad%current_level)
  CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,jfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected current density
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'J0') ! Add field to plotting file
  !---Cleanup objects used for projection
  CALL u%delete ! Destroy LHS vector
  CALL v%delete ! Destroy RHS vector
  CALL mop%delete ! Destroy mass matrix
  DEALLOCATE(u,v,mop) ! Deallocate objects
  CALL minv%pre%delete ! Destroy preconditioner
  DEALLOCATE(minv%pre)
  CALL minv%delete ! Destroy solver
  DEALLOCATE(minv)
  !---Finalize enviroment
  CALL oft_finalize
END IF
!!\subsection doc_mug_recon_ex_driver_run Run simulation
!!
!! Finally, the simulation can be run using either the driver routine for linear
!! (\ref xmhd::xmhd_lin_run "xmhd_lin_run") or non-linear MHD (\ref xmhd::xmhd_run "xmhd_run").
!! These routines both advance the solution in time with the physics specified in the input file,
!! see the documentation for \ref xmhd::xmhd_lin_run "xmhd_lin_run" and \ref xmhd::xmhd_run "xmhd_run",
!! and produces restart files that contain the solution at different times.
!!
!! Several quantities are also recorded to a history file `xmhd.hist` during
!! the simulation. The data in the history file may be plotted using the script
!! `plot_mug_hist.py`, which is located in `src/utilities/scripts` or `python` directories for the
!! repository and an installation respectively.
!!
!! \note OFT plotting scripts require the python packages `numpy` and `matplotlib` as well
!! as path access to the python modules provided in `python` for installations or `src/utilities` in the repo.
xmhd_minlev=minlev  ! Set minimum level for multigrid preconditioning
den_scale=N0        ! Set density scale
oft_env%pm=pm       ! Show linear iteration progress?
IF(linear)THEN
  CALL xmhd_ML_vlagrange%vec_create(pert_fields%V)
  CALL xmhd_ML_lagrange%vec_create(pert_fields%Ti)
  CALL xmhd_ML_lagrange%vec_create(pert_fields%Ne)
  CALL xmhd_lin_run(ic_fields,pert_fields)
ELSE
  CALL ic_fields%B%add(1.d0,1.d0,pert_fields%B)
  !---Run simulation
  temp_floor=T0*1.d-2 ! Set temperature floor
  den_floor=N0*1.d-2  ! Set density floor
  CALL xmhd_run(ic_fields)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM MUG_slab_recon
! STOP SOURCE
!!
!!\section doc_mug_recon_ex_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! As with \ref doc_mug_sph_ex1 and \ref doc_mug_sph_ex2 this example should only be run
!! with multiple processes. Some annotation of the options is provided inline below, for
!! more information on the available options in the `xmhd_options` group see
!! \ref xmhd::xmhd_run "xmhd_run".
!!
!!\verbatim
!!&runtime_options
!! ppn=1               ! Number of processors/tasks per node (used for heirarchical domain decomposition)
!! debug=0             ! Debug level (0-3)
!!/
!!
!!&mesh_options
!! meshname='slab'     ! Meshname (unused at present)
!! cad_type=1          ! Mesh format type (1 -> native)
!! nlevels=2           ! Number of total grid levels (see docs for more information)
!! nbase=1             ! Number grid levels before domain decomposition (see docs for more information)
!! part_meth=2         ! Partition "uniformly" along x-axis
!!/
!!
!!&t3d_options
!! filename='slab_gem.t3d'
!! inpname='slab_gem.inp'
!! reflect='xyz'       ! Reflect input grid in all directions
!! ref_per=T,T,F       ! Make grid periodic in the X,Y directions
!!/
!!
!!&slab_recon_options
!! order=2             ! FE order
!! minlev=2            ! Minimum level for MG preconditioning
!! linear=F            ! Perform linear simulation?
!! view_ic=F           ! View initial conditions but do not run simulation
!! plot_run=F          ! Run plotting instead of simulation
!! pm=F                ! View extended linear and non-linear iteration output?
!!/
!!
!!&xmhd_options
!! mu_ion=1.           ! Ion mass (atomic units)
!! xmhd_ohmic=T        ! Include Ohmic heating
!! xmhd_visc_heat=T    ! Include viscous heating
!! xmhd_hall=F         ! Include Hall terms?
!! bbc='bc'            ! Perfectly-conducting BC for B-field
!! vbc='all'           ! Zero-flow BC for velocity
!! nbc='n'             ! Neumann BC for density
!! tbc='n'             ! Neumann BC for temperature
!! dt=8.e-7            ! Maximum time step
!! eta=742.6           ! Constant resistivity
!!! eta_hyper=50.0      ! Hyper-resistivity
!!! me_factor=73.446    ! M_e/M_i = 25.0 (as per GEM challenge paper)
!! visc_type='iso'     ! Use isotropic viscosity tensor
!! nu_par=7425.9       ! Fluid viscosity
!! d_dens=50.          ! Density diffusion
!! kappa_par=2943.4    ! Parallel thermal conduction (fixed)
!! kappa_perp=2943.4   ! Perpendicular thermal conduction (fixed)
!! nsteps=1000         ! Number of time steps to take
!! rst_ind=0           ! Index of file to restart from (0 -> use subroutine arguments)
!! rst_freq=10         ! Restart file frequency
!! xmhd_mfnk=T         ! Use matrix-free method
!! lin_tol=1.E-8       ! Linear solver tolerance
!! nl_tol=1.E-6        ! Non-linear solver tolerance
!! ittarget=30         ! Target for # of linear iterations per time step
!! xmhd_prefreq=20     ! Preconditioner update frequency
!! xmhd_nparts=0,20,40 ! Number of parts for local matrix decomposition
!! nl_update=3         ! # of NL iterations that causes preconditioner update
!!/
!!\endverbatim
!!
!!\subsection doc_mug_recon_ex_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate multi-grid preconditioner
!! that uses GRMES+Block-Jacobi+LU on each level.
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./MUG_slab_recon oft.in oft_in.xml
!!\endverbatim
!!
!!```xml
!!<oft>
!!  <xmhd>
!!    <pre type="mg">
!!      <smoother direction="up">
!!        <solver type="gmres">
!!          <its>2</its>
!!          <nrits>2</nrits>
!!          <pre type="block_jacobi">
!!            <nlocal>2</nlocal>
!!            <solver type="lu"></solver>
!!          </pre>
!!        </solver>
!!      </smoother>
!!    </pre>
!!  </xmhd>
!!</oft>
!!```
!!
!!\section doc_mug_recon_ex_post Post-Processing options
!!
!! During the simulation the evolution of some global quantities is written to the history file `xmhd.hist`. A utility
!! script (`plot_mug_hist.py` located in `bin` after installation) is included as part of OFT for plotting the signals in this file.
!! This is useful for monitoring general progress of the simulation as well as numerical parameters like the iteration
!! count and solver time. Using this script we can see the instability growth and relaxation of the initial magnetic
!! equilibrium to the lowest energy state.
!!
!!\verbatim
!!~$ python plot_mug_hist.py xmhd.hist
!!\endverbatim
!!
!!\subsection doc_mug_recon_ex_post_plot Creating plot files
!!
!! To generate 3D plots, and perform additional diagnostic sampling (see \ref xmhd::oft_xmhd_probe "oft_xmhd_probe"), a plot run can
!! be performed by setting `plot_run=T` in the `slab_recon_options` input group, which calls \ref xmhd::xmhd_plot "xmhd_plot". With this option
!! additional run time options are available in the `xmhd_plot_options` group that control how restart files are sampled for plotting.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-6
!! rst_start=0
!! rst_end=1000
!!/
!!\endverbatim
!!
!! Once the post-processing run is complete `bin/build_xdmf.py` can be used to generate `*.xmf` files that can be loaded by
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!!
!! \image html MUG_gem_ex-J.png "Resulting current distribution at final time"
!! \image html MUG_gem_ex-V.png "Resulting velocity magnitude distribution at final time"
