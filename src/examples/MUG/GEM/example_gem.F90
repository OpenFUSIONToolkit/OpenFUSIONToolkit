!!Example: Slab Reconnection    {#doc_ex7}
!!============================
!!
!![TOC]
!!
!!This example demonstrates the use of the \ref xmhd "extended MHD" module in PSI-Tet. In
!!this example reconnection in a periodic-slab current sheet will be simulated.
!!
!!\section doc_ex7_code_helper Helper module
!!
!! Need docs
! START SOURCE
MODULE GEM_helpers
USE oft_base
USE oft_mesh_type, ONLY: mesh
USE fem_utils, ONLY: fem_interp
USE mhd_utils, ONLY: mu0
USE oft_la_base, ONLY: oft_vector
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: oft_lag_rinterp
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1curl_zerob, oft_h1_rinterp
USE diagnostic, ONLY: flux_probe
USE xmhd, ONLY: oft_xmhd_driver, oft_xmhd_probe, xmhd_sub_fields
IMPLICIT NONE
!---------------------------------------------------------------------------
! CLASS magnetospheric_interpolator
!---------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: GEM_interp
  REAL(8) :: dpsi = 1.d-1 !< Amplitude of flux perturbation
  REAL(8) :: den_inf = 2.d-1 !< Density at infinity
  REAL(8) :: lam = 0.5d0 !< Wavelength of perturbation
  REAL(8) :: Lx = 25.6 !< Length of domain in x-direction
  REAL(8) :: Lz = 12.8 !< Length of domain in z-direction
  CHARACTER(LEN=1) :: field = 'n' !< Field component to initialize
CONTAINS
  !> Reconstruct field
  PROCEDURE :: interp => GEM_interp_apply
END TYPE GEM_interp
!------------------------------------------------------------------------------
! CLASS GEM_probe
!------------------------------------------------------------------------------
!> Reduced MHD probe object for HIT-SI diagnostics
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_xmhd_probe) :: GEM_probe
  INTEGER(4) :: io_unit !< I/O unit for history file
  LOGICAL :: initialized = .FALSE. !< Flag to indicate setup has been called
  TYPE(oft_h1_rinterp), POINTER :: Bfield => NULL() !< Magnetic field interpolation class
  TYPE(flux_probe) :: flux_probe !< Synthetic flux probe
CONTAINS
  !> Extract probe signals
  PROCEDURE :: apply => GEM_probe_apply
END TYPE GEM_probe
CONTAINS
!!\subsection doc_ex7_code_helper_interp Initial condition field interpolator
!!
!! Need docs
SUBROUTINE GEM_interp_apply(self,cell,f,gop,val)
CLASS(GEM_interp), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: cell
REAL(8), INTENT(in) :: f(:)
REAL(8), INTENT(in) :: gop(3,4)
REAL(8), INTENT(out) :: val(:)
REAL(8) :: pt(3),Beq(3),Bper(3)
! Map logical positionto physical coordinates
pt = mesh%log2phys(cell,f)
! Return requested field evaluated at "(cell,f) -> pt"
SELECT CASE(self%field)
  CASE('n') ! Density
    val = 1.d0/(COSH(pt(3)/self%lam))**2 + self%den_inf
  CASE('b') ! Magnetic field
    Beq = (/TANH(pt(3)/self%lam), 0.d0, 0.d0/)
    Bper(1) = -pi*COS(2.d0*pi*pt(1)/self%Lx)*SIN(pi*pt(3)/self%Lz)/self%Lz
    Bper(2) = 0.d0
    Bper(3) = 2.d0*pi*SIN(2.d0*pi*pt(1)/self%Lx)*COS(pi*pt(3)/self%Lz)/self%Lx
    val = Beq + self%dpsi*Bper
  CASE DEFAULT
    CALL oft_abort('Unknown field component','GEM_interp_apply',__FILE__)
END SELECT
END SUBROUTINE GEM_interp_apply
!!\subsection doc_ex7_code_helper_probe Reconnected flux probe
!!
!! Need docs
SUBROUTINE GEM_probe_apply(self,sub_fields,t)
CLASS(GEM_probe), INTENT(inout) :: self
type(xmhd_sub_fields), intent(inout) :: sub_fields
REAL(8), INTENT(in) :: t
REAL(8) :: tflux
!---History file variables
INTEGER(4), PARAMETER :: nhist_fields=2
REAL(4), DIMENSION(nhist_fields) :: output
CHARACTER(LEN=20), DIMENSION(nhist_fields) :: hist_fields
CHARACTER(LEN=2), DIMENSION(nhist_fields) :: hist_forms
!---------------------------------------------------------------------------
! Setup if necessary
!---------------------------------------------------------------------------
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%Bfield)
  !---Setup internal flux probe
  self%flux_probe%hcpc=(/0.d0, 0.d0, 8.d0/)
  self%flux_probe%hcpv=(/0.125d0, 0.0d0, 0.0d0/)
  CALL self%flux_probe%setup
  self%flux_probe%B=>self%Bfield
  !---Setup history file I/O
  IF(oft_env%head_proc)THEN
    hist_fields(1)='time'
    hist_forms='r4'
    hist_fields(2)='flux'
    OPEN(NEWUNIT=self%io_unit,FILE='gem_flux.hist',FORM='UNFORMATTED')
    WRITE(self%io_unit)nhist_fields
    WRITE(self%io_unit)hist_fields
    WRITE(self%io_unit)hist_forms
  END IF
  self%initialized=.TRUE.
END IF
!---------------------------------------------------------------------------
! Sample signals and save to history file
!---------------------------------------------------------------------------
!---Setup interpolator
self%Bfield%u=>sub_fields%B
CALL self%Bfield%setup
!---Sample flux
CALL self%flux_probe%eval(tflux)
!---Save results
output=(/t,tflux/)
IF(oft_env%head_proc)THEN
  WRITE(self%io_unit)output
  FLUSH(self%io_unit)
END IF
END SUBROUTINE GEM_probe_apply
END MODULE GEM_helpers
!!\section doc_ex6_code_driver Driver program
!!
!! Need docs
PROGRAM example_gem
USE oft_base
!--Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: h1_setup_interp, h1_getmop, oft_h1_project, h1grad_zerop, &
  oft_h1_rinterp, oft_h1_cinterp
!---Physics
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, temp_floor, den_floor, den_scale, &
  xmhd_sub_fields
!---Self
USE GEM_helpers, ONLY: GEM_interp, GEM_probe
IMPLICIT NONE
!---Fixed scaling parameters
REAL(8), PARAMETER :: N0 = 5.196374d16
REAL(8), PARAMETER :: T0 = 1.d2
REAL(8), PARAMETER :: B0 = 0.002045692328575d0
!---Mass matrix solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Local variables
TYPE(xmhd_sub_fields) :: ic_fields
TYPE(oft_h1_rinterp), TARGET :: bfield
TYPE(oft_h1_cinterp), TARGET :: jfield
TYPE(GEM_interp), TARGET :: GEM_field
TYPE(GEM_probe) :: GEM_probes
REAL(8), POINTER :: vals(:),vec_vals(:,:)
INTEGER(4) :: i,ierr,io_unit
!---Input file options
INTEGER(4) :: order = 2
INTEGER(4) :: minlev = 1
LOGICAL :: plot_run = .FALSE.
LOGICAL :: view_ic = .FALSE.
NAMELIST/test_mhd_options/order,minlev,plot_run,view_ic
!!\subsection doc_ex6_code_driver_init Grid and FE setup
!!
!! Need docs
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_mhd_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
rgrnd=(/0.d0,0.d0,1.d0/)
CALL multigrid_construct
!---Setup I/0
IF(view_ic.OR.plot_run)CALL mesh%setup_io(order)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order,minlev)
CALL lag_setup_interp
!---H1(Curl) subspace
CALL oft_hcurl_setup(order,minlev)
CALL hcurl_setup_interp
!---H1(Grad) subspace
CALL oft_h0_setup(order+1,minlev)
CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order,minlev)
CALL h1_setup_interp
!---------------------------------------------------------------------------
! If plot run, just run and finish
!---------------------------------------------------------------------------
IF(plot_run)THEN
  CALL xmhd_plot(GEM_probes)
  CALL oft_finalize
END IF
!!\subsection doc_ex6_code_driver_ic Initial conditions
!!
!! Need docs
!---------------------------------------------------------------------------
! Set constant initial temperature
!---------------------------------------------------------------------------
CALL oft_lag_create(ic_fields%Ti)
CALL ic_fields%Ti%set(T0)
!---------------------------------------------------------------------------
! Set zero initial velocity
!---------------------------------------------------------------------------
CALL oft_lag_vcreate(ic_fields%V)
CALL ic_fields%V%set(0.d0)
!---------------------------------------------------------------------------
! Set intial density from analytic definition
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_getmop(mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL oft_lag_create(u)
CALL oft_lag_create(v)
!---Project onto scalar Lagrange basis
GEM_field%field='n'
CALL oft_lag_project(GEM_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create density field and set values
CALL oft_lag_create(ic_fields%Ne)
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
CALL h1_getmop(mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL oft_h1_create(u)
CALL oft_h1_create(v)
!---Project onto vector H(Curl) basis
GEM_field%field='b'
CALL oft_h1_project(GEM_field,v)
CALL h1grad_zerop(v) ! Zero out redundant vertex degrees of freedom
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Create magnetic field and set values
CALL oft_h1_create(ic_fields%B)
CALL ic_fields%B%add(0.d0,1.d0,u)
CALL ic_fields%B%scale(B0) ! Scale to desired value
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
! Save initial conditions if desired
!---------------------------------------------------------------------------
IF(view_ic)THEN
  !---Set up output arrays for plotting
  NULLIFY(vals)
  ALLOCATE(vec_vals(3,ic_fields%Ne%n))
  !---Plot density
  CALL ic_fields%Ne%get_local(vals) ! Fetch local values
  CALL mesh%save_vertex_scalar(vals,'N0') ! Add field to plotting file
  !---Plot temperature
  CALL ic_fields%Ti%get_local(vals) ! Fetch local values
  CALL mesh%save_vertex_scalar(vals,'T0') ! Add field to plotting file
  !---Plot velocity
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL ic_fields%V%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'V0') ! Add field to plotting file
  !---------------------------------------------------------------------------
  ! Project B and J for plotting
  !---------------------------------------------------------------------------
  !---Generate mass matrix
  NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
  CALL oft_lag_vgetmop(mop,"none") ! Construct mass matrix with "none" BC
  !---Setup linear solver
  CALL create_cg_solver(minv)
  minv%A=>mop ! Set matrix to be solved
  minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
  CALL create_diag_pre(minv%pre) ! Setup Preconditioner
  !---Create fields for solver
  CALL oft_lag_vcreate(u)
  CALL oft_lag_vcreate(v)
  !---Project B onto vector Lagrange basis
  bfield%u=>ic_fields%B
  CALL bfield%setup
  CALL oft_lag_vproject(bfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected magnetic field
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'B0') ! Add field to plotting file
  !---Project J onto vector Lagrange basis
  jfield%u=>ic_fields%B
  CALL jfield%setup
  CALL oft_lag_vproject(jfield,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  !---Retrieve and save projected current density
  DO i=1,3
    vals=>vec_vals(i,:)
    CALL u%get_local(vals,i)
  END DO
  CALL mesh%save_vertex_vector(vec_vals,'J0') ! Add field to plotting file
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
!!\subsection doc_ex6_code_driver_run Run simulation
!!
!! Need docs
xmhd_minlev=minlev  ! Set minimum level for multigrid preconditioning
temp_floor=T0*1.d-2 ! Set temperature floor
den_floor=N0*1.d-2  ! Set density floor
den_scale=N0        ! Set density scale
oft_env%pm=.FALSE.      ! Do not show linear iteration progress
!---Run simulation
CALL xmhd_run(ic_fields)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM example_gem
! STOP SOURCE
!!
!!\section doc_ex7_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! As with \ref doc_ex6 "Example 6" this example should only be run with multiple processes.
!! Some annotation of the options is provided inline below, for more information on the
!! available options in the \c xmhd_options group see \ref xmhd::xmhd_plot "xmhd_plot".
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='slab'
!! cad_type=1
!! nlevels=2
!! nbase=1
!!/
!!
!!&t3d_options
!! filename='slab_gem.t3d'
!! inpname='slab_gem.inp'
!! reflect='xyz'     ! Reflect input grid in all directions
!! ref_per=T,T,F     ! Make grid periodic in the X,Y directions
!!/
!!
!!&test_mhd_options
!! order=2           ! FE order
!! minlev=2          ! Minimum level for MG preconditioning
!! view_ic=F         ! View initial conditions but do not run simulation
!! plot_run=F        ! Run plotting instead of simulation
!!/
!!
!!&xmhd_options
!! mu_ion=1.         ! Ion mass (atomic units)
!! me_factor=73.45   ! Electron mass factor
!! xmhd_hall=T       ! Include Hall physics
!! xmhd_ohmic=T      ! Include Ohmic heating
!! xmhd_visc_heat=T  ! Include viscous heating
!! bbc='bc'          ! Perfectly-conducting BC for B-field
!! vbc='all'         ! Zero-flow BC for velocity
!! nbc='n'           ! Neumann BC for density
!! tbc='n'           ! Neumann BC for temperature
!! dt=2.e-7          ! Maximum time step
!! eta=978.7         ! Constant resistivity
!! visc_type='iso'   ! Use isotropic viscosity tensor
!! nu_par=9877.0     ! Fluid viscosity
!! d_dens=10.        ! Density diffusion
!! kappa_par=3914.9  ! Parallel thermal conduction (fixed)
!! kappa_perp=3914.9 ! Perpendicular thermal conduction (fixed)
!! nsteps=2000       ! Number of time steps to take
!! rst_ind=0         ! Index of file to restart from (0 -> use subroutine arguments)
!! rst_freq=10       ! Restart file frequency
!! xmhd_mfnk=T       ! Use matrix-free method
!! lin_tol=1.E-8     ! Linear solver tolerance
!! nl_tol=1.E-6      ! Non-linear solver tolerance
!! nu_xmhd=0,20,5    ! Number of smoother iterations for default preconditioner
!! ittarget=30       ! Target for # of linear iterations per time step
!! xmhd_prefreq=20   ! Preconditioner update frequency
!!/
!!\endverbatim
!!
!!\subsection doc_ex7_input_plot Post-Processing options
!!
!! When running the code for post-processing additional run time options are available.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-6
!! rst_start=0
!! rst_end=1000
!!/
!!
!! \image html example_gem-result.png "Resulting current distribution for the first eigenmode"
!!
!!\endverbatim
