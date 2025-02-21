!!MUG Example: Spheromak Heating    {#doc_mug_sph_ex2}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of \ref doc_mhd_main "MUG"
!! to model self-heating of a spheromak in a unit cylinder. This process
!! provides a simple example case illustrating the basic aspects of an MHD simulation,
!! including: 1) Temperature dependent resistivity, 2) Anisotropic thermal conduction and 3)
!! Ohmic and Viscous heating.
!!
!! The dynamics in this example will be prdominetly limited to heating. However, if the
!! simulation is run long enough evolution of the equilibrium profile will be observed and
!! eventual instability due to current peaking will occur.
!!
!!\section doc_mug_sph_ex2_code Code Walk Through
!!
!! The code consists of three basic sections, required imports and variable definitions,
!! finite element setup, and system creation and solution.
!!
!!\subsection doc_mug_sph_ex2_code_inc Module Includes
!! The first thing that we must do is include the OFT modules which contain
!! the required functions and variables. It is good practice to restrict the included elements
!! to only those needed. This is done using the `ONLY:` clause to specifically include only
!! certain definitions. The exceptions to this practice are the \ref oft_local and \ref oft_base
!! modules, which contain a controlled set of commonly used elements that can be safely imported
!! as a whole.
! START SOURCE
PROGRAM MUG_sph_heat
!---Runtime
USE oft_base
!---Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct, multigrid_add_quad
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, &
  create_ilu_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_level, oft_hcurl_nlevels
USE oft_hcurl_operators, ONLY: hcurl_zerob
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_divout, h1_zeroi, h1_mc, h1curl_zerob
!---Physics
USE taylor, ONLY: taylor_hmodes, taylor_hffa, taylor_hlam
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_taxis, vel_scale, den_scale, &
  den_floor, temp_floor, xmhd_sub_fields
IMPLICIT NONE
!!\subsection doc_mug_sph_ex2_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve and post-processing.
!---H1 divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
!---Local variables
INTEGER(i4) :: ierr,io_unit
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
TYPE(xmhd_sub_fields) :: ic_fields
!---Runtime options
INTEGER(i4) :: order = 2
REAL(r8) :: b0_scale = 1.E-1_r8
REAL(r8) :: n0 = 1.d19
REAL(r8) :: t0 = 6.d0
LOGICAL :: plot_run=.FALSE.
LOGICAL :: pm=.FALSE.
NAMELIST/sph_heat_options/order,b0_scale,n0,t0,plot_run,pm
!!\subsection doc_mug_sph_ex2_code_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,sph_heat_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
rgrnd=[2.d0,0.d0,0.d0]
CALL multigrid_construct
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order, -1)
!---H1(Curl) subspace
CALL oft_hcurl_setup(order, -1)
!---H1(Grad) subspace
CALL oft_h0_setup(order+1, -1)
!---H1 full space
CALL oft_h1_setup(order, -1)
!!\subsection doc_mug_sph_ex2_code_plot Perform post-processing
!!
!! To visualize the solution fields once a simulation has completed the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_run "xmhd_run" and generates plot files, and optionally
!! probe signals at evenly spaced points in time as specified in the input file, see
!! \ref xmhd::xmhd_plot "xmhd_plot" and \ref doc_mug_sph_ex2_post.
IF(plot_run)THEN
  !---Run post-processing routine
  CALL xmhd_plot
  !---Finalize enviroment and exit
  CALL oft_finalize()
END IF
!!\subsection doc_mug_sph_ex2_ic Computing Initial Conditions
!!
!! For this simulation we only need the spheromak mode, which is the lowest
!! force-free eignstate in this geometry. As a result the initial condition
!! is stable to all types of mode activity.
CALL taylor_hmodes(1)
CALL oft_lag_set_level(oft_lagrange_nlevels)
!! As in \ref doc_mug_sph_ex1 we must transform the gauge of the Taylor
!! state solution to the appropriate magnetic field BCs. For more information
!! on this see the description in \ref doc_mug_sph_ex1_ic of that example.
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop)
CALL oft_h0_getlop(lop,"grnd")
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
!---Setup Preconditioner
CALL create_bjacobi_pre(linv%pre,-1)
DEALLOCATE(linv%pre%pre)
CALL create_ilu_pre(linv%pre%pre)
divout%solver=>linv
divout%bc=>h0_zerogrnd
divout%keep_boundary=.TRUE.
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
CALL oft_h1_create(ic_fields%B)
CALL taylor_hffa(1,oft_hcurl_level)%f%get_local(tmp)
CALL ic_fields%B%restore_local(tmp,1)
CALL divout%apply(ic_fields%B)
NULLIFY(divout%solver)
CALL divout%delete()
CALL linv%pre%pre%delete()
DEALLOCATE(linv%pre%pre)
CALL linv%pre%delete()
DEALLOCATE(linv%pre)
CALL linv%delete()
DEALLOCATE(linv)
!! Now we set initial conditions for the simulation using the computed taylor
!! state, flat temperature and density profiles and zero initial velocity. The
!! lower bound for density (\ref xmhd::den_floor "den_floor") and temperature
!! (\ref xmhd::temp_floor "temp_floor") in the simulation are also set. These
!! variables act to prevent very low density and negative temperature regions
!! in the simulation. The scale factors for the velocity (\ref xmhd::vel_scale
!! "vel_scale") and density (\ref xmhd::den_scale "den_scale") evolution
!! equations are also set. These variables are used to scale the corresponding
!! rows in the non-linear and linear operators to provide even weighting in
!! the residual calculations. In general these scale factors should be set to
!! the order of magnitude expected for the corresponding variables, \f$ km/s \f$ and
!! \f$ 10^{19} m^{-3} \f$ in this simulation.
CALL ic_fields%B%scale(b0_scale*taylor_hlam(1,oft_hcurl_level))
!---Clean up temporary matrices and fields
CALL lop%delete
DEALLOCATE(tmp,lop)
!---Create velocity field
CALL oft_lag_vcreate(ic_fields%V)
vel_scale = 1.d3
!---Create density field
CALL oft_lag_create(ic_fields%Ne)
CALL ic_fields%Ne%set(n0)
den_scale = n0
den_floor = n0*1.d-2
!---Create temperature field
CALL oft_lag_create(ic_fields%Ti)
CALL ic_fields%Ti%set(t0)
temp_floor = t0*1.d-2
!!\subsection doc_mug_sph_ex2_code_run Run Simulation
!!
!! Finally, the simulation can be run using the driver routine for non-linear
!! extended MHD (\ref xmhd::xmhd_run "xmhd_run"). This routine advances the
!! solution in time with the physics specified in the input file, see the
!! documentation for \ref xmhd::xmhd_run "xmhd_run", and produces restart files
!! that contain the solution at different times.
!!
!! Several quantities are also recorded to a history file `xmhd.hist` during
!! the simulation, including the toroidal current (where the symmetry axis is specified by \ref xmhd::xmhd_taxis
!! "xmhd_taxis"). The data in the history file may be plotted using the script
!! `plot_mug_hist.py`, which is located in `src/utilities/scripts` or `python` directories for the
!! repository and an installation respectively.
!!
!! \note OFT plotting scripts require the python packages `numpy` and `matplotlib` as well
!! as path access to the python modules provided in `python` for installations or `src/utilities` in the repo.
!!
!! To visualize the solution fields once a simulation has completed the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_run "xmhd_run" and generates plot files, and optionally
!! probe signals at evenly spaced points in time as specified in the input file, see
!! \ref xmhd::xmhd_plot "xmhd_plot".
xmhd_taxis=3
oft_env%pm=pm
!---Run simulation
CALL xmhd_run(ic_fields)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM MUG_sph_heat
! STOP SOURCE
!!
!!\section doc_mug_sph_ex2_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! As with \ref doc_ex5 "Example 5" this example should only be run with multiple processes.
!! Some annotation of the options is provided inline below, for more information on the
!! available options in the `xmhd_options` group see \ref xmhd::xmhd_plot "xmhd_plot".
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='test'
!! cad_type=0
!! nlevels=2
!! nbase=1
!! grid_order=2
!! fix_boundary=T
!!/
!!
!!&native_mesh_options
!! filename='cyl_heat.h5'
!!/
!!
!!&sph_heat_options
!! order=3
!! b0_scale=1.e-1
!! n0=1.e19
!! t0=6.
!! plot_run=F
!!/
!!
!!&xmhd_options
!! xmhd_ohmic=T      ! Include Ohmic heating
!! xmhd_visc_heat=T  ! Include viscous heating
!! bbc='bc'          ! Perfectly-conducting BC for B-field
!! vbc='all'         ! Zero-flow BC for velocity
!! nbc='d'           ! Dirichlet BC for density
!! tbc='d'           ! Dirichlet BC for temperature
!! dt=2.e-7          ! Maximum time step
!! eta=25.           ! Resistivity at reference temperature (Spitzer-like)
!! eta_temp=6.       ! Reference temperature for resistivity
!! nu_par=400.       ! Fluid viscosity
!! d_dens=10.        ! Density diffusion
!! kappa_par=1.E4    ! Parallel thermal conduction (fixed)
!! kappa_perp=1.E2   ! Perpendicular thermal conduction (fixed)
!! nsteps=2000       ! Number of time steps to take
!! rst_freq=10       ! Restart file frequency
!! lin_tol=1.E-9     ! Linear solver tolerance
!! nl_tol=1.E-5      ! Non-linear solver tolerance
!! xmhd_mfnk=T       ! Use matrix-free Jacobian operator
!! rst_ind=0         ! Index of file to restart from (0 -> use subroutine arguments)
!! ittarget=40       ! Target for # of linear iterations per time step
!! mu_ion=2.         ! Ion mass (atomic units)
!! xmhd_prefreq=20   ! Preconditioner update frequency
!!/
!!\endverbatim
!!
!!\subsection doc_mug_sph_ex2_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate ILU(0) preconditioner.
!! Currently, this preconditioner method is the suggested starting preconditioner for all
!! time-dependent MHD solves.
!!
!!```xml
!!<oft>
!!  <xmhd>
!!    <pre type="gmres">
!!      <its>8</its>
!!      <nrits>8</nrits>
!!      <pre type="block_jacobi">
!!        <nlocal>-1</nlocal>
!!        <solver type="ilu"></solver>
!!      </pre>
!!    </pre>
!!  </xmhd>
!!</oft>
!!```
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./MUG_sph_heat oft.in oft_in.xml
!!\endverbatim
!!
!!\section doc_mug_sph_ex2_post Post-Processing options
!!
!! During the simulation the evolution of some global quantities is written to the history file `xmhd.hist`. A utility
!! script (`plot_mug_hist.py` located in `bin` after installation) is included as part of OFT for plotting the signals in this file.
!! This is useful for monitoring general progress of the simulation as well as numerical parameters like the iteration
!! count and solver time. Using this script we can see the increase in the average temperature in time and a brief velocity transient
!! early in time as equilibrium adjusts slightly to be consistent with the magnetic BCs used in the simulation.
!!
!!\verbatim
!!~$ python plot_mug_hist.py xmhd.hist
!!\endverbatim
!!
!! \image html MUG_heat_ex-Ti.png "Time evolution of volume averaged temperature, showing heating due to Ohmic dissipation of current."
!!
!! \image html MUG_heat_ex-Kin.png "Time evolution of kinetic energy throughout simulation, showing an initial transient that damps away."
!!
!!\subsection doc_mug_sph_ex2_post_plot Creating plot files
!!
!! To generate 3D plots, and perform additional diagnostic sampling (see \ref xmhd::oft_xmhd_probe "oft_xmhd_probe"), a plot run can
!! be performed by setting `plot_run=T` in the `sph_tilt_options` input group, which calls \ref xmhd::xmhd_plot "xmhd_plot". With this option
!! additional run time options are available in the `xmhd_plot_options` group that control how restart files are sampled for plotting.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-5
!! rst_start=0
!! rst_end=2000
!!/
!!\endverbatim
!!
!! Once the post-processing run is complete `bin/build_xdmf.py` can be used to generate `*.xmf` files that can be loaded by
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!!
!! \image html MUG_heat_ex-Fields.png "Resulting temperature distribution (shading) and magnetic field (vectors) at the end of the simulation, showing thermal confinement in the core of the torus."
!!
!!\section doc_mug_sph_ex2_mesh Mesh Creation
!! A mesh file `cyl_heat.h5` is provided with this example. Instructions to generate your
!! own mesh for the geometry using [CUBIT](https://cubit.sandia.gov/) and [GMSH](https://gmsh.info/).
!!
!!\subsection doc_mug_sph_ex2_cubit Meshing with CUBIT
!!
!! A suitable mesh for this example, with radius of 1m and height of 1m, can be created using
!! the CUBIT script below.
!!
!!\verbatim
!!reset
!!
!!create Cylinder height 1 radius 1
!!
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .2
!!mesh volume 1
!!
!!set duplicate block elements off
!!block 1 add volume 1 
!!block 1 element type tetra10
!!
!!set large exodus file on
!!export Genesis  "cyl_heat.g" overwrite block 1
!!\endverbatim
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_cubit.py --in_file=cyl_heat.g
!!\endverbatim
!!
!!\subsection doc_mug_sph_ex2_gmsh Meshing with Gmsh
!!
!! If the CUBIT mesh generation codes is not avilable the mesh can be created using the Gmsh code and the
!! geometry script below.
!!
!!\verbatim
!!Coherence;
!!Point(1) = {0, 0, 0, 1.0};
!!Point(2) = {1, 0, 0, 1.0};
!!Point(3) = {0, 1, 0, 1.0};
!!Point(4) = {-1, 0, 0, 1.0};
!!Point(5) = {0, -1, 0, 1.0};
!!Circle(1) = {2, 1, 3};
!!Circle(2) = {3, 1, 4};
!!Circle(3) = {4, 1, 5};
!!Circle(4) = {5, 1, 2};
!!Line Loop(5) = {2, 3, 4, 1};
!!Plane Surface(6) = {5};
!!Extrude {0, 0, 1} {
!!  Surface{6};
!!}
!!\endverbatim
!!
!! To generate a mesh, with resolution matching the Cubit example above, place the script contents in a file called
!! `cyl_heat.geo` and run the following command.
!!
!!\verbatim
!!~$ gmsh -3 -format mesh -optimize -clscale .2 -order 2 -o cyl_heat.mesh cyl_heat.geo
!!\endverbatim
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_gmsh.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_gmsh.py --in_file=cyl_heat.mesh
!!\endverbatim