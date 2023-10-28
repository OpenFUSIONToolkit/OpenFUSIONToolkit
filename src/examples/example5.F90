!!MUG Example 1: Spheromak Tilt    {#doc_mhd_ex1}
!!=========================
!!
!![TOC]
!!
!!This example demonstrates the use of the \ref xmhd "extended MHD" module in the Open FUSION Toolkit (OFT). It is well known
!!that a spheromak confined in a cylindrical flux conserver is unstable to an ideal tilt
!!mode when the ratio of the cylinder's height to radius is greater than 1.3. This result
!!was shown first by Bondeson and Marklin, where the instablity acts to dissipate energy
!!and drive the magnetic configuration toward the Taylor state.
!!
!!\section doc_ex5_cubit Mesh Creation with CUBIT
!!
!!A suitable mesh for this example, with radius of 1m and height of 2m, can be created using
!!the CUBIT script below.
!!
!!\verbatim
!!reset
!!create Cylinder height 2 radius 1
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .2
!!mesh volume 1
!!refine parallel fileroot 'cyl' overwrite no_execute
!!\endverbatim
!!
!!\section doc_ex5_code Code Walk Through
!!
!!The code consists of three basic sections, required imports and variable definitions,
!!finite element setup, and system creation and solution.
!!
!!\subsection doc_ex5_code_inc Module Includes
! START SOURCE
PROGRAM xmhd_cyl
!---Runtime
USE oft_base
!---Grid
USE oft_mesh_type, ONLY: mesh, rgrnd
USE multigrid_build, ONLY: multigrid_construct, multigrid_add_quad
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels, oft_lag_set_level
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_level, oft_hcurl_nlevels
USE oft_hcurl_operators, ONLY: hcurl_setup_interp, hcurl_mloptions, hcurl_zerob
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_divout, h1_zeroi, h1_mc, h1curl_zerob, h1_setup_interp
!---Physics
USE taylor, ONLY: taylor_hmodes, taylor_minlev, taylor_hffa, taylor_hlam
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_taxis, xmhd_sub_fields
IMPLICIT NONE
!!\section doc_ex5_code_vars Local Variables
!---H1 divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
CLASS(oft_matrix), pointer :: lop => NULL()
!---Local variables
INTEGER(i4) :: ierr,io_unit
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
CLASS(oft_vector), POINTER :: db => NULL()
TYPE(xmhd_sub_fields) :: ic_fields
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: b0_scale = 1.E-1_r8
REAL(r8) :: b1_scale = 1.E-5_r8
REAL(r8) :: n0 = 1.d19
REAL(r8) :: t0 = 6.d0
LOGICAL :: plot_run=.FALSE.
NAMELIST/cyl_options/order,minlev,b0_scale,b1_scale,plot_run,n0,t0
!!\section doc_ex5_code_setup OFT Initialization
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,cyl_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
rgrnd=(/2.d0,0.d0,0.d0/)
CALL multigrid_construct
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!---Lagrange
CALL oft_lag_setup(order, minlev)
CALL lag_setup_interp
CALL lag_mloptions
!---H1(Curl) subspace
CALL oft_hcurl_setup(order, minlev)
CALL hcurl_setup_interp
CALL hcurl_mloptions
!---H1(Grad) subspace
CALL oft_h0_setup(order+1, minlev)
CALL h0_setup_interp
!---H1 full space
CALL oft_h1_setup(order, minlev)
CALL h1_setup_interp
!!\section doc_ex5_code_taylor Computing Initial Conditions
!!
!! The spheromak mode corresponds to the thrid lowest force-free eignstate of the
!! flux conserver. This can be computed using the \ref taylor::taylor_hmodes
!! "taylor_hmodes" subroutine and with the desired number of modes set to three.
!! The fact that this mode is the third
!! eigenstate is indicative of the ideal instability we wish to simulate with this
!! example as there are two magentic configurations that contain less energy
!! for a given helicity. As a result, we expect instability to drive the configuration
!! toward one of these minimum energy states. This also leads us to our choice of an
!! intial perturbation to the equilibrium, which we will chose to match field of
!! the lowest eigenstate.
taylor_minlev=minlev
CALL taylor_hmodes(3)
CALL oft_lag_set_level(oft_lagrange_nlevels)
!!\subsection doc_ex5_code_taylor_gauge Setting Magnetic Boundary Conditions
!!
!! The \ref taylor::taylor_hmodes "taylor_hmodes" subroutine computes the vector
!! potential for each of the requested eignestates. However, the reduced MHD
!! solver uses magnetic field as the primary variable. With force-free eigenstate
!! solutions the vector potential and magentic field can be related by a gauge
!! transformation. Here this fact is used to produce an appropriate initial
!! condition from the vector potential without changing the representation
!! order or projecting the field.
!!
!! The gauge is transformed using a \ref oft_h1_operators::h1_divout "Divout"
!! procedure to remove divergence from the vector potential by adding a gradient.
!! For a gauge transformation all of the "divergence" is located at the boundary
!! and the \ref oft_h1_operators::h1_divout "Divout" solver only acts to change
!! the boundary condition of the field. This process is performed on both the
!! equilibrium and perturbing fields.
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop)
CALL oft_h0_getlop(lop,"grnd")
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
CALL create_diag_pre(linv%pre) ! Setup Preconditioner
divout%solver=>linv
divout%bc=>h0_zerogrnd
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
CALL oft_h1_create(ic_fields%B)
CALL taylor_hffa(1,oft_hcurl_level)%f%get_local(tmp)
CALL ic_fields%B%restore_local(tmp,1)
CALL divout%apply(ic_fields%B)
!---
CALL oft_h1_create(db)
CALL taylor_hffa(3,oft_hcurl_level)%f%get_local(tmp)
CALL db%restore_local(tmp,1)
CALL divout%apply(db)
!!\subsection doc_ex5_code_taylor_combine Set Initial Conditions
!!
!! With the required fields computed the full field can be initialized. The
!! full field is simply the equilibrium (3rd eigenmode) field plus a low
!! amplitude perturbation field (1st eigenmode). The velocity field is also
!! created but initialized to zero everywhere.
CALL ic_fields%B%scale(b1_scale*taylor_hlam(1,oft_hcurl_level))
CALL ic_fields%B%add(1.d0,b0_scale*taylor_hlam(3,oft_hcurl_level),db)
!---Clean up temporary matrices and fields
CALL lop%delete
CALL db%delete
DEALLOCATE(tmp,lop,db)
!!\section doc_ex5_code_run Run Simulation
!!
!! Finally, the simulation can be run using the driver routine for non-linear
!! MHD (\ref xmhd::xmhd_run "xmhd_run"). This routine advances the
!! solution in time with the physics specified in the input file, see the
!! documentation for \ref xmhd::xmhd_run "xmhd_run", and produces restart files
!! that contain the solution at different times.
!!
!! By default a MG preconditioner is used with the coarsest level specified by
!! \ref xmhd::xmhd_minlev "xmhd_minlev"  similar to the Taylor module. Several quantities
!! are also recorded to a history file \c "xmhd.hist" during the simulation, including the
!! toroidal current (where the symmetry axis is specified by \ref xmhd::xmhd_taxis
!! "xmhd_taxis"). The data in the history file may be plotted using the script
!! \c "src/utilities/scripts/plot_xmhd_hist.py"
!!
!! \note OFT plotting scripts require the python packages NUMPY and MATPLOTLIB as well
!! as path access to the python modules provided in "src/utilities".
!!
!! To visualize the solution fields once a simulation has completed the \ref xmhd::xmhd_plot
!! "xmhd_plot" subroutine is used. This subroutine steps through the restart files
!! produced by \ref xmhd::xmhd_run "xmhd_run" and generates plot files, and optionally
!! probe signals at evenly spaced points in time as specified in the input file, see
!! \ref xmhd::xmhd_plot "xmhd_plot".
xmhd_minlev=minlev
xmhd_taxis=3
oft_env%pm=.FALSE.
IF(plot_run)THEN
  !---Setup I/0
  CALL mesh%setup_io(order)
  !---Run post-processing routine
  CALL xmhd_plot
ELSE
  !---Create velocity field
  CALL oft_lag_vcreate(ic_fields%V)
  !---Create static density/temperature
  CALL oft_lag_create(ic_fields%Ne)
  CALL oft_lag_create(ic_fields%Ti)
  CALL ic_fields%Ne%set(n0)
  CALL ic_fields%Ti%set(t0)
  !---Run simulation
  CALL xmhd_run(ic_fields)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM xmhd_cyl
! STOP SOURCE
!!
!!\section doc_ex5_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment. This example
!! should not be run with only a single process as solving the time-depedent MHD equations is
!! significantly more challenging than previous examples. For more information on the options in the
!! \c xmhd_options group see \ref xmhd::xmhd_plot "xmhd_plot".
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='test'
!! cad_type=2
!! nlevels=2
!! nbase=1
!! grid_order=2
!! fix_boundary=true
!!/
!!
!!&cubit_options
!! filename='cyl.in.e'
!! inpname='cyl.3dm'
!! lf_file=T
!!/
!!
!!&hcurl_op_options
!! df_wop=0.,.65,.372,.324
!! nu_wop=0,64,2,1
!!/
!!
!!&lag_op_options
!! df_lop=0.,.9,0.86,0.64
!! nu_lop=0,64,2,1
!!/
!!
!!&cyl_options
!! order=3
!! minlev=2
!! b0_scale=1.e-1
!! b1_scale=1.e-5
!! n0=1.E19
!! t0=3.
!! plot_run=F
!!/
!!
!!&xmhd_options
!! bbc='bc'          ! Perfectly-conducting BC for B-field
!! vbc='all'         ! Zero-flow BC for velocity
!! dt=6.e-8          ! Maximum time step
!! eta=1.            ! Resistivity
!! nu_par=10.        ! Fluid viscosity
!! nsteps=2000       ! Number of time steps to take
!! rst_freq=10       ! Restart file frequency
!! lin_tol=1.E-10    ! Linear solver tolerance
!! nl_tol=1.E-5      ! Non-linear solver tolerance
!! nu_xmhd=0,1,8,2   ! Number of smoother iterations for default preconditioner
!! rst_ind=0         ! Index of file to restart from (0 -> use subroutine arguments)
!! ittarget=40       ! Target for # of linear iterations per time step
!! mu_ion=2.         ! Ion mass (atomic units)
!! xmhd_prefreq=20   ! Preconditioner update frequency
!!/
!!\endverbatim
!!
!!\subsection doc_ex5_input_plot Post-Processing options
!!
!! When running the code for post-processing additional run time options are available.
!!
!!\verbatim
!!&xmhd_plot_options
!! t0=1.E-8
!! dt=1.E-6
!! rst_start=0
!! rst_end=2000
!!/
!!\endverbatim
!!
!!\subsection doc_ex5_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate MG preconditioner.
!! Currently, this preconditioner method is the suggest preconditioner for all
!! time-dependent MHD solves.
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./example5 oft.in oft_in.xml
!!\endverbatim
!!
!! \warning Use of this preconditioner requires OFT be built with the PETSc and
!! FOX libraries.
!!
!!\verbatim
!!<oft>
!!  <xmhd>
!!    <pre type="mg">
!!      <smoother direction="both">
!!        <solver type="gmres">
!!          <its>0,0,2,2</its>
!!          <nrits>0,0,2,2</nrits>
!!          <pre type="block_jacobi">
!!            <nlocal>0,0,1,1</nlocal>
!!            <solver type="lu">
!!              <type>lu</type>
!!              <package>superd</package>
!!            </solver>
!!          </pre>
!!        </solver>
!!      </smoother>
!!      <coarse>
!!        <solver type="gmres">
!!          <its>8</its>
!!          <nrits>8</nrits>
!!          <pre type="block_jacobi">
!!            <solver type="lu">
!!              <type>lu</type>
!!              <package>superd</package>
!!            </solver>
!!          </pre>
!!        </solver>
!!      </coarse>
!!    </pre>
!!  </xmhd>
!!</oft>
!!\endverbatim
