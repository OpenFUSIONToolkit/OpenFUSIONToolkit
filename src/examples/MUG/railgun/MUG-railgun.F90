!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Plasma railgun    {#doc_mug_railgun}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of \ref doc_mhd_main "MUG" to model a simple pulsed plasma accelerator, in
!! either a parallel-plate ("railgun") or a coaxial (Marshall gun) geometry. The channel is initially filled with a
!! uniform, unmagnetized plasma. A bias voltage is then applied across the electrodes at the closed "breech" end of
!! the channel (minimum \f$ x \f$), which drives a current through the plasma and a magnetic field behind it. The
!! resulting \f$ J \times B \f$ force forms a current sheet that sweeps the plasma ahead of it down the channel,
!! like a snowplow.
!!
!!\section doc_mug_railgun_phys Physics background
!!
!! Accelerators of this type were studied extensively for space propulsion and plasma injection, starting with the
!! coaxial plasma gun of [Marshall (1960)](https://doi.org/10.1063/1.1705989). When the channel is filled with gas
!! before the discharge, the current sheet acts like a magnetic piston, gathering up the gas in front of it
!! (the "snowplow" mode). When the magnetic Reynolds number is low, field instead diffuses through the plasma
!! and a nearly stationary high-current region forms near the breech (the "deflagration" mode).
!! [Sitaraman and Raja (2014)](https://doi.org/10.1063/1.4861867) used resistive MHD simulations of a
!! coaxial accelerator to show that resistivity controls which mode forms. This example operates in the
!! snowplow mode, except in a deliberately resistive region at the breech (see \ref doc_mug_railgun_code_shape).
!!
!! A simple estimate of the sheet speed follows from treating the sheet as a thin piston that collects all of
!! the plasma it passes. If the sheet moves at a constant speed \f$ v \f$ and the field behind it is steady,
!! then \f$ \nabla \times E = 0 \f$ behind the sheet and the applied field \f$ E = V/h \f$ is carried unchanged
!! to the sheet, where \f$ E = v B \f$. The magnetic pressure balances the momentum flux of the swept-up mass,
!! \f$ B^2 / 2 \mu_0 = \rho_0 v^2 \f$. Combining these gives the snowplow speed for a voltage-driven accelerator
!!
!! \f[ v_{sp} = \left( \frac{E^2}{2 \mu_0 \rho_0} \right)^{1/4}, \qquad B = \frac{E}{v_{sp}}. \f]
!!
!! For the inputs below (\f$ V = 500 \f$ V across a gap of \f$ h = 10 \f$ cm, deuterium at
!! \f$ n_0 = 5 \times 10^{19} m^{-3} \f$), \f$ v_{sp} \approx 88 \f$ km/s and \f$ B \approx 0.057 \f$ T. This is
!! about Mach 2.8 relative to the initial 6 eV plasma, so the sheet drives a shock ahead of it.
!!
!! The snowplow model assumes the swept-up plasma is compressed into an infinitely thin layer. At this Mach number the
!! shock only compresses the plasma by a factor of about 3.3, so the swept-up plasma forms a layer that grows
!! between the shock and the sheet. Treating the sheet as a magnetic piston moving at \f$ v \f$, the pressure behind the
!! shock must balance the magnetic pressure, \f$ B^2 / 2 \mu_0 = p_0 + \rho_0 u_s v \f$, where the shock speed
!! \f$ u_s = \frac{\gamma+1}{4} v + \sqrt{\left(\frac{\gamma+1}{4} v\right)^2 + c_s^2} \f$ follows from the
!! Rankine-Hugoniot relations and \f$ c_s \f$ is the initial sound speed. With \f$ E = v B \f$ this gives a piston
!! speed of 79 km/s, about 10% below the snowplow estimate, with the shock running ahead at 114 km/s.
!!
!! In the coaxial case the bias field varies as \f$ E = V/(r \ln(r_o/r_i)) \f$. Applying the estimate above at
!! each radius gives \f$ v \propto r^{-1/2} \f$ (105 km/s at the inner electrode and 75 km/s at the outer
!! electrode for these inputs), so the sheet tilts, with the inner electrode leading. Pressure along the sheet
!! partially evens out this difference. Treating the sheet as a single flat piston instead gives
!! \f$ v^4 = V^2 / (\mu_0 \rho_0 (r_o^2 - r_i^2) \ln(r_o/r_i)) \f$, which is \f$ \approx 87 \f$ km/s.
!!
!!\section doc_mug_railgun_code Code Walk Through
!!
!! The code consists of a helper module, which defines the bias source and dissipation shaping, and the main
!! program, which sets up the finite element spaces and runs the simulation.
!!
!!\subsection doc_mug_railgun_code_helper Helper Module
! START SOURCE
MODULE railgun_helpers
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE oft_la_base, ONLY: oft_vector
USE fem_utils, ONLY: fem_interp
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_cross_sint, oft_hcurl_grad_gzerop
USE xmhd, ONLY: oft_xmhd_driver, oft_xmhd_pop, oft_xmhd_push, xmhd_sub_fields, &
  xmhd_ML_lagrange, xmhd_ML_hcurl_grad
IMPLICIT NONE
#include "local.h"
!! To apply the bias electric field in the MHD simulation, we define a class that inherits from
!! \ref xmhd::oft_xmhd_driver "oft_xmhd_driver". This class is called at the start of each nonlinear
!! solve to modify the RHS and initial guess vectors as needed for source terms.
!!
!! Here we construct a `source` vector, which holds the surface contribution of the bias E-field to the
!! B-field time-advance equation, and add it to the RHS with the appropriate voltage on every time step.
TYPE, EXTENDS(oft_xmhd_driver) :: bias_esource
  REAL(r8) :: voltage = 1.d0                     ! Bias voltage [V]
  REAL(r8) :: ramp_time = 1.d-6                  ! Voltage ramp time [s]
  CLASS(oft_vector), POINTER :: source => NULL() ! Bias E-field source (for 1 V)
CONTAINS
  ! Apply BC(s) to RHS and initial guess vectors
  PROCEDURE :: apply => bias_esource_apply
END TYPE bias_esource
!! To compute the `source` vector for the driver class, we must compute
!! \f[ \int \phi \cdot (E \times dS), \f]
!! where \f$ \phi \f$ are the basis/test functions for the combined \f$ H(Curl) \f$ + \f$ \nabla H^1 \f$
!! FE space and \f$ E \f$ is a pre-specified vector field. A subroutine (\ref oft_hcurl_grad_operators::oft_hcurl_grad_cross_sint "oft_hcurl_grad_cross_sint")
!! is already available to compute the desired integral for an arbitrary field, so here we define
!! a \ref fem_utils::fem_interp "fem_interp" class that interpolates the desired E-field during integration.
!!
!! In addition to the usual `interp` method we also define a `setup` method, which is used to detect
!! parallel-plate vs coaxial configurations and compute geometry-specific constants.
TYPE, EXTENDS(fem_interp) :: bias_einterp
  LOGICAL :: coax = .FALSE. ! Coaxial geometry?
  REAL(r8) :: h = 0.d0      ! Electrode gap (parallel-plate) or \f$ \ln(r_o/r_i) \f$ (coaxial)
CONTAINS
  ! Setup geometry-specific values
  PROCEDURE :: setup => bias_einterp_setup
  ! Interpolate field
  PROCEDURE :: interp => bias_einterp_apply
END TYPE bias_einterp
!! We must also define a few global variables for use across subroutines/functions.
TYPE(bias_einterp) :: bias_field ! Instantiation of `bias_einterp`
REAL(r8) :: xleft = 0.d0         ! "Left-most" extent of domain (breech, x-direction)
REAL(r8) :: xright = 0.d0        ! "Right-most" extent of domain (muzzle, x-direction)
CONTAINS
!! When an instance of \ref xmhd::oft_xmhd_driver "oft_xmhd_driver" is passed to \ref xmhd::xmhd_run "xmhd_run"
!! the `apply` method is called before each nonlinear solve, allowing the user to apply many generic BCs
!! without having to modify the code of the MHD module itself. In this case we apply the bias electric field
!! by adding its contribution to the RHS vector, with the voltage ramped linearly from zero over `ramp_time`.
SUBROUTINE bias_esource_apply(self,up,u,v,t,dt)
CLASS(bias_esource), INTENT(inout) :: self ! Class object
CLASS(oft_vector), INTENT(inout) :: up     ! Solution vector at start of time step
CLASS(oft_vector), INTENT(inout) :: u      ! Guess vector for nonlinear solve
CLASS(oft_vector), INTENT(inout) :: v      ! RHS vector for nonlinear solve
REAL(r8), INTENT(in) :: t                  ! Time at start of time step
REAL(r8), INTENT(in) :: dt                 ! Time step size
!---Local variables
REAL(r8) :: volt_curr
TYPE(xmhd_sub_fields) :: sub_fields
TYPE(oft_hcurl_grad_gzerop) :: grad_zerop
!---------------------------------------------------------------------------
! Create voltage source if not available
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(self%source))THEN
  CALL xmhd_ML_hcurl_grad%vec_create(self%source)
  CALL bias_field%setup(xmhd_ML_lagrange%ml_mesh%mesh)
  CALL oft_hcurl_grad_cross_sint(xmhd_ML_hcurl_grad%current_level,bias_field,self%source)
  grad_zerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
  CALL grad_zerop%apply(self%source)
END IF
!---------------------------------------------------------------------------
! Update RHS field with source BC
!---------------------------------------------------------------------------
CALL xmhd_ML_hcurl_grad%vec_create(sub_fields%B)         ! Allocate the field we want to retrieve
CALL oft_xmhd_pop(v,sub_fields)                          ! Grab desired field from RHS vector
volt_curr = self%voltage*MIN(1.d0,(t+dt)/self%ramp_time) ! Linearly ramp voltage
CALL sub_fields%B%add(1.d0,dt*volt_curr,self%source)     ! Add source to B-field RHS
CALL oft_xmhd_push(sub_fields,v)                         ! Return modified field to RHS vector
!---------------------------------------------------------------------------
! Delete temporaries
!---------------------------------------------------------------------------
CALL sub_fields%B%delete
DEALLOCATE(sub_fields%B)
END SUBROUTINE bias_esource_apply
!! This subroutine detects the geometry and computes the extents of the channel. A coaxial channel is
!! detected from the mesh, which in that case is periodic in angle about the x-axis
!! (`mesh%periodic%revolved`). The x-extents mark the breech, where the bias field is applied, and are also
!! used to locate the high-dissipation entrance region.
SUBROUTINE bias_einterp_setup(self,mesh)
CLASS(bias_einterp), INTENT(inout) :: self  ! Class object
CLASS(oft_mesh), TARGET, INTENT(in) :: mesh ! Mesh used for interpolation
REAL(r8) :: ymin,ymax
self%mesh=>mesh
IF(mesh%periodic%revolved)THEN
  !---Get radial extents
  ymin=SQRT(MINVAL(self%mesh%r(2,:)**2+self%mesh%r(3,:)**2))
  ymax=SQRT(MAXVAL(self%mesh%r(2,:)**2+self%mesh%r(3,:)**2))
  ymin=-oft_mpi_max(-ymin)
  ymax=oft_mpi_max(ymax)
  self%h=LOG(ymax/ymin)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'  Coaxial setup:'
    WRITE(*,'(A,ES10.2)')'    Inner radius  = ',ymin
    WRITE(*,'(A,ES10.2)')'    Outer radius  = ',ymax
  END IF
  self%coax=.TRUE.
ELSE
  !---Get electrode positions
  ymin=MINVAL(self%mesh%r(2,:))
  ymax=MAXVAL(self%mesh%r(2,:))
  ymin=-oft_mpi_max(-ymin)
  ymax=oft_mpi_max(ymax)
  self%h=(ymax-ymin)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'  Parallel plate setup:'
    WRITE(*,'(A,ES10.2)')'    Channel width = ',self%h
  END IF
END IF
xleft=MINVAL(self%mesh%r(1,:))
xright=MAXVAL(self%mesh%r(1,:))
xleft=-oft_mpi_max(-xleft)
xright=oft_mpi_max(xright)
IF(oft_env%head_proc)WRITE(*,'(A,ES10.2)')'    Length        = ',xright-xleft
END SUBROUTINE bias_einterp_setup
!! This subroutine evaluates the bias E-field for a 1 V bias, which is uniform (\f$ E_y = 1/h \f$) in the
!! parallel-plate case and \f$ E_r = 1/(r \ln(r_o/r_i)) \f$ in the coaxial case. The field is only
!! non-zero on the breech surface, as determined by the x-location of the integration points.
SUBROUTINE bias_einterp_apply(self,cell,f,gop,val)
CLASS(bias_einterp), INTENT(inout) :: self ! Class object
INTEGER(i4), INTENT(in) :: cell            ! Cell containing evaluation point
REAL(r8), INTENT(in) :: f(:)               ! Logical position of point in cell
REAL(r8), INTENT(in) :: gop(3,4)           ! Cell spatial Jacobian matrix
REAL(r8), INTENT(out) :: val(:)            ! Interpolated field
REAL(r8) :: pt(3)
val(1:3)=0.d0
pt=self%mesh%log2phys(cell,f)
IF(pt(1)>xleft+1.d-5)RETURN
IF(self%coax)THEN
  !---1/r electric field in r-direction
  val(1:3)=[0.d0,pt(2),pt(3)]
  val(1:3)=val(1:3)/(SUM(val(1:3)**2)*self%h)
ELSE
  !---Constant electric field in y-direction
  val(2)=1.d0/self%h
END IF
END SUBROUTINE bias_einterp_apply
!!\subsection doc_mug_railgun_code_shape Entrance dissipation
!! The breech is a closed, no-slip wall, so the plasma accelerated away from it would leave a near-vacuum
!! behind the sheet. To keep the numerics well behaved, density is held at \f$ n_0 \f$ on the breech
!! (see \ref doc_mug_railgun_input) and the resistivity (\f$ \eta \f$) and density diffusion
!! (\f$ D \f$) are raised by up to a factor of 101 near the inlet. This refills the region behind the sheet
!! from the breech, and also gives a resistive, deflagration-like region in which the current sheet forms.
!! The sheet leaves this region about 0.3 m from the breech, after which it propagates as a snowplow.
!!
!! As with the driver, "shaping" functions are used to modify the dissipation coefficients without having
!! to change the core MHD code. The function returns the factor applied to the base coefficient.
FUNCTION entrance_diff_shape(pt) RESULT(val)
REAL(r8), INTENT(in) :: pt(3) ! Location in physical space
REAL(r8) :: val               ! Scaling factor for dissipation coefficient
val = (1.d0 + 100.d0*(1.d0-TANH(10.d0*(pt(1)-(xleft*0.8d0+xright*0.2d0))))/2.d0)
END FUNCTION entrance_diff_shape
END MODULE railgun_helpers
!!\subsection doc_mug_railgun_code_prog Main Program
!! The main program reads the runtime options, builds the mesh and finite element spaces, sets the initial
!! conditions and runs the simulation (or the post-processing). Runtime options are set in the main input
!! file using the `railgun_options` group.
!!
!! **Option group:** `railgun_options`
!! |  Option         |  Description  | Type [dim] |
!! |-----------------|---------------|------------|
!! | `order=2`       | FE order | int |
!! | `minlev=-1`     | Lowest level for multigrid solvers (only polynomial MG is supported) | int |
!! | `plot_run=F`    | Perform plotting run instead of simulation | bool |
!! | `pm=F`          | Show linear solver iterations | bool |
!! | `n0=1E19`       | Initial density [m^-3] | real |
!! | `t0=5`          | Initial temperature [eV] | real |
!! | `voltage=1E4`   | Bias voltage [V] | real |
!! | `ramp_time=1E-6` | Time to linearly ramp the bias voltage [s] | real |
PROGRAM MUG_railgun
USE oft_base
!--- Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!--- Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp
!--- H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!--- Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp
!--- Physics
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, vel_scale, temp_floor, den_scale, &
  den_floor, xmhd_sub_fields, xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, &
  xmhd_ML_lagrange, xmhd_ML_vlagrange, res_shaping, diff_shaping
!--- Helper module
USE railgun_helpers, ONLY: bias_esource, entrance_diff_shape
IMPLICIT NONE
!--- Local variables
INTEGER(i4) :: ierr,io_unit
TYPE(multigrid_mesh) :: mg_mesh
TYPE(bias_esource) :: driver
TYPE(xmhd_sub_fields) :: ic_fields
!--- Input variables
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = -1
LOGICAL :: plot_run = .FALSE.
LOGICAL :: pm = .FALSE.
REAL(r8) :: n0 = 1.d19
REAL(r8) :: t0 = 5.d0
REAL(r8) :: voltage = 1.d4
REAL(r8) :: ramp_time = 1.d-6
NAMELIST/railgun_options/order,minlev,plot_run,n0,t0,voltage,ramp_time,pm
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!--- Read options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,railgun_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No "railgun_options" found in input file.','MUG_railgun',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "railgun_options" in input file.','MUG_railgun',__FILE__)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh)
IF((minlev<mg_mesh%mgdim).AND.(minlev/=-1))CALL oft_abort("Only polynomial MG currently supported for xMHD", &
  "MUG_railgun",__FILE__)
!---------------------------------------------------------------------------
! Build FE spaces
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
!---------------------------------------------------------------------------
! Run simulation or plotting
!---------------------------------------------------------------------------
xmhd_minlev=minlev       ! Lowest level for MG solvers
vel_scale=1.d3           ! Velocity equation scale factor
den_scale=n0             ! Density equation scale factor
den_floor=n0/5.d0        ! Set density floor
temp_floor=t0/3.d0       ! Set temperature floor
oft_env%pm=pm            ! Show linear iteration progress?
IF(plot_run)THEN
  !--- Run plotting routine
  CALL xmhd_plot
ELSE
  !--- Create fields for initial conditions
  CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
  CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
  CALL xmhd_ML_lagrange%vec_create(ic_fields%Ne)
  CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
  !--- Set initial conditions (uniform plasma at rest with no magnetic field)
  CALL ic_fields%Ne%set(n0)
  CALL ic_fields%Ti%set(t0)
  !--- Set custom dissipation shaping at entrance
  res_shaping=>entrance_diff_shape
  diff_shaping=>entrance_diff_shape
  !--- Run simulation
  driver%voltage=voltage
  driver%ramp_time=ramp_time
  CALL xmhd_run(ic_fields,driver=driver)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM MUG_railgun
! STOP SOURCE
!!
!!\section doc_mug_railgun_input Input files
!!
!! Below is an input file for the parallel-plate case, which can be used with this example in a parallel
!! environment. This example should not be run with only a single process as solving the time-dependent MHD
!! equations is significantly more challenging than previous examples. For more information on the options
!! in the `xmhd_options` group see \ref xmhd::xmhd_run "xmhd_run". The coaxial input file (`oft_coax.in`) is
!! the same apart from the mesh file and the number of partitions.
!!
!! The plasma model is single-fluid resistive MHD (`xmhd_hall=F`) with density and temperature advance.
!! The resistivity is close to the Spitzer value at 6 eV. Ohmic and viscous heating are disabled, so the plasma is
!! heated only by compression. The run covers 10.5 \f$ \mu s \f$, which ends shortly before the shock ahead of the
!! sheet reaches the end wall of the channel. With the tetrahedral meshes provided, each case takes roughly 15
!! (parallel plate) to 25 (coaxial) minutes on 3 MPI processes. The hexahedral meshes describe the same geometry
!! and also run with these inputs (set `filename` in `native_mesh_options`), but with the preconditioner below
!! each time step is about 5-15 times more expensive.
!!
! OFT_DOC_INCLUDE: oft_slab.in
!!
!!\subsection doc_mug_railgun_input_solver Solver and boundary condition specification
!!
!! The XML input file sets the per-surface boundary conditions and the preconditioner. The magnetic field
!! uses `bbc='bc'`, where all walls are perfect conductors (\f$ B \cdot \hat{n} = 0 \f$, tangential
!! \f$ E = 0 \f$) except where the bias source adds a tangential E-field on the breech. The velocity, density
!! and temperature BCs are set per nodeset (`vbc='xml'`, `nbc='x'`, `tbc='x'`), where each entry in
!! `<den_bc>`, `<temp_bc>` and `<vel_bc>` corresponds to the nodeset with that index in the mesh file.
!! For density and temperature, `1` fixes the value at the initial value and `0` gives a natural (zero-flux)
!! condition. For velocity, `1` fixes all components (no-slip), `2` fixes only the normal component (free-slip) and
!! `3` fixes only the tangential components.
!!
!! | Nodeset | Surface | Density | Temperature | Velocity |
!! |---------|---------|---------|-------------|----------|
!! | 1 | Breech (\f$ x = -0.5 \f$) | Fixed at \f$ n_0 \f$ | Fixed at \f$ T_0 \f$ | No-slip |
!! | 2 | Muzzle end wall (\f$ x = 0.5 \f$) | Zero flux | Zero flux | No-slip |
!! | 3 | Electrodes | Zero flux | Fixed at \f$ T_0 \f$ (cold wall) | Free-slip |
!!
!! The remaining boundaries (z-faces of the slab and angular faces of the coaxial wedge) are periodic.
!!
!! The preconditioner for the parallel-plate case is block Jacobi with an LU factorization on each MPI
!! process. The coaxial case (`oft_coax.xml`) wraps the same preconditioner in a short GMRES iteration.
!!
! OFT_DOC_INCLUDE: oft_slab.xml
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ mpirun -np 3 ./MUG-railgun oft_slab.in oft_slab.xml
!!~$ mpirun -np 3 ./MUG-railgun oft_coax.in oft_coax.xml
!!\endverbatim
!!
!!\section doc_mug_railgun_post Post-Processing options
!!
!!\subsection doc_mug_railgun_post_plot Creating plot files
!!
!! To generate 3D plots, and perform additional diagnostic sampling (see \ref xmhd::oft_xmhd_probe "oft_xmhd_probe"), a plot run can
!! be performed by setting `plot_run=T` in the `railgun_options` input group, which calls \ref xmhd::xmhd_plot "xmhd_plot". With this option
!! additional run time options are available in the `xmhd_plot_options` group that control how restart files are sampled for plotting.
!!
!! Once the post-processing run is complete `bin/OFT_build_xdmf.py` can be used to generate `*.xmf` files that can be loaded by
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!!
!!\subsection doc_mug_railgun_post_sheet Current sheet speed
!!
!! The script `plot_railgun.py`, included with this example, reads the plot files from both cases, tracks the
!! peak of the current density along the channel at several transverse locations and fits the sheet speed while
!! the sheet is between \f$ x = -0.2 \f$ and \f$ 0.35 \f$ m (after it leaves the resistive entrance region). It
!! prints the measured speeds alongside the piston and snowplow estimates from \ref doc_mug_railgun_phys and
!! creates the figures below. Run it from the directory above the two run directories (or pass their paths) after
!! a plotting run in each.
!!
!!\verbatim
!!~$ python plot_railgun.py --slab slab_run --coax coax_run --time 8.E-6
!!Sheet speed [km/s]
!!Location                   Measured     Piston   Snowplow
!!Coax, r = 0.11 m               88.3       91.1      100.6
!!Coax, r = 0.15 m               89.0       77.3       86.1
!!Coax, r = 0.18 m               88.4       70.1       78.6
!!Parallel plate                 80.1       78.9       87.8
!!Coax, flat sheet: snowplow 87.0 km/s
!!\endverbatim
!!
!! In the parallel-plate case the sheet forms in the resistive entrance region and then moves down the channel at a
!! constant 80 km/s. This is within 2% of the piston estimate and 9% below the snowplow estimate. Behind the sheet, the
!! magnetic field is about \f$ 1.1 B_{sp} \f$ and the plasma moves at about \f$ 0.9 v_{sp} \f$, consistent with
!! \f$ E = v B \f$ at the measured speed. The swept-up plasma forms a dense layer between the sheet and the shock. The peak compression (about 4.7) is larger than the adiabatic
!! shock value of 3.3, likely because thermal conduction to the fixed-temperature electrodes cools the layer.
!!
!! In the coaxial case the sheet tilts, with the inner electrode leading, as expected from the stronger field there.
!! However, the tilt does not keep growing as the local (\f$ v \propto r^{-1/2} \f$) estimate would suggest. Instead the
!! sheet moves as a single structure at about 88 km/s, close to the piston estimate at the inner electrode. Near
!! the outer electrode, where the drive is weakest, the current spreads out and part of it trails behind along
!! the electrode.
!!
!! \image html MUG_railgun_ex-Sheet.png "Left: Current sheet position against time at several transverse locations, with fitted speeds and piston estimates (in parentheses, dashed lines). The shaded band is the range used for the fits. Right: Magnetic field, density and axial velocity along the center of the parallel-plate channel at 8.1 us, normalized by the snowplow estimates and initial density."
!!
!! \image html MUG_railgun_ex-J.png "Transverse current density (radial for coaxial) in the coaxial (top) and parallel-plate (bottom) cases at 8.1 us"
!!
!! \image html MUG_railgun_ex-V.png "Axial flow velocity in the coaxial (top) and parallel-plate (bottom) cases at 8.1 us"
!!
!!\section doc_mug_railgun_mesh Mesh Generation
!! Mesh files for the parallel-plate (`MUG-railgun-tet.h5`, `MUG-railgun-hex.h5`) and coaxial
!! (`MUG-coax_railgun-tet.h5`, `MUG-coax_railgun-hex.h5`) cases are provided with this example, using
!! tetrahedral and hexahedral elements. Instructions to generate your own mesh for the geometry using
!! [CUBIT](https://cubit.sandia.gov/) are also provided below. Only a thin slice of each geometry is
!! meshed, which is then reflected across the \f$ z = 0 \f$ plane and made periodic in \f$ z \f$ (slab) or
!! angle (coaxial) when loaded (`reflect=T` and `ref_periodic=T` in `native_mesh_options`). The coaxial
!! geometry is detected automatically from the angular periodicity.
!!
!!\subsection doc_mug_railgun_mesh_par Parallel plate/Slab
!! A suitable mesh for this example in slab geometry, with an electrode gap of 10 cm and a length of 1 m,
!! can be created using the CUBIT script below.
!!
! OFT_DOC_INCLUDE: MUG-railgun.jou
!!
!!\subsection doc_mug_railgun_mesh_coax Coaxial/Cylindrical
!! A suitable mesh for this example in cylindrical geometry, with an inner radius of 10 cm, an outer radius of
!! 20 cm and a length of 1 m, can be created using the CUBIT script below. The 5 degree wedge becomes a 10 degree
!! periodic wedge after reflection.
!!
! OFT_DOC_INCLUDE: MUG-coax_railgun.jou
!!
!! Once complete the mesh should be converted into the native mesh format using the `OFT_convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python OFT_convert_cubit.py --in_file=MUG-railgun-tet.g --periodic_nodeset=-1
!!\endverbatim
!!
!! @note The option `--periodic_nodeset=-1` indicates that the surface corresponding to the last nodeset (python-style
!! reverse indexing) should be treated as periodic following reflection.
