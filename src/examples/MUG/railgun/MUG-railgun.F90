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
!! This example demonstrates the use of \ref doc_mhd_main "MUG" to model a simple plasma
!! railgun in either slab (parallel plate) or coaxial geometry. The system is initialized
!! with a uniform plasma throughout the domain and zero initial B-field. Then a suitable
!! E-field is applied on the "left-most" (x-direction) boundary to model biasing of the
!! top-bottom/inner-outer electrodes, driving a corresponding current and accelerating
!! the plasma across the simulation domain.
!!
!!\section doc_mug_railgun_code Code Walk Through
!!
!! The code consists of three basic sections, required imports and variable definitions,
!! finite element setup, and system creation and solution.
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
!! To apply the driving electric field in the MHD simulation we use define a
!! class that inherits from \ref xmhd::oft_xmhd_driver "oft_xmhd_driver". This class
!! is designed to be called at the start of the nonlinear solve to modify RHS and
!! initial guess vectors as needed for source terms.
!!
!! In this example we will construct a `source` vector, which corresponds to the
!! E-field surface contribution to the B-field time-advance equations.
TYPE, EXTENDS(oft_xmhd_driver) :: bias_esource
  REAL(8) :: voltage = 1.d0                      ! Parallel plate voltage
  REAL(r8) :: ramp_time = 1.E-6                  ! Voltage ramp time
  CLASS(oft_vector), POINTER :: source => NULL() ! Toroidal current source
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
!! In addition to the usual `interp` method we also define a `setup` method, which will be
!! used to determine geometry-specific constants and detect parallel plate vs coaxial configurations.
TYPE, EXTENDS(fem_interp) :: bias_einterp
  LOGICAL :: coax = .FALSE. ! Coaxial geometry
  REAL(8) :: xmin = 0.d0    ! Minimum x-coordinate
  REAL(8) :: h = 0.d0       ! Accelerator height (y-direction)
CONTAINS
  ! Setup geometry-specific values
  PROCEDURE :: setup => bias_einterp_setup
  ! Interpolate field
  PROCEDURE :: interp => bias_einterp_apply
END TYPE bias_einterp
!! We must also define a few global variables for use across subroutines/functions.
TYPE(bias_einterp) :: bias_field ! Instantiation of `bias_einterp`
REAL(r8) :: xleft = 0.d0         ! "Left-most" extent of domain (x-direction)
REAL(r8) :: xright = 0.d0        ! "Right-most" extent of domain (x-direction)
CONTAINS
!! When an instance of \ref xmhd::oft_xmhd_driver "oft_xmhd_driver" is passed to \ref xmhd::xmhd_run "xmhd_run"
!! the `apply` method is called before each nonlinear solve, allowing the user to apply many generic BCs
!! without having to modify the code of the MHD module itself. In this case we can apply the desired BC
!! for the bias electric field by adding its contribution to the RHS vector.
SUBROUTINE bias_esource_apply(self,up,u,v,t,dt)
CLASS(bias_esource), INTENT(inout) :: self ! Class object 
CLASS(oft_vector), INTENT(inout) :: up     ! Solution vector at start of time step
CLASS(oft_vector), INTENT(inout) :: u      ! Guess vector for nonlinear solve
CLASS(oft_vector), INTENT(inout) :: v      ! RHS vector for nonlinear solve
REAL(8), INTENT(in) :: t
REAL(8), INTENT(in) :: dt
!---Local variables
INTEGER(4) :: io_unit
REAL(8) :: volt_curr
TYPE(xmhd_sub_fields) :: sub_fields
TYPE(oft_hcurl_grad_gzerop) :: grad_zerop
type(oft_timer) :: mytimer
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
!! This subroutine computes the size and orientation (parallel or coaxial) of the
!! channel for the current grid and uses that to define the extents in the x-direction
!! as well as the characteristic length scale for the E-field interpolation.
subroutine bias_einterp_setup(self,mesh)
class(bias_einterp), intent(inout) :: self  ! Class object 
class(oft_mesh), target, intent(in) :: mesh ! Mesh used for interpolation
REAL(8) :: ymin,ymax
self%mesh=>mesh
IF(mesh%periodic%revolved)THEN
  !---Get mesh extents
  ymin=SQRT(MINVAL(self%mesh%r(2,:)**2+self%mesh%r(3,:)**2))
  ymax=SQRT(MAXVAL(self%mesh%r(2,:)**2+self%mesh%r(3,:)**2))
  ymin=oft_mpi_max(ymin)
  ymax=oft_mpi_max(ymax)
  self%h=LOG(ymax)-LOG(ymin)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'  Coaxial setup:'
    WRITE(*,'(A,3ES10.2)')'    Inner radius  = ',ymin
    WRITE(*,'(A,3ES10.2)')'    Outer radius  = ',ymax
  END IF
  self%coax=.TRUE.
ELSE
  !---Get mesh extents
  ymin=MINVAL(self%mesh%r(2,:))
  ymax=MAXVAL(self%mesh%r(2,:))
  ymin=-oft_mpi_max(-ymin)
  ymax=oft_mpi_max(ymax)
  self%h=(ymax-ymin)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(A)')'  Parallel plate setup:'
    WRITE(*,'(A,3ES10.2)')'    Channel width = ',self%h
  END IF
END IF
xleft=MINVAL(self%mesh%r(1,:))
xright=MAXVAL(self%mesh%r(1,:))
xleft=-oft_mpi_max(-xleft)
xright=oft_mpi_max(xright)
IF(oft_env%head_proc)WRITE(*,'(A,3ES10.2)')'    Length        = ',xright-xleft
end subroutine bias_einterp_setup
!! This subroutine evaluates the bias E-field for the given case, where
!! the field is constant and has a \f$ 1/r \f$ dependence for the parallel and
!! coaxial cases respectively. The field is only non-zero on the left-most surface
!! as determined by the x-location of the integration points.
subroutine bias_einterp_apply(self,cell,f,gop,val)
class(bias_einterp), intent(inout) :: self ! Class object 
integer(4), intent(in) :: cell             ! Cell containing evaluation point
real(8), intent(in) :: f(:)                ! Logical position of point in cell
real(8), intent(in) :: gop(3,4)            ! Cell spatial Jacobian matrix
real(8), intent(out) :: val(:)             ! Interpolated field
real(8) :: pt(3)
val(1:3)=0.d0
pt=self%mesh%log2phys(cell,f)
IF(self%coax)THEN
  !---1/r electric field in r-direction on MIN(x) boundary only
  IF(pt(1)<xleft+1.d-5)THEN
    val=[0.d0,pt(2),pt(3)]; val=val/(magnitude(val)**2)
    val=val/self%h
  END IF
ELSE
  !---Constant electric field in y-direction on MIN(x) boundary only
  IF(pt(1)<xleft+1.d-5)THEN
    val(2)=1.d0/self%h
  END IF
END IF
end subroutine bias_einterp_apply
!! For this example case we will fairly simple BCs at the inlet and outlet
!! of the channel. To ease the numerics somewhat, we will increase the resitivity (\f$ \eta \f$)
!! and diffusion (\f$ D \f$) near the inlet. This helps avoid very low density
!! at the inlet that would otherwise be caused by acceleration of the plasma without
!! an inflow or density source.
!!
!! As with the driver, "shaping" functions are defined to enable modifying dissipation
!! coefficients without having to change the core MHD code.
FUNCTION entrance_diff_shape(pt) result(val)
REAL(r8), INTENT(in) :: pt(3) ! Location in physical space
REAL(r8) :: val               ! Scaling factor for dissipation coefficient
val = (1.d0 + 100.d0*(1.d0-TANH(10.d0*(pt(1)-(xleft*0.8d0+xright*0.2d0))))/2.d0)
END FUNCTION entrance_diff_shape
END MODULE railgun_helpers
!!\subsection doc_mug_railgun_code_prog Main Program
!! Driver program for MHD simulation of parallel plate accelerator
!!
!! Runtime options are set in the main input file using the group
!! `railgun_options` group (see below).
PROGRAM mug_railgun
USE oft_base
!--- Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!--- Linear algebra
USE oft_la_base, ONLY: oft_vector
!--- Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp
!--- H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
!--- H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!--- Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp
!--- Physics
USE xmhd, ONLY: xmhd_run, xmhd_plot, xmhd_minlev, xmhd_bnorm_force, &
  vel_scale, temp_floor, den_scale, den_floor, xmhd_sub_fields, &
  xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, &
  xmhd_ML_lagrange, xmhd_ML_vlagrange, res_shaping, diff_shaping
!--- Helper module
USE railgun_helpers, ONLY: bias_esource, entrance_diff_shape
IMPLICIT NONE
!--- Local variables
INTEGER(4) :: ierr,io_unit
TYPE(multigrid_mesh) :: mg_mesh
TYPE(bias_esource) :: driver
TYPE(xmhd_sub_fields) :: ic_fields
!--- Input variables
INTEGER(4) :: order = 2
INTEGER(4) :: minlev = -1
LOGICAL :: plot_run=.FALSE.
LOGICAL :: pm = .FALSE.
REAL(8) :: n0=1.d19
REAL(8) :: T0=5.d0
REAL(8) :: voltage=10.d3
NAMELIST/railgun_options/order,minlev,plot_run,n0,T0,voltage,pm
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!--- Read options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,railgun_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh)
IF((minlev<mg_mesh%mgdim).AND.(minlev/=-1))CALL oft_abort("Only polynomial MG currently supported for xMHD","main",__FILE__)
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
! xmhd_bnorm_force=.FALSE. !
xmhd_minlev=minlev       ! Lowest level for MG solvers
vel_scale=1.d3           ! Velocity equation scale factor
den_scale=n0             ! Density equation scale factor
den_floor=n0/5.d0        ! Set temperature floor
temp_floor=T0/3.d0       ! Set density floor
oft_env%pm=pm            ! Show linear iteration progress?
IF(plot_run)THEN
  !--- Run plotting routine
  CALL xmhd_plot
ELSE
  !--- Create B fields
  CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
  CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
  CALL xmhd_ML_lagrange%vec_create(ic_fields%Ne)
  CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
  !--- Set initial conditions
  CALL ic_fields%Ne%set(n0)
  CALL ic_fields%Ti%set(T0)
  !--- Set custom diffusivity shaping at entrance
  res_shaping=>entrance_diff_shape
  diff_shaping=>entrance_diff_shape
  !--- Run simulation
  driver%voltage=voltage ! Set voltage
  CALL xmhd_run(ic_fields,driver=driver)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM mug_railgun
! STOP SOURCE
!!
!!\section doc_mug_railgun_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment. This example
!! should not be run with only a single process as solving the time-depedent MHD equations is
!! significantly more challenging than previous examples. For more information on the options in the
!! `xmhd_options` group see \ref xmhd::xmhd_run "xmhd_run".
!!
! OFT_DOC_INCLUDE: oft_slab.in
!!
!!\subsection doc_mug_railgun_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_slab.xml` file that constructs an appropriate ILU(0) preconditioner.
!! Currently, this preconditioner method is the suggested starting preconditioner for all
!! time-dependent MHD solves.
!!
! OFT_DOC_INCLUDE: oft_slab.xml
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./MUG-railgun oft_slab.in oft_slab.xml
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
!! Once the post-processing run is complete `bin/build_xdmf.py` can be used to generate `*.xmf` files that can be loaded by
!! [VisIt](https://visit-dav.github.io/visit-website/index.html), [ParaView](https://www.paraview.org/), or other visualization programs.
!!
!! \image html MUG_railgun_ex-J.png "Current density in coaxial (top) and parallel plate (bottom) cases"
!!
!! \image html MUG_railgun_ex-V.png "Axial flow velocity in coaxial (top) and parallel plate (bottom) cases"
!!
!!\section doc_mug_railgun_mesh Mesh Generation
!! Mesh files for the parallel plate (`MUG-railgun-tet.h5`) and coaxial (`MUG-coax_railgun-tet.h5`) are provided with this example.
!! Instructions to generate your own mesh for the geometry using [CUBIT](https://cubit.sandia.gov/) are also provided below.
!!
!!\subsection doc_mug_railgun_mesh_par Parallel plate/Slab
!! A suitable mesh for this example in slab geometry, with separation of 10 cm, can be created using
!! the CUBIT script below.
!!
! OFT_DOC_INCLUDE: MUG-railgun.jou 
!!
!!\subsection doc_mug_railgun_mesh_coax Coaxial/Cylindrical
!! A suitable mesh for this example in cylindrical geometry, with inner radius of 10 cm and outer radius of 20 cm, can be created using
!! the CUBIT script below.
!!
! OFT_DOC_INCLUDE: MUG-coax_railgun.jou 
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_cubit.py --in_file=MUG-railgun-tet.g --periodic_nodeset=-1
!!\endverbatim
!!
!! @note The option `--periodic_nodeset=-1` indicates that the surface corresponding to the last nodeset (python-style reverse indexing)
!! should be treated as periodic following reflection