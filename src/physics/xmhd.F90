!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file xmhd.F90
!
!> Module for modeling extended MHD with a mixed conforming/lagrange basis set.
!!
!! The equations currently solved are the equations of extended xMHD with
!! scalar temperature and pressure.
!!
!! \f[ \frac{\partial n}{\partial t} +  \nabla \cdot \left( n u \right) =
!! D \nabla^2 n \f]
!!
!! \f[ \rho \left[ \frac{\partial u}{\partial t} +  u \cdot \nabla u \right]
!! = J \times B - \nabla nk(T_e + T_i) - \nabla \cdot \Pi \f]
!!
!! \f[ \frac{\partial B}{\partial t} = - \nabla \times \left( -u \times B +
!! \eta J + \frac{1}{ne} \left( J \times B - \nabla nkT_e \right) + \frac{m_e}{ne^2}
!! \frac{\partial J}{\partial t} \right) \f]
!!
!! \f[ \frac{n}{\gamma-1} \left[ \frac{\partial T_i}{\partial t} + u \cdot \nabla T_i \right] =
!! -nkT_i \nabla \cdot u - \nabla \cdot q_i  + Q_i \f]
!!
!! For two temperature
!! \f[ \frac{n}{\gamma-1} \left[ \frac{\partial T_e}{\partial t} + u_e \cdot \nabla T_e \right] =
!! -nkT_e \nabla \cdot u_e - \nabla \cdot q_e  + Q_e \f]
!!
!! \f[ u_e = u - \frac{J}{ne} \f]
!!
!! **Thermal transport**
!!
!! \f[ q_s = - n \left[ \chi_{\parallel,s} \hat{b} \hat{b} + \chi_{\perp,s} \left( I
!! - \hat{b} \hat{b} \right) \right] \cdot \nabla T_s \f]
!!
!! With `xmhd_brag=T` and single temperature (\f$ T_e = T_i \f$)
!! \f[ \chi_{\parallel,i} = \chi_{\parallel,e}, \chi_{\perp,i} = min(\chi_{\perp,i},\chi_{\parallel,i}) \f]
!!
!! **Heat sources**
!!
!! For single temperature (\f$ T_e = T_i \f$)
!! \f[ Q_i = \left[ \eta J^2 - \left( \nabla u \right)^T : \Pi \right] / 2 \f]
!!
!! For two temperature
!! \f[ Q_i = - \left( \nabla u \right)^T : \Pi, Q_e = \eta J^2 \f]
!!
!! **Viscosity**
!!
!! \f[ W = \left( \nabla u + (\nabla u)^T - \frac{2}{3} I \nabla \cdot u \right) \f]
!!
!! With `visc_type='kin'`
!! \f[ \Pi = - \nu \nabla u \f]
!!
!! With `visc_type='iso'`
!! \f[ \Pi = - \nu W \f]
!!
!! With `visc_type='ani'`
!! \f[ \Pi = - \left[ \nu_{\parallel} \hat{b} \hat{b}
!! + \nu_{\perp} \left( I - \hat{b} \hat{b} \right) \right] \cdot W \f]
!!
!! Runtime options are available when using the driver (@ref xmhd::xmhd_run "xmhd_run")
!! and post-processing (@ref xmhd::xmhd_plot "xmhd_plot") routines, see their documentation
!! for available options.
!!
!! @note Electron mass is artificially enhanced by a factor using @ref xmhd::me_factor
!! "me_factor"
!!
!! @authors Chris Hansen
!! @date September 2013
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE xmhd
#define XMHD_ITCACHE 10
#if !defined(XMHD_RST_LEN)
#define XMHD_RST_LEN 5
#endif
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh, multigrid_level, multigrid_base_pushcc
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, map_list, vector_extrapolate, &
  oft_matrix, oft_matrix_ptr, oft_graph, oft_graph_ptr, oft_local_mat
USE oft_deriv_matrices, ONLY: oft_mf_matrix, oft_noop_matrix
USE oft_native_la, ONLY: oft_native_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_la_utils, ONLY: create_vector, create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_solver_xml, &
  create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: fem_max_levels, fem_common_linkage, oft_fem_type, oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type, oft_ml_fem_comp_type, oft_ml_fe_comp_vecspace
USE fem_utils, ONLY: fem_avg_bcc, fem_interp, cc_interp, cross_interp, &
  tensor_dot_interp, fem_partition, fem_dirichlet_diag, fem_dirichlet_vec
USE oft_lag_basis, ONLY: oft_lag_eval_all, oft_3D_lagrange_cast, &
  oft_lag_geval_all, oft_scalar_fem
USE oft_lag_operators, ONLY: oft_lag_vgetmop, oft_lag_vrinterp, oft_lag_vdinterp, &
  oft_lag_vproject, oft_lag_project_div, oft_lag_rinterp, oft_lag_ginterp, &
  lag_vbc_tensor, lag_vbc_diag
USE oft_hcurl_basis, ONLY: oft_hcurl_eval_all, oft_3D_hcurl_cast, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops, oft_hcurl_fem
USE oft_hcurl_operators, ONLY: oft_hcurl_rinterp
USE oft_h1_basis, ONLY: oft_h1_fem, oft_h1_geval_all, oft_3D_h1_cast
USE oft_h1_operators, ONLY: oft_h1_zeroi
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp, oft_hcurl_grad_cinterp, oft_hcurl_grad_dinterp, &
  oft_hcurl_grad_divout, hcurl_grad_jump_error, hcurl_grad_div
!---
USE diagnostic, ONLY: tfluxfun, vec_energy, weighted_vec_energy, scal_energy, &
  scal_int, field_probe
USE mhd_utils
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> Forcing object class for xMHD
!!
!! Used to apply arbitrary external fields or dirichlet boundary conditions
!! during time advance. If passed during a run the driver is called to modify
!! the time centered source vector and/or the initial field for the NL solve.
!!
!! @note If dirichlet boundary conditions are being used both the source vector
!! and initial field should be modified consistently. Otherwise these terms will
!! unnecessarily load the linear solver.
!---------------------------------------------------------------------------------
type, abstract, public :: oft_xmhd_driver
contains
  !> Modify vectors to apply forcing
  procedure(xmhd_driver_apply), deferred :: apply
  !> Setup driver from restart file
  procedure :: rst_load => xmhd_driver_rst_dummy
  !> Save driver info to restart file
  procedure :: rst_save => xmhd_driver_rst_dummy
end type oft_xmhd_driver
!---------------------------------------------------------------------------------
!> Probe object class for xMHD
!!
!! Used to extract probe signals from restart files during post-processing. If
!! passed to \ref xmhd::xmhd_plot "xmhd_plot" the object is called to sample and
!! output the desired probe signals.
!---------------------------------------------------------------------------------
type, abstract, public :: oft_xmhd_probe
contains
  !> Sample probe signals from solution
  procedure(xmhd_probe_apply), deferred :: apply
  !> Flush internal write buffer
  procedure :: flush => xmhd_probe_flush
end type oft_xmhd_probe
!---------------------------------------------------------------------------------
!> NL metric matrix for xMHD
!---------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_xmhd_errmatrix
  real(r8) :: dt = 1.E-3_r8 !< Time step
  real(r8) :: diag_vals(7) = 0.d0 !< Diagnostic values
  class(oft_vector), pointer :: up => NULL() !< Field at start of time step
contains
  !> Apply the matrix
  procedure :: apply_real => xmhd_errmatrix_apply
end type oft_xmhd_errmatrix
!---------------------------------------------------------------------------------
!> Mass metric matrix for Reduced MHD
!---------------------------------------------------------------------------------
type, public, extends(oft_noop_matrix) :: oft_xmhd_massmatrix
  real(r8) :: diag_vals(7) = 0.d0 !< Diagnostic values
  class(oft_vector), pointer :: u0 => NULL() !< Equilibrium field
contains
  !> Apply the matrix
  procedure :: apply_real => xmhd_massmatrix_apply
end type oft_xmhd_massmatrix
!------------------------------------------------------------------------------
!> xMHD operator container
!------------------------------------------------------------------------------
type :: xmhd_ops
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: b1_bc => NULL() !< B-field BC flag (Curl)
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: b2_bc => NULL() !< B-field BC flag (Grad)
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: v_bc => NULL() !< Velocity BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< Density BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: t_bc => NULL() !< Temperature BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: solid_node => NULL() !< Solid region flag
  CLASS(oft_matrix), POINTER :: J => NULL() !< Jacobian matrix
  CLASS(oft_matrix), POINTER :: interp => NULL() !< Interpolation matrix
  TYPE(oft_xmhd_errmatrix) :: A !< NL metric matrix
end type xmhd_ops
!------------------------------------------------------------------------------
!> Interpolate a xMHD operator fields
!------------------------------------------------------------------------------
type :: xmhd_interp_cache
  integer(i4) :: cell = -1 !< Cell index for values
  real(r8), contiguous, pointer, dimension(:) :: bcurl => NULL() !< HCurl values (Curl)
  real(r8), contiguous, pointer, dimension(:) :: bgrad => NULL() !< HCurl values (Grad)
  real(r8), contiguous, pointer, dimension(:) :: Te => NULL() !< Lagrange values
  real(r8), contiguous, pointer, dimension(:) :: N2 => NULL() !< Hyper-diff aux values
  real(r8), contiguous, pointer, dimension(:) :: J2 => NULL() !< Hyper-res aux values (HCurl)
  real(r8), contiguous, pointer, dimension(:,:) :: lf => NULL() !< Lagrange values (N,V,Ti)
end type xmhd_interp_cache
!------------------------------------------------------------------------------
!> Interpolate a xMHD operator fields
!------------------------------------------------------------------------------
type, extends(fem_interp) :: xmhd_interp
  real(r8), contiguous, pointer, dimension(:) :: bcurl_loc => NULL() !< Local HCurl values (Curl)
  real(r8), contiguous, pointer, dimension(:) :: bgrad_loc => NULL() !< Local HCurl values (Grad)
  real(r8), contiguous, pointer, dimension(:) :: Te_loc => NULL() !< Local Te values
  real(r8), contiguous, pointer, dimension(:) :: N2_loc => NULL() !< Local hyper-diff aux values
  real(r8), contiguous, pointer, dimension(:) :: J2_loc => NULL() !< Local hyper-res aux values (HCurl)
  real(r8), contiguous, pointer, dimension(:,:) :: lf_loc => NULL() !< Local Lagrange values (V,N,Ti)
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  class(oft_h1_fem), pointer :: grad_rep => NULL() !< Grad(H^1) FE representation
  class(oft_hcurl_fem), pointer :: curl_rep => NULL() !< H(Curl) FE representation
  class(oft_scalar_fem), pointer :: lag_rep => NULL() !< Lagrange FE representation
  type(xmhd_interp_cache), pointer :: cache(:) => NULL() !< Thread local field cache
contains
  !> Retrieve local values for interpolation
  procedure :: setup => xmhd_interp_setup
  !> Reconstruct field
  procedure :: interp => xmhd_interp_apply
  !> Destroy temporary internal storage
  procedure :: delete => xmhd_interp_delete
end type xmhd_interp
!------------------------------------------------------------------------------
!> Unpack object for local fields evaluated by @ref xmhd::xmhd_interp "xmhd_interp"
!------------------------------------------------------------------------------
type :: xmhd_loc_values
  REAL(r8) :: N = 0.d0 !< Density
  REAL(r8) :: Ti = 0.d0 !< Ion temperature (T in single fluid)
  REAL(r8) :: Te = 0.d0 !< Electron temperature
  REAL(r8) :: N2 = 0.d0 !< Hyper-diff aux field
  REAL(r8) :: B(3) = 0.d0 !< Magnetic field
  REAL(r8) :: V(3) = 0.d0 !< Velocity field
  REAL(r8) :: J2(3) = 0.d0 !< Hyper-res aux field
  REAL(r8) :: dN(3) = 0.d0 !< Density gradient
  REAL(r8) :: dTi(3) = 0.d0 !< Ion temperature gradient (T in single fluid)
  REAL(r8) :: dTe(3) = 0.d0 !< Electron temperature gradient
  REAL(r8) :: dN2(3) = 0.d0 !< Hyper-diff aux field gradient
  REAL(r8) :: J(3) = 0.d0 !< Current density (normalized, J = Curl(B))
  REAL(r8) :: J2c(3) = 0.d0 !< Hyper-res aux field curl
  REAL(r8) :: dV(3,3) = 0.d0 !< Velocity shear tensor
end type xmhd_loc_values
!------------------------------------------------------------------------------
!> Object for sub-fields used by @ref xmhd::oft_xmhd_push "oft_xmhd_push"
!! and @ref xmhd::oft_xmhd_pop "oft_xmhd_pop"
!------------------------------------------------------------------------------
type, public :: xmhd_sub_fields
  CLASS(oft_vector), POINTER :: B => NULL() !< Magnetic field "H(Curl) + Grad(H^1)"
  CLASS(oft_vector), POINTER :: V => NULL() !< Velocity "Vector Lagrange"
  CLASS(oft_vector), POINTER :: Ne => NULL() !< Electron density "Lagrange"
  CLASS(oft_vector), POINTER :: Ti => NULL() !< Ion temperature "Lagrange"
  CLASS(oft_vector), POINTER :: Te => NULL() !< Electron temperature "Lagrange"
  CLASS(oft_vector), POINTER :: N2 => NULL() !< Hyper-diff aux field "Lagrange"
  CLASS(oft_vector), POINTER :: J2 => NULL() !< Hyper-res aux field "H(Curl)"
end type xmhd_sub_fields
!---Interface definitions
abstract interface
!------------------------------------------------------------------------------
!> Modify solution vectors to apply forcing
!------------------------------------------------------------------------------
  subroutine xmhd_driver_apply(self,up,u,v,t,dt)
    import oft_xmhd_driver, oft_vector, r8
    class(oft_xmhd_driver), intent(inout) :: self !< Forcing object
    class(oft_vector), intent(inout) :: up !< Solution at start of time step
    class(oft_vector), intent(inout) :: u !< Initial guess for solution
    class(oft_vector), intent(inout) :: v !< RHS for non-linear solve
    real(r8), intent(in) :: t !< Current solution time
    real(r8), intent(in) :: dt !< Current solution timestep
  end subroutine xmhd_driver_apply
!------------------------------------------------------------------------------
!> Sample probe signals from solution
!------------------------------------------------------------------------------
  subroutine xmhd_probe_apply(self,sub_fields,t)
    import oft_xmhd_probe, xmhd_sub_fields, r8
    class(oft_xmhd_probe), intent(inout) :: self !< Probe object
    type(xmhd_sub_fields), intent(inout) :: sub_fields !< Unpacked fields at current time
    real(r8), intent(in) :: t !< Current time
  end subroutine xmhd_probe_apply
end interface
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC, extends(oft_ml_fe_comp_vecspace) :: ml_xmhd_vecspace
contains
  !> Needs docs
  PROCEDURE :: interp => oft_xmhd_interp
  !> Needs docs
  PROCEDURE :: inject => oft_xmhd_inject
end type ml_xmhd_vecspace
!------------------------------------------------------------------------------
! Global variables
!------------------------------------------------------------------------------
INTEGER(i4), PARAMETER :: xmhd_rst_version = 3 !< Restart file version number
TYPE(xml_node), POINTER :: xmhd_root_node => NULL() !< xMHD XML node
TYPE(xml_node), POINTER :: xmhd_pre_node => NULL() !< preconditioner XML node
!---Equation control
LOGICAL :: xmhd_jcb = .TRUE. !< Include JxB force on fluid
LOGICAL :: xmhd_advec = .TRUE. !< Include fluid advection
LOGICAL :: xmhd_hall = .FALSE. !< Include Hall terms
LOGICAL, PUBLIC :: xmhd_adv_b = .TRUE. !< Advance magnetic field
LOGICAL :: xmhd_adv_den = .TRUE. !< Advance density
LOGICAL :: xmhd_adv_temp = .TRUE. !< Advance temperature
LOGICAL :: xmhd_ohmic = .FALSE. !< Resistive heating
LOGICAL :: xmhd_visc_heat = .FALSE. !< Viscous heating
LOGICAL :: xmhd_brag = .FALSE. !< Braginskii thermal conduction
LOGICAL :: xmhd_linear = .FALSE. !< Compute linearized matrix
LOGICAL :: xmhd_upwind = .TRUE. !< Use upwinding
LOGICAL :: xmhd_rw = .FALSE. !< Resistive wall present
LOGICAL :: xmhd_monitor_div = .FALSE. !< Monitor div error?
LOGICAL :: xmhd_two_temp = .FALSE. !< Run 2 Temperature model
LOGICAL :: xmhd_therm_equil = .FALSE. !< Include cross-species thermalization
LOGICAL :: xmhd_diss_centered = .FALSE. !< Time-center dissipation terms
LOGICAL :: xmhd_skip_update = .FALSE. !< Internal flag for skipping Jac update
INTEGER(i4) :: te_ind = -1 !< Index of electron temperature field
INTEGER(i4) :: n2_ind = -1 !< Index of hyper-diff aux field
INTEGER(i4) :: j2_ind = -1 !< Index of hyper-res aux field
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: xmhd_mat_mask !< Matrix block mask
!---Boundary conditions
LOGICAL :: xmhd_bnorm_force = .TRUE. !< Force B-norm to zero in cleaning
CHARACTER(LEN=2) :: bbc = 'bc' !< Magnetic field BC ('bc','ic')
CHARACTER(LEN=4) :: vbc = 'all' !< Velocity BC ('all','norm','tang')
LOGICAL :: xmhd_vbcdir = .FALSE. !< Velocity BC type flag
INTEGER(i4) :: vbc_type = 1 !< Velocity BC type
CHARACTER(LEN=1) :: nbc = 'd' !< Density BC ('d','n')
CHARACTER(LEN=1) :: tbc = 'd' !< Temperature BC ('d','n')
CHARACTER(LEN=3) :: visc_type = 'kin' !< Viscosity type ('kin','iso','ani')
INTEGER(i4) :: visc_itype = 1 !< Viscosity type flag (1->kin,2->iso,3->ani)
LOGICAL, ALLOCATABLE, DIMENSION(:) :: solid_cell !< Solid region flag
!---Physical parameters
REAL(r8) :: jac_dt = 1.E-3_r8 !< Time step
REAL(r8) :: dt_scale = 1._r8 !< Time step scaling factor
REAL(r8) :: eta = 1.E-2_r8 !< Resistivity \f$ \eta / \mu_0 \f$
REAL(r8) :: eta_temp = -1._r8 !< Reference temperature for resistivity
REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_reg !< Resistivity scale factors
REAL(r8) :: eta_hyper = -1.d0 !< Hyper resistivity
REAL(r8) :: nu_par = 0.d0 !< Fluid (parallel) viscosity (note: dynamic viscosity is nu=nu_par*den_scale/n0)
REAL(r8) :: nu_perp = 0.d0 !< Fluid perpendicular viscosity if anisotropic, as nu_par
REAL(r8) :: m_ion = -1.d0 !< Ion mass in SI units
REAL(r8) :: mu_ion = 2.d0 !< Ion mass in atomic units
REAL(r8) :: d_dens = 1.E0_r8 !< Particle diffusivity
REAL(r8) :: d2_dens = -1.d0 !< Particle hyper-diffusivity
REAL(r8) :: te_factor = 1.d0 !< Electron temperature factor (single temperature)
REAL(r8) :: me_factor = 100.E0_r8 !< Artificial electron mass factor
REAL(r8) :: k_boltz = elec_charge !< Boltzmann constant
REAL(r8) :: temp_gamma = 5.d0/3.d0 !< Ratio of specific heats
REAL(r8) :: kappa_par = 1.d0 !< Parallel thermal conductivity factor
REAL(r8) :: kappa_perp = 1.d0 !< Perpendicular thermal conductivity factor
REAL(r8), PUBLIC :: vel_scale = 1.d3 !< Velocity scale factor
REAL(r8), PUBLIC :: den_scale = -1.d0 !< Density scale factor
REAL(r8), PUBLIC :: n2_scale = -1.d0 !< Hyper-diffusivity aux variable scale factor
REAL(r8), PUBLIC :: j2_scale = -1.d0  !< Hyper-resistivity aux variable scale factor
REAL(r8), PUBLIC :: den_floor = -1.d0 !< Density floor
REAL(r8), PUBLIC :: temp_floor = -1.d0 !< Temperature floor
REAL(r8), PUBLIC :: g_accel = 0.d0 !< Gravity (oriented in the -Z direction)
REAL(r8) :: xmhd_eps = 1.d-10 !< Epsilon for magnetic field normalization
INTEGER(i4) :: xmhd_taxis = 3 !< Axis for toroidal flux and current
!> Resistivity profile
PROCEDURE(oft_1d_func), PUBLIC, POINTER :: res_profile => NULL()
!---Multi-level environment
INTEGER(i4) :: xmhd_blevel = 0 !< Highest level on base meshes
INTEGER(i4) :: xmhd_lev = 1 !< Active FE level
INTEGER(i4) :: xmhd_level = 1 !< Active FE level
INTEGER(i4) :: xmhd_lin_level = 1 !< Highest linear element level
INTEGER(i4) :: xmhd_nlevels = 1 !< Number of total levels
INTEGER(i4) :: xmhd_minlev = -1 !< Lowest MG level
INTEGER(i4), DIMENSION(fem_max_levels) :: nu_xmhd = 1 !< Number of smoother iterations
!---Operators and preconditioning
TYPE(oft_fem_comp_type), POINTER :: xmhd_rep => NULL() !< Active field representation
TYPE(oft_ml_fem_comp_type), TARGET :: ML_xmhd_rep !< ML container for field representation
TYPE(ml_xmhd_vecspace), TARGET :: xmhd_ml_vecspace
TYPE(oft_mf_matrix), TARGET :: mfmat !< Matrix free operator
TYPE(xmhd_ops), POINTER :: oft_xmhd_ops => NULL() !< Operator container
TYPE(xmhd_ops), POINTER, DIMENSION(:) :: oft_xmhd_ops_ML => NULL() !< MG operator container
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_int => NULL() !< MG interpolation operators
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_J => NULL() !< MG Jacobian operators
CLASS(oft_solver), POINTER :: xmhd_pre => NULL() !< Preconditioner object
INTEGER(i4) :: xmhd_prefreq = 30 !< Desired update frequency for preconditioner
INTEGER(i4) :: xmhd_opcount = 1 !< Number of time steps since preconditioner update
INTEGER(i4), DIMENSION(fem_max_levels) :: xmhd_nparts = 1 !< Number of local partitions for preconditioning
REAL(r8) :: xmhd_opdt = 1.d0 !< Time step for current preconditioner object
LOGICAL :: xmhd_mfnk = .FALSE. !< Use Jacobian free non-linear advance
!---Interpolator cache variables
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: xmhd_lag_rop => NULL()
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: xmhd_lag_gop => NULL()
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: xmhd_hgrad_rop => NULL()
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: xmhd_hcurl_rop => NULL()
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: xmhd_hcurl_cop => NULL()
!$omp threadprivate(xmhd_lag_rop,xmhd_lag_gop,xmhd_hgrad_rop,xmhd_hcurl_rop,xmhd_hcurl_cop)
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: neg_source,neg_flag
CLASS(multigrid_mesh), POINTER :: mg_mesh
CLASS(oft_mesh), POINTER :: mesh
!
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: xmhd_ML_lagrange => NULL()
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: xmhd_ML_H1 => NULL()
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: xmhd_ML_H1grad => NULL()
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: xmhd_ML_hcurl => NULL()
TYPE(oft_ml_fem_comp_type), POINTER, PUBLIC :: xmhd_ML_hcurl_grad => NULL()
TYPE(oft_ml_fem_comp_type), POINTER, PUBLIC :: xmhd_ML_vlagrange => NULL()
!
CLASS(oft_scalar_fem), POINTER :: oft_lagrange => NULL()
CLASS(oft_hcurl_fem), POINTER :: oft_hcurl => NULL()
CLASS(oft_h1_fem), POINTER :: oft_hgrad => NULL()
!------------------------------------------------------------------------------
! Exported interfaces
!------------------------------------------------------------------------------
PUBLIC xmhd_run, oft_xmhd_push, oft_xmhd_pop, xmhd_driver_apply, xmhd_minlev
PUBLIC xmhd_taxis, xmhd_plot, xmhd_probe_apply, xmhd_lin_run, xmhd_bnorm_force
CONTAINS
!------------------------------------------------------------------------------
!> Read settings for extended MHD model from input files
!!
!! Runtime options are set in the main input file using the `xmhd_options` group.
!!
!! @warning Most default physical values in the namelist below will not be appropriate
!! for new simulations. They are merely placeholders to indicate the type of the input.
!!
!! **Option group:** `xmhd_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `xmhd_jcb=T`       | Apply magnetic force to fluid | bool |
!! | `xmhd_advec=T`     | Advect fluid | bool |
!! | `xmhd_adv_den=T`   | Advance density | bool |
!! | `xmhd_adv_temp=T`  | Advance temperature | bool |
!! | `xmhd_hall=F`      | Include Hall terms in Ohm's law | bool |
!! | `xmhd_ohmic=F`     | Include Ohmic heating | bool |
!! | `xmhd_visc_heat=F` | Include Viscous heating | bool |
!! | `xmhd_brag=F`      | Use Braginskii thermal conduction coefficients | bool |
!! | `xmhd_upwind=T`    | Use upwinding for advection terms | bool |
!! | `xmhd_therm_equil=F` | Include thermal equilibration between electrons/ions | bool |
!! | `xmhd_diss_centered=F` | Use time-centered dissipation terms for (V, Ti, Te) | bool |
!! | `bbc="bc"`         | Magnetic field boundary condition ("bc", "ic") | str(2) |
!! | `vbc="all"`        | Velocity boundary condition ("all") | str(3) |
!! | `nbc="d"`          | Density boundary condition ("d", "n") | str(1) |
!! | `tbc="d"`          | Temperature boundary condition ("d", "n") | str(1) |
!! | `visc_type="kin"`  | Fluid viscosity type ("kin", "iso", "ani") | str(3) |
!! | `dt=1`             | Maximum allowed time step | real |
!! | `eta=1`            | Plasma resistivity, in diffusivity units | real |
!! | `eta_temp=-1`      | Base resistivity temperature (neg to disable) | real |
!! | `eta_hyper=-1`     | Plasma hyper-resistivity, in diffusivity units (neg to disable) | real |
!! | `nu_par=0`         | Plasma (parallel) viscosity, in diffusivity units | real |
!! | `nu_perp=0`        | Plasma perpendicular viscosity, in diffusivity units | real |
!! | `d_dens=1`         | Plasma diffusivity, in diffusivity units | real |
!! | `d2_dens=-1`        | Plasma hyper-diffusivity, in diffusivity units (neg to disable) | real |
!! | `kappa_par=1`      | Plasma thermal conduction, in diffusivity units | real |
!! | `kappa_perp=1`     | Plasma thermal conduction, in diffusivity units | real |
!! | `nsteps=2000`      | Number of time steps to take | int |
!! | `rst_freq=10`      | Frequency of dump files | int |
!! | `lin_tol=1E-9`     | Linear solver tolerance | real |
!! | `nl_tol=1E-5`      | Non-Linear solver tolerance | real |
!! | `nu_xmhd=1`        | Number of smoother iterations for MG preconditioner | int [nlevels] |
!! | `xmhd_nparts=1`    | Number of local partitions for factor based preconditioning | int |
!! | `nclean=500`       | Frequency of divergence cleaning | int |
!! | `rst_ind=0`        | Index of restart file to initialize from | int |
!! | `maxextrap=2`      | Extrapolation order for NL initial guess | int |
!! | `ittarget=60`      | Target number of linear iterations per time step | int |
!! | `mu_ion=2`         | Ion mass in atomic units | real |
!! | `me_factor=100`    | Electron mass factor | real |
!! | `te_factor=1`      | Electron temperature factor (single temperature) | real |
!! | `xmhd_prefreq=30`  | Frequency of full preconditioner update | int |
!! | `xmhd_monitor_div=F` | Monitor divergence error | bool |
!! | `xmhd_mfnk=F`      | Use Matrix-Free NL solve | bool |
!
!------------------------------------------------------------------------------
SUBROUTINE xmhd_read_settings(dt,lin_tol,nl_tol,rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget,nl_update)
real(r8), intent(out) :: dt !< Maximum timestep
real(r8), intent(out) :: lin_tol !< Linear solver tolerance
real(r8), intent(out) :: nl_tol !< Nonlinear solver tolerance
integer(i4), intent(out) :: rst_ind !< Index of file for restart (rst_ind>0)
integer(i4), intent(out) :: nsteps !< Number of timesteps
integer(i4), intent(out) :: rst_freq !< Frequency to save restart files
integer(i4), intent(out) :: nclean !< Frequency to clean divergence
integer(i4), intent(out) :: maxextrap !< Extrapolation order for initial guess
integer(i4), intent(out) :: ittarget !< Maximum number of linear iterations
integer(i4), intent(out) :: nl_update !< Maximum number of linear iterations
integer(i4) :: io_unit,ierr
!---
namelist/xmhd_options/xmhd_jcb,xmhd_advec,xmhd_adv_den,xmhd_adv_temp,xmhd_hall,xmhd_ohmic, &
  xmhd_visc_heat,xmhd_brag,xmhd_upwind,xmhd_therm_equil,bbc,vbc,nbc,tbc,visc_type,dt,eta, &
  eta_temp,nu_par,nu_perp,d_dens,kappa_par,kappa_perp,nsteps,rst_freq,lin_tol,nl_tol,nu_xmhd, &
  xmhd_nparts,nclean,rst_ind,maxextrap,ittarget,nl_update,mu_ion,me_factor,te_factor,xmhd_prefreq, &
  xmhd_monitor_div,xmhd_mfnk,xmhd_diss_centered,eta_hyper,d2_dens
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Set defaults
!------------------------------------------------------------------------------
lin_tol=1.E-9_r8
nl_tol=1.E-5_r8
rst_ind=0
nsteps=2000
rst_freq=10
nclean=500
maxextrap=2
ittarget=60
nl_update=2
!------------------------------------------------------------------------------
! Read-in Parameters
!------------------------------------------------------------------------------
open(NEWUNIT=io_unit,FILE=oft_env%ifile)
read(io_unit,xmhd_options,IOSTAT=ierr)
close(io_unit)
if(ierr<0)call oft_abort('No MHD options found in input file.','xmhd_read_settings',__FILE__)
if(ierr>0)call oft_abort('Error parsing MHD options in input file.','xmhd_read_settings',__FILE__)
!---Look for xMHD node
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"xmhd",xmhd_root_node,ierr)
  IF(ierr==0)THEN
    !---Look for pre node
    CALL xml_get_element(xmhd_root_node,"pre",xmhd_pre_node,ierr)
  END IF
END IF
#endif
!------------------------------------------------------------------------------
! Check settings and setup
!------------------------------------------------------------------------------
m_ion=mu_ion*proton_mass ! Scale by proton mass
IF(m_ion<0.d0)CALL oft_abort('Invalid Ion mass','xmhd_read_settings',__FILE__)
!---Set velocity BC flag
IF(vbc(1:3)/='all')xmhd_vbcdir=.TRUE.
SELECT CASE(TRIM(vbc))
  CASE("norm")
    vbc_type=1
  CASE("tang")
    vbc_type=2
END SELECT
!---Check viscosity type
SELECT CASE(visc_type)
  CASE('kin')
    visc_itype=1
    nu_perp=nu_par
  CASE('iso')
    visc_itype=2
    nu_perp=nu_par
  CASE('ani')
    visc_itype=3
  CASE DEFAULT
    CALL oft_abort('Unknown viscosity type','xmhd_read_settings',__FILE__)
END SELECT
!---Disable JxB if no B advance
IF(.NOT.xmhd_adv_b)THEN
  xmhd_jcb=.FALSE.
  eta_hyper=-1.d0
END IF
IF(.NOT.xmhd_adv_den)d2_dens=-1.d0
IF(.NOT.xmhd_two_temp)xmhd_therm_equil=.FALSE.
CALL xmhd_setup_regions
DEBUG_STACK_POP
END SUBROUTINE xmhd_read_settings
!------------------------------------------------------------------------------
!> Main driver subroutine for extended MHD time advance
!!
!! Runtime options are set in the main input file using the group
!! `xmhd_options` group, see \ref xmhd_read_settings.
!------------------------------------------------------------------------------
subroutine xmhd_run(initial_fields,driver,probes,profile_only)
TYPE(xmhd_sub_fields), INTENT(inout) :: initial_fields !< Initial conditions
CLASS(oft_xmhd_driver), OPTIONAL, INTENT(inout) :: driver !< Forcing object
CLASS(oft_xmhd_probe), OPTIONAL, INTENT(inout) :: probes !< Probe object
LOGICAL, OPTIONAL, INTENT(in) :: profile_only !< Profile operator timing and stop?
!---H(Curl) + Grad(H^1) divout solver
TYPE(oft_hcurl_grad_divout) :: divout
!---Jacobian solver
type(oft_nksolver) :: nksolver
type(oft_native_gmres_solver), target :: solver
integer(i4) :: nlevels
integer(i4), allocatable, dimension(:) :: levels
!---Local variables
class(oft_vector), pointer :: u,v,up
type(xmhd_sub_fields) :: sub_fields
type(oft_hcurl_grad_rinterp), target :: Bfield
type(oft_hcurl_grad_dinterp) :: divfield
type(oft_timer) :: mytimer
real(r8), pointer :: vals(:)
integer(i4) :: i,j,k,ierr,io_unit,io_stat,rst_version,rst_tmp,baseit,nredo,nl_update
integer(i4) :: itcount(XMHD_ITCACHE)
real(r8) :: dt,dthist(XMHD_ITCACHE)
real(r8) :: mag,div_err,mer,merp,ver,verp,gerr,cerr,verr,elapsed_time
real(r8) :: fac,lramp,tflux,tcurr,t,dtin,div_error,jump_error,derror
real(r8) :: ndens,npart,temp_avg,tempe_avg,mesh_vol,tmpint(2)
character(LEN=XMHD_RST_LEN) :: rst_char
LOGICAL :: force_refactor,exists,rst
TYPE(oft_h1_zeroi), TARGET :: h1_zeroi
!---Extrapolation fields
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(11)
!---Input variables
real(r8) :: lin_tol,nl_tol,scale_tmp(4)
integer(i4) :: rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget
DEBUG_STACK_PUSH
mg_mesh=>xmhd_ML_hcurl%ml_mesh
IF(.NOT.oft_3D_hcurl_cast(oft_hcurl,xmhd_ML_hcurl%current_level))CALL oft_abort("Invalid Curl FE object","xmhd_run",__FILE__)
IF(.NOT.oft_3D_lagrange_cast(oft_lagrange,xmhd_ML_lagrange%current_level))CALL oft_abort("Invalid Lagrange FE object","xmhd_run",__FILE__)
IF(.NOT.oft_3D_h1_cast(oft_hgrad,xmhd_ML_H1grad%current_level))CALL oft_abort("Invalid Grad FE object","xmhd_run",__FILE__)
mesh=>oft_hcurl%mesh
!------------------------------------------------------------------------------
! Read-in Parameters
!------------------------------------------------------------------------------
IF(ASSOCIATED(initial_fields%Te))xmhd_two_temp=.TRUE.
CALL xmhd_read_settings(dt,lin_tol,nl_tol,rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget,nl_update)
IF(den_scale<0.d0)den_scale=SQRT(initial_fields%Ne%dot(initial_fields%Ne)/REAL(initial_fields%Ne%ng,8))
IF((d2_dens>0.d0).AND.(n2_scale<0.d0))n2_scale=den_scale*(REAL(oft_lagrange%order,8)/mesh%hrms)**2
IF((eta_hyper>0.d0).AND.(j2_scale<0.d0))j2_scale=(REAL(oft_lagrange%order,8)/mesh%hrms)**2
!------------------------------------------------------------------------------
! Setup ML environment
!------------------------------------------------------------------------------
CALL xmhd_setup_rep
!------------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------------
call oft_xmhd_create(u)
call oft_xmhd_create(up)
call oft_xmhd_create(v)
!------------------------------------------------------------------------------
! Create extrapolation fields
!------------------------------------------------------------------------------
IF(maxextrap>0)THEN
  ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
  DO i=1,maxextrap
    CALL oft_xmhd_create(extrap_fields(i)%f)
    extrapt(i)=0.d0
  END DO
  nextrap=0
END IF
!------------------------------------------------------------------------------
! Create divergence cleaner
!------------------------------------------------------------------------------
call xmhd_ML_hcurl_grad%vec_create(sub_fields%B)
CALL divout%setup(xmhd_ML_hcurl_grad,"grnd")
divout%pm=.TRUE.
IF(TRIM(bbc)=="ic")divout%keep_boundary=.TRUE.
!------------------------------------------------------------------------------
! Setup history file
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd.hist', desc="History file for non-linear xMHD run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('tflux','r4', desc="Toroidal flux [Wb]")
  CALL hist_file%add_field('tcurr','r4', desc="Toroidal current [A*mu0]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('jerr', 'r4', desc="Magnetic jump error")
  CALL hist_file%add_field('derr', 'r4', desc="Magnetic div error")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average ion temperature [eV]")
  CALL hist_file%add_field('te',   'r4', desc="Average electron temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
END IF
!------------------------------------------------------------------------------
! Import initial fields
!------------------------------------------------------------------------------
dtin=dt
104 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
WRITE(rst_char,104)rst_ind
READ(rst_char,104,IOSTAT=io_stat)rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=rst_ind))CALL oft_abort("Step count exceeds format width", "xmhd_run", __FILE__)
rst=oft_file_exist(TRIM('xMHD_'//rst_char//'.rst'))
IF(rst.AND.(rst_ind/=0))THEN
  CALL oft_xmhd_rst_load(u, 'xMHD_'//rst_char//'.rst', 'U', t, dt, rst_version_out=rst_version)
  IF(PRESENT(driver))CALL driver%rst_load('xMHD_'//rst_char//'.rst')
  IF(rst_version<2)THEN
    dt = dt*2.d0
  END IF
  IF(dt>dtin)dt=dtin
  !---Setup history file
  IF(oft_env%head_proc)THEN
    INQUIRE(FILE='xmhd.hist',EXIST=exists)
    IF(.NOT.exists)THEN
      CALL oft_warn('No "xmhd.hist" file present for restart, new file created.')
      CALL hist_file%write_header
    END IF
  END IF
ELSE
  IF(rst_ind>0)CALL oft_abort('Restart does not exist.','xmhd_run',__FILE__)
  CALL initial_fields%V%scale(1.d0/vel_scale)
  CALL initial_fields%Ne%scale(1.d0/den_scale)
  CALL oft_xmhd_push(initial_fields,u)
  t=0.d0
  !---Create initial conditions restart file
  WRITE(rst_char,104)0
  CALL oft_xmhd_rst_save(u, t, dt, 'xMHD_'//rst_char//'.rst', 'U')
  IF(PRESENT(driver))CALL driver%rst_save('xMHD_'//rst_char//'.rst')
  !---Setup history file
  IF(oft_env%head_proc)CALL hist_file%write_header
END IF
IF(oft_env%head_proc)CALL hist_file%open ! Open history file
itcount=ittarget
dthist=dt
xmhd_opdt=dt
!------------------------------------------------------------------------------
! Setup Operators
!------------------------------------------------------------------------------
call xmhd_alloc_ops
jac_dt=dt/2.d0
ALLOCATE(neg_flag(3,mesh%nc),neg_source(3,mesh%nc))
neg_flag=0.d0
neg_source=0.d0
call xmhd_set_ops(u)
!---Force velocity to zero in solid regions
IF(xmhd_rw)THEN
  NULLIFY(vals)
  DO j=1,3
    CALL initial_fields%V%get_local(vals,j)
    DO i=1,oft_lagrange%ne
      IF(oft_xmhd_ops%solid_node(i))vals(i)=0.d0
    END DO
    CALL initial_fields%V%restore_local(vals,j)
  END DO
  DEALLOCATE(vals)
END IF
!------------------------------------------------------------------------------
! Compute energies
!------------------------------------------------------------------------------
mesh_vol=SUM(mesh%cv)
mesh_vol=oft_mpi_sum(mesh_vol)
!
oft_xmhd_ops%A%dt=-ABS(dt/2.d0)
CALL oft_xmhd_ops%A%apply(u,v)
!---Get diagnostics
mer=oft_xmhd_ops%A%diag_vals(1)/2.d0
ver=oft_xmhd_ops%A%diag_vals(2)*m_ion/2.d0
tflux=oft_xmhd_ops%A%diag_vals(3)
tcurr=oft_xmhd_ops%A%diag_vals(4)
npart=oft_xmhd_ops%A%diag_vals(5)/mesh_vol
temp_avg=oft_xmhd_ops%A%diag_vals(6)/mesh_vol
IF(xmhd_two_temp)THEN
  tempe_avg=oft_xmhd_ops%A%diag_vals(7)/mesh_vol
ELSE
  tempe_avg=-1.d0
END IF
!
IF(xmhd_monitor_div)THEN
  CALL oft_xmhd_pop(u,sub_fields)
  Bfield%u=>sub_fields%B
  CALL Bfield%setup(oft_hcurl,oft_hgrad)
  !---Compute jump error
  jump_error=hcurl_grad_jump_error(xmhd_ML_hcurl_grad%current_level,sub_fields%B,oft_hcurl%quad%order)
  divfield%u=>sub_fields%B
  CALL divfield%setup(oft_hcurl,oft_hgrad)
  derror=scal_energy(mesh,divfield,oft_hcurl%quad%order)
ELSE
  jump_error=-1.d0
  derror=-1.d0
END IF
!---Setup B-norm
IF(.NOT.xmhd_bnorm_force)THEN
  IF(.NOT.xmhd_monitor_div)CALL oft_xmhd_pop(u,sub_fields)
  CALL xmhd_ML_H1grad%vec_create(divout%bnorm)
  CALL hcurl_grad_div(xmhd_ML_hcurl_grad%current_level,sub_fields%B,divout%bnorm)
  h1_zeroi%ML_H1_rep=>xmhd_ML_H1grad
  CALL h1_zeroi%apply(divout%bnorm)
END IF
!------------------------------------------------------------------------------
! Print run information
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
100 FORMAT(2X,A,L1)
101 FORMAT(2X,2A)
102 FORMAT(2X,A,ES11.3)
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting non-linear xMHD run'
  WRITE(*,100)'Lorentz    = ',xmhd_jcb
  WRITE(*,100)'V-Advec    = ',xmhd_advec
  WRITE(*,100)'B-Adv      = ',xmhd_adv_b
  WRITE(*,100)'N-Adv      = ',xmhd_adv_den
  WRITE(*,100)'T-Adv      = ',xmhd_adv_temp
  IF(xmhd_two_temp)THEN
    WRITE(*,100)'Two Temp   = ',xmhd_two_temp
  ELSE
    WRITE(*,102)'Te/Ti    = ',te_factor
  END IF
  WRITE(*,100)'Hall       = ',xmhd_hall
  WRITE(*,100)'Ohmic      = ',xmhd_ohmic
  WRITE(*,100)'Visc Heat  = ',xmhd_visc_heat
  WRITE(*,100)'Upwinding  = ',xmhd_upwind
  WRITE(*,101)'B-BC       = ',bbc
  WRITE(*,101)'V-BC       = ',vbc
  IF(xmhd_adv_den)WRITE(*,101) 'N-BC       = ',nbc
  IF(xmhd_adv_temp)WRITE(*,101)'T-BC       = ',tbc
  WRITE(*,102)'t0         = ',t
  WRITE(*,102)'dt         = ',dt
  WRITE(*,102)'eta0       = ',eta
  WRITE(*,102)'eta0_Temp  = ',eta_temp
  WRITE(*,101)'visc_type  = ',visc_type
  WRITE(*,102)'  nu_par   = ',nu_par
  WRITE(*,102)'  nu_perp  = ',nu_perp
  IF(xmhd_adv_den)WRITE(*,102)'D_dens     = ',d_dens
  IF(xmhd_adv_temp)THEN
    WRITE(*,102)'kappa_par  = ',kappa_par
    WRITE(*,102)'kappa_perp = ',kappa_perp
  END IF
  IF(temp_floor<0.d0)WRITE(*,102)'temp_floor = ',temp_floor
  IF(den_floor<0.d0)WRITE(*,102)'den_floor  = ',den_floor
  WRITE(*,102)'L-Tol      = ',lin_tol
  WRITE(*,102)'NL-Tol     = ',nl_tol
  WRITE(*,102)'Mag Energy = ',mer
  WRITE(*,102)'Kin Energy = ',ver
  IF(xmhd_monitor_div)WRITE(*,'(2X,A,2ES11.3)')'Div Error  = ',derror,jump_error
  WRITE(*,102)'# B-DOFs   = ',REAL(oft_hcurl%global%ne+oft_hgrad%global%ne,8)
  WRITE(*,102)'# V-DOFs   = ',REAL(3*oft_lagrange%global%ne,8)
  WRITE(*,102)'# N-DOFs   = ',REAL(oft_lagrange%global%ne,8)
  WRITE(*,102)'# T-DOFs   = ',REAL(oft_lagrange%global%ne,8)
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
END IF
IF(PRESENT(profile_only))THEN
  IF(profile_only)CALL xmhd_profile(u)
END IF
!------------------------------------------------------------------------------
! Setup linear solver
!------------------------------------------------------------------------------
IF(xmhd_mfnk)THEN
  CALL up%set(1.d0)
  CALL up%set(mesh%hrms, iblock=1)
  CALL up%set(mesh%hrms, iblock=2)
  CALL mfmat%setup(u,oft_xmhd_ops%A,up)
  mfmat%b0=1.d-5
  solver%A=>mfmat
ELSE
  solver%A=>oft_xmhd_ops%J
END IF
solver%its=200
solver%atol=lin_tol
solver%itplot=1
solver%nrits=20
solver%pm=oft_env%pm
!------------------------------------------------------------------------------
! Setup Preconditioner
!------------------------------------------------------------------------------
nlevels=xmhd_nlevels-xmhd_minlev+1
NULLIFY(solver%pre)
IF(nlevels==1)THEN
  IF(ASSOCIATED(xmhd_pre_node))THEN
    CALL create_solver_xml(solver%pre,xmhd_pre_node)
  ELSE
    CALL create_diag_pre(solver%pre)
  END IF
ELSE
  ALLOCATE(levels(nlevels))
  levels=(/(i,i=xmhd_minlev,xmhd_nlevels)/)
  xmhd_ml_vecspace%ML_FE_rep=>ML_xmhd_rep
  CALL create_mlpre(solver%pre,ml_J,levels,nlevels=nlevels, &
    ml_vecspace=xmhd_ml_vecspace, &
    stype=2,nu=nu_xmhd(xmhd_minlev:xmhd_nlevels),xml_root=xmhd_pre_node)
  DEALLOCATE(levels)
END IF
xmhd_pre=>solver%pre
xmhd_pre%A=>oft_xmhd_ops%J
!------------------------------------------------------------------------------
! Setup Newton solver
!------------------------------------------------------------------------------
nksolver%A=>oft_xmhd_ops%A
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
IF(xmhd_mfnk)THEN
  nksolver%J_update=>xmhd_mfnk_update
  nksolver%up_freq=1
ELSE
  nksolver%J_update=>xmhd_set_ops
  nksolver%up_freq=4
END IF
!------------------------------------------------------------------------------
! Begin time stepping
!------------------------------------------------------------------------------
baseit=-1
nredo=0
DO i=1,nsteps
    merp=mer
    verp=ver
    !---Clean accumulated divergence error
    IF(MOD(rst_ind+i,nclean)==0)THEN
      CALL oft_xmhd_pop(u,sub_fields)
      IF(oft_env%head_proc)WRITE(*,*)'  Cleaning Divergence'
      CALL divout%apply(sub_fields%B)
      CALL oft_xmhd_push(sub_fields,u)
    END IF
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Update extrapolation fields
    IF(maxextrap>0)THEN
      DO j=maxextrap,2,-1
        CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
        extrapt(j)=extrapt(j-1)
      END DO
      IF(nextrap<maxextrap)nextrap=nextrap+1
      CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
      extrapt(1)=t
    END IF
    !---Source terms for floor values
    neg_source=neg_source/2.d0
    DO j=1,mesh%nc
      IF(neg_source(1,j)<1.d-2)neg_source(1,j)=0.d0
      IF(neg_source(2,j)<1.d-2)neg_source(2,j)=0.d0
      IF(neg_source(3,j)<1.d-2)neg_source(3,j)=0.d0
      IF(neg_flag(1,j)>0.d0)neg_source(1,j)=MAX(neg_source(1,j),MIN(2.d0,neg_flag(1,j)))
      IF(neg_flag(2,j)>0.d0)neg_source(2,j)=MAX(neg_source(2,j),MIN(2.d0,neg_flag(2,j)))
      IF(neg_flag(3,j)>0.d0)neg_source(3,j)=MAX(neg_source(3,j),MIN(2.d0,neg_flag(3,j)))
    END DO
    k=COUNT(neg_source(1,:)>0.d0)
    IF(k>0)WRITE(*,*)oft_env%rank,'N-FLOOR cells  ',k
    k=COUNT(neg_source(2,:)>0.d0)
    IF(k>0)WRITE(*,*)oft_env%rank,'T-FLOOR cells  ',k
    IF(xmhd_two_temp)THEN
      k=COUNT(neg_source(3,:)>0.d0)
      IF(k>0)WRITE(*,*)oft_env%rank,'Te-FLOOR cells ',k
    END IF
    neg_flag=0.d0
!------------------------------------------------------------------------------
! Perform NL Advance
!------------------------------------------------------------------------------
    IF(i==2)xmhd_opcount=0
    IF(MOD(xmhd_opcount,xmhd_prefreq)==0.OR.ABS(xmhd_opdt-dt)/dt>.15d0)xmhd_opcount=0
    !---Update preconditioning matrices if necessary
    IF(xmhd_mfnk.AND.xmhd_opcount/=0.AND.(nksolver%lits>baseit+1.OR.nksolver%nlits>=nl_update))THEN
       IF(oft_env%head_proc)WRITE(*,*)'Update Jacobian',nksolver%lits,baseit,nredo
       IF(nredo<1)THEN
         jac_dt=dt/2.d0
         xmhd_skip_update=.FALSE.
         CALL xmhd_set_ops(u)
         xmhd_skip_update=.TRUE.
         baseit=-1
         nredo=nredo+1
       ELSE
         xmhd_opcount=0
       END IF
    END IF
    IF(xmhd_opcount==0)THEN
      jac_dt=dt/2.d0
      xmhd_skip_update=.FALSE.
      CALL xmhd_set_ops(u)
      xmhd_skip_update=.TRUE.
      IF(rst_ind==0.AND.i<3*xmhd_prefreq)force_refactor=.TRUE.
      CALL xmhd_pre%update(force_refactor)
      force_refactor=.FALSE.
      xmhd_opdt=dt
      xmhd_opcount=0
      baseit=-1
      nredo=0
    END IF
    xmhd_opcount=xmhd_opcount+1
    !---Retry loop for solver failure
    CALL up%add(0.d0,1.d0,u)
    DO j=1,4
      !---Compute metric to get RHS
      IF(xmhd_hall.OR.xmhd_upwind)NULLIFY(oft_xmhd_ops%A%up)
      oft_xmhd_ops%A%dt=-abs(dt/2.d0)
      CALL oft_xmhd_ops%A%apply(u,v)
      !---Get diagnostics
      mer=oft_xmhd_ops%A%diag_vals(1)/2.d0
      ver=oft_xmhd_ops%A%diag_vals(2)*m_ion/2.d0
      tflux=oft_xmhd_ops%A%diag_vals(3)
      tcurr=oft_xmhd_ops%A%diag_vals(4)
      npart=oft_xmhd_ops%A%diag_vals(5)/mesh_vol
      temp_avg=oft_xmhd_ops%A%diag_vals(6)/mesh_vol
      IF(xmhd_two_temp)tempe_avg=oft_xmhd_ops%A%diag_vals(7)/mesh_vol
      !---Extrapolate solution
      IF(maxextrap>0)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,t+dt,u)
      !---Apply driver BCs
      IF(PRESENT(driver))THEN
        CALL driver%apply(up,u,v,t,dt)
      END IF
      !---Solve Non-linear system
      IF(xmhd_hall.OR.xmhd_upwind)oft_xmhd_ops%A%up=>up
      oft_xmhd_ops%A%dt=abs(dt/2.d0)
      jac_dt=dt/2.d0
      CALL nksolver%apply(u,v)
      IF(nksolver%cits<0)THEN
        !---Apply new floor source
        DO k=1,mesh%nc
          IF(neg_flag(1,k)>0.d0)neg_source(1,k)=MAX(neg_source(1,k),MIN(2.d0,neg_flag(1,k)))
          IF(neg_flag(2,k)>0.d0)neg_source(2,k)=MAX(neg_source(2,k),MIN(2.d0,neg_flag(2,k)))
          IF(neg_flag(3,k)>0.d0)neg_source(3,k)=MAX(neg_source(3,k),MIN(2.d0,neg_flag(3,k)))
        END DO
        neg_flag=0.d0
        !
        CALL u%add(0.d0,1.d0,up)
        dt=dt/2.d0
        jac_dt=dt/2.d0
        xmhd_skip_update=.FALSE.
        CALL xmhd_set_ops(u)
        xmhd_skip_update=.TRUE.
        CALL xmhd_pre%update(force_refactor)
        force_refactor=.FALSE.
        xmhd_opdt=dt
        xmhd_opcount=0
        itcount=ittarget
        dthist=dt
      ELSE
        EXIT
      END IF
      IF(j==4)CALL oft_abort('NL Solver failed to converge.','xmhd_run',__FILE__)
    END DO
!------------------------------------------------------------------------------
! Compute energies and fluxes
!------------------------------------------------------------------------------
    !---Compute divergence error
    IF(xmhd_monitor_div)THEN
      CALL oft_xmhd_pop(u,sub_fields)
      Bfield%u=>sub_fields%B
      CALL Bfield%setup(oft_hcurl,oft_hgrad)
      !---Compute jump error
      jump_error=hcurl_grad_jump_error(xmhd_ML_hcurl_grad%current_level,sub_fields%B,oft_hcurl%quad%order)
      divfield%u=>sub_fields%B
      CALL divfield%setup(oft_hcurl,oft_hgrad)
      derror=scal_energy(mesh,divfield,oft_hcurl%quad%order)
    END IF
!------------------------------------------------------------------------------
! Write out initial solution progress
!------------------------------------------------------------------------------
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      hist_i4=(/rst_ind+i-1,nksolver%lits,nksolver%nlits/)
      hist_r4=REAL([t,tflux,tcurr,mer,jump_error,derror,ver,npart,temp_avg,tempe_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
      WRITE(*,103)rst_ind+i,t,nksolver%lits,nksolver%nlits,elapsed_time,dt
      IF(oft_debug_print(1))WRITE(*,*)
      CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
    END IF
    !---Process probes
    IF(PRESENT(probes))THEN
      CALL oft_xmhd_pop(up,initial_fields)
      CALL initial_fields%V%scale(vel_scale)
      CALL initial_fields%Ne%scale(den_scale)
      CALL probes%apply(initial_fields,t)
    END IF
!------------------------------------------------------------------------------
! Update time step and current time
!------------------------------------------------------------------------------
    IF(baseit<0)baseit=nksolver%lits
    itcount(MOD(i,XMHD_ITCACHE)+1)=nksolver%lits
    dthist(MOD(i,XMHD_ITCACHE)+1)=dt
    t=t+dt
    dt=ittarget*SUM(dthist)/SUM(itcount)
    IF(dt>dtin)dt=dtin
    IF(dt<1.d-10)CALL oft_abort('Time step dropped too low!','xmhd_run',__FILE__)
!------------------------------------------------------------------------------
! Write dump file
!------------------------------------------------------------------------------
    IF(MOD(i,rst_freq)==0)THEN
      IF(oft_env%head_proc)CALL mytimer%tick
      !---Create restart file
      WRITE(rst_char,104)rst_ind+i
      READ(rst_char,104,IOSTAT=io_stat)rst_tmp
      IF((io_stat/=0).OR.(rst_tmp/=rst_ind+i))CALL oft_abort("Step count exceeds format width", "xmhd_run", __FILE__)
      CALL oft_xmhd_rst_save(u, t, dt, 'xMHD_'//rst_char//'.rst', 'U')
      IF(PRESENT(driver))CALL driver%rst_save('xMHD_'//rst_char//'.rst')
      IF(oft_env%head_proc)THEN
        elapsed_time=mytimer%tock()
        WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
        CALL hist_file%flush
        IF(PRESENT(probes))CALL probes%flush
      END IF
    END IF
END DO
IF(oft_env%head_proc)CALL hist_file%close
CALL oft_xmhd_pop(u,initial_fields)
CALL initial_fields%V%scale(vel_scale)
CALL initial_fields%Ne%scale(den_scale)
!---Clear up memory
CALL nksolver%delete
CALL solver%delete
CALL divout%delete
!--- solver fields
CALL u%delete
CALL up%delete
CALL v%delete
DEALLOCATE(u,up,v)
!--- plotting fields
CALL sub_fields%B%delete
DEALLOCATE(sub_fields%B)
!--- interpolation objects
CALL Bfield%delete
CALL divfield%delete
!--- extrapolation fields
IF(maxextrap>0)THEN
  DO i=1,maxextrap
    CALL extrap_fields(i)%f%delete
    DEALLOCATE(extrap_fields(i)%f)
  END DO
  DEALLOCATE(extrap_fields,extrapt)
END IF
DEBUG_STACK_POP
end subroutine xmhd_run
!------------------------------------------------------------------------------
!> Main driver subroutine for extended MHD time advance
!!
!! Runtime options are set in the main input file using the group `xmhd_options` group,
!! see \ref xmhd_run.
!!
!! @note This method assumes that [B0, V0, Ne0, Ti0, Te0] constitute an
!! equilibrium.
!------------------------------------------------------------------------------
subroutine xmhd_lin_run(equil_fields,pert_fields,escale)
TYPE(xmhd_sub_fields), INTENT(inout) :: equil_fields !< Equilibrium fields
TYPE(xmhd_sub_fields), INTENT(inout) :: pert_fields !< Perurbed fields
REAL(r8), OPTIONAL, INTENT(in) :: escale !< Desired total energy to be used for rescaling (optional)
!---H(Curl) + Grad(H^1) divout solver
TYPE(oft_hcurl_grad_divout) :: divout
!---Jacobian solver
type(oft_native_gmres_solver), target :: solver
integer(i4) :: nlevels
integer(i4), allocatable, dimension(:) :: levels
!---Local variables
class(oft_vector), pointer :: du,u0,v,un1,un2
type(xmhd_sub_fields) :: sub_fields
type(oft_hcurl_grad_rinterp) :: Bfield
type(oft_hcurl_grad_dinterp) :: divfield
type(oft_xmhd_massmatrix) :: mop
type(oft_timer) :: mytimer
integer(i4) :: i,j,ierr,io_unit,io_stat,rst_tmp,nl_update
integer(i4) :: itcount(XMHD_ITCACHE)
real(r8) :: dt,dthist(XMHD_ITCACHE)
real(r8) :: mag,div_err,mer,merp,ver,verp,gerr,cerr,verr,elapsed_time
real(r8) :: fac,lramp,tflux,tcurr,t,dtin,div_error,jump_error,derror,de_scale
real(r8) :: ndens,npart,temp_avg,tempe_avg,mesh_vol
real(r8), pointer, dimension(:) :: vals => NULL()
character(LEN=XMHD_RST_LEN) :: rst_char
character(LEN=OFT_SLEN) :: comm_line
!---Extrapolation fields
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(11)
logical :: rst
!---Input variables
real(r8) :: lin_tol,nl_tol
integer(i4) :: rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(oft_hcurl,xmhd_ML_hcurl%current_level))CALL oft_abort("Invalid Curl FE object","xmhd_lin_run",__FILE__)
IF(.NOT.oft_3D_lagrange_cast(oft_lagrange,xmhd_ML_lagrange%current_level))CALL oft_abort("Invalid Lagrange FE object","xmhd_lin_run",__FILE__)
IF(.NOT.oft_3D_h1_cast(oft_hgrad,xmhd_ML_H1grad%current_level))CALL oft_abort("Invalid Grad FE object","xmhd_lin_run",__FILE__)
mg_mesh=>xmhd_ML_hcurl%ml_mesh
mesh=>oft_hcurl%mesh
!------------------------------------------------------------------------------
! Read-in Parameters
!------------------------------------------------------------------------------
IF(ASSOCIATED(equil_fields%Te).AND.ASSOCIATED(pert_fields%Te))xmhd_two_temp=.TRUE.
IF(XOR(ASSOCIATED(equil_fields%Te),ASSOCIATED(pert_fields%Te)))CALL oft_abort( &
  "Te0 and dTe ICs are required for two temp.", "xmhd_lin_run", __FILE__)
CALL xmhd_read_settings(dt,lin_tol,nl_tol,rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget,nl_update)
IF(den_scale<0.d0)den_scale=SQRT(equil_fields%Ne%dot(equil_fields%Ne)/REAL(equil_fields%Ne%ng,8))
IF((d2_dens>0.d0).AND.(n2_scale<0.d0))n2_scale=den_scale*(REAL(oft_lagrange%order,8)/mesh%hrms)**2
IF((eta_hyper>0.d0).AND.(j2_scale<0.d0))j2_scale=(REAL(oft_lagrange%order,8)/mesh%hrms)**2
!---Force heating and conduction terms to zero
xmhd_ohmic=.FALSE.
xmhd_visc_heat=.FALSE.
xmhd_brag=.FALSE.
xmhd_linear=.TRUE.
xmhd_diss_centered=.TRUE.
!------------------------------------------------------------------------------
! Setup ML environment
!------------------------------------------------------------------------------
CALL xmhd_setup_rep
!------------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------------
call oft_xmhd_create(du)
call oft_xmhd_create(u0)
call oft_xmhd_create(v)
CALL oft_xmhd_create(un1)
CALL oft_xmhd_create(un2)
!------------------------------------------------------------------------------
! Create extrapolation fields
!------------------------------------------------------------------------------
IF(maxextrap>0)THEN
  ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
  DO i=1,maxextrap
    CALL oft_xmhd_create(extrap_fields(i)%f)
    extrapt(i)=0.d0
  END DO
  nextrap=0
END IF
!------------------------------------------------------------------------------
! Create divergence cleaner
!------------------------------------------------------------------------------
call xmhd_ML_hcurl_grad%vec_create(sub_fields%B)
CALL divout%setup(xmhd_ML_hcurl_grad,"grnd")
divout%pm=.TRUE.
IF(TRIM(bbc)=="ic")divout%keep_boundary=.TRUE.
!------------------------------------------------------------------------------
! Setup history file
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd.hist', desc="History file for linear xMHD run")
104 FORMAT('E0 = ',ES11.3)
  IF(PRESENT(escale))THEN
    WRITE(comm_line,104)escale
    CALL hist_file%add_comm(comm_line)
  END IF
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('tflux','r4', desc="Toroidal flux [Wb]")
  CALL hist_file%add_field('tcurr','r4', desc="Toroidal current [A*mu0]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('jerr', 'r4', desc="Magnetic jump error")
  CALL hist_file%add_field('derr', 'r4', desc="Magnetic div error")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average ion temperature [eV]")
  CALL hist_file%add_field('te',   'r4', desc="Average electron temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
END IF
!------------------------------------------------------------------------------
! Import initial fields
!------------------------------------------------------------------------------
dtin=dt
105 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
WRITE(rst_char,105)rst_ind
READ(rst_char,105,IOSTAT=io_stat)rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=rst_ind))CALL oft_abort("Step count exceeds format width", "xmhd_lin_run", __FILE__)
rst=oft_file_exist(TRIM('xMHD_'//rst_char//'.rst'))
IF(rst.AND.(rst_ind/=0))THEN
  CALL oft_xmhd_rst_load(du, 'xMHD_'//rst_char//'.rst', 'U', t, dt)
  WRITE(rst_char,105)0
  rst=oft_file_exist(TRIM('xMHD_'//rst_char//'.rst'))
  IF(.NOT.rst)CALL oft_abort('Equilibrium restart does not exist.','xmhd_lin_run',__FILE__)
  CALL oft_xmhd_rst_load(u0, 'xMHD_'//rst_char//'.rst', 'U0')
  dt=dtin
ELSE
  IF(rst_ind>0)CALL oft_abort('Restart does not exist.','xmhd_run',__FILE__)
  CALL equil_fields%V%scale(1.d0/vel_scale)
  CALL equil_fields%Ne%scale(1.d0/den_scale)
  CALL oft_xmhd_push(equil_fields,u0)
  CALL equil_fields%V%scale(vel_scale)
  CALL equil_fields%Ne%scale(den_scale)
  !---
  CALL pert_fields%V%scale(1.d0/vel_scale)
  CALL pert_fields%Ne%scale(1.d0/den_scale)
  CALL oft_xmhd_push(pert_fields,du)
  t=0.d0
  !---Create initial conditions restart file
  WRITE(rst_char,105)0
  CALL oft_xmhd_rst_save(du, t, dt, 'xMHD_'//rst_char//'.rst', 'U')
  CALL xmhd_rep%vec_save(u0,'xMHD_'//rst_char//'.rst','U0', append=.TRUE.)
  !---Setup history file
  IF(oft_env%head_proc)CALL hist_file%write_header
END IF
IF(oft_env%head_proc)CALL hist_file%open ! Open history file
!------------------------------------------------------------------------------
! Setup Operators
!------------------------------------------------------------------------------
call xmhd_alloc_ops
jac_dt=dt
call xmhd_set_ops(u0)
mop%u0=>u0
!---Force perturbation to zero in solid regions
IF(xmhd_rw)THEN
  NULLIFY(vals)
  DO j=1,3
    CALL pert_fields%V%get_local(vals,j)
    DO i=1,oft_lagrange%ne
      IF(oft_xmhd_ops%solid_node(i))vals(i)=0.d0
    END DO
    CALL pert_fields%V%restore_local(vals,j)
  END DO
  DEALLOCATE(vals)
END IF
!------------------------------------------------------------------------------
! Compute energies
!------------------------------------------------------------------------------
mesh_vol=SUM(mesh%cv)
mesh_vol=oft_mpi_sum(mesh_vol)
CALL mop%apply(du,v)
!---Get diagnostics
mer=mop%diag_vals(1)/2.d0
ver=mop%diag_vals(2)*m_ion/2.d0
tflux=mop%diag_vals(3)
tcurr=mop%diag_vals(4)
npart=mop%diag_vals(5)/mesh_vol
temp_avg=mop%diag_vals(6)/mesh_vol
tempe_avg=mop%diag_vals(7)/mesh_vol
!---Compute divergence error
IF(xmhd_monitor_div)THEN
  CALL oft_xmhd_pop(du,sub_fields)
  Bfield%u=>sub_fields%B
  CALL Bfield%setup(oft_hcurl,oft_hgrad)
  !---Compute jump error
  jump_error=hcurl_grad_jump_error(xmhd_ML_hcurl_grad%current_level,sub_fields%B,oft_hcurl%quad%order)
  divfield%u=>sub_fields%B
  CALL divfield%setup(oft_hcurl,oft_hgrad)
  derror=scal_energy(mesh,divfield,oft_hcurl%quad%order)
ELSE
  jump_error=-1.d0
  derror=-1.d0
END IF
!------------------------------------------------------------------------------
! Print run information
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
100 FORMAT(2X,A,L1)
101 FORMAT(2X,2A)
102 FORMAT(2X,A,ES11.3)
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting linear xMHD run'
  WRITE(*,100)'Lorentz    = ',xmhd_jcb
  WRITE(*,100)'V-Advec    = ',xmhd_advec
  WRITE(*,100)'N-Adv      = ',xmhd_adv_den
  WRITE(*,100)'T-Adv      = ',xmhd_adv_temp
  WRITE(*,100)'Hall       = ',xmhd_hall
  WRITE(*,100)'Ohmic      = ',xmhd_ohmic
  WRITE(*,100)'Visc Heat  = ',xmhd_visc_heat
  WRITE(*,101)'B-BC       = ',bbc
  WRITE(*,101)'V-BC       = ',vbc
  IF(xmhd_adv_den)WRITE(*,101)'N-BC       = ',nbc
  IF(xmhd_adv_temp)WRITE(*,101)'T-BC       = ',tbc
  WRITE(*,102)'t0         = ',t
  WRITE(*,102)'dt         = ',dt
  WRITE(*,102)'eta0       = ',eta
  WRITE(*,102)'eta0_Temp  = ',eta_temp
  WRITE(*,101)'visc_type  = ',visc_type
  WRITE(*,102)'  nu_par   = ',nu_par
  WRITE(*,102)'  nu_perp  = ',nu_perp
  IF(xmhd_adv_den)WRITE(*,102)'D_dens     = ',d_dens
  IF(xmhd_adv_temp)THEN
    WRITE(*,102)'kappa_par  = ',kappa_par
    WRITE(*,102)'kappa_perp = ',kappa_perp
  END IF
  WRITE(*,102)'L-Tol      = ',lin_tol
  WRITE(*,102)'NL-Tol     = ',nl_tol
  WRITE(*,102)'Mag Energy = ',mer
  WRITE(*,102)'Kin Energy = ',ver
  IF(xmhd_monitor_div)WRITE(*,'(2X,A,2ES11.3)')'Div Error  = ',derror,jump_error
  WRITE(*,102)'# B-DOFs   = ',REAL(oft_hcurl%global%ne+oft_hgrad%global%ne,4)
  WRITE(*,102)'# V-DOFs   = ',REAL(3*oft_lagrange%global%ne,4)
  WRITE(*,102)'# N-DOFs   = ',REAL(oft_lagrange%global%ne,4)
  WRITE(*,102)'# T-DOFs   = ',REAL(oft_lagrange%global%ne,4)
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
END IF
!------------------------------------------------------------------------------
! Setup linear solver
!------------------------------------------------------------------------------
solver%A=>oft_xmhd_ops%J
solver%its=400
solver%atol=lin_tol
solver%itplot=1
solver%nrits=40
solver%pm=oft_env%pm
!------------------------------------------------------------------------------
! Setup Preconditioner
!------------------------------------------------------------------------------
nlevels=xmhd_nlevels-xmhd_minlev+1
NULLIFY(solver%pre)
IF(nlevels==1)THEN
  IF(ASSOCIATED(xmhd_pre_node))THEN
    CALL create_solver_xml(solver%pre,xmhd_pre_node)
  ELSE
    CALL create_diag_pre(solver%pre)
  END IF
ELSE
  ALLOCATE(levels(nlevels))
  levels=(/(i,i=xmhd_minlev,xmhd_nlevels)/)
  xmhd_ml_vecspace%ML_FE_rep=>ML_xmhd_rep
  CALL create_mlpre(solver%pre,ml_J,levels,nlevels=nlevels, &
    ml_vecspace=xmhd_ml_vecspace, &
    stype=2,nu=nu_xmhd(xmhd_minlev:xmhd_nlevels),xml_root=xmhd_pre_node)
  DEALLOCATE(levels)
END IF
xmhd_pre=>solver%pre
!------------------------------------------------------------------------------
! Begin time stepping
!------------------------------------------------------------------------------
DO i=1,nsteps
    merp=mer
    verp=ver
    !---Clean accumulated divergence error
    IF(MOD(rst_ind+i,nclean)==0)THEN
      CALL oft_xmhd_pop(du,sub_fields)
      IF(oft_env%head_proc)WRITE(*,*)'  Cleaning Divergence'
      CALL divout%apply(sub_fields%B)
      CALL oft_xmhd_push(sub_fields,du)
    END IF
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Update extrapolation fields
    IF(maxextrap>0)THEN
      DO j=maxextrap,2,-1
       CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
       extrapt(j)=extrapt(j-1)
      END DO
      IF(nextrap<maxextrap)nextrap=nextrap+1
      CALL extrap_fields(1)%f%add(0.d0,1.d0,du)
      extrapt(1)=t
    END IF
!------------------------------------------------------------------------------
! Perform Linear Advance
!------------------------------------------------------------------------------
    !---Compute metric to get RHS
    IF(i==2)THEN
      jac_dt=dt*2.d0/3.d0
      dt_scale=2.d0/3.d0
      CALL xmhd_set_ops(u0)
      CALL xmhd_pre%update(.TRUE.)
    END IF
    !---Compute metric to get RHS
    CALL un2%add(0.d0,1.d0,un1)
    CALL un1%add(0.d0,1.d0,du)
    IF(i>1)THEN
      CALL un2%add(-1.d0,4.d0,un1)
      CALL un2%scale(1.d0/3.d0)
      dt_scale=2.d0/3.d0
      CALL mop%apply(un2,v)
      CALL xmhd_diag(du,mop%diag_vals,u0)
    ELSE
      CALL mop%apply(du,v)
    END IF
    ! CALL mop%apply(du,v)
    !---Get diagnostics
    mer=mop%diag_vals(1)/2.d0
    ver=mop%diag_vals(2)*m_ion/2.d0
    tflux=mop%diag_vals(3)
    tcurr=mop%diag_vals(4)
    npart=mop%diag_vals(5)/mesh_vol
    temp_avg=mop%diag_vals(6)/mesh_vol
    tempe_avg=mop%diag_vals(7)/mesh_vol
    IF(maxextrap>0)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,t+dt,du)
!------------------------------------------------------------------------------
! Scale fields
!------------------------------------------------------------------------------
    IF(PRESENT(escale))THEN
      de_scale=SQRT(escale/(ver+mer/mu0))
      CALL v%scale(de_scale)
      CALL un1%scale(de_scale)
    END IF
    CALL solver%apply(du,v)
!------------------------------------------------------------------------------
! Compute energies and fluxes
!------------------------------------------------------------------------------
    !---Compute divergence error
    IF(xmhd_monitor_div)THEN
      CALL oft_xmhd_pop(du,sub_fields)
      Bfield%u=>sub_fields%B
      CALL Bfield%setup(oft_hcurl,oft_hgrad)
      !---Compute jump error
      jump_error=hcurl_grad_jump_error(xmhd_ML_hcurl_grad%current_level,sub_fields%B,oft_hcurl%quad%order)
      divfield%u=>sub_fields%B
      CALL divfield%setup(oft_hcurl,oft_hgrad)
      derror=scal_energy(mesh,divfield,oft_hcurl%quad%order)
    END IF
!------------------------------------------------------------------------------
! Write out initial solution progress
!------------------------------------------------------------------------------
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      hist_i4=(/rst_ind+i-1,solver%cits,0/)
      hist_r4=REAL([t,tflux,tcurr,mer,jump_error,derror,ver,npart,temp_avg,tempe_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,F12.3,ES12.2)
      WRITE(*,103)rst_ind+i,t,solver%cits,elapsed_time,dt
      IF(oft_debug_print(1))WRITE(*,*)
      CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
    END IF
    !---Update current time
    t=t+dt
!------------------------------------------------------------------------------
! Write dump file
!------------------------------------------------------------------------------
    IF(MOD(i,rst_freq)==0)THEN
      IF(oft_env%head_proc)CALL mytimer%tick
      !---Create restart file
      WRITE(rst_char,105)rst_ind+i
      READ(rst_char,105,IOSTAT=io_stat)rst_tmp
      IF((io_stat/=0).OR.(rst_tmp/=rst_ind+i))CALL oft_abort("Step count exceeds format width", "xmhd_lin_run", __FILE__)
      CALL oft_xmhd_rst_save(du, t, dt, 'xMHD_'//rst_char//'.rst', 'U')
      IF(oft_env%head_proc)THEN
        elapsed_time=mytimer%tock()
        WRITE(*,'(A,F12.3)')'  I/O Time = ',elapsed_time
        CALL hist_file%flush
      END IF
    END IF
END DO
IF(oft_env%head_proc)CALL hist_file%close
CALL oft_xmhd_pop(du,pert_fields)
CALL pert_fields%V%scale(vel_scale)
CALL pert_fields%Ne%scale(den_scale)
!---Clear up memory
CALL solver%delete
CALL divout%delete
!--- solver fields
CALL du%delete
CALL u0%delete
CALL v%delete
DEALLOCATE(du,u0,v)
CALL un1%delete
CALL un2%delete
DEALLOCATE(un1,un2)
!--- plotting fields
CALL sub_fields%B%delete
DEALLOCATE(sub_fields%B)
!--- interpolation objects
CALL Bfield%delete
CALL divfield%delete
!--- extrapolation fields
IF(maxextrap>0)THEN
  DO i=1,maxextrap
    CALL extrap_fields(i)%f%delete
    DEALLOCATE(extrap_fields(i)%f)
  END DO
  DEALLOCATE(extrap_fields,extrapt)
END IF
DEBUG_STACK_POP
end subroutine xmhd_lin_run
!------------------------------------------------------------------------------
!> Evalute upwinding parameter
!!
!! Upwinding parameter for inconsistent streamline upwinded Petrov-Galerkin method
!! (distorted bases applied to advective terms only).
!!
!! @result Upwinding weight
!------------------------------------------------------------------------------
PURE FUNCTION xmhd_upwind_weight(v_mag,he,nu) RESULT(w_up)
REAL(r8), INTENT(in) :: v_mag !< Advection velocity [m/s]
REAL(r8), INTENT(in) :: he !< Element size [m]
REAL(r8), INTENT(in) :: nu !< Dissipation amplitude [m^2/s^2]
REAL(r8) :: Pc,beta,w_up
Pc = MAX(1.d-8,v_mag*he/(2.d0*nu)) ! Peclet number
beta = 1.d0/TANH(Pc) - 1.d0/Pc ! Beta parameter
w_up = beta*he/(2.d0*MAX(1.d-10,v_mag)) ! Upwinding weight
END FUNCTION xmhd_upwind_weight
!------------------------------------------------------------------------------
!> Update Jacobian matrix with new fields
!------------------------------------------------------------------------------
subroutine xmhd_build_ops(fullin)
class(fem_interp), intent(inout) :: fullin !< Interpolator to evaluate fields
type(oft_local_mat), allocatable, dimension(:,:) :: mtmp
type(oft_1d_int), allocatable, dimension(:) :: iloc
integer(i4) :: ii,i,ip,j,m,n,jr,jc,jp,nrows,jtmp(1)
integer(i4), allocatable, target :: j_lag(:),j_hgrad(:),j_hcurl(:)
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
real(r8) :: f,diff,chi_temp(2),bmag,det,uvec(3),umat(3,3),nu_tmp(2),ten1(3,3)
real(r8), allocatable, target :: lag_rop(:),lag_gop(:,:),hgrad_rop(:,:)
real(r8), allocatable, target :: hcurl_rop(:,:),hcurl_cop(:,:)
real(r8), allocatable :: lag_gopt(:,:),hgrad_ropt(:,:)
real(r8), allocatable :: hcurl_ropt(:,:),hcurl_copt(:,:)
class(oft_vector), pointer :: tmp
real(r8) :: goptmp(3,4),vol,b(1,1),eta_curr,u,c1,c2,c3,c4,divv,tgam_fac,fulltmp(40)
real(r8) :: bhat(3),dvT(3,3)
real(r8) :: he,v0_mag,cgop(3,6),uten(3,3,3)
logical :: curved
class(oft_matrix), pointer :: Jac
type(oft_quad_type), pointer :: quad
type(xmhd_loc_values) :: u0
DEBUG_STACK_PUSH
IF(oft_debug_print(2))write(*,'(4X,A,I4)')'Building xMHD Jacobian: ',xmhd_level
!---
Jac=>oft_xmhd_ops%J
quad=>xmhd_ML_hcurl%levels(xmhd_ML_hcurl%nlevels)%fe%quad !oft_hcurl%quad
tgam_fac=temp_gamma - 1.d0
!------------------------------------------------------------------------------
! Reset matrices
!------------------------------------------------------------------------------
CALL Jac%zero
!--Setup thread locks
ALLOCATE(tlocks(xmhd_rep%nfields))
DO i=1,xmhd_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
!$omp parallel private(det,j_lag,j_hgrad,j_hcurl,lag_rop,lag_gop,hgrad_rop,hcurl_rop, &
!$omp hcurl_cop,curved,goptmp,vol,j,m,jr,jc,bhat,he,v0_mag,dvT,eta_curr,u,uvec,umat, &
!$omp mtmp,chi_temp,bmag,c1,c2,c3,c4,divv,nu_tmp,ten1,cgop,fulltmp,lag_gopt, &
!$omp hgrad_ropt,hcurl_ropt,hcurl_copt,uten,i,ii,iloc,u0)
!---Thread local arrays
allocate(j_lag(oft_lagrange%nce),j_hgrad(oft_hgrad%nce),j_hcurl(oft_hcurl%nce))
allocate(lag_rop(oft_lagrange%nce),lag_gop(3,oft_lagrange%nce))
allocate(hgrad_rop(3,oft_hgrad%nce))
allocate(hcurl_rop(3,oft_hcurl%nce),hcurl_cop(3,oft_hcurl%nce))
allocate(lag_gopt(oft_lagrange%nce,3),hgrad_ropt(oft_hgrad%nce,3))
allocate(hcurl_ropt(oft_hcurl%nce,3),hcurl_copt(oft_hcurl%nce,3))
allocate(mtmp(xmhd_rep%nfields,xmhd_rep%nfields))
CALL xmhd_rep%mat_setup_local(mtmp, xmhd_mat_mask)
allocate(iloc(xmhd_rep%nfields))
iloc(1)%v=>j_hcurl
iloc(2)%v=>j_hgrad
DO i=3,7
  iloc(i)%v=>j_lag
END DO
IF(xmhd_two_temp)iloc(te_ind)%v=>j_lag
IF(n2_ind>0)iloc(n2_ind)%v=>j_lag
IF(j2_ind>0)iloc(j2_ind)%v=>j_hcurl
!---Setup caching
IF(xmhd_level==xmhd_nlevels)THEN
  xmhd_lag_rop=>lag_rop
  xmhd_lag_gop=>lag_gop
  xmhd_hgrad_rop=>hgrad_rop
  xmhd_hcurl_rop=>hcurl_rop
  xmhd_hcurl_cop=>hcurl_cop
END IF
!--- Define viscosity temporary varible (currently static)
nu_tmp(1) = nu_par
nu_tmp(2) = nu_perp
!------------------------------------------------------------------------------
! Integration loop
!------------------------------------------------------------------------------
!$omp do schedule(static)
DO ip=1,mesh%nparts
do ii=1,mesh%tloc_c(ip)%n
  i=mesh%tloc_c(ip)%v(ii)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  he = (6.d0*SQRT(2.d0)*mesh%cv(i))**(1.d0/3.d0)/REAL(oft_lagrange%order,8) ! Cell node spacing
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j_lag)
  call oft_hgrad%ncdofs(i,j_hgrad)
  call oft_hcurl%ncdofs(i,j_hcurl)
  !---Zero local matrix
  CALL xmhd_rep%mat_zero_local(mtmp)
!------------------------------------------------------------------------------
! Quadrature Loop
!------------------------------------------------------------------------------
  IF(xmhd_hall.AND.xmhd_two_temp)mtmp(1,7)%m=>mtmp(1,te_ind)%m ! Point hall-term coupling to electrons
  DO m=1,quad%np
!------------------------------------------------------------------------------
! Get local reconstructed operators
!------------------------------------------------------------------------------
    if(curved.OR.m==1)then
      call mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*quad%wts(m)
    !---Lagrange basis functions
    CALL oft_lag_eval_all(oft_lagrange,i,quad%pts(:,m),lag_rop)
    CALL oft_lag_geval_all(oft_lagrange,i,quad%pts(:,m),lag_gop,goptmp)
    lag_gopt=TRANSPOSE(lag_gop)
    CALL oft_h1_geval_all(oft_hgrad,i,quad%pts(:,m),hgrad_rop,goptmp)
    hgrad_ropt=TRANSPOSE(hgrad_rop)
    CALL oft_hcurl_eval_all(oft_hcurl,i,quad%pts(:,m),hcurl_rop,goptmp)
    hcurl_ropt=TRANSPOSE(hcurl_rop)
    CALL oft_hcurl_ceval_all(oft_hcurl,i,quad%pts(:,m),hcurl_cop,cgop)
    hcurl_copt=TRANSPOSE(hcurl_cop)
!------------------------------------------------------------------------------
! Reconstruct fields
!------------------------------------------------------------------------------
    call fullin%interp(i,quad%pts(:,m),goptmp,fulltmp)
    CALL xmhd_interp_unpack(fulltmp,u0)
    v0_mag = magnitude(u0%V)
    dvT = TRANSPOSE(u0%dV)
    divv = u0%dV(1,1) + u0%dV(2,2) + u0%dV(3,3)
    bmag = SQRT(SUM(u0%B**2)) + xmhd_eps
    bhat = u0%B/bmag
    !---Handle density and temperature floors
    IF(temp_floor>0.d0)u0%Ti=MAX(temp_floor,u0%Ti)
    IF(temp_floor>0.d0)u0%Te=MAX(temp_floor,u0%Te)
    IF(den_floor>0.d0)u0%N=MAX(den_floor,u0%N)
    !---Transport coefficients
    eta_curr=eta*eta_reg(mesh%reg(i))
    IF(eta_temp>0.d0)THEN
      IF(ASSOCIATED(res_profile))THEN
        eta_curr=eta*res_profile(u0%Te)
      ELSE
        eta_curr=eta*(eta_temp/u0%Te)**(3.d0/2.d0)
      END IF
    END IF
    IF(xmhd_brag)THEN
      IF(xmhd_two_temp)THEN
        chi_temp=brag_ion_transport(u0%N,u0%Ti,mu_ion,bmag)/u0%N
      ELSE
        chi_temp=brag_comb_transport(u0%N,u0%Ti,mu_ion,bmag)/u0%N
      END IF
      chi_temp(1)=kappa_par*chi_temp(1)
      chi_temp(2)=kappa_perp*chi_temp(2)
    ELSE
      chi_temp(1)=kappa_par*den_scale/u0%N
      chi_temp(2)=kappa_perp*den_scale/u0%N
    END IF
    !------------------------------------------------------------------------------
    ! Handle solid region
    !------------------------------------------------------------------------------
    IF(solid_cell(i))THEN
      eta_curr=eta_reg(mesh%reg(i))
      do jc=1,oft_hcurl%nce
        uvec = eta_curr*hcurl_cop(:,jc)*jac_dt*det
        !$omp simd
        do jr=1,oft_hcurl%nce
          mtmp(1,1)%m(jr,jc) = mtmp(1,1)%m(jr,jc) &
            + (hcurl_ropt(jr,1)*hcurl_rop(1,jc) &
            +  hcurl_ropt(jr,2)*hcurl_rop(2,jc) &
            +  hcurl_ropt(jr,3)*hcurl_rop(3,jc))*det &
            + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
        end do
      end do
      do jc=1,oft_hgrad%nce
        uvec=hgrad_rop(:,jc)*det
        !$omp simd
        do jr=1,oft_hcurl%nce
          mtmp(1,2)%m(jr,jc) = mtmp(1,2)%m(jr,jc) &
            + DOT_PRODUCT(hcurl_ropt(jr,:),uvec)
        end do
      end do
      do jc=1,oft_hcurl%nce
        uvec=hcurl_rop(:,jc)*det
        !$omp simd
        do jr=1,oft_hgrad%nce
          mtmp(2,1)%m(jr,jc) = mtmp(2,1)%m(jr,jc) &
            + DOT_PRODUCT(hgrad_ropt(jr,:),uvec)
        end do
      end do
      do jc=1,oft_hgrad%nce
        uvec=hgrad_rop(:,jc)*det
        !$omp simd
        do jr=1,oft_hgrad%nce
          mtmp(2,2)%m(jr,jc) = mtmp(2,2)%m(jr,jc) &
            + DOT_PRODUCT(hgrad_ropt(jr,:),uvec)
        end do
      end do
      CYCLE
    END IF
!------------------------------------------------------------------------------
! Compute dB_c/dt matrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Compute dB_c/dt(B_c) matrix (1,1)
!------------------------------------------------------------------------------
    IF(xmhd_adv_b)THEN
      c1 = jac_dt*det
      c3 = eta_curr*c1
      umat(:,1) = u0%V*c1
      IF(xmhd_hall)THEN
        c2 = c1/(mu0*u0%N*elec_charge)
        c3 = c3 + dt_scale*me_factor*elec_mass*det/(mu0*u0%N*elec_charge**2)
        umat(:,1) = umat(:,1) - u0%J*c2
      END IF
      do jc=1,oft_hcurl%nce
        uvec = c3*hcurl_cop(:,jc) - cross_product(umat(:,1),hcurl_rop(:,jc))
        IF(xmhd_hall)uvec = uvec + cross_product(hcurl_cop(:,jc),u0%B)*c2
        !$omp simd
        do jr=1,oft_hcurl%nce
          mtmp(1,1)%m(jr,jc) = mtmp(1,1)%m(jr,jc) &
            + (hcurl_ropt(jr,1)*hcurl_rop(1,jc) &
            +  hcurl_ropt(jr,2)*hcurl_rop(2,jc) &
            +  hcurl_ropt(jr,3)*hcurl_rop(3,jc))*det &
            + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
        end do
      end do
      !---Hyper resistivity
      IF(j2_ind>0)THEN
        c1 = -eta_hyper*jac_dt*det*j2_scale
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0*c1
        do jc=1,oft_hcurl%nce
          uvec = hcurl_cop(:,jc)*c1
          !$omp simd
          do jr=1,oft_hcurl%nce
            mtmp(1,j2_ind)%m(jr,jc) = mtmp(1,j2_ind)%m(jr,jc) &
            + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
          end do
          uvec = hcurl_cop(:,jc)*det/j2_scale
          !$omp simd
          do jr=1,oft_hcurl%nce
            mtmp(j2_ind,1)%m(jr,jc) = mtmp(j2_ind,1)%m(jr,jc) &
            + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
          end do
          uvec = hcurl_rop(:,jc)*det
          !$omp simd
          do jr=1,oft_hcurl%nce
            mtmp(j2_ind,j2_ind)%m(jr,jc) = mtmp(j2_ind,j2_ind)%m(jr,jc) &
            + DOT_PRODUCT(hcurl_ropt(jr,:),uvec)
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dB_c/dt(B_g) matrix (1,2)
!------------------------------------------------------------------------------
      do jc=1,oft_hgrad%nce
        uvec = -cross_product(umat(:,1),hgrad_rop(:,jc))
        !$omp simd
        do jr=1,oft_hcurl%nce
          mtmp(1,2)%m(jr,jc) = mtmp(1,2)%m(jr,jc) &
            + (hcurl_ropt(jr,1)*hgrad_rop(1,jc) &
            +  hcurl_ropt(jr,2)*hgrad_rop(2,jc) &
            +  hcurl_ropt(jr,3)*hgrad_rop(3,jc))*det &
            + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dB_c/dt(V) matrix (1,3:5)
!------------------------------------------------------------------------------
      c1 = jac_dt*det*vel_scale
      do jc=1,oft_lagrange%nce
        uvec = u0%B*lag_rop(jc)*c1
        !$omp simd
        do jr=1,oft_hcurl%nce
          mtmp(1,3)%m(jr,jc) = mtmp(1,3)%m(jr,jc) &
            + (hcurl_copt(jr,2)*uvec(3) - hcurl_copt(jr,3)*uvec(2))
          mtmp(1,4)%m(jr,jc) = mtmp(1,4)%m(jr,jc) &
            + (hcurl_copt(jr,3)*uvec(1) - hcurl_copt(jr,1)*uvec(3))
          mtmp(1,5)%m(jr,jc) = mtmp(1,5)%m(jr,jc) &
            + (hcurl_copt(jr,1)*uvec(2) - hcurl_copt(jr,2)*uvec(1))
        end do
      end do
!------------------------------------------------------------------------------
! Compute dB_c/dt(n) matrix (1,6)
!------------------------------------------------------------------------------
      IF(xmhd_hall.AND.xmhd_adv_den)THEN
        c1 = -jac_dt*det*den_scale/(mu0*elec_charge*u0%N*u0%N)
        c2 = -k_boltz*u0%Te*jac_dt*det*den_scale/(elec_charge*u0%N)
        c3 = c2/u0%N
        IF(xmhd_linear)umat(:,1)=k_boltz*u0%dTe*jac_dt*det*den_scale/(elec_charge*u0%N)
        do jc=1,oft_lagrange%nce
          uvec = (cross_product(u0%J,u0%B)*c1 + u0%dN*c3)*lag_rop(jc) + lag_gop(:,jc)*c2
          IF(xmhd_linear)uvec = uvec + umat(:,1)*lag_rop(jc)
          !$omp simd
          do jr=1,oft_hcurl%nce
            mtmp(1,6)%m(jr,jc) = mtmp(1,6)%m(jr,jc) &
              + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dB_c/dt(T) matrix (1,7) or dB_c/dt(Te) matrix (1,8)
! If 'xmhd_two_temp=True', mtmp(1,7)%m => mtmp(1,8)%m
!------------------------------------------------------------------------------
      IF(xmhd_hall.AND.xmhd_adv_temp)THEN
        c1 = k_boltz*jac_dt*det/(elec_charge*u0%N)
        IF(.NOT.xmhd_two_temp)c1=c1*te_factor
        umat(:,1)=-c1*u0%dN
        do jc=1,oft_lagrange%nce
          uvec = umat(:,1)*lag_rop(jc)
          !$omp simd
          do jr=1,oft_hcurl%nce
            mtmp(1,7)%m(jr,jc) = mtmp(1,7)%m(jr,jc) &
              + DOT_PRODUCT(hcurl_copt(jr,:),uvec)
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dB_g/dt matrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Compute dB_g/dt(B_c) matrix (2,1)
!------------------------------------------------------------------------------
      do jc=1,oft_hcurl%nce
        uvec=hcurl_rop(:,jc)*det
        !$omp simd
        do jr=1,oft_hgrad%nce
          mtmp(2,1)%m(jr,jc) = mtmp(2,1)%m(jr,jc) &
            + DOT_PRODUCT(hgrad_ropt(jr,:),uvec)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dB_g/dt(B_g) matrix (2,2)
!------------------------------------------------------------------------------
      do jc=1,oft_hgrad%nce
        uvec=hgrad_rop(:,jc)*det
        !$omp simd
        do jr=1,oft_hgrad%nce
          mtmp(2,2)%m(jr,jc) = mtmp(2,2)%m(jr,jc) &
            + DOT_PRODUCT(hgrad_ropt(jr,:),uvec)
        end do
      end do
    END IF ! Advance B test
!------------------------------------------------------------------------------
! Compute dV/dt matrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Compute dV/dt(B_c) matrix (3:5,1)
!------------------------------------------------------------------------------
    IF(xmhd_jcb)THEN
      c1 = -jac_dt*det/(vel_scale*u0%N*m_ion*mu0)
      do jc=1,oft_hcurl%nce
        uvec = (cross_product(u0%J,hcurl_rop(:,jc)) + cross_product(hcurl_cop(:,jc),u0%B))*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(3,1)%m(jr,jc) = mtmp(3,1)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(4,1)%m(jr,jc) = mtmp(4,1)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(5,1)%m(jr,jc) = mtmp(5,1)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dV/dt(B_g) matrix (3:5,2)
!------------------------------------------------------------------------------
      c1 = -jac_dt*det/(vel_scale*u0%N*m_ion*mu0)
      do jc=1,oft_hgrad%nce
        uvec = cross_product(u0%J,hgrad_rop(:,jc))*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(3,2)%m(jr,jc) = mtmp(3,2)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(4,2)%m(jr,jc) = mtmp(4,2)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(5,2)%m(jr,jc) = mtmp(5,2)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
    END IF
!------------------------------------------------------------------------------
! Compute dV/dt(V) matrix (3:5,3:5)
!------------------------------------------------------------------------------
    !---Upwinding setup
    IF(xmhd_upwind)c1 = jac_dt*det*xmhd_upwind_weight(v0_mag,he,nu_tmp(2)*den_scale/u0%N)
    c2 = den_scale*jac_dt*det/u0%N
    IF(.NOT.xmhd_diss_centered)THEN
      c1 = 2.d0*c1; c2 = 2.d0*c2
    END IF
    IF(xmhd_linear)c2 = jac_dt*det
    c3=2.d0/3.d0
    IF(visc_itype==3)THEN
      !--- Anisotropic tensor (nu_par*bb + nu_perp*(I-bb))
      do j=1,3
        ten1(:,j) = (nu_tmp(1)-nu_tmp(2))*bhat(j)*bhat
        ten1(j,j) = ten1(j,j)+nu_tmp(2)
      end do
    ELSE
      c2 = nu_tmp(1)*c2
    END IF
    !---Advection and mass
    do jc=1,oft_lagrange%nce
      u =  lag_rop(jc)*det
      IF(xmhd_advec)THEN
        u = u + DOT_PRODUCT(u0%V,lag_gop(:,jc))*jac_dt*det
        umat = dvT*lag_rop(jc)*jac_dt*det
      ELSE
        umat = 0.d0
      END IF
      umat(1,1) = umat(1,1) + u
      umat(2,2) = umat(2,2) + u
      umat(3,3) = umat(3,3) + u
      IF(xmhd_upwind)THEN
        uvec = u0%V*(DOT_PRODUCT(u0%V,lag_gop(:,jc)))*c1
      ELSE
        uvec=0.d0
      END IF
      !---Viscous term
      SELECT CASE(visc_itype)
        CASE(1)
          uvec = uvec + c2*lag_gop(:,jc)
          !$omp simd private(u)
          do jr=1,oft_lagrange%nce
            u = DOT_PRODUCT(lag_gopt(jr,:),uvec)
            !---Vx submatrix (1,1:3)
            mtmp(3,3)%m(jr,jc) = mtmp(3,3)%m(jr,jc) + lag_rop(jr)*umat(1,1) + u
            mtmp(3,4)%m(jr,jc) = mtmp(3,4)%m(jr,jc) + lag_rop(jr)*umat(1,2)
            mtmp(3,5)%m(jr,jc) = mtmp(3,5)%m(jr,jc) + lag_rop(jr)*umat(1,3)
            !---Vy submatrix (2,1:3)
            mtmp(4,3)%m(jr,jc) = mtmp(4,3)%m(jr,jc) + lag_rop(jr)*umat(2,1)
            mtmp(4,4)%m(jr,jc) = mtmp(4,4)%m(jr,jc) + lag_rop(jr)*umat(2,2) + u
            mtmp(4,5)%m(jr,jc) = mtmp(4,5)%m(jr,jc) + lag_rop(jr)*umat(2,3)
            !---Vz submatrix (3,1:3)
            mtmp(5,3)%m(jr,jc) = mtmp(5,3)%m(jr,jc) + lag_rop(jr)*umat(3,1)
            mtmp(5,4)%m(jr,jc) = mtmp(5,4)%m(jr,jc) + lag_rop(jr)*umat(3,2)
            mtmp(5,5)%m(jr,jc) = mtmp(5,5)%m(jr,jc) + lag_rop(jr)*umat(3,3) + u
          end do
        CASE(2)
          uvec = uvec + c2*lag_gop(:,jc)
          !$omp simd private(u)
          do jr=1,oft_lagrange%nce
            u = DOT_PRODUCT(lag_gopt(jr,:),uvec)
            !---Vx submatrix (1,1:3)
            mtmp(3,3)%m(jr,jc) = mtmp(3,3)%m(jr,jc) + lag_rop(jr)*umat(1,1) &
              + lag_gopt(jr,1)*uvec(1) - c3*lag_gopt(jr,1)*uvec(1)  + u
            mtmp(3,4)%m(jr,jc) = mtmp(3,4)%m(jr,jc) + lag_rop(jr)*umat(1,2) &
              + lag_gopt(jr,2)*uvec(1) - c3*lag_gopt(jr,1)*uvec(2)
            mtmp(3,5)%m(jr,jc) = mtmp(3,5)%m(jr,jc) + lag_rop(jr)*umat(1,3) &
              + lag_gopt(jr,3)*uvec(1) - c3*lag_gopt(jr,1)*uvec(3)
            !---Vy submatrix (2,1:3)
            mtmp(4,3)%m(jr,jc) = mtmp(4,3)%m(jr,jc) + lag_rop(jr)*umat(2,1) &
              + lag_gopt(jr,1)*uvec(2) - c3*lag_gopt(jr,2)*uvec(1)
            mtmp(4,4)%m(jr,jc) = mtmp(4,4)%m(jr,jc) + lag_rop(jr)*umat(2,2) &
              + lag_gopt(jr,2)*uvec(2) - c3*lag_gopt(jr,2)*uvec(2) + u
            mtmp(4,5)%m(jr,jc) = mtmp(4,5)%m(jr,jc) + lag_rop(jr)*umat(2,3) &
              + lag_gopt(jr,3)*uvec(2) - c3*lag_gopt(jr,2)*uvec(3)
            !---Vz submatrix (3,1:3)
            mtmp(5,3)%m(jr,jc) = mtmp(5,3)%m(jr,jc) + lag_rop(jr)*umat(3,1) &
              + lag_gopt(jr,1)*uvec(3) - c3*lag_gopt(jr,3)*uvec(1)
            mtmp(5,4)%m(jr,jc) = mtmp(5,4)%m(jr,jc) + lag_rop(jr)*umat(3,2) &
              + lag_gopt(jr,2)*uvec(3) - c3*lag_gopt(jr,3)*uvec(2)
            mtmp(5,5)%m(jr,jc) = mtmp(5,5)%m(jr,jc) + lag_rop(jr)*umat(3,3) &
              + lag_gopt(jr,3)*uvec(3) - c3*lag_gopt(jr,3)*uvec(3) + u
          end do
        CASE(3) ! Case 'ani'
          DO jr=1,3
            DO j=1,3
              uten(:,j,jr)=c2*(ten1(:,jr)*lag_gop(j,jc) - c3*ten1(:,j)*lag_gop(jr,jc))
            END DO
          END DO
          uvec = uvec + c2*MATMUL(ten1,lag_gop(:,jc))
          !$omp simd private(u)
          do jr=1,oft_lagrange%nce
            u = DOT_PRODUCT(lag_gopt(jr,:),uvec)
            !---Vx submatrix (1,1:3)
            mtmp(3,3)%m(jr,jc) = mtmp(3,3)%m(jr,jc) + lag_rop(jr)*umat(1,1) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,1,1)) + u
            mtmp(3,4)%m(jr,jc) = mtmp(3,4)%m(jr,jc) + lag_rop(jr)*umat(1,2) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,1,2))
            mtmp(3,5)%m(jr,jc) = mtmp(3,5)%m(jr,jc) + lag_rop(jr)*umat(1,3) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,1,3))
            !---Vy submatrix (2,1:3)
            mtmp(4,3)%m(jr,jc) = mtmp(4,3)%m(jr,jc) + lag_rop(jr)*umat(2,1) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,2,1))
            mtmp(4,4)%m(jr,jc) = mtmp(4,4)%m(jr,jc) + lag_rop(jr)*umat(2,2) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,2,2)) + u
            mtmp(4,5)%m(jr,jc) = mtmp(4,5)%m(jr,jc) + lag_rop(jr)*umat(2,3) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,2,3))
            !---Vz submatrix (3,1:3)
            mtmp(5,3)%m(jr,jc) = mtmp(5,3)%m(jr,jc) + lag_rop(jr)*umat(3,1) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,3,1))
            mtmp(5,4)%m(jr,jc) = mtmp(5,4)%m(jr,jc) + lag_rop(jr)*umat(3,2) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,3,2))
            mtmp(5,5)%m(jr,jc) = mtmp(5,5)%m(jr,jc) + lag_rop(jr)*umat(3,3) &
              + DOT_PRODUCT(lag_gopt(jr,:),uten(:,3,3)) + u
          end do
      END SELECT
    end do
!------------------------------------------------------------------------------
! Compute dV/dt(n) matrix (3:5,6)
!------------------------------------------------------------------------------
    IF(xmhd_adv_den)THEN
      c1 = jac_dt*det*den_scale/(u0%N*u0%N*m_ion*mu0*vel_scale)
      c3 = (u0%Ti+u0%Te)*k_boltz*jac_dt*det*den_scale/(u0%N*m_ion*vel_scale)
      c2 = c3/u0%N
      IF(xmhd_linear)THEN
        c1 = 0.d0
        c2 = c3/(u0%Ti+u0%Te)
      END IF
      do jc=1,oft_lagrange%nce
        uvec = -lag_rop(jc)*u0%dN*c2 + lag_gop(:,jc)*c3
        IF(xmhd_linear)uvec = lag_rop(jc)*(u0%dTi+u0%dTe)*c2 + lag_gop(:,jc)*c3
        IF(xmhd_jcb)uvec = uvec + cross_product(u0%J,u0%B)*lag_rop(jc)*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(3,6)%m(jr,jc) = mtmp(3,6)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(4,6)%m(jr,jc) = mtmp(4,6)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(5,6)%m(jr,jc) = mtmp(5,6)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
    END IF
!------------------------------------------------------------------------------
! Compute dV/dt(T) matrix (3:5,7) and dV/dt(Te) matrix (3:5,8)
!------------------------------------------------------------------------------
    IF(xmhd_adv_temp)THEN
      c2 = k_boltz*jac_dt*det/(m_ion*vel_scale)
      IF(.NOT.xmhd_two_temp)c2=c2*(1.d0+te_factor)
      c1 = c2/u0%N
      do jc=1,oft_lagrange%nce
        uvec = lag_rop(jc)*u0%dN*c1 + lag_gop(:,jc)*c2
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(3,7)%m(jr,jc) = mtmp(3,7)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(4,7)%m(jr,jc) = mtmp(4,7)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(5,7)%m(jr,jc) = mtmp(5,7)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
      IF(xmhd_two_temp)THEN
        do jc=1,oft_lagrange%nce
          uvec = lag_rop(jc)*u0%dN*c1 + lag_gop(:,jc)*c2
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(3,te_ind)%m(jr,jc) = mtmp(3,te_ind)%m(jr,jc) + lag_rop(jr)*uvec(1)
            mtmp(4,te_ind)%m(jr,jc) = mtmp(4,te_ind)%m(jr,jc) + lag_rop(jr)*uvec(2)
            mtmp(5,te_ind)%m(jr,jc) = mtmp(5,te_ind)%m(jr,jc) + lag_rop(jr)*uvec(3)
          end do
        end do
      END IF
    END IF
!------------------------------------------------------------------------------
! Compute dn/dt matrix
!------------------------------------------------------------------------------
    IF(xmhd_adv_den)THEN
      c2 = jac_dt*det
      c1 = c2*vel_scale/den_scale
      !---Upwinding setup
      IF(xmhd_upwind)c3 = jac_dt*det*xmhd_upwind_weight(v0_mag,he,d_dens)
!------------------------------------------------------------------------------
! Compute dn/dt(V) matrix (6,3:5)
!------------------------------------------------------------------------------
      do jc=1,oft_lagrange%nce
        uvec = u0%dN*lag_rop(jc)*c1 + u0%N*lag_gop(:,jc)*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(6,3)%m(jr,jc) = mtmp(6,3)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(6,4)%m(jr,jc) = mtmp(6,4)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(6,5)%m(jr,jc) = mtmp(6,5)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dn/dt(n) matrix (6,6)
!------------------------------------------------------------------------------
      do jc=1,oft_lagrange%nce
        u = DOT_PRODUCT(u0%V,lag_gop(:,jc))*c2 + divv*lag_rop(jc)*c2 + lag_rop(jc)*det
        uvec = d_dens*c2*lag_gop(:,jc)
        IF(xmhd_upwind)uvec = uvec + u0%V*(DOT_PRODUCT(u0%V,lag_gop(:,jc)))*c3
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(6,6)%m(jr,jc) = mtmp(6,6)%m(jr,jc) + lag_rop(jr)*u &
            + DOT_PRODUCT(lag_gopt(jr,:),uvec)
        end do
      end do
      !---Hyper diffusivity
      IF(n2_ind>0)THEN
        c1 = -d2_dens*jac_dt*det*n2_scale/den_scale
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0*c1
        do jc=1,oft_lagrange%nce
          uvec = lag_gop(:,jc)*c1
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(6,n2_ind)%m(jr,jc) = mtmp(6,n2_ind)%m(jr,jc) &
              + DOT_PRODUCT(lag_gopt(jr,:),uvec)
          end do
          uvec = lag_gop(:,jc)*det*den_scale/n2_scale
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(n2_ind,6)%m(jr,jc) = mtmp(n2_ind,6)%m(jr,jc) &
              + DOT_PRODUCT(lag_gopt(jr,:),uvec)
          end do
          u = lag_rop(jc)*det
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(n2_ind,n2_ind)%m(jr,jc) = mtmp(n2_ind,n2_ind)%m(jr,jc) &
              + lag_rop(jr)*u
          end do
        end do
      END IF
    END IF
!------------------------------------------------------------------------------
! Compute dT/dt matrix
!------------------------------------------------------------------------------
    IF(xmhd_adv_temp)THEN
!------------------------------------------------------------------------------
! Compute dT/dt(B_c) matrix (7,1)
!------------------------------------------------------------------------------
      IF(xmhd_ohmic.AND.xmhd_adv_b.AND.(.NOT.xmhd_two_temp))THEN
        c1 = eta_curr*tgam_fac*jac_dt*det*2.d0/((1.d0+te_factor)*mu0*u0%N*k_boltz)
        do jc=1,oft_hcurl%nce
          u = -DOT_PRODUCT(u0%J,hcurl_cop(:,jc))*c1
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(7,1)%m(jr,jc) = mtmp(7,1)%m(jr,jc) + lag_rop(jr)*u
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dT/dt(V) matrix (7,3:5)
!------------------------------------------------------------------------------
      c3 = jac_dt*det*vel_scale ! advection
      c1 = u0%Ti*tgam_fac*c3 ! compression
      c2 = m_ion*tgam_fac*den_scale*c3/(k_boltz*u0%N)
      !---Viscous term
      IF(xmhd_visc_heat)THEN
        IF(.NOT.xmhd_two_temp)c2 = c2/(1.d0+te_factor)
        SELECT CASE(visc_itype)
          CASE(1)
            umat = -c2*nu_tmp(1)*u0%dV
          CASE(2)
            !--- Nabla(V) + (Nabla(V))^T
            umat = -c2*nu_tmp(1)*(u0%dV+dvT)
            !--- (- 2/3 I Div(V))
            umat(1,1) = umat(1,1) + c2*nu_tmp(1)*2.d0*divv/3.d0
            umat(2,2) = umat(2,2) + c2*nu_tmp(1)*2.d0*divv/3.d0
            umat(3,3) = umat(3,3) + c2*nu_tmp(1)*2.d0*divv/3.d0
          CASE(3)
            !--- Anisotropic tensor (nu_par*bb + nu_perp*(I-bb))
            do j=1,3
              ten1(:,j) = (nu_tmp(1)-nu_tmp(2))*bhat(j)*bhat
              ten1(j,j) = ten1(j,j)+nu_tmp(2)
            end do
            ten1=c2*ten1/2.d0
            umat = -MATMUL(ten1,dvT+2.d0*u0%dV) - MATMUL(dvT,ten1) &
              + 2.d0/3.d0*ten1*divv
            u = SUM(SUM(ten1*u0%dV,DIM=2))
            umat(1,1) = umat(1,1) + 2.d0/3.d0*u
            umat(2,2) = umat(2,2) + 2.d0/3.d0*u
            umat(3,3) = umat(3,3) + 2.d0/3.d0*u
        END SELECT
      ELSE
        umat=0.d0
      END IF
      do jc=1,oft_lagrange%nce
        uvec = lag_gop(:,jc)*c1 + u0%dTi*lag_rop(jc)*c3
        uvec = uvec + MATMUL(lag_gop(:,jc),umat)
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(7,3)%m(jr,jc) = mtmp(7,3)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(7,4)%m(jr,jc) = mtmp(7,4)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(7,5)%m(jr,jc) = mtmp(7,5)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dT/dt(n) matrix (7,6)
!------------------------------------------------------------------------------
      IF(xmhd_ohmic.AND.xmhd_adv_den.AND.(.NOT.xmhd_two_temp))THEN
        c1 = eta_curr*tgam_fac*DOT_PRODUCT(u0%J,u0%J)*den_scale*jac_dt*det/ &
          ((1.d0+te_factor)*mu0*u0%N*u0%N*k_boltz)
        do jc=1,oft_lagrange%nce
          u = lag_rop(jc)*c1
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(7,6)%m(jr,jc) = mtmp(7,6)%m(jr,jc) + lag_rop(jr)*u
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dT/dt(T) matrix (7,7)
!------------------------------------------------------------------------------
      !---Upwinding setup
      IF(xmhd_upwind)c1 = jac_dt*det*xmhd_upwind_weight(v0_mag,he,chi_temp(2)*tgam_fac)
      c2 = divv*tgam_fac ! compression
      IF(xmhd_therm_equil)c2 = c2 + 1.d0/elec_ion_therm_rate(u0%Te,u0%N,mu_ion)
      c3 = tgam_fac*jac_dt*det
      IF(.NOT.xmhd_diss_centered)THEN
        c1 = 2.d0*c1; c3 = 2.d0*c3
      END IF
      do jc=1,oft_lagrange%nce
        u = (lag_rop(jc)*c2 + DOT_PRODUCT(u0%V,lag_gop(:,jc)))*jac_dt*det
        u = u + lag_rop(jc)*det
        !---Thermal conduction
        uvec = ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(lag_gop(:,jc),bhat)*bhat &
          + chi_temp(2)*lag_gop(:,jc))*c3
        u = u - ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(u0%dN,bhat)*DOT_PRODUCT(bhat,lag_gop(:,jc)) &
          + chi_temp(2)*DOT_PRODUCT(u0%dN,lag_gop(:,jc)))*c3/u0%N
        !---Upwinding
        IF(xmhd_upwind)uvec = uvec + u0%V*(DOT_PRODUCT(u0%V,lag_gop(:,jc)))*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(7,7)%m(jr,jc) = mtmp(7,7)%m(jr,jc) + lag_rop(jr)*u &
            + DOT_PRODUCT(lag_gopt(jr,:),uvec)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dTe/dt(B_c) matrix (8,1)
!------------------------------------------------------------------------------
    IF(xmhd_two_temp)THEN
      IF(xmhd_adv_b)THEN
        c1 = 0.d0
        IF(xmhd_ohmic)c1 = -2.d0*eta_curr*tgam_fac*jac_dt*det/(mu0*u0%N*k_boltz)
        c2 = -jac_dt*det/(u0%N*elec_charge*mu0) ! advection
        do jc=1,oft_hcurl%nce
          u = DOT_PRODUCT(u0%J*c1 + u0%dTe*c2,hcurl_cop(:,jc))
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(te_ind,1)%m(jr,jc) = mtmp(te_ind,1)%m(jr,jc) + lag_rop(jr)*u
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dTe/dt(V) matrix (8,3:5)
!------------------------------------------------------------------------------
      c1 = jac_dt*det*vel_scale ! advection
      c2 = u0%Te*tgam_fac*c1 ! compression
      do jc=1,oft_lagrange%nce
        uvec = lag_gop(:,jc)*c2 + u0%dTe*lag_rop(jc)*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(te_ind,3)%m(jr,jc) = mtmp(te_ind,3)%m(jr,jc) + lag_rop(jr)*uvec(1)
          mtmp(te_ind,4)%m(jr,jc) = mtmp(te_ind,4)%m(jr,jc) + lag_rop(jr)*uvec(2)
          mtmp(te_ind,5)%m(jr,jc) = mtmp(te_ind,5)%m(jr,jc) + lag_rop(jr)*uvec(3)
        end do
      end do
!------------------------------------------------------------------------------
! Compute dTe/dt(n) matrix (8,6)
!------------------------------------------------------------------------------
      IF(xmhd_adv_den)THEN
        c1 = 0.d0
        IF(xmhd_ohmic)c1 = eta_curr*tgam_fac*DOT_PRODUCT(u0%J,u0%J)/(mu0*u0%N*u0%N*k_boltz)
        c1 = c1 + DOT_PRODUCT(u0%J,u0%dTe)/(mu0*u0%N*u0%N*elec_charge) ! advection
        c1 = c1*den_scale*jac_dt*det
        do jc=1,oft_lagrange%nce
          u = lag_rop(jc)*c1
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(te_ind,6)%m(jr,jc) = mtmp(te_ind,6)%m(jr,jc) + lag_rop(jr)*u
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dT/dt(Te) matrix (7,8) and dTe/dt(T) matrix (8,7)
!------------------------------------------------------------------------------
      IF(xmhd_therm_equil)THEN
        c1 = (1.d0 - 3.d0*u0%Ti/u0%Te)*jac_dt*det/(2.d0*elec_ion_therm_rate(u0%Te,u0%N,mu_ion))
        c2 = -jac_dt*det/elec_ion_therm_rate(u0%Te,u0%N,mu_ion)
        do jc=1,oft_lagrange%nce
          !$omp simd
          do jr=1,oft_lagrange%nce
            mtmp(7,te_ind)%m(jr,jc) = mtmp(7,te_ind)%m(jr,jc) + lag_rop(jr)*lag_rop(jc)*c1
            mtmp(te_ind,7)%m(jr,jc) = mtmp(te_ind,7)%m(jr,jc) + lag_rop(jr)*lag_rop(jc)*c2
          end do
        end do
      END IF
!------------------------------------------------------------------------------
! Compute dTe/dt(Te) matrix (8,8)
!------------------------------------------------------------------------------
      IF(xmhd_brag)THEN
        chi_temp = brag_elec_transport(u0%N,u0%Te,bmag)/u0%N
        chi_temp(1) = kappa_par*chi_temp(1)
        chi_temp(2) = MAX(1.d1, kappa_perp*chi_temp(2))
      END IF
      umat(:,1) = u0%V-u0%J/(mu0*u0%N*elec_charge) ! advection
      !---Upwinding setup
      IF(xmhd_upwind)THEN
        v0_mag = magnitude(umat(:,1))
        c1 = jac_dt*det*xmhd_upwind_weight(v0_mag,he,chi_temp(2)*tgam_fac)
      END IF
      c2 = divv*tgam_fac ! compression
      IF(xmhd_therm_equil)c2 = c2 - (1.d0 - 3.d0*u0%Ti/u0%Te)/(2.d0*elec_ion_therm_rate(u0%Te,u0%N,mu_ion))
      c3 = tgam_fac*jac_dt*det
      IF(.NOT.xmhd_diss_centered)c3 = 2.d0*c3
      do jc=1,oft_lagrange%nce
        u = (lag_rop(jc)*c2 + DOT_PRODUCT(umat(:,1),lag_gop(:,jc)))*jac_dt*det
        u = u + lag_rop(jc)*det
        !---Thermal conduction
        uvec = ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(lag_gop(:,jc),bhat)*bhat &
          + chi_temp(2)*lag_gop(:,jc))*c3
        u = u - ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(u0%dN,bhat)*DOT_PRODUCT(bhat,lag_gop(:,jc)) &
          + chi_temp(2)*DOT_PRODUCT(u0%dN,lag_gop(:,jc)))*c3/u0%N
        !---Upwinding
        IF(xmhd_upwind)uvec = uvec + 2.d0*umat(:,1)*(DOT_PRODUCT(umat(:,1),lag_gop(:,jc)))*c1
        !$omp simd
        do jr=1,oft_lagrange%nce
          mtmp(te_ind,te_ind)%m(jr,jc) = mtmp(te_ind,te_ind)%m(jr,jc) + lag_rop(jr)*u &
            + DOT_PRODUCT(lag_gopt(jr,:),uvec)
        end do
      end do
    END IF
    END IF
!------------------------------------------------------------------------------
! End quadrature loop
!------------------------------------------------------------------------------
  end do
  !------------------------------------------------------------------------------
  ! Apply BCs to local matrices
  !------------------------------------------------------------------------------
  CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%b1_bc(j_hcurl),1)
  CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%b2_bc(j_hgrad),2)
  CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%n_bc(j_lag),6)
  CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%t_bc(j_lag),7)
  IF(xmhd_two_temp)CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%t_bc(j_lag),te_ind)
  !---Velocity BC
  IF(xmhd_vbcdir)THEN
    DO jr=1,oft_lagrange%nce
      IF(oft_lagrange%global%gbe(j_lag(jr)))THEN
        CALL lag_vbc_tensor(oft_lagrange,j_lag(jr),vbc_type,umat)
        DO m=1,xmhd_rep%nfields
          IF(ASSOCIATED(mtmp(3,m)%m))THEN
            DO jc=1,SIZE(mtmp(3,m)%m,2)
              !---Unpack
              uvec(1)=mtmp(3,m)%m(jr,jc)
              uvec(2)=mtmp(4,m)%m(jr,jc)
              uvec(3)=mtmp(5,m)%m(jr,jc)
              !---Apply BC
              uvec = MATMUL(umat,uvec)
              !---Repack
              mtmp(3,m)%m(jr,jc)=uvec(1)
              mtmp(4,m)%m(jr,jc)=uvec(2)
              mtmp(5,m)%m(jr,jc)=uvec(3)
            END DO
          END IF
        END DO
      END IF
    END DO
  ELSE
    CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%v_bc(j_lag),3)
    CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%v_bc(j_lag),4)
    CALL xmhd_rep%mat_zero_local_rows(mtmp,oft_xmhd_ops%v_bc(j_lag),5)
  END IF
!------------------------------------------------------------------------------
! Add local contribution to full matrix
!------------------------------------------------------------------------------
  IF(xmhd_hall.AND.xmhd_two_temp)NULLIFY(mtmp(1,7)%m)
  CALL xmhd_rep%mat_add_local(Jac,mtmp,iloc,tlocks)
end do
END DO
!---Delete temporary matrices
deallocate(j_lag,j_hgrad,j_hcurl)
deallocate(lag_rop,lag_gop,hgrad_rop)
deallocate(hcurl_rop,hcurl_cop)
deallocate(lag_gopt,hgrad_ropt,hcurl_ropt,hcurl_copt)
NULLIFY(xmhd_lag_rop,xmhd_lag_gop,xmhd_hgrad_rop,xmhd_hcurl_rop,xmhd_hcurl_cop)
!---Destroy local matrix
CALL xmhd_rep%mat_destroy_local(mtmp)
deallocate(mtmp,iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,xmhd_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
!------------------------------------------------------------------------------
! Set BC diagnoals
!------------------------------------------------------------------------------
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_hcurl,Jac,oft_xmhd_ops%b1_bc,1)
CALL fem_dirichlet_diag(oft_hgrad,Jac,oft_xmhd_ops%b2_bc,2)
CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%n_bc,6)
CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%t_bc,7)
IF(xmhd_two_temp)CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%t_bc,te_ind)
!---Velocity BC
IF(xmhd_vbcdir)THEN
  DO i=1,oft_lagrange%nbe
    IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
    j=oft_lagrange%lbe(i)
    IF(.NOT.oft_lagrange%global%gbe(j))CYCLE
    jtmp=j
    CALL lag_vbc_diag(oft_lagrange,j,vbc_type,umat)
    DO jr=1,3
      DO jc=1,3
        CALL Jac%add_values(jtmp,jtmp,umat(jr,jc),1,1,jr+2,jc+2)
      END DO
    END DO
  END DO
ELSE
  CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%v_bc,3)
  CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%v_bc,4)
  CALL fem_dirichlet_diag(oft_lagrange,Jac,oft_xmhd_ops%v_bc,5)
END IF
!------------------------------------------------------------------------------
! Final assembly
!------------------------------------------------------------------------------
call oft_xmhd_create(tmp)
call Jac%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp)
DEBUG_STACK_POP
end subroutine xmhd_build_ops
!------------------------------------------------------------------------------
!> Setup and allocate operators used in xMHD advance
!------------------------------------------------------------------------------
subroutine xmhd_alloc_ops
integer(i4) :: i,j,level,levelin,offset
INTEGER(i4), ALLOCATABLE :: color_tmp(:)
integer(i4), ALLOCATABLE :: j_lag(:)
REAL(r8), POINTER :: node_flag(:)
CLASS(oft_vector), POINTER :: vectmp
type(xmhd_ops), pointer :: ops
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Setup Operators and ML environment
!------------------------------------------------------------------------------
allocate(ml_J(xmhd_ML_hcurl%nlevels-xmhd_minlev+1),ml_int(xmhd_ML_hcurl%nlevels-xmhd_minlev))
allocate(oft_xmhd_ops_ML(xmhd_ML_hcurl%nlevels))
xmhd_blevel=xmhd_ML_hcurl%blevel
xmhd_nlevels=xmhd_ML_hcurl%nlevels
xmhd_level=xmhd_ML_hcurl%level
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Allocating xMHD structures'
!------------------------------------------------------------------------------
! Setup matrix mask based on included physics
!------------------------------------------------------------------------------
ALLOCATE(xmhd_mat_mask(xmhd_rep%nfields,xmhd_rep%nfields))
xmhd_mat_mask=0
!---B_c rows
IF(xmhd_adv_b)THEN
  xmhd_mat_mask(1,1:5)=1
  IF(xmhd_hall)THEN
    IF(xmhd_adv_den)xmhd_mat_mask(1,6)=1
    IF(xmhd_adv_temp)THEN
      IF(xmhd_two_temp)THEN
        xmhd_mat_mask(1,te_ind)=1
      ELSE
        xmhd_mat_mask(1,7)=1
      END IF
    END IF
  END IF
  !---B_g rows
  xmhd_mat_mask(2,1:2)=1
  !---Hyper-resistivity
  IF(j2_ind>0)THEN
    xmhd_mat_mask(1,j2_ind)=1
    xmhd_mat_mask(j2_ind,1)=1
    xmhd_mat_mask(j2_ind,j2_ind)=1
  END IF
ELSE
  xmhd_mat_mask(1,1)=2
  xmhd_mat_mask(2,2)=2
END IF
!---V rows
DO j=0,2
  IF(xmhd_jcb)xmhd_mat_mask(3+j,1:2)=1
  xmhd_mat_mask(3+j,3:5)=1
  IF(xmhd_adv_den)xmhd_mat_mask(3+j,6)=1
  IF(xmhd_adv_temp)THEN
    xmhd_mat_mask(3+j,7)=1
    IF(xmhd_two_temp)xmhd_mat_mask(3+j,te_ind)=1
  END IF
END DO
!---n rows
IF(xmhd_adv_den)THEN
  xmhd_mat_mask(6,3:6)=1
  IF(n2_ind>0)THEN
    xmhd_mat_mask(6,n2_ind)=1
    xmhd_mat_mask(n2_ind,6)=1
    xmhd_mat_mask(n2_ind,n2_ind)=1
  END IF
ELSE
  xmhd_mat_mask(6,6)=2
END IF
!---T rows
IF(xmhd_adv_temp)THEN
  xmhd_mat_mask(7,3:5)=1
  xmhd_mat_mask(7,7)=1
  !---Te rows
  IF(xmhd_two_temp)THEN
    IF(xmhd_adv_b)xmhd_mat_mask(te_ind,1)=1
    xmhd_mat_mask(te_ind,3:5)=1
    IF(xmhd_adv_den)xmhd_mat_mask(te_ind,6)=1
    xmhd_mat_mask(te_ind,te_ind)=1
    IF(xmhd_therm_equil)THEN
      xmhd_mat_mask(7,te_ind)=1
      xmhd_mat_mask(te_ind,7)=1
    END IF
  ELSE
    IF(xmhd_ohmic.AND.xmhd_adv_b)xmhd_mat_mask(7,1)=1
    IF(xmhd_ohmic.AND.xmhd_adv_den)xmhd_mat_mask(7,6)=1
  END IF
ELSE
  xmhd_mat_mask(7,7)=2
  IF(xmhd_two_temp)xmhd_mat_mask(te_ind,te_ind)=2
END IF
!------------------------------------------------------------------------------
! Create matrix objects for each level
!------------------------------------------------------------------------------
levelin=xmhd_level
do level=levelin,xmhd_minlev,-1
  call xmhd_set_level(level)
  ops=>oft_xmhd_ops
  CALL xmhd_rep%mat_create(ops%J,xmhd_mat_mask)
  ml_J(level-xmhd_minlev+1)%m=>ops%J
  !---Set element coloring
  IF(xmhd_nparts(level)>1)THEN
    SELECT TYPE(this=>ops%J)
      CLASS IS(oft_native_matrix)
        ALLOCATE(this%color(this%nr))
        this%color=0
        offset=0
        DO i=1,xmhd_rep%nfields
          ALLOCATE(color_tmp(xmhd_rep%fields(i)%fe%ne))
          SELECT TYPE(this=>xmhd_rep%fields(i)%fe)
          CLASS IS(oft_fem_type)
            CALL fem_partition(this,color_tmp,xmhd_nparts(level))
          CLASS DEFAULT
            CALL oft_abort("Invalid FE representation for partitioning", &
              "xmhd_alloc_ops",__FILE__)
          END SELECT
          DO j=1,xmhd_rep%fields(i)%fe%ne
            this%color(offset+j)=color_tmp(j)
          END DO
          offset=offset+xmhd_rep%fields(i)%fe%ne
          DEALLOCATE(color_tmp)
        END DO
    END SELECT
  END IF
end do
CALL ML_xmhd_rep%build_interp(xmhd_minlev)
!---Alias
do level=levelin,xmhd_minlev,-1
  call xmhd_set_level(level)
  ops=>oft_xmhd_ops
  IF(level>xmhd_minlev)THEN
    ops%interp=>ML_xmhd_rep%interp_matrices(level)%m
    ml_int(level-xmhd_minlev)%M=>ML_xmhd_rep%interp_matrices(level)%m
  END IF
end do
!---Setup solid node flags
NULLIFY(node_flag,vectmp)
DO level=levelin,xmhd_minlev,-1
  call xmhd_set_level(level)
  !
  ALLOCATE(oft_xmhd_ops%solid_node(oft_lagrange%ne))
  oft_xmhd_ops%solid_node=.FALSE.
  CALL xmhd_ML_vlagrange%vec_create(vectmp)
  CALL vectmp%set(0.d0)
  CALL vectmp%get_local(node_flag)
  !
  ALLOCATE(j_lag(oft_lagrange%nce))
  DO i=1,mesh%nc
    IF(solid_cell(i))THEN
      CALL oft_lagrange%ncdofs(i,j_lag)
      node_flag(j_lag)=1.d0
    END IF
  END DO
  DEALLOCATE(j_lag)
  CALL vectmp%restore_local(node_flag,add=.TRUE.)
  CALL vectmp%get_local(node_flag)
  DO i=1,oft_lagrange%ne
    IF(node_flag(i)>0.d0)oft_xmhd_ops%solid_node(i)=.TRUE.
  END DO
  CALL vectmp%delete
  DEALLOCATE(node_flag,vectmp)
  CALL xmhd_set_bc()
END DO
call xmhd_set_level(levelin)
! IF(oft_debug_print(2))WRITE(*,*)'  Done'
DEBUG_STACK_POP
end subroutine xmhd_alloc_ops
!------------------------------------------------------------------------------
!> Setup material regions from XML input file
!------------------------------------------------------------------------------
subroutine xmhd_setup_regions()
!---XML solver fields
#ifdef HAVE_XML
integer(i4) :: nread_id,nread_eta,nread_type,ierr,i,j,reg_type(1)
real(r8) :: eta(1)
TYPE(xml_node), POINTER :: reg_node,inner_node
TYPE(xml_nodelist) :: reg_nodes
#endif
integer(i4), ALLOCATABLE :: regs(:),reg_types(:)
DEBUG_STACK_PUSH
ALLOCATE(regs(mesh%nreg),reg_types(mesh%nreg))
reg_types=1
IF(.NOT.ALLOCATED(eta_reg))THEN
  ALLOCATE(eta_reg(mesh%nreg),solid_cell(mesh%nc))
END IF
eta_reg=-1.d0
solid_cell=.FALSE.
#ifdef HAVE_XML
IF(ASSOCIATED(xmhd_root_node))THEN
  !---Look for pre node
  CALL xml_get_element(xmhd_root_node,"region",reg_nodes,ierr)
  IF(reg_nodes%n>0)THEN
    DO i=0,reg_nodes%n-1
      reg_node=>reg_nodes%nodes(i+1)%this
      !---
      CALL xml_get_element(reg_node,"id",inner_node,ierr)
      IF(ierr/=0)CALL oft_abort("Error reading regions IDs for group", &
        "xmhd_setup_regions",__FILE__)
      CALL xml_extractDataContent(inner_node,regs,num=nread_id,iostat=ierr)
      IF(nread_id==0)CALL oft_abort("Zero values given in id group", &
        "xmhd_setup_regions",__FILE__)
      IF(ierr>0)CALL oft_abort("Too many id values specified","xmhd_setup_regions", &
      __FILE__)
      IF(ANY(regs(1:nread_id)>mesh%nreg).OR.ANY(regs(1:nread_id)<=0))CALL oft_abort( &
        "Invalid region ID","xmhd_setup_regions",__FILE__)
      !---
      CALL xml_get_element(reg_node,"eta",inner_node,ierr)
      CALL xml_extractDataContent(inner_node,eta,num=nread_eta,iostat=ierr)
      IF(nread_eta==0)CALL oft_abort("Zero values given in eta group", &
        "xmhd_setup_regions",__FILE__)
      IF(ierr>0)CALL oft_abort("Too many eta values specified","xmhd_setup_regions", &
        __FILE__)
      IF(eta(1)<0.d0)CALL oft_abort("Invalid eta value specified","xmhd_setup_regions", &
        __FILE__)
      !---Get region type
      CALL xml_get_element(reg_node,"type",inner_node,ierr)
      IF(ierr/=0)THEN
        reg_type(1)=2.d0
      ELSE
        CALL xml_extractDataContent(inner_node,reg_type,num=nread_type,iostat=ierr)
        IF(nread_eta==0)CALL oft_abort("Zero values given in type group", &
          "xmhd_setup_regions",__FILE__)
        IF(ierr>0)CALL oft_abort("Too many type values specified","xmhd_setup_regions", &
          __FILE__)
        IF(reg_type(1)<1.OR.reg_type(1)>2)CALL oft_abort("Invalid type specified","xmhd_setup_regions", &
          __FILE__)
      END IF
      !---
      DO j=1,nread_id
        IF(eta_reg(regs(j))>0.d0)THEN
          CALL oft_abort("Region blocks overlap","xmhd_setup_regions",__FILE__)
        ELSE
          eta_reg(regs(j))=eta(1)
        END IF
        reg_types(regs(j))=reg_type(1)
      END DO
    END DO
  END IF
END IF
!
xmhd_rw=.FALSE.
DO i=1,mesh%nc
  IF(reg_types(mesh%reg(i))==2)THEN
    xmhd_rw=.TRUE.
    solid_cell(i)=.TRUE.
  END IF
END DO
#ifdef HAVE_MPI
call MPI_ALLREDUCE(MPI_IN_PLACE,xmhd_rw,1,OFT_MPI_LOGICAL,MPI_LOR,oft_env%COMM,ierr)
#endif
#endif
eta_reg=ABS(eta_reg)
DEBUG_STACK_POP
end subroutine xmhd_setup_regions
!------------------------------------------------------------------------------
!> Set BC flags
!------------------------------------------------------------------------------
subroutine xmhd_set_bc
integer(i4) :: i,j
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Set boundary condition for B
!------------------------------------------------------------------------------
ALLOCATE(oft_xmhd_ops%b1_bc(oft_hcurl%ne))
ALLOCATE(oft_xmhd_ops%b2_bc(oft_hgrad%ne))
IF(xmhd_adv_b)THEN
  oft_xmhd_ops%b1_bc=.FALSE.
  oft_xmhd_ops%b2_bc=.FALSE.
  IF(bbc(1:2)=='bc')THEN
    !$omp parallel do
    DO i=1,mesh%np
      oft_xmhd_ops%b2_bc(i)=.TRUE.
    END DO
  ELSE IF(bbc(1:2)=='ic')THEN
    !$omp parallel do
    DO i=1,mesh%np
      IF(mesh%global%gbp(i))CYCLE
      oft_xmhd_ops%b2_bc(i)=.TRUE.
    END DO
    !$omp parallel do private(j)
    DO i=1,oft_hcurl%nbe
      j=oft_hcurl%lbe(i)
      IF(oft_hcurl%global%gbe(j))oft_xmhd_ops%b1_bc(j)=.TRUE.
    END DO
  ELSE
    CALL oft_abort('Invalid B Boundary Condition.','xmhd_bc',__FILE__)
  END IF
ELSE
  oft_xmhd_ops%b1_bc=.TRUE.
  oft_xmhd_ops%b2_bc=.TRUE.
END IF
!------------------------------------------------------------------------------
! Apply boundary condition to V
!------------------------------------------------------------------------------
ALLOCATE(oft_xmhd_ops%v_bc(oft_lagrange%ne))
oft_xmhd_ops%v_bc=.FALSE.
IF(vbc(1:4)=='none'.OR.vbc(1:4)=='norm')THEN
  !---Do nothing
ELSE IF(vbc(1:3)=='all')THEN
  oft_xmhd_ops%v_bc=(oft_lagrange%global%gbe.OR.oft_xmhd_ops%solid_node)
ELSE
  CALL oft_abort('Invalid V Boundary Condition.','xmhd_bc',__FILE__)
END IF
!------------------------------------------------------------------------------
! Apply boundary condition to n
!------------------------------------------------------------------------------
ALLOCATE(oft_xmhd_ops%n_bc(oft_lagrange%ne))
IF(xmhd_adv_den)THEN
  oft_xmhd_ops%n_bc=.FALSE.
  IF(nbc=='n')THEN
    !---Do nothing
  ELSE IF(nbc=='d')THEN
    oft_xmhd_ops%n_bc=(oft_lagrange%global%gbe.OR.oft_xmhd_ops%solid_node)
  ELSE
    CALL oft_abort('Invalid N Boundary Condition.','xmhd_bc',__FILE__)
  END IF
ELSE
  oft_xmhd_ops%n_bc=.TRUE.
END IF
!------------------------------------------------------------------------------
! Apply boundary condition to T
!------------------------------------------------------------------------------
ALLOCATE(oft_xmhd_ops%t_bc(oft_lagrange%ne))
IF(xmhd_adv_temp)THEN
  oft_xmhd_ops%t_bc=.FALSE.
  IF(tbc=='n')THEN
    !---Do nothing
  ELSE IF(tbc=='d')THEN
    oft_xmhd_ops%t_bc=(oft_lagrange%global%gbe.OR.oft_xmhd_ops%solid_node)
  ELSE
    CALL oft_abort('Invalid T Boundary Condition.','xmhd_bc',__FILE__)
  END IF
ELSE
  oft_xmhd_ops%t_bc=.TRUE.
END IF
DEBUG_STACK_POP
end subroutine xmhd_set_bc
!------------------------------------------------------------------------------
!> Compute the NL metric for solution field
!!
!! b = F(a)
!------------------------------------------------------------------------------
subroutine xmhd_errmatrix_apply(self,a,b)
class(oft_xmhd_errmatrix), intent(inout) :: self !< Error matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
integer(i4) :: i,ii,ip,j,m,jr,jp,ncemax,nrows,ptind(2)
integer(i4), allocatable, dimension(:) :: j_lag,j_hgrad,j_hcurl
real(r8) :: f,diff,chi_temp(2),bmag,det
real(r8) :: mloc(3,3),goptmp(3,4),vol,eta_curr,neg_vols(3),divv
real(r8) :: bhat(3),jp0q(3),np0q,vad(3),vad_mag,he,c1,cgop(3,6),nu_tmp(3),ten1(3,3)
real(r8) :: diag_vals(7),pt(3),pt_r,floor_tmp(3),fulltmp(40),u,uvec(3),umat(3,3)
real(r8), allocatable, dimension(:) :: bc_loc,bg_loc,jpin_loc,te_loc,j2_loc,n2_loc
real(r8), allocatable, dimension(:,:) :: vpin_loc,lag_loc
real(r8), target, allocatable, dimension(:) :: lag_rop
real(r8), target, allocatable, dimension(:,:) :: lag_gop,hgrad_rop,hcurl_rop,hcurl_cop
real(r8), pointer, contiguous, dimension(:) :: jpin,npin
real(r8), pointer, contiguous, dimension(:) :: bcout,bgout,nout,tout,teout,j2out,n2out
real(r8), pointer, contiguous, dimension(:,:) :: vpin,vout
real(r8), pointer, dimension(:) :: vtmp
logical :: curved,back_step
type(oft_quad_type), pointer :: quad
TYPE(xmhd_interp) :: full_interp
type(xmhd_loc_values) :: u0
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2X,A)')'Apply xMHD non-linear function'
quad=>oft_hcurl%quad
neg_vols=0.d0
!------------------------------------------------------------------------------
! Get local field values from source
!------------------------------------------------------------------------------
full_interp%u=>a
CALL full_interp%setup(mesh)
!------------------------------------------------------------------------------
! Get local field values from previous step if necessary
!------------------------------------------------------------------------------
NULLIFY(jpin,npin,vpin)
back_step=.FALSE.
IF(ASSOCIATED(self%up))THEN
  CALL self%up%get_local(jpin,1)
  CALL self%up%get_local(npin,6)
  IF(xmhd_upwind)THEN
    ALLOCATE(vpin(oft_lagrange%ne,3))
    vtmp=>vpin(:,1)
    CALL self%up%get_local(vtmp,3)
    vtmp=>vpin(:,2)
    CALL self%up%get_local(vtmp,4)
    vtmp=>vpin(:,3)
    CALL self%up%get_local(vtmp,5)
  END IF
  back_step=.TRUE.
END IF
!------------------------------------------------------------------------------
! Get local field values for result
!------------------------------------------------------------------------------
NULLIFY(bcout,bgout,nout,tout,teout,j2out,n2out)
CALL b%set(0.d0)
CALL b%get_local(bcout,1)
CALL b%get_local(bgout,2)
ALLOCATE(vout(oft_lagrange%ne,3))
vtmp=>vout(:,1)
CALL b%get_local(vtmp,3)
vtmp=>vout(:,2)
CALL b%get_local(vtmp,4)
vtmp=>vout(:,3)
CALL b%get_local(vtmp,5)
CALL b%get_local(nout,6)
CALL b%get_local(tout,7)
IF(xmhd_two_temp)CALL b%get_local(teout,te_ind)
IF(j2_ind>0)CALL b%get_local(j2out,j2_ind)
IF(n2_ind>0)CALL b%get_local(n2out,n2_ind)
!---Zero diagnostic values
diag_vals=0.d0
ptind=(/1,2/)
IF(xmhd_taxis==1)ptind(1)=3
IF(xmhd_taxis==2)ptind(2)=3
!------------------------------------------------------------------------------
! Apply metric function
!------------------------------------------------------------------------------
!$omp parallel private(det,j_lag,j_hgrad,j_hcurl,lag_rop,lag_gop,hgrad_rop,hcurl_rop, &
!$omp hcurl_cop,curved,goptmp,vol,j,m,jr,bhat,jp0q,vad,vad_mag,he,c1,cgop,divv, &
!$omp np0q,eta_curr,u,bc_loc,bg_loc,lag_loc,pt_r,vpin_loc,jpin_loc,chi_temp,te_loc, &
!$omp j2_loc,n2_loc,bmag,uvec,umat,pt,nu_tmp,ten1,i,ii,floor_tmp,u0,fulltmp) &
!$omp reduction(+:neg_vols) reduction(+:diag_vals)
allocate(j_lag(oft_lagrange%nce),j_hgrad(oft_hgrad%nce),j_hcurl(oft_hcurl%nce))
allocate(lag_rop(oft_lagrange%nce),lag_gop(3,oft_lagrange%nce))
allocate(hgrad_rop(3,oft_hgrad%nce))
allocate(hcurl_rop(3,oft_hcurl%nce),hcurl_cop(3,oft_hcurl%nce))
allocate(bc_loc(oft_hcurl%nce),bg_loc(oft_hgrad%nce))
allocate(lag_loc(oft_lagrange%nce,5))
IF(xmhd_two_temp)allocate(te_loc(oft_lagrange%nce))
IF(j2_ind>0)allocate(j2_loc(oft_hcurl%nce))
IF(n2_ind>0)allocate(n2_loc(oft_lagrange%nce))
IF(back_step)allocate(vpin_loc(oft_lagrange%nce,4),jpin_loc(oft_hcurl%nce))
!--- Define viscosity temporary varible (currently static)
nu_tmp(1) = nu_par
nu_tmp(2) = nu_perp
!---Setup caching
xmhd_lag_rop=>lag_rop
xmhd_lag_gop=>lag_gop
xmhd_hgrad_rop=>hgrad_rop
xmhd_hcurl_rop=>hcurl_rop
xmhd_hcurl_cop=>hcurl_cop
!---Loop over cells
!$omp do schedule(static)
DO ip=1,mesh%nparts
do ii=1,mesh%tloc_c(ip)%n
  i=mesh%tloc_c(ip)%v(ii)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  he = (6.d0*SQRT(2.d0)*mesh%cv(i))**(1.d0/3.d0)/REAL(oft_lagrange%order,8) ! Cell node spacing
  call oft_lagrange%ncdofs(i,j_lag)
  call oft_hgrad%ncdofs(i,j_hgrad)
  call oft_hcurl%ncdofs(i,j_hcurl)
!------------------------------------------------------------------------------
! Init local entries
!------------------------------------------------------------------------------
  bc_loc=0.d0; bg_loc=0.d0; lag_loc=0.d0
  IF(xmhd_two_temp)te_loc=0.d0
  IF(j2_ind>0)j2_loc=0.d0
  IF(n2_ind>0)n2_loc=0.d0
  !---Lagrange elements
  IF(back_step)THEN
    do jr=1,oft_lagrange%nce
      vpin_loc(jr,1) = npin(j_lag(jr))
      IF(xmhd_upwind)vpin_loc(jr,2:4) = vpin(j_lag(jr),:)
    end do
    !---H(Curl) elements
    do jr=1,oft_hcurl%nce
      jpin_loc(jr) = jpin(j_hcurl(jr))
    end do
  END IF
!------------------------------------------------------------------------------
! Quadrature Loop
!------------------------------------------------------------------------------
  DO m=1,quad%np
!------------------------------------------------------------------------------
! Get local reconstructed operators
!------------------------------------------------------------------------------
    if(curved.OR.m==1)then
      call mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*quad%wts(m)
    !---Get local reconstruction operators
    CALL oft_lag_eval_all(oft_lagrange,i,quad%pts(:,m),lag_rop)
    CALL oft_lag_geval_all(oft_lagrange,i,quad%pts(:,m),lag_gop,goptmp)
    CALL oft_h1_geval_all(oft_hgrad,i,quad%pts(:,m),hgrad_rop,goptmp)
    CALL oft_hcurl_eval_all(oft_hcurl,i,quad%pts(:,m),hcurl_rop,goptmp)
    CALL oft_hcurl_ceval_all(oft_hcurl,i,quad%pts(:,m),hcurl_cop,cgop)
    !---Evaluate centering values
    call full_interp%interp(i,quad%pts(:,m),goptmp,fulltmp)
    CALL xmhd_interp_unpack(fulltmp,u0)
    divv = u0%dV(1,1)+u0%dV(2,2)+u0%dV(3,3)
    bmag = SQRT(SUM(u0%B**2)) + xmhd_eps
    bhat = u0%B/bmag
    !---Reconstruct previous step fields if necessary
    vad=u0%V; vad_mag = magnitude(vad)
    np0q=u0%N; jp0q=u0%J
    IF(back_step)THEN
      np0q=0.d0; jp0q=0.d0
      do jr=1,oft_lagrange%nce
        np0q = np0q + lag_rop(jr)*vpin_loc(jr,1)
      end do
      np0q = np0q*den_scale
      IF(xmhd_upwind)THEN
        vad=0.d0
        do jr=1,oft_lagrange%nce
          vad = vad + lag_rop(jr)*vpin_loc(jr,2:4)
        end do
        vad=vad*vel_scale; vad_mag = magnitude(vad)
      END IF
      do jr=1,oft_hcurl%nce
        jp0q = jp0q + jpin_loc(jr)*hcurl_cop(:,jr)
      end do
    END IF
    !---Handle floors
    IF(u0%N<0.d0)neg_vols(1)=neg_vols(1)+det
    IF(u0%Ti<0.d0)neg_vols(2)=neg_vols(2)+det
    IF(u0%Te<0.d0)neg_vols(3)=neg_vols(3)+det
    IF((den_floor>0.d0).AND.(u0%N<den_floor))neg_flag(1,i)=MAX(neg_flag(1,i),(den_floor-u0%N)/den_floor)
    IF((temp_floor>0.d0).AND.(u0%Ti<temp_floor))neg_flag(2,i)=MAX(neg_flag(2,i),(temp_floor-u0%Ti)/temp_floor)
    IF((temp_floor>0.d0).AND.(u0%Te<temp_floor))neg_flag(3,i)=MAX(neg_flag(3,i),(temp_floor-u0%Te)/temp_floor)
    floor_tmp=(/u0%N,u0%Ti,u0%Te/)
    IF(temp_floor>0.d0)u0%Ti=MAX(temp_floor,u0%Ti)
    IF(temp_floor>0.d0)u0%Te=MAX(temp_floor,u0%Te)
    IF(den_floor>0.d0)u0%N=MAX(den_floor,u0%N)
    IF(den_floor>0.d0)np0q=MAX(den_floor,np0q)
    !---Transport coefficients (if fixed)
    eta_curr=eta*eta_reg(mesh%reg(i))
    IF(.NOT.xmhd_brag)THEN
      chi_temp(1)=kappa_par*den_scale/u0%N
      chi_temp(2)=kappa_perp*den_scale/u0%N
    END IF
    !---Compute diagnostics
    IF(self%dt<0.d0)THEN
      pt=mesh%log2phys(i,quad%pts(:,m))
      pt_r = MAX(1.d-10,SUM(pt(ptind)**2))
      diag_vals(1) = diag_vals(1) + SUM(u0%B**2)*det
      diag_vals(2) = diag_vals(2) + u0%N*SUM(u0%V**2)*det
      uvec=cross_product(pt,u0%B)
      diag_vals(3) = diag_vals(3) + uvec(xmhd_taxis)*det/pt_r
      uvec=cross_product(pt,u0%J)
      diag_vals(4) = diag_vals(4) + uvec(xmhd_taxis)*det/pt_r
      diag_vals(5) = diag_vals(5) + u0%N*det
      diag_vals(6) = diag_vals(6) + u0%Ti*det
      diag_vals(7) = diag_vals(7) + u0%Te*det
    END IF
    !------------------------------------------------------------------------------
    ! Compute F (u) for solid region
    !------------------------------------------------------------------------------
    IF(solid_cell(i))THEN
      IF(xmhd_adv_b)THEN
        !---Eta J
        eta_curr=eta_reg(mesh%reg(i))
        uvec = eta_curr*u0%J*self%dt*det
        do jr=1,oft_hcurl%nce
          bc_loc(jr) = bc_loc(jr) + DOT_PRODUCT(hcurl_cop(:,jr),uvec) &
            + DOT_PRODUCT(hcurl_rop(:,jr),u0%B)*det
        end do
        do jr=1,oft_hgrad%nce
          bg_loc(jr) = bg_loc(jr) + DOT_PRODUCT(hgrad_rop(:,jr),u0%B)*det
        end do
      END IF
      CYCLE
    END IF
!------------------------------------------------------------------------------
! Compute F_bc (u)
!------------------------------------------------------------------------------
    IF(xmhd_adv_b)THEN
      IF(eta_temp>0.d0)THEN
        IF(ASSOCIATED(res_profile))THEN
          eta_curr=eta*res_profile(u0%Te)
        ELSE
          eta_curr=eta*(eta_temp/u0%Te)**(3.d0/2.d0)
        END IF
      END IF
      !---Eta J - VxB
      uvec = eta_curr*u0%J  + cross_product(u0%B,u0%V)
      IF((self%dt>0.d0).OR.xmhd_diss_centered)THEN
        c1 = 1.d0
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0
        uvec = uvec - eta_hyper*U0%J2c*c1
      END IF
      !---Hall (JxB - grad Pe)
      IF(xmhd_hall)THEN
        uvec = uvec + cross_product(u0%J,u0%B)/(mu0*u0%N*elec_charge) &
          - k_boltz*u0%Te*u0%dN/(u0%N*elec_charge)
      END IF
      uvec = uvec*self%dt
      !---Hall (dJ/dt)
      IF(xmhd_hall)THEN
        uvec = uvec + 0.5d0*u0%J*me_factor*elec_mass/(mu0*u0%N*elec_charge**2)
        IF(back_step)THEN
          uvec = uvec + 0.5d0*u0%J*me_factor*elec_mass/(mu0*np0q*elec_charge**2) &
            - 0.5d0*jp0q*me_factor*elec_mass/(mu0*u0%N*elec_charge**2)
        END IF
      END IF
      uvec = uvec*det
      !---Add contributions
      do jr=1,oft_hcurl%nce
        bc_loc(jr) = bc_loc(jr) + DOT_PRODUCT(hcurl_cop(:,jr),uvec) &
          + DOT_PRODUCT(hcurl_rop(:,jr),u0%B)*det
      end do
      !---Hyper-resistivity
      IF((self%dt>0.d0).AND.(j2_ind>0))THEN
        do jr=1,oft_hcurl%nce
          j2_loc(jr) = j2_loc(jr) + (DOT_PRODUCT(hcurl_rop(:,jr),u0%J2) &
            + DOT_PRODUCT(hcurl_cop(:,jr),u0%J))*det/j2_scale
        end do
      END IF
!------------------------------------------------------------------------------
! Compute F_bg (u)
!------------------------------------------------------------------------------
      uvec=u0%B*det
      do jr=1,oft_hgrad%nce
        bg_loc(jr) = bg_loc(jr) + DOT_PRODUCT(hgrad_rop(:,jr),uvec)
      end do
    END IF ! Advance B
!------------------------------------------------------------------------------
! Compute F_v (u)
!------------------------------------------------------------------------------
    uvec = 0.d0; umat = 0.d0
    !---JxB
    IF(xmhd_jcb)uvec = uvec - cross_product(u0%J,u0%B)/(u0%N*m_ion*mu0)
    !---Gravity
    uvec(3) = uvec(3) - g_accel
    !---Advection
    IF(xmhd_advec)uvec = uvec + MATMUL(u0%V,u0%dV)
    !---Pressure
    IF(xmhd_adv_den)uvec = uvec + u0%dN*(u0%Ti+u0%Te)*k_boltz/(u0%N*m_ion)
    IF(xmhd_adv_temp)uvec = uvec + (u0%dTi+u0%dTe)*k_boltz/m_ion
    !---Viscosity
    IF((self%dt>0.d0).OR.xmhd_diss_centered)THEN
      SELECT CASE(visc_itype)
        CASE(1)
          umat = nu_tmp(1)*u0%dV/(u0%N/den_scale)
        CASE(2)
          !--- Nabla(V) + (Nabla(V))^T
          umat = nu_tmp(1)*(u0%dV+TRANSPOSE(u0%dV))/(u0%N/den_scale)
          !--- - 2/3 I Div(V)
          u = 2.d0*nu_tmp(1)*divv/(3.d0*u0%N/den_scale)
          umat(1,1) = umat(1,1) - u
          umat(2,2) = umat(2,2) - u
          umat(3,3) = umat(3,3) - u
        CASE(3)
          !--- Nabla(V) + (Nabla(V))^T
          umat = (u0%dV+TRANSPOSE(u0%dV))/(u0%N/den_scale)
          !--- - 2/3 I Div(V)
          u = 2.d0*divv/(3.d0*u0%N/den_scale)
          umat(1,1) = umat(1,1) - u
          umat(2,2) = umat(2,2) - u
          umat(3,3) = umat(3,3) - u
          !--- Anisotropic tensor (nu_par*bb + nu_perp*(I-bb))
          do j=1,3
            ten1(:,j) = (nu_tmp(1)-nu_tmp(2))*bhat(j)*bhat
            ten1(j,j) = ten1(j,j)+nu_tmp(2)
          end do
          umat = MATMUL(ten1,umat)
      END SELECT
      IF(.NOT.xmhd_diss_centered)umat = 2.d0*umat
      IF(xmhd_adv_den)uvec = uvec - MATMUL(u0%dN,umat)/u0%N
      !---Upwinding
      IF(xmhd_upwind)THEN
        u = xmhd_upwind_weight(vad_mag,he,nu_tmp(2)*den_scale/u0%N)
        IF(.NOT.xmhd_diss_centered)u = 2.d0*u
        umat(:,1) = umat(:,1) + vad*u*DOT_PRODUCT(u0%V,u0%dV(:,1))
        umat(:,2) = umat(:,2) + vad*u*DOT_PRODUCT(u0%V,u0%dV(:,2))
        umat(:,3) = umat(:,3) + vad*u*DOT_PRODUCT(u0%V,u0%dV(:,3))
      END IF
    END IF
    !---Mass and scale
    umat = TRANSPOSE(umat)*self%dt*det/vel_scale
    uvec = (uvec*self%dt + u0%V)*det/vel_scale
    !---Add contributions
    do jr=1,oft_lagrange%nce
      lag_loc(jr,1:3) = lag_loc(jr,1:3) + lag_rop(jr)*uvec &
        + lag_gop(1,jr)*umat(:,1) + lag_gop(2,jr)*umat(:,2) &
        + lag_gop(3,jr)*umat(:,3)
    end do
!------------------------------------------------------------------------------
! Compute F_n (u)
!------------------------------------------------------------------------------
    IF(xmhd_adv_den)THEN
      u = 0.d0; uvec = 0.d0
      !---Compression and advection
      u = u + DOT_PRODUCT(u0%V,u0%dN)
      u = u + u0%N*divv
      !---Density diffusion
      uvec = d_dens*u0%dN
      IF((self%dt>0.d0).OR.xmhd_diss_centered)THEN
        c1 = 1.d0
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0
        uvec = uvec - d2_dens*u0%dN2
      END IF
      !---Upwinding
      IF(xmhd_upwind)THEN
        c1 = xmhd_upwind_weight(vad_mag,he,d_dens)
        uvec = uvec + vad*c1*DOT_PRODUCT(u0%V,u0%dN)
      END IF
      IF(neg_source(1,i)>0.d0)u = u - neg_source(1,i)*den_floor/ABS(4.d0*self%dt)
      !---Mass and scale
      uvec = uvec*self%dt*det/den_scale
      u = (u*self%dt + floor_tmp(1))*det/den_scale
      !---Add contributions
      do jr=1,oft_lagrange%nce
        lag_loc(jr,4) = lag_loc(jr,4) + lag_rop(jr)*u &
          + DOT_PRODUCT(lag_gop(:,jr),uvec)
      end do
      !---Hyper-diffusivity
      IF((self%dt>0.d0).AND.(n2_ind>0))THEN
        u = u0%N2*det/n2_scale
        uvec = u0%dN*det/n2_scale
        do jr=1,oft_lagrange%nce
          n2_loc(jr) = n2_loc(jr) + lag_rop(jr)*u &
            + DOT_PRODUCT(lag_gop(:,jr),uvec)
        end do
      END IF
    END IF
!------------------------------------------------------------------------------
! Compute F_T (u)
!------------------------------------------------------------------------------
    IF(xmhd_adv_temp)THEN
      u = 0.d0; uvec = 0.d0
      !---Viscous heating
      IF(xmhd_visc_heat)THEN
        c1 = m_ion*(temp_gamma-1.d0)/(k_boltz*u0%N/den_scale)
        IF(.NOT.xmhd_two_temp)c1=c1/(1.d0+te_factor)
        SELECT CASE(visc_itype)
          CASE(1)
            u = u - nu_tmp(1)*c1*SUM(SUM(u0%dV*u0%dV,DIM=1))
          CASE(2)
            !--- Nabla(V) + (Nabla(V))^T
            u = u - nu_tmp(1)*c1*SUM(SUM(u0%dV*(u0%dV + TRANSPOSE(u0%dV)),DIM=1))
            !--- - 2/3 I Div(V)
            u = u + 2.d0*nu_tmp(1)*c1*(divv**2)/3.d0
          CASE(3)
            !--- Nabla(V) + (Nabla(V))^T - 2/3 I Div(V)
            umat = TRANSPOSE(u0%dV)+u0%dV
            umat(1,1) = umat(1,1)-2.d0*divv/3.d0
            umat(2,2) = umat(2,2)-2.d0*divv/3.d0
            umat(3,3) = umat(3,3)-2.d0*divv/3.d0
            !--- Anisotropic tensor (nu_par*bb + nu_perp*(I-bb))
            do j=1,3
              ten1(:,j) = (nu_tmp(1)-nu_tmp(2))*bhat(j)*bhat
              ten1(j,j) = ten1(j,j)+nu_tmp(2)
            end do
            ten1 = MATMUL(ten1,umat)
            u = u - c1*SUM(SUM(u0%dV*ten1,DIM=1))
        END SELECT
      END IF
      !---Electron heat sources (Ohmic, Equilibration)
      IF(xmhd_two_temp)THEN
        IF(xmhd_therm_equil)u = u + (u0%Ti-u0%Te)/elec_ion_therm_rate(u0%Te,u0%N,mu_ion)
      ELSE
        IF(xmhd_ohmic)u = u &
          -eta_curr*DOT_PRODUCT(u0%J,u0%J)*(temp_gamma-1.d0)/((1.d0+te_factor)*mu0*u0%N*k_boltz)
      END IF
      !---Compression and advection
      u = u + DOT_PRODUCT(u0%V,u0%dTi)
      u = u + divv*u0%Ti*(temp_gamma-1.d0)
      !---Thermal conduction
      IF((self%dt>0.d0).OR.xmhd_diss_centered)THEN
        IF(xmhd_brag)THEN
          IF(xmhd_two_temp)THEN
            chi_temp=brag_ion_transport(u0%N,u0%Ti,mu_ion,bmag)/u0%N
          ELSE
            chi_temp=brag_comb_transport(u0%N,u0%Ti,mu_ion,bmag)/u0%N
          END IF
          chi_temp(1)=kappa_par*chi_temp(1)
          chi_temp(2)=kappa_perp*chi_temp(2)
        END IF
        c1 = 1.d0
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0
        uvec = uvec + ((chi_temp(1)-chi_temp(2))*bhat*DOT_PRODUCT(bhat,u0%dTi) &
          + chi_temp(2)*u0%dTi)*(temp_gamma-1.d0)*c1
        u = u - ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(u0%dN,bhat)*DOT_PRODUCT(bhat,u0%dTi) &
          + chi_temp(2)*DOT_PRODUCT(u0%dN,u0%dTi))*(temp_gamma-1.d0)*c1/u0%N
        !---Upwinding
        IF(xmhd_upwind)THEN
          c1 = xmhd_upwind_weight(vad_mag,he,chi_temp(2)*(temp_gamma-1.d0))
          IF(.NOT.xmhd_diss_centered)c1 = 2.d0*c1
          uvec = uvec + vad*c1*DOT_PRODUCT(u0%V,u0%dTi)
        END IF
      END IF
      IF(neg_source(2,i)>0.d0)u = u - neg_source(2,i)*temp_floor/ABS(4.d0*self%dt)
      !---Mass and scale
      uvec = uvec*self%dt*det
      u = (u*self%dt + floor_tmp(2))*det
      !---Add contributions
      do jr=1,oft_lagrange%nce
        lag_loc(jr,5) = lag_loc(jr,5) + lag_rop(jr)*u &
          + DOT_PRODUCT(lag_gop(:,jr),uvec)
      end do
    END IF
!------------------------------------------------------------------------------
! Compute F_Te (u)
!------------------------------------------------------------------------------
    IF(xmhd_two_temp.AND.xmhd_adv_temp)THEN
      u = 0.d0; uvec = 0.d0
      !---Collisional Heating
      IF(xmhd_therm_equil)u = u - (u0%Ti-u0%Te)/elec_ion_therm_rate(u0%Te,u0%N,mu_ion)
      !---Ohmic heating
      IF(xmhd_ohmic)u = u - eta_curr*DOT_PRODUCT(u0%J,u0%J)*(temp_gamma-1.d0)/(mu0*u0%N*k_boltz)
      !---Compression and advection
      u = u + DOT_PRODUCT(u0%V-u0%J/(mu0*u0%N*elec_charge),u0%dTe)
      u = u + (divv+DOT_PRODUCT(u0%J,u0%dN)/(mu0*elec_charge*u0%N**2))*u0%Te*(temp_gamma-1.d0)
      !---Thermal conduction
      IF((self%dt>0.d0).OR.xmhd_diss_centered)THEN
        IF(xmhd_brag)THEN
          chi_temp = brag_elec_transport(u0%N,u0%Te,bmag)/u0%N
          chi_temp(1) = kappa_par*chi_temp(1)
          chi_temp(2) = MAX(1.d1, kappa_perp*chi_temp(2))
        END IF
        c1 = 1.d0
        IF(.NOT.xmhd_diss_centered)c1 = 2.d0
        uvec = uvec + ((chi_temp(1)-chi_temp(2))*bhat*DOT_PRODUCT(bhat,u0%dTe) &
          + chi_temp(2)*u0%dTe)*(temp_gamma-1.d0)*c1
        u = u - ((chi_temp(1)-chi_temp(2))*DOT_PRODUCT(u0%dN,bhat)*DOT_PRODUCT(bhat,u0%dTe) &
          + chi_temp(2)*DOT_PRODUCT(u0%dN,u0%dTe))*(temp_gamma-1.d0)*c1/u0%N
        !---Upwinding
        IF(xmhd_upwind)THEN
          IF(back_step)THEN
            vad = vad-jp0q/(mu0*np0q*elec_charge)
          ELSE
            vad = vad-u0%J/(mu0*u0%N*elec_charge)
          END IF
          vad_mag = magnitude(vad)
          c1 = xmhd_upwind_weight(vad_mag,he,chi_temp(2)*(temp_gamma-1.d0))
          IF(.NOT.xmhd_diss_centered)c1 = 2.d0*c1
          uvec = uvec + vad*c1*DOT_PRODUCT(u0%V-u0%J/(mu0*u0%N*elec_charge),u0%dTe)
        END IF
      END IF
      IF(neg_source(3,i)>0.d0)u = u - neg_source(3,i)*temp_floor/ABS(4.d0*self%dt)
      !---Mass and scale
      uvec = uvec*self%dt*det
      u = (u*self%dt + floor_tmp(3))*det
      !---Add contributions
      do jr=1,oft_lagrange%nce
        te_loc(jr) = te_loc(jr) + lag_rop(jr)*u &
          + DOT_PRODUCT(lag_gop(:,jr),uvec)
      end do
    END IF
!------------------------------------------------------------------------------
! End quadrature loop
!------------------------------------------------------------------------------
  end do
!------------------------------------------------------------------------------
! Add local contributions
!------------------------------------------------------------------------------
  !$omp critical(xmhderr_red)
  !---B rows
  IF(xmhd_adv_b)THEN
    do jr=1,oft_hcurl%nce
      bcout(j_hcurl(jr)) = bcout(j_hcurl(jr)) + bc_loc(jr)
      IF(j2_ind>0)j2out(j_hcurl(jr)) = j2out(j_hcurl(jr)) + j2_loc(jr)
    end do
    do jr=1,oft_hgrad%nce
      bgout(j_hgrad(jr)) = bgout(j_hgrad(jr)) + bg_loc(jr)
    end do
  END IF
  !---Lagrange rows (V,n,T)
  IF(.NOT.solid_cell(i))THEN
    DO jr=1,oft_lagrange%nce
      vout(j_lag(jr),:) = vout(j_lag(jr),:) + lag_loc(jr,1:3)
      IF(xmhd_adv_den)nout(j_lag(jr)) = nout(j_lag(jr)) + lag_loc(jr,4)
      IF(xmhd_adv_temp)THEN
        tout(j_lag(jr)) = tout(j_lag(jr)) + lag_loc(jr,5)
        IF(xmhd_two_temp)teout(j_lag(jr)) = teout(j_lag(jr)) + te_loc(jr)
      END IF
      IF(n2_ind>0)n2out(j_lag(jr)) = n2out(j_lag(jr)) + n2_loc(jr)
    END DO
  END IF
  !$omp end critical(xmhderr_red)
end do
end do
!---Delete temporary matrices
deallocate(j_lag,j_hgrad,j_hcurl)
deallocate(lag_rop,lag_gop,hgrad_rop)
deallocate(hcurl_rop,hcurl_cop)
deallocate(bc_loc,bg_loc,lag_loc)
IF(xmhd_two_temp)deallocate(te_loc)
IF(n2_ind>0)DEALLOCATE(n2_loc)
IF(j2_ind>0)DEALLOCATE(j2_loc)
IF(back_step)deallocate(vpin_loc,jpin_loc)
NULLIFY(xmhd_lag_rop,xmhd_lag_gop,xmhd_hgrad_rop,xmhd_hcurl_rop,xmhd_hcurl_cop)
!$omp end parallel
!------------------------------------------------------------------------------
! Boundary conditions
!------------------------------------------------------------------------------
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_hcurl,full_interp%bcurl_loc,bcout,oft_xmhd_ops%b1_bc)
CALL fem_dirichlet_vec(oft_hgrad,full_interp%bgrad_loc,bgout,oft_xmhd_ops%b2_bc)
CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,4),nout,oft_xmhd_ops%n_bc)
CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,5),tout,oft_xmhd_ops%t_bc)
IF(xmhd_two_temp)CALL fem_dirichlet_vec(oft_lagrange,full_interp%Te_loc,teout,oft_xmhd_ops%t_bc)
!---Velocity BC
IF(xmhd_vbcdir)THEN
  !$omp parallel do private(j,mloc)
  DO i=1,oft_lagrange%nbe
    j=oft_lagrange%lbe(i)
    IF(.NOT.oft_lagrange%global%gbe(j))CYCLE
    CALL lag_vbc_tensor(oft_lagrange,j,vbc_type,mloc)
    vout(j,:)=MATMUL(mloc,vout(j,:))
    IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
    CALL lag_vbc_diag(oft_lagrange,j,vbc_type,mloc)
    vout(j,:)=vout(j,:)+MATMUL(mloc,full_interp%lf_loc(j,1:3))
  END DO
ELSE
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,1),vout(:,1),oft_xmhd_ops%v_bc)
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,2),vout(:,2),oft_xmhd_ops%v_bc)
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,3),vout(:,3),oft_xmhd_ops%v_bc)
END IF
!------------------------------------------------------------------------------
! Set result from local field values, summing contributions across seams
!------------------------------------------------------------------------------
CALL b%restore_local(bcout,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(bgout,2,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,1)
CALL b%restore_local(vtmp,3,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,2)
CALL b%restore_local(vtmp,4,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,3)
CALL b%restore_local(vtmp,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(nout,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(tout,7,add=.TRUE.,wait=(7/=xmhd_rep%nfields))
IF(xmhd_two_temp)CALL b%restore_local(teout,te_ind,add=.TRUE.,wait=(te_ind/=xmhd_rep%nfields))
IF(n2_ind>0)CALL b%restore_local(n2out,n2_ind,add=.TRUE.,wait=(n2_ind/=xmhd_rep%nfields))
IF(j2_ind>0)CALL b%restore_local(j2out,j2_ind,add=.TRUE.,wait=(j2_ind/=xmhd_rep%nfields))
!---Pack diagnostics
IF(self%dt<0.d0)THEN
  diag_vals(3:4)=diag_vals(3:4)/(2*pi)
  self%diag_vals=oft_mpi_sum(diag_vals,7)
END IF
!------------------------------------------------------------------------------
! Destroy worker arrays
!------------------------------------------------------------------------------
DEALLOCATE(bcout,bgout,vout,nout,tout)
IF(xmhd_two_temp)DEALLOCATE(teout)
IF(n2_ind>0)DEALLOCATE(n2out)
IF(j2_ind>0)DEALLOCATE(j2out)
IF(back_step)THEN
  DEALLOCATE(jpin,npin)
  IF(xmhd_upwind)DEALLOCATE(vpin)
END IF
CALL full_interp%delete
!----
IF(neg_vols(1)>0.d0)WRITE(*,'(A,I6,A,ES11.3)') '[',oft_env%rank,'] NEG Ne volume   = ',REAL(neg_vols(1),4)
IF(neg_vols(2)>0.d0)WRITE(*,'(A,I6,A,ES11.3)')'[',oft_env%rank,'] NEG Ti volume   = ',REAL(neg_vols(2),4)
IF(xmhd_two_temp)THEN
  IF(neg_vols(3)>0.d0)WRITE(*,'(A,I6,A,ES11.3)')'[',oft_env%rank,'] NEG Te volume   = ',REAL(neg_vols(3),4)
END IF
DEBUG_STACK_POP
end subroutine xmhd_errmatrix_apply
!------------------------------------------------------------------------------
!> Compute the mass matrix for solution field
!!
!! b = M(a)
!------------------------------------------------------------------------------
subroutine xmhd_massmatrix_apply(self,a,b)
class(oft_xmhd_massmatrix), intent(inout) :: self !< Mass matrix object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
real(r8), allocatable :: bc_loc(:),bg_loc(:),te_loc(:),lag_loc(:,:)
integer(i4) :: i,ii,ip,j,m,jr,jc,jp,ptind(2)
integer(i4), allocatable :: j_lag(:),j_hgrad(:),j_hcurl(:)
real(r8), target, allocatable :: lag_rop(:),hgrad_rop(:,:),hcurl_rop(:,:),hcurl_cop(:,:)
real(r8), pointer, contiguous, dimension(:) :: neqin,bcout,bgout,nout,tout,teout
real(r8), pointer, contiguous, dimension(:,:) :: vout
real(r8), pointer, dimension(:) :: vtmp
real(r8) :: mloc(3,3),goptmp(3,4),vol,f,det,cgop(3,6),fulltmp(40)
real(r8) :: neq,diag_vals(7),pt(3),pt_r
logical :: curved,neq_present
type(oft_quad_type), pointer :: quad
TYPE(xmhd_interp) :: full_interp
type(xmhd_loc_values) :: u0
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2X,A)')'Apply xMHD mass matrix'
quad=>oft_hcurl%quad
!------------------------------------------------------------------------------
! Get local field values from source
!------------------------------------------------------------------------------
full_interp%u=>a
CALL full_interp%setup(mesh)
NULLIFY(neqin)
IF(ASSOCIATED(self%u0))THEN
  neq_present=.TRUE.
  CALL self%u0%get_local(neqin,6)
ELSE
  neq_present=.FALSE.
END IF
!------------------------------------------------------------------------------
! Get local field values for result
!------------------------------------------------------------------------------
NULLIFY(bcout,bgout,nout,tout,teout)
CALL b%set(0.d0)
CALL b%get_local(bcout,1)
CALL b%get_local(bgout,2)
ALLOCATE(vout(oft_lagrange%ne,3))
vtmp=>vout(:,1)
CALL b%get_local(vtmp,3)
vtmp=>vout(:,2)
CALL b%get_local(vtmp,4)
vtmp=>vout(:,3)
CALL b%get_local(vtmp,5)
CALL b%get_local(nout,6)
CALL b%get_local(tout,7)
IF(xmhd_two_temp)CALL b%get_local(teout,te_ind)
!---Zero diagnostic values
diag_vals=0.d0
ptind=(/1,2/)
IF(xmhd_taxis==1)ptind(1)=3
IF(xmhd_taxis==2)ptind(2)=3
!------------------------------------------------------------------------------
! Apply metric function
!------------------------------------------------------------------------------
!$omp parallel private(det,j_lag,j_hgrad,j_hcurl,lag_rop,hgrad_rop,hcurl_rop, &
!$omp hcurl_cop,u0,curved,goptmp,vol,j,m,jr,neq,pt,pt_r,cgop,bc_loc,bg_loc,te_loc, &
!$omp lag_loc,i,ii,fulltmp) reduction(+:diag_vals)
allocate(j_lag(oft_lagrange%nce),j_hgrad(oft_hgrad%nce),j_hcurl(oft_hcurl%nce))
allocate(lag_rop(oft_lagrange%nce),hgrad_rop(3,oft_hgrad%nce),hcurl_rop(3,oft_hcurl%nce))
allocate(hcurl_cop(3,oft_hcurl%nce))
allocate(bc_loc(oft_hcurl%nce),bg_loc(oft_hgrad%nce))
allocate(lag_loc(oft_lagrange%nce,5))
IF(xmhd_two_temp)allocate(te_loc(oft_lagrange%nce))
!---Setup caching
xmhd_lag_rop=>lag_rop
ALLOCATE(xmhd_lag_gop(3,oft_lagrange%nce))
xmhd_lag_gop=0.d0
xmhd_hgrad_rop=>hgrad_rop
xmhd_hcurl_rop=>hcurl_rop
xmhd_hcurl_cop=>hcurl_cop
!---Loop over cells
!$omp do schedule(static)
DO ip=1,mesh%nparts
do ii=1,mesh%tloc_c(ip)%n
  i=mesh%tloc_c(ip)%v(ii)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j_lag)
  call oft_hgrad%ncdofs(i,j_hgrad)
  call oft_hcurl%ncdofs(i,j_hcurl)
  !---Init local entries
  bc_loc=0.d0; bg_loc=0.d0; lag_loc=0.d0
  IF(xmhd_two_temp)te_loc=0.d0
!------------------------------------------------------------------------------
! Quadrature Loop
!------------------------------------------------------------------------------
  DO m=1,quad%np
!------------------------------------------------------------------------------
! Get local reconstructed operators
!------------------------------------------------------------------------------
    if(curved.OR.m==1)then
      call mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    pt=mesh%log2phys(i,quad%pts(:,m))
    det=vol*quad%wts(m)
    !---Lagrange elements
    CALL oft_lag_eval_all(oft_lagrange,i,quad%pts(:,m),lag_rop)
    CALL oft_h1_geval_all(oft_hgrad,i,quad%pts(:,m),hgrad_rop,goptmp)
    CALL oft_hcurl_eval_all(oft_hcurl,i,quad%pts(:,m),hcurl_rop,goptmp)
    CALL oft_hcurl_ceval_all(oft_hcurl,i,quad%pts(:,m),hcurl_cop,cgop)
    !
    call full_interp%interp(i,quad%pts(:,m),goptmp,fulltmp)
    CALL xmhd_interp_unpack(fulltmp,u0)
    !---Compute diagnostics
    pt_r = MAX(1.d-10,SUM(pt(ptind)**2))
    diag_vals(1) = diag_vals(1) + SUM(u0%B**2)*det
    mloc(:,1)=cross_product(pt,u0%B)
    diag_vals(3) = diag_vals(3) + mloc(xmhd_taxis,1)*det/pt_r
    mloc(:,1)=cross_product(pt,u0%J)
    diag_vals(4) = diag_vals(4) + mloc(xmhd_taxis,1)*det/pt_r
    IF(neq_present)THEN
      neq=0.d0
      do jr=1,oft_lagrange%nce
        neq = neq + lag_rop(jr)*neqin(j_lag(jr))
      end do
      neq = neq*den_scale
      diag_vals(2) = diag_vals(2) + neq*SUM(u0%V**2)*det
    ELSE
      neq=u0%N
      diag_vals(2) = diag_vals(2) + u0%N*SUM(u0%V**2)*det
    END IF
    diag_vals(5) = diag_vals(5) + u0%N*det
    diag_vals(6) = diag_vals(6) + u0%Ti*det
    diag_vals(7) = diag_vals(7) + u0%Te*det
!------------------------------------------------------------------------------
! Compute mass matrix
!------------------------------------------------------------------------------
    !---Compute M_bc (u)
    do jr=1,oft_hcurl%nce
      bc_loc(jr) = bc_loc(jr) + DOT_PRODUCT(hcurl_rop(:,jr),u0%B)*det
    end do
    !---Electron inertia
    IF(xmhd_hall.AND.(.NOT.solid_cell(i)))THEN
      u0%J=u0%J*dt_scale*me_factor*elec_mass*det/(mu0*neq*elec_charge**2)
      do jr=1,oft_hcurl%nce
        bc_loc(jr) = bc_loc(jr) + DOT_PRODUCT(hcurl_cop(:,jr),u0%J)
      end do
    END IF
    !---Compute M_bg (u)
    do jr=1,oft_hgrad%nce
      bg_loc(jr) = bg_loc(jr) + DOT_PRODUCT(hgrad_rop(:,jr),u0%B)*det
    end do
    IF(solid_cell(i))CYCLE
    !---Compute F_{V,n,T) (u)
    do jr=1,oft_lagrange%nce
      lag_loc(jr,1:3) = lag_loc(jr,1:3) + lag_rop(jr)*u0%V*det
      IF(xmhd_adv_den)lag_loc(jr,4) = lag_loc(jr,4) + lag_rop(jr)*u0%N*det
      IF(xmhd_adv_temp)lag_loc(jr,5) = lag_loc(jr,5) + lag_rop(jr)*u0%Ti*det
      IF(xmhd_two_temp)te_loc(jr) = te_loc(jr) + lag_rop(jr)*u0%Te*det
    end do
!------------------------------------------------------------------------------
! End quadrature loop
!------------------------------------------------------------------------------
  end do
  !$omp critical(xmhderr_red)
!------------------------------------------------------------------------------
! Add local contributions
!------------------------------------------------------------------------------
  !---B rows
  do jr=1,oft_hcurl%nce
    bcout(j_hcurl(jr)) = bcout(j_hcurl(jr)) + bc_loc(jr)
  end do
  do jr=1,oft_hgrad%nce
    bgout(j_hgrad(jr)) = bgout(j_hgrad(jr)) + bg_loc(jr)
  end do
  !---Lagrange rows (V,n,T)
  do jr=1,oft_lagrange%nce
    vout(j_lag(jr),:) = vout(j_lag(jr),:) + lag_loc(jr,1:3)/vel_scale
    IF(xmhd_adv_den)nout(j_lag(jr)) = nout(j_lag(jr)) + lag_loc(jr,4)/den_scale
    IF(xmhd_adv_temp)tout(j_lag(jr)) = tout(j_lag(jr)) + lag_loc(jr,5)
    IF(xmhd_two_temp)teout(j_lag(jr)) = teout(j_lag(jr)) + te_loc(jr)
  end do
  !$omp end critical(xmhderr_red)
end do
end do
!---Delete temporary matrices
deallocate(j_lag,j_hgrad,j_hcurl)
deallocate(lag_rop,hgrad_rop,hcurl_rop)
deallocate(bc_loc,bg_loc,lag_loc)
IF(xmhd_two_temp)deallocate(te_loc)
NULLIFY(xmhd_lag_rop,xmhd_lag_gop,xmhd_hgrad_rop,xmhd_hcurl_rop,xmhd_hcurl_cop)
!$omp end parallel
!------------------------------------------------------------------------------
! Boundary conditions
!------------------------------------------------------------------------------
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_hcurl,full_interp%bcurl_loc,bcout,oft_xmhd_ops%b1_bc)
CALL fem_dirichlet_vec(oft_hgrad,full_interp%bgrad_loc,bgout,oft_xmhd_ops%b2_bc)
CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,4),nout,oft_xmhd_ops%n_bc)
CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,5),tout,oft_xmhd_ops%t_bc)
IF(xmhd_two_temp)CALL fem_dirichlet_vec(oft_lagrange,full_interp%Te_loc,teout,oft_xmhd_ops%t_bc)
!------------------------------------------------------------------------------
! Velocity BC
!------------------------------------------------------------------------------
IF(xmhd_vbcdir)THEN
  !---Copy values for boundary terms
  !$omp parallel do private(j,mloc)
  DO i=1,oft_lagrange%nbe
    j=oft_lagrange%lbe(i)
    IF(.NOT.oft_lagrange%global%gbe(j))CYCLE
    CALL lag_vbc_tensor(oft_lagrange,j,vbc_type,mloc)
    vout(j,:)=MATMUL(mloc,vout(j,:))
    IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
    CALL lag_vbc_diag(oft_lagrange,j,vbc_type,mloc)
    vout(j,:)=vout(j,:)+MATMUL(mloc,full_interp%lf_loc(j,1:3))
  END DO
ELSE
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,1),vout(:,1),oft_xmhd_ops%v_bc)
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,2),vout(:,2),oft_xmhd_ops%v_bc)
  CALL fem_dirichlet_vec(oft_lagrange,full_interp%lf_loc(:,3),vout(:,3),oft_xmhd_ops%v_bc)
END IF
!------------------------------------------------------------------------------
! Set result from local field values, summing contributions across seams
!------------------------------------------------------------------------------
CALL b%restore_local(bcout,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(bgout,2,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,1)
CALL b%restore_local(vtmp,3,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,2)
CALL b%restore_local(vtmp,4,add=.TRUE.,wait=.TRUE.)
vtmp=>vout(:,3)
CALL b%restore_local(vtmp,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(nout,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(tout,7,add=.TRUE.,wait=(7/=xmhd_rep%nfields))
IF(xmhd_two_temp)CALL b%restore_local(teout,te_ind,add=.TRUE.,wait=(te_ind/=xmhd_rep%nfields))
tout=0.d0
IF(n2_ind>0)CALL b%restore_local(tout,n2_ind,add=.TRUE.,wait=(n2_ind/=xmhd_rep%nfields))
bcout=0.d0
IF(j2_ind>0)CALL b%restore_local(bcout,j2_ind,add=.TRUE.,wait=(j2_ind/=xmhd_rep%nfields))
!---Pack diagnostics
diag_vals(2)=diag_vals(2)
diag_vals(3:4)=diag_vals(3:4)/(2*pi)
diag_vals(5)=diag_vals(5)
self%diag_vals=oft_mpi_sum(diag_vals,7)
!------------------------------------------------------------------------------
! Destroy worker arrays
!------------------------------------------------------------------------------
DEALLOCATE(bcout,bgout,vout,nout,tout)
IF(ASSOCIATED(neqin))DEALLOCATE(neqin)
IF(xmhd_two_temp)DEALLOCATE(teout)
CALL full_interp%delete
DEBUG_STACK_POP
end subroutine xmhd_massmatrix_apply
!------------------------------------------------------------------------------
!> Compute diagnostic values from fields
!------------------------------------------------------------------------------
subroutine xmhd_diag(u,diag_vals,ueq)
class(oft_vector), target, intent(inout) :: u !< Field for diagnostic evaluation
real(r8), intent(out) :: diag_vals(7) !< Resulting diagnostic values [7]
class(oft_vector), optional, intent(inout) :: ueq !< Equilibrium field (linear calculation)
integer(i4) :: i,ii,ip,j,m,jr,jc,jp,ptind(2)
integer(i4), allocatable :: j_lag(:),j_hgrad(:),j_hcurl(:)
real(r8), target, allocatable :: lag_rop(:),hgrad_rop(:,:),hcurl_rop(:,:),hcurl_cop(:,:)
real(r8), pointer, contiguous, dimension(:) :: neqin
real(r8), pointer, contiguous, dimension(:,:) :: vout
real(r8), pointer, dimension(:) :: vtmp
real(r8) :: mloc(3,3),goptmp(3,4),vol,f,det,cgop(3,6),fulltmp(40)
real(r8) :: neq,pt(3),pt_r
logical :: curved,neq_present
type(oft_quad_type), pointer :: quad
TYPE(xmhd_interp) :: full_interp
type(xmhd_loc_values) :: u0
DEBUG_STACK_PUSH
IF(oft_debug_print(1))write(*,'(2X,A)')'Apply xMHD mass matrix'
quad=>oft_hcurl%quad
!------------------------------------------------------------------------------
! Get local field values from source
!------------------------------------------------------------------------------
full_interp%u=>u
CALL full_interp%setup(mesh)
NULLIFY(neqin)
IF(PRESENT(ueq))THEN
  neq_present=.TRUE.
  CALL ueq%get_local(neqin,6)
ELSE
  neq_present=.FALSE.
END IF
!---Zero diagnostic values
diag_vals=0.d0
ptind=(/1,2/)
IF(xmhd_taxis==1)ptind(1)=3
IF(xmhd_taxis==2)ptind(2)=3
!------------------------------------------------------------------------------
! Apply metric function
!------------------------------------------------------------------------------
!$omp parallel private(det,j_lag,j_hgrad,j_hcurl,lag_rop,hgrad_rop,hcurl_rop, &
!$omp hcurl_cop,curved,goptmp,vol,j,m,jr,neq,pt,pt_r,cgop,u0,i,ii,fulltmp) &
!$omp reduction(+:diag_vals)
allocate(j_lag(oft_lagrange%nce),j_hgrad(oft_hgrad%nce),j_hcurl(oft_hcurl%nce))
allocate(lag_rop(oft_lagrange%nce),hgrad_rop(3,oft_hgrad%nce),hcurl_rop(3,oft_hcurl%nce))
allocate(hcurl_cop(3,oft_hcurl%nce))
!---Setup caching
xmhd_lag_rop=>lag_rop
ALLOCATE(xmhd_lag_gop(3,oft_lagrange%nce))
xmhd_lag_gop=0.d0
xmhd_hgrad_rop=>hgrad_rop
xmhd_hcurl_rop=>hcurl_rop
xmhd_hcurl_cop=>hcurl_cop
!---Loop over cells
!$omp do schedule(static)
DO ip=1,mesh%nparts
do ii=1,mesh%tloc_c(ip)%n
  i=mesh%tloc_c(ip)%v(ii)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j_lag)
  call oft_hgrad%ncdofs(i,j_hgrad)
  call oft_hcurl%ncdofs(i,j_hcurl)
!------------------------------------------------------------------------------
! Quadrature Loop
!------------------------------------------------------------------------------
  DO m=1,quad%np
!------------------------------------------------------------------------------
! Get local reconstructed operators
!------------------------------------------------------------------------------
    if(curved.OR.m==1)then
      call mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    pt=mesh%log2phys(i,quad%pts(:,m))
    det=vol*quad%wts(m)
    !---Lagrange elements
    CALL oft_lag_eval_all(oft_lagrange,i,quad%pts(:,m),lag_rop)
    CALL oft_h1_geval_all(oft_hgrad,i,quad%pts(:,m),hgrad_rop,goptmp)
    CALL oft_hcurl_eval_all(oft_hcurl,i,quad%pts(:,m),hcurl_rop,goptmp)
    CALL oft_hcurl_ceval_all(oft_hcurl,i,quad%pts(:,m),hcurl_cop,cgop)
    !
    call full_interp%interp(i,quad%pts(:,m),goptmp,fulltmp)
    CALL xmhd_interp_unpack(fulltmp,u0)
    !---Compute diagnostics
    pt_r = MAX(1.d-10,SUM(pt(ptind)**2))
    diag_vals(1) = diag_vals(1) + SUM(u0%B**2)*det
    mloc(:,1)=cross_product(pt,u0%B)
    diag_vals(3) = diag_vals(3) + mloc(xmhd_taxis,1)*det/pt_r
    mloc(:,1)=cross_product(pt,u0%J)
    diag_vals(4) = diag_vals(4) + mloc(xmhd_taxis,1)*det/pt_r
    IF(neq_present)THEN
      neq=0.d0
      do jr=1,oft_lagrange%nce
        neq = neq + lag_rop(jr)*neqin(j_lag(jr))
      end do
      neq = neq*den_scale
      diag_vals(2) = diag_vals(2) + neq*SUM(u0%V**2)*det
    ELSE
      neq = u0%N
      diag_vals(2) = diag_vals(2) + u0%N*SUM(u0%V**2)*det
    END IF
    diag_vals(5) = diag_vals(5) + u0%N*det
    diag_vals(6) = diag_vals(6) + u0%Ti*det
    diag_vals(7) = diag_vals(7) + u0%Te*det
  END DO
end do
end do
!---Delete temporary matrices
deallocate(j_lag,j_hgrad,j_hcurl)
deallocate(lag_rop,hgrad_rop,hcurl_rop)
NULLIFY(xmhd_lag_rop,xmhd_lag_gop,xmhd_hgrad_rop,xmhd_hcurl_rop,xmhd_hcurl_cop)
!$omp end parallel
!---Pack diagnostics
diag_vals(2)=diag_vals(2)
diag_vals(3:4)=diag_vals(3:4)/(2*pi)
diag_vals(5)=diag_vals(5)
diag_vals=oft_mpi_sum(diag_vals,7)
!------------------------------------------------------------------------------
! Destroy worker arrays
!------------------------------------------------------------------------------
IF(ASSOCIATED(neqin))DEALLOCATE(neqin)
CALL full_interp%delete
DEBUG_STACK_POP
end subroutine xmhd_diag
!------------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!------------------------------------------------------------------------------
subroutine xmhd_setup_rep
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating xMHD FE type'
!---Create FE representation
ML_xmhd_rep%nlevels=xmhd_ML_hcurl%nlevels
ML_xmhd_rep%nfields=7
IF(xmhd_two_temp)THEN
  ML_xmhd_rep%nfields=ML_xmhd_rep%nfields+1
  te_ind=ML_xmhd_rep%nfields
END IF
IF(d2_dens>0.d0)THEN
  ML_xmhd_rep%nfields=ML_xmhd_rep%nfields+1
  n2_ind=ML_xmhd_rep%nfields
END IF
IF(eta_hyper>0.d0)THEN
  ML_xmhd_rep%nfields=ML_xmhd_rep%nfields+1
  j2_ind=ML_xmhd_rep%nfields
END IF
ALLOCATE(ML_xmhd_rep%ml_fields(ML_xmhd_rep%nfields))
ALLOCATE(ML_xmhd_rep%field_tags(ML_xmhd_rep%nfields))
ML_xmhd_rep%ml_fields(1)%ml=>xmhd_ML_hcurl
ML_xmhd_rep%field_tags(1)='Bc'
ML_xmhd_rep%ml_fields(2)%ml=>xmhd_ML_H1grad
ML_xmhd_rep%field_tags(2)='Bg'
ML_xmhd_rep%ml_fields(3)%ml=>xmhd_ML_lagrange
ML_xmhd_rep%field_tags(3)='Vx'
ML_xmhd_rep%ml_fields(4)%ml=>xmhd_ML_lagrange
ML_xmhd_rep%field_tags(4)='Vy'
ML_xmhd_rep%ml_fields(5)%ml=>xmhd_ML_lagrange
ML_xmhd_rep%field_tags(5)='Vz'
ML_xmhd_rep%ml_fields(6)%ml=>xmhd_ML_lagrange
ML_xmhd_rep%field_tags(6)='n'
ML_xmhd_rep%ml_fields(7)%ml=>xmhd_ML_lagrange
ML_xmhd_rep%field_tags(7)='T'
IF(xmhd_two_temp)THEN
  ML_xmhd_rep%ml_fields(te_ind)%ml=>xmhd_ML_lagrange
  ML_xmhd_rep%field_tags(te_ind)='Te'
END IF
IF(n2_ind>0)THEN
  ML_xmhd_rep%ml_fields(n2_ind)%ml=>xmhd_ML_lagrange
  ML_xmhd_rep%field_tags(n2_ind)='n2'
END IF
IF(j2_ind>0)THEN
  ML_xmhd_rep%ml_fields(j2_ind)%ml=>xmhd_ML_hcurl
  ML_xmhd_rep%field_tags(j2_ind)='j2'
END IF
call ML_xmhd_rep%setup()
xmhd_rep=>ML_xmhd_rep%current_level
!---Declare legacy variables
xmhd_blevel=xmhd_ML_hcurl%blevel
xmhd_nlevels=xmhd_ML_hcurl%nlevels
xmhd_level=xmhd_ML_hcurl%level
IF(xmhd_minlev<0)xmhd_minlev=xmhd_nlevels
end subroutine xmhd_setup_rep
!------------------------------------------------------------------------------
!> Set the current level for xMHD model
!------------------------------------------------------------------------------
subroutine xmhd_set_level(level)
integer(i4), intent(in) :: level !< Desired level
if(level>xmhd_ML_hcurl%nlevels.OR.level<=0)then
  call oft_abort('Invalid FE level','xmhd_set_level',__FILE__)
end if
! if(level<mg_mesh%mgdim)then
!   call multigrid_level(mg_mesh,level)
! else
!   call multigrid_level(mg_mesh,mg_mesh%mgdim)
! end if
CALL ML_xmhd_rep%set_level(level)
xmhd_rep=>ML_xmhd_rep%current_level
!---
CALL xmhd_ML_lagrange%set_level(level)
IF(.NOT.oft_3D_lagrange_cast(oft_lagrange,xmhd_ML_lagrange%current_level))CALL oft_abort("Invalid FE object","xmhd_set_level",__FILE__)
CALL xmhd_ML_hcurl_grad%set_level(level,propogate=.TRUE.)
IF(.NOT.oft_3D_hcurl_cast(oft_hcurl,xmhd_ML_hcurl%current_level))CALL oft_abort("Invalid Curl FE object","xmhd_run",__FILE__)
IF(.NOT.oft_3D_h1_cast(oft_hgrad,xmhd_ML_H1grad%current_level))CALL oft_abort("Invalid Grad FE object","xmhd_run",__FILE__)
xmhd_level=level
! xmhd_lev=oft_hcurl_lev
oft_xmhd_ops=>oft_xmhd_ops_ML(xmhd_level)
end subroutine xmhd_set_level
!------------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new fields
!------------------------------------------------------------------------------
subroutine xmhd_mfnk_update(uin)
class(oft_vector), target, intent(inout) :: uin !< Current field
IF(oft_debug_print(1))write(*,*)'Updating xMHD MF-Jacobian'
CALL mfmat%update(uin)
END SUBROUTINE xmhd_mfnk_update
!------------------------------------------------------------------------------
!> Setup interpolator for xMHD solution fields
!!
!! - Fetches local vector values for interpolation
!------------------------------------------------------------------------------
subroutine xmhd_interp_setup(self,mesh)
CLASS(xmhd_interp), INTENT(inout) :: self !< Interpolation object
class(oft_mesh), target, intent(inout) :: mesh
INTEGER(i4) :: i
REAL(r8), POINTER, DIMENSION(:) :: vtmp
!---Set FEM links
self%grad_rep=>oft_hgrad
self%curl_rep=>oft_hcurl
self%lag_rep=>oft_lagrange
self%mesh=>mesh
!---Get local slice
CALL self%u%get_local(self%bcurl_loc,1)
CALL self%u%get_local(self%bgrad_loc,2)
IF(.NOT.ASSOCIATED(self%lf_loc))THEN
  ALLOCATE(self%lf_loc(self%lag_rep%ne,5))
  IF(xmhd_two_temp)ALLOCATE(self%Te_loc(self%lag_rep%ne))
  IF(j2_ind>0)ALLOCATE(self%J2_loc(self%curl_rep%ne))
  IF(n2_ind>0)ALLOCATE(self%N2_loc(self%lag_rep%ne))
END IF
vtmp=>self%lf_loc(:,1)
CALL self%u%get_local(vtmp,3)
vtmp=>self%lf_loc(:,2)
CALL self%u%get_local(vtmp,4)
vtmp=>self%lf_loc(:,3)
CALL self%u%get_local(vtmp,5)
vtmp=>self%lf_loc(:,4)
CALL self%u%get_local(vtmp,6)
vtmp=>self%lf_loc(:,5)
CALL self%u%get_local(vtmp,7)
IF(xmhd_two_temp)CALL self%u%get_local(self%Te_loc,te_ind)
IF(j2_ind>0)CALL self%u%get_local(self%J2_loc,j2_ind)
IF(n2_ind>0)CALL self%u%get_local(self%N2_loc,n2_ind)
!---Allocate thread cache
IF(.NOT.ASSOCIATED(self%cache))THEN
  ALLOCATE(self%cache(oft_env%nthreads))
  DO i=1,oft_env%nthreads
    ALLOCATE(self%cache(i)%bcurl(self%curl_rep%nce))
    ALLOCATE(self%cache(i)%bgrad(self%grad_rep%nce))
    ALLOCATE(self%cache(i)%lf(self%lag_rep%nce,5))
    IF(xmhd_two_temp)ALLOCATE(self%cache(i)%Te(self%lag_rep%nce))
    IF(j2_ind>0)ALLOCATE(self%cache(i)%J2(self%curl_rep%nce))
    IF(n2_ind>0)ALLOCATE(self%cache(i)%N2(self%lag_rep%nce))
    self%cache(i)%cell=-1
    self%cache(i)%bcurl=0.d0
    self%cache(i)%bgrad=0.d0
    self%cache(i)%lf=0.d0
    IF(xmhd_two_temp)self%cache(i)%Te=0.d0
    IF(j2_ind>0)self%cache(i)%J2=0.d0
    IF(n2_ind>0)self%cache(i)%N2=0.d0
  END DO
END IF
end subroutine xmhd_interp_setup
!------------------------------------------------------------------------------
!> Destroy interpolator for xMHD solution fields
!------------------------------------------------------------------------------
subroutine xmhd_interp_delete(self)
class(xmhd_interp), intent(inout) :: self !< Interpolation object
INTEGER(i4) :: i
!---Destroy locals
IF(ASSOCIATED(self%bgrad_loc))THEN
  DEALLOCATE(self%bgrad_loc,self%bcurl_loc,self%lf_loc)
  IF(xmhd_two_temp)DEALLOCATE(self%Te_loc)
  IF(j2_ind>0)DEALLOCATE(self%J2_loc)
  IF(n2_ind>0)DEALLOCATE(self%N2_loc)
END IF
IF(ASSOCIATED(self%cache))THEN
  DO i=1,oft_env%nthreads
    DEALLOCATE(self%cache(i)%bcurl)
    DEALLOCATE(self%cache(i)%bgrad)
    DEALLOCATE(self%cache(i)%lf)
    IF(xmhd_two_temp)DEALLOCATE(self%cache(i)%Te)
    IF(j2_ind>0)DEALLOCATE(self%cache(i)%J2)
    IF(n2_ind>0)DEALLOCATE(self%cache(i)%N2)
  END DO
  DEALLOCATE(self%cache)
END IF
NULLIFY(self%curl_rep,self%grad_rep,self%lag_rep,self%u)
end subroutine xmhd_interp_delete
!------------------------------------------------------------------------------
!> Reconstruct xMHD solution fields
!------------------------------------------------------------------------------
subroutine xmhd_interp_apply(self,cell,f,gop,val)
class(xmhd_interp), intent(inout) :: self !< Interpolation object
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [30]
integer(i4), allocatable :: j(:)
integer(i4) :: jc,k
real(r8) :: b0(3),j0(3),cgop(3,6),Te,dTe(3),tmp,vtmp(3)
real(r8) :: N2,dN2(3),J2(3),J2c(3),lf0(5),dlf0(3,5)
real(r8), allocatable :: lag_rop(:),lag_gop(:,:)
real(r8), allocatable :: hcurl_rop(:,:),hcurl_cop(:,:),hgrad_rop(:,:)
type(xmhd_interp_cache), pointer :: tloc
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%bgrad_loc))CALL oft_abort('Setup has not been called!', &
  'xmhd_interp_apply',__FILE__)
!---Pull cache
tloc=>self%cache(oft_tid+1)
IF(tloc%cell/=cell)THEN
  allocate(j(self%curl_rep%nce))
  call self%curl_rep%ncdofs(cell,j)
  do jc=1,self%curl_rep%nce
    tloc%bcurl(jc)=self%bcurl_loc(j(jc))
    IF(j2_ind>0)tloc%J2(jc)=self%J2_loc(j(jc))
  end do
  deallocate(j)
  allocate(j(self%grad_rep%nce))
  call self%grad_rep%ncdofs(cell,j)
  do jc=1,self%grad_rep%nce
    tloc%bgrad(jc)=self%bgrad_loc(j(jc))
  end do
  deallocate(j)
  allocate(j(self%lag_rep%nce))
  call self%lag_rep%ncdofs(cell,j)
  do jc=1,self%lag_rep%nce
    tloc%lf(jc,:)=self%lf_loc(j(jc),:)
    IF(xmhd_two_temp)tloc%Te(jc)=self%Te_loc(j(jc))
    IF(n2_ind>0)tloc%N2(jc)=self%N2_loc(j(jc))
  end do
  deallocate(j)
  tloc%cell=cell
END IF
!---
lf0=0.d0; dlf0=0.d0; b0=0.d0; j0=0.d0;
Te=0.d0; dTe=0.d0
N2=0.d0; dN2=0.d0
J2=0.d0; J2c=0.d0
IF(ASSOCIATED(xmhd_lag_rop))THEN
  !---Reconstruct velocity field and shear tensor
  DO k=1,5
    tmp=0.d0
    vtmp=0.d0
    !$omp simd reduction(+:tmp,vtmp)
    DO jc=1,self%lag_rep%nce
      tmp = tmp + tloc%lf(jc,k)*xmhd_lag_rop(jc)
      vtmp = vtmp + tloc%lf(jc,k)*xmhd_lag_gop(:,jc)
    END DO
    lf0(k)=tmp
    dlf0(:,k)=vtmp
  END DO
  IF(xmhd_two_temp)THEN
    !$omp simd reduction(+:Te,dTe)
    DO jc=1,self%lag_rep%nce
      Te = Te + tloc%Te(jc)*xmhd_lag_rop(jc)
      dTe = dTe + tloc%Te(jc)*xmhd_lag_gop(:,jc)
    END DO
  END IF
  IF(n2_ind>0)THEN
    !$omp simd reduction(+:N2,dN2)
    DO jc=1,self%lag_rep%nce
      N2 = N2 + tloc%N2(jc)*xmhd_lag_rop(jc)
      dN2 = dN2 + tloc%N2(jc)*xmhd_lag_gop(:,jc)
    END DO
  END IF
  !---Reconstruct magnetic field
  !$omp simd reduction(+:b0)
  do jc=1,self%grad_rep%nce
    b0 = b0 + tloc%bgrad(jc)*xmhd_hgrad_rop(:,jc)
  end do
  !---Reconstruct magnetic field and current density
  !$omp simd reduction(+:b0,j0)
  do jc=1,self%curl_rep%nce
    b0 = b0 + tloc%bcurl(jc)*xmhd_hcurl_rop(:,jc)
    j0 = j0 + tloc%bcurl(jc)*xmhd_hcurl_cop(:,jc)
  end do
  IF(j2_ind>0)THEN
    !$omp simd reduction(+:J2,J2c)
    do jc=1,self%curl_rep%nce
      J2 = J2 + tloc%J2(jc)*xmhd_hcurl_rop(:,jc)
      J2c = J2c + tloc%J2(jc)*xmhd_hcurl_cop(:,jc)
    end do
  END IF
ELSE
  !---Lagrange elements
  ALLOCATE(lag_rop(self%lag_rep%nce),lag_gop(3,self%lag_rep%nce))
  CALL oft_lag_eval_all(self%lag_rep,cell,f,lag_rop)
  CALL oft_lag_geval_all(self%lag_rep,cell,f,lag_gop,gop)
  DO k=1,5
    tmp=0.d0
    vtmp=0.d0
    !$omp simd reduction(+:tmp,vtmp)
    DO jc=1,self%lag_rep%nce
      tmp = tmp + tloc%lf(jc,k)*lag_rop(jc)
      vtmp = vtmp + tloc%lf(jc,k)*lag_gop(:,jc)
    END DO
    lf0(k)=tmp
    dlf0(:,k)=vtmp
  END DO
  IF(xmhd_two_temp)THEN
    !$omp simd reduction(+:Te,dTe)
    DO jc=1,self%lag_rep%nce
      Te = Te + tloc%Te(jc)*lag_rop(jc)
      dTe = dTe + tloc%Te(jc)*lag_gop(:,jc)
    END DO
  END IF
  IF(n2_ind>0)THEN
    !$omp simd reduction(+:N2,dN2)
    DO jc=1,self%lag_rep%nce
      N2 = N2 + tloc%N2(jc)*lag_rop(jc)
      dN2 = dN2 + tloc%N2(jc)*lag_gop(:,jc)
    END DO
  END IF
  DEALLOCATE(lag_rop,lag_gop)
  !---Grad(H^1) elements
  ALLOCATE(hgrad_rop(3,self%grad_rep%nce))
  CALL oft_h1_geval_all(self%grad_rep,cell,f,hgrad_rop,gop)
  !---Reconstruct magnetic field
  !$omp simd reduction(+:b0)
  do jc=1,self%grad_rep%nce
    b0 = b0 + tloc%bgrad(jc)*hgrad_rop(:,jc)
  end do
  DEALLOCATE(hgrad_rop)
  !---H(Curl) elements
  CALL oft_hcurl_get_cgops(gop,cgop)
  ALLOCATE(hcurl_rop(3,self%curl_rep%nce),hcurl_cop(3,self%curl_rep%nce))
  CALL oft_hcurl_eval_all(self%curl_rep,cell,f,hcurl_rop,gop)
  CALL oft_hcurl_ceval_all(self%curl_rep,cell,f,hcurl_cop,cgop)
  !---Reconstruct magnetic field and current density
  !$omp simd reduction(+:b0,j0)
  do jc=1,self%curl_rep%nce
    b0 = b0 + tloc%bcurl(jc)*hcurl_rop(:,jc)
    j0 = j0 + tloc%bcurl(jc)*hcurl_cop(:,jc)
  end do
  IF(j2_ind>0)THEN
    !$omp simd reduction(+:J2,J2c)
    do jc=1,self%curl_rep%nce
      J2 = J2 + tloc%J2(jc)*hcurl_rop(:,jc)
      J2c = J2c + tloc%J2(jc)*hcurl_cop(:,jc)
    end do
  END IF
  DEALLOCATE(hcurl_rop,hcurl_cop)
END IF
!---Pack
val(1:5)=lf0 ! (Vx, Vy, Vz, N, Ti)
val(6:14)=RESHAPE(dlf0(:,1:3),(/9/)) ! dV_i / dx_j
val(15:17)=dlf0(:,4) ! dN / dx_j
val(18:20)=dlf0(:,5) ! dTi /dx_j
val(21:23)=b0
val(24:26)=j0
IF(xmhd_two_temp)THEN
  val(27) = Te
  val(28:30) = dTe
END IF
IF(n2_ind>0)THEN
  val(31) = N2
  val(32:34) = dN2
END IF
IF(j2_ind>0)THEN
  val(35:37) = J2
  val(38:40) = J2c
END IF
DEBUG_STACK_POP
end subroutine xmhd_interp_apply
!------------------------------------------------------------------------------
!> Unpack interpolation fields into xMHD field structure
!------------------------------------------------------------------------------
subroutine xmhd_interp_unpack(vals,vloc)
real(r8), intent(in) :: vals(40) !< Packed fields [40]
type(xmhd_loc_values), intent(out) :: vloc !< Unpacked fields
!---Unpack
vloc%V = vals(1:3)*vel_scale
vloc%dV = RESHAPE(vals(6:14),(/3,3/))*vel_scale
vloc%N = vals(4)*den_scale; vloc%dN = vals(15:17)*den_scale
vloc%Ti = vals(5); vloc%dTi = vals(18:20)
vloc%B = vals(21:23); vloc%J = vals(24:26)
IF(xmhd_two_temp)THEN
  vloc%Te = vals(27); vloc%dTe = vals(28:30)
ELSE
  vloc%Te = te_factor*vloc%Ti; vloc%dTe = te_factor*vloc%dTi
END IF
IF(n2_ind>0)THEN
  vloc%N2 = vals(31)*n2_scale
  vloc%dN2 = vals(32:34)*n2_scale
END IF
IF(j2_ind>0)THEN
  vloc%J2 = vals(35:37)*j2_scale
  vloc%J2c = vals(38:40)*j2_scale
END IF
end subroutine xmhd_interp_unpack
!------------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new solution
!------------------------------------------------------------------------------
subroutine xmhd_set_ops(uin)
class(oft_vector), target, intent(inout) :: uin !< Current solution
TYPE(xmhd_interp) :: full_interp
TYPE(cc_interp) :: fullcc_interp
REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: fullcc,fullcctmp
integer(i4) :: i,l,k,levelin,minlev,qorder,nfields
type(oft_timer) :: mytimer
IF(xmhd_skip_update)THEN
  xmhd_skip_update=.FALSE.
  RETURN
END IF
IF(oft_debug_print(1))write(*,'(2X,A)')'Update xMHD Jacobians'
!------------------------------------------------------------------------------
! Setup interpolators
!------------------------------------------------------------------------------
NULLIFY(fullcc,fullcctmp)
full_interp%u=>uin
CALL full_interp%setup(mesh)
!------------------------------------------------------------------------------
! Construct Jacobian on each sublevel
!------------------------------------------------------------------------------
levelin=xmhd_level
minlev=xmhd_minlev
qorder=oft_hcurl%order
nfields=26
IF(xmhd_two_temp)nfields=30
do i=levelin,minlev,-1
  call xmhd_set_level(i)
!------------------------------------------------------------------------------
! Level is on finest mesh, use full field representations
!------------------------------------------------------------------------------
  if(oft_lagrange%mesh%nc==mg_mesh%meshes(mg_mesh%mgdim)%nc)then
    !---Construct Jacobian
    call xmhd_build_ops(full_interp)
    !---Average on lowest level to cell centers
    if(oft_lagrange%order==1.AND.i>xmhd_minlev)then
      allocate(fullcc(nfields,mesh%nc))
      call fem_avg_bcc(mesh,full_interp,fullcc,qorder,n=nfields)
    end if
!------------------------------------------------------------------------------
! Level is on coarse mesh, average fields to coarser grid
!------------------------------------------------------------------------------
  else
    !---Create cell centered storage
    if(ASSOCIATED(fullcctmp))DEALLOCATE(fullcctmp)
    allocate(fullcctmp(nfields,mesh%nc))
    !---Transfer from distributed to local grid
    if(xmhd_ML_lagrange%level==xmhd_ML_lagrange%blevel)then
      call multigrid_base_pushcc(mg_mesh,fullcc,fullcctmp,nfields)
    !---Average values to over child cells
    else
      !$omp parallel do private(k)
      do l=1,mesh%nc
        fullcctmp(:,l)=0.d0
        do k=1,8
          fullcctmp(:,l)=fullcctmp(:,l)+fullcc(:,8*(l-1)+k)/8.d0
        end do
      end do
    end if
    !---Switch storage pointers
    DEALLOCATE(fullcc)
    fullcc=>fullcctmp
    NULLIFY(fullcctmp)
    !---Setup interpolation
    fullcc_interp%bcc=>fullcc
    !---Construct Jacobian
    call xmhd_build_ops(fullcc_interp)
  end if
end do
!------------------------------------------------------------------------------
! Cleanup temporary storage
!------------------------------------------------------------------------------
if(ASSOCIATED(fullcc))deallocate(fullcc)
if(ASSOCIATED(fullcctmp))deallocate(fullcctmp)
call full_interp%delete()
call xmhd_set_level(levelin)
end subroutine xmhd_set_ops
!------------------------------------------------------------------------------
!> Create new xMHD solution vector
!------------------------------------------------------------------------------
subroutine oft_xmhd_create(new,level,cache,native)
class(oft_vector), pointer, intent(out) :: new !< Vector to create
integer(i4), optional, intent(in) :: level !< FE level (optional)
logical, optional, intent(in) :: cache !< Allow caching? (optional)
logical, optional, intent(in) :: native !< Force native representation? (optional)
DEBUG_STACK_PUSH
CALL ML_xmhd_rep%vec_create(new,level=level,cache=cache,native=native)
DEBUG_STACK_POP
end subroutine oft_xmhd_create
!------------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!------------------------------------------------------------------------------
subroutine oft_xmhd_rst_save(u,t,dt,filename,path)
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL xmhd_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
  CALL hdf5_write(xmhd_rst_version,filename,'xMHD_Version')
  !---Save scales
  CALL hdf5_write(vel_scale,filename,'vel_scale')
  CALL hdf5_write(den_scale,filename,'den_scale')
  IF(n2_ind>0)CALL hdf5_write(n2_scale,filename,'n2_scale')
  IF(j2_ind>0)CALL hdf5_write(j2_scale,filename,'j2_scale')
END IF
DEBUG_STACK_POP
end subroutine oft_xmhd_rst_save
!------------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!------------------------------------------------------------------------------
subroutine oft_xmhd_rst_load(u,filename,path,t,dt,rst_version_out)
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
integer(i4), optional, intent(out) :: rst_version_out !< Version number
integer(i4) :: i,rst_version
real(r8) :: scale_tmp
real(r8), pointer, dimension(:) :: stmp
DEBUG_STACK_PUSH
CALL xmhd_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
CALL hdf5_read(rst_version,filename,'xMHD_Version')
IF(PRESENT(rst_version_out))rst_version_out=rst_version
!---Rescale if necessary scales
IF(rst_version>2)THEN
  CALL hdf5_read(scale_tmp,filename,'vel_scale')
  IF(scale_tmp/=vel_scale)THEN
    NULLIFY(stmp)
    DO i=3,5
      CALL u%get_local(stmp, i)
      stmp=stmp*scale_tmp/vel_scale
      CALL u%restore_local(stmp, i)
    END DO
    DEALLOCATE(stmp)
  END IF
  CALL hdf5_read(scale_tmp,filename,'den_scale')
  IF(scale_tmp/=den_scale)THEN
    NULLIFY(stmp)
    CALL u%get_local(stmp, 6)
    stmp=stmp*scale_tmp/den_scale
    CALL u%restore_local(stmp, 6)
    DEALLOCATE(stmp)
  END IF
  IF(n2_ind>0)THEN
    CALL hdf5_read(scale_tmp,filename,'n2_scale')
    IF(scale_tmp/=n2_scale)THEN
      NULLIFY(stmp)
      CALL u%get_local(stmp, n2_ind)
      stmp=stmp*scale_tmp/n2_scale
      CALL u%restore_local(stmp, n2_ind)
      DEALLOCATE(stmp)
    END IF
  END IF
  IF(j2_ind>0)THEN
    CALL hdf5_read(scale_tmp,filename,'j2_scale')
    IF(scale_tmp/=j2_scale)THEN
      NULLIFY(stmp)
      CALL u%get_local(stmp, j2_ind)
      stmp=stmp*scale_tmp/j2_scale
      CALL u%restore_local(stmp, j2_ind)
      DEALLOCATE(stmp)
    END IF
  END IF
END IF
DEBUG_STACK_POP
end subroutine oft_xmhd_rst_load
!------------------------------------------------------------------------------
!> Replace subfields in xMHD solution vector with specified fields
!------------------------------------------------------------------------------
subroutine oft_xmhd_push(sub_fields,u)
type(xmhd_sub_fields), intent(inout) :: sub_fields !< Source fields
class(oft_vector), intent(inout) :: u !< Destination xMHD vector
real(r8), pointer, dimension(:) :: vals
DEBUG_STACK_PUSH
NULLIFY(vals)
!---Update magnetic field
IF(ASSOCIATED(sub_fields%B))THEN
  CALL sub_fields%B%get_slice(vals,1)
  CALL u%restore_slice(vals,1,wait=.TRUE.)
  DEALLOCATE(vals)
  CALL sub_fields%B%get_slice(vals,2)
  CALL u%restore_slice(vals,2)
  DEALLOCATE(vals)
END IF
!---Update velocity field
IF(ASSOCIATED(sub_fields%V))THEN
  CALL sub_fields%V%get_slice(vals,1)
  CALL u%restore_slice(vals,3,wait=.TRUE.)
  CALL sub_fields%V%get_slice(vals,2)
  CALL u%restore_slice(vals,4,wait=.TRUE.)
  CALL sub_fields%V%get_slice(vals,3)
  CALL u%restore_slice(vals,5)
END IF
!---Update density field
IF(ASSOCIATED(sub_fields%Ne))THEN
  CALL sub_fields%Ne%get_slice(vals)
  CALL u%restore_slice(vals,6)
END IF
!---Update temperature fields
IF(ASSOCIATED(sub_fields%Ti))THEN
  CALL sub_fields%Ti%get_slice(vals)
  CALL u%restore_slice(vals,7)
END IF
IF(ASSOCIATED(sub_fields%Te))THEN
  CALL sub_fields%Te%get_slice(vals)
  CALL u%restore_slice(vals,te_ind)
END IF
IF(ASSOCIATED(sub_fields%N2))THEN
  CALL sub_fields%N2%get_slice(vals)
  CALL u%restore_slice(vals,n2_ind)
END IF
IF(ASSOCIATED(vals))DEALLOCATE(vals)
IF(ASSOCIATED(sub_fields%J2))THEN
  CALL sub_fields%J2%get_slice(vals)
  CALL u%restore_slice(vals,j2_ind)
  DEALLOCATE(vals)
END IF
DEBUG_STACK_POP
end subroutine oft_xmhd_push
!------------------------------------------------------------------------------
!> Extract subfields from xMHD solution vector
!------------------------------------------------------------------------------
subroutine oft_xmhd_pop(u,sub_fields)
class(oft_vector), intent(inout) :: u !< Source xMHD vector
type(xmhd_sub_fields), intent(inout) :: sub_fields !< Destination fields
real(r8), pointer, dimension(:) :: vals
DEBUG_STACK_PUSH
NULLIFY(vals)
!---Extract magnetic field
IF(ASSOCIATED(sub_fields%B))THEN
  CALL u%get_slice(vals,1)
  CALL sub_fields%B%restore_slice(vals,1)
  DEALLOCATE(vals)
  CALL u%get_slice(vals,2)
  CALL sub_fields%B%restore_slice(vals,2)
  DEALLOCATE(vals)
END IF
!---Extract velocity field
IF(ASSOCIATED(sub_fields%V))THEN
  CALL u%get_slice(vals,3)
  CALL sub_fields%V%restore_slice(vals,1,wait=.TRUE.)
  CALL u%get_slice(vals,4)
  CALL sub_fields%V%restore_slice(vals,2,wait=.TRUE.)
  CALL u%get_slice(vals,5)
  CALL sub_fields%V%restore_slice(vals,3)
END IF
!---Extract density field
IF(ASSOCIATED(sub_fields%Ne))THEN
  CALL u%get_slice(vals,6)
  CALL sub_fields%Ne%restore_slice(vals)
END IF
!---Extract temperature field
IF(ASSOCIATED(sub_fields%Ti))THEN
  CALL u%get_slice(vals,7)
  CALL sub_fields%Ti%restore_slice(vals)
END IF
!---Extract electron temperature field
IF(ASSOCIATED(sub_fields%Te))THEN
  CALL u%get_slice(vals,te_ind)
  CALL sub_fields%Te%restore_slice(vals)
END IF
IF(ASSOCIATED(sub_fields%N2))THEN
  CALL u%get_slice(vals,n2_ind)
  CALL sub_fields%N2%restore_slice(vals)
END IF
IF(ASSOCIATED(vals))DEALLOCATE(vals)
IF(ASSOCIATED(sub_fields%J2))THEN
  CALL u%get_slice(vals,j2_ind)
  CALL sub_fields%J2%restore_slice(vals)
  DEALLOCATE(vals)
END IF
DEBUG_STACK_POP
end subroutine oft_xmhd_pop
!------------------------------------------------------------------------------
!> Interpolate a coarse level xMHD solution to the next finest level
!!
!! @note The global \ref xmhd_level is incremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine oft_xmhd_interp(self,acors,afine)
class(ml_xmhd_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Coarse level vector to interpolate
class(oft_vector), intent(inout) :: afine !< Interpolated solution on fine level
INTEGER(i4) :: i
REAL(r8), POINTER :: acurl(:),agrad(:),tmp(:)
DEBUG_STACK_PUSH
CALL xmhd_set_level(xmhd_level+1) ! Step one level up
!---Interpolate
CALL afine%set(0.d0)
CALL oft_xmhd_ops%interp%apply(acors,afine)
!---Correct gradient subspace following geometric interpolation
IF(oft_hcurl%order==1)THEN
  !---Get local values
  NULLIFY(agrad,acurl)
  CALL afine%get_local(acurl,1)
  CALL afine%get_local(agrad,2)
  ALLOCATE(tmp(mesh%np))
  !---
  IF(bbc(1:2)=='bc')THEN
    !$omp parallel do if(mesh%np>OFT_OMP_VTHRESH)
    DO i=1,mesh%np
      tmp(i)=agrad(i)
      agrad(i)=0.d0
    END DO
  ELSE IF(bbc(1:2)=='ic')THEN
    !$omp parallel do if(mesh%np>OFT_OMP_VTHRESH)
    DO i=1,mesh%np
      IF(mesh%global%gbp(i))THEN
        tmp(i)=0.d0
      ELSE
        tmp(i)=agrad(i)
        agrad(i)=0.d0
      END IF
    END DO
  END IF
  !---Off-diagonal part
  !$omp parallel do if(mesh%ne>OFT_OMP_VTHRESH)
  DO i=1,mesh%ne
    acurl(i) = acurl(i) + &
      (tmp(mesh%le(2,i))-tmp(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
  END DO
  !---
  CALL afine%restore_local(acurl,1,wait=.TRUE.)
  CALL afine%restore_local(agrad,2)
  DEALLOCATE(acurl,agrad,tmp)
END IF
DEBUG_STACK_POP
end subroutine oft_xmhd_interp
!------------------------------------------------------------------------------
!> Inject a fine level xMHD solution to the next coarsest level
!!
!! @note The global \ref xmhd_level is decremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine oft_xmhd_inject(self,afine,acors)
class(ml_xmhd_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Fine level vector to inject
class(oft_vector), intent(inout) :: acors !< Injected solution on coarse level
real(r8), pointer, dimension(:) :: vtmp
integer(i4) :: i
DEBUG_STACK_PUSH
call xmhd_set_level(xmhd_level-1) ! Step one level down
!---Restrict
CALL oft_xmhd_ops_ML(xmhd_level+1)%interp%applyT(afine,acors)
!---Apply boundary conditions
NULLIFY(vtmp)
CALL acors%get_local(vtmp,iblock=1)
WHERE(oft_xmhd_ops%b1_bc)
  vtmp = 0.d0
END WHERE
CALL acors%restore_local(vtmp,iblock=1,wait=.TRUE.)
DEALLOCATE(vtmp)
!
CALL acors%get_local(vtmp,iblock=2)
WHERE(oft_xmhd_ops%b2_bc)
  vtmp = 0.d0
END WHERE
CALL acors%restore_local(vtmp,iblock=2,wait=.TRUE.)
DEALLOCATE(vtmp)
!
DO i=1,3
  CALL acors%get_local(vtmp,iblock=2+i)
  WHERE(oft_xmhd_ops%v_bc)
    vtmp = 0.d0
  END WHERE
  CALL acors%restore_local(vtmp,iblock=2+i,wait=.TRUE.)
END DO
!
CALL acors%get_local(vtmp,iblock=6)
WHERE(oft_xmhd_ops%n_bc)
  vtmp = 0.d0
END WHERE
CALL acors%restore_local(vtmp,iblock=6,wait=.TRUE.)
!
CALL acors%get_local(vtmp,iblock=7)
WHERE(oft_xmhd_ops%t_bc)
  vtmp = 0.d0
END WHERE
CALL acors%restore_local(vtmp,iblock=7)
!
IF(te_ind>0)THEN
  CALL acors%get_local(vtmp,iblock=te_ind)
  WHERE(oft_xmhd_ops%t_bc)
    vtmp = 0.d0
  END WHERE
  CALL acors%restore_local(vtmp,iblock=te_ind)
END IF
DEALLOCATE(vtmp)
DEBUG_STACK_POP
end subroutine oft_xmhd_inject
!------------------------------------------------------------------------------
!> Save or load driver info from restart file
!------------------------------------------------------------------------------
subroutine xmhd_driver_rst_dummy(self,rst_file)
class(oft_xmhd_driver), intent(inout) :: self !< Forcing object
character(*), intent(in) :: rst_file !< Name of restart file
IF(oft_env%head_proc)CALL oft_warn('No driver rst specified')
end subroutine xmhd_driver_rst_dummy
!------------------------------------------------------------------------------
!> Flush internal write buffers for probe object
!------------------------------------------------------------------------------
subroutine xmhd_probe_flush(self)
class(oft_xmhd_probe), intent(inout) :: self !< Probe object
end subroutine xmhd_probe_flush
!------------------------------------------------------------------------------
!> Plotting subroutine for xMHD post-processing
!!
!! Runtime options are set in the main input file using the group `xmhd_plot_options`.
!! Additionally, options from the `xmhd_options` group are also read to get run parameters,
!! see @ref xmhd::xmhd_run "xmhd_run" for a list of options.
!!
!! **Option group:** `xmhd_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `probe_only=F`     | Sample probes but do not create plots | bool |
!! | `plot_div=F`       | Plot \f$ \nabla \cdot B \f$ field | bool |
!! | `file_list="none"` | List of history files to read during post-processing (w/o file ext) | str(40) |
!! | `t0=0`             | Starting time for sampling | real |
!! | `t1=1`             | End time for sampling | real |
!! | `dt=1E-6`          | Time between sampling points | real |
!! | `rst_start=0`      | First restart file to read | int |
!! | `rst_end=2000`     | Last restart file to read | int |
!------------------------------------------------------------------------------
subroutine xmhd_plot(probes)
CLASS(oft_xmhd_probe), OPTIONAL, INTENT(inout) :: probes !< Probe object (optional)
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
class(oft_vector), pointer :: u,uc,up,ap,bp,x
TYPE(xmhd_sub_fields) :: sub_fields
TYPE(oft_hcurl_grad_cinterp), TARGET :: Jfield
TYPE(oft_hcurl_grad_rinterp), TARGET :: Bfield
TYPE(oft_hcurl_rinterp), TARGET :: J2field
type(oft_timer) :: mytimer
real(r8), pointer :: bvout(:,:)
real(r8), pointer :: vals(:),hcvals(:)
integer(i4) :: i,j,ip,u_ip,ierr,io_flag,io_unit,io_stat,rst_cur,rst_tmp,nl_update
real(r8) :: mag,div_err,mer,merp,ver,verp,gerr,cerr,verr
real(r8) :: fac,lramp,tflux,tcurr,t,tp,td
real(r8) :: ndens
!---
TYPE(xdmf_plot_file) :: xdmf
LOGICAL :: rst,first
character(LEN=XMHD_RST_LEN) :: rst_char
CHARACTER(LEN=OFT_PATH_SLEN) :: file_tmp,file_prev
!---Input variables
real(r8) :: lin_tol,nl_tol,dt_run
integer(i4) :: rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget
!---Plot group
LOGICAL :: linear=.FALSE.
LOGICAL :: probe_only=.FALSE.
LOGICAL :: plot_div=.FALSE.
real(r8) :: t0=0._r8
real(r8) :: t1=1._r8
real(r8) :: dt=1.E-6_r8
CHARACTER(LEN=OFT_PATH_SLEN) :: file_list='none'
INTEGER(i4) :: rst_start=0
INTEGER(i4) :: rst_end=2000
!---
namelist/xmhd_plot_options/probe_only,plot_div,file_list,t0,t1,dt, &
  rst_start,rst_end
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(2X,A)')'Starting xMHD post-processing'
END IF
mg_mesh=>xmhd_ML_hcurl%ml_mesh
IF(.NOT.oft_3D_hcurl_cast(oft_hcurl,xmhd_ML_hcurl%current_level))CALL oft_abort("Invalid Curl FE object","xmhd_plot",__FILE__)
IF(.NOT.oft_3D_lagrange_cast(oft_lagrange,xmhd_ML_lagrange%current_level))CALL oft_abort("Invalid Lagrange FE object","xmhd_plot",__FILE__)
IF(.NOT.oft_3D_h1_cast(oft_hgrad,xmhd_ML_h1grad%current_level))CALL oft_abort("Invalid Grad FE object","xmhd_plot",__FILE__)
mesh=>oft_hcurl%mesh
file_list="none"
open(NEWUNIT=io_unit,FILE=oft_env%ifile)
read(io_unit,xmhd_plot_options,IOSTAT=ierr)
close(io_unit)
if(ierr<0)call oft_abort('No MHD plot options found in input file.','xmhd_plot',__FILE__)
if(ierr>0)call oft_abort('Error parsing MHD plot options in input file.','xmhd_plot',__FILE__)
!---
CALL xdmf%setup("mug")
CALL mesh%setup_io(xdmf,oft_hcurl%order)
!------------------------------------------------------------------------------
! Check run type and optional fields
!------------------------------------------------------------------------------
!---Check for linearization fields
rst_cur=0
100 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
WRITE(rst_char,100)rst_cur
READ(rst_char,100,IOSTAT=io_stat)rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
file_tmp='xMHD_'//rst_char//'.rst'
IF(oft_file_exist(file_tmp))linear=hdf5_field_exist(file_tmp,'U0_Bc')
!---Find first file to be processed
IF(TRIM(file_list)=="none")THEN
  rst_cur=rst_start
  WRITE(rst_char,100)rst_cur
  READ(rst_char,100,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
  file_tmp='xMHD_'//rst_char//'.rst'
ELSE
  OPEN(NEWUNIT=io_unit,FILE=TRIM(file_list))
  READ(io_unit,*,IOSTAT=io_flag)file_tmp
  IF(io_flag>0)CALL oft_abort('Error reading restart list','xmhd_plot',__FILE__)
  CLOSE(io_unit)
END IF
!---Check for optional fields
IF(oft_file_exist(file_tmp))THEN
  xmhd_two_temp=hdf5_field_exist(file_tmp,'U_Te')
  IF(hdf5_field_exist(file_tmp,'U_n2'))d2_dens=1.d0
  IF(hdf5_field_exist(file_tmp,'U_j2'))eta_hyper=1.d0
ELSE
  CALL oft_abort("First restart file cannot be found or is inaccessible",'xmhd_plot',__FILE__)
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,'(4X,A,L1)')'Linear   = ',linear
  WRITE(*,'(4X,A,L1)')'Te found = ',xmhd_two_temp
  WRITE(*,'(4X,A,L1)')'N2 found = ',(d2_dens>0.d0)
  WRITE(*,'(4X,A,L1)')'J2 found = ',(eta_hyper>0.d0)
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
END IF
!------------------------------------------------------------------------------
! Read-in run parameters (only `rst_freq` is used)
!------------------------------------------------------------------------------
CALL xmhd_read_settings(dt_run,lin_tol,nl_tol,rst_ind,nsteps,rst_freq,nclean,maxextrap,ittarget,nl_update)
rst_ind=rst_start
nsteps=rst_end
!------------------------------------------------------------------------------
! Setup ML environment
!------------------------------------------------------------------------------
!---Build composite representation if necessary
IF(ML_xmhd_rep%nlevels==0)THEN
  CALL xmhd_setup_rep
  ALLOCATE(oft_xmhd_ops_ML(xmhd_ML_hcurl%nlevels))
END IF
!------------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------------
call oft_xmhd_create(u)
call oft_xmhd_create(uc)
call oft_xmhd_create(up)
!------------------------------------------------------------------------------
! Create plotting fields
!------------------------------------------------------------------------------
NULLIFY(hcvals)
call xmhd_ML_lagrange%vec_create(x)
call xmhd_ML_vlagrange%vec_create(ap)
call xmhd_ML_vlagrange%vec_create(bp)
allocate(bvout(3,x%n))
!---Allocate sub-fields
call xmhd_ML_hcurl_grad%vec_create(sub_fields%B)
call xmhd_ML_vlagrange%vec_create(sub_fields%V)
call xmhd_ML_lagrange%vec_create(sub_fields%Ne)
call xmhd_ML_lagrange%vec_create(sub_fields%Ti)
IF(xmhd_two_temp)call xmhd_ML_lagrange%vec_create(sub_fields%Te)
IF(n2_ind>0)call xmhd_ML_lagrange%vec_create(sub_fields%N2)
IF(j2_ind>0)call xmhd_ML_hcurl%vec_create(sub_fields%J2)
!------------------------------------------------------------------------------
! Setup Lagrange mass solver
!------------------------------------------------------------------------------
NULLIFY(lmop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,lmop,'none')
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
!---Load linearization fields if available
IF(linear)THEN
  rst_cur=0
  WRITE(rst_char,100)rst_cur
  READ(rst_char,100,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
  CALL oft_xmhd_rst_load(uc, 'xMHD_'//rst_char//'.rst', 'U0')
!------------------------------------------------------------------------------
! Retrieve sub-fields and scale
!------------------------------------------------------------------------------
  CALL oft_xmhd_pop(uc,sub_fields)
  CALL sub_fields%V%scale(vel_scale)
  CALL sub_fields%Ne%scale(den_scale)
!------------------------------------------------------------------------------
! Extract density and plot
!------------------------------------------------------------------------------
  vals=>bvout(1,:)
  CALL sub_fields%Ne%get_local(vals)
  CALL mesh%save_vertex_scalar(vals,xdmf,'N0')
!------------------------------------------------------------------------------
! Extract temperatures and plot
!------------------------------------------------------------------------------
  vals=>bvout(1,:)
  CALL sub_fields%Ti%get_local(vals)
  CALL mesh%save_vertex_scalar(vals,xdmf,'T0')
  !---Electron temperature
  IF(xmhd_two_temp)THEN
    CALL sub_fields%Te%get_local(vals)
    CALL mesh%save_vertex_scalar(vals,xdmf,'Te0')
  END IF
!------------------------------------------------------------------------------
! Extract velocity and plot
!------------------------------------------------------------------------------
  vals=>bvout(1,:)
  CALL sub_fields%V%get_local(vals,1)
  vals=>bvout(2,:)
  CALL sub_fields%V%get_local(vals,2)
  vals=>bvout(3,:)
  CALL sub_fields%V%get_local(vals,3)
  CALL mesh%save_vertex_vector(bvout,xdmf,'V0')
!------------------------------------------------------------------------------
! Project magnetic field and plot
!------------------------------------------------------------------------------
  Bfield%u=>sub_fields%B
  CALL Bfield%setup(xmhd_ML_hcurl_grad%current_level)
  CALL oft_lag_vproject(oft_lagrange,Bfield,bp)
  CALL ap%set(0.d0)
  CALL lminv%apply(ap,bp)
  vals=>bvout(1,:)
  CALL ap%get_local(vals,1)
  vals=>bvout(2,:)
  CALL ap%get_local(vals,2)
  vals=>bvout(3,:)
  CALL ap%get_local(vals,3)
  CALL mesh%save_vertex_vector(bvout,xdmf,'B0')
  !---Project current density and plot
  Jfield%u=>sub_fields%B
  CALL Jfield%setup(xmhd_ML_hcurl_grad%current_level)
  CALL oft_lag_vproject(oft_lagrange,Jfield,bp)
  CALL ap%set(0.d0)
  CALL lminv%apply(ap,bp)
  vals=>bvout(1,:)
  CALL ap%get_local(vals,1)
  vals=>bvout(2,:)
  CALL ap%get_local(vals,2)
  vals=>bvout(3,:)
  CALL ap%get_local(vals,3)
  CALL mesh%save_vertex_vector(bvout,xdmf,'J0')
END IF
!------------------------------------------------------------------------------
! Count restart file list
!------------------------------------------------------------------------------
IF(file_list(1:4)/="none")THEN
  INQUIRE(FILE=TRIM(file_list),EXIST=rst)
  IF(.NOT.rst)CALL oft_abort('Restart list file does not exist.','xmhd_plot',__FILE__)
  !---
  nsteps=0
  rst_freq=1
  OPEN(NEWUNIT=io_unit,FILE=TRIM(file_list))
  DO
    READ(io_unit,*,IOSTAT=io_flag)file_tmp
    IF(io_flag>0)THEN
      CALL oft_abort('Error reading restart list','xmhd_plot',__FILE__)
    ELSE IF(io_flag<0)THEN
      EXIT
    END IF
    nsteps=nsteps+1
  END DO
  IF(oft_env%head_proc)WRITE(*,*)'Found ',INT(nsteps,2),' restart files in list'
  REWIND(io_unit)
END IF
!------------------------------------------------------------------------------
! Load fields
!------------------------------------------------------------------------------
IF(file_list(1:4)/="none")THEN
  READ(io_unit,*,IOSTAT=io_flag)file_tmp
  IF(io_flag>0)THEN
    CALL oft_abort('Error reading restart list','xmhd_plot',__FILE__)
  END IF
  rst_cur=0
ELSE
  rst_cur=rst_ind
  WRITE(rst_char,100)rst_cur
  READ(rst_char,100,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
  file_tmp='xMHD_'//rst_char//'.rst'
END IF
rst=oft_file_exist(TRIM(file_tmp))
CALL oft_mpi_barrier(ierr)
IF(.NOT.rst)THEN
  CALL oft_abort('Restart file does not exist.','xmhd_plot',__FILE__)
ELSE
  CALL oft_xmhd_rst_load(u, file_tmp, 'U', t=t)
END IF
!---Initialize time sampling
IF(t0<t.AND.dt>0.d0)CALL oft_abort('Requested start time is not in *.rst range', &
  'xmhd_plot',__FILE__)
tp=t
td=t0
IF(td==0.d0)td=t0+dt
ip=0
u_ip=0
file_prev=file_tmp
CALL up%add(0.d0,1.d0,u)
!------------------------------------------------------------------------------
! Loop over restart files
!------------------------------------------------------------------------------
first=.TRUE.
DO
  IF(td>t1)EXIT
!------------------------------------------------------------------------------
! Load fields
!------------------------------------------------------------------------------
  IF(file_list(1:4)/="none")THEN
    READ(io_unit,*,IOSTAT=io_flag)file_tmp
    IF(io_flag>0)THEN
      CALL oft_abort('Error reading restart list','xmhd_plot',__FILE__)
    ELSE IF(io_flag<0)THEN
      EXIT
    END IF
    rst_cur=rst_cur+1
  ELSE
    IF(dt<0.d0.AND.first)rst_cur=rst_cur-rst_freq
    rst_cur=rst_cur+rst_freq
    IF(rst_cur>nsteps)EXIT
    WRITE(rst_char,100)rst_cur
    READ(rst_char,100,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
    file_tmp='xMHD_'//rst_char//'.rst'
  END IF
  rst=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst)THEN
    CALL oft_abort('Restart file does not exist.','xmhd_plot',__FILE__)
  ELSE
    CALL hdf5_read(t,file_tmp,'t')
  END IF
  IF(dt<0.d0)td=t
  first=.FALSE.
!------------------------------------------------------------------------------
! Check if next sample point is in current time interval
!------------------------------------------------------------------------------
  IF(td>=tp.AND.td<=t)THEN
    CALL oft_xmhd_rst_load(u, file_tmp, 'U')
    IF(u_ip/=ip)THEN
      CALL oft_xmhd_rst_load(up, file_prev, 'U')
    END IF
    DO WHILE(td>=tp.AND.td<=t)
      IF(dt<0.d0)THEN
        CALL uc%add(0.d0,1.d0,u)
      ELSE
        CALL uc%add(0.d0,(td-t)/(tp-t),up,(td-tp)/(t-tp),u)
      END IF
      CALL xdmf%add_timestep(td)
!------------------------------------------------------------------------------
! Retrieve sub-fields and setup interpolations
!------------------------------------------------------------------------------
      CALL oft_xmhd_pop(uc,sub_fields)
      CALL sub_fields%V%scale(vel_scale)
      CALL sub_fields%Ne%scale(den_scale)
      IF(n2_ind>0)CALL sub_fields%N2%scale(n2_scale)
      IF(j2_ind>0)CALL sub_fields%J2%scale(j2_scale)
!------------------------------------------------------------------------------
! Save probe signals
!------------------------------------------------------------------------------
      IF(PRESENT(probes))CALL probes%apply(sub_fields,td)
      IF(probe_only)THEN
        !---Update sampling time
        td=td+dt
        IF(dt<0.d0)EXIT
        CYCLE
      END IF
!------------------------------------------------------------------------------
! Extract density and plot
!------------------------------------------------------------------------------
      vals=>bvout(1,:)
      CALL sub_fields%Ne%get_local(vals)
      CALL mesh%save_vertex_scalar(vals,xdmf,'N')
      !---Hyper-diff aux field
      IF(n2_ind>0)THEN
        CALL sub_fields%N2%get_local(vals)
        CALL mesh%save_vertex_scalar(vals,xdmf,'N2')
      END IF
!------------------------------------------------------------------------------
! Extract temperatures and plot
!------------------------------------------------------------------------------
      CALL sub_fields%Ti%get_local(vals)
      CALL mesh%save_vertex_scalar(vals,xdmf,'T')
      !---Electron temperature
      IF(xmhd_two_temp)THEN
        CALL sub_fields%Te%get_local(vals)
        CALL mesh%save_vertex_scalar(vals,xdmf,'Te')
      END IF
!------------------------------------------------------------------------------
! Extract velocity and plot
!------------------------------------------------------------------------------
      vals=>bvout(1,:)
      CALL sub_fields%V%get_local(vals,1)
      vals=>bvout(2,:)
      CALL sub_fields%V%get_local(vals,2)
      vals=>bvout(3,:)
      CALL sub_fields%V%get_local(vals,3)
      CALL mesh%save_vertex_vector(bvout,xdmf,'V')
!------------------------------------------------------------------------------
! Project magnetic field and plot
!------------------------------------------------------------------------------
      Bfield%u=>sub_fields%B
      CALL Bfield%setup(xmhd_ML_hcurl_grad%current_level)
      CALL oft_lag_vproject(oft_lagrange,Bfield,bp)
      CALL ap%set(0.d0)
      CALL lminv%apply(ap,bp)
      vals=>bvout(1,:)
      CALL ap%get_local(vals,1)
      vals=>bvout(2,:)
      CALL ap%get_local(vals,2)
      vals=>bvout(3,:)
      CALL ap%get_local(vals,3)
      CALL mesh%save_vertex_vector(bvout,xdmf,'B')
      !---Hyper-res aux field
      IF(j2_ind>0)THEN
        J2field%u=>sub_fields%J2
        CALL J2field%setup(oft_hcurl)
        CALL oft_lag_vproject(oft_lagrange,J2field,bp)
        CALL ap%set(0.d0)
        CALL lminv%apply(ap,bp)
        vals=>bvout(1,:)
        CALL ap%get_local(vals,1)
        vals=>bvout(2,:)
        CALL ap%get_local(vals,2)
        vals=>bvout(3,:)
        CALL ap%get_local(vals,3)
        CALL mesh%save_vertex_vector(bvout,xdmf,'J2')
      END IF
      !---Current density
      Jfield%u=>sub_fields%B
      CALL Jfield%setup(xmhd_ML_hcurl_grad%current_level)
      CALL oft_lag_vproject(oft_lagrange,Jfield,bp)
      CALL ap%set(0.d0)
      CALL lminv%apply(ap,bp)
      vals=>bvout(1,:)
      CALL ap%get_local(vals,1)
      vals=>bvout(2,:)
      CALL ap%get_local(vals,2)
      vals=>bvout(3,:)
      CALL ap%get_local(vals,3)
      CALL mesh%save_vertex_vector(bvout,xdmf,'J')
      !---Divergence error
      IF(plot_div)THEN
        CALL oft_lag_project_div(oft_lagrange,Bfield,x)
        vals=>bvout(1,:)
        CALL x%get_local(vals)
        CALL bp%set(0.d0)
        CALL bp%restore_local(vals,1)
        CALL ap%set(0.d0)
        CALL lminv%apply(ap,bp)
        vals=>bvout(1,:)
        CALL ap%get_local(vals,1)
        CALL mesh%save_vertex_scalar(vals,xdmf,'divB')
      END IF
!------------------------------------------------------------------------------
! Update sampling time
!------------------------------------------------------------------------------
      td=td+dt
      IF(dt<0.d0)EXIT
    END DO
    u_ip=rst_cur
    CALL up%add(0.d0,1.d0,u)
  END IF
  !---
  ip=rst_cur
  file_prev=file_tmp
  tp=t
END DO
IF(file_list(1:4)/="none")CLOSE(io_unit)
!---Clear up memory
!--- solver fields
CALL u%delete
CALL uc%delete
CALL up%delete
DEALLOCATE(u,uc,up)
!--- plotting fields
CALL x%delete
CALL ap%delete
CALL bp%delete
DEALLOCATE(x,ap,bp)
CALL sub_fields%B%delete
CALL sub_fields%V%delete
CALL sub_fields%Ne%delete
CALL sub_fields%Ti%delete
DEALLOCATE(sub_fields%B,sub_fields%V,sub_fields%Ne,sub_fields%Ti)
IF(xmhd_two_temp)THEN
  CALL sub_fields%Te%delete
  DEALLOCATE(sub_fields%Te)
END IF
IF(n2_ind>0)THEN
  CALL sub_fields%N2%delete
  DEALLOCATE(sub_fields%N2)
END IF
IF(j2_ind>0)THEN
  CALL sub_fields%J2%delete
  DEALLOCATE(sub_fields%J2)
END IF
!--- interpolation objects
CALL Bfield%delete
CALL Jfield%delete
DEALLOCATE(bvout)
DEBUG_STACK_POP
end subroutine xmhd_plot
!------------------------------------------------------------------------------
!> Simple subroutine to compute timing of different solve phases
!------------------------------------------------------------------------------
subroutine xmhd_profile(u)
class(oft_vector), intent(inout) :: u
INTEGER(4) :: i,j,k
REAL(8) :: elapsed_time,matvec_time,nlop_time,build_optime,dot_time
class(oft_vector), pointer :: v
type(oft_timer) :: mytimer
CALL xmhd_rep%vec_create(v)
matvec_time=0.d0
nlop_time=0.d0
build_optime=0.d0
dot_time=0.d0
IF(oft_env%head_proc)WRITE(*,'(A)')'Starting profiling'
comm_times=0
DO i=1,10
  IF(oft_env%head_proc)CALL mytimer%tick
  CALL xmhd_set_ops(u)
  IF(oft_env%head_proc)build_optime=build_optime+mytimer%tock()
  DO j=1,10
    IF(oft_env%head_proc)CALL mytimer%tick
    CALL oft_xmhd_ops%A%apply(u, v)
    IF(oft_env%head_proc)nlop_time=nlop_time+mytimer%tock()
    !
    IF(oft_env%head_proc)CALL mytimer%tick
    DO k=1,10
      CALL oft_xmhd_ops%J%apply(u, v)
    END DO
    IF(oft_env%head_proc)matvec_time=matvec_time+mytimer%tock()
    IF(oft_env%head_proc)CALL mytimer%tick
    DO k=1,100
      elapsed_time=u%dot(v)
    END DO
    IF(oft_env%head_proc)dot_time=dot_time+mytimer%tock()
  END DO
END DO
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,Es11.3)')'  Dot    = ',dot_time/REAL(10*10*100,8)
  WRITE(*,'(A,Es11.3)')'  MatVec = ',matvec_time/REAL(10*10*10,8)
  WRITE(*,'(A,Es11.3)')'  NL Op  = ',nlop_time/REAL(10*10,8)
  WRITE(*,'(A,Es11.3)')'  Build  = ',matvec_time/REAL(10,8)
END IF
CALL oft_finalize()
end subroutine xmhd_profile
end module xmhd
