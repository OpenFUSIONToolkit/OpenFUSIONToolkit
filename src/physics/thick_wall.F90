!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file thick_wall.F90
!
!> Module for modeling eddy currents in 3D volumetric ("thick") conducting
!! structures (ThickCurr).
!!
!! The magnetic field is represented by a vector potential \f$ B = \nabla \times A \f$,
!! with \f$ A \f$ in the mixed H(Curl) + Grad(H^1) (Nedelec) space used by the
!! @ref xmhd "xMHD" module. Using the temporal gauge
!! (\f$ E = -\partial A / \partial t \f$), Ohm's law and Ampere's law give
!!
!! \f[ \sigma \frac{\partial A}{\partial t} + \nabla \times \left( \frac{1}{\mu_0}
!! \nabla \times A \right) = 0 \f]
!!
!! where \f$ \sigma = 1 / \eta \f$ and \f$ \eta \f$ is the resistivity of the region in
!! \f$ \Omega \cdot m \f$ (SI units, matching ThinCurr). Non-conducting regions are
!! treated as a pseudo-vacuum with a large resistivity (`eta_vac`), which drives the
!! field in those regions toward \f$ \nabla \times B = 0 \f$.
!!
!! The higher order Grad(H^1) subspace is required to accurately represent the
!! electric field (\f$ -\partial A / \partial t \f$) in conducting regions. In vacuum
!! regions these DOFs only change the gauge and are removed (constrained to zero) for
!! all DOFs that do not touch a conducting cell, which reduces the number of active
!! DOFs and removes most of the gradient null space in vacuum. The vertex Grad(H^1)
!! DOFs duplicate the lowest order H(Curl) gradients and are always removed.
!!
!! The outer boundary of the domain is treated as a perfect conductor
!! (\f$ \hat{n} \times E = 0 \f$), which holds the tangential component of \f$ A \f$
!! (and therefore \f$ B \cdot \hat{n} \f$) fixed in time.
!!
!! Time advance uses a linear \f$ \theta \f$-scheme for the increment
!! \f$ \delta A = A^{n+1} - A^n \f$
!! \f[ \left( M_{\sigma} + \theta \Delta t K \right) \delta A = - \Delta t K A^n \f]
!! where \f$ M_{\sigma} \f$ is the mass matrix weighted by \f$ \mu_0 / \eta \f$ and
!! \f$ K_{ij} = \int \nabla \times w_i \cdot \nabla \times w_j dV \f$. The default
!! \f$ \theta = 1 \f$ (backward Euler) is L-stable and damps the stiff pseudo-vacuum
!! modes.
!!
!! The decay times (L/R time constants) of the slowest eigenmodes may also be
!! computed with @ref thick_wall::thickcurr_eigs "thickcurr_eigs".
!!
!! Runtime options are available when using the setup (@ref thick_wall::thickcurr_setup
!! "thickcurr_setup") and post-processing (@ref thick_wall::thickcurr_plot "thickcurr_plot")
!! routines, see their documentation for available options.
!!
!! @date November 2026
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE thick_wall
#define THICKCURR_RST_LEN 5
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_stitching, ONLY: oft_seam
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr
USE oft_deriv_matrices, ONLY: oft_noop_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_mlpre, create_pre_xml, create_diag_pre, &
  create_cg_solver
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
!---
USE fem_base, ONLY: fem_max_levels, oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE fem_utils, ONLY: fem_avg_bcc
USE oft_lag_basis, ONLY: oft_scalar_fem, oft_3D_lagrange_cast
USE oft_lag_operators, ONLY: oft_lag_vgetmop, oft_lag_vproject
USE oft_hcurl_basis, ONLY: oft_hcurl_fem, oft_3D_hcurl_cast, oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops
USE oft_h1_basis, ONLY: oft_h1_fem, oft_3D_h1_cast, oft_h1_geval_all
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp, oft_hcurl_grad_cinterp, &
  oft_ml_hcurl_grad_vecspace, hcurl_grad_getmop, hcurl_grad_grad
!---
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> Forcing object class for ThickCurr
!!
!! Used to apply external sources or boundary conditions during the time advance.
!! If passed to @ref thick_wall::thickcurr_run "thickcurr_run" the driver is called
!! each time step, after the RHS has been computed and before the linear solve,
!! and may modify the RHS and/or the initial guess for the increment.
!---------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_thickcurr_driver
CONTAINS
  !> Modify vectors to apply forcing
  PROCEDURE(thickcurr_driver_apply), DEFERRED :: apply
END TYPE oft_thickcurr_driver
!---------------------------------------------------------------------------------
!> Probe object class for ThickCurr
!!
!! Used to sample signals from the solution. If passed to
!! @ref thick_wall::thickcurr_run "thickcurr_run" or
!! @ref thick_wall::thickcurr_plot "thickcurr_plot" the object is called to
!! sample the desired signals.
!---------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_thickcurr_probe
CONTAINS
  !> Sample probe signals from solution
  PROCEDURE(thickcurr_probe_apply), DEFERRED :: apply
  !> Flush internal write buffer
  PROCEDURE :: flush => thickcurr_probe_flush
END TYPE oft_thickcurr_probe
!---Interface definitions
ABSTRACT INTERFACE
!------------------------------------------------------------------------------
!> Modify solution vectors to apply forcing
!------------------------------------------------------------------------------
  SUBROUTINE thickcurr_driver_apply(self,ap,da,rhs,t,dt)
    IMPORT oft_thickcurr_driver, oft_vector, r8
    CLASS(oft_thickcurr_driver), INTENT(inout) :: self !< Forcing object
    CLASS(oft_vector), INTENT(inout) :: ap !< Vector potential at start of time step
    CLASS(oft_vector), INTENT(inout) :: da !< Initial guess for increment in vector potential
    CLASS(oft_vector), INTENT(inout) :: rhs !< RHS for linear solve
    REAL(r8), INTENT(in) :: t !< Time at start of time step
    REAL(r8), INTENT(in) :: dt !< Current timestep
  END SUBROUTINE thickcurr_driver_apply
!------------------------------------------------------------------------------
!> Sample probe signals from solution
!------------------------------------------------------------------------------
  SUBROUTINE thickcurr_probe_apply(self,a,t)
    IMPORT oft_thickcurr_probe, oft_vector, r8
    CLASS(oft_thickcurr_probe), INTENT(inout) :: self !< Probe object
    CLASS(oft_vector), INTENT(inout) :: a !< Vector potential at current time
    REAL(r8), INTENT(in) :: t !< Current time
  END SUBROUTINE thickcurr_probe_apply
END INTERFACE
!------------------------------------------------------------------------------
!> Constrained (removed) DOF flags for a given FE level
!------------------------------------------------------------------------------
TYPE :: thickcurr_bc_flags
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: curl => NULL() !< Fixed H(Curl) DOFs (boundary)
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: grad => NULL() !< Fixed Grad(H^1) DOFs (vertex, boundary, or vacuum)
  LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: grad_vac => NULL() !< Grad(H^1) DOFs removed in vacuum (interior, high order)
END TYPE thickcurr_bc_flags
!------------------------------------------------------------------------------
!> Interpolate a H(Curl) + Grad(H^1) field scaled by \f$ \alpha / \eta \f$ in conducting
!! regions and zero in vacuum regions (eg. \f$ J = -\frac{1}{\eta} \frac{\partial A}{\partial t} \f$)
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_hcurl_grad_rinterp) :: thickcurr_cond_interp
  REAL(r8) :: alpha = 1.d0 !< Scale factor
  LOGICAL :: use_eta = .TRUE. !< Divide by resistivity (and zero in vacuum)?
CONTAINS
  !> Reconstruct field
  PROCEDURE :: interp => cond_interp_apply
END TYPE thickcurr_cond_interp
!------------------------------------------------------------------------------
!> Operator for eigenvalue solves, \f$ M_{\sigma} P \f$, where \f$ P \f$ is the
!! \f$ M_{\sigma} \f$-orthogonal projection that removes the gradient null space
!------------------------------------------------------------------------------
TYPE, EXTENDS(oft_noop_matrix) :: thickcurr_eig_amat
CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => eig_amat_apply
END TYPE thickcurr_eig_amat
!------------------------------------------------------------------------------
! Global variables
!------------------------------------------------------------------------------
!---FE representation (must be set by the user before calling setup)
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: thickcurr_ML_hcurl => NULL() !< H(Curl) FE space
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: thickcurr_ML_H1grad => NULL() !< Grad(H^1) FE space
TYPE(oft_ml_fem_comp_type), POINTER, PUBLIC :: thickcurr_ML_hcurl_grad => NULL() !< Full H(Curl) + Grad(H^1) FE space
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: thickcurr_ML_lagrange => NULL() !< Lagrange FE space (optional, for plotting)
TYPE(oft_ml_fem_comp_type), POINTER, PUBLIC :: thickcurr_ML_vlagrange => NULL() !< Vector Lagrange FE space (optional, for plotting)
INTEGER(i4), PUBLIC :: thickcurr_minlev = -1 !< Lowest MG level (negative to use top level only)
!---Run parameters
REAL(r8) :: tc_dt = 1.E-4_r8 !< Time step [s]
REAL(r8) :: tc_theta = 1._r8 !< Time centering parameter
REAL(r8) :: tc_eta_vac = -1._r8 !< Pseudo-vacuum resistivity [Ohm-m]
REAL(r8) :: tc_lin_tol = 1.E-9_r8 !< Linear solver tolerance
INTEGER(i4) :: tc_nsteps = 100 !< Number of time steps
INTEGER(i4) :: tc_rst_freq = 10 !< Frequency of restart files
INTEGER(i4) :: tc_rst_ind = 0 !< Index of restart file to start from
INTEGER(i4), DIMENSION(fem_max_levels) :: tc_nu = 20 !< Number of smoother iterations for MG
!---Region information
REAL(r8), ALLOCATABLE, DIMENSION(:) :: tc_eta_reg !< Resistivity for each region [Ohm-m]
LOGICAL, ALLOCATABLE, DIMENSION(:) :: tc_cond_reg !< Is region a conductor?
!---Operators
LOGICAL :: tc_initialized = .FALSE. !< Has setup been completed?
INTEGER(i4) :: tc_level = 1 !< Top FE level
INTEGER(i4) :: tc_minlev = 1 !< Lowest FE level for MG
TYPE(thickcurr_bc_flags), DIMENSION(fem_max_levels) :: tc_bcs !< Constrained DOFs on each level
REAL(r8) :: tc_ndofs(3) = 0.d0 !< Number of DOFs (total, active, removed in vacuum)
CLASS(oft_matrix), POINTER :: tc_kop => NULL() !< Curl-curl operator (no BC)
CLASS(oft_matrix), POINTER :: tc_mcop => NULL() !< Conductivity weighted mass matrix (conductors only, no BC)
CLASS(oft_matrix), POINTER :: tc_mop => NULL() !< Conductivity weighted mass matrix (boundary DOFs removed)
TYPE(oft_matrix_ptr), DIMENSION(fem_max_levels) :: tc_lhs !< LHS operator on each level
CLASS(oft_solver), POINTER :: tc_solver => NULL() !< Linear solver for LHS
TYPE(oft_ml_hcurl_grad_vecspace), POINTER :: tc_vecspace => NULL() !< ML vector space for MG
TYPE(xml_node) :: tc_root_node !< ThickCurr XML node
TYPE(xml_node) :: tc_pre_node !< Preconditioner XML node
!---Gradient projection (eigenvalue solves)
LOGICAL :: tc_proj_initialized = .FALSE. !< Has gradient projection been setup?
LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: tc_lop_fixed => NULL() !< Fixed DOFs for weighted Laplacian
CLASS(oft_matrix), POINTER :: tc_lop => NULL() !< Conductivity weighted Laplacian
CLASS(oft_solver), POINTER :: tc_lop_solver => NULL() !< Solver for weighted Laplacian
!---Plotting
CLASS(oft_matrix), POINTER :: tc_hcurl_mop => NULL() !< Unweighted H(Curl) + Grad(H^1) mass matrix
CLASS(oft_solver), POINTER :: tc_hcurl_minv => NULL() !< Unweighted H(Curl) + Grad(H^1) mass solver
CLASS(oft_matrix), POINTER :: tc_lag_mop => NULL() !< Vector Lagrange mass matrix
CLASS(oft_solver), POINTER :: tc_lag_minv => NULL() !< Vector Lagrange mass solver
!------------------------------------------------------------------------------
! Exported interfaces
!------------------------------------------------------------------------------
PUBLIC thickcurr_setup, thickcurr_run, thickcurr_eigs, thickcurr_plot, thickcurr_save_field
PUBLIC thickcurr_diag, thickcurr_plot_order, thickcurr_rst_save, thickcurr_rst_load, thickcurr_delete
PUBLIC thickcurr_zerob, thickcurr_clean_field, thickcurr_driver_apply, thickcurr_probe_apply
CONTAINS
!------------------------------------------------------------------------------
!> Read run settings for ThickCurr from input files
!!
!! Runtime options are set in the main input file using the `thickcurr_options` group.
!!
!! **Option group:** `thickcurr_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `dt=1E-4`        | Time step [s] | real |
!! | `theta=1`        | Time centering (1 = backward Euler, 0.5 = Crank-Nicolson) | real |
!! | `eta_vac=-1`     | Pseudo-vacuum resistivity [Ohm-m] (negative for 1E4*max(eta)) | real |
!! | `lin_tol=1E-9`   | Linear solver tolerance | real |
!! | `nsteps=100`     | Number of time steps to take | int |
!! | `rst_freq=10`    | Frequency of restart files | int |
!! | `rst_ind=0`      | Index of restart file to initialize from | int |
!! | `nu_tc=20`       | Number of smoother iterations for MG preconditioner | int [nlevels] |
!!
!! Conducting regions are specified in the XML input file using `region` elements
!! within the `thickcurr` element, where `id` is a list of mesh region IDs and `eta`
!! is the resistivity in Ohm-m. Regions not specified are treated as a pseudo-vacuum.
!! An optional `pre` element may be used to specify the preconditioner used in the
!! linear solver.
!!
!! \verbatim
!! <thickcurr>
!!   <region>
!!     <id>1</id>
!!     <eta>1.7E-8</eta>
!!   </region>
!!   <pre type="block_jacobi">
!!     <nlocal>1</nlocal>
!!     <solver type="lu"></solver>
!!   </pre>
!! </thickcurr>
!! \endverbatim
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_read_settings()
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: dt,theta,eta_vac,lin_tol
INTEGER(i4) :: nsteps,rst_freq,rst_ind
INTEGER(i4), DIMENSION(fem_max_levels) :: nu_tc
NAMELIST/thickcurr_options/dt,theta,eta_vac,lin_tol,nsteps,rst_freq,rst_ind,nu_tc
DEBUG_STACK_PUSH
!---Set defaults
dt=tc_dt
theta=tc_theta
eta_vac=tc_eta_vac
lin_tol=tc_lin_tol
nsteps=tc_nsteps
rst_freq=tc_rst_freq
rst_ind=tc_rst_ind
nu_tc=tc_nu
!---Read-in parameters
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thickcurr_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error parsing ThickCurr options in input file.','thickcurr_read_settings',__FILE__)
IF(ierr<0.AND.oft_env%head_proc)CALL oft_warn('No "thickcurr_options" found in input file, using defaults')
!---Check settings
IF(dt<=0.d0)CALL oft_abort('Invalid time step, must be > 0','thickcurr_read_settings',__FILE__)
IF(theta<0.5d0.OR.theta>1.d0)CALL oft_abort('Invalid "theta", must be in [0.5,1]','thickcurr_read_settings',__FILE__)
IF(rst_freq<=0)CALL oft_abort('Invalid "rst_freq", must be > 0','thickcurr_read_settings',__FILE__)
tc_dt=dt
tc_theta=theta
tc_eta_vac=eta_vac
tc_lin_tol=lin_tol
tc_nsteps=nsteps
tc_rst_freq=rst_freq
tc_rst_ind=rst_ind
tc_nu=nu_tc
!---Look for XML nodes
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"thickcurr",tc_root_node,ierr)
  IF(ierr==0)CALL xml_get_element(tc_root_node,"pre",tc_pre_node,ierr)
END IF
DEBUG_STACK_POP
END SUBROUTINE thickcurr_read_settings
!------------------------------------------------------------------------------
!> Setup resistivity for each mesh region from XML input or a passed list
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_setup_regions(mesh,eta_reg)
CLASS(oft_mesh), INTENT(in) :: mesh !< Mesh
REAL(r8), OPTIONAL, INTENT(in) :: eta_reg(:) !< Resistivity for each region [Ohm-m] (<=0 for vacuum)
INTEGER(i4) :: i,j,ierr
INTEGER(i4), POINTER :: regs(:)
REAL(r8) :: eta
TYPE(xml_node) :: inner_node
TYPE(xml_node), POINTER :: reg_node
TYPE(xml_nodelist) :: reg_nodes
DEBUG_STACK_PUSH
IF(ALLOCATED(tc_eta_reg))DEALLOCATE(tc_eta_reg,tc_cond_reg)
ALLOCATE(tc_eta_reg(mesh%nreg),tc_cond_reg(mesh%nreg))
tc_eta_reg=-1.d0
IF(PRESENT(eta_reg))THEN
  IF(SIZE(eta_reg)/=mesh%nreg)CALL oft_abort("Size of resistivity list does not match number of regions", &
    "thickcurr_setup_regions",__FILE__)
  tc_eta_reg=eta_reg
ELSE IF(tc_root_node%associated())THEN
  CALL xml_get_element(tc_root_node,"region",reg_nodes,ierr)
  DO i=1,reg_nodes%n
    reg_node=>reg_nodes%nodes(i)
    !---Get region IDs
    CALL xml_get_element(reg_node,"id",inner_node,ierr)
    IF(ierr/=0)CALL oft_abort("Error reading region IDs for group","thickcurr_setup_regions",__FILE__)
    NULLIFY(regs)
    CALL xml_read_content(inner_node,regs,iostat=ierr)
    IF(ierr/=0)CALL oft_abort("Error reading values in id group","thickcurr_setup_regions",__FILE__)
    IF(ANY(regs>mesh%nreg).OR.ANY(regs<=0))CALL oft_abort("Invalid region ID","thickcurr_setup_regions",__FILE__)
    !---Get resistivity
    CALL xml_get_element(reg_node,"eta",inner_node,ierr)
    IF(ierr/=0)CALL oft_abort("Error reading eta for group","thickcurr_setup_regions",__FILE__)
    CALL xml_read_content(inner_node,eta,iostat=ierr)
    IF(ierr/=0)CALL oft_xml_abort("Error reading values in eta group","thickcurr_setup_regions",__FILE__)
    IF(eta<=0.d0)CALL oft_abort("Invalid eta value specified, must be > 0","thickcurr_setup_regions",__FILE__)
    !---Set values
    DO j=1,SIZE(regs)
      IF(tc_eta_reg(regs(j))>0.d0)CALL oft_abort("Region blocks overlap","thickcurr_setup_regions",__FILE__)
      tc_eta_reg(regs(j))=eta
    END DO
    DEALLOCATE(regs)
  END DO
END IF
tc_cond_reg=(tc_eta_reg>0.d0)
!---Set pseudo-vacuum resistivity
IF(tc_eta_vac<=0.d0)THEN
  IF(.NOT.ANY(tc_cond_reg))CALL oft_abort("No conducting regions specified and no vacuum resistivity set", &
    "thickcurr_setup_regions",__FILE__)
  tc_eta_vac=1.d4*MAXVAL(tc_eta_reg)
END IF
DO i=1,mesh%nreg
  IF(.NOT.tc_cond_reg(i))tc_eta_reg(i)=tc_eta_vac
END DO
DEBUG_STACK_POP
END SUBROUTINE thickcurr_setup_regions
!------------------------------------------------------------------------------
!> Setup ThickCurr model and construct operators and solvers
!!
!! The FE representation (@ref thick_wall::thickcurr_ml_hcurl "thickcurr_ML_hcurl",
!! @ref thick_wall::thickcurr_ml_h1grad "thickcurr_ML_H1grad", and
!! @ref thick_wall::thickcurr_ml_hcurl_grad "thickcurr_ML_hcurl_grad") must be
!! set before calling this subroutine. See @ref thick_wall::thickcurr_read_settings
!! "thickcurr_read_settings" for available runtime options.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_setup(eta_reg)
REAL(r8), OPTIONAL, INTENT(in) :: eta_reg(:) !< Resistivity for each region [Ohm-m] (<=0 for vacuum, overrides XML)
INTEGER(i4) :: i,level,nlevels
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels
TYPE(oft_matrix_ptr), ALLOCATABLE, DIMENSION(:) :: ml_mats
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
TYPE(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(tc_initialized)CALL thickcurr_delete
IF(.NOT.ASSOCIATED(thickcurr_ML_hcurl_grad))CALL oft_abort("FE representation not set", &
  "thickcurr_setup",__FILE__)
IF(oft_env%head_proc)CALL mytimer%tick()
!---Read settings
CALL thickcurr_read_settings
tc_level=thickcurr_ML_hcurl_grad%level
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Invalid Curl FE object","thickcurr_setup",__FILE__)
CALL thickcurr_setup_regions(hcurl_rep%mesh,eta_reg)
!---Determine MG levels
tc_minlev=thickcurr_minlev
IF(tc_minlev<0)tc_minlev=tc_level
IF(tc_minlev>tc_level.OR.tc_minlev<thickcurr_ML_hcurl_grad%minlev) &
  CALL oft_abort("Invalid minimum MG level","thickcurr_setup",__FILE__)
nlevels=tc_level-tc_minlev+1
IF(nlevels>1)THEN
  IF(.NOT.ASSOCIATED(thickcurr_ML_hcurl_grad%interp_matrices(tc_level)%m))THEN
    CALL thickcurr_ML_hcurl_grad%build_interp(tc_minlev)
    CALL thickcurr_ML_hcurl_grad%set_level(tc_level,propogate=.TRUE.)
  END IF
END IF
!------------------------------------------------------------------------------
! Build constraint flags and operators
!------------------------------------------------------------------------------
DO level=tc_level,tc_minlev,-1
  CALL thickcurr_ML_hcurl_grad%set_level(level,propogate=.TRUE.)
  CALL thickcurr_set_bcs(level)
  NULLIFY(tc_lhs(level)%m)
  CALL thickcurr_build_mat(tc_lhs(level)%m,1.d0,tc_theta*tc_dt,.TRUE.,1.d0,.TRUE.)
END DO
CALL thickcurr_ML_hcurl_grad%set_level(tc_level,propogate=.TRUE.)
CALL thickcurr_build_mat(tc_kop,0.d0,1.d0,.FALSE.,0.d0,.TRUE.)
CALL thickcurr_build_mat(tc_mcop,1.d0,0.d0,.FALSE.,0.d0,.FALSE.)
CALL thickcurr_build_mat(tc_mop,1.d0,0.d0,.TRUE.,0.d0,.TRUE.)
!------------------------------------------------------------------------------
! Setup linear solver
!------------------------------------------------------------------------------
ALLOCATE(oft_native_gmres_solver::tc_solver)
SELECT TYPE(this=>tc_solver)
  CLASS IS(oft_native_gmres_solver)
    this%nrits=20
END SELECT
tc_solver%A=>tc_lhs(tc_level)%m
tc_solver%its=500
tc_solver%atol=tc_lin_tol
tc_solver%itplot=1
!---Preconditioner
NULLIFY(tc_solver%pre)
IF(nlevels==1)THEN
  IF(tc_pre_node%associated())THEN
    CALL create_pre_xml(tc_solver%pre,tc_pre_node,.TRUE.)
  ELSE
    CALL create_diag_pre(tc_solver%pre)
  END IF
ELSE
  ALLOCATE(levels(nlevels),ml_mats(nlevels))
  DO i=1,nlevels
    levels(i)=tc_minlev+i-1
    ml_mats(i)%m=>tc_lhs(levels(i))%m
  END DO
  ALLOCATE(tc_vecspace)
  tc_vecspace%ML_FE_rep=>thickcurr_ML_hcurl_grad
  CALL create_mlpre(tc_solver%pre,ml_mats,levels,nlevels=nlevels,ml_vecspace=tc_vecspace, &
    stype=2,nu=tc_nu(tc_minlev:tc_level),xml_root=tc_pre_node)
  DEALLOCATE(levels,ml_mats)
END IF
tc_solver%pre%A=>tc_lhs(tc_level)%m
tc_initialized=.TRUE.
!------------------------------------------------------------------------------
! Print setup information
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'ThickCurr setup'
  WRITE(*,'(2X,A,ES11.3)')'dt         = ',tc_dt
  WRITE(*,'(2X,A,ES11.3)')'theta      = ',tc_theta
  WRITE(*,'(2X,A,ES11.3)')'eta_vac    = ',tc_eta_vac
  WRITE(*,'(2X,A,ES11.3)')'L-Tol      = ',tc_lin_tol
  WRITE(*,'(2X,A,I4,A,I4)')'MG levels  = ',tc_minlev,' -> ',tc_level
  WRITE(*,'(2X,A)')'Regions:'
  DO i=1,SIZE(tc_eta_reg)
    IF(tc_cond_reg(i))THEN
      WRITE(*,'(4X,I4,A,ES11.3)')i,': eta [Ohm-m] = ',tc_eta_reg(i)
    ELSE
      WRITE(*,'(4X,I4,A)')i,': vacuum'
    END IF
  END DO
  WRITE(*,'(2X,A,ES11.3)')'# A-DOFs   = ',tc_ndofs(1)
  WRITE(*,'(2X,A,ES11.3)')'  Active   = ',tc_ndofs(2)
  WRITE(*,'(2X,A,ES11.3)')'  Vac-grad = ',tc_ndofs(3)
  WRITE(*,'(2X,A,ES11.3)')'Setup time = ',mytimer%tock()
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
END IF
DEBUG_STACK_POP
END SUBROUTINE thickcurr_setup
!------------------------------------------------------------------------------
!> Destroy ThickCurr operators and solvers
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_delete()
INTEGER(i4) :: level
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(ASSOCIATED(tc_solver))THEN
  CALL tc_solver%pre%delete
  DEALLOCATE(tc_solver%pre)
  CALL tc_solver%delete
  DEALLOCATE(tc_solver)
END IF
IF(ASSOCIATED(tc_vecspace))DEALLOCATE(tc_vecspace)
DO level=tc_level,tc_minlev,-1
  IF(ASSOCIATED(tc_lhs(level)%m))THEN
    CALL tc_lhs(level)%m%delete
    DEALLOCATE(tc_lhs(level)%m)
  END IF
  IF(ASSOCIATED(tc_bcs(level)%curl))DEALLOCATE(tc_bcs(level)%curl,tc_bcs(level)%grad,tc_bcs(level)%grad_vac)
END DO
CALL tc_kop%delete
CALL tc_mcop%delete
CALL tc_mop%delete
DEALLOCATE(tc_kop,tc_mcop,tc_mop)
IF(tc_proj_initialized)THEN
  CALL tc_lop_solver%pre%delete
  DEALLOCATE(tc_lop_solver%pre)
  CALL tc_lop_solver%delete
  DEALLOCATE(tc_lop_solver)
  CALL tc_lop%delete
  DEALLOCATE(tc_lop,tc_lop_fixed)
  tc_proj_initialized=.FALSE.
END IF
tc_initialized=.FALSE.
DEBUG_STACK_POP
END SUBROUTINE thickcurr_delete
!------------------------------------------------------------------------------
!> Determine constrained (removed) DOFs on the current FE level
!!
!! - H(Curl): DOFs on the boundary (perfectly conducting wall)
!! - Grad(H^1): vertex DOFs (gauge), DOFs on the boundary, and higher order DOFs
!!   that do not touch a conducting cell (vacuum gauge)
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_set_bcs(level)
INTEGER(i4), INTENT(in) :: level !< FE level
INTEGER(i4) :: i,j
INTEGER(i4), ALLOCATABLE :: j_grad(:)
REAL(r8), POINTER, DIMENSION(:) :: flag
CLASS(oft_vector), POINTER :: vtmp
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_set_bcs",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,thickcurr_ML_hcurl_grad%current_level%fields(2)%fe)) &
  CALL oft_abort("Incorrect Grad FE type","thickcurr_set_bcs",__FILE__)
mesh=>hcurl_rep%mesh
IF(ASSOCIATED(tc_bcs(level)%curl))DEALLOCATE(tc_bcs(level)%curl,tc_bcs(level)%grad,tc_bcs(level)%grad_vac)
ALLOCATE(tc_bcs(level)%curl(hcurl_rep%ne))
ALLOCATE(tc_bcs(level)%grad(hgrad_rep%ne),tc_bcs(level)%grad_vac(hgrad_rep%ne))
tc_bcs(level)%curl=hcurl_rep%global%gbe
!---Flag Grad(H^1) DOFs touching a conducting cell (summed over processors)
NULLIFY(flag)
CALL hgrad_rep%vec_create(vtmp)
CALL vtmp%set(0.d0)
CALL vtmp%get_local(flag)
ALLOCATE(j_grad(hgrad_rep%nce))
DO i=1,mesh%nc
  IF(.NOT.tc_cond_reg(mesh%reg(i)))CYCLE
  CALL hgrad_rep%ncdofs(i,j_grad)
  flag(j_grad)=1.d0
END DO
DEALLOCATE(j_grad)
CALL vtmp%restore_local(flag,add=.TRUE.)
CALL vtmp%get_local(flag)
DO j=1,hgrad_rep%ne
  tc_bcs(level)%grad_vac(j)=(flag(j)<0.5d0).AND.(j>mesh%np).AND.(.NOT.hgrad_rep%global%gbe(j))
  tc_bcs(level)%grad(j)=(j<=mesh%np).OR.hgrad_rep%global%gbe(j).OR.tc_bcs(level)%grad_vac(j)
END DO
!---Count DOFs (top level only)
IF(level==tc_level)THEN
  tc_ndofs(1)=REAL(hcurl_rep%global%ne+hgrad_rep%global%ne,8)
  flag=0.d0
  WHERE(tc_bcs(level)%grad_vac)flag=1.d0
  CALL vtmp%restore_local(flag)
  tc_ndofs(3)=vtmp%dot(vtmp)
  flag=0.d0
  WHERE(tc_bcs(level)%grad)flag=1.d0
  CALL vtmp%restore_local(flag)
  tc_ndofs(2)=vtmp%dot(vtmp)
  DEALLOCATE(flag)
  CALL vtmp%delete
  DEALLOCATE(vtmp)
  CALL hcurl_rep%vec_create(vtmp)
  NULLIFY(flag)
  CALL vtmp%get_local(flag)
  flag=0.d0
  WHERE(tc_bcs(level)%curl)flag=1.d0
  CALL vtmp%restore_local(flag)
  tc_ndofs(2)=tc_ndofs(1)-tc_ndofs(2)-vtmp%dot(vtmp)
END IF
CALL vtmp%delete
DEALLOCATE(vtmp,flag)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_set_bcs
!------------------------------------------------------------------------------
!> Add a diagonal value for constrained DOFs in one block of a matrix
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_add_bc_diag(mat,ne,linkage,fixed,diag_val,iblock)
CLASS(oft_matrix), INTENT(inout) :: mat !< Matrix object
INTEGER(i4), INTENT(in) :: ne !< Number of DOFs in block
TYPE(oft_seam), INTENT(in) :: linkage !< Linkage information for block
LOGICAL, INTENT(in) :: fixed(:) !< Constrained DOF flags [ne]
REAL(r8), INTENT(in) :: diag_val !< Diagonal value
INTEGER(i4), OPTIONAL, INTENT(in) :: iblock !< Block index (optional)
INTEGER(i4) :: i,j_tmp(1)
REAL(r8) :: val_tmp(1,1)
DEBUG_STACK_PUSH
val_tmp=diag_val
!---DOFs not shared with other processors
DO i=1,ne
  IF(.NOT.fixed(i))CYCLE
  IF(linkage%be(i))CYCLE
  j_tmp=i
  IF(PRESENT(iblock))THEN
    CALL mat%add_values(j_tmp,j_tmp,val_tmp,1,1,iblock,iblock)
  ELSE
    CALL mat%add_values(j_tmp,j_tmp,val_tmp,1,1)
  END IF
END DO
!---Shared DOFs (owner only)
DO i=1,linkage%nbe
  j_tmp=linkage%lbe(i)
  IF(.NOT.fixed(j_tmp(1)))CYCLE
  IF(.NOT.linkage%leo(i))CYCLE
  IF(PRESENT(iblock))THEN
    CALL mat%add_values(j_tmp,j_tmp,val_tmp,1,1,iblock,iblock)
  ELSE
    CALL mat%add_values(j_tmp,j_tmp,val_tmp,1,1)
  END IF
END DO
DEBUG_STACK_POP
END SUBROUTINE thickcurr_add_bc_diag
!------------------------------------------------------------------------------
!> Construct a ThickCurr operator on the current FE level
!!
!! The matrix is given by \f$ \alpha_M M_{\sigma} + \alpha_K K \f$ where
!! \f$ M_{\sigma} \f$ is the H(Curl) + Grad(H^1) mass matrix weighted by
!! \f$ \mu_0 / \eta \f$ and \f$ K \f$ is the curl-curl operator,
!! \f$ K_{ij} = \int \nabla \times w_i \cdot \nabla \times w_j dV \f$.
!!
!! If `apply_bc` is true, the rows and columns corresponding to constrained DOFs
!! (see @ref thick_wall::thickcurr_set_bcs "thickcurr_set_bcs") are removed and a
!! value of `bc_diag` is placed on the diagonal of these rows.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_build_mat(mat,alpha_m,alpha_k,apply_bc,bc_diag,include_vac)
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat !< Matrix object
REAL(r8), INTENT(in) :: alpha_m !< Scale factor for mass matrix
REAL(r8), INTENT(in) :: alpha_k !< Scale factor for curl-curl operator
LOGICAL, INTENT(in) :: apply_bc !< Remove constrained DOFs?
REAL(r8), INTENT(in) :: bc_diag !< Diagonal value for constrained DOFs
LOGICAL, INTENT(in) :: include_vac !< Include pseudo-vacuum regions in mass matrix?
INTEGER(i4) :: i,m,jr,jc,level
INTEGER(i4), ALLOCATABLE :: j_curl(:),j_grad(:)
REAL(r8) :: vol,det,goptmp(3,4),cgop(3,6),sigma_m,elapsed_time
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop_curl,cop_curl,rop_grad
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: mop11,mop12,mop21,mop22
LOGICAL :: curved
LOGICAL, POINTER, CONTIGUOUS, DIMENSION(:) :: bc_curl,bc_grad
CLASS(oft_vector), POINTER :: vtmp
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
CLASS(oft_mesh), POINTER :: mesh
TYPE(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing ThickCurr operator'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_build_mat",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,thickcurr_ML_hcurl_grad%current_level%fields(2)%fe)) &
  CALL oft_abort("Incorrect Grad FE type","thickcurr_build_mat",__FILE__)
mesh=>hcurl_rep%mesh
level=thickcurr_ML_hcurl_grad%level
bc_curl=>tc_bcs(level)%curl
bc_grad=>tc_bcs(level)%grad
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL thickcurr_ML_hcurl_grad%current_level%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
! Operator integration
!------------------------------------------------------------------------------
!$omp parallel private(j_curl,j_grad,rop_curl,cop_curl,rop_grad,det,mop11,mop12,mop21,mop22, &
!$omp curved,goptmp,cgop,m,vol,jc,jr,sigma_m)
ALLOCATE(j_curl(hcurl_rep%nce),j_grad(hgrad_rep%nce))
ALLOCATE(rop_curl(3,hcurl_rep%nce),cop_curl(3,hcurl_rep%nce),rop_grad(3,hgrad_rep%nce))
ALLOCATE(mop11(hcurl_rep%nce,hcurl_rep%nce),mop12(hcurl_rep%nce,hgrad_rep%nce))
ALLOCATE(mop21(hgrad_rep%nce,hcurl_rep%nce),mop22(hgrad_rep%nce,hgrad_rep%nce))
!$omp do schedule(guided)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Conductivity weight for cell (mu0/eta)
  sigma_m=alpha_m*mu0/tc_eta_reg(mesh%reg(i))
  IF((.NOT.include_vac).AND.(.NOT.tc_cond_reg(mesh%reg(i))))sigma_m=0.d0
  !---Get local reconstructed operators
  mop11=0.d0; mop12=0.d0; mop22=0.d0
  DO m=1,hcurl_rep%quad%np ! Loop over quadrature points
    IF(curved.OR.m==1)THEN
      CALL mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
      CALL oft_hcurl_get_cgops(goptmp,cgop)
    END IF
    det=vol*hcurl_rep%quad%wts(m)
    CALL oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop_curl,goptmp)
    CALL oft_hcurl_ceval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),cop_curl,cgop)
    CALL oft_h1_geval_all(hgrad_rep,i,hcurl_rep%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local matrix contributions
    DO jr=1,hcurl_rep%nce
      DO jc=1,hcurl_rep%nce
        mop11(jr,jc) = mop11(jr,jc) + (sigma_m*DOT_PRODUCT(rop_curl(:,jr),rop_curl(:,jc)) &
          + alpha_k*DOT_PRODUCT(cop_curl(:,jr),cop_curl(:,jc)))*det
      END DO
      DO jc=1,hgrad_rep%nce
        mop12(jr,jc) = mop12(jr,jc) + sigma_m*DOT_PRODUCT(rop_curl(:,jr),rop_grad(:,jc))*det
      END DO
    END DO
    DO jr=1,hgrad_rep%nce
      DO jc=1,hgrad_rep%nce
        mop22(jr,jc) = mop22(jr,jc) + sigma_m*DOT_PRODUCT(rop_grad(:,jr),rop_grad(:,jc))*det
      END DO
    END DO
  END DO
  !---Get local to global DOF mapping
  CALL hcurl_rep%ncdofs(i,j_curl)
  CALL hgrad_rep%ncdofs(i,j_grad)
  !---Remove constrained DOFs
  IF(apply_bc)THEN
    DO jr=1,hcurl_rep%nce
      IF(bc_curl(j_curl(jr)))THEN
        mop11(jr,:)=0.d0
        mop11(:,jr)=0.d0
        mop12(jr,:)=0.d0
      END IF
    END DO
    DO jr=1,hgrad_rep%nce
      IF(bc_grad(j_grad(jr)))THEN
        mop12(:,jr)=0.d0
        mop22(jr,:)=0.d0
        mop22(:,jr)=0.d0
      END IF
    END DO
  END IF
  mop21=TRANSPOSE(mop12)
  !---Add local values to global matrix
  CALL mat%atomic_add_values(j_curl,j_curl,mop11,hcurl_rep%nce,hcurl_rep%nce,1,1)
  CALL mat%atomic_add_values(j_curl,j_grad,mop12,hcurl_rep%nce,hgrad_rep%nce,1,2)
  CALL mat%atomic_add_values(j_grad,j_curl,mop21,hgrad_rep%nce,hcurl_rep%nce,2,1)
  CALL mat%atomic_add_values(j_grad,j_grad,mop22,hgrad_rep%nce,hgrad_rep%nce,2,2)
END DO
DEALLOCATE(j_curl,j_grad,rop_curl,cop_curl,rop_grad)
DEALLOCATE(mop11,mop12,mop21,mop22)
!$omp end parallel
!---Set diagonal entries for constrained rows
IF(apply_bc.AND.(bc_diag/=0.d0))THEN
  CALL thickcurr_add_bc_diag(mat,hcurl_rep%ne,hcurl_rep%linkage,bc_curl,bc_diag,1)
  CALL thickcurr_add_bc_diag(mat,hgrad_rep%ne,hgrad_rep%linkage,bc_grad,bc_diag,2)
END IF
CALL thickcurr_ML_hcurl_grad%current_level%vec_create(vtmp)
CALL mat%assemble(vtmp)
CALL vtmp%delete
DEALLOCATE(vtmp)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
END SUBROUTINE thickcurr_build_mat
!------------------------------------------------------------------------------
!> Zero constrained components of a vector (boundary, gauge, and vacuum DOFs)
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_zerob(a)
CLASS(oft_vector), INTENT(inout) :: a !< Vector
REAL(r8), POINTER, DIMENSION(:) :: vals
INTEGER(i4) :: level
DEBUG_STACK_PUSH
level=thickcurr_ML_hcurl_grad%level
NULLIFY(vals)
CALL a%get_local(vals,1)
WHERE(tc_bcs(level)%curl)vals=0.d0
CALL a%restore_local(vals,1,wait=.TRUE.)
DEALLOCATE(vals)
CALL a%get_local(vals,2)
WHERE(tc_bcs(level)%grad)vals=0.d0
CALL a%restore_local(vals,2)
DEALLOCATE(vals)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_zerob
!------------------------------------------------------------------------------
!> Prepare a vector potential for use in ThickCurr
!!
!! Vertex Grad(H^1) components are moved to the equivalent lowest order H(Curl)
!! DOFs and higher order Grad(H^1) components in vacuum are removed. Neither
!! operation changes the magnetic field.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_clean_field(a)
CLASS(oft_vector), INTENT(inout) :: a !< Vector potential
INTEGER(i4) :: i,level
REAL(r8), POINTER, DIMENSION(:) :: acurl,agrad
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)CALL oft_abort("ThickCurr not initialized","thickcurr_clean_field",__FILE__)
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_clean_field",__FILE__)
mesh=>hcurl_rep%mesh
level=thickcurr_ML_hcurl_grad%level
NULLIFY(acurl,agrad)
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
DO i=1,mesh%ne
  acurl((i-1)*hcurl_rep%gstruct(2)+1)=acurl((i-1)*hcurl_rep%gstruct(2)+1) + &
    (agrad(mesh%le(2,i))-agrad(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
END DO
agrad(1:mesh%np)=0.d0
WHERE(tc_bcs(level)%grad_vac)agrad=0.d0
CALL a%restore_local(acurl,1,wait=.TRUE.)
CALL a%restore_local(agrad,2)
DEALLOCATE(acurl,agrad)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_clean_field
!------------------------------------------------------------------------------
!> Compute diagnostic quantities for a given vector potential
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_diag(a,energy,dadt,power)
CLASS(oft_vector), INTENT(inout) :: a !< Vector potential [T-m]
REAL(r8), INTENT(out) :: energy !< Magnetic energy [J]
CLASS(oft_vector), OPTIONAL, INTENT(inout) :: dadt !< Time derivative of vector potential [V/m] (optional)
REAL(r8), OPTIONAL, INTENT(out) :: power !< Ohmic power dissipated in conducting regions [W] (requires `dadt`)
CLASS(oft_vector), POINTER :: tmp
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)CALL oft_abort("ThickCurr not initialized","thickcurr_diag",__FILE__)
CALL a%new(tmp)
CALL tc_kop%apply(a,tmp)
energy=a%dot(tmp)/(2.d0*mu0)
IF(PRESENT(power))THEN
  IF(.NOT.PRESENT(dadt))CALL oft_abort("Power requires time derivative","thickcurr_diag",__FILE__)
  CALL tc_mcop%apply(dadt,tmp)
  power=dadt%dot(tmp)/mu0
END IF
CALL tmp%delete
DEALLOCATE(tmp)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_diag
!------------------------------------------------------------------------------
!> Save ThickCurr solution to a restart file
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_rst_save(a,t,dt,filename,dadt)
CLASS(oft_vector), POINTER, INTENT(inout) :: a !< Vector potential to save
REAL(r8), INTENT(in) :: t !< Current solution time
REAL(r8), INTENT(in) :: dt !< Current timestep
CHARACTER(LEN=*), INTENT(in) :: filename !< Name of restart file
CLASS(oft_vector), POINTER, OPTIONAL, INTENT(inout) :: dadt !< Time derivative of vector potential (optional)
DEBUG_STACK_PUSH
CALL thickcurr_ML_hcurl_grad%current_level%vec_save(a,filename,'A')
IF(PRESENT(dadt))CALL thickcurr_ML_hcurl_grad%current_level%vec_save(dadt,filename,'dAdt',append=.TRUE.)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
END SUBROUTINE thickcurr_rst_save
!------------------------------------------------------------------------------
!> Load ThickCurr solution from a restart file
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_rst_load(a,filename,t,dt,dadt,dadt_found)
CLASS(oft_vector), POINTER, INTENT(inout) :: a !< Vector potential
CHARACTER(LEN=*), INTENT(in) :: filename !< Name of restart file
REAL(r8), OPTIONAL, INTENT(out) :: t !< Time of loaded solution
REAL(r8), OPTIONAL, INTENT(out) :: dt !< Timestep at loaded time
CLASS(oft_vector), POINTER, OPTIONAL, INTENT(inout) :: dadt !< Time derivative of vector potential (optional)
LOGICAL, OPTIONAL, INTENT(out) :: dadt_found !< Was the time derivative found in the file? (optional)
LOGICAL :: found
DEBUG_STACK_PUSH
CALL thickcurr_ML_hcurl_grad%current_level%vec_load(a,filename,'A')
IF(PRESENT(dadt))THEN
  found=hdf5_field_exist(filename,'dAdt_c')
  IF(found)CALL thickcurr_ML_hcurl_grad%current_level%vec_load(dadt,filename,'dAdt')
  IF(PRESENT(dadt_found))dadt_found=found
END IF
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
END SUBROUTINE thickcurr_rst_load
!------------------------------------------------------------------------------
!> Get restart filename for a given index
!------------------------------------------------------------------------------
FUNCTION thickcurr_rst_name(ind) RESULT(filename)
INTEGER(i4), INTENT(in) :: ind !< Restart index
CHARACTER(LEN=OFT_PATH_SLEN) :: filename !< Restart filename
INTEGER(i4) :: io_stat,rst_tmp
CHARACTER(LEN=THICKCURR_RST_LEN) :: rst_char
100 FORMAT (I THICKCURR_RST_LEN.THICKCURR_RST_LEN)
WRITE(rst_char,100)ind
READ(rst_char,100,IOSTAT=io_stat)rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=ind))CALL oft_abort("Step count exceeds format width", &
  "thickcurr_rst_name",__FILE__)
filename='thickcurr_'//rst_char//'.rst'
END FUNCTION thickcurr_rst_name
!------------------------------------------------------------------------------
!> Advance the vector potential in time
!!
!! Runtime options are set in the main input file using the `thickcurr_options`
!! group, see @ref thick_wall::thickcurr_read_settings "thickcurr_read_settings".
!!
!! The tangential components of the vector potential on the boundary are held
!! fixed at their initial values (perfectly conducting wall). The initial field is
!! first processed with @ref thick_wall::thickcurr_clean_field "thickcurr_clean_field".
!!
!! Restart files are written to `thickcurr_XXXXX.rst` every `rst_freq` steps and
!! time history information is written to `thickcurr.hist`. Restart files also contain
!! the time derivative of the vector potential (\f$ -E \f$) over the preceding step,
!! which is used to compute the electric field and current density for plotting. For
!! the initial restart file the time derivative from the first step is used.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_run(a,driver,probes,t_final)
CLASS(oft_vector), INTENT(inout) :: a !< Initial vector potential on input, final potential on output [T-m]
CLASS(oft_thickcurr_driver), OPTIONAL, INTENT(inout) :: driver !< Forcing object (optional)
CLASS(oft_thickcurr_probe), OPTIONAL, INTENT(inout) :: probes !< Probe object (optional)
REAL(r8), OPTIONAL, INTENT(out) :: t_final !< Time at end of run [s] (optional)
INTEGER(i4) :: i,lits
REAL(r8) :: t,dt,energy,power,elapsed_time
CLASS(oft_vector), POINTER :: u,up,du,rhs
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_file
TYPE(oft_bin_file) :: hist_file
TYPE(oft_timer) :: mytimer
LOGICAL :: exists,pm_save,save_rst0
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)CALL thickcurr_setup
!------------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------------
CALL thickcurr_ML_hcurl_grad%vec_create(u)
CALL thickcurr_ML_hcurl_grad%vec_create(up)
CALL thickcurr_ML_hcurl_grad%vec_create(du)
CALL thickcurr_ML_hcurl_grad%vec_create(rhs)
dt=tc_dt
!------------------------------------------------------------------------------
! Setup history file
!------------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('thickcurr.hist', desc="History file for ThickCurr time-dependent run")
  CALL hist_file%add_field('ts',    'i4', desc="Time step index")
  CALL hist_file%add_field('lits',  'i4', desc="Linear iteration count")
  CALL hist_file%add_field('time',  'r8', desc="Simulation time [s]")
  CALL hist_file%add_field('energy','r8', desc="Magnetic energy [J]")
  CALL hist_file%add_field('power', 'r8', desc="Ohmic power in conductors (averaged over step) [W]")
  CALL hist_file%add_field('stime', 'r8', desc="Walltime [s]")
END IF
!------------------------------------------------------------------------------
! Import initial fields
!------------------------------------------------------------------------------
rst_file=thickcurr_rst_name(tc_rst_ind)
save_rst0=.FALSE.
IF((tc_rst_ind>0).AND.oft_file_exist(TRIM(rst_file)))THEN
  CALL thickcurr_rst_load(u,TRIM(rst_file),t=t)
  IF(oft_env%head_proc)THEN
    INQUIRE(FILE='thickcurr.hist',EXIST=exists)
    IF(.NOT.exists)THEN
      CALL oft_warn('No "thickcurr.hist" file present for restart, new file created.')
      CALL hist_file%write_header
    END IF
  END IF
ELSE
  IF(tc_rst_ind>0)CALL oft_abort('Restart does not exist.','thickcurr_run',__FILE__)
  CALL u%add(0.d0,1.d0,a)
  CALL thickcurr_clean_field(u)
  t=0.d0
  save_rst0=.TRUE. ! Initial restart file is saved after first step (with time derivative)
  IF(oft_env%head_proc)CALL hist_file%write_header
END IF
IF(oft_env%head_proc)CALL hist_file%open
CALL thickcurr_diag(u,energy)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting ThickCurr time advance'
  WRITE(*,'(2X,A,ES11.3)')'t0         = ',t
  WRITE(*,'(2X,A,ES11.3)')'dt         = ',dt
  WRITE(*,'(2X,A,I8)')    'nsteps     = ',tc_nsteps
  WRITE(*,'(2X,A,ES11.3)')'Mag Energy = ',energy
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
  WRITE(*,'(A6,A6,4A13)')'Step','Lits','Time','Energy','Power','Walltime'
  CALL hist_file%write(data_i4=(/tc_rst_ind,0/),data_r8=(/t,energy,0.d0,0.d0/))
END IF
IF(PRESENT(probes))CALL probes%apply(u,t)
!------------------------------------------------------------------------------
! Begin time stepping
!------------------------------------------------------------------------------
IF(oft_env%head_proc)CALL mytimer%tick()
CALL du%set(0.d0)
DO i=1,tc_nsteps
  CALL up%add(0.d0,1.d0,u)
  !---Compute RHS, -dt*K*A^n (boundary values held fixed)
  CALL tc_kop%apply(up,rhs)
  CALL rhs%scale(-dt)
  CALL thickcurr_zerob(rhs)
  !---Apply driver (initial guess is the previous increment)
  IF(PRESENT(driver))CALL driver%apply(up,du,rhs,t,dt)
  !---Solve linear system for increment
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL tc_solver%apply(du,rhs)
  oft_env%pm=pm_save
  lits=tc_solver%cits
  IF(lits<0)CALL oft_warn("Linear solver failed to converge")
  CALL u%add(1.d0,1.d0,du)
  t=t+dt
  !---Diagnostics (du temporarily holds dA/dt)
  CALL du%scale(1.d0/dt)
  CALL thickcurr_diag(u,energy,du,power)
  IF(save_rst0)THEN
    CALL thickcurr_rst_save(up,t-dt,dt,TRIM(thickcurr_rst_name(0)),dadt=du)
    save_rst0=.FALSE.
  END IF
  IF(MOD(tc_rst_ind+i,tc_rst_freq)==0)CALL thickcurr_rst_save(u,t,dt,TRIM(thickcurr_rst_name(tc_rst_ind+i)),dadt=du)
  CALL du%scale(dt)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    WRITE(*,'(I6,I6,4ES13.5)')tc_rst_ind+i,lits,t,energy,power,elapsed_time
    CALL hist_file%write(data_i4=(/tc_rst_ind+i,lits/),data_r8=(/t,energy,power,elapsed_time/))
  END IF
  IF(PRESENT(probes))CALL probes%apply(u,t)
  !---Flush output with restart files
  IF(MOD(tc_rst_ind+i,tc_rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL hist_file%flush
    IF(PRESENT(probes))CALL probes%flush
  END IF
END DO
IF(save_rst0)CALL thickcurr_rst_save(u,t,dt,TRIM(thickcurr_rst_name(0))) ! No steps taken
IF(oft_env%head_proc)CALL hist_file%close
IF(PRESENT(probes))CALL probes%flush
!---Copy final solution
CALL a%add(0.d0,1.d0,u)
IF(PRESENT(t_final))t_final=t
!---Cleanup
CALL u%delete
CALL up%delete
CALL du%delete
CALL rhs%delete
DEALLOCATE(u,up,du,rhs)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_run
!------------------------------------------------------------------------------
!> Setup projection used to remove the gradient null space in eigenvalue solves
!!
!! The null space of the curl-curl operator is spanned by gradients of Grad(H^1)
!! functions consistent with the constrained DOFs (zero on the boundary and on
!! higher order DOFs removed in vacuum). The projection
!! \f$ P = I - G L^{-1} G^T M_{\sigma} \f$, where \f$ G \f$ is the discrete gradient
!! and \f$ L = G^T M_{\sigma} G \f$ is the conductivity weighted Laplacian, removes
!! these components in an \f$ M_{\sigma} \f$-orthogonal manner.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_grad_proj_setup()
INTEGER(i4) :: i,j,m,jr,jc,level
INTEGER(i4), ALLOCATABLE :: j_grad(:)
REAL(r8) :: vol,det,goptmp(3,4),sigma_m
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gop,lop
LOGICAL :: curved
CLASS(oft_vector), POINTER :: vtmp
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(tc_proj_initialized)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_grad_proj_setup",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,thickcurr_ML_hcurl_grad%current_level%fields(2)%fe)) &
  CALL oft_abort("Incorrect Grad FE type","thickcurr_grad_proj_setup",__FILE__)
mesh=>hgrad_rep%mesh
level=thickcurr_ML_hcurl_grad%level
!---Fixed scalar DOFs: boundary and higher order DOFs removed in vacuum
ALLOCATE(tc_lop_fixed(hgrad_rep%ne))
DO j=1,hgrad_rep%ne
  tc_lop_fixed(j)=hgrad_rep%global%gbe(j).OR.tc_bcs(level)%grad_vac(j)
END DO
!---Construct conductivity weighted Laplacian (same quadrature as mass matrix)
NULLIFY(tc_lop)
CALL hgrad_rep%mat_create(tc_lop)
!$omp parallel private(j_grad,gop,lop,det,curved,goptmp,m,vol,jc,jr,sigma_m)
ALLOCATE(j_grad(hgrad_rep%nce),gop(3,hgrad_rep%nce),lop(hgrad_rep%nce,hgrad_rep%nce))
!$omp do schedule(guided)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  sigma_m=mu0/tc_eta_reg(mesh%reg(i))
  lop=0.d0
  DO m=1,hcurl_rep%quad%np
    IF(curved.OR.m==1)CALL mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    CALL oft_h1_geval_all(hgrad_rep,i,hcurl_rep%quad%pts(:,m),gop,goptmp)
    DO jr=1,hgrad_rep%nce
      DO jc=1,hgrad_rep%nce
        lop(jr,jc) = lop(jr,jc) + sigma_m*DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      END DO
    END DO
  END DO
  CALL hgrad_rep%ncdofs(i,j_grad)
  DO jr=1,hgrad_rep%nce
    IF(tc_lop_fixed(j_grad(jr)))THEN
      lop(jr,:)=0.d0
      lop(:,jr)=0.d0
    END IF
  END DO
  CALL tc_lop%atomic_add_values(j_grad,j_grad,lop,hgrad_rep%nce,hgrad_rep%nce)
END DO
DEALLOCATE(j_grad,gop,lop)
!$omp end parallel
CALL thickcurr_add_bc_diag(tc_lop,hgrad_rep%ne,hgrad_rep%linkage,tc_lop_fixed,1.d0)
CALL hgrad_rep%vec_create(vtmp)
CALL tc_lop%assemble(vtmp)
CALL vtmp%delete
DEALLOCATE(vtmp)
!---Setup solver
CALL create_cg_solver(tc_lop_solver)
tc_lop_solver%A=>tc_lop
tc_lop_solver%its=5000
CALL create_diag_pre(tc_lop_solver%pre)
tc_proj_initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE thickcurr_grad_proj_setup
!------------------------------------------------------------------------------
!> Apply the transposed discrete gradient, \f$ g = G^T y \f$
!!
!! Vertex components are computed from the lowest order H(Curl) DOFs and higher
!! order components from the Grad(H^1) DOFs. Shared DOFs are only included by
!! their owning processor before the result is summed.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_gradtp(y,g)
CLASS(oft_vector), INTENT(inout) :: y !< H(Curl) + Grad(H^1) dual vector
CLASS(oft_vector), INTENT(inout) :: g !< Grad(H^1) result
INTEGER(i4) :: i,j,k
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: emap,gmap
REAL(r8), POINTER, DIMENSION(:) :: ycurl,ygrad,gloc
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_gradtp",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,thickcurr_ML_hcurl_grad%current_level%fields(2)%fe)) &
  CALL oft_abort("Incorrect Grad FE type","thickcurr_gradtp",__FILE__)
mesh=>hgrad_rep%mesh
NULLIFY(ycurl,ygrad,gloc)
CALL y%get_local(ycurl,1)
CALL y%get_local(ygrad,2)
CALL g%set(0.d0)
CALL g%get_local(gloc)
gloc=0.d0
!---Vertex part from lowest order edge DOFs
ALLOCATE(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
!$omp parallel do private(j,k)
DO i=1,mesh%np
  DO j=mesh%kpe(i),mesh%kpe(i+1)-1
    k=mesh%lpe(j)
    IF(.NOT.mesh%fullmesh)THEN
      IF(emap(k)>0)THEN
        IF(.NOT.mesh%estitch%leo(emap(k)))CYCLE
      END IF
    END IF
    IF(mesh%global%le(k)<0)THEN
      gloc(i)=gloc(i)+ycurl((k-1)*hcurl_rep%gstruct(2)+1)*(1-2*MIN(1,ABS(mesh%le(1,k)-i)))
    ELSE
      gloc(i)=gloc(i)-ycurl((k-1)*hcurl_rep%gstruct(2)+1)*(1-2*MIN(1,ABS(mesh%le(1,k)-i)))
    END IF
  END DO
END DO
DEALLOCATE(emap)
!---Higher order part (identity)
ALLOCATE(gmap(hgrad_rep%ne))
CALL get_inverse_map(hgrad_rep%linkage%lbe,hgrad_rep%linkage%nbe,gmap,hgrad_rep%ne)
!$omp parallel do
DO i=mesh%np+1,hgrad_rep%ne
  IF(gmap(i)>0)THEN
    IF(.NOT.hgrad_rep%linkage%leo(gmap(i)))CYCLE
  END IF
  gloc(i)=gloc(i)+ygrad(i)
END DO
DEALLOCATE(gmap)
CALL g%restore_local(gloc,add=.TRUE.)
DEALLOCATE(ycurl,ygrad,gloc)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_gradtp
!------------------------------------------------------------------------------
!> Remove gradient components from a vector potential, \f$ x = P x \f$
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_grad_proj(x)
CLASS(oft_vector), INTENT(inout) :: x !< Vector potential (constrained DOFs must be zero)
REAL(r8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: tmp,g,u
LOGICAL :: pm_save
DEBUG_STACK_PUSH
CALL x%new(tmp)
CALL thickcurr_ML_H1grad%vec_create(g)
CALL thickcurr_ML_H1grad%vec_create(u)
!---Weighted weak divergence, G^T M x
CALL tc_mop%apply(x,tmp)
CALL thickcurr_gradtp(tmp,g)
NULLIFY(vals)
CALL g%get_local(vals)
WHERE(tc_lop_fixed)vals=0.d0
CALL g%restore_local(vals)
DEALLOCATE(vals)
!---Solve for gradient correction
CALL u%set(0.d0)
tc_lop_solver%atol=MAX(1.d-14,1.d-12*SQRT(g%dot(g)))
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL tc_lop_solver%apply(u,g)
oft_env%pm=pm_save
!---Remove gradient, x = x - G u
CALL u%scale(-1.d0)
CALL hcurl_grad_grad(thickcurr_ML_hcurl_grad%current_level,u,x)
CALL thickcurr_zerob(x)
CALL tmp%delete
CALL g%delete
CALL u%delete
DEALLOCATE(tmp,g,u)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_grad_proj
!------------------------------------------------------------------------------
!> Apply eigenvalue operator \f$ b = M_{\sigma} P a \f$
!------------------------------------------------------------------------------
SUBROUTINE eig_amat_apply(self,a,b)
CLASS(thickcurr_eig_amat), INTENT(inout) :: self !< Matrix object
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Source field
CLASS(oft_vector), INTENT(inout) :: b !< Result of matrix product
CLASS(oft_vector), POINTER :: tmp
DEBUG_STACK_PUSH
CALL a%new(tmp)
CALL tmp%add(0.d0,1.d0,a)
CALL thickcurr_zerob(tmp)
CALL thickcurr_grad_proj(tmp)
CALL tc_mop%apply(tmp,b)
CALL tmp%delete
DEALLOCATE(tmp)
DEBUG_STACK_POP
END SUBROUTINE eig_amat_apply
!------------------------------------------------------------------------------
!> Compute the slowest decaying eigenmodes (L/R time constants) of the system
!!
!! Eigenmodes satisfy \f$ K x = \lambda M_{\sigma} x \f$ (with constrained DOFs removed), with decay time \f$ \tau = 1/\lambda \f$. To exclude the gradient
!! null space and target the slowest modes, the equivalent problem
!! \f$ M_{\sigma} P x = \mu (M_{\sigma} + \theta \Delta t K) x \f$ is solved for the
!! largest \f$ \mu = 1 / (1 + \theta \Delta t \lambda) \f$ using ARPACK, where \f$ P \f$
!! is the \f$ M_{\sigma} \f$-orthogonal projection that removes gradients. This reuses
!! the LHS matrix and linear solver from the time-dependent advance, so the time step
!! `dt` acts as a shift and should be chosen to be comparable to the time constants
!! of interest.
!!
!! @note Fields with no resistive decay (eg. toroidal flux in a perfectly conducting
!! torus) have \f$ \mu = 1 \f$ and are returned with \f$ \tau = -1 \f$.
!!
!! If `save_plot` is true, the modes are saved to the plot group `thickcurr_eigs` in the
!! `eigs/` subdirectory (so that they are not overwritten by other plot files), with the
!! mode index used as the time.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_eigs(nev,tau,vecs,save_plot,tol)
INTEGER(i4), INTENT(in) :: nev !< Number of eigenvalues to compute
REAL(r8), INTENT(out) :: tau(nev) !< Decay times for each mode (largest first) [s]
TYPE(oft_vector_ptr), OPTIONAL, INTENT(inout) :: vecs(nev) !< Eigenvectors (optional, created if not associated)
LOGICAL, OPTIONAL, INTENT(in) :: save_plot !< Save modes to plot file? (optional)
REAL(r8), OPTIONAL, INTENT(in) :: tol !< Eigensolver tolerance (optional)
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
INTEGER(i4), ALLOCATABLE :: sort_inds(:)
REAL(r8) :: alam,mu,dt_eff,energy
REAL(r8), POINTER, DIMENSION(:) :: vslice
CLASS(oft_vector), POINTER :: u,dadt
TYPE(thickcurr_eig_amat), TARGET :: amat
TYPE(oft_irlm_eigsolver) :: arsolver
TYPE(xdmf_plot_file) :: xdmf
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
LOGICAL :: do_plot
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)CALL thickcurr_setup
do_plot=.FALSE.
IF(PRESENT(save_plot))do_plot=save_plot
dt_eff=tc_theta*tc_dt
CALL thickcurr_grad_proj_setup
!------------------------------------------------------------------------------
! Setup eigensolver
!------------------------------------------------------------------------------
arsolver%A=>amat
arsolver%M=>tc_lhs(tc_level)%m
arsolver%Minv=>tc_solver
arsolver%mode=2
arsolver%which='LA'
arsolver%nev=nev
arsolver%ncv=MAX(2*nev+1,20)
arsolver%its=1000
arsolver%tol=1.d-8
IF(PRESENT(tol))arsolver%tol=tol
CALL thickcurr_ML_hcurl_grad%vec_create(u)
CALL u%set(1.d0)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting ThickCurr eigenvalue solve'
  WRITE(*,'(2X,A,I6)')    'nev        = ',nev
  WRITE(*,'(2X,A,ES11.3)')'Shift (dt) = ',dt_eff
END IF
CALL arsolver%apply(u,alam)
IF(arsolver%info<0)CALL oft_abort("Eigenvalue solve failed","thickcurr_eigs",__FILE__)
!------------------------------------------------------------------------------
! Compute time constants (largest first)
!------------------------------------------------------------------------------
ALLOCATE(sort_inds(nev))
DO i=1,nev
  sort_inds(i)=nev-i+1 ! ARPACK returns eigenvalues in ascending order
  mu=arsolver%eig_val(1,sort_inds(i))
  IF(mu>=1.d0-1.d-10)THEN
    tau(i)=-1.d0
  ELSE IF(mu<=0.d0)THEN
    tau(i)=0.d0
  ELSE
    tau(i)=dt_eff*mu/(1.d0-mu)
  END IF
END DO
IF(oft_env%head_proc)THEN
  WRITE(*,'(2X,A)')'Decay times [s]:'
  DO i=1,nev
    WRITE(*,'(4X,I4,ES16.7)')i,tau(i)
  END DO
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
END IF
!------------------------------------------------------------------------------
! Save eigenvectors
!------------------------------------------------------------------------------
IF(PRESENT(vecs).OR.do_plot)THEN
  IF(do_plot)THEN
    IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
      CALL oft_abort("Invalid Curl FE object","thickcurr_eigs",__FILE__)
    CALL xdmf%setup("thickcurr_eigs","eigs/")
    CALL hcurl_rep%mesh%setup_io(xdmf,thickcurr_plot_order())
  END IF
  DO i=1,nev
    vslice=>arsolver%eig_vec(:,sort_inds(i))
    CALL u%restore_slice(vslice)
    CALL thickcurr_grad_proj(u)
    !---Normalize to unit magnetic energy
    CALL thickcurr_diag(u,energy)
    IF(energy>0.d0)CALL u%scale(1.d0/SQRT(energy))
    IF(PRESENT(vecs))THEN
      IF(.NOT.ASSOCIATED(vecs(i)%f))CALL thickcurr_ML_hcurl_grad%vec_create(vecs(i)%f)
      CALL vecs(i)%f%add(0.d0,1.d0,u)
    END IF
    IF(do_plot)THEN
      CALL xdmf%add_timestep(REAL(i,8))
      IF(tau(i)>0.d0)THEN
        CALL u%new(dadt)
        CALL dadt%add(0.d0,-1.d0/tau(i),u) ! Mode decays as exp(-t/tau)
        CALL thickcurr_save_field(u,xdmf,dadt=dadt)
        CALL dadt%delete
        DEALLOCATE(dadt)
      ELSE
        CALL thickcurr_save_field(u,xdmf)
      END IF
    END IF
  END DO
END IF
!---Cleanup
CALL u%delete
DEALLOCATE(u,sort_inds)
CALL arsolver%delete
DEBUG_STACK_POP
#else
tau=0.d0
CALL oft_abort("Subroutine requires ARPACK","thickcurr_eigs",__FILE__)
#endif
END SUBROUTINE thickcurr_eigs
!------------------------------------------------------------------------------
!> Reconstruct a H(Curl) + Grad(H^1) field scaled by \f$ \alpha / \eta \f$
!! (zero in vacuum regions) or by \f$ \alpha \f$ only if `use_eta` is false
!------------------------------------------------------------------------------
SUBROUTINE cond_interp_apply(self,cell,f,gop,val)
CLASS(thickcurr_cond_interp), INTENT(inout) :: self !< Interpolation object
INTEGER(i4), INTENT(in) :: cell !< Cell for interpolation
REAL(r8), INTENT(in) :: f(:) !< Position in cell in logical coordinates
REAL(r8), INTENT(in) :: gop(3,4) !< Logical gradient vectors at f
REAL(r8), INTENT(out) :: val(:) !< Reconstructed field at f
INTEGER(i4) :: reg
CALL self%oft_hcurl_grad_rinterp%interp(cell,f,gop,val)
IF(self%use_eta)THEN
  reg=self%mesh%reg(cell)
  IF(tc_cond_reg(reg))THEN
    val=val*self%alpha/tc_eta_reg(reg)
  ELSE
    val=0.d0
  END IF
ELSE
  val=val*self%alpha
END IF
END SUBROUTINE cond_interp_apply
!------------------------------------------------------------------------------
!> Get tessellation order for plotting (Lagrange order if available, otherwise 1)
!------------------------------------------------------------------------------
FUNCTION thickcurr_plot_order() RESULT(order)
INTEGER(i4) :: order !< Tessellation order for plot files
CLASS(oft_scalar_fem), POINTER :: lag_rep
order=1
IF(ASSOCIATED(thickcurr_ML_vlagrange))THEN
  IF(.NOT.oft_3D_lagrange_cast(lag_rep,thickcurr_ML_lagrange%current_level)) &
    CALL oft_abort("Invalid Lagrange FE object","thickcurr_plot_order",__FILE__)
  order=lag_rep%order
END IF
END FUNCTION thickcurr_plot_order
!------------------------------------------------------------------------------
!> Save vector potential, magnetic field, and current density to a plot file
!!
!! The vector potential (`A`) and magnetic field (`B`), \f$ B = \nabla \times A \f$, are
!! always saved. If the time derivative of the vector potential is provided, the electric
!! field (`E`), \f$ E = -\partial A / \partial t \f$, and current density (`J`),
!! \f$ J = E / \eta \f$ in conducting regions and zero in vacuum, are also saved.
!! Otherwise the current density is computed from the weak curl of \f$ B \f$,
!! \f$ \mu_0 J = M^{-1} K A \f$, where \f$ M \f$ is the (unweighted) H(Curl) + Grad(H^1)
!! mass matrix. Note that this form includes surface currents on the perfectly conducting
!! walls, which appear in the layer of cells adjacent to the boundary.
!!
!! If the Lagrange FE representation (@ref thick_wall::thickcurr_ml_lagrange
!! "thickcurr_ML_lagrange" and @ref thick_wall::thickcurr_ml_vlagrange
!! "thickcurr_ML_vlagrange") is set, fields are projected to vertices.
!! Otherwise cell averaged values are saved. The plot file must be setup with a
!! tessellation order given by @ref thick_wall::thickcurr_plot_order "thickcurr_plot_order".
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_save_field(a,xdmf,suffix,dadt)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Vector potential [T-m]
TYPE(xdmf_plot_file), INTENT(in) :: xdmf !< Plot file object
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: suffix !< Suffix for field names (optional)
CLASS(oft_vector), TARGET, OPTIONAL, INTENT(inout) :: dadt !< Time derivative of vector potential [V/m] (optional)
CHARACTER(LEN=:), ALLOCATABLE :: tag
CLASS(oft_vector), POINTER :: jvec,tmp
TYPE(oft_hcurl_grad_rinterp), TARGET :: Afield,Jfield
TYPE(oft_hcurl_grad_cinterp), TARGET :: Bfield
TYPE(thickcurr_cond_interp), TARGET :: Efield,Jcfield
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
LOGICAL :: pm_save
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Invalid Curl FE object","thickcurr_save_field",__FILE__)
tag=''
IF(PRESENT(suffix))tag=suffix
!---Save potential and magnetic field
Afield%u=>a
CALL Afield%setup(thickcurr_ML_hcurl_grad%current_level)
Bfield%u=>a
CALL Bfield%setup(thickcurr_ML_hcurl_grad%current_level)
CALL save_interp(Afield,'A'//tag)
CALL save_interp(Bfield,'B'//tag)
CALL Afield%delete
CALL Bfield%delete
!---Save electric field and current density from time derivative
IF(PRESENT(dadt))THEN
  Efield%u=>dadt
  CALL Efield%setup(thickcurr_ML_hcurl_grad%current_level)
  Efield%mesh=>hcurl_rep%mesh
  Efield%alpha=-1.d0
  Efield%use_eta=.FALSE.
  CALL save_interp(Efield,'E'//tag)
  Jcfield%u=>dadt
  CALL Jcfield%setup(thickcurr_ML_hcurl_grad%current_level)
  Jcfield%mesh=>hcurl_rep%mesh
  Jcfield%alpha=-1.d0
  CALL save_interp(Jcfield,'J'//tag)
  CALL Efield%delete
  CALL Jcfield%delete
  DEBUG_STACK_POP
  RETURN
END IF
!---Compute current density from weak curl of B
IF(.NOT.ASSOCIATED(tc_hcurl_minv))THEN
  NULLIFY(tc_hcurl_mop)
  CALL hcurl_grad_getmop(thickcurr_ML_hcurl_grad%current_level,tc_hcurl_mop,'none')
  CALL create_cg_solver(tc_hcurl_minv)
  tc_hcurl_minv%A=>tc_hcurl_mop
  tc_hcurl_minv%its=-2
  CALL create_diag_pre(tc_hcurl_minv%pre)
END IF
CALL a%new(jvec)
CALL a%new(tmp)
CALL tc_kop%apply(a,tmp)
CALL tmp%scale(1.d0/mu0)
CALL jvec%set(0.d0)
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL tc_hcurl_minv%apply(jvec,tmp)
oft_env%pm=pm_save
Jfield%u=>jvec
CALL Jfield%setup(thickcurr_ML_hcurl_grad%current_level)
CALL save_interp(Jfield,'J'//tag)
CALL Jfield%delete
CALL jvec%delete
CALL tmp%delete
DEALLOCATE(jvec,tmp)
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Save an interpolated vector field at vertices or cell centers
!------------------------------------------------------------------------------
SUBROUTINE save_interp(field,name)
CLASS(oft_hcurl_grad_rinterp), INTENT(inout) :: field !< Field interpolator
CHARACTER(LEN=*), INTENT(in) :: name !< Name of field in plot file
INTEGER(i4) :: i
REAL(r8), POINTER, DIMENSION(:,:) :: bvout
REAL(r8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: ap,bp
CLASS(oft_scalar_fem), POINTER :: lag_rep
IF(ASSOCIATED(thickcurr_ML_vlagrange))THEN
  IF(.NOT.oft_3D_lagrange_cast(lag_rep,thickcurr_ML_lagrange%current_level)) &
    CALL oft_abort("Invalid Lagrange FE object","thickcurr_save_field",__FILE__)
  IF(.NOT.ASSOCIATED(tc_lag_minv))THEN
    NULLIFY(tc_lag_mop)
    CALL oft_lag_vgetmop(thickcurr_ML_vlagrange%current_level,tc_lag_mop,'none')
    CALL create_cg_solver(tc_lag_minv)
    tc_lag_minv%A=>tc_lag_mop
    tc_lag_minv%its=-2
    CALL create_diag_pre(tc_lag_minv%pre)
  END IF
  CALL thickcurr_ML_vlagrange%vec_create(ap)
  CALL thickcurr_ML_vlagrange%vec_create(bp)
  ALLOCATE(bvout(3,lag_rep%ne))
  CALL oft_lag_vproject(lag_rep,field,bp)
  CALL ap%set(0.d0)
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL tc_lag_minv%apply(ap,bp)
  oft_env%pm=pm_save
  DO i=1,3
    vals=>bvout(i,:)
    CALL ap%get_local(vals,i)
  END DO
  CALL hcurl_rep%mesh%save_vertex_vector(bvout,xdmf,name)
  CALL ap%delete
  CALL bp%delete
  DEALLOCATE(ap,bp,bvout)
ELSE
  ALLOCATE(bvout(3,hcurl_rep%mesh%nc))
  CALL fem_avg_bcc(hcurl_rep%mesh,field,bvout,2*hcurl_rep%order,n=3)
  CALL hcurl_rep%mesh%save_cell_vector(bvout,xdmf,name)
  DEALLOCATE(bvout)
END IF
END SUBROUTINE save_interp
END SUBROUTINE thickcurr_save_field
!------------------------------------------------------------------------------
!> Plot ThickCurr solution from restart files
!!
!! Runtime options are set in the main input file using the group `thickcurr_plot_options`.
!! Additionally, options from the `thickcurr_options` group are also read to get run
!! parameters (only `rst_freq` is used).
!!
!! **Option group:** `thickcurr_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `rst_start=0`      | First restart file to read | int |
!! | `rst_end=2000`     | Last restart file to read | int |
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_plot(probes)
CLASS(oft_thickcurr_probe), OPTIONAL, INTENT(inout) :: probes !< Probe object (optional)
INTEGER(i4) :: i,io_unit,ierr
INTEGER(i4) :: rst_start=0
INTEGER(i4) :: rst_end=2000
REAL(r8) :: t
REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_cell
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_file
CLASS(oft_vector), POINTER :: u,dadt
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
TYPE(xdmf_plot_file) :: xdmf
LOGICAL :: dadt_found
NAMELIST/thickcurr_plot_options/rst_start,rst_end
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(thickcurr_ML_hcurl_grad))CALL oft_abort("FE representation not set", &
  "thickcurr_plot",__FILE__)
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thickcurr_plot_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error parsing ThickCurr plot options in input file.','thickcurr_plot',__FILE__)
IF(.NOT.tc_initialized)CALL thickcurr_setup
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl_grad%current_level%fields(1)%fe)) &
  CALL oft_abort("Invalid Curl FE object","thickcurr_plot",__FILE__)
mesh=>hcurl_rep%mesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(2X,A)')'Starting ThickCurr post-processing'
END IF
CALL xdmf%setup("thickcurr")
CALL mesh%setup_io(xdmf,thickcurr_plot_order())
!---Save resistivity (only available for untessellated output)
IF(thickcurr_plot_order()==1)THEN
  ALLOCATE(eta_cell(mesh%nc))
  DO i=1,mesh%nc
    eta_cell(i)=tc_eta_reg(mesh%reg(i))
    IF(.NOT.tc_cond_reg(mesh%reg(i)))eta_cell(i)=-1.d0
  END DO
  CALL mesh%save_cell_scalar(eta_cell,xdmf,'eta')
  DEALLOCATE(eta_cell)
END IF
!---Loop over restart files
CALL thickcurr_ML_hcurl_grad%vec_create(u)
CALL thickcurr_ML_hcurl_grad%vec_create(dadt)
DO i=rst_start,rst_end,tc_rst_freq
  rst_file=thickcurr_rst_name(i)
  IF(.NOT.oft_file_exist(TRIM(rst_file)))EXIT
  CALL thickcurr_rst_load(u,TRIM(rst_file),t=t,dadt=dadt,dadt_found=dadt_found)
  IF(oft_env%head_proc)WRITE(*,'(4X,2A,ES13.5)')TRIM(rst_file),': t = ',t
  IF(PRESENT(probes))CALL probes%apply(u,t)
  CALL xdmf%add_timestep(t)
  IF(dadt_found)THEN
    CALL thickcurr_save_field(u,xdmf,dadt=dadt)
  ELSE
    CALL thickcurr_save_field(u,xdmf)
  END IF
END DO
CALL dadt%delete
DEALLOCATE(dadt)
IF(PRESENT(probes))CALL probes%flush
IF(oft_env%head_proc)WRITE(*,'(A)')'============================'
CALL u%delete
DEALLOCATE(u)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_plot
!------------------------------------------------------------------------------
!> Flush internal write buffers for probe object
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_probe_flush(self)
CLASS(oft_thickcurr_probe), INTENT(inout) :: self !< Probe object
! Do nothing
END SUBROUTINE thickcurr_probe_flush
END MODULE thick_wall
