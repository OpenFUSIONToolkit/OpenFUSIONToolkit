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
!! with \f$ A \f$ in the H(Curl) (Nedelec) space used by the @ref xmhd "xMHD" module,
!! but without the higher order Grad(H^1) subspace. Using the temporal gauge
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
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
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
USE fem_base, ONLY: fem_max_levels, oft_ml_fem_type, oft_ml_fe_vecspace
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE fem_utils, ONLY: fem_avg_bcc
USE oft_lag_basis, ONLY: oft_lag_setup, oft_scalar_fem, oft_3D_lagrange_cast, &
  oft_lag_geval_all
USE oft_lag_operators, ONLY: oft_lag_vgetmop, oft_lag_vproject, oft_lag_zerob
USE oft_hcurl_basis, ONLY: oft_hcurl_fem, oft_3D_hcurl_cast, oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops
USE oft_hcurl_operators, ONLY: oft_hcurl_rinterp, oft_hcurl_cinterp, oft_hcurl_zerob, &
  oft_hcurl_getmop, hcurl_grad, hcurl_gradtp
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
CLASS(oft_matrix), POINTER :: tc_kop => NULL() !< Curl-curl operator (no BC)
CLASS(oft_matrix), POINTER :: tc_mcop => NULL() !< Conductivity weighted mass matrix (conductors only, no BC)
CLASS(oft_matrix), POINTER :: tc_mop => NULL() !< Conductivity weighted mass matrix (boundary DOFs removed)
TYPE(oft_matrix_ptr), DIMENSION(fem_max_levels) :: tc_lhs !< LHS operator on each level
CLASS(oft_solver), POINTER :: tc_solver => NULL() !< Linear solver for LHS
TYPE(oft_ml_fe_vecspace), POINTER :: tc_vecspace => NULL() !< ML vector space for MG
TYPE(oft_hcurl_zerob), TARGET :: tc_zerob !< Boundary condition for H(Curl) vectors
TYPE(xml_node) :: tc_root_node !< ThickCurr XML node
TYPE(xml_node) :: tc_pre_node !< Preconditioner XML node
!---Gradient projection (eigenvalue solves)
LOGICAL :: tc_proj_initialized = .FALSE. !< Has gradient projection been setup?
TYPE(oft_ml_fem_type), TARGET :: tc_ML_lag1 !< Linear Lagrange FE space for gradient projection
TYPE(oft_lag_zerob), TARGET :: tc_lag_zerob !< Boundary condition for linear Lagrange vectors
CLASS(oft_matrix), POINTER :: tc_lop => NULL() !< Conductivity weighted Laplacian
CLASS(oft_solver), POINTER :: tc_lop_solver => NULL() !< Solver for weighted Laplacian
!---Plotting
CLASS(oft_matrix), POINTER :: tc_hcurl_mop => NULL() !< Unweighted H(Curl) mass matrix
CLASS(oft_solver), POINTER :: tc_hcurl_minv => NULL() !< Unweighted H(Curl) mass solver
CLASS(oft_matrix), POINTER :: tc_lag_mop => NULL() !< Vector Lagrange mass matrix
CLASS(oft_solver), POINTER :: tc_lag_minv => NULL() !< Vector Lagrange mass solver
!------------------------------------------------------------------------------
! Exported interfaces
!------------------------------------------------------------------------------
PUBLIC thickcurr_setup, thickcurr_run, thickcurr_eigs, thickcurr_plot, thickcurr_save_field
PUBLIC thickcurr_diag, thickcurr_plot_order, thickcurr_rst_save, thickcurr_rst_load, thickcurr_delete
PUBLIC thickcurr_zerob, thickcurr_driver_apply, thickcurr_probe_apply
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
!! The FE representation (@ref thick_wall::thickcurr_ml_hcurl "thickcurr_ML_hcurl")
!! must be set before calling this subroutine. See @ref thick_wall::thickcurr_read_settings
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
IF(.NOT.ASSOCIATED(thickcurr_ML_hcurl))CALL oft_abort("FE representation not set", &
  "thickcurr_setup",__FILE__)
IF(oft_env%head_proc)CALL mytimer%tick()
!---Read settings
CALL thickcurr_read_settings
tc_level=thickcurr_ML_hcurl%level
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl%current_level)) &
  CALL oft_abort("Invalid Curl FE object","thickcurr_setup",__FILE__)
CALL thickcurr_setup_regions(hcurl_rep%mesh,eta_reg)
tc_zerob%ML_hcurl_rep=>thickcurr_ML_hcurl
!---Determine MG levels
tc_minlev=thickcurr_minlev
IF(tc_minlev<0)tc_minlev=tc_level
IF(tc_minlev>tc_level.OR.tc_minlev<thickcurr_ML_hcurl%minlev) &
  CALL oft_abort("Invalid minimum MG level","thickcurr_setup",__FILE__)
nlevels=tc_level-tc_minlev+1
!------------------------------------------------------------------------------
! Build operators
!------------------------------------------------------------------------------
CALL thickcurr_build_mat(tc_kop,0.d0,1.d0,.FALSE.,0.d0,.TRUE.)
CALL thickcurr_build_mat(tc_mcop,1.d0,0.d0,.FALSE.,0.d0,.FALSE.)
CALL thickcurr_build_mat(tc_mop,1.d0,0.d0,.TRUE.,0.d0,.TRUE.)
DO level=tc_level,tc_minlev,-1
  CALL thickcurr_ML_hcurl%set_level(level)
  NULLIFY(tc_lhs(level)%m)
  CALL thickcurr_build_mat(tc_lhs(level)%m,1.d0,tc_theta*tc_dt,.TRUE.,1.d0,.TRUE.)
END DO
CALL thickcurr_ML_hcurl%set_level(tc_level)
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
  tc_vecspace%ML_FE_rep=>thickcurr_ML_hcurl
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
  WRITE(*,'(2X,A,ES11.3)')'# A-DOFs   = ',REAL(hcurl_rep%global%ne,8)
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
  DEALLOCATE(tc_lop)
  tc_proj_initialized=.FALSE.
END IF
tc_initialized=.FALSE.
DEBUG_STACK_POP
END SUBROUTINE thickcurr_delete
!------------------------------------------------------------------------------
!> Construct a ThickCurr operator on the current FE level
!!
!! The matrix is given by \f$ \alpha_M M_{\sigma} + \alpha_K K \f$ where
!! \f$ M_{\sigma} \f$ is the H(Curl) mass matrix weighted by \f$ \mu_0 / \eta \f$
!! and \f$ K \f$ is the curl-curl operator,
!! \f$ K_{ij} = \int \nabla \times w_i \cdot \nabla \times w_j dV \f$.
!!
!! If `apply_bc` is true, the rows and columns corresponding to boundary DOFs are
!! removed and a value of `bc_diag` is placed on the diagonal of these rows.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_build_mat(mat,alpha_m,alpha_k,apply_bc,bc_diag,include_vac)
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat !< Matrix object
REAL(r8), INTENT(in) :: alpha_m !< Scale factor for mass matrix
REAL(r8), INTENT(in) :: alpha_k !< Scale factor for curl-curl operator
LOGICAL, INTENT(in) :: apply_bc !< Remove boundary DOFs?
REAL(r8), INTENT(in) :: bc_diag !< Diagonal value for boundary DOFs
LOGICAL, INTENT(in) :: include_vac !< Include pseudo-vacuum regions in mass matrix?
INTEGER(i4) :: i,m,jr,jc
INTEGER(i4), ALLOCATABLE :: j_curl(:)
REAL(r8) :: vol,det,goptmp(3,4),cgop(3,6),sigma_m,elapsed_time
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop,cop,mop
LOGICAL :: curved
LOGICAL, ALLOCATABLE, DIMENSION(:) :: bc_flag
CLASS(oft_vector), POINTER :: vtmp
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_mesh), POINTER :: mesh
TYPE(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing ThickCurr operator'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl%current_level)) &
  CALL oft_abort("Incorrect Curl FE type","thickcurr_build_mat",__FILE__)
mesh=>hcurl_rep%mesh
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
! Operator integration
!------------------------------------------------------------------------------
!$omp parallel private(j_curl,rop,cop,det,mop,curved,goptmp,cgop,m,vol,jc,jr,sigma_m,bc_flag)
ALLOCATE(j_curl(hcurl_rep%nce),bc_flag(hcurl_rep%nce))
ALLOCATE(rop(3,hcurl_rep%nce),cop(3,hcurl_rep%nce))
ALLOCATE(mop(hcurl_rep%nce,hcurl_rep%nce))
!$omp do schedule(guided)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Conductivity weight for cell (mu0/eta)
  sigma_m=alpha_m*mu0/tc_eta_reg(mesh%reg(i))
  IF((.NOT.include_vac).AND.(.NOT.tc_cond_reg(mesh%reg(i))))sigma_m=0.d0
  !---Get local reconstructed operators
  mop=0.d0
  DO m=1,hcurl_rep%quad%np ! Loop over quadrature points
    IF(curved.OR.m==1)THEN
      CALL mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
      CALL oft_hcurl_get_cgops(goptmp,cgop)
    END IF
    det=vol*hcurl_rep%quad%wts(m)
    CALL oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop,goptmp)
    CALL oft_hcurl_ceval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    DO jr=1,hcurl_rep%nce
      DO jc=1,hcurl_rep%nce
        mop(jr,jc) = mop(jr,jc) + (sigma_m*DOT_PRODUCT(rop(:,jr),rop(:,jc)) &
          + alpha_k*DOT_PRODUCT(cop(:,jr),cop(:,jc)))*det
      END DO
    END DO
  END DO
  !---Get local to global DOF mapping
  CALL hcurl_rep%ncdofs(i,j_curl)
  !---Remove boundary DOFs
  IF(apply_bc)THEN
    DO jr=1,hcurl_rep%nce
      bc_flag(jr)=hcurl_rep%global%gbe(j_curl(jr))
    END DO
    DO jr=1,hcurl_rep%nce
      IF(bc_flag(jr))THEN
        mop(jr,:)=0.d0
        mop(:,jr)=0.d0
      END IF
    END DO
  END IF
  !---Add local values to global matrix
  CALL mat%atomic_add_values(j_curl,j_curl,mop,hcurl_rep%nce,hcurl_rep%nce)
END DO
DEALLOCATE(j_curl,bc_flag,rop,cop,mop)
!$omp end parallel
!---Set diagonal entries for boundary rows
IF(apply_bc.AND.(bc_diag/=0.d0))THEN
  ALLOCATE(mop(1,1),j_curl(1))
  mop(1,1)=bc_diag
  DO i=1,hcurl_rep%nbe
    jr=hcurl_rep%lbe(i)
    IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
    IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
    j_curl=jr
    CALL mat%add_values(j_curl,j_curl,mop,1,1)
  END DO
  DEALLOCATE(j_curl,mop)
END IF
CALL hcurl_rep%vec_create(vtmp)
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
!> Zero boundary (tangential) components of a vector potential
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_zerob(a)
CLASS(oft_vector), INTENT(inout) :: a !< Vector
DEBUG_STACK_PUSH
tc_zerob%ML_hcurl_rep=>thickcurr_ML_hcurl
CALL tc_zerob%apply(a)
DEBUG_STACK_POP
END SUBROUTINE thickcurr_zerob
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
SUBROUTINE thickcurr_rst_save(a,t,dt,filename)
CLASS(oft_vector), POINTER, INTENT(inout) :: a !< Vector potential to save
REAL(r8), INTENT(in) :: t !< Current solution time
REAL(r8), INTENT(in) :: dt !< Current timestep
CHARACTER(LEN=*), INTENT(in) :: filename !< Name of restart file
DEBUG_STACK_PUSH
CALL thickcurr_ML_hcurl%current_level%vec_save(a,filename,'A')
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
END SUBROUTINE thickcurr_rst_save
!------------------------------------------------------------------------------
!> Load ThickCurr solution from a restart file
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_rst_load(a,filename,t,dt)
CLASS(oft_vector), POINTER, INTENT(inout) :: a !< Vector potential
CHARACTER(LEN=*), INTENT(in) :: filename !< Name of restart file
REAL(r8), OPTIONAL, INTENT(out) :: t !< Time of loaded solution
REAL(r8), OPTIONAL, INTENT(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL thickcurr_ML_hcurl%current_level%vec_load(a,filename,'A')
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
!! fixed at their initial values (perfectly conducting wall).
!!
!! Restart files are written to `thickcurr_XXXXX.rst` every `rst_freq` steps and
!! time history information is written to `thickcurr.hist`.
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
LOGICAL :: exists,pm_save
DEBUG_STACK_PUSH
IF(.NOT.tc_initialized)CALL thickcurr_setup
!------------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------------
CALL thickcurr_ML_hcurl%vec_create(u)
CALL thickcurr_ML_hcurl%vec_create(up)
CALL thickcurr_ML_hcurl%vec_create(du)
CALL thickcurr_ML_hcurl%vec_create(rhs)
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
  t=0.d0
  CALL thickcurr_rst_save(u,t,dt,TRIM(thickcurr_rst_name(0)))
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
  !---Diagnostics
  CALL du%scale(1.d0/dt)
  CALL thickcurr_diag(u,energy,du,power)
  CALL du%scale(dt)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    WRITE(*,'(I6,I6,4ES13.5)')tc_rst_ind+i,lits,t,energy,power,elapsed_time
    CALL hist_file%write(data_i4=(/tc_rst_ind+i,lits/),data_r8=(/t,energy,power,elapsed_time/))
  END IF
  IF(PRESENT(probes))CALL probes%apply(u,t)
  !---Save restart file
  IF(MOD(tc_rst_ind+i,tc_rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL hist_file%flush
    IF(PRESENT(probes))CALL probes%flush
    CALL thickcurr_rst_save(u,t,dt,TRIM(thickcurr_rst_name(tc_rst_ind+i)))
  END IF
END DO
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
!! The null space of the curl-curl operator is spanned by gradients of linear
!! Lagrange functions (that vanish on the boundary). The projection
!! \f$ P = I - G L^{-1} G^T M_{\sigma} \f$, where \f$ G \f$ is the discrete gradient
!! and \f$ L = G^T M_{\sigma} G \f$ is the conductivity weighted Laplacian, removes
!! these components in an \f$ M_{\sigma} \f$-orthogonal manner.
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_grad_proj_setup()
INTEGER(i4) :: i,m,jr,jc
INTEGER(i4), ALLOCATABLE :: j_lag(:)
REAL(r8) :: vol,det,goptmp(3,4),sigma_m
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gop,lop
LOGICAL :: curved
CLASS(oft_vector), POINTER :: vtmp
CLASS(oft_scalar_fem), POINTER :: lag_rep
CLASS(multigrid_mesh), POINTER :: mg_mesh
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(tc_proj_initialized)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---Create linear Lagrange space on the finest mesh
IF(tc_ML_lag1%nlevels==0)THEN
  mg_mesh=>thickcurr_ML_hcurl%ml_mesh
  CALL oft_lag_setup(mg_mesh,1,tc_ML_lag1,minlev=mg_mesh%mgdim)
END IF
tc_lag_zerob%ML_lag_rep=>tc_ML_lag1
IF(.NOT.oft_3D_lagrange_cast(lag_rep,tc_ML_lag1%current_level)) &
  CALL oft_abort("Invalid Lagrange FE object","thickcurr_grad_proj_setup",__FILE__)
mesh=>lag_rep%mesh
!---Construct conductivity weighted Laplacian (zero on boundary)
NULLIFY(tc_lop)
CALL lag_rep%mat_create(tc_lop)
!$omp parallel private(j_lag,gop,lop,det,curved,goptmp,m,vol,jc,jr,sigma_m)
ALLOCATE(j_lag(lag_rep%nce),gop(3,lag_rep%nce),lop(lag_rep%nce,lag_rep%nce))
!$omp do schedule(guided)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  sigma_m=mu0/tc_eta_reg(mesh%reg(i))
  lop=0.d0
  DO m=1,lag_rep%quad%np
    IF(curved.OR.m==1)CALL mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    CALL oft_lag_geval_all(lag_rep,i,lag_rep%quad%pts(:,m),gop,goptmp)
    DO jr=1,lag_rep%nce
      DO jc=1,lag_rep%nce
        lop(jr,jc) = lop(jr,jc) + sigma_m*DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      END DO
    END DO
  END DO
  CALL lag_rep%ncdofs(i,j_lag)
  DO jr=1,lag_rep%nce
    IF(lag_rep%global%gbe(j_lag(jr)))THEN
      lop(jr,:)=0.d0
      lop(:,jr)=0.d0
    END IF
  END DO
  CALL tc_lop%atomic_add_values(j_lag,j_lag,lop,lag_rep%nce,lag_rep%nce)
END DO
DEALLOCATE(j_lag,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j_lag(1))
lop(1,1)=1.d0
DO i=1,lag_rep%nbe
  jr=lag_rep%lbe(i)
  IF(.NOT.lag_rep%global%gbe(jr))CYCLE
  IF(.NOT.lag_rep%linkage%leo(i))CYCLE
  j_lag=jr
  CALL tc_lop%add_values(j_lag,j_lag,lop,1,1)
END DO
DEALLOCATE(j_lag,lop)
CALL lag_rep%vec_create(vtmp)
CALL tc_lop%assemble(vtmp)
CALL vtmp%delete
DEALLOCATE(vtmp)
!---Setup solver
CALL create_cg_solver(tc_lop_solver)
tc_lop_solver%A=>tc_lop
tc_lop_solver%its=2000
CALL create_diag_pre(tc_lop_solver%pre)
tc_proj_initialized=.TRUE.
DEBUG_STACK_POP
END SUBROUTINE thickcurr_grad_proj_setup
!------------------------------------------------------------------------------
!> Remove gradient components from a vector potential, \f$ x = P x \f$
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_grad_proj(x)
CLASS(oft_vector), INTENT(inout) :: x !< Vector potential (boundary DOFs must be zero)
CLASS(oft_vector), POINTER :: tmp,g,u
LOGICAL :: pm_save
DEBUG_STACK_PUSH
CALL x%new(tmp)
CALL tc_ML_lag1%vec_create(g)
CALL tc_ML_lag1%vec_create(u)
!---Weighted weak divergence, G^T M x
CALL tc_mop%apply(x,tmp)
CALL g%set(0.d0)
CALL hcurl_gradtp(thickcurr_ML_hcurl%current_level,tmp,g)
CALL tc_lag_zerob%apply(g)
!---Solve for gradient correction
CALL u%set(0.d0)
tc_lop_solver%atol=MAX(1.d-14,1.d-12*SQRT(g%dot(g)))
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL tc_lop_solver%apply(u,g)
oft_env%pm=pm_save
!---Remove gradient, x = x - G u
CALL u%scale(-1.d0)
CALL hcurl_grad(thickcurr_ML_hcurl%current_level,u,x)
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
!! Eigenmodes satisfy \f$ K x = \lambda M_{\sigma} x \f$ (with \f$ \hat{n} \times x = 0 \f$
!! on the boundary), with decay time \f$ \tau = 1/\lambda \f$. To exclude the gradient
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
CLASS(oft_vector), POINTER :: u
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
CALL thickcurr_ML_hcurl%vec_create(u)
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
    IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl%current_level)) &
      CALL oft_abort("Invalid Curl FE object","thickcurr_eigs",__FILE__)
    CALL xdmf%setup("thickcurr_eigs")
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
      IF(.NOT.ASSOCIATED(vecs(i)%f))CALL thickcurr_ML_hcurl%vec_create(vecs(i)%f)
      CALL vecs(i)%f%add(0.d0,1.d0,u)
    END IF
    IF(do_plot)THEN
      CALL xdmf%add_timestep(REAL(i,8))
      CALL thickcurr_save_field(u,xdmf)
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
!! The magnetic field is \f$ B = \nabla \times A \f$ and the current density is
!! computed from the weak curl of \f$ B \f$, \f$ \mu_0 J = M^{-1} K A \f$, where
!! \f$ M \f$ is the (unweighted) H(Curl) mass matrix.
!!
!! If the Lagrange FE representation (@ref thick_wall::thickcurr_ml_lagrange
!! "thickcurr_ML_lagrange" and @ref thick_wall::thickcurr_ml_vlagrange
!! "thickcurr_ML_vlagrange") is set, fields are projected to vertices.
!! Otherwise cell averaged values are saved. The plot file must be setup with a
!! tessellation order given by @ref thick_wall::thickcurr_plot_order "thickcurr_plot_order".
!------------------------------------------------------------------------------
SUBROUTINE thickcurr_save_field(a,xdmf,suffix)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Vector potential [T-m]
TYPE(xdmf_plot_file), INTENT(in) :: xdmf !< Plot file object
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: suffix !< Suffix for field names (optional)
CHARACTER(LEN=:), ALLOCATABLE :: tag
CLASS(oft_vector), POINTER :: jvec,tmp
TYPE(oft_hcurl_rinterp), TARGET :: Afield,Jfield
TYPE(oft_hcurl_cinterp), TARGET :: Bfield
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
LOGICAL :: pm_save
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl%current_level)) &
  CALL oft_abort("Invalid Curl FE object","thickcurr_save_field",__FILE__)
tag=''
IF(PRESENT(suffix))tag=suffix
!---Compute current density from weak curl of B
IF(.NOT.ASSOCIATED(tc_hcurl_minv))THEN
  NULLIFY(tc_hcurl_mop)
  CALL oft_hcurl_getmop(hcurl_rep,tc_hcurl_mop,'none')
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
!---Save fields
Afield%u=>a
CALL Afield%setup(hcurl_rep)
Bfield%u=>a
CALL Bfield%setup(hcurl_rep)
Jfield%u=>jvec
CALL Jfield%setup(hcurl_rep)
CALL save_interp(Afield,'A'//tag)
CALL save_interp(Bfield,'B'//tag)
CALL save_interp(Jfield,'J'//tag)
CALL Afield%delete
CALL Bfield%delete
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
CLASS(oft_hcurl_rinterp), INTENT(inout) :: field !< Field interpolator
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
CLASS(oft_vector), POINTER :: u
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
TYPE(xdmf_plot_file) :: xdmf
NAMELIST/thickcurr_plot_options/rst_start,rst_end
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(thickcurr_ML_hcurl))CALL oft_abort("FE representation not set", &
  "thickcurr_plot",__FILE__)
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thickcurr_plot_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error parsing ThickCurr plot options in input file.','thickcurr_plot',__FILE__)
IF(.NOT.tc_initialized)CALL thickcurr_setup
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,thickcurr_ML_hcurl%current_level)) &
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
CALL thickcurr_ML_hcurl%vec_create(u)
DO i=rst_start,rst_end,tc_rst_freq
  rst_file=thickcurr_rst_name(i)
  IF(.NOT.oft_file_exist(TRIM(rst_file)))EXIT
  CALL thickcurr_rst_load(u,TRIM(rst_file),t=t)
  IF(oft_env%head_proc)WRITE(*,'(4X,2A,ES13.5)')TRIM(rst_file),': t = ',t
  IF(PRESENT(probes))CALL probes%apply(u,t)
  CALL xdmf%add_timestep(t)
  CALL thickcurr_save_field(u,xdmf)
END DO
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
